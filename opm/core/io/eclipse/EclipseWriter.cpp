/*
  Copyright (c) 2013-2014 Andreas Lauser
  Copyright (c) 2013 SINTEF ICT, Applied Mathematics.
  Copyright (c) 2013 Uni Research AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "config.h"

#include "EclipseWriter.hpp"

#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/props/phaseUsageFromDeck.hpp>
#include <opm/core/simulator/SimulatorState.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/parameters/Parameter.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/wells.h> // WellType

#include <opm/parser/eclipse/Deck/DeckKeyword.hpp>
#include <opm/parser/eclipse/Utility/SpecgridWrapper.hpp>
#include <opm/parser/eclipse/Utility/WelspecsWrapper.hpp>

#include <boost/algorithm/string/case_conv.hpp> // to_upper_copy
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp> // path

#include <ctime>      // mktime
#include <forward_list>
#include <memory>     // unique_ptr
#include <utility>    // move

using namespace Opm;
using namespace Opm::parameter;

#ifdef HAVE_ERT
#include <ert/ecl/fortio.h>
#include <ert/ecl/ecl_endian_flip.h>
#include <ert/ecl/ecl_grid.h>
#include <ert/ecl/ecl_kw_magic.h>
#include <ert/ecl/ecl_kw.h>
#include <ert/ecl/ecl_sum.h>
#include <ert/ecl/ecl_util.h>
#include <ert/ecl/ecl_init_file.h>
#include <ert/ecl/ecl_file.h>
#include <ert/ecl/ecl_rst_file.h>

// namespace start here since we don't want the ERT headers in it
namespace {

namespace {
/// Helper function when we don't really want any transformation
/// (The C++ committee removed std::identity because it was "troublesome" (!?!)
static double noConversion (const double& u) { return u; }


/// Helper method that can be used in keyword transformation (must carry
/// the barsa argument)
static double toBar (const double& pressure) {
    return Opm::unit::convert::to (pressure, Opm::unit::barsa);
}

/// Helper method that can be used in keyword transformation (must carry
/// the milliDarcy argument)
static double toMilliDarcy (const double& permeability) {
    return Opm::unit::convert::to (permeability, Opm::prefix::milli * Opm::unit::darcy);
}

/// Names of the saturation property for each phase. The order of these
/// names are critical; they must be the same as the BlackoilPhases enum
static const char* SAT_NAMES[] = { "SWAT", "SOIL", "SGAS" };

} // anonymous namespace

/// Smart pointer/handle class for ERT opaque types, such as ecl_kw_type*.
///
/// \tparam T Type of handle being wrapper
template <typename T>
struct EclipseHandle {
    /// Instead of inheriting std::unique_ptr and letting the compiler
    /// provide a default implementation which calls the base class, we
    /// define the move constructor and assignment operator ourselves
    /// and call and aggregated pointer, because of bugs in GCC 4.4
    EclipseHandle <T> (EclipseHandle <T>&& rhs)
        : h_ (std::move (rhs.h_)) { }

    EclipseHandle <T>& operator= (EclipseHandle <T>&& rhs) {
        h_ = std::move (rhs.h_);
        return *this;
    }

    /// Prevent GCC 4.4 from the urge of generating a copy constructor
    EclipseHandle (const EclipseHandle&) = delete;
    EclipseHandle <T>& operator= (const EclipseHandle <T>&) = delete;

    /// Construct a new smart handle based on the returned value of
    /// an allocation function and the corresponding destroyer function.
    EclipseHandle <T> (T* t, void (*destroy)(T*))
        : h_ (t, destroy) { }

    /// Construct an object whose memory is freed as part of another
    /// structure. This constructor is dangerous! Make sure that you
    /// have the lifetime management correct before using it.
    EclipseHandle <T> (T* t) : h_ (t, no_delete) { }

    /// Convenience operator that lets us use this type as if
    /// it was a handle directly.
    operator T* () const { return h_.get (); }

private:
    std::unique_ptr <T, void (*)(T*) throw()> h_; // handle

    // helper function to pass to the second pointer constructor, since
    // the runtime library does not like this construct
    static void no_delete (T*) { }
};

// retrieve all data fields in SI units of a deck keyword
std::vector<double> getAllSiDoubles_(Opm::DeckKeywordConstPtr keywordPtr)
{
    std::vector<double> retBuff;
    for (unsigned i = 0; i < keywordPtr->size(); ++i) {
        Opm::DeckRecordConstPtr recordPtr(keywordPtr->getRecord(i));
        for (unsigned j = 0; j < recordPtr->size(); ++j) {
            Opm::DeckItemConstPtr itemPtr(recordPtr->getItem(j));
            for (unsigned k = 0; k < itemPtr->size(); ++k) {
                retBuff.push_back(itemPtr->getSIDouble(k));
            }
        }
    }
    return retBuff;
}

// retrieve all integer data fields of a deck keyword
std::vector<int> getAllIntegers_(Opm::DeckKeywordConstPtr keywordPtr)
{
    std::vector<int> retBuff;
    for (unsigned i = 0; i < keywordPtr->size(); ++i) {
        Opm::DeckRecordConstPtr recordPtr(keywordPtr->getRecord(i));
        for (unsigned j = 0; j < recordPtr->size(); ++j) {
            Opm::DeckItemConstPtr itemPtr(recordPtr->getItem(j));
            for (unsigned k = 0; k < itemPtr->size(); ++k) {
                retBuff.push_back(itemPtr->getInt(k));
            }
        }
    }
    return retBuff;
}

// throw away the data for all non-active cells in an array
void restrictToActiveCells_(std::vector<double> &data, const std::vector<int> &actnumData)
{
    assert(actnumData.size() == data.size());

    size_t curActiveIdx = 0;
    for (size_t curIdx = 0; curIdx < data.size(); ++curIdx) {
        if (!actnumData[curIdx])
            continue; // ignore non-active cells

        assert(curActiveIdx <= curIdx);
        data[curActiveIdx] = data[curIdx];
        ++ curActiveIdx;
    }

    data.resize(curActiveIdx);
}

// throw away the data for all non-active cells in an array. (this is
// the variant of the function which takes an UnstructuredGrid object.)
void restrictToActiveCells_(std::vector<double> &data, const UnstructuredGrid &grid)
{
    if (!grid.global_cell)
        // if there is no active -> global mapping, all cells
        // are considered active
        return;

    // activate those cells that are actually there
    for (int i = 0; i < grid.number_of_cells; ++i) {
        // make sure that global cell indices are always at least as
        // large as the active one and that the global cell indices
        // are in increasing order. the latter might become
        // problematic if cells are extensively re-ordered, but that
        // does not seem to be the case so far
        assert(grid.global_cell[i] >= i);
        assert(i == 0 || grid.global_cell[i - 1] < grid.global_cell[i]);

        data[i] = data[grid.global_cell[i]];
    }
    data.resize(grid.number_of_cells);
}

// convert the units of an array
template <class TransferFunction>
void convertUnit_(std::vector<double> &data, TransferFunction &transferFn)
{
    for (size_t curIdx = 0; curIdx < data.size(); ++curIdx) {
        data[curIdx] = transferFn(data[curIdx]);
    }
}

// extract a sub-array of a larger one which represents multiple
// striped ones
void extractFromStripedData_(std::vector<double> &data,
                             int offset,
                             int stride)
{
    size_t tmpIdx = 0;
    for (size_t curIdx = offset; curIdx < data.size(); curIdx += stride) {
        assert(tmpIdx <= curIdx);
        data[tmpIdx] = data[curIdx];
        ++tmpIdx;
    }
    // shirk the result
    data.resize(tmpIdx);
}

// enclosure of the current grid in a Cartesian space
int getCartesianSize_(const UnstructuredGrid& grid) {
    const int nx = grid.cartdims[0];
    const int ny = grid.cartdims[1];
    const int nz = grid.cartdims[2];
    return nx * ny * nz;
}

void getActiveCells_(const UnstructuredGrid& grid,
                     std::vector <int>& actnum)
{
    // we must fill the Cartesian grid with flags
    const int size = getCartesianSize_(grid);

    // if we don't have a global_cells field, then assume that all
    // grid cells is active
    if (!grid.global_cell) {
        if (grid.number_of_cells != size) {
            OPM_THROW (std::runtime_error,
                       "No ACTNUM map but grid size != Cartesian size");
        }
        actnum.assign (size, 1);
    }
    else {
        // start out with entire map being inactive
        actnum.assign (size, 0);

        // activate those cells that are actually there
        for (int i = 0; i < grid.number_of_cells; ++i) {
            actnum[grid.global_cell[i]] = 1;
        }
    }
}

/**
 * Eclipse "keyword" (i.e. named data) for a vector. (This class is
 * different from EclKW in the constructors it provide).
 */
template <typename T>
struct EclipseKeyword : public EclipseHandle <ecl_kw_type> {
    /// Special initialization from double-precision array.
    EclipseKeyword (const std::string& name,
                    const std::vector<double>& data)
        : EclipseHandle<ecl_kw_type>(ecl_kw_alloc(name.c_str(),
                                                  data.size(),
                                                  type()),
                                     ecl_kw_free)
    { copyData (data, &noConversion, /*offset=*/0, /*stride=*/1); }

    /// Initialization from integer array.
    EclipseKeyword (const std::string& name,
                    const std::vector<int>& data)
        : EclipseHandle<ecl_kw_type>(ecl_kw_alloc(name.c_str(),
                                                  data.size(),
                                                  type()),
                                     ecl_kw_free)
    { copyData (data, &noConversion, /*offset=*/0, /*stride=*/1); }

    /// Constructor for optional fields
    EclipseKeyword (const std::string& name)
        : EclipseHandle <ecl_kw_type> (0, ecl_kw_free) {
        static_cast<void> (name);
    }

    // constructor to keep compatibility with the old eclipse parser
    EclipseKeyword (const std::string& name,
                    const EclipseGridParser& parser);

    // GCC 4.4 doesn't generate these constructors for us; provide the
    // default implementation explicitly here instead
    EclipseKeyword (EclipseKeyword&& rhs)
        : EclipseHandle <ecl_kw_type> (std::move (rhs))
    { }

    EclipseKeyword& operator= (EclipseKeyword&& rhs)
    {
        EclipseHandle <ecl_kw_type>::operator= (std::move(rhs));
        return *this;
    }
    EclipseKeyword (const EclipseKeyword&) = delete;
    EclipseKeyword& operator= (const EclipseKeyword&) = delete;

private:
    /// Map the C++ data type (given by T) to an Eclipse type enum
    static ecl_type_enum type ();

    /// Helper function that is the meat of the constructor
    template <typename U>
    void copyData (const std::vector <U>& data,
                   double (* const transf)(const double&),
                   const int offset,
                   const int stride) {
        // number of elements to take
        const int num = dataSize (data, offset, stride);

        // fill it with values
        T* target = static_cast <T*> (ecl_kw_get_ptr (*this));
        for (int i = 0; i < num; ++i) {
            target[i] = static_cast <T> (transf (data[i * stride + offset]));
        }
    }

    // Compute the number of outputs this dataset will give
    template <typename U>
    int dataSize (const std::vector <U>& data,
                  const int offset,
                  const int stride) {
        // number of elements we can possibly take from the vector
        const int num = data.size ();

        // range cannot start outside of data set
        assert(offset >= 0 && offset < num);

        // don't jump out of the set when trying to
        assert(stride > 0 && stride < num - offset);

        // number of (strided) entries it will provide. the last item
        // in the array is num - 1. the last strided item we can pick
        // (from recs number of records) is (recs - 1) * stride + offset,
        // which must be <= num - 1. we are interested in the maximum
        // case where it holds to equals. rearranging the above gives us:
        const int recs = (num - 1 - offset) / stride + 1;
        return recs;
    }

    int dataSize (Opm::DeckKeywordConstPtr keyword,
                  const int offset = 0,
                  const int stride = 1)
    {
        int numFlatItems = 0;
        Opm::DeckRecordConstPtr record = keyword->getRecord(0);
        for (unsigned itemIdx = 0; itemIdx < record->size(); ++itemIdx) {
            numFlatItems += record->getItem(itemIdx)->size();
        }

        // range cannot start outside of data set
        assert(offset >= 0 && offset < numFlatItems);

        // don't jump out of the set when trying to
        assert(stride > 0 && stride < numFlatItems - offset);

        // number of (strided) entries it will provide. the last item
        // in the array is num - 1. the last strided item we can pick
        // (from recs number of records) is (recs - 1) * stride + offset,
        // which must be <= num - 1. we are interested in the maximum
        // case where it holds to equals. rearranging the above gives us:
        const int recs = (numFlatItems - 1 - offset) / stride + 1;
        return recs;
    }

};

// specializations for known keyword types
template <> ecl_type_enum EclipseKeyword<int   >::type () { return ECL_INT_TYPE   ; }
template <> ecl_type_enum EclipseKeyword<float >::type () { return ECL_FLOAT_TYPE ; }
template <> ecl_type_enum EclipseKeyword<double>::type () { return ECL_DOUBLE_TYPE; }

/// keywords in ERT requires single-precision type, but OPM have them
/// stored as double-precision. this template specialization instantiates
/// a copy function that downcast the data to the required type.
template <>
EclipseKeyword <float>::EclipseKeyword (
        const std::string& name,
        const EclipseGridParser& parser)
    // allocate handle and put in smart pointer base class
    : EclipseHandle <ecl_kw_type> (
          ecl_kw_alloc (name.c_str(),
                        // we can safely use the *size* of the original
                        dataSize (parser.getValue <double> (name), 0, 1),
                        type ()),
          ecl_kw_free) {
    const std::vector <double>& data = parser.getValue <double> (name);
    copyData (data, &noConversion, 0, 1);
}

/**
 * Pointer to memory that holds the name to an Eclipse output file.
 */
struct EclipseFileName : public EclipseHandle <const char> {
    EclipseFileName (const std::string& outputDir,
                     const std::string& baseName,
                     ecl_file_enum type,
                     int outputStepIdx)

        // filename formatting function returns a pointer to allocated
        // memory that must be released with the free() function
        : EclipseHandle <const char> (
              ecl_util_alloc_filename (outputDir.c_str(),
                                       baseName.c_str(),
                                       type,
                                       false, // formatted?
                                       outputStepIdx),
              freestr) { }
private:
    /// Facade which allows us to free a const char*
    static void freestr (const char* ptr) {
        ::free (const_cast<char*>(ptr));
    }
};

/// Get dimensions of the grid from the parse of the input file
std::vector <int> parserDim (const EclipseGridParser& parser) {
    std::vector<int> dim(/* n = */ 3);
    // dimensions explicitly given
    if (parser.hasField("SPECGRID")) {
        dim = parser.getSPECGRID ().dimensions;
    }
    // dimensions implicitly given by number of deltas
    else if (parser.hasField("DXV")) {
        assert(parser.hasField("DYV"));
        assert(parser.hasField("DZV"));
        dim[0] = parser.getFloatingPointValue("DXV").size();
        dim[1] = parser.getFloatingPointValue("DYV").size();
        dim[2] = parser.getFloatingPointValue("DZV").size();
    }
    else {
        OPM_THROW(std::runtime_error,
                  "Only decks featureing either the SPECGRID or the D[XYZ]V keywords "
                  "are currently supported");
    }
    return dim;
}

/// Get dimensions of the grid from the parse of the input file
std::vector <int> parserDim (Opm::DeckConstPtr newParserDeck) {
    std::vector<int> dim(/* n = */ 3);
    // dimensions explicitly given
    if (newParserDeck->hasKeyword("SPECGRID")) {
        SpecgridWrapper specgrid(newParserDeck->getKeyword("SPECGRID"));
        dim = specgrid.numBlocksVector();
    }
    // dimensions implicitly given by number of deltas
    else if (newParserDeck->hasKeyword("DXV")) {
        assert(newParserDeck->hasKeyword("DYV"));
        assert(newParserDeck->hasKeyword("DZV"));
        dim[0] = newParserDeck->getKeyword("DXV")->getRawDoubleData().size();
        dim[1] = newParserDeck->getKeyword("DYV")->getRawDoubleData().size();
        dim[2] = newParserDeck->getKeyword("DZV")->getRawDoubleData().size();
    }
    else {
        OPM_THROW(std::runtime_error,
                  "Only decks featureing either the SPECGRID or the D[XYZ]V keywords "
                  "are currently supported");
    }
    return dim;
}

/// Convert OPM phase usage to ERT bitmask
static int phaseMask (const PhaseUsage uses) {
    return (uses.phase_used [BlackoilPhases::Liquid] ? ECL_OIL_PHASE   : 0)
         | (uses.phase_used [BlackoilPhases::Aqua]   ? ECL_WATER_PHASE : 0)
         | (uses.phase_used [BlackoilPhases::Vapour] ? ECL_GAS_PHASE   : 0);
}

struct EclipseRestart : public EclipseHandle <ecl_rst_file_type> {
    EclipseRestart (const std::string& outputDir,
                    const std::string& baseName,
                    const SimulatorTimer& timer,
                    int outputStepIdx)
        // notice the poor man's polymorphism of the allocation function
        : EclipseHandle <ecl_rst_file_type> (
              (timer.currentStepNum () > 0 ? ecl_rst_file_open_append
                                           : ecl_rst_file_open_write)(
                  EclipseFileName (outputDir,
                                   baseName,
                                   ECL_UNIFIED_RESTART_FILE,
                                   outputStepIdx)),
              ecl_rst_file_close) { }

    void writeHeader (const SimulatorTimer& timer,
                      int outputStepIdx,
                      const PhaseUsage uses,
                      const EclipseGridParser parser,
                      const int num_active_cells) {
        const std::vector <int> dim = parserDim (parser);
        ecl_rst_file_fwrite_header (*this,
                                    outputStepIdx,
                                    timer.currentPosixTime(),
                                    Opm::unit::convert::to (timer.simulationTimeElapsed (),
                                                            Opm::unit::day),
                                    dim[0],
                                    dim[1],
                                    dim[2],
                                    num_active_cells,
                                    phaseMask (uses));
    }

    void writeHeader (const SimulatorTimer& timer,
                      int outputStepIdx,
                      const PhaseUsage uses,
                      Opm::DeckConstPtr newParserDeck,
                      const int num_active_cells) {
        const std::vector <int> dim = parserDim (newParserDeck);
        ecl_rst_file_fwrite_header (*this,
                                    outputStepIdx,
                                    timer.currentPosixTime(),
                                    Opm::unit::convert::to (timer.simulationTimeElapsed (),
                                                            Opm::unit::day),
                                    dim[0],
                                    dim[1],
                                    dim[2],
                                    num_active_cells,
                                    phaseMask (uses));
    }
};

/**
 * The EclipseSolution class wraps the actions that must be done to the
 * restart file while writing solution variables; it is not a handle on
 * its own.
 */
struct EclipseSolution : public EclipseHandle <ecl_rst_file_type> {
    EclipseSolution (EclipseRestart& rst_file)
        : EclipseHandle <ecl_rst_file_type> (start_solution (rst_file),
                                             ecl_rst_file_end_solution) { }

    template <typename T>
    void add (const EclipseKeyword<T>& kw) {
        ecl_rst_file_add_kw (*this, kw);
    }

private:
    /// Helper method to call function *and* return the handle
    static ecl_rst_file_type* start_solution (EclipseRestart& rst_file) {
        ecl_rst_file_start_solution (rst_file);
        return rst_file;
    }
};

/**
 * Representation of an Eclipse grid.
 */
struct EclipseGrid : public EclipseHandle <ecl_grid_type> {
    /// Create a grid based on the keywords available in input file
    static EclipseGrid make (const EclipseGridParser& parser,
                             const UnstructuredGrid& grid) {
        if (parser.hasField("DXV")) {
            // make sure that the DYV and DZV keywords are present if the
            // DXV keyword is used in the deck...
            assert(parser.hasField("DYV"));
            assert(parser.hasField("DZV"));

            const auto& dxv = parser.getFloatingPointValue("DXV");
            const auto& dyv = parser.getFloatingPointValue("DYV");
            const auto& dzv = parser.getFloatingPointValue("DZV");

            return EclipseGrid (dxv, dyv, dzv);
        }
        else if (parser.hasField("ZCORN")) {
            struct grdecl g = parser.get_grdecl ();

            auto coordData = parser.getFloatingPointValue(COORD_KW);
            auto zcornData = parser.getFloatingPointValue(ZCORN_KW);

            EclipseKeyword<float> coord_kw (COORD_KW, coordData);
            EclipseKeyword<float> zcorn_kw (ZCORN_KW, zcornData);

            // get the actually active cells, after processing
            std::vector <int> actnum;
            getActiveCells_(grid, actnum);
            EclipseKeyword<int> actnum_kw (ACTNUM_KW, actnum);

            EclipseKeyword<float> mapaxes_kw (MAPAXES_KW);
            if (g.mapaxes) {
                auto mapaxesData = parser.getFloatingPointValue(MAPAXES_KW);
                mapaxes_kw = std::move (EclipseKeyword<float> (MAPAXES_KW, mapaxesData));
            }

            return EclipseGrid (g.dims, zcorn_kw, coord_kw, actnum_kw, mapaxes_kw);
        }
        else {
            OPM_THROW(std::runtime_error,
                  "Can't create an ERT grid (no supported keywords found in deck)");
        }
    }

    /// Create a grid based on the keywords available in input file
    static EclipseGrid make (Opm::DeckConstPtr newParserDeck,
                             const UnstructuredGrid& grid)
    {
        if (newParserDeck->hasKeyword("DXV")) {
            // make sure that the DYV and DZV keywords are present if the
            // DXV keyword is used in the deck...
            assert(newParserDeck->hasKeyword("DYV"));
            assert(newParserDeck->hasKeyword("DZV"));

            const auto& dxv = newParserDeck->getKeyword("DXV")->getSIDoubleData();
            const auto& dyv = newParserDeck->getKeyword("DYV")->getSIDoubleData();
            const auto& dzv = newParserDeck->getKeyword("DZV")->getSIDoubleData();

            return EclipseGrid (dxv, dyv, dzv);
        }
        else if (newParserDeck->hasKeyword("ZCORN")) {
            struct grdecl g;
            GridManager::createGrdecl(newParserDeck, g);

            auto coordData = getAllSiDoubles_(newParserDeck->getKeyword(COORD_KW));
            auto zcornData = getAllSiDoubles_(newParserDeck->getKeyword(ZCORN_KW));
            EclipseKeyword<float> coord_kw (COORD_KW, coordData);
            EclipseKeyword<float> zcorn_kw (ZCORN_KW, zcornData);

            // get the actually active cells, after processing
            std::vector <int> actnum;
            getActiveCells_(grid, actnum);
            EclipseKeyword<int> actnum_kw (ACTNUM_KW, actnum);

            EclipseKeyword<float> mapaxes_kw (MAPAXES_KW);
            if (g.mapaxes) {
                auto mapaxesData = getAllSiDoubles_(newParserDeck->getKeyword(MAPAXES_KW));
                mapaxes_kw = std::move (EclipseKeyword<float> (MAPAXES_KW, mapaxesData));
            }

            return EclipseGrid (g.dims, zcorn_kw, coord_kw, actnum_kw, mapaxes_kw);
        }
        else {
            OPM_THROW(std::runtime_error,
                  "Can't create an ERT grid (no supported keywords found in deck)");
        }
    }

    /**
     * Save the grid in an .EGRID file.
     */
    void write (const std::string& outputDir,
                const std::string& baseName,
                int outputStepIdx) {
        ecl_grid_fwrite_EGRID (*this,
                               EclipseFileName (outputDir,
                                                baseName,
                                                ECL_EGRID_FILE,
                                                outputStepIdx));
    }

    // GCC 4.4 doesn't generate these constructors for us; provide the
    // default implementation explicitly here instead
    EclipseGrid (EclipseGrid&& rhs)
        : EclipseHandle <ecl_grid_type> (std::move (rhs)) { }
    EclipseGrid& operator= (EclipseGrid&& rhs) {
        EclipseHandle <ecl_grid_type>::operator= (std::move(rhs));
        return *this;
    }
    EclipseGrid (const EclipseGrid&) = delete;
    EclipseGrid& operator= (const EclipseGrid&) = delete;

private:
    // each of these cases could have been their respective subclass,
    // but there is not any polymorphism on each of these grid types
    // once we have the handle

    // setup smart pointer for Cartesian grid
    EclipseGrid (const std::vector<double>& dxv,
                 const std::vector<double>& dyv,
                 const std::vector<double>& dzv)
        : EclipseHandle <ecl_grid_type> (
              ecl_grid_alloc_dxv_dyv_dzv (dxv.size (),
                                          dyv.size (),
                                          dzv.size (),
                                          &dxv[0],
                                          &dyv[0],
                                          &dzv[0],
                                          NULL),
              ecl_grid_free) { }

    // setup smart pointer for cornerpoint grid
    EclipseGrid (const int dims[],
                 const EclipseKeyword<float>& zcorn,
                 const EclipseKeyword<float>& coord,
                 const EclipseKeyword<int>& actnum,
                 const EclipseKeyword<float>& mapaxes)
        : EclipseHandle <ecl_grid_type> (
              ecl_grid_alloc_GRDECL_kw(dims[0],
                                       dims[1],
                                       dims[2],
                                       zcorn,
                                       coord,
                                       actnum,
                                       mapaxes),
              ecl_grid_free) { }
};

/**
 * Initialization file which contains static properties (such as
 * porosity and permeability) for the simulation field.
 */
struct EclipseInit : public EclipseHandle <fortio_type> {
    // contrary to the grid, the location of the file goes here because
    // there is only one construction method but several write methods
    // (but we need to do a bit of logic before we can call the actual
    // constructor, so we'll have to do with a static wrapper)
    static EclipseInit make (const std::string& outputDir,
                             const std::string& baseName,
                             int outputStepIdx) {
        EclipseFileName initFileName (outputDir,
                                      baseName,
                                      ECL_INIT_FILE,
                                      outputStepIdx);
        bool fmt_file;
        if (!ecl_util_fmt_file(initFileName, &fmt_file)) {
            OPM_THROW(std::runtime_error,
                      "Could not determine formatted/unformatted status of file:" << initFileName << " non-standard name?" << std::endl);
        }
        return EclipseInit (initFileName, fmt_file);
    }

    void writeHeader (const EclipseGrid& grid,
                      const SimulatorTimer& timer,
                      const EclipseGridParser& parser,
                      const PhaseUsage uses) {
        EclipseKeyword<float> poro (PORO_KW, parser);
        ecl_init_file_fwrite_header (*this,
                                     grid,
                                     poro,
                                     phaseMask (uses),
                                     timer.currentPosixTime());
    }

    void writeHeader (const UnstructuredGrid& grid,
                      const SimulatorTimer& timer,
                      Opm::DeckConstPtr newParserDeck,
                      const PhaseUsage uses)
    {
        auto dataField = getAllSiDoubles_(newParserDeck->getKeyword(PORO_KW));
        restrictToActiveCells_(dataField, grid);

        EclipseGrid eclGrid = EclipseGrid::make (newParserDeck, grid);

        EclipseKeyword<float> poro (PORO_KW, dataField);
        ecl_init_file_fwrite_header (*this,
                                     eclGrid,
                                     poro,
                                     phaseMask (uses),
                                     timer.currentPosixTime ());
    }

    void writeKeyword (const std::string& keywordName,
                       const EclipseGridParser& parser,
                       double (* const transf)(const double&)) {
        auto data = parser.getValue <double> (keywordName);
        convertUnit_(data, transf);
        EclipseKeyword <float> kw (keywordName, data);
        ecl_kw_fwrite (kw, *this);
    }

    void writeKeyword (const std::string& keywordName, const std::vector<double> &data)
    {
        EclipseKeyword <float> kw (keywordName, data);
        ecl_kw_fwrite(kw, *this);
    }

    // GCC 4.4 doesn't generate these constructors for us; provide the
    // default implementation explicitly here instead
    EclipseInit (EclipseInit&& rhs)
        : EclipseHandle <fortio_type> (std::move (rhs)) { }
    EclipseInit& operator= (EclipseInit& rhs) {
        EclipseHandle <fortio_type>::operator= (std::move(rhs));
        return *this;
    }
    EclipseInit (const EclipseInit&) = delete;
    EclipseInit& operator= (const EclipseInit&) = delete;
private:
    EclipseInit (const EclipseFileName& fname, const bool formatted)
        : EclipseHandle <fortio_type> (
              fortio_open_writer (fname, formatted, ECL_ENDIAN_FLIP),
              fortio_fclose) { }
};


// in order to get RTTI for this "class" (which is just a typedef), we must
// ask the compiler to explicitly instantiate it.
template struct EclipseHandle<ecl_sum_tstep_struct>;


} // anonymous namespace

// Note: the following parts were taken out of the anonymous
// namespace, since EclipseSummary is now used as a pointer member in
// EclipseWriter and forward declared in EclipseWriter.hpp.

// forward decl. of mutually dependent type
struct EclipseWellReport;

struct EclipseSummary : public EclipseHandle <ecl_sum_type> {
    EclipseSummary (const std::string& outputDir,
                    const std::string& baseName,
                    const SimulatorTimer& timer,
                    const EclipseGridParser parser)
        : EclipseHandle <ecl_sum_type> (
              alloc_writer (outputDir, baseName, timer, parser),
              ecl_sum_free) { }

    EclipseSummary (const std::string& outputDir,
                    const std::string& baseName,
                    const SimulatorTimer& timer,
                    Opm::DeckConstPtr newParserDeck)
        : EclipseHandle <ecl_sum_type> (
              alloc_writer (outputDir, baseName, timer, newParserDeck),
              ecl_sum_free) { }

    typedef std::unique_ptr <EclipseWellReport> var_t;
    typedef std::vector <var_t> vars_t;

    EclipseSummary& add (var_t var) {
        vars_.push_back (std::move (var));
        return *this;
    }

    // make sure the summary section is flushed before it goes away
    // (this will happen before all the timesteps are individually
    // destroyed, so their memory is still valid at this point)
    ~EclipseSummary () {
        ecl_sum_fwrite (*this);
    }

    // add rate variables for each of the well in the input file
    void addWells (const EclipseGridParser& parser,
                   const PhaseUsage& uses);

    // add rate variables for each of the well in the input file
    void addWells (Opm::DeckConstPtr newParserDeck,
                   const PhaseUsage& uses);

    // no inline implementation of this since it depends on the
    // EclipseWellReport type being completed first
    void writeTimeStep (const SimulatorTimer& timer,
                         const WellState& wellState);

private:
    vars_t vars_;

    // don't define a new type for timesteps (since they should all
    // be created with makeTimeStep anyway), just use the basic handle
    // type and a typedef.
    typedef EclipseHandle <ecl_sum_tstep_type> EclipseTimeStep;

    /// Create a new time step and add it to this summary. The summary
    /// will take care of memory management, the object returned is a
    /// "view" into it. Make sure that that view does not outlive the
    /// summary object! Notice that there is no deleter in the constructor.
    std::unique_ptr <EclipseTimeStep> makeTimeStep (const SimulatorTimer& timer) {
        EclipseTimeStep* tstep = new EclipseTimeStep (
                    ecl_sum_add_tstep (*this,
                                       timer.currentStepNum (),
                                       Opm::unit::convert::to (timer.simulationTimeElapsed (),
                                                               Opm::unit::day)));
        return std::unique_ptr <EclipseTimeStep> (tstep);
    }

    /// Helper routine that lets us use local variables to hold
    /// intermediate results while filling out the allocations function's
    /// argument list.
    static ecl_sum_type* alloc_writer (const std::string& outputDir,
                                        const std::string& baseName,
                                        const SimulatorTimer& timer,
                                        const EclipseGridParser& parser) {
        boost::filesystem::path casePath (outputDir);
        casePath /= boost::to_upper_copy (baseName);

        const std::vector <int> dim = parserDim (parser);
        return ecl_sum_alloc_writer (casePath.string ().c_str (),
                                     false, /* formatted   */
                                     true,  /* unified     */
                                     ":",    /* join string */
                                     timer.simulationTimeElapsed (),
                                     dim[0],
                                     dim[1],
                                     dim[2]);
    }

    /// Helper routine that lets us use local variables to hold
    /// intermediate results while filling out the allocations function's
    /// argument list.
    static ecl_sum_type* alloc_writer (const std::string& outputDir,
                                       const std::string& baseName,
                                       const SimulatorTimer& timer,
                                       Opm::DeckConstPtr newParserDeck) {
        boost::filesystem::path casePath (outputDir);
        casePath /= boost::to_upper_copy (baseName);

        const std::vector <int> dim = parserDim (newParserDeck);
        return ecl_sum_alloc_writer (casePath.string ().c_str (),
                                     false, /* formatted   */
                                     true,  /* unified     */
                                     ":",    /* join string */
                                     timer.simulationTimeElapsed (),
                                     dim[0],
                                     dim[1],
                                     dim[2]);
    }
};


/**
 * Summary variable that reports a characteristics of a well.
 */
struct EclipseWellReport : public EclipseHandle <smspec_node_type> {
protected:

    EclipseWellReport (const EclipseSummary& summary,    /* section to add to  */
                       Opm::DeckConstPtr newParserDeck,  /* well names         */
                       int whichWell,                    /* index of well line */
                       PhaseUsage uses,                  /* phases present     */
                       BlackoilPhases::PhaseIndex phase, /* oil, water or gas  */
                       WellType type,                    /* prod. or inj.      */
                       char aggregation,                 /* rate or total      */
                       std::string unit)
        : EclipseHandle <smspec_node_type> (
              ecl_sum_add_var (summary,
                               varName (phase,
                                        type,
                                        aggregation).c_str (),
                               wellName (newParserDeck, whichWell).c_str (),
                               /* num = */ 0,
                               unit.c_str(),
                               /* defaultValue = */ 0.))
        // save these for when we update the value in a timestep
        , index_ (whichWell * uses.num_phases + uses.phase_pos [phase])

        // producers can be seen as negative injectors
        , sign_ (type == INJECTOR ? +1. : -1.) { }

public:
    /// Allows us to pass this type to ecl_sum_tstep_iset
    operator int () {
        return smspec_node_get_params_index (*this);
    }

    /// Update the monitor according to the new state of the well, and
    /// get the reported value. Note: Only call this once for each timestep.
    virtual double update (const SimulatorTimer& timer,
                             const WellState& wellState) = 0;

private:
    /// index into a (flattened) wells*phases matrix
    const int index_;

    /// natural sign of the rate
    const double sign_;

    /// Get the name associated with this well
    /*std::string wellName (const EclipseGridParser& parser,
                          int whichWell) {
        return parser.getWELSPECS().welspecs[whichWell].name_;
    }
    */

    /// Get the name associated with this well
    std::string wellName (Opm::DeckConstPtr newParserDeck,
                          int whichWell)
    {
        Opm::WelspecsWrapper welspecs(newParserDeck->getKeyword("WELSPECS"));
        return welspecs.wellName(whichWell);
    }

    /// Compose the name of the summary variable, e.g. "WOPR" for
    /// well oil production rate.
    std::string varName (BlackoilPhases::PhaseIndex phase,
                         WellType type,
                         char aggregation) {
        std::string name;
        name += 'W'; // well
        if (aggregation == 'B') {
            name += "BHP";
        } else {
            switch (phase) {
            case BlackoilPhases::Aqua:   name += 'W'; break; /* water */
            case BlackoilPhases::Vapour: name += 'G'; break; /* gas */
            case BlackoilPhases::Liquid: name += 'O'; break; /* oil */
            default:
                OPM_THROW(std::runtime_error,
                          "Unknown phase used in blackoil reporting");
            }
            switch (type) {
            case WellType::INJECTOR: name += 'I'; break;
            case WellType::PRODUCER: name += 'P'; break;
            default:
                OPM_THROW(std::runtime_error,
                          "Unknown well type used in blackoil reporting");
            }
            name += aggregation; /* rate ('R') or total ('T') */
        }
        return name;
    }
protected:
    double rate (const WellState& wellState) {
        // convert m^3/s of injected fluid to m^3/d of produced fluid
        const double convFactor = Opm::unit::convert::to (1., Opm::unit::day);
        const double value = sign_ * wellState.wellRates () [index_] * convFactor;
        return value;
    }

    double bhp (const WellState& wellstate) {
        // Note that 'index_' is used here even though it is meant
        // to give a (well,phase) pair.
        const int num_phases = wellstate.wellRates().size() / wellstate.bhp().size();
        return wellstate.bhp()[index_/num_phases];
    }
};

/// Monitors the rate given by a well.
struct EclipseWellRate : public EclipseWellReport {
    EclipseWellRate (const EclipseSummary& summary,
                     Opm::DeckConstPtr newParserDeck,
                     int whichWell,
                     PhaseUsage uses,
                     BlackoilPhases::PhaseIndex phase,
                     WellType type)
        : EclipseWellReport (summary,
                             newParserDeck,
                             whichWell,
                             uses,
                             phase,
                             type,
                             'R',
                             "SM3/DAY" /* surf. cub. m. per day */ ) { }

    virtual double update (const SimulatorTimer& /*timer*/,
                             const WellState& wellState) {
        // TODO: Why only positive rates?
        return std::max (0., rate (wellState));
    }
};

/// Monitors the total production in a well.
struct EclipseWellTotal : public EclipseWellReport {

    EclipseWellTotal (const EclipseSummary& summary,
                      Opm::DeckConstPtr newParserDeck,
                      int whichWell,
                      PhaseUsage uses,
                      BlackoilPhases::PhaseIndex phase,
                      WellType type)
        : EclipseWellReport (summary,
                             newParserDeck,
                             whichWell,
                             uses,
                             phase,
                             type,
                             'T',
                             "SM3" /* surface cubic meter */ )

        // nothing produced when the reporting starts
        , total_ (0.) { }

    virtual double update (const SimulatorTimer& timer,
                             const WellState& wellState) {
        if (timer.currentStepNum() == 0) {
            // We are at the initial state.
            // No step has been taken yet.
            return 0.0;
        }
        // TODO: Is the rate average for the timestep, or is in
        // instantaneous (in which case trapezoidal or Simpson integration
        // would probably be better)
        const double intg = timer.stepLengthTaken () * rate (wellState);
        // add this timesteps production to the total
        total_ += intg;
        // report the new production total
        return total_;
    }

private:
    /// Aggregated value of the course of the simulation
    double total_;
};

/// Monitors the bottom hole pressure in a well.
struct EclipseWellBhp : public EclipseWellReport {
    EclipseWellBhp   (const EclipseSummary& summary,
                      Opm::DeckConstPtr newParserDeck,
                      int whichWell,
                      PhaseUsage uses,
                      BlackoilPhases::PhaseIndex phase,
                      WellType type)
        : EclipseWellReport (summary,
                             newParserDeck,
                             whichWell,
                             uses,
                             phase,
                             type,
                             'B',
                             "Pascal")
    { }

    virtual double update (const SimulatorTimer& /*timer*/,
                           const WellState& wellState)
    {
        return bhp(wellState);
    }
};

inline void
EclipseSummary::writeTimeStep (const SimulatorTimer& timer,
                               const WellState& wellState)
{
    // internal view; do not move this code out of EclipseSummary!
    std::unique_ptr <EclipseTimeStep> tstep = makeTimeStep (timer);
    // write all the variables
    for (vars_t::iterator v = vars_.begin(); v != vars_.end(); ++v) {
        const double value = (*v)->update (timer, wellState);
        ecl_sum_tstep_iset(*tstep, *(*v).get (), value);
    }
}

/// Supported well types. Enumeration doesn't let us get all the members,
/// so we must have an explicit array.
static WellType WELL_TYPES[] = { INJECTOR, PRODUCER };


inline void
EclipseSummary::addWells (Opm::DeckConstPtr newParserDeck,
                          const PhaseUsage& uses) {
    // TODO: Only create report variables that are requested with keywords
    // (e.g. "WOPR") in the input files, and only for those wells that are
    // mentioned in those keywords
    Opm::DeckKeywordConstPtr welspecsKeyword = newParserDeck->getKeyword("WELSPECS");
    const int numWells = welspecsKeyword->size();
    for (int phaseCounter = 0;
          phaseCounter != BlackoilPhases::MaxNumPhases;
          ++phaseCounter) {
        const BlackoilPhases::PhaseIndex phase =
                static_cast <BlackoilPhases::PhaseIndex> (phaseCounter);
        // don't bother with reporting for phases that aren't there
        if (!uses.phase_used [phaseCounter]) {
            continue;
        }
        for (size_t typeIndex = 0;
             typeIndex < sizeof (WELL_TYPES) / sizeof (WELL_TYPES[0]);
             ++typeIndex) {
            const WellType type = WELL_TYPES[typeIndex];
            for (int whichWell = 0; whichWell != numWells; ++whichWell) {
                // W{O,G,W}{I,P}R
                add (std::unique_ptr <EclipseWellReport> (
                              new EclipseWellRate (*this,
                                                   newParserDeck,
                                                   whichWell,
                                                   uses,
                                                   phase,
                                                   type)));
                // W{O,G,W}{I,P}T
                add (std::unique_ptr <EclipseWellReport> (
                              new EclipseWellTotal (*this,
                                                    newParserDeck,
                                                    whichWell,
                                                    uses,
                                                    phase,
                                                    type)));
            }
        }
    }

    // Add BHP monitors
    for (int whichWell = 0; whichWell != numWells; ++whichWell) {
        // In the call below: uses, phase and the well type arguments
        // are not used, except to set up an index that stores the
        // well indirectly. For details see the implementation of the
        // EclipseWellReport constructor, and the method
        // EclipseWellReport::bhp().
        BlackoilPhases::PhaseIndex phase = BlackoilPhases::Liquid;
        if (!uses.phase_used[BlackoilPhases::Liquid]) {
            phase = BlackoilPhases::Vapour;
        }
        add (std::unique_ptr <EclipseWellReport> (
                        new EclipseWellBhp (*this,
                                            newParserDeck,
                                            whichWell,
                                            uses,
                                            phase,
                                            WELL_TYPES[0])));
    }

}

namespace Opm {

void EclipseWriter::writeInit(const SimulatorTimer &timer,
                              const SimulatorState& reservoirState,
                              const WellState& wellState)
{
    // if we don't want to write anything, this method becomes a
    // no-op...
    if (!enableOutput_) {
        return;
    }
    if (newParserDeck_) {
        /* Grid files */
        EclipseGrid eclGrid = EclipseGrid::make (newParserDeck_, *grid_);
        eclGrid.write (outputDir_, baseName_, /*stepIdx=*/0);

        EclipseInit fortio = EclipseInit::make (outputDir_, baseName_, /*stepIdx=*/0);
        fortio.writeHeader (*grid_,
                            timer,
                            newParserDeck_,
                            uses_);

        if (newParserDeck_->hasKeyword("PERM")) {
            auto data = getAllSiDoubles_(newParserDeck_->getKeyword("PERM"));
            convertUnit_(data, toMilliDarcy);
            fortio.writeKeyword ("PERM", data);
        }

        if (newParserDeck_->hasKeyword("PERMX")) {
            auto data = getAllSiDoubles_(newParserDeck_->getKeyword("PERMX"));
            convertUnit_(data, toMilliDarcy);
            fortio.writeKeyword ("PERMX", data);
        }
        if (newParserDeck_->hasKeyword("PERMY")) {
            auto data = getAllSiDoubles_(newParserDeck_->getKeyword("PERMY"));
            convertUnit_(data, toMilliDarcy);
            fortio.writeKeyword ("PERMY", data);
        }
        if (newParserDeck_->hasKeyword("PERMZ")) {
            auto data = getAllSiDoubles_(newParserDeck_->getKeyword("PERMZ"));
            convertUnit_(data, toMilliDarcy);
            fortio.writeKeyword ("PERMZ", data);
        }

        /* Initial solution (pressure and saturation) */
        writeSolution_(timer, reservoirState);

        /* Create summary object (could not do it at construction time,
           since it requires knowledge of the start time). */
        summary_.reset(new EclipseSummary(outputDir_, baseName_, timer, newParserDeck_));
        summary_->addWells (newParserDeck_, uses_);
    }
    else {
        /* Grid files */
        EclipseGrid ecl_grid = EclipseGrid::make (*parser_, *grid_);
        ecl_grid.write (outputDir_, baseName_, /*stepIdx=*/0);

        EclipseInit fortio = EclipseInit::make (outputDir_, baseName_, /*stepIdx=*/0);
        fortio.writeHeader (ecl_grid,
                            timer,
                            *parser_,
                            uses_);

        fortio.writeKeyword ("PERMX", *parser_, &toMilliDarcy);
        fortio.writeKeyword ("PERMY", *parser_, &toMilliDarcy);
        fortio.writeKeyword ("PERMZ", *parser_, &toMilliDarcy);

        /* Initial solution (pressure and saturation) */
        writeSolution_(timer, reservoirState);

        /* Create summary object (could not do it at construction time,
           since it requires knowledge of the start time). */
        summary_.reset(new EclipseSummary(outputDir_, baseName_, timer, *parser_));
        summary_->addWells (*parser_, uses_);
    }
}

void EclipseWriter::activeToGlobalCellData_(std::vector<double> &globalCellsBuf,
                                              const std::vector<double> &activeCellsBuf,
                                              const std::vector<double> &inactiveCellsBuf) const
{
    globalCellsBuf = inactiveCellsBuf;

    // overwrite the values of active cells
    for (int activeCellIdx = 0;
         activeCellIdx < grid_->number_of_cells;
         ++activeCellIdx)
    {
        int globalCellIdx = grid_->global_cell[activeCellIdx];
        globalCellsBuf[globalCellIdx] = activeCellsBuf[activeCellIdx];
    }
}

void EclipseWriter::writeSolution_(const SimulatorTimer& timer,
                                   const SimulatorState& reservoirState)
{
    if (newParserDeck_) {
        // start writing to files
        EclipseRestart rst(outputDir_, baseName_, timer, outputTimeStepIdx_);
        rst.writeHeader (timer, outputTimeStepIdx_, uses_, newParserDeck_, reservoirState.pressure().size ());
        EclipseSolution sol (rst);

        // write out the pressure of the reference phase (whatever
        // phase that is...). this is not the most performant solution
        // thinkable, but this is also not in the most performance
        // critical code path!
        std::vector<double> tmp = reservoirState.pressure();
        convertUnit_(tmp, toBar);

        sol.add(EclipseKeyword<float>("PRESSURE", tmp));

        for (int phase = 0; phase != BlackoilPhases::MaxNumPhases; ++phase) {
            // Eclipse never writes the oil saturation, so all post-processors
            // must calculate this from the other saturations anyway
            if (phase == BlackoilPhases::PhaseIndex::Liquid) {
                continue;
            }
            if (uses_.phase_used [phase]) {
                tmp = reservoirState.saturation();
                extractFromStripedData_(tmp,
                                        /*offset=*/uses_.phase_pos[phase],
                                        /*stride=*/uses_.num_phases);
                sol.add(EclipseKeyword<float>(SAT_NAMES[phase], tmp));
            }
        }
    }
    else {
        // start writing to files
        EclipseRestart rst (outputDir_,
                            baseName_,
                            timer,
                            outputTimeStepIdx_);
        rst.writeHeader (timer,
                         outputTimeStepIdx_,
                         uses_,
                         *parser_,
                         reservoirState.pressure ().size ());
        EclipseSolution sol (rst);

        // write pressure and saturation fields (same as DataMap holds)
        // convert the pressures from Pascals to bar because Eclipse
        // seems to write bars
        auto data = reservoirState.pressure();
        convertUnit_(data, toBar);
        sol.add(EclipseKeyword<float>("PRESSURE", data));

        for (int phase = 0; phase != BlackoilPhases::MaxNumPhases; ++phase) {
            // Eclipse never writes the oil saturation, so all post-processors
            // must calculate this from the other saturations anyway
            if (phase == BlackoilPhases::PhaseIndex::Liquid) {
                continue;
            }
            if (uses_.phase_used [phase]) {
                auto tmp = reservoirState.saturation();
                extractFromStripedData_(tmp,
                                        /*offset=*/uses_.phase_pos[phase],
                                        /*stride=*/uses_.num_phases);
                sol.add (EclipseKeyword<float>(SAT_NAMES[phase], tmp));
            }
        }
    }

    ++outputTimeStepIdx_;
}

void EclipseWriter::writeTimeStep(const SimulatorTimer& timer,
                                  const SimulatorState& reservoirState,
                                  const WellState& wellState)
{
    // if we don't want to write anything, this method becomes a
    // no-op...
    if (!enableOutput_) {
        return;
    }

    // respected the output_interval parameter
    if (timer.currentStepNum() % outputInterval_ != 0) {
        return;
    }

    /* Field variables (pressure, saturation) */
    writeSolution_(timer, reservoirState);

    /* Summary variables (well reporting) */
    // TODO: instead of writing the header (smspec) every time, it should
    // only be written when there is a change in the well configuration
    // (first timestep, in practice), and reused later. but how to do this
    // without keeping the complete summary in memory (which will then
    // accumulate all the timesteps)?
    //
    // Note: The answer to the question above is still not settled,
    // but now we do keep the complete summary in memory, as a member
    // variable in the EclipseWriter class, instead of creating a
    // temporary EclipseSummary in this function every time it is
    // called.  This has been changed so that the final summary file
    // will contain data from the whole simulation, instead of just
    // the last step.
    summary_->writeTimeStep(timer, wellState);
}

#else
namespace Opm {

void EclipseWriter::writeInit(const SimulatorTimer&,
                              const SimulatorState&,
                              const WellState&)
 {
    // if we don't want to write anything, this method becomes a
    // no-op...
     if (!enableOutput_) {
        return;
     }

    OPM_THROW(std::runtime_error,
              "The ERT libraries are required to write ECLIPSE output files.");
}

void EclipseWriter::writeTimeStep(
        const SimulatorTimer&,
        const SimulatorState&,
        const WellState&)
{
    // if we don't want to write anything, this method becomes a
    // no-op...
    if (!enableOutput_) {
        return;
    }

    OPM_THROW(std::runtime_error,
              "The ERT libraries are required to write ECLIPSE output files.");
}

#endif // HAVE_ERT

EclipseWriter::EclipseWriter (
        const ParameterGroup& params,
        std::shared_ptr <EclipseGridParser> parser,
        std::shared_ptr <const UnstructuredGrid> grid)
    : parser_ (parser)
    , newParserDeck_()
    , grid_ (grid)
    , uses_ (phaseUsageFromDeck (*parser))
{
    // set the index of the first time step written to 0...
    outputTimeStepIdx_ = 0;


    // get the base name from the name of the deck
    using boost::filesystem::path;
    path deck (params.get <std::string> ("deck_filename"));
    if (boost::to_upper_copy (path (deck.extension ()).string ()) == ".DATA") {
        baseName_ = path (deck.stem ()).string ();
    }
    else {
        baseName_ = path (deck.filename ()).string ();
    }

    // make uppercase of everything (or otherwise we'll get uppercase
    // of some of the files (.SMSPEC, .UNSMRY) and not others
    baseName_ = boost::to_upper_copy (baseName_);

    // retrieve the value of the "output" parameter
    enableOutput_ = params.getDefault<bool>("output", /*defaultValue=*/true);

    // retrieve the interval at which something should get written to
    // disk (once every N timesteps)
    outputInterval_ = params.getDefault<int>("output_interval", /*defaultValue=*/1);

    // store in current directory if not explicitly set
    outputDir_ = params.getDefault<std::string>("output_dir", ".");

    // set the index of the first time step written to 0...
    outputTimeStepIdx_ = 0;

    if (enableOutput_) {
        // make sure that the output directory exists, if not try to create it
        if (!boost::filesystem::exists(outputDir_)) {
            std::cout << "Trying to create directory \"" << outputDir_ << "\" for the simulation output\n";
            boost::filesystem::create_directories(outputDir_);
        }

        if (!boost::filesystem::is_directory(outputDir_)) {
            OPM_THROW(std::runtime_error,
                      "The path specified as output directory '" << outputDir_
                      << "' is not a directory");
        }
    }
}

EclipseWriter::EclipseWriter (
        const ParameterGroup& params,
        Opm::DeckConstPtr newParserDeck,
        std::shared_ptr <const UnstructuredGrid> grid)
    : parser_ ()
    , newParserDeck_(newParserDeck)
    , grid_ (grid)
    , uses_ (phaseUsageFromDeck (newParserDeck))
{
    // get the base name from the name of the deck
    using boost::filesystem::path;
    path deck (params.get <std::string> ("deck_filename"));
    if (boost::to_upper_copy (path (deck.extension ()).string ()) == ".DATA") {
        baseName_ = path (deck.stem ()).string ();
    }
    else {
        baseName_ = path (deck.filename ()).string ();
    }

    // make uppercase of everything (or otherwise we'll get uppercase
    // of some of the files (.SMSPEC, .UNSMRY) and not others
    baseName_ = boost::to_upper_copy (baseName_);

    // retrieve the value of the "output" parameter
    enableOutput_ = params.getDefault<bool>("output", /*defaultValue=*/true);

    // retrieve the interval at which something should get written to
    // disk (once every N timesteps)
    outputInterval_ = params.getDefault<int>("output_interval", /*defaultValue=*/1);

    // store in current directory if not explicitly set
    outputDir_ = params.getDefault<std::string>("output_dir", ".");

    // set the index of the first time step written to 0...
    outputTimeStepIdx_ = 0;

    if (enableOutput_) {
        // make sure that the output directory exists, if not try to create it
        if (!boost::filesystem::exists(outputDir_)) {
            std::cout << "Trying to create directory \"" << outputDir_ << "\" for the simulation output\n";
            boost::filesystem::create_directories(outputDir_);
        }

        if (!boost::filesystem::is_directory(outputDir_)) {
            OPM_THROW(std::runtime_error,
                      "The path specified as output directory '" << outputDir_
                      << "' is not a directory");
        }
    }
}

// default destructor is OK, just need to be defined
EclipseWriter::~EclipseWriter() { }

} // namespace Opm
