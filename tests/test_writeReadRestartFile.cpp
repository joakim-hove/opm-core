
/*
  Copyright 2014 Statoil IT
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

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define BOOST_TEST_MODULE EclipseWriter
#include <boost/test/unit_test.hpp>

#include <opm/core/io/eclipse/EclipseWriter.hpp>
#include <opm/core/io/eclipse/EclipseReader.hpp>
#include <opm/core/io/eclipse/EclipseIOUtil.hpp>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/grid/GridHelpers.hpp>
#include <opm/core/props/phaseUsageFromDeck.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/wells.h>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseMode.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/props/IncompPropertiesFromDeck.hpp>
#include <opm/core/simulator/TwophaseState.hpp>

// ERT stuff
#include <ert/ecl/ecl_kw.h>
#include <ert/ecl/ecl_file.h>
#include <ert/ecl/ecl_kw_magic.h>
#include <ert/ecl_well/well_info.h>
#include <ert/ecl_well/well_state.h>
#include <ert/util/test_work_area.h>

#include <string.h>


std::string input =
           "RUNSPEC\n"
           "OIL\n"
           "GAS\n"
           "WATER\n"
           "DISGAS\n"
           "VAPOIL\n"
           "UNIFOUT\n"
           "UNIFIN\n"
           "DIMENS\n"
           " 10 10 10 /\n"

           "GRID\n"
           "DXV\n"
           "10*0.25 /\n"
           "DYV\n"
           "10*0.25 /\n"
           "DZV\n"
           "10*0.25 /\n"
           "TOPS\n"
           "100*0.25 /\n"
           "\n"

           "SOLUTION\n"
           "RESTART\n"
           "TESTWELLSTATE 1/\n"
           "\n"

           "START             -- 0 \n"
           "1 NOV 1979 / \n"

           "SCHEDULE\n"
           "SKIPREST\n"
           "RPTRST\n"
           "BASIC=1\n"
           "/\n"
           "DATES             -- 1\n"
           " 10  OKT 2008 / \n"
           "/\n"
           "WELSPECS\n"
           "    'OP_1'       'OP'   9   9 1*     'OIL' 1*      1*  1*   1*  1*   1*  1*  / \n"
           "    'OP_2'       'OP'   9   9 1*     'OIL' 1*      1*  1*   1*  1*   1*  1*  / \n"
           "/\n"
           "COMPDAT\n"
           " 'OP_1'  9  9   1   1 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
           " 'OP_2'  9  9   2   2 'OPEN' 1*   46.825   0.311  4332.346 1*  1*  'X'  22.123 / \n"
           " 'OP_1'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
           "/\n"
           "WCONPROD\n"
               "'OP_1' 'OPEN' 'ORAT' 20000  4* 1000 /\n"
           "/\n"
           "WCONINJE\n"
               "'OP_2' 'GAS' 'OPEN' 'RATE' 100 200 400 /\n"
           "/\n"

           "DATES             -- 2\n"
           " 20  JAN 2011 / \n"
           "/\n"
           "WELSPECS\n"
           "    'OP_3'       'OP'   9   9 1*     'OIL' 1*      1*  1*   1*  1*   1*  1*  / \n"
           "/\n"
           "COMPDAT\n"
           " 'OP_3'  9  9   1   1 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
           "/\n"
           "WCONPROD\n"
               "'OP_3' 'OPEN' 'ORAT' 20000  4* 1000 /\n"
           "/\n"

           "DATES             -- 3\n"
           " 15  JUN 2013 / \n"
           "/\n"
           "COMPDAT\n"
           " 'OP_2'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
           " 'OP_1'  9  9   7  7 'SHUT' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
           "/\n"

           "DATES             -- 4\n"
           " 22  APR 2014 / \n"
           "/\n"
           "WELSPECS\n"
           "    'OP_4'       'OP'   9   9 1*     'OIL' 1*      1*  1*   1*  1*   1*  1*  / \n"
           "/\n"
           "COMPDAT\n"
           " 'OP_4'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
           " 'OP_3'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
           "/\n"
           "WCONPROD\n"
               "'OP_4' 'OPEN' 'ORAT' 20000  4* 1000 /\n"
           "/\n"

           "DATES             -- 5\n"
           " 30  AUG 2014 / \n"
           "/\n"
           "WELSPECS\n"
           "    'OP_5'       'OP'   9   9 1*     'OIL' 1*      1*  1*   1*  1*   1*  1*  / \n"
           "/\n"
           "COMPDAT\n"
           " 'OP_5'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
           "/\n"
           "WCONPROD\n"
               "'OP_5' 'OPEN' 'ORAT' 20000  4* 1000 /\n"
           "/\n"

           "DATES             -- 6\n"
           " 15  SEP 2014 / \n"
           "/\n"
           "WCONPROD\n"
               "'OP_3' 'SHUT' 'ORAT' 20000  4* 1000 /\n"
           "/\n"

           "DATES             -- 7\n"
           " 9  OCT 2014 / \n"
           "/\n"
           "WELSPECS\n"
           "    'OP_6'       'OP'   9   9 1*     'OIL' 1*      1*  1*   1*  1*   1*  1*  / \n"
           "/\n"
           "COMPDAT\n"
           " 'OP_6'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
           "/\n"
           "WCONPROD\n"
               "'OP_6' 'OPEN' 'ORAT' 20000  4* 1000 /\n"
           "/\n"
           "TSTEP            -- 8\n"
           "10 /"
           "/\n";


std::shared_ptr<Opm::BlackoilState> createBlackOilState(Opm::EclipseGridConstPtr eclGrid , const Opm::PhaseUsage& phaseUsage) {

  std::shared_ptr<Opm::GridManager> grid(new Opm::GridManager(eclGrid));
  const UnstructuredGrid& ug_grid = *(grid->c_grid());
  std::shared_ptr<Opm::BlackoilState> blackoilState(new Opm::BlackoilState( Opm::UgGridHelpers::numCells(ug_grid) , Opm::UgGridHelpers::numFaces(ug_grid) , phaseUsage.num_phases) );

  return blackoilState;
}

Opm::EclipseWriterPtr createEclipseWriter(Opm::DeckConstPtr deck,
                                          Opm::EclipseStatePtr eclipseState,
                                          std::string& eclipse_data_filename) {

  Opm::parameter::ParameterGroup params;
  params.insertParameter("deck_filename", eclipse_data_filename);

  const Opm::PhaseUsage phaseUsage = Opm::phaseUsageFromDeck(deck);

  Opm::EclipseWriterPtr eclWriter(new Opm::EclipseWriter(params,
                                                         eclipseState,
                                                         phaseUsage,
                                                         eclipseState->getEclipseGrid()->getCartesianSize(),
                                                         0));
  return eclWriter;
}

void setValuesInWellState(std::shared_ptr<Opm::WellState> wellState){
    wellState->bhp()[0] = 1.23;
    wellState->bhp()[1] = 2.34;

    wellState->temperature()[0] = 3.45;
    wellState->temperature()[1] = 4.56;

    wellState->wellRates()[0] = 5.67;
    wellState->wellRates()[1] = 6.78;
    wellState->wellRates()[2] = 7.89;
    wellState->wellRates()[3] = 8.90;
    wellState->wellRates()[4] = 9.01;
    wellState->wellRates()[5] = 10.12;

    wellState->perfPress()[0] = 20.41;
    wellState->perfPress()[1] = 21.19;
    wellState->perfPress()[2] = 22.41;
    wellState->perfPress()[3] = 23.19;
    wellState->perfPress()[4] = 24.41;
    wellState->perfPress()[5] = 25.19;
    wellState->perfPress()[6] = 26.41;
    wellState->perfPress()[7] = 27.19;
    wellState->perfPress()[8] = 28.41;

    wellState->perfRates()[0] = 30.45;
    wellState->perfRates()[1] = 31.19;
    wellState->perfRates()[2] = 32.45;
    wellState->perfRates()[3] = 33.19;
    wellState->perfRates()[4] = 34.45;
    wellState->perfRates()[5] = 35.19;
    wellState->perfRates()[6] = 36.45;
    wellState->perfRates()[7] = 37.19;
    wellState->perfRates()[8] = 38.45;
}

BOOST_AUTO_TEST_CASE(EclipseReadWriteWellStateData)
{
    std::string eclipse_data_filename    = "TestWellState.DATA";
    test_work_area_type * test_area = test_work_area_alloc("EclipseReadWriteWellStateData");

    Opm::Parser parser;
    Opm::ParseMode parseMode;
    Opm::DeckConstPtr deck = parser.parseString(input, parseMode);
    Opm::EclipseStatePtr  eclipseState(new Opm::EclipseState(deck , parseMode));
    Opm::EclipseWriterPtr eclipseWriter = createEclipseWriter(deck, eclipseState, eclipse_data_filename);

    std::shared_ptr<Opm::SimulatorTimer> simTimer( new Opm::SimulatorTimer() );
    simTimer->init(eclipseState->getSchedule()->getTimeMap());
    eclipseWriter->writeInit(*simTimer);
    std::shared_ptr<Opm::WellState> wellState(new Opm::WellState());
    Opm::PhaseUsage phaseUsage = Opm::phaseUsageFromDeck(deck);

    Opm::GridManager gridManager(deck);
    Opm::WellsManager wellsManager(eclipseState, 1, *gridManager.c_grid(), NULL);
    const Wells* wells = wellsManager.c_wells();
    std::shared_ptr<Opm::BlackoilState> blackoilState = createBlackOilState(eclipseState->getEclipseGrid(), phaseUsage);
    wellState->init(wells, *blackoilState);

    //Set test data for pressure
    std::vector<double>& pressure = blackoilState->getCellData( Opm::BlackoilState::PRESSURE );
    std::fill( pressure.begin() , pressure.end() , 6.0 );

    //Set test data for temperature
    std::vector<double>& temperature = blackoilState->getCellData( Opm::BlackoilState::TEMPERATURE );
    std::fill( temperature.begin() , temperature.end() , 7 );

    std::vector<double>& saturation = blackoilState->getCellData( Opm::BlackoilState::SATURATION );
    //Set test data for saturation water
    std::vector<double> swatdata(1000, 8);
    Opm::EclipseIOUtil::addToStripedData(swatdata, saturation, phaseUsage.phase_pos[Opm::BlackoilPhases::Aqua], phaseUsage.num_phases);

    //Set test data for saturation gas
    std::vector<double> sgasdata(1000, 9);
    Opm::EclipseIOUtil::addToStripedData(sgasdata, saturation, phaseUsage.phase_pos[Opm::BlackoilPhases::Vapour], phaseUsage.num_phases);

    // Set test data for rs
    std::vector<double>& rs = blackoilState->getCellData( Opm::BlackoilState::GASOILRATIO );
    {    
	double current_rs = 300.0;
	std::for_each(rs.begin(), rs.end(), 
		      [&current_rs](double &rs)
		      { 
			  rs = current_rs;
			  current_rs++;
		      });
    }

    
    // Set testdata for rv
    std::vector<double>& rv = blackoilState->getCellData( Opm::BlackoilState::RV );
    {    
	double current_rv = 400.0;
	std::for_each(rv.begin(), rv.end(), 
		      [&current_rv](double &rv)
		      { 
			  rv = current_rv;
			  current_rv++;
		      });
    }

    setValuesInWellState(wellState);
    simTimer->setCurrentStepNum(1);
    eclipseWriter->writeTimeStep(*simTimer, *blackoilState, *wellState , false);

    std::shared_ptr<Opm::WellState> wellStateRestored(new Opm::WellState());
    wellStateRestored->init(wells, *blackoilState);

    //Read and verify OPM XWEL data, and solution data: pressure, temperature, saturation data, rs and rv
    std::shared_ptr<Opm::BlackoilState> blackoilStateRestored = createBlackOilState(eclipseState->getEclipseGrid(), phaseUsage);
    Opm::init_from_restart_file(eclipseState, Opm::UgGridHelpers::numCells(*gridManager.c_grid()), phaseUsage, *blackoilStateRestored, *wellStateRestored);
    
    BOOST_CHECK_EQUAL_COLLECTIONS(wellState->bhp().begin(), wellState->bhp().end(), wellStateRestored->bhp().begin(), wellStateRestored->bhp().end());
    BOOST_CHECK_EQUAL_COLLECTIONS(wellState->temperature().begin(), wellState->temperature().end(), wellStateRestored->temperature().begin(), wellStateRestored->temperature().end());
    BOOST_CHECK_EQUAL_COLLECTIONS(wellState->wellRates().begin(), wellState->wellRates().end(), wellStateRestored->wellRates().begin(), wellStateRestored->wellRates().end());
    BOOST_CHECK_EQUAL_COLLECTIONS(wellState->perfRates().begin(), wellState->perfRates().end(), wellStateRestored->perfRates().begin(), wellStateRestored->perfRates().end());
    BOOST_CHECK_EQUAL_COLLECTIONS(wellState->perfPress().begin(), wellState->perfPress().end(), wellStateRestored->perfPress().begin(), wellStateRestored->perfPress().end());

    std::vector<double>& saturation_restored = blackoilStateRestored->getCellData( Opm::BlackoilState::SATURATION );
    std::vector<double> swat_restored;
    std::vector<double> swat;
    std::vector<double> sgas_restored;
    std::vector<double> sgas;
    Opm::EclipseIOUtil::extractFromStripedData(saturation_restored, swat_restored, phaseUsage.phase_pos[Opm::BlackoilPhases::Aqua], phaseUsage.num_phases);
    Opm::EclipseIOUtil::extractFromStripedData(saturation, swat, phaseUsage.phase_pos[Opm::BlackoilPhases::Aqua], phaseUsage.num_phases);
    Opm::EclipseIOUtil::extractFromStripedData(saturation_restored, sgas_restored, phaseUsage.phase_pos[Opm::BlackoilPhases::Vapour], phaseUsage.num_phases);
    Opm::EclipseIOUtil::extractFromStripedData(saturation, sgas, phaseUsage.phase_pos[Opm::BlackoilPhases::Vapour], phaseUsage.num_phases);


    std::vector<double>& pressure_restored = blackoilStateRestored->getCellData( Opm::BlackoilState::PRESSURE );
    std::vector<double>& temperature_restored = blackoilStateRestored->getCellData( Opm::BlackoilState::TEMPERATURE );
    std::vector<double>& rs_restored = blackoilStateRestored->getCellData( Opm::BlackoilState::GASOILRATIO );
    std::vector<double>& rv_restored = blackoilStateRestored->getCellData( Opm::BlackoilState::RV );


    
    for (size_t cellindex = 0; cellindex < 10; ++cellindex) {
	BOOST_CHECK_CLOSE(pressure[cellindex], pressure_restored[cellindex], 0.00001);
        BOOST_CHECK_CLOSE(temperature[cellindex], temperature_restored[cellindex], 0.00001);
	BOOST_CHECK_CLOSE(rs[cellindex], rs_restored[cellindex], 0.0000001);
        BOOST_CHECK_CLOSE(rv[cellindex], rv_restored[cellindex], 0.0000001);

        BOOST_CHECK_CLOSE(swat[cellindex], swat_restored[cellindex], 0.00001);
        BOOST_CHECK_CLOSE(sgas[cellindex], sgas_restored[cellindex], 0.00001);
    }

    test_work_area_free(test_area);
}
