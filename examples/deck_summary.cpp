/*
  Copyright 2014 Statoil ASA

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


#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well.hpp>

#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>

#include <iostream>
#include <memory>
#include <vector>


void summarizeSchedule(const std::shared_ptr<const Opm::EclipseState> eclipseState , const Opm::GridManager& gridManager) {
    Opm::SimulatorTimer simTimer;
    std::shared_ptr<const Opm::Schedule> schedule = eclipseState->getSchedule();
    simTimer.init(  schedule->getTimeMap() );
    
    for (size_t timeStep = 0; timeStep < simTimer.numSteps(); timeStep++) {
        Opm::WellsManager wells( eclipseState , timeStep , *gridManager.c_grid() , NULL /* Perm */ ); 
    }
    std::cout << "Timesteps: " << simTimer.numSteps() << std::endl;
    
    {
        size_t wellIndex;
        std::vector<Opm::WellConstPtr> wells = schedule->getWells();
        
        std::cout << "Wells: ";
        for (wellIndex = 0; wellIndex  < schedule->numWells(); wellIndex++) {
            Opm::WellConstPtr well = wells[wellIndex];
            std::cout << well->name() << " ";
            if ((wellIndex % 8) == 7)
                std::cout << std::endl << "       ";
            
        }
        std::cout << std::endl;
    }
}




int main(int argc, char ** argv) {
    Opm::ParserPtr parser(new Opm::Parser());
    std::string file = argv[1];
    std::string gridFile = argv[2];
    Opm::DeckConstPtr deck = parser->parseFile(file, false);
    std::shared_ptr<const Opm::EclipseState> state( new Opm::EclipseState( deck) );

    Opm::EclipseGridParser oldDeck(gridFile);
    Opm::GridManager gridManager(oldDeck);

    summarizeSchedule( state , gridManager );
    
    return 0;
}
