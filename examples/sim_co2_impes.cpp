/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.

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

#include <opm/common/utility/platform_dependent/disable_warnings.h>

#include <dune/common/version.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/porsol/common/SimulatorUtilities.hpp>
#include <opm/porsol/common/Rock.hpp>
#include <opm/porsol/blackoil/co2fluid/BlackoilCo2PVT.hpp>
#include <opm/porsol/blackoil/BlackoilSimulator.hpp>
#include <opm/porsol/blackoil/BlackoilFluid.hpp>
#include <opm/porsol/blackoil/BlackoilWells.hpp>
#include <opm/porsol/mimetic/TpfaCompressible.hpp>
#include <opm/porsol/blackoil/ComponentTransport.hpp>

#include <opm/core/utility/StopWatch.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#include <dune/grid/CpGrid.hpp>

#include <iostream>

typedef Dune::CpGrid Grid;
typedef Opm::Rock<Grid::dimension> Rock;
typedef Opm::BlackoilFluid Fluid;
typedef Opm::BlackoilWells Wells;
typedef Opm::BasicBoundaryConditions<true, false>  FBC;
typedef Opm::TpfaCompressible<Grid, Rock, Fluid, Wells, FBC> FlowSolver;
typedef Opm::ExplicitCompositionalTransport<Grid, Rock, Fluid, Wells> TransportSolver;

typedef Opm::BlackoilSimulator<Grid, Rock, Fluid, Wells, FlowSolver, TransportSolver> Simulator;


int main(int argc, char** argv)
try
{
    Opm::ParameterGroup param(argc, argv);
    Dune::MPIHelper::instance(argc,argv);

    // Initialize.
    Simulator sim;
    sim.init(param);

    // Run simulation.
    Opm::time::StopWatch clock;
    clock.start();
    sim.simulate();
    clock.stop();
    std::cout << "\n\nSimulation clock time (secs): " << clock.secsSinceStart() << std::endl;
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}


