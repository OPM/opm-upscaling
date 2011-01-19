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



#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <dune/porsol/blackoil/BlackoilSimulator.hpp>
#include <dune/common/mpihelper.hh>
#include <dune/porsol/common/SimulatorUtilities.hpp>
#include <dune/grid/CpGrid.hpp>
#include <dune/porsol/common/Rock.hpp>
#include <dune/porsol/mimetic/TpfaCompressible.hpp>
#include <dune/common/StopWatch.hpp>
#include <dune/porsol/common/Wells.hpp>
#include <dune/porsol/blackoil/BlackoilFluid.hpp>
#include <dune/porsol/blackoil/ComponentTransport.hpp>


typedef Dune::CpGrid Grid;
typedef Dune::Rock<Grid::dimension> Rock;
typedef Opm::BlackoilFluid Fluid;
// typedef Opm::BlackoilWells Wells; // WELLS
typedef Opm::Wells Wells;
typedef Dune::BasicBoundaryConditions<true, false>  FBC;
typedef Dune::TpfaCompressible<Grid, Rock, Fluid, FBC> FlowSolver;
typedef Opm::ExplicitCompositionalTransport<Grid, Rock, Fluid, Wells> TransportSolver;


typedef Opm::BlackoilSimulator<Grid, Rock, Fluid, Wells, FlowSolver, TransportSolver> Simulator;


int main(int argc, char** argv)
{
    Dune::parameter::ParameterGroup param(argc, argv);
    Dune::MPIHelper::instance(argc,argv);

    // Make a grid and props.
    Grid grid;
    Rock rock;
    Fluid fluid;
    Wells wells;
    FlowSolver flow_solver;
    TransportSolver transport_solver;

    using namespace Dune;

    // Initialization.
    std::string fileformat = param.getDefault<std::string>("fileformat", "cartesian");
    if (fileformat == "eclipse") {
        Dune::EclipseGridParser parser(param.get<std::string>("filename"));
        double z_tolerance = param.getDefault<double>("z_tolerance", 0.0);
        bool periodic_extension = param.getDefault<bool>("periodic_extension", false);
        bool turn_normals = param.getDefault<bool>("turn_normals", false);
        grid.processEclipseFormat(parser, z_tolerance, periodic_extension, turn_normals);
        double perm_threshold_md = param.getDefault("perm_threshold_md", 0.0);
        double perm_threshold = Dune::unit::convert::from(perm_threshold_md, Dune::prefix::milli*Dune::unit::darcy);
        rock.init(parser, grid.globalCell(), perm_threshold);
        fluid.init(parser);
        wells.init(parser);
    } else if (fileformat == "cartesian") {
        Dune::array<int, 3> dims = {{ param.getDefault<int>("nx", 1),
                                      param.getDefault<int>("ny", 1),
                                      param.getDefault<int>("nz", 1) }};
        Dune::array<double, 3> cellsz = {{ param.getDefault<double>("dx", 1.0),
                                           param.getDefault<double>("dy", 1.0),
                                           param.getDefault<double>("dz", 1.0) }};
        grid.createCartesian(dims, cellsz);
        double default_poro = param.getDefault("default_poro", 1.0);
        double default_perm_md = param.getDefault("default_perm_md", 100.0);
        double default_perm = unit::convert::from(default_perm_md, prefix::milli*unit::darcy);
        MESSAGE("Warning: For generated cartesian grids, we use uniform rock properties.");
        rock.init(grid.size(0), default_poro, default_perm);
        EclipseGridParser parser(param.get<std::string>("filename")); // Need a parser for the fluids anyway.
        fluid.init(parser);
        wells.init(parser);
    } else {
        THROW("Unknown file format string: " << fileformat);
    }
    flow_solver.init(param);
    transport_solver.init(param);
    double total_time = param.getDefault("total_time", 30*unit::day);
    double initial_stepsize = param.getDefault("initial_stepsize", 1.0*unit::day);
    bool do_impes = param.getDefault("do_impes", false);
    std::string output_dir = param.getDefault<std::string>("output_dir", "output");
    bool gravity_test = param.getDefault("gravity_test", false);
    bool newcode = param.getDefault("newcode", true);

    // Run simulation.
    Dune::time::StopWatch clock;
    clock.start();
    Simulator sim;
    sim.simulate(grid, rock, fluid, wells, flow_solver, transport_solver,
                 total_time, initial_stepsize, do_impes, output_dir, gravity_test, newcode);
    clock.stop();
    std::cout << "\n\nSimulation clock time (secs): " << clock.secsSinceStart() << std::endl;
}

