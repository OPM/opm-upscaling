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

#include <algorithm>
#include <iostream>
#include <iomanip>

#include <boost/static_assert.hpp>

#include <dune/common/array.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/Units.hpp>

#include <dune/porsol/common/SimulatorUtilities.hpp>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/CpGrid.hpp>
#include <dune/common/EclipseGridParser.hpp>
#include <dune/common/EclipseGridInspector.hpp>

#include <dune/porsol/common/fortran.hpp>
#include <dune/porsol/common/blas_lapack.hpp>
#include <dune/porsol/common/Matrix.hpp>
#include <dune/porsol/common/GridInterfaceEuler.hpp>
#include <dune/porsol/common/Rock.hpp>
#include <dune/porsol/common/BoundaryConditions.hpp>

#include <dune/porsol/mimetic/TpfaCompressible.hpp>
#include <dune/common/param/ParameterGroup.hpp>
#include <dune/porsol/common/setupGridAndProps.hpp>
#include <dune/porsol/blackoil/fluid/FluidMatrixInteractionBlackoil.hpp>
#include <dune/porsol/blackoil/BlackoilFluid.hpp>


template<int dim, class Grid, class Rock, class Fluid>
void test_flowsolver(const Grid& grid,
                     const Rock& rock,
                     const Fluid& fluid,
                     const double dt,
                     const int num_iter)
{
    // Boundary conditions.
    typedef Dune::BasicBoundaryConditions<true, false>  FBC;
    typedef Dune::FlowBC BC;
    FBC flow_bc(7);
    flow_bc.flowCond(1) = BC(BC::Dirichlet, 300.0*Dune::unit::barsa);
    flow_bc.flowCond(2) = BC(BC::Dirichlet, 100.0*Dune::unit::barsa);

    // Gravity.
    typename Grid::Vector gravity(0.0);
//     gravity[2] = Dune::unit::gravity;

    // Flow solver.
    typedef Dune::TpfaCompressible<Grid, Rock, FBC> FlowSolver;
    FlowSolver solver;
    solver.init(grid, rock, gravity, flow_bc);

    // Source terms.
    std::vector<double> src(grid.numCells(), 0.0);
//     if (g.numberOfCells() > 1) {
//         src[0]     = 1.0;
//         src.back() = -1.0;
//     }


    // Initial state.
    typedef typename Fluid::CompVec CompVec;
    typedef typename Fluid::PhaseVec PhaseVec;
    CompVec init_z(0.0);
    init_z[0] = 1.0;
    std::vector<CompVec> z(grid.numCells(), init_z);
    MESSAGE("******* Assuming zero capillary pressures *******");
    PhaseVec init_p(100.0*Dune::unit::barsa);
    std::vector<PhaseVec> phase_pressure(grid.numCells(), init_p);
    // Rescale z values so that pore volume is filled exactly
    // (to get zero initial volume discrepancy).
    for (int cell = 0; cell < grid.numCells(); ++cell) {
        double pore_vol = grid.cellVolume(cell)*rock.porosity(cell);
        typename Fluid::FluidState state = fluid.computeState(phase_pressure[cell], z[cell]);
        double fluid_vol = state.total_phase_volume_;
        z[cell] *= pore_vol/fluid_vol;
    }

    // Solve flow system.
    solver.solve(fluid, phase_pressure, z, flow_bc, src, dt, num_iter, 1e-8, 3, 1);

    // Output to VTK.
    typedef typename FlowSolver::SolutionType FlowSolution;
    FlowSolution soln = solver.getSolution();
    std::vector<typename Grid::Vector> cell_velocity;
    estimateCellVelocitySimpleInterface(cell_velocity, grid, soln);
    // Dune's vtk writer wants multi-component data to be flattened.
    std::vector<double> cell_velocity_flat(&*cell_velocity.front().begin(),
                                           &*cell_velocity.back().end());
//     getCellPressure(cell_pressure, grid, soln);
//    output_pressure = soln.cellPressure();
    Dune::VTKWriter<typename Grid::LeafGridView> vtkwriter(grid.leafView());
    vtkwriter.addCellData(soln.cellPressure(), "pressure");
    vtkwriter.addCellData(cell_velocity_flat, "velocity", dim);
    vtkwriter.write("testsolution-" + boost::lexical_cast<std::string>(0),
                    Dune::VTKOptions::ascii);

    // Dump pressures to Matlab.
    std::ofstream dump("pressure");
    std::copy(soln.cellPressure().begin(), soln.cellPressure().end(),
              std::ostream_iterator<double>(dump, " "));
}


using namespace Dune;

int main(int argc, char** argv)
{
    Dune::parameter::ParameterGroup param(argc, argv);
    Dune::MPIHelper::instance(argc,argv);

    // Make a grid and props.
    Dune::CpGrid grid;
    Dune::Rock<3> rock;
    Opm::BlackoilFluid fluid;

    // Initialization.
    std::string fileformat = param.getDefault<std::string>("fileformat", "cartesian");
    if (fileformat == "eclipse") {
        EclipseGridParser parser(param.get<std::string>("filename"));
        double z_tolerance = param.getDefault<double>("z_tolerance", 0.0);
        bool periodic_extension = param.getDefault<bool>("periodic_extension", false);
        bool turn_normals = param.getDefault<bool>("turn_normals", false);
        grid.processEclipseFormat(parser, z_tolerance, periodic_extension, turn_normals);
        double perm_threshold_md = param.getDefault("perm_threshold_md", 0.0);
        double perm_threshold = unit::convert::from(perm_threshold_md, prefix::milli*unit::darcy);
        rock.init(parser, grid.globalCell(), perm_threshold);
        fluid.init(parser);
    } else if (fileformat == "cartesian") {
        array<int, 3> dims = {{ param.getDefault<int>("nx", 1),
                                param.getDefault<int>("ny", 1),
                                param.getDefault<int>("nz", 1) }};
        array<double, 3> cellsz = {{ param.getDefault<double>("dx", 1.0),
                                     param.getDefault<double>("dy", 1.0),
                                     param.getDefault<double>("dz", 1.0) }};
        grid.createCartesian(dims, cellsz);
        double default_poro = param.getDefault("default_poro", 0.2);
        double default_perm_md = param.getDefault("default_perm_md", 100.0);
        double default_perm = unit::convert::from(default_perm_md, prefix::milli*unit::darcy);
        MESSAGE("Warning: For generated cartesian grids, we use uniform rock properties.");
        rock.init(grid.size(0), default_poro, default_perm);
        EclipseGridParser parser(param.get<std::string>("filename")); // Need a parser for the fluids anyway.
        fluid.init(parser);
    } else {
        THROW("Unknown file format string: " << fileformat);
    }

    double dt = param.getDefault("dt", 1.0);
    int num_iter = param.getDefault("num_iter", 5);

    // Run test.
    test_flowsolver<3>(grid, rock, fluid, dt, num_iter);
}

