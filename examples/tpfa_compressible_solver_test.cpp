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

#include <opm/core/utility/have_boost_redef.hpp>

#include <algorithm>
#include <iostream>
#include <iomanip>

#include <boost/static_assert.hpp>

#include <dune/common/array.hh>
#include <dune/common/mpihelper.hh>
#include <opm/core/utility/Units.hpp>

#include <dune/porsol/common/SimulatorUtilities.hpp>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/CpGrid.hpp>
#include <opm/core/eclipse/EclipseGridParser.hpp>
#include <opm/core/eclipse/EclipseGridInspector.hpp>

#include <dune/porsol/common/fortran.hpp>
#include <dune/porsol/common/blas_lapack.hpp>
#include <dune/porsol/common/Matrix.hpp>
#include <dune/porsol/common/GridInterfaceEuler.hpp>
#include <dune/porsol/common/Rock.hpp>
#include <dune/porsol/common/BoundaryConditions.hpp>

#include <dune/porsol/blackoil/fluid/BlackoilPVT.hpp>
#include <dune/porsol/blackoil/BlackoilFluid.hpp>

#include <dune/porsol/mimetic/TpfaCompressible.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <dune/porsol/common/setupGridAndProps.hpp>
#include <dune/porsol/common/Wells.hpp>
#include <dune/porsol/blackoil/fluid/FluidMatrixInteractionBlackoil.hpp>


template<int dim, class Grid, class Rock, class Fluid, class FlowSolver>
void test_flowsolver(const Grid& grid,
                     const Rock& rock,
                     const Fluid& fluid,
                     FlowSolver& solver,
                     const double dt)
{
    // Boundary conditions.
    typedef Dune::FlowBC BC;
    typedef Dune::BasicBoundaryConditions<true, false>  FBC;
    FBC flow_bc(7);
    flow_bc.flowCond(1) = BC(BC::Dirichlet, 300.0*Opm::unit::barsa);
    flow_bc.flowCond(2) = BC(BC::Dirichlet, 100.0*Opm::unit::barsa);

    // Gravity.
    typename Grid::Vector gravity(0.0);
//     gravity[2] = Dune::unit::gravity;

    Opm::Wells wells;

    // Flow solver setup.
    solver.setup(grid, rock, fluid, wells, gravity, flow_bc);

    // Source terms.
    std::vector<double> src(grid.numCells(), 0.0);
//     if (g.numberOfCells() > 1) {
//         src[0]     = 1.0;
//         src.back() = -1.0;
//     }

    int num_cells = grid.numCells();
    int num_faces = grid.numFaces();

    // Initial state.
    typedef typename Fluid::CompVec CompVec;
    typedef typename Fluid::PhaseVec PhaseVec;
    CompVec init_z(0.0);
    init_z[Fluid::Oil] = 1.0;
    std::vector<CompVec> z(grid.numCells(), init_z);
    MESSAGE("******* Assuming zero capillary pressures *******");
    PhaseVec init_p(100.0*Opm::unit::barsa);
    std::vector<PhaseVec> cell_pressure(grid.numCells(), init_p);
    // Rescale z values so that pore volume is filled exactly
    // (to get zero initial volume discrepancy).
    for (int cell = 0; cell < grid.numCells(); ++cell) {
        typename Fluid::FluidState state = fluid.computeState(cell_pressure[cell], z[cell]);
        double fluid_vol = state.total_phase_volume_density_;
        z[cell] *= 1.0/fluid_vol;
    }
    std::vector<PhaseVec> face_pressure(num_faces);
    for (int face = 0; face < num_faces; ++face) {
        int bid = grid.boundaryId(face);
        if (flow_bc.flowCond(bid).isDirichlet()) {
            face_pressure[face] = flow_bc.flowCond(bid).pressure();
        } else {
            int c[2] = { grid.faceCell(face, 0), grid.faceCell(face, 1) };
            face_pressure[face] = 0.0;
            int num = 0;
            for (int j = 0; j < 2; ++j) {
                if (c[j] >= 0) {
                    face_pressure[face] += cell_pressure[c[j]];
                    ++num;
                }
            }
            face_pressure[face] /= double(num);
        }
    }
    std::vector<double> face_flux, well_bhp_pressure, well_perf_pressure, well_flux;

    // Solve flow system.
    solver.solve(cell_pressure, face_pressure, z, face_flux, well_bhp_pressure, well_perf_pressure, well_flux, src, dt);

    // Output to VTK.
    std::vector<typename Grid::Vector> cell_velocity;
    estimateCellVelocitySimpleInterface(cell_velocity, grid, face_flux);
    // Dune's vtk writer wants multi-component data to be flattened.
    std::vector<double> cell_pressure_flat(&*cell_pressure.front().begin(),
                                           &*cell_pressure.back().end());
    std::vector<double> cell_velocity_flat(&*cell_velocity.front().begin(),
                                           &*cell_velocity.back().end());
    Dune::VTKWriter<typename Grid::LeafGridView> vtkwriter(grid.leafView());
    vtkwriter.addCellData(cell_pressure_flat, "pressure", Fluid::numPhases);
    vtkwriter.addCellData(cell_velocity_flat, "velocity", Grid::dimension);
    vtkwriter.write("testsolution", Dune::VTKOptions::ascii);

    // Dump data for Matlab.
    std::ofstream dump("celldump");
    dump.precision(15);
    std::vector<double> liq_press(num_cells);
    for (int cell = 0; cell < num_cells; ++cell) {
        liq_press[cell] = cell_pressure[cell][Fluid::Liquid];
    }
    std::copy(liq_press.begin(), liq_press.end(),
              std::ostream_iterator<double>(dump, " "));
    dump << '\n';
}


typedef Dune::CpGrid Grid;
typedef Dune::Rock<Grid::dimension> Rock;
typedef Opm::BlackoilFluid Fluid;
typedef Dune::BasicBoundaryConditions<true, false>  FBC;
typedef Dune::TpfaCompressible<Grid, Rock, Fluid, Opm::Wells, FBC> FlowSolver;

int main(int argc, char** argv)
{
    Opm::parameter::ParameterGroup param(argc, argv);
    Dune::MPIHelper::instance(argc,argv);

    // Make a grid and props.
    Grid grid;
    Rock rock;
    Fluid fluid;
    FlowSolver solver;

    using namespace Dune;

    // Initialization.
    std::string fileformat = param.getDefault<std::string>("fileformat", "cartesian");
    if (fileformat == "eclipse") {
        Opm::EclipseGridParser parser(param.get<std::string>("filename"));
        double z_tolerance = param.getDefault<double>("z_tolerance", 0.0);
        bool periodic_extension = param.getDefault<bool>("periodic_extension", false);
        bool turn_normals = param.getDefault<bool>("turn_normals", false);
        grid.processEclipseFormat(parser, z_tolerance, periodic_extension, turn_normals);
        double perm_threshold_md = param.getDefault("perm_threshold_md", 0.0);
        double perm_threshold = Opm::unit::convert::from(perm_threshold_md, Opm::prefix::milli*Opm::unit::darcy);
        rock.init(parser, grid.globalCell(), perm_threshold);
        fluid.init(parser);
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
        double default_perm = Opm::unit::convert::from(default_perm_md, Opm::prefix::milli*Opm::unit::darcy);
        MESSAGE("Warning: For generated cartesian grids, we use uniform rock properties.");
        rock.init(grid.size(0), default_poro, default_perm);
	Opm::EclipseGridParser parser(param.get<std::string>("filename")); // Need a parser for the fluids anyway.
        fluid.init(parser);
    } else {
        THROW("Unknown file format string: " << fileformat);
    }
    solver.init(param);
    double dt = param.getDefault("dt", 1.0);

    // Run test.
    test_flowsolver<3>(grid, rock, fluid, solver, dt);
}
