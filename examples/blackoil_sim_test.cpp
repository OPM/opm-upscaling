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
#include <dune/common/StopWatch.hpp>
#include <dune/porsol/common/setupGridAndProps.hpp>
#include <dune/porsol/common/Wells.hpp>
#include <dune/porsol/blackoil/BlackoilWells.hpp>
#include <dune/porsol/blackoil/fluid/FluidMatrixInteractionBlackoil.hpp>
#include <dune/porsol/blackoil/BlackoilFluid.hpp>

#include <dune/porsol/blackoil/ComponentTransport.hpp>

#include <boost/filesystem/convenience.hpp>


template<class Grid, class Fluid>
void output(const Grid& grid,
            const std::vector<typename Fluid::PhaseVec>& cell_pressure,
            const std::vector<typename Fluid::CompVec>& z,
            const std::vector<double>& face_flux,
            const int step,
            const std::string& filebase)
{
    // Ensure directory exists.
    boost::filesystem::path fpath(filebase);
    if (fpath.has_branch_path()) {
        create_directories(fpath.branch_path());
    }

    // Output to VTK.
    std::vector<typename Grid::Vector> cell_velocity;
    estimateCellVelocitySimpleInterface(cell_velocity, grid, face_flux);
    // Dune's vtk writer wants multi-component data to be flattened.
    std::vector<double> cell_pressure_flat(&*cell_pressure.front().begin(),
                                           &*cell_pressure.back().end());
    std::vector<double> cell_velocity_flat(&*cell_velocity.front().begin(),
                                           &*cell_velocity.back().end());
    std::vector<double> z_flat(&*z.front().begin(),
                               &*z.back().end());
    Dune::VTKWriter<typename Grid::LeafGridView> vtkwriter(grid.leafView());
    vtkwriter.addCellData(cell_pressure_flat, "pressure", Fluid::numPhases);
    vtkwriter.addCellData(cell_velocity_flat, "velocity", Grid::dimension);
    vtkwriter.addCellData(z_flat, "z", Fluid::numComponents);
    vtkwriter.write(filebase + '-' + boost::lexical_cast<std::string>(step),
                    Dune::VTKOptions::ascii);

    // Dump data for Matlab.
    std::vector<double> zv[Fluid::numComponents];
    for (int comp = 0; comp < Fluid::numComponents; ++comp) {
        zv[comp].resize(grid.numCells());
        for (int cell = 0; cell < grid.numCells(); ++cell) {
            zv[comp][cell] = z[cell][comp];
        }
    }
    std::string matlabdumpname(filebase + "-");
    matlabdumpname += boost::lexical_cast<std::string>(step);
    matlabdumpname += ".dat";
    std::ofstream dump(matlabdumpname.c_str());
    dump.precision(15);
    int num_cells = cell_pressure.size();
    std::vector<double> liq_press(num_cells);
    for (int cell = 0; cell < num_cells; ++cell) {
        liq_press[cell] = cell_pressure[cell][Fluid::Liquid];
    }
    std::copy(liq_press.begin(), liq_press.end(),
              std::ostream_iterator<double>(dump, " "));
    dump << '\n';
    for (int comp = 0; comp < Fluid::numComponents; ++comp) {
        std::copy(zv[comp].begin(), zv[comp].end(),
                  std::ostream_iterator<double>(dump, " "));
        dump << '\n';
    }
}




template<class Grid, class Rock, class Fluid, class Wells, class FlowSolver, class TransportSolver>
void simulate(const Grid& grid,
              const Rock& rock,
              const Fluid& fluid,
              Wells& wells,  // Will change after pressure updates.
              FlowSolver& flow_solver,
              TransportSolver& transport_solver,
              const double total_time,
              const double initial_stepsize,
              const bool do_impes,
              const std::string& output_dir)
{
    // Boundary conditions.
    typedef Dune::FlowBC BC;
    typedef Dune::BasicBoundaryConditions<true, false>  FBC;
    FBC flow_bc(7);
    flow_bc.flowCond(1) = BC(BC::Dirichlet, 300.0*Dune::unit::barsa);
    flow_bc.flowCond(2) = BC(BC::Dirichlet, 100.0*Dune::unit::barsa); // WELLS

    // Gravity.
    typename Grid::Vector gravity(0.0);
//     gravity[2] = Dune::unit::gravity;

    // Flow solver setup.
    flow_solver.setup(grid, rock, fluid, wells, gravity, flow_bc);

    // Transport solver setup.
    transport_solver.setup(grid, rock, fluid, wells);

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
    init_z[Fluid::Oil] = 1.0;
    CompVec bdy_z = flow_solver.inflowMixture();
    std::vector<CompVec> cell_z(grid.numCells(), init_z);
    MESSAGE("******* Assuming zero capillary pressures *******");
    PhaseVec init_p(100.0*Dune::unit::barsa);
    std::vector<PhaseVec> cell_pressure(grid.numCells(), init_p);
    PhaseVec bdy_p(300.0*Dune::unit::barsa);
    // PhaseVec bdy_p(100.0*Dune::unit::barsa); // WELLS
    // Rescale z values so that pore volume is filled exactly
    // (to get zero initial volume discrepancy).
    for (int cell = 0; cell < grid.numCells(); ++cell) {
        double pore_vol = grid.cellVolume(cell)*rock.porosity(cell);
        typename Fluid::FluidState state = fluid.computeState(cell_pressure[cell], cell_z[cell]);
        double fluid_vol = state.total_phase_volume_;
        cell_z[cell] *= pore_vol/fluid_vol;
    }
    int num_faces = grid.numFaces();
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
    double voldisclimit = flow_solver.volumeDiscrepancyLimit();
    double stepsize = initial_stepsize;
    double current_time = 0.0;
    int step = -1;
    std::vector<double> face_flux;
    std::vector<double> well_pressure;
    std::vector<double> well_flux;
    std::vector<PhaseVec> cell_pressure_start;
    std::vector<PhaseVec> face_pressure_start;
    std::vector<CompVec> cell_z_start;
    while (current_time < total_time) {
        cell_pressure_start = cell_pressure;
        face_pressure_start = face_pressure;
        cell_z_start = cell_z;

        // Do not run past total_time.
        if (current_time + stepsize > total_time) {
            stepsize = total_time - current_time;
        }
        ++step;
        std::cout << "\n\n================    Simulation step number " << step
                  << "    ==============="
                  << "\n      Current time (days)     " << Dune::unit::convert::to(current_time, Dune::unit::day)
                  << "\n      Current stepsize (days) " << Dune::unit::convert::to(stepsize, Dune::unit::day)
                  << "\n      Total time (days)       " << Dune::unit::convert::to(total_time, Dune::unit::day)
                  << "\n" << std::endl;

        // Solve flow system.
        enum FlowSolver::ReturnCode result
            = flow_solver.solve(cell_pressure, face_pressure, cell_z, face_flux, well_pressure, well_flux, src, stepsize, do_impes);

        // Check if the flow solver succeeded.
        if (result != FlowSolver::SolveOk) {
            THROW("Flow solver refused to run due to too large volume discrepancy.");
        }

        // Update wells with new perforation pressures and fluxes.
        wells.update(grid.numCells(), well_pressure, well_flux);

        // Transport and check volume discrepancy.
        bool voldisc_ok = true;
        if (!do_impes) {
            double actual_computed_time
                = transport_solver.transport(bdy_p, bdy_z,
                                             face_flux, cell_pressure, face_pressure,
                                             stepsize, voldisclimit, cell_z);
            voldisc_ok = (actual_computed_time == stepsize);
        } else {
            voldisc_ok = flow_solver.volumeDiscrepancyAcceptable(cell_pressure, face_pressure, cell_z, stepsize);
        }

        // If discrepancy too large, redo entire pressure step.
        if (!voldisc_ok) {
            std::cout << "********* Shortening (pressure) stepsize, redoing step number " << step <<" **********" << std::endl;
            stepsize *= 0.5;
            --step;
            cell_pressure = cell_pressure_start;
            face_pressure = face_pressure_start;
            cell_z = cell_z_start;
            continue;
        }

        // Output.
        std::string output_name = output_dir + "/" + "blackoil-output";
        output<Grid, Fluid>(grid, cell_pressure, cell_z, face_flux, step, output_name);

        // Adjust time.
        current_time += stepsize;
        if (voldisc_ok) {
            // stepsize *= 1.5;
        }
    }
}


typedef Dune::CpGrid Grid;
typedef Dune::Rock<Grid::dimension> Rock;
typedef Opm::BlackoilFluid Fluid;
// typedef Opm::BlackoilWells Wells; // WELLS
typedef Opm::Wells Wells;
typedef Dune::BasicBoundaryConditions<true, false>  FBC;
typedef Dune::TpfaCompressible<Grid, Rock, Fluid, FBC> FlowSolver;
typedef Opm::ExplicitCompositionalTransport<Grid, Rock, Fluid, Wells> TransportSolver;

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

    // Run simulation.
    Dune::time::StopWatch clock;
    clock.start();
    simulate(grid, rock, fluid, wells, flow_solver, transport_solver, total_time, initial_stepsize, do_impes, output_dir);
    clock.stop();
    std::cout << "\n\nSimulation clock time (secs): " << clock.secsSinceStart() << std::endl;
}

