//===========================================================================
//
// File: mimetic_solver_test.cpp
//
// Created: Tue Jul  7 15:35:21 2009
//
// Author(s): B�rd Skaflestad     <bard.skaflestad@sintef.no>
//            Atgeirr F Rasmussen <atgeirr@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

  This file is part of The Open Reservoir Simulator Project (OpenRS).

  OpenRS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OpenRS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "config.h"

#include <opm/common/utility/platform_dependent/disable_warnings.h>

#include <dune/common/version.hh>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/grid/yaspgrid.hh>

#if HAVE_ALUGRID
#include <dune/common/shared_ptr.hh>
#include <dune/grid/io/file/gmshreader.hh>
// dune-grid 2.2.0 tests for this define instead of HAVE_ALUGRID
#define ENABLE_ALUGRID 1
#include <dune/grid/alugrid.hh>
#endif

#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/parser/eclipse/Units/Units.hpp>

#include <dune/grid/CpGrid.hpp>

#include <opm/porsol/common/BoundaryConditions.hpp>
#include <opm/porsol/common/GridInterfaceEuler.hpp>
#include <opm/porsol/common/Matrix.hpp>
#include <opm/porsol/common/ReservoirPropertyCapillary.hpp>
#include <opm/porsol/common/SimulatorUtilities.hpp>
#include <opm/porsol/common/blas_lapack.hpp>
#include <opm/porsol/common/fortran.hpp>
#include <opm/porsol/common/setupGridAndProps.hpp>
#include <opm/porsol/mimetic/IncompFlowSolverHybrid.hpp>
#include <opm/porsol/mimetic/MimeticIPEvaluator.hpp>

#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>

#if USE_ALUGRID
template<class GType>
Dune::shared_ptr<GType>
make_gmsh(const std::string& msh_file)
{
    return Dune::shared_ptr<GType>(Dune::GmshReader<GType>::read(msh_file));
}
#endif



template<int dim, class GI, class RI>
void test_flowsolver(const GI& g, const RI& r)
{
    typedef typename GI::CellIterator                   CI;
    typedef Opm::BasicBoundaryConditions<true, false>  FBC;
    typedef Opm::IncompFlowSolverHybrid<GI, RI, FBC,
                                         Opm::MimeticIPEvaluator> FlowSolver;

    FlowSolver solver;

    typedef Opm::FlowBC BC;
    FBC flow_bc(7);

#if !USE_ALUGRID
    flow_bc.flowCond(5) = BC(BC::Dirichlet, 100.0*Opm::unit::barsa);
    flow_bc.flowCond(6) = BC(BC::Dirichlet, 0.0*Opm::unit::barsa);
#endif

    typename CI::Vector gravity(0.0);
    // gravity[2] = Dune::unit::gravity;

    solver.init(g, r, gravity, flow_bc);

    std::vector<double> src(g.numberOfCells(), 0.0);
    std::vector<double> sat(g.numberOfCells(), 0.0);
//     if (g.numberOfCells() > 1) {
//         src[0]     = 1.0;
//         src.back() = -1.0;
//     }

    solver.solve(r, sat, flow_bc, src, 5e-9, 3, 1);

#if 1
    typedef typename FlowSolver::SolutionType FlowSolution;
    FlowSolution soln = solver.getSolution();

    std::vector<typename GI::Vector> cell_velocity;
    estimateCellVelocity(cell_velocity, g, soln);
    // Dune's vtk writer wants multi-component data to be flattened.
    std::vector<double> cell_velocity_flat(&*cell_velocity.front().begin(),
                                           &*cell_velocity.back().end());
    std::vector<double> cell_pressure;
    getCellPressure(cell_pressure, g, soln);

    Dune::VTKWriter<typename GI::GridType::LeafGridView> vtkwriter(g.grid().leafGridView());
    vtkwriter.addCellData(cell_velocity_flat, "velocity", dim);
    vtkwriter.addCellData(cell_pressure, "pressure");
    vtkwriter.write("testsolution-" + boost::lexical_cast<std::string>(0),
                    Dune::VTK::ascii);
#else    
    solver.printSystem("system");
    typedef typename FlowSolver::SolutionType FlowSolution;
    FlowSolution soln = solver.getSolution();

    std::cout << "Cell Pressure:\n" << std::scientific << std::setprecision(15);
    for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
        std::cout << '\t' << soln.pressure(c) << '\n';
    }

    std::cout << "Cell (Out) Fluxes:\n";
    std::cout << "flux = [\n";
    for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
        for (FI f = c->facebegin(); f != c->faceend(); ++f) {
            std::cout << soln.outflux(f) << ' ';
        }
        std::cout << "\b\n";
    }
    std::cout << "]\n";
#endif
}


using namespace Opm;

int main(int argc, char** argv)
try
{
    Opm::ParameterGroup param(argc, argv);
    Dune::MPIHelper::instance(argc,argv);

    // Make a grid
#if !USE_ALUGRID
    // Make a grid and props.
    Dune::CpGrid grid;
    ReservoirPropertyCapillary<3> res_prop;
    setupGridAndProps(param, grid, res_prop);

    // Make the grid interface
    Opm::GridInterfaceEuler<Dune::CpGrid> g(grid);
#else
    typedef Dune::ALUSimplexGrid<3,3> GType;
    Dune::shared_ptr<GType> pgrid = make_gmsh<GType>(param.get<std::string>("filename"));
    Opm::GridInterfaceEuler<GType> g(*pgrid);

    ReservoirPropertyCapillary<3> res_prop;
    res_prop.init(g.numberOfCells());
#endif

    test_flowsolver<3>(g, res_prop);
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}
