//===========================================================================
//
// File: known_answer_test.cpp
//
// Created: Thu Mar 25 13:57:12 2010
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Jostein R Natvig    <jostein.r.natvig@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2010 Statoil ASA.

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

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2,3)
#include <dune/common/parallel/mpihelper.hh>
#else
#include <dune/common/mpihelper.hh>
#endif

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/yaspgrid.hh>

#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/utility/Units.hpp>

#include <dune/grid/CpGrid.hpp>

#include <opm/porsol/common/BoundaryConditions.hpp>
#include <opm/porsol/common/GridInterfaceEuler.hpp>
#include <opm/porsol/common/Matrix.hpp>
#include <opm/porsol/common/ReservoirPropertyCapillary.hpp>
#include <opm/porsol/common/SimulatorUtilities.hpp>
#include <opm/porsol/common/blas_lapack.hpp>
#include <opm/porsol/common/fortran.hpp>
#include <opm/porsol/mimetic/IncompFlowSolverHybrid.hpp>
#include <opm/porsol/mimetic/MimeticIPEvaluator.hpp>

#include <opm/upscaling/initCPGrid.hpp>



#include <algorithm>
#include <array>
#include <cassert>
#include <iomanip>
#include <iostream>

#include <boost/lexical_cast.hpp>



// ------------ Specifying the solution ------------

typedef Dune::FieldVector<double, 3> Vec;

namespace {

    double u(const Vec& x)
    {
        const double pi = 3.14159265358979323846264338327950288;
        return std::sin(2*pi*x[0]) * std::cos(2*pi*x[1]) * x[2];
    }

    double Lu(const Vec& x)
    {
        const double pi = 3.14159265358979323846264338327950288;
        return -2 * 2*pi * 2*pi * std::sin(2*pi*x[0]) * std::cos(2*pi*x[1]) * x[2];
    }

} // namespace anonymous

namespace Opm
{
    template <class BoundaryFunc>
    class FunctionBoundaryConditions : public PeriodicConditionHandler
    {
    public:
        FunctionBoundaryConditions(BoundaryFunc bfunc)
            : bfunc_(bfunc)
        {
        }

        template <class BoundaryFace>
        FlowBC flowCond(const BoundaryFace& bf) const
        {
            assert(bf.boundary());
            return FlowBC(FlowBC::Dirichlet, bfunc_(bf.centroid()));
        }

    private:
        BoundaryFunc bfunc_;
    };

}

template<class GI>
void assign_src(const GI& g, std::vector<double>& src)
{
    typedef typename GI::CellIterator CI;

    int count = 0;
    for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
        src[count++] = -Lu(c->centroid()) * c->volume();
    }
}

template<class GI, class BCS>
void assign_bc(const GI& g, BCS& bcs)
{
    typedef Opm::FlowBC BC;
    typedef typename GI::CellIterator CI;
    typedef typename CI::FaceIterator FI;

    int max_bid = 0;
    for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
        for (FI f = c->facebegin(); f != c->faceend(); ++f) {
            int bid = f->boundaryId();

            if (bid > max_bid) {
                max_bid = bid;
                bcs.resize(bid + 1);
            }

            bcs.flowCond(bid) = BC(BC::Dirichlet, u(f->centroid()));
        }
    }
}

template<class GI>
void compare_pressure(const GI& g, const std::vector<double>& p)
{
    typedef typename GI::CellIterator CI;
    int count = 0;
    double l1err = 0.0;
    double l2err = 0.0;
    double linferr = 0.0;
    double totv = 0.0;

    for (CI c = g.cellbegin(); c != g.cellend(); ++c, ++count) {
        Vec    cen  = c->centroid();
        double uval = u(cen);
        double diff = uval - p[count];
        double v    = c->volume();

        l1err   += std::fabs(diff * v);
        l2err   += diff * diff * v;
        linferr  = std::max(std::fabs(diff), linferr);

        totv += v;
    }

    l2err = std::sqrt(l2err);
    std::cout << "\n\n"
              << "\n     L1 error density: " << l1err / totv
              << "\n     L2 error density: " << l2err / totv
              << "\n     Linf error:       " << linferr << "\n\n\n";
}


template<class GI, class RI>
void test_flowsolver(const GI& g, const RI& r, double tol, int kind)
{
    typedef typename GI::CellIterator                   CI;
    typedef double (*SolutionFuncPtr)(const Vec&);

    typedef Opm::FunctionBoundaryConditions<SolutionFuncPtr> FBC;
    typedef Opm::IncompFlowSolverHybrid<GI, RI, FBC,
                                        Opm::MimeticIPEvaluator> FlowSolver;

    FlowSolver solver;

    // FBC flow_bc;
    // assign_bc(g, flow_bc);
    FBC flow_bc(&u);

    typename CI::Vector gravity(0.0);

    std::cout << "========== Init pressure solver =============" << std::endl;
    Opm::time::StopWatch rolex;
    rolex.start();
    solver.init(g, r, gravity, flow_bc);
    rolex.stop();
    std::cout << "========== Time in seconds: " << rolex.secsSinceStart() << " =============" << std::endl;

    std::vector<double> src(g.numberOfCells(), 0.0);
    assign_src(g, src);
    std::vector<double> sat(g.numberOfCells(), 0.0);


    std::cout << "========== Starting pressure solve =============" << std::endl;
    rolex.start();
    solver.solve(r, sat, flow_bc, src, tol, 3, kind);
    rolex.stop();
    std::cout << "========== Time in seconds: " << rolex.secsSinceStart() << " =============" << std::endl;

    typedef typename FlowSolver::SolutionType FlowSolution;
    FlowSolution soln = solver.getSolution();

    std::vector<typename GI::Vector> cell_velocity;
    estimateCellVelocity(cell_velocity, g, solver.getSolution());
    // Dune's vtk writer wants multi-component data to be flattened.
    std::vector<double> cell_velocity_flat(&*cell_velocity.front().begin(),
                                           &*cell_velocity.back().end());
    std::vector<double> cell_pressure;
    getCellPressure(cell_pressure, g, soln);

    compare_pressure(g, cell_pressure);

    Dune::VTKWriter<typename GI::GridType::LeafGridView> vtkwriter(g.grid().leafGridView());
    vtkwriter.addCellData(cell_velocity_flat, "velocity", GI::GridType::dimension);
    vtkwriter.addCellData(cell_pressure, "pressure");
    vtkwriter.write("testsolution-" + boost::lexical_cast<std::string>(0),
                    Dune::VTK::ascii);
}



int main(int argc, char** argv)
try
{
    Opm::parameter::ParameterGroup param(argc, argv);
    Dune::MPIHelper::instance(argc,argv);

    // Make a Dune::CpGrid.
    typedef Dune::CpGrid Grid;
    Grid grid;
    Opm::initCPGrid(grid , param);
    grid.setUniqueBoundaryIds(true);

    // Make the grid interface
    Opm::GridInterfaceEuler<Grid> g(grid);

    // Reservoir properties.
    Opm::ReservoirPropertyCapillary<Grid::dimension> res_prop;
    res_prop.init(g.numberOfCells(), 1.0, 1.0);
    res_prop.setViscosities(1.0, 1.0);

    test_flowsolver(g, res_prop,
                    param.getDefault("tolerance", 1e-8),
                    param.getDefault("linear_solver_type", 1));
}
catch (const std::exception& e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}
