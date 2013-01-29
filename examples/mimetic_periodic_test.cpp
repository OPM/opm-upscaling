//===========================================================================
//
// File: mimetic_periodic_test.cpp
//
// Created: Wed Aug 19 13:42:04 2009
//
// Author(s): Bård Skaflestad     <bard.skaflestad@sintef.no>
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

#include <opm/core/utility/have_boost_redef.hpp>

#include <iostream>

#include <boost/array.hpp>

#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <dune/grid/CpGrid.hpp>

#include <opm/porsol/common/PeriodicHelpers.hpp>
#include <opm/porsol/common/BoundaryConditions.hpp>
#include <opm/porsol/common/GridInterfaceEuler.hpp>
#include <opm/porsol/common/ReservoirPropertyCapillary.hpp>

#include <opm/porsol/mimetic/MimeticIPEvaluator.hpp>
#include <opm/porsol/mimetic/IncompFlowSolverHybrid.hpp>

using namespace Opm;

int main(int argc, char** argv)
{
    typedef Opm::GridInterfaceEuler<Dune::CpGrid>                       GI;
    typedef GI  ::CellIterator                                     CI;
    typedef CI  ::FaceIterator                                     FI;
    typedef Opm::BasicBoundaryConditions<true, false>                  BCs;
    typedef Opm::ReservoirPropertyCapillary<3>                    RI;
    typedef Opm::IncompFlowSolverHybrid<GI, RI, BCs,
                                         Opm::MimeticIPEvaluator> FlowSolver;

    Opm::parameter::ParameterGroup param(argc, argv);
    Dune::CpGrid grid;
    grid.init(param);
    grid.setUniqueBoundaryIds(true);
    GridInterfaceEuler<Dune::CpGrid> g(grid);
    typedef FlowBC FBC;
    Dune::array<FBC, 6> cond = {{ FBC(FBC::Periodic,  1.0*Opm::unit::barsa),
                            FBC(FBC::Periodic, -1.0*Opm::unit::barsa),
                            FBC(FBC::Periodic,  0.0),
                            FBC(FBC::Periodic,  0.0),
                            FBC(FBC::Neumann,   0.0),
                            FBC(FBC::Neumann,   0.0) }};
    BCs fbc;
    createPeriodic(fbc, g, cond);

    RI r;
    r.init(g.numberOfCells());

    CI::Vector gravity;
    gravity[0] = gravity[1] = gravity[2] = 0.0;
#if 0
    gravity[2] = Dune::unit::gravity;
#endif

    FlowSolver solver;
    solver.init(g, r, gravity, fbc);

    std::vector<double> src(g.numberOfCells(), 0.0);
    std::vector<double> sat(g.numberOfCells(), 0.0);

    solver.solve(r, sat, fbc, src);

#if 0
    FlowSolver::SolutionType soln = solver.getSolution();
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
    
    return 0;
}

