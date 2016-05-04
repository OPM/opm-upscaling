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

#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <dune/grid/CpGrid.hpp>

#include <opm/porsol/common/PeriodicHelpers.hpp>
#include <opm/porsol/common/BoundaryConditions.hpp>
#include <opm/porsol/common/GridInterfaceEuler.hpp>
#include <opm/porsol/common/ReservoirPropertyCapillary.hpp>
#include <opm/porsol/mimetic/MimeticIPEvaluator.hpp>
#include <opm/porsol/mimetic/IncompFlowSolverHybrid.hpp>

#include <opm/upscaling/initCPGrid.hpp>



#include <array>
#include <iostream>




int main(int argc, char** argv)
try
{
    typedef Opm::GridInterfaceEuler<Dune::CpGrid>      GI;
    typedef GI ::CellIterator                          CI;
    typedef Opm::BasicBoundaryConditions<true, false>  BCs;
    typedef Opm::ReservoirPropertyCapillary<3>         RI;
    typedef Opm::IncompFlowSolverHybrid<GI, RI, BCs,
                                        Opm::MimeticIPEvaluator> FlowSolver;

    Dune::CpGrid grid;
    {
        Opm::parameter::ParameterGroup param(argc, argv);
        Opm::initCPGrid(grid , param);
        grid.setUniqueBoundaryIds(true);
    }

    const Opm::GridInterfaceEuler<Dune::CpGrid> g(grid);

    using FBC = Opm::FlowBC;
    const auto cond =
        std::array<FBC, 6>{{ FBC(FBC::Periodic,  1.0*Opm::unit::barsa),
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

    FlowSolver solver{};
    solver.init(g, r, gravity, fbc);

    std::vector<double> src(g.numberOfCells(), 0.0);
    std::vector<double> sat(g.numberOfCells(), 0.0);

    solver.solve(r, sat, fbc, src);
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}

