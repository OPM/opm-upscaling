//===========================================================================
//
// File: mimetic_solver_test.cpp
//
// Created: Tue Jul  7 15:35:21 2009
//
// Author(s): B�rd Skaflestad     <bard.skaflestad@sintef.no>
//            Atgeirr F Rasmussen <atgeirr@sintef.no>
//
// $Date: 2009-09-01 20:55:38 +0200 (Tue, 01 Sep 2009) $
//
// $Revision: 441 $
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

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <opm/common/utility/platform_dependent/disable_warnings.h>

#include <dune/common/version.hh>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/grid/CpGrid.hpp>

#include <opm/parser/eclipse/Units/Units.hpp>
#include <opm/common/utility/parameters/ParameterGroup.hpp>

#include <opm/output/eclipse/EclipseGridInspector.hpp>

#include <opm/porsol/common/BoundaryConditions.hpp>
#include <opm/porsol/common/GridInterfaceEuler.hpp>
#include <opm/porsol/common/Matrix.hpp>
#include <opm/porsol/common/ReservoirPropertyCapillaryAnisotropicRelperm.hpp>
#include <opm/porsol/common/blas_lapack.hpp>
#include <opm/porsol/common/fortran.hpp>
#include <opm/porsol/mimetic/IncompFlowSolverHybrid.hpp>
#include <opm/porsol/mimetic/MimeticIPAnisoRelpermEvaluator.hpp>

#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <memory>

namespace
{
    void build_grid(const Opm::Deck&  deck,
                    Dune::CpGrid&      grid,
                    std::array<int,3>& cartDims)
    {
        auto g = grdecl{};
        {
            const auto insp = Opm::EclipseGridInspector(deck);

            cartDims[0] = g.dims[0] = insp.gridSize()[0];
            cartDims[1] = g.dims[1] = insp.gridSize()[1];
            cartDims[2] = g.dims[2] = insp.gridSize()[2];
        }

        g.coord = &deck["COORD"].back().getSIDoubleData()[0];
        g.zcorn = &deck["ZCORN"].back().getSIDoubleData()[0];

        if (deck.hasKeyword("ACTNUM")) {
            g.actnum = &deck["ACTNUM"].back().getIntData()[0];
            grid.processEclipseFormat(g, false, false);
        }
        else {
            const auto dflt_actnum =
                std::vector<int>(g.dims[0] * g.dims[1] * g.dims[2], 1);

            g.actnum = dflt_actnum.data();
            grid.processEclipseFormat(g, false, false);
        }
    }
} // namespace anonymous


template <int dim, class Interface>
void test_evaluator(const Interface& g)
{
    typedef typename Interface::CellIterator CI;
    typedef typename CI       ::FaceIterator FI;
    typedef typename CI       ::Scalar       Scalar;

    typedef Opm::SharedFortranMatrix FMat;

    std::cout << "Called test_evaluator()" << std::endl;

    std::vector<int> numf; numf.reserve(g.numberOfCells());
    int max_nf = -1;
    for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
        numf.push_back(0);
        int& nf = numf.back();

        for (FI f = c->facebegin(); f != c->faceend(); ++f) {
            ++nf;
        }

        max_nf = std::max(max_nf, nf);
    }

    typedef int DummyClass;
    Opm::MimeticIPAnisoRelpermEvaluator<Interface, DummyClass> ip(max_nf);

    // Set dummy permeability K=diag(10,1,...,1,0.1).
    std::vector<Scalar> perm(dim * dim, Scalar(0.0));
    Opm::SharedCMatrix K(dim, dim, &perm[0]);
    for (int i = 0; i < dim; ++i)
        K(i,i) = 1.0;
    K(0    ,0    ) *= 10.0;
    K(dim-1,dim-1) /= 10.0;

    // Storage for inverse ip.
    std::vector<Scalar> ip_store(max_nf * max_nf, Scalar(0.0));

    // Loop grid whilst building (and outputing) the inverse IP matrix.
    int count = 0;
    for (CI c = g.cellbegin(); c != g.cellend(); ++c, ++count) {
        FMat Binv(numf[count], numf[count], &ip_store[0]);

        ip.evaluate(c, K, Binv);

        std::cout << count << " -> Binv = [\n" << Binv << "]\n";
    }
}


template<int dim, class RI>
void assign_permeability(RI& r, int nc, double k)
{
    typedef typename RI::SharedPermTensor Tensor;

    for (int c = 0; c < nc; ++c) {
        Tensor K = r.permeabilityModifiable(c);
        for (int i = 0; i < dim; ++i) {
            K(i,i) = k;
        }
    }
}


template<int dim, class GI, class RI>
void test_flowsolver(const GI& g, const RI& r)
{
    typedef typename GI::CellIterator                 CI;
    typedef Opm::BasicBoundaryConditions<true, false> FBC;

    using FlowSolver =
        Opm::IncompFlowSolverHybrid<GI, RI, FBC,
                                    Opm::MimeticIPAnisoRelpermEvaluator>;

    FlowSolver solver;

    typedef Opm::FlowBC BC;
    FBC flow_bc(7);
    flow_bc.flowCond(5) = BC(BC::Dirichlet, 100.0*Opm::unit::barsa);

    typename CI::Vector gravity;
    gravity[0] = gravity[1] = 0.0;
    gravity[2] = Opm::unit::gravity;

    solver.init(g, r, gravity, flow_bc);

    std::vector<double> src(g.numberOfCells(), 0.0);
    std::vector<double> sat(g.numberOfCells(), 0.0);

    solver.solve(r, sat, flow_bc, src);
}


int main(int argc, char** argv)
try
{
    Opm::ParameterGroup param(argc, argv);
    Dune::MPIHelper::instance(argc, argv);

    auto parser = std::make_shared<Opm::Parser>();

    const Opm::Deck& deck(parser->parseFile(param.get<std::string>("filename")));

    // Make a grid
    Dune::CpGrid grid;
    {
        std::array<int,3> cartDims;
        build_grid(deck, grid, cartDims);
    }

    // Make the grid interface
    Opm::GridInterfaceEuler<Dune::CpGrid> g(grid);

    // Reservoir properties.
    auto res_prop = Opm::ReservoirPropertyCapillaryAnisotropicRelperm<3>{};
    res_prop.init(deck, grid.globalCell());

    assign_permeability<3>(res_prop, g.numberOfCells(), 0.1*Opm::unit::darcy);

    test_flowsolver<3>(g, res_prop);
}
catch (const std::exception& e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}
