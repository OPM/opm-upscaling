/*===========================================================================
//
// File: spu_2p.cpp
//
// Created: 2011-10-10 10:35:13+0200
//
// Authors: Ingeborg S. Ligaarden <Ingeborg.Ligaarden@sintef.no>
//          Jostein R. Natvig     <Jostein.R.Natvig@sintef.no>
//          Halvor M. Nilsen      <HalvorMoll.Nilsen@sintef.no>
//          Atgeirr F. Rasmussen  <atgeirr@sintef.no>
//          Bård Skaflestad       <Bard.Skaflestad@sintef.no>
//
//==========================================================================*/


/*
  Copyright 2011 SINTEF ICT, Applied Mathematics.
  Copyright 2011 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

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

#include <opm/core/utility/have_boost_redef.hpp>

#include <cstddef>
#include <cassert>

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <vector>

#include <dune/common/array.hh>

#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>

#include <dune/grid/CpGrid.hpp>
#include <dune/porsol/common/ImplicitTransportDefs.hpp>
#include <dune/porsol/common/BCRSMatrixBlockAssembler.hpp>
#include <dune/porsol/common/ReservoirPropertyCapillary.hpp>
#include <dune/porsol/common/setupGridAndProps.hpp>
#include <dune/porsol/common/SimulatorUtilities.hpp>
#include <dune/porsol/common/LinearSolverISTL.hpp>

#include <opm/core/GridAdapter.hpp>

#include <opm/core/linalg/sparse_sys.h>
#include <opm/core/pressure/tpfa/ifs_tpfa.h>
#include <opm/core/pressure/tpfa/trans_tpfa.h>

#include <opm/core/transport/ImplicitAssembly.hpp>
#include <opm/core/transport/ImplicitTransport.hpp>
#include <opm/core/transport/JacobianSystem.hpp>
#include <opm/core/transport/SinglePointUpwindTwoPhase.hpp>

class PressureSolver {
public:
    PressureSolver(UnstructuredGrid* g, const Opm::Rock& rock)
        : htrans_(g->cell_facepos[ g->number_of_cells ]),
          trans_ (g->number_of_faces),
          gpress_(g->cell_facepos[ g->number_of_cells ])
    {
        tpfa_htrans_compute(g, &rock.perm()[0], &htrans_[0]);

        h_ = ifs_tpfa_construct(g, 0);
    }

    ~PressureSolver() {
        ifs_tpfa_destroy(h_);
    }

    template <class State>
    void
    solve(UnstructuredGrid*                      g     ,
          const ::std::vector<double>& totmob,
          const ::std::vector<double>& src   ,
          State&                       state ) {

        tpfa_eff_trans_compute(g, &totmob[0], &htrans_[0], &trans_[0]);

        // No gravity
        ::std::fill(gpress_.begin(), gpress_.end(), double(0.0));

        ifs_tpfa_forces f = { 0 };
        f.src = &src[0];

        ifs_tpfa_assemble(g, &f, &trans_[0], &gpress_[0], h_);

        Opm::LinearSolverISTLAMG linsolve;
        linsolve.solve(h_->A, h_->b, h_->x);

        ifs_tpfa_solution soln = { 0 };
        soln.cell_press = &state.pressure()[0];
        soln.face_flux  = &state.faceflux()[0];

        ifs_tpfa_press_flux(g, &f, &trans_[0], h_, &soln);
    }

private:
    ::std::vector<double> htrans_;
    ::std::vector<double> trans_ ;
    ::std::vector<double> gpress_;

    struct ifs_tpfa_data* h_;
};


typedef Opm::SinglePointUpwindTwoPhase<Opm::TwophaseFluidWrapper> TransportModel;

using namespace Opm::ImplicitTransportDefault;

typedef Dune::FieldVector<double, 1>    ScalarVectorBlockType;
typedef Dune::FieldMatrix<double, 1, 1> ScalarMatrixBlockType;

typedef Dune::BlockVector<ScalarVectorBlockType> ScalarBlockVector;
typedef Dune::BCRSMatrix <ScalarMatrixBlockType> ScalarBCRSMatrix;

typedef NewtonVectorCollection< ScalarBlockVector >          NVecColl;
typedef JacobianSystem        < ScalarBCRSMatrix, NVecColl > JacSys;

typedef Opm::ImplicitTransport<TransportModel  ,
                               JacSys          ,
                               Opm::MaxNormDune,
                               VectorNegater   ,
                               VectorZero      ,
                               MatrixZero      ,
                               VectorAssign    > TransportSolver;





int
main(int argc, char** argv)
{
    Opm::parameter::ParameterGroup param(argc, argv);

    Dune::CpGrid                        cp_grid;
    Dune::ReservoirPropertyCapillary<3> res_prop;

    setupGridAndProps(param, cp_grid, res_prop);

    res_prop.init(cp_grid.size(0), 1, 1);
    res_prop.setViscosities(1.0, 1.0);
    res_prop.setDensities  (0.0, 0.0);

    GridAdapter grid;
    grid.init(cp_grid);

    Opm::Rock rock(grid.c_grid()->number_of_cells, grid.c_grid()->dimensions);

    rock.perm_homogeneous(1);
    rock.poro_homogeneous(1);

    PressureSolver psolver(grid.c_grid(), rock);

    std::vector<double> totmob(grid.c_grid()->number_of_cells, 1.0);
    std::vector<double> src   (grid.c_grid()->number_of_cells, 0.0);

    src[0]                                  =  1.0;
    src[grid.c_grid()->number_of_cells - 1] = -1.0;

    Opm::ReservoirState<> state(grid.c_grid());
    std::vector<double>& sat=state.saturation();
    std::vector<double> tmp_sat(sat.size(),0.3);
    sat =  tmp_sat;

    psolver.solve(grid.c_grid(), totmob, src, state);

    Opm::TransportSource tsrc;
    ::std::vector<double> ssrc (2, 0.0);  ssrc[0] = 1.0;
    ::std::vector<double> ssink(2, 0.0);
    append_transport_source(0, 1.0, src[0]    , ssrc , tsrc);
    append_transport_source(grid.c_grid()->number_of_cells - 1,
                            1.0, src.back(), ssink, tsrc);

    Opm::ImplicitTransportDetails::NRReport  rpt;
    Opm::ImplicitTransportDetails::NRControl ctrl;

    std::vector<double> porevol;
    UnstructuredGrid& cgrid=*grid.c_grid();
    std::vector<double> htrans(cgrid.cell_facepos[ cgrid.number_of_cells ],1);
    Dune::array<double,3> gravity;
    gravity[2]=10.0;
    compute_porevolume(grid.c_grid(), rock, porevol);

    Opm::TwophaseFluidWrapper   fluid  (res_prop);
    TransportModel  model  (fluid, *grid.c_grid(), porevol, &gravity[0]);//, &trans[0]);
    model.initGravityTrans(*grid.c_grid(),htrans);
    TransportSolver tsolver(model);

    double dt   = 1e2;
    ctrl.max_it = 20 ;
    ctrl.rtol = 1e-40;
    ctrl.atol = 1e-12;

    Opm::LinearSolverBICGSTAB linsolve;
    tsolver.solve(*grid.c_grid(), &tsrc, dt, ctrl, state, linsolve, rpt);

    std::cerr << "Number of linear solves: " << rpt.nit        << '\n'
              << "Process converged:       " << (rpt.flag > 0) << '\n'
              << "Convergence flag:        " << rpt.flag       << '\n'
              << "Final residual norm:     " << rpt.norm_res   << '\n'
              << "Final increment norm:    " << rpt.norm_dx    << '\n';

    ::std::ofstream sfile("saturation-00.txt");
    sfile.setf(::std::ios::showpos | ::std::ios::scientific);
    sfile.precision(15);
    ::std::copy(state.saturation().begin(),
                state.saturation().end  (),
                ::std::ostream_iterator<double>(sfile, "\n"));
}
