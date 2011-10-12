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

#include <cstddef>
#include <cassert>

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <vector>

#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>

#include <dune/grid/CpGrid.hpp>

#include <dune/porsol/common/BCRSMatrixBlockAssembler.hpp>
#include <dune/porsol/common/ReservoirPropertyCapillary.hpp>
#include <dune/porsol/common/setupGridAndProps.hpp>
#include <dune/porsol/common/SimulatorUtilities.hpp>
#include <dune/porsol/common/LinearSolverISTL.hpp>

#include <dune/porsol/opmpressure/src/GridAdapter.hpp>

#include <dune/porsol/opmpressure/src/sparse_sys.h>
#include <dune/porsol/opmpressure/src/ifs_tpfa.h>
#include <dune/porsol/opmpressure/src/trans_tpfa.h>

#include <dune/porsol/opmtransport/src/ImplicitAssembly.hpp>
#include <dune/porsol/opmtransport/src/ImplicitTransport.hpp>
#include <dune/porsol/opmtransport/src/JacobianSystem.hpp>

#include <dune/porsol/opmtransport/src/SinglePointUpwindTwoPhase.hpp>

class Rock {
public:
    Rock(::std::size_t nc, ::std::size_t dim)
        : dim_ (dim           ),
          perm_(nc * dim * dim),
          poro_(nc            ) {}

    const ::std::vector<double>& perm() const { return perm_; }
    const ::std::vector<double>& poro() const { return poro_; }

    void
    perm_homogeneous(double k) {
        setVector(0.0, perm_);

        const ::std::size_t d2 = dim_ * dim_;

        for (::std::size_t c = 0, nc = poro_.size(); c < nc; ++c) {
            for (::std::size_t i = 0; i < dim_; ++i) {
                perm_[c*d2 + i*(dim_ + 1)] = k;
            }
        }
    }

    void
    poro_homogeneous(double phi) {
        setVector(phi, poro_);
    }

private:
    void
    setVector(double x, ::std::vector<double>& v) {
        ::std::fill(v.begin(), v.end(), x);
    }

    ::std::size_t         dim_ ;
    ::std::vector<double> perm_;
    ::std::vector<double> poro_;
};

template <int np = 2>
class ReservoirState {
public:
    ReservoirState(const grid_t* g)
        : press_ (g->number_of_cells),
          fpress_(g->number_of_faces),
          flux_  (g->number_of_faces),
          sat_   (np * g->number_of_cells)
    {}

    ::std::vector<double>& pressure    () { return press_ ; }
    ::std::vector<double>& facepressure() { return fpress_; }
    ::std::vector<double>& faceflux    () { return flux_  ; }
    ::std::vector<double>& saturation  () { return sat_   ; }

    const ::std::vector<double>& faceflux    () const { return flux_; }
    const ::std::vector<double>& saturation  () const { return sat_ ; }

private:
    ::std::vector<double> press_ ;
    ::std::vector<double> fpress_;
    ::std::vector<double> flux_  ;
    ::std::vector<double> sat_   ;
};

class PressureLinearSolver {
public:
    PressureLinearSolver()
    {
        Dune::parameter::ParameterGroup params;

        params.insertParameter("linsolver_tolerance",
                               boost::lexical_cast<std::string>(5.0e-9));

        params.insertParameter("linsoler_verbosity",
                               boost::lexical_cast<std::string>(1));

        params.insertParameter("linsolver_type",
                               boost::lexical_cast<std::string>(1));

        ls_.init(params);
    }

    void
    solve(struct CSRMatrix* A,
          const double*     b,
          double*           x)
    {
        Dune::LinearSolverISTL::LinearSolverResults res =
            ls_.solve(A->m, A->nnz, A->ia, A->ja, A->sa, b, x);
    }

private:
    Dune::LinearSolverISTL ls_;
};

class PressureSolver {
public:
    PressureSolver(grid_t* g, const Rock& rock)
        : htrans_(g->cell_facepos[ g->number_of_cells ]),
          trans_ (g->number_of_faces),
          gpress_(g->cell_facepos[ g->number_of_cells ])
    {
        tpfa_htrans_compute(g, &rock.perm()[0], &htrans_[0]);

        h_ = ifs_tpfa_construct(g);
    }

    ~PressureSolver() {
        ifs_tpfa_destroy(h_);
    }

    template <class State>
    void
    solve(grid_t*                      g     ,
          const ::std::vector<double>& totmob,
          const ::std::vector<double>& src   ,
          State&                       state ) {

        tpfa_eff_trans_compute(g, &totmob[0], &htrans_[0], &trans_[0]);

        // No gravity
        ::std::fill(gpress_.begin(), gpress_.end(), double(0.0));

        ifs_tpfa_assemble(g, &trans_[0], &src[0], &gpress_[0], h_);

        PressureLinearSolver linsolve;
        linsolve.solve(h_->A, h_->b, h_->x);
        
        ifs_tpfa_press_flux(g, &trans_[0], h_,
                            &state.pressure()[0],
                            &state.faceflux()[0]);
    }

private:
    ::std::vector<double> htrans_;
    ::std::vector<double> trans_ ;
    ::std::vector<double> gpress_;

    struct ifs_tpfa_data* h_;
};

class TwophaseFluid {
public:
    TwophaseFluid(const Dune::ReservoirPropertyCapillary<3>& r)
        : r_(r)
    {}
    void init(const Dune::ReservoirPropertyCapillary<3>& r)
    {
        r_ = r;
    }

    template <class Sat ,
              class Mob ,
              class DMob>
    void
    mobility(int c, const Sat& s, Mob& mob, DMob& dmob) const {
        const double s1 = s[0];

        r_.phaseMobilities     (c, s1, mob );
        r_.phaseMobilitiesDeriv(c, s1, dmob);
    }

    double density(int p) const {
        if (p == 0) {
            return r_.densityFirstPhase();
        } else {
            return r_.densitySecondPhase();
        }
    }
    

private:
    Dune::ReservoirPropertyCapillary<3> r_;
};

typedef Opm::SinglePointUpwindTwoPhase<TwophaseFluid> TransportModel;

using namespace Opm::ImplicitTransportDefault;

typedef Dune::FieldVector<double, 1>    ScalarVectorBlockType;
typedef Dune::FieldMatrix<double, 1, 1> ScalarMatrixBlockType;

typedef Dune::BlockVector<ScalarVectorBlockType> ScalarBlockVector;
typedef Dune::BCRSMatrix <ScalarMatrixBlockType> ScalarBCRSMatrix;

typedef NewtonVectorCollection< ScalarBlockVector >          NVecColl;
typedef JacobianSystem        < ScalarBCRSMatrix, NVecColl > JacSys;

class TransportLinearSolver {
public:
    void
    solve(const ScalarBCRSMatrix&  A,
          const ScalarBlockVector& b,
          ScalarBlockVector&       x) {

        Dune::MatrixAdapter<ScalarBCRSMatrix,
                            ScalarBlockVector,
                            ScalarBlockVector> opA(A);

        Dune::SeqILU0<ScalarBCRSMatrix,
                      ScalarBlockVector,
                      ScalarBlockVector> precond(A, 1.0);

        int maxit  = A.N();
        double tol = 5.0e-7;
        int verb   = 1;

        Dune::BiCGSTABSolver<ScalarBlockVector>
            solver(opA, precond, tol, maxit, verb);

        ScalarBlockVector           bcpy(b);
        Dune::InverseOperatorResult res;
        solver.apply(x, bcpy, res);
    }
};

template <class Vector>
class MaxNorm {
public:
    static double
    norm(const Vector& v) {
        return v.infinity_norm();
    }
};

typedef Opm::ImplicitTransport<TransportModel,
                               JacSys        ,
                               MaxNorm       ,
                               VectorNegater ,
                               VectorZero    ,
                               MatrixZero    > TransportSolver;

void
compute_porevolume(const grid_t*        g,
                   const Rock&          rock,
                   std::vector<double>& porevol)
{
    const ::std::vector<double>& poro = rock.poro();

    assert (poro.size() == (::std::size_t)(g->number_of_cells));

    porevol.resize(rock.poro().size());

    ::std::transform(poro.begin(), poro.end(),
                     g->cell_volumes,
                     porevol.begin(),
                     ::std::multiplies<double>());
}

class TransportSource {
public:
    TransportSource() : nsrc(0) {}

    int                   nsrc      ;
    ::std::vector< int  > cell      ;
    ::std::vector<double> pressure  ;
    ::std::vector<double> flux      ;
    ::std::vector<double> saturation;
};

template <class Arr>
void
append_transport_source(int c, double p, double v, const Arr& s,
                        TransportSource& src)
{
    src.cell      .push_back(c);
    src.pressure  .push_back(p);
    src.flux      .push_back(v);
    src.saturation.insert(src.saturation.end(),
                          s.begin(), s.end());
    ++src.nsrc;
}

int
main(int argc, char** argv)
{
    Dune::parameter::ParameterGroup param(argc, argv);

    Dune::CpGrid                        cp_grid;
    Dune::ReservoirPropertyCapillary<3> res_prop;

    setupGridAndProps(param, cp_grid, res_prop);

    res_prop.init(cp_grid.size(0), 1, 1);
    res_prop.setViscosities(1.0, 1.0);
    res_prop.setDensities  (0.0, 0.0);

    GridAdapter grid;
    grid.init(cp_grid);

    Rock rock(grid.c_grid()->number_of_cells, grid.c_grid()->dimensions);

    rock.perm_homogeneous(1);
    rock.poro_homogeneous(1);

    PressureSolver psolver(grid.c_grid(), rock);

    std::vector<double> totmob(grid.c_grid()->number_of_cells, 1.0);
    std::vector<double> src   (grid.c_grid()->number_of_cells, 0.0);

    src[0]                                  =  1.0;
    src[grid.c_grid()->number_of_cells - 1] = -1.0;

    ReservoirState<> state(grid.c_grid());

    psolver.solve(grid.c_grid(), totmob, src, state);

    TransportSource tsrc;
    ::std::vector<double> ssrc (2, 0.0);  ssrc[0] = 1.0;
    ::std::vector<double> ssink(2, 0.0);
    append_transport_source(0, 1.0, src[0]    , ssrc , tsrc);
    append_transport_source(grid.c_grid()->number_of_cells - 1,
                            1.0, src.back(), ssink, tsrc);

    Opm::ImplicitTransportDetails::NRReport  rpt;
    Opm::ImplicitTransportDetails::NRControl ctrl;

    std::vector<double> porevol;
    compute_porevolume(grid.c_grid(), rock, porevol);

    TwophaseFluid   fluid  (res_prop);
    TransportModel  model  (fluid, *grid.c_grid(), porevol);
    TransportSolver tsolver(model);

    double dt   = 1e2;
    ctrl.max_it = 20 ;

    TransportLinearSolver linsolve;
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
