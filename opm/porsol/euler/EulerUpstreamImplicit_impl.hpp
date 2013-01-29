//===========================================================================
//
// File: EulerUpstreamImplicit_impl.hpp
//
// Created: Tue Jun 16 14:25:24 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            B�rd Skaflestad     <bard.skaflestad@sintef.no>
//            Halvor M Nilsen     <hnil@sintef.no>
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

#ifndef OPENRS_EULERUPSTREAMIMPLICIT_IMPL_HEADER
#define OPENRS_EULERUPSTREAMIMPLICIT_IMPL_HEADER

#include <cmath>
#include <algorithm>
#include <limits>
#include <vector>

#include <dune/common/array.hh>

#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/Average.hpp>
#include <opm/core/utility/Units.hpp>
#include <dune/grid/common/Volumes.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/porsol/common/ImplicitTransportDefs.hpp>
#include <opm/core/pressure/tpfa/trans_tpfa.h>

namespace Opm
{


    template <class GI, class RP, class BC>
    inline EulerUpstreamImplicit<GI, RP, BC>::EulerUpstreamImplicit()
    {
        check_sat_ = true;
        clamp_sat_ = false;
    }
    template <class GI, class RP, class BC>
    inline EulerUpstreamImplicit<GI, RP, BC>::EulerUpstreamImplicit(const GI& g, const RP& r, const BC& b)
    {
        check_sat_ = true;
        clamp_sat_ = false;
        initObj(g, r, b);
    }

    template <class GI, class RP, class BC>
    inline void EulerUpstreamImplicit<GI, RP, BC>::init(const Opm::parameter::ParameterGroup& param)
    {
        check_sat_ = param.getDefault("check_sat", check_sat_);
        clamp_sat_ = param.getDefault("clamp_sat", clamp_sat_);
        //Opm::ImplicitTransportDetails::NRControl ctrl_;
        ctrl_.max_it = param.getDefault("transport_nr_max_it", 10);
        max_repeats_ = param.getDefault("transport_max_rep", 10);
        ctrl_.atol  = param.getDefault("transport_atol", 1.0e-6);
        ctrl_.rtol  = param.getDefault("transport_rtol", 5.0e-7);
        ctrl_.max_it_ls = param.getDefault("transport_max_it_ls", 20);
        ctrl_.dxtol = param.getDefault("transport_dxtol", 1e-6);
        ctrl_.verbosity = param.getDefault("transport_verbosity", 0);
    }

    template <class GI, class RP, class BC>
    inline void EulerUpstreamImplicit<GI, RP, BC>::init(const Opm::parameter::ParameterGroup& param,
                                                        const GI& g, const RP& r, const BC& b)
    {
        init(param);
        initObj(g, r, b);
    }


    template <class GI, class RP, class BC>
    inline void EulerUpstreamImplicit<GI, RP, BC>::initObj(const GI& g, const RP& r, const BC& b)
    {
        //residual_computer_.initObj(g, r, b);

        mygrid_.init(g.grid());
        porevol_.resize(mygrid_.numCells());
        for (int i = 0; i < mygrid_.numCells(); ++i){
            porevol_[i]= mygrid_.cellVolume(i)*r.porosity(i);
        }
        // int numf=mygrid_.numFaces();
        int num_cells = mygrid_.numCells();
        int ngconn  = mygrid_.c_grid()->cell_facepos[num_cells];
        //std::vector<double> htrans_(ngconn);
        htrans_.resize(ngconn);
        const double* perm = &(r.permeability(0)(0,0));
        tpfa_htrans_compute(mygrid_.c_grid(), perm, &htrans_[0]);
        // int count = 0;

        myrp_= r;

        typedef typename GI::CellIterator CIt;
        typedef typename CIt::FaceIterator FIt;
        std::vector<FIt> bid_to_face;
        int maxbid = 0;
        for (CIt c = g.cellbegin(); c != g.cellend(); ++c) {
            for (FIt f = c->facebegin(); f != c->faceend(); ++f) {
                int bid = f->boundaryId();
                maxbid = std::max(maxbid, bid);
            }
        }

        bid_to_face.resize(maxbid + 1);
        std::vector<int> egf_cf(mygrid_.numFaces());
        int cix=0;
        for (CIt c = g.cellbegin(); c != g.cellend(); ++c) {
            int loc_fix=0;
            for (FIt f = c->facebegin(); f != c->faceend(); ++f) {
                if (f->boundary() && b.satCond(*f).isPeriodic()) {
                    bid_to_face[f->boundaryId()] = f;
                }
                int egf=f->index();
                int cf=mygrid_.cellFace(cix,loc_fix);
                egf_cf[egf]=cf;
                loc_fix+=1;
            }
            cix+=1;
        }

        const UnstructuredGrid& c_grid=*mygrid_.c_grid();
        int hf_ind=0;
        int bf_ind=0;
        periodic_cells_.resize(0);
        periodic_faces_.resize(0);
        periodic_hfaces_.resize(0);
        periodic_nbfaces_.resize(0);
        //cell1 = cell0;
        direclet_cells_.resize(0);
        direclet_sat_.resize(0);
        direclet_sat_.resize(0);
        direclet_hfaces_.resize(0);

        ASSERT(periodic_cells_.size()==0);
        for (CIt c = g.cellbegin(); c != g.cellend(); ++c) {
            int cell0 = c->index();
            for (FIt f = c->facebegin(); f != c->faceend(); ++f) {
                // Neighbour face, will be changed if on a periodic boundary.
                // Compute cell[1], cell_sat[1]
                FIt nbface = f;
                if (f->boundary()) {
                    bf_ind+=1;
                    if (b.satCond(*f).isPeriodic()) {
                        nbface = bid_to_face[b.getPeriodicPartner(f->boundaryId())];
                        ASSERT(nbface != f);
                        int cell1 = nbface->cellIndex();
                        ASSERT(cell0 != cell1);

                        int f_ind=f->index();

                        int fn_ind=nbface->index();
                        // mapping face indices
                        f_ind=egf_cf[f_ind];
                        fn_ind=egf_cf[fn_ind];
                        ASSERT((c_grid.face_cells[2*f_ind]==-1) || (c_grid.face_cells[2*f_ind+1]==-1));
                        ASSERT((c_grid.face_cells[2*fn_ind]==-1) || (c_grid.face_cells[2*fn_ind+1]==-1));
                        ASSERT((c_grid.face_cells[2*f_ind]==cell0) || (c_grid.face_cells[2*f_ind+1]==cell0));
                        ASSERT((c_grid.face_cells[2*fn_ind]==cell1) || (c_grid.face_cells[2*fn_ind+1]==cell1));
                        periodic_cells_.push_back(cell0);
                        periodic_cells_.push_back(cell1);
                        periodic_faces_.push_back(f_ind);
                        periodic_hfaces_.push_back(hf_ind);
                        periodic_nbfaces_.push_back(fn_ind);
                    } else if (!( b.flowCond(*f).isNeumann() && b.flowCond(*f).outflux() == 0.0)) {
                        //cell1 = cell0;
                        direclet_cells_.push_back(cell0);
                        direclet_sat_.push_back(b.satCond(*f).saturation());
                        direclet_sat_.push_back(1-b.satCond(*f).saturation());//only work for 2 phases
                        direclet_hfaces_.push_back(hf_ind);
                    }
                }
                hf_ind+=1;
            }
        }

        mygrid_.makeQPeriodic(periodic_hfaces_,periodic_cells_);
        // use fractional flow instead of saturation as src
        TwophaseFluid myfluid(myrp_);
        int num_b=direclet_cells_.size();
        for(int i=0; i <num_b; ++i){
            Dune::array<double,2> sat = {{direclet_sat_[2*i] ,direclet_sat_[2*i+1] }};
            Dune::array<double,2> mob;
            Dune::array<double,2*2> dmob;
            myfluid.mobility(direclet_cells_[i], sat, mob, dmob);
            double fl = mob[0]/(mob[0]+mob[1]);
            direclet_sat_[2*i] = fl;
            direclet_sat_[2*i+1] = 1-fl;
        }
    }

    template <class GI, class RP, class BC>
    inline void EulerUpstreamImplicit<GI, RP, BC>::display()
    {
        using namespace std;
        cout << endl;
        cout <<"Displaying some members of EulerUpstreamImplicit" << endl;
        cout << endl;
    }

    template <class GI, class RP, class BC>
    template <class PressureSolution>
    bool EulerUpstreamImplicit<GI, RP, BC>::transportSolve(std::vector<double>& saturation,
                                                           const double time,
                                                           const typename GI::Vector& gravity,
                                                           const PressureSolution& pressure_sol,
                                                           const Opm::SparseVector<double>& injection_rates) const
    {

        Opm::ReservoirState<2> state(mygrid_.c_grid());
        {
            std::vector<double>& sat = state.saturation();
            for (int i=0; i < mygrid_.numCells(); ++i){
                sat[2*i] = saturation[i];
                sat[2*i+1] = 1-saturation[i];
            }
        }

        //int count=0;
        const UnstructuredGrid* cgrid = mygrid_.c_grid();
        int numhf = cgrid->cell_facepos[cgrid->number_of_cells];

        std::vector<double>     faceflux(numhf);

        for (int c = 0, i = 0; c < cgrid->number_of_cells; ++c){
            for (; i < cgrid->cell_facepos[c + 1]; ++i) {
                int f= cgrid->cell_faces[i];
                double outflux = pressure_sol.outflux(i);
                double sgn = 2.0*(cgrid->face_cells[2*f + 0] == c) - 1;
                faceflux[f] = sgn * outflux;
            }
        }
        int num_db=direclet_hfaces_.size();
        std::vector<double> sflux(num_db);
        for (int i=0; i < num_db;++i){
            sflux[i]=-pressure_sol.outflux(direclet_hfaces_[i]);
        }
        state.faceflux()=faceflux;

        double dt_transport = time;
        int nr_transport_steps = 1;
	Opm::time::StopWatch clock;
        int repeats = 0;
        bool finished = false;
        clock.start();

        TwophaseFluid myfluid(myrp_);
        double* tmp_grav=0;
        const UnstructuredGrid& c_grid=*mygrid_.c_grid();
        TransportModel model(myfluid,c_grid,porevol_,tmp_grav);
        model.makefhfQPeriodic(periodic_faces_,periodic_hfaces_, periodic_nbfaces_);
        model.initGravityTrans(*mygrid_.c_grid(),htrans_);
        TransportSolver tsolver(model);
        LinearSolver linsolve_;
        Opm::ImplicitTransportDetails::NRReport  rpt_;

        Opm::TransportSource tsrc;//create_transport_source(0, 2);
        // the input flux is assumed to be the saturation times the flux in the transport solver

        tsrc.nsrc =direclet_cells_.size();
        tsrc.saturation = direclet_sat_;
        tsrc.cell = direclet_cells_;
        tsrc.flux = sflux;

        while (!finished) {
            for (int q = 0; q < nr_transport_steps; ++q) {
                tsolver.solve(*mygrid_.c_grid(), &tsrc, dt_transport, ctrl_, state, linsolve_, rpt_);
                if(rpt_.flag<0){
                    break;
                }
            }
            if(!(rpt_.flag<0) ){
                finished =true;
            }else{
                if(repeats >max_repeats_){
                    finished=true;
                }else{
                    MESSAGE("Warning: Transport failed, retrying with more steps.");
                    nr_transport_steps *= 2;
                    dt_transport = time/nr_transport_steps;
                    if (ctrl_.verbosity){
                        std::cout << "Warning: Transport failed, retrying with more steps. dt = "
                                  << dt_transport/Opm::unit::year << " year.\n";
                    }

                    std::vector<double>& sat = state.saturation();
                    for (int i=0; i < mygrid_.numCells(); ++i){
                        sat[2*i] = saturation[i];
                        sat[2*i+1] = 1-saturation[i];
                    }
                }
            }
            repeats +=1;
        }
        clock.stop();
        std::cout << "EulerUpstreamImplicite used  " << repeats
                  << " repeats and " << nr_transport_steps <<" steps"<< std::endl;
#ifdef VERBOSE
        std::cout << "Seconds taken by transport solver: " << clock.secsSinceStart() << std::endl;
#endif // VERBOSE
        {
            std::vector<double>& sat = state.saturation();
            for (int i=0; i < mygrid_.numCells(); ++i){
                saturation[i] = sat[2*i];
            }
        }
        if((rpt_.flag<0)){
            std::cerr << "EulerUpstreamImplicit did not converge" << std::endl;
            return false;
        }else{
            return true;
        }
    }

    template <class GI, class RP, class BC>
    inline void EulerUpstreamImplicit<GI, RP, BC>::checkAndPossiblyClampSat(std::vector<double>& s) const
    {
        int num_cells = s.size();
        for (int cell = 0; cell < num_cells; ++cell) {
            if (s[cell] > 1.0 || s[cell] < 0.0) {
                if (clamp_sat_) {
                    s[cell] = std::max(std::min(s[cell], 1.0), 0.0);
                } else if (s[cell] > 1.001 || s[cell] < -0.001) {
                    THROW("Saturation out of range in EulerUpstreamImplicit: Cell " << cell << "   sat " << s[cell]);
                }
            }
        }
    }

}



#endif // OPENRS_EULERUPSTREAM_IMPL_HEADER
