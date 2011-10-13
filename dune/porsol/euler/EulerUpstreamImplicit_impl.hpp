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



#include <cassert>
#include <cmath>
#include <algorithm>
#include <limits>

#include <dune/common/ErrorMacros.hpp>
#include <dune/common/Average.hpp>
#include <dune/common/Units.hpp>
#include <dune/grid/common/Volumes.hpp>
//#include <dune/porsol/euler/CflCalculator.hpp>
#include <dune/common/StopWatch.hpp>
#include <dune/porsol/opmtransport/examples/ImplicitTransportDefs.hpp>
#include <dune/porsol/opmtransport/examples/transport_source.h>
#include <vector>
#include <array>

namespace Dune
{


    template <class GI, class RP, class BC>
    inline EulerUpstreamImplicit<GI, RP, BC>::EulerUpstreamImplicit()
    {
	//std::array<double, 2> mu  = {{ 1.0, 1.0 }};
	//std::array<double, 2> rho = {{ 0.0, 0.0 }};
	// myfluid_ = TwophaseFluid(mu,rho);
	check_sat_ = true;
	clamp_sat_ = false;
    }	
    template <class GI, class RP, class BC>
    inline EulerUpstreamImplicit<GI, RP, BC>::EulerUpstreamImplicit(const GI& g, const RP& r, const BC& b)
    {
        //residual_computer_.initObj(g, r, b);
	//std::array<double, 2> mu  = {{ 1.0, 1.0 }};
	//std::array<double, 2> rho = {{ 0.0, 0.0 }};
	//myfluid_  = Opm::SimpleFluid2pWrapper<RP>(r);
	check_sat_ = true;
	clamp_sat_ = false;
	//residual_computer_.initObj(g, r, b);
    initObj(g, r, b);
    }



    template <class GI, class RP, class BC>
    inline void EulerUpstreamImplicit<GI, RP, BC>::init(const parameter::ParameterGroup& param)
    {
	check_sat_ = param.getDefault("check_sat", check_sat_);
	clamp_sat_ = param.getDefault("clamp_sat", clamp_sat_);
    }

    template <class GI, class RP, class BC>
    inline void EulerUpstreamImplicit<GI, RP, BC>::init(const parameter::ParameterGroup& param,
						const GI& g, const RP& r, const BC& b)
    {
	init(param);
	initObj(g, r, b);
    }


    template <class GI, class RP, class BC>
    inline void EulerUpstreamImplicit<GI, RP, BC>::initObj(const GI& g, const RP& r, const BC& b)
    {
    	//        residual_computer_.initObj(g, r, b);

    	mygrid_.init(g);
    	porevol_.resize(mygrid_.numCells());
    	for (int i = 0; i < mygrid_.numCells(); ++i){
    		porevol_[i]= mygrid_.cellVolume(i)*r.porosity(i);
    	}


    	myrp_= r;
    	/*
    	typedef typename GI::CellIterator Cit;
    	typedef typename GI::CellIterator Cit;
    	for (FIt f = c->facebegin(); f != c->faceend(); ++f) {
    	  // Neighbour face, will be changed if on a periodic boundary.
    	   FIt nbface = f;
    	  // Compute cell[1], cell_sat[1]
    	   if (f->boundary()) {
    		   if (s.pboundary_->satCond(*f).isPeriodic()) {
    			   nbface = s.bid_to_face_[s.pboundary_->getPeriodicPartner(f->boundaryId())];
    			   ASSERT(nbface != f);
    			   cell[1] = nbface->cellIndex();
    			   ASSERT(cell[0] != cell[1]);
    			   // Periodic faces will be visited twice, but only once
    			   // should they contribute. We make sure that we skip the
    			   // periodic faces half the time.
    			   if (cell[0] > cell[1]) {
    				   // We skip this face.
    				   continue;
    			   }
    		   } else {
    			   ASSERT(s.pboundary_->satCond(*f).isDirichlet());
    			   cell[1] = cell[0];
    			   cell_sat[1] = s.pboundary_->satCond(*f).saturation();
    		   }
    	   }
    	}
    	 */
		//model_ = TransportModel(myfluid_,mygrid_.c_grid(),porevol_,&gravity[0],&trans[0]);
		//myfluid_.init(r);
		//model_.init(myfluid_,mygrid_.c_grid(),porevol_,&gravity[0],&trans[0]);
		//model.init(myfluid_,mygrid_.c_grid(), &porevol_, gravity, &trans);
		//tsolver_ = TransportSolver(model_);
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
    void EulerUpstreamImplicit<GI, RP, BC>::transportSolve(std::vector<double>& saturation,
						   const double time,
						   const typename GI::Vector& gravity,
						   const PressureSolution& pressure_sol,
						   const SparseVector<double>& injection_rates) const
    {

    	ReservoirState<2> state(mygrid_.c_grid());
    	std::vector<double>& sat = state.saturation();
    	for (int i=0; i < mygrid_.numCells(); ++i){
    			sat[2*i] = saturation[i];
    			sat[2*i+1] = 1-saturation[i];
    	}


    	//int count=0;
    	const grid_t* cgrid = mygrid_.c_grid();
    	int numhf = cgrid->cell_facepos[cgrid->number_of_cells];
    	//for (int i=0; i < numhf); ++i){
    	std::vector<double>	faceflux(numhf);
    	for (int c = 0, i = 0; c < cgrid->number_of_cells; ++c){
    		for (; i < cgrid->cell_facepos[c + 1]; ++i) {
    		   	int f= cgrid->cell_faces[i];
    		   	double outflux = pressure_sol.outflux(i);
    		   	double sgn = 2.0*(cgrid->face_cells[2*f + 0] == c) - 1;
    		   	faceflux[f] = sgn * outflux;
    		}
    	}
    	state.faceflux()=faceflux;

		double dt_transport = time;
		int nr_transport_steps = 1;
		time::StopWatch clock;
		int repeats = 0;
		const int max_repeats = 10;
		bool finished = false;
		clock.start();
		std::vector< double > saturation_initial(saturation);
		//std::array< double, 3 >    mygravity;
		std::vector< double >		trans;
		trans.resize(mygrid_.numFaces());

		TwophaseFluid myfluid(myrp_);
		TransportModel model(myfluid,*mygrid_.c_grid(),porevol_,&gravity[0],&trans[0]);
		TransportSolver		tsolver(model);
		Opm::ImplicitTransportDetails::NRReport  rpt_;
		Opm::ImplicitTransportDetails::NRControl ctrl_;
		//Opm::ImplicitTransportLinAlgSupport::CSRMatrixUmfpackSolver linsolve_;
		TransportLinearSolver linsolve_;
		//std::vector<double> totmob(mygrid_.numCells(), 1.0);
	    //std::vector<double> src   (mygrid_.numCells(), 0.0);
		//src[0]                         =  1.0;
		//    src[grid->number_of_cells - 1] = -1.0;
	    TransportSource* tsrc = 0;//create_transport_source(0, 2);
		while (!finished) {
	 	   		for (int q = 0; q < nr_transport_steps; ++q) {
	 	   			tsolver.solve(*mygrid_.c_grid(), tsrc, dt_transport, ctrl_, state, linsolve_, rpt_);
	 	   			if(rpt_.flag<0){
	 	   			  break;
	 	   			}
	    		}
	 	   		if(~(rpt_.flag<0) ){
	 	   			finished =true;
	 	   		}else{
	 	   			if(repeats >max_repeats){
	 	   				finished=true;
	 	   			}else{
	 	   				MESSAGE("Warning: Transport failed, retrying with more steps.");
	 	   				nr_transport_steps *= 2;
	 	   				dt_transport = time/nr_transport_steps;
	 	   				saturation = saturation_initial;
	 	   			}
	 	   		}
	 	   		repeats +=1;
		}
        clock.stop();
        if((rpt_.flag<0)){
        	std::cerr << "EulerUpstreamImplicit did not converge" << std::endl;
        }
#ifdef VERBOSE
        std::cout << "Seconds taken by transport solver: " << clock.secsSinceStart() << std::endl;
#endif // VERBOSE
        for (int i=0; i < mygrid_.numCells(); ++i){
	 		   saturation[i] = sat[2*i];
		}
    }




    /*


    template <class GI, class RP, class BC>
    inline void EulerUpstreamImplicit<GI, RP, BC>::initFinal()
    {
	// Build bid_to_face_ mapping for handling periodic conditions.
	int maxbid = 0;
	for (typename GI::CellIterator c = pgrid_->cellbegin(); c != pgrid_->cellend(); ++c) {
	    for (typename GI::CellIterator::FaceIterator f = c->facebegin(); f != c->faceend(); ++f) {
		int bid = f->boundaryId();
		maxbid = std::max(maxbid, bid);
	    }
	}
	bid_to_face_.clear();
	bid_to_face_.resize(maxbid + 1);
	for (typename GI::CellIterator c = pgrid_->cellbegin(); c != pgrid_->cellend(); ++c) {
	    for (typename GI::CellIterator::FaceIterator f = c->facebegin(); f != c->faceend(); ++f) {
		if (f->boundary() && pboundary_->satCond(*f).isPeriodic()) {
		    bid_to_face_[f->boundaryId()] = f;
		}
	    }
	}

        // Build cell_iters_.
        const int num_cells_per_iter = std::min(50, pgrid_->numberOfCells());
        int counter = 0;
	for (typename GI::CellIterator c = pgrid_->cellbegin(); c != pgrid_->cellend(); ++c, ++counter) {
            if (counter % num_cells_per_iter == 0) {
                cell_iters_.push_back(c);
            }
        }
        cell_iters_.push_back(pgrid_->cellend());
    }

    */




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
