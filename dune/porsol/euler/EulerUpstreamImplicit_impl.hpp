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
#include <vector>
#include <array>

namespace Dune
{


    template <class GI, class RP, class BC>
    inline EulerUpstreamImplicit<GI, RP, BC>::EulerUpstreamImplicit()
    {
	std::array<double, 2> mu  = {{ 1.0, 1.0 }};
	std::array<double, 2> rho = {{ 0.0, 0.0 }};
	myfluid_ = TwophaseFluid(mu,rho);
	check_sat_ = true;
	clamp_sat_ = false;
    }	
    template <class GI, class RP, class BC>
    inline EulerUpstreamImplicit<GI, RP, BC>::EulerUpstreamImplicit(const GI& g, const RP& r, const BC& b)
    {
        //residual_computer_.initObj(g, r, b);
	std::array<double, 2> mu  = {{ 1.0, 1.0 }};
	std::array<double, 2> rho = {{ 0.0, 0.0 }};
	myfluid_ = TwophaseFluid(mu,rho);
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
    	porevol_.resize(g.numberOfCells());
    	for (int i = 0; i < mygrid_.numCells(); ++i){
    		porevol_[i]= g.cell_volume_[i]*r.porosity(i);
    	}
    	typedef typename GI::CellIterator CIt;
    	typedef typename CIt::FaceIterator FIt;
    	typedef typename FIt::Vector Vector;

    	dunefaceind_.resize(2*g.numberOfFaces());
    	//state.faceflux() = pressure_sol.flux;
    	//std::vector<double>& flux = state.faceflux();
    	int count=0;
    	for (CIt c = c->cellbegin(); c != c->cellend(); ++c) {
    		for (FIt f = c->facebegin(); f != c->faceend(); ++f) {
    			dunefaceind_[f->index()] = f->index();
    			count=count+1;
    		}
    	}
    	/*
        for (CIt c = g.cellbegin(); c != g.cellend(); ++c) {
            porevol_[c->index()] = c->volume()*r.porosity(c->index());
       	 }
     	*/
    	std::array< double, 3 >    gravity;
    	std::vector< double >		trans;
    	trans.resize(mygrid_.numCells());
		model_ = TransportModel(myfluid_,mygrid_.c_grid(),porevol_,gravity);
		tsolver_ = TransportSolver(model_);
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
    	typedef typename GI::CellIterator CIt;
    	typedef typename CIt::FaceIterator FIt;
    	typedef typename FIt::Vector Vector;
    	int count=0;
    	for (CIt c = c->cellbegin(); c != c->cellend(); ++c) {
    		for (FIt f = c->facebegin(); f != c->faceend(); ++f) {
    			faceflux_[dunefaceind_[count]] = pressure_sol.outflux(f);
    			count=count+1;
    		}
    	}

		double dt_transport = time;
		int nr_transport_steps = 1;
		time::StopWatch clock;
		int repeats = 0;
		const int max_repeats = 10;
		bool finished = false;
		clock.start();
		std::vector< double > saturation_initial(saturation);

		while (!finished) {
			try {
	 	   		for (int q = 0; q < nr_transport_steps; ++q) {
	    			//grid_t* pp=mygrid_.c_grid();
	    			//tsolver_.solve(*pp, NULL, dt_transport, ctrl_, state, linsolve_, rpt_);
	    		}
#ifdef VERBOSE
		    	std::cout << "Doing " << 1
				  << " steps for saturation equation with stepsize "
					  << dt_transport << " in seconds." << std::endl;
#endif // VERBOSE
	   	 	}
	    	catch (...) {
				if (repeats > max_repeats) {
		    		throw;
				}
				MESSAGE("Warning: Transport failed, retrying with more steps.");
				nr_transport_steps *= 2;
				dt_transport = time/nr_transport_steps;
				saturation = saturation_initial;
	    	}
		}

        clock.stop();
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
