//===========================================================================
//
// File: EulerUpstream_impl.hpp
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

#ifndef OPENRS_EULERUPSTREAM_IMPL_HEADER
#define OPENRS_EULERUPSTREAM_IMPL_HEADER



#include <opm/common/ErrorMacros.hpp>
#include <opm/core/utility/Average.hpp>
#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/grid/common/Volumes.hpp>
#include <opm/porsol/euler/CflCalculator.hpp>
#include <opm/grid/utility/StopWatch.hpp>

#include <cassert>
#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>

namespace Opm
{


    template <class GI, class RP, class BC>
    inline EulerUpstream<GI, RP, BC>::EulerUpstream()
	: method_viscous_(true),
	  method_gravity_(true),
	  method_capillary_(true),
	  use_cfl_viscous_(true),
	  use_cfl_gravity_(true),
	  use_cfl_capillary_(true),
	  courant_number_(0.5),
	  minimum_small_steps_(1),
          maximum_small_steps_(10000),
	  check_sat_(true),
	  clamp_sat_(false)
    {
    }

    template <class GI, class RP, class BC>
    inline EulerUpstream<GI, RP, BC>::EulerUpstream(const GI& g, const RP& r, const BC& b)
	: method_viscous_(true),
	  method_gravity_(true),
	  method_capillary_(true),
	  use_cfl_viscous_(true),
	  use_cfl_gravity_(true),
	  use_cfl_capillary_(true),
	  courant_number_(0.5),
	  minimum_small_steps_(1),
          maximum_small_steps_(10000),
	  check_sat_(true),
	  clamp_sat_(false)
    {
        residual_computer_.initObj(g, r, b);
    }



    template <class GI, class RP, class BC>
    inline void EulerUpstream<GI, RP, BC>::init(const Opm::ParameterGroup& param)
    {
	courant_number_ = param.getDefault("courant_number", courant_number_);
	method_viscous_ = param.getDefault("method_viscous", method_viscous_);
	method_gravity_ = param.getDefault("method_gravity", method_gravity_);
	method_capillary_ = param.getDefault("method_capillary", method_capillary_);
	use_cfl_viscous_ = param.getDefault("use_cfl_viscous", use_cfl_viscous_);
	use_cfl_gravity_ = param.getDefault("use_cfl_gravity", use_cfl_gravity_);
	use_cfl_capillary_ = param.getDefault("use_cfl_capillary", use_cfl_capillary_);
	minimum_small_steps_ = param.getDefault("minimum_small_steps", minimum_small_steps_);
	maximum_small_steps_ = param.getDefault("maximum_small_steps", maximum_small_steps_);
	check_sat_ = param.getDefault("check_sat", check_sat_);
	clamp_sat_ = param.getDefault("clamp_sat", clamp_sat_);
    }

    template <class GI, class RP, class BC>
    inline void EulerUpstream<GI, RP, BC>::init(const Opm::ParameterGroup& param,
						const GI& g, const RP& r, const BC& b)
    {
	init(param);
	initObj(g, r, b);
    }


    template <class GI, class RP, class BC>
    inline void EulerUpstream<GI, RP, BC>::initObj(const GI& g, const RP& r, const BC& b)
    {
        residual_computer_.initObj(g, r, b);
        porevol_.resize(g.numberOfCells());
        for (CIt c = g.cellbegin(); c != g.cellend(); ++c) {
            porevol_[c->index()] = c->volume()*r.porosity(c->index());
        }
    }



    template <class GI, class RP, class BC>
    inline void EulerUpstream<GI, RP, BC>::display()
    {
	using namespace std;
	cout << endl;
	cout <<"Displaying some members of EulerUpstream" << endl;
	cout << endl;
	cout << "courant_number = " << courant_number_ << endl;
    }



    template <class GI, class RP, class BC>
    inline void EulerUpstream<GI, RP, BC>::setCourantNumber(double cn)
    {
	courant_number_ = cn;
    }



    template <class GI, class RP, class BC>
    template <class PressureSolution>
    void EulerUpstream<GI, RP, BC>::transportSolve(std::vector<double>& saturation,
						   const double time,
						   const typename GI::Vector& gravity,
						   const PressureSolution& pressure_sol,
						   const Opm::SparseVector<double>& injection_rates) const
    {
	// Compute the cfl time-step.
	double cfl_dt = computeCflTime(saturation, time, gravity, pressure_sol);

	// Compute the number of small steps to take, and the actual small timestep.
	int nr_transport_steps;
	if (cfl_dt > time){
	    nr_transport_steps = minimum_small_steps_;
	} else {
            double steps = std::min<double>(std::ceil(time/cfl_dt), std::numeric_limits<int>::max());
            nr_transport_steps = int(steps);
	    nr_transport_steps = std::max(nr_transport_steps, minimum_small_steps_);
            nr_transport_steps = std::min(nr_transport_steps, maximum_small_steps_);
	}
	double dt_transport = time/nr_transport_steps;

	// Do the timestepping. The try-catch blocks are there to handle
	// the situation that smallTimeStep throws, which may happen due
	// to saturation out of bounds (if check_sat_ is true).
	// We cannot guarantee that this does not happen, since we do not
	// (yet) compute a capillary cfl condition.
	// Using exception for "alternate control flow" like this is bad
	// design, should rather use error return values for this.
	std::vector<double> saturation_initial(saturation);
	bool finished = false;
	int repeats = 0;
	const int max_repeats = 10;
	Opm::time::StopWatch clock;
        clock.start();
	while (!finished) {
	    try {
#ifdef VERBOSE
		std::cout << "Doing " << nr_transport_steps
			  << " steps for saturation equation with stepsize "
			  << dt_transport << " in seconds." << std::endl;
#endif // VERBOSE
		for (int q = 0; q < nr_transport_steps; ++q) {
		    smallTimeStep(saturation,
				  dt_transport,
				  gravity,
				  pressure_sol,
                                  injection_rates);
		}
		finished = true;
	    }
	    catch (...) {
		++repeats;
		if (repeats > max_repeats) {
		    throw;
		}
		OPM_MESSAGE("Warning: Transport failed, retrying with more steps.");
		nr_transport_steps *= 2;
		dt_transport = time/nr_transport_steps;
		saturation = saturation_initial;
	    }
	}
        clock.stop();
#ifdef VERBOSE
        std::cout << "Seconds taken by transport solver: " << clock.secsSinceStart() << std::endl;
#endif // VERBOSE
    }





    /*


    template <class GI, class RP, class BC>
    inline void EulerUpstream<GI, RP, BC>::initFinal()
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
    template <class PressureSolution>
    inline double EulerUpstream<GI, RP, BC>::computeCflTime(const std::vector<double>& /*saturation*/,
#ifdef VERBOSE
							    const double time,
#else
							    const double,
#endif
							    const typename GI::Vector& gravity,
							    const PressureSolution& pressure_sol) const
    {
	// Deal with cfl stuff, compute the necessary number of small time steps.
	double cfl_dt_v = 1e99;
	double cfl_dt_g = 1e99;
	double cfl_dt_c = 1e99;

	// Viscous cfl.
	if (method_viscous_ && use_cfl_viscous_) {
	    cfl_dt_v = cfl_calculator::findCFLtimeVelocity(residual_computer_.grid(),
							   residual_computer_.reservoirProperties(),
							   pressure_sol);
#ifdef VERBOSE
	    std::cout << "CFL dt for velocity is  "
                      << cfl_dt_v << " seconds   ("
                      << Opm::unit::convert::to(cfl_dt_v, Opm::unit::day)
                      << " days)." << std::endl;
#endif // VERBOSE
	}

	// Gravity cfl.
	if (method_gravity_ && use_cfl_gravity_) {
	    cfl_dt_g = cfl_calculator::findCFLtimeGravity(residual_computer_.grid(),
                                                          residual_computer_.reservoirProperties(),
                                                          gravity);
#ifdef VERBOSE
	    std::cout << "CFL dt for gravity is   "
                      << cfl_dt_g << " seconds   ("
                      << Opm::unit::convert::to(cfl_dt_g, Opm::unit::day)
                      << " days)." << std::endl;
#endif // VERBOSE
	}

	// Capillary cfl.
	if (method_capillary_ && use_cfl_capillary_) {
            cfl_dt_c = cfl_calculator::findCFLtimeCapillary(residual_computer_.grid(),
                                                            residual_computer_.reservoirProperties());
#ifdef VERBOSE
	    std::cout << "CFL dt for capillary term is "
                      << cfl_dt_c << " seconds   ("
                      << Opm::unit::convert::to(cfl_dt_c, Opm::unit::day)
                      << " days)." << std::endl;
#endif // VERBOSE
	}

	double cfl_dt = std::min(std::min(cfl_dt_v, cfl_dt_g), cfl_dt_c);
	cfl_dt *= courant_number_;
#ifdef VERBOSE
        std::cout << "Total impes time is          "
                  << time << " seconds   ("
                  << Opm::unit::convert::to(time, Opm::unit::day)
                  << " days)." << std::endl;

	std::cout << "Final modified CFL dt is     "
                  << cfl_dt << " seconds   ("
                  << Opm::unit::convert::to(cfl_dt, Opm::unit::day)
                  << " days)." << std::endl;
#endif // VERBOSE
	return cfl_dt;
    }




    template <class GI, class RP, class BC>
    inline void EulerUpstream<GI, RP, BC>::checkAndPossiblyClampSat(std::vector<double>& s) const
    {
	int num_cells = s.size();
	for (int cell = 0; cell < num_cells; ++cell) {
	    if (s[cell] > 1.0 || s[cell] < 0.0) {
		if (clamp_sat_) {
		    s[cell] = std::max(std::min(s[cell], 1.0), 0.0);
		} else if (s[cell] > 1.001 || s[cell] < -0.001) {
            OPM_THROW(std::runtime_error,
                      "Saturation out of range in EulerUpstream: "
                      "Cell: " + std::to_string(cell) +
                      "   sat " + std::to_string(s[cell]));
		}
	    }
	}
    }




	
    template <class GI, class RP, class BC>
    template <class PressureSolution>
    inline void EulerUpstream<GI, RP, BC>::smallTimeStep(std::vector<double>& saturation,
							 const double dt,
							 const typename GI::Vector& gravity,
							 const PressureSolution& pressure_sol,
                                                         const Opm::SparseVector<double>& injection_rates) const
    {
        if (method_capillary_) {
            residual_computer_.computeCapPressures(saturation);
        }
	residual_computer_.computeResidual(saturation, gravity, pressure_sol, injection_rates,
                                           method_viscous_, method_gravity_, method_capillary_,
                                           residual_);
// 	double max_ok_dt = 1e100;
//         const double tol = 1e-10;
	int num_cells = saturation.size();
	for (int i = 0; i < num_cells; ++i) {
	    const double sat_change = dt*residual_[i]/porevol_[i];
	    saturation[i] += sat_change;
// 	    if (sat_change > tol) {
// 		max_ok_dt = std::min(max_ok_dt, (1.0 - saturation[i])/sc);
// 	    } else if (sat_change < -tol) {
// 		max_ok_dt = std::min(max_ok_dt, -saturation[i]/sc);
// 	    }
	}
// 	std::cout << "Maximum nonviolating timestep is " << max_ok_dt << " seconds\n";
	if (check_sat_ || clamp_sat_) {
	    checkAndPossiblyClampSat(saturation);
	}
    }


} // end namespace Opm


#endif // OPENRS_EULERUPSTREAM_IMPL_HEADER
