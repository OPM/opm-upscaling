//===========================================================================
//
// File: SteadyStateUpscaler_impl.hpp
//
// Created: Fri Aug 28 14:07:51 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Brd Skaflestad     <bard.skaflestad@sintef.no>
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

#ifndef OPM_STEADYSTATEUPSCALER_IMPL_HEADER
#define OPM_STEADYSTATEUPSCALER_IMPL_HEADER


#include <boost/lexical_cast.hpp>
#include <opm/porsol/common/MatrixInverse.hpp>
#include <opm/porsol/common/SimulatorUtilities.hpp>
#include <opm/porsol/common/ReservoirPropertyFixedMobility.hpp>
#include <opm/core/utility/Units.hpp>
#include <algorithm>
#include <iostream>

namespace Opm
{



    template <class Traits>
    inline SteadyStateUpscaler<Traits>::SteadyStateUpscaler()
	: Super(),
          use_gravity_(false),
	  output_vtk_(false),
          print_inoutflows_(false),
	  simulation_steps_(10),
	  stepsize_(0.1),
	  relperm_threshold_(1.0e-8),
          maximum_mobility_contrast_(1.0e9),
          sat_change_threshold_(0.0)
    {
    }




    template <class Traits>
    inline void SteadyStateUpscaler<Traits>::initImpl(const Opm::parameter::ParameterGroup& param)
    {
	Super::initImpl(param);
        use_gravity_ =  param.getDefault("use_gravity", use_gravity_);        
	output_vtk_ = param.getDefault("output_vtk", output_vtk_);
	print_inoutflows_ = param.getDefault("print_inoutflows", print_inoutflows_);
	simulation_steps_ = param.getDefault("simulation_steps", simulation_steps_);
	stepsize_ = Opm::unit::convert::from(param.getDefault("stepsize", stepsize_),
					      Opm::unit::day);
	relperm_threshold_ = param.getDefault("relperm_threshold", relperm_threshold_);
        maximum_mobility_contrast_ = param.getDefault("maximum_mobility_contrast", maximum_mobility_contrast_);
        sat_change_threshold_ = param.getDefault("sat_change_threshold", sat_change_threshold_);

	transport_solver_.init(param);
        // Set viscosities and densities if given.
        double v1_default = this->res_prop_.viscosityFirstPhase();
        double v2_default = this->res_prop_.viscositySecondPhase();
        this->res_prop_.setViscosities(param.getDefault("viscosity1", v1_default), param.getDefault("viscosity2", v2_default));
        double d1_default = this->res_prop_.densityFirstPhase();
        double d2_default = this->res_prop_.densitySecondPhase();
        this->res_prop_.setDensities(param.getDefault("density1", d1_default), param.getDefault("density2", d2_default));
    }




    namespace {
        double maxMobility(double m1, double m2)
        {
            return std::max(m1, m2);
        }
        // The matrix variant expects diagonal mobilities.
        template <class SomeMatrixType>
        double maxMobility(double m1, SomeMatrixType& m2)
        {
            double m = m1;
            for (int i = 0; i < std::min(m2.numRows(), m2.numCols()); ++i) {
                m = std::max(m, m2(i,i));
            }
            return m;
        }
        void thresholdMobility(double& m, double threshold)
        {
            m = std::max(m, threshold);
        }
        // The matrix variant expects diagonal mobilities.
        template <class SomeMatrixType>
        void thresholdMobility(SomeMatrixType& m, double threshold)
        {
            for (int i = 0; i < std::min(m.numRows(), m.numCols()); ++i) {
                m(i, i) = std::max(m(i, i), threshold);
            }
        }
    } // anon namespace




    template <class Traits>
    inline std::pair<typename SteadyStateUpscaler<Traits>::permtensor_t,
                     typename SteadyStateUpscaler<Traits>::permtensor_t>
    SteadyStateUpscaler<Traits>::
    upscaleSteadyState(const int flow_direction,
                       const std::vector<double>& initial_saturation,
		       const double boundary_saturation,
		       const double pressure_drop,
		       const permtensor_t& upscaled_perm)
    {
	static int count = 0;
	++count;
	int num_cells = this->ginterf_.numberOfCells();
	// No source or sink.
	std::vector<double> src(num_cells, 0.0);
	Opm::SparseVector<double> injection(num_cells);
	// Gravity.
	Dune::FieldVector<double, 3> gravity(0.0);
        if (use_gravity_) {
            gravity[2] = Opm::unit::gravity;
        }
	if (gravity.two_norm() > 0.0) {
	    OPM_MESSAGE("Warning: Gravity is experimental for flow solver.");
	}

        // Set up initial saturation profile.
        std::vector<double> saturation = initial_saturation;

        // Set up boundary conditions.
        setupUpscalingConditions(this->ginterf_, this->bctype_, flow_direction,
                                 pressure_drop, boundary_saturation, this->twodim_hack_, this->bcond_);

        // Set up solvers.
        if (flow_direction == 0) {
            this->flow_solver_.init(this->ginterf_, this->res_prop_, gravity, this->bcond_);
        }
        transport_solver_.initObj(this->ginterf_, this->res_prop_, this->bcond_);

        // Run pressure solver.
        this->flow_solver_.solve(this->res_prop_, saturation, this->bcond_, src,
                                 this->residual_tolerance_, this->linsolver_verbosity_, 
                                 this->linsolver_type_, false,
                                 this->linsolver_maxit_, this->linsolver_prolongate_factor_,
                                 this->linsolver_smooth_steps_);
        double max_mod = this->flow_solver_.postProcessFluxes();
        std::cout << "Max mod = " << max_mod << std::endl;

        // Do a run till steady state. For now, we just do some pressure and transport steps...
        std::vector<double> saturation_old = saturation;
        for (int iter = 0; iter < simulation_steps_; ++iter) {
            // Run transport solver.
            transport_solver_.transportSolve(saturation, stepsize_, gravity, this->flow_solver_.getSolution(), injection);

            // Run pressure solver.
            this->flow_solver_.solve(this->res_prop_, saturation, this->bcond_, src,
                                     this->residual_tolerance_, this->linsolver_verbosity_,
                                     this->linsolver_type_, false,
                                     this->linsolver_maxit_, this->linsolver_prolongate_factor_,
                                     this->linsolver_smooth_steps_);
            max_mod = this->flow_solver_.postProcessFluxes();
            std::cout << "Max mod = " << max_mod << std::endl;

            // Print in-out flows if requested.
            if (print_inoutflows_) {
                std::pair<double, double> w_io, o_io;
                computeInOutFlows(w_io, o_io, this->flow_solver_.getSolution(), saturation);
                std::cout << "Pressure step " << iter
                          << "\nWater flow [in] " << w_io.first
                          << "  [out] " << w_io.second
                          << "\nOil flow   [in] " << o_io.first
                          << "  [out] " << o_io.second
                          << std::endl;
            }

            // Output.
            if (output_vtk_) {
                writeVtkOutput(this->ginterf_,
                               this->res_prop_,
                               this->flow_solver_.getSolution(),
                               saturation,
                               std::string("output-steadystate")
                               + '-' + boost::lexical_cast<std::string>(count)
                               + '-' + boost::lexical_cast<std::string>(flow_direction)
                               + '-' + boost::lexical_cast<std::string>(iter));
            }

            // Comparing old to new.
            int num_cells = saturation.size();
            double maxdiff = 0.0;
            for (int i = 0; i < num_cells; ++i) {
                maxdiff = std::max(maxdiff, std::fabs(saturation[i] - saturation_old[i]));
            }
#ifdef VERBOSE
            std::cout << "Maximum saturation change: " << maxdiff << std::endl;
#endif
            if (maxdiff < sat_change_threshold_) {
#ifdef VERBOSE
                std::cout << "Maximum saturation change is under steady state threshold." << std::endl;
#endif
                break;
            }

            // Copy to old.
            saturation_old = saturation;
        }

        // Compute phase mobilities.
        // First: compute maximal mobilities.
        typedef typename Super::ResProp::Mobility Mob;
        Mob m;
        double m1max = 0;
        double m2max = 0;
        for (int c = 0; c < num_cells; ++c) {
            this->res_prop_.phaseMobility(0, c, saturation[c], m.mob);
            m1max = maxMobility(m1max, m.mob);
            this->res_prop_.phaseMobility(1, c, saturation[c], m.mob);
            m2max = maxMobility(m2max, m.mob);
        }
        // Second: set thresholds.
        const double mob1_abs_thres = relperm_threshold_ / this->res_prop_.viscosityFirstPhase();
        const double mob1_rel_thres = m1max / maximum_mobility_contrast_;
        const double mob1_threshold = std::max(mob1_abs_thres, mob1_rel_thres);
        const double mob2_abs_thres = relperm_threshold_ / this->res_prop_.viscositySecondPhase();
        const double mob2_rel_thres = m2max / maximum_mobility_contrast_;
        const double mob2_threshold = std::max(mob2_abs_thres, mob2_rel_thres);
        // Third: extract and threshold.
        std::vector<Mob> mob1(num_cells);
        std::vector<Mob> mob2(num_cells);
        for (int c = 0; c < num_cells; ++c) {
            this->res_prop_.phaseMobility(0, c, saturation[c], mob1[c].mob);
            thresholdMobility(mob1[c].mob, mob1_threshold);
            this->res_prop_.phaseMobility(1, c, saturation[c], mob2[c].mob);
            thresholdMobility(mob2[c].mob, mob2_threshold);
        }

        // Compute upscaled relperm for each phase.
        ReservoirPropertyFixedMobility<Mob> fluid_first(mob1);
        permtensor_t eff_Kw = Super::upscaleEffectivePerm(fluid_first);
        ReservoirPropertyFixedMobility<Mob> fluid_second(mob2);
        permtensor_t eff_Ko = Super::upscaleEffectivePerm(fluid_second);

        // Set the steady state saturation fields for eventual outside access.
        last_saturation_state_.swap(saturation);

	// Compute the (anisotropic) upscaled mobilities.
        // eff_Kw := lambda_w*K
        //  =>  lambda_w = eff_Kw*inv(K); 
	permtensor_t lambda_w(matprod(eff_Kw, inverse3x3(upscaled_perm)));
	permtensor_t lambda_o(matprod(eff_Ko, inverse3x3(upscaled_perm)));

        // Compute (anisotropic) upscaled relative permeabilities.
        // lambda = k_r/mu
        permtensor_t k_rw(lambda_w);
        k_rw *= this->res_prop_.viscosityFirstPhase();
        permtensor_t k_ro(lambda_o);
        k_ro *= this->res_prop_.viscositySecondPhase();
	return std::make_pair(k_rw, k_ro);
    }




    template <class Traits>
    inline const std::vector<double>&
    SteadyStateUpscaler<Traits>::lastSaturationState() const
    {
	return last_saturation_state_;
    }




    template <class Traits>
    double SteadyStateUpscaler<Traits>::lastSaturationUpscaled() const
    {
        typedef typename GridInterface::CellIterator CellIter;
        double pore_vol = 0.0;
        double sat_vol = 0.0;
        for (CellIter c = this->ginterf_.cellbegin(); c != this->ginterf_.cellend(); ++c) {
            double cell_pore_vol = c->volume()*this->res_prop_.porosity(c->index());
            pore_vol += cell_pore_vol;
            sat_vol += cell_pore_vol*last_saturation_state_[c->index()];
        }
        // Dividing by pore volume gives average saturations.
        return sat_vol/pore_vol;
    }




    template <class Traits>
    template <class FlowSol>
    void SteadyStateUpscaler<Traits>::computeInOutFlows(std::pair<double, double>& water_inout,
                                                        std::pair<double, double>& oil_inout,
                                                        const FlowSol& flow_solution,
                                                        const std::vector<double>& saturations) const
    {
        typedef typename GridInterface::CellIterator CellIter;
        typedef typename CellIter::FaceIterator FaceIter;

	double side1_flux = 0.0;
	double side2_flux = 0.0;
	double side1_flux_oil = 0.0;
	double side2_flux_oil = 0.0;
        std::map<int, double> frac_flow_by_bid;
        int num_cells = this->ginterf_.numberOfCells();
        std::vector<double> cell_inflows_w(num_cells, 0.0);
        std::vector<double> cell_outflows_w(num_cells, 0.0);

        // Two passes: First pass, deal with outflow, second pass, deal with inflow.
        // This is for the periodic case, so that we are sure all fractional flows have
        // been set in frac_flow_by_bid.
        for (int pass = 0; pass < 2; ++pass) {
            for (CellIter c = this->ginterf_.cellbegin(); c != this->ginterf_.cellend(); ++c) {
                for (FaceIter f = c->facebegin(); f != c->faceend(); ++f) {
                    if (f->boundary()) {
                        double flux = flow_solution.outflux(f);
                        const SatBC& sc = this->bcond_.satCond(f);
                        if (flux < 0.0 && pass == 1) {
                            // This is an inflow face.
                            double frac_flow = 0.0;
                            if (sc.isPeriodic()) {
                                assert(sc.saturationDifference() == 0.0);
                                int partner_bid = this->bcond_.getPeriodicPartner(f->boundaryId());
                                std::map<int, double>::const_iterator it = frac_flow_by_bid.find(partner_bid);
                                if (it == frac_flow_by_bid.end()) {
                                    OPM_THROW(std::runtime_error, "Could not find periodic partner fractional flow. Face bid = " << f->boundaryId()
                                          << " and partner bid = " << partner_bid);
                                }
                                frac_flow = it->second;
                            } else {
                                assert(sc.isDirichlet());
                                frac_flow = this->res_prop_.fractionalFlow(c->index(), sc.saturation());
                            }
                            cell_inflows_w[c->index()] += flux*frac_flow;
                            side1_flux += flux*frac_flow;
                            side1_flux_oil += flux*(1.0 - frac_flow);
                        } else if (flux >= 0.0 && pass == 0) {
                            // This is an outflow face.
                            double frac_flow = this->res_prop_.fractionalFlow(c->index(), saturations[c->index()]);
                            if (sc.isPeriodic()) {
                                frac_flow_by_bid[f->boundaryId()] = frac_flow;
//                                 std::cout << "Inserted bid " << f->boundaryId() << std::endl;
                            }
                            cell_outflows_w[c->index()] += flux*frac_flow;
                            side2_flux += flux*frac_flow;
                            side2_flux_oil += flux*(1.0 - frac_flow);
                        }
                    }
                }
            }
        }
	water_inout = std::make_pair(side1_flux, side2_flux);
	oil_inout = std::make_pair(side1_flux_oil, side2_flux_oil);
    }




} // namespace Opm


#endif // OPM_STEADYSTATEUPSCALER_IMPL_HEADER
