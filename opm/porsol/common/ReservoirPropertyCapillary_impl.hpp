//===========================================================================
//
// File: ReservoirPropertyCapillary_impl.hpp
//
// Created: Thu Oct 22 20:16:15 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Bård Skaflestad     <bard.skaflestad@sintef.no>
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

#ifndef OPENRS_RESERVOIRPROPERTYCAPILLARY_IMPL_HEADER
#define OPENRS_RESERVOIRPROPERTYCAPILLARY_IMPL_HEADER


namespace Opm
{


    template <int dim>
    double ReservoirPropertyCapillary<dim>::mobilityFirstPhase(int cell_index, double saturation) const
    {
        return relPermFirstPhase(cell_index, saturation) / Super::viscosity1_;
    }


    template <int dim>
    double ReservoirPropertyCapillary<dim>::mobilitySecondPhase(int cell_index, double saturation) const
    {
        return relPermSecondPhase(cell_index, saturation) / Super::viscosity2_;
    }


    template <int dim>
    void ReservoirPropertyCapillary<dim>::phaseMobility(int phase_index,
							int cell_index,
							double saturation,
							double& phase_mob) const
    {
	if (phase_index == 0) {
	    phase_mob = mobilityFirstPhase(cell_index, saturation);
	} else {
	    assert(phase_index == 1);
	    phase_mob = mobilitySecondPhase(cell_index, saturation);
	}
    }


    template <int dim>
    double ReservoirPropertyCapillary<dim>::totalMobility(int cell_index, double saturation) const
    {
        double l1 = mobilityFirstPhase(cell_index, saturation);
        double l2 = mobilitySecondPhase(cell_index, saturation);
        return l1 + l2;
    }


    template <int dim>
    double ReservoirPropertyCapillary<dim>::fractionalFlow(int cell_index, double saturation) const
    {
        double l1 = mobilityFirstPhase(cell_index, saturation);
        double l2 = mobilitySecondPhase(cell_index, saturation);
        return l1/(l1 + l2);
    }


    template <int dim>
    template<class Vector>
    void ReservoirPropertyCapillary<dim>::phaseMobilities(int cell_index, double saturation, Vector& mobility) const
    {
        //assert (mobility.size() >= Super::NumberOfPhases);
        mobility[0] = mobilityFirstPhase(cell_index, saturation);
        mobility[1] = mobilitySecondPhase(cell_index, saturation);
    }


    template <int dim>
    template <class Vector>
    void
    ReservoirPropertyCapillary<dim>::phaseMobilitiesDeriv(int c, double s,
                                                          Vector& dmob) const {

        dmob[0] =   relPermFirstPhaseDeriv (c, s) / Super::viscosity1_;
        dmob[3] = - relPermSecondPhaseDeriv(c, s) / Super::viscosity2_;
        dmob[1] = dmob[2] = 0;
    }



    // ------ Private methods ------


    template <int dim>
    double ReservoirPropertyCapillary<dim>::relPermFirstPhase(int cell_index, double saturation) const
    {
        if (Super::rock_.size() > 0) {
            const int region = Super::cell_to_rock_[cell_index];
            assert (region < int(Super::rock_.size()));
	    double res;
	    Super::rock_[region].krw(saturation, res);
            return res;
        } else {
            // HACK ALERT!
            // Use quadratic rel-perm if no known rock table exists.
            return saturation * saturation;
        }
    }

    template <int dim>
    double
    ReservoirPropertyCapillary<dim>::
    relPermFirstPhaseDeriv(int cell_index, double saturation) const
    {
        if (Super::rock_.size() > 0) {
            const int region = Super::cell_to_rock_[cell_index];
            assert (region < int(Super::rock_.size()));
            double res;
            Super::rock_[region].dkrw(saturation, res);
            return res;
        } else {
            // HACK ALERT!
            // Use quadratic rel-perm if no known rock table exists.
            return 2 * saturation;
        }
    }




    template <int dim>
    double ReservoirPropertyCapillary<dim>::relPermSecondPhase(int cell_index, double saturation) const
    {
        if (Super::rock_.size() > 0) {
            const int region = Super::cell_to_rock_[cell_index];
            assert (region < int(Super::rock_.size()));
	    double res;
	    Super::rock_[region].kro(saturation, res);
            return res;
        } else {
            // HACK ALERT!
            // Use quadratic rel-perm if no known rock table exists.
            return (1 - saturation) * (1 - saturation);
        }
    }

    template <int dim>
    double
    ReservoirPropertyCapillary<dim>::
    relPermSecondPhaseDeriv(int cell_index, double saturation) const
    {
        if (Super::rock_.size() > 0) {
            const int region = Super::cell_to_rock_[cell_index];
            assert (region < int(Super::rock_.size()));
            double res;
            Super::rock_[region].dkro(saturation, res);
            return res;
        } else {
            // HACK ALERT!
            // Use quadratic rel-perm if no known rock table exists.
            return - 2 * (1 - saturation);
        }
    }



    template <int dim>
    void ReservoirPropertyCapillary<dim>::cflFracFlows(int rock, double s, double& ff_first, double& ff_gravity) const
    {
        if (rock == -1) {
            // No rock dependency, we might just as well use the first cell.
            const int cell_index = 0;
            double l1 = mobilityFirstPhase(cell_index, s);
            double l2 = mobilitySecondPhase(cell_index, s);
            ff_first = l1/(l1 + l2);
            ff_gravity = l1*l2/(l1 + l2);
        } else {
            double krw, kro;
	    Super::rock_[rock].krw(s, krw);
            Super::rock_[rock].kro(s, kro);
            double l1 = krw/Super::viscosity1_;
            double l2 = kro/Super::viscosity2_;
            ff_first = l1/(l1 + l2);
            ff_gravity = l1*l2/(l1 + l2);
        }
    }




    template <int dim>
    std::array<double, 3> ReservoirPropertyCapillary<dim>::computeSingleRockCflFactors(int rock, double min_perm, double max_poro) const
    {
        // Make min_perm matrix.
        OwnCMatrix min_perm_matrix(dim, dim, (double*)0);
        eye(min_perm_matrix);
        min_perm_matrix *= min_perm;

        // Sample values at many saturation points.
        const int N = 257;
        double delta = 1.0/double(N - 1);
        double last_ff1, last_ffg;
        double max_der1 = -1e100;
        double max_derg = -1e100;
        cflFracFlows(rock, 0.0, last_ff1, last_ffg);
        double max_ffg = last_ffg;
        double max_derpc = rock == -1 ? 0.0 :
            std::fabs(Super::rock_[rock].capPressDeriv(min_perm_matrix, max_poro, 0.0));
        for (int i = 1; i < N; ++i) {
            double s = double(i)*delta;
            double ff1, ffg;
            cflFracFlows(rock, s, ff1, ffg);
            double est_deriv_ff1 = std::fabs(ff1 - last_ff1)/delta;
            double est_deriv_ffg = std::fabs(ffg - last_ffg)/delta;
            max_der1 = std::max(max_der1, est_deriv_ff1);
            max_derg = std::max(max_derg, est_deriv_ffg);
            max_ffg = std::max(max_ffg, ffg);
            max_derpc = rock == -1 ? 0.0 :
                std::max(max_derpc, std::fabs(Super::rock_[rock].capPressDeriv(min_perm_matrix, max_poro, s)));
            last_ff1 = ff1;
            last_ffg = ffg;
        }
        std::array<double, 3> retval = {{ 1.0/max_der1, 1.0/max_derg, max_ffg*max_derpc }};
        return retval;
    }




    template <int dim>
    void ReservoirPropertyCapillary<dim>::computeCflFactors()
    {
        if (Super::rock_.empty()) {
            std::array<double, 3> fac = computeSingleRockCflFactors(-1, 0.0, 0.0);
            Super::cfl_factor_ = fac[0];
            Super::cfl_factor_gravity_ = fac[1];
            Super::cfl_factor_capillary_ = fac[2];
        } else {
            // Compute min perm and max poro per rock (for J-scaling cap pressure funcs).
            std::vector<double> min_perm(Super::rock_.size(), 1e100);
            std::vector<double> max_poro(Super::rock_.size(), 0.0);
            int num_cells = Super::porosity_.size();
            for (int c = 0; c < num_cells; ++c) {
                int r = Super::cell_to_rock_[c];
                min_perm[r] = std::min(min_perm[r], trace(Super::permeability(c))/double(dim));
                max_poro[r] = std::max(max_poro[r], Super::porosity(c));
            }
            Super::cfl_factor_ = 1e100;
            Super::cfl_factor_gravity_ = 1e100;
            Super::cfl_factor_capillary_ = 0.0;
            for (int r = 0; r < int(Super::rock_.size()); ++r) {
                std::array<double, 3> fac = computeSingleRockCflFactors(r, min_perm[r], max_poro[r]);
                Super::cfl_factor_ = std::min(Super::cfl_factor_, fac[0]);
                Super::cfl_factor_gravity_ = std::min(Super::cfl_factor_gravity_, fac[1]);
                Super::cfl_factor_capillary_ = std::max(Super::cfl_factor_capillary_, fac[2]);
            }
        }
    }





} // namespace Opm


#endif // OPENRS_RESERVOIRPROPERTYCAPILLARY_IMPL_HEADER
