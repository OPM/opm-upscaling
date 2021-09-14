//===========================================================================
//
// File: ReservoirPropertyCapillary.hpp
//
// Created: Fri Jul  3 12:28:48 2009
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

#ifndef OPENRS_RESERVOIRPROPERTYCAPILLARY_HEADER
#define OPENRS_RESERVOIRPROPERTYCAPILLARY_HEADER

#include <opm/porsol/common/RockJfunc.hpp>
#include <opm/porsol/common/ReservoirPropertyCommon.hpp>
#include <array>

namespace Opm
{

    /// @brief A wrapper for a scalar.
    struct ScalarMobility
    {
	double mob;
	void setToAverage(const ScalarMobility& m1, const ScalarMobility& m2)
	{
	    mob = 0.5*(m1.mob + m2.mob);
	}
	void setToSum(const ScalarMobility& m1, const ScalarMobility& m2)
	{
	    mob = m1.mob + m2.mob;
	}
	void setToInverse(const ScalarMobility& m)
	{
	    mob = 1.0/m.mob;
	}
	template <class Vec>
	Vec multiply(const Vec& v)
	{
	    Vec ret(v);
	    ret *= mob;
	    return ret;
	}
        ScalarMobility& operator *=(const ScalarMobility& other)
        {
            mob *= other.mob;
            return *this;
        }
    };

    /// @brief A property class for incompressible two-phase flow.
    /// @tparam dim the dimension of the space, used for giving permeability tensors the right size.
    template <int dim>
    class ReservoirPropertyCapillary : public ReservoirPropertyCommon<dim, ReservoirPropertyCapillary<dim>, RockJfunc>
    {
    public:
	/// @brief The (scalar) mobility type.
	typedef ScalarMobility Mobility;

	/// @brief Mobility of first (water) phase.
        /// @param cell_index index of a grid cell.
	/// @param saturation a saturation value.
        /// @return mobility value of first phase at the given cell and saturation.
        double mobilityFirstPhase(int cell_index, double saturation) const;

	/// @brief Mobility of second (oil) phase.
        /// @param cell_index index of a grid cell.
	/// @param saturation a saturation value.
        /// @return mobility value of second phase at the given cell and saturation.
        double mobilitySecondPhase(int cell_index, double saturation) const;

	/// @brief Phase mobility.
	/// @param phase_index phase for which to compute mobility.
	/// @param cell_index index of a grid cell.
	/// @param saturation a saturation value.
	/// @param[out] phase_mob phase mobility at the given cell and saturation.
        void phaseMobility(int phase_index, int cell_index, double saturation, double& phase_mob) const;

	/// @brief Total mobility.
        /// @param cell_index index of a grid cell.
	/// @param saturation a saturation value.
        /// @return total mobility value at the given cell and saturation.
        double totalMobility(int cell_index, double saturation) const;

	/// @brief Fractional flow (of the first phase).
        /// @param cell_index index of a grid cell.
	/// @param saturation a saturation value.
        /// @return fractional flow value at the given cell and saturation.
        double fractionalFlow(int cell_index, double saturation) const;

        /// @brief Mobilities for both phases.
	/// @tparam Vector a class with size() and operator[].
        /// @param cell_index index of a grid cell.
	/// @param saturation a saturation value.
        /// @param[out] mobility the phase mobilities at the given cell and saturation.
	///                      Expected to be of size 2 before (and after) the call.
        template<class Vector>
        void phaseMobilities(int cell_index, double saturation, Vector& mobility) const;

        template<class Vector>
        void phaseMobilitiesDeriv(int c, double s, Vector& dmob) const;

	/// @brief Computes cfl factors. Called from ReservoirPropertyCommon::init().
        void computeCflFactors();
    private:
	typedef ReservoirPropertyCommon<dim, ReservoirPropertyCapillary<dim>, RockJfunc> Super;
	// Methods
        double relPermFirstPhase(int cell_index, double saturation) const;
        double relPermSecondPhase(int cell_index, double saturation) const;
        double relPermFirstPhaseDeriv(int cell_index, double saturation) const;
        double relPermSecondPhaseDeriv(int cell_index, double saturation) const;
        void cflFracFlows(int rock, double s, double& ff_first, double& ff_gravity) const;
        std::array<double, 3> computeSingleRockCflFactors(int rock, double min_perm, double max_poro) const;
    };


} // namespace Opm

#include "ReservoirPropertyCapillary_impl.hpp"

#endif // OPENRS_RESERVOIRPROPERTYCAPILLARY_HEADER
