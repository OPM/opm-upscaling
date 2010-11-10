/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

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

#ifndef OPM_BLACKOILUPSTREAMEULER_HEADER_INCLUDED
#define OPM_BLACKOILUPSTREAMEULER_HEADER_INCLUDED


#include <dune/common/param/ParameterGroup.hpp>
#include <dune/common/SparseVector.hpp>

namespace Opm
{

    template <class SimpleGrid, class Fluid>
    class BlackoilUpstreamEuler
    {
    public:
	/// @brief
	/// @todo Doc me
	/// @param
        void init(const Dune::parameter::ParameterGroup& param)
        {
        }

	/// @brief
	/// @todo Doc me
	/// @param
	void initObj(const SimpleGrid& grid)
        {
        }

	/// \brief Solve transport equation, evolving \param fluid_state
	/// for \param time seconds.
	/// Cfl type conditions may force many explicit timesteps to
	/// be taken, before the function returns.
	/// @tparam
	/// @param
	template <class PressureSolution>
	void transportSolve(std::vector<typename Fluid::ComponentVec>& fluid_state,
			    const double time,
			    const typename SimpleGrid::Vector& gravity,
			    const PressureSolution& pressure_sol,
			    const Dune::SparseVector<double>& injection_rates) const;
    };

} // namespace Opm

#endif // OPM_BLACKOILUPSTREAMEULER_HEADER_INCLUDED
