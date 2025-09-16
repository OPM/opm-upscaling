//===========================================================================
//
// File: ImplicitCapillarity.hpp
//
// Created: Thu May  6 15:29:51 2010
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Jostein R Natvig    <jostein.r.natvig@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2010 Statoil ASA.

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

#ifndef OPENRS_IMPLICITCAPILLARITY_HEADER
#define OPENRS_IMPLICITCAPILLARITY_HEADER


#include <opm/porsol/euler/EulerUpstreamResidual.hpp>
#include <opm/common/utility/parameters/ParameterGroup.hpp>
#include <opm/common/utility/numeric/SparseVector.hpp>
#include <opm/porsol/mimetic/IncompFlowSolverHybrid.hpp>
#include <opm/porsol/mimetic/MimeticIPEvaluator.hpp>

namespace Opm {

    /// Class for doing simple transport by explicit Euler upstream method for general grid.
    /// @tparam
    template <class GridInterface, class ReservoirProperties, class BoundaryConditions,
              template <class, class> class InnerProd = MimeticIPEvaluator>
    class ImplicitCapillarity
    {
    public:
        typedef IncompFlowSolverHybrid<GridInterface,
                                       ReservoirProperties, 
                                       BoundaryConditions,
                                       InnerProd> PressureSolver;

	/// @brief
	/// @todo Doc me
	ImplicitCapillarity();
	/// @brief
	/// @todo Doc me
	/// @param
 	ImplicitCapillarity(const GridInterface& grid,
                            const ReservoirProperties& resprop,
                            const BoundaryConditions& boundary);
	/// @brief
	/// @todo Doc me
	/// @param
	void init(const Opm::ParameterGroup& param);
	/// @brief
	/// @todo Doc me
	/// @param
	void init(const Opm::ParameterGroup& param,
		  const GridInterface& grid,
		  const ReservoirProperties& resprop,
		  const BoundaryConditions& boundary);
	/// @brief
	/// @todo Doc me
	/// @param
	void initObj(const GridInterface& grid,
		     const ReservoirProperties& resprop,
		     const BoundaryConditions& boundary);


    /// @brief Solve transport equation.
    /// @param saturation the evolving saturation
    /// @param time Time in seconds.
    /// @param gravity Gravity
    /// @param pressure_sol Pressure solution
    /// @param injection_rates Injection ratees
    /// @details Cfl type conditions may force many explicit timesteps to
    ///          be taken, before the function returns.
    template <class PressureSolution>
	void transportSolve(std::vector<double>& saturation,
			    const double time,
			    const typename GridInterface::Vector& gravity,
			    const PressureSolution& pressure_sol,
			    const Opm::SparseVector<double>& injection_rates) const;

    protected:
	typedef typename GridInterface::CellIterator CIt;
	typedef typename CIt::FaceIterator FIt;
	typedef typename FIt::Vector Vector;

        mutable PressureSolver psolver_;

	void checkAndPossiblyClampSat(std::vector<double>& s) const;

        EulerUpstreamResidual<GridInterface,
                              ReservoirProperties,
                              BoundaryConditions> residual_;

	bool method_viscous_;
	bool method_gravity_;
	bool check_sat_;
	bool clamp_sat_;
        double residual_tolerance_;
        int linsolver_verbosity_;
        int linsolver_type_;
        double update_relaxation_;

    };

} // namespace Opm


#include "ImplicitCapillarity_impl.hpp"



#endif // OPENRS_IMPLICITCAPILLARITY_HEADER
