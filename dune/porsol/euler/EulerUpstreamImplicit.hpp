//===========================================================================
//
// File: EulerUpstreamImplicit	.hpp
//
// Created: Tue Jun 16 14:07:53 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            B�rd Skaflestad     <bard.skaflestad@sintef.no>
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

#ifndef OPENRS_EULERUPSTREAMIMPLICIT_HEADER
#define OPENRS_EULERUPSTREAMIMPLICIT_HEADER

//#include <dune/porsol/euler/EulerUpstreamImplicit	Residual.hpp>
#include <dune/porsol/opmpressure/src/GridAdapter.hpp>
#include <tr1/unordered_map>
#include <dune/common/param/ParameterGroup.hpp>
#include <dune/porsol/opmtransport/src/CSRMatrixUmfpackSolver.hpp>
#include <dune/porsol/opmtransport/src/NormSupport.hpp>
#include <dune/porsol/opmtransport/src/ImplicitAssembly.hpp>
#include <dune/porsol/opmtransport/src/ImplicitTransport.hpp>
#include <dune/porsol/opmtransport/src/JacobianSystem.hpp>

#include <dune/common/SparseVector.hpp>
#include <dune/porsol/opmtransport/examples/ImplicitTransportDefs.hpp>



#include <dune/porsol/opmtransport/src/CSRMatrixBlockAssembler.hpp>

#include <dune/porsol/opmtransport/src/SimpleFluid2pWrapper.hpp>
#include <dune/porsol/opmtransport/src/SinglePointUpwindTwoPhase.hpp>

// temporatly
#include <array>

namespace Dune {

    /// Class for doing simple transport by explicit Euler upstream method for general grid.
    /// @tparam
    template <class GridInterface, class ReservoirProperties, class BoundaryConditions>
    class EulerUpstreamImplicit	
    {
    public:
	/// @brief
	/// @todo Doc me
	EulerUpstreamImplicit	();
	/// @brief
	/// @todo Doc me
	/// @param
 	EulerUpstreamImplicit	(const GridInterface& grid,
 		      const ReservoirProperties& resprop,
		      const BoundaryConditions& boundary);
	/// @brief
	/// @todo Doc me
	/// @param
	void init(const parameter::ParameterGroup& param);
	/// @brief
	/// @todo Doc me
	/// @param
	void init(const parameter::ParameterGroup& param,
		  const GridInterface& grid,
		  const ReservoirProperties& resprop,
		  const BoundaryConditions& boundary);
	/// @brief
	/// @todo Doc me
	/// @param
	void initObj(const GridInterface& grid,
		     const ReservoirProperties& resprop,
		     const BoundaryConditions& boundary);
	/// @brief
	/// @todo Doc me
	/// @param
	void display();

	/// \brief Set the Courant number.
	/// That is dt = dt_cfl*courant_number.
	/// For this explicit method it should be < 1.
	void setCourantNumber(double cn);

	/// \brief Solve transport equation, evolving \param saturation
	/// for \param time seconds.
	/// Cfl type conditions may force many explicit timesteps to
	/// be taken, before the function returns.
	/// @tparam
	/// @param
	template <class PressureSolution>
	void transportSolve(std::vector<double>& saturation,
			    const double time,
			    const typename GridInterface::Vector& gravity,
			    const PressureSolution& pressure_sol,
			    const SparseVector<double>& injection_rates) const;

    protected:
	// typedef typename GridInterface::CellIterator CIt;
	// typedef typename CIt::FaceIterator FIt;
	// typedef typename FIt::Vector Vector;
     typedef ReservoirProperties RP;

	/*
	template <class PressureSolution>
	double computeCflTime(const std::vector<double>& saturation,
			      const double time,
			      const typename GridInterface::Vector& gravity,
			      const PressureSolution& pressure_sol) const;
	*/
	template <class PressureSolution>
	void smallTimeStep(std::vector<double>& saturation,
			   const double time,
			   const typename GridInterface::Vector& gravity,
			   const PressureSolution& pressure_sol,
                           const SparseVector<double>& injection_rates) const;

	void checkAndPossiblyClampSat(std::vector<double>& s) const;


	
	// try to make bard ingerfaces
	typedef Opm::SimpleFluid2pWrapper<ReservoirProperties>                          TwophaseFluid;
	typedef Opm::SinglePointUpwindTwoPhase<TwophaseFluid> TransportModel;
	// using namespace Opm::ImplicitTransportDefault
	typedef Opm::ImplicitTransportDefault::NewtonVectorCollection< ::std::vector<double> >      NVecColl;
	typedef Opm::ImplicitTransportDefault::JacobianSystem        < struct CSRMatrix, NVecColl > JacSys;
	
	typedef Opm::ImplicitTransport<TransportModel,
                               JacSys        ,
                               MaxNorm       ,
                               Opm::ImplicitTransportDefault::VectorNegater ,
                               Opm::ImplicitTransportDefault::VectorZero    ,
                               Opm::ImplicitTransportDefault::MatrixZero    > TransportSolver;
	// should be initialized by param
	Opm::ImplicitTransportDetails::NRReport  rpt_;
	Opm::ImplicitTransportDetails::NRControl ctrl_;
	Opm::ImplicitTransportLinAlgSupport::CSRMatrixUmfpackSolver linsolve_;
	
	TransportModel  model_;
	TransportSolver tsolver_;
	GridAdapter mygrid_;

	Opm::SimpleFluid2pWrapper< ReservoirProperties > myfluid_;

	bool check_sat_;
	bool clamp_sat_;
    std::vector<double> porevol_;
    std::vector<double> faceflux_;
    std::vector<double> dunefaceind_;
    std::vector<double> dunehfacesign_;

	// Storing residual so that we won't have to reallocate it for every step.
	//mutable std::vector<double> residual_;
    };

} // namespace Dune

#include "EulerUpstreamImplicit_impl.hpp"

#endif // OPENRS_EULERUPSTREAM_HEADER
