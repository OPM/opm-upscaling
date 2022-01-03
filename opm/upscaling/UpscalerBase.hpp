//===========================================================================
//
// File: UpscalerBase.hpp
//
// Created: Thu Apr 29 10:20:22 2010
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
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

#ifndef OPM_UPSCALERBASE_HEADER
#define OPM_UPSCALERBASE_HEADER

#include <opm/input/eclipse/Deck/Deck.hpp>

#include <opm/common/utility/parameters/ParameterGroup.hpp>

#include <opm/grid/CpGrid.hpp>

#include <opm/porsol/common/GridInterfaceEuler.hpp>
#include <opm/porsol/common/BoundaryConditions.hpp>


namespace Opm
{
    /**
       @brief A base class for upscaling.
       @author Atgeirr F. Rasmussen <atgeirr@sintef.no>
    */
    template <class Traits>
    class UpscalerBase
    {
    protected:
    public:
	// ------- Typedefs  -------
	typedef Dune::CpGrid GridType;
 	enum { Dimension = GridType::dimension };
	typedef GridInterfaceEuler<GridType> GridInterface;
    typedef typename Traits::template ResProp<Dimension>::Type ResProp;

	/// A type for the upscaled permeability.
	typedef typename ResProp::MutablePermTensor permtensor_t;

	enum BoundaryConditionType { Fixed = 0, Linear = 1, Periodic = 2 };

	// ------- Methods -------

	/// Default constructor.
	UpscalerBase();

        virtual ~UpscalerBase() {;} ;

	/// Initializes the upscaler from parameters.
	void init(const Opm::ParameterGroup& param);

	/// Initializes the upscaler from given arguments.
	void init(const Opm::Deck& deck,
                  BoundaryConditionType bctype,
                  double perm_threshold,
                  double residual_tolerance = 1e-8,
                  int linsolver_verbosity = 0,
                  int linsolver_type = 3,
                  bool twodim_hack = false,
                  int linsolver_maxit = 0,
                  double linsolver_prolongate_factor = 1.0,
                  int linsolver_smooth_steps = 1,
                  const double gravity = 0.0);

	/// Access the grid.
	const GridType& grid() const;

        /// Set boundary condition type. This may not be used to swicth
        /// between Periodic and the other types, since the grid is
        /// modified for Periodic conditions.
        void setBoundaryConditionType(BoundaryConditionType type);

        /// Set the permeability of a cell directly. This will override
        /// the permeability that was read from the eclipse file.
        void setPermeability(const int cell_index, const permtensor_t& k);

	/// Does a single-phase upscaling.
	/// @return an upscaled permeability tensor.
	permtensor_t upscaleSinglePhase();

        /// Compute upscaled porosity.
        /// @return total pore volume of all cells divided by total volume.
        double upscalePorosity() const;

        /// Compute upscaled net porosity.
        /// @return total pore volume (with NTG) of all cells divided by total volume.
        double upscaleNetPorosity() const;

        /// Compute upscaled NTG.
        /// @return total net of all cells divided by total volume.
        double upscaleNTG() const;

        /// Compute upscaled SWCR.
        /// @return total irreducible water volume divided by total pore volume
        double upscaleSWCR(const bool NTG) const;

        /// Compute upscaled SOWCR.
        /// @return total irreducible oil volume divided by total pore volume
        double upscaleSOWCR(const bool NTG) const;

    protected:
	// ------- Typedefs and enums -------
	typedef GridInterface::CellIterator                CellIter;
	typedef CellIter::FaceIterator                     FaceIter;
	typedef BasicBoundaryConditions<true, true>             BCs;
        typedef typename Traits::template FlowSolver<GridInterface, BCs>::Type FlowSolver;

	// ------- Methods -------
	template <class FlowSol>
	double computeAverageVelocity(const FlowSol& flow_solution,
				      const int flow_dir,
				      const int pdrop_dir) const;

	double computeDelta(const int flow_dir) const;

        template <class FluidInterface>
        permtensor_t upscaleEffectivePerm(const FluidInterface& fluid);

	virtual void initImpl(const Opm::ParameterGroup& param);

	virtual void initFinal(const Opm::ParameterGroup& param);

	// ------- Data members -------
	BoundaryConditionType bctype_;
	bool twodim_hack_;
	double residual_tolerance_;
        int linsolver_maxit_;
        double linsolver_prolongate_factor_;
	int linsolver_verbosity_;
        int linsolver_type_;
        int linsolver_smooth_steps_;
        double gravity_;

	GridType grid_;
	GridInterface ginterf_;
	ResProp res_prop_;
	BCs bcond_;
	FlowSolver flow_solver_;
    };

} // namespace Opm

#include "UpscalerBase_impl.hpp"




#endif // OPM_UPSCALERBASE_HEADER
