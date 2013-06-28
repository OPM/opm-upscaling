//===========================================================================
//
// File: SteadyStateUpscaler.hpp
//
// Created: Fri Aug 28 14:01:19 2009
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

#ifndef OPM_STEADYSTATEUPSCALERIMPLICIT_HEADER
#define OPM_STEADYSTATEUPSCALERIMPLICIT_HEADER

#include <opm/upscaling/UpscalerBase.hpp>
#include <opm/porsol/euler/EulerUpstream.hpp>
#include <opm/porsol/euler/ImplicitCapillarity.hpp>
#include <dune/grid/common/GridAdapter.hpp>
#include <boost/array.hpp>

namespace Opm
{
    /**
       @brief A class for doing steady state upscaling.
       @author Atgeirr F. Rasmussen <atgeirr@sintef.no>
    */
    template <class Traits>
    class SteadyStateUpscalerImplicit : public UpscalerBase<Traits>
    {
    public:
        // ------- Typedefs and enums -------

        typedef UpscalerBase<Traits> Super;
        typedef typename Super::permtensor_t permtensor_t;
        typedef typename UpscalerBase<Traits>::GridInterface GridInterface;
        typedef typename UpscalerBase<Traits>::GridType GridType;
        enum { Dimension = UpscalerBase<Traits>::Dimension };

        // ------- Methods -------

        /// Default constructor.
        SteadyStateUpscalerImplicit();

        /// Does a steady-state upscaling.
        /// @param flow_direction The cardinal direction in which to impose a pressure gradient for the purpose of converging to steady state.
        /// @param initial_saturation the initial saturation profile for the steady-state computation.
        /// The vector must have size equal to the number of cells in the grid.
        /// @param boundary_saturation the saturation of fluid flowing in across the boundary,
        /// only needed for nonperiodic upscaling.
        /// @param pressure_drop the pressure drop in Pascal over the domain.
        /// @param upscaled_perm typically the output of upscaleSinglePhase().
        /// @return the upscaled relative permeability matrices of both phases.
        /// The relative permeability matrix, call it k, is such that if K_w is the phase
        /// permeability and K the absolute permeability, K_w = k*K.
        std::pair<permtensor_t, permtensor_t> upscaleSteadyState(const int flow_direction,
                                                                 const std::vector<double>& initial_saturation,
                                                                 const double boundary_saturation,
                                                                 const double pressure_drop,
                                                                 const permtensor_t& upscaled_perm,bool& success);

        /// Accessor for the steady state saturation field. This is empty until
        /// upscaleSteadyState() is called, at which point it will
        /// contain the last computed (steady) saturation state.
        const std::vector<double>& lastSaturationState() const;

        /// Computes the upscaled saturation corresponding to the saturation field
        /// returned by lastSaturationState(). Does this by computing total saturated
        /// volume divided by total pore volume.
        double lastSaturationUpscaled() const;

        /// Ensure saturations are not outside table
        void initSatLimits(std::vector<double>& s) const;

        void setToCapillaryLimit(double average_s, std::vector<double>& s) const;


    protected:
        // ------- Typedefs -------
   typedef typename Traits::template TransportSolver<GridInterface, typename Super::BCs>::Type TransportSolver;

        // ------- Methods -------
        template <class FlowSol>
        void computeInOutFlows(std::pair<double, double>& water_inout,
                               std::pair<double, double>& oil_inout,
                               const FlowSol& flow_solution,
                               const std::vector<double>& saturations) const;
        /// Override from superclass.
        virtual void initImpl(const Opm::parameter::ParameterGroup& param);


        // ------- Data members -------
        std::vector<double> last_saturation_state_;
        bool use_gravity_;
        bool output_vtk_;
        bool output_ecl_;
        bool print_inoutflows_;
        int simulation_steps_;
        double init_stepsize_;
        double relperm_threshold_;
        double maximum_mobility_contrast_;
        double sat_change_year_;
        int    max_it_;
        double  max_stepsize_;
        double dt_sat_tol_;
        bool use_maxdiff_;
        TransportSolver transport_solver_;
        GridAdapter grid_adapter_;
    };

} // namespace Opm

#include "SteadyStateUpscalerImplicit_impl.hpp"


#endif // OPM_STEADYSTATEUPSCALERIMPLICIT_HEADER
