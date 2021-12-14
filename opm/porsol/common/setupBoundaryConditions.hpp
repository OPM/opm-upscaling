//===========================================================================
//
// File: setupBoundaryConditions.hpp
//
// Created: Fri Aug 21 09:07:09 2009
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

#ifndef OPENRS_SETUPBOUNDARYCONDITIONS_HEADER
#define OPENRS_SETUPBOUNDARYCONDITIONS_HEADER

#include <opm/common/utility/parameters/ParameterGroup.hpp>
#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/porsol/common/BoundaryConditions.hpp>
#include <opm/porsol/common/PeriodicHelpers.hpp>

namespace Opm
{

    /// @brief Setup boundary conditions for a simulation.
    /// It is assumed that the boundary ids are 1-6, similar to cartesian case/Yaspgrid,
    /// unless periodic, in which case we assume unique boundary ids.
    template <class GridInterface, class BCs>
    inline void setupBoundaryConditions(const Opm::ParameterGroup& param,
					const GridInterface& g,
					BCs& bcs)
    {
	if (param.getDefault("upscaling", false)) {
	    int bct = param.get<int>("boundary_condition_type");
	    int pddir = param.getDefault("pressure_drop_direction", 0);
	    double pdrop = param.getDefault("boundary_pressuredrop", 1.0e5);
	    double bdy_sat = param.getDefault("boundary_saturation", 1.0);
	    bool twodim_hack = param.getDefault("2d_hack", false);
	    setupUpscalingConditions(g, bct, pddir, pdrop, bdy_sat, twodim_hack, bcs);
	    return;
	}
        if (param.getDefault("region_based_bcs", false)) {
            setupRegionBasedConditions(param, g, bcs);
            return;
        }
	// Make flow equation boundary conditions.
	// Default is pressure 1.0e5 on the left, 0.0 on the right.
	// Recall that the boundary ids range from 1 to 6 for the cartesian edges,
	// and that boundary id 0 means interiour face/intersection.
	std::string flow_bc_type = param.getDefault<std::string>("flow_bc_type", "dirichlet");
	FlowBC::BCType bct = FlowBC::Dirichlet;
	double leftval = 1.0*Opm::unit::barsa;
	double rightval = 0.0;
	if (flow_bc_type == "neumann") {
	    bct = FlowBC::Neumann;
	    leftval = param.get<double>("left_flux");
	    rightval = param.getDefault<double>("right_flux", -leftval);
	} else if (flow_bc_type == "dirichlet") {
	    leftval = param.getDefault<double>("left_pressure", leftval);
	    rightval = param.getDefault<double>("right_pressure", rightval);
	} else if (flow_bc_type == "periodic") {
	    OPM_THROW(std::runtime_error, "Periodic conditions not here yet.");
	} else {
	    OPM_THROW(std::runtime_error, "Unknown flow boundary condition type " << flow_bc_type);
	}
	bcs.resize(7);
	bcs.flowCond(1) = FlowBC(bct, leftval);
	bcs.flowCond(2) = FlowBC(bct, rightval);

	// Default transport boundary conditions are used.
    }

    /// @brief
    /// @todo Doc me!
    /// @param
    template <class GridInterface, class BCs>
    inline void setupUpscalingConditions(const GridInterface& g,
					 int bct,
					 int pddir,
					 double pdrop,
					 double bdy_sat,
					 bool twodim_hack,
					 BCs& bcs)
    {
	// Caution: This enum is copied from Upscaler.hpp.
	enum BoundaryConditionType { Fixed = 0, Linear = 1, Periodic = 2, PeriodicSingleDirection = 3, Noflow = 4 };
        if (bct < 0 || bct > 2) {
            OPM_THROW(std::runtime_error, "Illegal boundary condition type (0-2 are legal): " << bct); // Later on, we may allow 3 and 4.
        }
	BoundaryConditionType bctype = static_cast<BoundaryConditionType>(bct);
        assert(pddir >=0 && pddir <= 2);

	// Flow conditions.
	switch (bctype) {
	case Fixed:
	    {
		// assert(!g.uniqueBoundaryIds());
		bcs.clear();
		bcs.resize(7);
		bcs.flowCond(2*pddir + 1) = FlowBC(FlowBC::Dirichlet, pdrop);
		bcs.flowCond(2*pddir + 2) = FlowBC(FlowBC::Dirichlet, 0.0);
		bcs.satCond(2*pddir + 1) = SatBC(SatBC::Dirichlet, bdy_sat); // The only possible inflow location.
		for (int i = 0; i < 7; ++i) {
		    bcs.setCanonicalBoundaryId(i, i);
		}
		break;
	    }
	case Linear:
	    {
		// assert(g.uniqueBoundaryIds());
		createLinear(bcs, g, pdrop, pddir, bdy_sat, twodim_hack);
		break;
	    }
	case Periodic:
	    {
		// assert(g.uniqueBoundaryIds());
		FlowBC fb(FlowBC::Periodic, 0.0);
		std::array<FlowBC, 6> fcond = {{ fb, fb, fb, fb, fb, fb }};
		fcond[2*pddir] = FlowBC(FlowBC::Periodic, pdrop);
		fcond[2*pddir + 1] = FlowBC(FlowBC::Periodic, -pdrop);
		SatBC sb(SatBC::Periodic, 0.0);
		std::array<SatBC, 6> scond = {{ sb, sb, sb, sb, sb, sb }};
		if (twodim_hack) {
// 		    fcond[2] = FlowBC(FlowBC::Neumann, 0.0);
// 		    fcond[3] = FlowBC(FlowBC::Neumann, 0.0);
		    fcond[4] = FlowBC(FlowBC::Neumann, 0.0);
		    fcond[5] = FlowBC(FlowBC::Neumann, 0.0);
// 		    scond[2] = SatBC(SatBC::Dirichlet, 1.0);
// 		    scond[3] = SatBC(SatBC::Dirichlet, 1.0);
		    scond[4] = SatBC(SatBC::Dirichlet, 1.0);
		    scond[5] = SatBC(SatBC::Dirichlet, 1.0);
		}
		createPeriodic(bcs, g, fcond, scond);
		break;
	    }
	default:
	    OPM_THROW(std::runtime_error, "Error in switch statement, should never be here.");
	}

	// Default transport boundary conditions are used.
    }



    namespace
    {
        template <class Vector>
        bool isInside(const Vector& low, const Vector& high, const Vector& pt)
        {
            return low[0] < pt[0]
                && low[1] < pt[1]
                && low[2] < pt[2]
                && high[0] > pt[0]
                && high[1] > pt[1]
                && high[2] > pt[2];
        }
    } // anon namespace


    template <class GridInterface, class BCs>
    inline void setupRegionBasedConditions(const Opm::ParameterGroup& param,
                                           const GridInterface& g,
                                           BCs& bcs)
    {
        // Extract region and pressure value for Dirichlet bcs.
        typedef typename GridInterface::Vector Vector;
        Vector low;
        low[0] = param.getDefault("dir_block_low_x", 0.0);
        low[1] = param.getDefault("dir_block_low_y", 0.0);
        low[2] = param.getDefault("dir_block_low_z", 0.0);
        Vector high;
        high[0] = param.getDefault("dir_block_high_x", 1.0);
        high[1] = param.getDefault("dir_block_high_y", 1.0);
        high[2] = param.getDefault("dir_block_high_z", 1.0);
        double dir_block_pressure = param.get<double>("dir_block_pressure");

        // Set flow conditions for that region.
        // For this to work correctly, unique boundary ids should be used,
        // otherwise conditions may spread outside the given region, to all
        // faces with the same bid as faces inside the region.
        typedef typename GridInterface::CellIterator CI;
        typedef typename CI::FaceIterator FI;
        int max_bid = 0;
        std::vector<int> dir_bids;
        for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
            for (FI f = c->facebegin(); f != c->faceend(); ++f) {
                int bid = f->boundaryId();
                max_bid = std::max(bid, max_bid);
                if (bid != 0 && isInside(low, high, f->centroid())) {
                    dir_bids.push_back(bid);
                }
            }
        }
        bcs.resize(max_bid + 1);
        for (std::vector<int>::const_iterator it = dir_bids.begin(); it != dir_bids.end(); ++it) {
            bcs.flowCond(*it) = FlowBC(FlowBC::Dirichlet, dir_block_pressure);
        }

        // Transport BCs are defaulted.
    }



} // namespace Opm


#endif // OPENRS_SETUPBOUNDARYCONDITIONS_HEADER
