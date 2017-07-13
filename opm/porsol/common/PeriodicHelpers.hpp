//===========================================================================
//
// File: PeriodicHelpers.hpp
//
// Created: Fri Aug 14 10:59:57 2009
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

#ifndef OPENRS_PERIODICHELPERS_HEADER
#define OPENRS_PERIODICHELPERS_HEADER

#include "BoundaryPeriodicity.hpp"
#include "BoundaryConditions.hpp"
#include <dune/common/fvector.hh>
#include <array>
#include <algorithm>
#include <iostream>

namespace Opm
{





    template <class BCs>
    void storeFlowCond(BCs& bcs,
		       const std::vector<BoundaryFaceInfo>& bfinfo,
		       const std::array<FlowBC, 6>& fconditions,
		       const std::array<double, 6>& side_areas)
    {
	int num_bdy = bfinfo.size();
	for (int i = 0; i < num_bdy; ++i) {
	    FlowBC face_bc = fconditions[bfinfo[i].canon_pos];
	    if (face_bc.isNeumann()) {
		face_bc = FlowBC(FlowBC::Neumann, face_bc.outflux()*bfinfo[i].area/side_areas[bfinfo[i].canon_pos]);
	    }
	    bcs.flowCond(bfinfo[i].bid) = face_bc;
	}
    }





    template <class BCs>
    void storeSatCond(BCs& bcs,
		      const std::vector<BoundaryFaceInfo>& bfinfo,
		      const std::array<SatBC, 6>& sconditions)
    {
	int num_bdy = bfinfo.size();
	for (int i = 0; i < num_bdy; ++i) {
	    bcs.satCond(bfinfo[i].bid) = sconditions[bfinfo[i].canon_pos];
	}
    }


    template <class BC>
    std::array<bool, 6> extractPeriodic(const std::array<BC, 6>& bcs)
    {
	std::array<bool, 6> retval = {{ bcs[0].isPeriodic(),
					  bcs[1].isPeriodic(),
					  bcs[2].isPeriodic(),
					  bcs[3].isPeriodic(),
					  bcs[4].isPeriodic(),
					  bcs[5].isPeriodic() }};
	return retval;
    }


    /// @brief Makes a boundary condition object representing
    /// periodic boundary conditions in any cartesian directions.
    /// The grid interface needs to export boundary ids that are
    /// unique for each boundary face for this to be possible.
    /// @tparam BCs the boundary condition class
    /// @tparam GridInterface grid interface class
    template <class BCs, class GridInterface>
    void createPeriodic(BCs& fbcs,
			const GridInterface& g,
			const std::array<FlowBC, 2*GridInterface::Dimension>& fconditions,
			const std::array<SatBC, 2*GridInterface::Dimension>& sconditions,
			double spatial_tolerance = 1e-6)
    {
	static_assert(BCs::HasFlowConds, "");
	static_assert(BCs::HasSatConds, "");
	// Check the conditions given.
	for (int i = 0; i < GridInterface::Dimension; ++i) {
	    if (fconditions[2*i].isPeriodic()) {
		assert(fconditions[2*i + 1].isPeriodic());
		assert(fconditions[2*i].pressureDifference() == -fconditions[2*i + 1].pressureDifference());
		assert(sconditions[2*i].isPeriodic());
		assert(sconditions[2*i + 1].isPeriodic());
		assert(sconditions[2*i].saturationDifference() == 0.0);
		assert(sconditions[2*i + 1].saturationDifference() == 0.0);
	    }
	}
	std::vector<BoundaryFaceInfo> bfinfo;
	std::array<double, 6> side_areas;
	createPeriodicImpl(fbcs, bfinfo, side_areas, g, extractPeriodic(fconditions), spatial_tolerance);
	storeFlowCond(fbcs, bfinfo, fconditions, side_areas);
	storeSatCond(fbcs, bfinfo, sconditions);
    }




    template <class BCs, class GridInterface>
    void createPeriodic(BCs& fbcs,
			const GridInterface& g,
			const std::array<FlowBC, 2*GridInterface::Dimension>& fconditions,
			double spatial_tolerance = 1e-6)
    {
	static_assert(BCs::HasFlowConds, "");
	static_assert(!BCs::HasSatConds, "");
	// Check the conditions given.
	for (int i = 0; i < GridInterface::Dimension; ++i) {
	    if (fconditions[2*i].isPeriodic()) {
		assert(fconditions[2*i + 1].isPeriodic());
		assert(fconditions[2*i].pressureDifference() == -fconditions[2*i + 1].pressureDifference());
	    }
	}
	std::vector<BoundaryFaceInfo> bfinfo;
	std::array<double, 6> side_areas;
	createPeriodicImpl(fbcs, bfinfo, side_areas, g, extractPeriodic(fconditions), spatial_tolerance);
	storeFlowCond(fbcs, bfinfo, fconditions, side_areas);
    }




    template <class BCs, class GridInterface>
    void createPeriodic(BCs& fbcs,
			const GridInterface& g,
			const std::array<SatBC, 2*GridInterface::Dimension>& sconditions,
			double spatial_tolerance = 1e-6)
    {
	static_assert(!BCs::HasFlowConds, "");
	static_assert(BCs::HasSatConds, "");
	// Check the conditions given.
	for (int i = 0; i < GridInterface::Dimension; ++i) {
	    if (sconditions[2*i].isPeriodic()) {
		assert(sconditions[2*i + 1].isPeriodic());
		assert(sconditions[2*i].saturationDifference() == -sconditions[2*i + 1].saturationDifference());
	    }
	}
	std::vector<BoundaryFaceInfo> bfinfo;
	std::array<double, 6> side_areas;
	createPeriodicImpl(fbcs, bfinfo, side_areas, g, extractPeriodic(sconditions), spatial_tolerance);
	storeSatCond(fbcs, bfinfo, sconditions);
    }



    /// @brief Common implementation for the various createPeriodic functions.
    template <class BCs, class GridInterface>
    void createPeriodicImpl(BCs& fbcs,
			    std::vector<BoundaryFaceInfo>& bfinfo,
			    std::array<double, 6>& side_areas,
			    const GridInterface& g,
			    const std::array<bool, 2*GridInterface::Dimension>& is_periodic,
			    double spatial_tolerance = 1e-6)
    {
        findPeriodicPartners(bfinfo, side_areas, g.grid().leafGridView(), is_periodic, spatial_tolerance);
        int num_bdy = bfinfo.size();
        // This will likely change with boundarySegmentIndex() instead of boundaryId():
        int max_bid = num_bdy;

	// Resize the conditions object. We clear it first to make sure it's all defaults.
	fbcs.clear();
	fbcs.resize(max_bid + 1);

	// Now we have all the info, and the object to store it in. So we store it...
	for (int i = 0; i < num_bdy; ++i) {
	    int bid1 = bfinfo[i].bid;
	    int bid2 = bfinfo[i].partner_bid;
	    if (bid1 < bid2) {
		// If there is no periodic partner, bid2 will be zero, so we will not end up here.
		fbcs.setPeriodicPartners(bid1, bid2);
	    }
	    fbcs.setCanonicalBoundaryId(bid1, bfinfo[i].canon_pos + 1);
	}
    }



    /// @brief Makes a boundary condition object representing
    /// linear boundary conditions in any cartesian direction.
    /// The grid interface needs to export boundary ids that are
    /// unique for each boundary face for this to be possible.
    /// @tparam BCs the boundary condition class
    /// @tparam GridInterface grid interface class
    template <class BCs, class GridInterface>
    void createLinear(BCs& fbcs,
                      const GridInterface& g,
                      const double pdrop,
                      const int pddir,
                      const double bdy_sat,
                      const bool twodim_hack = false,
                      const double spatial_tolerance = 1e-6)
    {
        // NOTE: Section copied from createPeriodicImpl().
	// Pick out all boundary faces, simultaneously find the
	// bounding box of their centroids, and the max id.
	typedef typename GridInterface::CellIterator CI;
	typedef typename CI::FaceIterator FI;
	typedef typename GridInterface::Vector Vector;
        std::vector<FI> bface_iters;
	Vector low(1e100);
	Vector hi(-1e100);
	int max_bid = 0;
	for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
	    for (FI f = c->facebegin(); f != c->faceend(); ++f) {
		if (f->boundaryId()) {
                    bface_iters.push_back(f);
		    Vector fcent = f->centroid();
		    for (int dd = 0; dd < GridInterface::Dimension; ++dd) {
			low[dd] = std::min(low[dd], fcent[dd]);
			hi[dd] = std::max(hi[dd], fcent[dd]);
		    }
		    max_bid = std::max(max_bid, f->boundaryId());
		}
	    }
	}
        int num_bdy = bface_iters.size();
	if (max_bid != num_bdy) {
	    OPM_THROW(std::runtime_error, "createLinear() assumes that every boundary face has a unique boundary id. That seems to be violated.");
	}
        fbcs.resize(max_bid + 1);

        // Iterate over boundary faces, setting boundary conditions for their corresponding ids.
        double cmin = low[pddir];
        double cmax = hi[pddir];
        double cdelta = cmax - cmin;
        
        for (int i = 0; i < num_bdy; ++i) {
            Vector fcent = bface_iters[i]->centroid();
            int canon_pos = -1;
	    for (int dd = 0; dd < GridInterface::Dimension; ++dd) {
		double coord = fcent[dd];
		if (fabs(coord - low[dd]) <= spatial_tolerance) {
		    canon_pos = 2*dd;
		    break;
		} else if (fabs(coord - hi[dd]) <= spatial_tolerance) {
		    canon_pos = 2*dd + 1;
		    break;
		}
	    }
	    if (canon_pos == -1) {
		std::cerr << "Centroid: " << fcent << "\n";
		std::cerr << "Bounding box min: " << low << "\n";
		std::cerr << "Bounding box max: " << hi << "\n";
		OPM_THROW(std::runtime_error, "Boundary face centroid not on bounding box. Maybe the grid is not an axis-aligned shoe-box?");
	    }
            double relevant_coord = fcent[pddir];
            double pressure = pdrop*(1.0 - (relevant_coord - cmin)/cdelta);
            int bid = bface_iters[i]->boundaryId();
            fbcs.setCanonicalBoundaryId(bid, canon_pos + 1);
            fbcs.satCond(bid) = SatBC(SatBC::Dirichlet, bdy_sat);
            fbcs.flowCond(bid) = FlowBC(FlowBC::Dirichlet, pressure);
            if (twodim_hack && canon_pos >= 4) {
                // Noflow for Z- and Z+ boundaries in the 3d-grid-that-is-really-2d case.
                fbcs.flowCond(bid) = FlowBC(FlowBC::Neumann, 0.0);
            }
        }
    }


} // namespace Opm

#endif // OPENRS_PERIODICHELPERS_HEADER
