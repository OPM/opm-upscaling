/*
  Copyright 2011 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_BOUNDARYPERIODICITY_HEADER_INCLUDED
#define OPM_BOUNDARYPERIODICITY_HEADER_INCLUDED


#include <dune/common/fvector.hh>
#include <vector>

namespace Dune
{

    /// @brief
    /// @todo Doc me!
    struct BoundaryFaceInfo
    {
	/// Face index in [0, ..., #faces - 1]
	int face_index;
        /// Boundary id of this face.
	int bid;
        /// Canonical position of face:
        ///  0 -> xmin
        ///  1 -> xmax
        ///  2 -> ymin
        ///  3 -> ymax
        /// ...
	int canon_pos;
        /// Face index of periodic partner, or -1 if no partner.
	int partner_face_index;
        /// Boundary id of periodic partner face, or 0 if no parner.
	int partner_bid;
        /// Face area.
	double area;
        /// Face centroid.
	FieldVector<double,3> centroid;

        /// Comparison operator.
        /// Intended to make periodic partners appear close in
        /// sorted order, but this is only a heuristic.
	bool operator<(const BoundaryFaceInfo& other) const
	{
	    return cmpval() < other.cmpval();
	}

    private:
	double cmpval() const
	{
            const double pi = 3.14159265358979323846264338327950288;
	    return centroid[(canon_pos/2 + 1)%3] + pi*centroid[(canon_pos/2 + 2)%3];
	}
    };


    /// @brief Find a match (periodic partner) for the given face.
    /// @param[inout] bfaces the boundary face info list.
    /// @param[in] face the face for which we seek a periodic partner
    /// @param[in] lower lower end of search interval [lower, upper)
    /// @param[in] upper upper end of search interval [lower, upper)
    /// @return true if a match was found, otherwise false
    bool match(std::vector<BoundaryFaceInfo>& bfaces, int face, int lower, int upper)
    {
	const double area_tol = 1e-6;
	const double centroid_tol = 1e-6;
	int cp = bfaces[face].canon_pos;
	int target_cp = (cp%2 == 0) ? cp + 1 : cp - 1;
	FieldVector<double, 3> cent_this = bfaces[face].centroid;
	for (int j = lower; j < upper; ++j) {
	    if (bfaces[j].canon_pos == target_cp) {
		if (fabs(bfaces[face].area - bfaces[j].area) <= area_tol) {
		    FieldVector<double, 3> cent_other = bfaces[j].centroid;
		    cent_other -= cent_this;
		    double dist = cent_other.two_norm();
		    if (dist <= centroid_tol) {
			bfaces[face].partner_face_index = bfaces[j].face_index;
			bfaces[face].partner_bid = bfaces[j].bid;
			bfaces[j].partner_face_index = bfaces[face].face_index;
			bfaces[j].partner_bid = bfaces[face].bid;
			break;
		    }
		}
	    }
	}
	return (bfaces[face].partner_face_index != -1);
    }

} // namespace Dune

#endif // OPM_BOUNDARYPERIODICITY_HEADER_INCLUDED
