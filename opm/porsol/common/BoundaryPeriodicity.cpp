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


#include <opm/porsol/common/BoundaryPeriodicity.hpp>

namespace Opm {

bool match(std::vector<BoundaryFaceInfo>& bfaces, int face, int lower, int upper)
{
    const double area_tol = 1e-6;
    const double centroid_tol = 1e-6;
    int cp = bfaces[face].canon_pos;
    int target_cp = (cp%2 == 0) ? cp + 1 : cp - 1;
    Dune::FieldVector<double, 3> cent_this = bfaces[face].centroid;
    for (int j = lower; j < upper; ++j) {
        if (bfaces[j].canon_pos == target_cp) {
            if (fabs(bfaces[face].area - bfaces[j].area) <= area_tol) {
                Dune::FieldVector<double, 3> cent_other = bfaces[j].centroid;
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

}
