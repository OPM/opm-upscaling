/*****************************************************************************
 *   Copyright (C) 2009-2010 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief A generic template for tabulated material laws that depend
 *        on two parameters.
 *
 * This class is used for rectengular areas of higher resolution and
 * is required to make sure that the material law is continuous.
 */
#ifndef OPM_TABULATED_MATERIAL2_HIRES_HH
#define OPM_TABULATED_MATERIAL2_HIRES_HH

#include "tabulatedmaterial2.hh"

namespace Opm {
/*!
 * \brief A generic template for tabulated material laws that depend
 *        on two parameters.
 *
 * This class is used for rectengular areas of higher resolution and
 * is required to make sure that the material law is continuous.
 *
 * In order to achieve this, we use "transition zones" at the
 * edges of the high resolution area. If we need a value inside
 * such a transition zone, we interpolate between the value of the
 * low-resolution table and the one of the high-resolution table,
 * depending on how far the point is away from the edge.
 */
template <class Traits>
class TabulatedMaterial2HiRes : public TabulatedMaterial2<Traits>
{
    typedef typename Traits::Scalar Scalar;
    typedef TabulatedMaterial2<Traits> ParentType;
    enum { numX = Traits::numX, numY = Traits::numY };

public:
    TabulatedMaterial2HiRes()
    {
    };

    /*!
     * \brief Returns the total weight of the value of the hires
     *        area compared to the one of the lowres area.
     */
    Scalar hiresWeight(Scalar x, Scalar y) const
    {
        return southernWeight_(y)*northernWeight_(y)*
            westernWeight_(x)*easternWeight_(x);
    };

protected:
    /*!
     * \brief Returns the weighting factor of the high resolution
     *        area relative to the low resolution area on the
     *        southern border.
     */
    Scalar southernWeight_(Scalar y) const
    {
        assert(Traits::transitionSouth >= 0);

        // make sure that we are in the transition zone
        if (y <= Traits::yMin) {
            return 0.0;
        }
        if (y >= Traits::yMin + Traits::transitionSouth)
            return 1.0;

        // we need the tWidth variable to avoid a bogus
        // division by zero warning in gcc
        Scalar tWidth = Traits::transitionSouth;
        Scalar relPos = (y - Traits::yMin)/tWidth;

        return transitionFn_(relPos);
    }

    /*!
     * \brief Returns the weighting factor of the high resolution
     *        area relative to the low resolution area on the
     *        northern border.
     */
    Scalar northernWeight_(Scalar y) const
    {
        assert(Traits::transitionNorth >= 0);

        // make sure that we are in the transition zone
        if (y >= Traits::yMax) {
            return 0.0;
        }
        if (y <= Traits::yMax - Traits::transitionNorth)
            return 1.0;

        // we need the tWidth variable to avoid a bogus
        // division by zero warning in gcc
        Scalar tWidth = Traits::transitionNorth;
        Scalar relPos = (Traits::yMax - y)/tWidth;

        return transitionFn_(relPos);
    }

    /*!
     * \brief Returns the weighting factor of the high resolution
     *        area relative to the low resolution area on the
     *        western border.
     */
    Scalar westernWeight_(Scalar x) const
    {
        assert(Traits::transitionWest >= 0);

        // make sure that we are in the transition zone
        if (x <= Traits::xMin) {
            return 0.0;
        }
        if (x >= Traits::xMin + Traits::transitionWest)
            return 1.0;

        // we need the tWidth variable to avoid a bogus
        // division by zero warning in gcc
        Scalar tWidth = Traits::transitionWest;
        Scalar relPos = (x - Traits::xMin)/tWidth;

        return transitionFn_(relPos);
    }

    /*!
     * \brief Returns the weighting factor of the high resolution
     *        area relative to the low resolution area on the
     *        eastern border.
     */
    Scalar easternWeight_(Scalar x) const
    {
        assert(Traits::transitionEast >= 0);

        // make sure that we are in the transition zone
        if (x >= Traits::xMax) {
            return 0.0;
        }
        if (x <= Traits::xMax - Traits::transitionEast)
            return 1.0;

        // we need the tWidth variable to avoid a bogus
        // division by zero warning in gcc
        Scalar tWidth = Traits::transitionEast;
        Scalar relPos = (Traits::xMax - x)/tWidth;

        return transitionFn_(relPos);
    }

    /*!
     * \brief Returns the highres weight given a relative positon in
     *        the transition zone.
     *
     *  We use a quadratic function which quickly decreases the weight
     *  of the lowres table. TODO (?): a spline would probably be a
     *  good idea since it is C1 continous...
     */
    Scalar transitionFn_(Scalar relPos) const
    {  return std::min(1.0,  1 - (1 - relPos)*(1 - relPos)); }

};
}

#endif
