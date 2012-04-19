/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser                                    *
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
 */
#ifndef OPM_TABULATED_MATERIAL2_HH
#define OPM_TABULATED_MATERIAL2_HH

#include <dune/porsol/blackoil/co2fluid/opm/common/exceptions.hh>

#include <assert.h>

namespace Opm {
/*!
 * \brief A generic template for tabulated material laws that depend
 *        on two parameters.
 */
template <class Traits>
class TabulatedMaterial2
{
    typedef typename Traits::Scalar Scalar;
    enum { numX = Traits::numX, numY = Traits::numY };

public:
    TabulatedMaterial2()
    {
    };

    Scalar xMin() const
    { return Traits::xMin; }

    Scalar xMax() const
    { return Traits::xMax; }

    Scalar yMin() const
    { return Traits::yMin; }

    Scalar yMax() const
    { return Traits::yMax; }

    bool applies(Scalar x, Scalar y) const
    {
        return xMin() <= x && x <= xMax() &&
            yMin() <= y && y <= yMax();
    }

    Scalar at(Scalar x, Scalar y) const
    {
#if !defined NDEBUG || defined CHECK_TABLE_RANGE
        if (!applies(x,y))
        {
            std::cout << "*** TabulatedMaterial2():  outside table range " << x << " , " << y << " ****" << std::endl;
            x = std::min(std::max(x,xMin()),xMax());
            y = std::min(std::max(y,yMin()),yMax());
            /*
            OPM_THROW(Opm::NumericalProblem,
                       "Attempt to get tabulated " << Traits::name << " value for ("
                       << x << ", " << y
                       << ") on a table of extend "
                       << xMin() << " to " << xMax() << " times "
                       << yMin() << " to " << yMax());
*/
        };
#endif

        // use higher resolution if available
        Scalar hiresValue = 0.0;
        Scalar hiresWeight = 0.0;
        if (Traits::hires.applies(x,y)) {
            hiresValue = Traits::hires.at(x, y);
            hiresWeight = Traits::hires.hiresWeight(x, y);
            assert(hiresWeight >= 0.0 && hiresWeight <= 1.0);

            if (hiresWeight == 1.0)
                return hiresValue;
        }

        int i = findXIdx_(x);
        int j = findYIdx_(y);

        Scalar xAtI = xAt_(i);
        Scalar xAtI1 = xAt_(i + 1);
        Scalar yAtJ = yAt_(j);
        Scalar yAtJ1 = yAt_(j + 1);

        Scalar alpha = (x - xAtI)/(xAtI1 - xAtI);
        Scalar beta = (y - yAtJ)/(yAtJ1 - yAtJ);

        // bi-linear interpolation
        Scalar lowresValue =
            (1-alpha)*(1-beta)*val(i, j) +
            (1-alpha)*(  beta)*val(i, j + 1) +
            (  alpha)*(1-beta)*val(i + 1, j) +
            (  alpha)*(  beta)*val(i + 1, j + 1);

        // return the weighted sum of the low- and high-resolution
        // values
        return hiresValue*hiresWeight + lowresValue*(1 - hiresWeight);
    };

    Scalar val(int i, int j) const
    {
#if !defined NDEBUG || defined CHECK_TABLE_RANGE
        if (i < 0 || i >= Traits::numX ||
            j < 0 || j >= Traits::numY) {
            OPM_THROW(Opm::NumericalProblem,
                       "Attempt to access element ("
                       << i << ", " << j
                       << ") on a " << Traits::name << " table of size ("
                       << Traits::numX << ", " << Traits::numY
                       << ")\n");
        };
#endif
        return Traits::vals[i][j];
    };

protected:
    int findXIdx_(Scalar x) const
    {
        if (x == xMax())
            return numX - 2;
        const int result = static_cast<int>((x - xMin())/(xMax() - xMin())*(numX - 1));
        return std::max(0, std::min(result, numX - 2));
    };

    int findYIdx_(Scalar y) const
    {
        if (y == yMax())
            return numY - 2;
        const int result = static_cast<int>((y - yMin())/(yMax() - yMin())*(numY - 1));
        return std::max(0, std::min(result, numY - 2));
    };

    Scalar xAt_(int i) const
    { return i*(xMax() - xMin())/(numX - 1) + xMin(); }
    Scalar yAt_(int j) const
    { return j*(yMax() - yMin())/(numY - 1) + yMin(); }
};
}

#endif
