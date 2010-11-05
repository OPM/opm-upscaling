/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_OILGASFLUID_HEADER_INCLUDED
#define OPM_OILGASFLUID_HEADER_INCLUDED


#include <dune/common/param/ParameterGroup.hpp>
#include <dune/common/fvector.hh>

namespace Opm
{

    class OilGasFluid
    {
    public:
        enum { NumPhases = 2 };
        enum { NumComponents = 2 };
        typedef Dune::FieldVector<double, NumPhases> PhaseVec;
        typedef Dune::FieldVector<double, NumPhases> ComponentVec;

        void init(const Dune::parameter::ParameterGroup& param)
        {
        }
    };


} // namespace Opm


#endif // OPM_OILGASFLUID_HEADER_INCLUDED
