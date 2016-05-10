/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.

  This file is part of The Open Porous Media Project (OPM).

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

#ifndef INIT_CPGRID_HPP
#define INIT_CPGRID_HPP

#include <dune/grid/CpGrid.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

namespace Opm {

    void initCPGrid(Dune::CpGrid& grid, const parameter::ParameterGroup& param);

}
#endif
