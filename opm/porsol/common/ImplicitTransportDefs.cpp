/*
  Copyright 2010 Statoil ASA.

  This file is part of The Open Porous Media project (OPM).

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

#include <config.h>

#include <opm/porsol/common/ImplicitTransportDefs.hpp>

#include <opm/grid/UnstructuredGrid.h>

#include <algorithm>
#include <cassert>

void
Opm::compute_porevolume(const UnstructuredGrid* g,
                        const Rock&             rock,
                        std::vector<double>&    porevol)
{
    const auto& poro = rock.poro();

    assert (poro.size() == (::std::size_t)(g->number_of_cells));

    porevol.resize(rock.poro().size());

    ::std::transform(poro.begin(), poro.end(),
                     g->cell_volumes,
                     porevol.begin(),
                     ::std::multiplies<double>());
}
