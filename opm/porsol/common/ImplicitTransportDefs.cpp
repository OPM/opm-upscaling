#include <config.h>

#include <opm/porsol/common/ImplicitTransportDefs.hpp>

#include <opm/core/grid.h>

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
