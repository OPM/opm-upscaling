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

#include <dune/porsol/blackoil/BlackoilUpstreamEuler.hpp>
#include <dune/porsol/common/SimpleRock.hpp>
#include <dune/porsol/blackoil/OilGasFluid.hpp>
#include <dune/grid/CpGrid.hpp>
#include <vector>


// Initialize the variables to some valid state.
template <class Grid, class Rock, class Fluid>
void initializeBlackoilState(const Grid& grid,
                             const Rock& rock,
                             const Fluid& fluid,
                             std::vector<typename Fluid::ComponentVec>& z,
                             std::vector<double>& cell_pressure)
{
}


int main(int argc, char** argv)
{
    // Parameters.
    Dune::parameter::ParameterGroup param(argc, argv);

    // Initialization.
    typedef Dune::CpGrid Grid;
    typedef Opm::SimpleRock Rock;
    typedef Opm::OilGasFluid Fluid;
    typedef Opm::BlackoilUpstreamEuler<Grid, Fluid> Transport;
    Grid grid;
    grid.init(param);
    Rock rock;
    rock.init(param);
    Fluid fluid;
    fluid.init(param);
    Transport transport;
    transport.init(param);

    // State variables: z and pressure.
    int num_cells = grid.size(0);
    std::vector<Fluid::ComponentVec> z(num_cells);
    std::vector<double> cell_pressure(num_cells);
    initializeBlackoilState(grid, rock, fluid, z, cell_pressure);
}
