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

#include <dune/common/Units.hpp>
#include <dune/porsol/blackoil/BlackoilWells.hpp>
#include <dune/porsol/common/Rock.hpp>


int main(int argc, char** argv)
{
    // Parser.
    //const std::string ecl_file("SPE9.DATA");
    const std::string ecl_file("SPJ8.DATA");
    Dune::EclipseGridParser parser(ecl_file);
    if (parser.hasField("WELSPECS")) {
	parser.getWELSPECS().write(std::cout);    
    }
    if (parser.hasField("COMPDAT")) {
	parser.getCOMPDAT().write(std::cout);    
    }
    if (parser.hasField("WCONINJE")) {
	parser.getWCONINJE().write(std::cout);    
    }
    if (parser.hasField("WCONPROD")) {
	parser.getWCONPROD().write(std::cout);
    }
    if (parser.hasField("WELTARG")) {
	parser.getWELTARG().write(std::cout);    
    }

    // Grid
    Dune::CpGrid grid;
    Dune::array<int, 3> dims;
    Dune::array<double, 3> cellsize;
    std::fill(dims.begin(), dims.end(), 3);
    cellsize[0] = 2.0;
    cellsize[1] = 3.0;
    cellsize[2] = 4.0;
    grid.createCartesian(dims, cellsize);
    
    // Rock
    Dune::Rock<3> rock;
    rock.init(grid.numCells(), 1.0,
	      100.0*Dune::prefix::milli*Dune::unit::darcy);

    // Test the BlackoilWells class.
    Opm::BlackoilWells wells;
    wells.init(parser, grid, rock);
}
