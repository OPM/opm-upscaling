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



#include <array>
#include <string>

#include "config.h"
#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <dune/common/version.hh>

#include <dune/common/parallel/mpihelper.hh>

#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/common/ErrorMacros.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <opm/grid/CpGrid.hpp>

#include <opm/common/utility/parameters/ParameterGroup.hpp>

#include <opm/upscaling/initCPGrid.hpp>


void Opm::initCPGrid(Dune::CpGrid& grid, const Opm::ParameterGroup& param) {
    std::string fileformat = param.get<std::string>("fileformat");
    if (fileformat == "eclipse") {
        std::string filename = param.get<std::string>("filename");
        if (param.has("z_tolerance")) {
            std::cerr << "****** Warning: z_tolerance parameter is obsolete, use PINCH in deck input instead\n";
        }
        bool periodic_extension = param.getDefault<bool>("periodic_extension", false);
        bool turn_normals = param.getDefault<bool>("turn_normals", false);
        Opm::Parser parser;
        auto deck = parser.parseFile(filename);
        Opm::EclipseGrid inputGrid(deck);
        grid.processEclipseFormat(&inputGrid,
                                  /* ecl_state = */ nullptr,
                                  periodic_extension,
                                  turn_normals,
                                  /* clip_z = */ false,
                                  /* pinchActive = */ true,
                                  /* edge_conformal = */ false);
    } else if (fileformat == "cartesian") {
        std::array<int, 3> dims = {{ param.getDefault<int>("nx", 1),
                                     param.getDefault<int>("ny", 1),
                                     param.getDefault<int>("nz", 1) }};
        std::array<double, 3> cellsz = {{ param.getDefault<double>("dx", 1.0),
                                          param.getDefault<double>("dy", 1.0),
                                          param.getDefault<double>("dz", 1.0) }};
        grid.createCartesian(dims, cellsz);
    } else {
        OPM_THROW(std::runtime_error, "Unknown file format string: " + fileformat);
    }
}
