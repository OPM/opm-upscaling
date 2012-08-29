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



#include "config.h"

#include <opm/core/utility/have_boost_redef.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <dune/porsol/blackoil/co2fluid/BlackoilCo2PVT.hpp>

int main(int argc, char** argv)
{
    double temperature = 300.;
    if (argc == 2) temperature = std::atof(argv[1]);

    Opm::BlackoilCo2PVT boPvt;
    Opm::EclipseGridParser ep;
    Opm::time::StopWatch clock;
  clock.start();
    boPvt.init(ep);

    boPvt.generateBlackOilTables(temperature);
  clock.stop();
  std::cout << "\n\nInitialisation and table generation - clock time (secs): " << clock.secsSinceStart() << std::endl;
}

