/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

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

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <opm/common/utility/platform_dependent/disable_warnings.h>

#include <dune/common/version.hh>

#include <dune/common/parallel/mpihelper.hh>

#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/parser/eclipse/Units/Units.hpp>

#include <opm/upscaling/SinglePhaseUpscaler.hpp>

#include <iostream>

using namespace Opm;
using namespace Opm::prefix;
using namespace Opm::unit;

int main(int argc, char** argv)
try
{
    Dune::MPIHelper::instance(argc, argv);

    SinglePhaseUpscaler upscaler;
    {
        auto param = Opm::ParameterGroup(argc, argv);

        upscaler.init(param);
    }

    auto upscaled_K = upscaler.upscaleSinglePhase();
    {
        const auto fact = convert::to(1.0*square(meter), milli*darcy);

        upscaled_K *= fact;
    }

    std::cout.precision(15);
    std::cout << "Upscaled K in millidarcy:\n" << upscaled_K << std::endl;
}
catch (const std::exception& e) {
    std::cerr << "Program threw an exception: " << e.what() << '\n';
    throw;
}
