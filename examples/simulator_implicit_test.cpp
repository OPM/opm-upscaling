/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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



#define VERBOSE

#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/core/utility/have_boost_redef.hpp>

#include "SimulatorTester.hpp"
#include "SimulatorTesterFlexibleBC.hpp"
#include <dune/porsol/euler/EulerUpstreamImplicit.hpp>
#include <dune/porsol/common/SimulatorTraits.hpp>

namespace Dune
{
    template <class IsotropyPolicy>
    struct Implicit
    {
        template <class GridInterface, class BoundaryConditions>
        struct TransportSolver
        {
            //enum { Dimension = GridInterface::Dimension };
        	enum { Dimension = GridInterface::Dimension };
            typedef typename IsotropyPolicy::template ResProp<Dimension>::Type RP;

            typedef EulerUpstreamImplicit<GridInterface,
                                  RP,
                                  BoundaryConditions> Type;

        };
    };
}

using namespace Dune;

typedef SimulatorTraits<Isotropic, Implicit> SimTraits;
typedef SimulatorTesterFlexibleBC<SimTraits> Simulator;

int main(int argc, char** argv)
{
    Opm::parameter::ParameterGroup param(argc, argv);
    MPIHelper::instance(argc,argv);
    Simulator sim;
    sim.init(param);
    sim.run();
}

