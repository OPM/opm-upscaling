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


#define VERBOSE
//#define USE_TBB

#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/core/utility/have_boost_redef.hpp>

#include "SimulatorTester.hpp"
#include "SimulatorTesterFlexibleBC.hpp"
#include <opm/porsol/common/SimulatorTraits.hpp>
#include <dune/common/mpihelper.hh>

#ifdef USE_TBB
#include <tbb/task_scheduler_init.h>
#endif // USE_TBB

using namespace Opm;

typedef SimulatorTraits<Isotropic, Explicit> SimTraits;
typedef SimulatorTesterFlexibleBC<SimTraits> Simulator;

int main(int argc, char** argv)
{
    Opm::parameter::ParameterGroup param(argc, argv);
    Dune::MPIHelper::instance(argc,argv);
#ifdef USE_TBB
    int num_threads = param.getDefault("num_threads", tbb::task_scheduler_init::default_num_threads());
    tbb::task_scheduler_init init(num_threads);
#endif
    Simulator sim;
    sim.init(param);
    sim.run();
}

