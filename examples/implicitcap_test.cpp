//===========================================================================
//
// File: simulator_test.cpp
//
// Created: Fri Aug  7 10:08:17 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Bård Skaflestad     <bard.skaflestad@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

  This file is part of The Open Reservoir Simulator Project (OpenRS).

  OpenRS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OpenRS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
*/


#define VERBOSE
//#define USE_TBB

#include "config.h"

#include "SimulatorTester.hpp"
#include "SimulatorTesterFlexibleBC.hpp"
#include <opm/porsol/common/SimulatorTraits.hpp>
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
#include <dune/common/parallel/mpihelper.hh>
#else
#include <dune/common/mpihelper.hh>
#endif

#ifdef USE_TBB
#include <tbb/task_scheduler_init.h>
#endif // USE_TBB

#include <iostream>

using namespace Opm;

typedef SimulatorTraits<Isotropic, ImplicitCap> SimTraits;
typedef SimulatorTesterFlexibleBC<SimTraits> Simulator;

int main(int argc, char** argv)
try
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
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}


