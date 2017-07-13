//===========================================================================
//
// File: aniso_simulator_test.cpp
//
// Created: Wed Oct 28 15:10:02 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            B�rd Skaflestad     <bard.skaflestad@sintef.no>
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

#include <opm/common/utility/platform_dependent/disable_warnings.h>

#include <dune/common/version.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include "SimulatorTester.hpp"
#include "SimulatorTesterFlexibleBC.hpp"

#include <opm/porsol/common/SimulatorTraits.hpp>

#ifdef USE_TBB
#include <tbb/task_scheduler_init.h>
#endif // USE_TBB

#include <iostream>

using namespace Opm;

typedef SimulatorTraits<Anisotropic, Explicit> SimTraits;
typedef SimulatorTesterFlexibleBC<SimTraits> Simulator;

int main(int argc, char** argv)
try
{
    Opm::ParameterGroup param(argc, argv);
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
