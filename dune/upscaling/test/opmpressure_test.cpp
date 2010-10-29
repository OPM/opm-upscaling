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



#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/upscaling/UpscalerBase.hpp>
#include <dune/common/Units.hpp>
#include <dune/porsol/common/SimulatorTraits.hpp>
#include <dune/porsol/mimetic/IfshInterface.hpp>

using namespace Dune;
using namespace Dune::prefix;
using namespace Dune::unit;


/// Combines the component traits into a single, parametrized type.
template <class RelpermPolicy, template <class> class TransportPolicy>
struct MyTraits : public RelpermPolicy, TransportPolicy<RelpermPolicy>
{
    /// The pressure/flow solver type.
    template <class GridInterface, class BoundaryConditions>
    struct FlowSolver
    {
        typedef IfshInterface<GridInterface,
                              typename RelpermPolicy::template ResProp<GridInterface::Dimension>::Type,
                              BoundaryConditions> Type;
    };
};

typedef UpscalerBase<MyTraits<Isotropic, Explicit> > Upscaler;

int main(int argc, char** argv)
{
    parameter::ParameterGroup param(argc, argv);
    // MPIHelper::instance(argc,argv);
    Upscaler upscaler;
    upscaler.init(param);
    Upscaler::permtensor_t upscaled_K = upscaler.upscaleSinglePhase();
    upscaled_K *= (1.0/(milli*darcy));
    std::cout.precision(15);
    std::cout << "Upscaled K in millidarcy:\n" << upscaled_K << std::endl;
}
