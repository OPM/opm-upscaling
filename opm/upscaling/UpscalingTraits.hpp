//===========================================================================
//
// File: UpscalingTraits.hpp
//
// Created: Wed Apr 28 10:36:42 2010
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2010 Statoil ASA.

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

#ifndef OPM_UPSCALINGTRAITS_HEADER
#define OPM_UPSCALINGTRAITS_HEADER

#include <opm/porsol/common/SimulatorTraits.hpp>

namespace Opm
{
    typedef SimulatorTraits<Isotropic, Explicit> UpscalingTraitsBasic;
    //typedef SimulatorTraits<Isotropic, Implicit> UpscalingTraitsBasicImplicit;
    typedef SimulatorTraits<Anisotropic, Explicit> UpscalingTraitsAnisoRelperm;

} // namespace Opm


#endif // OPM_UPSCALINGTRAITS_HEADER
