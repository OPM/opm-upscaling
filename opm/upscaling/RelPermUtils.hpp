/*
  Copyright 2010 Statoil ASA.

  This file is part of The Open Porous Media project (OPM).

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

/** @file RelPermutils.hpp
    @brief Helper utilities for relperm upscaling applications
 */

#ifndef OPM_UPSCALING_RELPERM_UTILS_HPP
#define OPM_UPSCALING_RELPERM_UTILS_HPP

#include <opm/upscaling/SinglePhaseUpscaler.hpp>

namespace Opm {
  //! \brief Get value from tensor
  //! \param[in] K The tensor to extract the value from
  //! \param[in] voigt_idx The voigt index for value to extract (0..8)
  //! \return The requested value
  double getVoigtValue(const SinglePhaseUpscaler::permtensor_t& K, int voigt_idx);

  //! \brief Set value in tensor
  //! \param[out] K The tensor to set value in
  //! \param[in] voigt_idx The voigt index for value to set (0..8)
  //! \param[in] val Value to set in tensor
  void setVoigtValue(SinglePhaseUpscaler::permtensor_t& K, int voigt_idx, double val);
}

#endif
