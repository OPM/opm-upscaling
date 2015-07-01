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

#include <config.h>
#include <opm/upscaling/RelPermUtils.hpp>
#include <iostream>
#include <exception>

namespace Opm {

static const std::vector<size_t> voigt_idx_tab = {0,4,8,5,2,1,7,6,3}; //!< Voigt-to-C index table

// Assumes that permtensor_t use C ordering.
double getVoigtValue(const SinglePhaseUpscaler::permtensor_t& K, int voigt_idx)
{
#if !defined(NDEBUG)
    OPM_ERROR_IF(not ((K.numRows() == 3) && (K.numCols() == 3)),
                 "Function getVoigtValue() is only supported "
                 "for 3-by-3 tensors");
#endif

    if (voigt_idx < 0 || voigt_idx > 8) {
        std::cerr << "Voigt index out of bounds (only 0-8 allowed)" << std::endl;
        throw std::exception();
    }

    return K.data()[voigt_idx_tab[voigt_idx]];
}

// Assumes that permtensor_t use C ordering.
void setVoigtValue(SinglePhaseUpscaler::permtensor_t& K, int voigt_idx, double val)
{
#if !defined(NDEBUG)
    OPM_ERROR_IF(not ((K.numRows() == 3) && (K.numCols() == 3)),
                 "Function setVoigtValue() is only supported "
                 "for 3-by-3 tensors.");
#endif

    if (voigt_idx < 0 || voigt_idx > 8) {
        std::cerr << "Voigt index out of bounds (only 0-8 allowed)" << std::endl;
        throw std::exception();
    }

    K.data()[voigt_idx_tab[voigt_idx]] = val;
}

}
