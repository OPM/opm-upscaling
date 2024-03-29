/*
  Copyright 2014 by Andreas Lauser

  This file is part of The Open Porous Media Project (OPM).

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

#ifndef OPM_UPSCALING_PARSER_ADDITIONS_HPP
#define OPM_UPSCALING_PARSER_ADDITIONS_HPP

#include <opm/input/eclipse/Parser/Parser.hpp>

namespace Opm {

    /*!
     * \brief This function registers non-standard keywords used by
     *        opm-upscaling in a parser object.
     *
     * The name of this function is intentionally long and awkward. This is
     * to discourage its use unless it is *really* necessary!
     */
    void addNonStandardUpscalingKeywords(Parser& parser);
}

#endif
