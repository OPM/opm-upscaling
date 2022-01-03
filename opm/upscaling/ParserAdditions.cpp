/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.

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

#include <opm/upscaling/ParserAdditions.hpp>

#include <opm/json/JsonObject.hpp>

#include <opm/input/eclipse/Parser/ParserKeyword.hpp>

namespace {
    void add_RHO(Opm::Parser& parser)
    {
        const auto rhoJsonData =
R"~~({
        "name"     : "RHO",
        "sections" : [],
        "data"     : {"value_type" : "DOUBLE" }
}
)~~";

        parser.addParserKeyword(Json::JsonObject(rhoJsonData));
    }
} // Anonymous

void Opm::addNonStandardUpscalingKeywords(Parser& parser)
{
    add_RHO(parser);
}
