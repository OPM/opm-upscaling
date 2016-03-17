#include <opm/upscaling/ParserAdditions.hpp>

#include <opm/json/JsonObject.hpp>

#include <opm/parser/eclipse/Parser/ParserKeyword.hpp>

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
