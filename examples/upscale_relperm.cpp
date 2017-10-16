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

/**
   @file upscale_relperm.C
   @brief Upscales relative permeability as a fuction of water saturation assuming capillary equilibrium.

   Description:

   Reads in a lithofacies geometry in Eclipse format, reads in J(S_w)
   and relpermcurve(S_w) for each stone type, and calculates upscaled
   (three directions) relative permeability curves as a function of Sw.

   The relative permeability computation is based on
     - Capillary equilibrium, p_c is spatially invariant.
     - Optional gravitational effects. If gravity is not specified,
       gravity will be assumed to be zero.
   Units handling:
     - Assumes cornerpoint file reports lengths in cm.
     - Input surface tension is in dynes/cm
     - Input density is in g/cm^3
     - The denominator \sigma * cos(\phi) in J-function scaling
       is what we call "surface tension". If angle dependency is to be
       included, calculate the "surface tension" yourself.
     - Outputted capillary pressure is in Pascals.

   Steps in the code:

   1: Process command line options.
   2: Read Eclipse file
   3: Read relperm- and J-function for each stone-type.
   4: Tesselate the grid (Sintef code)
   5: Find minimum and maximum capillary pressure from the
      J-functions in each cell.
   6: Upscale water saturation as a function of capillary pressure
   7: Upscale single phase permeability.
   8: Upscale phase permeability for capillary pressures
      that corresponds to a uniform saturation grid, and
      compute relative permeability.
   9: Print output to screen and optionally to file.

 */
#include <config.h>

#include <opm/common/utility/platform_dependent/disable_warnings.h>

#include <dune/common/version.hh>

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
#include <dune/common/parallel/mpihelper.hh>
#else
#include <dune/common/mpihelper.hh>
#endif

#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/parser/eclipse/Units/Units.hpp>

#include <opm/upscaling/RelPermUtils.hpp>
#include <opm/upscaling/SinglePhaseUpscaler.hpp>

#include <algorithm>
#include <cfloat> // for DBL_MAX/DBL_MIN
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <sys/utsname.h>

using namespace Opm;
using namespace std;


static void usage()
{
    cerr << "Usage: upscale_relperm <options> <eclipsefile> stoneA.txt stoneB.txt ..." << endl <<
        "where the options are:" << endl <<
        "  -bc <string>                 -- which boundary conditions to use." << endl <<
        "                                  Possible values are f (fixed), l (linear)" << endl <<
        "                                  and p (periodic). Default f (fixed)." << endl <<
        "  -points <integer>            -- Number of saturation points to upscale for." << endl <<
        "                                  Uniformly distributed within saturation endpoints." << endl <<
        "                                  Default 30." << endl <<
        "  -relPermCurve <integer>      -- For isotropic input, the column number in the stone-files" << endl <<
        "                                  that represents the phase to be upscaled," << endl <<
        "                                  typically 2 (default) for water and 3 for oil." << endl <<
        "  -jFunctionCurve <integer>    -- the column number in the stone-files that" << endl <<
        "                                  represent the Leverett J-function. Default 4." << endl <<
        "  -upscaleBothPhases <bool>    -- If this is true, relperm for both phases will be upscaled" << endl <<
        "                                  and both will be outputted to Eclipse format. Default true." << endl <<
        "                                  For isotropic input, relPermCurves is assumed to be 2 and 3," << endl <<
        "                                  for anisotropic input, relPermCurves are assumed to be 3-5" << endl <<
        "                                  and 6-8 respectively for the two phases" << endl <<
        "  -gravity <float>             -- use 9.81 for standard gravity. Default zero. Unit m/s^2." << endl <<
        "  -surfaceTension <float>      -- Surface tension to use in J-function/Pc conversion." << endl <<
        "                                  Default 11 dynes/cm (oil-water systems). In absence of" << endl <<
        "                                  a correct value, the surface tension for gas-oil systems " << endl <<
        "                                  could be 22.5 dynes/cm." << endl <<
        "  -waterDensity <float>        -- density of water, only applicable to non-zero" << endl <<
        "                                  gravity, g/cm³. Default 1" << endl <<
        "  -oilDensity <float>          -- density of oil, only applicable to non-zero" << endl <<
        "                                  gravity, g/cm³. Default 0.6" << endl <<
        "  -output <string>             -- filename for where to write upscaled values." << endl <<
        "                                  If not supplied, output will only go to " << endl <<
        "                                  the terminal (standard out)." << endl <<
        "  -interpolate <integer>       -- If supplied, the output data points will be" << endl <<
        "                                  interpolated using monotone cubic interpolation" << endl <<
        "                                  on a uniform grid with the specified number of" << endl <<
        "                                  points. Suggested value: 1000." << endl <<
        "  -maxPermContrast <float>     -- maximal permeability contrast in model." << endl <<
        "                                  Default 10^7" << endl <<
        "  -minPerm <float>             -- Minimum floating point value allowed for" << endl <<
        "                                  phase permeability in computations. If set to zero," << endl <<
        "                                  some models can end up singular. Default 10^-12" << endl <<
        "  -maxPerm <float>             -- Maximum floating point value allowed for" << endl <<
        "                                  permeability. " << endl <<
        "                                  Default 100000. Unit Millidarcy." << endl <<
        "  -fluids <string>             -- Either ow for oil/water systems or go for gas/oil systems. Default ow." << endl <<
        "                                  In case of go, the waterDensity option should be set to gas density" << endl <<
        "                                  Also remember to specify the correct surface tension" << endl <<
        "  -krowxswirr <float>          -- Oil relative permeability in x-direction at Swirr(from SWOF table)." << endl <<
        "                                  In case of oil/gas, this value is needed to ensure consistensy" << endl <<
        "                                  between SWOF and SGOF tables. Only has affect if fluids is set to go" << endl <<
        "                                  and upscaleBothPhases is true." << endl <<
        "                                  If not set, the point is not inserted into the final table." << endl <<
        "  -krowyswirr <float>          -- Oil relative permeability in y-direction at Swirr(from SWOF table). See krowxswirr." << endl <<
        "  -krowzswirr <float>          -- Oil relative permeability in z-direction at Swirr(from SWOF table). See krowxswirr." << endl <<
        "  -doEclipseCheck <bool>       -- Default true. Check that input relperm curves includes relperms at critical" << endl <<
        "                                  saturation points, i.e. that krw(swcrit)=0 and krow(swmax) = 0 and similar for oil/gas." << endl <<
        "  -critRelpermThresh <float>   -- If minimum relperm values are less than this threshold, they are set to zero" << endl <<
        "                                  and will pass the EclipseCheck. Default 10^-6" << endl <<
        "If only one stone-file is supplied, it is used for all stone-types defined" << endl <<
        "in the geometry. If more than one, it corresponds to the SATNUM-values." << endl;
    // "minPoro" intentionally left undocumented
    // "saturationThreshold"  also
}


static void usageandexit() {
    usage();
    exit(1);
}

//! \brief Parse command line arguments into string map.
//! \param[in,out] options The map of options. Should be filled with default values on entry.
//! \param[in] varnum Number of arguments
//! \param[in] vararg The arguments
//! \param[in] verbose Whether or not to print parsed command line arguments
//! \returns index of input file if positive, negated index to offending argument on failure.
static int parseCommandLine(std::map<std::string,std::string>& options,
                            int varnum, char** vararg, bool verbose)
{
    int argeclindex = 0;
    for (int argidx = 1; argidx < varnum; argidx += 2)  {
        if (string(vararg[argidx]).substr(0,1) == "-")    {
            string searchfor = string(vararg[argidx]).substr(1); // Chop off leading '-'
            /* Check if it is a match */
            if (options.count(searchfor) == 1) {
                options[searchfor] = string(vararg[argidx+1]);
                if (verbose)
                    cout << "Parsed command line option: "
                         << searchfor << " := " << vararg[argidx+1] << endl;
                argeclindex = argidx + 2;
            }
            else
                return -argidx;
        }
        else {
            // if vararg[argidx] does not start in '-',
            // assume we have found the position of the Eclipse-file.
            argeclindex = argidx;
            break; // out of for-loop,
        }
    }

    return argeclindex;
}

//! \brief Return eclipse-style output filename.
//! \param[in] opfname Base output file name.
//! \param[in] comp Component (X, Y, Z).
//! \param[in] sat Fluid system type.
//! \return Eclipse-style filename for requested component/fluid system combination.
static std::string getEclipseOutputFile(const std::string& opfname, char comp, char sat)
{
    string fnbase = opfname.substr(0,opfname.find_last_of('.'));
    return fnbase + "-" +comp + ".S" + sat + "OF";
}


//! \brief Write eclipse-style output file.
//! \param[in] RelPermValues RelPerm values to write.
//! \param[in] Satvalues Saturation values to write.
//! \param[in] Pvalues Pressure values to write.
//! \param[in] options Option structure.
//! \param[in] component Component to write (0..2).
//! \param[in] owsystem Fluid system type.
  template<class Lazy>
static void writeEclipseOutput(Lazy& RelPermValues,
                               const std::vector<double>& Satvalues,
                               const std::vector<double>& Pvalues,
                               std::map<std::string,std::string>& options,
                               int component, bool owsystem)
{
    std::stringstream swof;
    char sat                  = (owsystem?'W':'G');
    char comp                 = 'x'+component;
    std::string krowstring = std::string("krow") + comp + "swirr";
    double krowswirr = atof(options[krowstring].c_str());
    const int outputprecision = atoi(options["outputprecision"].c_str());
    const int fieldwidth      = outputprecision + 8;

    // x-direction
    swof << "-- This file is based on the results in " << endl
         << "-- " << options["output"] << endl
         << "-- for relperm in " << comp << "-direction." << endl
         << "-- Pressure values (Pc) given in bars." << endl
         << "--        S" << (char)std::tolower(sat) << "       Kr"
                          << (char)std::tolower(sat) << comp << comp
                          << "      Kro" << (char)std::tolower(sat) << comp << comp
                          << "      Pc(bar)" << endl
         << "--S" << sat << "OF" << endl;
    if (krowswirr > 0) {
        swof << showpoint << setw(fieldwidth) << setprecision(outputprecision) << 0
             << showpoint << setw(fieldwidth) << setprecision(outputprecision) << 0
             << showpoint << setw(fieldwidth) << setprecision(outputprecision) << krowswirr
             << showpoint << setw(fieldwidth) << setprecision(outputprecision) << 0 << endl;
    }
    for (size_t i=0; i < Satvalues.size(); ++i) {
        swof << showpoint << setw(fieldwidth) << setprecision(outputprecision) << Satvalues[i]
             << showpoint << setw(fieldwidth) << setprecision(outputprecision) << RelPermValues[0][component][i]
             << showpoint << setw(fieldwidth) << setprecision(outputprecision) << RelPermValues[1][component][i]
             << showpoint << setw(fieldwidth) << setprecision(outputprecision) << Pvalues[i]/100000.0 << endl;
    }
    swof << "/" << endl;
    std::ofstream file;
    file.open(getEclipseOutputFile(options["output"], std::toupper(comp), sat),
                                   std::ios::out | std::ios::trunc);
    file << swof.str();
    file.close();
}

namespace {
    std::map<std::string, std::string> defineOptions()
    {
        return
        {
            {"bc",                            "f"}, // Fixed boundary conditions
            {"points",                       "30"}, // Number of saturation points (uniformly distributed within saturation endpoints)
            {"relPermCurve",                  "2"}, // Which column in the rock types are upscaled
            {"upscaleBothPhases",          "true"}, // Whether to upscale for both phases in the same run. Default true.
            {"jFunctionCurve",                "4"}, // Which column in the rock type file is the J-function curve
            {"surfaceTension",               "11"}, // Surface tension given in dynes/cm
            {"output",                         ""}, // If this is set, output goes to screen and to this file.
            {"gravity",                     "0.0"}, // default is no gravitational effects
            {"waterDensity",                "1.0"}, // default density of water, only applicable to gravity
            {"oilDensity",                  "0.6"}, // ditto
            {"interpolate",                   "0"}, // default is not to interpolate
            {"maxpoints",                  "1000"}, // maximal number of saturation points.
            {"outputprecision",               "4"}, // number of significant numbers to print
            {"maxPermContrast",             "1e7"}, // maximum allowed contrast in each single-phase computation
            {"minPerm",                   "1e-12"}, // absolute minimum for allowed cell permeability
            {"maxPerm",                  "100000"}, // maximal allowed cell permeability
            {"minPoro",                  "0.0001"}, // this limit is necessary for pcmin/max computation
            {"saturationThreshold",     "0.00001"}, // accuracy threshold for saturation, we ignore Pc values that
                                                    // give so small contributions near endpoints.
            {"linsolver_tolerance",       "1e-12"}, // residual tolerance for linear solver
            {"linsolver_verbosity",           "0"}, // verbosity level for linear solver
            {"linsolver_max_iterations",      "0"}, // Maximum number of iterations allow, specify 0 for default
            {"linsolver_type",                "3"}, // Type of linear solver: 0 = ILU0/CG, 1 = AMG/CG, 2 KAMG/CG, 3 FAST_AMG/CG
            {"linsolver_prolongate_factor", "1.0"}, // Prolongation factor in AMG
            {"linsolver_smooth_steps",        "1"}, // Number of smoothing steps in AMG
            {"fluids",                       "ow"}, // Whether upscaling for oil/water (ow) or gas/oil (go)
            {"krowxswirr",                   "-1"}, // Relative permeability in x direction of oil in corresponding oil/water system
            {"krowyswirr",                   "-1"}, // Relative permeability in y direction of oil in corresponding oil/water system
            {"krowzswirr",                   "-1"}, // Relative permeability in z direction of oil in corresponding oil/water system
            {"doEclipseCheck",             "true"}, // Check if minimum relpermvalues in input are zero (specify critical saturations)
            {"critRelpermThresh",          "1e-6"}, // Threshold for setting minimum relperm to 0 (thus specify critical saturations)
        };
    }

    void invertCapillaryPressureIsotropic(const char*                ROCKFILENAME,
                                          const int                  i,
                                          const int                  jFunctionCurve,
                                          const int                  relPermCurve,
                                          std::vector<std::string>&  JfunctionNames,
                                          Opm::RelPermUpscaleHelper& helper)
    {
        MonotCubicInterpolator Jtmp;

        try {
            Jtmp = MonotCubicInterpolator(ROCKFILENAME, 1, jFunctionCurve);
        }
        catch (const char* errormessage) {
            std::stringstream str;

            str << "Error: " << errormessage << endl
                << "Check filename and -jFunctionCurve" << endl;

            throw std::runtime_error(str.str());
        }

        // Invert J-function, now we get saturation as a function of pressure:
        if (Jtmp.isStrictlyMonotone()) {
            helper.InvJfunctions.emplace_back(Jtmp.get_fVector(),
                                              Jtmp.get_xVector());

            JfunctionNames.push_back(ROCKFILENAME);

            auto& krfunc = helper.Krfunctions[0];

            if (helper.upscaleBothPhases) {
                krfunc[0].emplace_back(ROCKFILENAME, 1, 2);
                krfunc[1].emplace_back(ROCKFILENAME, 1, 3);
            }
            else {
                krfunc[0].emplace_back(ROCKFILENAME, 1, relPermCurve);
            }
        }
        else {
            std::stringstream str;

            str << "Error: Jfunction " << (i + 1)
                << " in rock file "    << ROCKFILENAME
                << " was not invertible.";

            throw std::runtime_error(str.str());
        }
    }

    void invertCapillaryPressureAnIsotropic(const char*                ROCKFILENAME,
                                            const int                  i,
                                            const int                  /* jFunctionCurve */,
                                            const int                  /* relPermCurve */,
                                            std::vector<std::string>&  JfunctionNames,
                                            Opm::RelPermUpscaleHelper& helper)
    {
        MonotCubicInterpolator Pctmp;

        try {
            Pctmp = MonotCubicInterpolator(ROCKFILENAME, 2, 1);
        }
        catch (const char* errormessage) {
            std::stringstream str;

            str << "Error: " << errormessage << '\n'
                << "Check filename and columns 1 and 2 (Pc and "
                << helper.saturationstring << ')';

            throw std::domain_error(str.str());
        }

        // Invert Pc(Sw) curve into Sw(Pc):
        if (Pctmp.isStrictlyMonotone()) {
            helper.SwPcfunctions.emplace_back(Pctmp.get_fVector(),
                                              Pctmp.get_xVector());

            JfunctionNames.push_back(ROCKFILENAME);

            auto& krfunc = helper.Krfunctions;

            krfunc[0][0].emplace_back(ROCKFILENAME, 2, 3);
            krfunc[1][0].emplace_back(ROCKFILENAME, 2, 4);
            krfunc[2][0].emplace_back(ROCKFILENAME, 2, 5);

            if (helper.upscaleBothPhases) {
                krfunc[0][1].emplace_back(ROCKFILENAME, 2, 6);
                krfunc[1][1].emplace_back(ROCKFILENAME, 2, 7);
                krfunc[2][1].emplace_back(ROCKFILENAME, 2, 8);
            }
        }
        else {
            std::stringstream str;

            str << "Error: Pc(" << helper.saturationstring << ") curve "
                << (i + 1) << " in rock file "
                << ROCKFILENAME << " was not invertible.";

            throw std::runtime_error(str.str());
        }
    }

    void extractSaturationFunctions(const int                  varnum,
                                    char*                      vararg[],
                                    const int                  rockfileindex,
                                    const int                  stone_types,
                                    const int                  jFunctionCurve,
                                    const int                  relPermCurve,
                                    std::vector<std::string>&  JfunctionNames,
                                    Opm::RelPermUpscaleHelper& helper)
    {
        // If input is anisotropic, then we are in second mode with
        // different input file format

        auto invertPc = helper.anisotropic_input
            ? invertCapillaryPressureAnIsotropic
            : invertCapillaryPressureIsotropic;

        for (int i = 0 ; i < stone_types; ++i) {
            const char* ROCKFILENAME =
                vararg[(rockfileindex + stone_types == varnum)
                       ? rockfileindex + i
                       : rockfileindex];

            // Check if rock file exists and is readable:
            {
                std::ifstream rockfile(ROCKFILENAME, std::ios::in);

                if (!rockfile) {
                    std::stringstream str;

                    str << "Error: Filename " << ROCKFILENAME
                        << " not found or not readable.";

                    throw std::runtime_error(str.str());
                }
            }

            invertPc(ROCKFILENAME, i, jFunctionCurve,
                     relPermCurve, JfunctionNames, helper);
        }
    }
}

int main(int varnum, char** vararg)
try
{
   // Variables used for timing/profiling:
   clock_t start, finish;
   double timeused = 0.0, timeused_tesselation = 0.0;
   double timeused_upscale_wallclock = 0.0;

   /******************************************************************************
    * Step 1:
    * Process command line options
    */

   Dune::MPIHelper& mpi=Dune::MPIHelper::instance(varnum, vararg);
   const int mpi_rank = mpi.rank();
#ifdef HAVE_MPI
   const int mpi_nodecount = mpi.size();
#endif

   if (varnum == 1) { /* If no arguments supplied ("upscale_relperm" is the first "argument") */
      usage();
      exit(1);
   }

   /*
     Populate options-map with default values
   */
   auto options = defineOptions();

   /* Check first if there is anything on the command line to look for */
   if (varnum == 1) {
      if (mpi_rank == 0)
        cout << "Error: No Eclipsefile or stonefiles found on command line." << endl;
      usageandexit();
   }

   /*
      'argeclindex' is so that vararg[argeclindex] = the eclipse filename.
   */
   int argeclindex = parseCommandLine(options, varnum, vararg, mpi_rank == 0);
   if (argeclindex < 0) {
       if (mpi_rank == 0)
           cout << "Option -" << vararg[-argeclindex] << " unrecognized." << endl;
        usageandexit();
   }

   RelPermUpscaleHelper helper(mpi_rank, options);
   const auto owsystem = helper.saturationstring == "Sw";

   // argeclindex should now point to the eclipse file
   static char* ECLIPSEFILENAME(vararg[argeclindex]);
   argeclindex += 1; // argeclindex jumps to next input argument, now it points to the stone files.

   // argeclindex now points to the first J-function. This index is not
   // to be touched now.
   static int rockfileindex = argeclindex;


   /* Check if at least one J-function is supplied on command line */
   if (varnum <= rockfileindex)
        throw std::runtime_error("Error: No J-functions found on command line.");

   /* Check validity of boundary conditions chosen, and make booleans
      for boundary conditions, this allows more readable code later. */
   helper.setupBoundaryConditions();

   bool isFixed    = helper.boundaryCondition == SinglePhaseUpscaler::Fixed,
        isLinear   = helper.boundaryCondition == SinglePhaseUpscaler::Linear,
        isPeriodic = helper.boundaryCondition == SinglePhaseUpscaler::Periodic;

   // If this number is 1 or higher, the output will be interpolated, if not
   // the computed data is untouched.
   const int interpolationPoints = atoi(options["interpolate"].c_str());

   const auto doInterpolate = interpolationPoints > 1;

   /***********************************************************************
    * Step 2:
    * Load geometry and data from Eclipse file
    */


   // Read data from the Eclipse file and
   // populate our vectors with data from the file

   // Test if filename exists and is readable
   start = clock();

   auto deck = RelPermUpscaleHelper::parseEclipseFile(ECLIPSEFILENAME);

   finish   = clock();
   timeused = (double(finish) - double(start)) / CLOCKS_PER_SEC;
   if (helper.isMaster) {
       cout << " (" << timeused << " secs)" << endl;
   }

   {
       const double minPerm = atof(options["minPerm"].c_str());
       const double maxPerm = atof(options["maxPerm"].c_str());
       const double minPoro = atof(options["minPoro"].c_str());

       helper.sanityCheckInput(deck, minPerm, maxPerm, minPoro);
   }

   /***************************************************************************
    * Step 3:
    * Load relperm- and J-function-curves for the stone types.
    * We read columns from text-files, syntax allowed is determined
    * by MonotCubicInterpolator which actually opens and parses the
    * text files.
    *
    * If a standard eclipse data file is given as input, the data columns
    * should be:
    *    Sw    Krw   Kro   J-func
    * (In this case, the option -relPermCurve determines which of Krw or Kro is used)
    *
    * If output from this very program is given as input, then the data columns read
    *   Pc    Sw    Krx   Kry   Krz
    *
    * (and the option -relPermCurve and -jFunctionCurve are ignored)
    *
    * How do we determine which mode of operation?
    *  - If PERMY and PERMZ are present in grdecl-file, we are in the anisotropic mode
    *
    */

   // Number of stone-types is max(satnums):

   // If there is only one J-function supplied on the command line,
   // use that for all stone types.

   const auto stone_types =
       static_cast<int>(*(max_element(helper.satnums.begin(),
                                      helper.satnums.end())));

   // This decides whether we are upscaling both phases in this run or only one
   helper.upscaleBothPhases = (options["upscaleBothPhases"] == "true");
   helper.points            = atoi(options["points"].c_str());

   // Input for surfaceTension is dynes/cm
   // SI units are Joules/square metre => multiply with 10^-3 to obtain SI units.
   const double surfaceTension = atof(options["surfaceTension"].c_str()) * 1e-3;

   const double gravity      = atof(options["gravity"].c_str());
   const bool includeGravity = (std::fabs(gravity) > DBL_MIN); // true for non-zero gravity

   const int outputprecision = atoi(options["outputprecision"].c_str());

   // Handle two command line input formats, either one J-function for all stone types
   // or one each. If there is only one stone type, both code blocks below are equivalent.

   if (! ((varnum == rockfileindex + stone_types) ||
          (varnum == rockfileindex + 1          )))
   {
       throw std::runtime_error("Error:  Wrong number of stone-functions provided.");
   }

   // Placeholder for the names of the loaded J-functions.
   std::vector<string> JfunctionNames;
   {
       // This decides whether we are upscaling water or oil relative permeability
       const int relPermCurve   = atoi(options["relPermCurve"].c_str());
       const int jFunctionCurve = atoi(options["jFunctionCurve"].c_str());

       extractSaturationFunctions(varnum, vararg,
                                  rockfileindex, stone_types,
                                  jFunctionCurve, relPermCurve,
                                  JfunctionNames, helper);
   }

   // Check if input relperm curves satisfy ECLIPSE requirement of
   // specifying critical saturations
   if (helper.doEclipseCheck) {
       helper.checkCriticalSaturations();
   }

   /*****************************************************************************
    * Step 4:
    * Generate tesselated grid:
    * This is a step needed for the later discretization code to figure out which
    * cells are connected to which. Each cornerpoint-cell is tesselated into 8 tetrahedrons.
    *
    * In case of non-zero gravity, calculate z-values of every cell:
    *   1) Compute height of model by averaging z-values of the top layer corners.
    *   2) Calculate density difference between phases in SI-units
    *   3) Go through each cell and find the z-values of the eight corners of the cell.
    *      Set height of cell equal to average of z-values of the corners minus half of
    *      model height. Now the cell height is relative to model centre.
    *      Set pressure difference for the cell equal to density difference times gravity
    *      constant times cell height times factor 10^-7 to obtain bars (same as p_c)
    */

   timeused_tesselation = helper.tesselateGrid(deck);

   /* If gravity is to be included, calculate z-values of every cell: */
   if (includeGravity) {
       helper.calculateCellPressureGradients();
   }

   /******************************************************************************
    * Step 5:
    * Go through each cell and calculate the minimum and
    * maximum capillary pressure possible in the cell, given poro,
    * perm and the J-function for the cell.  This depends on the
    * J-function in that they represent all possible saturations,
    * ie. we do not want to extrapolate the J-functions (but we might
    * have to do that later in the computations).
    *
    * The user-supplied surface tension is ignored until
    * the final output of results.
    */

   helper.calculateMinMaxCapillaryPressure();

   /***************************************************************************
    * Step 6:
    * Upscale capillary pressure function.
    *
    * This is upscaled in advance in order to be able to have uniformly distributed
    * saturation points for which upscaling is performed.
    *
    * Capillary pressure points are chosen heuristically in order to
    * ensure largest saturation interval between two saturation points
    * is 1/500 of the saturation interval. Monotone cubic interpolation
    * will be used afterwards for accessing the tabulated values.
    */

   helper.upscaleCapillaryPressure();

   /*****************************************************************************
    * Step 7:
    * Upscale single phase permeability
    * This uses the PERMX in the eclipse file as data, and upscales using
    * fixed boundary (no-flow) conditions
    *
    * In an MPI-environment, this is only done on the master node.
    */

   helper.upscaleSinglePhasePermeability();

   /*****************************************************************
    * Step 8:
    *
    * Loop through a given number of uniformly distributed saturation points
    * and upscale relative permeability for each of them.
    *    a: Make vector of capillary pressure points corresponding to uniformly
    *       distributed water saturation points between saturation endpoints.
    *    b: Loop over capillary pressure points
    *       1)  Loop over all cells to find the saturation value given the
    *           capillary pressure found in (a). Given the saturation value, find the
    *           phase permeability in the cell given input relperm curve and input
    *           permeability values.
    *       2)  Upscale phase permeability for the geometry.
    *    c: Calculate relperm tensors from all the phase perm tensors.
    */

    double avg_upscaling_time_pr_point;
    std::tie(timeused_upscale_wallclock,
             avg_upscaling_time_pr_point) =
        helper.upscalePermeability(mpi_rank);

   /*
    * Step 8c: Make relperm values from phaseperms
    *          (only master node can do this)
    */
   std::array<vector<vector<double>>,2> RelPermValues;
   RelPermValues[0] = helper.getRelPerm(0);
   if (helper.upscaleBothPhases) {
       RelPermValues[1] = helper.getRelPerm(1);
   }

   /*********************************************************************************
    *  Step 9
    *
    * Output results to stdout and optionally to file. Note, we only output to
    * file if the '-outputWater'-option and/or '-outputOil' has been set, as this option is an
    * empty string by default.
    */
   if (helper.isMaster) {
       stringstream outputtmp;

       // Print a table of all computed values:
       outputtmp << "######################################################################\n"
                 << "# Results from upscaling relative permeability.\n"
                 << "#\n";

#if defined(HAVE_MPI) && HAVE_MPI
       outputtmp << "#          (MPI-version)\n";
#endif  // HAVE_MPI

       {
           time_t now = std::time(NULL);
           outputtmp << "# Finished: " << asctime(localtime(&now));
       }

       {
           utsname hostname;   uname(&hostname);
           outputtmp << "# Hostname: " << hostname.nodename << endl;
       }

       outputtmp << "#" << '\n'
                 << "# ECLIPSE file: " << ECLIPSEFILENAME << '\n'
                 << "#        cells: " << helper.tesselatedCells << '\n'
                 << "#  Pore volume: " << helper.poreVolume << '\n'
                 << "#       volume: " << helper.volume << '\n'
                 << "#     Porosity: " << (helper.poreVolume / helper.volume) << '\n'
                 << "#\n";

       if (! helper.anisotropic_input) {
           const auto& invJ = helper.InvJfunctions;

           for (int i = 0; i < stone_types ; ++i) {
               outputtmp << "# Stone " << (i + 1) << ": "
                         << JfunctionNames[i]
                         << " (" << invJ[i].getSize() << " points)\n";
           }

           outputtmp << "#         jFunctionCurve: "
                     << options["jFunctionCurve"] << '\n';

           if (! helper.upscaleBothPhases) {
               outputtmp << "#           relPermCurve: "
                         << options["relPermCurve"] << '\n';
           }
       }
       else {
           // anisotropic input, not J-functions that are supplied on
           // command line (but vector JfunctionNames is still used)

           const auto& kr = helper.Krfunctions[0][0];

           for (int i = 0; i < stone_types ; ++i) {
               outputtmp << "# Stone " << (i + 1) << ": "
                         << JfunctionNames[i]
                         << " (" << kr[i].getSize() << " points)\n";
           }
       }

       outputtmp << "#" << '\n'
                 << "# Timings:   Tesselation: "
                 << timeused_tesselation << " secs" << '\n'
                 << "#              Upscaling: "
                 << timeused_upscale_wallclock << " secs";

#ifdef HAVE_MPI
       outputtmp << " (wallclock time)\n"
                 << "#                         "
                 << avg_upscaling_time_pr_point
                 << " secs pr. saturation point" << '\n'
                 << "#              MPI-nodes: "
                 << mpi_nodecount << '\n';

       // Single phase upscaling time is included here, in possibly a hairy
       // way.
       const auto speedup =
           (avg_upscaling_time_pr_point * (helper.points + 1) + timeused_tesselation)
           / (timeused_upscale_wallclock +
              avg_upscaling_time_pr_point +
              timeused_tesselation);

       outputtmp << "#                Speedup: " << speedup
                 << ", efficiency: " << (speedup / mpi_nodecount) << endl;

#else // !HAVE_MPI

       outputtmp << ", " << avg_upscaling_time_pr_point
                 << " secs avg for " << helper.points
                 << " runs" << endl;

#endif // HAVE_MPI

       outputtmp << "#\n"
                 << "# Options used:\n"
                 << "#     Boundary conditions: ";

       if (isFixed)    { outputtmp << "Fixed (no-flow)\n"; }
       if (isPeriodic) { outputtmp << "Periodic\n"; }
       if (isLinear)   { outputtmp << "Linear\n"; }

       outputtmp << "#                  points: " << options["points"] << '\n'
                 << "#         maxPermContrast: " << options["maxPermContrast"] << '\n'
                 << "#                 minPerm: " << options["minPerm"] << " [mD]\n"
                 << "#                 minPoro: " << options["minPoro"] << '\n'
                 << "#          surfaceTension: " << options["surfaceTension"] << " dynes/cm\n";

       if (includeGravity) {
           outputtmp << "#                 gravity: "
                     << options["gravity"] << " m/s²\n";

           if (owsystem) {
               outputtmp << "#            waterDensity: "
                         << options["waterDensity"] << " g/cm³\n";
           }
           else {
               outputtmp << "#              gasDensity: "
                         << options["waterDensity"] << " g/cm³\n";
           }

           outputtmp << "#              oilDensity: "
                     << options["oilDensity"] << " g/cm³\n";
       }
       else {
           outputtmp << "#                 gravity: 0\n";
       }

       if (doInterpolate) {
           outputtmp << "#             interpolate: "
                     << options["interpolate"] << " points\n";
       }

       {
           auto mD = [](const double k)
           {
               return Opm::unit::convert::to(k, Opm::prefix::milli*Opm::unit::darcy);
           };

           outputtmp << "#\n"
                     << "# Single phase permeability (mD)\n"
                     << "#  |Kxx  Kxy  Kxz| = "
                     << mD(helper.permTensor(0,0)) << "  "
                     << mD(helper.permTensor(0,1)) << "  "
                     << mD(helper.permTensor(0,2)) << '\n';

           outputtmp << "#  |Kyx  Kyy  Kyz| = "
                     << mD(helper.permTensor(1,0)) << "  "
                     << mD(helper.permTensor(1,1)) << "  "
                     << mD(helper.permTensor(1,2)) << '\n';

           outputtmp << "#  |Kzx  Kzy  Kzz| = "
                     << mD(helper.permTensor(2,0)) << "  "
                     << mD(helper.permTensor(2,1)) << "  "
                     << mD(helper.permTensor(2,2)) << '\n';

           outputtmp << "#\n";
       }

       if (doInterpolate) {
           outputtmp << "# NB: Data points shown are interpolated.\n";
       }

       outputtmp << "######################################################################\n";

       if (helper.upscaleBothPhases) {
           string phase1, phase2;

           if (owsystem) {
               phase1 = "w";
           }
           else {
               phase1 = "g";
           }

           phase2 = "o";

           if (isFixed) {
               outputtmp << "#  Pc (Pa)        " << helper.saturationstring
                         << "           Kr" << phase1
                         << "xx       Kr" << phase1
                         << "yy       Kr" << phase1
                         << "zz       Kr" << phase2
                         << "xx       Kr" << phase2
                         << "yy       Kr" << phase2
                         << "zz\n";
           }
           else if (isPeriodic || isLinear) {
               outputtmp << "#  Pc (Pa)        " << helper.saturationstring
                         << "           Kr" << phase1
                         << "xx       Kr" << phase1
                         << "yy       Kr" << phase1
                         << "zz       Kr" << phase1
                         << "yz       Kr" << phase1
                         << "xz       Kr" << phase1
                         << "xy       Kr" << phase1
                         << "zy       Kr" << phase1
                         << "zx       Kr" << phase1
                         << "yx       Kr" << phase2
                         << "xx       Kr" << phase2
                         << "yy       Kr" << phase2
                         << "zz       Kr" << phase2
                         << "yz       Kr" << phase2
                         << "xz       Kr" << phase2
                         << "xy       Kr" << phase2
                         << "zy       Kr" << phase2
                         << "zx       Kr" << phase2
                         << "yx\n";
           }
       }
       else {
           if (isFixed) {
               outputtmp << "#  Pc (Pa)        "
                         << helper.saturationstring
                         << "            Krxx        Kryy        Krzz\n";
           }
           else if (isPeriodic || isLinear) {
               outputtmp << "#  Pc (Pa)        "
                         << helper.saturationstring
                         << "            Krxx        Kryy        Krzz        "
                         << "Kryz        Krxz        Krxy        Krzy        "
                         << "Krzx        Kryx\n";
           }
       }

       auto Pvalues = helper.pressurePoints;

       // Multiply all pressures with the surface tension (potentially)
       // supplied at the command line. This multiplication has been
       // postponed to here to avoid division by zero and to avoid special
       // handling of negative capillary pressure in the code above.
       for (auto& press : Pvalues) {
           press *= surfaceTension;
       }

       auto Satvalues = helper.WaterSaturation; //.get_fVector();

       // If user wants interpolated output, do monotone cubic interpolation
       // by modifying the data vectors that are to be printed.
       if (doInterpolate) {

           // Split interval between minumum and maximum saturation into
           // (interpolationPoints - 1) equally sized smaller intervals.
           auto SatvaluesInterp = std::vector<double>{};
           {
               const auto m = std::minmax_element(Satvalues.begin(), Satvalues.end());

               const auto xmin = *m.first;
               const auto xmax = *m.second;

               // Make uniform grid in saturation axis
               SatvaluesInterp.reserve(interpolationPoints);

               auto interp = [xmin, xmax, interpolationPoints](const int i) -> double {
                   const auto t = static_cast<double>(i) / (interpolationPoints - 1);

                   return xmin + t*(xmax - xmin);
               };

               for (int i = 0; i < interpolationPoints; ++i) {
                   SatvaluesInterp.push_back(interp(i));
               }
           }

           // Now capillary pressure and computed relperm-values must be
           // viewed as functions of saturation, and then interpolated on
           // the uniform saturation grid.

           // Now overwrite existing Pvalues and relperm-data with
           // interpolated data:
           {
               const auto PvaluesVsSaturation =
                   MonotCubicInterpolator(Satvalues, Pvalues);

               Pvalues.clear();  Pvalues.reserve(SatvaluesInterp.size());

               for (const auto& s : SatvaluesInterp) {
                   Pvalues.push_back(PvaluesVsSaturation.evaluate(s));
               }
           }

           for (int voigtIdx = 0; voigtIdx < helper.tensorElementCount; ++voigtIdx) {
               auto& kr = RelPermValues[0][voigtIdx];

               const auto RelPermVsSaturation =
                   MonotCubicInterpolator(Satvalues, kr);

               kr.clear();  kr.reserve(SatvaluesInterp.size());

               for (const auto& s : SatvaluesInterp) {
                   kr.push_back(RelPermVsSaturation.evaluate(s));
               }
           }

           if (helper.upscaleBothPhases) {
               for (int voigtIdx = 0; voigtIdx < helper.tensorElementCount; ++voigtIdx) {
                   auto& kr = RelPermValues[1][voigtIdx];

                   const auto RelPermVsSaturation =
                       MonotCubicInterpolator(Satvalues, kr);

                   kr.clear();  kr.reserve(SatvaluesInterp.size());

                   for (const auto& s : SatvaluesInterp) {
                       kr.push_back(RelPermVsSaturation.evaluate(s));
                   }
               }
           }

           // Now also overwrite Satvalues
           Satvalues = std::move(SatvaluesInterp);
       }

       // The code below does not care whether the data is interpolated or not.
       const auto fieldwidth = outputprecision + 8;
       for (decltype(Satvalues.size())
                i = 0, n = Satvalues.size(); i < n; ++i)
       {
           outputtmp << showpoint << setw(fieldwidth)
                     << setprecision(outputprecision)
                     << Pvalues[i];

           outputtmp << showpoint << setw(fieldwidth)
                     << setprecision(outputprecision)
                     << Satvalues[i];

           for (int voigtIdx = 0; voigtIdx < helper.tensorElementCount; ++voigtIdx) {
               outputtmp << showpoint << setw(fieldwidth)
                         << setprecision(outputprecision)
                         << RelPermValues[0][voigtIdx][i];
           }

           if (helper.upscaleBothPhases) {
               for (int voigtIdx = 0; voigtIdx < helper.tensorElementCount; ++voigtIdx) {
                   outputtmp << showpoint << setw(fieldwidth)
                             << setprecision(outputprecision)
                             << RelPermValues[1][voigtIdx][i];
               }
           }

           outputtmp << endl;
       }

       cout << outputtmp.str();

       if (options["output"] != "") {
           cout << "Writing results to " << options["output"] << endl;

           ofstream outfile;

           outfile.open(options["output"].c_str(), ios::out | ios::trunc);
           outfile << outputtmp.str();
           outfile.close();
       }

       // If both phases are upscaled and output is specified, create SWOF
       // or SGOF files for ECLIPSE
       if (options["output"] != "" && helper.upscaleBothPhases) {
            cout << "Writing ECLIPSE compatible files to "
                 << getEclipseOutputFile(options["output"], 'X', owsystem ? 'W' : 'G')
                 << ", "
                 << getEclipseOutputFile(options["output"], 'Y', owsystem ? 'W' : 'G')
                 << " and "
                 << getEclipseOutputFile(options["output"], 'Z', owsystem ? 'W' : 'G')
                 << endl;

            for (int comp = 0; comp < 3; ++comp) {
                writeEclipseOutput(RelPermValues, Satvalues, Pvalues,
                                   options, comp, owsystem);
            }
       }
   }
}
catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    usageandexit();
}
