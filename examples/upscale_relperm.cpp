
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <cfloat> // for DBL_MAX/DBL_MIN
#include <map>
#include <sys/utsname.h>

#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
#include <dune/common/parallel/mpihelper.hh>
#else
#include <dune/common/mpihelper.hh>
#endif

#include <opm/core/utility/Units.hpp>
#include <opm/upscaling/SinglePhaseUpscaler.hpp>
#include <opm/upscaling/ParserAdditions.hpp>
#include <opm/upscaling/RelPermUtils.hpp>

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
   RelPermUpscaleHelper helper;

   helper.isMaster = (mpi_rank == 0);
   if (varnum == 1) { /* If no arguments supplied ("upscale_relperm" is the first "argument") */
      usage();
      exit(1);
   }

   /*
     Populate options-map with default values
   */
   map<string,string> options =
       {{"bc",                            "f"}, // Fixed boundary conditions
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
        {"critRelpermThresh",          "1e-6"}};// Threshold for setting minimum relperm to 0 (thus specify critical saturations)

   // Conversion factor, multiply mD numbers with this to get m² numbers
   const double milliDarcyToSqMetre =
       Opm::unit::convert::to(1.0*Opm::prefix::milli*Opm::unit::darcy,
                              Opm::unit::square(Opm::unit::meter));
   // Reference: http://www.spe.org/spe-site/spe/spe/papers/authors/Metric_Standard.pdf

   /* Check first if there is anything on the command line to look for */
   if (varnum == 1) {
      if (helper.isMaster) cout << "Error: No Eclipsefile or stonefiles found on command line." << endl;
      usageandexit();
   }


   /*
      'argeclindex' is so that vararg[argeclindex] = the eclipse filename.
   */
   int argeclindex = parseCommandLine(options, varnum, vararg, helper.isMaster);
   if (argeclindex < 0) {
       if (helper.isMaster)
           cout << "Option -" << vararg[-argeclindex] << " unrecognized." << endl;
        usageandexit();
   }
     
   // What fluid system are we dealing with? (oil/water or gas/oil)
   bool owsystem;
   if (options["fluids"] == "ow" || options["fluids"] == "wo") {
       owsystem=true;
       helper.saturationstring = "Sw";
   }
   else if (options["fluids"] == "go" || options["fluids"] == "og") {
       owsystem=false;
       helper.saturationstring = "Sg";
   }
   else {
       if (helper.isMaster) cerr << "Fluidsystem " << options["fluids"] << " not valid (-fluids option). Should be ow or go" << endl << endl;
       usageandexit();
   }

   // argeclindex should now point to the eclipse file
   static char* ECLIPSEFILENAME(vararg[argeclindex]);
   argeclindex += 1; // argeclindex jumps to next input argument, now it points to the stone files.

   // Boolean set to true if input permeability in eclipse-file has diagonal anisotropy.
   // (full-tensor anisotropy will be ignored)
   helper.anisotropic_input = false;
   
   // argeclindex now points to the first J-function. This index is not
   // to be touched now.
   static int rockfileindex = argeclindex;
   

   /* Check if at least one J-function is supplied on command line */
   if (varnum <= rockfileindex) {
       if (helper.isMaster) cerr << "Error: No J-functions found on command line." << endl;
       usageandexit();
   }
    
   /* Check validity of boundary conditions chosen, and make booleans 
      for boundary conditions, this allows more readable code later. */
   helper.setupBoundaryConditions(options);

   bool isFixed    = helper.boundaryCondition == SinglePhaseUpscaler::Fixed,
        isLinear   = helper.boundaryCondition == SinglePhaseUpscaler::Linear,
        isPeriodic = helper.boundaryCondition == SinglePhaseUpscaler::Periodic;

   // If this number is 1 or higher, the output will be interpolated, if not
   // the computed data is untouched.
   const int interpolationPoints = atoi(options["interpolate"].c_str());
   bool doInterpolate = false;
   if (interpolationPoints > 1) {
       doInterpolate = true;
   }
   
   /***********************************************************************
    * Step 2:
    * Load geometry and data from Eclipse file
    */

   
   // Read data from the Eclipse file and 
   // populate our vectors with data from the file
  
   // Test if filename exists and is readable
   ifstream eclipsefile(ECLIPSEFILENAME, ios::in);
   if (eclipsefile.fail()) {
       if (helper.isMaster) cerr << "Error: Filename " << ECLIPSEFILENAME << " not found or not readable." << endl;
       usageandexit();
   }
   eclipsefile.close(); 

   if (helper.isMaster) cout << "Parsing Eclipse file <" << ECLIPSEFILENAME << "> ... ";
   flush(cout);   start = clock();
   Opm::ParseMode parseMode;
   Opm::ParserPtr parser(new Opm::Parser());
   Opm::addNonStandardUpscalingKeywords(parser);
   Opm::DeckConstPtr deck(parser->parseFile(ECLIPSEFILENAME , parseMode));
   finish = clock();   timeused = (double(finish)-double(start))/CLOCKS_PER_SEC;
   if (helper.isMaster) cout << " (" << timeused <<" secs)" << endl;

   Opm::DeckRecordConstPtr specgridRecord = deck->getKeyword("SPECGRID")->getRecord(0);
   std::array<int,3> res;
   res[0] = specgridRecord->getItem("NX")->getInt(0);
   res[1] = specgridRecord->getItem("NY")->getInt(0);
   res[2] = specgridRecord->getItem("NZ")->getInt(0);

   const double maxPermContrast = atof(options["maxPermContrast"].c_str());
   const double minPerm = atof(options["minPerm"].c_str());
   const double maxPerm = atof(options["maxPerm"].c_str());
   const double minPoro = atof(options["minPoro"].c_str());
   const double saturationThreshold = atof(options["saturationThreshold"].c_str());

   helper.sanityCheckInput(deck, minPerm, maxPerm, minPoro);

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

   int stone_types = int(*(max_element(helper.satnums.begin(), helper.satnums.end())));
   
   std::vector<string> JfunctionNames; // Placeholder for the names of the loaded J-functions.

   // This decides whether we are upscaling water or oil relative permeability
   const int relPermCurve = atoi(options["relPermCurve"].c_str());
   // This decides whether we are upscaling both phases in this run or only one
   helper.upscaleBothPhases = (options["upscaleBothPhases"] == "true");

   const int jFunctionCurve        = atoi(options["jFunctionCurve"].c_str());
   helper.points                   = atoi(options["points"].c_str());
   const double gravity            = atof(options["gravity"].c_str());

   // Input for surfaceTension is dynes/cm
   // SI units are Joules/square metre
   const double surfaceTension     = atof(options["surfaceTension"].c_str()) * 1e-3; // multiply with 10^-3 to obtain SI units 
   const double waterDensity       = atof(options["waterDensity"].c_str());
   const double oilDensity         = atof(options["oilDensity"].c_str());
   const bool includeGravity       = (fabs(gravity) > DBL_MIN); // true for non-zero gravity
   const int outputprecision       = atoi(options["outputprecision"].c_str());

   // Handle two command line input formats, either one J-function for all stone types
   // or one each. If there is only one stone type, both code blocks below are equivalent.
   
   if (varnum == rockfileindex + stone_types) {
      for (int i=0 ; i < stone_types; ++i) {
         const char* ROCKFILENAME = vararg[rockfileindex+i];
         // Check if rock file exists and is readable:
         ifstream rockfile(ROCKFILENAME, ios::in);
         if (rockfile.fail()) {
            if (helper.isMaster) cerr << "Error: Filename " << ROCKFILENAME << " not found or not readable." << endl;
            usageandexit();
         }
         rockfile.close(); 
         
         if (! helper.anisotropic_input) {
             
             MonotCubicInterpolator Jtmp;
             try {
                 Jtmp = MonotCubicInterpolator(ROCKFILENAME, 1, jFunctionCurve); 
             }
             catch (const char * errormessage) {
                 if (helper.isMaster) cerr << "Error: " << errormessage << endl;
                 if (helper.isMaster) cerr << "Check filename and -jFunctionCurve" << endl;
                 usageandexit();
             }
             
             // Invert J-function, now we get saturation as a function of pressure:
             if (Jtmp.isStrictlyMonotone()) {
                 helper.InvJfunctions.push_back(MonotCubicInterpolator(Jtmp.get_fVector(), Jtmp.get_xVector()));
                 JfunctionNames.push_back(ROCKFILENAME);
                 if (helper.upscaleBothPhases) {
                     helper.Krfunctions[0][0].push_back(MonotCubicInterpolator(ROCKFILENAME, 1, 2));
                     helper.Krfunctions[0][1].push_back(MonotCubicInterpolator(ROCKFILENAME, 1, 3));
                 }
                 else {
                     helper.Krfunctions[0][0].push_back(MonotCubicInterpolator(ROCKFILENAME, 1, relPermCurve));
                 }
             }
             else {
                 if (helper.isMaster) cerr << "Error: Jfunction " << i+1 << " in rock file " << ROCKFILENAME << " was not invertible." << endl;
                 usageandexit();
             }
         }
         else {  // If input is anisotropic, then we are in second mode with different input file format
             MonotCubicInterpolator Pctmp;
             try {
                 Pctmp = MonotCubicInterpolator(ROCKFILENAME, 2, 1);
             }
             catch (const char * errormessage) {
                 if (helper.isMaster) cerr << "Error: " << errormessage << endl;
                 if (helper.isMaster) cerr << "Check filename and columns 1 and 2 (Pc and " << helper.saturationstring <<")" << endl;
                 usageandexit();
             }
             
             // Invert Pc(Sw) curve into Sw(Pc):
              if (Pctmp.isStrictlyMonotone()) {
                 helper.SwPcfunctions.push_back(MonotCubicInterpolator(Pctmp.get_fVector(), Pctmp.get_xVector()));
                 JfunctionNames.push_back(ROCKFILENAME);
                 helper.Krfunctions[0][0].push_back(MonotCubicInterpolator(ROCKFILENAME, 2, 3));
                 helper.Krfunctions[1][0].push_back(MonotCubicInterpolator(ROCKFILENAME, 2, 4));
                 helper.Krfunctions[2][0].push_back(MonotCubicInterpolator(ROCKFILENAME, 2, 5));
                 if (helper.upscaleBothPhases) {
                     helper.Krfunctions[0][1].push_back(MonotCubicInterpolator(ROCKFILENAME, 2, 6));
                     helper.Krfunctions[1][1].push_back(MonotCubicInterpolator(ROCKFILENAME, 2, 7));
                     helper.Krfunctions[2][1].push_back(MonotCubicInterpolator(ROCKFILENAME, 2, 8));
                 }
              }
             else {
                 if (helper.isMaster) cerr << "Error: Pc(" << helper.saturationstring << ") curve " << i+1 << " in rock file " << ROCKFILENAME << " was not invertible." << endl;
                 usageandexit();
             }
         }
      } 
   }
   // The code below loads the same file once for every rock type in
   // the file. This is stone_types-1 more than strictly necessary, so
   // it could have been simplified.
   else if (varnum == rockfileindex + 1) {
       const char* ROCKFILENAME = vararg[rockfileindex];
       // Check if rock file exists and is readable:
       ifstream rockfile(ROCKFILENAME, ios::in);
       if (rockfile.fail()) {
           if (helper.isMaster) cerr << "Error: Filename " << ROCKFILENAME << " not found or not readable." << endl;
           usageandexit();
       }
       rockfile.close(); 
       if (! helper.anisotropic_input) {
           MonotCubicInterpolator Jtmp; 
           try {
               Jtmp = MonotCubicInterpolator(ROCKFILENAME, 1, jFunctionCurve);
           }
           catch (const char * errormessage) {
               if (helper.isMaster) cerr << "Error: " << errormessage << endl;
               if (helper.isMaster) cerr << "Check filename and -jFunctionCurve" << endl;
               usageandexit();
           }
           if (Jtmp.isStrictlyMonotone()) {
               for (int i=0; i < stone_types; ++i) {
                   // Invert J-function, now we get saturation as a function of pressure:
                   helper.InvJfunctions.push_back(MonotCubicInterpolator(Jtmp.get_fVector(), Jtmp.get_xVector()));
                   JfunctionNames.push_back(vararg[rockfileindex]);
                   if (helper.upscaleBothPhases) {
                       helper.Krfunctions[0][0].push_back(MonotCubicInterpolator(vararg[rockfileindex], 1, 2));
                       helper.Krfunctions[0][1].push_back(MonotCubicInterpolator(vararg[rockfileindex], 1, 3));
                   }
                   else {
                       helper.Krfunctions[0][0].push_back(MonotCubicInterpolator(vararg[rockfileindex], 1, relPermCurve));
                   }
               }
           }
           else {
               if (helper.isMaster) cerr << "Error: Jfunction " << 1 << " in rock file " << ROCKFILENAME << " was not invertible." << endl;
               usageandexit();
           }
       }
       else {
           MonotCubicInterpolator Pctmp;
           try {
               Pctmp = MonotCubicInterpolator(ROCKFILENAME, 2, 1);
           }
           catch (const char * errormessage) {
               if (helper.isMaster) cerr << "Error: " << errormessage << endl;
               if (helper.isMaster) cerr << "Check filename and columns 1 and 2 (Pc and " << helper.saturationstring <<")" << endl;
               usageandexit();
           }
           // Invert Pc(Sw) curve into Sw(Pc):
           if (Pctmp.isStrictlyMonotone()) {
               for (int i=0; i < stone_types; ++i) { 
                   helper.SwPcfunctions.push_back(MonotCubicInterpolator(Pctmp.get_fVector(), Pctmp.get_xVector()));
                   JfunctionNames.push_back(ROCKFILENAME);
                   helper.Krfunctions[0][0].push_back(MonotCubicInterpolator(ROCKFILENAME, 2, 3));
                   helper.Krfunctions[1][0].push_back(MonotCubicInterpolator(ROCKFILENAME, 2, 4));
                   helper.Krfunctions[2][0].push_back(MonotCubicInterpolator(ROCKFILENAME, 2, 5));
                   if (helper.upscaleBothPhases) {
                       helper.Krfunctions[0][1].push_back(MonotCubicInterpolator(ROCKFILENAME, 2, 6));
                       helper.Krfunctions[1][1].push_back(MonotCubicInterpolator(ROCKFILENAME, 2, 7));
                       helper.Krfunctions[2][1].push_back(MonotCubicInterpolator(ROCKFILENAME, 2, 8));
                   }
               }
           }
           else {
               if (helper.isMaster) cerr << "Error: Pc(" << helper.saturationstring << ") curve " << 1 << " in rock file " << ROCKFILENAME << " was not invertible." << endl;
               usageandexit();
           }           
       }
   }
   else {
       if (helper.isMaster) cerr << "Error:  Wrong number of stone-functions provided. " << endl;
       usageandexit();
   }
   
   // Check if input relperm curves satisfy Eclipse requirement of specifying critical saturations
   helper.doEclipseCheck = (options["doEclipseCheck"] == "true");
   helper.critRelpThresh = atof(options["critRelpermThresh"].c_str());
   if (helper.doEclipseCheck)
       helper.checkCriticalSaturations();

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

   double linsolver_tolerance = atof(options["linsolver_tolerance"].c_str());
   timeused_tesselation = helper.tesselateGrid(deck, options);

   vector<double> dP;
   double dPmin = +DBL_MAX;
   double dPmax = -DBL_MAX;

   /* If gravity is to be included, calculate z-values of every cell: */
   if (includeGravity) {
       // height of model is calculated as the average of the z-values at the top layer
       // This calculation makes assumption on the indexing of cells in the grid, going from bottom to top.
       double modelHeight = 0;
       for (unsigned int zIdx = (4 * res[0] * res[1] * (2*res[2]-1)); zIdx < helper.zcorns.size(); ++zIdx) {
           modelHeight += helper.zcorns[zIdx] / (4*res[0]*res[1]);
       }
       
       // We assume that the spatial units in the grid file is in centimetres, 
       // so we divide by 100 to get to metres.
       modelHeight = modelHeight/100.0; 

       // Input water and oil density is given in g/cm3, we convert it to kg/m3 (SI) 
       // by multiplying with 1000.
       double dRho = (waterDensity-oilDensity) * 1000; // SI unit (kg/m3)

       // Calculating difference in capillary pressure for all cells
       dP = vector<double>(helper.satnums.size(), 0);
       for (unsigned int cellIdx = 0; cellIdx < helper.satnums.size(); ++cellIdx) {
           int i,j,k; // Position of cell in cell hierarchy
           vector<int> zIndices(8,0); // 8 corners with 8 heights
           int horIdx = (cellIdx+1) - int(std::floor(((double)(cellIdx+1))/((double)(res[0]*res[1]))))*res[0]*res[1]; // index in the corresponding horizon
           if (horIdx == 0) { 
               horIdx = res[0]*res[1];
           }
           i = horIdx - int(std::floor(((double)horIdx)/((double)res[0])))*res[0];
           if (i == 0) { 
               i = res[0];
           }
           j = (horIdx-i)/res[0]+1;
           k = ((cellIdx+1)-res[0]*(j-1)-1)/(res[0]*res[1])+1;
           int zBegin = 8*res[0]*res[1]*(k-1); // indices of Z-values of bottom
           int level2 = 4*res[0]*res[1]; // number of z-values in one horizon
           zIndices[0] = zBegin + 4*res[0]*(j-1)+2*i-1;
           zIndices[1] = zBegin + 4*res[0]*(j-1)+2*i;
           zIndices[2] = zBegin + 2*res[0]*(2*j-1)+2*i;
           zIndices[3] = zBegin + 2*res[0]*(2*j-1)+2*i-1;
           zIndices[4] = zBegin + level2 + 4*res[0]*(j-1)+2*i-1;
           zIndices[5] = zBegin + level2 + 4*res[0]*(j-1)+2*i;
           zIndices[6] = zBegin + level2 + 2*res[0]*(2*j-1)+2*i;
           zIndices[7] = zBegin + level2 + 2*res[0]*(2*j-1)+2*i-1;
           
           double cellDepth = 0;
           for (unsigned int corner = 0; corner < 8; ++corner) {
               cellDepth += helper.zcorns[zIndices[corner]-1] / 8.0;
           }
           // cellDepth is in cm, convert to m by dividing by 100
           cellDepth = cellDepth / 100.0;
           dP[cellIdx] = dRho * gravity * (cellDepth-modelHeight/2.0);

           // assume distances in grid are given in cm. 
           dPmin = min(dPmin,dP[cellIdx]); 
           dPmax = max(dPmax,dP[cellIdx]);
           
       }
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

   if (maxPermContrast == 0) {
       if (helper.isMaster) cout << "Illegal contrast value" << endl;
       usageandexit();
   }

   vector<double> cellVolumes;
   cellVolumes.resize(helper.satnums.size(), 0.0);
   helper.cellPoreVolumes.resize(helper.satnums.size(), 0.0);


   /* Find minimium and maximum capillary pressure values in each
      cell, and use the global min/max as the two initial pressure
      points for computations.
   
      Also find max single-phase permeability, used to obey the 
      maxPermContrast option.

      Also find properly upscaled saturation endpoints, these are
      printed out to stdout for reference during computations, but will 
      automatically appear as the lowest and highest saturation points
      in finished output.
   */
   helper.tesselatedCells = 0; // for counting "active" cells (Sintef interpretation of "active")
   helper.Pcmax = -DBL_MAX, helper.Pcmin = DBL_MAX;
   double maxSinglePhasePerm = 0;
   double Swirvolume = 0;
   double Sworvolume = 0;
   // cell_idx is the eclipse index.
   const std::vector<int>& ecl_idx = helper.upscaler.grid().globalCell();
   Dune::CpGrid::Codim<0>::LeafIterator c = helper.upscaler.grid().leafbegin<0>();
   for (; c != helper.upscaler.grid().leafend<0>(); ++c) {
       unsigned int cell_idx = ecl_idx[c->index()];
       if (helper.satnums[cell_idx] > 0) { // Satnum zero is "no rock"

	   cellVolumes[cell_idx] = c->geometry().volume();
	   helper.cellPoreVolumes[cell_idx] = cellVolumes[cell_idx] * helper.poros[cell_idx];
	   
	   double Pcmincandidate, Pcmaxcandidate, minSw, maxSw;
           
	   if (! helper.anisotropic_input) {
	       Pcmincandidate = helper.InvJfunctions[int(helper.satnums[cell_idx])-1].getMinimumX().first
		   / sqrt(helper.perms[0][cell_idx] * milliDarcyToSqMetre / helper.poros[cell_idx]);
	       Pcmaxcandidate = helper.InvJfunctions[int(helper.satnums[cell_idx])-1].getMaximumX().first
		   / sqrt(helper.perms[0][cell_idx] * milliDarcyToSqMetre/helper.poros[cell_idx]);
	       minSw = helper.InvJfunctions[int(helper.satnums[cell_idx])-1].getMinimumF().second;
	       maxSw = helper.InvJfunctions[int(helper.satnums[cell_idx])-1].getMaximumF().second;
	   }
	   else { // anisotropic input, we do not to J-function scaling
	       Pcmincandidate = helper.SwPcfunctions[int(helper.satnums[cell_idx])-1].getMinimumX().first;
	       Pcmaxcandidate = helper.SwPcfunctions[int(helper.satnums[cell_idx])-1].getMaximumX().first;
               
	       minSw = helper.SwPcfunctions[int(helper.satnums[cell_idx])-1].getMinimumF().second;
	       maxSw = helper.SwPcfunctions[int(helper.satnums[cell_idx])-1].getMaximumF().second;
	   }
	   helper.Pcmin = min(Pcmincandidate, helper.Pcmin);
	   helper.Pcmax = max(Pcmaxcandidate, helper.Pcmax);
           
	   maxSinglePhasePerm = max( maxSinglePhasePerm, helper.perms[0][cell_idx]);
           
	   //cout << "minSwc: " << minSw << endl;
	   //cout << "maxSwc: " << maxSw << endl;
           
	   // Add irreducible water saturation volume
	   Swirvolume += minSw * helper.cellPoreVolumes[cell_idx];
	   Sworvolume += maxSw * helper.cellPoreVolumes[cell_idx];
       }
       ++helper.tesselatedCells; // keep count.
   }
   helper.minSinglePhasePerm = max(maxSinglePhasePerm/maxPermContrast, minPerm);
   
   
   if (includeGravity) {
       helper.Pcmin = helper.Pcmin - dPmax;
       helper.Pcmax = helper.Pcmax - dPmin;
   }

   if (helper.isMaster) cout << "Pcmin:    " << helper.Pcmin << endl;
   if (helper.isMaster) cout << "Pcmax:    " << helper.Pcmax << endl;

   if (helper.Pcmin > helper.Pcmax) {
       if (helper.isMaster) cerr << "ERROR: No legal capillary pressures found for this system. Exiting..." << endl;
       usageandexit();
   }

   // Total porevolume and total volume -> upscaled porosity:
   helper.poreVolume = std::accumulate(helper.cellPoreVolumes.begin(),
                                       helper.cellPoreVolumes.end(),
                                       0.0);
   helper.volume = std::accumulate(cellVolumes.begin(),
                                   cellVolumes.end(),
                                   0.0);

   helper.Swir = Swirvolume/helper.poreVolume;
   helper.Swor = Sworvolume/helper.poreVolume;

   if (helper.isMaster) {
       cout << "LF Pore volume:    " << helper.poreVolume << endl;
       cout << "LF Volume:         " << helper.volume << endl;
       cout << "Upscaled porosity: " << helper.poreVolume/helper.volume << endl;
       cout << "Upscaled " << helper.saturationstring << "ir:     " << helper.Swir << endl;
       cout << "Upscaled " << helper.saturationstring << "max:    " << helper.Swor << endl; //Swor=1-Swmax
       cout << "Saturation points to be computed: " << helper.points << endl;
   }

   // Sometimes, if Swmax=1 or Swir=0 in the input tables, the upscaled 
   // values can be a little bit larger (within machine precision) and
   // the check below fails. Hence, check if these values are within the 
   // the [0 1] interval within some precision (use linsolver_precision)
   if (helper.Swor > 1.0 && helper.Swor - linsolver_tolerance < 1.0) {
       helper.Swor = 1.0;
   }
   if (helper.Swir < 0.0 && helper.Swir + linsolver_tolerance > 0.0) {
       helper.Swir = 0.0;
   }
   if (helper.Swir < 0 || helper.Swir > 1 || helper.Swor < 0 || helper.Swor > 1) {
       if (helper.isMaster) cerr << "ERROR: " << helper.saturationstring << "ir/" << helper.saturationstring << "or unsensible. Check your input. Exiting";
       usageandexit();
   }      
   

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

   double largestSaturationInterval = helper.Swor-helper.Swir;
   
   double Ptestvalue = helper.Pcmax;
   
   while (largestSaturationInterval > (helper.Swor-helper.Swir)/500.0) {
       //       cout << Ptestvalue << endl;
       if (helper.Pcmax == helper.Pcmin) {
           // This is a dummy situation, we go through once and then 
           // we are finished (this will be triggered by zero permeability)
           Ptestvalue = helper.Pcmin;
           largestSaturationInterval = 0;
       }
       else if (helper.WaterSaturationVsCapPressure.getSize() == 0) {
           /* No data values previously computed */
           Ptestvalue = helper.Pcmax;
       }
       else if (helper.WaterSaturationVsCapPressure.getSize() == 1) {
           /* If only one point has been computed, it was for Pcmax. So now
              do Pcmin */
           Ptestvalue = helper.Pcmin;
       }
       else {
           /* Search for largest saturation interval in which there are no
              computed saturation points (and estimate the capillary pressure
              that will fall in the center of this saturation interval)
           */
           pair<double,double> SatDiff = helper.WaterSaturationVsCapPressure.getMissingX();
           Ptestvalue = SatDiff.first;
           largestSaturationInterval = SatDiff.second;
       }
       
       // Check for saneness of Ptestvalue:
       if (std::isnan(Ptestvalue) || std::isinf(Ptestvalue)) {
           if (helper.isMaster) cerr << "ERROR: Ptestvalue was inf or nan" << endl;
           break; // Jump out of while-loop, just print out the results
           // up to now and exit the program
       }
       
       double waterVolume = 0.0;
       for (unsigned int i = 0; i < ecl_idx.size(); ++i) {
           unsigned int cell_idx = ecl_idx[i];
               double waterSaturationCell = 0.0;
               if (helper.satnums[cell_idx] > 0) { // handle "no rock" cells with satnum zero
                   double PtestvalueCell;
                   if (includeGravity) {
                       PtestvalueCell = Ptestvalue - dP[cell_idx];
                   }
                   else {
                       PtestvalueCell = Ptestvalue;
                   }
                   if (! helper.anisotropic_input ) {   
                       double Jvalue = sqrt(helper.perms[0][cell_idx] * milliDarcyToSqMetre /helper.poros[cell_idx]) * PtestvalueCell;
                       //cout << "JvalueCell: " << Jvalue << endl;
                       waterSaturationCell 
                           = helper.InvJfunctions[int(helper.satnums[cell_idx])-1].evaluate(Jvalue);
                   }
                   else { // anisotropic_input, then we do not do J-function-scaling
                       waterSaturationCell = helper.SwPcfunctions[int(helper.satnums[cell_idx])-1].evaluate(PtestvalueCell);
                       //cout << Ptestvalue << "\t" <<  waterSaturationCell << endl;
                   }
               }
               waterVolume += waterSaturationCell  * helper.cellPoreVolumes[cell_idx];
       }
       helper.WaterSaturationVsCapPressure.addPair(Ptestvalue, waterVolume/helper.poreVolume);
   }
   //   cout << WaterSaturationVsCapPressure.toString();

   // Now, it may happen that we have a large number of cells, and
   // some cells with near zero poro and perm. This may cause that
   // Pcmax has been estimated so high that it does not affect Sw
   // within machine precision, and then we need to truncate the
   // largest Pc values:
   helper.WaterSaturationVsCapPressure.chopFlatEndpoints(saturationThreshold);

   // Now we can also invert the upscaled water saturation
   // (it should be monotonic)
   if (!helper.WaterSaturationVsCapPressure.isStrictlyMonotone()) {
       if (helper.isMaster) {
           cerr << "Error: Upscaled water saturation not strictly monotone in capillary pressure." << endl;
           cerr << "       Unphysical input data, exiting." << endl;
           cerr << "       Trying to dump " << helper.saturationstring << " vs Pc to file swvspc_debug.txt for inspection" << endl;
           ofstream outfile; 
           outfile.open("swvspc_debug.txt", ios::out | ios::trunc); 
           outfile << "# Pc      " << helper.saturationstring << endl;
           outfile << helper.WaterSaturationVsCapPressure.toString();
           outfile.close();   
       }
       usageandexit();
   }
   MonotCubicInterpolator CapPressureVsWaterSaturation(helper.WaterSaturationVsCapPressure.get_fVector(),
                                                       helper.WaterSaturationVsCapPressure.get_xVector());
   
   /*****************************************************************************
    * Step 7:
    * Upscale single phase permeability
    * This uses the PERMX in the eclipse file as data, and upscales using
    * fixed boundary (no-flow) conditions
    *
    * In an MPI-environment, this is only done on the master node.
    */

   helper.upscaleSinglePhasePermeability();
   typedef SinglePhaseUpscaler::permtensor_t Matrix;

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

   // Put correct number of zeros in, just to be able to access RelPerm[index] later 
   for (int idx=0; idx < helper.points; ++idx) {
       helper.WaterSaturation.push_back(0.0);
       vector<double> tmp;
       helper.PhasePerm[0].push_back(tmp);
       for (int voigtIdx=0; voigtIdx < helper.tensorElementCount; ++voigtIdx) { 
           helper.PhasePerm[0][idx].push_back(0.0);
       } 
       if (helper.upscaleBothPhases){
           helper.PhasePerm[1].push_back(tmp);
           for (int voigtIdx=0; voigtIdx < helper.tensorElementCount; ++voigtIdx) { 
               helper.PhasePerm[1][idx].push_back(0.0);
           }            
       }
   }

   // Make vector of capillary pressure points corresponding to uniformly distribued
   // saturation points between Swor and Swir.
   
   for (int pointidx = 1; pointidx <= helper.points; ++pointidx) {
       // pointidx=1 corresponds to Swir, pointidx=points to Swor.
       double saturation = helper.Swir + (helper.Swor-helper.Swir)/(helper.points-1)*(pointidx-1);
       helper.pressurePoints.push_back(CapPressureVsWaterSaturation.evaluate(saturation));
   }
   // Preserve max and min pressures
   helper.pressurePoints[0]=helper.Pcmax;
   helper.pressurePoints[helper.pressurePoints.size()-1]=helper.Pcmin;

   // Fill with zeros initially (in case of non-mpi)
   for (int idx=0; idx < helper.points; ++idx) {
       helper.node_vs_pressurepoint.push_back(0);
   }
   
#if HAVE_MPI
   // Distribute work load over mpi nodes.
   for (int idx=0; idx < points; ++idx) {
       // Ensure master node gets equal or less work than the other nodes, since
       // master node also computes single phase perm.
       node_vs_pressurepoint[idx] = (mpi_nodecount-1) - idx % mpi_nodecount;
       /*if (isMaster) {
           cout << "Pressure point " << idx << " assigned to node " << node_vs_pressurepoint[idx] << endl;
       }*/
   }   
#endif


   clock_t start_upscale_wallclock = clock();

   double waterVolumeLF;
   // Now loop through the vector of capillary pressure points that
   // this node should compute.
   for (int pointidx = 0; pointidx < helper.points; ++pointidx) {

       // Should "I" (mpi-wise) compute this pressure point?
       if (helper.node_vs_pressurepoint[pointidx] == mpi_rank) {
           
           Ptestvalue = helper.pressurePoints[pointidx];
           
           double accPhasePerm = 0.0;
           double accPhase2Perm = 0.0;
           
           double maxPhasePerm = 0.0;
           double maxPhase2Perm = 0.0;
           
           vector<double> phasePermValues, phase2PermValues;
           vector<vector<double> > phasePermValuesDiag, phase2PermValuesDiag;
           phasePermValues.resize(helper.satnums.size());
           phasePermValuesDiag.resize(helper.satnums.size());
           if (helper.upscaleBothPhases) {
               phase2PermValues.resize(helper.satnums.size());
               phase2PermValuesDiag.resize(helper.satnums.size());
           }
           waterVolumeLF = 0.0;
           for (unsigned int i = 0; i < ecl_idx.size(); ++i) {
               unsigned int cell_idx = ecl_idx[i];
               double cellPhasePerm = minPerm;
               double cellPhase2Perm = minPerm;
               vector<double>  cellPhasePermDiag, cellPhase2PermDiag;
               cellPhasePermDiag.push_back(minPerm);
               cellPhasePermDiag.push_back(minPerm);
               cellPhasePermDiag.push_back(minPerm);
               if (helper.upscaleBothPhases) {
                   cellPhase2PermDiag.push_back(minPerm);
                   cellPhase2PermDiag.push_back(minPerm);
                   cellPhase2PermDiag.push_back(minPerm);
               }
               
               if (helper.satnums[cell_idx] > 0) { // handle "no rock" cells with satnum zero
                   //            cout << endl << "Cell no. " << cell_idx << endl;
                   double PtestvalueCell;
                   if (includeGravity) {
                       PtestvalueCell = Ptestvalue - dP[cell_idx];
                   }
                   else {
                       PtestvalueCell = Ptestvalue;
                   }
                   
                   if (! helper.anisotropic_input) {
                       
                       double Jvalue = sqrt(helper.perms[0][cell_idx] * milliDarcyToSqMetre/helper.poros[cell_idx]) * PtestvalueCell;
                       //cout << "JvalueCell: " << Jvalue << endl;
                       double waterSaturationCell 
                           = helper.InvJfunctions[int(helper.satnums[cell_idx])-1].evaluate(Jvalue);
                       waterVolumeLF += waterSaturationCell * helper.cellPoreVolumes[cell_idx];
                       
                       // Compute cell relative permeability. We use a lower cutoff-value as we
                       // easily divide by zero here.  When water saturation is
                       // zero, we get 'inf', which is circumvented by the cutoff value.
                       cellPhasePerm = 
                           helper.Krfunctions[0][0][int(helper.satnums[cell_idx])-1].evaluate(waterSaturationCell) *
                           helper.perms[0][cell_idx];
                       if (helper.upscaleBothPhases) {
                           cellPhase2Perm = 
                               helper.Krfunctions[0][1][int(helper.satnums[cell_idx])-1].evaluate(waterSaturationCell) *
                               helper.perms[0][cell_idx];
                       }
                   }
                   else {
                       double waterSaturationCell = helper.SwPcfunctions[int(helper.satnums[cell_idx])-1].evaluate(PtestvalueCell);
                       //cout << PtestvalueCell << "\t" << waterSaturationCell << endl;
                       waterVolumeLF += waterSaturationCell * helper.cellPoreVolumes[cell_idx];
                       
                       cellPhasePermDiag[0] = helper.Krfunctions[0][0][int(helper.satnums[cell_idx])-1].evaluate(waterSaturationCell) *
                           helper.perms[0][cell_idx];
                       cellPhasePermDiag[1] = helper.Krfunctions[1][0][int(helper.satnums[cell_idx])-1].evaluate(waterSaturationCell) *
                           helper.perms[1][cell_idx];
                       cellPhasePermDiag[2] = helper.Krfunctions[2][0][int(helper.satnums[cell_idx])-1].evaluate(waterSaturationCell) *
                           helper.perms[2][cell_idx];
                       if (helper.upscaleBothPhases) {
                           cellPhase2PermDiag[0] = helper.Krfunctions[0][1][int(helper.satnums[cell_idx])-1].evaluate(waterSaturationCell) *
                               helper.perms[0][cell_idx];
                           cellPhase2PermDiag[1] = helper.Krfunctions[1][1][int(helper.satnums[cell_idx])-1].evaluate(waterSaturationCell) *
                               helper.perms[1][cell_idx];
                           cellPhase2PermDiag[2] = helper.Krfunctions[2][1][int(helper.satnums[cell_idx])-1].evaluate(waterSaturationCell) *
                               helper.perms[2][cell_idx];
                       }
                   }
                   
                   phasePermValues[cell_idx] = cellPhasePerm;
                   phasePermValuesDiag[cell_idx] = cellPhasePermDiag;
                   maxPhasePerm = max(maxPhasePerm, cellPhasePerm);
                   maxPhasePerm = max(maxPhasePerm, *max_element(cellPhasePermDiag.begin(),
                                                                 cellPhasePermDiag.end()));
                   if (helper.upscaleBothPhases) {
                       phase2PermValues[cell_idx] = cellPhase2Perm;
                       phase2PermValuesDiag[cell_idx] = cellPhase2PermDiag;
                       maxPhase2Perm = max(maxPhase2Perm, cellPhase2Perm);
                       maxPhase2Perm = max(maxPhase2Perm, *max_element(cellPhase2PermDiag.begin(),
                                                                     cellPhase2PermDiag.end()));
                   }
               }
           }
           // Now we can determine the smallest permitted permeability we can calculate for

           // We have both a fixed bottom limit, as well as a possible higher limit determined
           // by a maximum allowable permeability.
           double minPhasePerm = max(maxPhasePerm/maxPermContrast, minPerm);
           double minPhase2Perm;
           if (helper.upscaleBothPhases) minPhase2Perm = max(maxPhase2Perm/maxPermContrast, minPerm);

           // Now remodel the phase permeabilities obeying minPhasePerm
           Matrix cellperm(3,3,nullptr);
           zero(cellperm);
           for (unsigned int i = 0; i < ecl_idx.size(); ++i) {
               unsigned int cell_idx = ecl_idx[i];
               zero(cellperm);
               if (! helper.anisotropic_input) {
                   double cellPhasePerm = max(minPhasePerm, phasePermValues[cell_idx]);
                   accPhasePerm += cellPhasePerm;
                   double kval = max(minPhasePerm, cellPhasePerm);
                   cellperm(0,0) = kval;
                   cellperm(1,1) = kval;
                   cellperm(2,2) = kval;
               }
               else { // anisotropic_input
                   // Truncate values lower than minPhasePerm upwards.
                   phasePermValuesDiag[cell_idx][0] = max(minPhasePerm, phasePermValuesDiag[cell_idx][0]);
                   phasePermValuesDiag[cell_idx][1] = max(minPhasePerm, phasePermValuesDiag[cell_idx][1]);
                   phasePermValuesDiag[cell_idx][2] = max(minPhasePerm, phasePermValuesDiag[cell_idx][2]);
                   accPhasePerm += phasePermValuesDiag[cell_idx][0]; // not correct anyway                   
                   cellperm(0,0) = phasePermValuesDiag[cell_idx][0];
                   cellperm(1,1) = phasePermValuesDiag[cell_idx][1];
                   cellperm(2,2) = phasePermValuesDiag[cell_idx][2];
               }
               helper.upscaler.setPermeability(i, cellperm);
           }
           
           // Output average phase perm, this is just a reality check so that we are not way off.
           //cout << ", Arith. mean phase perm = " << accPhasePerm/float(tesselatedCells) << " mD, ";
           
           //  Call single-phase upscaling code 
           Matrix phasePermTensor = helper.upscaler.upscaleSinglePhase();
           
           // Now upscale phase permeability for phase 2
           Matrix phase2PermTensor;
           if (helper.upscaleBothPhases) {
               zero(cellperm);
               for (unsigned int i = 0; i < ecl_idx.size(); ++i) {
                   unsigned int cell_idx = ecl_idx[i];
                   zero(cellperm);
                   if (! helper.anisotropic_input) {
                       double cellPhase2Perm = max(minPhase2Perm, phase2PermValues[cell_idx]);
                       accPhase2Perm += cellPhase2Perm;
                       double kval = max(minPhase2Perm, cellPhase2Perm);
                       cellperm(0,0) = kval;
                       cellperm(1,1) = kval;
                       cellperm(2,2) = kval;
                   }
                   else { // anisotropic_input
                       // Truncate values lower than minPhasePerm upwards.
                       phase2PermValuesDiag[cell_idx][0] = max(minPhase2Perm, phase2PermValuesDiag[cell_idx][0]);
                       phase2PermValuesDiag[cell_idx][1] = max(minPhase2Perm, phase2PermValuesDiag[cell_idx][1]);
                       phase2PermValuesDiag[cell_idx][2] = max(minPhase2Perm, phase2PermValuesDiag[cell_idx][2]);
                       accPhase2Perm += phase2PermValuesDiag[cell_idx][0]; // not correct anyway                   
                       cellperm(0,0) = phase2PermValuesDiag[cell_idx][0];
                       cellperm(1,1) = phase2PermValuesDiag[cell_idx][1];
                       cellperm(2,2) = phase2PermValuesDiag[cell_idx][2];
                   }
                   helper.upscaler.setPermeability(i, cellperm);
               }
               phase2PermTensor = helper.upscaler.upscaleSinglePhase();
           }

           //cout << phasePermTensor << endl;


           // Here we recalculate the upscaled water saturation,
           // although it is already known when we asked for the
           // pressure point to compute for. Nonetheless, we
           // recalculate here to avoid any minor roundoff-error and
           // interpolation error (this means that the saturation
           // points are not perfectly uniformly distributed)
           helper.WaterSaturation[pointidx] =  waterVolumeLF/helper.poreVolume;

           
#ifdef HAVE_MPI
           cout << "Rank " << mpi_rank << ": ";
#endif
           cout << Ptestvalue << "\t" << helper.WaterSaturation[pointidx];
           // Store and print phase-perm-result
           for (int voigtIdx=0; voigtIdx < helper.tensorElementCount; ++voigtIdx) {
               helper.PhasePerm[0][pointidx][voigtIdx] = getVoigtValue(phasePermTensor, voigtIdx);
               cout << "\t" << getVoigtValue(phasePermTensor, voigtIdx);
               if (helper.upscaleBothPhases){
                   helper.PhasePerm[1][pointidx][voigtIdx] = getVoigtValue(phase2PermTensor, voigtIdx);
                   cout << "\t" << getVoigtValue(phase2PermTensor, voigtIdx);
               }
           } 
           cout << endl;
       }
   }
   
   clock_t finish_upscale_wallclock = clock();
   timeused_upscale_wallclock = (double(finish_upscale_wallclock)-double(start_upscale_wallclock))/CLOCKS_PER_SEC;

   helper.collectResults();

   // Average time pr. upscaling point:
#ifdef HAVE_MPI
   // Sum the upscaling time used by all processes
   double timeused_total;
   MPI_Reduce(&timeused_upscale_wallclock, &timeused_total, 1, MPI_DOUBLE, 
              MPI_SUM, 0, MPI_COMM_WORLD);
   double avg_upscaling_time_pr_point = timeused_total/(double)helper.points;
       
#else
   double avg_upscaling_time_pr_point = timeused_upscale_wallclock / (double)helper.points;
#endif


   /* 
    * Step 8c: Make relperm values from phaseperms
    *          (only master node can do this)
    */
   std::array<vector<vector<double>>,2> RelPermValues;
   RelPermValues[0] = helper.getRelPerm(0);
   if (helper.upscaleBothPhases)
       RelPermValues[1] = helper.getRelPerm(1);

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
       outputtmp << "######################################################################" << endl;
       outputtmp << "# Results from upscaling relative permeability."<< endl;
       outputtmp << "#" << endl;
#if HAVE_MPI
       outputtmp << "#          (MPI-version)" << endl;
#endif
       time_t now = std::time(NULL);
       outputtmp << "# Finished: " << asctime(localtime(&now));
       
       utsname hostname;   uname(&hostname);
       outputtmp << "# Hostname: " << hostname.nodename << endl;
       
       outputtmp << "#" << endl;
       outputtmp << "# Eclipse file: " << ECLIPSEFILENAME << endl;
       outputtmp << "#        cells: " << helper.tesselatedCells << endl;
       outputtmp << "#  Pore volume: " << helper.poreVolume << endl;
       outputtmp << "#       volume: " << helper.volume << endl;
       outputtmp << "#     Porosity: " << helper.poreVolume/helper.volume << endl;
       outputtmp << "#" << endl;
       if (! helper.anisotropic_input) {
           for (int i=0; i < stone_types ; ++i) {
               outputtmp << "# Stone " << i+1 << ": " << JfunctionNames[i] << " (" << helper.InvJfunctions[i].getSize() << " points)" <<  endl;
           }
           outputtmp << "#         jFunctionCurve: " << options["jFunctionCurve"] << endl;
           if (!helper.upscaleBothPhases) outputtmp << "#           relPermCurve: " << options["relPermCurve"] << endl;
       }
       else { // anisotropic input, not J-functions that are supplied on command line (but vector JfunctionNames is still used)
           for (int i=0; i < stone_types ; ++i) {
               outputtmp << "# Stone " << i+1 << ": " << JfunctionNames[i] << " (" << helper.Krfunctions[0][0][i].getSize() << " points)" <<  endl;
           }
       }
           
       outputtmp << "#" << endl;
       outputtmp << "# Timings:   Tesselation: " << timeused_tesselation << " secs" << endl;
       outputtmp << "#              Upscaling: " << timeused_upscale_wallclock << " secs";
#ifdef HAVE_MPI
       outputtmp << " (wallclock time)" << endl;
       outputtmp << "#                         " << avg_upscaling_time_pr_point << " secs pr. saturation point" << endl;
       outputtmp << "#              MPI-nodes: " << mpi_nodecount << endl;

       // Single phase upscaling time is included here, in possibly a hairy way.
       double speedup = (avg_upscaling_time_pr_point * (points + 1) + timeused_tesselation)/(timeused_upscale_wallclock + avg_upscaling_time_pr_point + timeused_tesselation);
       outputtmp << "#                Speedup: " << speedup << ", efficiency: " << speedup/mpi_nodecount << endl;
#else
       outputtmp << ", " << avg_upscaling_time_pr_point << " secs avg for " << helper.points << " runs" << endl;
#endif
       outputtmp << "# " << endl;
       outputtmp << "# Options used:" << endl;
       outputtmp << "#     Boundary conditions: ";
       if (isFixed)    outputtmp << "Fixed (no-flow)" << endl;
       if (isPeriodic) outputtmp << "Periodic" << endl;
       if (isLinear)   outputtmp << "Linear" << endl;
       outputtmp << "#                  points: " << options["points"] << endl;
       outputtmp << "#         maxPermContrast: " << options["maxPermContrast"] << endl;
       outputtmp << "#                 minPerm: " << options["minPerm"] << endl;
       outputtmp << "#                 minPoro: " << options["minPoro"] << endl;
       outputtmp << "#          surfaceTension: " << options["surfaceTension"] << " dynes/cm" << endl;
       if (includeGravity) {
           outputtmp << "#                 gravity: " << options["gravity"] << " m/s²" << endl;
           if (owsystem) outputtmp << "#            waterDensity: " << options["waterDensity"] << " g/cm³" << endl;
           else outputtmp << "#              gasDensity: " << options["waterDensity"] << " g/cm³" << endl;
           outputtmp << "#              oilDensity: " << options["oilDensity"] << " g/cm³" << endl;
       }
       else {
           outputtmp << "#                 gravity: 0" << endl;
       }
       if (doInterpolate) {
           outputtmp << "#             interpolate: " << options["interpolate"] << " points" << endl;
       }
       outputtmp << "# " << endl;
       outputtmp << "# Single phase permeability" << endl; 
       outputtmp << "#  |Kxx  Kxy  Kxz| = " << helper.permTensor(0,0) << "  " << helper.permTensor(0,1) << "  " << helper.permTensor(0,2) << endl; 
       outputtmp << "#  |Kyx  Kyy  Kyz| = " << helper.permTensor(1,0) << "  " << helper.permTensor(1,1) << "  " << helper.permTensor(1,2) << endl; 
       outputtmp << "#  |Kzx  Kzy  Kzz| = " << helper.permTensor(2,0) << "  " << helper.permTensor(2,1) << "  " << helper.permTensor(2,2) << endl; 
       outputtmp << "# " << endl;
       if (doInterpolate) {
           outputtmp << "# NB: Data points shown are interpolated." << endl;
       }
       outputtmp << "######################################################################" << endl;
       if (helper.upscaleBothPhases) {
           string phase1, phase2;
           if (owsystem) phase1="w"; else phase1="g";
           phase2="o";
           if (isFixed) { 
               outputtmp << "#  Pc (Pa)        " << helper.saturationstring << "           Kr" << phase1 << "xx       Kr" << phase1 << "yy       Kr" << phase1 << "zz"
                         <<  "       Kr" << phase2 << "xx       Kr" << phase2 << "yy       Kr" << phase2 << "zz" <<  endl; 
           } 
           else if (isPeriodic || isLinear) { 
               outputtmp << "#  Pc (Pa)        " << helper.saturationstring << "           Kr" << phase1 << "xx       Kr" << phase1 << "yy       Kr" << phase1 << "zz       Kr"
                         << phase1 << "yz       Kr" << phase1 << "xz       Kr" << phase1 << "xy       Kr" << phase1 << "zy       Kr" << phase1 << "zx       Kr" << phase1 << "yx" 
                         << "       Kr" << phase2 << "xx       Kr" << phase2 << "yy       Kr" << phase2 << "zz       Kr" 
                         << phase2 << "yz       Kr" << phase2 << "xz       Kr" << phase2 << "xy       Kr" << phase2 << "zy       Kr" << phase2 << "zx       Kr" << phase2 << "yx" << endl;
           }
       }
       else {
           if (isFixed) { 
               outputtmp << "#  Pc (Pa)        " << helper.saturationstring << "            Krxx        Kryy        Krzz" << endl;
           } 
           else if (isPeriodic || isLinear) { 
               outputtmp << "#  Pc (Pa)        " << helper.saturationstring << "            Krxx        Kryy        Krzz        Kryz        Krxz        Krxy        Krzy        Krzx        Kryx" << endl;
           }
       }
       
       vector<double> Pvalues = helper.pressurePoints;

       // Multiply all pressures with the surface tension (potentially) supplied
       // at the command line. This multiplication has been postponed to here
       // to avoid division by zero and to avoid special handling of negative
       // capillary pressure in the code above.
       std::transform(Pvalues.begin(), Pvalues.end(), Pvalues.begin(), 
		      std::bind1st(std::multiplies<double>(), surfaceTension));
       vector<double> Satvalues = helper.WaterSaturation; //.get_fVector(); 
       
       // If user wants interpolated output, do monotone cubic interpolation
       // by modifying the data vectors that are to be printed
       if (doInterpolate) {
           // Find min and max for saturation values
           double xmin = +DBL_MAX;
           double xmax = -DBL_MAX;
           for (unsigned int i = 0; i < Satvalues.size(); ++i) {
               if (Satvalues[i] < xmin) {
                   xmin = Satvalues[i];
               }
               if (Satvalues[i] > xmax) {
                   xmax = Satvalues[i];
               }
           }
           // Make uniform grid in saturation axis
           vector<double> SatvaluesInterp;
           for (int i = 0; i < interpolationPoints; ++i) {
               SatvaluesInterp.push_back(xmin + ((double)i)/((double)interpolationPoints-1)*(xmax-xmin));
           }
           // Now capillary pressure and computed relperm-values must be viewed as functions
           // of saturation, and then interpolated on the uniform saturation grid.
           
           // Now overwrite existing Pvalues and relperm-data with interpolated data:
           MonotCubicInterpolator PvaluesVsSaturation(Satvalues, Pvalues);
           Pvalues.clear();
           for (int i = 0; i < interpolationPoints; ++i) {
               Pvalues.push_back(PvaluesVsSaturation.evaluate(SatvaluesInterp[i]));
           }
           for (int voigtIdx = 0; voigtIdx < helper.tensorElementCount; ++voigtIdx) {
               MonotCubicInterpolator RelPermVsSaturation(Satvalues, RelPermValues[0][voigtIdx]);
               RelPermValues[0][voigtIdx].clear();
               for (int i=0; i < interpolationPoints; ++i) {
                   RelPermValues[0][voigtIdx].push_back(RelPermVsSaturation.evaluate(SatvaluesInterp[i]));
               }
           }
           if (helper.upscaleBothPhases) {
               for (int voigtIdx = 0; voigtIdx < helper.tensorElementCount; ++voigtIdx) {
                   MonotCubicInterpolator RelPermVsSaturation(Satvalues, RelPermValues[1][voigtIdx]);
                   RelPermValues[1][voigtIdx].clear();
                   for (int i=0; i < interpolationPoints; ++i) {
                       RelPermValues[1][voigtIdx].push_back(RelPermVsSaturation.evaluate(SatvaluesInterp[i]));
                   }
               }
           }
           
           // Now also overwrite Satvalues
           Satvalues.clear();
           Satvalues = SatvaluesInterp;
       }
       
       // The code below does not care whether the data is interpolated or not.
       const int fieldwidth = outputprecision + 8;
       for (unsigned int i=0; i < Satvalues.size(); ++i) {
           outputtmp << showpoint << setw(fieldwidth) << setprecision(outputprecision) << Pvalues[i]; 
           outputtmp << showpoint << setw(fieldwidth) << setprecision(outputprecision) << Satvalues[i]; 
           
           for (int voigtIdx = 0; voigtIdx < helper.tensorElementCount; ++voigtIdx) {
               outputtmp << showpoint << setw(fieldwidth) << setprecision(outputprecision) 
                         << RelPermValues[0][voigtIdx][i]; 
           } 
           if (helper.upscaleBothPhases) {
               for (int voigtIdx = 0; voigtIdx < helper.tensorElementCount; ++voigtIdx) {
                   outputtmp << showpoint << setw(fieldwidth) << setprecision(outputprecision) 
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
       
       // If both phases are upscaled and output is specyfied, create SWOF or SGOF files for Eclipse
       if (options["output"] != "" && helper.upscaleBothPhases) {
           // krow(swirr)-values if given
           double krowxswirr = atof(options["krowxswirr"].c_str());
           double krowyswirr = atof(options["krowyswirr"].c_str());
           double krowzswirr = atof(options["krowzswirr"].c_str());

           stringstream swofx, swofy, swofz;
           string satstringCap = ""; if (owsystem) satstringCap = "W"; else satstringCap = "G"; 
           string satstring = ""; if (owsystem) satstring = "w"; else satstring = "g"; 
           // x-direction
           swofx << "-- This file is based on the results in " << endl 
                 << "-- " << options["output"] << endl
                 << "-- for relperm in x-direction." << endl
                 << "-- Pressure values (Pc) given in bars." << endl
                 << "--        S" << satstring << "       Kr" << satstring << "xx      Kro" << satstring << "xx      Pc(bar)" << endl
                 << "--S" << satstringCap << "OF" << endl;
           if (krowxswirr > 0){
               swofx << showpoint << setw(fieldwidth) << setprecision(outputprecision) << 0
                     << showpoint << setw(fieldwidth) << setprecision(outputprecision) << 0
                     << showpoint << setw(fieldwidth) << setprecision(outputprecision) << krowxswirr
                     << showpoint << setw(fieldwidth) << setprecision(outputprecision) << 0 << endl;
           }
           for (unsigned int i=0; i < Satvalues.size(); ++i) {
               swofx << showpoint << setw(fieldwidth) << setprecision(outputprecision) << Satvalues[i]
                     << showpoint << setw(fieldwidth) << setprecision(outputprecision) << RelPermValues[0][0][i]
                     << showpoint << setw(fieldwidth) << setprecision(outputprecision) << RelPermValues[1][0][i]
                     << showpoint << setw(fieldwidth) << setprecision(outputprecision) << Pvalues[i]/100000.0 << endl;
           }
           swofx << "/" << endl;
           // y-direction
           swofy << "-- This file is based on the results in " << endl 
                 << "-- " << options["output"] << endl
                 << "-- for relperm in y-direction." << endl
                 << "-- Pressure values (Pc) given in bars." << endl
                 << "--        S" << satstring << "       Kr" << satstring << "yy      Kro" << satstring << "yy      Pc(bar)" << endl
                 << "--S" << satstringCap << "OF" << endl;
           if (krowyswirr > 0){
               swofy << showpoint << setw(fieldwidth) << setprecision(outputprecision) << 0
                     << showpoint << setw(fieldwidth) << setprecision(outputprecision) << 0
                     << showpoint << setw(fieldwidth) << setprecision(outputprecision) << krowyswirr
                     << showpoint << setw(fieldwidth) << setprecision(outputprecision) << 0 << endl;
           }
           for (unsigned int i=0; i < Satvalues.size(); ++i) {
               swofy << showpoint << setw(fieldwidth) << setprecision(outputprecision) << Satvalues[i]
                     << showpoint << setw(fieldwidth) << setprecision(outputprecision) << RelPermValues[0][1][i]
                     << showpoint << setw(fieldwidth) << setprecision(outputprecision) << RelPermValues[1][1][i]
                     << showpoint << setw(fieldwidth) << setprecision(outputprecision) << Pvalues[i]/100000.0 << endl;
           }
           swofy << "/" << endl;
           // z-direction
           swofz << "-- This file is based on the results in " << endl 
                 << "-- " << options["output"] << endl
                 << "-- for relperm in z-direction." << endl
                 << "-- Pressure values (Pc) given in bars." << endl
                 << "--        S" << satstring << "       Kr" << satstring << "zz      Kro" << satstring << "zz      Pc(bar)" << endl
                 << "--S" << satstringCap << "OF" << endl;
           if (krowzswirr > 0){
               swofz << showpoint << setw(fieldwidth) << setprecision(outputprecision) << 0
                     << showpoint << setw(fieldwidth) << setprecision(outputprecision) << 0
                     << showpoint << setw(fieldwidth) << setprecision(outputprecision) << krowzswirr
                     << showpoint << setw(fieldwidth) << setprecision(outputprecision) << 0 << endl;
           }
           for (unsigned int i=0; i < Satvalues.size(); ++i) {
               swofz << showpoint << setw(fieldwidth) << setprecision(outputprecision) << Satvalues[i]
                     << showpoint << setw(fieldwidth) << setprecision(outputprecision) << RelPermValues[0][2][i]
                     << showpoint << setw(fieldwidth) << setprecision(outputprecision) << RelPermValues[1][2][i]
                     << showpoint << setw(fieldwidth) << setprecision(outputprecision) << Pvalues[i]/100000.0 << endl;
           }
           swofz << "/" << endl;
           //cout << swofx.str() << endl;
           //cout << swofy.str() << endl;
           //cout << swofz.str() << endl;
           ofstream xfile, yfile, zfile;
           string opfname = options["output"];
           string fnbase = opfname.substr(0,opfname.find_first_of('.'));
           string xfilename = fnbase + "-x.S" + satstringCap + "OF";
           string yfilename = fnbase + "-y.S" + satstringCap + "OF";
           string zfilename = fnbase + "-z.S" + satstringCap + "OF";
           
           cout << "Writing Eclipse compatible files to " << xfilename << ", " << yfilename << " and " << zfilename << endl;
           xfile.open(xfilename.c_str(), ios::out | ios::trunc);
           xfile << swofx.str();
           xfile.close();
           yfile.open(yfilename.c_str(), ios::out | ios::trunc);
           yfile << swofy.str();
           yfile.close();
           zfile.open(zfilename.c_str(), ios::out | ios::trunc);
           zfile << swofz.str();
           zfile.close();
       }

   }

   return 0;
}
catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    usageandexit();
}
