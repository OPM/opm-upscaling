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
   @file upscale_cond.C
   @brief Upscales resistivity as a function of water saturation
   
   Description:  
   
   Reads in a lithofacies geometry in Eclipse format, reads in J(S_w)
   for each stone type, and calculates upscaled 
   resistivity values for values of S_w.
   
   Steps in the code:
   
   1: Process command line options.
   2: Read and parse Eclipse file 
   3: Read a J-function for each stone-type.
   4: Tesselate the grid (Sintef code)
   5: Find minimum and maximum capillary pressure from the J-functions in each cell.
   6: Upscale water saturation as a function of capillary pressure
   7: Upscale conductivity ( = 1/resistivity ) for each saturation point
   8: Print output to screen and optionally to file.
   
   The relative permeability computation is based on 
     - Capillary equilibrium, p_c is spatially invariant.

   Units handling:
     - Input surface tension is in dynes/cm
     - The denominator \sigma * cos(\phi) in J-function scaling
       is what we call "surface tension". If angle dependency is to be
       included, calculate the "surface tension" yourself.
     - Outputted capillary pressure is in Pascals.

 */
#include <config.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <cfloat> // for DBL_MAX
#include <map>
#include <sys/utsname.h>

#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
#include <dune/common/parallel/mpihelper.hh>
#else
#include <dune/common/mpihelper.hh>
#endif

#include <opm/core/utility/MonotCubicInterpolator.hpp>
#include <opm/upscaling/SinglePhaseUpscaler.hpp>
 
using namespace Opm;
using namespace std;



/**
   Explains how to use the file. Shows possible option parameters.
 */
void usage()
{
    cout << "Usage: upscale_cond [<options>] <eclipsefile> Jfunc1.data [Jfunc2.data [...]]" << endl << 
        "where the options are:" << endl <<
        "  -bc <string>                    -- which boundary conditions to use." << endl <<
        "                                     Possible values are f (fixed), " << endl <<
        "                                     l (linear) and p (periodic). Default f." << endl <<
        "  -points <integer>               -- Number of saturation points to compute for. " <<
        "                                     Default 30" << endl <<
        "  -resistivityCutoff <float>      -- Default 10000. Not necessary to touch" << endl <<
        "  -lithologyCoeff <float>         -- Archie parameter. Default 1.0" << endl <<
        "  -cementationExponent <float>    -- Archie parameter. Default 1.8" << endl <<
        "  -saturationExponent <float>     -- Archie parameter. Default 2.3" << endl <<
        "  -waterresistivity               -- Resistivity of formation water." << endl <<
        "  -output <string>                -- Output to file as well as screen if" << endl <<
        "                                     provided." << endl <<
        "                                     Default none" << endl <<
        "  -mudresistivity <float>         -- Constant resistivity for mud. Default 1.4" << endl <<
        "  -mud1rocktype <int>             -- First mud rock type." << endl << 
        "                                     0 if there are no mud types." << endl <<
        "                                     Default 0" << endl <<
        "  -mud2rocktype <int>             -- Second mud rock type. Default 0" << endl <<
        "  -jFunctionCurve <int>           -- Which column in the supplied J-files" << endl <<
        "                                     corresponds to the J-function values." << endl <<
        "                                     Default 2." << endl <<
        "  -surfaceTension <float>         -- Surface tension to use in J-function/Pc conversion." << endl << 
        "                                     Default 11 dynes/cm (oil-water systems). In absence of" << endl <<  
        "                                     a correct value, the surface tension for gas-oil systems " << endl << 
        "                                     could be 22.5 dynes/cm." << endl <<                 
        "  -interpolate <integer>          -- If supplied and > 1, the output data points will be" << endl << 
        "                                     interpolated using monotone cubic interpolation" << endl << 
        "                                     on a uniform grid with the specified number of" << endl << 
        "                                     points. Suggested value: 1000." << endl <<         "" << endl <<
        "  -rock<int>cemexp <float>        -- Cementation exponent can be set on a per rocktype basis" << endl <<
        "  -rock<int>satexp <float>        -- Saturation exponent can be set on a per rocktype basis" << endl <<
        "Jfunctions are data files with two colums of numbers. The first column is water " << endl <<
        "saturation, the second is the J-function. The first Jfunction, Jfunc1.data" << endl <<
        "corresponds to the first rock type defined in the eclipsefile's SATNUM. The" << endl <<
        "second correspond to the second rock type and so on. If just one Jfunc is" << endl <<
        "given, this is used for all rock types" << endl;
    // "minPoro" intentionally left undocumented

}

void usageandexit() {
    usage();
    exit(1);
}
// Assumes that permtensor_t use C ordering.
double getVoigtValue(const SinglePhaseUpscaler::permtensor_t& K, int voigt_idx)
{
    assert(K.numRows() == 3 && K.numCols() == 3);
    switch (voigt_idx) {
    case 0: return K.data()[0];
    case 1: return K.data()[4];
    case 2: return K.data()[8];
    case 3: return K.data()[5];
    case 4: return K.data()[2];
    case 5: return K.data()[1];
    case 6: return K.data()[7];
    case 7: return K.data()[6];
    case 8: return K.data()[3];
    default:
        std::cout << "Voigt index out of bounds (only 0-8 allowed)" << std::endl;
        throw std::exception();
    }
}


// Assumes that permtensor_t use C ordering.
void setVoigtValue(SinglePhaseUpscaler::permtensor_t& K, int voigt_idx, double val)
{
    assert(K.numRows() == 3 && K.numCols() == 3);
    switch (voigt_idx) {
    case 0: K.data()[0] = val; break;
    case 1: K.data()[4] = val; break;
    case 2: K.data()[8] = val; break;
    case 3: K.data()[5] = val; break;
    case 4: K.data()[2] = val; break;
    case 5: K.data()[1] = val; break;
    case 6: K.data()[7] = val; break;
    case 7: K.data()[6] = val; break;
    case 8: K.data()[3] = val; break;
    default:
        std::cout << "Voigt index out of bounds (only 0-8 allowed)" << std::endl;
        throw std::exception();
    }
}

int main(int varnum, char** vararg)
try
{ 
   // Variables used for timing/profiling:
   clock_t start, finish;
   double timeused = 0, timeused_tesselation = 0;
   double timeused_upscale_wallclock = 0.0;

   Dune::MPIHelper& mpi=Dune::MPIHelper::instance(varnum, vararg);
   const int mpi_rank = mpi.rank();
#ifdef HAVE_MPI
   const int mpi_nodecount = mpi.size();
#endif
   bool isMaster = (mpi_rank == 0);
   if (varnum == 1) { /* If no arguments supplied ("upscale_cond" is the first "argument") */
       usageandexit();
   }


   
   /*
     Populate options-map with default values
   */
   map<string,string> options;
   options.insert(make_pair("bc",                     "f"     )); // Fixed boundary conditions
   options.insert(make_pair("points",                 "30"    )); // Number of saturation points
   options.insert(make_pair("resistivityCutoff",      "10000" )); // we do not need higher resistivity
   options.insert(make_pair("lithologyCoeff",         "1.0"   )); // Archie parameter
   options.insert(make_pair("cementationExponent",    "1.8"   )); // Archie parameter
   options.insert(make_pair("saturationExponent",     "2.3"   )); // Archie parameter
   options.insert(make_pair("output",                 ""      )); // Output to file as well as screen if provided
   options.insert(make_pair("mudresistivity",         "1.4"   )); // Constant resistivity for mud
   options.insert(make_pair("mud1rocktype",           "0"     )); // First mud rock type. 0 if there is no mud types
   options.insert(make_pair("mud2rocktype",           "0"     )); // Second mud rock type. 
   options.insert(make_pair("jFunctionCurve",         "4"     )); // column number of J-function values in input files
   options.insert(make_pair("surfaceTension",         "11"    )); // Surface tension given in dynes/cm
   options.insert(make_pair("interpolate",            "0"     )); // default is not to interpolate  
   options.insert(make_pair("minPerm",                "1e-12" )); // minimum modelled permeability (for saturation distr)
   options.insert(make_pair("minPoro",            "0.0001")); // this limit is necessary for pcmin/max computation

   // dune-cornerpoint specific options
   options.insert(make_pair("linsolver_tolerance", "1e-12"));  // residual tolerance for linear solver
   options.insert(make_pair("linsolver_verbosity", "0"));     // verbosity level for linear solver
   options.insert(make_pair("linsolver_type",      "3"));     // type of linear solver: 0 = ILU/BiCGStab, 1 = AMG/CG, 2 = KAMG/CG, 3 = FastAMG/CG
   options.insert(make_pair("linsolver_prolongate_factor", "1.0")); // Factor to scale the prolongate coarse grid correction
   options.insert(make_pair("linsolver_smooth_steps", "1")); // Number of pre and postsmoothing steps for AMG

   /*
     Extra options for CT-experiments. If you encounter a rock with more than 6
     rocktypes, additional lines only need to be added here, the rest of the code adapts.
   */
   options.insert(make_pair("waterresistivity", "1.32")); // Formation water resistivity in Ohm meters.
   options.insert(make_pair("rock1cemexp", "0")); // Cementation exponent for SATNUM 1
   options.insert(make_pair("rock1satexp", "0")); // Saturation exponent for SATNUM 1
   options.insert(make_pair("rock2cemexp", "0")); // etc..
   options.insert(make_pair("rock2satexp", "0")); // A value of zero means it that these
   options.insert(make_pair("rock3cemexp", "0")); // options is not in use, and then the global option
   options.insert(make_pair("rock3satexp", "0")); // above will be used.
   options.insert(make_pair("rock4cemexp", "0"));
   options.insert(make_pair("rock4satexp", "0"));
   options.insert(make_pair("rock5cemexp", "0"));
   options.insert(make_pair("rock5satexp", "0"));
   options.insert(make_pair("rock6cemexp", "0"));
   options.insert(make_pair("rock6satexp", "0"));

   // Conversion factor, multiply mD numbers with this to get m² numbers
   const double milliDarcyToSqMetre = 9.869233e-16;
   // Reference: http://www.spe.org/spe-site/spe/spe/papers/authors/Metric_Standard.pdf
 
   /*
     Look for strings in args matching the entries in the options map, 
     if found, replace default values with command line values.
   */

   /* Check first if there is anything on the command line to look for */
   if (varnum == 1) {
       if (isMaster) cout << "Error: No Eclipsefile or J-functions found on command line." << endl;
      usageandexit();
   }

   /* Loop over all command line options in order to look 
      for options. 

      argidx loops over all the arguments here, and updates the
      variable 'argeclindex' *if* it finds any legal options,
      'argeclindex' is so that vararg[argeclindex] = the eclipse
      filename. If options are illegal, argeclindex will be wrong, 
      
   */
   int argeclindex = 0;
   for (int argidx = 1; argidx < varnum; argidx += 2)  {
       if (string(vararg[argidx]).substr(0,1) == "-")    {
           string searchfor = string(vararg[argidx]).substr(1); // Chop off leading '-'
           /* Check if it is a match */
           if (options.count(searchfor) == 1) {
               options[searchfor] = string(vararg[argidx+1]);
               if (isMaster) cout << "Parsed command line option: " << searchfor << " := " << vararg[argidx+1] << endl;
               argeclindex = argidx + 2;
           }
           else {
               if (isMaster) cout << "Option -" << searchfor << " unrecognized." << endl;
               usageandexit();
           }
       }
       else { 
           // if vararg[argidx] does not start in '-', 
           // assume we have found the position of the Eclipse-file.
           argeclindex = argidx;
           break; // out of for-loop, 
       }
   }
   
   
   // argeclindex should now point at the eclipse file
   static char* ECLIPSEFILENAME(vararg[argeclindex]);
   argeclindex += 1; // areglcindex jumps to next input argument, now it points to the first J-file
   
   // argcindex now points to the first J-function. This index is not
   // to be touched now.
   static int JFindex = argeclindex;
   
   /* Check if at least one J-function is supplied on command line */
   if (varnum <= JFindex) {
      if (isMaster) cerr << "Error: No J-functions found on command line." << endl;
      usageandexit();
   }

   // Check that boundary conditions are valid , and make booleans 
   // for boundary conditions. This allows more readable code later.
   bool isFixed = false, isLinear = false, isPeriodic = false;
   SinglePhaseUpscaler::BoundaryConditionType boundaryCondition = SinglePhaseUpscaler::Fixed ; 

   int tensorElementCount = 0; // Number of independent elements in resulting tensor
   if (options["bc"].substr(0,1) == "f") {
       isFixed = true; isLinear = isPeriodic = false;
       boundaryCondition = SinglePhaseUpscaler::Fixed; // This refers to the mimetic namespace (Sintef)
       tensorElementCount = 3; // Diagonal
   }
   else if (options["bc"].substr(0,1) == "l") {
       isLinear = true; isFixed = isPeriodic = false;
       boundaryCondition = SinglePhaseUpscaler::Linear;
       tensorElementCount = 9;
   }
   else if (options["bc"].substr(0,1) == "p") {
       isPeriodic = true; isFixed = isLinear = false;
       boundaryCondition = SinglePhaseUpscaler::Periodic;
       tensorElementCount = 9; // This should be symmetric, so 6 should suffice, but
       // the code calculates 9 elements that are supposed to be symmetric.
   }
   else {
       if (isMaster) cout << "Invalid boundary condition: " << options["bc"] << endl;
       usageandexit();
   }
   
   // Read data from the Eclipse file and 
   // populate our vectors with data from the file
   // Test if filename exists and is readable
   ifstream eclipsefile(ECLIPSEFILENAME, ios::in);
   if (eclipsefile.fail()) {
       if (isMaster) cerr << "Error: Filename " << ECLIPSEFILENAME << " not found or not readable." << endl;
       usageandexit();
   }
   eclipsefile.close(); 

   if (isMaster) cout << "Parsing Eclipse file <" << ECLIPSEFILENAME << "> ... ";
   flush(cout);   start = clock();
   Opm::EclipseGridParser eclParser(ECLIPSEFILENAME, false);
   finish = clock();   timeused = (double(finish)-double(start))/CLOCKS_PER_SEC;
   if (isMaster) cout << " (" << timeused <<" secs)" << endl;

   // Check that we have the information we need from the eclipse file:
   if (! (eclParser.hasField("SPECGRID") && eclParser.hasField("COORD") && eclParser.hasField("ZCORN")
          && eclParser.hasField("PORO") && eclParser.hasField("PERMX") && eclParser.hasField("SATNUM"))) {
       if (isMaster) cerr << "Error: Did not find SPECGRID, COORD and ZCORN in Eclipse file " << ECLIPSEFILENAME << endl;
       usageandexit();
   }

   const int points                = atoi(options["points"].c_str());

   vector<int> satnums   = eclParser.getIntegerValue("SATNUM");
   vector<double> poros  = eclParser.getFloatingPointValue("PORO");
   vector<double> permxs = eclParser.getFloatingPointValue("PERMX");
   const double minPerm = atof(options["minPerm"].c_str()); 
   const double minPoro = atof(options["minPoro"].c_str());


   // Read in J-functions for all stone-types.
   // Number of stone-types is max(satnums):
   
   // If there is only one J-function supplied on the command line,
   // use that for all stone types.

   int stone_types = int(*(max_element(satnums.begin(), satnums.end())));
   std::vector<MonotCubicInterpolator> InvJfunctions; // Holds the inverse of the loaded J-functions.
   std::vector<string> JfunctionNames; // Placeholder for the names of the loaded J-functions.

   // Handle two command line input formats, either one J-function for all stone types
   // or one each. If there is only one stone type, both code blocks below are equivalent.
   const int jFunctionCurve = atoi(options["jFunctionCurve"].c_str()); // default 4
   // Input for surfaceTension is dynes/cm
   // SI units are Joules/square metre
   const double surfaceTension     = atof(options["surfaceTension"].c_str()) * 1e-3; // multiply 
   if (varnum == JFindex + stone_types) {
      for (int i=0 ; i < stone_types; ++i) {
         const char* ROCKFILENAME = vararg[JFindex+i];
         // Check if rock file exists and is readable:
         ifstream rockfile(ROCKFILENAME, ios::in);
         if (rockfile.fail()) {
             if (isMaster) cerr << "Error: Filename " << ROCKFILENAME << " not found or not readable." << endl;
            usageandexit();
         }
         rockfile.close(); 
         
         MonotCubicInterpolator Jtmp;
         try {
             Jtmp = MonotCubicInterpolator(vararg[JFindex + i], 1, jFunctionCurve); 
         }
         catch (const char * errormessage) {
             if (isMaster) cerr << "Error: " << errormessage << endl;
             if (isMaster) cerr << "Check filename and -jFunctionCurve" << endl;
             usageandexit();
         }
         
         // Invert J-function, now we get saturation as a function of pressure:
         if (Jtmp.isStrictlyMonotone()) {
            InvJfunctions.push_back(MonotCubicInterpolator(Jtmp.get_fVector(), Jtmp.get_xVector()));
            JfunctionNames.push_back(vararg[JFindex + i]);
            //cout << vararg[JFindex + i] << endl;
            //cout << InvJfunctions[i].toString();
         }
         else {
             if (isMaster) cerr << "Error: Jfunction " << i+1 << " in rock file " << ROCKFILENAME << " was not invertible." << endl;
             usageandexit();
         }
      } 
   }
   else if (varnum == JFindex + 1) {
      for (int i=0; i < stone_types; ++i) {
        const char* ROCKFILENAME = vararg[JFindex];
         // Check if rock file exists and is readable:
         ifstream rockfile(ROCKFILENAME, ios::in);
         if (rockfile.fail()) {
             if (isMaster) cerr << "Error: Filename " << ROCKFILENAME << " not found or not readable." << endl;
	     usageandexit();
         }
         rockfile.close(); 

         MonotCubicInterpolator Jtmp;
         try {
             Jtmp = MonotCubicInterpolator(vararg[JFindex], 1, jFunctionCurve);
         }
         catch (const char * errormessage) {
             if (isMaster) cerr << "Error: " << errormessage << endl;
             if (isMaster) cerr << "Check filename and -jFunctionCurve" << endl;
             usageandexit();
         }

         // Invert J-function, now we get saturation as a function of pressure:
         if (Jtmp.isStrictlyMonotone()) {
            InvJfunctions.push_back(MonotCubicInterpolator(Jtmp.get_fVector(), Jtmp.get_xVector()));
            JfunctionNames.push_back(vararg[JFindex]);
         }
         else {
             if (isMaster) cerr << "Error: Jfunction " << i+1 << " in rock file " << ROCKFILENAME << " was not invertible." << endl;
            usageandexit();
         }
      }
   }
   else {
       if (isMaster) cerr << "Error:  Wrong number of J-functions provided. " << endl << 
               "Note that all input arguments after the eclipsefile " << endl << 
               "are interpreted as J-functions." << endl;
       usageandexit();
   }


   /*
     Step 4
     
     Generate tesselated grid. Tesselation depends on boundary conditions.
     For periodic boundary conditions, the grid needs to be massaged slightly
     (crop top and bottom). These modifications ruin the computations for 
     fixed and linear boundary conditions.
   */
   vector<int>  griddims = eclParser.getSPECGRID().dimensions;
   SinglePhaseUpscaler upscaler;
   double ztol = 0.0;
   double linsolver_tolerance = atof(options["linsolver_tolerance"].c_str());
   int linsolver_verbosity    = atoi(options["linsolver_verbosity"].c_str());
   int linsolver_type         = atoi(options["linsolver_type"].c_str());
   bool twodim_hack = false;
   eclParser.convertToSI();
   upscaler.init(eclParser, boundaryCondition,
                 Opm::unit::convert::from(minPerm, Opm::prefix::milli*Opm::unit::darcy),
                 ztol, linsolver_tolerance, linsolver_verbosity, linsolver_type, twodim_hack);

   finish = clock();   timeused_tesselation = (double(finish)-double(start))/CLOCKS_PER_SEC;
   if (isMaster) cout << " (" << timeused_tesselation <<" secs)" << endl;    
   
  
   int maxSatnum = 0;
   int tesselatedCells = 0;

   /* Sanity check/fix on input for each cell:
      - Check that SATNUM are set sensibly, that is => 0 and < 1000, error if not.
      - Check that porosity is between 0 and 1, error if not.
        Set to minPoro if zero or less than minPoro (due to pcmin/max computation)
      - Check that permeability is zero or positive. Error if negative. Set to minPerm if
        zero and less than minPerm.
      - Check maximum number of SATNUM values (can be number of rock types present)
   */
   for (size_t i = 0; i < satnums.size(); ++i) {
       if (satnums[i] < 0 || satnums[i] > 1000) { 
           if (isMaster) cerr << "satnums[" << i << "] = " << satnums[i] << ", not sane, quitting." << endl;
           usageandexit();
       }
       if (satnums[i] > maxSatnum) {
           maxSatnum = satnums[i];
       }
       if ((poros[i] >= 0) && (poros[i] < minPoro)) { // Truncate porosity from below
           poros[i] = minPoro;
       }
       if (poros[i] < 0 || poros[i] > 1) {
           if (isMaster) cerr << "poros[" << i <<"] = " << poros[i] << ", not sane, quitting." << endl;
           usageandexit();
       }
       if (permxs[i] == 0) {
           permxs[i] = minPerm;
       }
       if (permxs[i] < 0) {
           if (isMaster) cerr << "permx[" << i <<"] = " << permxs[i] << ", not sane, quitting." << endl;
           usageandexit();
       }
       // Explicitly handle "no rock" cells, set them to minimum perm and zero porosity.
       if (satnums[i] == 0) {
           permxs[i] = minPerm;
           poros[i] = 0;
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
    */


   vector<double> cellVolumes, cellPoreVolumes; 
   cellVolumes.resize(satnums.size(), 0.0);
   cellPoreVolumes.resize(satnums.size(), 0.0);

   double Pcmax = -DBL_MAX, Pcmin = DBL_MAX;
   double Swirvolume = 0;
   double Sworvolume = 0;
   const std::vector<int>& ecl_idx = upscaler.grid().globalCell();
   Dune::CpGrid::Codim<0>::LeafIterator c = upscaler.grid().leafbegin<0>();
   for (; c != upscaler.grid().leafend<0>(); ++c) {
       size_t cell_idx = ecl_idx[c->index()];
       if (satnums[cell_idx] > 0) { // Satnum zero is "no rock"
           cellVolumes[cell_idx] = c->geometry().volume();
           cellPoreVolumes[cell_idx] = cellVolumes[cell_idx] * poros[cell_idx];
           
           double Pcmincandidate, Pcmaxcandidate, minSw, maxSw;
           
           Pcmincandidate = InvJfunctions[int(satnums[cell_idx])-1].getMinimumX().first
               / sqrt(permxs[cell_idx] * milliDarcyToSqMetre / poros[cell_idx]) * surfaceTension;
           Pcmaxcandidate = InvJfunctions[int(satnums[cell_idx])-1].getMaximumX().first
               / sqrt(permxs[cell_idx] * milliDarcyToSqMetre/poros[cell_idx]) * surfaceTension;
           minSw = InvJfunctions[int(satnums[cell_idx])-1].getMinimumF().second;
           maxSw = InvJfunctions[int(satnums[cell_idx])-1].getMaximumF().second;
           
           Pcmin = min(Pcmincandidate, Pcmin);
           Pcmax = max(Pcmaxcandidate, Pcmax);
           
           // Add irreducible water saturation volume
           Swirvolume += minSw * cellPoreVolumes[cell_idx];
           Sworvolume += maxSw * cellPoreVolumes[cell_idx];
       }
       ++tesselatedCells; // keep count.
   }

   if (isMaster) cout << "Pcmin:    " << Pcmin << endl;
   if (isMaster) cout << "Pcmax:    " << Pcmax << endl;

   if (Pcmin > Pcmax) {
       if (isMaster) cerr << "ERROR: No legal capillary pressures found for this system. Exiting..." << endl;
       usageandexit();
   }

   // Total porevolume and total volume -> upscaled porosity:
   double poreVolume = std::accumulate(cellPoreVolumes.begin(), 
                                       cellPoreVolumes.end(),
                                       0.0);
   double volume = std::accumulate(cellVolumes.begin(),
                                   cellVolumes.end(),
                                   0.0);
   double Swir = Swirvolume/poreVolume;
   double Swor = Sworvolume/poreVolume;

   if (isMaster) {
       cout << "LF Pore volume:    " << poreVolume << endl;
       cout << "LF Volume:         " << volume << endl;
       cout << "Upscaled porosity: " << poreVolume/volume << endl;
       cout << "Upscaled Swir:     " << Swir << endl;
       cout << "Upscaled Swmax:    " << Swor << endl; //Swor=1-Swmax
       cout << "Saturation points to be computed: " << points << endl;
   }

   // Sometimes, if Swmax=1 or Swir=0 in the input tables, the upscaled 
   // values can be a little bit larger (within machine precision) and
   // the check below fails. Hence, check if these values are within the 
   // the [0 1] interval within some precision (use linsolver_precision)
   if (Swor > 1.0 && Swor - linsolver_tolerance < 1.0) {
       Swor = 1.0;
   }
   if (Swir < 0.0 && Swir + linsolver_tolerance > 0.0) {
       Swir = 0.0;
   }
   if (Swir < 0.0 || Swir > 1.0 || Swor < 0.0 || Swor > 1.0) {
       if (isMaster) cerr << "ERROR: Swir/Swor unsensible. Check your input. Exiting";
       usageandexit();
   }      
   

   /************************************************************
    * Handle petrophysical constants/parameters
    */
   const double resistivityCutoff = atof(options["resistivityCutoff"].c_str()); // 10000 ohmm
   const double mudresistivity = atof(options["mudresistivity"].c_str()); // 1.4 ohmm
   const double mud1rocktype = atoi(options["mud1rocktype"].c_str()); // rocktype 0 (DEFAULT IS WRONG)
   const double mud2rocktype = atoi(options["mud2rocktype"].c_str()); // rocktype 0 (DEFAULT IS WRONG)
   const double a_lithologycoeff = atof(options["lithologyCoeff"].c_str()); // 1.0
   const double m_cementation = atof(options["cementationExponent"].c_str()); //1.8
   const double n_saturationexponent = atof(options["saturationExponent"].c_str()); // 2.3
   const double waterresistivity = atof(options["waterresistivity"].c_str());

   vector<double> cementationexponents, saturationexponents;
   cementationexponents.resize(maxSatnum + 1); // index 0 unused!!
   saturationexponents.resize(maxSatnum + 1); // index 0 unused!!
   
   /* The saturation and cementation exponents may be set globally for
      all rock types using the options cementationExponent and
      saturationExponent. However, these can also be set on a rock by
      rock basis, the rock types that are not set explicitly, will be
      given the general value.
   */
   for (size_t i = 1; i <= (size_t)maxSatnum; ++i) {

       stringstream rocktypestring;
       rocktypestring << i; // integer to string conversion.
       
       //cout << "rock" + rocktypestring.str() + "cemexp" << endl;
       double cemoptionvalue = atof(options["rock" + rocktypestring.str() + "cemexp"].c_str());
       double satoptionvalue = atof(options["rock" + rocktypestring.str() + "satexp"].c_str());
       if (cemoptionvalue > 0) {
           cementationexponents[i] = cemoptionvalue;
       }
       else {
           cementationexponents[i] = m_cementation;
       }
       if (satoptionvalue > 0) {
           saturationexponents[i] = satoptionvalue;
       }
       else {
           saturationexponents[i] = n_saturationexponent;
       }
       //cout << "rocktype " << i << endl;
       //cout << " cem: " << cementationexponents[i] << endl;
       //cout << " sat: " << saturationexponents[i] << endl;
       
   }

   if (waterresistivity <= 0) {
       if (isMaster) cout << "Error: Water resistivity must be positive." << endl;
       usageandexit();
   }
                               

   // If this number is 1 or higher, the output will be interpolated, if not 
   // the computed data is untouched. 
   const int interpolationPoints = atoi(options["interpolate"].c_str()); 
   bool doInterpolate = false; 
   if (interpolationPoints > 1) {
       doInterpolate = true;
   }


   /* We skip TVD and let the user supply the temperature via an option instead. */
   // const double trueVerticalDepth = 3000; /* unit meters */
   // const double temperature = 0.0378*trueVerticalDepth - 7.2; /* unit deg C */

   if (isMaster) cout << "Rw:  " << waterresistivity << endl;



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

   MonotCubicInterpolator WaterSaturationVsCapPressure;
   
   double largestSaturationInterval = Swor-Swir;
   
   double Ptestvalue = Pcmax;
   
   while (largestSaturationInterval > (Swor-Swir)/500.0) {
       //       cout << Ptestvalue << endl;
       if (Pcmax == Pcmin) {
           // This is a dummy situation, we go through once and then 
           // we are finished (this will be triggered by zero permeability)
           Ptestvalue = Pcmin;
           largestSaturationInterval = 0;
       }
       else if (WaterSaturationVsCapPressure.getSize() == 0) {
           /* No data values previously computed */
           Ptestvalue = Pcmax;
       }
       else if (WaterSaturationVsCapPressure.getSize() == 1) {
           /* If only one point has been computed, it was for Pcmax. So now
              do Pcmin */
           Ptestvalue = Pcmin;
       }
       else {
           /* Search for largest saturation interval in which there are no
              computed saturation points (and estimate the capillary pressure
              that will fall in the center of this saturation interval)
           */
           pair<double,double> SatDiff = WaterSaturationVsCapPressure.getMissingX();
           Ptestvalue = SatDiff.first;
           largestSaturationInterval = SatDiff.second;
       }
       
       // Check for saneness of Ptestvalue:
       if (std::isnan(Ptestvalue) || std::isinf(Ptestvalue)) {
           if (isMaster) cerr << "ERROR: Ptestvalue was inf or nan" << endl;
           break; // Jump out of while-loop, just print out the results
           // up to now and exit the program
       }
       
       double waterVolume = 0.0;
       Dune::CpGrid::Codim<0>::LeafIterator c = upscaler.grid().leafbegin<0>();
       for (; c != upscaler.grid().leafend<0>(); ++c) {
           size_t cell_idx = ecl_idx[c->index()];
           //           for (size_t cell_idx = 0; cell_idx < satnums.size(); ++cell_idx) {
           //if (LFgrid.getCellIndex(cell_idx) != EMPTY) {
           double waterSaturationCell = 0.0;
           if (satnums[cell_idx] > 0) { // handle "no rock" cells with satnum zero
               double PtestvalueCell;
               PtestvalueCell = Ptestvalue;
               double Jvalue = sqrt(permxs[cell_idx] * milliDarcyToSqMetre /poros[cell_idx]) * PtestvalueCell / surfaceTension;
               //cout << "JvalueCell: " << Jvalue << endl;
               waterSaturationCell 
                   = InvJfunctions[int(satnums[cell_idx])-1].evaluate(Jvalue);
           }
           waterVolume += waterSaturationCell  * cellPoreVolumes[cell_idx];
       }
       
       WaterSaturationVsCapPressure.addPair(Ptestvalue, waterVolume/poreVolume);
   }
   
   //cout << WaterSaturationVsCapPressure.toString();
   
   // Now, it may happen that we have a large number of cells, and
   // some cells with near zero poro and perm. This may cause that
   // Pcmax has been estimated so high that it does not affect Sw
   // within machine precision, and then we need to truncate the
   // largest Pc values:
   WaterSaturationVsCapPressure.chopFlatEndpoints();
   WaterSaturationVsCapPressure.shrinkFlatAreas();
   // Now we can also invert the upscaled water saturation
   // (it should be monotonic)
   if (!WaterSaturationVsCapPressure.isStrictlyMonotone()) {
       if (isMaster) {
           cerr << "Error: Upscaled water saturation not strictly monotone in capillary pressure." << endl;
           cerr << "       Unphysical input data, exiting." << endl;
       }
       usageandexit();
   }
   MonotCubicInterpolator CapPressureVsWaterSaturation(WaterSaturationVsCapPressure.get_fVector(), 
                                                       WaterSaturationVsCapPressure.get_xVector());

   /*****************************************************************
    * Step 7:
    * 
    * Loop through a given number of uniformly distributed saturation points
    * and upscale conductivity for each of them.
    *    a: Make vector of capillary pressure points corresponding to uniformly
    *       distributed water saturation points between saturation endpoints.
    *    b: Loop over capillary pressure points
    *       1)  Loop over all cells to find the saturation value given the 
    *           capillary pressure found in (a). Given the saturation value, find the
    *           conductivity in the cell given its physical properties
    *       2)  Upscale conductivity for the geometry.
    */

   // Empty vectors for computed data. Will be null for some of the data in an MPI-setting

   typedef SinglePhaseUpscaler::permtensor_t Matrix;

   Matrix zeroMatrix(3,3,(double*)0);
   zero(zeroMatrix);
 
   vector<double> WaterSaturation; // This will hold re-upscaled water saturation for the computed pressure points.
   vector<vector<double> > UpscaledConductivity;  // 'tensorElementCount' phaseperm values per pressurepoint.


   // Put correct number of zeros in
   for (int idx=0; idx < points; ++idx) {
       WaterSaturation.push_back(0.0);
       vector<double> tmp;
       UpscaledConductivity.push_back(tmp);
       for (int voigtIdx=0; voigtIdx < tensorElementCount; ++voigtIdx) { 
           UpscaledConductivity[idx].push_back(0.0); 
       } 
   }

   // Make vector of capillary pressure points corresponding to uniformly distributed
   // saturation points between Swor and Swir.
   
   vector<double> pressurePoints;
   for (int pointidx = 1; pointidx <= points; ++pointidx) {
       // pointidx=1 corresponds to Swir, pointidx=points to Swor.
       double saturation = Swir + (Swor-Swir)/(points-1)*(pointidx-1);
       pressurePoints.push_back(CapPressureVsWaterSaturation.evaluate(saturation));
   }


   // Construct a vector that tells for each pressure point which mpi-node (rank) should compute for that
   // particular pressure point
   vector<int> node_vs_pressurepoint;
   // Fill with zeros initially (in case of non-mpi)
   for (int idx=0; idx < points; ++idx) {
       node_vs_pressurepoint.push_back(0);
   }
   
#ifdef HAVE_MPI
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
   
   double waterVolumeLF = 0.0;
   // Now loop through the vector of capillary pressure points that
   // this node should compute.
   for (int pointidx = 0; pointidx < points; ++pointidx) {
       
       double accRes = 0.0;
       // Should "I" (mpi-wise) compute this pressure point?
       if (node_vs_pressurepoint[pointidx] == mpi_rank) {
           
           Ptestvalue = pressurePoints[pointidx];
           
           cout << "Upscaling resistivity for Pc = " << Ptestvalue;
           flush(cout);
           
           start = clock();
           
           // Loop over each cell again to find saturations given this particular
           // capillary pressure:
           waterVolumeLF = 0.0;

           Dune::CpGrid::Codim<0>::LeafIterator c = upscaler.grid().leafbegin<0>();
           for (; c != upscaler.grid().leafend<0>(); ++c) {
               size_t cell_idx = ecl_idx[c->index()];
               //               for (size_t cell_idx = 0; cell_idx < satnums.size(); ++cell_idx) {
               //  if (LFgrid.getCellIndex(cell_idx) != EMPTY) {
               double resistivityCell = resistivityCutoff;
               if (satnums[cell_idx] > 0) { // SATNUM  index == 0 is legal, means "not rock"
                   //cout << endl << "Cell no. " << cell_idx << " satnum: " <<  satnums[cell_idx] << endl;
                   double Jvalue = sqrt(permxs[cell_idx] * milliDarcyToSqMetre / poros[cell_idx]) * Ptestvalue / surfaceTension;
                   //cout << "JvalueCell: " << Jvalue << endl;
                   double waterSaturationCell 
                       = InvJfunctions[int(satnums[cell_idx])-1].evaluate(Jvalue);
                   //cout << "WatersaturationCell: " << waterSaturationCell << endl;
                   waterVolumeLF += waterSaturationCell  * cellPoreVolumes[cell_idx];
                   
                   // Compute cell resistivity. We use a cutoff-value as we
                   // easily divide by zero here.  When water saturation is
                   // zero, we get 'inf', which is circumvented by the cutoff value.
                   
                   if ((satnums[cell_idx] == mud1rocktype) || (satnums[cell_idx] == mud2rocktype)) {
                       // Handle mud specially
                       resistivityCell = mudresistivity;
                   }
                   else {
                       resistivityCell 
                           = min(a_lithologycoeff * waterresistivity / pow(poros[cell_idx], cementationexponents[satnums[cell_idx]]) 
                                 / pow(waterSaturationCell, saturationexponents[satnums[cell_idx]]), 
                                 resistivityCutoff);
                   }
               }              
               // Insert conductivity (reciprocal of resistivity) into
               // the grid for upscaling:
               Matrix cellCond = zeroMatrix;
               double condval = 1.0/resistivityCell;
               cellCond(0,0) = condval;
               cellCond(1,1) = condval;
               cellCond(2,2) = condval;
               upscaler.setPermeability(c->index(), cellCond);
               accRes += resistivityCell;
               //cout << "cell " << cell_idx << " resistivity" << resistivityCell << endl;
           }
       }
       
       
       // Output average resistity
       cout << ", Arith. mean res = " << accRes/float(tesselatedCells) << " ohms, ";
       
       //  Call upscaling code (SINTEF/FRAUNHOFER)
       
       Matrix condTensor;
       condTensor = upscaler.upscaleSinglePhase();
      
       
       
      // Here we recalculate the upscaled water saturation,
      // although it is already known when we asked for the
      // pressure point to compute for. Nonetheless, we
      // recalculate here to avoid any minor roundoff-error and
      // interpolation error (this means that the saturation
      // points are not perfectly uniformly distributed)
      WaterSaturation[pointidx] = waterVolumeLF/poreVolume;
      cout << WaterSaturation[pointidx] << endl;
      
      invert(condTensor);
      Matrix resTensor(condTensor);


      // Store result
      for (int voigtIdx = 0; voigtIdx < tensorElementCount; ++voigtIdx) {
          UpscaledConductivity[pointidx][voigtIdx] = getVoigtValue(resTensor, voigtIdx);
      }

      //      finish = clock();    timeused = (double(finish)-double(start))/CLOCKS_PER_SEC;
      //cout << " (" << timeused <<" secs)" << endl;
      //timeused_upscale_acc += timeused;

      // Output computed values for impatient users..
#ifdef HAVE_MPI
           cout << "Rank " << mpi_rank << ": ";
#endif
      cout << Ptestvalue << "\t" << waterVolumeLF/poreVolume;
      for (int voigtIdx = 0; voigtIdx < tensorElementCount; ++voigtIdx) {
          cout << "\t" << getVoigtValue(resTensor, voigtIdx);
      }
      cout << endl;
   }
   clock_t finish_upscale_wallclock = clock();
   timeused_upscale_wallclock = (double(finish_upscale_wallclock)-double(start_upscale_wallclock))/CLOCKS_PER_SEC;

   //   double timeused_upscale_total = timeused_upscale_wallclock;

      /****** FINISHED WITH MAIN COMPUTATION ******/

#ifdef HAVE_MPI
   /* Step 7b: Transfer all computed data to master node.
      Master node should post a receive for all values missing,
      other nodes should post a send for all the values they have.
    */
   MPI_Barrier(MPI_COMM_WORLD); // Not strictly necessary.
   if (isMaster) {
       // Loop over all values, receive data and put into local data structure
       for (int idx=0; idx < points; ++idx) {
           if (node_vs_pressurepoint[idx] != 0) {
               // Receive data
               double recvbuffer[2+tensorElementCount];
               MPI_Recv(recvbuffer, 2+tensorElementCount, MPI_DOUBLE, 
                        node_vs_pressurepoint[idx], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               // Put received data into correct place.
               WaterSaturation[(int)recvbuffer[0]] = recvbuffer[1];
               for (int voigtIdx=0; voigtIdx < tensorElementCount; ++voigtIdx) {
                   UpscaledConductivity[(int)recvbuffer[0]][voigtIdx] = recvbuffer[2+voigtIdx];
               }
           }
       }
   }
   else {
       for (int idx=0; idx < points; ++idx) {
           if (node_vs_pressurepoint[idx] == mpi_rank) {
               // Pack and send data. C-style.
               double sendbuffer[2+tensorElementCount];
               sendbuffer[0] = (double)idx;
               sendbuffer[1] = WaterSaturation[idx];
               for (int voigtIdx=0; voigtIdx < tensorElementCount; ++voigtIdx) {
                   sendbuffer[2+voigtIdx] = UpscaledConductivity[idx][voigtIdx];
               }
               MPI_Send(sendbuffer, 2+tensorElementCount, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
           }
       }
   }
#endif


   // Average time pr. upscaling point:
#ifdef HAVE_MPI
   // Sum the upscaling time used by all processes
   double timeused_total;
   MPI_Reduce(&timeused_upscale_wallclock, &timeused_total, 1, MPI_DOUBLE, 
              MPI_SUM, 0, MPI_COMM_WORLD);
   double avg_upscaling_time_pr_point = timeused_total/(double)points;
       
#else
   double avg_upscaling_time_pr_point = timeused_upscale_wallclock / (double)points;
#endif

   /*********************************************************************************
    *  Step 8
    *
    * Output results to stdout and optionally to file.
    */
   // If no data computed, we do not have more to do:
   if (WaterSaturation.size() == 0) {
      return(1); // non-zero return value, this means something wrong with input data.
   }
   
   if (isMaster) {
       stringstream outputtmp;
       
       // Print a table of all computed values:
       outputtmp << "######################################################################" << endl;
       outputtmp << "# Results from upscaling resistivity."<< endl;
       outputtmp << "#" << endl;
#ifdef HAVE_MPI
       outputtmp << "#          (MPI-version)" << endl;
#endif
       time_t now = std::time(NULL);
       outputtmp << "# Finished: " << asctime(localtime(&now));
       
       utsname hostname;   uname(&hostname);
       outputtmp << "# Hostname: " << hostname.nodename << endl;
       
       outputtmp << "#" << endl;
       outputtmp << "# Eclipse file: " << ECLIPSEFILENAME << endl;
       outputtmp << "#        cells: " << tesselatedCells << endl;
       outputtmp << "#  Pore volume: " << poreVolume << endl;
       outputtmp << "#       volume: " << volume << endl;
       outputtmp << "#     Porosity: " << poreVolume/volume << endl;
       outputtmp << "#" << endl;
       for (int i=0; i < stone_types ; ++i) {
           outputtmp << "# Stone " << i+1 << ": " << JfunctionNames[i] << " (" << InvJfunctions[i].getSize() << " points)" <<  endl;
       }
       outputtmp << "#" << endl;
             outputtmp << "#" << endl;
       outputtmp << "# Timings:   Tesselation: " << timeused_tesselation << " secs" << endl;
       outputtmp << "#              Upscaling: " << timeused_upscale_wallclock << " secs";
#ifdef HAVE_MPI
       outputtmp << " (wallclock time)" << endl;
       outputtmp << "#                         " << avg_upscaling_time_pr_point << " secs pr. saturation point" << endl;
       outputtmp << "#              MPI-nodes: " << mpi_nodecount << endl;

       double speedup = (avg_upscaling_time_pr_point * (points) + timeused_tesselation)/(timeused_upscale_wallclock + avg_upscaling_time_pr_point + timeused_tesselation);
       outputtmp << "#                Speedup: " << speedup << ", efficiency: " << speedup/mpi_nodecount << endl;
#else
       outputtmp << ", " << avg_upscaling_time_pr_point << " secs avg for " << points << " runs" << endl;
#endif
       outputtmp << "# " << endl;
       outputtmp << "# Options used:" << endl;
       outputtmp << "#      Boundary conditions: ";
       if (isFixed)    outputtmp << "Fixed (no-flow)" << endl;
       if (isPeriodic) outputtmp << "Periodic" << endl;
       if (isLinear)   outputtmp << "Linear" << endl;
       outputtmp << "#        resistivityCutoff: " << options["resistivityCutoff"] << " ohms" << endl;
       outputtmp << "#           mudresistivity: " << options["mudresistivity"] << " ohms" << endl;
       outputtmp << "#             mud1rocktype: " << options["mud1rocktype"] << endl;
       outputtmp << "#             mud2rocktype: " << options["mud2rocktype"] << endl;
       outputtmp << "#           jFunctionCurve: " << options["jFunctionCurve"] << endl;
       outputtmp << "#           surfaceTension: " << options["surfaceTension"] << " dynes/cm" << endl;
       outputtmp << "#                 minPerm: " << options["minPerm"] << endl;
       outputtmp << "#                 minPoro: " << options["minPoro"] << endl;

       if (doInterpolate) { 
           outputtmp << "#              interpolate: " << options["interpolate"] << " points" << endl; 
       }
       outputtmp << "#" << endl;
       outputtmp << "# Archie options (global):" << endl;
       outputtmp << "#           lithologyCoeff: " << options["lithologyCoeff"] << endl;
       outputtmp << "#      cementationExponent: " << options["cementationExponent"] << endl;
       outputtmp << "#       saturationExponent: " << options["saturationExponent"] << endl;
       outputtmp << "# Archie options pr. rocktype:" << endl;
       for (size_t i = 1; i <= (size_t)maxSatnum; ++i) {
           stringstream rocktypestring;
           rocktypestring << i; // integer to string conversion.
           double cemoptionvalue = atof(options["rock" + rocktypestring.str() + "cemexp"].c_str());
           double satoptionvalue = atof(options["rock" + rocktypestring.str() + "satexp"].c_str());
           if (cemoptionvalue > 0) {
               outputtmp << "# CementationExp, rocktype " << i << ": " << cemoptionvalue << endl;
           }
           if (satoptionvalue > 0) {
               outputtmp << "# SaturationExp, rocktype  " << i << ": " << satoptionvalue << endl;
           }
       }
       outputtmp << "# " << endl;
       if (doInterpolate) { 
           outputtmp << "# NB: Data points shown are interpolated." << endl; 
       }
       outputtmp << "######################################################################" << endl;
       if (isFixed) {
           outputtmp << "#       Pc (Pa)       Sw            Rxx           Ryy           Rzz" << endl;
       }
       else if (isLinear) {
           outputtmp << "#       Pc (Pa)       Sw            Rxx           Ryy           Rzz           Ryz           Rxz           Rxy           Rzy           Rzx           Ryx " << endl;             
       }
       else if (isPeriodic) {
           outputtmp << "#       Pc (Pa)       Sw            Rxx           Ryy           Rzz           Ryz           Rxz           Rxy           Rzy           Rzx           Ryx " << endl;             
       }
       

       vector<double> Pvalues = pressurePoints; // WaterSaturation.get_xVector();
       vector<double> Satvalues = WaterSaturation; //.get_fVector();

       
       /* Rearrange the UpscaledConductivity array so that voigtIdx is the first index
          (to facilitate interpolation in the then last variable) 

          (results from this is only to be trusted on the master node)
       */
       
       vector<vector <double> > ResDirValues; // voigtIdx is first index.
       for (int voigtIdx=0; voigtIdx < tensorElementCount; ++voigtIdx) {
           vector<double> tmp;
           ResDirValues.push_back(tmp);
       }
         // Loop over all pressure points 
       for (int idx=0; idx < points; ++idx) {
           Matrix condTensor(zeroMatrix);
           for (int voigtIdx = 0; voigtIdx < tensorElementCount; ++voigtIdx) {
               setVoigtValue(condTensor, voigtIdx, UpscaledConductivity[idx][voigtIdx]);
           }
           //cout << phasePermTensor << endl;
           for (int voigtIdx = 0; voigtIdx < tensorElementCount; ++voigtIdx) {
               ResDirValues[voigtIdx].push_back(getVoigtValue(condTensor, voigtIdx));
           }
           //cout << relPermTensor << endl;
       }
       /*   for (int voigtIdx = 0; voigtIdx < tensorElementCount; ++voigtIdx) {
            ResDirValues.push_back(Res[voigtIdx].get_fVector());
            }*/
       
       /* If user wants interpolated output, do monotone cubic interpolation
          by modifying the data vectors that are to be printed */
       if (doInterpolate) {
           // Find min and max for saturation values
           double xmin = DBL_MAX;
           double xmax = -DBL_MAX;
           for (size_t i = 0; i < Satvalues.size(); ++i) {
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
           // Now capillary pressure and computed conductivity-values must be viewed as functions
           // of saturation, and then interpolated on the uniform saturation grid.
           
           // Now overwrite existing Pvalues and ResDirValues-data with interpolated data:
           MonotCubicInterpolator PvaluesVsSaturation(Satvalues, Pvalues);
           Pvalues.clear();
           for (int i = 0; i < interpolationPoints; ++i) {
               Pvalues.push_back(PvaluesVsSaturation.evaluate(SatvaluesInterp[i]));
           }
           for (int voigtIdx = 0; voigtIdx < tensorElementCount; ++voigtIdx) {
               MonotCubicInterpolator ResDirVsSaturation(Satvalues, ResDirValues[voigtIdx]);
               ResDirValues[voigtIdx].clear();
               for (int i=0; i < interpolationPoints; ++i) {
                   ResDirValues[voigtIdx].push_back(ResDirVsSaturation.evaluate(SatvaluesInterp[i]));
               }
           }
           // Now also overwrite Satvalues
           Satvalues.clear();
           Satvalues = SatvaluesInterp;
           
       }
       

       /* Output computed resistivity data */
       for (size_t i=0; i < Satvalues.size(); ++i) {
           // Note: The Interpolator-object's values contain the log10 of the real values.
           outputtmp << std::showpoint << std::setw(14) << Pvalues[i];
           outputtmp << std::showpoint << std::setw(14) << Satvalues[i];   
           
           for (int voigtIdx = 0; voigtIdx < tensorElementCount; ++voigtIdx) {
               outputtmp << showpoint << setw(14) << ResDirValues[voigtIdx][i];
           }
           outputtmp << endl;
       }
       
       cout << outputtmp.str();
       
       /* Possibly write to output file */
       if (options["output"] != "") {
           cout << "Writing results to " << options["output"] << endl;
           ofstream outfile;
           outfile.open(options["output"].c_str(), ios::out | ios::trunc);
           outfile << outputtmp.str();
           outfile.close();      
       }
   }

    return 0;
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}

