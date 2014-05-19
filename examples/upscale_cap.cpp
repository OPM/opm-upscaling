/*
  Copyright 2008-2011 Statoil ASA.

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

   @file upscale_cap.C
   @brief Upscales capillary pressure function, water saturation, and water saturation pr. rocktype.

   Description: 
 
   Reads in a lithofacies geometry in Eclipse format, reads in J(S_w)
   permeability and porosity, and upscale the capillary pressure function
   and water saturation, both for the whole geometry, and on a per rocktype basis

   Assumption for capillary pressure function
     - Capillary equilibrium, p_c is spatially invariant.
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
   3: Read and J-function for each stone-type.
   4: Find minimum and maximum capillary pressure from the 
      J-functions in each cell.
   5: Upscale water saturation as a function of capillary pressure
   6: Print output to screen and optionally to file.
 
    @author: Håvard Berland <havb <at> statoil.com>
 */
#include <config.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <cfloat> // for DB_MAX/DBL_MIN
#include <map>
#include <sys/utsname.h>

#include <opm/core/utility/MonotCubicInterpolator.hpp>
#include <opm/upscaling/SinglePhaseUpscaler.hpp>
#include <opm/upscaling/ParserAdditions.hpp>

#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
#include <dune/common/parallel/mpihelper.hh>
#else
#include <dune/common/mpihelper.hh>
#endif

using namespace Opm;
using namespace std;

/**
   Explains how to use the file. Shows possible option parameters.

 */
void usage()
{
    cout << "Usage: upscale_cap <options> <eclipsefile> stoneA.txt stoneB.txt ..." << endl << 
        "where the options are:" << endl <<
        "  -points <integer>            -- Number of saturation points to report for." << endl <<
        "                                  Uniformly distributed within saturation endpoints." << endl <<
        "                                  Default 50." << endl << 
        "  -jFunctionCurve <integer>    -- the column number in the stone-files that" << endl << 
        "                                  represent the Leverett J-function. Default 4." << endl <<
        "  -output <string>             -- filename for where to write upscaled values." << endl <<
        "                                  If not supplied, output will only go to " << endl <<
        "                                  the terminal (standard out)." << endl <<
        "  -surfaceTension <float>      -- Surface tension to use in J-function/Pc conversion." << endl << 
        "                                  Default 11 dynes/cm (oil-water systems). In absence of" << endl <<  
        "                                  a correct value, the surface tension for gas-oil systems " << endl << 
        "                                  could be 22.5 dynes/cm." << endl << 
        "  -maxPermContrast <float>     -- maximal permeability contrast in model." << endl <<
        "                                  Default 10^7" << endl <<
        "  -minPerm <float>             -- Minimum floating point value allowed for" << endl << 
        "                                  phase permeability in computations. If set to zero," << endl << 
        "                                  some models can end up singular for permeability" << endl <<
        "                                  upscaling. Default 10^-12" << endl << 
        "" << endl <<
        "If only one stone-file is supplied, it is used for all stone-types defined" << endl <<
        "in the geometry. If more than one, it corresponds to the SATNUM-values." << endl;
    // "minPoro" intentionally left undocumented

}


void usageandexit() {
    usage();
    exit(1);
}

int main(int varnum, char** vararg)
try
{ 

   /******************************************************************************
    * Step 1:
    * Process command line options
    */

    Dune::MPIHelper::instance(varnum, vararg);

   if (varnum == 1) { /* If no arguments supplied ("upscale_cap" is the first ('zero')  "argument") */
      usage();
      exit(1);
   }

   /*
     Populate options-map with default values
   */
   map<string,string> options;
   options.insert(make_pair("points",             "50"   )); // Number of saturation points (uniformly distributed within saturation endpoints)
   options.insert(make_pair("jFunctionCurve",     "4")); // Which column in the rock type file is the J-function curve
   options.insert(make_pair("output",             "")); // If this is set, output goes to screen and to this file. 
   options.insert(make_pair("outputprecision",    "8")); // number of decimals to print
   options.insert(make_pair("surfaceTension",     "11")); // Surface tension given in dynes/cm 
   options.insert(make_pair("maxPermContrast",    "1e7")); // maximum allowed contrast in each single-phase computation
   options.insert(make_pair("minPerm",            "1e-12")); // absoluted minimum allowed minimal cell permeability
   options.insert(make_pair("minPoro",            "0.0001")); // this limit is necessary for pcmin/max computation
   options.insert(make_pair("linsolver_tolerance", "1e-12"));  // used for swir/swmax check in upscale_cap

   // Conversion factor, multiply mD numbers with this to get m² numbers
   const double milliDarcyToSqMetre = 9.869233e-16;
   // Reference: http://www.spe.org/spe-site/spe/spe/papers/authors/Metric_Standard.pdf


   /* Check first if there is anything on the command line to look for */
   if (varnum == 1) {
      cout << "Error: No Eclipsefile or stonefiles found on command line." << endl;
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
               cout << "Parsed command line option: " << searchfor << " := " << vararg[argidx+1] << endl;
               argeclindex = argidx + 2;
           }
           else {
               cout << "Option -" << searchfor << " unrecognized." << endl;
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
     
   // argeclindex should now point to the eclipse file
   static char* ECLIPSEFILENAME(vararg[argeclindex]);
   argeclindex += 1; // argeclindex jumps to next input argument, now it points to the stone files.
   
   // argeclindex now points to the first J-function. This index is not
   // to be touched now.
   static int JFindex = argeclindex;
   

   /* Check if at least one J-function is supplied on command line */
   if (varnum <= JFindex) {
       cerr << "Error: No J-functions found on command line." << endl;
       usageandexit();
   }
    
   
   /***********************************************************************
    * Step 2:
    * Load geometry and data from Eclipse file
    */
   
   // Read data from the Eclipse file and 
   // populate our vectors with data from the file

   const double emptycellvolumecutoff = 1e-10;
      
   // Test if filename exists and is readable
   ifstream eclipsefile(ECLIPSEFILENAME, ios::in);
   if (eclipsefile.fail()) {
       cerr << "Error: Filename " << ECLIPSEFILENAME << " not found or not readable." << endl;
       usageandexit();
   }
   eclipsefile.close(); 

   cout << "Parsing Eclipse file <" << ECLIPSEFILENAME << "> ... " << endl;
   Opm::ParserPtr parser(new Opm::Parser());
   Opm::addNonStandardUpscalingKeywords(parser);
   Opm::DeckConstPtr deck(parser->parseFile(ECLIPSEFILENAME));
   
   // Check that we have the information we need from the eclipse file:  
   if (! (deck->hasKeyword("SPECGRID") && deck->hasKeyword("COORD") && deck->hasKeyword("ZCORN")  
          && deck->hasKeyword("PORO") && deck->hasKeyword("PERMX") && deck->hasKeyword("SATNUM"))) {  
       cerr << "Error: Did not find SPECGRID, COORD and ZCORN in Eclipse file " << ECLIPSEFILENAME << endl;  
       usageandexit();  
   }  
   
   vector<int>   satnums = deck->getKeyword("SATNUM")->getIntData();  
   vector<double>  poros = deck->getKeyword("PORO")->getRawDoubleData();
   vector<double> permxs = deck->getKeyword("PERMX")->getRawDoubleData();
   vector<int>  griddims(3);
   Opm::DeckRecordConstPtr specgridRecord(deck->getKeyword("SPECGRID")->getRecord(0));
   griddims[0] = specgridRecord->getItem("NX")->getInt(0);
   griddims[1] = specgridRecord->getItem("NY")->getInt(0);
   griddims[2] = specgridRecord->getItem("NZ")->getInt(0);

   unsigned int maxSatnum = 0;
   const double maxPermContrast = atof(options["maxPermContrast"].c_str());
   const double minPerm = atof(options["minPerm"].c_str());
   const double minPoro = atof(options["minPoro"].c_str());

   /* Sanity check/fix on input for each cell:
      - Check that SATNUM are set sensibly, that is => 0 and < 1000, error if not.
      - Check that porosity is between 0 and 1, error if not.
        Set to minPoro if zero or less than minPoro (due to pcmin/max computation)
      - Check that permeability is zero or positive. Error if negative. 
        Set to minPerm if zero or less than minPerm.
      - Check maximum number of SATNUM values (can be number of rock types present)
   */
   for (unsigned int i = 0; i < satnums.size(); ++i) {
       if (satnums[i] < 0 || satnums[i] > 1000) { 
           cerr << "satnums[" << i << "] = " << satnums[i] << ", not sane, quitting." << endl;
           usageandexit();
       }
       if (satnums[i] > (int)maxSatnum) {
           maxSatnum = satnums[i];
       }
       if ((poros[i] >= 0) && (poros[i] < minPoro)) { // Truncate porosity from below
           poros[i] = minPoro;
       }
       if (poros[i] < 0 || poros[i] > 1) {
           cerr << "poros[" << i <<"] = " << poros[i] << ", not sane, quitting." << endl;
           usageandexit();
       }
       if ((permxs[i] >= 0) && (permxs[i] < minPerm)) { // Truncate permeability from below
           permxs[i] = minPerm;
       }
       if (permxs[i] < 0) {
           cerr << "permx[" << i <<"] = " << permxs[i] << ", not sane, quitting." << endl;
           usageandexit();
       }
       // Explicitly handle "no rock" cells, set them to minimum perm and zero porosity.
       if (satnums[i] == 0) {
           permxs[i] = minPerm;
           poros[i] = 0; // zero poro is fine for these cells, as they are not 
                         // used in pcmin/max computation.
       }
   }  


   /***************************************************************************
    * Step 3:
    * Load relperm- and J-function-curves for the stone types.
    * We read columns from text-files, syntax allowed is determined 
    * by MonotCubicInterpolator which actually opens and parses the 
    * text files.
    */

   // Number of stone-types is max(satnums):
   
   // If there is only one J-function supplied on the command line,
   // use that for all stone types.

   int stone_types = int(*(max_element(satnums.begin(), satnums.end())));
   std::vector<MonotCubicInterpolator> InvJfunctions; // Holds the inverse of the loaded J-functions.
   
   std::vector<string> JfunctionNames; // Placeholder for the names of the loaded J-functions.

   // Input for surfaceTension is dynes/cm
   // SI units are Joules/square metre
   const double surfaceTension     = atof(options["surfaceTension"].c_str()) * 1e-3; // multiply with 10^-3 to obtain SI units 
   const int jFunctionCurve        = atoi(options["jFunctionCurve"].c_str());
   const int interpolationPoints   = atoi(options["points"].c_str());
   const int outputprecision       = atoi(options["outputprecision"].c_str());

   // Handle two command line input formats, either one J-function for all stone types
   // or one each. If there is only one stone type, both code blocks below are equivalent.
   
   if (varnum == JFindex + stone_types) {
      for (int i=0 ; i < stone_types; ++i) {
         const char* ROCKFILENAME = vararg[JFindex+i];
         // Check if rock file exists and is readable:
         ifstream rockfile(ROCKFILENAME, ios::in);
         if (rockfile.fail()) {
            cerr << "Error: Filename " << ROCKFILENAME << " not found or not readable." << endl;
            usageandexit();
         }
         rockfile.close(); 
         MonotCubicInterpolator Jtmp;
         try {
             Jtmp = MonotCubicInterpolator(ROCKFILENAME, 1, jFunctionCurve); 
         }
         catch (const char * errormessage) {
             cerr << "Error: " << errormessage << endl;
             cerr << "Check filename and -jFunctionCurve" << endl;
             usageandexit();
         }
         // Invert J-function, now we get saturation as a function of pressure:
         if (Jtmp.isStrictlyMonotone()) {
            InvJfunctions.push_back(MonotCubicInterpolator(Jtmp.get_fVector(), Jtmp.get_xVector()));
            JfunctionNames.push_back(vararg[JFindex + i]);
         }
         else {
             cerr << "Error: Jfunction " << i+1 << " in rock file " << ROCKFILENAME << " was not invertible." << endl;
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
            cerr << "Error: Filename " << ROCKFILENAME << " not found or not readable." << endl;
            usageandexit();
         }
         rockfile.close(); 
         MonotCubicInterpolator Jtmp;
         try {
             Jtmp = MonotCubicInterpolator(ROCKFILENAME, 1, jFunctionCurve);
         }
         catch(const char * errormessage) {
             cerr << "Error: " << errormessage << endl;
             cerr << "Check filename and -jFunctionCurve" << endl;
             usageandexit();
         }
         // Invert J-function, now we get saturation as a function of pressure:
         if (Jtmp.isStrictlyMonotone()) {
            InvJfunctions.push_back(MonotCubicInterpolator(Jtmp.get_fVector(), Jtmp.get_xVector()));
            JfunctionNames.push_back(vararg[JFindex]);
         }
         else {
            cerr << "Error: Jfunction " << i+1 << " in rock file " << ROCKFILENAME << " was not invertible." << endl;
            usageandexit();
         }
      }
   }
   else {
      cerr << "Error:  Wrong number of stone-functions provided. " << endl;
      usageandexit();
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

   if (maxPermContrast == 0) {
       cout << "Illegal contrast value" << endl;
       usageandexit();
   }
   

   // Construct an object for single-phase upscaling, since we need to get some
   // information from the grid.
   SinglePhaseUpscaler upscaler;
   upscaler.init(deck, SinglePhaseUpscaler::Fixed,
                 Opm::unit::convert::from(minPerm, Opm::prefix::milli*Opm::unit::darcy),
                 0.0, 1e-8, 0, 1, false);  // options on this line are noops for upscale_cap

   vector<double> cellVolumes, cellPoreVolumes; 
   cellVolumes.resize(satnums.size(), 0.0);
   cellPoreVolumes.resize(satnums.size(), 0.0);

   /* Volumes/saturations pr. rocktype */
   vector<double> cellporevolume_rocktype;
   cellporevolume_rocktype.resize(maxSatnum + 1, 0.0);

   vector<double> watervolume_rocktype;
   watervolume_rocktype.resize(maxSatnum + 1, 0.0);


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
   int tesselatedCells = 0; // for counting "active" cells (Sintef interpretation of "active")
   double Pcmax = -DBL_MAX, Pcmin = DBL_MAX;
   double maxSinglePhasePerm = 0;
   double Swirvolume = 0;
   double Sworvolume = 0;

   const std::vector<int>& ecl_idx = upscaler.grid().globalCell();
   Dune::CpGrid::Codim<0>::LeafIterator c = upscaler.grid().leafbegin<0>();
   for (; c != upscaler.grid().leafend<0>(); ++c) {
       unsigned int cell_idx = ecl_idx[c->index()];
       if (satnums[cell_idx] > 0) { // Satnum zero is "no rock"

	   cellVolumes[cell_idx] = c->geometry().volume();
	   cellPoreVolumes[cell_idx] = cellVolumes[cell_idx] * poros[cell_idx];

	   
	   double Pcmincandidate = InvJfunctions[int(satnums[cell_idx])-1].getMinimumX().first
	       / sqrt(permxs[cell_idx] * milliDarcyToSqMetre/poros[cell_idx]) * surfaceTension;
	   Pcmin = min(Pcmincandidate, Pcmin);
           
	   double Pcmaxcandidate = InvJfunctions[int(satnums[cell_idx])-1].getMaximumX().first
	       / sqrt(permxs[cell_idx] * milliDarcyToSqMetre/poros[cell_idx]) * surfaceTension;
	   Pcmax = max(Pcmaxcandidate, Pcmax);
           
	   maxSinglePhasePerm = max( maxSinglePhasePerm, permxs[cell_idx]);
           
	   cellporevolume_rocktype[satnums[cell_idx]] += cellPoreVolumes[cell_idx];
	   double minSw = InvJfunctions[int(satnums[cell_idx])-1].getMinimumF().second;
	   double maxSw = InvJfunctions[int(satnums[cell_idx])-1].getMaximumF().second;
           
	   // cout << "minSwc: " << minSw << endl;
	   // cout << "maxSwc: " << maxSw << endl;
           
	   // Add irreducible water saturation volume
	   Swirvolume += minSw * cellPoreVolumes[cell_idx];
	   Sworvolume += maxSw * cellPoreVolumes[cell_idx];
           
       }
       ++tesselatedCells; // keep count.
   }

   //double minSinglePhasePerm = max(maxSinglePhasePerm/maxPermContrast, minPerm);
   
   cout << "Pcmin:    " << Pcmin << endl;
   cout << "Pcmax:    " << Pcmax << endl;

   if (Pcmin > Pcmax) {
       cerr << "ERROR: No legal capillary pressures found for this system. Exiting..." << endl;
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

   cout << "LF Pore volume:    " << poreVolume << endl;
   cout << "LF Volume:         " << volume << endl;
   cout << "Upscaled porosity: " << poreVolume/volume << endl;
   cout << "Upscaled Swir:     " << Swir << endl;
   cout << "Upscaled Swmax:    " << Swor << endl; //Swor=1-Swmax

   // Sometimes, if Swmax=1 or Swir=0 in the input tables, the upscaled 
   // values can be a little bit larger (within machine precision) and
   // the check below fails. Hence, check if these values are within the 
   // the [0 1] interval within some precision
   double linsolver_tolerance = atof(options["linsolver_tolerance"].c_str());
   if (Swor > 1.0 && Swor - linsolver_tolerance < 1.0) {
       Swor = 1.0;
   }
   if (Swir < 0.0 && Swir + linsolver_tolerance > 0.0) {
       Swir = 0.0;
   }
   if (Swir < 0.0 || Swir > 1.0 || Swor < 0.0 || Swor > 1.0) {
       cerr << "ERROR: Swir/Swor unsensible. Check your input. Exiting";
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
    * is 1/1000 of the saturation interval. Monotone cubic interpolation
    * will be used afterwards for accessing the tabulated values.
    */

   MonotCubicInterpolator WaterSaturationVsCapPressure;
   
   double largestSaturationInterval = Swor-Swir;
   
   double Ptestvalue = Pcmax;

   vector<MonotCubicInterpolator> watersaturation_rocktype;
   for (unsigned int satidx=0; satidx <= maxSatnum; ++satidx) {
       MonotCubicInterpolator tmp;
       watersaturation_rocktype.push_back(tmp);
   }
   
   while (largestSaturationInterval > (Swor-Swir)/200.0) {
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
           cerr << "ERROR: Ptestvalue was inf or nan" << endl;
           break; // Jump out of while-loop, just print out the results
           // up to now and exit the program
       }

       // Reset array to zero.
       for (unsigned int satidx = 0; satidx <= maxSatnum; ++satidx) {
           watervolume_rocktype[satidx] = 0.0;
       }
     
       double waterVolume = 0.0;
       for (unsigned int cell_idx = 0; cell_idx < satnums.size(); ++cell_idx) {
           if (cellVolumes[cell_idx] > emptycellvolumecutoff) {
               double waterSaturationCell = 0.0;
               if (satnums[cell_idx] > 0) { // handle "no rock" cells with satnum zero
                   double PtestvalueCell;
                   PtestvalueCell = Ptestvalue;
                   double Jvalue = sqrt(permxs[cell_idx] * milliDarcyToSqMetre/poros[cell_idx]) 
                       * PtestvalueCell / surfaceTension;
                   //cout << "JvalueCell: " << Jvalue << endl;
                   waterSaturationCell 
                       = InvJfunctions[int(satnums[cell_idx])-1].evaluate(Jvalue);
               }
               waterVolume += waterSaturationCell  * cellPoreVolumes[cell_idx];
               watervolume_rocktype[satnums[cell_idx]] += waterSaturationCell * cellPoreVolumes[cell_idx];
               
           }
       }
       WaterSaturationVsCapPressure.addPair(Ptestvalue, waterVolume/poreVolume);
       //       cout << "Ptestvalue " << Ptestvalue << " sat: " << waterVolume/poreVolume << endl;
       for (unsigned int satidx = 1; satidx <= maxSatnum; ++satidx) {
           // cout << "satidx "<< satidx << " " << watervolume_rocktype[satidx]/cellporevolume_rocktype[satidx] << endl;
           //cout << "watvol: " << watervolume_rocktype[satidx] << " porevol " << cellporevolume_rocktype[satidx] << endl;
           if (cellporevolume_rocktype[satidx] > 0) {
               watersaturation_rocktype[satidx].addPair(Ptestvalue, watervolume_rocktype[satidx]/cellporevolume_rocktype[satidx]);
           }
           else {
               watersaturation_rocktype[satidx].addPair(Ptestvalue, 0.0);
           }
       }

   }

   // Check if the saturation vs cap pressure curve is monotone
   // If not, it would have been a problem for upscale_relperm, but
   // it is not as critical here, so we only issue a warning
   // (upscale_relperm solves this by issung chopFlatEndpoints and possibly shrinkFlatAreas, 
   // but this is trickier to implement in this code due to watersaturation_rocktype[satidx])
   if (!WaterSaturationVsCapPressure.isStrictlyMonotone()) {
       {
           cerr << "Warning: Upscaled water saturation not strictly monotone in capillary pressure." << endl;
           cerr << "         Unphysical input data?." << endl;
       }
   }
   MonotCubicInterpolator CapPressureVsWaterSaturation(WaterSaturationVsCapPressure.get_fVector(), 
                                                       WaterSaturationVsCapPressure.get_xVector());

   
   /*********************************************************************************
    *  Step 9
    *
    * Output results to stdout and optionally to file. Note, we only output to
    * file if the '-outputWater'-option and/or '-outputOil' has been set, as this option is an
    * empty string by default.
    */
   vector<double> Pvalues = WaterSaturationVsCapPressure.get_xVector(); 
   vector<double> Satvalues = WaterSaturationVsCapPressure.get_fVector(); 
   
   vector<vector<double> > watersaturation_rocktype_values;
   vector<double> tmp;
   watersaturation_rocktype_values.push_back(tmp); // dummy zero index element
   for (unsigned int satidx=1; satidx <= maxSatnum; ++satidx) {
       watersaturation_rocktype_values.push_back(watersaturation_rocktype[satidx].get_fVector());
   }
   stringstream outputtmp;
   
   // Print a table of all computed values:
   outputtmp << "######################################################################" << endl;
   outputtmp << "# Results from upscaling capillary pressure and water saturations."<< endl;
   outputtmp << "#" << endl;
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
   outputtmp << "# " << endl;
   outputtmp << "# Options used:" << endl;
   outputtmp << "#          jFunctionCurve: " << options["jFunctionCurve"] << endl;
   outputtmp << "#                  points: " << options["points"] << endl;
   outputtmp << "#         maxPermContrast: " << options["maxPermContrast"] << endl;
   outputtmp << "#          surfaceTension: " << options["surfaceTension"] << endl;   
   outputtmp << "######################################################################" << endl;
   outputtmp << "#         Pc (Pa)         Sw              Sw1           Sw2       Sw3 etc.." << endl; 
   
   
  // If user wants interpolated output, do monotone cubic interpolation
   // by modifying the data vectors that are to be printed
   if (interpolationPoints > 0) {
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

       // Now overwrite existing Pvalues and saturation-data with interpolated data:
       MonotCubicInterpolator PvaluesVsSaturation(Satvalues, Pvalues);
       Pvalues.clear();
       for (int i = 0; i < interpolationPoints; ++i) {
           Pvalues.push_back(PvaluesVsSaturation.evaluate(SatvaluesInterp[i]));
       }
       for (unsigned int satidx = 1; satidx <= maxSatnum; ++satidx) {
           MonotCubicInterpolator WaterSaturationRocktypeVsSaturation(Satvalues, watersaturation_rocktype_values[satidx]);
           watersaturation_rocktype_values[satidx].clear();
           for (int i=0; i < interpolationPoints; ++i) {
               watersaturation_rocktype_values[satidx].push_back(WaterSaturationRocktypeVsSaturation.evaluate(SatvaluesInterp[i]));
           }
       }
       // Now also overwrite Satvalues
       Satvalues.clear();
       Satvalues = SatvaluesInterp;
   }


   const int fieldwidth = outputprecision + 8;
   for (unsigned int i=0; i < Satvalues.size(); ++i) {
       outputtmp << showpoint << setw(fieldwidth) << setprecision(outputprecision) << Pvalues[i]; 
       outputtmp << showpoint << setw(fieldwidth) << setprecision(outputprecision) << Satvalues[i]; 
       for (unsigned int satidx = 1; satidx <= maxSatnum; ++satidx) { 
           outputtmp << showpoint << setw(fieldwidth) << setprecision(outputprecision) 
                     << watersaturation_rocktype_values[satidx][i]; 
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


   return 0;
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}

