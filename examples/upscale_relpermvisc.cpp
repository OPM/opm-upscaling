
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
  @file upscale_relpermvisc.C
  @brief Upscales relative permeability as a funtion of watersaturation in the viscous limit.
  
  Description: 
 
  Upscaling of relative permeability (two-phase) using steady-state
  in the viscous limit, on Eclipse corner-point geometries in shoe-box format
  
  Reads in a lithofacies geometry in Eclipse format, reads in
  viscosities for oil and water from command line, relpermcurve(S_w)
  for each stone type from file, and calculates upscaled (three
  directions) relative permeability curves as a function of Sw.

  Fixed, linear and periodic boundary conditions are supported,
  yielding 3, 9 or 9 values respectively (but only 6 significant for
  periodic) for relative permeability at each saturation point.
  
  The relative permeability computation is based on 
    - Steady-state, viscous limit. v_w/v_o = constant, === lambda_w/lamda_o = constant
    - No gravitational effects.
    - No capillary effects.
 
  Steps in the code:
 
  1: Process command line options.
  2: Read Eclipse file 
  3: Read relperm-function for each stone-type.
  4: Tesselate the grid (Sintef code)
  5: Generate simple statistics by looping over the cells.
  6: Upscale fractional flow ratio vs. water saturation
  6: Find upper and lower bounds for fractional flow ratio v_w/v_o 
  7: Upscale single phase permeability (in order to make relative perm later)
  8: For uniformly distributed water saturation values between upscaled Swir and Swor, 
      a: Find fractional flow ratio corresponding to wanted upscaled water saturation
      b: Model water saturation in each cell given fractional flow ratio
      c: Model phase permeability in each cell given cell water saturation and inputted
         relative permeability curves.
      d: Upscale phase permeability
   9: Repeat step 8 for the oil phase. 
  10: Print output to screen and optionally to files. 
 
  The relperm-functions must be defined in text files, with water saturation in 
  the first column, and then relative permeability for water in column 2, and relative 
  permeability for oil in column 3. Lines starting with -- or # are ignored.
 
  Units for (dynamic) viscosity are assumed to be in Pascal seconds. (Pa s)
  1000 Pa s = 1 centiPoise.
 */
#include <config.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <cfloat>  // FOR DBL_MAX/DBL_MIN
#include <map>
#include <sys/utsname.h>

#ifdef USEMPI
#include <mpi.h>
#endif

#include <opm/core/utility/MonotCubicInterpolator.hpp>
#include <opm/upscaling/SinglePhaseUpscaler.hpp>

using namespace Opm;
using namespace std;

/**
   The usage() function displays a message with syntax and available 
   options for the program, and is typically called whenever the code 
   encounters any errors in the options or in the command line syntax 
   (missing filenames etc.) 
 */
void usage()
{
    cout << "Usage: upscale_relpermvisc <options> <eclipsefile> rock1.txt rock2.txt... (isotropic case)" << endl <<
        "       upscale_relpermvisc <options> <eclipsefile> rock1_water.txt rock1_oil.txt rock2_water.txt rock2_oil.txt ... (anisotropic case)" << endl <<
        " where the options are:" << endl <<
        "-bc <string>                 -- which boundary conditions to use." << endl << 
        "                                Possible values are f (fixed), " << endl << 
        "                                l (linear) and p (periodic). Default f." << endl <<
        "-points <integer>            -- Number of saturation points to upscale for." << endl <<
        "                                Uniformly distributed within saturation endpoints." << endl <<
        "                                Default 30." << endl <<
        "-waterViscosity <float >     -- viscosity of water given in Pascal seconds" << endl <<
        "                                (1000 cP = 1 Pa s)" << endl <<
        "-oilViscosity <float>        -- ditto" << endl <<
        "-waterCurve <integer>        -- the column number in the stone files that represent" << endl <<
        "                                relative permeability for water. Default 2" << endl <<
        "-oilCurve <integer>          -- the column number in the stone files that represents " << endl <<
        "                                relative permeability for oil. Default 3." << endl << 
        "-outputWater <string>        -- filename for where to write upscaled values for" << endl <<
        "                                water relperm. If not supplied, output will only" << endl <<
        "                                go to the terminal (standard out)." << endl <<
        "-outputOil <string>          -- ditto" << endl <<
        "-interpolate <integer>       -- If supplied and > 1, the output data points will be" << endl << 
 	"                                interpolated using monotone cubic interpolation" << endl << 
 	"                                on a uniform grid with the specified number of" << endl << 
 	"                                points. Suggested value: 1000." << endl <<         "" << endl <<
        "-minPerm <float>             -- Minimum floating point value allowed for" << endl <<
        "                                phase permeability in computations. If set to zero," << endl <<
        "                                some models can end up singular. Default 1e-12" << endl <<
        "If only one stone-file (only relperm-values are used) is supplied, it" << endl <<
        "is used for all stone-types defined in the geometry. If more than one," << endl <<
        "it corresponds to the SATNUM-values." << endl;
    cout << "$Rev: 491 $" << endl;
    // Intentionally undocumented features:
    //  -outputprecision
    //  -maxPermContrast
    //  -minPoro
    //  -saturationThreshold
                    
}

void usageandexit() {
#ifdef USEMPI
    MPI_Finalize();
#endif
    usage();
    exit(1);
}

// Assumes that permtensor_t use C ordering.
double getVoigtValue(const SinglePhaseUpscaler::permtensor_t& K, int voigt_idx)
{
    ASSERT(K.numRows() == 3 && K.numCols() == 3);
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
    ASSERT(K.numRows() == 3 && K.numCols() == 3);
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
{ 
    // Variables used for timing/profiling:
    clock_t start, finish;
    double timeused = 0, timeused_tesselation = 0; //, timeused_upscale_acc_water = 0, timeused_upscale_acc_oil = 0; 
    double timeused_upscale_wallclock = 0.0;
    
    /******************************************************************************
     * Step 1:
     * Process command line options
     */
    
    int mpi_rank = 0;
#ifdef USEMPI
    int mpi_nodecount = 1;
    MPI_Init(&varnum, &vararg);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_nodecount);
#endif 
    bool isMaster = (mpi_rank == 0);
    if (varnum == 1) { /* If no arguments supplied ("upscale_relpermvisc" is the first "argument") */
        usage();
        exit(1);
    }
    
    /*
      Populate options-map with default values
    */
    map<string,string> options;
    options.insert(make_pair("bc",                 "f"     )); // Boundary conditions, default is fixed. 
    options.insert(make_pair("points",             "30"   )); // Number of saturation points (uniformly distributed within saturation endpoints)
    options.insert(make_pair("outputWater",        ""));       // If set, where to write upscaled water permeabilities
    options.insert(make_pair("outputOil",          ""));       // If set, where to write upscaled oil permeabilites
    options.insert(make_pair("waterCurve",         "2"));      // column in stonefiles that represents rel.perm. for water
    options.insert(make_pair("oilCurve",           "3"));      // column in stonefiles that represents rel.perm. for oil
    options.insert(make_pair("waterViscosity",     "0.001"));  // viscosity of water in Pascal seconds.
    options.insert(make_pair("oilViscosity",       "0.1"));    // viscosity of oil in Pascal seconds.
    options.insert(make_pair("interpolate",        "0"));    // default is not to interpolate 
    //options.insert(make_pair("outputprecision",    "8")); // number of decimals to print
    options.insert(make_pair("minPerm",            "1e-12")); // perm values below this will be truncated upwards
    options.insert(make_pair("maxPerm",            "100000")); // perm values below this will be truncated upwards
    options.insert(make_pair("minPoro",            "0.0001")); // poro values below this will be truncated upwards
    options.insert(make_pair("maxPermContrast",    "1e7")); // maximum allowed contrast in each single-phase computation
    options.insert(make_pair("saturationThreshold", "0.00001")); // accuracy threshold for saturation, we ignore Pc values
    // that give so small contributions to Sw near endpoints.
    options.insert(make_pair("linsolver_tolerance", "1e-12"));  // residual tolerance for linear solver
    options.insert(make_pair("linsolver_verbosity", "0"));     // verbosity level for linear solver
    options.insert(make_pair("linsolver_type",      "1"));     // type of linear solver: 0 = ILU/BiCGStab, 1 = AMG/CG
    
    /* Check first if there is anything on the command line to look for */
    if (varnum == 1) {
        if (isMaster) cerr << "Error: No Eclipsefile or stonefiles found on command line." << endl;
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
    
    // argeclindex should now point to the eclipse file
    static char* ECLIPSEFILENAME(vararg[argeclindex]);
    argeclindex += 1; // argeclindex jumps to next input argument, now it points to the stone files.
    
    // Boolean set to true of input permeability in eclipse-file has diagonal anisotropy.
    // (full-tensor anisotropy will be ignored)
    bool anisotropic_input = false;
    
    // argcindex now points to the first J-function. This index is not
    // to be touched now.
    static int rockfileindex = argeclindex;
    

    /* Check if at least one J-function is supplied on command line */
    if (varnum <= rockfileindex) {
        if (isMaster) cerr << "Error: No J-functions found on command line." << endl;
        usageandexit();
    }
    
    /* Check validity of boundary conditions chosen, and make booleans 
       for boundary conditions, this allows more readable code later. */
    bool isFixed, isLinear, isPeriodic; 
    SinglePhaseUpscaler::BoundaryConditionType boundaryCondition; 
    int tensorElementCount; // Number of independent elements in resulting tensor. 
    if (options["bc"].substr(0,1) == "f") { 
        isFixed = true; isLinear = false; isPeriodic = false; 
        boundaryCondition = SinglePhaseUpscaler::Fixed; 
        tensorElementCount = 3; // Diagonal 
    } 
    else if (options["bc"].substr(0,1) == "l") { 
        isLinear = true; isFixed = false; isPeriodic = false; 
        boundaryCondition = SinglePhaseUpscaler::Linear; 
        tensorElementCount = 9; // Full-tensor 
    } 
    else if (options["bc"].substr(0,1) == "p") { 
        isPeriodic = true; isLinear = false; isFixed = false; 
        boundaryCondition = SinglePhaseUpscaler::Periodic; 
        tensorElementCount = 9; // Symmetric. 
    } 
    else { 
        if (isMaster) cout << "Invalid boundary condition. Only one of the letters f, l or p are allowed." << endl; 
        usageandexit();
    } 
 
    // Index values for the phases. We are going to use the same code
    // twice to avoid source code duplication in this file, first for
    // water, and then for oil. To accomplish this, we have some data
    // in a vector where the first element is the waterphase and the
    // second is the oilphase, defined by these constants.
    const int waterPhaseIndex  = 0; 
    const int oilPhaseIndex    = 1; 
  
    const double points = atof(options["points"].c_str());
    
    const int waterCurveColumn = atoi(options["waterCurve"].c_str());
    const int oilCurveColumn   = atoi(options["oilCurve"].c_str());
    //const int outputprecision  = atoi(options["outputprecision"].c_str());
 
 
    // Vector of column index for input rel perm curves 
    vector<int> relPermCurves; 
    relPermCurves.push_back(waterCurveColumn); 
    relPermCurves.push_back(oilCurveColumn); 
    // Now this vector can be accessed with 
    // relPermCurves[waterPhaseIndex] etc. 
    
    // Vector of viscosities     
    vector<double> viscosities; 
    viscosities.push_back(atof(options["waterViscosity"].c_str())); 
    viscosities.push_back(atof(options["oilViscosity"].c_str()));    
    // Now this vector can be accessed with  
    // viscosities[waterPhaseIndex] to give the viscosity for water. 
    
 
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
        if (isMaster) cerr << "Error: Filename " << ECLIPSEFILENAME << " not found or not readable." << endl;
        usageandexit();
    }
    eclipsefile.close(); 
 
    if (isMaster) cout << "Parsing Eclipse file <" << ECLIPSEFILENAME << "> ... ";
    flush(cout);   start = clock();
 
    Opm::EclipseGridParser eclParser(ECLIPSEFILENAME,false);
    finish = clock();   timeused = (double(finish)-double(start))/CLOCKS_PER_SEC;
    if (isMaster) cout << " (" << timeused <<" secs)" << endl;
  
    // Check that we have the information we need from the eclipse file: 
    if (! (eclParser.hasField("SPECGRID") && eclParser.hasField("COORD") && eclParser.hasField("ZCORN") 
           && eclParser.hasField("PORO") && eclParser.hasField("PERMX") && eclParser.hasField("SATNUM"))) { 
        cerr << "Error: Did not find SPECGRID, COORD and ZCORN in Eclipse file " << ECLIPSEFILENAME << endl; 
        usage(); 
        exit(1); 
    } 
    
    vector<double> poros  = eclParser.getFloatingPointValue("PORO"); 
    vector<double> permxs = eclParser.getFloatingPointValue("PERMX"); 
 
    // Load anisotropic (only diagonal supported) input if present in grid
    vector<double> permys, permzs;
    
    if (eclParser.hasField("PERMY") && eclParser.hasField("PERMZ")) {
        anisotropic_input = true;
        permys = eclParser.getFloatingPointValue("PERMY");
        permzs = eclParser.getFloatingPointValue("PERMZ");
        if (isMaster) cout << "Info: PERMY and PERMZ present, going into anisotropic input mode, no J-functions\n"; 
        if (isMaster) cout << "      Options -relPermCurve and -jFunctionCurve is meaningless.\n"; 
    } 
    
    
    /* Initialize a default satnums-vector with only "ones" (meaning only one rocktype) */ 
    vector<int> satnums(poros.size(), 1); 
    
    if (eclParser.hasField("SATNUM")) { 
        satnums = eclParser.getIntegerValue("SATNUM"); 
    } 
    else if (eclParser.hasField("ROCKTYPE")) { 
        satnums = eclParser.getIntegerValue("ROCKTYPE"); 
    } 
    else { 
        if (isMaster) cout << "Warning: SATNUM or ROCKTYPE not found in input file, assuming only one rocktype" << endl; 
    } 
 

    /* Sanity check/fix on input for each cell:
       - Check that SATNUM are set sensibly, that is => 0 and < 1000, error if not.
       - Check that porosity is between 0 and 1, error if not.
         Set to minPoro if zero or less than minPoro (due to pcmin/max computation)
       - Check that permeability is zero or positive. Error if negative. 
         Set to minPerm if zero or less than minPerm.
       - Check maximum number of SATNUM values (can be number of rock types present)
    */
    const double maxPermContrast = atof(options["maxPermContrast"].c_str());
    const double minPerm = atof(options["minPerm"].c_str());
    const double maxPerm = atof(options["maxPerm"].c_str());
    const double minPoro = atof(options["minPoro"].c_str());   
    const double saturationThreshold = atof(options["saturationThreshold"].c_str());
    double maxPermInInputFile = 0.0;
    int cells_truncated_from_below_poro = 0;
    int cells_truncated_from_below_permx = 0;
    int cells_truncated_from_above_permx = 0;
    if (minPerm == 0) {
        if (isMaster) cout << "Warning: minPerm set to zero." << endl;
    }
    if (maxPermContrast <= 0) {
        if (isMaster) cout << "Illegal maxPermContrast value" << endl;
        usageandexit();
    }
 
    int maxSatnum = 0;
    for (unsigned int i = 0; i < satnums.size(); ++i) {
        if (satnums[i] < 0 || satnums[i] > 1000) { 
            if (isMaster) cerr << "satnums[" << i << "] = " << satnums[i] << ", not sane, quitting." << endl;
            usageandexit();
        }
        if (satnums[i] > maxSatnum) {
            maxSatnum = satnums[i];
        }
        if ((poros[i] >= 0) && (poros[i] < minPoro)) {
            poros[i] = minPoro;
            ++cells_truncated_from_below_poro;
        }
        if (poros[i] < 0 || poros[i] > 1) {
            if (isMaster) cerr << "poros[" << i <<"] = " << poros[i] << ", not sane, quitting." << endl;
            usageandexit();
        }
        if (permxs[i] > maxPermInInputFile) {
            maxPermInInputFile = permxs[i];
        }
        if ((permxs[i] >= 0) && (permxs[i] < minPerm)) {
            permxs[i] = minPerm;
            ++cells_truncated_from_below_permx;
        }
        if (permxs[i] > maxPerm) { // Truncate permeability from above
            permxs[i] = maxPerm;
            ++cells_truncated_from_above_permx;
        }
        if (permxs[i] < 0) {
            if (isMaster) cerr << "permx[" << i <<"] = " << permxs[i] << ", not sane, quitting." << endl;
            usageandexit();
        }
        if ((permxs[i] >= 0) && (permxs[i] > maxPerm)) {
          permxs[i] = maxPerm;
        }
        if (anisotropic_input) {
            if (permys[i] < 0) {
                if (isMaster) cerr << "permy[" << i <<"] = " << permys[i] << ", not sane, quitting." << endl;
                usageandexit();
            }
            if (permzs[i] < 0) {
                if (isMaster) cerr << "permz[" << i <<"] = " << permzs[i] << ", not sane, quitting." << endl;
                usageandexit();
            }
        }
        // Explicitly handle "no rock" cells, set them to minimum perm and zero porosity.
        if (satnums[i] == 0) {
            permxs[i] = minPerm;
            poros[i] = 0; // zero poro is fine for these cells, as they are not 
            // used in pcmin/max computation.
            if (anisotropic_input) {
                permys[i] = minPerm;
                permzs[i] = minPerm;
            }
        }
    }  
    if (cells_truncated_from_below_poro > 0) {
        cout << "Cells with truncated porosity: " << cells_truncated_from_below_poro << endl;
    }
    if (cells_truncated_from_below_permx > 0) {
        cout << "Cells with permx truncated from below: " << cells_truncated_from_below_permx << endl;
    }
    if (cells_truncated_from_above_permx > 0) {
        cout << "Cells with permx truncated from above: " << cells_truncated_from_above_permx << endl;
    }

    /***************************************************************************
     * Step 3:
     * Load relperm-curves for the stone types.
     * We read columns from text-files, syntax allowed is determined 
     * by MonotCubicInterpolator which actually opens and parses the 
     * text files.
     * 
     * How do we determine which mode of operation?
     *  - If PERMY and PERMZ are present in grdecl-file, we are in the anisotropic mode
     *  
     */
 
    int stone_types = int(*(max_element(satnums.begin(), satnums.end())));
    std::vector<MonotCubicInterpolator> Krw; // Holds relperm-water-curves for each stone type
    std::vector<MonotCubicInterpolator> Kro; // Holds relperm-oil-curves for each stone type
    // If anisotropic input
    std::vector<MonotCubicInterpolator> Krwx, Krwy, Krwz, Krox, Kroy, Kroz;
    std::vector<string> rockTypeNames; // Placeholder for the names of the loaded rock perm files,
                                       // to be used in final output.
 
    // Handle two command line input formats, either stone.txt for all stone types
    // *or* one each. If there is only one stone type, both code blocks below are equivalent.
    if (! anisotropic_input) {
        if (varnum == rockfileindex + stone_types) { // one .txt for each stone type in isotropic case
            for (int i=0 ; i < stone_types; ++i) { 
                const char* ROCKFILENAME = vararg[rockfileindex+i];
                // Check if rock file exists and is readable:
                ifstream rockfile(ROCKFILENAME, ios::in);
                if (rockfile.fail()) {
                    if (isMaster) cerr << "Error: Filename " << ROCKFILENAME << " not found or not readable." << endl;
                    usageandexit();
                }
                rockfile.close(); 
                MonotCubicInterpolator Krwtmp;
                try {
                    Krwtmp = MonotCubicInterpolator(ROCKFILENAME, 1, relPermCurves[waterPhaseIndex]); 
                }
                catch (const char * errormessage) {
                    if (isMaster) cerr << "Error: " << errormessage << " Filename: " << ROCKFILENAME << endl;
                    if (isMaster) cerr << "Check filename and -waterCurve" << endl;
                    usageandexit();
                }
                if (!Krwtmp.isStrictlyMonotone()) {
                    if (isMaster) cerr << "Error: Data in column " << relPermCurves[waterPhaseIndex] << " of file " << ROCKFILENAME << endl <<
                        "       was not strictly monotone. Exiting." << endl;
                    usageandexit();
                }
                MonotCubicInterpolator Krotmp;
                try {
                    Krotmp = MonotCubicInterpolator(ROCKFILENAME, 1, relPermCurves[oilPhaseIndex]); 
                }
                catch (const char * errormessage) {
                    if (isMaster) cerr << "Error: " << errormessage << " Filename: " << ROCKFILENAME <<  endl;
                    if (isMaster) cerr << "Check filename and -oilCurve" << endl;
                    usageandexit();
                }
                
                if (!Krotmp.isStrictlyMonotone()) {
                    if (isMaster) cerr << "Error: Data in column " << relPermCurves[oilPhaseIndex] << " of file " << ROCKFILENAME << endl <<
                        "       was not strictly monotone. Exiting." << endl;
                    usageandexit();
                }
                
                if ((Krwtmp.isStrictlyIncreasing() && Krotmp.isStrictlyIncreasing()) || (Krwtmp.isStrictlyDecreasing() && Krotmp.isStrictlyDecreasing())) {
                    if (isMaster) cerr << "Error: Input relperm curves are both increasing or decreasing in file " << ROCKFILENAME << endl;
                    usageandexit();
                }
                
                Krw.push_back(Krwtmp); 
                Kro.push_back(Krotmp); 
                
                if (isMaster) cout << "Loaded rock file: " << ROCKFILENAME  
                                   << ", for stone type " << i+1 << endl; 
                
                rockTypeNames.push_back(ROCKFILENAME); 
            }
        } 
        else if (varnum == rockfileindex + 1) { // one .txt for all stone types 
            const char* ROCKFILENAME = vararg[rockfileindex];
            // Check if rock files exists and is readable:
            ifstream rockfile(ROCKFILENAME, ios::in);
            if (rockfile.fail()) {
                if (isMaster) cerr << "Error: Filename " << ROCKFILENAME << " not found or not readable." << endl;
                usageandexit();
            }
            rockfile.close();
            MonotCubicInterpolator Krwtmp;
            try {
                Krwtmp = MonotCubicInterpolator(ROCKFILENAME, 1, relPermCurves[waterPhaseIndex]); 
            }
            catch (const char * errormessage) {
                if (isMaster) cerr << "Error: " << errormessage << " Filename: " << ROCKFILENAME << endl;
                if (isMaster) cerr << "Check filename and -waterCurve" << endl;
                usageandexit();
            }
            if (!Krwtmp.isStrictlyMonotone()) {
                if (isMaster) cerr << "Error: Data in column " << relPermCurves[waterPhaseIndex] << " of file " << ROCKFILENAME << endl <<
                    "       was not strictly monotone. Exiting." << endl;
                usageandexit();
            }
            MonotCubicInterpolator Krotmp;
            try {
                Krotmp = MonotCubicInterpolator(ROCKFILENAME, 1, relPermCurves[oilPhaseIndex]); 
            }
            catch (const char * errormessage) {
                if (isMaster) cerr << "Error: " << errormessage << " Filename: " << ROCKFILENAME <<  endl;
                if (isMaster) cerr << "Check filename and -oilCurve" << endl;
                usageandexit();
            }
            if (!Krotmp.isStrictlyMonotone()) {
                if (isMaster) cerr << "Error: Data in column " << relPermCurves[oilPhaseIndex] << " of file " << ROCKFILENAME << endl <<
                    "       was not strictly monotone. Exiting." << endl;
                usageandexit();
            }
            if ((Krwtmp.isStrictlyIncreasing() && Krotmp.isStrictlyIncreasing()) || (Krwtmp.isStrictlyDecreasing() && Krotmp.isStrictlyDecreasing())) {
                if (isMaster) cerr << "Error: Input relperm curves are both increasing or decreasing in file " << ROCKFILENAME << endl;
                usageandexit();
            }
            for (int i=0 ; i < stone_types; ++i) { //Insert the same input curves for all rock types
                Krw.push_back(Krwtmp); 
                Kro.push_back(Krotmp); 
                rockTypeNames.push_back(ROCKFILENAME); 
            }
                       
            if (isMaster) cout << "Loaded rock file: " << ROCKFILENAME  
                               << ", for all stone types" << endl; 
        }
        else { 
            cerr << "Error:  Wrong number of stone-functions provided. " << endl 
                 << "Note that all input arguments after eclipse file are " << endl  
                 << "interpreted as stone functions." << endl; 
            return 1; 
        } 
    }
    else {   // Anisotropic case here! (double set of input curves needed
        cout << "rock types: " << stone_types << endl;
        cout << "varnum: " << varnum << endl;
        cout << "rockfileindex: " << rockfileindex << endl;
        if (varnum == rockfileindex + 2*stone_types) { // two .txt for each stone type in anisotropic case
            int rockidx=0;
            for (int i=0 ; i < 2*stone_types; i+=2) {  
                rockidx++;
                const char* WATERFILENAME = vararg[rockfileindex+i];
                const char* OILFILENAME = vararg[rockfileindex+i+1];
                // Check if rock files exists and is readable:
                ifstream waterfile(WATERFILENAME, ios::in);
                ifstream oilfile(OILFILENAME, ios::in);
                if (waterfile.fail()) {
                    if (isMaster) cerr << "Error: Filename " << WATERFILENAME << " not found or not readable." << endl;
                    usageandexit();
                }
                waterfile.close(); 
                if (oilfile.fail()) {
                    if (isMaster) cerr << "Error: Filename " << OILFILENAME << " not found or not readable." << endl;
                    usageandexit();
                }
                oilfile.close(); 
 
                MonotCubicInterpolator Krwxtmp, Krwytmp, Krwztmp, Kroxtmp, Kroytmp, Kroztmp;
                try {
                    Krwxtmp = MonotCubicInterpolator(WATERFILENAME, 2, 3); 
                    Krwytmp = MonotCubicInterpolator(WATERFILENAME, 2, 4); 
                    Krwztmp = MonotCubicInterpolator(WATERFILENAME, 2, 5); 
                }
                catch (const char * errormessage) {
                    if (isMaster) cerr << "Error: " << errormessage << " Filename: " << WATERFILENAME << endl;
                    if (isMaster) cerr << "Check filename and relpermcurves" << endl;
                    usageandexit();
                }
                if (!Krwxtmp.isStrictlyMonotone()) {
                    if (isMaster) cerr << "Error: Data in column " << 3 << " of file " << WATERFILENAME << endl <<
                        "       was not strictly monotone. Exiting." << endl;
                    usageandexit();
                }
                if (!Krwytmp.isStrictlyMonotone()) {
                    if (isMaster) cerr << "Error: Data in column " << 4 << " of file " << WATERFILENAME << endl <<
                        "       was not strictly monotone. Exiting." << endl;
                    usageandexit();
                }
                if (!Krwztmp.isStrictlyMonotone()) {
                    if (isMaster) cerr << "Error: Data in column " << 5 << " of file " << WATERFILENAME << endl <<
                        "       was not strictly monotone. Exiting." << endl;
                    usageandexit();
                }
                try {
                    Kroxtmp = MonotCubicInterpolator(OILFILENAME, 2, 3); 
                    Kroytmp = MonotCubicInterpolator(OILFILENAME, 2, 4); 
                    Kroztmp = MonotCubicInterpolator(OILFILENAME, 2, 5); 
                }
                catch (const char * errormessage) {
                    if (isMaster) cerr << "Error: " << errormessage << " Filename: " << OILFILENAME << endl;
                    if (isMaster) cerr << "Check filename and relpermcurves" << endl;
                    usageandexit();
                }
                if (!Kroxtmp.isStrictlyMonotone()) {
                    if (isMaster) cerr << "Error: Data in column " << 3 << " of file " << OILFILENAME << endl <<
                        "       was not strictly monotone. Exiting." << endl;
                    usageandexit();
                }
                if (!Kroytmp.isStrictlyMonotone()) {
                    if (isMaster) cerr << "Error: Data in column " << 4 << " of file " << OILFILENAME << endl <<
                        "       was not strictly monotone. Exiting." << endl;
                    usageandexit();
                }
                if (!Kroztmp.isStrictlyMonotone()) {
                    if (isMaster) cerr << "Error: Data in column " << 5 << " of file " << OILFILENAME << endl <<
                        "       was not strictly monotone. Exiting." << endl;
                    usageandexit();
                }
                Krwx.push_back(Krwxtmp); Krwy.push_back(Krwytmp); Krwz.push_back(Krwztmp); 
                Krox.push_back(Kroxtmp); Kroy.push_back(Kroytmp); Kroz.push_back(Kroztmp); 
                if (isMaster) cout << "Loaded rock files: " << WATERFILENAME << " and " 
                                   << OILFILENAME  << ", for stone type " << rockidx << endl; 
                rockTypeNames.push_back(WATERFILENAME);       
                rockTypeNames.push_back(OILFILENAME);       
            }
        }
        else if (varnum == rockfileindex + 2) { // one waterfile and one oilfile for all stone types 
            const char* WATERFILENAME = vararg[rockfileindex];
            const char* OILFILENAME = vararg[rockfileindex+1];
            // Check if rock files exists and is readable:
            ifstream waterfile(WATERFILENAME, ios::in);
            ifstream oilfile(OILFILENAME, ios::in);
            if (waterfile.fail()) {
                if (isMaster) cerr << "Error: Filename " << WATERFILENAME << " not found or not readable." << endl;
                usageandexit();
            }
            waterfile.close(); 
            if (oilfile.fail()) {
                if (isMaster) cerr << "Error: Filename " << OILFILENAME << " not found or not readable." << endl;
                usageandexit();
            }
            oilfile.close(); 
            MonotCubicInterpolator Krwxtmp, Krwytmp, Krwztmp;
            try {
                Krwxtmp = MonotCubicInterpolator(WATERFILENAME, 2, 3); 
                Krwytmp = MonotCubicInterpolator(WATERFILENAME, 2, 4); 
                Krwztmp = MonotCubicInterpolator(WATERFILENAME, 2, 5); 
            }
            catch (const char * errormessage) {
                if (isMaster) cerr << "Error: " << errormessage << " Filename: " << WATERFILENAME << endl;
                if (isMaster) cerr << "Check filename and relpermcurves" << endl;
                usageandexit();
            }
            if (!Krwxtmp.isStrictlyMonotone()) {
                if (isMaster) cerr << "Error: Data in column " << 3 << " of file " << WATERFILENAME << endl <<
                    "       was not strictly monotone. Exiting." << endl;
                usageandexit();
            }
            if (!Krwytmp.isStrictlyMonotone()) {
                if (isMaster) cerr << "Error: Data in column " << 4 << " of file " << WATERFILENAME << endl <<
                    "       was not strictly monotone. Exiting." << endl;
                usageandexit();
            }
            if (!Krwztmp.isStrictlyMonotone()) {
                if (isMaster) cerr << "Error: Data in column " << 5 << " of file " << WATERFILENAME << endl <<
                    "       was not strictly monotone. Exiting." << endl;
                usageandexit();
            }
            MonotCubicInterpolator Kroxtmp, Kroytmp, Kroztmp;
            try {
                Kroxtmp = MonotCubicInterpolator(OILFILENAME, 2, 3); 
                Kroytmp = MonotCubicInterpolator(OILFILENAME, 2, 4); 
                Kroztmp = MonotCubicInterpolator(OILFILENAME, 2, 5); 
            }
            catch (const char * errormessage) {
                if (isMaster) cerr << "Error: " << errormessage << " Filename: " << OILFILENAME << endl;
                if (isMaster) cerr << "Check filename and relpermcurves" << endl;
                usageandexit();
            }
            if (!Kroxtmp.isStrictlyMonotone()) {
                if (isMaster) cerr << "Error: Data in column " << 3 << " of file " << OILFILENAME << endl <<
                    "       was not strictly monotone. Exiting." << endl;
                usageandexit();
            }
            if (!Kroytmp.isStrictlyMonotone()) {
                if (isMaster) cerr << "Error: Data in column " << 4 << " of file " << OILFILENAME << endl <<
                    "       was not strictly monotone. Exiting." << endl;
                usageandexit();
            }
            if (!Kroztmp.isStrictlyMonotone()) {
                if (isMaster) cerr << "Error: Data in column " << 5 << " of file " << OILFILENAME << endl <<
                    "       was not strictly monotone. Exiting." << endl;
                usageandexit();
            }
            for (int i=0 ; i < 2*stone_types; i+=2) {  
                Krwx.push_back(Krwxtmp); Krwy.push_back(Krwytmp); Krwz.push_back(Krwztmp); 
                Krox.push_back(Kroxtmp); Kroy.push_back(Kroytmp); Kroz.push_back(Kroztmp); 
                rockTypeNames.push_back(WATERFILENAME);       
                rockTypeNames.push_back(OILFILENAME);       
            }
            if (isMaster) cout << "Loaded rock files: " << WATERFILENAME << " and " 
                               << OILFILENAME  << ", for all stone types" << endl; 
 
        }
        else { 
            cerr << "Error:  Wrong number of stone-functions provided. " << endl 
                 << "Note that all input arguments after eclipse file are " << endl  
                 << "interpreted as input functions." << endl; 
            return 1; 
        } 
    }


    /*****************************************************************************
     * Step 4:
     * Generate tesselated grid:
     * This is a step needed for the later discretization code to figure out which 
     * cells are connected to which. Each cornerpoint-cell is tesselated into 8 tetrahedrons.
     */
    if (isMaster) cout << "Tesselating grid... ";
    flush(cout);   start = clock();
    SinglePhaseUpscaler upscaler;
    double ztol = 0.0;
    double linsolver_tolerance = atof(options["linsolver_tolerance"].c_str());
    int linsolver_verbosity = atoi(options["linsolver_verbosity"].c_str());
    int linsolver_type = atoi(options["linsolver_type"].c_str());
    bool twodim_hack = false;
    eclParser.convertToSI();
    upscaler.init(eclParser, boundaryCondition,
                  Opm::unit::convert::from(minPerm, Opm::prefix::milli*Opm::unit::darcy),
                  ztol, linsolver_tolerance, linsolver_verbosity, linsolver_type, twodim_hack);
 
    finish = clock();   timeused_tesselation = (double(finish)-double(start))/CLOCKS_PER_SEC;
    if (isMaster) cout << " (" << timeused_tesselation <<" secs)" << endl;
 
    /******************************************************************************
     * Step 5:
     * Loop over cells to calculate
     *  - total volume for each rock type
     *  - Upscaled Swir
     *  - Upscaled Swor
     *  - maximum single phase perm, to be used for permeability contrast control 
     * 
     * Find total volume for each stone type
     * This is information to be used later on as it allows for faster computation
     * (because in the viscous limit, only the rock type matters for each cells
     * saturation distribution, not its porosity or permeability)
     */

    vector<double> cellVolumes, cellPoreVolumes; // values for each cell (check if really needed)
    cellVolumes.resize(satnums.size(), 0.0);
    cellPoreVolumes.resize(satnums.size(), 0.0);

    vector<double> rocktypeVolume; // rocktypeVolume[rockIdx] == aggregated pore-volume for that type
    rocktypeVolume.resize(stone_types);
 
    int tesselatedCells = 0;
    
    double maxSinglePhasePerm = 0;
    double Swirvolume = 0;
    double Sworvolume = 0;
 
    // Loop through all cells. Add the cells porevolume to the corresponding rock type volume 
    // Also determine bounds for fractional flow ratio 
    const std::vector<int>& ecl_idx = upscaler.grid().globalCell();
    Dune::CpGrid::Codim<0>::LeafIterator c = upscaler.grid().leafbegin<0>();
    for (; c != upscaler.grid().leafend<0>(); ++c) {
        unsigned int cell_idx = ecl_idx[c->index()];
        if (satnums[cell_idx] > 0) { // Satnum zero is "no rock"
            cellVolumes[cell_idx]     = c->geometry().volume();
            cellPoreVolumes[cell_idx] = cellVolumes[cell_idx] * poros[cell_idx]; 
            
            rocktypeVolume[int(satnums[cell_idx])-1] += cellPoreVolumes[cell_idx]; 
            
            // Also find max single-phase perm in input file:
            maxSinglePhasePerm = max( maxSinglePhasePerm, permxs[cell_idx]);
            double minSw, maxSw;
            if (! anisotropic_input) {
                minSw = Krw[int(satnums[cell_idx])-1].getMinimumX().first;
                maxSw = Krw[int(satnums[cell_idx])-1].getMaximumX().first;
                //cout << "minSwc: " << minSw << endl;
                //cout << "maxSwc: " << maxSw << endl;
            }
            else {
                minSw = Krwx[int(satnums[cell_idx])-1].getMinimumX().first;
                maxSw = Krwx[int(satnums[cell_idx])-1].getMaximumX().first;               
            }
            // Add irreducible water saturation volume
            Swirvolume += minSw * cellPoreVolumes[cell_idx];
            Sworvolume += maxSw * cellPoreVolumes[cell_idx];
        }
        ++tesselatedCells; // keep count  (also counts non-rock-cells)
    } 
    
    double minSinglePhasePerm = max(maxSinglePhasePerm/maxPermContrast, minPerm);
 
    // Total porevolume and total volume -> upscaled porosity: 
    double poreVolume = accumulate(cellPoreVolumes.begin(),  
                                   cellPoreVolumes.end(), 
                                   0.0); 
    double volume = accumulate(cellVolumes.begin(), 
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

    /*********************************************************************** 
     * Step 6 
     * Upscale fractional flow ratio vs water saturation
     *
     * This is upscaled in advance in order to be able to have uniformly distributed
     * saturation points for which upscaling is performed.
     *
     * Fractional flow ratio points are chosen heuristically in order to
     * ensure largest saturation interval between two saturation points
     * is 1/500 of the saturation interval. Monotone cubic
     * interpolation will be used afterwards for accessing the
     * tabulated values.
     */ 
    
    MonotCubicInterpolator WaterSaturationVsFractionalFlow;
    double largestSaturationInterval = Swor - Swir;
 
    // Make fractional flow rato vs watersaturation curve for each stone type 
    const int interpPoints = 1000; 
    const double lowerRelPermLimit = 1e-12; 
    vector<MonotCubicInterpolator> FracFlowRatioInv; // We only need the inverse of this 
    vector<MonotCubicInterpolator> FracFlowRatioInvY, FracFlowRatioInvZ;
    for (int stone_idx = 0; stone_idx < stone_types; ++stone_idx) { 
        MonotCubicInterpolator waterCurve, oilCurve, 
            waterCurveY, oilCurveY, waterCurveYInv, oilCurveYInv,
            waterCurveZ, oilCurveZ, waterCurveZInv, oilCurveZInv; 
       // NOTE: This is a huge hack for anisotropic_input, we only look at
       // x-directions, and hope this is good enough for establishing
       // which fracflows to choose. The only side effect is probably slightly 
       // variations in distance between each saturation point.
        if (! anisotropic_input) {
            waterCurve = MonotCubicInterpolator(Krw[stone_idx]);
            oilCurve = MonotCubicInterpolator(Kro[stone_idx]); 
        }
        else {
            waterCurve = MonotCubicInterpolator(Krwx[stone_idx]);
            oilCurve = MonotCubicInterpolator(Krox[stone_idx]); 
            waterCurveY   = MonotCubicInterpolator(Krwy[stone_idx]);
            oilCurveY     = MonotCubicInterpolator(Kroy[stone_idx]);
            waterCurveZ   = MonotCubicInterpolator(Krwz[stone_idx]);
            oilCurveZ     = MonotCubicInterpolator(Kroz[stone_idx]);            
        }
        MonotCubicInterpolator waterCurveInv(waterCurve.get_fVector(), waterCurve.get_xVector()); 
        MonotCubicInterpolator oilCurveInv(oilCurve.get_fVector(), oilCurve.get_xVector()); 
        // Finding min and max bounds for watersaturation 
        double waterSwMin = waterCurve.getMinimumX().first; 
        double waterSwMax = waterCurve.getMaximumX().first; 
        double oilSwMin = oilCurve.getMinimumX().first; 
        double oilSwMax = oilCurve.getMaximumX().first; 
        if (anisotropic_input) {
            waterSwMin = min(waterSwMin, waterCurveY.getMinimumX().first);
            waterSwMax = max(waterSwMax, waterCurveY.getMaximumX().first);
            waterSwMin = min(waterSwMin, waterCurveZ.getMinimumX().first);
            waterSwMax = max(waterSwMax, waterCurveZ.getMaximumX().first);
            oilSwMin = min(oilSwMin, oilCurveY.getMinimumX().first);
            oilSwMax = max(oilSwMax, oilCurveY.getMaximumX().first);
            oilSwMin = min(oilSwMin, oilCurveZ.getMinimumX().first);
            oilSwMax = max(oilSwMax, oilCurveZ.getMaximumX().first);
        }        
        // Next two lines limit the Sw-span. CHECK!!
        double leftBound = min(oilSwMin, waterSwMin); // Maximum of the two curves left bounds 
        double rightBound = max(oilSwMax, waterSwMax); // Minimum of the two curves right bounds 
        // Make sure the rel perm values are not too small (avoid infinity in the fracflow) 
        if (waterCurve.evaluate(leftBound) < lowerRelPermLimit) { 
            leftBound = waterCurveInv.evaluate(lowerRelPermLimit); 
        } 
        if (oilCurve.evaluate(rightBound) < lowerRelPermLimit) { 
            rightBound = oilCurveInv.evaluate(lowerRelPermLimit); 
            // Assume that also holds for aniso input
        } 
        
        // Uniformly distributed watersaturation points for each rocktype
        vector<double> Sw(interpPoints, leftBound); 
        vector<double> f(interpPoints, 0); 
        vector<double> fY(interpPoints, 0); 
        vector<double> fZ(interpPoints, 0); 
        double stepSize = (rightBound-leftBound)/(interpPoints-1); 
        for (int i = 0; i < interpPoints; ++i) { 
            Sw[i] += i*stepSize; 
            f[i] = max(waterCurve.evaluate(Sw[i]), lowerRelPermLimit) / 
                max(oilCurve.evaluate(Sw[i]), lowerRelPermLimit)
                * viscosities[oilPhaseIndex]/viscosities[waterPhaseIndex]; 
            if (anisotropic_input) {
                fY[i] = max(waterCurveY.evaluate(Sw[i]), lowerRelPermLimit) / 
                    max(oilCurveY.evaluate(Sw[i]), lowerRelPermLimit)
                    * viscosities[oilPhaseIndex]/viscosities[waterPhaseIndex]; 
                fZ[i] = max(waterCurveZ.evaluate(Sw[i]), lowerRelPermLimit) / 
                    max(oilCurveZ.evaluate(Sw[i]), lowerRelPermLimit)
                    * viscosities[oilPhaseIndex]/viscosities[waterPhaseIndex]; 
            }
        }     
        FracFlowRatioInv.push_back(MonotCubicInterpolator(f,Sw));   
        if (anisotropic_input) {
            FracFlowRatioInvY.push_back(MonotCubicInterpolator(fY,Sw));
            FracFlowRatioInvZ.push_back(MonotCubicInterpolator(fY,Sw));
        }
    } 
    
    // Find max/min fracflowratio over all rock types.
    double fracflowratioMin = numeric_limits<double>().max(); 
    double fracflowratioMax = 0; 
    for (int rockIdx = 0; rockIdx < stone_types; ++rockIdx) { 
        double fracflowratioMinRock = FracFlowRatioInv[rockIdx].getMinimumX().first; 
        double fracflowratioMaxRock = FracFlowRatioInv[rockIdx].getMaximumX().first; 
        fracflowratioMin = min(fracflowratioMin, fracflowratioMinRock); 
        fracflowratioMax = max(fracflowratioMax, fracflowratioMaxRock); 
        if (anisotropic_input) {
            double fracflowratioMinRockY = FracFlowRatioInvY[rockIdx].getMinimumX().first; 
            double fracflowratioMaxRockY = FracFlowRatioInvY[rockIdx].getMaximumX().first; 
            fracflowratioMin = min(fracflowratioMin, fracflowratioMinRockY); 
            fracflowratioMax = max(fracflowratioMax, fracflowratioMaxRockY); 
            double fracflowratioMinRockZ = FracFlowRatioInvZ[rockIdx].getMinimumX().first; 
            double fracflowratioMaxRockZ = FracFlowRatioInvZ[rockIdx].getMaximumX().first; 
            fracflowratioMin = min(fracflowratioMin, fracflowratioMinRockZ); 
            fracflowratioMax = max(fracflowratioMax, fracflowratioMaxRockZ); 
        }  
    } 
    
    if (isMaster) cout << endl << "Lower fracflowratio: " << fracflowratioMin << ", Upper fracflowratio: " << fracflowratioMax << endl; 
    // Now upscale fractional flow vs water saturation 
    // (i.e., populate the vector WaterSaturationVsFractionalFlow)
 
    double fracFlowRatioTestvalue;
    
    while (largestSaturationInterval > (Swor-Swir)/500.0) {
        if (fracflowratioMax == fracflowratioMin) {
            // This is a dummy situation, we go through once and then 
            // we are finished (this will be triggered by zero permeability)
            // (not checked for applicatiblity in viscous limit, this is copied from cap limit)
            fracFlowRatioTestvalue = fracflowratioMin;
            largestSaturationInterval = 0;
        }
        else if (WaterSaturationVsFractionalFlow.getSize() == 0) {
            /* No data previously computed */
            fracFlowRatioTestvalue = fracflowratioMax;
        }
        else if (WaterSaturationVsFractionalFlow.getSize() == 1) {
            /* Means that this is second pass through this while-loop */
            fracFlowRatioTestvalue = fracflowratioMin;
        }
        else {
            /* Search for largest saturation interval in which there are no
               computed saturation points (and estimate the fracflow ratio
               that will fall in the center of this saturation interval)
            */
            pair<double,double> SatDiff = WaterSaturationVsFractionalFlow.getMissingX();
            fracFlowRatioTestvalue = SatDiff.first;
            largestSaturationInterval = SatDiff.second;
        }
            
        // Check for saneness of fracFlowRatioTestvalue
        if (std::isnan(fracFlowRatioTestvalue) || std::isinf(fracFlowRatioTestvalue)) {
            if (isMaster) cerr << "ERROR: fracFlowRatioTestvalue was inf or nan." << endl;
            break; // Jump out out while-loop, just print the results up to now and exit
        }
 
        // Do the saturation modelling with the current fracFlowRatioTestvalue
        double waterVolume = 0.0;
        vector<double> waterSaturationRockType;
        waterSaturationRockType.resize(stone_types);
        
        for (int rockIdx = 0; rockIdx < stone_types; ++rockIdx) {
            waterSaturationRockType[rockIdx] = FracFlowRatioInv[rockIdx].evaluate(fracFlowRatioTestvalue);
            if (anisotropic_input) {
                waterSaturationRockType[rockIdx] += FracFlowRatioInvY[rockIdx].evaluate(fracFlowRatioTestvalue);
                waterSaturationRockType[rockIdx] += FracFlowRatioInvZ[rockIdx].evaluate(fracFlowRatioTestvalue);
                waterSaturationRockType[rockIdx] /= 3.0; // arithmetic average of three directions.
            }
            waterVolume += waterSaturationRockType[rockIdx] * rocktypeVolume[rockIdx];
        }
        WaterSaturationVsFractionalFlow.addPair(fracFlowRatioTestvalue, waterVolume/poreVolume);
    }
    // In case we have all flat endpoints in this curve, we chop them off:
    // (We still preserve Swir and Swor)
    WaterSaturationVsFractionalFlow.chopFlatEndpoints(saturationThreshold);
 
    // In case we have a monotone function, but not strictly monotone,
    // we can remove some data points in order to make it strictly
    // monotone.  This has been seen to happen on one occasion, quite
    // close to the saturation endpoint. Either we have to do this, or
    // increase lowerRelPermLimit from 1e-12 to 1e-8. We might as well
    // do both.
    WaterSaturationVsFractionalFlow.shrinkFlatAreas();
    
    // Now we can also invert the upscaled water saturation
    // (it should be monotonic)
    if (!WaterSaturationVsFractionalFlow.isStrictlyMonotone()) {
        if (isMaster) {
            cout << WaterSaturationVsFractionalFlow.toString();
            cerr << "Error: Upscaled water saturation not strictly monotone in fractional flow ratio." << endl;
            cerr << "       Unphysical input data, exiting." << endl;
        }
        usageandexit();
    }
    MonotCubicInterpolator FractionalFlowVsWaterSaturation(WaterSaturationVsFractionalFlow.get_fVector(), 
                                                           WaterSaturationVsFractionalFlow.get_xVector());
    
   
    /*****************************************************************************
     * Step 7:
     * Upscale single phase permeability
     * This uses the PERMX in the eclipse file as data, and upscales using
     * fixed boundary (no-flow) conditions
     *
     * In an MPI-environment, this is only done on the master node.
     */
    typedef SinglePhaseUpscaler::permtensor_t Matrix;
    Matrix zeroMatrix(3,3,(double*)0);
    zero(zeroMatrix);
    Matrix permTensor = zeroMatrix;
    Matrix permTensorInv = zeroMatrix;
 
    if (isMaster) {
        //cout << "Rank " << mpi_rank << " upscaling single-phase permeability..."; flush(cout);
        Matrix cellperm = zeroMatrix;
        for (unsigned int i = 0; i < ecl_idx.size(); ++i) {
            unsigned int cell_idx = ecl_idx[i];
            zero(cellperm);
            if (! anisotropic_input) {
                double kval = max(permxs[cell_idx], minSinglePhasePerm);
                cellperm(0,0) = kval;
                cellperm(1,1) = kval;
                cellperm(2,2) = kval;
            }
            else {
                cellperm(0,0) = max(minSinglePhasePerm, permxs[cell_idx]);
                cellperm(1,1) = max(minSinglePhasePerm, permys[cell_idx]);
                cellperm(2,2) = max(minSinglePhasePerm, permzs[cell_idx]);
            }
            upscaler.setPermeability(i, cellperm);
        }
        permTensor = upscaler.upscaleSinglePhase();
        permTensorInv = permTensor;
        invert(permTensorInv);
    }

   /*****************************************************************
    * Step 8 and 9  
    * (step 8 is for water, 9 is for oil, the code is identical, we just 
    * swap some variables) 
    * 
    * For uniformly distributed water saturation values between upscaled Swir and Swor, 
    *    a: Find fractional flow ratio corresponding to wanted upscaled water saturation
    *    b: Model water saturation in each cell given fractional flow ratio
    *    c: Model phase permeability in each cell given cell water saturation and inputted
    *       relative permeability curves.
    *    d: Upscale phase permeability
    */

    vector<double> WaterSaturation;
    
    vector<vector<vector<double> > > PhasePerm; // 'phases' * 'tensorElementCount' phaseperm values per fracflowpoint
 
    vector<vector<double> > tmp1;
    vector<vector<double> > tmp2;
    PhasePerm.push_back(tmp1);
    PhasePerm.push_back(tmp2);
 
    // Put empty interpolator objects for the upscaled results, number of 
    // interpolators depending on the boundary condition. 
    for (int idx=0; idx < points; ++idx) {
        WaterSaturation.push_back(0.0); // Pad with zeroes
        vector<double> foo1, foo2;
        PhasePerm[waterPhaseIndex].push_back(foo1);
        PhasePerm[oilPhaseIndex].push_back(foo2);
        for (int voigtIdx=0; voigtIdx < tensorElementCount; ++voigtIdx) { 
            PhasePerm[waterPhaseIndex][idx].push_back(0.0); // Pad with zeroes. 
            PhasePerm[oilPhaseIndex][idx].push_back(0.0);  
        }    
    }

    // Make vector of frac flow ratio points corresponding to uniformly distributed 
    // saturation points between Swor and Swir
    vector <double> fracFlowRatioPoints;
    for (int pointidx=1; pointidx <= points; ++pointidx) {
        // pointidx=1 corresponds to Swir, pointidx=points to Swor.
        double saturation = Swir + (Swor-Swir)/(points-1)*(pointidx-1);
        fracFlowRatioPoints.push_back(FractionalFlowVsWaterSaturation.evaluate(saturation));
    }
 
    // Construct a vector that tells for each fracflowratio point which mpi-node (rank) should compute for that
    // particular fracflowratio point
    vector<int> node_vs_fracflowratiopoint;
    // Fill with zeros initially (in case of non-mpi)
    for (int idx=0; idx < points; ++idx) {
        node_vs_fracflowratiopoint.push_back(0);
    }
    
#if USEMPI
    // Distribute work load over mpi nodes.
    for (int idx=0; idx < points; ++idx) {
        // Ensure master node gets equal or less work than the other nodes, since
        // master node also computes single phase perm.
        node_vs_fracflowratiopoint[idx] = (mpi_nodecount-1) - idx % mpi_nodecount;
        /*if (isMaster) {
          cout << "Fracflowratio point " << idx << " assigned to node " << node_vs_fracflowratiopoint[idx] << endl;
          }*/
    }   
#endif

    clock_t start_upscale_wallclock = clock();

    double waterVolumeLF; // water volume for the whole model
    // Now loop through the vector of fractional flow ratios that
    // this node should compute.
    for (int phase = waterPhaseIndex; phase <= oilPhaseIndex; ++phase) {
        
        string phaseName;
        if (phase == waterPhaseIndex) {
            phaseName = string("water");
        }
        else if (phase == oilPhaseIndex) {
            phaseName = string("oil");
        }
        if (isMaster) cout << endl << "Upscaling relative permeability for " << phaseName << "... " << endl;
        for (int pointidx = 0; pointidx < points; ++pointidx) {
            
            // Should "I" (mpi-wise) compute this fracflowratio point?
            if (node_vs_fracflowratiopoint[pointidx] == mpi_rank) {
                
                fracFlowRatioTestvalue = fracFlowRatioPoints[pointidx];
                
                double accPhasePerm = 0.0; // accumulated, can be used for debugging
                
                double maxPhasePerm = 0.0;
                
                vector<double> phasePermValues;
                vector<vector<double> > phasePermValuesDiag;
                phasePermValues.resize(satnums.size());
                phasePermValuesDiag.resize(satnums.size());
                waterVolumeLF = 0.0;
 
                vector<double> waterSaturationRockType;
                waterSaturationRockType.resize(stone_types);
                
                for (int rockIdx = 0; rockIdx < stone_types; ++rockIdx) {
                    waterSaturationRockType[rockIdx] = FracFlowRatioInv[rockIdx].evaluate(fracFlowRatioTestvalue);
                    if (anisotropic_input) {
                        waterSaturationRockType[rockIdx] += FracFlowRatioInvY[rockIdx].evaluate(fracFlowRatioTestvalue);
                        waterSaturationRockType[rockIdx] += FracFlowRatioInvZ[rockIdx].evaluate(fracFlowRatioTestvalue);
                        waterSaturationRockType[rockIdx] /= 3.0; // arithmetic average of three directions.
                    }
                    waterVolumeLF += waterSaturationRockType[rockIdx] * rocktypeVolume[rockIdx];
                }   
                for (unsigned int i = 0; i < ecl_idx.size(); ++i) {
                    unsigned int cell_idx = ecl_idx[i];
                    double cellPhasePerm = minPerm;
                    vector<double>  cellPhasePermDiag;
                    cellPhasePermDiag.push_back(minPerm);
                    cellPhasePermDiag.push_back(minPerm);
                    cellPhasePermDiag.push_back(minPerm);

                    if (satnums[cell_idx] > 0) { // Satnum zero is "no rock", model those with minPerm.
                        
                        // Water saturation is only a function of the rock type
                        double saturationCell 
                            = waterSaturationRockType[int(satnums[cell_idx])-1];
                        if (! anisotropic_input) {
                            double cellRelPerm; 
                            if (phase == waterPhaseIndex) {
                                cellRelPerm = Krw[int(satnums[cell_idx])-1].evaluate(saturationCell);
                            }
                            else { // if (phase == oilPhaseIndex) {
                                cellRelPerm = Kro[int(satnums[cell_idx])-1].evaluate(saturationCell);
                            }
                            cellPhasePerm = cellRelPerm * permxs[cell_idx];
                        }
                        else { //anisotropic
                            if (phase == waterPhaseIndex) {
                                cellPhasePermDiag[0] = Krwx[int(satnums[cell_idx])-1].evaluate(saturationCell) * 
                                    permxs[cell_idx];
                                cellPhasePermDiag[1] = Krwy[int(satnums[cell_idx])-1].evaluate(saturationCell) * 
                                    permys[cell_idx];
                                cellPhasePermDiag[2] = Krwz[int(satnums[cell_idx])-1].evaluate(saturationCell) * 
                                    permzs[cell_idx];
                            }
                            else { // oil
                                cellPhasePermDiag[0] = Krox[int(satnums[cell_idx])-1].evaluate(saturationCell) * 
                                    permxs[cell_idx];
                                cellPhasePermDiag[1] = Kroy[int(satnums[cell_idx])-1].evaluate(saturationCell) * 
                                    permys[cell_idx];
                                cellPhasePermDiag[2] = Kroz[int(satnums[cell_idx])-1].evaluate(saturationCell) * 
                                    permzs[cell_idx];
                            }
                        }
                    }
                    phasePermValues[cell_idx] = cellPhasePerm;
                    phasePermValuesDiag[cell_idx] = cellPhasePermDiag;
                    maxPhasePerm = max(maxPhasePerm, cellPhasePerm);
                    maxPhasePerm = max(maxPhasePerm, *max_element(cellPhasePermDiag.begin(),
                                                                  cellPhasePermDiag.end()));
                }
                
                // Now we can determine the smallest permitted permeability we can calculate for
                // We have both a fixed bottom limit, as well as a possible higher limit determined
                // by a maximum allowable permeability.
                double minPhasePerm = max(maxPhasePerm/maxPermContrast, minPerm);
 
                // Now remodel the phase permeabilities obeying minPhasePerm.
                Matrix cellperm = zeroMatrix;
                for (unsigned int i = 0; i < ecl_idx.size(); ++i) {
                    unsigned int cell_idx = ecl_idx[i];
                    zero(cellperm);
                    if (! anisotropic_input) {
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
                    upscaler.setPermeability(i, cellperm);
                }
                Matrix phasePermTensor = upscaler.upscaleSinglePhase();
                
                // Here we recalculate the upscaled water saturation,
                // although it is already known when we asked for the
                // fracflowratio point to compute for. Nonetheless, we
                // recalculate here to avoid any minor roundoff-error and
                // interpolation error (this means that the saturation
                // points are not perfectly uniformly distributed)
                //cout << waterVolumeLF/poreVolume;
                WaterSaturation[pointidx] =  waterVolumeLF/poreVolume;
                
#ifdef USEMPI
                cout << "Rank " << mpi_rank << ": " << endl;;
#endif
                cout << fracFlowRatioTestvalue << "\t" << WaterSaturation[pointidx];
                // Store and print phase-perm-result
                for (int voigtIdx=0; voigtIdx < tensorElementCount; ++voigtIdx) { 
                    PhasePerm[phase][pointidx][voigtIdx] = getVoigtValue(phasePermTensor,voigtIdx); 
                    cout << "\t" << getVoigtValue(phasePermTensor,voigtIdx); 
                } 
                cout << endl; 
            }
        }  // end loop over saturation points
    } // end phase loop
    clock_t finish_upscale_wallclock = clock();
    timeused_upscale_wallclock = (double(finish_upscale_wallclock)-double(start_upscale_wallclock))/CLOCKS_PER_SEC;
    //double timeused_upscale_total = timeused_upscale_wallclock;
#ifdef USEMPI   
    /* Step Xb: Transfer all computed data to master node.
       Master node should post a receive for all values missing,
       other nodes should post a send for all the values they have.
    */
    MPI_Barrier(MPI_COMM_WORLD); // Not strictly necessary.
    if (isMaster) {
        // Loop over all values, receive data and put into local data structure
        for (int phase = waterPhaseIndex; phase <= oilPhaseIndex; ++phase) {
            for (int idx=0; idx < points; ++idx) {
                if (node_vs_fracflowratiopoint[idx] != 0) {
                    // Receive data
                    double recvbuffer[2+tensorElementCount];
                    MPI_Recv(recvbuffer, 2+tensorElementCount, MPI_DOUBLE, 
                             node_vs_fracflowratiopoint[idx], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    // Put received data into correct place.
                    WaterSaturation[(int)recvbuffer[0]] = recvbuffer[1];
                    for (int voigtIdx=0; voigtIdx < tensorElementCount; ++voigtIdx) {
                        PhasePerm[phase][(int)recvbuffer[0]][voigtIdx] = recvbuffer[2+voigtIdx];
                    }
                }
            }
        }
    }
    else {
        for (int phase = waterPhaseIndex; phase <= oilPhaseIndex; ++phase) {
            for (int idx=0; idx < points; ++idx) {
                if (node_vs_fracflowratiopoint[idx] == mpi_rank) {
                    // Pack and send data. C-style.
                    double sendbuffer[2+tensorElementCount];
                    sendbuffer[0] = (double)idx;
                    sendbuffer[1] = WaterSaturation[idx];
                    for (int voigtIdx=0; voigtIdx < tensorElementCount; ++voigtIdx) {
                        sendbuffer[2+voigtIdx] = PhasePerm[phase][idx][voigtIdx];
                    }
                    MPI_Send(sendbuffer, 2+tensorElementCount, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
                }
            }
        }
    }
#endif

    // Average time pr. upscaling point:
#ifdef USEMPI
    // Sum the upscaling time used by all processes
    double timeused_total;
    MPI_Reduce(&timeused_upscale_wallclock, &timeused_total, 1, MPI_DOUBLE, 
               MPI_SUM, 0, MPI_COMM_WORLD);
    double avg_upscaling_time_pr_point = timeused_total/(2.0*(double)points);
       
#else
    double avg_upscaling_time_pr_point = timeused_upscale_wallclock / (2.0*(double)points);
#endif

    /* 
     * Step Xc: Make relperm values from phaseperms
     *          (only master node can do this)
     */
   
    vector<vector<vector <double> > > RelPermValues; // phase is first index voigtIdx is second index.
    vector<vector <double> > tmp10, tmp20;
    RelPermValues.push_back(tmp10);
    RelPermValues.push_back(tmp20);
    for (int voigtIdx=0; voigtIdx < tensorElementCount; ++voigtIdx) {
        vector<double> tmpfoo,tmpbar;
        RelPermValues[waterPhaseIndex].push_back(tmpfoo);
        RelPermValues[oilPhaseIndex].push_back(tmpbar);
    }
    if (isMaster) {
        for (int phase = waterPhaseIndex; phase <= oilPhaseIndex; ++phase) {
            // Loop over all fracflowratio points 
            for (int idx=0; idx < points; ++idx) {
                Matrix phasePermTensor = zeroMatrix;
                zero(phasePermTensor);
                for (int voigtIdx = 0; voigtIdx < tensorElementCount; ++voigtIdx) {
                    setVoigtValue(phasePermTensor, voigtIdx, PhasePerm[phase][idx][voigtIdx]);
                }
                //cout << phasePermTensor << endl;
                Matrix relPermTensor = zeroMatrix;
                // relPermTensor = phasePermTensor;
                // relPermTensor *= permTensorInv;
                prod(phasePermTensor, permTensorInv, relPermTensor);
                for (int voigtIdx = 0; voigtIdx < tensorElementCount; ++voigtIdx) {
                    RelPermValues[phase][voigtIdx].push_back(getVoigtValue(relPermTensor, voigtIdx));
                }
                //cout << relPermTensor << endl;
            }
        }
    }

   /*********************************************************************************
    * Step 10
    *
    * Output results to stdout and optionally to file. Note, we only output to
    * file if the '-outputWater'-option and/or '-outputOil' has been set, as this option is an
    * empty string by default.
    */
   
   if (isMaster) {

       // If no data computed, we do not have more to do:
       if (WaterSaturation.size() == 0) {
           return(1); // non-zero return value, this means something wrong with input data.
       }
   
       stringstream outheadtmp;
       stringstream outwatertmp;
       stringstream outoiltmp;
       
       // Print a table of all computed values:
       outheadtmp << "######################################################################" << endl;
       outheadtmp << "# Results from upscaling relative permeability."<< endl;
       outheadtmp << "#" << endl;
       time_t now = std::time(NULL);
       outheadtmp << "# Finished: " << asctime(localtime(&now));
       
       utsname hostname;   uname(&hostname);
       outheadtmp << "# Hostname: " << hostname.nodename << endl;
       
       outheadtmp << "#" << endl;
       outheadtmp << "# Eclipse file: " << ECLIPSEFILENAME << endl;
       outheadtmp << "#        cells: " << tesselatedCells << endl;
       outheadtmp << "#  Pore volume: " << poreVolume << endl;
       outheadtmp << "#       volume: " << volume << endl;
       outheadtmp << "#     Porosity: " << poreVolume/volume << endl;
       outheadtmp << "#" << endl;
       for (int i=0; i < stone_types ; ++i) {
           if (! anisotropic_input) {
               outheadtmp << "# Stone " << i+1 << ": " << rockTypeNames[i] << " (" << Krw[i].getSize() << " points)" <<  endl;
           }
           else {
               outheadtmp << "# Stone " << i+1 << ": " << rockTypeNames[i] << " (" << Krwx[i].getSize() << " points)" <<  endl;
           }
       }
       outheadtmp << "#" << endl;
       outheadtmp << "# Timings:   Tesselation: " << timeused_tesselation << " secs" << endl;
       outheadtmp << "#              Upscaling: " << timeused_upscale_wallclock << " secs";
#ifdef USEMPI
       outheadtmp << " (wallclock time)" << endl;
       outheadtmp << "#                         " << avg_upscaling_time_pr_point << " secs pr. saturation point" << endl;
       outheadtmp << "#              MPI-nodes: " << mpi_nodecount << endl;

       // Single phase upscaling time is included here, in possibly a hairy way.
       double speedup = (avg_upscaling_time_pr_point * (2*points + 1) + timeused_tesselation)/(timeused_upscale_wallclock + avg_upscaling_time_pr_point + timeused_tesselation);
       outheadtmp << "#                Speedup: " << speedup << ", efficiency: " << speedup/mpi_nodecount << endl;
#else
       outheadtmp << ", " << avg_upscaling_time_pr_point << " secs avg for " << (2*points) << " runs" << endl;
#endif
       outheadtmp << "# " << endl;
       outheadtmp << "# Options used:" << endl;
       outheadtmp << "#     Boundary conditions: ";
       if (isFixed)    outheadtmp << "Fixed (no-flow)" << endl;
       if (isPeriodic) outheadtmp << "Periodic" << endl;
       if (isLinear)   outheadtmp << "Linear" << endl;
       outheadtmp << "#                  points: " << options["points"] <<  endl;
       outheadtmp << "#         maxPermContrast: " << options["maxPermContrast"] << endl;
       outheadtmp << "#                 minPerm: " << options["minPerm"] << endl;
       outheadtmp << "#                 minPoro: " << options["minPoro"] << endl;      
       outheadtmp << "#          waterViscosity: " << options["waterViscosity"] << " Pa s" << endl;
       outheadtmp << "#            oilViscosity: " << options["oilViscosity"] << " Pa s" << endl;
       if (doInterpolate) { 
           outheadtmp << "#             interpolate: " << options["interpolate"] << " points" << endl; 
       }
       outheadtmp << "# " << endl;
       outheadtmp << "# Single phase permeability" << endl;
       outheadtmp << "#  |Kxx  Kxy  Kxz| = " << permTensor(0,0) << "  " << permTensor(0,1) << "  " << permTensor(0,2) << endl;
       outheadtmp << "#  |Kyx  Kyy  Kyz| = " << permTensor(1,0) << "  " << permTensor(1,1) << "  " << permTensor(1,2) << endl;
       outheadtmp << "#  |Kzx  Kzy  Kzz| = " << permTensor(2,0) << "  " << permTensor(2,1) << "  " << permTensor(2,2) << endl;
       outheadtmp << "# " << endl;
       if (doInterpolate) { 
           outheadtmp << "# NB: Data points shown are interpolated." << endl; 
       } outheadtmp << "######################################################################" << endl;
       
       if (isFixed) {
           outwatertmp << "#     v_w/v_o          Sw          Krwxx         Krwyy         Krwzz" << endl;
           outoiltmp   << "#     v_w/v_o          Sw          Kroxx         Kroyy         Krozz" << endl;
       }
       else if (isLinear) {
           outwatertmp << "#     v_w/v_o          Sw          Krwxx         Krwyy         Krwzz         Krwyz         Krwxz         Krwxy         Krwzy         Krwzx         Krwyx" << endl;
           outoiltmp   << "#     v_w/v_o          Sw          Kroxx         Kroyy         Krozz         Kroyz         Kroxz         Kroxy         Krozy         Krozx         Kroyx" << endl;
       }
       else if (isPeriodic) {
           outwatertmp << "#     v_w/v_o          Sw          Krwxx         Krwyy         Krwzz         Krwyz         Krwxz         Krwxy         Krwzy         Krwzx         Krwyx" << endl;
           outoiltmp   << "#     v_w/v_o          Sw          Kroxx         Kroyy         Krozz         Kroyz         Kroxz         Kroxy         Krozy         Krozx         Kroyx" << endl;
       }
       
       // If user wants interpolated output, do monotone cubic interpolation 
       // by modifying the data vectors that are to be printed 
       if (doInterpolate) {
           // Find min and max for saturation values 
           double satmin = +DBL_MAX;
           double satmax = -DBL_MAX;
           for (unsigned int i = 0; i < WaterSaturation.size(); ++i) { 
               if (WaterSaturation[i] < satmin) { 
                   satmin = WaterSaturation[i]; 
               } 
               if (WaterSaturation[i] > satmax) { 
                   satmax = WaterSaturation[i]; 
               } 
           } 

           // Make uniform grid in saturation axis (for both water and oil)
           vector<double> SatvaluesInterp, fracFlowRatioPointsInterp;
           for (int i = 0; i < interpolationPoints; ++i) { 
               SatvaluesInterp.push_back(satmin + ((double)i)/((double)interpolationPoints-1)*(satmax-satmin)); 
               fracFlowRatioPointsInterp.push_back(WaterSaturationVsFractionalFlow.evaluate(SatvaluesInterp[i]));
           } 

           // Now fracflowratio and computed relperm-values must be viewed as functions 
           // of saturation, and then interpolated on the uniform saturation grid. 
           
           // Now overwrite existing FlowRatioValues and relperm-data with interpolated data: 
           MonotCubicInterpolator FracFlowRatioValuesVsSaturationWater, FracFlowRatioValuesVsSaturationOil;
           
           for (int voigtIdx = 0; voigtIdx < tensorElementCount; ++voigtIdx) { 
               MonotCubicInterpolator RelPermWaterVsSaturation(WaterSaturation, RelPermValues[waterPhaseIndex][voigtIdx]); 
               MonotCubicInterpolator RelPermOilVsSaturation(  WaterSaturation, RelPermValues[oilPhaseIndex][voigtIdx]); 
               RelPermValues[waterPhaseIndex][voigtIdx].clear(); 
               RelPermValues[oilPhaseIndex][voigtIdx].clear();
               for (int i=0; i < interpolationPoints; ++i) { 
                   RelPermValues[waterPhaseIndex][voigtIdx].push_back(RelPermWaterVsSaturation.evaluate(SatvaluesInterp[i])); 
                   RelPermValues[oilPhaseIndex][voigtIdx].push_back(RelPermOilVsSaturation.evaluate(SatvaluesInterp[i]));                
               } 
           } 
           

           // Now also overwrite Satvalues and fracFlowRatioPoints
           WaterSaturation.clear(); 
           WaterSaturation = SatvaluesInterp; 
 
           fracFlowRatioPoints.clear();
           fracFlowRatioPoints = fracFlowRatioPointsInterp;
       }
       
       
       for (unsigned int i=0; i < WaterSaturation.size(); ++i) {
           outwatertmp << showpoint << setw(14) << fracFlowRatioPoints[i];
           outwatertmp << showpoint << setw(14) << WaterSaturation[i];
           
           for (int voigtIdx=0; voigtIdx < tensorElementCount; ++voigtIdx) {
               outwatertmp << showpoint << setw(14)  << RelPermValues[waterPhaseIndex][voigtIdx][i];              
           }
           outwatertmp << endl;
           
           // Ignore further output if we have reached no water saturation
           if (WaterSaturation[i] == 0.0) break; /* maybe == on floats will fail at some point */ 
       }
       
       for (unsigned int i=0; i < WaterSaturation.size(); ++i) {
           outoiltmp << showpoint << setw(14) << fracFlowRatioPoints[i];
           outoiltmp << showpoint << setw(14) << WaterSaturation[i];     
           
           for (int voigtIdx=0; voigtIdx < tensorElementCount; ++voigtIdx) {
               outoiltmp << showpoint << setw(14) << RelPermValues[oilPhaseIndex][voigtIdx][i];              
           }
           outoiltmp << endl;
           
           // Ignore further output if we have reached no water saturation
           if (WaterSaturation[i] == 0.0) break; /* maybe == on floats will fail at some point */ 
       }
       
       cout << outheadtmp.str();
       cout << "# Relperm curve for water: " << endl;
       cout << outwatertmp.str();
       cout << "####################################################################" << endl;
       cout << "# Relperm curve for oil: " << endl;
       cout << outoiltmp.str();
       
       
       if (options["outputWater"] != "") {
           cout << "Writing (water) results to " << options["outputWater"] << endl;
           ofstream outfile;
           outfile.open(options["outputWater"].c_str(), ios::out | ios::trunc);
           outfile << outheadtmp.str();
           outfile << outwatertmp.str();
           outfile.close();      
       }
       
       if (options["outputOil"] != "") {
           cout << "Writing (oil) results to " << options["outputOil"] << endl;
           ofstream outfile;
           outfile.open(options["outputOil"].c_str(), ios::out | ios::trunc);
           outfile << outheadtmp.str();
           outfile << outoiltmp.str();
           outfile.close();      
       }
   }

#if USEMPI
   MPI_Finalize();
#endif

   return 0;
};
