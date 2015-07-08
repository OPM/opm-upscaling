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
   @file upscale_relperm_benchmark.C
   @brief Benchmark version of upscale_relperm.

   BENCHMARK VERSION

   This is a benchmark version of upscale_relperm, whose ordinary description is given below.
   The main machinery is the same, but with some changes:

   - Input data (grid file, rock file and reference solution) is not provided from command line,
     but is built in at compiler time by embedding hexadecimal (1 byte) input data files. See 
     README for further documentation.
   - Other command line options are not supported.
   - The construction of deck and stone data is changed due to change in input routine.
   - All checks on number of input from command line are removed (since this no longer is valid).
   - Some option defaults are changed:
       points = 20
       upscaleBothPhases = 0 (false)
       jFunctionCurve = 3
       outputprecision = 20
   - All output is surpressed, except for a start-up message and the final output in step 9.
   - A test for verification of computed solution is implemented in Step 9. The comparison is done
     within an absolute tolerance that can be changed below.
   - Output provided in Step 9 is changed to fit the benchmark suite

   There are two models available with different model sizes. It is important that the model size
   don't fit into the cache. The model type can be chosen below by changing the macro 'MODEL_TYPE'.
   The tolerance to be used when comparing results can also be changed below. To re-build simply
   run 'make' inside the root directory (opm-benchmarks).

   Assumptions:
   - Only one stone type
   - Isotropic input data
   - upscales only one phase
*/

/**
   Original description of upscale_relperm:

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

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <opm/upscaling/RelPermUtils.hpp>

#include <opm/core/utility/Units.hpp>

// Choose model:
//   - Small: MODEL_TYPE 1  (35751 active cells, ~5 MB)
#define MODEL_TYPE 1

// Define tolerance to be used when comparing results.
double tolerance = 1e-4;

// Include eclipse grid file and reference solution by embedding hexadecimal input data files.
#if MODEL_TYPE == 1
char model_name[] = "Small";
unsigned char eclipseInput[] = {
    #include <benchmarks/input/benchmark20_grid.grdecl.gz.hex>
};
unsigned char resultString[] = {
    #include <benchmarks/input/benchmark20_upscaled_relperm.out.gz.hex>
};
#else
#error The macro 'MODEL_TYPE' is invalid. Possible values are 1-1.
#endif

// Include rock file by embedding hexadecimal input data file.
unsigned char rockString[] = {
    #include <benchmarks/input/stonefile_benchmark.txt.gz.hex>
};


using namespace Opm;
using namespace std;
namespace io = boost::iostreams;

// Function for displaying a vector. Used for testing purposes
void dispVec(string name, vector<double> vec) {
    cout.precision(10);
    cout << name << ": ";
    for (vector<double>::iterator iter = vec.begin(); iter < vec.end(); ++iter) {
        cout << *iter << " ";
    }
    cout << "(size = " << vec.size() << endl;
}

static void usage() // Benchmark version
{
    cerr << "Usage: This is a benchmark version of upscale_relperm. This executable takes no input," << endl <<
            "the model and stone data is included at compiler time." << endl;
}

static void usageandexit() {
    usage();
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    exit(1);
}

static void inflate (const unsigned char* input, const int size, stringstream& output) {
    // read compressed data raw from .data segment into buffer
    stringstream compressed;
    std::copy (input, input+size, ostream_iterator <char> (compressed));

    // setup a filter which decompress the memory stream
    io::filtering_istream filter;
    filter.push (io::gzip_decompressor ());
    filter.push (compressed);

    // return the decompressed copy
    io::copy (filter, output);
}

int main(int varnum, char** vararg)
try
{

    // Variables used for timing/profiling:
    clock_t start, finish;
    double timeused = 0.0, timeused_tesselation = 0.0;
    double timeused_upscale_wallclock = 0.0;

    clock_t global_start = clock(); // Timing used for benchmarking

    /******************************************************************************
     * Step 1:
     * Process command line options
     */

    int mpi_rank = 0;
#ifdef HAVE_MPI
    int mpi_nodecount = 1;
    MPI_Init(&varnum, &vararg);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_nodecount);
#endif
    RelPermUpscaleHelper helper;
    helper.isMaster = (mpi_rank == 0);

    // Print start-up message
    if (helper.isMaster)
        cout << "Running benchmark version of upscale_relperm (model type: " << model_name << ") ..." << endl;

    // Suppress output in benchmark version (both cout and cerr):
    std::streambuf* cout_sbuf = std::cout.rdbuf(); // save original streambuf
    stringstream ss_null; // Stringstream to redirect all unwanted output
    std::cout.rdbuf(ss_null.rdbuf()); // redirect 'cout' to ss_null

    std::streambuf* cerr_sbuf = std::cerr.rdbuf();
    std::cerr.rdbuf(ss_null.rdbuf());

    if (varnum > 1) { // If arguments are supplied, show error ("upscale_relperm_benchmark" is the first (and only) "argument")
        usageandexit();
    }

    /*
      Populate options-map with default values
    */
    map<string,string> options =
        {{"bc",                       "f"}, // Fixed boundary conditions
        {"points",                   "20"}, // Number of saturation points (uniformly distributed within saturation endpoints) [benchmark]
        {"relPermCurve",              "2"}, // Which column in the rock types are upscaled
        {"upscaleBothPhases",     "false"}, // Whether to upscale for both phases in the same run. Default true. [benchmark]
        {"jFunctionCurve",            "3"}, // Which column in the rock type file is the J-function curve [benchmark]
        {"surfaceTension",           "11"}, // Surface tension given in dynes/cm
        {"output",                     ""}, // If this is set, output goes to screen and to this file.
        {"gravity",                 "0.0"}, // default is no gravitational effects
        {"waterDensity",            "1.0"}, // default density of water, only applicable to gravity
        {"oilDensity",              "0.6"}, // ditto
        {"interpolate",               "0"}, // default is not to interpolate
        {"maxpoints",              "1000"}, // maximal number of saturation points.
        {"outputprecision",          "20"}, // number of significant numbers to print [benchmark]
        {"maxPermContrast",         "1e7"}, // maximum allowed contrast in each single-phase computation
        {"minPerm",               "1e-12"}, // absolute minimum for allowed cell permeability
        {"maxPerm",              "100000"}, // maximal allowed cell permeability
        {"minPoro",              "0.0001"}, // this limit is necessary for pcmin/max computation
        {"saturationThreshold", "0.00001"}, // accuracy threshold for saturation, we ignore Pc values that
                                            // give so small contributions near endpoints.
        {"linsolver_tolerance",   "1e-12"}, // residual tolerance for linear solver
        {"linsolver_verbosity",       "0"}, // verbosity level for linear solver
        {"linsolver_type",            "3"}, // type of linear solver: 0 = ILU0/CG, 1 = AMG/CG, 2 KAMG/CG, 3 FAST_AMG/CG
        {"fluids",                   "ow"}, // wheater upscaling for oil/water (ow) or gas/oil (go)
        {"krowxswirr",               "-1"}, // relative permeability in x direction of oil in corresponding oil/water system
        {"krowyswirr",               "-1"}, // relative permeability in y direction of oil in corresponding oil/water system
        {"krowzswirr",               "-1"}, // relative permeability in z direction of oil in corresponding oil/water system
        {"doEclipseCheck",         "true"}, // Check if minimum relpermvalues in input are zero (specify critical saturations)
        {"critRelpermThresh",      "1e-6"}};// Threshold for setting minimum relperm to 0 (thus specify critical saturations)

    // Conversion factor, multiply mD numbers with this to get m² numbers
    const double milliDarcyToSqMetre =
        Opm::unit::convert::to(1.0*Opm::prefix::milli*Opm::unit::darcy,
                               Opm::unit::square(Opm::unit::meter));
    // Reference: http://www.spe.org/spe-site/spe/spe/papers/authors/Metric_Standard.pdf


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

    // Boolean set to true if input permeability in eclipse-file has diagonal anisotropy.
    // (full-tensor anisotropy will be ignored)
    helper.anisotropic_input = false;


    /* Check validity of boundary conditions chosen, and make booleans
       for boundary conditions, this allows more readable code later. */
    bool isFixed, isLinear, isPeriodic;
    if (options["bc"].substr(0,1) == "f") {
        isFixed = true; isLinear = false; isPeriodic = false;
        helper.boundaryCondition = SinglePhaseUpscaler::Fixed; // This refers to the mimetic namespace (Sintef)
        helper.tensorElementCount = 3; // Diagonal
    }
    else if (options["bc"].substr(0,1) == "l") {
        isLinear = true; isFixed = false; isPeriodic = false;
        helper.boundaryCondition = SinglePhaseUpscaler::Linear;
        helper.tensorElementCount = 9; // Full-tensor
    }
    else if (options["bc"].substr(0,1) == "p") {
        isPeriodic = true; isLinear = false; isFixed = false;
        helper.boundaryCondition = SinglePhaseUpscaler::Periodic;
        helper.tensorElementCount = 9; // Symmetric.
    }
    else {
        if (helper.isMaster) cout << "Invalid boundary condition. Only one of the letters f, l or p are allowed." << endl;
        usageandexit();
    }

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

    if (helper.isMaster) cout << "Parsing Eclipse file... ";
    flush(cout);   start = clock();

    // create the parser
    Opm::ParserPtr parser(new Opm::Parser);
    stringstream *gridstringstream(new stringstream(stringstream::in | stringstream::out));
    std::shared_ptr<std::istream> gridstream(gridstringstream);
    inflate (eclipseInput, sizeof (eclipseInput) / sizeof (eclipseInput[0]), *gridstringstream);
    Opm::DeckConstPtr deck = parser->parseStream(gridstream);

    finish = clock();   timeused = (double(finish)-double(start))/CLOCKS_PER_SEC;
    if (helper.isMaster) cout << " (" << timeused <<" secs)" << endl;

    const double maxPermContrast = atof(options["maxPermContrast"].c_str());
    const double minPerm = atof(options["minPerm"].c_str());
    const double maxPerm = atof(options["maxPerm"].c_str());
    const double minPoro = atof(options["minPoro"].c_str());
    const double saturationThreshold = atof(options["saturationThreshold"].c_str());
    start = clock();
    helper.sanityCheckInput(deck, minPerm, maxPerm, minPoro);

    Opm::DeckRecordConstPtr specgridRecord(deck->getKeyword("SPECGRID")->getRecord(0));
    int x_res = specgridRecord->getItem("NX")->getInt(0);
    int y_res = specgridRecord->getItem("NY")->getInt(0);
    int z_res = specgridRecord->getItem("NZ")->getInt(0);

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
    //const int outputprecision       = atoi(options["outputprecision"].c_str());


    // Benchmark version: (assumes only one phase to be upscaled, and only one stone type)
    stringstream stonestream(stringstream::in | stringstream::out);
    inflate (rockString, sizeof (rockString) / sizeof (rockString[0]), stonestream);
    vector<double> inputWaterSaturation;
    vector<double> inputRelPerm;
    vector<double> inputJfunction;

    string nextStoneLine;
    while (getline(stonestream, nextStoneLine)) {
        if (nextStoneLine[0] == '#') continue;
        double nextStoneValue;
        stringstream line(stringstream::in | stringstream::out);
        line << nextStoneLine;
        int colNr = 1;
        while (line >> nextStoneValue) {
            if (colNr == 1) inputWaterSaturation.push_back(nextStoneValue);
            else if (colNr == relPermCurve) inputRelPerm.push_back(nextStoneValue);
            else if (colNr == jFunctionCurve) inputJfunction.push_back(nextStoneValue);
            colNr++;
        }
    }

    MonotCubicInterpolator Jtmp(inputWaterSaturation, inputJfunction);
    for (int i=0; i < stone_types; ++i) {
        helper.Krfunctions[0][0].push_back(MonotCubicInterpolator(inputWaterSaturation, inputRelPerm));
        // Invert J-function
        if (Jtmp.isStrictlyMonotone()) {
            helper.InvJfunctions.push_back(MonotCubicInterpolator(Jtmp.get_fVector(), Jtmp.get_xVector()));
        }
        else {
            if (helper.isMaster) cout << "Error: J-function curve not strictly monotone" << endl;
            usageandexit();
        }
        JfunctionNames.push_back("stonefile_benchmark.txt");
    }

    // Check if input relperm curves satisfy Eclipse requirement of specifying critical saturations
    helper.doEclipseCheck = (options["doEclipseCheck"] == "true");
    helper.critRelpThresh = atof(options["critRelpermThresh"].c_str());
    int numberofrockstocheck = 1;
    // if (varnum == rockfileindex + stone_types) numberofrockstocheck = stone_types;
    // else numberofrockstocheck = 1;
    if (helper.doEclipseCheck) {
        for (int i=0 ; i < numberofrockstocheck; ++i) {
            if (helper.anisotropic_input) {
                double minrelpx = Krxfunctions[i].getMinimumF().second;
                double minrelpy = Kryfunctions[i].getMinimumF().second;
                double minrelpz = Krzfunctions[i].getMinimumF().second;
                if (minrelpx == 0) ; // Do nothing
                else if (minrelpx < helper.critRelpThresh) {
                    // set to 0
                    vector<double> svec, kvec;
                    svec = Krxfunctions[i].get_xVector();
                    kvec = Krxfunctions[i].get_fVector();
                    if (kvec[0] < helper.critRelpThresh) {
                        kvec[0] = 0.0;
                    }
                    else if (kvec[kvec.size()-1] < helper.critRelpThresh) {
                        kvec[kvec.size()-1] = 0.0;
                    }
                    Krxfunctions[i] = MonotCubicInterpolator(svec, kvec);
                }
                else {
                    // Error message
                    cerr << "Relperm curve for rock " << i << " does not specify critical saturation." << endl
                         << "Minimum relperm value is " << minrelpx << ", critRelpermThresh is " << helper.critRelpThresh << endl;
                    usageandexit();
                }
                //
                if (minrelpy == 0) ; // Do nothing
                else if (minrelpy < helper.critRelpThresh) {
                    // set to 0
                    vector<double> svec, kvec;
                    svec = Kryfunctions[i].get_xVector();
                    kvec = Kryfunctions[i].get_fVector();
                    if (kvec[0] < helper.critRelpThresh) {
                        kvec[0] = 0.0;
                    }
                    else if (kvec[kvec.size()-1] < helper.critRelpThresh) {
                        kvec[kvec.size()-1] = 0.0;
                    }
                    Kryfunctions[i] = MonotCubicInterpolator(svec, kvec);
                }
                else {
                    // Error message
                    cerr << "Relperm curve for rock " << i << " does not specify critical saturation." << endl
                         << "Minimum relperm value is " << minrelpy << ", critRelpermThresh is " << helper.critRelpThresh << endl;
                    usageandexit();
                }
                //
                if (minrelpz == 0) ; // Do nothing
                else if (minrelpz < helper.critRelpThresh) {
                    // set to 0
                    vector<double> svec, kvec;
                    svec = Krzfunctions[i].get_xVector();
                    kvec = Krzfunctions[i].get_fVector();
                    if (kvec[0] < helper.critRelpThresh) {
                        kvec[0] = 0.0;
                    }
                    else if (kvec[kvec.size()-1] < helper.critRelpThresh) {
                        kvec[kvec.size()-1] = 0.0;
                    }
                    Krzfunctions[i] = MonotCubicInterpolator(svec, kvec);
                }
                else {
                    // Error message
                    cerr << "Relperm curve for rock " << i << " does not specify critical saturation." << endl
                         << "Minimum relperm value is " << minrelpz << ", critRelpermThresh is " << helper.critRelpThresh << endl;
                    usageandexit();
                }
                //
                if (helper.upscaleBothPhases) {
                    minrelpx = Krxfunctions2[i].getMinimumF().second;
                    minrelpy = Kryfunctions2[i].getMinimumF().second;
                    minrelpz = Krzfunctions2[i].getMinimumF().second;
                    if (minrelpx == 0) ; // Do nothing
                    else if (minrelpx < helper.critRelpThresh) {
                        // set to 0
                        vector<double> svec, kvec;
                        svec = Krxfunctions2[i].get_xVector();
                        kvec = Krxfunctions2[i].get_fVector();
                        if (kvec[0] < helper.critRelpThresh) {
                            kvec[0] = 0.0;
                        }
                        else if (kvec[kvec.size()-1] < helper.critRelpThresh) {
                            kvec[kvec.size()-1] = 0.0;
                        }
                        Krxfunctions2[i] = MonotCubicInterpolator(svec, kvec);
                    }
                    else {
                        // Error message
                        cerr << "Relperm curve for rock " << i << " does not specify critical saturation." << endl
                             << "Minimum relperm value is " << minrelpx << ", critRelpermThresh is " << helper.critRelpThresh << endl;
                        usageandexit();
                    }
                    //
                    if (minrelpy == 0) ; // Do nothing
                    else if (minrelpy < helper.critRelpThresh) {
                        // set to 0
                        vector<double> svec, kvec;
                        svec = Kryfunctions2[i].get_xVector();
                        kvec = Kryfunctions2[i].get_fVector();
                        if (kvec[0] < helper.critRelpThresh) {
                            kvec[0] = 0.0;
                        }
                        else if (kvec[kvec.size()-1] < helper.critRelpThresh) {
                            kvec[kvec.size()-1] = 0.0;
                        }
                        Kryfunctions2[i] = MonotCubicInterpolator(svec, kvec);
                    }
                    else {
                        // Error message
                        cerr << "Relperm curve for rock " << i << " does not specify critical saturation." << endl
                             << "Minimum relperm value is " << minrelpy << ", critRelpermThresh is " << helper.critRelpThresh << endl;
                        usageandexit();
                    }
                    //
                    if (minrelpz == 0) ; // Do nothing
                    else if (minrelpz < helper.critRelpThresh) {
                        // set to 0
                        vector<double> svec, kvec;
                        svec = Krzfunctions2[i].get_xVector();
                        kvec = Krzfunctions2[i].get_fVector();
                        if (kvec[0] < helper.critRelpThresh) {
                            kvec[0] = 0.0;
                        }
                        else if (kvec[kvec.size()-1] < helper.critRelpThresh) {
                            kvec[kvec.size()-1] = 0.0;
                        }
                        Krzfunctions2[i] = MonotCubicInterpolator(svec, kvec);
                    }
                    else {
                        // Error message
                        cerr << "Relperm curve for rock " << i << " does not specify critical saturation." << endl
                             << "Minimum relperm value is " << minrelpz << ", critRelpermThresh is " << helper.critRelpThresh << endl;
                        usageandexit();
                    }
                    //
                }
            }
            else {
                double minrelp = helper.Krfunctions[0][i].getMinimumF().second;
                if (minrelp == 0) ; // Do nothing
                else if (minrelp < helper.critRelpThresh) {
                    // set to 0
                    vector<double> svec, kvec;
                    svec = helper.Krfunctions[0][i].get_xVector();
                    kvec = helper.Krfunctions[0][i].get_fVector();
                    if (kvec[0] < helper.critRelpThresh) {
                        kvec[0] = 0.0;
                    }
                    else if (kvec[kvec.size()-1] < helper.critRelpThresh) {
                        kvec[kvec.size()-1] = 0.0;
                    }
                    helper.Krfunctions[0][i] = MonotCubicInterpolator(svec, kvec);
                }
                else {
                    // Error message
                    cerr << "Relperm curve for rock " << i << " does not specify critical saturation." << endl
                         << "Minimum relperm value is " << minrelp << ", critRelpermThresh is " << helper.critRelpThresh << endl;
                    usageandexit();
                }
                if (helper.upscaleBothPhases) {
                    minrelp = helper.Krfunctions[1][i].getMinimumF().second;
                    if (minrelp == 0) ;
                    else if (minrelp < helper.critRelpThresh) {
                        // set to 0
                        vector<double> svec, kvec;
                        svec = helper.Krfunctions[1][i].get_xVector();
                        kvec = helper.Krfunctions[1][i].get_fVector();
                        if (kvec[0] < helper.critRelpThresh) kvec[0] = 0.0;
                        else if (kvec[kvec.size()-1] < helper.critRelpThresh) kvec[kvec.size()-1] = 0.0;
                        helper.Krfunctions[1][i] = MonotCubicInterpolator(svec, kvec);
                    }
                    else {
                        // Error message
                        cerr << "Relperm curve for rock " << i << " does not specify critical saturation."
                             << "Minimum relperm value is " << minrelp << ", critRelpermThresh is " << helper.critRelpThresh << endl;
                        usageandexit();
                    }
                }
            }
        }
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

    if (helper.isMaster) cout << "Tesselating grid... ";
    flush(cout);   start = clock();
    double linsolver_tolerance = atof(options["linsolver_tolerance"].c_str());
    int linsolver_verbosity = atoi(options["linsolver_verbosity"].c_str());
    int linsolver_type = atoi(options["linsolver_type"].c_str());
    bool twodim_hack = false;
    helper.upscaler.init(deck, helper.boundaryCondition,
                         Opm::unit::convert::from(minPerm, Opm::prefix::milli*Opm::unit::darcy),
                         linsolver_tolerance, linsolver_verbosity, linsolver_type, twodim_hack);

    finish = clock();   timeused_tesselation = (double(finish)-double(start))/CLOCKS_PER_SEC;
    if (helper.isMaster) cout << " (" << timeused_tesselation <<" secs)" << endl;

    vector<double> dP;
    double dPmin = +DBL_MAX;
    double dPmax = -DBL_MAX;

    /* If gravity is to be included, calculate z-values of every cell: */
    if (includeGravity) {
        // height of model is calculated as the average of the z-values at the top layer
        // This calculation makes assumption on the indexing of cells in the grid, going from bottom to top.
        double modelHeight = 0;
        for (unsigned int zIdx = (4 * x_res * y_res * (2*z_res-1)); zIdx < helper.zcorns.size(); ++zIdx) {
            modelHeight += helper.zcorns[zIdx] / (4*x_res*y_res);
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
            int horIdx = (cellIdx+1) - int(std::floor(((double)(cellIdx+1))/((double)(x_res*y_res))))*x_res*y_res; // index in the corresponding horizon
            if (horIdx == 0) {
                horIdx = x_res*y_res;
            }
            i = horIdx - int(std::floor(((double)horIdx)/((double)x_res)))*x_res;
            if (i == 0) {
                i = x_res;
            }
            j = (horIdx-i)/x_res+1;
            k = ((cellIdx+1)-x_res*(j-1)-1)/(x_res*y_res)+1;
            int zBegin = 8*x_res*y_res*(k-1); // indices of Z-values of bottom
            int level2 = 4*x_res*y_res; // number of z-values in one horizon
            zIndices[0] = zBegin + 4*x_res*(j-1)+2*i-1;
            zIndices[1] = zBegin + 4*x_res*(j-1)+2*i;
            zIndices[2] = zBegin + 2*x_res*(2*j-1)+2*i;
            zIndices[3] = zBegin + 2*x_res*(2*j-1)+2*i-1;
            zIndices[4] = zBegin + level2 + 4*x_res*(j-1)+2*i-1;
            zIndices[5] = zBegin + level2 + 4*x_res*(j-1)+2*i;
            zIndices[6] = zBegin + level2 + 2*x_res*(2*j-1)+2*i;
            zIndices[7] = zBegin + level2 + 2*x_res*(2*j-1)+2*i-1;

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
            double WaterSaturationCell = 0.0;
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
                    WaterSaturationCell
                        = helper.InvJfunctions[int(helper.satnums[cell_idx])-1].evaluate(Jvalue);
                }
                else { // anisotropic_input, then we do not do J-function-scaling
                    WaterSaturationCell = helper.SwPcfunctions[int(helper.satnums[cell_idx])-1].evaluate(PtestvalueCell);
                    //cout << Ptestvalue << "\t" <<  helper.WaterSaturationCell << endl;
                }
            }
            waterVolume += WaterSaturationCell  * helper.cellPoreVolumes[cell_idx];
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


    clock_t start_upscaling = clock();

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
        /*if (helper.isMaster) {
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
                        double WaterSaturationCell
                            = helper.InvJfunctions[int(helper.satnums[cell_idx])-1].evaluate(Jvalue);
                        waterVolumeLF += WaterSaturationCell * helper.cellPoreVolumes[cell_idx];

                        // Compute cell relative permeability. We use a lower cutoff-value as we
                        // easily divide by zero here.  When water saturation is
                        // zero, we get 'inf', which is circumvented by the cutoff value.
                        cellPhasePerm =
                            helper.Krfunctions[0][0][int(helper.satnums[cell_idx])-1].evaluate(WaterSaturationCell) *
                            helper.perms[0][cell_idx];
                        if (helper.upscaleBothPhases) {
                            cellPhase2Perm =
                                helper.Krfunctions[0][1][int(helper.satnums[cell_idx])-1].evaluate(WaterSaturationCell) *
                                helper.perms[0][cell_idx];
                        }
                    }
                    else {
                        double WaterSaturationCell = helper.SwPcfunctions[int(helper.satnums[cell_idx])-1].evaluate(PtestvalueCell);
                        //cout << PtestvalueCell << "\t" << helper.WaterSaturationCell << endl;
                        waterVolumeLF += WaterSaturationCell * helper.cellPoreVolumes[cell_idx];

                        cellPhasePermDiag[0] = helper.Krfunctions[0][0][int(helper.satnums[cell_idx])-1].evaluate(WaterSaturationCell) *
                            helper.perms[0][cell_idx];
                        cellPhasePermDiag[1] = helper.Krfunctions[1][0][int(helper.satnums[cell_idx])-1].evaluate(WaterSaturationCell) *
                            helper.perms[1][cell_idx];
                        cellPhasePermDiag[2] = helper.Krfunctions[2][0][int(helper.satnums[cell_idx])-1].evaluate(WaterSaturationCell) *
                            helper.perms[2][cell_idx];
                        if (helper.upscaleBothPhases) {
                            cellPhase2PermDiag[0] = helper.Krfunctions[0][1][int(helper.satnums[cell_idx])-1].evaluate(WaterSaturationCell) *
                                helper.perms[0][cell_idx];
                            cellPhase2PermDiag[1] = helper.Krfunctions[1][1][int(helper.satnums[cell_idx])-1].evaluate(WaterSaturationCell) *
                                helper.perms[1][cell_idx];
                            cellPhase2PermDiag[2] = helper.Krfunctions[2][1][int(helper.satnums[cell_idx])-1].evaluate(WaterSaturationCell) *
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
    double avg_upscaling_time_pr_point = timeused_total/(double)points;

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
     *  Step 9 - Benchmark version
     *
     * Assuming fixed BCs, one phase. Everything done only by master node.
     *
     * a: verify results
     * b: output
     */

     vector<double> Pvalues = helper.pressurePoints;
     // Multiply all pressures with the surface tension (potentially) supplied
     // at the command line. This multiplication has been postponed to here
     // to avoid division by zero and to avoid special handling of negative
     // capillary pressure in the code above.
     std::transform(Pvalues.begin(), Pvalues.end(), Pvalues.begin(),
                    std::bind1st(std::multiplies<double>(), surfaceTension));

    /*
     * Step 9a: Verify results
     */

    if (helper.isMaster) {

        // Define reference solutions to compare with
        vector<double> PvaluesReference;
        vector<double> WaterSaturationReference;
        vector<vector <double> > RelPermValuesReference;
        bool pressuresEqual = true;
        bool saturationsEqual = true;
        bool relpermsEqual = true;
        bool referenceResultsMatch = true;

        vector<double> tempVec;
        RelPermValuesReference.push_back(tempVec);
        RelPermValuesReference.push_back(tempVec);
        RelPermValuesReference.push_back(tempVec);

        // Read reference solution
        stringstream referencestream(stringstream::in | stringstream::out);
        inflate (resultString, sizeof (resultString) / sizeof (resultString[0]), referencestream);
        string nextReferenceLine;
        while (getline(referencestream, nextReferenceLine)) {
            if (nextReferenceLine[0] == '#') continue;
            stringstream line(stringstream::in | stringstream::out);
            double nextReferenceValue;
            line << nextReferenceLine;
            line >> nextReferenceValue;
            PvaluesReference.push_back(nextReferenceValue);
            line >> nextReferenceValue;
            WaterSaturationReference.push_back(nextReferenceValue);
            line >> nextReferenceValue;
            RelPermValuesReference[0].push_back(nextReferenceValue);
            line >> nextReferenceValue;
            RelPermValuesReference[1].push_back(nextReferenceValue);
            line >> nextReferenceValue;
            RelPermValuesReference[2].push_back(nextReferenceValue);
        }

        // Check that number of upscaled points matches
        size_t npoints = helper.points;
        if (WaterSaturationReference.size() != npoints) referenceResultsMatch = false;
        if (PvaluesReference.size() != npoints) referenceResultsMatch = false;
        if (RelPermValuesReference[0].size() != npoints) referenceResultsMatch = false;
        if (RelPermValuesReference[1].size() != npoints) referenceResultsMatch = false;
        if (RelPermValuesReference[2].size() != npoints) referenceResultsMatch = false;

        // Verify results
        if (referenceResultsMatch) {
            for (int i=0; i<helper.points; ++i) {
                if (fabs(WaterSaturationReference[i] - helper.WaterSaturation[i]) > tolerance) {
                    saturationsEqual = false;
                    break;
                }
                if (fabs(PvaluesReference[i] - Pvalues[i]) > tolerance*100) {
                    pressuresEqual = false;
                    break;
                }
            }
            for (int voigtIdx=0; voigtIdx<helper.tensorElementCount; ++voigtIdx) {
                for (int i=0; i<helper.points; ++i) {
                    if (fabs(RelPermValuesReference[voigtIdx][i] - RelPermValues[voigtIdx][i]) > tolerance) {
                        relpermsEqual = false;
                        break;
                    }
                }
            }
        }


        /*
         * Step 9b: Output run-time results
         */

        // Allow output again
        std::cout.rdbuf(cout_sbuf); // restore the original stream buffer
        std::cerr.rdbuf(cerr_sbuf);

        clock_t global_finish = clock();
        double processing_time = (double(start_upscaling) - double(global_start))/CLOCKS_PER_SEC;
        double upscaling_time = (double(global_finish) - double(start_upscaling))/CLOCKS_PER_SEC;
        double benchmark_time = (double(global_finish) - double(global_start))/CLOCKS_PER_SEC;
        double benchmark_time_min = floor(benchmark_time/60.0);
        double benchmark_time_sec = benchmark_time - benchmark_time_min*60;
        stringstream outputtmp;
        string versiondate = "17.07.2012";
        string dashed_line = "----------------------------------------------------------------------\n";

        outputtmp << endl;
        outputtmp << dashed_line;
#ifdef HAVE_MPI
        outputtmp << "upscale_relperm, MPI version (" << versiondate << ")" << endl;
#else
        outputtmp << "upscale_relperm, no-MPI version (" << versiondate << ")" << endl;
#endif
        outputtmp << "Part of the OPM project, http://www.opm-project.org\n";

        // Calculate approx model size
        int nCellsTotal = x_res*y_res*z_res;
        int model_size = (8*nCellsTotal + 2*nCellsTotal + (x_res+1)*(y_res+1)*2)*sizeof(double) + 2*nCellsTotal*sizeof(int);

        outputtmp << dashed_line;
        outputtmp << "Model type :      " << model_name << endl;
        outputtmp << "Active cells:     " << helper.tesselatedCells << endl;
        outputtmp << "Model size:       ~" << model_size/1000000 << " MB" << endl;
        outputtmp << "Upscaling points: " << helper.points << endl;
        // outputtmp << "Stone file:       " << JfunctionNames[1] << endl;
        outputtmp << "Different model sizes are available, change model type to increase the size.\n"
            "Model type can be changed by editing macro 'MODEL_TYPE' in source file." << endl;

#ifdef HAVE_MPI
        outputtmp << dashed_line;
        outputtmp << "MPI-nodes: " << mpi_nodecount << endl;
        double speedup = (avg_upscaling_time_pr_point * (points + 1) + timeused_tesselation)/(timeused_upscale_wallclock + avg_upscaling_time_pr_point + timeused_tesselation);
        outputtmp << "Speedup:   " << speedup << endl;
#endif

        outputtmp << dashed_line << "Verification of results:" << endl;
        if (!referenceResultsMatch) {
            outputtmp << "The number of upscaled points does not match the number of \n"
                      << "upscaled points in reference solution. Validation not possible.\n";
        }
        else if (pressuresEqual && saturationsEqual && relpermsEqual) {
            outputtmp << "Computed results are verified to be equal to reference\n"
                      << "solution within an absolute tolerance of " << tolerance << ".\n"
                      << "The tolerance can be changed in the source file.\n";
        }
        else {
            outputtmp << "Computed results are not equal to reference solution\n"
                      << "within an absolute tolerance of " << tolerance << ".\n"
                      << "The tolerance can be changed in the source file.\n";
        }

        outputtmp << dashed_line;
        outputtmp << "Wallclock timing:\n";
        outputtmp << "Input- and grid processing: " << processing_time << " sec" << endl;
        outputtmp << "Upscaling:                  " << upscaling_time << " sec" << endl;
        outputtmp << "Total wallclock time:       " << benchmark_time << " sec";

        if (benchmark_time > 60.0) {
            outputtmp << "  (" << int(benchmark_time_min) << " min " << benchmark_time_sec << " sec)\n";
        }
        else {
            outputtmp << endl;
        }
        outputtmp << dashed_line;

        cout << outputtmp.str();

    }


#if HAVE_MPI
    MPI_Finalize();
#endif

    return 0;
}
catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    usageandexit();
}
