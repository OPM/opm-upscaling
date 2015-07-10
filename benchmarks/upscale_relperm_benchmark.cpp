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
    // Print start-up message
    if (mpi_rank == 0)
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
        {{"bc",                          "f"}, // Fixed boundary conditions
        {"points",                      "20"}, // Number of saturation points (uniformly distributed within saturation endpoints) [benchmark]
        {"relPermCurve",                 "2"}, // Which column in the rock types are upscaled
        {"upscaleBothPhases",        "false"}, // Whether to upscale for both phases in the same run. Default true. [benchmark]
        {"jFunctionCurve",               "3"}, // Which column in the rock type file is the J-function curve [benchmark]
        {"surfaceTension",              "11"}, // Surface tension given in dynes/cm
        {"output",                        ""}, // If this is set, output goes to screen and to this file.
        {"gravity",                    "0.0"}, // default is no gravitational effects
        {"waterDensity",               "1.0"}, // default density of water, only applicable to gravity
        {"oilDensity",                 "0.6"}, // ditto
        {"interpolate",                  "0"}, // default is not to interpolate
        {"maxpoints",                 "1000"}, // maximal number of saturation points.
        {"outputprecision",             "20"}, // number of significant numbers to print [benchmark]
        {"maxPermContrast",            "1e7"}, // maximum allowed contrast in each single-phase computation
        {"minPerm",                  "1e-12"}, // absolute minimum for allowed cell permeability
        {"maxPerm",                 "100000"}, // maximal allowed cell permeability
        {"minPoro",                 "0.0001"}, // this limit is necessary for pcmin/max computation
        {"saturationThreshold",    "0.00001"}, // accuracy threshold for saturation, we ignore Pc values that
                                               // give so small contributions near endpoints.
        {"linsolver_tolerance",      "1e-12"}, // residual tolerance for linear solver
        {"linsolver_verbosity",          "0"}, // verbosity level for linear solver
        {"linsolver_type",               "3"}, // type of linear solver: 0 = ILU0/CG, 1 = AMG/CG, 2 KAMG/CG, 3 FAST_AMG/CG
        {"linsolver_max_iterations",     "0"}, // Maximum number of iterations allow, specify 0 for default
        {"linsolver_prolongate_factor","1.0"}, // Prolongation factor in AMG
        {"linsolver_smooth_steps",       "1"}, // Number of smoothing steps in AMG
        {"fluids",                      "ow"}, // wheater upscaling for oil/water (ow) or gas/oil (go)
        {"krowxswirr",                  "-1"}, // relative permeability in x direction of oil in corresponding oil/water system
        {"krowyswirr",                  "-1"}, // relative permeability in y direction of oil in corresponding oil/water system
        {"krowzswirr",                  "-1"}, // relative permeability in z direction of oil in corresponding oil/water system
        {"doEclipseCheck",            "true"}, // Check if minimum relpermvalues in input are zero (specify critical saturations)
        {"critRelpermThresh",         "1e-6"}};// Threshold for setting minimum relperm to 0 (thus specify critical saturations)

    RelPermUpscaleHelper helper(mpi_rank, options);
    helper.setupBoundaryConditions(options);

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

    const double minPerm = atof(options["minPerm"].c_str());
    const double maxPerm = atof(options["maxPerm"].c_str());
    const double minPoro = atof(options["minPoro"].c_str());
    start = clock();
    helper.sanityCheckInput(deck, minPerm, maxPerm, minPoro);

    Opm::DeckRecordConstPtr specgridRecord(deck->getKeyword("SPECGRID")->getRecord(0));
    std::array<int,3> res;
    res[0] = specgridRecord->getItem("NX")->getInt(0);
    res[1] = specgridRecord->getItem("NY")->getInt(0);
    res[2] = specgridRecord->getItem("NZ")->getInt(0);

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
        else
            throw std::runtime_error("Error: J-function curve not strictly monotone");

        JfunctionNames.push_back("stonefile_benchmark.txt");
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

    timeused_tesselation = helper.tesselateGrid(deck, options);

    /* If gravity is to be included, calculate z-values of every cell: */
    if (includeGravity)
        helper.calculateCellPressureGradients(res, options);

    /******************************************************************************
     * Step 5:
     * Go through each cell and calculate the minimum and
     * maximum capillary pressure possible in the cell, given poro,
     * perm and the J-function for the cell.  This depends on the
     * J-function in that they represent all possible saturations,
     * ie. we do not want to extrapolate the J-functions (but we might
     * have to do that later in the computations).
     */
    helper.calculateMinMaxCapillaryPressure(options);
    const std::vector<int>& ecl_idx = helper.upscaler.grid().globalCell();

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

    helper.upscaleCapillaryPressure(options);

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

    double timeused_upscale_wallclock, avg_upscaling_time_pr_point;
    std::tie(timeused_upscale_wallclock, avg_upscaling_time_pr_point) =
                helper.upscalePermeability(options, mpi_rank);

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
                    if (fabs(RelPermValuesReference[voigtIdx][i] - RelPermValues[0][voigtIdx][i]) > tolerance) {
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
        int nCellsTotal = res[0]*res[1]*res[2];
        int model_size = (8*nCellsTotal + 2*nCellsTotal + (res[0]+1)*(res[1]+1)*2)*sizeof(double) + 2*nCellsTotal*sizeof(int);

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
