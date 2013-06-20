// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:

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

/** @file upscale_perm.C
 *  @brief Upscales permeability
 *  
 *  Upscales permeability for a given Eclipse-model. 
 *  Upscaling is performed for periodic, linear and
 *  fixed boundary conditions and output is ASCII to standard out.
 *  
 *  Input is Eclipse grid format specifying the corner-point
 *  grid (must be of shoebox-shape, but this condition is slightly relaxed
 *  on top and bottom surfaces).
 * 
 *  The input eclipse file must specify the permeability properties
 *  for each cell.
 * 
 *  If only PERMX is supplied, isotropic permeability is assumed in each cell
 *
 *  If PERMX, PERMY and PERMZ are supplied, diagonal anisotropic permeability
 *  is assumed in each cell.
 *
 *  If PERMXX, PERMYY, PERMZZ, PERMXY, PERMYZ and PERMZX are supplied,
 *  full-tensor (symmetric) anisotropic permeability is assumed in each cell.
 *
 */
#include <config.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <sys/utsname.h>

#include <opm/upscaling/SinglePhaseUpscaler.hpp>
#include <dune/common/mpihelper.hh>

using namespace std;


void usage() {
    cout << endl <<
        "Usage: upscale_perm <options> <eclipsefile>" << endl <<
        "where the options are:" << endl <<
        "-output <string>  -- filename for where to write upscaled values." << endl <<
        "                     If not supplied, output will only go to " << endl <<
        "                     the terminal (standard out)." << endl <<
        "-bc <string>      -- which boundary conditions to compute for. " << endl <<
        "                     <string> may contain a combination of the" << endl <<
        "                     letters lfp to compute linear (l), fixed (f) " << endl << 
        "                     and periodic (p) boundary conditions." << endl <<
        "                     Default: f (fixed boundary conditions)" << endl <<
        "-minPerm <float>  -- Minimum floating point value allowed for" << endl <<
        "                     permeability. If zero, the problem is singular" << endl <<
        "                     Default 1e-9. Unit Millidarcy." << endl;
}

/**
   @brief Upscales permeability
   
   @param varnum Number of input arguments
   @param vararg Input arguments
   @return int
*/
int upscale(int varnum, char** vararg) {
    Dune::MPIHelper& mpi=Dune::MPIHelper::instance(varnum, vararg);
    mpi.rank();
    if (varnum ==  1) { // If no arguments supplied ("upscale_perm" is the first argument)
        cout << "Error: No eclipsefile provided" << endl;
        usage();
        exit(1);
    } 
    map<string,string> options;
    options.insert(make_pair("output", "")); // If this is set, output goes to screen and to this file 
    options.insert(make_pair("bc",     "f")); // Fixed boundary conditions are default    
    options.insert(make_pair("minPerm", "1e-9")); // Minimum allowable permeability value (for diagonal tensor entries)
    
    options.insert(make_pair("linsolver_tolerance", "1e-8"));  // residual tolerance for linear solver
    options.insert(make_pair("linsolver_verbosity", "0"));     // verbosity level for linear solver
    options.insert(make_pair("linsolver_max_iterations", "0"));         // Maximum number of iterations allow, specify 0 for default
    options.insert(make_pair("linsolver_prolongate_factor", "1.6")); // Factor to scale the prolongate coarse grid correction
    options.insert(make_pair("linsolver_type",      "1"));     // type of linear solver: 0 = ILU/BiCGStab, 1 = AMG/CG
    options.insert(make_pair("linsolver_smooth_steps", "2")); // Number of pre and postsmoothing steps for AMG

    // Parse options from command line
    int eclipseindex = 1; // Index for the eclipsefile in the command line options

    for (int argidx = 1; argidx < varnum - 1; argidx += 2)   {
        string searchfor = string(vararg[argidx]).substr(1); //Chop off leading '-'
        if (string(vararg[argidx]).substr(0,1) != "-") {
            cout << "Error: Options must come before eclipse file" << endl;
            usage();
            exit(1);
        }
        // Check if it is a match
        if (options.count(searchfor) == 1) {
            options[searchfor] = string(vararg[argidx+1]);
            cout << "Parsed command line option: " << searchfor << " := " << vararg[argidx+1] << endl;
            eclipseindex += 2;
        }
    }               

    const char* ECLIPSEFILENAME(vararg[eclipseindex]);

   
    // Test if filename exists and is readable
    ifstream eclipsefile(ECLIPSEFILENAME, ios::in);
    if (eclipsefile.fail()) {
        cerr << "Error: Filename " << ECLIPSEFILENAME << " not found or not readable." << endl;
        usage();
        exit(1);
    }
    eclipsefile.close();
    
    // Check validity of boundary conditions chosen, and make booleans 
    // for boundary conditions, this allows more readable code later
    bool isFixed, isLinear, isPeriodic;
    isFixed = isLinear = isPeriodic = false;
    
    // Read in default or user-specified boundary conditions:
    string boundcond(options["bc"]);
    
    // Length of string must be between 1 and 3:
    if (boundcond.length()>= 1 && boundcond.length() <= 3) {
        if (boundcond.find(string("p")) < 3) {
            isPeriodic = true;
        }
        if (boundcond.find(string("f")) < 3) {
            isFixed = true;
        }
        if (boundcond.find(string("l")) < 3) {
            isLinear = true;
        }
        
        // If no boundary conditions are set now, issue error:
        if (!isFixed && !isLinear && !isPeriodic) {
            cerr << "Error: No boundary conditions specified: " << boundcond << endl;
            usage();
            exit(1);                                                                        
        }
    }
    else {
        cerr << "Error: Syntax error in specifying boundary conditions: " << boundcond << endl;
        usage();
        exit(1);
    }

   
    // Variables for timing/profiling
    clock_t start, finish;
    double timeused = 0; // reusable variable
    double timeused_periodic_tesselation = 0, timeused_nonperiodic_tesselation = 0;
    double timeused_periodic = 0,  timeused_fixed = 0, timeused_linear = 0;
        
    cout << endl;
    
    // Storage for upscaled results:
    using Opm::SinglePhaseUpscaler;
    typedef SinglePhaseUpscaler::permtensor_t Matrix;
    Matrix Kfixed, Klinear, Kperiodic;
   

    /***********************************************************************
     * Load geometry and data from Eclipse file
     */
    cout << "Parsing Eclipse file <" << ECLIPSEFILENAME << "> ... ";
    flush(cout);   start = clock();
    Opm::EclipseGridParser * eclParser_p;
    try {
        eclParser_p = new Opm::EclipseGridParser(ECLIPSEFILENAME);
    }
    catch (...) {
        cout << "Error: Filename " << ECLIPSEFILENAME << " does not look like an eclipse grid file." << endl;
        usage();
        exit(1);
    }
    Opm::EclipseGridParser& eclParser = *eclParser_p;

    finish = clock();   timeused = (double(finish)-double(start))/CLOCKS_PER_SEC;
    cout << " (" << timeused <<" secs)" << endl;
 
    // Check that we have the information we need from the eclipse file, we will check PERM-fields later
    if (! (eclParser.hasField("SPECGRID") && eclParser.hasField("COORD") && eclParser.hasField("ZCORN"))) {  
        cerr << "Error: Did not find SPECGRID, COORD and ZCORN in Eclipse file " << ECLIPSEFILENAME << endl;  
        usage();  
        exit(1);  
    }


    /*****************************************************************
     * Tesselate grid 
     * 
     * Possibly twice because, the grid must be massaged slightly
     * (crop top and bottom) for periodic boundary conditions. These
     * modifications ruin the computations for linear and fixed
     * boundary conditions, so we must tesselate twice, thus use more
     * memory and more time for processing.
     */


    double ztol = 0.0; 
    double linsolver_tolerance = atof(options["linsolver_tolerance"].c_str());
    int linsolver_verbosity = atoi(options["linsolver_verbosity"].c_str());
    int linsolver_type = atoi(options["linsolver_type"].c_str());
    int linsolver_maxit = atoi(options["linsolver_max_iterations"].c_str());
    int smooth_steps = atoi(options["linsolver_smooth_steps"].c_str());
    double linsolver_prolongate_factor = atof(options["linsolver_prolongate_factor"].c_str());
    bool twodim_hack = false;

    SinglePhaseUpscaler upscaler_nonperiodic;
    SinglePhaseUpscaler upscaler_periodic;

    const double minPerm = Opm::unit::convert::from(atof(options["minPerm"].c_str()),
                                                    Opm::prefix::milli*Opm::unit::darcy);

    if (isFixed || isLinear)  {
        cout << "Tesselating non-periodic grid ...";
        start = clock();
        upscaler_nonperiodic.init(eclParser, 
                                  isFixed ? SinglePhaseUpscaler::Fixed : SinglePhaseUpscaler::Linear,
                                  minPerm, ztol,  linsolver_tolerance, linsolver_verbosity, linsolver_type, 
                                  twodim_hack, linsolver_maxit, linsolver_prolongate_factor, smooth_steps);
        finish = clock();
        timeused_nonperiodic_tesselation = (double(finish)-double(start))/CLOCKS_PER_SEC;
        cout << " (" << timeused_nonperiodic_tesselation << " secs)" << endl << endl;
    }
    if (isPeriodic) {
        cout << "Tesselating periodic grid ...  ";
        start = clock();
        upscaler_periodic.init(eclParser, SinglePhaseUpscaler::Periodic, minPerm,
                               ztol,  linsolver_tolerance, linsolver_verbosity, linsolver_type, twodim_hack,
                               linsolver_maxit, linsolver_prolongate_factor, smooth_steps);
        finish = clock();
        timeused_periodic_tesselation = (double(finish)-double(start))/CLOCKS_PER_SEC;
        cout << " (" << timeused_periodic_tesselation << " secs)" << endl << endl;
    }


    
    
    /*********************************************************************
     * Do porosity upscaling
     *
     * This is an added feature. It is done since does not cost anything
     * in terms of cpu-resources (compared to permeability upscaling).
     */
    double upscaledPorosity = 0.0;
    if (eclParser.hasField("PORO")) {
        if (isPeriodic) {
            upscaledPorosity = upscaler_periodic.upscalePorosity();
        } else {
            upscaledPorosity = upscaler_nonperiodic.upscalePorosity();
        }
    }
 
    /*********************************************************************
     * Do single-phase permeability upscaling 
     */
    
    if (isFixed)  {
        cout << "Compute for fixed boundary conditions: ...  ";
        start = clock();
        upscaler_nonperiodic.setBoundaryConditionType(SinglePhaseUpscaler::Fixed);
        Kfixed = upscaler_nonperiodic.upscaleSinglePhase();
        Kfixed *= 1.0/(Opm::prefix::milli*Opm::unit::darcy);
        finish = clock();
        timeused_fixed = (double(finish)-double(start))/CLOCKS_PER_SEC;
        cout << " ( " << timeused_fixed << " secs)" << endl;
        cout << Kfixed << endl;
        cout << endl;
    }
    
    if (isLinear)  {
        cout << "Compute for linear boundary conditions: ... " << endl;
        start = clock();
        upscaler_nonperiodic.setBoundaryConditionType(SinglePhaseUpscaler::Linear);
        Klinear = upscaler_nonperiodic.upscaleSinglePhase();
        Klinear *= 1.0/(Opm::prefix::milli*Opm::unit::darcy);
        finish = clock();
        timeused_linear = (double(finish)-double(start))/CLOCKS_PER_SEC;
        cout << Klinear << endl;
        cout << " ( " << timeused_linear << " secs)" << endl;
        cout << endl << endl;
    }
    
    if (isPeriodic)  {                
        cout << "Compute for periodic boundary conditions: ... ";
        start = clock();
        upscaler_periodic.setBoundaryConditionType(SinglePhaseUpscaler::Periodic);
        Kperiodic = upscaler_periodic.upscaleSinglePhase();
        Kperiodic *= 1.0/(Opm::prefix::milli*Opm::unit::darcy);
        finish = clock();
        timeused_periodic =  (double(finish)-double(start))/CLOCKS_PER_SEC;     
        cout << " (" << timeused_periodic << " secs)" << endl;           
        cout << Kperiodic << endl;
        cout << endl;
    }     
    
    /***********************************************************************
     * Output results to stdout or optionally to file
     */
    
    stringstream outputtmp;
    
    // Print a table of all computed values:
    outputtmp << "###############################################################################" << endl;
    outputtmp << "# Results from upscaling permeability."<< endl;
    outputtmp << "#" << endl;
    time_t now = time(NULL);
    outputtmp << "# Finished: " << asctime(localtime(&now));
    utsname hostname; uname(&hostname);
    outputtmp << "# Hostname: " << hostname.nodename << endl;

    outputtmp << "#" << endl;
    outputtmp << "# Eclipse file: " << ECLIPSEFILENAME << endl;
    outputtmp << "# Porosity : "  << upscaledPorosity << endl;
    outputtmp << "#" << endl;
    outputtmp << "# Options used:" << endl;
    outputtmp << "#     Boundary conditions: ";
    if (isFixed)    outputtmp << "Fixed (no-flow)  "; 
    if (isPeriodic) outputtmp << "Periodic  ";
    if (isLinear)   outputtmp << "Linear  ";
    outputtmp << endl;
    outputtmp << "#                 minPerm: " << options["minPerm"] << endl;
    outputtmp << "#" << endl;
    outputtmp << "# If both linear and fixed boundary conditions are calculated, " << endl <<
        "# the nonperiodic tesselation is done only once" << endl <<  "# " << endl  << "#" << endl;
    
    if (isFixed) {
        outputtmp << "# Upscaled permeability for fixed boundary conditions:" << endl;
        outputtmp << "# Tesselation time: " << timeused_nonperiodic_tesselation << " s" << endl;
        outputtmp << "# Computation time: " << timeused_fixed << " s" << endl;
        outputtmp << Kfixed;
    }
    if (isLinear) {
        outputtmp << "# Upscaled permeability for linear boundary conditions:" << endl;
        if  (!isFixed) {
            // Only display this if fixed BC was not used, as the tesselation time is shared.
            outputtmp << "# Tesselation time: " << timeused_nonperiodic_tesselation << " s" << endl;
        }
        outputtmp << "# Computation time: " << timeused_linear << " s" << endl;
        outputtmp << Klinear;
    }
    if (isPeriodic) {
        outputtmp << "# Upscaled permeability for periodic boundary conditions:" << endl;
        outputtmp << "# Tesselation time: " << timeused_periodic_tesselation << " s" << endl;
        outputtmp << "# Computation time: " << timeused_periodic << " s" << endl;
        outputtmp << Kperiodic;
    }    
    cout << endl << outputtmp.str();
    
    if (options["output"] != "")  {
        cout << "Writing results to " << options["output"] << endl;
        ofstream outfile;
        outfile.open(options["output"].c_str(), ios::out | ios::trunc);
        outfile << outputtmp.str();
        outfile.close();
    }
    return 0;   
}

/**
   @brief Upscales permeability
   
   @param varnum Number of input arguments
   @param vararg Input arguments
   @return int
*/
int main(int varnum, char** vararg) {
    return upscale(varnum, vararg);
}
