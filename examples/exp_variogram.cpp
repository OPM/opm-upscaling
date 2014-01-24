/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2010 Statoil ASA.

  This file is part of the Open Porous Media project (OPM).

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

/*
  This program computes data for an experimental 
  variogram from a cornerpoint geometry with properties.

  This works by choosing pairs of volumes chosen
  randomly, and comparing their distance with the
  difference in their properties (poro or perm)

  Direction for pairing is set from the command line

*/
#include <config.h>

#include <opm/core/io/eclipse/CornerpointChopper.hpp>
#include <opm/upscaling/SinglePhaseUpscaler.hpp>
#include <opm/porsol/common/setupBoundaryConditions.hpp>
#include <opm/core/utility/Units.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include <ios>
#include <iomanip>
#include <sys/utsname.h>
#include <ctime>
#include <cmath>
#include <sstream>
#include <fstream>
#include <iostream>

#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
#include <dune/common/parallel/mpihelper.hh>
#else
#include <dune/common/mpihelper.hh>
#endif

int main(int argc, char** argv)
try
{
    Dune::MPIHelper::instance(argc, argv);

    Opm::parameter::ParameterGroup param(argc, argv);
    std::string gridfilename = param.get<std::string>("gridfilename");
    Opm::CornerPointChopper ch(gridfilename);

    // The cells with i coordinate in [imin, imax) are included, similar for j.
    // The z limits may be changed inside the chopper to match actual min/max z.
    const int* dims = ch.dimensions();
    int imin = param.getDefault("imin", 0);
    int imax = param.getDefault("imax", dims[0]);
    int jmin = param.getDefault("jmin", 0);
    int jmax = param.getDefault("jmax", dims[1]);
    double zmin = param.getDefault("zmin", ch.zLimits().first);
    double zmax = param.getDefault("zmax", ch.zLimits().second);
    int pairs = param.getDefault("pairs", 100);
    std::string direction = param.get<std::string>("direction");
    int ilen = param.getDefault("ilen", 0);
    int jlen = param.getDefault("jlen", 0);
    double zlen = param.getDefault("zlen", 0.0);
    boost::mt19937::result_type userseed = param.getDefault("seed", 0);

    int outputprecision = param.getDefault("outputprecision", 8);
    std::string resultfile = param.getDefault<std::string>("resultfile", "");

    double minperm = param.getDefault("minperm", 1e-9);
    double minpermSI = Opm::unit::convert::from(minperm, Opm::prefix::milli*Opm::unit::darcy);
    double z_tolerance = param.getDefault("z_tolerance", 0.0);
    double residual_tolerance = param.getDefault("residual_tolerance", 1e-8);
    int linsolver_verbosity = param.getDefault("linsolver_verbosity", 0);
    int linsolver_type = param.getDefault("linsolver_type", 1);
    
    // Check for unused parameters (potential typos).
    if (param.anyUnused()) {
	std::cout << "*****     WARNING: Unused parameters:     *****\n";
	param.displayUsage();
    }


    if (ilen <= 0) {
        std::cerr << "Error: ilen (" << ilen << ") must be greater than zero\n";
        exit(1);
    }
    if (jlen <= 0) {
        std::cerr << "Error: jlen (" << jlen <<") must be greater than zero\n";
        exit(1);
    }
    if (zlen <= 0.0) {
        std::cerr << "Eror: zlen (" << zlen <<") must be greater than zero\n";
        exit(1);
    }

    // Check user supplied variogram direction, either horizontal or vertical 
    enum variogram_directions { undefined, horizontal, vertical };
    variogram_directions variogram_direction = undefined;
    std::string distancemetric;
    if ( direction == "vertical" || direction == "vert" ) {
	variogram_direction = vertical;
	distancemetric = "zcorn";
    }
    else { // if (direction == "horizontal" || direction == "horiz") {
	variogram_direction = horizontal;
	distancemetric = "cells";
    }


    if (variogram_direction == undefined) {
        std::cerr << "Error: variogram direction is undefined, user supplied '" << direction << "'.\n";
        exit(1);
    }
    
    
    // Check that we do not have any user input 
    // that goes outside the coordinates described in
    // the cornerpoint file (runtime-exception will be thrown in case of error)
    ch.verifyInscribedShoebox(imin, ilen, imax,
                              jmin, jlen, jmax,
                              zmin, zlen, zmax);
    
    // Random number generator from boost.
    boost::mt19937 gen;
    
    // Seed the random number generators with the current time, unless specified on command line
    // Warning: Current code does not allow 0 for the seed!!
    if (userseed == 0) {
        gen.seed(time(NULL));
    }
    else {
        gen.seed(userseed);
    }
        

    // Note that end is included in interval for uniform_int.
    boost::uniform_int<> disti(imin, imax - ilen);
    boost::uniform_int<> distj(jmin, jmax - jlen);
    boost::uniform_real<> distz(zmin, std::max(zmax - zlen, zmin));
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > ri(gen, disti);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > rj(gen, distj);
    boost::variate_generator<boost::mt19937&, boost::uniform_real<> > rz(gen, distz);
    
    // Storage for results
    std::vector<double> distances;
    std::vector<double> porodiffs;
    std::vector<double> permxdiffs;
    std::vector<double> permydiffs;
    std::vector<double> permzdiffs;
    
    for (int pair = 1; pair <= pairs; ++pair) {
        int istart_1 = ri();
        int jstart_1 = rj();
        double zstart_1 = rz();
        Opm::CornerPointChopper::ChopContext context;
        ch.chop(istart_1, istart_1 + ilen, jstart_1, jstart_1 + jlen, zstart_1, zstart_1 + zlen, context, false);
	
        Opm::EclipseGridParser subparser_1 = ch.subparser(context);
        subparser_1.convertToSI();
        Opm::SinglePhaseUpscaler upscaler_1;
        upscaler_1.init(subparser_1, Opm::SinglePhaseUpscaler::Fixed, minpermSI, z_tolerance,
			residual_tolerance, linsolver_verbosity, linsolver_type, false);
        Opm::SinglePhaseUpscaler::permtensor_t upscaled_K_1 = upscaler_1.upscaleSinglePhase();
        upscaled_K_1 *= (1.0/(Opm::prefix::milli*Opm::unit::darcy));
	double porosity_1 = upscaler_1.upscalePorosity();

        // Pick another location to form a location-pair to be compared
        int istart_2 = 0, jstart_2 = 0;
        double zstart_2 = 0;
        if (variogram_direction == horizontal) {
            istart_2 = ri();
            jstart_2 = rj();
            zstart_2 = zstart_1;
        }
        else if (variogram_direction == vertical) {
            istart_2 = istart_1;
            jstart_2 = jstart_1;
            zstart_2 = rz();
        }   
        Opm::CornerPointChopper::ChopContext context2;
        ch.chop(istart_2, istart_2 + ilen, jstart_2, jstart_2 + jlen, zstart_2, zstart_2 + zlen, context2, false);
	
        Opm::EclipseGridParser subparser_2 = ch.subparser(context2);
	subparser_2.convertToSI();
        Opm::SinglePhaseUpscaler upscaler_2;
        upscaler_2.init(subparser_2, Opm::SinglePhaseUpscaler::Fixed, minpermSI, z_tolerance,
			residual_tolerance, linsolver_verbosity, linsolver_type, false);
        Opm::SinglePhaseUpscaler::permtensor_t upscaled_K_2 = upscaler_2.upscaleSinglePhase();
        upscaled_K_2 *= (1.0/(Opm::prefix::milli*Opm::unit::darcy));
	double porosity_2 = upscaler_2.upscalePorosity();

        if (variogram_direction == horizontal) {
            distances.push_back(sqrt(pow(istart_2 - istart_1,2) + pow(jstart_2 - jstart_1,2)));
	}
        else if (variogram_direction == vertical) {
            distances.push_back(fabs(zstart_2 - zstart_1));
        }
        porodiffs.push_back(fabs(porosity_2 - porosity_1));
        permxdiffs.push_back(fabs(upscaled_K_2(0,0) - upscaled_K_1(0,0)));
	permydiffs.push_back(fabs(upscaled_K_2(1,1) - upscaled_K_1(1,1)));
        permzdiffs.push_back(fabs(upscaled_K_2(2,2) - upscaled_K_1(2,2)));
    }
    
    // Make stream of output data, to be outputted to screen and optionally to file
    std::stringstream outputtmp;
    
    outputtmp << "################################################################################################" << std::endl;
    outputtmp << "# Data for experimental variogram" << std::endl;
    outputtmp << "#" << std::endl;
    time_t now = time(NULL);
    outputtmp << "# Finished: " << asctime(localtime(&now));
        
    utsname hostname;   uname(&hostname);
    outputtmp << "# Hostname: " << hostname.nodename << std::endl;
    outputtmp << "#" << std::endl;
    outputtmp << "# Options used:" << std::endl;
    outputtmp << "#     gridfilename: " << gridfilename << std::endl;
    outputtmp << "#   i; min,len,max: " << imin << " " << ilen << " " << imax << std::endl;
    outputtmp << "#   j; min,len,max: " << jmin << " " << jlen << " " << jmax << std::endl;
    outputtmp << "#   z; min,len,max: " << zmin << " " << zlen << " " << zmax << std::endl;
    outputtmp << "#        direction: " << direction << std::endl;
    outputtmp << "#            pairs: " << pairs << std::endl;
    outputtmp << "################################################################################################" << std::endl;
    outputtmp << "# distance (" << distancemetric << ")      porositydiff                permxdiff               permydiff                permzdiff" << std::endl;
    
    const int fieldwidth = outputprecision + 8;
    for (int pair = 1; pair <= pairs; ++pair) {
	outputtmp << std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << distances[pair-1] << '\t' <<
	    std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << porodiffs[pair-1] << '\t' <<
	    std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << permxdiffs[pair-1] << '\t' <<
	    std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << permydiffs[pair-1] << '\t' <<
	    std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << permzdiffs[pair-1] << '\t' <<
	    std::endl;
    }
    
    if (resultfile != "") {
	std::cout << "Writing results to " << resultfile << std::endl;
	std::ofstream outfile;
	outfile.open(resultfile.c_str(), std::ios::out | std::ios::trunc);
	outfile << outputtmp.str();
	outfile.close();      
    }
    
        
    
    std::cout << outputtmp.str();
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}

