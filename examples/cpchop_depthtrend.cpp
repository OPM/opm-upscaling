/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.

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

#include <config.h>

#include <opm/common/utility/platform_dependent/disable_warnings.h>

#include <dune/common/version.hh>

#include <dune/common/parallel/mpihelper.hh>

#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/upscaling/CornerpointChopper.hpp>

#include <opm/porsol/common/setupBoundaryConditions.hpp>

#include <opm/upscaling/SinglePhaseUpscaler.hpp>

#include <cmath> // for min()
#include <ctime>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <sstream>

#include <sys/utsname.h>

#include <random>

/**
   This program is a variant of cpchop. Instead of subsampling randomly, 
   it picks subsamples downwards in a model. It is specifically designed
   for extracting laterally extensive subsamples in specified intervals
   downwards through a model, in order to analyze depth trends in 
   porosity in particular
   
   Default is to pick a new subsample every meter (zresolution=1)
   The subsample height might be smaller or larger than zresolution.
*/


int main(int argc, char** argv)
try
{
   if (argc == 1) {
        std::cout << "Usage: cpchop_depthtrend gridfilename=filename.grdecl [zresolution=1] [zlen=1] [ilen=5] [jlen=5] " << std::endl;
        std::cout << "       [zlen=5] [imin=] [imax=] [jmin=] [jmax=] [upscale=true] [resettoorigin=true]" << std::endl;
        std::cout << "       [seed=111] [minperm=1e-9] " << std::endl;
        exit(1);
    }

    Dune::MPIHelper::instance(argc, argv);

    Opm::ParameterGroup param(argc, argv);
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
    int ilen = param.getDefault("ilen", imax - imin);
    int jlen = param.getDefault("jlen", jmax - jmin);
    double zresolution = param.getDefault("zresolution", 1.0);
    double zlen = param.getDefault("zlen", zresolution);
    bool upscale = param.getDefault("upscale", true);
    bool resettoorigin = param.getDefault("resettoorigin", true);
    std::mt19937::result_type userseed = param.getDefault("seed", 0);

    int outputprecision = param.getDefault("outputprecision", 8);
    std::string filebase = param.getDefault<std::string>("filebase", "");
    std::string resultfile = param.getDefault<std::string>("resultfile", "");

    double minperm = param.getDefault("minperm", 1e-9);
    double minpermSI = Opm::unit::convert::from(minperm, Opm::prefix::milli*Opm::unit::darcy);

    if (param.has("z_tolerance")) {
        std::cerr << "****** Warning: z_tolerance parameter is obsolete, use PINCH in deck input instead\n";
    }
    double residual_tolerance = param.getDefault("residual_tolerance", 1e-8);
    double linsolver_verbosity = param.getDefault("linsolver_verbosity", 0);
    double linsolver_type = param.getDefault("linsolver_type", 1);
        
    // Check for unused parameters (potential typos).
    if (param.anyUnused()) {
	std::cout << "*****     WARNING: Unused parameters:     *****\n";
	param.displayUsage();
    }
    
    // Check that we do not have any user input 
    // that goes outside the coordinates described in
    // the cornerpoint file (runtime-exception will be thrown in case of error)
    ch.verifyInscribedShoebox(imin, ilen, imax, 
			      jmin, jlen, jmax,
			      zmin, zlen, zmax);

    std::mt19937 gen;
    
    // Seed the random number generators with the current time, unless specified on command line
    // Warning: Current code does not allow 0 for the seed!!
    if (userseed == 0) {
        gen.seed(time(NULL));
    }
    else {
        gen.seed(userseed);
    }
        

    // Note that end is included in interval for uniform_int.
    std::uniform_int_distribution<> disti(imin, imax - ilen);
    std::uniform_int_distribution<> distj(jmin, jmax - jlen);
    auto ri = [&disti, &gen] { return disti(gen); };
    auto rj = [&distj, &gen] { return distj(gen); };
    
    // Storage for results
    std::vector<double> zstarts;
    std::vector<double> porosities;
    std::vector<double> permxs;
    std::vector<double> permys;
    std::vector<double> permzs;

    
    /* z_start is the topmost point of the subsample to extract */
    for (double zstart = 0.0; zstart  <= zmax-zlen; zstart += zresolution) {
        /* Horizontally, we pick by random, even though default behaviour is
           to have ilen=imax-min so that there is no randomness */
        int istart = ri();
        int jstart = rj();
        ch.chop(istart, istart + ilen, jstart, jstart + jlen, zstart, std::min(zstart + zlen,zmax), resettoorigin);
        std::string subsampledgrdecl = filebase;

        // Output grdecl-data to file if a filebase is supplied.
        if (filebase != "") {
            std::ostringstream oss;
            oss << 'Z' << std::setw(4) << std::setfill('0') << zstart;
            subsampledgrdecl += oss.str();
            subsampledgrdecl += ".grdecl";
            ch.writeGrdecl(subsampledgrdecl);
        }

        try { /* The upscaling may fail to converge on icky grids, lets just pass by those */
            if (upscale) {
                auto subdeck = ch.subDeck();
                Opm::SinglePhaseUpscaler upscaler;
                upscaler.init(subdeck, Opm::SinglePhaseUpscaler::Fixed, minpermSI,
                              residual_tolerance, linsolver_verbosity, linsolver_type, false);
                
                Opm::SinglePhaseUpscaler::permtensor_t upscaled_K = upscaler.upscaleSinglePhase();
                upscaled_K *= (1.0/(Opm::prefix::milli*Opm::unit::darcy));
                
                
                zstarts.push_back(zstart);
                porosities.push_back(upscaler.upscalePorosity());
                permxs.push_back(upscaled_K(0,0));
                permys.push_back(upscaled_K(1,1));
                permzs.push_back(upscaled_K(2,2));
                
            }   
        }
        catch (...) {
            std::cerr << "Warning: Upscaling chopped subsample at z=" << zstart << "failed, proceeding to next subsample\n";
        }
    }
     
    if (upscale) {
        
        // Make stream of output data, to be outputted to screen and optionally to file
        std::stringstream outputtmp;
        
        outputtmp << "################################################################################################" << std::endl;
        outputtmp << "# Results from depth trend analysis on subsamples" << std::endl;
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
        outputtmp << "#      zresolution: " << zresolution << std::endl;
        outputtmp << "################################################################################################" << std::endl;
        outputtmp << "# zstart          porosity                 permx                   permy                   permz" << std::endl;
        
        const int fieldwidth = outputprecision + 8;
        for (size_t sample = 1; sample <= porosities.size(); ++sample) {
            outputtmp << 
                std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << zstarts[sample-1] << '\t' <<
                std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << porosities[sample-1] << '\t' <<
                std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << permxs[sample-1] << '\t' <<
                std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << permys[sample-1] << '\t' <<
                std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << permzs[sample-1] << '\t' <<
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
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}

