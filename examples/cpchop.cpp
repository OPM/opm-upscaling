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


#include <dune/common/CornerpointChopper.hpp>
#include <dune/upscaling/SinglePhaseUpscaler.hpp>
#include <dune/porsol/common/setupBoundaryConditions.hpp>
#include <dune/common/Units.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include <ios>
#include <iomanip>
#include <sys/utsname.h>
#include <ctime>
#include <sstream>
#include <fstream>
#include <iostream>


int main(int argc, char** argv)
{
    Dune::parameter::ParameterGroup param(argc, argv);
    std::string gridfilename = param.get<std::string>("gridfilename");
    Dune::CornerPointChopper ch(gridfilename);

    // The cells with i coordinate in [imin, imax) are included, similar for j.
    // The z limits may be changed inside the chopper to match actual min/max z.
    const int* dims = ch.dimensions();
    int imin = param.getDefault("imin", 0);
    int imax = param.getDefault("imax", dims[0]);
    int jmin = param.getDefault("jmin", 0);
    int jmax = param.getDefault("jmax", dims[1]);
    double zmin = param.getDefault("zmin", ch.zLimits().first);
    double zmax = param.getDefault("zmax", ch.zLimits().second);
    int subsamples = param.getDefault("subsamples", 1);
    int ilen = param.getDefault("ilen", imax - imin);
    int jlen = param.getDefault("jlen", jmax - jmin);
    double zlen = param.getDefault("zlen", zmax - zmin);
    bool upscale = param.getDefault("upscale", true);
    bool resettoorigin = param.getDefault("resettoorigin", true);
    boost::mt19937::result_type userseed = param.getDefault("seed", 0);

    int outputprecision = param.getDefault("outputprecision", 8);
    std::string filebase = param.getDefault<std::string>("filebase", "");
    std::string resultfile = param.getDefault<std::string>("resultfile", "");

    double z_tolerance = param.getDefault("z_tolerance", 0.0);
    double residual_tolerance = param.getDefault("residual_tolerance", 1e-8);
    double linsolver_verbosity = param.getDefault("linsolver_verbosity", 0);
    double linsolver_type = param.getDefault("linsolver_type", 1);

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
    std::vector<double> porosities;
    std::vector<double> permxs;
    std::vector<double> permys;
    std::vector<double> permzs;

    for (int sample = 1; sample <= subsamples; ++sample) {
        int istart = ri();
        int jstart = rj();
        double zstart = rz();
        ch.chop(istart, istart + ilen, jstart, jstart + jlen, zstart, zstart + zlen, resettoorigin);
        std::string subsampledgrdecl = filebase;

        // Output grdecl-data to file if a filebase is supplied.
        if (filebase != "") {
            std::ostringstream oss;
            oss << 'R' << std::setw(4) << std::setfill('0') << sample;
            subsampledgrdecl += oss.str();
            subsampledgrdecl += ".grdecl";
            ch.writeGrdecl(subsampledgrdecl);
        }

        if (upscale) {
            Dune::EclipseGridParser subparser = ch.subparser();
            
            Dune::SinglePhaseUpscaler upscaler;
            upscaler.init(subparser, Dune::SinglePhaseUpscaler::Fixed, 0.0, z_tolerance,
                          residual_tolerance, linsolver_verbosity, linsolver_type, false);
            
            Dune::SinglePhaseUpscaler::permtensor_t upscaled_K = upscaler.upscaleSinglePhase();
            upscaled_K *= (1.0/(Dune::prefix::milli*Dune::unit::darcy));
            
            
            porosities.push_back(upscaler.upscalePorosity());
            permxs.push_back(upscaled_K(0,0));
            permys.push_back(upscaled_K(1,1));
            permzs.push_back(upscaled_K(2,2));
        
        }   
    }
     
    if (upscale) {
        
        // Make stream of output data, to be outputted to screen and optionally to file
        std::stringstream outputtmp;
        
        outputtmp << "################################################################################################" << std::endl;
        outputtmp << "# Results from property analysis on subsamples" << std::endl;
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
        outputtmp << "#       subsamples: " << subsamples << std::endl;
        outputtmp << "################################################################################################" << std::endl;
        outputtmp << "# id          porosity                 permx                   permy                   permz" << std::endl;
        
        const int fieldwidth = outputprecision + 8;
        for (int sample = 1; sample <= subsamples; ++sample) {
            outputtmp << sample << '\t' <<
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
