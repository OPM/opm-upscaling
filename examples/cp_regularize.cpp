/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2010 Statoil ASA

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
    int ires = param.getDefault("ires", 10);
    int jres = param.getDefault("jres", 10);
    int zres = param.getDefault("zres", 10);

    std::string resultgrid = param.getDefault<std::string>("resultgrid", "regularizedgrid.grdecl");

    double z_tolerance = param.getDefault("z_tolerance", 0.0);
    double residual_tolerance = param.getDefault("residual_tolerance", 1e-8);
    double linsolver_verbosity = param.getDefault("linsolver_verbosity", 0);
    double linsolver_type = param.getDefault("linsolver_type", 1);

    // Check that we do not have any user input 
    // that goes outside the coordinates described in
    // the cornerpoint file (runtime-exception will be thrown in case of error)
    // (ilen, jlen and zlen set to zero, does not apply here)
    ch.verifyInscribedShoebox(imin, 0, imax, 
			      jmin, 0, jmax,
			      zmin, 0, zmax);

    
    // Storage for properties for regularized cells
    std::vector<double> poro;
    std::vector<double> permx;
    std::vector<double> permy;
    std::vector<double> permz;


    // Construct mapping from coarse i and j indices to fine
    std::vector<int> iidx_f, jidx_f;
    int finesprcoarse_i = floor(dims[0] / ires);
    int remainder_i = dims[0] - ires*finesprcoarse;
    for (int iidx_c=0; iidx_c < remainder_i; ++iidx_c) {
        iidx_f.push_back(finesprcoarse_i + 1);  
    }
    for (int iidx_c=remainder_i; iidx_c < ires; ++iidx_c) {
        iidx_f.push_back(finesprcoarse_i);
    }
    int finesprcoarse_j = floor(dims[1] / jres);
    int remainder_j = dims[1] - jres*finesprcoarse_j;
    for (int jidx_c=0; iidx_c < remainder; ++iidx_c) {
        jidx_f.push_back(finesprcoarse_j + 1);
    }
    for (int jidx_c=remainder; jidx_c < jres; ++jidx_c) {
        jidx_f.push_back(finesprcoarse_j);
    }



    // Run through the new regular grid to find its properties
    for (zidx_c=0; zidx_c < zres; ++zidx_c) {
        for (jidx_c=0; jidx_c < jres; ++jidx_c) {
            for (iidx_c=0; iidx_c < ires; ++iidx_c) {
                cp.chop(iidx_f(iidx_c), iidx_f(iidx_c+1),
                    jidx_f(jidx_c), jidx_f(jidx_c+1),
                    zcorn_c(zidc_c), zcorn_c(zidx_c+1),
                    false);
                Dune::EclipseGridParser subparser = ch.subparser();
                Dune::SinglePhaseUpscaler upscaler;
                upscaler.init(subparser, Dune::SinglePhaseUpscaler::Fixed, 0.0, z_tolerance,
                          residual_tolerance, linsolver_verbosity, linsolver_type, false);
            
                Dune::SinglePhaseUpscaler::permtensor_t upscaled_K = upscaler.upscaleSinglePhase();
                upscaled_K *= (1.0/(Dune::prefix::milli*Dune::unit::darcy));
                poro.push_back(upscaler.upscalePorosity());
                permx.push_back(upscaled_K(0,0));
                permy.push_back(upscaled_K(1,1));
                permz.push_back(upscaled_K(2,2)); 
            }
        }
    }
    // Write regularized grid to outputfile
    std::ofstream out(resultgrid.c_str());
    if (!out) {
        std::cerr << "Could not open file " << filename << "\n";
        throw std::runtime_error("Could not open output file.");
    }
    out << "SPECGRID\n" << ires << ' ' << jres << ' ' << zres 
                << " 1 F\n/\n\n";

    out << "COORD\n";
    for (int j = 0; j < jres; ++j) {
        for (int i = 0; i < ires; ++i) {
            // Calculate these on the fly instead!
            out << x(i) << " " << y(j) << " " << zmin << " "
                << x(i) << " " << y(j) << " " << zmax << "\n";
        }
    }

    /*
     Write ZCORN, that is the Z-coordinates along the pillars, specifying
     the eight corners of each cell. Each corner is specified for each
     cell, even though it is the same corner that is used in other
     cells. 

     We loop over corners in each grid cell, directions: z, y, x (x innermost).
     The code here *IS* redundant, but the grid is also very redundant
     for a grid that is really uniform..
   */ 
  out << "ZCORN\n";
 //  for zidx=1:numel(z)-1,
 //    for yidx=1:numel(y)-1,
  //       for xidx=1:numel(x)-1,
   //         fprintf(outfile, ' %g %g', ...
    //                z(zidx), z(zidx));
     //    end
      //   fprintf(outfile, '\n');
       //  for xidx=1:numel(x)-1,
        //    fprintf(outfile, ' %g %g', ...
         //           z(zidx), z(zidx));
   //      end
    //  end
     // for yidx=1:numel(y)-1,
      //   for xidx=1:numel(x)-1,
       //     fprintf(outfile, ' %g %g', ...
        //            z(zidx+1), z(zidx+1));
      //   end
       //  fprintf(outfile, '\n');
        // for xidx=1:numel(x)-1,
         //   fprintf(outfile, ' %g %g', ...
          //          z(zidx+1), z(zidx+1));
     //    end
     // end

  // end
   //fprintf(outfile, '\n/\n\n');
    out.close();
}

