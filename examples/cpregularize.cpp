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

/**
   Program to regularize cornerpoint grids

   Caveats:
   - Only grids with vertical pillars

   - CornerPointChopper can only chop along existing pillars. In case
     your asked-for horizontal resolution does not divide the initial
     number of pillars in x/y, you will not obtain a fully regular
     grid, but still easier numerically.
     
   - Be careful with non-flat top and bottom boundary.

*/
#include <config.h>

#include <opm/core/io/eclipse/CornerpointChopper.hpp>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/io/eclipse/EclipseGridInspector.hpp>
#include <opm/upscaling/SinglePhaseUpscaler.hpp>
#include <opm/porsol/common/setupBoundaryConditions.hpp>
#include <opm/core/utility/Units.hpp>

#include <ios>
#include <iomanip>
#include <sys/utsname.h>
#include <ctime>
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
   if (argc == 1) {
        std::cout << "Usage: cpregularize gridfilename=filename.grdecl [ires=5] [jres=5] [zres=5] " << std::endl;
        std::cout << "       [imin=] [imax=] [jmin=] [jmax=] [zmin=] [zmax=] " << std::endl;
        std::cout << "       [z_tolerance=0.0] [minperm=1e-9] " << std::endl;
        std::cout << "       [resultgrid=regularizedgrid.grdecl]" << std::endl;
        exit(1);
    }

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
    int ires = param.getDefault("ires", 1);
    int jres = param.getDefault("jres", 1);
    int zres = param.getDefault("zres", 1);

    std::string resultgrid = param.getDefault<std::string>("resultgrid", "regularizedgrid.grdecl");

    double minperm = param.getDefault("minperm", 1e-9);
    double minpermSI = Opm::unit::convert::from(minperm, Opm::prefix::milli*Opm::unit::darcy);
    double z_tolerance = param.getDefault("z_tolerance", 1e-8);
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
    // (ilen, jlen and zlen set to zero, does not apply here)
    ch.verifyInscribedShoebox(imin, 0, imax, 
			      jmin, 0, jmax,
			      zmin, 0, zmax);

    
    // Storage for properties for regularized cells
    std::vector<double> poro;
    std::vector<double> permx;
    std::vector<double> permy;
    std::vector<double> permz;


    // Original x/y resolution in terms of coordinate values (not indices)
    Opm::EclipseGridParser gridparser(gridfilename); // TODO: REFACTOR!!!! it is stupid to parse this again
    Opm::EclipseGridInspector gridinspector(gridparser);
    std::array<double, 6> gridlimits=gridinspector.getGridLimits();
    double finegridxresolution = (gridlimits[1]-gridlimits[0])/dims[0];
    double finegridyresolution = (gridlimits[3]-gridlimits[2])/dims[1];

    // Construct mapping from coarse i and j indices to fine
    // and COORDS values for regularized pillars.
    std::vector<int> iidx_f, jidx_f;
    std::vector<double> newcoords_x;
    int finesprcoarse_i = floor(dims[0] / ires);
    int remainder_i = dims[0] - ires*finesprcoarse_i;
    for (int iidx_c=0; iidx_c < remainder_i+1; ++iidx_c) {
        iidx_f.push_back(iidx_c*(finesprcoarse_i + 1));  // Spread remainder evenly
    }
    for (int iidx_c=remainder_i + 1; iidx_c < ires; ++iidx_c) {
        iidx_f.push_back(iidx_c*finesprcoarse_i + remainder_i);
    }
    iidx_f.push_back(imax); // endpoint needed below

    int finesprcoarse_j = floor(dims[1] / jres);
    int remainder_j = dims[1] - jres*finesprcoarse_j;
    for (int jidx_c=0; jidx_c < remainder_j+1; ++jidx_c) {
        jidx_f.push_back(jidx_c*(finesprcoarse_j + 1)); // Spread remainder evenly
    }
    for (int jidx_c=remainder_j + 1; jidx_c < jres; ++jidx_c) {
        jidx_f.push_back(jidx_c*finesprcoarse_j + remainder_j);
    }
    jidx_f.push_back(jmax); // endpoint needed below

    // Construct new ZCORN for regular grid
    std::vector<double> zcorn_c;
    for (int zidx_c=0; zidx_c < zres; ++zidx_c) {
	zcorn_c.push_back(zmin + zidx_c * (zmax-zmin)/zres);
    }
    zcorn_c.push_back(zmax);
    


    // Run through the new regular grid to find its properties
    for (int zidx_c=0; zidx_c < zres; ++zidx_c) {
        for (int jidx_c=0; jidx_c < jres; ++jidx_c) {
            for (int iidx_c=0; iidx_c < ires; ++iidx_c) {
                Opm::CornerPointChopper::ChopContext context;
                ch.chop(iidx_f[iidx_c], iidx_f[iidx_c+1],
			jidx_f[jidx_c], jidx_f[jidx_c+1],
			zcorn_c[zidx_c], zcorn_c[zidx_c+1],
			context, false);
		try {
		    Opm::EclipseGridParser subparser = ch.subparser(context);
                    subparser.convertToSI(); // Because the upscaler expects SI units.
		    Opm::SinglePhaseUpscaler upscaler;
		    upscaler.init(subparser, Opm::SinglePhaseUpscaler::Fixed, minpermSI, z_tolerance,
				  residual_tolerance, linsolver_verbosity, linsolver_type, false);
            
		    Opm::SinglePhaseUpscaler::permtensor_t upscaled_K = upscaler.upscaleSinglePhase();
		    upscaled_K *= (1.0/(Opm::prefix::milli*Opm::unit::darcy));
		    poro.push_back(upscaler.upscalePorosity());
		    permx.push_back(upscaled_K(0,0));
		    permy.push_back(upscaled_K(1,1));
		    permz.push_back(upscaled_K(2,2)); 
		}
		catch (...) {
		    std::cout << "Warning: Upscaling for cell failed to convert, values set to zero\n";
		    poro.push_back(0.0);
		    permx.push_back(0.0);
		    permy.push_back(0.0);
		    permz.push_back(0.0);
		}
	    }
        }
    }
    // Write regularized grid to outputfile
    std::ofstream out(resultgrid.c_str());
    if (!out) {
        std::cerr << "Could not open file " << resultgrid << "\n";
        throw std::runtime_error("Could not open output file.");
    }
    out << "SPECGRID\n" << ires << ' ' << jres << ' ' << zres 
                << " 1 F\n/\n\n";

    out << "COORD\n";
    for (int j = 0; j <= jres; ++j) {
        for (int i = 0; i <= ires; ++i) {
	    out << finegridxresolution*iidx_f[i] << " " << finegridyresolution*jidx_f[j] << " " << zmin << " "
                << finegridxresolution*iidx_f[i] << " " << finegridyresolution*jidx_f[j] << " " << zmax << "\n";
        }
    }
    out << "/\n\n";

    /*
     Write ZCORN, that is the Z-coordinates along the pillars, specifying
     the eight corners of each cell. Each corner is specified for each
     cell, even though it is the same corner that is used in other
     cells. 

     We loop over corners in each grid cell, directions: z, y, x (x innermost).
     The code here *IS* redundant, but the grid is also very redundant
     for a grid that is really regular..
   */ 
    out << "ZCORN\n";
    double zlen = zmax-zmin;
    for (int zidx=0; zidx < zres; ++zidx) {
	for (int j = 0; j < jres; ++j) {
	    for (int i = 0; i < ires; ++i) {
		out << zlen/zres*zidx << "  " << zlen/zres*zidx << "  ";
	    }
	    out << "\n";
	    for (int i = 0; i < ires; ++i) {
		out << zlen/zres*zidx << "  " << zlen/zres*zidx << "  ";
	    }
	}
	for (int j = 0; j < jres; ++j) {
	    for (int i = 0; i < ires; ++i) {
		out << zlen/zres*(zidx+1) << "  " << zlen/zres*(zidx+1) << "  ";
	    }
	    out << "\n";
	    for (int i = 0; i < ires; ++i) {
		out << zlen/zres*(zidx+1) << "  " << zlen/zres*(zidx+1) << "  ";
	    }
	}
    }
    out << "/\n\n";
    
    out << "PORO\n";
    for (size_t idx=0; idx < (size_t)poro.size(); ++idx) {
	out << poro[idx] << std::endl;
    }
    out << "/\n\n";
    
    out << "PERMX\n";
    for (size_t idx=0; idx < (size_t)permx.size(); ++idx) {
	out << permx[idx] << std::endl;
    }
    out << "/\n\n";
    
    out << "PERMY\n\n";
    for (size_t idx=0; idx < (size_t)permy.size(); ++idx) {
	out << permy[idx] << std::endl;
    }
    out << "/\n\n";
    
    out << "PERMZ\n\n";
    for (size_t idx=0; idx < (size_t)permz.size(); ++idx) {
	out << permz[idx] << std::endl;
    }
    out << "/\n";
    
    out.close();
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}


