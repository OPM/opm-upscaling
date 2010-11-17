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

    // Run through the new regular grid to find its properties

    // Write regularized grid to outputfile
}
