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


#include "CornerpointChopper.hpp"
#include "dune/upscaling/SinglePhaseUpscaler.hpp"
#include <dune/porsol/common/setupBoundaryConditions.hpp>
#include <dune/common/Units.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <iomanip>


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
    std::string filebase = param.getDefault<std::string>("filebase", "");




    // Random number generator from boost.
    boost::mt19937 gen;
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
        ch.chop(istart, istart + ilen, jstart, jstart + jlen, zstart, zstart + zlen);
        std::string outname = filebase;
        // Output to file if a filebase is supplied.
        if (outname != "") {
            std::ostringstream oss;
            oss << 'R' << std::setw(4) << std::setfill('0') << sample;
            outname += oss.str();
            outname += ".grdecl";
            ch.writeGrdecl(outname);
        }

        Dune::EclipseGridParser subparser = ch.subparser();

	Dune::SinglePhaseUpscaler upscaler;
	upscaler.init(subparser, Dune::SinglePhaseUpscaler::Fixed, 1e-9);

	Dune::SinglePhaseUpscaler::permtensor_t upscaled_K = upscaler.upscaleSinglePhase();
        upscaled_K *= (1.0/(Dune::prefix::milli*Dune::unit::darcy));


	porosities.push_back(upscaler.upscalePorosity());
        permxs.push_back(upscaled_K(0,0));
        permys.push_back(upscaled_K(1,1));
        permzs.push_back(upscaled_K(2,2));

    }
    
    // Output results
    // (TODO: Redirect to user-supplied filename)
    
    for (int sample = 1; sample <= subsamples; ++sample) {
        std::cout << sample << '\t' << porosities[sample-1] << 
            permxs[sample-1] << '\t' <<
            permys[sample-1] << '\t' <<
            permzs[sample-1] << '\t' <<
            std::endl;
    }
    
}
