/*
  Copyright 2011 IRIS - International Research Institute of Stavanger.

  This file is part of The Open Reservoir Simulator Project (OpenRS).

  OpenRS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OpenRS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <config.h>

#include <opm/core/utility/have_boost_redef.hpp>

#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <Appleyard.hpp>

using namespace Opm;

bool getSaturationUpdate(double time0,
                         double dt,
                         const std::vector<double> & saturation, 
                         const std::vector<double> & saturation0,
                         std::vector<double> & dSat)
{
	for (unsigned int i=0; i<dSat.size(); ++i) {
		double s=saturation[i];
		double s0=saturation0[i];
		const double a=0.4;
		
		dSat[i] = -(s-s0-a*(s*s*(1.-s*s))*dt) / (1.-2.*a*s*(1.-2.*s*s)*dt);
	}
	
	return true;
}

int main(int argc, char** argv)
{
	parameter::ParameterGroup param(argc, argv);
	
	std::vector<double> sat0(9, 0.5);
	double timeStart=0.0;
	double timeEnd=10.0;

	AppleyardControl ac(sat0, timeStart, timeEnd);

	// Specifies min and max allowed saturation for each phase: <satMin,satMax>
	//   (Black oil: Sw_co <= Sw <= 1-So_r     0 <= Sg <= 1-Sw_co     So <= 1-Sw_co)
	// Layout of the saturation vector: <begin,pastEnd,stride>
	//   (Two options: ...,s1_i,s2_i,s3_i,...   versus  s1_1,...,s1_N,s2_1,...,s2_N,s3_1,...,s3_N) 
	ac.phaseLayout(0.0, 1.0, 0, sat0.size()-2, 3); 
	ac.phaseLayout(0.1, 0.9, 1, sat0.size()-1, 3);
	ac.phaseLayout(0.2, 0.8, 2, sat0.size(), 3);
	
	do {
		// Compute saturation increment:
		// Input:
		//   current time 
		//   current (micro) timestep
		//   current saturation approximation
		//   initial saturation
		// Output:
		//   saturation increment
		// Return value:
		//   true if saturation increment successfully computed
		
		ac.isOk() = getSaturationUpdate(ac.time(), ac.dt(), ac.saturation(), ac.saturationIntermediate(), ac.dSat());
	
	} while(ac.continueUpdate());
	
	if (ac.hasConverged()) {
		ac.finalizeUpdate(sat0);
	}
	else {
		// Time for last successful update:
		ac.time();
		// Last successful saturation update:
		ac.saturationIntermediate();
	}
	
}
