//===========================================================================
//                                                                           
// File: MiscibilityLiveGas.cpp                                               
//                                                                           
// Created: Wed Feb 10 09:21:53 2010                                         
//                                                                           
// Author: Bjørn Spjelkavik <bsp@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
//===========================================================================
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

#include "config.h"
#include <algorithm>
#include "MiscibilityLiveGas.hpp"
#include <opm/common/ErrorMacros.hpp>
#include <opm/core/utility/linearInterpolation.hpp>

using namespace std;
using namespace Opm;

namespace Opm
{

    //------------------------------------------------------------------------
    // Member functions
    //-------------------------------------------------------------------------

    /// Constructor
    MiscibilityLiveGas::MiscibilityLiveGas(const Opm::PvtgTable& pvtgTable)
    {
        // GAS, PVTG
        const auto &saturatedPvtgTable = pvtgTable.getSaturatedTable();

        saturated_gas_table_.resize(4);
        saturated_gas_table_[0] = saturatedPvtgTable.getColumn("PG").vectorCopy();
        saturated_gas_table_[1] = saturatedPvtgTable.getColumn("BG").vectorCopy();
        saturated_gas_table_[2] = saturatedPvtgTable.getColumn("MUG").vectorCopy();
        saturated_gas_table_[3] = saturatedPvtgTable.getColumn("RV").vectorCopy();

        int sz = saturated_gas_table_[0].size();
        undersat_gas_tables_.resize(sz);
        for (int i=0; i<sz; ++i) {
            const auto &undersaturatedPvtgTable = pvtgTable.getUnderSaturatedTable(i);

            undersat_gas_tables_[i][0] = undersaturatedPvtgTable.getColumn("RV").vectorCopy();
            undersat_gas_tables_[i][1] = undersaturatedPvtgTable.getColumn("BG").vectorCopy();
            undersat_gas_tables_[i][2] = undersaturatedPvtgTable.getColumn("MUG").vectorCopy();
        }
    }

    // Destructor
     MiscibilityLiveGas::~MiscibilityLiveGas()
    {
    }

    double MiscibilityLiveGas::getViscosity(int /*region*/, double press, const surfvol_t& surfvol) const
    {
	return miscible_gas(press, surfvol, 2, false);
    }

    void MiscibilityLiveGas::getViscosity(const std::vector<PhaseVec>& pressures,
                                          const std::vector<CompVec>& surfvol,
                                          int phase,
                                          std::vector<double>& output) const
    {
        assert(pressures.size() == surfvol.size());
        int num = pressures.size();
        output.resize(num);
        for (int i = 0; i < num; ++i) {
            output[i] = miscible_gas(pressures[i][phase], surfvol[i], 2, false);
        }
    }

    // Vaporised oil-gas ratio.
    double MiscibilityLiveGas::R(int /*region*/, double press, const surfvol_t& surfvol) const
    {
        if (surfvol[Liquid] == 0.0) {
            return 0.0;
        }
	double R = linearInterpolation(saturated_gas_table_[0],
					     saturated_gas_table_[3], press);
	double maxR = surfvol[Liquid]/surfvol[Vapour];
	if (R < maxR ) {  // Saturated case
	    return R;
	} else {
	    return maxR;  // Undersaturated case
	}
    }

    void MiscibilityLiveGas::R(const std::vector<PhaseVec>& pressures,
                               const std::vector<CompVec>& surfvol,
                               int phase,
                               std::vector<double>& output) const
    {
        assert(pressures.size() == surfvol.size());
        int num = pressures.size();
        output.resize(num);
        for (int i = 0; i < num; ++i) {
            output[i] = R(0, pressures[i][phase], surfvol[i]);
        }
    }

    // Vaporised oil-gas ratio derivative
    double MiscibilityLiveGas::dRdp(int /*region*/, double press, const surfvol_t& surfvol) const
    {
	double R = linearInterpolation(saturated_gas_table_[0],
					     saturated_gas_table_[3], press);
	double maxR = surfvol[Liquid]/surfvol[Vapour];
	if (R < maxR ) {  // Saturated case
	    return linearInterpolationDerivative(saturated_gas_table_[0],
					    saturated_gas_table_[3],
					    press);
	} else {
	    return 0.0;  // Undersaturated case
	}	
    }

    void MiscibilityLiveGas::dRdp(const std::vector<PhaseVec>& pressures,
                                  const std::vector<CompVec>& surfvol,
                                  int phase,
                                  std::vector<double>& output_R,
                                  std::vector<double>& output_dRdp) const
    {
        assert(pressures.size() == surfvol.size());
        R(pressures, surfvol, phase, output_R);
        int num = pressures.size();
        output_dRdp.resize(num);
        for (int i = 0; i < num; ++i) {
            output_dRdp[i] = dRdp(0, pressures[i][phase], surfvol[i]); // \TODO Speedup here by using already evaluated R.
        }
    }

    double MiscibilityLiveGas::B(int /*region*/, double press, const surfvol_t& surfvol) const
    {
        if (surfvol[Vapour] == 0.0) return 1.0; // To handle no-gas case.
        return  miscible_gas(press, surfvol, 1, false);
    }

    void MiscibilityLiveGas::B(const std::vector<PhaseVec>& pressures,
                               const std::vector<CompVec>& surfvol,
                               int phase,
                               std::vector<double>& output) const
    {
        assert(pressures.size() == surfvol.size());
        int num = pressures.size();
        output.resize(num);
        for (int i = 0; i < num; ++i) {
            output[i] = B(0, pressures[i][phase], surfvol[i]);
        }
    }

    double MiscibilityLiveGas::dBdp(int /*region*/, double press, const surfvol_t& surfvol) const
    {	
        if (surfvol[Vapour] == 0.0) return 0.0; // To handle no-gas case.
        return miscible_gas(press, surfvol, 1, true);
    }

    void MiscibilityLiveGas::dBdp(const std::vector<PhaseVec>& pressures,
                                  const std::vector<CompVec>& surfvol,
                                  int phase,
                                  std::vector<double>& output_B,
                                  std::vector<double>& output_dBdp) const
    {
        assert(pressures.size() == surfvol.size());
        B(pressures, surfvol, phase, output_B);
        int num = pressures.size();
        output_dBdp.resize(num);
        for (int i = 0; i < num; ++i) {
            output_dBdp[i] = dBdp(0, pressures[i][phase], surfvol[i]); // \TODO Speedup here by using already evaluated B.
        }
    }

    double MiscibilityLiveGas::miscible_gas(double press, const surfvol_t& surfvol, int item,
					    bool deriv) const
    {
	int section;
	double R = linearInterpolation(saturated_gas_table_[0],
					     saturated_gas_table_[3], press,
					     section);
	double maxR = surfvol[Liquid]/surfvol[Vapour];
	if (deriv) {
	    if (R < maxR ) {  // Saturated case
		return linearInterpolationDerivative(saturated_gas_table_[0],
						saturated_gas_table_[item],
						press);
	    } else {  // Undersaturated case
		int is = section;
		if (undersat_gas_tables_[is][0].size() < 2) {
		    double val = (saturated_gas_table_[item][is+1]
				  - saturated_gas_table_[item][is]) /
			(saturated_gas_table_[0][is+1] -
			 saturated_gas_table_[0][is]);
		    return val;
		}
		double val1 =
		    linearInterpolation(undersat_gas_tables_[is][0],
					      undersat_gas_tables_[is][item],
					      maxR);
		double val2 = 
		    linearInterpolation(undersat_gas_tables_[is+1][0],
					      undersat_gas_tables_[is+1][item],
					      maxR);
		double val = (val2 - val1)/
		    (saturated_gas_table_[0][is+1] - saturated_gas_table_[0][is]);
		return val;
	    }
	} else {
	    if (R < maxR ) {  // Saturated case
		return linearInterpolation(saturated_gas_table_[0],
						 saturated_gas_table_[item],
						 press);
	    } else {  // Undersaturated case
		int is = section;
		// Extrapolate from first table section
		if (is == 0 && press < saturated_gas_table_[0][0]) {
		    return linearInterpolation(undersat_gas_tables_[0][0],
						     undersat_gas_tables_[0][item],
						     maxR);
		}

		// Extrapolate from last table section
		int ltp = saturated_gas_table_[0].size() - 1;
		if (is+1 == ltp && press > saturated_gas_table_[0][ltp]) {
		    return linearInterpolation(undersat_gas_tables_[ltp][0],
						     undersat_gas_tables_[ltp][item],
						     maxR);
		}

		// Interpolate between table sections
		double w = (press - saturated_gas_table_[0][is]) /
		    (saturated_gas_table_[0][is+1] - 
		     saturated_gas_table_[0][is]);
		if (undersat_gas_tables_[is][0].size() < 2) {
		    double val = saturated_gas_table_[item][is] +
			w*(saturated_gas_table_[item][is+1] -
			   saturated_gas_table_[item][is]);
		    return val;
		}
		double val1 =
		    linearInterpolation(undersat_gas_tables_[is][0],
					      undersat_gas_tables_[is][item],
					      maxR);
		double val2 = 
		    linearInterpolation(undersat_gas_tables_[is+1][0],
					      undersat_gas_tables_[is+1][item],
					      maxR);
		double val = val1 + w*(val2 - val1);
		return val;
	    }
	}
    }


} // namespace Opm
