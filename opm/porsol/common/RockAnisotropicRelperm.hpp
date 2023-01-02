//===========================================================================
//
// File: RockAnisotropicRelperm.hpp
//
// Created: Fri Oct 23 14:11:24 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Bård Skaflestad     <bard.skaflestad@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

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

#ifndef OPENRS_ROCKANISOTROPICRELPERM_HEADER
#define OPENRS_ROCKANISOTROPICRELPERM_HEADER


#include <dune/common/fvector.hh>
#include <opm/common/utility/numeric/NonuniformTableLinear.hpp>
#include <opm/porsol/common/Matrix.hpp>

#include <fstream>
#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>

namespace Opm
{

    class RockAnisotropicRelperm
    {
    public:

        void setUseJfunctionScaling(const bool use_j)
        {
            if (use_j) {
                OPM_THROW(std::runtime_error, "RockAnisotropicRelperm cannot use J-scaling.");
            }
        }

        void setSigmaAndTheta(const double, const double)
        {
            OPM_THROW(std::runtime_error, "RockAnisotropicRelperm cannot accept sigma and theta arguments.");
        }


	template <template <class> class SP, class OP>
	void kr(const int phase_index, const double saturation, FullMatrix<double, SP, OP>& kr_value) const
	{
	    zero(kr_value);
	    kr_value(0,0) = krxx_[phase_index](saturation);
	    kr_value(1,1) = kryy_[phase_index](saturation);
	    kr_value(2,2) = krzz_[phase_index](saturation);
	}

	template <template <class> class SP, class OP>
	double capPress(const FullMatrix<double, SP, OP>& /*perm*/, const double /*poro*/, const double saturation) const
	{
            return cap_press_(saturation);
	}

	template <template <class> class SP, class OP>
	double satFromCapPress(const FullMatrix<double, SP, OP>& /*perm*/, const double /*poro*/, const double cp) const
	{
            return cap_press_.inverse(cp);
	}

	void read(const std::string& directory, const std::string& specification)
	{
	    // For this type of rock, the specification is a line with two file names.
	    std::istringstream specstream(specification);
	    std::string rockname[2];
        specstream >> rockname[0] >> rockname[1];
        for (int i = 0; i < 2; ++i) {
        std::string rockfilename = directory + rockname[i];
        std::ifstream rock_stream(rockfilename.c_str());
        if (!rock_stream) {
            OPM_THROW(std::runtime_error, "Could not open file " + rockfilename);
        }
        readAnisoFormat(rock_stream, i);
        }
	}

    private:
	void readAnisoFormat(std::istream& is, int phase)
	{
	    // Ignore comments.
	    while (is.peek() == '#') {
		is.ignore(std::numeric_limits<int>::max(), '\n');
	    }
	    typedef Dune::FieldVector<double, 5> Data;
	    std::istream_iterator<Data> start(is);
	    std::istream_iterator<Data> end;
	    std::vector<Data> data(start, end);
	    if (!is.eof()) {
		OPM_THROW(std::runtime_error, "Reading stopped but we're not at eof: something went wrong reading data.");
	    }
            // Making sure that the tables start/end at zero relperm.
            bool all_ok = true;
            for (int i = 2; i < 5; ++i) {
                if (phase == 0 && data[0][i] != 0.0) { // Water table
                    data[0][i] = 0.0;
                    all_ok = false;
                } else if (phase == 1 && data.back()[i] != 0.0) { // Oil table
                    data.back()[i] = 0.0;
                    all_ok = false;
                }
            }
            if (!all_ok) {
                std::cout << "Warning: Relperm tables were modified to go to zero at the beginning/end." << std::endl;
            }
	    std::vector<double> pcow, svals, krxx, kryy, krzz;
	    for (int i = 0; i < int(data.size()); ++i) {
		pcow.push_back(data[i][0]);
		svals.push_back(data[i][1]);
		krxx.push_back(data[i][2]);
		kryy.push_back(data[i][3]);
		krzz.push_back(data[i][4]);
	    }
	    krxx_[phase] = TabFunc(svals, krxx);
	    kryy_[phase] = TabFunc(svals, kryy);
	    krzz_[phase] = TabFunc(svals, krzz);
	    if (phase == 0) {
		cap_press_ = TabFunc(svals, pcow);
	    } else {
		assert(phase == 1);
		assert(cap_press_ == TabFunc(svals, pcow));
	    }
	}

	typedef NonuniformTableLinear<double> TabFunc;
	TabFunc krxx_[2];
	TabFunc kryy_[2];
	TabFunc krzz_[2];
	TabFunc cap_press_;
    };

} // namespace Opm



#endif // OPENRS_ROCKANISOTROPICRELPERM_HEADER
