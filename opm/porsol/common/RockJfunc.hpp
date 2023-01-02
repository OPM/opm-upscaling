//===========================================================================
//
// File: RockJfunc.hpp
//
// Created: Fri Oct 23 08:59:52 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            B�rd Skaflestad     <bard.skaflestad@sintef.no>
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

#ifndef OPENRS_ROCKJFUNC_HEADER
#define OPENRS_ROCKJFUNC_HEADER

#include <dune/common/fvector.hh>
#include <opm/common/utility/numeric/NonuniformTableLinear.hpp>
#include <opm/porsol/common/Matrix.hpp>

#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>

namespace Opm
{

    class RockJfunc
        {
    public:
        RockJfunc()
            : use_jfunction_scaling_(true), sigma_cos_theta_(1.0)
        {
        }

        void setUseJfunctionScaling(const bool use_j)
        {
            use_jfunction_scaling_ = use_j;
        }

        void setSigmaAndTheta(const double sigma, const double theta)
        {
            sigma_cos_theta_ = sigma*std::cos(theta);
        }

	void krw(const double saturation, double& krw_value) const
	{
	    krw_value = krw_(saturation);
	}

	void kro(const double saturation, double& kro_value) const
	{
	    kro_value = kro_(saturation);
	}

	void dkrw(const double saturation, double& dkrw_value) const
	{
	    dkrw_value = krw_.derivative(saturation);
	}

	void dkro(const double saturation, double& dkro_value) const
	{
	    dkro_value = kro_.derivative(saturation);
	}

	double s_min() const
	{
		return s_min_;
	}
	double s_max() const
	{
		return s_max_;
	}

	template <template <class> class SP, class OP>
	double capPress(const FullMatrix<double, SP, OP>& perm, const double poro, const double saturation) const
	{
            if (use_jfunction_scaling_) {
                // p_{cow} = J\frac{\sigma \cos \theta}{\sqrt{k/\phi}}
                // \sigma \cos \theta is by default approximated by 1.0;
                // k is approximated by the average of the diagonal terms.
                double sqrt_k_phi = std::sqrt(trace(perm)/(perm.numRows()*poro));
                return Jfunc_(saturation)*sigma_cos_theta_/sqrt_k_phi;
            } else {
                // The Jfunc_ table actually contains the pressure directly.
                return Jfunc_(saturation);
            }
	}

	template <template <class> class SP, class OP>
	double capPressDeriv(const FullMatrix<double, SP, OP>& perm, const double poro, const double saturation) const
	{
            if (use_jfunction_scaling_) {
                // p_{cow} = J\frac{\sigma \cos \theta}{\sqrt{k/\phi}}
                // \sigma \cos \theta is by default approximated by 1.0;
                // k is approximated by the average of the diagonal terms.
                double sqrt_k_phi = std::sqrt(trace(perm)/(perm.numRows()*poro));
                return Jfunc_.derivative(saturation)*sigma_cos_theta_/sqrt_k_phi;
            } else {
                // The Jfunc_ table actually contains the pressure directly.
                return Jfunc_.derivative(saturation);
            }
	}

	template <template <class> class SP, class OP>
	double satFromCapPress(const FullMatrix<double, SP, OP>& perm, const double poro, const double cap_press) const
	{
			double s = 0;
            if (use_jfunction_scaling_) {
                // p_{cow} = J\frac{\sigma \cos \theta}{\sqrt{k/\phi}}
                // \sigma \cos \theta is by default approximated by 1.0;
                // k is approximated by the average of the diagonal terms.
                double sqrt_k_phi = std::sqrt(trace(perm)/(perm.numRows()*poro));
                s = Jfunc_.inverse(cap_press*sqrt_k_phi/sigma_cos_theta_);
            } else {
                // The Jfunc_ table actually contains the pressure directly.
                s = Jfunc_.inverse(cap_press);
            }
            s = std::min(s_max_, std::max(s_min_, s));
            return s;
	}

	void read(const std::string& directory, const std::string& specification)
	{
	    // For this type of rock, the specification is simply a line with the file name.
	    std::istringstream specstream(specification);
	    std::string rockname;
	    specstream >> rockname;
            std::string rockfilename = directory + rockname;
            std::ifstream rock_stream(rockfilename.c_str());
            if (!rock_stream) {
                OPM_THROW(std::runtime_error, "Could not open file " + rockfilename);
            }
	    readStatoilFormat(rock_stream);
	}

    private:
	void readStatoilFormat(std::istream& is)
	{
            /* Skip lines at the top of the file starting with '#' or '--' */
            char c = is.peek();
            while ((c == '-' || c == '#') && is) {
                if (c == '#') {
                    std::string commentline;
                    std::getline(is, commentline);
                    c = is.peek();
                }
                else { /* Check if the initial '-' was first byte of the sequence '--' */
                    is.get(c);
                    if (c == '-') {
                        std::string commentline;
                        std::getline(is, commentline);
                        c = is.peek();
                    } else {
                        is.putback(c);
                    }
                }
            }
            if (!is) {
                OPM_THROW(std::runtime_error, "Something went wrong while skipping optional comment header of rock file.");
            }
	    typedef Dune::FieldVector<double, 4> Data;
	    std::istream_iterator<Data> start(is);
	    std::istream_iterator<Data> end;
	    std::vector<Data> data(start, end);
	    if (!is.eof()) {
		OPM_THROW(std::runtime_error, "Reading stopped but we're not at eof: something went wrong reading data of rock file.");
	    }
	    std::vector<double> svals, krw, kro, Jfunc;
	    for (int i = 0; i < int(data.size()); ++i) {
		svals.push_back(data[i][0]);
		krw.push_back(data[i][1]);
		kro.push_back(data[i][2]);
		Jfunc.push_back(data[i][3]);
	    }
            if (krw[0] != 0.0) {
                krw[0] = 0.0;
                std::cout << "Warning: krw table were modified to go to zero at the beginning." << std::endl;
            }
            if (kro.back() != 0.0) {
                kro.back() = 0.0;
                std::cout << "Warning: kro table were modified to go to zero at the end." << std::endl;
            }
            s_min_ = svals.front();
	    s_max_ = svals.back();
	    krw_ = TabFunc(svals, krw);
	    kro_ = TabFunc(svals, kro);
	    Jfunc_ = TabFunc(svals, Jfunc);
	    std::vector<double> invJfunc(Jfunc);
	    std::reverse(invJfunc.begin(), invJfunc.end());
	    std::vector<double> invsvals(svals);
	    std::reverse(invsvals.begin(), invsvals.end());
	    invJfunc_ = TabFunc(invJfunc, invsvals);
	}

	typedef NonuniformTableLinear<double> TabFunc;
	TabFunc krw_;
	TabFunc kro_;
	TabFunc Jfunc_;
	TabFunc invJfunc_;
	bool use_jfunction_scaling_;
	double sigma_cos_theta_;
	double s_min_;
	double s_max_;


    };

} // namespace Opm


#endif // OPENRS_ROCKJFUNC_HEADER
