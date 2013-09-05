//===========================================================================
//
// File: SteadyStateUpscalerManager.hpp
//
// Created: Thu Apr 29 11:12:56 2010
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Jostein R Natvig    <jostein.r.natvig@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2010 Statoil ASA.

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

#ifndef OPM_STEADYSTATEUPSCALERMANAGER_HEADER
#define OPM_STEADYSTATEUPSCALERMANAGER_HEADER


#include <opm/upscaling/SteadyStateUpscaler.hpp>
#include <opm/upscaling/UpscalingTraits.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/SparseTable.hpp>
#include <cmath>
#include <fstream>
#include <iostream>


namespace Opm
{


    /// Reads saturation and pressure drop data from an input stream.
    template <class Istream>
    void readControl(Istream& is, std::vector<double>& saturations, Opm::SparseTable<double>& all_pdrops)
    {
        int num_sats;
        is >> num_sats;
        std::vector<double> sat(num_sats);
        Opm::SparseTable<double> ap;
        std::vector<double> tmp_pd;
        for (int i = 0; i < num_sats; ++i) {
            is >> sat[i];
            int num_pdrops;
            is >> num_pdrops;
            tmp_pd.resize(num_pdrops);
            for (int j = 0; j < num_pdrops; ++j) {
                is >> tmp_pd[j];
            }
            ap.appendRow(tmp_pd.begin(), tmp_pd.end());
        }
        all_pdrops.swap(ap);
        saturations.swap(sat);
    }




    /// Writes saturation and pressure drop data to an output stream.
    template <class Ostream>
    void writeControl(Ostream& os, const std::vector<double>& saturations, const Opm::SparseTable<double>& all_pdrops)
    {
        int num_sats = saturations.size();
        os << num_sats << '\n';
        for (int i = 0; i < num_sats; ++i) {
            os << saturations[i] << '\n';
            int num_pdrops = all_pdrops[i].size();
            os << num_pdrops;
            for (int j = 0; j < num_pdrops; ++j) {
                os << ' ' << all_pdrops[i][j];
            }
            os << '\n';
        }
    }




    template<class Ostream, class Tensor>
    void writeRelPerm(Ostream& os, const Tensor& K, double sat, double pdrop)
    {
        /* const int num_rows = K.numRows(); */
        /* const int num_cols = K.numCols(); */
	
        /* We write tensor output in Voigt notation, 
           (but we also output the remainding three terms)
           http://en.wikipedia.org/wiki/Voigt_notation
        */

	/*
	  TODO:
	  If fixed boundary condition, only output diagonal elements
	  If linear or periodic bc's, output all 9 elements (even though 6 is strictly necessary for periodic bc's)
	  Use the Tensor class to order output into Voigt notation, so that it works for more than 3x3 matrices
	*/
        
        os << pdrop << '\t';
        os << sat << '\t';
        os << K(0,0) << '\t'; /* xx */
        os << K(1,1) << '\t'; /* yy */
        os << K(2,2) << '\t'; /* zz */
        os << K(1,2) << '\t'; /* yz */
        os << K(0,2) << '\t'; /* xz */
        os << K(0,1) << '\t'; /* xy */
        os << K(2,1) << '\t'; /* zy */
        os << K(2,0) << '\t'; /* zx */
        os << K(1,2);         /* yz */
        

	os << std::endl;

    }



    template <class Traits>
    class SteadyStateUpscalerManager
    {
    public:
        void upscale(const Opm::parameter::ParameterGroup& param)
        {
            // Control structure.
            std::vector<double> saturations;
            Opm::SparseTable<double> all_pdrops;
            bool from_file = param.has("sat_pdrop_filename");
            if (from_file) {
                std::string filename = param.get<std::string>("sat_pdrop_filename");
                std::ifstream file(filename.c_str());
                if (!file) {
                    OPM_THROW(std::runtime_error, "Could not open file " << filename);
                }
                readControl(file, saturations, all_pdrops);
            } else {
                // Get a linear range of saturations.
                int num_sats = param.getDefault("num_sats", 4);
                double min_sat = param.getDefault("min_sat", 0.2);
                double max_sat = param.getDefault("max_sat", 0.8);
                saturations.resize(num_sats);
                for (int i = 0; i < num_sats; ++i) {
                    double factor = num_sats == 1 ? 0 : double(i)/double(num_sats - 1);
                    saturations[i] = (1.0 - factor)*min_sat + factor*max_sat;
                }
                // Get a logarithmic range of pressure drops.
                int num_pdrops = param.getDefault("num_pdrops", 5);
                double log_min_pdrop = std::log(param.getDefault("min_pdrop", 1e2));
                double log_max_pdrop = std::log(param.getDefault("max_pdrop", 1e6));
                std::vector<double> pdrops;
                pdrops.resize(num_pdrops);
                for (int i = 0; i < num_pdrops; ++i) {
                    double factor = num_pdrops == 1 ? 0 : double(i)/double(num_pdrops - 1);
                    pdrops[i] = std::exp((1.0 - factor)*log_min_pdrop + factor*log_max_pdrop);
                }
                // Assign the same pressure drops to all saturations.
                for (int i = 0; i < num_sats; ++i) {
                    all_pdrops.appendRow(pdrops.begin(), pdrops.end());
                }
            }
            int flow_direction = param.getDefault("flow_direction", 0);

            // Print the saturations and pressure drops.
            // writeControl(std::cout, saturations, all_pdrops);

            // Initialize upscaler.
            typedef SteadyStateUpscaler<Traits> Upscaler;
            typedef typename Upscaler::permtensor_t permtensor_t;
            Upscaler upscaler;
            upscaler.init(param);

            // First, compute an upscaled permeability.
            permtensor_t upscaled_K = upscaler.upscaleSinglePhase();
            permtensor_t upscaled_K_copy = upscaled_K;
            upscaled_K_copy *= (1.0/(Opm::prefix::milli*Opm::unit::darcy));
            std::cout.precision(15);
            std::cout << "Upscaled K in millidarcy:\n" << upscaled_K_copy << std::endl;
            std::cout << "Upscaled porosity: " << upscaler.upscalePorosity() << std::endl;

            // Create output streams for upscaled relative permeabilities
            std::string kr_filename = param.getDefault<std::string>("kr_filename", "upscaled_relperm");
            std::string krw_filename = kr_filename + "_water";
            std::string kro_filename = kr_filename + "_oil";
            std::ofstream krw_out(krw_filename.c_str());
            std::ofstream kro_out(kro_filename.c_str());
            krw_out << "# Result from steady state upscaling" << std::endl;
            krw_out << "# Pressuredrop  Sw  Krxx  Kryy  Krzz" << std::endl;
            kro_out << "# Result from steady state upscaling" << std::endl;
            kro_out << "# Pressuredrop  Sw  Krxx  Kryy  Krzz" << std::endl;


            krw_out.precision(15);  krw_out.setf(std::ios::scientific | std::ios::showpoint);
            kro_out.precision(15);  kro_out.setf(std::ios::scientific | std::ios::showpoint);
	    //#endif

            // Then, compute some upscaled relative permeabilities.
            int num_cells = upscaler.grid().size(0);
            int num_sats = saturations.size();
            for (int i = 0; i < num_sats; ++i) {
                // Starting every computation with a trio of uniform profiles.
                std::vector<double> init_sat(num_cells, saturations[i]);
                const Opm::SparseTable<double>::row_type pdrops = all_pdrops[i];
                int num_pdrops = pdrops.size();
                for (int j = 0; j < num_pdrops; ++j) {
                    double pdrop = pdrops[j];
                    std::pair<permtensor_t, permtensor_t> lambda
                        = upscaler.upscaleSteadyState(flow_direction, init_sat, saturations[i], pdrop, upscaled_K);
                    double usat = upscaler.lastSaturationUpscaled();
                    std::cout << "\n\nTensor of upscaled relperms for initial saturation " << saturations[i]
                              << ", real steady-state saturation " << usat
                              << " and pressure drop " << pdrop
                              << ":\n\n[water]\n" << lambda.first
                              << "\n[oil]\n" << lambda.second << std::endl;
                    // Changing initial saturations for next pressure drop to equal the steady state of the last
                    init_sat = upscaler.lastSaturationState();

		    
                    writeRelPerm(krw_out, lambda.first , usat, pdrop);
                    writeRelPerm(kro_out, lambda.second, usat, pdrop);
		    
                }
            }
        }
    };






} // namespace Opm

#endif // OPM_STEADYSTATEUPSCALERMANAGER_HEADER
