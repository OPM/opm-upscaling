//===========================================================================
//
// File: setupGridAndProps.hpp
//
// Created: Tue Aug 11 14:47:35 2009
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

#ifndef OPM_SETUPGRIDANDPROPS_HEADER
#define OPM_SETUPGRIDANDPROPS_HEADER

#include <opm/common/utility/parameters/ParameterGroup.hpp>
#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/grid/CpGrid.hpp>
#include <opm/porsol/common/ReservoirPropertyCapillary.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <opm/common/utility/platform_dependent/disable_warnings.h>

#include <opm/common/utility/FileSystem.hpp>

#include <opm/common/utility/platform_dependent/reenable_warnings.h>

namespace Opm
{

    /// Helper for determining whether we should
    template <class RP>
    bool useJ()
    {
        return false;
    }

    template<>
    bool useJ< ReservoirPropertyCapillary<3> >();

    /// @brief
    /// @todo Doc me!
    /// @param
    template <template <int> class ResProp>
    inline void setupGridAndProps(const Opm::ParameterGroup& param,
                                  Dune::CpGrid& grid,
                                  ResProp<3>& res_prop)
    {
        // Initialize grid and reservoir properties.
        // Parts copied from Dune::CpGrid::init().
        std::string fileformat = param.getDefault<std::string>("fileformat", "cartesian");
        if (fileformat == "eclipse") {
            std::string ecl_file = param.get<std::string>("filename");

            Opm::Parser parser;
            auto deck = parser.parseFile(ecl_file);
            if (param.has("z_tolerance")) {
                std::cerr << "****** Warning: z_tolerance parameter is obsolete, use PINCH in deck input instead\n";
            }
            bool periodic_extension = param.getDefault<bool>("periodic_extension", false);
            bool clip_z = param.getDefault<bool>("clip_z", false);
            bool turn_normals = param.getDefault<bool>("turn_normals", false);
            {
                Opm::EclipseGrid inputGrid(deck);
                grid.processEclipseFormat(&inputGrid, /* ecl_state = */ nullptr,
                                          periodic_extension,
                                          turn_normals,
                                          clip_z,
                                          /* pinchActive = */ true,
                                          /* edge_conformal = */ false);
            }
            // Save EGRID file in case we are writing ECL output.
            if (param.getDefault("output_ecl", false)) {
                OPM_THROW(std::runtime_error, "Saving to EGRID files is not yet implemented");
                /*
                Opm::filesystem::path ecl_path(ecl_file);
                const std::vector<int>& globalCell = grid.globalCell();
                ecl_path.replace_extension(".EGRID");
                parser.saveEGRID(ecl_path.string() , (int) globalCell.size() , &globalCell[0]);
                */
            }
            double perm_threshold_md = param.getDefault("perm_threshold_md", 0.0);
            double perm_threshold = Opm::unit::convert::from(perm_threshold_md, Opm::prefix::milli*Opm::unit::darcy);
            std::string rock_list = param.getDefault<std::string>("rock_list", "no_list");
            std::string* rl_ptr = (rock_list == "no_list") ? 0 : &rock_list;
            bool use_j = param.getDefault("use_jfunction_scaling", useJ<ResProp<3> >());
            double sigma = 1.0;
            double theta = 0.0;
            if (use_j) {
                sigma = param.getDefault("sigma", sigma);
                theta = param.getDefault("theta", theta);
            }
            if (param.has("viscosity1") || param.has("viscosity2")) {
                double v1 = param.getDefault("viscosity1", 0.001);
                double v2 = param.getDefault("viscosity2", 0.003);
                res_prop.setViscosities(v1, v2);
            }
            res_prop.init(deck, grid.globalCell(), perm_threshold, rl_ptr,
                          use_j, sigma, theta);
        } else if (fileformat == "cartesian") {
            std::array<int, 3> dims = {{ param.getDefault<int>("nx", 1),
                                    param.getDefault<int>("ny", 1),
                                    param.getDefault<int>("nz", 1) }};
            std::array<double, 3> cellsz = {{ param.getDefault<double>("dx", 1.0),
                                         param.getDefault<double>("dy", 1.0),
                                         param.getDefault<double>("dz", 1.0) }};
            grid.createCartesian(dims, cellsz);
            double default_poro = param.getDefault("default_poro", 0.2);
            double default_perm_md = param.getDefault("default_perm_md", 100.0);
            double default_perm = Opm::unit::convert::from(default_perm_md, Opm::prefix::milli*Opm::unit::darcy);
            OPM_MESSAGE("Warning: For generated cartesian grids, we use uniform reservoir properties.");
            res_prop.init(grid.size(0), default_poro, default_perm);
        } else {
            OPM_THROW(std::runtime_error,
                      "Unknown file format string: " + fileformat);
        }
        if (param.getDefault("use_unique_boundary_ids", false)) {
            grid.setUniqueBoundaryIds(true);
        }
    }

    /// @brief
    /// @todo Doc me!
    /// @param
    template <template <int> class ResProp>
    inline void setupGridAndPropsEclipse(const Opm::Deck& deck,
                                         bool periodic_extension,
                                         bool turn_normals,
                                         bool clip_z,
                                         bool unique_bids,
                                         double perm_threshold,
                                         const std::string& rock_list,
                                         bool use_jfunction_scaling,
                                         double sigma,
                                         double theta,
                                         Dune::CpGrid& grid,
                                         ResProp<3>& res_prop)
    {
        const Opm::EclipseGrid eg(deck);
        const std::string* rl_ptr = (rock_list == "no_list") ? 0 : &rock_list;

        grid.processEclipseFormat(&eg, /* ecl_state = */ nullptr,
                                  periodic_extension,
                                  turn_normals,
                                  clip_z,
                                  /* pinchActive = */ true,
                                  /* edge_conformal = */ false);

        res_prop.init(deck, grid.globalCell(), perm_threshold, rl_ptr, use_jfunction_scaling, sigma, theta);

        if (unique_bids) {
            grid.setUniqueBoundaryIds(true);
        }
    }

} // namespace Opm


#endif // OPENRS_SETUPGRIDANDPROPS_HEADER
