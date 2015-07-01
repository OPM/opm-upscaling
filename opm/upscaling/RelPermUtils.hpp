/*
  Copyright 2010 Statoil ASA.

  This file is part of The Open Porous Media project (OPM).

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

/** @file RelPermutils.hpp
    @brief Helper utilities for relperm upscaling applications
 */

#ifndef OPM_UPSCALING_RELPERM_UTILS_HPP
#define OPM_UPSCALING_RELPERM_UTILS_HPP

#include <opm/core/utility/MonotCubicInterpolator.hpp>
#include <opm/upscaling/SinglePhaseUpscaler.hpp>
#include <array>
#include <vector>

namespace Opm {
  //! \brief Get value from tensor
  //! \param[in] K The tensor to extract the value from
  //! \param[in] voigt_idx The voigt index for value to extract (0..8)
  //! \return The requested value
  double getVoigtValue(const SinglePhaseUpscaler::permtensor_t& K, int voigt_idx);

  //! \brief Set value in tensor
  //! \param[out] K The tensor to set value in
  //! \param[in] voigt_idx The voigt index for value to set (0..8)
  //! \param[in] val Value to set in tensor
  void setVoigtValue(SinglePhaseUpscaler::permtensor_t& K, int voigt_idx, double val);

  //! \brief Helper class for relperm upscaling applications.
  class RelPermUpscaleHelper {
    public:
      bool isMaster;                               //!< Whether this is the master MPI node or not.
      bool doEclipseCheck;                         //!< Whether to check that input relperm curves include relperm at critical saturation points.
      bool anisotropic_input;                      //!< Whether input eclipse file has diagonal anisotrophy.
      double critRelpThresh;                       //!< Threshold for eclipse check of relative permeabilities.
      int points;                                  //!< Number of saturation points to upscale for.
      int tensorElementCount;                      //!< Number of independent elements in resulting tensor.
      std::vector<int> node_vs_pressurepoint;      //!< Distribution of pressure points to MPI nodes.
      bool upscaleBothPhases;                      //!< Whether to upscale both phases.
      std::vector<double> WaterSaturation;         //!< Re-upscaled water saturation for the computed pressure points.
      std::array<std::vector<std::vector<double>>,2> PhasePerm; //!< Permeability values per pressure point for each phase.
      SinglePhaseUpscaler::permtensor_t permTensor; //!< Tensor of upscaled results.
      SinglePhaseUpscaler::permtensor_t permTensorInv; //!< Inverted tensor of upscaled results.
      SinglePhaseUpscaler upscaler;                //!< The upscaler class.
      std::array<std::vector<double>,3> perms;     //!< 'PERM' values from input file.
      std::vector<double> poros;                   //!< Cell porosities
      std::vector<double> zcorns;                  //!< Cell heights.
      std::vector<int> satnums;                    //!< Cell satnums.
      double minSinglePhasePerm;                   //!< Minimum single phase permability value.
      std::vector<MonotCubicInterpolator> InvJfunctions; //!< Inverse of the loaded J-functions.
      std::vector<MonotCubicInterpolator> Krfunctions;  //!< Relperm-curves for phase 1 for each stone type
      std::vector<MonotCubicInterpolator> Krfunctions2; //!< Relperm-curves for phase 2 for each stone type
      SinglePhaseUpscaler::BoundaryConditionType boundaryCondition; //!< Boundary conditions to use.
      std::vector<MonotCubicInterpolator> SwPcfunctions; //!< Holds Sw(Pc) for each rocktype.
      std::string saturationstring; //!< Fluid system type.
      double Swir; //!< Upscaled saturation ir.
      double Swor; //!< Upscaled saturation max.
      double Pcmin; //!< Minimum capillary pressure.
      double Pcmax; //!< Maximum capillary pressure.
      std::vector<double> cellPoreVolumes; //!< Pore volume for each grid cell.
      size_t tesselatedCells; //!< Number of "active" cells (Sintef interpretation of "active")
      double volume; //!< Total volume.
      double poreVolume; //!< Total pore volume.
      MonotCubicInterpolator WaterSaturationVsCapPressure; //!< Water saturation as a function of capillary pressure.
      std::vector<double> pressurePoints; //!< Vector of capillary pressure points between Swor and Swir.

      RelPermUpscaleHelper() : permTensor(3,3,nullptr) {}

      //! \brief Collect results from all MPI nodes.
      void collectResults();
  };
}

#endif
