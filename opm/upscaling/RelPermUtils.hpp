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

#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <opm/core/utility/MonotCubicInterpolator.hpp>

#include <opm/upscaling/ParserAdditions.hpp>
#include <opm/upscaling/SinglePhaseUpscaler.hpp>

#include <array>
#include <map>
#include <memory>
#include <tuple>
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
      int points;                                  //!< Number of saturation points to upscale for.
      int tensorElementCount;                      //!< Number of independent elements in resulting tensor.
      bool upscaleBothPhases;                      //!< Whether to upscale both phases.
      std::vector<double> WaterSaturation;         //!< Re-upscaled water saturation for the computed pressure points.
      SinglePhaseUpscaler::permtensor_t permTensor; //!< Tensor of upscaled results.
      std::vector<int> satnums;                    //!< Cell satnums.
      std::vector<MonotCubicInterpolator> InvJfunctions; //!< Inverse of the loaded J-functions.
      //! \brief Relperm-curves for each (component->phase->stone type)
      std::array<std::array<std::vector<MonotCubicInterpolator>,2>,3> Krfunctions;
      SinglePhaseUpscaler::BoundaryConditionType boundaryCondition; //!< Boundary conditions to use.
      std::vector<MonotCubicInterpolator> SwPcfunctions; //!< Holds Sw(Pc) for each rocktype.
      std::string saturationstring; //!< Fluid system type.
      size_t tesselatedCells; //!< Number of "active" cells (Sintef interpretation of "active")
      double volume; //!< Total volume.
      double poreVolume; //!< Total pore volume.
      std::vector<double> pressurePoints; //!< Vector of capillary pressure points between Swor and Swir.

      //! \brief Form Deck object from input file (ECLIPSE format)
      //!
      //! \tparam String String type for representing pathnames.  Typically
      //! \code std::string \endcode or \code const char* \endcode.
      //!
      //! \param[in] eclipseFileName Name of input ECLIPSE file.
      //!
      //! \return Deck object resulting from parsing input file.
      template <class String>
      static std::shared_ptr<Deck>
      parseEclipseFile(const String& eclipseFileName);

      //! \brief Default constructor.
      //! \param[in] mpi_rank Rank of this process (for parallel simulations).
      //! \param[in] options Options structure.
      //! \details Uses the following options: fluids
      RelPermUpscaleHelper(int mpi_rank, std::map<std::string,std::string>& options_);

      //! \brief Collect results from all MPI nodes.
      void collectResults();

      //! \brief Calculate relperm values from phase permeabilities.
      //! \param[in] phase The phase to calculate values for (0-indexed).
      //! \return The phase permeability tensor values.
      //! \details First index is voigt index, second index is pressure point.
      std::vector<std::vector<double>> getRelPerm(int phase) const;

      //! \brief Upscale single phase permeability
      void upscaleSinglePhasePermeability();

      //! \brief Do sanity checks for input file.
      //! \param[in] deck The deck to sanity check.
      //! \param[in] minPerm Minimum permeability.
      //! \param[in] maxPerm Maximum permeability.
      //! \param[in] minPoro Minimum porosity.
      //! \details Throws error string.
      void sanityCheckInput(Opm::DeckConstPtr deck,
                            const double      minPerm,
                            const double      maxPerm,
                            const double      minPoro);

      //! \brief Check that input relperm curevs specify critical saturations.
      void checkCriticalSaturations();

      //! \brief Setup requested boundary conditions.
      //! \details Uses the following options: bc
      void setupBoundaryConditions();

      //! \brief Tesselate grid
      //! \param[in] deck The grid to tesselate.
      //! \param[in] options Option structure.
      //! \details Uses the following options: linsolver_tolerance,
      //!          linsolver_verbosity, linsolver_type, linsolver_max_iterations,
      //!          linsolver_smooth_steps, linsolver_prolongate_factor, minPerm
      //! \return Time used for tesselation.
      double tesselateGrid(Opm::DeckConstPtr deck);

      //! \brief Find cell center pressure gradient for every cell.
      //! \details Uses the following options: gravity, waterDensity, oilDensity
      void calculateCellPressureGradients();

      //! \brief Calculate minimum and maximum capillary pressures.
      //! \details Uses the following options: maxPermContrast, minPerm,
      //!                                      gravity, linsolver_tolerance
      void calculateMinMaxCapillaryPressure();

      //! \brief Upscale capillary pressure.
      //! \details Uses the following options: saturationThreshold
      void upscaleCapillaryPressure();

      //! \brief Upscale permeabilities.
      //! \param[in] mpi_rank MPI rank of this process.
      //! \details Uses the following options:  minPerm, maxPermContrast
      //! \return Tuple with (total time, time per point).
      std::tuple<double,double> upscalePermeability(int mpi_rank);

    private:
      //! \brief Perform critical saturation check for a single curve.
      //! \param[in,out] func Function to check for.
      //! \details Checks that minimum relperm is larger than the threshold.
      //!          If not func is set to a constant value of 0.
      //! \return True if test passed, false otherwise.
      bool checkCurve(MonotCubicInterpolator& func);

      double Swir; //!< Upscaled saturation ir.
      double Swor; //!< Upscaled saturation max.
      double Pcmin; //!< Minimum capillary pressure.
      double Pcmax; //!< Maximum capillary pressure.
      double critRelpThresh; //!< Threshold for eclipse check of relative permeabilities.
      std::vector<int> node_vs_pressurepoint; //!< Distribution of pressure points to MPI nodes.
      std::array<std::vector<std::vector<double>>,2> PhasePerm; //!< Permeability values per pressure point for each phase.
      SinglePhaseUpscaler::permtensor_t permTensorInv; //!< Inverted tensor of upscaled results.
      SinglePhaseUpscaler upscaler; //!< The upscaler class.
      std::array<std::vector<double>,3> perms; //!< 'PERM' values from input file.
      std::vector<double> poros; //!< Cell porosities.
      std::vector<double> zcorns; //!< Cell heights.
      double minSinglePhasePerm; //!< Minimum single phase permability value.
      std::vector<double> cellPoreVolumes; //!< Pore volume for each grid cell.
      MonotCubicInterpolator WaterSaturationVsCapPressure; //!< Water saturation as a function of capillary pressure.

#if defined(UNITTEST_TRESPASS_PRIVATE_PROPERTY_DP)
    public: // Intrusive unit testing of calculcateCellPressureGradients()
#endif // UNITTEST_TRESPASS_PRIVATE_PROPERTY_DP
      std::vector<double> dP; //!<  Cell center pressure gradients due to gravity effects.

#if defined(UNITTEST_TRESPASS_PRIVATE_PROPERTY_DP)
    private:
#endif // UNITTEST_TRESPASS_PRIVATE_PROPERTY_DP

      std::map<std::string,std::string>& options; //!< Reference to options structure.
  };

  template <class String>
  std::shared_ptr<Deck>
  RelPermUpscaleHelper::parseEclipseFile(const String& eclipseFileName)
  {
    auto parser = Parser{};
    addNonStandardUpscalingKeywords(parser);

    return parser.parseFile(eclipseFileName, ParseContext{});
  }
}

#endif  // OPM_UPSCALING_RELPERM_UTILS_HPP
