//===========================================================================
//
// File: ReservoirPropertyCommon.hpp
//
// Created: Mon Oct 26 08:23:31 2009
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

#ifndef OPENRS_RESERVOIRPROPERTYCOMMON_HEADER
#define OPENRS_RESERVOIRPROPERTYCOMMON_HEADER

#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/porsol/common/Matrix.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>

namespace Opm
{



    /// @brief Enum for the kind of permeability field originally retrieved.
    enum PermeabilityKind { ScalarPerm, DiagonalPerm, TensorPerm, None, Invalid };




    /// @brief A property class for incompressible two-phase flow.
    /// @tparam dim the dimension of the space, used for giving permeability tensors the right size.
    template <int dim, class RPImpl, class RockType>
    class ReservoirPropertyCommon
    {
    public:
        /// @brief Tensor type for read-only access to permeability.
        typedef ImmutableCMatrix PermTensor;
        /// @brief Tensor type to be used for holding copies of permeability tensors.
        typedef OwnCMatrix       MutablePermTensor;
        /// @brief Tensor type for read and write access to permeability.
        typedef SharedCMatrix    SharedPermTensor;

        /// @brief The number of phases 
        enum { NumberOfPhases = 2 };

        /// @brief Default constructor.
        ReservoirPropertyCommon();

        /// @brief Initialize from a grdecl file.
        /// @param deck the deck holding the grdecl data.
        /// @param global_cell the mapping from cell indices to the logical
        ///                    cartesian indices of the grdecl file.
        /// @param perm_threshold lower threshold for permeability.
        /// @param rock_list_filename if non-null, the referred string gives
        ///                           the filename for the rock list.
        /// @param use_jfunction_scaling if true, use j-function scaling of capillary
        ///                              pressure, if applicable.
        /// @param sigma interface tension for j-scaling, if applicable.
        /// @param theta angle for j-scaling, if applicable.
        void init(const Opm::Deck&,
                  const std::vector<int>& global_cell,
                  const double perm_threshold = 0.0,
                  const std::string* rock_list_filename = 0,
                  const bool use_jfunction_scaling = true,
                  const double sigma = 1.0,
                  const double theta = 0.0);

        /// @brief Initialize a uniform reservoir.
        /// @param num_cells number of cells in the grid.
        /// @param uniform_poro the uniform porosity.
        /// @param uniform_perm the uniform (scalar) permeability.
        void init(const int num_cells,
                  const double uniform_poro = 0.2,
                  const double uniform_perm = 100.0*Opm::prefix::milli*Opm::unit::darcy);

        /// @brief Set viscosities of both faces.
        /// @param v1 the viscosity of the first (water) phase.
        /// @param v2 the viscosity of the second (oil) phase.
        void setViscosities(double v1, double v2);

        /// @brief Set densitities of both faces.
        /// @param d1 the densitity of the first (water) phase.
        /// @param d2 the densitity of the second (oil) phase.
        void setDensities(double d1, double d2);

	/// @brief Viscosity of first (water) phase.
	/// @return the viscosity value.
        double viscosityFirstPhase() const;

	/// @brief Viscosity of second (oil) phase.
	/// @return the viscosity value.
        double viscositySecondPhase() const;

	/// @brief Density of first (water) phase.
	/// @return the density value.
        double densityFirstPhase() const;

	/// @brief Density of second (oil) phase.
	/// @return the density value.
        double densitySecondPhase() const;

        /// @brief Read-access to porosity.
        /// @param cell_index index of a grid cell.
        /// @return porosity value of the cell.
        double porosity(int cell_index) const;

        /// @brief Read-access to ntg.
        /// @param cell_index index of a grid cell.
        /// @return ntg value of the cell.
        double ntg(int cell_index) const;

        /// @brief Read-access to swcr.
        /// @param cell_index index of a grid cell.
        /// @return swcr value of the cell.
        double swcr(int cell_index) const;

        /// @brief Read-access to sowcr.
        /// @param cell_index index of a grid cell.
        /// @return sowcr value of the cell.
        double sowcr(int cell_index) const;

        /// @brief Read-access to permeability.
        /// @param cell_index index of a grid cell.
        /// @return permeability value of the cell.
        PermTensor permeability(int cell_index) const;

        /// @brief Read- and write-access to permeability. Use with caution.
        /// @param cell_index index of a grid cell.
        /// @return permeability value of the cell.
	SharedPermTensor permeabilityModifiable(int cell_index);

        /// @brief Densities for both phases.
	/// @tparam Vector a class with size() and operator[].
        /// @param cell_index index of a grid cell (not used).
        /// @param[out] density the phase densities.
	///                     Expected to be of size 2 before (and after) the call.
        template<class Vector>
        void phaseDensities(int /*cell_index*/, Vector& density) const;

        /// @brief Difference of densities.
        /// @return densityFirstPhase() - densitySecondPhase()
        double densityDifference() const;

        /// @brief A factor useful in cfl computations.
        /// @return the cfl factor.
        double cflFactor() const;

        /// @brief A factor useful in gravity cfl computations.
        /// @return the gravity cfl factor.
        double cflFactorGravity() const;

        /// @brief A factor useful in gravity cfl computations.
        /// @return the capillary cfl factor.
        double cflFactorCapillary() const;

        /// @brief Capillary pressure.
        /// @param cell_index index of a grid cell.
	    /// @param saturation a saturation value.
        /// @return capillary pressure at the given cell and saturation.
        double capillaryPressure(int cell_index, double saturation) const;
        /// @brief Derivative of Capillary pressure.
        /// @param c index of a grid cell.
        /// @param s saturation value.
        /// @return capillary pressure at the given cell and saturation.
        double capillaryPressureDeriv(int c, double s) const;

        // @brief Minimum of saturation in rock table.
        /// @param c index of a grid cell.
        /// @return minimum saturation in given cell.
        double s_min(int c) const;

        // @brief Maximum of saturation in rock table.
        /// @param c index of a grid cell.
        /// @return maximum saturation in given cell.
        double s_max(int c) const;

        /// @brief Inverse of the capillary pressure function.
        /// @param cell_index index of a grid cell.
        /// @param cap_press a capillary pressure value.
        /// @return the saturation at the given cell has the given capillary pressure.
        double saturationFromCapillaryPressure(int cell_index, double cap_press) const;

        /// @brief Write permeability and porosity in the Sintef legacy format.
        /// @param grid_prefix the prefix of all files output by this function.
        void writeSintefLegacyFormat(const std::string& grid_prefix) const;

    protected:
	// Methods
        void assignPorosity(const Opm::Deck& deck,
                            const std::vector<int>& global_cell);
        void assignNTG(const Opm::Deck& deck,
                       const std::vector<int>& global_cell);
        void assignSWCR(const Opm::Deck& deck,
                        const std::vector<int>& global_cell);
        void assignSOWCR(const Opm::Deck& deck,
                         const std::vector<int>& global_cell);
        void assignPermeability(const Opm::Deck& deck,
                                const std::vector<int>& global_cell,
                                const double perm_threshold);
        void assignRockTable(const Opm::Deck& deck,
                             const std::vector<int>& global_cell);
        void readRocks(const std::string& rock_list_file);

	// Supporting Barton/Nackman trick (also known as the curiously recurring template pattern).
	RPImpl& asImpl();

	// Data members.
        std::vector<double>        porosity_;
        std::vector<double>        ntg_;
        std::vector<double>        swcr_;
        std::vector<double>        sowcr_;
        std::vector<double>        permeability_;
        std::vector<unsigned char> permfield_valid_;
        double density1_;
        double density2_;
        double viscosity1_;
        double viscosity2_;
        double cfl_factor_;
        double cfl_factor_gravity_;
        double cfl_factor_capillary_;
        std::vector<RockType> rock_;
        std::vector<int> cell_to_rock_;
        PermeabilityKind permeability_kind_;
    };


} // namespace Opm

#include "ReservoirPropertyCommon_impl.hpp"

#endif // OPENRS_RESERVOIRPROPERTYCOMMON_HEADER
