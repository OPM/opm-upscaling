/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

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

#ifndef OPM_ROCK_HEADER_INCLUDED
#define OPM_ROCK_HEADER_INCLUDED

#include <opm/core/eclipse/EclipseGridParser.hpp>
#include <opm/porsol/common/Matrix.hpp>

#include <opm/porsol/common/ReservoirPropertyCommon.hpp>

namespace Opm
{

    /// @brief A property class for porous media rock.
    /// @tparam dim the dimension of the space, used for giving permeability tensors the right size.
    template <int dim>
    class Rock
    {
    public:
        /// @brief Tensor type for read-only access to permeability.
        typedef ImmutableCMatrix PermTensor;
        /// @brief Tensor type to be used for holding copies of permeability tensors.
        typedef OwnCMatrix       MutablePermTensor;
        /// @brief Tensor type for read and write access to permeability.
        typedef SharedCMatrix    SharedPermTensor;


        /// @brief Default constructor.
        Rock();

        /// @brief Initialize from a grdecl file.
        /// @param parser the parser holding the grdecl data.
	/// @param global_cell the mapping from cell indices to the logical
	///                    cartesian indices of the grdecl file.
 	/// @param perm_threshold lower threshold for permeability.
 	/// @param rock_list_filename if non-null, the referred string gives
        ///                           the filename for the rock list.
        /// @param use_jfunction_scaling if true, use j-function scaling of capillary
        ///                              pressure, if applicable.
        /// @param sigma interface tension for j-scaling, if applicable.
        /// @param theta angle for j-scaling, if applicable.
        void init(const Opm::EclipseGridParser& parser,
                  const std::vector<int>& global_cell,
                  const double perm_threshold = 0.0);

        /// @brief Initialize a uniform reservoir.
        /// @param num_cells number of cells in the grid.
        /// @param uniform_poro the uniform porosity.
        /// @param uniform_perm the uniform (scalar) permeability.
        void init(const int num_cells,
                  const double uniform_poro,
                  const double uniform_perm);


        /// @brief Read-access to porosity.
        /// @param cell_index index of a grid cell.
        /// @return porosity value of the cell.
        double porosity(int cell_index) const;

        /// @brief Read-access to permeability.
        /// @param cell_index index of a grid cell.
        /// @return permeability value of the cell.
        PermTensor permeability(int cell_index) const;

        /// @brief Read- and write-access to permeability. Use with caution.
        /// @param cell_index index of a grid cell.
        /// @return permeability value of the cell.
	SharedPermTensor permeabilityModifiable(int cell_index);

    protected:
	// Methods
        void assignPorosity(const Opm::EclipseGridParser& parser,
                            const std::vector<int>& global_cell);
        void assignPermeability(const Opm::EclipseGridParser& parser,
                                const std::vector<int>& global_cell,
                                const double perm_threshold);

	// Data members.
        std::vector<double>        porosity_;
        std::vector<double>        permeability_;
        std::vector<unsigned char> permfield_valid_;
        PermeabilityKind permeability_kind_;
    };


} // namespace Opm

#include "Rock_impl.hpp"


#endif // OPM_ROCK_HEADER_INCLUDED
