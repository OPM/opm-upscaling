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

#ifndef OPM_ROCK_IMPL_HEADER_INCLUDED
#define OPM_ROCK_IMPL_HEADER_INCLUDED


#include <fstream>
#include <boost/static_assert.hpp>
#include <boost/array.hpp>
#include <opm/core/io/eclipse/EclipseGridInspector.hpp>

namespace Opm
{

    // ----- Methods of Rock -----




    template <int dim>
    Rock<dim>::Rock()
        : permeability_kind_(Invalid)
    {
    }


    template <int dim>
    void Rock<dim>::init(const Opm::EclipseGridParser& parser,
                         const std::vector<int>& global_cell,
                         const double perm_threshold)
    {
        // This code is mostly copied from ReservoirPropertyCommon::init(...).
        BOOST_STATIC_ASSERT(dim == 3);

        permfield_valid_.assign(global_cell.size(), false);

        assignPorosity    (parser, global_cell);
        assignPermeability(parser, global_cell, perm_threshold);
    }


    template <int dim>
    void Rock<dim>::init(const int num_cells,
                         const double uniform_poro,
                         const double uniform_perm)
    {
        permfield_valid_.assign(num_cells, true);
        porosity_.assign(num_cells, uniform_poro);
        permeability_.assign(dim*dim*num_cells, 0.0);
        for (int i = 0; i < num_cells; ++i) {
            SharedPermTensor K = permeabilityModifiable(i);
            for (int dd = 0; dd < dim; ++dd) {
                K(dd, dd) = uniform_perm;
            }
        }
    }


    template <int dim>
    double Rock<dim>::porosity(int cell_index) const
    {
        return porosity_[cell_index];
    }


    template <int dim>
    typename Rock<dim>::PermTensor
    Rock<dim>::permeability(int cell_index) const
    {
        ASSERT (permfield_valid_[cell_index]);

        const PermTensor K(dim, dim, &permeability_[dim*dim*cell_index]);
        return K;
    }


    template <int dim>
    typename Rock<dim>::SharedPermTensor
    Rock<dim>::permeabilityModifiable(int cell_index)
    {
        // Typically only used for assigning synthetic perm values.
        SharedPermTensor K(dim, dim, &permeability_[dim*dim*cell_index]);

        // Trust caller!
        permfield_valid_[cell_index] = std::vector<unsigned char>::value_type(1);

        return K;
    }




    // ------ Private methods ------




    template <int dim>
    void Rock<dim>::assignPorosity(const Opm::EclipseGridParser& parser,
                                   const std::vector<int>& global_cell)
    {
        porosity_.assign(global_cell.size(), 1.0);

        if (parser.hasField("PORO")) {
            const std::vector<double>& poro = parser.getFloatingPointValue("PORO");

            for (int c = 0; c < int(porosity_.size()); ++c) {
                porosity_[c] = poro[global_cell[c]];
            }
        }
    }



    template <int dim>
    void Rock<dim>::assignPermeability(const Opm::EclipseGridParser& parser,
                                       const std::vector<int>& global_cell,
                                       double perm_threshold)
    {
	Opm::EclipseGridInspector insp(parser);
        std::tr1::array<int, 3> dims = insp.gridSize();
        int num_global_cells = dims[0]*dims[1]*dims[2];
        ASSERT (num_global_cells > 0);

        permeability_.assign(dim * dim * global_cell.size(), 0.0);

        std::vector<const std::vector<double>*> tensor;
        tensor.reserve(10);

        const std::vector<double> zero(num_global_cells, 0.0);
        tensor.push_back(&zero);

        BOOST_STATIC_ASSERT(dim == 3);
        boost::array<int,9> kmap;
        permeability_kind_ = fillTensor(parser, tensor, kmap);

        // Assign permeability values only if such values are
        // given in the input deck represented by 'parser'.  In
        // other words: Don't set any (arbitrary) default values.
        // It is infinitely better to experience a reproducible
        // crash than subtle errors resulting from a (poorly
        // chosen) default value...
        //
        if (tensor.size() > 1) {
            const int nc  = global_cell.size();
            int       off = 0;

            for (int c = 0; c < nc; ++c, off += dim*dim) {
                SharedPermTensor K(dim, dim, &permeability_[off]);
                int       kix  = 0;
                const int glob = global_cell[c];

                for (int i = 0; i < dim; ++i) {
                    for (int j = 0; j < dim; ++j, ++kix) {
                        K(i,j) = (*tensor[kmap[kix]])[glob];
                    }
                    K(i,i) = std::max(K(i,i), perm_threshold);
                }

                permfield_valid_[c] = std::vector<unsigned char>::value_type(1);
            }
        }
    }





} // namespace Opm


#endif // OPM_ROCK_IMPL_HEADER_INCLUDED
