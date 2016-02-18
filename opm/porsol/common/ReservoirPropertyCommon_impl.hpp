//===========================================================================
//
// File: ReservoirPropertyCommon_impl.hpp
//
// Created: Mon Oct 26 08:29:09 2009
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

#ifndef OPENRS_RESERVOIRPROPERTYCOMMON_IMPL_HEADER
#define OPENRS_RESERVOIRPROPERTYCOMMON_IMPL_HEADER

#include <opm/core/io/eclipse/EclipseGridInspector.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Deck/DeckKeyword.hpp>

#include <fstream>
#include <array>

namespace Opm
{

    namespace {

        /// @brief
        ///    Classify and verify a given permeability specification
        ///    from a structural point of view.  In particular, we
        ///    verify that there are no off-diagonal permeability
        ///    components such as @f$k_{xy}@f$ unless the
        ///    corresponding diagonal components are known as well.
        ///
        /// @param deck [in]
        ///    An Eclipse deck capable of answering which
        ///    permeability components are present in a given input
        ///    deck.
        ///
        /// @return
        ///    An enum value with the following possible values:
        ///        ScalarPerm     only one component was given.
        ///        DiagonalPerm   more than one component given.
        ///        TensorPerm     at least one cross-component given.
        ///        None           no components given.
        ///        Invalid        invalid set of components given.
        PermeabilityKind classifyPermeability(Opm::DeckConstPtr deck)
        {
            const bool xx = deck->hasKeyword("PERMX" );
            const bool xy = deck->hasKeyword("PERMXY");
            const bool xz = deck->hasKeyword("PERMXZ");

            const bool yx = deck->hasKeyword("PERMYX");
            const bool yy = deck->hasKeyword("PERMY" );
            const bool yz = deck->hasKeyword("PERMYZ");

            const bool zx = deck->hasKeyword("PERMZX");
            const bool zy = deck->hasKeyword("PERMZY");
            const bool zz = deck->hasKeyword("PERMZ" );

            int num_cross_comp = xy + xz + yx + yz + zx + zy;
            int num_comp       = xx + yy + zz + num_cross_comp;
            PermeabilityKind retval = None;
            if (num_cross_comp > 0) {
                retval = TensorPerm;
            } else {
                if (num_comp == 1) {
                    retval = ScalarPerm;
                } else if (num_comp >= 2) {
                    retval = DiagonalPerm;
                }
            }

            bool ok = true;
            if (num_comp > 0) {
                // At least one tensor component specified on input.
                // Verify that any remaining components are OK from a
                // structural point of view.  In particular, there
                // must not be any cross-components (e.g., k_{xy})
                // unless the corresponding diagonal component (e.g.,
                // k_{xx}) is present as well...
                //
                ok =        xx || !(xy || xz || yx || zx) ;
                ok = ok && (yy || !(yx || yz || xy || zy));
                ok = ok && (zz || !(zx || zy || xz || yz));
            }
            if (!ok) {
                retval = Invalid;
            }

            return retval;
        }


        /// @brief
        ///    Copy isotropic (scalar) permeability to other diagonal
        ///    components if the latter have not (yet) been assigned a
        ///    separate value.  Specifically, this function assigns
        ///    copies of the @f$i@f$ permeability component (e.g.,
        ///    'PERMX') to the @f$j@f$ and @f$k@f$ permeability (e.g.,
        ///    'PERMY' and 'PERMZ') components if these have not
        ///    previously been assigned.
        ///
        /// @param kmap
        ///    Permeability indirection map.  In particular @code
        ///    kmap[i] @endcode is the index (an integral number in
        ///    the set [1..9]) into the permeability tensor
        ///    representation of function @code fillTensor @endcode
        ///    which represents permeability component @code i
        ///    @endcode.
        ///
        /// @param [in] i
        /// @param [in] j
        /// @param [in] k
        void setScalarPermIfNeeded(std::array<int,9>& kmap,
                                   int i, int j, int k)
        {
            if (kmap[j] == 0) { kmap[j] = kmap[i]; }
            if (kmap[k] == 0) { kmap[k] = kmap[i]; }
        }


        /// @brief
        ///   Extract pointers to appropriate tensor components from
        ///   input deck.  The permeability tensor is, generally,
        ///   @code
        ///        [ kxx  kxy  kxz ]
        ///    K = [ kyx  kyy  kyz ]
        ///        [ kzx  kzy  kzz ]
        ///   @endcode
        ///   We store these values in a linear std::array using natural
        ///   ordering with the column index cycling the most rapidly.
        ///   In particular we use the representation
        ///   @code
        ///        [  0    1    2    3    4    5    6    7    8  ]
        ///    K = [ kxx, kxy, kxz, kyx, kyy, kyz, kzx, kzy, kzz ]
        ///   @endcode
        ///   Moreover, we explicitly enforce symmetric tensors by
        ///   assigning
        ///   @code
        ///     3     1       6     2       7     5
        ///    kyx = kxy,    kzx = kxz,    kzy = kyz
        ///   @endcode
        ///   However, we make no attempt at enforcing positive
        ///   definite tensors.
        ///
        /// @param [in]  deck
        ///    An Eclipse deck capable of answering which
        ///    permeability components are present in a given input
        ///    deck as well as retrieving the numerical value of each
        ///    permeability component in each grid cell.
        ///
        /// @param [out] tensor
        /// @param [out] kmap
        PermeabilityKind fillTensor(Opm::DeckConstPtr deck,
                                    std::vector<const std::vector<double>*>& tensor,
                                    std::array<int,9>&                     kmap)
        {
            PermeabilityKind kind = classifyPermeability(deck);
            if (kind == Invalid) {
                OPM_THROW(std::runtime_error, "Invalid set of permeability fields given.");
            }
            assert (tensor.size() == 1);
            for (int i = 0; i < 9; ++i) { kmap[i] = 0; }

            enum { xx, xy, xz,    // 0, 1, 2
                   yx, yy, yz,    // 3, 4, 5
                   zx, zy, zz };  // 6, 7, 8

            // -----------------------------------------------------------
            // 1st row: [kxx, kxy, kxz]
            if (deck->hasKeyword("PERMX" )) {
                kmap[xx] = tensor.size();
                tensor.push_back(&deck->getKeyword("PERMX").getSIDoubleData());

                setScalarPermIfNeeded(kmap, xx, yy, zz);
            }
            if (deck->hasKeyword("PERMXY")) {
                kmap[xy] = kmap[yx] = tensor.size();  // Enforce symmetry.
                tensor.push_back(&deck->getKeyword("PERMXY").getSIDoubleData());
            }
            if (deck->hasKeyword("PERMXZ")) {
                kmap[xz] = kmap[zx] = tensor.size();  // Enforce symmetry.
                tensor.push_back(&deck->getKeyword("PERMXZ").getSIDoubleData());
            }

            // -----------------------------------------------------------
            // 2nd row: [kyx, kyy, kyz]
            if (deck->hasKeyword("PERMYX")) {
                kmap[yx] = kmap[xy] = tensor.size();  // Enforce symmetry.
                tensor.push_back(&deck->getKeyword("PERMYX").getSIDoubleData());
            }
            if (deck->hasKeyword("PERMY" )) {
                kmap[yy] = tensor.size();
                tensor.push_back(&deck->getKeyword("PERMY").getSIDoubleData());

                setScalarPermIfNeeded(kmap, yy, zz, xx);
            }
            if (deck->hasKeyword("PERMYZ")) {
                kmap[yz] = kmap[zy] = tensor.size();  // Enforce symmetry.
                tensor.push_back(&deck->getKeyword("PERMYZ").getSIDoubleData());
            }

            // -----------------------------------------------------------
            // 3rd row: [kzx, kzy, kzz]
            if (deck->hasKeyword("PERMZX")) {
                kmap[zx] = kmap[xz] = tensor.size();  // Enforce symmetry.
                tensor.push_back(&deck->getKeyword("PERMZX").getSIDoubleData());
            }
            if (deck->hasKeyword("PERMZY")) {
                kmap[zy] = kmap[yz] = tensor.size();  // Enforce symmetry.
                tensor.push_back(&deck->getKeyword("PERMZY").getSIDoubleData());
            }
            if (deck->hasKeyword("PERMZ" )) {
                kmap[zz] = tensor.size();
                tensor.push_back(&deck->getKeyword("PERMZ").getSIDoubleData());

                setScalarPermIfNeeded(kmap, zz, xx, yy);
            }
            return kind;
        }

    } // anonymous namespace




    // ----- Methods of ReservoirPropertyCommon -----




    template <int dim, class RPImpl, class RockType>
    ReservoirPropertyCommon<dim, RPImpl, RockType>::ReservoirPropertyCommon()
#if 1
        : density1_  (1013.9*Opm::unit::kilogram/Opm::unit::cubic(Opm::unit::meter)),
          density2_  ( 834.7*Opm::unit::kilogram/Opm::unit::cubic(Opm::unit::meter)),
          viscosity1_(   1.0*Opm::prefix::centi*Opm::unit::Poise),
          viscosity2_(   3.0*Opm::prefix::centi*Opm::unit::Poise),
#else
        : density1_  (1000.0*Opm::unit::kilogram/Opm::unit::cubic(Opm::unit::meter)),
          density2_  (1000.0*Opm::unit::kilogram/Opm::unit::cubic(Opm::unit::meter)),
          viscosity1_(   1000.0*Opm::prefix::centi*Opm::unit::Poise),
          viscosity2_(   1000.0*Opm::prefix::centi*Opm::unit::Poise),
#endif
          permeability_kind_(Invalid)
    {
    }


    template <int dim, class RPImpl, class RockType>
    void ReservoirPropertyCommon<dim, RPImpl, RockType>::init(Opm::DeckConstPtr deck,
                                                              const std::vector<int>& global_cell,
                                                              const double perm_threshold,
                                                              const std::string* rock_list_filename,
                                                              const bool use_jfunction_scaling,
                                                              const double sigma,
                                                              const double theta)
    {
        // This code is mostly copied from ReservoirPropertyCommon::init(...).
        static_assert(dim == 3, "");

        permfield_valid_.assign(global_cell.size(),
                                std::vector<unsigned char>::value_type(0));

        assignPorosity    (deck, global_cell);
        assignNTG         (deck, global_cell);
        assignSWCR        (deck, global_cell);
        assignSOWCR       (deck, global_cell);
        assignPermeability(deck, global_cell, perm_threshold);
        assignRockTable   (deck, global_cell);

        if (rock_list_filename) {
            readRocks(*rock_list_filename);
        }

        // Added section. This is a hack, because not all rock classes
        // may care about J-scaling. They still have to implement
        // setUseJfunctionScaling() and setSigmaAndTheta(), though the
        // latter may throw if called for a rock where it does not make sense.
        int num_rocks = rock_.size();
        for (int i = 0; i < num_rocks; ++i) {
            rock_[i].setUseJfunctionScaling(use_jfunction_scaling);
            if (use_jfunction_scaling) {
                rock_[i].setSigmaAndTheta(sigma, theta);
            }
        }
        // End of added section.

        asImpl().computeCflFactors();
    }



    template <int dim, class RPImpl, class RockType>
    void ReservoirPropertyCommon<dim, RPImpl, RockType>::init(const int num_cells,
                                                              const double uniform_poro,
                                                              const double uniform_perm)
    {
        permfield_valid_.assign(num_cells, std::vector<unsigned char>::value_type(1));
        porosity_.assign(num_cells, uniform_poro);
        permeability_.assign(dim*dim*num_cells, 0.0);
        for (int i = 0; i < num_cells; ++i) {
            SharedPermTensor K = permeabilityModifiable(i);
            for (int dd = 0; dd < dim; ++dd) {
                K(dd, dd) = uniform_perm;
            }
        }
        cell_to_rock_.assign(num_cells, 0);
        asImpl().computeCflFactors();
    }

    template <int dim, class RPImpl, class RockType>
    void ReservoirPropertyCommon<dim, RPImpl, RockType>::setViscosities(double v1, double v2)
    {
        viscosity1_ = v1;
        viscosity2_ = v2;
    }

    template <int dim, class RPImpl, class RockType>
    void ReservoirPropertyCommon<dim, RPImpl, RockType>::setDensities(double d1, double d2)
    {
        density1_ = d1;
        density2_ = d2;
    }

    template <int dim, class RPImpl, class RockType>
    double ReservoirPropertyCommon<dim, RPImpl, RockType>::viscosityFirstPhase() const
    {
        return viscosity1_;
    }


    template <int dim, class RPImpl, class RockType>
    double ReservoirPropertyCommon<dim, RPImpl, RockType>::viscositySecondPhase() const
    {
        return viscosity2_;
    }


    template <int dim, class RPImpl, class RockType>
    double ReservoirPropertyCommon<dim, RPImpl, RockType>::densityFirstPhase() const
    {
        return density1_;
    }


    template <int dim, class RPImpl, class RockType>
    double ReservoirPropertyCommon<dim, RPImpl, RockType>::densitySecondPhase() const
    {
        return density2_;
    }


    template <int dim, class RPImpl, class RockType>
    double ReservoirPropertyCommon<dim, RPImpl, RockType>::porosity(int cell_index) const
    {
        return porosity_[cell_index];
    }

    template <int dim, class RPImpl, class RockType>
    double ReservoirPropertyCommon<dim, RPImpl, RockType>::ntg(int cell_index) const
    {
        return ntg_[cell_index];
    }

    template <int dim, class RPImpl, class RockType>
    double ReservoirPropertyCommon<dim, RPImpl, RockType>::swcr(int cell_index) const
    {
        return swcr_[cell_index];
    }

    template <int dim, class RPImpl, class RockType>
    double ReservoirPropertyCommon<dim, RPImpl, RockType>::sowcr(int cell_index) const
    {
        return sowcr_[cell_index];
    }


    template <int dim, class RPImpl, class RockType>
    typename ReservoirPropertyCommon<dim, RPImpl, RockType>::PermTensor
    ReservoirPropertyCommon<dim, RPImpl, RockType>::permeability(int cell_index) const
    {
        assert (permfield_valid_[cell_index]);

        const PermTensor K(dim, dim, &permeability_[dim*dim*cell_index]);
        return K;
    }


    template <int dim, class RPImpl, class RockType>
    typename ReservoirPropertyCommon<dim, RPImpl, RockType>::SharedPermTensor
    ReservoirPropertyCommon<dim, RPImpl, RockType>::permeabilityModifiable(int cell_index)
    {
        // Typically only used for assigning synthetic perm values.
        SharedPermTensor K(dim, dim, &permeability_[dim*dim*cell_index]);

        // Trust caller!
        permfield_valid_[cell_index] = std::vector<unsigned char>::value_type(1);

        return K;
    }


    template <int dim, class RPImpl, class RockType>
    template<class Vector>
    void ReservoirPropertyCommon<dim, RPImpl, RockType>::phaseDensities(int /*cell_index*/, Vector& density) const
    {
        assert (density.size() >= NumberOfPhases);
        density[0] = densityFirstPhase();
        density[1] = densitySecondPhase();
    }


    template <int dim, class RPImpl, class RockType>
    double ReservoirPropertyCommon<dim, RPImpl, RockType>::densityDifference() const
    {
        return density1_ - density2_;
    }


    template <int dim, class RPImpl, class RockType>
    double ReservoirPropertyCommon<dim, RPImpl, RockType>::cflFactor() const
    {
        return cfl_factor_;
    }


    template <int dim, class RPImpl, class RockType>
    double ReservoirPropertyCommon<dim, RPImpl, RockType>::cflFactorGravity() const
    {
        return cfl_factor_gravity_;
    }


    template <int dim, class RPImpl, class RockType>
    double ReservoirPropertyCommon<dim, RPImpl, RockType>::cflFactorCapillary() const
    {
        return cfl_factor_capillary_;
    }

    template <int dim, class RPImpl, class RockType>
    double ReservoirPropertyCommon<dim, RPImpl, RockType>::capillaryPressure(int cell_index, double saturation) const
        {
            if (rock_.size() > 0) {
                int r = cell_to_rock_[cell_index];
                return rock_[r].capPress(permeability(cell_index), porosity(cell_index), saturation);
            } else {
                // HACK ALERT!
                // Use zero capillary pressure if no known rock table exists.
                return 1e5*(1-saturation);
            }
    }

    template <int dim, class RPImpl, class RockType>
    double ReservoirPropertyCommon<dim, RPImpl, RockType>::capillaryPressureDeriv(int cell_index, double saturation) const
    {
        if (rock_.size() > 0) {
            int r = cell_to_rock_[cell_index];
            double dpc = rock_[r].capPressDeriv(permeability(cell_index), porosity(cell_index), saturation);
            return dpc;
        } else {
            // HACK ALERT!
            // Use zero capillary pressure if no known rock table exists.
            return -1.0e5;
        }
    }
    template <int dim, class RPImpl, class RockType>
    double ReservoirPropertyCommon<dim, RPImpl, RockType>::s_min(int cell_index) const
    {
        if (rock_.size() > 0) {
            int r = cell_to_rock_[cell_index];
            return rock_[r].s_min();
        } else {
            // HACK ALERT!
            // Use zero as minimum saturation if no known rock table exists.
            return 0;
        }
     }

    template <int dim, class RPImpl, class RockType>
    double ReservoirPropertyCommon<dim, RPImpl, RockType>::s_max(int cell_index) const
    {
        if (rock_.size() > 0) {
            int r = cell_to_rock_[cell_index];
            return rock_[r].s_max();
        } else {
            // HACK ALERT!
            // Use 1 as maximum saturation if no known rock table exists.
            return 1;
        }
    }

    template <int dim, class RPImpl, class RockType>
    double ReservoirPropertyCommon<dim, RPImpl, RockType>::saturationFromCapillaryPressure(int cell_index, double cap_press) const
    {
        if (rock_.size() > 0) {
            int r = cell_to_rock_[cell_index];
            return rock_[r].satFromCapPress(permeability(cell_index), porosity(cell_index), cap_press);
        } else {
            // HACK ALERT!
            // Use a zero saturation if no known rock table exists.
            return 0.0;
        }
    }


    template <int dim, class RPImpl, class RockType>
    void ReservoirPropertyCommon<dim, RPImpl, RockType>::writeSintefLegacyFormat(const std::string& grid_prefix) const
    {
        int num_cells = porosity_.size();
        // Write porosity.
        {
            std::string filename = grid_prefix + "-poro.dat";
            std::ofstream file(filename.c_str());
            if (!file) {
                OPM_THROW(std::runtime_error, "Could not open file " << filename);
            }
            file << num_cells << '\n';
            std::copy(porosity_.begin(), porosity_.end(), std::ostream_iterator<double>(file, "\n"));
        }
        // Write permeability.
        {
            std::string filename = grid_prefix + "-perm.dat";
            std::ofstream file(filename.c_str());
            if (!file) {
                OPM_THROW(std::runtime_error, "Could not open file " << filename);
            }
            file << num_cells << '\n';
            switch (permeability_kind_) {
            case TensorPerm:
                std::copy(permeability_.begin(), permeability_.end(), std::ostream_iterator<double>(file, "\n"));
                break;
            case DiagonalPerm:
                for (int c = 0; c < num_cells; ++c) {
                    int index = c*dim*dim;
                    for (int dd = 0; dd < dim; ++dd) {
                        file << permeability_[index + (dim + 1)*dd] << ' ';
                    }
                    file << '\n';
                }
                break;
            case ScalarPerm:
            case None: // Treated like a scalar permeability.
                for (int c = 0; c < num_cells; ++c) {
                    file << permeability_[c*dim*dim] << '\n';
                }
                break;
            default:
                OPM_THROW(std::runtime_error, "Cannot write invalid permeability.");
            }
        }
    }


    // ------ Private methods ------


    template <int dim, class RPImpl, class RockType>
    RPImpl& ReservoirPropertyCommon<dim, RPImpl, RockType>::asImpl()
    {
        return static_cast<RPImpl&>(*this);
    }




    template <int dim, class RPImpl, class RockType>
    void ReservoirPropertyCommon<dim, RPImpl, RockType>::assignPorosity(Opm::DeckConstPtr deck,
                                                                        const std::vector<int>& global_cell)
    {
        porosity_.assign(global_cell.size(), 1.0);

        if (deck->hasKeyword("PORO")) {
            Opm::EclipseGridInspector insp(deck);
            std::array<int, 3> dims = insp.gridSize();
            int num_global_cells = dims[0]*dims[1]*dims[2];
            const std::vector<double>& poro = deck->getKeyword("PORO").getSIDoubleData();
            if (int(poro.size()) != num_global_cells) {
                OPM_THROW(std::runtime_error, "PORO field must have the same size as the "
                      "logical cartesian size of the grid: "
                      << poro.size() << " != " << num_global_cells);
            }
            for (int c = 0; c < int(porosity_.size()); ++c) {
                porosity_[c] = poro[global_cell[c]];
            }
        }
    }

    template <int dim, class RPImpl, class RockType>
    void ReservoirPropertyCommon<dim, RPImpl, RockType>::assignNTG(Opm::DeckConstPtr deck,
                                                                   const std::vector<int>& global_cell)
    {
        ntg_.assign(global_cell.size(), 1.0);

        if (deck->hasKeyword("NTG")) {
            Opm::EclipseGridInspector insp(deck);
            std::array<int, 3> dims = insp.gridSize();
            int num_global_cells = dims[0]*dims[1]*dims[2];
            const std::vector<double>& ntg = deck->getKeyword("NTG").getSIDoubleData();
            if (int(ntg.size()) != num_global_cells) {
                OPM_THROW(std::runtime_error, "NTG field must have the same size as the "
                      "logical cartesian size of the grid: "
                      << ntg.size() << " != " << num_global_cells);
            }
            for (int c = 0; c < int(ntg_.size()); ++c) {
                ntg_[c] = ntg[global_cell[c]];
            }
        }
    }

    template <int dim, class RPImpl, class RockType>
    void ReservoirPropertyCommon<dim, RPImpl, RockType>::assignSWCR(Opm::DeckConstPtr deck,
                                                                    const std::vector<int>& global_cell)
    {
        swcr_.assign(global_cell.size(), 0.0);

        if (deck->hasKeyword("SWCR")) {
            Opm::EclipseGridInspector insp(deck);
            std::array<int, 3> dims = insp.gridSize();
            int num_global_cells = dims[0]*dims[1]*dims[2];
            const std::vector<double>& swcr =
                deck->getKeyword("SWCR").getSIDoubleData();
            if (int(swcr.size()) != num_global_cells) {
                OPM_THROW(std::runtime_error, "SWCR field must have the same size as the "
                      "logical cartesian size of the grid: "
                      << swcr.size() << " != " << num_global_cells);
            }
            for (int c = 0; c < int(swcr_.size()); ++c) {
                swcr_[c] = swcr[global_cell[c]];
            }
        }
    }

    template <int dim, class RPImpl, class RockType>
    void ReservoirPropertyCommon<dim, RPImpl, RockType>::assignSOWCR(Opm::DeckConstPtr deck,
                                                                     const std::vector<int>& global_cell)
    {
        sowcr_.assign(global_cell.size(), 0.0);

        if (deck->hasKeyword("SOWCR")) {
            Opm::EclipseGridInspector insp(deck);
            std::array<int, 3> dims = insp.gridSize();
            int num_global_cells = dims[0]*dims[1]*dims[2];
            const std::vector<double>& sowcr =
                deck->getKeyword("SOWCR").getSIDoubleData();
            if (int(sowcr.size()) != num_global_cells) {
                OPM_THROW(std::runtime_error, "SOWCR field must have the same size as the "
                      "logical cartesian size of the grid: "
                      << sowcr.size() << " != " << num_global_cells);
            }
            for (int c = 0; c < int(sowcr_.size()); ++c) {
                sowcr_[c] = sowcr[global_cell[c]];
            }
        }
    }

    template <int dim, class RPImpl, class RockType>
    void ReservoirPropertyCommon<dim, RPImpl, RockType>::assignPermeability(Opm::DeckConstPtr deck,
                                                                            const std::vector<int>& global_cell,
                                                                            double perm_threshold)
    {
        Opm::EclipseGridInspector insp(deck);
        std::array<int, 3> dims = insp.gridSize();
        int num_global_cells = dims[0]*dims[1]*dims[2];
        assert (num_global_cells > 0);

        permeability_.assign(dim * dim * global_cell.size(), 0.0);

        std::vector<const std::vector<double>*> tensor;
        tensor.reserve(10);

        const std::vector<double> zero(num_global_cells, 0.0);
        tensor.push_back(&zero);

        static_assert(dim == 3, "");
        std::array<int,9> kmap;
        permeability_kind_ = fillTensor(deck, tensor, kmap);
        for (int i = 1; i < int(tensor.size()); ++i) {
            if (int(tensor[i]->size()) != num_global_cells) {
                OPM_THROW(std::runtime_error, "All permeability fields must have the same size as the "
                      "logical cartesian size of the grid: "
                      << (tensor[i]->size()) << " != " << num_global_cells);
            }
        }

        // Assign permeability values only if such values are
        // given in the input deck represented by 'deck'.  In
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




    template <int dim, class RPImpl, class RockType>
    void ReservoirPropertyCommon<dim, RPImpl, RockType>::assignRockTable(Opm::DeckConstPtr deck,
                                                          const std::vector<int>& global_cell)
    {
        const int nc = global_cell.size();

        cell_to_rock_.assign(nc, 0);

        if (deck->hasKeyword("SATNUM")) {
            Opm::EclipseGridInspector insp(deck);
            std::array<int, 3> dims = insp.gridSize();
            int num_global_cells = dims[0]*dims[1]*dims[2];
            const std::vector<int>& satnum = deck->getKeyword("SATNUM").getIntData();
            if (int(satnum.size()) != num_global_cells) {
                OPM_THROW(std::runtime_error, "SATNUM field must have the same size as the "
                      "logical cartesian size of the grid: "
                      << satnum.size() << " != " << num_global_cells);
            }
            for (int c = 0; c < nc; ++c) {
                // Note: SATNUM is FORTRANish, ranging from 1 to n, therefore we subtract one.
                cell_to_rock_[c] = satnum[global_cell[c]] - 1;
            }
        }
        else if (deck->hasKeyword("ROCKTYPE")) {
            Opm::EclipseGridInspector insp(deck);
            std::array<int, 3> dims = insp.gridSize();
            int num_global_cells = dims[0]*dims[1]*dims[2];
            const std::vector<int>& satnum = deck->getKeyword("ROCKTYPE").getIntData();
            if (int(satnum.size()) != num_global_cells) {
                OPM_THROW(std::runtime_error, "ROCKTYPE field must have the same size as the "
                      "logical cartesian size of the grid: "
                      << satnum.size() << " != " << num_global_cells);
            }
            for (int c = 0; c < nc; ++c) {
                // Note: ROCKTYPE is FORTRANish, ranging from 1 to n, therefore we subtract one.
                cell_to_rock_[c] = satnum[global_cell[c]] - 1;
            }
        }
    }




    template <int dim, class RPImpl, class RockType>
    void ReservoirPropertyCommon<dim, RPImpl, RockType>::readRocks(const std::string& rock_list_file)
    {
        std::ifstream rl(rock_list_file.c_str());
        if (!rl) {
            OPM_THROW(std::runtime_error, "Could not open file " << rock_list_file);
        }
        int num_rocks = -1;
        rl >> num_rocks;
        assert(num_rocks >= 1);
        rock_.resize(num_rocks);
        std::string dir(rock_list_file.begin(), rock_list_file.begin() + rock_list_file.find_last_of('/') + 1);
        for (int i = 0; i < num_rocks; ++i) {
            std::string spec;
            while (spec.empty()) {
                std::getline(rl, spec);
            }
            rock_[i].read(dir, spec);
        }
    }




} // namespace Opm


#endif // OPENRS_RESERVOIRPROPERTYCOMMON_IMPL_HEADER
