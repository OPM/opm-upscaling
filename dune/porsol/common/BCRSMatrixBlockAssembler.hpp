/*===========================================================================
//
// File: BCRSMatrixBlockAssembler.hpp
//
// Created: 2011-10-02 14:04:38+0200
//
// Authors: Ingeborg S. Ligaarden <Ingeborg.Ligaarden@sintef.no>
//          Jostein R. Natvig     <Jostein.R.Natvig@sintef.no>
//          Halvor M. Nilsen      <HalvorMoll.Nilsen@sintef.no>
//          Atgeirr F. Rasmussen  <atgeirr@sintef.no>
//          Bård Skaflestad       <Bard.Skaflestad@sintef.no>
//
//==========================================================================*/


 /*
   Copyright 2011 SINTEF ICT, Applied Mathematics.
   Copyright 2011 Statoil ASA.

   This file is part of the Open Porous Media Project (OPM).

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

#ifndef OPM_BCRSMATRIXBLOCKASSEMBLER_HPP_HEADER
#define OPM_BCRSMATRIXBLOCKASSEMBLER_HPP_HEADER

#include <cassert>
#include <cstddef>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <dune/porsol/opmtransport/src/JacobianSystem.hpp>

namespace Opm {
    namespace ImplicitTransportDefault {
        typedef Dune::FieldMatrix<double, 1, 1> ScalarBlock;
        typedef Dune::BCRSMatrix<ScalarBlock>   ScalarBCRSMatrix;

        template <>
        class MatrixBlockAssembler<ScalarBCRSMatrix>
        {
        public:
            template <class Block>
            void
            assembleBlock(size_t ndof, size_t i, size_t j, const Block& b) {
                assert (ndof == 1);  (void) ndof;

                ScalarBCRSMatrix::block_type blk(b[0]);

                mat_[i][j] += blk;
            }

            template <class Connections>
            void
            createBlockRow(size_t i, const Connections& conn, size_t ndof) {
                assert (ndof == 1);  (void) ndof;

                ScalarBCRSMatrix::CreateIterator ci(mat_, i);

                for (typename Connections::const_iterator
                         c = conn.begin(), e = conn.end(); c != e; ++c) {
                    ci.insert(*c);
                }
            }

            void
            finalizeStructure() {}

            void
            setSize(size_t nrow, size_t ncol, size_t nnz = 0) {
                mat_.setSize     (nrow, ncol, nnz);
                mat_.setBuildMode(ScalarBCRSMatrix::row_wise);
            }

            const ScalarBCRSMatrix&
            Matrix() const { return mat_; }

        private:
            ScalarBCRSMatrix mat_;
        };
    }
 }

#endif  /* OPM_BCRSMATRIXBLOCKASSEMBLER_HPP_HEADER */
