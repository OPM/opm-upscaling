/*===========================================================================
//
// File: Dune::BCRSMatrixBlockAssembler.hpp
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

#include <opm/core/transport/JacobianSystem.hpp>

namespace Opm {
    namespace ImplicitTransportDefault {
        namespace ISTLTypeDetails {
            typedef Dune::FieldVector<double, 1>    ScalarVectorBlockType;
            typedef Dune::FieldMatrix<double, 1, 1> ScalarMatrixBlockType;

            typedef Dune::BlockVector<ScalarVectorBlockType> ScalarBlockVector;
            typedef Dune::BCRSMatrix <ScalarMatrixBlockType> ScalarBCRSMatrix;
        }

        template <>
        class VectorAdder<ISTLTypeDetails::ScalarBlockVector> {
        public:
            static void
            add(const ISTLTypeDetails::ScalarBlockVector& x,
                ISTLTypeDetails::ScalarBlockVector&       y) {
                y += x;
            }
        };

        template <>
        class VectorNegater<ISTLTypeDetails::ScalarBlockVector> {
        public:
            static void
            negate(ISTLTypeDetails::ScalarBlockVector& x) {
                x *= -1.0;
            }
        };

        template <>
        class VectorZero<ISTLTypeDetails::ScalarBlockVector> {
        public:
            static void
            zero(ISTLTypeDetails::ScalarBlockVector& x) {
                x = 0.0;
            }
        };

        template <>
        class VectorBlockAssembler<ISTLTypeDetails::ScalarBlockVector> {
        public:
            template <class Block>
            static void
            assemble(::std::size_t                       ndof,
                     ::std::size_t                       i   ,
                     const Block&                        b   ,
                     ISTLTypeDetails::ScalarBlockVector& vec ) {
                assert (ndof == 1);  (void) ndof;

                ISTLTypeDetails::ScalarBlockVector::block_type blk(b[0]);

                vec[i] += blk;
            }
        };

        template <>
        class VectorAssign<ISTLTypeDetails::ScalarBlockVector> {
        public:
            // y <- x
            static void
            assign(const ISTLTypeDetails::ScalarBlockVector& x,
                   ISTLTypeDetails::ScalarBlockVector&       y) {
                y = x;
            }

            // y <- a*x
            template <class Scalar>
            static void
            assign(const Scalar& a,
                   const ISTLTypeDetails::ScalarBlockVector& x,
                   ISTLTypeDetails::ScalarBlockVector&       y) {
                y  = x;
                y *= a;
            }
        };

        template <>
        class MatrixZero<ISTLTypeDetails::ScalarBCRSMatrix> {
        public:
            static void
            zero(ISTLTypeDetails::ScalarBCRSMatrix& A) {
                A = 0.0;
            }
        };

        template <>
        class MatrixBlockAssembler<ISTLTypeDetails::ScalarBCRSMatrix>
        {
        public:
            template <class Block>
            void
            assembleBlock(size_t ndof, size_t i, size_t j, const Block& b) {
                assert (ndof == 1);  (void) ndof;

                ISTLTypeDetails::ScalarBCRSMatrix::block_type blk(b[0]);

                mat_[i][j] += blk;
            }

            template <class Connections>
            void
            createBlockRow(size_t i, const Connections& conn, size_t ndof) {
                assert (ndof == 1);            (void) ndof;
                assert (i    == i_prev_ + 1);  (void) i   ;

                ISTLTypeDetails::ScalarBCRSMatrix::CreateIterator ci(mat_, i);

                for (typename Connections::const_iterator
                         c = conn.begin(), e = conn.end(); c != e; ++c) {
                    ci.insert(*c);
                }

                ++i_prev_;
                ++ci;
            }

            void
            finalizeStructure() {}

            void
            setSize(::std::size_t ndof,
                    ::std::size_t nrow,
                    ::std::size_t ncol,
                    ::std::size_t nnz = 0) {

                assert (ndof == 1);  (void) ndof;
                (void) nnz;

                mat_.setSize     (nrow, ncol);
                mat_.setBuildMode(ISTLTypeDetails::ScalarBCRSMatrix::row_wise);

                i_prev_ = -1;
            }

            const ISTLTypeDetails::ScalarBCRSMatrix&
            matrix() const { return mat_; }

            ISTLTypeDetails::ScalarBCRSMatrix&
            matrix() { return mat_; }

        private:
            ::std::size_t                     i_prev_;
            ISTLTypeDetails::ScalarBCRSMatrix mat_   ;
        };
    }
}

#endif  /* OPM_BCRSMATRIXBLOCKASSEMBLER_HPP_HEADER */
