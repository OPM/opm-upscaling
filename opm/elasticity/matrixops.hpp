//==============================================================================
//!
//! \file matrixops.hpp
//!
//! \date Nov 9 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Helper class with some matrix operations
//!
//==============================================================================
#ifndef MATRIXOPS_HPP_
#define MATRIXOPS_HPP_

#include <dune/common/fmatrix.hh> 
#include <dune/common/dynmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>

namespace Opm {
namespace Elasticity {

//! \brief A sparse matrix holding our operator
typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > Matrix;

//! \brief For storing matrix adjacency/sparsity patterns
typedef std::vector< std::set<int> > AdjacencyPattern;

//! \brief A vector holding our RHS
typedef Dune::BlockVector<Dune::FieldVector<double,1> > Vector;

//! \brief Helper class with some matrix operations
class MatrixOps {
  public:
    //! \brief Create a sparse matrix from a given adjacency pattern
    //! \param[in] adj The adjacency pattern
    //! \param[in] rows The number of rows in the matrix
    //! \param[in] cols The number of columns in the matrix
    //! \param[out] A The created matrix
    static void fromAdjacency(Matrix& A, const AdjacencyPattern& adj,
                              int rows, int cols);

    //! \brief Create a sparse matrix from a dense matrix
    //! \param[in] T The dense matrix
    //! \returns The sparse matrix
    static Matrix fromDense(const Dune::DynamicMatrix<double>& T);

    //! \brief Print a matrix to stdout
    //! \param[in] A The matrix to print
    static void print(const Matrix& A);

    //! \brief axpy like operation - returns A+alpha*B
    //! \param[in] A The matrix to subtract from
    //! \param[in] B The matrix to subtract
    //! \param[in] alpha The constant in front of B
    //! \returns A+alpha*B
    static Matrix Axpy(const Matrix& A, const Matrix& B, double alpha);

    //! \brief Augment a matrix with another
    //! \param[in] A The matrix to be augmented
    //! \param[in] B The matrix to augment with
    //! \param[in] r0 The starting row of the augment matrix
    //! \param[in] c0 The starting column of the augment matrix
    //! \param[in] symmetric If true, augment symmetrically
    static Matrix augment(const Matrix& A, const Matrix& B,
                          size_t r0, size_t c0, bool symmetric);

    //! \brief Extract the diagonal of a matrix into a new matrix
    //! \param[in] The matrix to extract the diagonal from
    //! \returns M = diag(A)
    static Matrix extractDiagonal(const Matrix& A);

    //! \brief Returns a diagonal matrix
    //! \param[in] N The dimension of the matrix
    static Matrix diagonal(size_t N);

    //! \brief Extract a subblock of a matrix into a new matrix
    //! \param[in] The matrix to extract from
    //! \returns The subblock
    static Matrix extractBlock(const Matrix& A,
                               size_t r0, size_t N, size_t c0, size_t M);

    //! \brief Save a matrix as a dense asc file
    //! \param[in] A The matrix to save
    //! \param[in] file File name
    //! \details This is only useful for debugging as the files grow very big
    static void saveAsc(const Matrix& A, const std::string& file);
};

}
}

#endif
