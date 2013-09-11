//==============================================================================
//!
//! \file mortar_evaluator.hpp
//!
//! \date Nov 15 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Linear operator for a Mortar block
//!
//==============================================================================
#ifndef MORTAR_EVALUATOR_HPP_
#define MORTAR_EVALUATOR_HPP_

#include <dune/istl/matrixmatrix.hh>
#include <dune/istl/solvers.hh>
#include <opm/elasticity/matrixops.hpp>
#include <opm/elasticity/mortar_utils.hpp>

namespace Opm {
namespace Elasticity {

/*! This implements a operator evaluation for the
 *  schur mortar-block S = B^T*A^-1*B
!*/

class MortarEvaluator : public Dune::LinearOperator<Vector, Vector> {
  public:
    // define the category
    enum {
      //! \brief The category the preconditioner is part of.
      category=Dune::SolverCategory::sequential
    };

    //! \brief Constructor
    //! \param[in] Ai Evaluator for A^-1
    //! \param[in] B The mortar coupling matrix
    MortarEvaluator(const Matrix& A_,
                    const Matrix& B_) :
      A(A_), B(B_)
    {
    }

    //! \brief Apply the multiplier block
    //! \param[in] x The vector to apply the operator to
    //! \param[out] y The result of the operator evaluation
    void apply(const Vector& x, Vector& y) const
    {
      Vector lambda, l2;

      MortarUtils::extractBlock(lambda, x, B.M(), A.N());
      A.mv(x, y);
      B.umv(lambda, y);
      l2.resize(lambda.size());
      B.mtv(x, l2);
      MortarUtils::injectBlock(y, l2, B.M(), A.N());
    }

    //! \brief Apply the multiplier block with an embedded axpy
    //! \param[in] alpha The scalar to scale with
    //! \param[in] x The vector to apply the operator to
    //! \param[out] y The result of the operator evaluation
    void applyscaleadd(field_type alpha, const Vector& x, Vector& y) const
    {
      Vector lambda, l2;

      A.usmv(alpha, x, y);
      MortarUtils::extractBlock(lambda, x, B.M(), A.N());
      B.umv(lambda, y);
      l2.resize(lambda.size());
      B.mtv(x, l2);
      for (size_t i=A.N(); i < y.size(); ++i)
        y[i] += alpha*l2[i-A.N()];
    }
  protected:
    const Matrix& A; //!< Reference to A matrix
    const Matrix& B; //!< Reference to the mortar coupling matrix
};

}
}

#endif
