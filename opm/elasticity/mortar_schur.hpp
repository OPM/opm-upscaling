//==============================================================================
//!
//! \file mortar_schur.hpp
//!
//! \date Nov 15 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Schur based linear operator for a Mortar block
//!
//==============================================================================
#ifndef MORTAR_SCHUR_HPP_
#define MORTAR_SCHUR_HPP_

#include <dune/istl/matrixmatrix.hh>
#include <opm/elasticity/applier.hpp>
#include <opm/elasticity/matrixops.hpp>

namespace Opm {
namespace Elasticity {

/*! This implements a operator evaluation for the
 *  schur mortar-block S = B^T*A^-1*B
!*/

  template<class T>
class MortarBlockEvaluator : public Dune::LinearOperator<Vector, Vector> {
  public:
    // define the category
    enum {
      //! \brief The category the preconditioner is part of.
      category=Dune::SolverCategory::sequential
    };

    //! \brief Constructor
    //! \param[in] Ai Solver or preconditioner for A^-1
    //! \param[in] B The mortar coupling matrix
    MortarBlockEvaluator(T& Ai_,
                         const Matrix& B_) :
      Ai(Ai_), B(B_), op(Ai)
    {
    }

    //! \brief Apply the multiplier block
    //! \param[in] x The vector to apply the operator to
    //! \param[out] y The result of the operator evaluation
    void apply(const Vector& x, Vector& y) const
    {
      y = 0;
      applyscaleadd(1.0, x, y);
    }

    //! \brief Apply the multiplier block with an embedded axpy
    //! \param[in] alpha The scalar to scale with
    //! \param[in] x The vector to apply the operator to
    //! \param[out] y The result of the operator evaluation
    void applyscaleadd(field_type alpha, const Vector& x, Vector& y) const
    {
      Vector temp1(B.N());
      B.mv(x, temp1);
      Vector temp(temp1.size());
      op.pre(temp, temp1);
      Dune::InverseOperatorResult res;
      temp = 0;
      op.pre(temp, temp1);
      op.apply(temp, temp1);
      op.post(temp);
      B.usmtv(alpha, temp, y);
    }
  protected:
    T& Ai;            //!< Reference to solver or evaluator for inverse operator
    const Matrix& B;  //!< Reference to the mortar coupling matrix
    mutable OperatorApplier<T> op; //!< Applier for the preconditioner / inverse solver
};

}
}

#endif
