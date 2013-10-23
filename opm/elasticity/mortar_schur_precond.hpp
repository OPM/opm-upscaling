//==============================================================================
//!
//! \file schur_precond.hpp
//!
//! \date Nov 15 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Schur based preconditioner for saddle point systems
//!
//==============================================================================
#ifndef SCHUR_PRECOND_HPP_
#define SCHUR_PRECOND_HPP_

#include <dune/istl/preconditioners.hh>
#include <dune/istl/matrixmatrix.hh>
#include <opm/elasticity/matrixops.hpp>
#include <opm/elasticity/mortar_utils.hpp>
#include <opm/elasticity/elasticity_preconditioners.hpp>

namespace Opm {
namespace Elasticity {

/*! This implements a Schur-decomposition based preconditioner for the
 *  mortar-elasticity system
 *  [A   B]
 *  [B'   ]
 *  
 *  The preconditioner is
 *  [Apre B]
 *  [     P]
 *  Here Apre is some preconditioner for A and P some preconditioner for
 *  S = B^TA^-1B
!*/

  template<class PrecondElasticityBlock>
class MortarSchurPre : public Dune::Preconditioner<Vector,Vector> {
  public:
    // define the category
    enum {
      //! \brief The category the preconditioner is part of.
      category=Dune::SolverCategory::sequential
    };

    //! \brief Constructor
    //! \param[in] P The multiplier block with diagonal A approximation
    //! \param[in] B The mortar coupling matrix
    //! \param[in] Apre_ A preconfigured preconditioner for A
    //! \param[in] symmetric If true, use symmetric preconditioning
    MortarSchurPre(const Matrix& P_, const Matrix& B_,
                   PrecondElasticityBlock& Apre_, bool symmetric_=false) :
      Apre(Apre_), B(B_), N(B.N()), M(B.M()),
      Lpre(P_, false), symmetric(symmetric_)
    {
    }

    //! \brief Destructor
    virtual ~MortarSchurPre()
    {
    }

    //! \brief Preprocess preconditioner
    virtual void pre(Vector& x, Vector& b)
    {
      // extract displacement DOFs
      Vector tempx, tempb;
      MortarUtils::extractBlock(tempx, x, N);
      MortarUtils::extractBlock(tempb, b, N);
      Apre.pre(tempx, tempb);
      MortarUtils::injectBlock(x, tempx, N);
      MortarUtils::injectBlock(b, tempb, N);
    }

    //! \brief Applies the preconditioner
    //! \param[out] v The resulting vector
    //! \param[in] d The vector to apply the preconditioner to
    virtual void apply(Vector& v, const Vector& d)
    {
      // multiplier block (second row)
      Vector temp;
      MortarUtils::extractBlock(temp, d, M, N);
      Dune::InverseOperatorResult r;
      Vector temp2;
      temp2.resize(temp.size(), false);
      Lpre.apply(temp2, temp, r);
      MortarUtils::injectBlock(v, temp2, M, N);

      // elasticity block (first row)
      MortarUtils::extractBlock(temp, d, N);
      if (!symmetric)
        B.mmv(temp2, temp);

      temp2.resize(N, false);
      temp2 = 0;
      Apre.apply(temp2, temp);
      MortarUtils::injectBlock(v, temp2, N);
    }

    //! \brief Dummy post-process function
    virtual void post(Vector& x)
    {
      Apre.post(x);
    }
  protected:
    //! \brief The preconditioner for the elasticity operator
    PrecondElasticityBlock& Apre;

    //! \brief The mortar coupling matrix
    const Matrix& B;

    //! \brief Number of displacement DOFs
    int N;

    //! \brief Number of multiplier DOFs
    int M;

    //! \brief Linear solver for the multiplier block
    LUSolver Lpre;

    //! \brief Whether or not to use a symmetric preconditioner
    bool symmetric; 
};

}
}

#endif
