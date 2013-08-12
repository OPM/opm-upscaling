#ifndef UZAWA_SOLVER_HPP_
#define UZAWA_SOLVER_HPP_
//==============================================================================
//!
//! \file uzawa_solver.hpp
//!
//! \date Nov 9 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Uzawa scheme helper class
//!
//==============================================================================

#include <dune/istl/solvers.hh>

namespace Opm {
  namespace Elasticity {

  template<class X, class Y>
class UzawaSolver : public Dune::InverseOperator<X,Y>
{
  public:
    typedef std::shared_ptr<Dune::InverseOperator<X,Y> > OperatorPtr;
    UzawaSolver(OperatorPtr& innersolver_,
                OperatorPtr& outersolver_,
                const Matrix& B_) :
      innersolver(innersolver_), outersolver(outersolver_), B(B_)
    {
    }

    void apply(X& x, Y& b, double reduction, 
               Dune::InverseOperatorResult& res)
    {
      apply(x, b, res);
    }

    void apply(X& x, Y& b, Dune::InverseOperatorResult& res)
    {
      Vector lambda, lambda2, u, u2;
      MortarUtils::extractBlock(u, b, B.N());
      Dune::InverseOperatorResult res2;
      u2.resize(u.size());
      u2 = 0;
      innersolver->apply(u2, u, res2);
      lambda.resize(B.M());
      B.mtv(u2, lambda);
      lambda2.resize(lambda.size());
      lambda2 = 0;
      outersolver->apply(lambda2, lambda, res2);
      MortarUtils::extractBlock(u, b, B.N());
      B.usmv(-1.0, lambda2, u);
      u2 = 0;
      innersolver->apply(u2, u, res);
      MortarUtils::injectBlock(x, u2, B.N());
      MortarUtils::injectBlock(x, lambda2, B.M(), B.N());
    }
  protected:
    OperatorPtr innersolver;
    OperatorPtr outersolver;
    const Matrix& B;
};

}
}

#endif
