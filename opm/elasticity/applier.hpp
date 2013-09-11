//==============================================================================
//!
//! \file applier.hpp
//!
//! \date Dec 22 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Helper class for abstracting differences between inverse operators
//!        and preconditioners in DUNE
//!
//==============================================================================
#ifndef APPLIER_H_
#define APPLIER_H_

#include <opm/elasticity/matrixops.hpp>

namespace Opm {
  namespace Elasticity {

  template<class T>
struct OperatorApplier
{
  OperatorApplier(T& t) : A(t)
  {
  }

  void apply(Vector& v, Vector& d);

  void pre(Vector& x, Vector& b);

  void post(Vector& x);

  T& A;
};

typedef OperatorApplier<Dune::InverseOperator<Vector, Vector> > InverseApplier;
typedef OperatorApplier<Dune::Preconditioner<Vector, Vector> > PreApplier;

  template<>
void InverseApplier::apply(Vector& v, Vector& d)
{
  Dune::InverseOperatorResult r;
  A.apply(v, d, r);
}

  template<>
void PreApplier::apply(Vector& v, Vector& d)
{
  A.apply(v, d);
}

  template<>
void InverseApplier::pre(Vector& x, Vector& b)
{
}

  template<>
void PreApplier::pre(Vector& x, Vector& b)
{
  A.pre(x,b);
}

  template<>
void InverseApplier::post(Vector& x)
{
}

  template<>
void PreApplier::post(Vector& x)
{
  A.post(x);
}

}
}

#endif
