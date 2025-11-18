//==============================================================================
//!
//! \file elasticity.cpp
//!
//! \date Nov 9 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Elasticity helper class - template implementations
//!
//==============================================================================

#include <opm/elasticity/elasticity.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>

Dune::FieldVector<double,3>
Opm::Elasticity::waveSpeeds(const Dune::FieldMatrix<double,6,6>& C,
                            const double phi,
                            const double theta,
                            const double density)
{
  const double r = 1;
  Dune::FieldVector<double, 3> x;
  x[0] = r*std::cos(theta)*std::cos(phi);
  x[1] = r*std::sin(theta)*std::cos(phi);
  x[2] = r*std::sin(phi);

  Dune::FieldMatrix<double, 3, 6> D;
  D[0][0] = x[0];
  D[0][4] = x[2];
  D[0][5] = x[1];
  D[1][1] = x[1];
  D[1][3] = x[2];
  D[1][5] = x[0];
  D[2][2] = x[2];
  D[2][3] = x[1];
  D[2][4] = x[0];
  D /= x.two_norm();
  Dune::FieldMatrix<double, 6, 3> Dt;
  for (std::size_t i=0;i<6;++i)
    for (std::size_t j=0;j<3;++j)
      Dt[i][j] = D[j][i];

  Dune::FieldMatrix<double, 3, 6> T;
  Dune::FieldMatrix<double, 3, 3> E;
  Dune::FMatrixHelp::multMatrix(D, C, T);
  Dune::FMatrixHelp::multMatrix(T, Dt, E);
  Dune::FieldVector<double, 3> eigenvalues;
  Dune::FMatrixHelp::eigenValues(E, eigenvalues);
  std::sort(eigenvalues.begin(), eigenvalues.end(), std::greater<>{});
  for (std::size_t i=0;i<3;++i)
    eigenvalues[i] = std::sqrt(eigenvalues[i]/density);

  return eigenvalues;
}
