//==============================================================================
//!
//! \file materials.cpp
//!
//! \date Oct 1 2007
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Material properties.
//!
//==============================================================================

#include "config.h"
#include "materials.hh"
#include <string.h>

namespace Opm {
namespace Elasticity {

std::ostream& Isotropic::write (std::ostream& os) const
{
  os <<"Isotropic material "<< this->num()
     <<": E = "<< E <<" nu = "<< nu;
  return os << std::endl;
}


/*!
  For edim = 3: \f[
  [C] = \frac{E}{(1+\nu)(1-2\nu)} \left[\begin{array}{cccccc}
  1-\nu & \nu & \nu \\
  \nu & 1-\nu & \nu \\
  \nu & \nu & 1-\nu \\
  & & & \frac{1}{2}-\nu \\
  & & & & \frac{1}{2}-\nu \\
  & & & & & \frac{1}{2}-\nu
  \end{array}\right] \f]

  For edim = 2 (plain strain): \f[
  [C] = \frac{E}{(1+\nu)(1-2\nu)} \left[\begin{array}{ccc}
  1-\nu & \nu \\
  \nu & 1-\nu \\
  & & \frac{1}{2}-\nu
  \end{array}\right] \f]
*/

bool Isotropic::getConstitutiveMatrix (Dune::FieldMatrix<double,6,6>& C,
                                       bool invers) const
{
  C = 0;
  const double one  = 1.f;
  const double two  = 2.f;
  const double half = 0.5f;
  if (nu < 0 || nu >= half) return false;

  const double fact = E / ((one + nu) * (one - nu - nu));

  C[0][0] = invers ? one / E : (one - nu) * fact;
  C[1][0] = invers ? -nu / E : nu * fact;
  C[2][0] = C[1][0];

  C[0][1] = C[1][0];
  C[1][1] = C[0][0];
  C[2][1] = C[1][0];

  C[0][2] = C[1][0];
  C[1][2] = C[1][0];
  C[2][2] = C[0][0];

  C[3][3] = invers ? (two + nu + nu) / E : (half - nu) * fact;
  C[4][4] = C[3][3];
  C[5][5] = C[3][3];

  return true;
}

bool Isotropic::getConstitutiveMatrix (Dune::FieldMatrix<double,3,3>& C,
                                       bool invers) const
{
  const double one  = 1.f;
  const double two  = 2.f;
  const double half = 0.5f;
  if (nu < 0 || nu >= half) return false;

  const double fact = E / ((one + nu) * (one - nu - nu));

  C[0][0] = invers ? (one - nu*nu) / E : (one - nu) * fact;
  C[1][0] = invers ? (-nu - nu*nu) / E : nu * fact;

  C[0][1] = C[1][0];
  C[1][1] = C[0][0];

  C[2][2] = invers ? (two + nu + nu) / E : (half - nu) * fact;

  return true;
}


OrthotropicD::OrthotropicD (int ID, double Ex, double Ey, double Ez,
                            double Gxy, double Gxz, double Gyz) : Material(ID)
{
  E[0] = Ex;
  E[1] = Ey;
  E[2] = Ez;
  E[3] = Gxy;
  E[4] = Gxz > double(0) ? Gxz : E[3];
  E[5] = Gyz > double(0) ? Gyz : E[4];
}


std::ostream& OrthotropicD::write (std::ostream& os) const
{
  os <<"Orthotropic material "<< this->num()
     <<": Ex = "<< E[0] <<" Ey = "<< E[1] <<" Ez = "<< E[2];
  if (E[4] == E[3] && E[5] == E[3])
    os <<" G = "<< E[3];
  else
    os <<" Gxy = "<< E[3] <<" Gxz = "<< E[4] <<" Gyz = "<< E[5];
  return os << std::endl;
}


double OrthotropicD::getPar (int ipar) const
{
  if (ipar > 0 && ipar <= 6)
    return E[ipar-1];
  else
    return double(0);
}


/*!
  \f[ [C] = \left[\begin{array}{cccccc}
  E_x \\ & E_y \\ & & E_z \\
  & & & G_{xy} \\ & & & & G_{xz} \\ & & & & & G_{yz}
  \end{array}\right] \f]
*/

bool OrthotropicD::getConstitutiveMatrix (Dune::FieldMatrix<double,6,6>& C,
                                          bool invers) const
{
  for (size_t i = 0; i < 6; i++)
    C[i][i] = invers ? double(1)/E[i] : E[i];

  return true;
}

bool OrthotropicD::getConstitutiveMatrix (Dune::FieldMatrix<double,3,3>& /* C */,
                                          bool /* invers */) const
{
  return false;
}


OrthotropicSym::OrthotropicSym (int ID,
                        const Dune::DynamicVector<double>& Cu) : Material(ID)
{
  // The input matrix is assumed with Voigt ordering (xx,yy,zz,yz,zx,xy)
  // Have to reorder to the internal order required  (xx,yy,zz,xy,xz,yz)
  memcpy(Cupper,&Cu[0],21*sizeof(double));
  std::swap(Cupper[ 3],Cupper[ 5]);
  std::swap(Cupper[ 8],Cupper[10]);
  std::swap(Cupper[12],Cupper[14]);
  std::swap(Cupper[15],Cupper[20]);
  std::swap(Cupper[16],Cupper[19]);
}


std::ostream& OrthotropicSym::write (std::ostream& os) const
{
  //TODO: Print matrix
  return os;
}


double OrthotropicSym::getPar (int ipar) const
{
  if (ipar > 0 && ipar <= 21)
    return Cupper[ipar-1];
  else
    return double(0);
}


/*!
  \f[ [C] = \left[\begin{array}{cccccc}
    C_{11} & C_{12} & C_{13} & C_{14} & C_{15} & C_{16} \\
           & C_{22} & C_{23} & C_{24} & C_{25} & C_{26} \\
           &        & C_{33} & C_{34} & C_{35} & C_{36} \\
           &        &        & C_{44} & C_{45} & C_{46} \\
    \multicolumn{3}{c}{\rm symm.} &   & C_{55} & C_{56} \\
           &        &        &        &        & C_{66}
  \end{array}\right] \f]
*/

bool OrthotropicSym::getConstitutiveMatrix (Dune::FieldMatrix<double,6,6>& C,
                                            bool invers) const
{
  size_t i, j, k = 0;
  for (i = 0; i < 6; i++)
  {
    C[i][i] = Cupper[k++];
    for (j = i+1; j < 6; j++)
      C[i][j] = C[j][i] = Cupper[k++];
  }

  if (invers) {
    //TODO!
    assert(0);
  }

  return true;
}

bool OrthotropicSym::getConstitutiveMatrix (Dune::FieldMatrix<double,3,3>& /* C */,
                                            bool /* invers */) const
{
  return false;
}

}
}
