//==============================================================================
//!
//! \file elasticity_impl.hpp
//!
//! \date Nov 9 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Elasticity helper class - template implementations
//!
//==============================================================================
#ifndef OPM_ELASTICITY_IMPL_HPP
#define OPM_ELASTICITY_IMPL_HPP

#include <opm/elasticity/shapefunctions.hpp>

namespace Opm {
namespace Elasticity {

  template<class GridType>
    template<int components, int funcdim>
void Elasticity<GridType>::getBmatrix(Dune::FieldMatrix<ctype,components,funcdim>& B,
                               const Dune::FieldVector<ctype,dim>& point,
                               const Dune::FieldMatrix<ctype,dim,dim>& Jinv)
{
    P1ShapeFunctionSet<ctype,ctype,dim> basis 
                             = P1ShapeFunctionSet<ctype,ctype,dim>::instance();
  int funcs = funcdim/dim;

  Dune::FieldMatrix<ctype,funcdim/dim,dim> dNdX;
  for (int i=0;i < funcs; i++)
    Jinv.mv(basis[i].evaluateGradient(point),dNdX[i]);
  static const int rows = 6+(dim-2)*12;

  Dune::FieldMatrix<ctype,rows,funcdim/dim> temp;
  temp = 0;

  if (dim == 3) {
#define INDEX(i,j) i+6*j
    for (int i=0; i < funcs; ++i) {
      // normal strain part
      temp[INDEX(0,0)][i] = dNdX[i][0];
      temp[INDEX(1,1)][i] = dNdX[i][1];
      temp[INDEX(2,2)][i] = dNdX[i][2];

      // shear strain part
      temp[INDEX(3,0)][i] = dNdX[i][1];
      temp[INDEX(3,1)][i] = dNdX[i][0];

      temp[INDEX(4,0)][i] = dNdX[i][2];
      temp[INDEX(4,2)][i] = dNdX[i][0];
      temp[INDEX(5,1)][i] = dNdX[i][2];
      temp[INDEX(5,2)][i] = dNdX[i][1];
    }
  } else if (dim == 2) {
#undef INDEX
#define INDEX(i,j) i+3*j
    for (int i=0; i < funcs; ++i) {
      // normal strain part
      temp[INDEX(0,0)][i] = dNdX[i][0];
      temp[INDEX(1,1)][i] = dNdX[i][1];

      // shear strain part
      temp[INDEX(2,0)][i] = dNdX[i][1];
      temp[INDEX(2,1)][i] = dNdX[i][0];
    }
  }
  size_t k=0;
  for (int j=0;j<funcs*dim;++j)
    for (size_t i=0;i<components;++i,++k)
      B[i][j] = temp[k%rows][k/rows];
}


    template<class GridType>
    template<int funcdims>
void Elasticity<GridType>::getBVector(Dune::FieldVector<ctype,funcdims>& Bvector,
                                      const Dune::FieldVector<ctype,dim>& point)
{
    P1ShapeFunctionSet<ctype,ctype,dim> basis 
        = P1ShapeFunctionSet<ctype,ctype,dim>::instance();

    Dune::FieldMatrix<ctype,funcdims,dim> N;
    for (int i=0;i < funcdims; i++){
        Bvector[i] = basis[i].evaluateFunction(point);
    }
}
    
  template <class GridType>
    template<int comp, int funcdim>
void Elasticity<GridType>::getStiffnessMatrix(
                                Dune::FieldMatrix<ctype,funcdim,funcdim>& A,
                                const Dune::FieldMatrix<ctype,comp,funcdim>& B,
                                const Dune::FieldMatrix<ctype,comp,comp>& C,
                                ctype detJW)
{
  Dune::FieldMatrix<ctype,funcdim,comp> B2;
  for (int i=0;i<comp;++i)
    for (int j=0;j<funcdim;++j)
      B2[j][i] = B[i][j];
  A = B.leftmultiplyany(C).leftmultiplyany(B2);
  A *= detJW;
}

  template<class GridType>
    template<int comp, int funcdim>
void Elasticity<GridType>::getStressVector(Dune::FieldVector<ctype,comp>& sigma,
                             const Dune::FieldVector<ctype,funcdim>& v,
                             const Dune::FieldVector<ctype,comp>& eps0,
                             const Dune::FieldMatrix<ctype,comp,funcdim>& B,
                             const Dune::FieldMatrix<ctype,comp,comp>& C)
{
  sigma = Dune::FMatrixHelp::mult(C,Dune::FMatrixHelp::mult(B,v)+eps0);
}

Dune::FieldVector<double,3> waveSpeeds(const Dune::FieldMatrix<double,6,6>& C, double phi, double theta, double density)
{
  const double r = 1;
  Dune::FieldVector<double, 3> x;
  x[0] = r*cos(theta)*cos(phi);
  x[1] = r*sin(theta)*cos(phi);
  x[2] = r*sin(phi);
  
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
  for (size_t i=0;i<6;++i)
    for (size_t j=0;j<3;++j)
      Dt[i][j] = D[j][i];

  Dune::FieldMatrix<double, 3, 6> T;
  Dune::FieldMatrix<double, 3, 3> E;
  Dune::FMatrixHelp::multMatrix(D, C, T);
  Dune::FMatrixHelp::multMatrix(T, Dt, E);
  Dune::FieldVector<double, 3> eigenvalues;
  Dune::FMatrixHelp::eigenValues(E, eigenvalues);
  std::sort(eigenvalues.begin(), eigenvalues.end(), std::greater<double>());
  for (size_t i=0;i<3;++i)
    eigenvalues[i] = sqrt(eigenvalues[i]/density);

  return eigenvalues;
}
}}

#endif
