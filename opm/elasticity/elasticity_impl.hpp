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
}}
