//==============================================================================
//!
//! \file elasticity.hpp
//!
//! \date Nov 9 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Elasticity helper class
//!
//==============================================================================
#ifndef ELASTICITY_HPP_
#define ELASTICITY_HPP_

#include <dune/common/fmatrix.hh>

namespace Opm {
namespace Elasticity {

//! \brief Elasticity helper class
  template<class GridType>
class Elasticity {
  public:
    //! \brief The dimension of our grid
    static const int dim = GridType::dimension;

    //! \brief A basic number
    typedef typename GridType::LeafGridView::ctype ctype;

    //! \brief Default constructor
    //! \param[in] gv_ The grid we are doing the calculations on
    Elasticity(const GridType& gv_) : gv(gv_) {}

    //! \brief Returns the B matrix in a quadrature point
    //! \param[in] point (Reference) coordinates of quadrature point 
    //! \param[in] Jinv Jacobian matrix in quadrature point
    //! \param[out] B The B matrix
      template<int funcdim>
    void getBVector(Dune::FieldVector<ctype,funcdim>& BVector,
                    const Dune::FieldVector<ctype,dim>& point);
  
      template<int components, int funcdim>
    void getBmatrix(Dune::FieldMatrix<ctype,components,funcdim>& B,
                    const Dune::FieldVector<ctype,dim>& point,
                    const Dune::FieldMatrix<ctype,dim,dim>& Jinv);

    //! \brief Return the stiffness matrix contributions in a quadrature point
    //! \param[in] B The B matrix in the quadrature point
    //! \param[in] C The constitutive matrix for the cell material
    //! \param[in] detJW Det J times quadrature weight
    //! \param[out] A The stiffness matrix contributions in the quadrature point
      template<int comp, int funcdim>
    void getStiffnessMatrix(Dune::FieldMatrix<ctype,funcdim,funcdim>& A,
                            const Dune::FieldMatrix<ctype,comp,funcdim>& B,
                            const Dune::FieldMatrix<ctype,comp,comp>& C,
                            ctype detJW);

    //! \brief Return the stress vector in a quadrature point
    //! \param[in] v The displacements in the quadrature point
    //! \param[in] eps0 The load case vector
    //! \param[in] B The B matrix in the quadrature point
    //! \param[in] C The constitutive matrix for the cell material
    //! \param[out] sigma The stress vector in the given quadrature point
      template<int comp, int funcdim>
    void getStressVector(Dune::FieldVector<ctype,comp>& sigma,
                         const Dune::FieldVector<ctype,funcdim>& v, 
                         const Dune::FieldVector<ctype,comp>& eps0,
                         const Dune::FieldMatrix<ctype,comp,funcdim>& B,
                         const Dune::FieldMatrix<ctype,comp,comp>& C);
  protected:
     //! \brief Const reference to our grid
    const GridType& gv;
};

//! \brief Compute the elastic wave velocities
//! \brief C The elastic tensor
//! \brief phi dip angle
//! \brief theta Azimuth angle
//! \brief density Density of material
Dune::FieldVector<double,3> waveSpeeds(const Dune::FieldMatrix<double,6,6>& C, double phi,
                                       double theta, double density);

}
}

#include "elasticity_impl.hpp"

#endif
