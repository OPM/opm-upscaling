//==============================================================================
//!
//! \file shapefunctions.hpp
//!
//! \date Nov 9 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for linear shape functions. Loosely based on code in dune-grid-howto 
//!
//==============================================================================
#pragma once

#include <dune/common/fvector.hh>

//! \brief Represents a linear shape function on a Q4/Q8 element
  template<class ctype, class rtype, int dim>
class LinearShapeFunction
{
  public:
    //! \brief The dimension of the shape function
    enum { dimension = dim };

    //! \brief Default constructor
    LinearShapeFunction() : coeff0(0.0), coeff1(0.0) {}

    //! \brief Construct a shape function with the given coefficients
    //! \param[in] coeff0_ The constant coefficients
    //! \param[in] coeff1_ The linear coefficients
    LinearShapeFunction(const Dune::FieldVector<rtype,dim>& coeff0_,
                        const Dune::FieldVector<rtype,dim>& coeff1_)
      : coeff0(coeff0_), coeff1(coeff1_) {}

    //! \brief Set the given conefficients
    //! \param[in] coeff0_ The constant coefficients
    //! \param[in] coeff1_ The linear coefficients
    void setCoeff(const Dune::FieldVector<rtype,dim>& coeff0_, const Dune::FieldVector<rtype,dim>& coeff1_)
    {
      coeff0 = coeff0_;
      coeff1 = coeff1_;
    }

    //! \brief Evaluate the shape function
    //! \param[in] local The local coordinates
    rtype evaluateFunction(const Dune::FieldVector<ctype,dim>& local) const
    {
      rtype result = 1;
      for (int i = 0; i < dim; ++i)
        result *= (coeff0[i]+coeff1[i] * local[i]);
      return result;
    }

    //! \brief Evaluate the gradient of the shape function
    //! \param[in] local The local coordinates
    Dune::FieldVector<rtype,dim>
    evaluateGradient(const Dune::FieldVector<ctype,dim>& local) const
    {
      Dune::FieldVector<rtype, dim> result;
      for (int i=0;i<dim;++i) {
        result[i] = 1;
        for (int j=0;j<dim;++j) {
          if (i == j)
            result[i] *= coeff1[j];
          else
            result[i] *= (coeff0[j]+coeff1[j]*local[j]);
        }
      }

      return result;
    }

  private:
    //! \brief The constant coefficients
    Dune::FieldVector<rtype,dim> coeff0;

    //! \brief The linear coefficients
    Dune::FieldVector<rtype,dim> coeff1;
};

//! \brief Singleton handler for the set of LinearShapeFunction
  template<class ctype, class rtype, int dim>
class P1ShapeFunctionSet
{
public:
    //! \brief The number of shape functions in the set
    enum { n = dim==2?4:8 };

    //! \brief A single shape function
    typedef LinearShapeFunction<ctype,rtype,dim> ShapeFunction;

    //! \brief The type of the return value from a shape function
    typedef rtype resulttype;

    //! \brief Get the only instance of this class
    static const P1ShapeFunctionSet& instance()
    {
      static const P1ShapeFunctionSet sfs;
      return sfs;
    }

    //! \brief Obtain a given shape function
    //! \param[in] i The requested shape function
    const ShapeFunction& operator[](int i) const
    {
      return f[i];
    }

private:
    //! \brief Private constructor prevents additional instances
    P1ShapeFunctionSet()
    {
        static rtype coeffs11[] = {0,
                                   1};
        static rtype coeffs12[] = {1,
                                   -1};

        static rtype coeffs21[] = { 1, 1,
                                    0, 1,
                                    1, 0,
                                    0, 0};
        static rtype coeffs22[] = {-1,-1, 
                                    1,-1, 
                                   -1, 1,
                                    1, 1};

        static rtype coeffs31[] = { 1, 1, 1,
                                    0, 1, 1,
                                    1, 0, 1,
                                    0, 0, 1,
                                    1, 1, 0,
                                    0, 1, 0,
                                    1, 0, 0,
                                    0, 0, 0};
        static rtype coeffs32[] = {-1,-1,-1, 
                                    1,-1,-1, 
                                   -1, 1,-1,
                                    1, 1,-1,
                                   -1,-1, 1,
                                    1,-1, 1,
                                   -1, 1, 1,
                                    1, 1, 1};

        rtype* coeffs1;
        rtype* coeffs2;
        if (dim == 1) {
          coeffs1 = coeffs11;
          coeffs2 = coeffs12;
        }
        if (dim == 2) {
          coeffs1 = coeffs21;
          coeffs2 = coeffs22;
        }
        if (dim == 3) {
          coeffs1 = coeffs31;
          coeffs2 = coeffs32;
        }

        Dune::FieldVector<rtype,dim> c1;
        Dune::FieldVector<rtype,dim> c2;
        for (int i=0;i<n;++i) {
          for (int j=0;j<dim;++j) {
            c1[j] = coeffs1[i*dim+j];
            c2[j] = coeffs2[i*dim+j];
          }
          f[i].setCoeff(c1,c2);
        }
    }

    //! \brief Our shape functions
    ShapeFunction f[n];
};
