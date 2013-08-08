//==============================================================================
//!
//! \file shapefunctions.hpp
//!
//! \date Nov 9 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Classes for shape functions. Loosely based on code in dune-grid-howto
//!
//==============================================================================
#ifndef SHAPEFUNCTIONS_HPP_
#define SHAPEFUNCTIONS_HPP_

#include <dune/common/fvector.hh>
#include <dune/common/dynmatrixev.hh>

#include <complex>

namespace Opm {
namespace Elasticity {

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

//! \brief Represents a cardinal function on a line
  template<class ctype, class rtype>
class LagrangeCardinalFunction
{
  public:
    //! \brief Empty default constructor
    LagrangeCardinalFunction() {}

    //! \brief Construct a cardinal function with the given nodes
    //! \param[in] nodes_ The nodes
    //! \param[in] i The node this function is associated with
    LagrangeCardinalFunction(const std::vector<rtype>& nodes_,
                             size_t i)
      : nodes(nodes_), node(i) {}

    //! \brief Evaluate the shape function
    //! \param[in] local The local coordinates
    rtype evaluateFunction(const ctype& local) const
    {
      rtype result = 1;
      for (size_t i=0; i < nodes.size(); ++i) {
        if (i != node)
          result *= (local-nodes[i])/(nodes[node]-nodes[i]);
      }

      return result;
    }

    //! \brief Evaluate the derivative of the cardinal function
    //! \param[in] local The local coordinates
    rtype evaluateGradient(const ctype& local) const
    {
      rtype result = 0;
      for (size_t i=0; i < nodes.size(); ++i) {
        rtype f = 1;
        for (int j=0; j < nodes.size(); ++j) {
          if (i != j && j != node)
            f *= (local-nodes[j])/(nodes[node]-nodes[j]);
        }
        result += f/(nodes[node]-nodes[i]);
      }

      return result;
    }

  private:
    //! \brief The nodes
    std::vector<rtype> nodes;

    size_t node;
};

//! \brief Represents a tensor-product of 1D functions
  template<class rtype, class ctype, class ftype, int dim>
class TensorProductFunction
{
  public:
    //! \brief The dimension of the function
    enum { dimension = dim };

    //! \brief Empty default constructor
    TensorProductFunction() {}

    //! \brief Construct a tensor-product function
    //! \param[in] funcs_ The functions
    TensorProductFunction(const Dune::FieldVector<ftype, dim>& funcs_)
      : funcs(funcs_) {}

    //! \brief Evaluate the function
    //! \param[in] local The local coordinates
    rtype evaluateFunction(const Dune::FieldVector<ctype,dim>& local) const
    {
      rtype result = 1;
      for (int i=0; i < dim; ++i)
        result *= funcs[i].evaluateFunction(local[i]);

      return result;
    }

    Dune::FieldVector<rtype, dim>
    evaluateGradient(const Dune::FieldVector<ctype,dim>& local) const
    {
      Dune::FieldVector<rtype, dim> result;
      for (int i=0; i < dim; ++i)
        result[i] = funcs[i].evaluateGradient(local[i]);
    }
  private:
    Dune::FieldVector<ftype, dim> funcs;
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

  template<int dim>
class PNShapeFunctionSet
{
public:
    typedef LagrangeCardinalFunction<double, double> CardinalFunction;

    typedef TensorProductFunction<double, double, CardinalFunction, dim>
            ShapeFunction;

    PNShapeFunctionSet(int n1, int n2, int n3=0)
    {
      int dims[3] = {n1, n2, n3};
      cfuncs.resize(dim);
      for (int i=0; i < dim; ++i) {
        std::vector<double> grid;
        grid = gaussLobattoLegendreGrid(dims[i]);
        for (int j=0;j<dims[i];++j)
          cfuncs[i].push_back(CardinalFunction(grid,j));
      }
      int l=0;
      Dune::FieldVector<CardinalFunction,dim> fs;
      if (dim == 3) {
        f.resize(n1*n2*n3);
        for (int k=0; k < n3; ++k) {
          for (int j=0; j < n2; ++j)
            for (int i=0; i< n1; ++i) {
              fs[0] = cfuncs[0][i];
              fs[1] = cfuncs[1][j];
              fs[2] = cfuncs[2][k];
              f[l++] = ShapeFunction(fs);
            }
        }
      } else {
        f.resize(n1*n2);
        for (int j=0; j < n2; ++j) {
          for (int i=0; i< n1; ++i) {
            fs[0] = cfuncs[0][i];
            fs[1] = cfuncs[1][j];
            f[l++] = ShapeFunction(fs);
          }
        }
      }
    }

    //! \brief Obtain a given shape function
    //! \param[in] i The requested shape function
    const ShapeFunction& operator[](int i) const
    {
      return f[i];
    }

    int size()
    {
      return f.size();
    }
protected:
    std::vector< std::vector<CardinalFunction> > cfuncs;
    std::vector<ShapeFunction> f;

    double legendre(double x, int n)
    {
      std::vector<double> Ln;
      Ln.resize(n+1);
      Ln[0] = 1.f;
      Ln[1] = x;
      if( n > 1 ) {
        for( int i=1;i<n;i++ )
          Ln[i+1] = (2*i+1.0)/(i+1.0)*x*Ln[i]-i/(i+1.0)*Ln[i-1];
      }

      return Ln[n];
    }

    double legendreDerivative(double x, int n)
    {
      std::vector<double> Ln;
      Ln.resize(n+1);

      Ln[0] = 1.0; Ln[1] = x;

      if( (x == 1.0) || (x == -1.0) )
        return( pow(x,n-1)*n*(n+1.0)/2.0 );
      else {
        for( int i=1;i<n;i++ )
          Ln[i+1] = (2.0*i+1.0)/(i+1.0)*x*Ln[i]-(double)i/(i+1.0)*Ln[i-1];
        return( (double)n/(1.0-x*x)*Ln[n-1]-n*x/(1-x*x)*Ln[n] );
      }
    }

    std::vector<double> gaussLegendreGrid(int n)
    {
      Dune::DynamicMatrix<double> A(n,n,0.0);

      A[0][1] = 1.f;
      for (int i=1;i<n-1;++i) {
        A[i][i-1] = (double)i/(2.0*(i+1.0)-1.0);
        A[i][i+1] = (double)(i+1.0)/(2*(i+1.0)-1.0);
      }
      A[n-1][n-2] = (n-1.0)/(2.0*n-1.0);

      Dune::DynamicVector<std::complex<double> > eigenValues(n);
      Dune::DynamicMatrixHelp::eigenValuesNonSym(A, eigenValues);

      std::vector<double> result(n);
      for (int i=0; i < n; ++i)
        result[i] = std::real(eigenValues[i]);
      std::sort(result.begin(),result.begin()+n);

      return result;
    }

    std::vector<double> gaussLobattoLegendreGrid(int n)
    {
      assert(n > 1);
      const double tolerance = 1.e-15;

      std::vector<double> result(n);
      result[0] = 0.0;
      result[n-1] = 1.0;
      if (n == 3)
        result[1] = 0.5;

      if (n < 4)
        return result;

      std::vector<double> glgrid = gaussLegendreGrid(n-1);
      for (int i=1;i<n-1;++i) {
        result[i] = (glgrid[i-1]+glgrid[i])/2.0;
        double old = 0.0;
        while (std::abs(old-result[i]) > tolerance) {
          old = result[i];
          double L = legendre(old, n-1);
          double Ld = legendreDerivative(old, n-1);
          result[i] += (1.0-old*old)*Ld/((n-1.0)*n*L);
        }
        result[i] = (result[i]+1.0)/2.0;
      }

      return result;
    }
};

}
}

#endif
