//==============================================================================
//!
//! \file seqlu.hpp
//!
//! \date Dec 22 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief SuperLU preconditioner interface
//!
//==============================================================================
#pragma once

#if HAVE_SUPERLU
#include <dune/istl/superlu.hh>
#endif

#if HAVE_SUPERLU
  template<class M, class X, class Y, int l=1>
class SeqLU : public Dune::Preconditioner<X,Y> {
  public:
    //! \brief The matrix type the preconditioner is for.
    typedef typename Dune::remove_const<M>::type matrix_type;
    //! \brief The domain type of the preconditioner.
    typedef X domain_type;
    //! \brief The range type of the preconditioner.
    typedef Y range_type;
    //! \brief The field type of the preconditioner.
    typedef typename X::field_type field_type;

    // define the category
    enum {
      //! \brief The category the preconditioner is part of.
      category=Dune::SolverCategory::sequential
    };

    /*! \brief Constructor.
      
    Constructor gets all parameters to operate the prec.
    \param A The matrix to operate on.
    */
    SeqLU (const M& A) :
      slu(A, false)
    {
    }

    /*!
      \brief Prepare the preconditioner.
      
      \copydoc Preconditioner::pre(X&,Y&)
    */
    virtual void pre (X& x, Y& b) {}

    /*!
      \brief Apply the precondioner.
      
      \copydoc Preconditioner::apply(X&,const Y&)
    */
    virtual void apply (X& v, const Y& d)
    {
      Dune::InverseOperatorResult res;
      Y t(d);
      slu.apply(v, t, res);
    }

    /*!
      \brief Clean up.
      
      \copydoc Preconditioner::post(X&)
    */
    virtual void post (X& x) {}

  private:
    Dune::SuperLU<M> slu;
  };
#endif

