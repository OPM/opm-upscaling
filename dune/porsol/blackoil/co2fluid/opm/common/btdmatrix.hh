#ifndef OPM_BLOCK_TRIDIAGONAL_MATRIX_HH
#define OPM_BLOCK_TRIDIAGONAL_MATRIX_HH

#include <opm/common/fmatrix.hh>
#include <opm/common/bcrsmatrix.hh>

/** \file
    \author Oliver Sander
    \brief Implementation of the BTDMatrix class
*/

namespace Opm {
  /** 
   * @addtogroup ISTL_SPMV 
   * @{
   */
    /** \brief A block-tridiagonal matrix 

    \todo It would be safer and more efficient to have a real implementation of
    a block-tridiagonal matrix and not just subclassing from BCRSMatrix.  But that's
    quite a lot of work for that little advantage.*/
    template <class B, class A=std::allocator<B> >
class BTDMatrix : public BCRSMatrix<B,A>
{
public:

    //===== type definitions and constants
    
    //! export the type representing the field
    typedef typename B::field_type field_type;
    
    //! export the type representing the components
    typedef B block_type;
    
    //! export the allocator type
    typedef A allocator_type;
    
    //! implement row_type with compressed vector
    //typedef BCRSMatrix<B,A>::row_type row_type;

    //! The type for the index access and the size
    typedef typename A::size_type size_type;

    //! increment block level counter
    enum {blocklevel = B::blocklevel+1};

    /** \brief Default constructor */
    BTDMatrix() : BCRSMatrix<B,A>() {}

    explicit BTDMatrix(int size) 
        : BCRSMatrix<B,A>(size, size, BCRSMatrix<B,A>::random) 
    {
        // special handling for 1x1 matrices
        if (size==1) {

            this->BCRSMatrix<B,A>::setrowsize(0, 1);
            this->BCRSMatrix<B,A>::endrowsizes();
            
            this->BCRSMatrix<B,A>::addindex(0, 0);
            this->BCRSMatrix<B,A>::endindices();
            
            return;
        }

        // Set number of entries for each row
        this->BCRSMatrix<B,A>::setrowsize(0, 2);

        for (int i=1; i<size-1; i++)
            this->BCRSMatrix<B,A>::setrowsize(i, 3);

        this->BCRSMatrix<B,A>::setrowsize(size-1, 2);

        this->BCRSMatrix<B,A>::endrowsizes();

        // The actual entries for each row
        this->BCRSMatrix<B,A>::addindex(0, 0);
        this->BCRSMatrix<B,A>::addindex(0, 1);

        for (int i=1; i<size-1; i++) {
            this->BCRSMatrix<B,A>::addindex(i, i-1);
            this->BCRSMatrix<B,A>::addindex(i, i  );
            this->BCRSMatrix<B,A>::addindex(i, i+1);
        }

        this->BCRSMatrix<B,A>::addindex(size-1, size-2);
        this->BCRSMatrix<B,A>::addindex(size-1, size-1);

        this->BCRSMatrix<B,A>::endindices();

    }

    //! assignment
    BTDMatrix& operator= (const BTDMatrix& other) {
        this->BCRSMatrix<B,A>::operator=(other);
        return *this;
    }

    //! assignment from scalar
    BTDMatrix& operator= (const field_type& k) {
        this->BCRSMatrix<B,A>::operator=(k);
        return *this;
    }

    /** \brief Use the Thomas algorithm to solve the system Ax=b in O(n) time
     *
     * \exception ISTLError if the matrix is singular
     *
     */
    template <class V>
    void solve (V& x, const V& rhs) const {

        // special handling for 1x1 matrices.  The generic algorithm doesn't work for them
        if (this->N()==1) {
            (*this)[0][0].solve(x[0],rhs[0]);
            return;
        }

        // Make copies of the rhs and the right matrix band
        V d = rhs;
        std::vector<block_type> c(this->N()-1);
        for (size_t i=0; i<this->N()-1; i++)
            c[i] = (*this)[i][i+1];

	/* Modify the coefficients. */
        block_type a_00_inv;
        FMatrixHelp::invertMatrix((*this)[0][0], a_00_inv);

        //c[0] /= (*this)[0][0];	/* Division by zero risk. */
        block_type c_0_tmp = c[0];
        FMatrixHelp::multMatrix(a_00_inv, c_0_tmp, c[0]);
        
        // d = a^{-1} d        /* Division by zero would imply a singular matrix. */
        typename V::block_type d_0_tmp = d[0];
        (*this)[0][0].solve(d[0], d_0_tmp);
        
	for (unsigned int i = 1; i < this->N(); i++) {

            // id = ( a_ii - c_{i-1} a_{i, i-1} ) ^{-1}
            block_type tmp;
            FMatrixHelp::multMatrix(c[i-1], (*this)[i][i-1], tmp);
            block_type id = (*this)[i][i];
            id -= tmp;
            id.invert();     /* Division by zero risk. */

            if (i<c.size()) {
                // c[i] *= id
                tmp = c[i];
                FMatrixHelp::multMatrix(tmp, id, c[i]);                /* Last value calculated is redundant. */
            }
            
            // d[i] = (d[i] - d[i-1] * (*this)[i][i-1]) * id;
            (*this)[i][i-1].mmv(d[i-1], d[i]);
            typename V::block_type tmpVec = d[i];
            id.mv(tmpVec, d[i]);
            //d[i] *= id;

	}
 
	/* Now back substitute. */
	x[this->N() - 1] = d[this->N() - 1];
	for (int i = this->N() - 2; i >= 0; i--) {
            //x[i] = d[i] - c[i] * x[i + 1];
            x[i] = d[i];
            c[i].mmv(x[i+1], x[i]);
        }

    }

private:	

    // ////////////////////////////////////////////////////////////////////////////
    //   The following methods from the base class should now actually be called
    // ////////////////////////////////////////////////////////////////////////////

    // createbegin and createend should be in there, too, but I can't get it to compile
    //     BCRSMatrix<B,A>::CreateIterator createbegin () {}
    //     BCRSMatrix<B,A>::CreateIterator createend () {}
    void setrowsize (size_type i, size_type s) {}
    void incrementrowsize (size_type i) {}
    void endrowsizes () {}
    void addindex (size_type row, size_type col) {}
    void endindices () {}
};
  /** @}*/

}  // end namespace Opm

#endif
