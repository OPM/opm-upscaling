// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef OPM_BASEARRAY_HH
#define OPM_BASEARRAY_HH

#include "assert.h"
#include<cmath>
#include<complex>
#include<cstddef>
#include<memory>

#include "istlexception.hh"
#include <opm/common/iteratorfacades.hh>

/** \file 
 \brief Implements several basic array containers.
*/

namespace Opm {
   
  /**  \brief A simple array container for objects of type B 

  Implement.

       -  iterator access 
       -  const_iterator access
       -  random access

	   This container has *NO* memory management at all,
	   copy constuctor, assignment and destructor are all default.

	   The constructor is made protected to emphasize that objects
       are only usably in derived classes.

	   Error checking: no error checking is provided normally.
	   Setting the compile time switch OPM_ISTL_WITH_CHECKING
	   enables error checking.
  */
  template<class B, class A=std::allocator<B> >
  class base_array_unmanaged
  {
  public:

	//===== type definitions and constants

	//! export the type representing the components
	typedef B member_type;

	//! export the allocator type
	typedef A allocator_type;

    //! the type for the index access
    typedef typename A::size_type size_type;
    

	//===== access to components

	//! random access to blocks
	B& operator[] (size_type i)
	{
#ifdef OPM_ISTL_WITH_CHECKING
	  if (i>=n) OPM_THROW(ISTLError,"index out of range");
#endif
	  return p[i];
	}

	//! same for read only access
	const B& operator[] (size_type i) const
	{
#ifdef OPM_ISTL_WITH_CHECKING
	  if (i>=n) OPM_THROW(ISTLError,"index out of range");
#endif
	  return p[i];
	}

      /** \brief Iterator implementation class  */
    template<class T>
    class RealIterator 
      :  public RandomAccessIteratorFacade<RealIterator<T>, T>
    {
    public:
      //! \brief The unqualified value type
      typedef typename remove_const<T>::type ValueType;
      
      friend class RandomAccessIteratorFacade<RealIterator<const ValueType>, const ValueType>;
      friend class RandomAccessIteratorFacade<RealIterator<ValueType>, ValueType>;
      friend class RealIterator<const ValueType>;
      friend class RealIterator<ValueType>;
      
      //! constructor
      RealIterator ()
	: p(0), i(0)
      {}
      
      RealIterator (const B* _p, B* _i) : p(_p), i(_i)
      {	  }
      
      RealIterator(const RealIterator<ValueType>& it)
	: p(it.p), i(it.i)
      {}
           
      //! return index
      size_type index () const
      {
	return i-p;
      }
            
      //! equality
      bool equals (const RealIterator<ValueType>& other) const
      {
	assert(other.p==p);
	return i==other.i;
      }
      
      //! equality
      bool equals (const RealIterator<const ValueType>& other) const
      {
	assert(other.p==p);
	return i==other.i;
      }

      std::ptrdiff_t distanceTo(const RealIterator& o) const
      {
	return o.i-i;
      }

    private:
      //! prefix increment
      void increment()
      {
	++i;
      }
      
      //! prefix decrement
      void decrement()
      {
	--i;
      }

      //! dereferencing
      B& dereference () const
      {
	return *i;
      }

      void advance(std::ptrdiff_t d)
      {
	i+=d;
      }
      
      const B* p;
      B* i;
    };

    //! iterator type for sequential access
    typedef RealIterator<B> iterator;
    
    
	//! begin iterator
	iterator begin ()
	{
	  return iterator(p,p);
	}

	//! end iterator
	iterator end ()
	{
	  return iterator(p,p+n);
	}

    //! @deprecated This method was renamed to make
    //! it distinct from the STL version which returns
    //! a reverse iterator. Use the new method beforeEnd
    //! instead.
    iterator rbegin() OPM_DEPRECATED
    {
      return beforeEnd();
    }


    //! @returns an iterator that is positioned before
    //! the end iterator of the vector, i.e. at the last entry.
	iterator beforeEnd ()
	{
	  return iterator(p,p+n-1);
	}

    //! @deprecated This method was renamed to make
    //! it distinct from the STL version which returns
    //! a reverse iterator. Use the new method beforeBegin
    //! instead.
    iterator rend ()  OPM_DEPRECATED
    {
      return beforeBegin();
    }

    //! @returns an iterator that is positioned before
    //! the first entry of the vector.
	iterator beforeBegin ()
	{
	  return iterator(p,p-1);
	}

	//! random access returning iterator (end if not contained)
	iterator find (size_type i)
	{
	  if (i<n)
		return iterator(p,p+i);
	  else
		return iterator(p,p+n);
	}

    //! iterator class for sequential access
    typedef RealIterator<const B> const_iterator;

	//! begin const_iterator
	const_iterator begin () const
	{
	  return const_iterator(p,p+0);
	}

	//! end const_iterator
	const_iterator end () const
	{
	  return const_iterator(p,p+n);
	}

    //! @deprecated This method was renamed to make
    //! it distinct from the STL version which returns
    //! a reverse iterator. Use the new method beforeEnd
    //! instead.
    const_iterator rbegin() const OPM_DEPRECATED
    {
      return beforeEnd();
    }

    //! @returns an iterator that is positioned before
    //! the end iterator of the vector. i.e. at the last element.
	const_iterator beforeEnd () const
	{
	  return const_iterator(p,p+n-1);
	}

    //! @deprecated This method was renamed to make
    //! it distinct from the STL version which returns
    //! a reverse iterator. Use the new method beforeBegin
    //! instead.
    const_iterator rend () const OPM_DEPRECATED
    {
      return beforeBegin();
    }
    
    //! @returns an iterator that is positioned before
    //! the first entry of the vector.
	const_iterator beforeBegin () const
	{
	  return const_iterator(p,p-1);
	}

	//! random access returning iterator (end if not contained)
	const_iterator find (size_type i) const
	{
	  if (i<n)
		return const_iterator(p,p+i);
	  else
		return const_iterator(p,p+n);
	}


	//===== sizes

	//! number of blocks in the array (are of size 1 here)
	size_type size () const
	{
	  return n;
	}

  protected:
	//! makes empty array
    base_array_unmanaged ()
      : n(0), p(0)
    {}
    //! make an initialized array
    base_array_unmanaged (size_type n_, B* p_)
      : n(n_), p(p_)
    {}
	size_type n; // number of elements in array
	B *p;  // pointer to dynamically allocated built-in array
  };



  /**  \brief Extend base_array_unmanaged by functions to manipulate

	   This container has *NO* memory management at all,
	   copy constuctor, assignment and destructor are all default.
	   A container can be constructed as empty or from a given pointer
       and size. This class is used to implement a view into a larger
       array.

	   You can copy such an object to a base_array to make a real copy.

	   Error checking: no error checking is provided normally.
	   Setting the compile time switch OPM_ISTL_WITH_CHECKING
	   enables error checking.
  */
  template<class B, class A=std::allocator<B> >
  class base_array_window : public base_array_unmanaged<B,A>
  {
  public:

	//===== type definitions and constants

	//! export the type representing the components
	typedef B member_type;

	//! export the allocator type
	typedef A allocator_type;

	//! make iterators available as types
	typedef typename base_array_unmanaged<B,A>::iterator iterator;

	//! make iterators available as types
	typedef typename base_array_unmanaged<B,A>::const_iterator const_iterator;

    //! The type used for the index access
    typedef typename base_array_unmanaged<B,A>::size_type size_type;

    //! The type used for the difference between two iterator positions
    typedef typename A::difference_type difference_type;
    
	//===== constructors and such

	//! makes empty array
    base_array_window ()
      : base_array_unmanaged<B,A>()
    {	}

    //! make array from given pointer and size
    base_array_window (B* _p, size_type _n)
      : base_array_unmanaged<B,A>(_n ,_p)
    {}

	//===== window manipulation methods

	//! set pointer and length
	void set (size_type _n, B* _p)
	{
	  this->n = _n;
	  this->p = _p;
	}

	//! advance pointer by newsize elements and then set size to new size
	void advance (difference_type newsize)
	{
	  this->p += this->n;
	  this->n = newsize;
	}

	//! increment pointer by offset and set size
	void move (difference_type offset, size_type newsize)
	{
	  this->p += offset;
	  this->n = newsize;
	}

	//! increment pointer by offset, leave size
	void move (difference_type offset)
	{
	  this->p += offset;
	}

	//! return the pointer
	B* getptr ()
	{
	  return this->p;
	}
  };



 /**  \brief  This container extends base_array_unmanaged by memory management
       with the usual copy semantics providing the full range of 
       copy constructor, destructor and assignment operators.

	   You can make

       - empty array
       - array with n components dynamically allocated
       - resize an array with complete loss of data
       - assign/construct from a base_array_window to make a copy of the data

	   Error checking: no error checking is provided normally.
	   Setting the compile time switch OPM_ISTL_WITH_CHECKING
	   enables error checking.
  */
  template<class B, class A=std::allocator<B> >
  class base_array : public base_array_unmanaged<B,A>
  {
  public:

	//===== type definitions and constants

	//! export the type representing the components
	typedef B member_type;

	//! export the allocator type
	typedef A allocator_type;

	//! make iterators available as types
	typedef typename base_array_unmanaged<B,A>::iterator iterator;

	//! make iterators available as types
	typedef typename base_array_unmanaged<B,A>::const_iterator const_iterator;
   
    //! The type used for the index access
    typedef typename base_array_unmanaged<B,A>::size_type size_type;

    //! The type used for the difference between two iterator positions
    typedef typename A::difference_type difference_type;
    
	//===== constructors and such

	//! makes empty array
    base_array ()
      : base_array_unmanaged<B,A>()
    {}

	//! make array with _n components
	base_array (size_type _n)
	  : base_array_unmanaged<B,A>(_n, 0)
	{
            if (this->n>0) {
		this->p = allocator_.allocate(this->n);
                new (this->p) B[this->n];
            } else
		{
		  this->n = 0;
		  this->p = 0;
		}
	}

	//! copy constructor
	base_array (const base_array& a)
	{
	  // allocate memory with same size as a
	  this->n = a.n;

	  if (this->n>0) {
              this->p = allocator_.allocate(this->n);
              new (this->p) B[this->n];
          } else
		{
		  this->n = 0;
		  this->p = 0;
		}

	  // and copy elements
	  for (size_type i=0; i<this->n; i++) this->p[i]=a.p[i];
	}

	//! construct from base class object 
	base_array (const base_array_unmanaged<B,A>& _a) 
	{
	  const base_array& a = static_cast<const base_array&>(_a);

	  // allocate memory with same size as a
	  this->n = a.n;
	  if (this->n>0) {
		this->p = allocator_.allocate(this->n);
                new (this->p) B[this->n];
          } else
		{
		  this->n = 0;
		  this->p = 0;
		}

	  // and copy elements
	  for (size_type i=0; i<this->n; i++) this->p[i]=a.p[i];
	}


	//! free dynamic memory
	~base_array () 
	{ 
            if (this->n>0) {
                int i=this->n;
                while (i)
                    this->p[--i].~B();
                allocator_.deallocate(this->p,this->n);
            }
	}

	//! reallocate array to given size, any data is lost
	void resize (size_type _n)
	{
	  if (this->n==_n) return;

	  if (this->n>0) {
              int i=this->n;
              while (i)
                  this->p[--i].~B();
              allocator_.deallocate(this->p,this->n); 
          }
	  this->n = _n;
	  if (this->n>0) {
		this->p = allocator_.allocate(this->n);
                new (this->p) B[this->n];
	  } else
		{
		  this->n = 0;
		  this->p = 0;
		}
	}

	//! assignment
	base_array& operator= (const base_array& a)
	{
	  if (&a!=this) // check if this and a are different objects
		{
		  // adjust size of array
		  if (this->n!=a.n) // check if size is different
			{
                            if (this->n>0) {
                                int i=this->n;
                                while (i)
                                    this->p[--i].~B();
                                allocator_.deallocate(this->p,this->n); // delete old memory
                            }
			  this->n = a.n;
			  if (this->n>0) {
				this->p = allocator_.allocate(this->n);
                                new (this->p) B[this->n];
			  } else
				{
				  this->n = 0;
				  this->p = 0;
				}
			}
		  // copy data
		  for (size_type i=0; i<this->n; i++) this->p[i]=a.p[i];
		}
	  return *this;
	}

	//! assign from base class object
	base_array& operator= (const base_array_unmanaged<B,A>& a)
	{
	  return this->operator=(static_cast<const base_array&>(a));
	}

        protected:

        A allocator_;
  };




  /** \brief A simple array container with non-consecutive index set.

       Elements in the array are of type B. This class provides

       -  iterator access 
       -  const_iterator access
       -  random access in log(n) steps using binary search
	   -  find returning iterator

	   This container has *NO* memory management at all,
	   copy constuctor, assignment and destructor are all default.

	   The constructor is made protected to emphasize that objects
       are only usably in derived classes.

	   Error checking: no error checking is provided normally.
	   Setting the compile time switch OPM_ISTL_WITH_CHECKING
	   enables error checking.
  */
  template<class B, class A=std::allocator<B> >
  class compressed_base_array_unmanaged
  {
  public:

	//===== type definitions and constants

	//! export the type representing the components
	typedef B member_type;

	//! export the allocator type
	typedef A allocator_type;

       //! The type used for the index access
    typedef typename A::size_type size_type;

	//===== access to components

	//! random access to blocks, assumes ascending ordering
	B& operator[] (size_type i)
	{
	  size_type l=0, r=n-1;
	  while (l<r)
		{
		  size_type q = (l+r)/2;
		  if (i <= j[q]) r=q;
		  else l = q+1;
		}
	  if (j[l]!=i){
	     OPM_THROW(ISTLError,"index "<<i<<" not in compressed array");
	  }	  
	  return p[l];
	}

	//! same for read only access, assumes ascending ordering
	const B& operator[] (size_type i) const
	{
	  size_type l=0, r=n-1;
	  while (l<r)
		{
		  size_type q = (l+r)/2;
		  if (i <= j[q]) r=q;
		  else l = q+1;
		}
	  if (j[l]!=i){
	     OPM_THROW(ISTLError,"index "<<i<<" not in compressed array");
	  }
	  return p[l];
	}

	//! iterator class for sequential access
    template<class T>
    class RealIterator
      : public BidirectionalIteratorFacade<RealIterator<T>, T>
	{
	public:
	  //! \brief The unqualified value type
	  typedef typename remove_const<T>::type ValueType;

	  friend class BidirectionalIteratorFacade<RealIterator<const ValueType>, const ValueType>;
	  friend class BidirectionalIteratorFacade<RealIterator<ValueType>, ValueType>;
	  friend class RealIterator<const ValueType>;
	  friend class RealIterator<ValueType>;

	  //! constructor
	  RealIterator ()
	    : p(0), j(0), i(0)
	  {}

	  //! constructor
	  RealIterator (B* _p, size_type* _j, size_type _i) 
	    : p(_p), j(_j), i(_i)
	  {	  }

	  /**
	   * @brief Copy constructor from mutable iterator
	   */
	  RealIterator(const RealIterator<ValueType>& it)
	    : p(it.p), j(it.j), i(it.i)
	  {}
	  
	  
	  //! equality
	  bool equals (const RealIterator<ValueType>& it) const
	  {
		assert(p==it.p);
		return (i)==(it.i);
	  }

	  //! equality
	  bool equals (const RealIterator<const ValueType>& it) const
	  {
		assert(p==it.p);
		return (i)==(it.i);
	  }


	  //! return index corresponding to pointer
	  size_type index () const
	  {
		return j[i];
	  }

	  //! Set index corresponding to pointer
	  void setindex (size_type k)
	  {
		return j[i] = k;
	  }
	  
	  /**
	   * @brief offset from the first entry.
	   *
	   * An iterator positioned at the beginning
	   * has to be increment this amount of times to
	   * the same position.
	   */
	  size_type offset () const
	  {
	    return i;
	  }
	  
	private:
	  //! prefix increment
	  void increment()
	  {
	    ++i;
	  }
	  
	  //! prefix decrement
	  void decrement()
	  {
	    --i;
	  }

	  //! dereferencing
	  B& dereference () const
	  {
	    return p[i];
	  }

	  B* p;
	  size_type* j;
	  size_type i;
	};

    /** @brief The iterator type. */
    typedef RealIterator<B> iterator;
    
	//! begin iterator
	iterator begin ()
	{
	  return iterator(p,j,0);
	}

	//! end iterator
	iterator end ()
	{
	  return iterator(p,j,n);
	}

    //! @deprecated This method was renamed to make
    //! it distinct from the STL version which returns
    //! a reverse iterator. Use the new method beforeEnd
    //! instead.
    iterator rbegin() OPM_DEPRECATED
    {
      return beforeEnd();
    }

    //! @returns an iterator that is positioned before
    //! the end iterator of the vector, i.e. at the last entry.
	iterator beforeEnd ()
	{
	  return iterator(p,j,n-1);
	}

    //! @deprecated This method was renamed to make
    //! it distinct from the STL version which returns
    //! a reverse iterator. Use the new method beforeBegin
    //! instead.
    iterator rend ()  OPM_DEPRECATED
    {
      return beforeBegin();
    }

    //! @returns an iterator that is positioned before
    //! the first entry of the vector.
	iterator beforeBegin ()
	{
	  return iterator(p,j,-1);
	}

	//! random access returning iterator (end if not contained)
	iterator find (size_type i)
	{
	  size_type l=0, r=n-1;
	  while (l<r)
		{
		  size_type q = (l+r)/2;
		  if (i <= j[q]) r=q;
		  else l = q+1;
		}

	  if (n && i==j[l])
		return iterator(p,j,l);
	  else
		return iterator(p,j,n);
	}

	//! const_iterator class for sequential access
    typedef RealIterator<const B> const_iterator;
    
	//! begin const_iterator
	const_iterator begin () const
	{
	  return const_iterator(p,j,0);
	}

	//! end const_iterator
	const_iterator end () const
	{
	  return const_iterator(p,j,n);
	}

    //! @deprecated This method was renamed to make
    //! it distinct from the STL version which returns
    //! a reverse iterator. Use the new method beforeEnd
    //! instead.
    const_iterator rbegin() const OPM_DEPRECATED
    {
      return beforeEnd();
    }

    //! @returns an iterator that is positioned before
    //! the end iterator of the vector. i.e. at the last element.
	const_iterator beforeEnd () const
	{
	  return const_iterator(p,j,n-1);
	}

    //! @deprecated This method was renamed to make
    //! it distinct from the STL version which returns
    //! a reverse iterator. Use the new method beforeBegin
    //! instead.
    const_iterator rend () const OPM_DEPRECATED
    {
      return beforeBegin();
    }
    
    //! @returns an iterator that is positioned before
    //! the first entry of the vector.
	const_iterator beforeBegin () const
	{
	  return const_iterator(p,j,-1);
	}

	//! random access returning iterator (end if not contained)
	const_iterator find (size_type i) const
	{
	  size_type l=0, r=n-1;
	  while (l<r)
		{
		  size_type q = (l+r)/2;
		  if (i <= j[q]) r=q;
		  else l = q+1;
		}

	  if (n && i==j[l])
		return const_iterator(p,j,l);
	  else
		return const_iterator(p,j,n);
	}

	//===== sizes

	//! number of blocks in the array (are of size 1 here)
	size_type size () const
	{
	  return n;
	}

  protected:
	//! makes empty array
	compressed_base_array_unmanaged ()
	  : n(0), p(0), j(0)
	{}

	size_type n;  // number of elements in array
	B *p;   // pointer to dynamically allocated built-in array
	size_type* j; // the index set
  };

} // end namespace

#endif
