/*
  Copyright 2011 - Magne Haveraaen, Helmer André Friis and Hans Munthe-Kaas.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/





#ifndef ArrayH
#define ArrayH

//---------------------- Include Files ----------------------------------------
#include <cstdlib>
#include <iostream>

#include <opm/porsol/twophase2/OPMKvasiSophusCode/GlobConst.h>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/GlobType.h>



//---------------------  Constants --------------------------------------------


//---------------------- Types ------------------------------------------------






#ifdef DEBUG_Array_h
#define DEBUG_NOW
#endif


template <class T>
class Array {
  template<class R> friend class Array; 
  public:
/*---------------------------------------------------*/
/*              Data invariant                       */
/*---------------------------------------------------*/
    Boolean u_DI_0u() const;
/*---------------------------------------------------*/
/*              Constructors                         */
/*---------------------------------------------------*/
    Array(); // Array0
    Array(const Array<T>&); // Array1
    Array(const int&); // Array2
    Array(const int&, const T&); // Array3
/*---------------------------------------------------*/
/*              Observers                            */
/*---------------------------------------------------*/
    int u_getSize_0u() const;
    Boolean u_indexOK_0u(const int& p) const;
    void u_getElement_0g(const int&, T& e) const;
    //    const T& operator[](const int) const;
    T& operator[](const int) const;
    int u_findLastPosition_0u(const T& elem) const;
/*---------------------------------------------------*/
/*              Generators                           */
/*---------------------------------------------------*/
    void u_changeElement_1u(const int&, const T&);
    void u_insert_1u(const int&, const T&);
    void u_prepend_1u(const T&);
    void u_append_1u(const T&);
    void u_deleteElement_1u(const int& p);
/*--------------------------------------------------- */
/*              Composite member functions           */
/*---------------------------------------------------*/
/*              map functions                        */
    template<class S> void u_fold_0g(S(*f)(const S&, const T&), const S& z, S& R) const;	// fold1
    template<class S> void u_fold_0g(S(*f)(const S&, const T&), S& R) const;	// fold2
/*              list functions                       */
    void u_concat_1u(const Array<T>&);
    void u_reverse_1u ();
    void u_shift_1u(const int&);
    void u_substring_1u(const int& , const int&);
/*              multiset functions                   */
    Boolean u_appearin_0u(const T& ) const;
    int u_count_0u (const T&  ) const;
    void u_contract_1u ();
    void u_arUnion_1u (const Array<T>&);
    void u_intersect_1u( const Array<T>& B);
    void u_removeOne_1u( const T& t);
    void u_removeAll_1u( const T& t);
    void u_subtractOne_1u( const Array<T>& B);
    void u_subtractAll_1u( const Array<T>& B);
    Boolean u_equivalentMultisets_0u( const Array<T>& B) const;

/*---------------------------------------------------*/
/*              Elementwise operators                */
/*---------------------------------------------------*/
    void operator+=(const Array<T>& b); // default assignment operator
    void operator-=(const Array<T>& b); // default assignment operator
    void operator*=(const Array<T>& b); // default assignment operator

    void operator+=(const real& b); // default assignment operator
    void operator-=(const real& b); // default assignment operator
    void operator*=(const real& b); // default assignment operator

/*---------------------------------------------------*/
/*              SophusStandard conventions           */
/*---------------------------------------------------*/
    Boolean operator==(const Array<T>&) const;
    void    operator=(const Array<T>&); // default assignment operator
    
    ~Array(); // destructor
/*              SophusStandard copy operations       */
    void u_copy_1u(); // copy0.
    void u_copy_1u(const int&); // copy2
    void u_copy_1u(const int& , const T&); // copy3

/*---------------------------------------------------*/
/*              OrderRelations                       */
/*---------------------------------------------------*/
    Boolean operator<=( const Array<T>& a ) const;

/*---------------------------------------------------*/
/*              Implementation details               */
/*---------------------------------------------------*/
  private :
/*---------------------------------------------------*/
/*              Friend with ArrayArray               */
/*---------------------------------------------------*/
    template <class R> friend int getSize(const Array<Array <R> > & A, const int& p);
    template <class R> friend Boolean indexOK(const Array<Array <R> > & A, const int& p, const int& q);
    template <class R> friend R getElement(const Array<Array <R> > & A, const int& p, const int& q);
    template <class R> friend Array<Array <R> > changeElement(const Array<Array <R> > & A, const int& p, const int& q, const R& el);
    template <class R> friend void changeElement_upd(Array<Array <R> > & A, const int& p, const int& q, const R& el);
    template <class R> friend Array<Array <R> > insert(const Array<Array <R> > & A, const int& p, const int& q, const R& el);
    template <class R> friend Array<Array <R> > prepend(const Array<Array <R> > & A, const int& p, const R& el);
    template <class R> friend Array<Array <R> > append(const Array<Array <R> > & A, const int& p, const R& el);
    template <class R> friend Array<Array <R> > deleteElement(const Array<Array <R> > & A, const int& p, const int& q);
    template <class R> friend Boolean appearin (const Array<Array<R> >& A, const int& p, const R& t );
    template <class R> friend Boolean appearin (const Array<Array<R> >& A, const R& t );
    template <class R> friend int count (const Array<Array<R> > & A, const int& p, const R& t );
    template <class R> friend int count (const Array<Array<R> > & A, const R& t );
/*---------------------------------------------------*/
/*              Local functions                      */
/*---------------------------------------------------*/
    void ul_setNewData(T* newdata, const int& newsize, const Boolean& newconst);

    void ul_setEmpty();
    void ul_makeCopy();
    void ul_expand();
    int ul_Find1Rep(const int& i, const int& k) const;

/*---------------------------------------------------*/
/*              Data Structure                       */
/*---------------------------------------------------*/
    int        arrSize ;
    T*         data    ;
    Boolean    isConst ;
    unsigned*  refCnt  ;

//DECLARE_SOPHUS(3.21)

};


//#pragma SophusCode Functional
//#pragma SophusCode __LINE__ __FILE__


/*---------------------------------------------------*/
/*              Data invariant                       */
/*---------------------------------------------------*/

template <class T>
Boolean Array<T>::u_DI_0u() const
{
  Boolean OK = true;
  if (arrSize < 0) OK = false;
  if (arrSize == 0)
  {
    if (!(data == NULL && refCnt == NULL )) OK = false;
    if (isConst) OK = false;
  };
  if (arrSize > 0)
  {
    if (data == NULL ) OK = false;
    if (refCnt == NULL ) 
    { OK = false ; 
    }
    else
    { OK = OK && (*refCnt >= 1) ;
    }
  }
  
  //TEST_DEEPINV(for(int i=0;i<arrSize;i++){OK=OK&&DI(data[isConst?0:i]);});

  return OK;
}


/*---------------------------------------------------*/
/*              Constructors                         */
/*---------------------------------------------------*/

template <class T>
Array<T>::Array ( ) //Array0
: arrSize ( 0 ) ,
  data    ( NULL ) ,
  isConst ( false ) ,
  refCnt  ( NULL )
{
}


/*---------------------------------------------------*/
template <class T>
Array<T>::Array ( const Array<T> & a ) // Array1

: arrSize ( a.arrSize ),
  data    ( a.data    ),
  isConst ( a.isConst ),
  refCnt  ( a.refCnt  )
{
  if ( refCnt != NULL )
    ++ ( *refCnt ) ;

}


/*---------------------------------------------------*/
template <class T>
Array<T>::Array ( const int& s ) //Array2
: arrSize ( s ),
  data    ( NULL ),
  isConst ( false ),
  refCnt  ( NULL )
{
   if ( arrSize > 0 ) {
      isConst = true ;
      data    = new T[1] ;
      //TEST_ASSERT(data != NULL, "unable to allocate memory",FATAL);
      //c //TEST_ASSERT(data[0] == T(), " default Dune::array element "<<data[0]<<" "<<T(),ERROR);
      data[0] = T();
      //c //TEST_ASSERT(data[0] == T(), " default Dune::array element assigned",ERROR);
      refCnt  = new unsigned (1) ;
      //TEST_ASSERT(refCnt != NULL, "unable to allocate memory",FATAL);
   }
}

/*---------------------------------------------------*/

template <class T>
Array<T>::Array ( const int& s, const T & el )  //Array3
: arrSize ( s ),
  data    ( NULL ),
  isConst ( false ),
  refCnt  ( NULL )
{
   if ( arrSize > 0 )
   {
      isConst = true ;
      data    = new T[1];
      //TEST_ASSERT(data != NULL, "unable to allocate memory",FATAL);
      data[0] = el;
      refCnt  = new unsigned(1);
      //TEST_ASSERT(refCnt != NULL, "unable to allocate memory",FATAL);
   }
}

/*---------------------------------------------------*/
/*              Observers                            */
/*---------------------------------------------------*/

template <class T>
int Array<T>::u_getSize_0u() const
{
   return arrSize;
}

template <class T>
Boolean Array<T>::u_indexOK_0u(const int& p) const
{

   Boolean ret = (p >= 0) && (arrSize > p);

   return ret;
}








//-----------------------------------------------------------------------------
//
//  getElement
//
//  DESCRIPTION:
//    Get one single element from Dune::array. Indexing from 0 - n-1 for Dune::array
//    size of n elements.
//
//  RETURNS:
//
//  NOTE:
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------
template <class T>
void Array<T>::u_getElement_0g(const int& i, T& e) const
{
   if(isConst){
      e = data[0];
   }
   else{
      e = data [i] ;
   }
}



//-----------------------------------------------------------------------------
//
//  operator[] - Subscripting operator
//
//  DESCRIPTION:
//    Get one single element from Dune::array. Indexing from 0 - n-1 for Dune::array
//    size of n elements.
//
//  RETURNS:
//
//  NOTE:
//    Due to special implementation dependent matters, this indexing operator 
//    may only be used in rvalues, not as lvalue.
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------
template <class T>
T& Array<T>::operator[](const int i) const
{
   return data [ isConst ? 0 : i ];
}



//-----------------------------------------------------------------------------
//
//  findLastPosition
//
//  DESCRIPTION:
//    The routine requires that elem appears in the Dune::array.  
//     In that case it returns the last postion in the Dune::array where elem occurs
//  RETURNS: the last postion in the Dune::array where elem occurs.
//
//  NOTE:
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------
template <class T>
int Array<T>::u_findLastPosition_0u(const T& elem) const
{

   int lastPos = -1; //to detect possible ERRORS in the input !!!
   int totNumb = u_count_0u(elem);
   int numb = 0;
   int i = 0;

   while (numb < totNumb)
   {
     if (data[isConst ? 0 : i] == elem)
     {
       numb++;
       lastPos = i;
     }
     i++;
   }


   return lastPos;
}



/*---------------------------------------------------*/
/*              Generators                           */
/*---------------------------------------------------*/



//-----------------------------------------------------------------------------
//
//  changeElement
//
//  DESCRIPTION:
//    Assign element "el" to position "i".
//    If this object is currently utilizing the compact form, expand it
//    and sequre that it does not share the data block and reference counter
//    with any other objects. Then assign the data element in the actual
//    Dune::array position.
//
//  RETURNS:
//
//  NOTE:
//    It is the users responsibility to check that position is within 
//    actual Dune::array size.
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------
template <class T>
void Array<T>::u_changeElement_1u( const int& i, const T& el)
  // i  - the index of the element to change, 0 - n-1
  // el - the new element value
{
   ul_expand();
   data [i] = el;
}







//-----------------------------------------------------------------------------
//
//  insert
//
//  DESCRIPTION:
//    Insert element "el" into position "p".
//    If object is compact and insert element equals existing element, 
//    simply increment size. If insert element differs, expand data block to
//    new size and copy data including the insert element.
//    If object is non compact, expand and copy data including insert.
//
//  RETURNS:
//
//  NOTE:
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------
template <class T>
void Array<T>::u_insert_1u( const int& p, const T& el)
{
   if(isConst)
   {
     if(data[0] == el)
     {
       arrSize = arrSize + 1;
     }
     else 
     {
       T*   newdata = new T[arrSize + 1];
       //TEST_ASSERT(newdata != NULL, "unable to allocate memory",FATAL);
       { for ( int i = 0; i < p; i++ )
         { newdata [i] = data [0] ;
         }
       }
       newdata [p] = el;
       { for ( int i = p; i < arrSize; i++ )
         { newdata [i+1] = data [0] ;
         }
       }
       ul_setNewData(newdata,arrSize+1,false);
     }
   }
   else
   {
     T*	newdata = new T[arrSize + 1];
     //TEST_ASSERT(newdata != NULL, "unable to allocate memory",FATAL);
     { for ( int i = 0; i < p; i++ )
       { newdata [i] = data [i] ;
         //TEST_ASSERT(DI(newdata[i]), "invariant u_insert_1u " << i,FATAL);
       }
     }
     newdata [p] = el;
     { for ( int i = p; i < arrSize; i++ )
       { newdata [i+1] = data [i] ;
       }
     }
     ul_setNewData(newdata,arrSize+1,false);
   }

}








//-----------------------------------------------------------------------------
//
//  prepend
//
//  DESCRIPTION:
//    Insert element at the beginning of the data block. This is just a 
//    particular case of insert.
//
//  RETURNS:
//
//  NOTE:
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------
template <class T>
void Array<T>::u_prepend_1u( const T& el)
{
   u_insert_1u(0,el);
}







//-----------------------------------------------------------------------------
//
//  append
//
//  DESCRIPTION:
//    Append the element "el" after the last element of this Dune::array.
//
//  RETURNS:
//
//  NOTE:
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------
template <class T>
void Array<T>::u_append_1u( const T& el)
{
   u_insert_1u(arrSize,el);
}






/*---------------------------------------------------*/
template <class T>
void Array<T>::u_deleteElement_1u(const int& p)
{
   if(arrSize == 1)
   {
     ul_setEmpty();
   } 
   else if(isConst)
   {
     arrSize = arrSize - 1;
   }
   else
   {
      T*	newdata = new T[arrSize - 1];
      //TEST_ASSERT(newdata != NULL, "unable to allocate memory",FATAL);
      { for ( int i = 0; i < p; i++ )
        { newdata [i] = data [i] ;
        }
      }
      { for ( int i = p+1; i < arrSize; i++ )
        {  newdata [i-1] = data [i] ;
        }
      }
      ul_setNewData(newdata,arrSize-1,false);
   }

}

/*---------------------------------------------------*/
/*              Composite member functions           */
/*---------------------------------------------------*/

/*              map functions                        */
/*---------------------------------------------------*/
/*---------------------------------------------------*/


template <class T> template <class S>
void Array<T>::u_fold_0g(S(*f)(const S&, const T&), const S& z, S& R) const
{
   // //TEST_TRACE("u_fold_0g"<<z);

   R = z ;

   if ( arrSize > 0){
      if ( isConst ) {
         for ( register int i = 0; i < arrSize; i++ )
            R = f ( R, data[0] ) ;
      }
      else
         for ( register int i = 0; i < arrSize; i++ )
            R = f ( R, data[i] ) ;
   }

}

/*---------------------------------------------------*/


template <class T> template <class S>
void Array<T>::u_fold_0g( S(*f)(const S&, const T&), S& R) const
{
   // //TEST_TRACE("u_fold_0g");

   R = data[0] ;

   if ( isConst ) {
      for ( register int i = 1; i < arrSize; i++ )
         R = f ( R, data[0] ) ;
   }
   else
      for ( register int i = 1; i < arrSize; i++ )
         R = f ( R, data[i] ) ;

}


/*---------------------------------------------------*/
/*              list functions                       */
/*---------------------------------------------------*/







//-----------------------------------------------------------------------------
//
//  concat
//
//  DESCRIPTION:
//    Concats another Dune::array at the end of this object. This object will
//    become a non constant object.
//
//  RETURNS:
//
//  NOTE:
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------
template <class T>
void Array<T>::u_concat_1u(const Array<T>& B)
{

   int  newsize = arrSize + B.arrSize;
   if ( newsize != 0 )
   {
     T*   newdata = new T[newsize];
     //TEST_ASSERT(newdata != NULL, "unable to allocate memory",FATAL);
     if ( isConst )
     { for ( int i = 0 ; i < arrSize ; i++ ) newdata[i] = data[0] ;
     }
     else
     { for ( int i = 0 ; i < arrSize ; i++ ) newdata[i] = data[i] ;
     }
  
     if ( B.isConst )
     { for ( int i = 0 ; i < B.arrSize ; i++ ) newdata[arrSize+i] = B.data[0] ;
     }
     else
     { for ( int i = 0 ; i < B.arrSize ; i++ ) newdata[arrSize+i] = B.data[i] ;
     }
  
     ul_setNewData (newdata, newsize, false);
   }

}








//-----------------------------------------------------------------------------
//
//  reverse
//
//  DESCRIPTION:
//    Reverses the element order of the object.
//
//  RETURNS:
//
//  NOTE:
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------
template <class T>
void Array<T>::u_reverse_1u ()
{

   if ( arrSize>0 && ! isConst )	// if there is a need to move data
   {
      T* newdata = new T[arrSize];
      //TEST_ASSERT(newdata != NULL, "unable to allocate memory",FATAL);

      for (register int i=0; i<arrSize; i++)
	  newdata[arrSize-1-i] = data[i];

      ul_setNewData(newdata,arrSize,isConst);

   }
}


/*---------------------------------------------------*/

template <class T>
void Array<T>::u_shift_1u ( const int& offset )
{
   if ( arrSize>1 && ! isConst )	// if there is a need to shift data
   {
      T*   newdata = new T[arrSize];
      //TEST_ASSERT(newdata != NULL, "unable to allocate memory",FATAL);
      int  p = offset<0 ? offset+arrSize : offset%arrSize ;      
      while (p<0) { p = p + arrSize; }

      for (int i=0; i<arrSize; i++)
	{ newdata[i] = data[(i+p)%arrSize]; }

      ul_setNewData(newdata,arrSize,isConst);
   }
}


/*--------------------------------------------------*/

template <class T>
void Array<T>::u_substring_1u (const int& m , const int& n)
{
   if ( m == n )
   {
     ul_setEmpty();
   }
   else if(isConst )
   {
     arrSize = n - m;
   }
   else // if there is a need to move data
   {
     T*   newdata = new T[n-m];
     //TEST_ASSERT(newdata != NULL, "unable to allocate memory",FATAL);
     for ( register int i=0; i<n-m; i++)
     { newdata[i] = data[i+m];
     }
     ul_setNewData(newdata,n-m,false);
   }
}




/*-------------------------------------------------- */
/*              multiset functions                   */
/*-------------------------------------------------- */







//-----------------------------------------------------------------------------
//
//  appearin
//
//  DESCRIPTION:
//    Checks whether element "t" appears in this object.
//
//  RETURNS:
//
//  NOTE:
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------
template <class T>
Boolean Array<T>::u_appearin_0u (const T& t ) const
{
   Boolean ret = false;
   if (isConst)
   { ret = data[0] == t;
   }
   else
   {
     for (register int i=0; i < arrSize; i++)
     {
       if ( t == data[i] )
       {
         ret = true;
         break;
       }
     }
   }
   return ret;
}








//-----------------------------------------------------------------------------
//
//  count
//
//  DESCRIPTION:
//    Counts how many object elements which equals the element "t".
//
//  RETURNS:
//
//  NOTE:
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------
template <class T>
int Array<T>::u_count_0u (const T& t ) const
{

   int Icount = 0;
   if ( isConst )
   { Icount = data[0] == t ? arrSize : 0;
   }
   else
   {
     for ( register int i=0; i < arrSize; i++)
     {
         if ( t == data[i] ){
            Icount = Icount + 1 ;
         }
     }
   }
   return Icount;
}

/*----------------------------------------------------*/

template <class T>
void Array<T>::u_arUnion_1u( const Array<T>& B)
{

   u_concat_1u(B);

}

/*----------------------------------------------------*/

template <class T>
void Array<T>::u_intersect_1u( const Array<T>& B)
{

   Boolean Found = true;
   while ( Found == true ) {
      Found = false;
      for ( int i = 1; i <= arrSize; i++ )
         {
         if ( B.u_count_0u(data[isConst ? 0 : arrSize-i]) 
              < u_count_0u(data[isConst ? 0 : arrSize-i]) )
            {
              u_deleteElement_1u(arrSize-i);
              Found = true ;
              break;
            }
         }
   }

}

/*----------------------------------------------------*/

template <class T>
void Array<T>::u_removeOne_1u( const T& t)
{

   if ( u_appearin_0u (t) )
   { if ( arrSize == 1 )
     { ul_setEmpty();
     }
     else if ( isConst )
     { arrSize = arrSize - 1;
     }
     else
     { int size = arrSize;
       for ( register int i = 1; i <= size; i++ )
        {
          if ( data[size-i] == t)
          {
            u_deleteElement_1u(size-i);
            break;
          }
        }
     }
   }

}

/*----------------------------------------------------*/

template <class T>
void Array<T>::u_removeAll_1u( const T& t)
{

   if ( (arrSize == 1) || isConst  )
   { if ( data[0] == t )
     { 
       ul_setEmpty();
     }
   }
   else
   { int size = arrSize;
     for ( register int i = 1; i <= size; i++ )
      {
        if ( data[size-i] == t)
        { 
          u_deleteElement_1u(size-i);
        }
      }
   }

}

/*----------------------------------------------------*/

template <class T>
void Array<T>::u_subtractOne_1u( const Array<T>& B)
{

   for ( register int i=0; i < B.arrSize; i++)
     {
     u_removeOne_1u(B.data[B.isConst?0:i]);
     }

}

/*----------------------------------------------------*/

template <class T>
void Array<T>::u_subtractAll_1u( const Array<T>& B)
{

   for ( register int i=0; i < B.arrSize; i++)
     {
     u_removeAll_1u(B.data[B.isConst?0:i]);
     }

}

/*----------------------------------------------------*/

template <class T>
Boolean Array<T>::u_equivalentMultisets_0u( const Array<T>& B) const
{

   Boolean ret = arrSize == B.arrSize;

   if(arrSize==0 || B.arrSize==0)
   {
   }
   else if(isConst)	// arrSize == B.arrSize
   {
     ret = arrSize == B.u_count_0u(data[0]);
   }
   else if(B.isConst)	// arrSize == B.arrSize
   {
     ret = B.arrSize == u_count_0u(B.data[0]);
   }
   else
   {
     for ( register int i=0; i < arrSize; i++)
     {
       if( u_count_0u(data[i]) != B.u_count_0u(data[i]) )
       {
         ret = false;
         break;
       }
     }
   }
   return ret;
}




/*---------------------------------------------------*/
/*              Elementwise operators                */
/*---------------------------------------------------*/

//-----------------------------------------------------------------------------
//
//  operator+= - Elementwise addition operator
//
//  DESCRIPTION:
//    Adds each element of b to the corresponding Array class elemement.
//
//  RETURNS:
//
//  NOTE:
//    We MUST have that arrSize == b.arrSize to use thie function properly !!!
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------
template <class T>
void Array<T>::operator+=( const Array<T> & b )
{

   ul_expand();

   for (int i=0; i < arrSize; i++)
   {
     data[i] += b.data[i];
   }

}


//-----------------------------------------------------------------------------
//
//  operator-= - Elementwise subtraction operator
//
//  DESCRIPTION:
//    Subtracts each element of b from the corresponding Array class elemement.
//
//  RETURNS:
//
//  NOTE:
//    We MUST have that arrSize == b.arrSize to use thie function properly !!!
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------
template <class T>
void Array<T>::operator-=( const Array<T> & b )
{

   ul_expand();

   for (int i=0; i < arrSize; i++)
   {
     data[i] -= b.data[i];
   }

}



//-----------------------------------------------------------------------------
//
//  operator+= - Elementwise multiplication operator
//
//  DESCRIPTION:
//    Multiplicates each element of b to the corresponding Array class elemement.
//
//  RETURNS:
//
//  NOTE:
//    We MUST have that arrSize == b.arrSize to use thie function properly !!!
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------
template <class T>
void Array<T>::operator*=( const Array<T> & b )
{

   ul_expand();

   for (int i=0; i < arrSize; i++)
   {
     data[i] *= b.data[i];
   }

}


//-----------------------------------------------------------------------------
//
//  operator+= - Elementwise addition operator
//
//  DESCRIPTION:
//    Adds b to the corresponding Array class elemement.
//
//  RETURNS:
//
//  NOTE:
//    We MUST have that arrSize == b.arrSize to use thie function properly !!!
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------
template <class T>
void Array<T>::operator+=( const real& b )
{

   ul_expand();

   for (int i=0; i < arrSize; i++)
   {
     data[i] += b;
   }

}


//-----------------------------------------------------------------------------
//
//  operator-= - Elementwise subtraction operator
//
//  DESCRIPTION:
//    Subtracts b from the corresponding Array class elemement.
//
//  RETURNS:
//
//  NOTE:
//    We MUST have that arrSize == b.arrSize to use thie function properly !!!
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------
template <class T>
void Array<T>::operator-=( const real& b )
{

   ul_expand();

   for (int i=0; i < arrSize; i++)
   {
     data[i] -= b;
   }

}



//-----------------------------------------------------------------------------
//
//  operator+= - Elementwise multiplication operator
//
//  DESCRIPTION:
//    Multiplicates b to the corresponding Array class elemement.
//
//  RETURNS:
//
//  NOTE:
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------
template <class T>
void Array<T>::operator*=( const real& b )
{

   ul_expand();

   for (int i=0; i < arrSize; i++)
   {
     data[i] *= b;
   }

}




/*---------------------------------------------------*/
/*              SophusStandard conventions           */
/*---------------------------------------------------*/

template <class T>
Boolean Array<T>::operator==(const Array<T> & b ) const
{
  Boolean res  =  ( arrSize == b.arrSize ) ;
   if ( res && (arrSize > 0) )
   {  
     if ( data == b.data )
     { // pointers to the same data
     }
     else 
     if ( isConst && b.isConst )
     {
       res = ( data[0] == b.data[0] ) ;
     }
     else if (isConst)
     {
       int i = 0;
       while ( res && i < b.arrSize ) 
       {
         if ( !(data[0] == b.data[i]) )
         { res = false ;
         }
         i++ ;
       }
     }
     else if (b.isConst)
     { 
       int i = 0;
       while ( res && i < arrSize ) 
       {
         if ( !(data[i] == b.data[0]) )
         { res = false ;
         }
         i++ ;
       }
     }
     else 
     { 
       int i = 0;
       while ( res && i < arrSize ) 
       {
         if ( !(data[i] == b.data[i]) )
         { res = false ;
         }
         i++ ;
       }
     }
   }
   return res;
}


/*---------------------------------------------------*/

template <class T>
void Array<T>::operator=( const Array<T> & a )
{

      ul_setEmpty();

      arrSize = a.arrSize ;
      data    = a.data    ;
      isConst = a.isConst ;
      refCnt  = a.refCnt  ;
      if ( refCnt != NULL ) (*refCnt)++ ;

}




/*---------------------------------------------------*/

template <class T>
Array<T>::~Array ( )
{

   ul_setEmpty();

}



/*              SophusStandard copy operations       */





//-----------------------------------------------------------------------------
//
//  copy(0) - Erase object data
//
//  DESCRIPTION:
//    Erases the object data.
//    The object will appear as if it was constructed by constructor(0).
//
//  RETURNS:
//
//  NOTE:
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------
template <class T>
void Array<T>::u_copy_1u() // copy0
{

  ul_setEmpty();

}







//-----------------------------------------------------------------------------
//
//  copy(2) - Preset size
//
//  DESCRIPTION:
//    Presets size of object to new size. This imply cleaning the object, then
//    initializing it as a "single constant" object with actual size.
//    The data block is created but contains pure garbage.
//    The object will appear as if it was constructed by constructor(2).
//
//  RETURNS:
//
//  NOTE:
//
//  EXAMPLE:
//    Array<int> a1;       // Instanciate empty Dune::array
//    a1 = copy(a1, 10);   // Preset the Dune::array size to 10 int
//
//
//-----------------------------------------------------------------------------
template <class T>
void Array<T>::u_copy_1u(const int& newSize) // copy2
{

  ul_setEmpty();

  if (newSize > 0)
  {
      arrSize = newSize;
      isConst = true;
      data = new T[1];
      //TEST_ASSERT(data != NULL, "unable to allocate memory",FATAL);
      data[0] = T();
      refCnt = new unsigned(1);
      //TEST_ASSERT(refCnt != NULL, "unable to allocate memory",FATAL);
  }

}








//-----------------------------------------------------------------------------
//
//  copy(3) - Preset size/data
//
//  DESCRIPTION:
//    Presets size of object to new size. This imply cleaning the object, then
//    initializing it as a "single constant" object with actual size.
//    The data block is created and the constant element assigned.
//    The object will appear as if it was constructed by constructor(3).
//
//  RETURNS:
//
//  NOTE:
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------
template <class T>
void Array<T>::u_copy_1u(const int& newSize, const T& t) // copy3
{

  ul_setEmpty();

  if (newSize > 0)
  {
      arrSize = newSize;
      isConst = true;
      data = new T[1];
      //TEST_ASSERT(data != NULL, "unable to allocate memory",FATAL);
      data[0] = t;
      refCnt = new unsigned (1);
      //TEST_ASSERT(refCnt != NULL, "unable to allocate memory",FATAL);
  }


}




/*---------------------------------------------------*/
/*              Local functions                      */
/*---------------------------------------------------*/




//-----------------------------------------------------------------------------
//
//  ul_setNewData
//
//  DESCRIPTION:
//    Assigns new object header data to object. Delete any existing data block and 
//    sequre uncoupling if any common data block exist.
//
//  RETURNS:
//
//  NOTE:
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------
template <class T>
void Array<T>::ul_setNewData(T* newdata, const int& newsize, const Boolean& newconst)

{

   if ( refCnt == NULL)
   {
      arrSize = newsize;
      data    = newdata;
      isConst = newconst;
      refCnt  = new unsigned (1) ;
      //TEST_ASSERT(refCnt != NULL, "unable to allocate memory",FATAL);
   }
      else
   {
      if ( *refCnt == 1)
      {
         delete [] data;		// special test for isConst not needed
         delete refCnt;

         arrSize = newsize;
         data    = newdata;
         isConst = newconst;
         refCnt  = new unsigned (1);
         //TEST_ASSERT(refCnt != NULL, "unable to allocate memory",FATAL);
      }
      if ( *refCnt >  1)
      {
         --(*refCnt) ;
         arrSize = newsize;
         data    = newdata;
         isConst = newconst;
         refCnt  = new unsigned (1);
         //TEST_ASSERT(refCnt != NULL, "unable to allocate memory",FATAL);
      }
   }

}





//-----------------------------------------------------------------------------
//
//  ul_setEmpty
//
//  DESCRIPTION:
//    Erase object data block and reference counter.
//    Delete any data if object is "single" object, and uncouple any 
//    connection to common data block if object is a "common" object. 
//    The object will be as created by constructor 0 (zero size, 
//    no data block, no ref. counter).
//
//  RETURNS:
//
//  NOTE:
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------
template <class T>
void Array<T>::ul_setEmpty()
{

  if ( arrSize > 0 )
  {
     if ((*refCnt) > 1 )
     {
       (*refCnt)--;
     }
     else // if ((*refCnt) == 1)
     {
       delete [] data;		// special test for isConst not needed
       delete refCnt;
     }
     data = NULL;       // Set the data-pointer (from "this object" to data block) to NULL
     refCnt = NULL;     // Set the refCnt-pointer (from "this object" to refCnt) to NULL
     arrSize = 0;
     isConst = false;
  }


}






//-----------------------------------------------------------------------------
//
//  ul_makeCopy
//
//  DESCRIPTION:
//    If this object refer to a common data block with other objects,
//    uncouple it by making a copy of the data block and reference counter.
//
//  RETURNS:
//
//  NOTE:
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------
template <class T>
void  Array<T>::ul_makeCopy ()
{

   if (refCnt != NULL) 
   { if ( (*refCnt) != 1 )
     {
       (*refCnt)-- ;
       refCnt = new unsigned (1) ;
       //TEST_ASSERT(refCnt != NULL, "unable to allocate memory",FATAL);

       T   * arg = data    ;
       T   * res = new T [ isConst ? 1 : arrSize ]  ;
       //TEST_ASSERT(res != NULL, "unable to allocate memory",FATAL);

       data   = res ;

       if ( isConst )
       { data[0] = arg[0] ;
       }
       else 
       {
         for ( register int i = 0 ; i < arrSize ; i++ )
            data [ i ] = arg [ i ] ;
       }
     }
   }

}





//-----------------------------------------------------------------------------
//
//  ul_expand
//
//  DESCRIPTION:
//    If the object is a constant object, it is changed to non constant and
//    expanded to the actual object size. If the data block is common with
//    other objects, this object uncouples from the common datablock and 
//    reference counter.
//
//    If the object is non constant, copy the data block if several
//    objects refers to it.
//
//
//  RETURNS:
//
//  NOTE:
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------
template <class T>
void  Array<T>::ul_expand ()
{

   if ( isConst && arrSize > 1) 
   {
      isConst = false ;

      T* res = new T [ arrSize ] ;
      //TEST_ASSERT(res != NULL, "unable to allocate memory",FATAL);
      for (int i = 0 ; i < arrSize ; i++ )
         res[i] = data[0] ;

      if ((*refCnt) != 1) // data is shared
      { (*refCnt)--;
	refCnt = new unsigned(1);
        //TEST_ASSERT(refCnt != NULL, "unable to allocate memory",FATAL);
      }
      else		 // private copy of old data
      { 
        delete [] data; 
      }

     data   = res ;
   }
   else
   {
     ul_makeCopy();
   }

}

/*---------------------------------------------------*/

template <class T>
int Array<T>::ul_Find1Rep(const int& i, const int& k) const
{


   int ret = -1;
   for ( register int j = k ; j < arrSize ; ++ j )
   {
      if(data[j] == data[k] )
      {
        ret = j;
        break;
      }
   }
   return ret;
}


/*---------------------------------------------------*/
/*              OrderRelations                       */
/*---------------------------------------------------*/


template <class T>
Boolean Array<T>::operator<=( const Array<T>& a ) const
{

   Boolean ret = true;
   Boolean allequal = true;
   int size = arrSize <= a.arrSize ? arrSize : a.arrSize;
   if (isConst && a.isConst)
   { allequal = data[0] == a.data[0];
     ret = allequal ? true : data[0] <= a.data[0];
   }
   else if (isConst)
   {
     for ( int i = 0 ; i < size ; i++ )
     {
       if( allequal && data[0] == a.data[i] )
       {
       }
       else
       { allequal = false;
         if ( data[0] > a.data[i] )
         {
           ret = false;
         }
         break;
       }
     }
   }
   else if (a.isConst)
   {
     for ( register int i = 0 ; i < size ; i++ )
     {
       if( allequal && data[i] == a.data[0] )
       {
       }
       else
       { allequal = false;
         if ( data[i] > a.data[0] )
         {
           ret = false;
         }
         break;
       }
     }
   }
   else 
   {
     for ( register int i = 0 ; i < size ; i++ )
     {
       if( allequal && data[i] == a.data[i] )
       {
       }
       else
       { allequal = false;
         if ( data[i] > a.data[i] )
         {
           ret = false;
         }
         break;
       }
     }
   }
   ret = allequal ? arrSize <= a.arrSize : ret;
   return ret;
}




#endif //ArrayH




/*
===============================================================================
    ------------------------ CLASS DESCRIPTION --------------------------
===============================================================================

  DESCRIPTION:
  ------------------
    <Public description for the class users>


  USAGE/EXAMPLES:
  ---------------
    <Examples how to use the class>
    
    
  MISCELLANEOUS:
  --------------
  


===============================================================================
    --------------------- END OF CLASS DESCRIPTION ----------------------
===============================================================================
*/
