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



#ifndef PairH
#define PairH

//---------------------- Include Files ----------------------------------------


#include <opm/porsol/twophase2/OPMKvasiSophusCode/GlobConst.h>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/GlobType.h>


//---------------------  Constants --------------------------------------------


//---------------------- Types ------------------------------------------------




#ifdef DEBUG_Pair_h
#define DEBUG_NOW
#endif


//-----------------------------------------------------------------------------
template<class S, class T>
class Pair
{
public:
  // Constructors
  Pair(); // default constructor
  Pair(const Pair & PP); // copy constructor
  Pair(const S & sc, const T & tc); // constructor

  // Observers
 
  Boolean u_DI_0u() const; // Data invariant
  S u_first_0u() const; // returns first element of pair
  T u_second_0u() const; // returns second element of pair

  // Generators
  void u_updateFirst_1u(const S& sc); // updates first element of pair with input
  void u_updateSecond_1u(const T& tc); // updates second element of pair with input

  // Map operations
  template<class U> Pair<U,T> u_mapFirst_0u(U (*F)(const S&)) const;
  template<class V> Pair<S,V> u_mapSecond_0u(V (*G)(const T&)) const;
  template<class U, class V> Pair<U,V> u_mapPair_0u(U (*F)(const S&), V (*G)(const T&)) const;
  template<class U, class V> Pair<U,V> u_mapPair_0u(Pair<U,V> (*FG)(const S&, const T&)) const;
  template<class U> Pair<U,U> u_map_0u(U (*FG)(const T&)) const;


  // Composite operations
  template<class U, class V> Pair<U,V> u_both_0u(U(*F)(const S&, const T&), V(*G)(const S&, const T&)) const;


/*              SophusStandard conventions           */
    Boolean operator==(const Pair& PP) const;
    void operator=(const Pair& PP);
    ~Pair(); // destructor
/*              SophusStandard copy operations       */
    void u_copy_1u(); // copy0.
    void u_copy_1u(const S & sc, const T & tc); // copy2

private:
  S m_s; // the first template argument
  T m_t; // the second template argument

//  DECLARE_SOPHUS(1.14)

};




// Constructors

// Default constructor
//-----------------------------------------------------------------------------
template<class S, class T>
Pair<S,T>::Pair():
m_s (S()),  
m_t (T())
//----------------------------------------------------------------------------
{ // creates an empty Pair
}

// Copy constructor
//-----------------------------------------------------------------------------
template<class S, class T>
Pair<S,T>::Pair(const Pair & PP):
m_s (PP.u_first_0u()),  
m_t (PP.u_second_0u())
//----------------------------------------------------------------------------
{
}

// Builder
//-----------------------------------------------------------------------------
template<class S, class T>
Pair<S,T>::Pair(const S & sc, const T & tc):
m_s (sc),
m_t (tc)
//----------------------------------------------------------------------------
{
}

// Data Invariant
//-----------------------------------------------------------------------------
template<class S, class T>
Boolean Pair<S,T>::u_DI_0u() const
//-----------------------------------------------------------------------------
{
  bool OK = DI(m_s) && DI(m_t);
  return OK;
}

// Observers

// returns first component
//-----------------------------------------------------------------------------
template<class S, class T>
S Pair<S,T>::u_first_0u() const
//----------------------------------------------------------------------------
{
  return m_s;
}

// returns second component
//-----------------------------------------------------------------------------
template<class S, class T>
T Pair<S,T>::u_second_0u() const
//----------------------------------------------------------------------------
{
 return m_t;
}

// Generators

// Sets m_s equal to sc
//-----------------------------------------------------------------------------
template<class S, class T>
void Pair<S,T>::u_updateFirst_1u(const S& sc)
//----------------------------------------------------------------------------
{
  m_s = sc;
}

// Sets m_t equal to tc
//-----------------------------------------------------------------------------
template<class S, class T>
void Pair<S,T>::u_updateSecond_1u(const T& tc)
//----------------------------------------------------------------------------
{
  m_t = tc;
}

// Map Operations

// Mapping of a function on one element of the pair, returns a new Pair
//-----------------------------------------------------------------------------
template<class S, class T> template<class U>
Pair<U,T> Pair<S,T>::u_mapFirst_0u(U (*F)(const S&)) const
//----------------------------------------------------------------------------
{
  U retSP = (*F)(m_s);
  Pair<U, T> newPair = Pair(retSP, m_t);
  return newPair;
}

// Mapping of a function on one element of the pair, returns a new Pair
//-----------------------------------------------------------------------------
template<class S, class T> template<class V>
Pair<S,V> Pair<S,T>::u_mapSecond_0u(V (*G)(const T&)) const
//----------------------------------------------------------------------------
{
  V retTP = (*G)(m_t);
  Pair<S, V> newPair = Pair(m_s, retTP);
  return newPair;

}

// Mapping of two functions on each element of the pair, returns a new Pair
//-----------------------------------------------------------------------------
template<class S, class T> template<class U, class V>
Pair<U,V> Pair<S,T>::u_mapPair_0u(U (*F)(const S&), V (*G)(const T&)) const
//----------------------------------------------------------------------------
{
  U retSP = (*F)(m_s);
  V retTP = (*G)(m_t);
  Pair<U, V> newPair = Pair(retSP, retTP);
  return newPair;
}








#ifdef COMPLETE_SOPHUS
//---------------------- Removed from -----------------------------------------


//-----------------------------------------------------------------------------
template<class S, class T> template<class U, class V>
Pair<U,V> Pair<S,T>::u_mapPair_0u(Pair<U,V> (*FG)(const S&, const T&)) const
//----------------------------------------------------------------------------
{
  Pair<u, v> newPair = Pair(FG(m_s,m_t));
  return newPair;
}

// Mapping of a function on both element of the pair if they are of the same type, returns a new Pair
//-----------------------------------------------------------------------------
template<class S, class T> template<class U>
Pair<U,U> Pair<S,T>::u_map_0u(U (*FG)(const T&)) const
//----------------------------------------------------------------------------
// This may not work: then put comments around entire function
{
  U retSP = (*FG)(m_s);
  U retTP = (*FG)(m_t);
  Pair<U, U> newPair = Pair(retSP, retTP);
  return newPair;
}


// Composite operations

template<class S, class T> template<class U, class V>
Pair<U,V> Pair<S,T>::u_both_0u(U(*F)(const S&, const T&), V(*G)(const S&, const T&)) const
//----------------------------------------------------------------------------
{
  U retSP = (*F)(m_s, m_t);
  V retTP = (*G)(m_s, m_t);
  Pair<U, V> newPair = Pair(retSP, retTP);
  return newPair;
}


//---------------------- Removed to -------------------------------------------
#endif //COMPLETE_SOPHUS




/*---------------------------------------------------*/
/*              SophusStandard conventions           */
/*---------------------------------------------------*/


template<class S, class T>
Boolean Pair<S,T>::operator==(const Pair& PP) const
{
  bool ok = false;
  if(m_s == PP.m_s && m_t == PP.m_t) ok = true;
  return ok;
}


template<class S, class T>
void Pair<S,T>::operator=(const Pair& PP) 
{
  m_s = PP.m_s;
  m_t = PP.m_t;
}

template <class S, class T>
Pair<S,T>::~Pair()
{
}


/*              SophusStandard copy operations       */


template <class S, class T>
void Pair<S,T>::u_copy_1u() // copy0
{

  m_s = copy(m_s);
  m_t = copy(m_t);

}


template <class S, class T>
void Pair<S,T>::u_copy_1u(const S & sc, const T & tc) // copy2
{

  m_s = sc;
  m_t = tc;

}





#endif //PairH



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

