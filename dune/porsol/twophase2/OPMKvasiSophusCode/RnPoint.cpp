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


//---------------------- Include Files ----------------------------------------
#include "RnPoint.h"

#include <cmath>


//---------------------  Constants --------------------------------------------
#define TOLERANCE 10e-9


//---------------------- Types ------------------------------------------------








//---------------------- Public Functions -------------------------------------





// Constructors

//-----------------------------------------------------------------------------
RnPoint::RnPoint()
: m_rns(RnShape()), m_nsd(0)
//-----------------------------------------------------------------------------
{

}

//-----------------------------------------------------------------------------
RnPoint::RnPoint(const RnPoint& tp) 
: m_rns(tp.m_rns), m_nsd(tp.m_nsd)
//-----------------------------------------------------------------------------
{

    for (int i = 0; i < m_nsd; i++) m_coor[i] = tp.m_coor[i];

}


//-----------------------------------------------------------------------------
RnPoint::RnPoint(const RnShape& ts) 
: m_rns(ts), m_nsd(ts.u_domainDimension_0u())
//-----------------------------------------------------------------------------
{

    for (int i = 0; i < m_nsd; i++) m_coor[i] = 0.0;

}

  
// Observers 
//----------------------------------------------------------------------------
int RnPoint::u_domainDimension_0u() const
//----------------------------------------------------------------------------
{
  return m_nsd;
}


//----------------------------------------------------------------------------
RnShape RnPoint::u_getShape_0u() const
//----------------------------------------------------------------------------
{
  return m_rns;
}

//----------------------------------------------------------------------------
Boolean RnPoint::u_sameShape_0u(const RnPoint& Q) const
//----------------------------------------------------------------------------
{
  Boolean OK;
  if (m_rns == Q.m_rns) {
    OK = true;
  }
  else {
    OK = false;
  }
  return OK;

}

//----------------------------------------------------------------------------
Boolean RnPoint::u_legalDirection_0u(const int dir) const
//----------------------------------------------------------------------------
{

  Boolean is_dir = false;
  if (dir >= 1 && dir <= m_nsd)
    is_dir = true;
  return is_dir;
}


//----------------------------------------------------------------------------
Boolean RnPoint::u_indexOK_0u(const int ind) const
//----------------------------------------------------------------------------
{

   Boolean ret = (ind >= 0) && (ind < m_nsd);

   return ret;
}

//-----------------------------------------------------------------------------
real RnPoint::u_getElement_0u(int dir) const
//-----------------------------------------------------------------------------
{
  return m_coor[dir];
}

//-----------------------------------------------------------------------------
real RnPoint::u_getComponent_0u(const int dir) const
//-----------------------------------------------------------------------------
{
    return m_coor[dir-1];
}

//-----------------------------------------------------------------------------
Boolean RnPoint::u_isOrigin_0u() const
//-----------------------------------------------------------------------------
{
  Boolean is_origin = true;
  int i;

  for ( i= 0; i < m_nsd; i++) {
    if ( m_coor[i] != 0.0) { 
      is_origin = false;
      break;
    }
  }
  return is_origin;
}


// Generators

//-----------------------------------------------------------------------------
void RnPoint::u_setShape_1u(const RnShape& S) 
//-----------------------------------------------------------------------------
{
  m_rns = S;
  m_nsd = S.u_domainDimension_0u();
  for (int i=0; i<m_nsd; i++) {
    m_coor[i] = 0.0;
  }
}

//-----------------------------------------------------------------------------
void RnPoint::u_setPoint_1u(const RnShape& S, const real v[MAXDIRECTIONSIZE]) 
//-----------------------------------------------------------------------------
{
  m_rns = S;
  m_nsd = S.u_domainDimension_0u();
  for (int i=0; i<m_nsd; i++) {
    m_coor[i] = v[i];
  }
}

//-----------------------------------------------------------------------------
void RnPoint::u_setPoint_1u(const RnShape& S, const real v) 
//-----------------------------------------------------------------------------
{
  m_rns = S;
  m_nsd = S.u_domainDimension_0u();
  for (int i=0; i<m_nsd; i++) {
    m_coor[i] = v;
  }
}

//-----------------------------------------------------------------------------
void RnPoint::u_changeElement_1u(int dir, real index) 
//-----------------------------------------------------------------------------
{
  
  m_coor[dir] = index;

}

//-----------------------------------------------------------------------------
void RnPoint::u_setComponentsAll_1u(const real val) 
//-----------------------------------------------------------------------------
{
  for (int i=0; i<m_nsd; i++) {
    m_coor[i] = val;
  }
}

//-----------------------------------------------------------------------------
void RnPoint::u_setComponent_1u(const int dir, real val) 
//-----------------------------------------------------------------------------
{
  m_coor[dir-1] = val;
}

//-----------------------------------------------------------------------------
void RnPoint::u_setComponents_1u(const real v[MAXDIRECTIONSIZE]) 
//-----------------------------------------------------------------------------
{
 
  for ( int i = 0; i < m_nsd; i++) 
  { 
    m_coor[i] = v[i];
  }
}


// Vector arithmetic

//-----------------------------------------------------------------------------
void RnPoint::operator+= (const RnPoint& tp) 
//-----------------------------------------------------------------------------
{

  for ( int i = 0; i < m_nsd; i++)
  { 
    m_coor[i] = m_coor[i] + tp.m_coor[i];
  }
}   
 
//-----------------------------------------------------------------------------
void RnPoint::operator-= (const RnPoint& tp)
//-----------------------------------------------------------------------------
{

  for ( int i = 0; i < m_nsd; i++)
  { 
    m_coor[i] = m_coor[i] - tp.m_coor[i];
  }
}

//-----------------------------------------------------------------------------
void RnPoint::operator*= (const real v)
//-----------------------------------------------------------------------------
{

  for ( int i = 0; i < m_nsd; i++)
  { 
    m_coor[i] = m_coor[i] * v;
  }
}

//-----------------------------------------------------------------------------
void RnPoint::operator*=(const Array<int>& a)
//-----------------------------------------------------------------------------
// Componentwise multiplication with corresponding Array elements
{

  for ( int i = 0; i < m_nsd; i++)
  { 
    int d;
    a.u_getElement_0g(i,d);
    m_coor[i] = m_coor[i] * (real)d;
  }
}

//-----------------------------------------------------------------------------
void RnPoint::u_zero_1u()
//-----------------------------------------------------------------------------
{
  for (int i=0; i<m_nsd; i++) {
    m_coor[i] = 0.0;
  }
}

//-----------------------------------------------------------------------------
void RnPoint::u_one_1u()
//-----------------------------------------------------------------------------
{
  for (int i=0; i<m_nsd; i++) {
    m_coor[i] = 1.0;
  }
}



// Vector and ring arithmetic



//-----------------------------------------------------------------------------
//
//  u_norm2_0u - Calculate the vector length
//
//  DESCRIPTION:
//    Calculate the length of the vector.
//
//  RETURNS:
//    Returns the vector length.
//
//  NOTE:
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------
real RnPoint::u_norm2_0u() const
{
  real norm2=0;
  for (int i=0; i<m_nsd; i++) 
  {
    norm2 += m_coor[i]*m_coor[i];
  }
  norm2 = sqrt(norm2);
  return norm2;
}







//-----------------------------------------------------------------------------
//
//  normalize - Normalize the vector
//
//  DESCRIPTION:
//    Normalize the vector situated in this object.
//
//  RETURNS:
//    No return value.
//
//  NOTE:
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------
void RnPoint::normalize(void)
{
  
  real length = u_norm2_0u();
  for (int i=0; i < m_nsd; i++) 
  {
    m_coor[i] = m_coor[i] / length;
  }
  
  return;
}//normalize






//-----------------------------------------------------------------------------
real RnPoint::u_dot_0u(const RnPoint& tp) const
//-----------------------------------------------------------------------------
{
  real dot = 0;
  for ( int i = 0; i < m_nsd; i++)
  { 
    dot += m_coor[i] * tp.m_coor[i];
  }
  return dot;
}   







//-----------------------------------------------------------------------------
//
//  cross - Cross product of the vector
//
//  DESCRIPTION:
//    Make the cross product of this vector and the input vector. The 
//    result is kept in this object.
//
//              C
//             /\                    .
//            /  \                   .
//           /    \                  .
//          /      \                 .
//         /________\                .
//        A          B
//
//    vecAB x vecAC = [x1,y1,z1] x [x2,y2,z2] =
//                  = [y1z2 - z1y2, z1x2 - x1z2, x1y2 - y1x2]
//
//    vecAB is this object,
//    vecAC is the entered object
//
//  RETURNS:
//    No return value.
//
//  NOTE:
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------
void RnPoint::cross(const RnPoint& mp)
{
  
  real tmpCoord[MAXDIRECTIONSIZE];  // Temporary storage of this vector

  for (int i=0; i < m_nsd; i++) tmpCoord[i] = m_coor[i];

  m_coor[0] = tmpCoord[1] * mp.m_coor[2] - tmpCoord[2] * mp.m_coor[1];
  m_coor[1] = tmpCoord[2] * mp.m_coor[0] - tmpCoord[0] * mp.m_coor[2];
  m_coor[2] = tmpCoord[0] * mp.m_coor[1] - tmpCoord[1] * mp.m_coor[0];
  
  return;
}//cross


//**********************************************************
//NOTE: This function may only be suitable to use in 2D.
//**********************************************************
bool RnPoint::isCounterClockwise(const RnPoint& x1, const RnPoint& x2) const
{

  //*************************************************************
  //This routine checks wheter we have a counter-clock wise rotation
  // of the points with coords. x1 and x2 (NOTE:coordP and coordC just below)
  //with respect to the point with coord. coordQ,
  //i.e. the coordinates which ownes the method (see just below).
  //*************************************************************

  RnPoint coordQ = *this;
  RnPoint coordP = x1;
  RnPoint coordC = x2;

  bool CCW = false;

  //Must compute "the circle" radius:
  RnPoint help = coordP;
  help -= coordQ;
  real radius = help.u_norm2_0u();

  //Must locate and find the radius of the point S (see documentation):
  help = coordC;
  help -= coordQ;
  real LQC = help.u_norm2_0u();
  real t = radius/LQC;
  RnPoint coordS = coordC;
  coordS.u_changeElement_1u(0,(1.0-t)*coordQ.u_getElement_0u(0) + t*coordC.u_getElement_0u(0));
  coordS.u_changeElement_1u(1,(1.0-t)*coordQ.u_getElement_0u(1) + t*coordC.u_getElement_0u(1));

  //We perform a simple tranformation such that origo becomes the center of the circle:
  coordP -= coordQ;
  coordS -= coordQ;

  //We then decide which quadrants of the circel coordP and coordS belong to:
  int kvadrantP, kvadrantS;
  real cpx = coordP.u_getElement_0u(0);
  real cpy = coordP.u_getElement_0u(1);
  real csx = coordS.u_getElement_0u(0);
  real csy = coordS.u_getElement_0u(1);

  if ((cpx >= 0.0) && (cpy >= 0.0))
  {
    kvadrantP = 1;
  }
  else if ((cpx < 0.0) && (cpy >= 0.0))
  {
    kvadrantP = 2;
  }
  else if ((cpx < 0.0) && (cpy < 0.0))
  {
    kvadrantP = 3;
  }
  else if ((cpx >= 0.0) && (cpy < 0.0))
  {
    kvadrantP = 4;
  }

  if ((csx >= 0.0) && (csy >= 0.0))
  {
    kvadrantS = 1;
  }
  else if ((csx < 0.0) && (csy >= 0.0))
  {
    kvadrantS = 2;
  }
  else if ((csx < 0.0) && (csy < 0.0))
  {
    kvadrantS = 3;
  }
  else if ((csx >= 0.0) && (csy < 0.0))
  {
    kvadrantS = 4;
  }

  //Finally we can decide whether the rotation is counterclockwise or not:
  if (kvadrantP == kvadrantS)
  {
    if ((kvadrantP == 1) || (kvadrantP == 2))
    {
      if (cpx > csx) CCW = true;
    }
    else if ((kvadrantP == 3) || (kvadrantP == 4))
    {
      if (cpx < csx) CCW = true;
    }
  }
  else
  {
    //Here it is enough the consider those cases which makes CCW = true, since it is initialized to false:
    if ((kvadrantP == 1) && ((kvadrantS == 2) || (kvadrantS == 3)))
    {
      if ((kvadrantP == 1) && (kvadrantS == 2))
      {
	CCW = true;
      }
      else
      {//Special check must be performed: (S are in quadrant 3)
	if (fabs(cpy) >= fabs(csy)) CCW = true;//OK???
      }
    }
    else if ((kvadrantP == 2) && ((kvadrantS == 3) || (kvadrantS == 4)))
    {
      if ((kvadrantP == 2) && (kvadrantS == 3))
      {
	CCW = true;
      }
      else
      {//Special check must be performed: (S are in quadrant 4)
	if (fabs(cpx) >= fabs(csx)) CCW = true;//OK???
      }
    }
    else if ((kvadrantP == 3) && ((kvadrantS == 4) || (kvadrantS == 1)))
    {
      if ((kvadrantP == 3) && (kvadrantS == 4))
      {
	CCW = true;
      }
      else
      {//Special check must be performed: (S are in quadrant 1)
	if (fabs(cpy) >= fabs(csy)) CCW = true;//OK???
      }
    }
    else if ((kvadrantP == 4) && ((kvadrantS == 1) || (kvadrantS == 2)))
    {
      if ((kvadrantP == 4) && (kvadrantS == 1))
      {
	CCW = true;
      }
      else
      {//Special check must be performed: (S are in quadrant 2)
	if (fabs(cpx) >= fabs(csx)) CCW = true;//OK???
      }
    }
  }

  return CCW;
}



#ifdef COMPLETE_SOPHUS
//---------------------- Removed from -----------------------------------------


//-----------------------------------------------------------------------------
real RnPoint::u_angle_0u(const RnPoint& tp) const
//-----------------------------------------------------------------------------
{
 // Johannes Mykkeltveit, RF: blyant 26/2 2002.
 real ret;
 if(u_dot_0u(tp) <= TOLERANCE)
 {
   ret = M_PI_2;
 }
 else
 {
   real t = tp.u_norm2_0u() / u_dot_0u(tp);
    // atan2 have 2 double arguments, and returns double.
    // double atan2(double y, double x); (unix help)
    // The intension is that the return value is converted to real by '='.
   ret = atan2( norm2(tp * t - (*this)), u_norm2_0u());
 }
 return ret;
}


//---------------------- Removed to -------------------------------------------
#endif //COMPLETE_SOPHUS






//-----------------------------------------------------------------------------
void RnPoint::operator*= (const RnPoint& tp)
//-----------------------------------------------------------------------------
{

  for ( int i = 0; i < m_nsd; i++)
  { 
    m_coor[i] = m_coor[i] * tp.m_coor[i];
  }
}

//-----------------------------------------------------------------------------
void RnPoint::operator/= (const RnPoint& tp)
//-----------------------------------------------------------------------------
{

  for ( int i = 0; i < m_nsd; i++)
  { 
    m_coor[i] = m_coor[i] / tp.m_coor[i];
  }
}



// Embedding
//-----------------------------------------------------------------------------
void RnPoint::u_embed_1u(const RnShape& S)
//-----------------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
void RnPoint::u_embed_1u(const RnPoint& Q)
//-----------------------------------------------------------------------------
{
 
  for ( int i = 0; i < m_nsd; i++)
  { m_coor[i] = m_coor[i] + Q.m_coor[i];
  }
}

// String Operations

//-----------------------------------------------------------------------------
void RnPoint::u_concatenate_1u (const RnPoint& Q)
//-----------------------------------------------------------------------------
{
    m_rns.u_concatenate_1u(Q.m_rns);
    int max = m_nsd+Q.m_nsd;
    for ( int i = m_nsd; i < max; i++) m_coor[i] = Q.m_coor[i];
    m_nsd = m_rns.u_domainDimension_0u();
} 

//----------------------------------------------------------------------------
void RnPoint::u_extract_1u(const int& l,const int& r)
//----------------------------------------------------------------------------
{

  for (int i=l; i<r; i++) 
    m_coor[i-l] = m_coor[i];
  
  m_nsd = r-l;
  m_rns.u_extract_1u(l,r);
}

// Conventions

//-----------------------------------------------------------------------------
void RnPoint::operator= (const RnPoint& tp)
//-----------------------------------------------------------------------------
{
    m_rns = tp.m_rns;
    m_nsd = tp.m_nsd;
    for ( int i = 0; i < m_nsd; i++) m_coor[i] = tp.m_coor[i];
} 

//-----------------------------------------------------------------------------
RnPoint::~RnPoint()
//-----------------------------------------------------------------------------
{
}
 
//-----------------------------------------------------------------------------
Boolean RnPoint::operator==(const RnPoint& rp) const
//-----------------------------------------------------------------------------
{
  Boolean is_equal = (m_rns == rp.m_rns);
  for ( int i = 0; i < m_nsd; i++) 
  { 
    if (fabs(m_coor[i]-rp.m_coor[i]) > TOLERANCE)
    { 
      is_equal = false;
      break;
    }
  }
  return is_equal;
}       



//-----------------------------------------------------------------------------
Boolean RnPoint::u_DI_0u() const
//-----------------------------------------------------------------------------
{

  Boolean ok = true;
 
  if ( m_nsd != m_rns.u_domainDimension_0u() )
    ok = false;
  else
  { if ( (m_nsd < 0) || (MAXDIRECTIONSIZE < m_nsd) )
      ok = false;
  }
  return ok;
}





//-----------------------------------------------------------------------------
void RnPoint::u_copy_1u()
//-----------------------------------------------------------------------------
{

  m_rns = RnShape();
  m_nsd = m_rns.u_domainDimension_0u();

}

//-----------------------------------------------------------------------------
void RnPoint::u_copy_1u(const RnShape& ts)
//-----------------------------------------------------------------------------
{
  m_rns = ts;
  m_nsd = ts.u_domainDimension_0u();
  for (int i = 0; i < m_nsd; i++) m_coor[i] = 0.0;
}







//---------------------- Protected Functions ----------------------------------



//---------------------- Private Functions ------------------------------------








/*
===============================================================================
    ------------------------ CLASS DESCRIPTION --------------------------
===============================================================================

  DESCRIPTION:
  ------------------
    <Description applicable regarding implementation, modifications and
     corrections etc>



  ENHANCEMENTS:
  -------------

    
    
  MISCELLANEOUS:
  --------------
  

===============================================================================
    --------------------- END OF CLASS DESCRIPTION ----------------------
===============================================================================
*/


