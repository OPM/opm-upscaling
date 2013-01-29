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
#include "MeshPoint.h"


//---------------------  Constants --------------------------------------------


//---------------------- Types ------------------------------------------------












//---------------------- Public Functions -------------------------------------



// Data invariant

//-----------------------------------------------------------------------------
Boolean MeshPoint::u_DI_0u() const
//-----------------------------------------------------------------------------
{
  Boolean ok = true;
  int i;
 
  if ( nsd != shape.u_domainDimension_0u())
    ok = false;
  else
  {
    if ( nsd < 0)
      ok = false;
    else {
      for ( i = 0; i < nsd; i++) 
	{
	  if (shape.u_getDirectionSize_0u(i+1) == 0 ) {
              if ( ind[i] != 0)
		{
		  ok = false;
		  break;
		}
	    }
	  else if ( ind[i] < 0  || ind[i] >= shape.u_getDirectionSize_0u(i+1) )
            {
              ok = false;
              break;
	    }
	}
    }
  }
  return ok;
} 


// Constructors

//-----------------------------------------------------------------------------
MeshPoint::MeshPoint()
//-----------------------------------------------------------------------------
{
  nsd = shape.u_domainDimension_0u();
}

//-----------------------------------------------------------------------------
MeshPoint::MeshPoint(const MeshPoint& mp)   // MeshPoint1
: nsd(mp.nsd), shape(mp.shape)
//-----------------------------------------------------------------------------
{
    for (int i = 0; i < nsd; i++) ind[i] = mp.ind[i];
}


//-----------------------------------------------------------------------------
MeshPoint::MeshPoint(const MeshShape& ms)   // MeshPoint2
: nsd(ms.u_domainDimension_0u()), shape(ms)
//-----------------------------------------------------------------------------
{
    for (int i = 0; i < nsd; i++) ind[i] = 0;
}


//-----------------------------------------------------------------------------
MeshPoint::MeshPoint(const MeshShape& ms, const int& def)   // MeshPoint3
: nsd(ms.u_domainDimension_0u()), shape(ms)
//-----------------------------------------------------------------------------
{
    for (int i = 0; i < nsd; i++) ind[i] = def;
}


//----------------------------------------------------------------------------
void MeshPoint::u_copy_1u()     // copy0
//----------------------------------------------------------------------------
{
  nsd = shape.u_domainDimension_0u();;
}

//-----------------------------------------------------------------------------
void MeshPoint::u_copy_1u(const MeshShape& ms)   // copy2
//-----------------------------------------------------------------------------
{
    shape = ms;
    nsd = ms.u_domainDimension_0u();
    for (int i = 0; i < nsd; i++) ind[i] = 0;
}


//-----------------------------------------------------------------------------
void MeshPoint::u_copy_1u(const MeshShape& ms, const int& def)   // copy3
//-----------------------------------------------------------------------------
{
    shape = ms;
    nsd = ms.u_domainDimension_0u();
    for (int i = 0; i < nsd; i++) ind[i] = def;
}



// Observers 
//----------------------------------------------------------------------------
int MeshPoint::u_domainDimension_0u() const
//----------------------------------------------------------------------------
{
  return nsd;
}
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
MeshShape MeshPoint::u_getShape_0u() const
//----------------------------------------------------------------------------
{
  return shape;
}

//----------------------------------------------------------------------------
Boolean MeshPoint::u_sameShape_0u(const MeshPoint& Q) const
//---------------------------------------------------------------------
{
  Boolean result = shape == Q.shape;
  return result;

}

//-----------------------------------------------------------------------------
Boolean MeshPoint::u_isDivisible_0u(const MeshPoint& mp) const
//-----------------------------------------------------------------------------

{

  Boolean is_div = true;
  int a, b, i, inv_ind, limit;
  if(mp.u_isOrigin_0u() == true){
    is_div = false;
  } else {
    for ( i = 0; i < nsd; i++) {
      a = ind[i];
      b = mp.ind[i];
      limit = shape.u_getDirectionSize_0u(i+1);
      inv_ind = ul_isDiv(a,b,limit);
      if ( inv_ind < 0) {
        is_div = false;
        break;
      } 
    }
  }
  return is_div;
}


//----------------------------------------------------------------------------
Boolean MeshPoint::u_legalDirection_0u(int dir) const
//----------------------------------------------------------------------------
{
  Boolean result = dir > 0 && dir <= nsd;
  return result;
}

//----------------------------------------------------------------------------
Boolean MeshPoint::u_indexOK_0u(const int ind) const
//----------------------------------------------------------------------------
{
   Boolean ret = (ind >= 0) && (ind < nsd);
   return ret;
}

//-----------------------------------------------------------------------------
int MeshPoint::u_getSize_0u() const
//-----------------------------------------------------------------------------
{
  return nsd;
}

//-----------------------------------------------------------------------------
void MeshPoint::u_getElement_0g(const int dir, int& e) const
//-----------------------------------------------------------------------------
{
  e = ind[dir];
}

//-----------------------------------------------------------------------------
int MeshPoint::u_getComponent_0u(int dir) const
//-----------------------------------------------------------------------------
{
  return ind[dir-1];
}

//-----------------------------------------------------------------------------
void MeshPoint::u_getComponents_0g(int cp[MAXDIRECTIONSIZE]) const
//-----------------------------------------------------------------------------
{
  for(int i = 0; i < nsd; i++)
    cp[i] = ind[i];
}

//-----------------------------------------------------------------------------
void MeshPoint::u_concatenate_1u(const MeshPoint& Q) 
//-----------------------------------------------------------------------------
{
  MeshShape qms = Q.u_getShape_0u();
  shape.u_concatenate_1u(qms);
  
  int qnsd = qms.u_domainDimension_0u();
  int oldnsd = nsd;
  nsd = shape.u_domainDimension_0u();    
  for(int i=0;i<qnsd;i++) 
    ind[i+oldnsd] = Q.u_getComponent_0u(i+1);
}

//----------------------------------------------------------------------------
void MeshPoint::u_extract_1u(const int& l,const int& r)
//----------------------------------------------------------------------------
{
  for (int i=l; i<r; i++) 
    ind[i-l] = ind[i];
  
  nsd = r-l;
  shape.u_extract_1u(l,r);
}


//-----------------------------------------------------------------------------
Boolean MeshPoint::u_isOrigin_0u() const
//-----------------------------------------------------------------------------
{
  Boolean is_origin = true;
  int i;

  for ( i= 0; i < nsd; i++) {
    if ( ind[i] != 0) {
      is_origin = false;
      break;
    }
  }
  return is_origin;
}


// Generators

//-----------------------------------------------------------------------------
void MeshPoint::u_setShape_1u(const MeshShape& S) 
//-----------------------------------------------------------------------------
{
  nsd = S.u_domainDimension_0u();
  shape = S; 
  for (int i=0; i<nsd; i++) {
    ind[i] = 0;
  }
}

//-----------------------------------------------------------------------------
void MeshPoint::u_setPoint_1u(const MeshShape& S,const int v[MAXDIRECTIONSIZE]) 
//-----------------------------------------------------------------------------
{
  nsd = S.u_domainDimension_0u();
  shape = S;
  int newind,limit;
  for (int i=0; i<nsd; i++) {
    newind = v[i];
    limit = shape.u_getDirectionSize_0u(i+1);
    while(newind<0)
      newind+=limit;
    if(newind>=limit && limit!=0)
      newind %= limit;
    ind[i] = newind;
  }
}

//-----------------------------------------------------------------------------
void MeshPoint::u_setPoint_1u(const MeshShape& S,int v) 
//-----------------------------------------------------------------------------
{
  nsd = S.u_domainDimension_0u();
  shape = S;
  int newind,limit;
  for (int i=0; i<nsd; i++) {
    newind = v;
    limit = shape.u_getDirectionSize_0u(i+1);
    while(newind<0)
      newind+=limit;
    if(newind>=limit && limit!=0)
      newind %= limit;
    ind[i] = newind;
  }
}

//-----------------------------------------------------------------------------
void MeshPoint::u_changeElement_1u(const int dir, int e)
//-----------------------------------------------------------------------------
{
//  27.01.2004 GIE/HAF: Original Sophus code
//  int limit = shape[dir];  
  int limit;
  shape.u_getElement_0g(dir, limit);

  while(e < 0)
    e+=limit;
  if(e >= limit && limit!=0)
    e %= limit;
  ind[dir] = e;
}

//-----------------------------------------------------------------------------
void MeshPoint::u_setComponentsAll_1u(int val) 
//-----------------------------------------------------------------------------
{
  int i, limit;
  for ( i = 0; i < nsd; i++) {
    limit = shape.u_getDirectionSize_0u(i+1);
    ind[i] = val;
    while ( ind[i] < 0 )
      ind[i] += limit;
    if (ind[i] >= limit && limit!=0 )
      ind[i] = ind[i] % limit;
  }
}


//-----------------------------------------------------------------------------
void MeshPoint::u_setComponent_1u(int dir, int index) 
//-----------------------------------------------------------------------------
{
  int limit = shape.u_getDirectionSize_0u(dir);

  while(index < 0)
    index+=limit;
  if(index >= limit && limit!=0)
    index %= limit;
  ind[dir-1] = index;
}

//-----------------------------------------------------------------------------
void MeshPoint::u_setComponents_1u(const int v[MAXDIRECTIONSIZE]) 
//-----------------------------------------------------------------------------
{
  int i, limit,newind;
  for ( i = 0; i < nsd; i++) {
    limit = shape.u_getDirectionSize_0u(i+1);
    newind = v[i];
    while(newind < 0)
      newind += limit;
    if(newind>= limit && limit!=0)
      newind %= limit;
    ind[i] = newind;
  }
}

// Composite Operations

//-----------------------------------------------------------------------------
int MeshPoint::u_getLexicographic_0u() const 
//-----------------------------------------------------------------------------
{
  int factor[MAXDIRECTIONSIZE];
  int i,k;

  for (i=1; i<=nsd; i++) {
    factor[i-1]=1;
    for (k=i+1; k<=nsd; k++) {
      factor[i-1] *=shape.u_getDirectionSize_0u(k);
    }
  }
  int result=0;
  for (i=1; i<=nsd; i++) {
    result += u_getComponent_0u(i)*factor[i-1];
  }
  return result;
}
//-----------------------------------------------------------------------------
void MeshPoint::u_setLexicographic_1u(int cp) 
//-----------------------------------------------------------------------------
{
  int i;  
  int b=cp;
  int limit;
  
  if(shape.u_domainDimension_0u()!=0) {

    for (i=shape.u_domainDimension_0u();i>0;i--) {
      limit = shape.u_getDirectionSize_0u(i);
      int icoord = b % limit;
      u_setComponent_1u(i,icoord);
      b /=limit;
    }
  }
}
// Vector arithmetic

//-----------------------------------------------------------------------------
void MeshPoint::operator += (const MeshPoint& mp) 
//-----------------------------------------------------------------------------
{
  int i, limit;
  for ( i = 0; i < nsd; i++) {
    limit = shape.u_getDirectionSize_0u(i+1);
    ind[i] += mp.ind[i];
    ind[i] %= limit;
  }
}   
 
//-----------------------------------------------------------------------------
void MeshPoint::operator -= (const MeshPoint& mp)
//-----------------------------------------------------------------------------
{
  int i, limit;
    for ( i = 0; i < nsd; i++)
    {
      limit = shape.u_getDirectionSize_0u(i+1);
      ind[i] -= mp.ind[i];
      while (ind[i] < 0)
	ind[i] += limit;
      if(ind[i]>=limit && limit!=0)
	ind[i] %= limit;
    }
}

//-----------------------------------------------------------------------------
void MeshPoint::operator *= (int v)
//-----------------------------------------------------------------------------
{
  int limit;
  int i;
  for (i=0; i<nsd; i++) {
    limit = shape.u_getDirectionSize_0u(i+1);
    ind[i] = ind[i]*v ; //+ limit;
    while (ind[i]<0) 
      ind[i] += limit;
    if(ind[i]>=limit && limit!=0)
      ind[i] %= limit;
  }
}

//-----------------------------------------------------------------------------
void MeshPoint::u_zero_1u()
//-----------------------------------------------------------------------------
{
  for (int i=0; i<nsd; i++) {
    ind[i] = 0;
  }
}

// Vector and ring arithmetic

//-----------------------------------------------------------------------------
void MeshPoint::operator *= (const MeshPoint& mp)
//-----------------------------------------------------------------------------
{
  int i, limit;
    for ( i = 0; i < nsd; i++)
    {
      limit = shape.u_getDirectionSize_0u(i+1);
      ind[i] *= mp.ind[i];
      if(ind[i]>=limit && limit!=0)
	ind[i] %= limit;
    }
}

//-----------------------------------------------------------------------------
void MeshPoint::operator /= (const MeshPoint& mp)
//-----------------------------------------------------------------------------
{
  for (int i = 0; i < nsd; i++) {  //TO BE CHANGED
    ind[i] /= mp.ind[i];
  }
}

// Embedding

//-----------------------------------------------------------------------------
void MeshPoint::u_embed_1u(const MeshShape& S)
//---This is u_embed_1u, overloading
//-----------------------------------------------------------------------------

{
  shape = S;

  int i;
  int limit;
  
  for (i=0;i<nsd;i++){
    limit = shape.u_getDirectionSize_0u(i+1);
    if (ind[i]>=limit && limit!=0)
      ind[i] %= limit;
  }
}

//-----------------------------------------------------------------------------
void MeshPoint::u_borderEmbed_1u(const MeshShape& S,const MeshShape& SP)
//-----------------------------------------------------------------------------
{
  int n=nsd-1;
  for(int i=S.u_domainDimension_0u();i>0;i--) {
    if(SP.u_getDirectionSize_0u(i)==0)
      ind[i-1]=0;
    else 
      ind[i-1]=ind[n--];
  }

  shape = S;
  nsd = S.u_domainDimension_0u();
}

//-----------------------------------------------------------------------------
void MeshPoint::u_embed_1u(const MeshPoint& Q)
//---This is uembed_2
//-----------------------------------------------------------------------------

{
  for (int i=0; i<nsd; i++) {
    ind[i] += Q.ind[i]; 

    int limit = shape.u_getDirectionSize_0u(i+1);
    if (ind[i] >= limit && limit!=0) {
      ind[i] = ind[i] % limit;
    }
  }
}

// Scaling

//-----------------------------------------------------------------------------
void MeshPoint::u_scale_1u(const MeshShape& S)
//-----------------------------------------------------------------------------
{
  shape.u_scale_1u(S);		
  for (int i=0; i<nsd; i++) {
    int limit = shape.u_getDirectionSize_0u(i+1);
    ind[i] = ind[i]*S.u_getDirectionSize_0u(i+1);
    if (ind[i] >= limit && limit!=0) {
      ind[i] = ind[i] % limit;
    }
  }
}

// Mirroring

//-----------------------------------------------------------------------------
void MeshPoint::u_mirror_1u(const MeshShape& S)
//-----------------------------------------------------------------------------
{
  for (int i=0; i<nsd; i++) {
    if (S.u_getDirectionSize_0u(i+1) == 0) {
      int limit = shape.u_getDirectionSize_0u(i+1);
      ind[i] = -ind[i];
      if(ind[i]<0)
        ind[i] += limit;
    }
  }
}

//-----------------------------------------------------------------------------
void MeshPoint::u_mirror_1u(const MeshPoint& Q, const MeshShape& S)
//-----------------------------------------------------------------------------
{
  for (int i=0; i<nsd; i++) {
    if (S.u_getDirectionSize_0u(i+1) == 0) {
      int limit = shape.u_getDirectionSize_0u(i+1);
      ind[i] = ind[i] - Q.ind[i];
      ind[i] = -ind[i]; // Se s.16 i spesifikasjonen
      ind[i] = ind[i] + Q.ind[i];
      while (ind[i] < 0) ind[i] += limit;
      if (ind[i] >= limit) 
	ind[i] %= limit;
    }
  }
}


// Conventions

//-----------------------------------------------------------------------------
void MeshPoint::operator = (const MeshPoint& mp)
//-----------------------------------------------------------------------------
{
  int i;
  
  shape = mp.shape;
  nsd = mp.nsd;

  for ( i = 0; i < nsd; i++) ind[i] = mp.ind[i];
} 

//-----------------------------------------------------------------------------
MeshPoint::~MeshPoint()
//-----------------------------------------------------------------------------
{
}
 
//-----------------------------------------------------------------------------
Boolean MeshPoint::operator==(const MeshPoint& mp) const
//-----------------------------------------------------------------------------
{
  Boolean is_equal = true;
  int i; 

  if (!(shape == mp.shape))
    is_equal = false;
  else
    for ( i = 0; i < nsd; i++) {
      if (ind[i] != mp.ind[i]) {
	is_equal = false;
	break;
      }
    }
  return is_equal;
}




//-----------------------------------------------------------------------------
int MeshPoint::operator[](const int& dir) const
//-----------------------------------------------------------------------------
{
  return ind[dir];
}       






// private            

//-----------------------------------------------------------------------------
void MeshPoint::ul_gcd( int x, int n, int& g, int& c, int& d) const
//-----------------------------------------------------------------------------
{
  int u[3], v[3], t[3];
  int q, i;
  if ( n > 0 &&  x > 0)
  {
    u[0] = 1;
    u[1] = 0;
    u[2] = x;
    v[0] = 0;
    v[1] = 1;
    v[2] = n;
    for (i = 0; i < 3; i++)
      t[i] = 0;
    while (v[2] != 0)
    {
      q = u[2]/v[2];
      if (u[2] * v[2] < 0) 
        q -= 1; // fake floor function
      for ( i = 0; i < 3; i++) 
      {
        t[i] = u[i] - q * v[i];
        u[i] = v[i];
        v[i] = t[i];
      }
    }
    c = u[0];
    d = u[1];
    g = u[2];
  }
  else 
  {
    if ( x == 0 && n == 0) 
    {
      g = 0;
      c = 0;
      d = 0;
    } 
    else
    {
      if ( x == 0) 
      {
        g = n;
        c = 0;
        d = 1;
      }
      else
      {
        g = x;
        c = 1;
        d = 0;
      }
    }
  }
}


//-----------------------------------------------------------------------------
int MeshPoint::ul_isDiv( int a, int b, int limit) const
//-----------------------------------------------------------------------------
{
  int c, d, g;
 
  ul_gcd(a,b,g,c,d);
  if ( g >= 1) 
  {
    if ( g > 1) 
    {
      a /= g;
      b /= g;
    }
    ul_gcd(b, limit, g, c, d);
    if ( g == 1)
    {
      if ( c >= limit)
        c %= limit;
      while ( c < 0)
        c+= limit;
      c *= a;
      c %= limit;
    }
    else
      c = -1;
  }
  else
    c = -1;
  return c;
}

//template<class codeboost>
//void MeshPoint::rules()
//{
//  MeshPoint p;
//  int ii;
//  topdown: getLexicographic(setLexicographic(p,ii)) = ii;
//}



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
