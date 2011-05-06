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
#include "MeshShape.h"

#include <climits>



//---------------------  Constants --------------------------------------------


//---------------------- Types ------------------------------------------------





//---------------------- Public Functions -------------------------------------






//----------------------------------------------------------------------------
Boolean MeshShape::u_DI_0u() const
//----------------------------------------------------------------------------
{
  Boolean is_ok = true;

  if ( nsd < 0 || nsd > MAXDIRECTIONSIZE ) { 
    is_ok = false;
  }
  else {
    for ( int i = 0; i < nsd; i++) {
      if ( v[i] < 0 ) {
          is_ok = false;
         break;
      }
    }
  }
  return is_ok;
}


//----------------------------------------------------------------------------
int MeshShape::u_domainDimension_0u() const
//----------------------------------------------------------------------------
{
  return nsd;
}

//----------------------------------------------------------------------------
void MeshShape::u_getElement_0g(const int dir, int& e) const
//----------------------------------------------------------------------------
{
  e = v[dir];
}

//----------------------------------------------------------------------------
int MeshShape::u_getDirectionSize_0u(int dir) const
//----------------------------------------------------------------------------
{
  return v[dir-1];
}

//----------------------------------------------------------------------------
Boolean MeshShape::u_indexOK_0u(const int dir) const
//----------------------------------------------------------------------------
{
   Boolean ret = (dir >= 0) && (dir < nsd);
   return ret;
}

//------------------------------------------------------------------------------
int MeshShape::u_getMaxDirectionVolume_0u() const
//------------------------------------------------------------------------------
{
  int max_val = 0;
    for ( int i = 0; i < nsd; i++)
    {
      if ( max_val < v[i])
        max_val = v[i];
    }
  return max_val;
}

//------------------------------------------------------------------------------
int MeshShape::u_getMinDirectionVolume_0u() const
//------------------------------------------------------------------------------
{
    int min_val=INT_MAX;
    for ( int i = 0; i < nsd; i++)
    {
      if ( min_val > v[i])
        min_val = v[i];
    }
  return min_val;
}
    

// Constructors

//----------------------------------------------------------------------------
MeshShape::MeshShape()     // MeshShape0
  : nsd(0)
//----------------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
MeshShape::MeshShape(const MeshShape& ms) // MeshShape1
  : nsd(ms.nsd)
//-----------------------------------------------------------------------------
{
  for ( int i = 0; i < nsd; i++) {
    v[i] = ms.v[i];
  }
}


//-----------------------------------------------------------------------------
MeshShape::MeshShape(const int& size)   // MeshShape2
  : nsd(size)
//-----------------------------------------------------------------------------
{
  for ( int i = 0; i < nsd; i++) {
    v[i] = 0;
  }
}


//-----------------------------------------------------------------------------
MeshShape::MeshShape(const int& size, const int& def)   // MeshShape3
  : nsd(size)
//-----------------------------------------------------------------------------
{
  for ( int i = 0; i < nsd; i++) {
    v[i] = def;
  }
}


//----------------------------------------------------------------------------
void MeshShape::u_copy_1u()     // copy0
//----------------------------------------------------------------------------
{
  nsd = 0;
}

//-----------------------------------------------------------------------------
void MeshShape::u_copy_1u(const int& size)   // copy2
//-----------------------------------------------------------------------------
{
  nsd = size;
  for ( int i = 0; i < nsd; i++) {
    v[i] = 0;
  }
}


//-----------------------------------------------------------------------------
void MeshShape::u_copy_1u(const int& size, const int& def)   // copy3
//-----------------------------------------------------------------------------
{
  nsd = size;
  for ( int i = 0; i < nsd; i++) {
    v[i] = def;
  }
}


// Observers

//-----------------------------------------------------------------------------
Boolean MeshShape::u_sameShape_0u(const MeshShape& msm) const
//-----------------------------------------------------------------------------
{
  Boolean result = nsd == msm.nsd;
  return result;

}

//-----------------------------------------------------------------------------
Boolean MeshShape::u_legalDirection_0u(int dir) const
//-----------------------------------------------------------------------------
{
  Boolean result = (dir >= 1 && dir <= nsd);
  return result;

}

//-----------------------------------------------------------------------------
void MeshShape::u_getDirectionSizes_0g(int cp[MAXDIRECTIONSIZE]) const
//-----------------------------------------------------------------------------
{
  //rettet av Helge i samråd med Magne & Steinar
  for (int i = 0; i < nsd; i++)
    cp[i] = v[i];
}


//-----------------------------------------------------------------------------
Boolean MeshShape::u_inShape_0u(int vn) const
//-----------------------------------------------------------------------------
{
  Boolean inshape=true;
  if(nsd==0)
    inshape = false;
  if(vn<0)
    inshape = false;
  int i=0;
  while(i<nsd && inshape) {
    if(vn> v[i])
      inshape = false;
    i++;
  }
  return inshape;
}

//-----------------------------------------------------------------------------
Boolean MeshShape::u_inShape_0u(const int vn[MAXDIRECTIONSIZE]) const
//-----------------------------------------------------------------------------
{
  Boolean inshape=true;
  int i=0;
  while(i<nsd && inshape) {
    if(vn[i] > v[i] || vn[i] < 0)
      inshape = false;
    i++;
  }
  return inshape;
}


// Observers

//-----------------------------------------------------------------------------
int MeshShape::u_getSize_0u() const
//-----------------------------------------------------------------------------
{
  int result = u_volume_0u();
  return result;
}


// Generators

//-----------------------------------------------------------------------------
void MeshShape::u_setShape_1u(int n)
//-----------------------------------------------------------------------------
{
  nsd = n;
  for ( int i = 0; i < nsd; i++) {
    v[i] = 0; //Litt i tvil her hvis n=0
  }
}

//-----------------------------------------------------------------------------
void MeshShape::u_setShape_1u(int dim,int sizes)
//-----------------------------------------------------------------------------
{
  nsd = dim;
  for ( int i = 0; i < nsd; i++) {
    v[i] = sizes;
  }
}

//-----------------------------------------------------------------------------
void MeshShape::u_setShape_1u(int dim,const int sizes[MAXDIRECTIONSIZE])
//-----------------------------------------------------------------------------
{
  nsd = dim;
  for ( int i = 0; i < nsd; i++) {
    v[i] = sizes[i];
  }
}

//-----------------------------------------------------------------------------
void MeshShape::u_changeElement_1u(const int & dir, const int & n_i)
//-----------------------------------------------------------------------------
{
  v[dir] = n_i;
} 

//-----------------------------------------------------------------------------
void MeshShape::u_setDirectionSizeAll_1u(int I)
//-----------------------------------------------------------------------------
{
  for (int i=0;i<nsd;i++){
    v[i] = I;
  }
}

//-----------------------------------------------------------------------------
void MeshShape::u_setDirectionSize_1u(int dir, int n_i)
//-----------------------------------------------------------------------------
{
  v[dir-1] = n_i;
} 
 
//-----------------------------------------------------------------------------
void MeshShape::u_setDirectionSizes_1u(const int vp[MAXDIRECTIONSIZE])
//-----------------------------------------------------------------------------
{
  for (int i=0;i<nsd;i++){
    v[i] = vp[i];
  }
}

// Composite

//----------------------------------------------------------------------------
void MeshShape::u_concatenate_1u(const MeshShape& msm)
//----------------------------------------------------------------------------
{ 
  for (int i=0; i<msm.nsd;i++)
    v[i+nsd]=msm.v[i]; 
  nsd=nsd+msm.nsd;
}


//----------------------------------------------------------------------------
void MeshShape::u_extract_1u(const int& l,const int& r)
//----------------------------------------------------------------------------
{
  for (int i=l; i<r; i++) 
    v[i-l] = v[i];
  
  nsd = r-l;
}

// Volumes

//----------------------------------------------------------------------------
int MeshShape:: u_getDirectionVolume_0u(int d) const
//----------------------------------------------------------------------------
// direction d
{
  return v[d-1];
}

//----------------------------------------------------------------------------
int MeshShape::u_volume_0u() const
//----------------------------------------------------------------------------
{
  int nno = 1;
  for ( int i = 0; i < nsd; i++) 
    nno *= v[i];
  return nno;
}

// Subshapes

//----------------------------------------------------------------------------
Boolean MeshShape::u_subShape_0u(const MeshShape& msm) const
//----------------------------------------------------------------------------
{
  Boolean issubshape;
  if (u_sameShape_0u(msm)){
    issubshape = true;
    for (int i=0;i<nsd;i++){ 
      issubshape = (issubshape && (v[i] >= msm.v[i])); 
      if (!issubshape) break;
    }
  } else {			// Not a subshape if the number
    issubshape = false;         // of directions are not equal
  }

  return issubshape;
}


// Scale factors

//----------------------------------------------------------------------------
Boolean MeshShape::u_isScale_0u() const
//----------------------------------------------------------------------------
{
  Boolean notempty;
  notempty = true;
  for (int i=0;i<nsd;i++) {
    notempty = notempty && (v[i] != 0);  
    if (!notempty) break;
  }
  return notempty;
}

//----------------------------------------------------------------------------
Boolean MeshShape::u_isScalable_0u(const MeshShape& SP) const
//----------------------------------------------------------------------------
{
  Boolean scalable=true;

  // Steinar: I changed some here - returning before the end makes the RETTRACE
  //          not work.

  if (u_sameShape_0u(SP) && SP.u_isScale_0u() && u_subShape_0u(SP)) {

    // finds out if v[i]/SP.v[i] is an integer
    int i=0;
    while(i<nsd && scalable) {
      if (v[i] % SP.v[i] != 0) 
	scalable = false;
      i++;
    }
  }
  else
    scalable = false;
  return scalable;
}

//----------------------------------------------------------------------------
void MeshShape::u_scale_1u(const MeshShape& SP)
//----------------------------------------------------------------------------
{
  for (int i=0;i<nsd;i++) {
    v[i] = v[i]*SP.v[i];
  }
}

//----------------------------------------------------------------------------
void MeshShape::u_getScaleFactor_1u(const MeshShape& SP)
//----------------------------------------------------------------------------
{
  for (int i=0;i<nsd;i++) {
    v[i] = v[i]/SP.v[i];
  }
}

// Mirrors

//----------------------------------------------------------------------------
Boolean MeshShape::u_isMirror_0u(const MeshShape& msm) const
//----------------------------------------------------------------------------
{
  Boolean result = u_sameShape_0u(msm);
  return result;

}

//----------------------------------------------------------------------------
int MeshShape::u_mirrorShape_0u() const
//----------------------------------------------------------------------------
{
  int result = 0;
  for (int i=0;i<nsd;i++) {
    if (v[i] != 0){
	 result++;
    }
  }
  return result;

}

// Borders

//----------------------------------------------------------------------------
Boolean MeshShape::u_isBorder_0u(const MeshShape& msm) const
//----------------------------------------------------------------------------
{
  Boolean result = u_sameShape_0u(msm);
  return result;

}


//----------------------------------------------------------------------------
Boolean MeshShape::u_isEmbeddableBorder_0u(const MeshShape& msm) const
//----------------------------------------------------------------------------
{
  Boolean result = u_sameShape_0u(msm) && u_subShape_0u(msm);
  return result;

}


//----------------------------------------------------------------------------
int MeshShape::u_borderShape_0u() const
//----------------------------------------------------------------------------
{
  int result = 0;
  for (int i=0;i<nsd;i++) {
    if (v[i] != 0) result++;
  }
  return result;

}

//----------------------------------------------------------------------------
void MeshShape::u_getBorder_1u(const MeshShape& msm)
//----------------------------------------------------------------------------
{
  // finds nsd for the boundary shape
  nsd = msm.u_borderShape_0u();
  int i,k=0;
  for (i=0;i<msm.nsd;i++) {
    if (msm.v[i] != 0) {
      v[k] = v[i];
      k++;
    }
  }
}





// Conventions

//---------------------------------------------------------------------------
void MeshShape::operator=(const MeshShape& ms)
//---------------------------------------------------------------------------
{
  if ( this != &ms) 
  {
    nsd = ms.nsd;
    if (nsd > 0) {
      for (int i = 0; i < nsd; i++)
	v[i] = ms.v[i];
    }
  }
}

//----------------------------------------------------------------------------
Boolean MeshShape::operator==(const MeshShape& ms) const
//----------------------------------------------------------------------------
{
  Boolean correct = true;
  if (nsd != ms.nsd)
    correct = false;
  else 
  {
    if ( nsd > 0 )
    {
      for ( int i = 0; i < nsd; i++)
      {
        if ( v[i] != ms.v[i]) 
        {
          correct = false;
          break;
        }
      }
    }
  }
  return correct;
}

//----------------------------------------------------------------------------
MeshShape::~MeshShape()
//----------------------------------------------------------------------------
{
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
  
  NOTE: MAXDIRECTIONSIZE must be defined before including this file
  --- End ---

  
  


===============================================================================
    --------------------- END OF CLASS DESCRIPTION ----------------------
===============================================================================
*/

