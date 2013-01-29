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
#include "RnShape.h"



//---------------------  Constants --------------------------------------------


//---------------------- Types ------------------------------------------------




//---------------------- Public Functions -------------------------------------





//----------------------------------------------------------------------------
Boolean RnShape::u_DI_0u() const
//----------------------------------------------------------------------------
{

  Boolean is_ok = true;

  if ( (nsd < 0) || (nsd > MAXDIRECTIONSIZE) ) { 
    is_ok = false;
  }

  return is_ok;
}


// Constructors

//-----------------------------------------------------------------------------
RnShape::RnShape()     // RnShape0
  : nsd(0)
//-----------------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
RnShape::RnShape(const RnShape & S) // RnShape1
  : nsd(S.nsd)
//-----------------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
RnShape::RnShape(const int n) // RnShape2
  : nsd(n)
//-----------------------------------------------------------------------------
{
}


// Observers

//----------------------------------------------------------------------------
int RnShape::u_domainDimension_0u() const
//----------------------------------------------------------------------------
{ 
  return nsd;
}


//-----------------------------------------------------------------------------
Boolean RnShape::u_sameShape_0u(const RnShape& tsm) const
//-----------------------------------------------------------------------------
{
  Boolean OK;
  if (nsd == tsm.nsd) {
    OK = true;
  }
  else {
    OK = false;
  }
  return OK;
}

//-----------------------------------------------------------------------------
Boolean RnShape::u_legalDirection_0u(const int dir) const
//-----------------------------------------------------------------------------
{
  Boolean legal;
  if (dir >= 1 && dir <= nsd) {
    legal = true;
  }
  else {
    legal = false;
  }
  return legal;
}

//----------------------------------------------------------------------------
Boolean RnShape::u_indexOK_0u(const int dir) const
//----------------------------------------------------------------------------
{

   Boolean ret = (dir >= 0) && (dir < nsd);

   return ret;
}


//-----------------------------------------------------------------------------
Boolean RnShape::u_inShape_0u(const real v) const
//-----------------------------------------------------------------------------
{
  return true;
}

//-----------------------------------------------------------------------------
Boolean RnShape::u_inShape_0u(const real v[MAXDIRECTIONSIZE]) const
//-----------------------------------------------------------------------------
{
  return true;
}




// Generators

//-----------------------------------------------------------------------------
void RnShape::u_setShape_1u(const int n)
//-----------------------------------------------------------------------------
{
  nsd = n;
}

// Composite

//----------------------------------------------------------------------------
void RnShape::u_concatenate_1u(const RnShape& tsm)
//----------------------------------------------------------------------------
{

  nsd = nsd + tsm.nsd;
}

//----------------------------------------------------------------------------
void RnShape::u_extract_1u(const int& l, const int& r)
//----------------------------------------------------------------------------
{

  nsd = r - l;
}


// Indexable
//----------------------------------------------------------------------------
real RnShape::u_getElement_0u(int dir) const
//----------------------------------------------------------------------------
{
  return REAL_MAX;
}

//-----------------------------------------------------------------------------
void RnShape::u_changeElement_1u(int dir, real n_i)
//-----------------------------------------------------------------------------
{
} 
 


// Conventions

//---------------------------------------------------------------------------
void RnShape::operator= (const RnShape& ts)
//---------------------------------------------------------------------------
{
    nsd = ts.nsd;
}

//----------------------------------------------------------------------------
Boolean RnShape::operator== (const RnShape& ts) const
//----------------------------------------------------------------------------
{
  Boolean correct = true;
  if (nsd != ts.nsd)
  { correct = false;
  }
  return correct;
}

//----------------------------------------------------------------------------
RnShape::~RnShape()
//----------------------------------------------------------------------------
{
}






//-----------------------------------------------------------------------------
void RnShape::u_copy_1u() // copy0
//-----------------------------------------------------------------------------
{

  nsd = 0;
  
}

//-----------------------------------------------------------------------------
void RnShape::u_copy_1u(const int n) //copy2
//-----------------------------------------------------------------------------
{
  nsd = n;
}




/*---------------------------------------------------*/
/*              OrderRelations                       */
/*---------------------------------------------------*/
Boolean RnShape::operator<=( const RnShape & ts ) const
{
  return nsd <= ts.nsd;
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
