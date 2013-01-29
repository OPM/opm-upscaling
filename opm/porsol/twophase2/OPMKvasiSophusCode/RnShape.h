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



#ifndef RnShapeH
#define RnShapeH

//---------------------- Include Files ----------------------------------------

#include <opm/porsol/twophase2/OPMKvasiSophusCode/GlobConst.h>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/GlobType.h>



//---------------------  Constants --------------------------------------------


//---------------------- Types ------------------------------------------------









#ifdef DEBUG_RnShape_h
#define DEBUG_NOW
#endif


//------------------------------------------------------------------------------
class RnShape
//------------------------------------------------------------------------------
{
public:
  Boolean u_DI_0u() const;

  // Constructors
  RnShape(); // RnShape0
  RnShape(const RnShape & S); // RnShape1
  RnShape(const int n); // RnShape2
 
  // Observers
  int u_domainDimension_0u() const;
  Boolean u_sameShape_0u(const RnShape& tsm) const;
  Boolean u_legalDirection_0u(const int dir) const;
  Boolean u_indexOK_0u(const int ind) const;
  Boolean u_inShape_0u(const real v) const;
  Boolean u_inShape_0u(const real v[MAXDIRECTIONSIZE]) const;

  // Generators
  void u_setShape_1u(const int nsd);

  // Composite
  void u_concatenate_1u(const RnShape& tsm);     
  void u_extract_1u(const int& l, const int& r);

  // Indexable
  real u_getElement_0u(int dir) const;
  void u_changeElement_1u(int dir, real n_i);

  // Conventions
  void operator= (const RnShape& ts);
  ~RnShape();
  Boolean operator== (const RnShape& ts) const; 

  void u_copy_1u(); // copy0
  void u_copy_1u(const int n); // copy2



  /*---------------------------------------------------*/
  /*              OrderRelations                       */
  /*---------------------------------------------------*/
  Boolean operator<=( const RnShape & ts ) const;

private:  

  int nsd;   // number of space dimensions
  /* There are no bounds in the values in each direction, no wraparound*/

  //DECLARE_SOPHUS(1.12)
}; //RnShape



//#pragma SophusCode Functional
//#pragma SophusCode __LINE__ __FILE__



#endif //RnShapeH



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
