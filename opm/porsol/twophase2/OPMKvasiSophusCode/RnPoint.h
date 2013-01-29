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


#ifndef RnPointH
#define RnPointH


//---------------------- Include Files ----------------------------------------
#include <opm/porsol/twophase2/OPMKvasiSophusCode/GlobConst.h>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/RnShape.h>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/Array.h>


#ifdef DEBUG_RnPoint_h
#define DEBUG_NOW
#endif


//----------------------------------------------------------------------------
class RnPoint
//----------------------------------------------------------------------------
{
  public:
     // Constructors
     RnPoint(); // RnPoint0
     RnPoint(const RnPoint& tp); // RnPoint1
     RnPoint(const RnShape& ts); // RnPoint2
  
     // Observers
     int u_domainDimension_0u() const;
     RnShape u_getShape_0u() const; 
     Boolean u_sameShape_0u(const RnPoint& Q) const;
     Boolean u_legalDirection_0u(const int dir) const;
     Boolean u_indexOK_0u(const int ind) const;
     real u_getElement_0u(int dir) const;     
     real u_getComponent_0u(const int dir) const;
     // real[] u_getComponents_0u() const;
     Boolean u_isOrigin_0u() const;

     // Generators
     void u_setShape_1u(const RnShape& S);
     void u_setPoint_1u(const RnShape& S,const real v[MAXDIRECTIONSIZE]);
     void u_setPoint_1u(const RnShape& S,const real v);
     void u_changeElement_1u(int dir, real index);
     void u_setComponentsAll_1u(const real val);
     void u_setComponent_1u(const int dir, real val);
     void u_setComponents_1u(const real v[MAXDIRECTIONSIZE]);

     // Vector arithmetic
     void operator+= (const RnPoint& mp); 
     void operator-= (const RnPoint& mp); 
     void operator*= (const real v); 
     void operator*= (const Array<int>& a);
     void u_zero_1u(); 
     void u_one_1u(); 


     // Vector and ring arithmetic
     real u_norm2_0u() const;
     void normalize();    // Normalize the vector
     real u_dot_0u(const RnPoint& mp) const;
     void cross(const RnPoint& mp);
     
     // Rotation
     //**********************************************************
     //NOTE: This function may only be suitable to use in 2D.
     //**********************************************************
     bool isCounterClockwise(const RnPoint& x1, const RnPoint& x2) const;

#ifdef COMPLETE_SOPHUS
//---------------------- Removed from -----------------------------------------
     real u_angle_0u(const RnPoint& mp) const;
//---------------------- Removed to -------------------------------------------
#endif //COMPLETE_SOPHUS

     void operator*= (const RnPoint& mp);
     void operator/= (const RnPoint& mp);  
     
     // Embedding
     void u_embed_1u(const RnPoint& Q);
     void u_embed_1u(const RnShape& S);

     // Scaling
     void u_scale_1u(const RnShape& S);

     // Mirroring
     void u_mirror_1u(const RnShape& S);
     void u_mirror_1u(const RnPoint& Q, const RnShape& S);

     // String Operations
     void u_concatenate_1u(const RnPoint& Q);
     void u_extract_1u(const int & l, const int & r);

     // Conventions
     void operator= (const RnPoint& tp);   
     ~RnPoint();
     Boolean operator== (const RnPoint& rp) const;
     Boolean u_DI_0u() const;
     void u_copy_1u(); // copy0
     void u_copy_1u(const RnShape& ts); // copy2

  private:

     // data members
     RnShape m_rns; // the underlying RnShape for the RnPoint
     int m_nsd; // number of space dimensions
     real m_coor[MAXDIRECTIONSIZE]; // the coordinates of the RnPoint

//   DECLARE_SOPHUS(2.8)

}; //RnPoint







#endif //RnPointH

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












