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


#ifndef MeshShapeH
#define MeshShapeH

//---------------------- Include Files ----------------------------------------

#include <opm/porsol/twophase2/OPMKvasiSophusCode/GlobConst.h>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/GlobType.h>


//---------------------  Constants --------------------------------------------


//---------------------- Types ------------------------------------------------





#ifdef DEBUG_MeshShape_h
#define DEBUG_NOW
#endif



class MeshShape
{
public:
  //--------------------
  // SophusStandard
  //--------------------

  // Constructors
  MeshShape();                      // MeshShape0
  MeshShape(const MeshShape& ms);   // MeshShape1
  MeshShape(const int& size);   // MeshShape2
  MeshShape(const int& size, const int& def);   // MeshShape3
 
  void u_copy_1u();                     // copy0
  void u_copy_1u(const int& size);   // copy2
  void u_copy_1u(const int& size, const int& def);   // copy3

  // Observers
  Boolean u_DI_0u() const;
  Boolean operator == (const MeshShape& ms) const; 


  // Conventions
  void operator = (const MeshShape& ms);
  ~MeshShape();


  //--------------------
  // CartShape
  //--------------------

  // Observers
  int u_domainDimension_0u() const;
  Boolean u_sameShape_0u(const MeshShape& msm) const;
  Boolean u_legalDirection_0u(int dir) const;
  Boolean u_indexOK_0u(const int ind) const;
  void u_getElement_0g(const int dir, int& e) const;
  int u_getDirectionSize_0u(int dir) const;
  void u_getDirectionSizes_0g(int cp[MAXDIRECTIONSIZE]) const; 
  Boolean u_inShape_0u(int v) const;
  Boolean u_inShape_0u(const int cp[MAXDIRECTIONSIZE]) const;

  // Generators
  void u_setShape_1u(int n,int v);
  void u_setShape_1u(int n,const int vp[MAXDIRECTIONSIZE]);
  void u_changeElement_1u(const int & dir, const int & e);
  void u_setDirectionSizeAll_1u(int n_i); 
  void u_setDirectionSize_1u(int dir, int n_i); 
  void u_setDirectionSizes_1u(const int vp[MAXDIRECTIONSIZE]); 

  // String operations
  void u_concatenate_1u(const MeshShape& msm);
  void u_extract_1u(const int & l, const int & r); 

  // Volumes
  int u_volume_0u() const;
  int u_getDirectionVolume_0u(int d) const; 
  int u_getMaxDirectionVolume_0u() const; 
  int u_getMinDirectionVolume_0u() const; 

  // SubShapes
  Boolean u_subShape_0u(const MeshShape& msm) const;

  // Scale factors
  Boolean u_isScale_0u() const;
  void u_scale_1u(const MeshShape& SP);

  // Mirrors
  Boolean u_isMirror_0u(const MeshShape& msm) const;
  int u_mirrorShape_0u() const;

  // Borders
  Boolean u_isBorder_0u(const MeshShape& msm) const;
  Boolean u_isEmbeddableBorder_0u(const MeshShape& msm) const;
  int u_borderShape_0u() const;
  void u_getBorder_1u(const MeshShape& msm);


  //--------------------
  // MeshShape
  //--------------------

  // Observers
  int u_getSize_0u() const;

  // Generators
  void u_setShape_1u(int nsd);

  // Scale factors
  Boolean u_isScalable_0u(const MeshShape& SP) const;
  void u_getScaleFactor_1u(const MeshShape& SP); 

private:  

  int nsd;   // number of space dimensions
  int v[MAXDIRECTIONSIZE]; // the number of nodes in each dimension 


//  DECLARE_SOPHUS(2.8)

}; //MeshShape - Class








#endif //MeshShapeH



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
  
  NOTE: MAXDIRECTIONSIZE must be defined before including this file
  --- End ---
  
  
  


===============================================================================
    --------------------- END OF CLASS DESCRIPTION ----------------------
===============================================================================
*/

