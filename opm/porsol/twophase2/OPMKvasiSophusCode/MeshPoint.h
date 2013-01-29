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



#ifndef MeshPointH
#define MeshPointH

//---------------------- Include Files ----------------------------------------

#include <opm/porsol/twophase2/OPMKvasiSophusCode/GlobConst.h>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/GlobType.h>

#include <opm/porsol/twophase2/OPMKvasiSophusCode/MeshShape.h>

//---------------------  Constants --------------------------------------------


//---------------------- Types ------------------------------------------------







#ifdef DEBUG_MeshPoint_h
#define DEBUG_NOW
#endif



//----------------------------------------------------------------------------
class MeshPoint
//----------------------------------------------------------------------------
{
  public:
    //--------------------
    // SophusStandard
    //--------------------

    // Constructors
    MeshPoint();                      // MeshPoint0
    MeshPoint(const MeshPoint& mp);   // MeshPoint1
    MeshPoint(const MeshShape& ms);   // MeshPoint2
    MeshPoint(const MeshShape& ms, const int& def);   // MeshPoint3

    void u_copy_1u();                 // copy0
    void u_copy_1u(const MeshShape& ms);   // copy2
    void u_copy_1u(const MeshShape& ms, const int& def);   // copy3

    // Observers
    Boolean u_DI_0u() const;
    Boolean operator == (const MeshPoint& ms) const; 
    int operator[] (const int& dir) const; 



    // Conventions
    void operator = (const MeshPoint& ms);
    ~MeshPoint();

    //--------------------
    // Indexable
    //--------------------

    int u_getSize_0u() const;
    void u_getElement_0g(const int dir, int& e) const;
    void u_changeElement_1u(int dir, int e);

    //--------------------
    // CartPoint
    //--------------------

    // Observers

    int u_domainDimension_0u() const;
    MeshShape u_getShape_0u() const; 
    Boolean u_sameShape_0u(const MeshPoint& Q) const;
    Boolean u_legalDirection_0u(int dir) const;
    Boolean u_indexOK_0u(const int ind) const;
    int u_getComponent_0u(int dir) const;
    void u_getComponents_0g(int cp[MAXDIRECTIONSIZE]) const;

    // Generators
    void u_setPoint_1u(const MeshShape& S, const int coords[MAXDIRECTIONSIZE]);
    void u_setPoint_1u(const MeshShape& S, int coords);
    void u_setComponentsAll_1u(int val);
    void u_setComponent_1u(int dir, int index);
    void u_setComponents_1u(const int v[MAXDIRECTIONSIZE]);

    // Embeddings and borders
    void u_embed_1u(const MeshShape& S);
    void u_borderEmbed_1u(const MeshShape& S,const MeshShape& SP);

    // String operations
    void u_concatenate_1u(const MeshPoint& Q);
    void u_extract_1u(const int & l, const int & r); 

    //--------------------
    // CartPointGroup
    //--------------------

    // Observers
    Boolean u_isOrigin_0u() const;

    // Generators
    void u_setShape_1u(const MeshShape& S);

    // Vector arithmetic
    void operator += (const MeshPoint& mp); 
    void operator -= (const MeshPoint& mp); 
    void operator *= (int v); 
    void u_zero_1u(); 

    // Embedding
    void u_embed_1u(const MeshPoint& Q);

    // Scaling
    void u_scale_1u(const MeshShape& S);

    // Mirroring
    void u_mirror_1u(const MeshShape& S);
    void u_mirror_1u(const MeshPoint& Q, const MeshShape& S);


    //--------------------
    // MeshPoint
    //--------------------

    // Observers
    Boolean u_isDivisible_0u(const MeshPoint& mp) const;

    // Composite
    int u_getLexicographic_0u() const;
    void u_setLexicographic_1u(int ind);

    // Vector and ring arithmetic
    void operator *= (const MeshPoint& mp);
    void operator /= (const MeshPoint& mp);  



  private:

    void ul_gcd( int x, int n, int& g, int& c, int& d) const;
    int ul_isDiv(int a, int b, int limit) const;

    // data members
    int nsd; // number of space dimensions
    int ind[MAXDIRECTIONSIZE]; // the coordinates of the MeshPoint
    MeshShape shape; // the underlying MeshShape for the MeshPoint


    //template<class codeboost> void rules();

//DECLARE_SOPHUS(2.9)


}; //MeshPoint - Class






#endif //MeshPointH



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

