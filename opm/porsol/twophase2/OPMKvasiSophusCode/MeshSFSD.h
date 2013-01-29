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


#ifndef MeshSFSD_H
#define MeshSFSD_H

//---------------------- Include Files ----------------------------------------
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <opm/porsol/twophase2/OPMKvasiSophusCode/GlobConst.h>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/GlobType.h>

#include <opm/porsol/twophase2/OPMKvasiSophusCode/MeshShape.h>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/MeshPoint.h>


//---------------------  Constants --------------------------------------------


//---------------------- Types ------------------------------------------------





#ifdef DEBUG_MeshSFSD_h
#define DEBUG_NOW
#endif


extern "C"
void matinv_(double *matrix, int *size) ;

extern "C"
void eqsolve_(double *matrix, double *rightHandSides, int *size, int *noOfRightHandSides) ;

extern "C"
void eqsolvb_(double *matrix, double *rightHandSides, int *size, int *noOfRightHandSides, int *numbSubDiag, int *numbSuperDiag, int *bandStorageSize);

class MeshSFSD
{
  public :

    //**********************************************************************************
    //NBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNB
    //NBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNB
    //Har lagt inn muligheter for bånd matriser i datastrukturen. IMIDLERTID er IKKE
    //alle operasjoner oppgradert med tanke på dette ennå !!!
    //SPESIELT gjelder dette rutinene:
    //u_mMatMult_ 2u, u_mVecMult_2u, u_getData_0u og u_mdivide_1u samt ul_Allocate
    //NBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNB
    //NBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNB
    //**********************************************************************************

    Boolean u_DI_0u() const;
    MeshSFSD() ;
    MeshSFSD(const MeshSFSD &V) ;
    MeshSFSD(const MeshShape &S) ;
    MeshSFSD(const MeshShape &S, const double &v) ;
    void u_copy_1u() ;
    void u_copy_1u(const MeshShape &S) ;
    void u_copy_1u(const MeshShape &S, const double &v) ;
    MeshShape u_getShape_0u() const ;
    int u_getNumberOfRows_0u() const ;
    int u_getNumberOfColumns_0u() const ;
    Boolean u_sameShape_0u(const MeshSFSD &CP) const ;
    Boolean u_legalPoint_0u(const MeshPoint &P) const ;
    double u_getData_0u(const MeshPoint &P) const ;
    void u_print_0u() const ;
    int u_domainDimension_0u() const ;
    double u_getElementFromMatMultOfNonQuadraticMatrices_0u(const int& indRowObj, const MeshSFSD& MatB, const int& indColB) const ;

    void u_setupForNonQuadraticMatrices_1u(const int& numbOfRows, const int& numbOfColumns, const double &v);

    void u_setupForBandedMatrices_1u(const MeshShape &S, const double &v, const int& numbSubDiag, const int& numbSuperDiag);

    void u_setData_1u(const MeshPoint &P, const double &v) ;
    void u_mMatMult_2u(MeshSFSD &M) const ;
    void u_mVecMult_2u(MeshSFSD &V) const ;
    void u_mdivide_1u(const MeshSFSD &V) ;
    void u_mdivideBanded_1u(MeshSFSD &V) ;
    void u_invert_1u() ;

    /* Vector operations */
    void operator+=(const MeshSFSD& CP);
    void operator*=(const double& v);

    Boolean operator==(const MeshSFSD &CP) const ;
    void operator=(const MeshSFSD &CP) ;
    ~MeshSFSD();


  private : 

    double *ul_GetItemPtr(double *p, const int r, const int c) const ;
  double *ul_GetItemPtrGeneral(double *p, const int r, const int c, const int StepLength) const ;
    double *ul_GetItemPtrBanded(double *p, const int r, const int c) const ;
    void ul_Allocate() ;
    void ul_Preset(const double initVal) ;


    MeshShape m_shape;
    int m_domainDimension;
    int m_size;
    int m_nRows;
    int m_nCols;
    double *m_p;
    //The rest of the data structure pertains to banded MATRICES ONLY !!!
    bool m_isMatBanded;
    int m_numbSubDiag;
    int m_numbSuperDiag;

//  DECLARE_SOPHUS(1.11)
};



#endif //MeshSFSD_H



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







