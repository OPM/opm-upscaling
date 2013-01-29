/*
  Copyright 2011 IRIS - International Research Institute of Stavanger.

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



#ifndef SFMpfaFpsH
#define SFMpfaFpsH



//*******************************************************************************

//---------------------- Include Files ----------------------------------------
#include <opm/porsol/twophase2/OPMIRISCode/IRISDuneGridInterface.hpp>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/GlobType.h>
#include <opm/porsol/twophase2/OPMIRISCode/IrisOpmTdefs.h>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/Array.h>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/Pair.h>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/RnShape.h>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/RnPoint.h>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/MeshSFSD.h>
#include <opm/porsol/twophase2/OPMIRISCode/SFPMPhysicalFieldProperties.h>
#include <opm/porsol/twophase2/OPMIRISCode/SFRelPerm.h>
#include <opm/porsol/twophase2/OPMIRISCode/SFCapPress.h>

//---------------------  Constants --------------------------------------------




//---------------------- typedef's and enumerations ----------------------------


//---------------------- Directives -------------------------------------------
using namespace std;  // This leads to that all names in the namespace std may be 
                      // used without the qualifier "std::".


//*******************************************************************************



//----------------------------------------------------------------------------
class SFMpfaFps
//----------------------------------------------------------------------------
{
  public:
     // Constructor
     SFMpfaFps();

     // Copy constructor
     SFMpfaFps(const SFMpfaFps& Fps);

     // Destructor
     ~SFMpfaFps();//@HAF: Do we need this???

     // Generators
     double computeOnePhaseTransmissibility(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const SFPMPhysicalFieldProperties& PMpp, const std::vector<int>& outRelk, const int& itilde, const IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt, const int& ihat, const int& ihatTilde, const IRISDuneGridInterface<DuneGridType>::EdgePointer khatPt);//@HAF: OK???

     double computeOnePhaseBoundaryTransmissibility(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const SFPMPhysicalFieldProperties& PMpp, const std::vector<int>& outRelk, const IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt, const int& ihat, const IRISDuneGridInterface<DuneGridType>::EdgePointer khatPt);

     void generateLocalInverseMatrixAndRHS(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const SFPMPhysicalFieldProperties& PMpp, const int& pressureProblemType);//@HAF: OK???

     void generateLocalBoundaryInverseMatrixAndRHS(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const SFPMPhysicalFieldProperties& PMpp, const SFCapPress& cPress, const int& pressureProblemType);

  private:

     // local functions
     double ml_matMultMH(const int& kBar, const int& itilde) const;

     void ml_generateLocalInverseMatrixAndRHS(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const Array<double>& K11, const Array<double>& K22, const Array<double>& K12, const Array< double >& alpha1, const Array< double >& alpha2);
     void ml_generateLocalBoundaryInverseMatrixAndRHS(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const Array<double>& K11, const Array<double>& K22, const Array<double>& K12, const Array<double>& alpha1, const Array<double>& alpha2, const Array<double>& g);

     double ml_FE1(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt, const int &ihat, const int& j, const int& normalVector, const Array<double>& K11, const Array<double>& K22, const Array<double>& K12);
     double ml_FE1B(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt, const int &ihat, const int& j, const IRISDuneGridInterface<DuneGridType>::EdgePointer khatPt, const int& normalVector, const Array<double>& K11, const Array<double>& K22, const Array<double>& K12, const Array< double > &alpha1, const Array< double > &alpha2);
     double ml_FE2(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt, const int &ihat, const int& normalVector, const Array<double>& K11, const Array<double>& K22, const Array<double>& K12);
     double ml_FE2B(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt, const int &ihat, const IRISDuneGridInterface<DuneGridType>::EdgePointer khatPt, const int& normalVector, const Array<double>& K11, const Array<double>& K22, const Array<double>& K12, const Array< double > &alpha1, const Array< double > &alpha2);
     double ml_FE3(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt, const int &ihat, const int& j, const int& normalVector, const Array<double>& K11, const Array<double>& K22, const Array<double>& K12);
     double ml_FE4(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt, const int &ihat, const int& normalVector, const Array<double>& K11, const Array<double>& K22, const Array<double>& K12);

     void ml_createConstBoundaryConditionData(const int& BCConstNR, const IRISDuneGridInterface<DuneGridType>& Igrid,  const Array<double>& valuesExternalSurfaces,  Array<double>& BCData) const;
     void ml_changeBoundaryConditionDataFor_g(const SFPMPhysicalFieldProperties& PMpp, const SFCapPress& cPress, const IRISDuneGridInterface<DuneGridType>& Igrid, Array<double>& BCData) const;
     int ml_findRockTypeAtBoundary(const SFPMPhysicalFieldProperties& PMpp, const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::FacePointer khatPt) const;
     int ml_getTransformSpaceNormalVectorLabel(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt, const IRISDuneGridInterface<DuneGridType>::EdgePointer khatPt) const;
     double ml_computeTransformSpaceIntegralOfTensorComponents(const IRISDuneGridInterface<DuneGridType>& Igrid, const bool& BC, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt, const int& ihat, const Array<double>& K11, const Array<double>& K22, const Array<double>& K12, const int& integralType) const;

     double ml_x_xi_Forxi_ANDeta_eq1_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_x_xi_Forxi_ANDeta_eq1_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_x_xi_Forxi_eq1ANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_x_xi_Forxi_eq1ANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     
     double ml_y_xi_Forxi_ANDeta_eq1_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_y_xi_Forxi_ANDeta_eq1_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_y_xi_Forxi_eq1ANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_y_xi_Forxi_eq1ANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     
     double ml_x_eta_Forxi_ANDeta_eq1_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_x_eta_Forxi_ANDeta_eq1_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_x_eta_Forxi_eq1ANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_x_eta_Forxi_eq1ANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;

     double ml_y_eta_Forxi_ANDeta_eq1_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_y_eta_Forxi_ANDeta_eq1_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_y_eta_Forxi_eq1ANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_y_eta_Forxi_eq1ANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     
     double ml_determinant_Forxi_ANDeta_eq1_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_determinant_Forxi_ANDeta_eq1_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_determinant_Forxi_eq1ANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_determinant_Forxi_eq1ANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     
     //NOTE: The following is only needed for bilinear expansions:

     double ml_x_xi_Forxi_ANDeta_eq0_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_x_xi_Forxi_ANDeta_eq0_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_x_xi_Forxi_eq0ANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_x_xi_Forxi_eq0ANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     
     double ml_y_xi_Forxi_ANDeta_eq0_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_y_xi_Forxi_ANDeta_eq0_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_y_xi_Forxi_eq0ANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_y_xi_Forxi_eq0ANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     
     double ml_x_eta_Forxi_ANDeta_eq0_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_x_eta_Forxi_ANDeta_eq0_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_x_eta_Forxi_eq0ANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_x_eta_Forxi_eq0ANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     
     double ml_y_eta_Forxi_ANDeta_eq0_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_y_eta_Forxi_ANDeta_eq0_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_y_eta_Forxi_eq0ANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_y_eta_Forxi_eq0ANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     
     double ml_determinant_Forxi_ANDeta_eq0_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_determinant_Forxi_ANDeta_eq0_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_determinant_Forxi_eq0ANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_determinant_Forxi_eq0ANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     
     //NOTE: The following is only needed for bilinear expansions in the case of 
     //auxillary dual-cells

     double ml_x_xi_Forxi_ANDeta_eqPARAM_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4, const double& eta) const ;
     double ml_x_xi_Forxi_ANDeta_eqPARAM_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_x_xi_Forxi_eqPARAMANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_x_xi_Forxi_eqPARAMANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     
     double ml_y_xi_Forxi_ANDeta_eqPARAM_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4, const double& eta) const ;
     double ml_y_xi_Forxi_ANDeta_eqPARAM_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_y_xi_Forxi_eqPARAMANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_y_xi_Forxi_eqPARAMANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     
     double ml_x_eta_Forxi_ANDeta_eqPARAM_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_x_eta_Forxi_ANDeta_eqPARAM_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_x_eta_Forxi_eqPARAMANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4, const double& xi) const ;
     double ml_x_eta_Forxi_eqPARAMANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     
     double ml_y_eta_Forxi_ANDeta_eqPARAM_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_y_eta_Forxi_ANDeta_eqPARAM_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     double ml_y_eta_Forxi_eqPARAMANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4, const double& xi) const ;
     double ml_y_eta_Forxi_eqPARAMANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const ;
     
     double ml_determinant_Forxi_ANDeta_eqPARAM_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4, const double& eta) const ;
     double ml_determinant_Forxi_ANDeta_eqPARAM_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4, const double& eta) const ;
     double ml_determinant_Forxi_eqPARAMANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4, const double& xi) const ;
     double ml_determinant_Forxi_eqPARAMANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4, const double& xi) const ;
     

     // data members
     MeshSFSD invMmatrix_, hmatrix_;

}; //SFMpfaFps


#endif //SFMpfaFpsH

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
