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


#ifndef SFCentralUpwindH
#define SFCentralUpwindH



//*******************************************************************************

//---------------------- Include Files ----------------------------------------
#include <opm/porsol/twophase2/OPMIRISCode/IRISDuneGridInterface.hpp>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/GlobType.h>
#include <opm/porsol/twophase2/OPMIRISCode/IrisOpmTdefs.h>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/Array.h>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/RnShape.h>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/RnPoint.h>
#include <opm/porsol/twophase2/OPMIRISCode/SFRelPerm.h>
#include <opm/porsol/twophase2/OPMIRISCode/SFPMPhysicalFieldProperties.h>


//---------------------  Constants --------------------------------------------




//---------------------- typedef's and enumerations ----------------------------


//---------------------- Directives -------------------------------------------
using namespace std;  // This leads to that all names in the namespace std may be 
                      // used without the qualifier "std::".


//*******************************************************************************



//----------------------------------------------------------------------------
class SFCentralUpwind
//----------------------------------------------------------------------------
{
  public:
     // Constructor
     SFCentralUpwind();

     // Copy constructor
     SFCentralUpwind(const SFCentralUpwind& SFC);

     // Destructor
     ~SFCentralUpwind();

     // Generators
     void setup(const IRISDuneGridInterface<DuneGridType>* IgridPtr);
  
     // Observers
     double getMinGridSize() const;
     double computeTimeStep(const double& minGridSize, const double& maxVelocity) const;
     void computeFluxOutOfAllCells(const SFRelPerm& relPerm, const SFPMPhysicalFieldProperties& PMPhysProp, SFSaturation& saturation, const SFVelocity& velocity_a, Array<double>& rhs);//@HAF: OK med ref. paa SFSaturation??? 

  private:

     // local functions
     double ml_computeVelocity_Kurganov(const double& sPluss, const double& sMinus, const int& RTPluss, const int& RTMinus, const double& totalVelocityN, const SFRelPerm& relPerm, const SFPMPhysicalFieldProperties& PMPhysProp);
     void ml_computeSaturationGradientArminjon(const SFSaturation& saturation);
     void ml_MINMODLimitingOfSaturationGradients();
     void ml_generateLimitersForStructuredGrid(const SFSaturation& saturation, const int& NX, const int& NY, const double& DX, const double& DY);
     double ml_MINMODLimiterStructuredGrid(const SFSaturation& saturation, const int& ID, const int& IDdirP1, const int& IDdirM1) const;
     double ml_MINMOD_2ARG(const double& arg1, const double& arg2);
     double ml_MINMOD_3ARG(const double& arg1, const double& arg2, const double& arg3);
     double ml_MINMOD_4ARG(const double& arg1, const double& arg2, const double& arg3, const double& arg4);

     int ml_getNonPeriodicBCType(const SFPMPhysicalFieldProperties& PMPhysProp, const IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt, const int& k);
     double ml_getBCDirichletCond(const SFPMPhysicalFieldProperties& PMPhysProp, const IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt, const int& k);

     // data members
     const IRISDuneGridInterface<DuneGridType>* IGIPtr_; // the grid
     Array<RnPoint> saturationGrad_, saturationGradNew_;//Important saturation gradients for second order scheme.
     Array<bool> fluxIsComputed_;//Help-Dune::array in order to avoid dublicate flux evaluation at edges.


}; //SFCentralUpwind


#endif //SFCentralUpwindH

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
