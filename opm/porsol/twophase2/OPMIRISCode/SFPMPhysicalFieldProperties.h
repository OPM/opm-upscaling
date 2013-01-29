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



#ifndef SFPMPhysicalFieldProperties_H
#define SFPMPhysicalFieldProperties_H

/*
#ifdef DEBUG_SFPMPhysicalFieldProperties_h

#define DEBUG_NOW

#endif
*/
//---------------------- Include Files ----------------------------------------

#include <opm/porsol/twophase2/OPMIRISCode/IRISDuneGridInterface.hpp>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/GlobType.h>
#include <opm/porsol/twophase2/OPMIRISCode/IrisOpmTdefs.h>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/Array.h>

#include <opm/porsol/twophase2/OPMIRISCode/PMTwoPhaseSimulationData.h>

//---------------------  Constants --------------------------------------------


//---------------------- Types ------------------------------------------------


class SFPMPhysicalFieldProperties
{
  public : 
  SFPMPhysicalFieldProperties();

  void setup(const IRISDuneGridInterface<DuneGridType>* IGIPtr, const PMTwoPhaseSimulationData* simulationDataPtr);//@HAF: OK???

  const PMTwoPhaseSimulationData* getSimulationDataPtr() const ;

  const Array<int>& getRockType() const ;
  const Array<double>& getPermeability(const int tensorComp) const ;
  const Array<double>& getPorosity() const ;
  double getViscosity(const int phase) const ;
  double getDensityGravity(const int phase) const ;

  void operator=(const SFPMPhysicalFieldProperties& PMpp);

  ~SFPMPhysicalFieldProperties() ;


  private : 

    void ml_fillInForPermeabilityRockTypeAndPorosity();
  void ml_generatePhysicalFieldProperties_DisplacementProb_StructGrid();
  void ml_readPhysicalFieldPropertiesForTriangularGridFromFile();

  const IRISDuneGridInterface<DuneGridType>* IGIPtr_;
  const PMTwoPhaseSimulationData* simulationDataPtr_;

  Array<int> rockType_;
  Array<double> Kxx_, Kyy_, Kzz_, Kxy_, Kxz_, Kyz_;
  Array<double> porosity_;
  /*
    Note that:
    phase=0 : oil
    phase=1 : water
    phase=3 : gas
  */
  double viscosity_[3];
  double density_gravity_[3]; //NB: This is meant to be phase density mulitiplied with accel. of gravity.
  //****************************************************************************************
};


#endif //SFPMPhysicalFieldProperties_H



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

