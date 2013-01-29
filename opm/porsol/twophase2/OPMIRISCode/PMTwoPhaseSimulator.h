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


#ifndef PMTwoPhaseSimulatorH
#define PMTwoPhaseSimulatorH



//*****************************************************************************************

//---------------------- Include Files ----------------------------------------
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

#include <opm/porsol/twophase2/OPMIRISCode/IRISDuneGridInterface.hpp>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/GlobType.h>
#include <opm/porsol/twophase2/OPMIRISCode/IrisOpmTdefs.h>

#include <opm/porsol/twophase2/OPMIRISCode/VTKHeaderWriter.h>
#include <opm/porsol/twophase2/OPMIRISCode/VTKStructuredGridWriter.h>
#include <opm/porsol/twophase2/OPMIRISCode/SFRelPerm.h>
#include <opm/porsol/twophase2/OPMIRISCode/SFPMPhysicalFieldProperties.h>
#include <opm/porsol/twophase2/OPMIRISCode/SFCapPress.h>
#include <opm/porsol/twophase2/OPMIRISCode/PMTwoPhaseSimulationData.h>
#include <opm/porsol/twophase2/OPMIRISCode/PMTwoPhaseElliptic.h>
#include <opm/porsol/twophase2/OPMIRISCode/PMTwoPhaseSaturation.h>


//---------------------  Constants --------------------------------------------




//---------------------- typedef's and enumerations ----------------------------

//---------------------- Directives -------------------------------------------
using namespace std;  // This leads to that all names in the namespace std may be 
                      // used without the qualifier "std::".


//*********************************************************************************************



template <class FVMImplicit, class FVMExplicit> 
class PMTwoPhaseSimulator
{
 public:
  //--------------------------------------------------------------------------
  //---------------------- Constructors --------------------------------------
  //--------------------------------------------------------------------------
  PMTwoPhaseSimulator(const int& numbElements);
  PMTwoPhaseSimulator(const FVMImplicit& FVMIMeth, const FVMExplicit& FVMEMeth, const int& numbElements);

  //--------------------------------------------------------------------------
  //---------------------- Generators-----------------------------------------
  //--------------------------------------------------------------------------
  void setup(const IRISDuneGridInterface<DuneGridType>& Igrid, const PMTwoPhaseSimulationData& simulationData, const SFPMPhysicalFieldProperties& PMpp);
  void performSimulation(const IRISDuneGridInterface<DuneGridType>& Igrid, const PMTwoPhaseSimulationData& simulationData, const TimeIntegrationType& integrationType);

  //--------------------------------------------------------------------------
  //---------------------- Observers -----------------------------------------
  //--------------------------------------------------------------------------


  //@HAF: Trenger vi en destructor???


  //-----------------------------------------------------------------------------------

  private : 
  //-----------------Useful local functions--------------------------------------------
  void ml_initializeSaturation(const IRISDuneGridInterface<DuneGridType>& Igrid, const PMTwoPhaseSimulationData& simulationData);
  void ml_generateInitializeSaturationStructuredGridPaperEx2(const IRISDuneGridInterface<DuneGridType>& Igrid, const PMTwoPhaseSimulationData& simulationData);
  void ml_generateInitializeSaturationUnstructuredGridPaperEx3(const IRISDuneGridInterface<DuneGridType>& Igrid, const PMTwoPhaseSimulationData& simulationData);
  void ml_writePressureToFileInVTKFormat(const IRISDuneGridInterface<DuneGridType>& Igrid, const PMTwoPhaseSimulationData& simulationData, const int& numbElements, const int& fileID);
  void ml_writeSaturationToFileInVTKFormat(const IRISDuneGridInterface<DuneGridType>& Igrid, const PMTwoPhaseSimulationData& simulationData, const int& fileID);


  //-----------------Data structure----------------------------------------------------

  double deltat_;//Time-step in IMPES algorithm.
  SFSaturation saturation_;
  PMTwoPhaseElliptic<FVMImplicit> pressureEq_;
  PMTwoPhaseSaturation<FVMExplicit> saturationEq_;
  //@HAF: What about visualization stuff etc.????????????????????????????
};




/*---------------------------------------------------*/
/*              Constructors                         */
/*---------------------------------------------------*/

template <class FVMImplicit, class FVMExplicit>
PMTwoPhaseSimulator<FVMImplicit, FVMExplicit>::PMTwoPhaseSimulator(const int& numbElements)
  : deltat_(0.0), saturation_(), pressureEq_(numbElements), saturationEq_()
{

}


template <class FVMImplicit, class FVMExplicit>
PMTwoPhaseSimulator<FVMImplicit, FVMExplicit>::PMTwoPhaseSimulator(const FVMImplicit& FVMIMeth, const FVMExplicit& FVMEMeth, const int& numbElements)
  : deltat_(0.0), saturation_(), pressureEq_(FVMIMeth, numbElements), saturationEq_(FVMEMeth)
{

}


/*---------------------------------------------------*/
/*              Generators                           */
/*---------------------------------------------------*/
template <class FVMImplicit, class FVMExplicit>
void PMTwoPhaseSimulator<FVMImplicit, FVMExplicit>::setup(const IRISDuneGridInterface<DuneGridType>& Igrid, const PMTwoPhaseSimulationData& simulationData, const SFPMPhysicalFieldProperties& PMpp)
{
  deltat_ = simulationData.IMPESd_.delta_T;
  ml_initializeSaturation(Igrid, simulationData);

  pressureEq_.setup(Igrid, simulationData, PMpp);
  saturationEq_.setup(Igrid, saturation_, simulationData, PMpp);
}


template <class FVMImplicit,class FVMExplicit>
void PMTwoPhaseSimulator<FVMImplicit,FVMExplicit>::performSimulation(const IRISDuneGridInterface<DuneGridType>& Igrid, const PMTwoPhaseSimulationData& simulationData, const TimeIntegrationType& integrationType)
{
  //The IMPES algorithm (fractional flow formulation of Hoteit & Firoozabadi (2008)) is implemented below.
  double actualT = simulationData.IMPESd_.T_Start;
  double TOL = 0.000000001;//@HAF: Hm...

  int step = 0;
  while (actualT < simulationData.IMPESd_.T_End - TOL)
  {
    pressureEq_.setupTwoPhasePressureEq(Igrid, saturation_);//@HAF: If this should be OK, BC's etc. must already have beem transmitted to pressureEq_ through the setup routine...+ computation of relevant mobilities and Psi are inside this routine!!!
    pressureEq_.solveLinearEquationSystemForPressure();//@HAF: In this case it is NOT precisely the pressure, but nevertheless...
    pressureEq_.computeVelocity(Igrid, saturation_);//@HAF: In this case it is actually NOT the total velocity, but the so called v_a

    saturationEq_.computeForwardInTime(deltat_, integrationType, pressureEq_.getVelocity());
    saturation_ = saturationEq_.getSaturation();
    //***************DEBUGGING***START****************************
    cout << "Saturation= " << endl;
    cout << saturation_[0] << endl;
    cout << saturation_[1] << endl;
    cout << saturation_[2] << endl;
    cout << saturation_[3] << endl;
    cout << "***********************************************" << endl;
    //***************DEBUGGING***END******************************
    //exit(0);
    //Post processing of Data:
    if (((step+1) % simulationData.IMPESd_.OutputFrequency) == 0)//@HAF: OK???
    {
      ml_writePressureToFileInVTKFormat(Igrid, simulationData, Igrid.cellCount(0,0), step+1);
      ml_writeSaturationToFileInVTKFormat(Igrid, simulationData, step+1);
    }

    actualT += deltat_;
    step++;
    pressureEq_.reset();//Makes sure that all elements of velocity, internal matrix, right hand side and solution vector are put to zero
  }
}



/*---------------------------------------------------*/
/*             Local functions                       */
/*---------------------------------------------------*/

template <class FVMImplicit,class FVMExplicit>
void PMTwoPhaseSimulator<FVMImplicit,FVMExplicit>::ml_initializeSaturation(const IRISDuneGridInterface<DuneGridType>& Igrid, const PMTwoPhaseSimulationData& simulationData)
{
  IRISDuneGridInterface<DuneGridType>::GridType gridType = Igrid.getGridType();//@HAF: gridType MUST of course be an enum (of type GridType see IRISDuneGridInterface)!!!

  saturation_.u_copy_1u(Igrid.cellCount(0,0), 0.0);

  if (gridType == IRISDuneGridInterface<DuneGridType>::_Triangular)
  {
    ml_generateInitializeSaturationUnstructuredGridPaperEx3(Igrid, simulationData);
  }
  else if (gridType == IRISDuneGridInterface<DuneGridType>::_Structured2D)
  {
    ml_generateInitializeSaturationStructuredGridPaperEx2(Igrid, simulationData);
  }
}


template <class FVMImplicit,class FVMExplicit>
void PMTwoPhaseSimulator<FVMImplicit,FVMExplicit>::ml_generateInitializeSaturationStructuredGridPaperEx2(const IRISDuneGridInterface<DuneGridType>& Igrid, const PMTwoPhaseSimulationData& simulationData)
{
  //***********************************************************************************
  //We generate the initial saturation for the PAPER_EX2 example on a structured grid:
  //***********************************************************************************

  //***********************************************************************************
  //@HAF: NOTE: We think this is correct according to the numbering scheme used in Dune,
  //but this MUST be thoroughly checked in the debuging process!!!
  //***********************************************************************************

  int NX = simulationData.td_.NoGridCells_X;
  int NY = simulationData.td_.NoGridCells_Y;

  for (int j=0; j < NY; j++)
  {
    int k = 0;
    for (int i=0; i < NX; i++)
    {
      k = i + j*NX;
      if (i == 0)
      {
	saturation_.u_changeElement_1u(k, 1.0);
      }
      else
      {
	saturation_.u_changeElement_1u(k, 0.0);
      }
    }  
  }
}


template <class FVMImplicit,class FVMExplicit>
void PMTwoPhaseSimulator<FVMImplicit,FVMExplicit>::ml_generateInitializeSaturationUnstructuredGridPaperEx3(const IRISDuneGridInterface<DuneGridType>& Igrid, const PMTwoPhaseSimulationData& simulationData)
{
  //***********************************************************************************
  //We generate the initial saturation for the PAPER_EX3 example on a unstructured grid:;
  //but, in fact, there is nothing to do in this case :-),
  //This is because the initial saturation field should be simply zero.
  //***********************************************************************************

}


template <class FVMImplicit,class FVMExplicit>
void PMTwoPhaseSimulator<FVMImplicit,FVMExplicit>::ml_writePressureToFileInVTKFormat(const IRISDuneGridInterface<DuneGridType>& Igrid, const PMTwoPhaseSimulationData& simulationData, const int& numbElements, const int& fileID)
{
  const SFPressure& pressure = pressureEq_.getPressure(numbElements);//@HAF: OK eller feil pga. referanse overføring som output fra funksjonen???

  IRISDuneGridInterface<DuneGridType>::GridType gridType = Igrid.getGridType();//@HAF: gridType MUST of course be an enum (of type GridType see IRISDuneGridInterface)!!!
  if (gridType == IRISDuneGridInterface<DuneGridType>::_Triangular)
  {
    //NOTE: NB: Remember to allocate memory for Kxx_ etc.!!!!!!!
    cout << "Nothing implemented yet." << endl; exit(0);
  }
  else if (gridType == IRISDuneGridInterface<DuneGridType>::_Structured2D)
  {
    VTKStructuredGridWriter::SurfaceOutput surface_output;
    surface_output = VTKStructuredGridWriter::planar2D;
    VTKStructuredGridWriter CellwriterP(&Igrid, pressure, &simulationData, 
					surface_output,
					VTKStructuredGridWriter::on_cells,
					"Pressure output");
  
    if (!CellwriterP.inputVerifiedOK()) {
      cerr << "Input grid not OK, aborting.";
    };

    char vtkFileNameP[100];
    sprintf(vtkFileNameP, "PressureField_%d.vtk", fileID);
  
    ofstream ofileCellsP(vtkFileNameP);
    if (!ofileCellsP) {
      cerr << "writetest.cpp: Cannot open output file (Pressure).  Aborting." << endl;
      exit(0);
    }
  
    CellwriterP.write(ofileCellsP);
  }
}
 
template <class FVMImplicit,class FVMExplicit>
void PMTwoPhaseSimulator<FVMImplicit,FVMExplicit>::ml_writeSaturationToFileInVTKFormat(const IRISDuneGridInterface<DuneGridType>& Igrid, const PMTwoPhaseSimulationData& simulationData, const int& fileID)
{
  IRISDuneGridInterface<DuneGridType>::GridType gridType = Igrid.getGridType();//@HAF: gridType MUST of course be an enum (of type GridType see IRISDuneGridInterface)!!!
  if (gridType == IRISDuneGridInterface<DuneGridType>::_Triangular)
  {
    //NOTE: NB: Remember to allocate memory for Kxx_ etc.!!!!!!!
    cout << "Nothing implemented yet." << endl; exit(0);
  }
  else if (gridType == IRISDuneGridInterface<DuneGridType>::_Structured2D)
  {
    VTKStructuredGridWriter::SurfaceOutput surface_output;
    surface_output = VTKStructuredGridWriter::planar2D;
    VTKStructuredGridWriter CellwriterS(&Igrid, saturation_, &simulationData, 
					surface_output,
					VTKStructuredGridWriter::on_cells,
					"Saturation output");
  
    if (!CellwriterS.inputVerifiedOK()) {
      cerr << "Input grid not OK, aborting.";
    };

    char vtkFileNameS[100];
    sprintf(vtkFileNameS, "SaturationField_%d.vtk", fileID);
  
    ofstream ofileCellsS(vtkFileNameS);
    if (!ofileCellsS) {
      cerr << "writetest.cpp: Cannot open output file (Saturation).  Aborting." << endl;
      exit(0);
    }
    
    CellwriterS.write(ofileCellsS);
  }
}


#endif //PMTwoPhaseSimulatorH

/*
===============================================================================
    ------------------------ CLASS DESCRIPTION --------------------------
===============================================================================

  DESCRIPTION:
  ------------------
  
    


  USAGE/EXAMPLES:
  ---------------
    <Examples how to use the class>
    
    
  MISCELLANEOUS:
  --------------
  
  
  
  


===============================================================================
    --------------------- END OF CLASS DESCRIPTION ----------------------
===============================================================================
*/




