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


#ifndef PMTwoPhaseSimulationData_H
#define PMTwoPhaseSimulationData_H

//---------------------- Include Files ----------------------------------------

#include <iostream>
//#include <strstream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <opm/porsol/twophase2/OPMKvasiSophusCode/GlobType.h>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/Array.h>

//---------------------- Directives -------------------------------------------
//using namespace std;  // This leads to that all names in the namespace std may be 
                      // used without the qualifier "std::".


//---------------------  Constants --------------------------------------------
#define FILE_INPUT_STR_LEN  100

//---------------------- Types ------------------------------------------------


typedef struct
{
  char NameOfGridType[FILE_INPUT_STR_LEN];
  char NameOfDGFFile[FILE_INPUT_STR_LEN];
  int  DiscretizationType;
  int  UnStructuredGrid;
  int  ReadStructuredGrid;
  int  PeriodicBoundaryConditions;
  int  ProblemDimension;
  double X_Extent;
  double Y_Extent; 
  double Z_Extent;
  int  NoGridCells_X;
  int  NoGridCells_Y;
  int  NoGridCells_Z;
  int  NoTensorComponents;

  double BoundaryCondLeft_g;
  double BoundaryCondLeft_alpha1;
  double BoundaryCondLeft_alpha2;

  double BoundaryCondRight_g;
  double BoundaryCondRight_alpha1;
  double BoundaryCondRight_alpha2;

  double BoundaryCondTop_g;
  double BoundaryCondTop_alpha1;
  double BoundaryCondTop_alpha2;

  double BoundaryCondBottom_g;
  double BoundaryCondBottom_alpha1;
  double BoundaryCondBottom_alpha2;

  double BoundaryCondFront_g;
  double BoundaryCondFront_alpha1;
  double BoundaryCondFront_alpha2;

  double BoundaryCondBack_g;
  double BoundaryCondBack_alpha1;
  double BoundaryCondBack_alpha2;

  int  NumberOfWells;

  double BoundaryCondWell_g;
  double BoundaryCondWell_alpha1;
  double BoundaryCondWell_alpha2;

} TempSimInputData;


typedef struct
{
  int OutputFrequency;
  int NumbRockTypes;

  double T_Start;
  double T_End; 
  double delta_T;

  double DS_Max;

  double Viscosity_p1;
  double Viscosity_p2;
  double Viscosity_p3;
  double DensityGravity_p1;
  double DensityGravity_p2;
  double DensityGravity_p3;

  //Now comes the boundary condition info. connected to the saturation.

  double BoundaryCondLeft_g;
  double BoundaryCondLeft_alpha1;
  double BoundaryCondLeft_alpha2;
  int RockTypeLeft;

  double BoundaryCondRight_g;
  double BoundaryCondRight_alpha1;
  double BoundaryCondRight_alpha2;
  int RockTypeRight;

  double BoundaryCondTop_g;
  double BoundaryCondTop_alpha1;
  double BoundaryCondTop_alpha2;
  int RockTypeTop;

  double BoundaryCondBottom_g;
  double BoundaryCondBottom_alpha1;
  double BoundaryCondBottom_alpha2;
  int RockTypeBottom;

  double BoundaryCondFront_g;
  double BoundaryCondFront_alpha1;
  double BoundaryCondFront_alpha2;
  int RockTypeFront;

  double BoundaryCondBack_g;
  double BoundaryCondBack_alpha1;
  double BoundaryCondBack_alpha2;
  int RockTypeBack;

} IMPESInputData;


class PMTwoPhaseSimulationData
{
  public : 

  void readPMTwoPhaseSimulationData();

  //Data://@HAF: We prefer to let the data be public (at least for the time being).
  TempSimInputData td_;
  IMPESInputData IMPESd_;

 private:
  //local functions:
  bool ml_ReadSimInputData(const char *fn, TempSimInputData &tempSimData);
  bool ml_ReadIMPESInputData(const char *fn, IMPESInputData &IMPESData);
  
};


#endif //PMTwoPhaseSimulationData_H



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

