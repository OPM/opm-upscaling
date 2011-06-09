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


//---------------------- Include Files ----------------------------------------
#include "PMTwoPhaseSimulationData.h"


//---------------------- Directives -------------------------------------------
using namespace std;  // This leads to that all names in the namespace std may be 
                      // used without the qualifier "std::".


void PMTwoPhaseSimulationData::readPMTwoPhaseSimulationData()
{
  //@HAF: Obviously this can be done much more general!!!
  char tempSimInputDataFile[]     = "TempSimInputData.dta";
  ml_ReadSimInputData(tempSimInputDataFile, td_);
  cout << endl << "------------------------------------"
       << endl << "File \"" << tempSimInputDataFile << "\" is read....";

  char IMPESInputDataFile[]     = "IMPESInputData.dta"; 
  ml_ReadIMPESInputData(IMPESInputDataFile, IMPESd_);
  cout << endl << "------------------------------------"
       << endl << "File \"" << IMPESInputDataFile << "\" is read....";

}

//------------------------------------------------------------
//Local fuctions
//------------------------------------------------------------

//-----------------------------------------------------------------------------
//
//  ReadSimInputData
//
//  DESCRIPTION:
//    Read temperature simulation data from file.
//
//  RETURNS:
//    True if succeeded, else false.
//
//  NOTE:
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------

bool PMTwoPhaseSimulationData::ml_ReadSimInputData(const char *fn, TempSimInputData &tempSimData)
{
  FILE  *fp = fopen(fn,"r");
  char  versionStr[FILE_INPUT_STR_LEN];  // The file version string
  char  fileRecStr[FILE_INPUT_STR_LEN];  // File recognition string
  char  inStr[FILE_INPUT_STR_LEN];       // Input string
  char  parName[FILE_INPUT_STR_LEN];     // Parameter name string
  char  dummy[FILE_INPUT_STR_LEN];       // Dummy string
  int   version;                         // Version number
  int   revision;                        // Revision number

  if (fp)
  {
    if ( fscanf(fp, "%s%d.%d\n%s\n\n", versionStr, &version, &revision, 
                   fileRecStr) == 4 )
    {
      while ( fgets (inStr, FILE_INPUT_STR_LEN, fp) )
      {
        parName[0] = '\0';
        sscanf(inStr, "%s", parName);

        if ( strcmp(parName, "NameOfGridType") == 0)
        {
          sscanf(inStr, "%s%s%s", parName, dummy, tempSimData.NameOfGridType);
        }

        else if ( strcmp(parName, "NameOfDGFFile") == 0)
        {
          sscanf(inStr, "%s%s%s", parName, dummy, tempSimData.NameOfDGFFile);
        }

        else if ( strcmp(parName, "DiscretizationType") == 0)
        {
          sscanf(inStr, "%s%s%d", parName, dummy, &tempSimData.DiscretizationType);
        }

        else if ( strcmp(parName, "UnStructuredGrid") == 0)
        {
          sscanf(inStr, "%s%s%d", parName, dummy, &tempSimData.UnStructuredGrid);
        }

        else if ( strcmp(parName, "ReadStructuredGrid") == 0)
        {
          sscanf(inStr, "%s%s%d", parName, dummy, &tempSimData.ReadStructuredGrid);
        }

        else if ( strcmp(parName, "PeriodicBoundaryConditions") == 0)
        {
          sscanf(inStr, "%s%s%d", parName, dummy, &tempSimData.PeriodicBoundaryConditions);
        }

        else if ( strcmp(parName, "ProblemDimension") == 0)
        {
          sscanf(inStr, "%s%s%d", parName, dummy, &tempSimData.ProblemDimension);
        }

        else if ( strcmp(parName, "X_Extent") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &tempSimData.X_Extent);
        }

        else if ( strcmp(parName, "Y_Extent") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &tempSimData.Y_Extent);
        }

        else if ( strcmp(parName, "Z_Extent") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &tempSimData.Z_Extent);
        }

        else if ( strcmp(parName, "NoGridCells_X") == 0)
        {
          sscanf(inStr, "%s%s%d", parName, dummy, &tempSimData.NoGridCells_X);
        }

        else if ( strcmp(parName, "NoGridCells_Y") == 0)
        {
          sscanf(inStr, "%s%s%d", parName, dummy, &tempSimData.NoGridCells_Y);
        }

        else if ( strcmp(parName, "NoGridCells_Z") == 0)
        {
          sscanf(inStr, "%s%s%d", parName, dummy, &tempSimData.NoGridCells_Z);
        }

        else if ( strcmp(parName, "NoTensorComponents") == 0)
        {
          sscanf(inStr, "%s%s%d", parName, dummy, &tempSimData.NoTensorComponents);
        }


        else if ( strcmp(parName, "BoundaryCondLeft_g") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &tempSimData.BoundaryCondLeft_g);
        }

        else if ( strcmp(parName, "BoundaryCondLeft_alpha1") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &tempSimData.BoundaryCondLeft_alpha1);
        }

        else if ( strcmp(parName, "BoundaryCondLeft_alpha2") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &tempSimData.BoundaryCondLeft_alpha2);
        }


        else if ( strcmp(parName, "BoundaryCondRight_g") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &tempSimData.BoundaryCondRight_g);
        }

        else if ( strcmp(parName, "BoundaryCondRight_alpha1") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &tempSimData.BoundaryCondRight_alpha1);
        }

        else if ( strcmp(parName, "BoundaryCondRight_alpha2") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &tempSimData.BoundaryCondRight_alpha2);
        }


        else if ( strcmp(parName, "BoundaryCondTop_g") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &tempSimData.BoundaryCondTop_g);
        }

        else if ( strcmp(parName, "BoundaryCondTop_alpha1") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &tempSimData.BoundaryCondTop_alpha1);
        }

        else if ( strcmp(parName, "BoundaryCondTop_alpha2") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &tempSimData.BoundaryCondTop_alpha2);
        }


        else if ( strcmp(parName, "BoundaryCondBottom_g") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &tempSimData.BoundaryCondBottom_g);
        }

        else if ( strcmp(parName, "BoundaryCondBottom_alpha1") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &tempSimData.BoundaryCondBottom_alpha1);
        }

        else if ( strcmp(parName, "BoundaryCondBottom_alpha2") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &tempSimData.BoundaryCondBottom_alpha2);
        }


        else if ( strcmp(parName, "BoundaryCondFront_g") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &tempSimData.BoundaryCondFront_g);
        }

        else if ( strcmp(parName, "BoundaryCondFront_alpha1") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &tempSimData.BoundaryCondFront_alpha1);
        }

        else if ( strcmp(parName, "BoundaryCondFront_alpha2") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &tempSimData.BoundaryCondFront_alpha2);
        }


        else if ( strcmp(parName, "BoundaryCondBack_g") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &tempSimData.BoundaryCondBack_g);
        }

        else if ( strcmp(parName, "BoundaryCondBack_alpha1") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &tempSimData.BoundaryCondBack_alpha1);
        }

        else if ( strcmp(parName, "BoundaryCondBack_alpha2") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &tempSimData.BoundaryCondBack_alpha2);
        }

        else if ( strcmp(parName, "NumberOfWells") == 0)
        {
          sscanf(inStr, "%s%s%d", parName, dummy, &tempSimData.NumberOfWells);
	}

	else if ( strcmp(parName, "BoundaryCondWell_g") == 0)
	{
	  sscanf(inStr, "%s%s%lf", parName, dummy, &tempSimData.BoundaryCondWell_g);
	}

	else if ( strcmp(parName, "BoundaryCondWell_alpha1") == 0)
	{
	  sscanf(inStr, "%s%s%lf", parName, dummy, &tempSimData.BoundaryCondWell_alpha1);
	}

	else if ( strcmp(parName, "BoundaryCondWell_alpha2") == 0)
	{
	  sscanf(inStr, "%s%s%lf", parName, dummy, &tempSimData.BoundaryCondWell_alpha2);
	}

      }
    }
    else
    {
      fprintf(stderr, "\n\nReadSimInputData : " 
        "Temperature simulation input file has wrong format...");
      assert(0);
    }
  }
  fclose(fp);
  return(true);

}//ReadSimInputData




//-----------------------------------------------------------------------------
//
//  ReadIMPESInputData
//
//  DESCRIPTION:
//    Read IMPES simulation data from file.
//
//  RETURNS:
//    True if succeeded, else false.
//
//  NOTE:
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------

bool PMTwoPhaseSimulationData::ml_ReadIMPESInputData(const char *fn, IMPESInputData &IMPESData)
{
  FILE  *fp = fopen(fn,"r");
  char  versionStr[FILE_INPUT_STR_LEN];  // The file version string
  char  fileRecStr[FILE_INPUT_STR_LEN];  // File recognition string
  char  inStr[FILE_INPUT_STR_LEN];       // Input string
  char  parName[FILE_INPUT_STR_LEN];     // Parameter name string
  char  dummy[FILE_INPUT_STR_LEN];       // Dummy string
  int   version;                         // Version number
  int   revision;                        // Revision number

  if (fp)
  {
    if ( fscanf(fp, "%s%d.%d\n%s\n\n", versionStr, &version, &revision, 
                   fileRecStr) == 4 )
    {
      while ( fgets (inStr, FILE_INPUT_STR_LEN, fp) )
      {
        parName[0] = '\0';
        sscanf(inStr, "%s", parName);

        if ( strcmp(parName, "OutputFrequency") == 0)
        {
          sscanf(inStr, "%s%s%d", parName, dummy, &IMPESData.OutputFrequency);
        }

        else if ( strcmp(parName, "NumbRockTypes") == 0)
        {
          sscanf(inStr, "%s%s%d", parName, dummy, &IMPESData.NumbRockTypes);
        }

        else if ( strcmp(parName, "T_Start") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.T_Start);
        }

        else if ( strcmp(parName, "T_End") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.T_End);
        }

        else if ( strcmp(parName, "delta_T") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.delta_T);
        }

        else if ( strcmp(parName, "DS_Max") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.DS_Max);
        }

       else if ( strcmp(parName, "Viscosity_p1") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.Viscosity_p1);
        }

        else if ( strcmp(parName, "Viscosity_p2") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.Viscosity_p2);
        }

        else if ( strcmp(parName, "Viscosity_p3") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.Viscosity_p3);
        }

        else if ( strcmp(parName, "DensityGravity_p1") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.DensityGravity_p1);
        }

        else if ( strcmp(parName, "DensityGravity_p2") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.DensityGravity_p2);
        }

        else if ( strcmp(parName, "DensityGravity_p3") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.DensityGravity_p3);
        }

        else if ( strcmp(parName, "BoundaryCondLeft_g") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.BoundaryCondLeft_g);
        }

        else if ( strcmp(parName, "BoundaryCondLeft_alpha1") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.BoundaryCondLeft_alpha1);
        }

        else if ( strcmp(parName, "BoundaryCondLeft_alpha2") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.BoundaryCondLeft_alpha2);
        }

        else if ( strcmp(parName, "RockTypeLeft") == 0)
        {
          sscanf(inStr, "%s%s%d", parName, dummy, &IMPESData.RockTypeLeft);
        }

        else if ( strcmp(parName, "BoundaryCondRight_g") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.BoundaryCondRight_g);
        }

        else if ( strcmp(parName, "BoundaryCondRight_alpha1") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.BoundaryCondRight_alpha1);
        }

        else if ( strcmp(parName, "BoundaryCondRight_alpha2") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.BoundaryCondRight_alpha2);
        }

        else if ( strcmp(parName, "RockTypeRight") == 0)
        {
          sscanf(inStr, "%s%s%d", parName, dummy, &IMPESData.RockTypeRight);
        }

        else if ( strcmp(parName, "BoundaryCondTop_g") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.BoundaryCondTop_g);
        }

        else if ( strcmp(parName, "BoundaryCondTop_alpha1") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.BoundaryCondTop_alpha1);
        }

        else if ( strcmp(parName, "BoundaryCondTop_alpha2") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.BoundaryCondTop_alpha2);
        }

        else if ( strcmp(parName, "RockTypeTop") == 0)
        {
          sscanf(inStr, "%s%s%d", parName, dummy, &IMPESData.RockTypeTop);
        }

        else if ( strcmp(parName, "BoundaryCondBottom_g") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.BoundaryCondBottom_g);
        }

        else if ( strcmp(parName, "BoundaryCondBottom_alpha1") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.BoundaryCondBottom_alpha1);
        }

        else if ( strcmp(parName, "BoundaryCondBottom_alpha2") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.BoundaryCondBottom_alpha2);
        }

        else if ( strcmp(parName, "RockTypeBottom") == 0)
        {
          sscanf(inStr, "%s%s%d", parName, dummy, &IMPESData.RockTypeBottom);
        }

        else if ( strcmp(parName, "BoundaryCondFront_g") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.BoundaryCondFront_g);
        }

        else if ( strcmp(parName, "BoundaryCondFront_alpha1") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.BoundaryCondFront_alpha1);
        }

        else if ( strcmp(parName, "BoundaryCondFront_alpha2") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.BoundaryCondFront_alpha2);
        }

        else if ( strcmp(parName, "RockTypeFront") == 0)
        {
          sscanf(inStr, "%s%s%d", parName, dummy, &IMPESData.RockTypeFront);
        }

        else if ( strcmp(parName, "BoundaryCondBack_g") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.BoundaryCondBack_g);
        }

        else if ( strcmp(parName, "BoundaryCondBack_alpha1") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.BoundaryCondBack_alpha1);
        }

        else if ( strcmp(parName, "BoundaryCondBack_alpha2") == 0)
        {
          sscanf(inStr, "%s%s%lf", parName, dummy, &IMPESData.BoundaryCondBack_alpha2);
        }

        else if ( strcmp(parName, "RockTypeBack") == 0)
        {
          sscanf(inStr, "%s%s%d", parName, dummy, &IMPESData.RockTypeBack);
        }

      }
    }
    else
    {
      fprintf(stderr, "\n\nReadIMPESInputData : " 
        "IMPES simulation input file has wrong format...");
      assert(0);
    }
  }
  fclose(fp);
  return(true);

}//ReadIMPESInputData
