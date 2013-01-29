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
#include "SFPMPhysicalFieldProperties.h"

#include <iostream>
#include <iomanip>
#include <stdexcept>

//---------------------- Directives -------------------------------------------
using namespace std;  // This leads to that all names in the namespace std may be 
                      // used without the qualifier "std::".





SFPMPhysicalFieldProperties::SFPMPhysicalFieldProperties()  : rockType_(), Kxx_(), Kyy_(), Kzz_(), Kxy_(), Kxz_(), Kyz_(), porosity_()
{
  for (int i=0; i < 3; i++)
  {
    viscosity_[i] = 0.0;
    density_gravity_[i] = 0.0;
  }
}


void SFPMPhysicalFieldProperties::setup(const IRISDuneGridInterface<DuneGridType>* IGIPtr, const PMTwoPhaseSimulationData* simulationDataPtr)
{
  IGIPtr_ = IGIPtr;
  simulationDataPtr_ = simulationDataPtr;//@HAF: OK???

  viscosity_[0] = simulationDataPtr_->IMPESd_.Viscosity_p1;//@HAF: OK???
  viscosity_[1] = simulationDataPtr_->IMPESd_.Viscosity_p2;//@HAF: OK???
  viscosity_[2] = simulationDataPtr_->IMPESd_.Viscosity_p3;//@HAF: OK???
  density_gravity_[0] = simulationDataPtr_->IMPESd_.DensityGravity_p1;//@HAF: OK???
  density_gravity_[1] = simulationDataPtr_->IMPESd_.DensityGravity_p2;//@HAF: OK???
  density_gravity_[2] = simulationDataPtr_->IMPESd_.DensityGravity_p3;//@HAF: OK???

  //Here we must handle permeability, rockType and porosity:
  ml_fillInForPermeabilityRockTypeAndPorosity();
}



const PMTwoPhaseSimulationData* SFPMPhysicalFieldProperties::getSimulationDataPtr() const
{
  return simulationDataPtr_;
}


const Array<int>& SFPMPhysicalFieldProperties::getRockType() const
{
  return rockType_;
}


const Array<double>& SFPMPhysicalFieldProperties::getPermeability(const int tensorComp) const
{
  //@HAF: Is the following implementation OK??? i.e. usable???
  if (IGIPtr_->domainDimension() == 2)
  {
    if (tensorComp == 0)
    {
      return Kxx_;
    }
    else if (tensorComp == 1)
    {
      return Kyy_;
    }
    else if (tensorComp == 2)
    {
      return Kxy_;
    }
  }
  else if (IGIPtr_->domainDimension() == 3)
  {
    if (tensorComp == 0)
    {
      return Kxx_;
    }
    else if (tensorComp == 1)
    {
      return Kyy_;
    }
    else if (tensorComp == 2)
    {
      return Kzz_;
    }
    else if (tensorComp == 3)
    {
      return Kxy_;
    }
    else if (tensorComp == 4)
    {
      return Kxz_;
    }
    else if (tensorComp == 5)
    {
      return Kyz_;
    }
  }
  throw std::runtime_error("Error in call to getPermeability().");
}


const Array<double>& SFPMPhysicalFieldProperties::getPorosity() const
{
  return porosity_;
}


double SFPMPhysicalFieldProperties::getViscosity(const int phase) const
{
  /*
    sw is the water saturation.
    phase=0 : oil
    phase=1 : water
    phase=3 : gas
  */

  return viscosity_[phase];
}


double SFPMPhysicalFieldProperties::getDensityGravity(const int phase) const
{
  /*
    sw is the water saturation.
    phase=0 : oil
    phase=1 : water
    phase=3 : gas
  */

  return density_gravity_[phase];
}



void SFPMPhysicalFieldProperties::operator=(const SFPMPhysicalFieldProperties& PMpp) 
{
  if (this != &PMpp)
  {
    IGIPtr_ = PMpp.IGIPtr_;
    simulationDataPtr_ = PMpp.simulationDataPtr_;

    rockType_ = PMpp.rockType_;
    Kxx_ = PMpp.Kxx_; 
    Kyy_ = PMpp.Kyy_;
    Kzz_ = PMpp.Kzz_;
    Kxy_ = PMpp.Kxy_;
    Kxz_ = PMpp.Kxz_;
    Kyz_ = PMpp.Kyz_;
    porosity_ = PMpp.porosity_;
 
    for (int i=0; i < 3; i++)
    {
      viscosity_[i] = PMpp.viscosity_[i];
      density_gravity_[i] = PMpp.density_gravity_[i];
    }
  }
}


SFPMPhysicalFieldProperties::~SFPMPhysicalFieldProperties() 
{

}





//*************************************************************
//*************************************************************
// Local functions ********************************************
//*************************************************************
//*************************************************************


void SFPMPhysicalFieldProperties::ml_fillInForPermeabilityRockTypeAndPorosity()
{
  //@HAF: The implementation here will be highly dependent on the type of grid used.
  //E.g. a triangular grid may (or may not) need certain quantities to be read etc.
  //In the first place we only provide a simple implemetation for a structured grid.

  IRISDuneGridInterface<DuneGridType>::GridType gridType = IGIPtr_->getGridType();//@HAF: gridType MUST of course be an enum (of type GridType see IRISDuneGridInterface)!!!
  if (gridType == IRISDuneGridInterface<DuneGridType>::_Triangular)
  {
    int L = IGIPtr_->cellCount(0,0);
    rockType_.u_copy_1u(L, 0);
    Kxx_.u_copy_1u(L, 0.0);
    Kyy_.u_copy_1u(L, 0.0); 
    Kxy_.u_copy_1u(L, 0.0);
    porosity_.u_copy_1u(L, 0.0);
    ml_readPhysicalFieldPropertiesForTriangularGridFromFile();
  }
  else if (gridType == IRISDuneGridInterface<DuneGridType>::_Structured2D)
  {
    int L = (simulationDataPtr_->td_.NoGridCells_X)*(simulationDataPtr_->td_.NoGridCells_Y);
    rockType_.u_copy_1u(L, 0);
    Kxx_.u_copy_1u(L, 0.0);
    Kyy_.u_copy_1u(L, 0.0); 
    Kxy_.u_copy_1u(L, 0.0);
    porosity_.u_copy_1u(L, 0.0);
    ml_generatePhysicalFieldProperties_DisplacementProb_StructGrid();//@HAF: NOTE: This is just an example...
  }
}


void SFPMPhysicalFieldProperties::ml_generatePhysicalFieldProperties_DisplacementProb_StructGrid()
{
  //***********************************************************************************
  //We generate the permeability, rock type and porosity distribution
  //for the (isotropic) displacement example on a structured grid (in a unit rectangular domain):
  //***********************************************************************************

  int NX = simulationDataPtr_->td_.NoGridCells_X;
  int NY = simulationDataPtr_->td_.NoGridCells_Y;
  double DX = simulationDataPtr_->td_.X_Extent/NX;
  double DY = simulationDataPtr_->td_.Y_Extent/NY;

  double length = simulationDataPtr_->td_.Y_Extent;
  double perm1 = 1.0;
  double perm2 = 0.01;
  //double perm2 = 1.0;

  double posX = 0.0;
  double posY = 0.0;
  int k=0;
  for (int j=0; j < NY; j++)
  {
    posY = 0.5*DY + j*DY;
    for (int i=0; i < NX; i++)
    {
      posX = 0.5*DX + i*DX;
      k = i + j*NX;
      if ((posY < (length - (1.0/3.0)*length)) && (posY > (1.0/3.0)*length)) //A 2D example.
	//if ((posX >= 0.25) && (posX <= 0.5) && (posY >= 0.25) && (posY <= 0.75)) //A 2D unit square example.
	//if ((posX >= 0.25) && (posX <= 0.5)) //A 1D example
      {
	Kxx_.u_changeElement_1u(k,perm2);
	Kyy_.u_changeElement_1u(k,perm2);
	porosity_.u_changeElement_1u(k,1.0);
	rockType_.u_changeElement_1u(k,1);
      }
      else
      {
	Kxx_.u_changeElement_1u(k,perm1);
	Kyy_.u_changeElement_1u(k,perm1);
	porosity_.u_changeElement_1u(k,1.0);
	rockType_.u_changeElement_1u(k,0);//@HAF: MEGET VIKTIG: Er dette OK???
      }
    }
  }
}


void SFPMPhysicalFieldProperties::ml_readPhysicalFieldPropertiesForTriangularGridFromFile()
{
  //***********************************************************************************
  //Here we read the permeability, rock type and porosity distribution
  //from file.
  //We note that these files should be generated by external software.
  //***********************************************************************************

  //***********************************************************************************
  //The way this is presently done, there are some prerequisites: 
  //a) The name of the files MUST be as below (i.e. rockType.dta, Kxx.dta, Kyy.dta, 
  //Kxy.dta and rockType.dta).
  //b) Each files consits of L values (where L is number of elements in the corresponding
  //triangular grid), seperated by whitespaces, but (NOTE nevertheless) in one
  //single line. This may seem akward, but steems from the fact that these files are generated
  //by external software on so called sgml-format. However, the sgml-tags (also
  //including the number of elements in the file) in both ends of the files MUST have been
  //removed before the files are read in this routine!!!
  //***********************************************************************************
  int fidx;

  ifstream rockTypeInStream("rockType.dta");
  fidx = 0;
  while (!rockTypeInStream.eof())
  {
    rockTypeInStream >> rockType_[fidx];
    fidx++;
  }
  rockTypeInStream.close();

 ifstream permInStream("Kxx.dta");
  fidx = 0;
  while (!permInStream.eof())
  {
    permInStream >> Kxx_[fidx];
    fidx++;
  }
  permInStream.close();

  permInStream.open("Kyy.dta");
  fidx = 0;
  while (!permInStream.eof())
  {
    permInStream >> Kyy_[fidx];
    fidx++;
  }
  permInStream.close();

  permInStream.open("Kxy.dta");
  fidx = 0;
  while (!permInStream.eof())
  {
    permInStream >> Kxy_[fidx];
    fidx++;
  }
  permInStream.close();

 ifstream porosityInStream("porosity.dta");
  fidx = 0;
  while (!porosityInStream.eof())
  {
    porosityInStream >> porosity_[fidx];
    fidx++;
  }
  porosityInStream.close();
}

