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


#ifndef PMTwoPhaseEllipticH
#define PMTwoPhaseEllipticH



//*****************************************************************************************

//---------------------- Include Files ----------------------------------------
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

#include <opm/porsol/twophase2/OPMIRISCode/IRISDuneGridInterface.hpp>
#include <opm/porsol/twophase2/OPMIRISCode/IRISDuneISTLInterface.hpp>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/GlobType.h>
#include <opm/porsol/twophase2/OPMIRISCode/IrisOpmTdefs.h>

#include <opm/porsol/twophase2/OPMIRISCode/SFRelPerm.h>
#include <opm/porsol/twophase2/OPMIRISCode/SFCapPress.h>
#include <opm/porsol/twophase2/OPMIRISCode/SFPMPhysicalFieldProperties.h>
#include <opm/porsol/twophase2/OPMIRISCode/SFMpfaFps.h>


//---------------------  Constants --------------------------------------------




//---------------------- typedef's and enumerations ----------------------------



//---------------------- Directives -------------------------------------------
using namespace std;  // This leads to that all names in the namespace std may be 
                      // used without the qualifier "std::".


//*********************************************************************************************



template <class FVMImplicit>   
class PMTwoPhaseElliptic
{
 public:
  //--------------------------------------------------------------------------
  //---------------------- Constructors --------------------------------------
  //--------------------------------------------------------------------------
  PMTwoPhaseElliptic(const int& numbElements);
  PMTwoPhaseElliptic(const FVMImplicit& FVMIMeth, const int& numbElements);

  //--------------------------------------------------------------------------
  //---------------------- Generators-----------------------------------------
  //--------------------------------------------------------------------------
  void setup(const IRISDuneGridInterface<DuneGridType>& Igrid, const PMTwoPhaseSimulationData& simulationData, const SFPMPhysicalFieldProperties& PMpp);
  void setupTwoPhasePressureEq(const IRISDuneGridInterface<DuneGridType>& Igrid, const SFSaturation& saturation);//@HAF: OK med referanse paa SFSaturation???
  void solveLinearEquationSystemForPressure();
  void reset();
  void computeVelocity(const IRISDuneGridInterface<DuneGridType>& Igrid, const SFSaturation& saturation);

  //--------------------------------------------------------------------------
  //---------------------- Observers -----------------------------------------
  //--------------------------------------------------------------------------
  const SFVelocity& getVelocity() const;//@HAF: OK med aa levere en referanse for SFVelocity???
  const SFPressure getPressure(const int& numbElements);//@HAF: OK med aa levere en referanse for SFPressure???

  //@HAF: Trenger vi en destructor???


  //-----------------------------------------------------------------------------------

  private : 
  //-----------------Useful local functions--------------------------------------------
  void ml_initializeVelocity(const IRISDuneGridInterface<DuneGridType>& Igrid);
  double ml_computePsiCAtCV(const int& ihat, const SFSaturation& saturation);
  double ml_computeMobilitiesAtCV(const int& ihat, const SFSaturation& saturation, const int& mobilityType);
  double ml_computeMobilitiesAtVertices(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer& qhatPt, const SFSaturation& saturation, const int& mobilityType);
  void ml_computeOnePhaseTransmissibilities(const IRISDuneGridInterface<DuneGridType>& Igrid);
  int ml_computeVelocityStorageIndex(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::ElementPointer& ihatPt, const int& ihat, const int& khat);

  //-----------------Data structure----------------------------------------------------

  FVMImplicit FVMMeth_;
  SFRelPerm rPerm_;
  SFCapPress cPress_;
  SFPMPhysicalFieldProperties PMpp_;
  IRISDuneISTLInterface IISTLI_;//The dune-framework for doing the linear eq. solve etc.
  SFVelocity velocity_;//@HAF: HUSK aat SFVelocity bare er en "typdefet" Array<double> eller lignende...

  //@HAF: Hva med lagring av transmissibiliteter???
  Array<double> transmissOnePhase_, transmissOnePhaseBC_;
  Array<double> transmissOnePhasePC_, transmissOnePhasePCBC_;
};




/*---------------------------------------------------*/
/*              Constructors                         */
/*---------------------------------------------------*/

template <class FVMImplicit>
PMTwoPhaseElliptic<FVMImplicit>::PMTwoPhaseElliptic(const int& numbElements)
: FVMMeth_(), rPerm_(), cPress_(), PMpp_(), IISTLI_(numbElements), velocity_(), transmissOnePhase_(), transmissOnePhaseBC_(), transmissOnePhasePC_(), transmissOnePhasePCBC_()
{

}

template <class FVMImplicit>
PMTwoPhaseElliptic<FVMImplicit>::PMTwoPhaseElliptic(const FVMImplicit& FVMIMeth, const int& numbElements)
: FVMMeth_(FVMIMeth), rPerm_(), cPress_(), PMpp_(), IISTLI_(numbElements), velocity_(), transmissOnePhase_(), transmissOnePhaseBC_(), transmissOnePhasePC_(), transmissOnePhasePCBC_()
{

}


/*---------------------------------------------------*/
/*              Generators                           */
/*---------------------------------------------------*/
template <class FVMImplicit>
void PMTwoPhaseElliptic<FVMImplicit>::setup(const IRISDuneGridInterface<DuneGridType>& Igrid, const PMTwoPhaseSimulationData& simulationData, const SFPMPhysicalFieldProperties& PMpp)//@HAF: Burde simulationData komme som en pointer???
{
  //The rel. perm. data
  //ifstream relPermInStream(relPermFormula_In_file);"relPermFormula.dta";
  ifstream relPermInStream("relPermFormula.dta");//@HAF: Should be with the file name as a parameter to the setup routine???
  rPerm_.readRelPermDataInFormulaFormat(relPermInStream,simulationData.IMPESd_.NumbRockTypes);
  relPermInStream.close();
  cout << endl << "File \"" << "relPermFormula.dta" << "\" is read...." << endl;

  //The cap. press. data
  //ifstream capPressInStream(capPressFormula_In_file);"capPressFormula.dta";
  ifstream capPressInStream("capPressFormula.dta");//@HAF: Should be with the file name as a parameter to the setup routine???
  cPress_.readCapPressDataInFormulaFormat(capPressInStream,simulationData.IMPESd_.NumbRockTypes);
  capPressInStream.close();
  cout << endl << "File \"" << "capPressFormula.dta" << "\" is read...." << endl;

  PMpp_ = PMpp;

  IISTLI_.setupMatrix(Igrid, 0);

  ml_initializeVelocity(Igrid);

  //@HAF:Hva med beregning av one-phase transmissibiliteter etc.???????
  ml_computeOnePhaseTransmissibilities(Igrid);//@HAF: Skal randen i det hele tatt involveres her??? Trenger vi egentlig noe eget for PC???
}


template <class FVMImplicit>
void PMTwoPhaseElliptic<FVMImplicit>::setupTwoPhasePressureEq(const IRISDuneGridInterface<DuneGridType>& Igrid, const SFSaturation& saturation)
{
  //@HAF: NOTE: It seems like mobilities will only be needed at the vertices. Could it also be necessary to compute them at the CV-centre as well i.e. see paper...???
  //@HAF: NOTE: From the point of view of structuring the software, this algorithm should possibly be placed inside FVMMeth_. However, for the time being it is more convenient to put it here. This can be changed later, at a relatively small cost.

  //@HAF: In the first place, this code is only for 2D.
  if (Igrid.cellDimension() == 2)
  {
    int idT = 0;
    int idTPC = 0;
    int idTBC = 0;
    int idTPCBC = 0;
    double integrationValue, integrationValue_pc, integrationValuePC;
    for (IRISDuneGridInterface<DuneGridType>::VertexIterator qhatIt = Igrid.setEntityPointerToFirst<IRISDuneGridInterface<DuneGridType>::VertexIterator>(0); qhatIt!=Igrid.entityPointerIsAtEnd<IRISDuneGridInterface<DuneGridType>::VertexIterator>(0); ++qhatIt)
    {
      IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt = qhatIt;
      bool qhatAtBoundary = Igrid.boundaryIndex(qhatPt);

      //-----Some mobilities at the vertex------------------------------------
      double lambda = ml_computeMobilitiesAtVertices(Igrid, qhatPt, saturation, 1);
      double lambda_0 = ml_computeMobilitiesAtVertices(Igrid, qhatPt, saturation, 0);
      //----------------------------------------------------------------------

      //Find the global index from the entity pointer:
      // int qhat = Igrid.getGlobalIndexFromEntityPointer(qhatPt);//@HAF: Kan vi droppe denne???
      int numbElemNeighs = Igrid.getNumberOfNeighbours(qhatPt, 0);

      for (int i=0; i < numbElemNeighs; i++)
      {
	IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt = Igrid.getGlobalEntityPointerFromRelativeIndex<IRISDuneGridInterface<DuneGridType>::ElementPointer>(qhatPt, i);//@HAF: OK???
	int ihat = Igrid.getGlobalIndexFromEntityPointer(ihatPt);

	//int outRelk[2];
	std::vector<int> outRelk(2,0);
	Igrid.getIntersectionRelativeEdgeIndexSetRightHandCSSorted(qhatPt, ihat, outRelk);//@HAF: First part of the old "findTopPath..."
	//int outRelz[2];
	//Igrid.getIntersectionRelativeElementIndexSetRightHandCSSorted<2>(qhatPt, ihat, outRelz);//@HAF: Second part of the old "findTopPath..."//NOTE: TRENGER VI EGENTLIG DENNE, OG I SAAFALL TIL HVA?????????????????????????????????????????????

	for (int iv = 0;iv < 2; iv++)//Since we are in 2D...
	{
	  IRISDuneGridInterface<DuneGridType>::EdgePointer khatPt = Igrid.getGlobalEntityPointerFromRelativeIndex<IRISDuneGridInterface<DuneGridType>::EdgePointer>(qhatPt, outRelk[iv]);

	  if (qhatAtBoundary)
	  {//The special boundary condition stuff.
	    integrationValue = lambda*transmissOnePhaseBC_[idTBC];//@HAF: Important: lambda at vertex or CV???
	    idTBC += 1;
	    integrationValuePC = lambda_0*transmissOnePhasePCBC_[idTPCBC];//@HAF: Important: lambda at vertex or CV???
	    idTPCBC += 1;
	    IISTLI_.rhsElement(ihat) += -integrationValue; 
	    IISTLI_.rhsElement(ihat) += -integrationValuePC;
	  }

	  integrationValuePC = 0.0;
	  for (int itilde = 0; itilde < numbElemNeighs; itilde++)
	  {
	    IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPtt = Igrid.getGlobalEntityPointerFromRelativeIndex<IRISDuneGridInterface<DuneGridType>::ElementPointer>(qhatPt, itilde);//@HAF: OK???
	    int ihatTilde = Igrid.getGlobalIndexFromEntityPointer(ihatPtt);
	    if (qhatAtBoundary)
	    {
	      integrationValue = lambda*transmissOnePhase_[idT];
	      idT += 1;
	      //FVMMeth_.generateLocalInverseMatrixAndRHS(Igrid, qhatPt, PMpp_, 1);//@HAF: NOTE this additional setup for PC!!! OK????????
	      //double valPC = FVMMeth_.computeOnePhaseTransmissibility(Igrid, qhatPt, PMpp_, outRelk, itilde, ihatPt, ihat, ihatTilde, khatPt);
	      //integrationValue_pc = lambda_0*valPC;//@HAF: Important: lambda at vertex or CV???
	      integrationValue_pc = lambda_0*transmissOnePhasePC_[idTPC];;//@HAF: Important: lambda at vertex or CV???
	      idTPC += 1;
	    }
	    else
	    {//Here we employ the precomputed one-phase transmissibilities:
	      //@HAF: NOTE: We use the same transmissibilites for "pressure" and PC, since the
	      //precomputed transmissibilities does NOT involve the boundary. OK???
	      integrationValue = lambda*transmissOnePhase_[idT];
	      integrationValue_pc = lambda_0*transmissOnePhase_[idT];
	      idT += 1;
	    }
	    IISTLI_.matrixElement(ihat,ihatTilde) += integrationValue;//@HAF: OK??? 
	    integrationValuePC += integrationValue_pc*ml_computePsiCAtCV(ihatTilde, saturation);
	  }
	  IISTLI_.rhsElement(ihat) += -integrationValuePC; 
	}
      }
    }
    //***************DEBUGGING***********START***************************
    //cout << "******************************************************************" << endl;
    //cout << "mat row 0= " << IISTLI_.matrixElement(0,0) << " , " << IISTLI_.matrixElement(0,1) << " , " << IISTLI_.matrixElement(0,2) << " , " << IISTLI_.matrixElement(0,3) << endl;
    //cout << "mat row 1= " << IISTLI_.matrixElement(1,0) << " , " << IISTLI_.matrixElement(1,1) << " , " << IISTLI_.matrixElement(1,2) << " , " << IISTLI_.matrixElement(1,3) << endl;
    //cout << "mat row 2= " << IISTLI_.matrixElement(2,0) << " , " << IISTLI_.matrixElement(2,1) << " , " << IISTLI_.matrixElement(2,2) << " , " << IISTLI_.matrixElement(2,3) << endl;
    //cout << "mat row 3= " << IISTLI_.matrixElement(3,0) << " , " << IISTLI_.matrixElement(3,1) << " , " << IISTLI_.matrixElement(3,2) << " , " << IISTLI_.matrixElement(3,3) << endl;
    //cout << "******************************************************************" << endl;
    //cout << "rhs= " << IISTLI_.rhsElement(0) << " , " << IISTLI_.rhsElement(1) << " , " << IISTLI_.rhsElement(2) << " , " << IISTLI_.rhsElement(3) << endl;
    //cout << "******************************************************************" << endl;
    //***************DEBUGGING***********END*****************************
  }
  else if (Igrid.cellDimension() == 3)
  {
    //3D code yet to come...
  }
}


template <class FVMImplicit>
void PMTwoPhaseElliptic<FVMImplicit>::solveLinearEquationSystemForPressure()
{
  IISTLI_.solveLinearSystem();
}


template <class FVMImplicit>
void PMTwoPhaseElliptic<FVMImplicit>::reset()
{
  velocity_ *= 0.0;
  IISTLI_.reset();
}


template <class FVMImplicit>
void PMTwoPhaseElliptic<FVMImplicit>::computeVelocity(const IRISDuneGridInterface<DuneGridType>& Igrid, const SFSaturation& saturation)
{
  //@HAF: NOTE: From the point of view of structuring the software, this algorithm should possibly be placed inside FVMMeth_. However, for the time being it is more convenient to put it here. This can be changed later, at a relatively small cost.

  //@HAF: In the first place, this code is only for 2D.
  if (Igrid.cellDimension() == 2)
  {
    int idT = 0;
    int idTBC = 0;
    double Fluxh, Fluxf;
    double integrationValue;//, integrationValue_pc, integrationValuePC;
    for (IRISDuneGridInterface<DuneGridType>::VertexIterator qhatIt = Igrid.setEntityPointerToFirst<IRISDuneGridInterface<DuneGridType>::VertexIterator>(0); qhatIt!=Igrid.entityPointerIsAtEnd<IRISDuneGridInterface<DuneGridType>::VertexIterator>(0); ++qhatIt)
    {
      IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt = qhatIt;
      bool qhatAtBoundary = Igrid.boundaryIndex(qhatPt);

      //-----A mobility at the vertex------------------------------------
      double lambda = ml_computeMobilitiesAtVertices(Igrid, qhatPt, saturation, 1);
      //----------------------------------------------------------------------

      //Find the global index from the entity pointer:
      // int qhat = Igrid.getGlobalIndexFromEntityPointer(qhatPt);//@HAF: Kan vi droppe denne???
      int numbElemNeighs = Igrid.getNumberOfNeighbours(qhatPt, 0);
      for (int i=0; i < numbElemNeighs; i++)
      {
	IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt = Igrid.getGlobalEntityPointerFromRelativeIndex<IRISDuneGridInterface<DuneGridType>::ElementPointer>(qhatPt, i);//@HAF: OK???
	int ihat = Igrid.getGlobalIndexFromEntityPointer(ihatPt);

	//int outRelk[2];
	std::vector<int> outRelk(2,0);
	Igrid.getIntersectionRelativeEdgeIndexSetRightHandCSSorted(qhatPt, ihat, outRelk);//@HAF: First part of the old "findTopPath..."
	//int outRelz[2];
	//Igrid.getIntersectionRelativeElementIndexSetRightHandCSSorted<2>(qhatPt, ihat, outRelz);//@HAF: Second part of the old "findTopPath..."//NOTE: TRENGER VI EGENTLIG DENNE, OG I SAAFALL TIL HVA?????????????????????????????????????????????

	for (int iv = 0;iv < 2; iv++)//Since we are in 2D...
	{
	  Fluxh = 0.0;
	  IRISDuneGridInterface<DuneGridType>::EdgePointer khatPt = Igrid.getGlobalEntityPointerFromRelativeIndex<IRISDuneGridInterface<DuneGridType>::EdgePointer>(qhatPt, outRelk[iv]);//@HAF: OK???
	  int khat = Igrid.getGlobalIndexFromEntityPointer(khatPt);
	  if (qhatAtBoundary)
	  {//The special boundary condition stuff.
	    integrationValue = lambda*transmissOnePhaseBC_[idTBC];//@HAF: Important: lambda at vertex or CV???
	    idTBC += 1;
	    Fluxh = integrationValue; 
	  }

	  Fluxf = 0.0;
	  for (int itilde = 0; itilde < numbElemNeighs; itilde++)
	  {
	    IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPtt = Igrid.getGlobalEntityPointerFromRelativeIndex<IRISDuneGridInterface<DuneGridType>::ElementPointer>(qhatPt, itilde);//@HAF: OK???
	    int ihatTilde = Igrid.getGlobalIndexFromEntityPointer(ihatPtt);
	   
	    integrationValue = lambda*transmissOnePhase_[idT];
	    idT += 1;
      
	    Fluxf += integrationValue*IISTLI_.solutionElement(ihatTilde);
	  }
	  int idx = ml_computeVelocityStorageIndex(Igrid, ihatPt, ihat, khat);
	  //velocity_.u_changeElement_1u(idx, velocity_[idx]+(Fluxh+Fluxf));
	  velocity_.u_changeElement_1u(idx, velocity_[idx]-(Fluxh+Fluxf));//@HAF: Change per. 02.03.2010.Correct???
	  //Eller skulle vi hatt velocity_.u_changeElement_1u(idx, velocity_[idx]-(-Fluxh+Fluxf))?????
	}
      }
    }
  }
  else if (Igrid.cellDimension() == 3)
  {
    //3D code yet to come...
  }
}


/*---------------------------------------------------*/
/*              Observers                            */
/*---------------------------------------------------*/

template <class FVMImplicit>
const SFVelocity& PMTwoPhaseElliptic<FVMImplicit>::getVelocity() const
{
  return velocity_;
}


template <class FVMImplicit>
const SFPressure PMTwoPhaseElliptic<FVMImplicit>::getPressure(const int& numbElements)
{
  SFPressure pressure(numbElements, 0.0);

  for (int i=0; i < numbElements; i++)
  {
    double val = IISTLI_.solutionElement(i);
    pressure.u_changeElement_1u(i, val);
  } 

  return pressure;
}


/*---------------------------------------------------*/
/*             Local functions                       */
/*---------------------------------------------------*/

template <class FVMImplicit>
void PMTwoPhaseElliptic<FVMImplicit>::ml_initializeVelocity(const IRISDuneGridInterface<DuneGridType>& Igrid)
{
  //@HAF: Maa sjekk gridtype (e.g. structured 2D or triangular)
  IRISDuneGridInterface<DuneGridType>::GridType gridType = Igrid.getGridType();//@HAF: gridType MUST of course be an enum (of type GridType see IRISDuneGridInterface)!!!

  //The "layout" of the velocity will be dependent on the type of the grid:
  //NOTE: Code only for 2D yet
  if (Igrid.cellDimension() == 2)
  {
    if (gridType == IRISDuneGridInterface<DuneGridType>::_Triangular)
    {
      velocity_.u_copy_1u(3*Igrid.cellCount(0,0), 0.0);
    }
    else if (gridType == IRISDuneGridInterface<DuneGridType>::_Structured2D)
    {
      velocity_.u_copy_1u(4*Igrid.cellCount(0,0), 0.0);
    }
  }
  else if (Igrid.cellDimension() == 3)
  {
    //3D code yet to come...
  }
}


template <class FVMImplicit>
void PMTwoPhaseElliptic<FVMImplicit>::ml_computeOnePhaseTransmissibilities(const IRISDuneGridInterface<DuneGridType>& Igrid)
{
  //@HAF: NOTE: From the point of view of structuring the software, this algorithm should possibly be placed inside FVMMeth_. However, for the time being it is more convenient to put it here. This can be changed later, at a relatively small cost.

  //@HAF: In the first place, this code is only for 2D.
  if (Igrid.cellDimension() == 2)
  {
    int idT = 0;
    int idTPC = 0;
    int idTBC = 0;
    int idTPCBC = 0;
    //Must first find the correct length of transmissOnePhase_ etc.
    for (IRISDuneGridInterface<DuneGridType>::VertexIterator qhatIt = Igrid.setEntityPointerToFirst<IRISDuneGridInterface<DuneGridType>::VertexIterator>(0); qhatIt!=Igrid.entityPointerIsAtEnd<IRISDuneGridInterface<DuneGridType>::VertexIterator>(0); ++qhatIt)
    {
      IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt = qhatIt;
      bool qhatAtBoundary = Igrid.boundaryIndex(qhatPt);
      int numbElemNeighs = Igrid.getNumberOfNeighbours(qhatPt, 0);
      for (int i=0; i < numbElemNeighs; i++)
      {
	for (int iv = 0;iv < 2; iv++)//Since we are in 2D...
	{
	  if (qhatAtBoundary)
	  {
	    idTBC += 1;
	    idTPCBC += 1;
	  }
	  for (int itilde = 0; itilde < numbElemNeighs; itilde++)
	  {
	    if (qhatAtBoundary)
	    {
	      idT += 1;
	      idTPC += 1;
	    }
	    else
	    {
	      idT += 1;
	    }
	  }
	}
      }
    }

    cout << "idT= " << idT << endl;
    cout << "idTBC= " << idTBC << endl;
    cout << "idTPC= " << idTPC << endl;
    cout << "idTPCBC= " << idTPCBC << endl;
    transmissOnePhase_.u_copy_1u(idT, 0.0); //@HAF: OK???
    transmissOnePhasePC_.u_copy_1u(idTPC, 0.0); //@HAF: OK???
    transmissOnePhaseBC_.u_copy_1u(idTBC, 0.0); //@HAF: OK???
    transmissOnePhasePCBC_.u_copy_1u(idTPCBC, 0.0); //@HAF: OK???
    //Must reset the indices to zero
    idT = 0;
    idTPC = 0;
    idTBC = 0;
    idTPCBC = 0;

    for (IRISDuneGridInterface<DuneGridType>::VertexIterator qhatIt = Igrid.setEntityPointerToFirst<IRISDuneGridInterface<DuneGridType>::VertexIterator>(0); qhatIt!=Igrid.entityPointerIsAtEnd<IRISDuneGridInterface<DuneGridType>::VertexIterator>(0); ++qhatIt)
    {
      IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt = qhatIt;
      bool qhatAtBoundary = Igrid.boundaryIndex(qhatPt);

      //------------------------------------------------------------------------------------
      //Setup for the local flux matrix and rhs.
      //However, NOTE: If qhat is on the boundary additional setups (in relevant routines) are necessary....
      
      FVMMeth_.generateLocalInverseMatrixAndRHS(Igrid, qhatPt, PMpp_, 0);
      //------------------------------------------------------------------------------------

      int numbElemNeighs = Igrid.getNumberOfNeighbours(qhatPt, 0);
      for (int i=0; i < numbElemNeighs; i++)
      {
	IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt = Igrid.getGlobalEntityPointerFromRelativeIndex<IRISDuneGridInterface<DuneGridType>::ElementPointer>(qhatPt, i);//@HAF: OK???
	int ihat = Igrid.getGlobalIndexFromEntityPointer(ihatPt);

	//int outRelk[2];
	std::vector<int> outRelk(2,0);
	Igrid.getIntersectionRelativeEdgeIndexSetRightHandCSSorted(qhatPt, ihat, outRelk);//@HAF: First part of the old "findTopPath..."
	//int outRelz[2];
	//Igrid.getIntersectionRelativeElementIndexSetRightHandCSSorted<2>(qhatPt, ihat, outRelz);//@HAF: Second part of the old "findTopPath..."//NOTE: TRENGER VI EGENTLIG DENNE, OG I SAAFALL TIL HVA?????????????????????????????????????????????
	for (int iv = 0;iv < 2; iv++)//Since we are in 2D...
	{
	  IRISDuneGridInterface<DuneGridType>::EdgePointer khatPt = Igrid.getGlobalEntityPointerFromRelativeIndex<IRISDuneGridInterface<DuneGridType>::EdgePointer>(qhatPt, outRelk[iv]);//@HAF: OK???

	  if (qhatAtBoundary)
	  {//The special boundary condition stuff.
	    FVMImplicit FVMMethBC_;
	    FVMMethBC_.generateLocalBoundaryInverseMatrixAndRHS(Igrid, qhatPt, PMpp_, cPress_, 0);
	    transmissOnePhaseBC_.u_changeElement_1u(idTBC, FVMMethBC_.computeOnePhaseBoundaryTransmissibility(Igrid, qhatPt, PMpp_, outRelk, ihatPt, ihat, khatPt));
	    idTBC += 1;
	  }

	  for (int itilde = 0; itilde < numbElemNeighs; itilde++)
	  {
	    IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPtt = Igrid.getGlobalEntityPointerFromRelativeIndex<IRISDuneGridInterface<DuneGridType>::ElementPointer>(qhatPt, itilde);//@HAF: OK???
	    int ihatTilde = Igrid.getGlobalIndexFromEntityPointer(ihatPtt);
	    
	    //Here we compute the one-phase transmissibilities, BUT only if qhat is NOT at the boundary.
	    //------------------------------------------------------------------------------------
	    transmissOnePhase_.u_changeElement_1u(idT, FVMMeth_.computeOnePhaseTransmissibility(Igrid, qhatPt, PMpp_, outRelk, itilde, ihatPt, ihat, ihatTilde, khatPt));
	    //cout << "idT= " << idT << " Transmissibility= " << transmissOnePhase_[idT] << endl;
	    idT += 1;
	  }
	}
      }
    }

    for (IRISDuneGridInterface<DuneGridType>::VertexIterator qhatIt = Igrid.setEntityPointerToFirst<IRISDuneGridInterface<DuneGridType>::VertexIterator>(0); qhatIt!=Igrid.entityPointerIsAtEnd<IRISDuneGridInterface<DuneGridType>::VertexIterator>(0); ++qhatIt)
    {
      IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt = qhatIt;
      bool qhatAtBoundary = Igrid.boundaryIndex(qhatPt);

      //------------------------------------------------------------------------------------
      //Setup for the local flux matrix and rhs.
      //However, NOTE: If qhat is on the boundary additional setups (in relevant routines) are necessary....
      if (qhatAtBoundary)
      {
	FVMMeth_.generateLocalInverseMatrixAndRHS(Igrid, qhatPt, PMpp_, 1);
      }
      //------------------------------------------------------------------------------------

      int numbElemNeighs = Igrid.getNumberOfNeighbours(qhatPt, 0);
      for (int i=0; i < numbElemNeighs; i++)
      {
	IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt = Igrid.getGlobalEntityPointerFromRelativeIndex<IRISDuneGridInterface<DuneGridType>::ElementPointer>(qhatPt, i);//@HAF: OK???
	int ihat = Igrid.getGlobalIndexFromEntityPointer(ihatPt);

	//int outRelk[2];
	std::vector<int> outRelk(2,0);
	Igrid.getIntersectionRelativeEdgeIndexSetRightHandCSSorted(qhatPt, ihat, outRelk);//@HAF: First part of the old "findTopPath..."
	//int outRelz[2];
	//Igrid.getIntersectionRelativeElementIndexSetRightHandCSSorted<2>(qhatPt, ihat, outRelz);//@HAF: Second part of the old "findTopPath..."//NOTE: TRENGER VI EGENTLIG DENNE, OG I SAAFALL TIL HVA?????????????????????????????????????????????
	for (int iv = 0;iv < 2; iv++)//Since we are in 2D...
	{
	  IRISDuneGridInterface<DuneGridType>::EdgePointer khatPt = Igrid.getGlobalEntityPointerFromRelativeIndex<IRISDuneGridInterface<DuneGridType>::EdgePointer>(qhatPt, outRelk[iv]);//@HAF: OK???

	  if (qhatAtBoundary)
	  {//The special boundary condition stuff.
	    FVMImplicit FVMMethBC_;
	    FVMMethBC_.generateLocalBoundaryInverseMatrixAndRHS(Igrid, qhatPt, PMpp_, cPress_, 1);
	    transmissOnePhasePCBC_.u_changeElement_1u(idTPCBC, FVMMethBC_.computeOnePhaseBoundaryTransmissibility(Igrid, qhatPt, PMpp_, outRelk, ihatPt, ihat, khatPt));
	    idTPCBC += 1;
	  }

	  for (int itilde = 0; itilde < numbElemNeighs; itilde++)
	  {
	    IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPtt = Igrid.getGlobalEntityPointerFromRelativeIndex<IRISDuneGridInterface<DuneGridType>::ElementPointer>(qhatPt, itilde);//@HAF: OK???
	    int ihatTilde = Igrid.getGlobalIndexFromEntityPointer(ihatPtt);
	    if (qhatAtBoundary)
	    {
	      //Here we compute the one-phase transmissibilities, BUT only if qhat IS at the boundary.
	      transmissOnePhasePC_.u_changeElement_1u(idTPC, FVMMeth_.computeOnePhaseTransmissibility(Igrid, qhatPt, PMpp_, outRelk, itilde, ihatPt, ihat, ihatTilde, khatPt));
	      //cout << "idT= " << idT << " Transmissibility= " << transmissOnePhase_[idT] << endl;
	      idTPC += 1;
	    }
	  }
	}
      }
    }

  }
  else if (Igrid.cellDimension() == 3)
  {
    //3D code yet to come...
  }
  //exit(0);
}


template <class FVMImplicit>
double PMTwoPhaseElliptic<FVMImplicit>::ml_computePsiCAtCV(const int& ihat, const SFSaturation& saturation)
{
  //@HAF: gravity is for the time being ignored in the computation  of PsiC
  return cPress_.getCapPressData(0, saturation[ihat], PMpp_.getRockType()[ihat]);
}


template <class FVMImplicit>
double PMTwoPhaseElliptic<FVMImplicit>::ml_computeMobilitiesAtCV(const int& ihat, const SFSaturation& saturation, const int& mobilityType)
{
  //NOTE: mobilityType==0 delivers lambda0

  double mu_o = PMpp_.getViscosity(0);
  double mu_w = PMpp_.getViscosity(1);
  double CVValue = saturation[ihat];
  double K_ro = rPerm_.getRelPermData(0, CVValue, PMpp_.getRockType()[ihat]);
  double K_rw = rPerm_.getRelPermData(1, CVValue, PMpp_.getRockType()[ihat]);
  
  return (mobilityType==0) ? K_ro/mu_o : (K_ro/mu_o)+(K_rw/mu_w);
}


template <class FVMImplicit>
double PMTwoPhaseElliptic<FVMImplicit>::ml_computeMobilitiesAtVertices(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer& qhatPt, const SFSaturation& saturation, const int& mobilityType)
{
  double volume;
  int numbNeigh = Igrid.getNumberOfNeighbours(qhatPt, 0);
  std::vector<IRISDuneGridInterface<DuneGridType>::ElementPointer> ihatPtElemNeigh; //@HAF: OK???
  //ihatPtElemNeigh.reserve(numbNeigh);
  Igrid.getNeighbours(qhatPt, ihatPtElemNeigh);//@HAF: OK???
  
  double sum = 0.0;
  double sum_N = 0.0;
  for (int j=0; j < numbNeigh; j++)
  {
    volume = Igrid.getVolume(ihatPtElemNeigh[j]);
    int ihat = Igrid.getGlobalIndexFromEntityPointer(ihatPtElemNeigh[j]);
    sum += ml_computeMobilitiesAtCV(ihat, saturation, mobilityType)*volume;  
    sum_N += volume;
  }

  return sum/sum_N;
}


template <class FVMImplicit>
int PMTwoPhaseElliptic<FVMImplicit>::ml_computeVelocityStorageIndex(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::ElementPointer& ihatPt, const int& ihat, const int& khat)
{
  int index = -1;
  int neighMax = -1;

  //@HAF: Maa sjekk gridtype (e.g. structured 2D or triangular)
  IRISDuneGridInterface<DuneGridType>::GridType gridType = Igrid.getGridType();//@HAF: gridType MUST of course be an enum (of type GridType see IRISDuneGridInterface)!!!

  //NOTE: Code only for 2D yet
  if (Igrid.cellDimension() == 2)
  {
    if (gridType == IRISDuneGridInterface<DuneGridType>::_Triangular)
    {
      neighMax = 3;
    }
    else if (gridType == IRISDuneGridInterface<DuneGridType>::_Structured2D)
    {
      neighMax = 4;
    }
    bool OK = false;
    for (int k=0; (k < neighMax) && !OK; k++)
    {
      IRISDuneGridInterface<DuneGridType>::EdgePointer khatPt = Igrid.getGlobalEntityPointerFromRelativeIndex<IRISDuneGridInterface<DuneGridType>::EdgePointer>(ihatPt, k);//@HAF: OK???
      int khatLoc = Igrid.getGlobalIndexFromEntityPointer(khatPt);
      if (khatLoc == khat) 
      {
	index = neighMax*ihat + k;
	OK = true;
      }
    }
  }
  else if (Igrid.cellDimension() == 3)
  {
    //3D code yet to come...
  }

  if (index == -1)
  {
    cout << "Something is very wrong in ml_computeVelocityStorageIndex. Program is stopped." << endl;
    exit(0);
  }

  return index;
}


#endif //PMTwoPhaseEllipticH

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




