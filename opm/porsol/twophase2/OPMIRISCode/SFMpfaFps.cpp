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
#include "SFMpfaFps.h"

#include <cmath>


//---------------------  Constants --------------------------------------------
#define TOLERANCE 10e-9


//---------------------- Types ------------------------------------------------








//---------------------- Public Functions -------------------------------------





// Constructors

//-----------------------------------------------------------------------------
SFMpfaFps::SFMpfaFps()
  : invMmatrix_(), hmatrix_()
//-----------------------------------------------------------------------------
{

}


// Copy Constructors

//-----------------------------------------------------------------------------
SFMpfaFps::SFMpfaFps(const SFMpfaFps& Fps)
  : invMmatrix_(Fps.invMmatrix_), hmatrix_(Fps.hmatrix_)
//-----------------------------------------------------------------------------
{

}

// Destructors

//-----------------------------------------------------------------------------
SFMpfaFps::~SFMpfaFps()
//-----------------------------------------------------------------------------
{

}



// Generators

double SFMpfaFps::computeOnePhaseTransmissibility(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const SFPMPhysicalFieldProperties& PMpp, const std::vector<int>& outRelk, const int& itilde, const IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt, const int& ihat, const int& ihatTilde, const IRISDuneGridInterface<DuneGridType>::EdgePointer khatPt)//@HAF: OK???
{
  double transRes = 0.0;

  const Array<double>& K11 = PMpp.getPermeability(0);//@HAF: OK with reference here???
  const Array<double>& K22 = PMpp.getPermeability(1);
  const Array<double>& K12 = PMpp.getPermeability(2); 

  //The rest of the transmissibility algorithm:

  int normalVectorLabel = ml_getTransformSpaceNormalVectorLabel(Igrid, qhatPt, ihatPt, khatPt);

  double deltaF = (ihat == ihatTilde) ? 1.0 : 0.0;//@HAF: OK???

  double biLinearIntegral = ml_FE2(Igrid, qhatPt, ihatPt, ihat, normalVectorLabel, K11, K22, K12);
  for (int j=0; j < 2; j++)
  {
    double shapeFuncCoeff = ml_matMultMH(outRelk[j],itilde);
    transRes += (shapeFuncCoeff - deltaF)*(ml_FE1(Igrid, qhatPt, ihatPt, ihat, j, normalVectorLabel, K11, K22, K12)) - shapeFuncCoeff*biLinearIntegral;
  }
		  
  //The bilinear contribution:
  int kMax = Igrid.getNumberOfNeighbours(qhatPt, 1);//@HAF: OK???
  transRes += (ml_matMultMH(kMax,itilde) + deltaF)*biLinearIntegral;

  return transRes;
}


double SFMpfaFps::computeOnePhaseBoundaryTransmissibility(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const SFPMPhysicalFieldProperties& PMpp, const std::vector<int>& outRelk, const IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt, const int& ihat, const IRISDuneGridInterface<DuneGridType>::EdgePointer khatPt)//@HAF: OK???
{
  double transRes = 0.0;

  const Array<double>& K11 = PMpp.getPermeability(0);//@HAF: OK with reference here???
  const Array<double>& K22 = PMpp.getPermeability(1);
  const Array<double>& K12 = PMpp.getPermeability(2); 

  //The rest of the transmissibility algorithm:

  int normalVectorLabel = ml_getTransformSpaceNormalVectorLabel(Igrid, qhatPt, ihatPt, khatPt);

  double biLinearIntegral = ml_FE2(Igrid, qhatPt, ihatPt, ihat, normalVectorLabel, K11, K22, K12);
  for (int j=0; j < 2; j++)
  {
    transRes += ml_matMultMH(outRelk[j],0)*(ml_FE1(Igrid, qhatPt, ihatPt, ihat, j, normalVectorLabel, K11, K22, K12) - biLinearIntegral);
  }
		  
  //The bilinear contribution:
  int kMax = Igrid.getNumberOfNeighbours(qhatPt, 1);//@HAF: OK???
  transRes += ml_matMultMH(kMax,0)*biLinearIntegral;

  return transRes;
}


void SFMpfaFps::generateLocalInverseMatrixAndRHS(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const SFPMPhysicalFieldProperties& PMpp, const int& pressureProblemType)//@HAF: OK???
{
  /*
    ------------IMPRTANT INFO.-------------------------------------------
    ---------------------------------------------------------------------
    pressureProblemType == 0 -> Pressure problem for psi.
    pressureProblemType == 1 -> Pressure problem for PC.
    ---------------------------------------------------------------------
  */
  
  //Fix BC-stuff etc.
  const PMTwoPhaseSimulationData* simDataPtr = PMpp.getSimulationDataPtr();
  Array<double> alpha1;
  Array<double> alpha2;
  int BCLL = 2*(*simDataPtr).td_.ProblemDimension; 
  Array<double> alpha1ExternalSurfaces(BCLL); 
  Array<double> alpha2ExternalSurfaces(BCLL); 

  if (Igrid.boundaryIndex(qhatPt))//@HAF: Is it OK to only allocate alpha1 and alpha2 in THIS case???
  {
    if (pressureProblemType == 0)
    {
      //Here we treat the problem for psi:
      //creates boundary data for Robin boundary conditions.

      //----------------------------------------------------------------------
      //NB: Vi maa vite hvilken av sidene nummerne som representerer hva 
      //geometrisk her.
      //Below we ASSUME the following to be correct:
      //In 2D:
      //0 and 2: verical lines (0 is left and 2 is right).
      //1 and 3: horizontal lines (1 is the bottom and 3 is the top).
      //In 3D:
      //0 and 3: verical lines (0 is left and 3 is right).
      //1 and 4: horizontal lines (1 is the front and 4 is the back).
      //2 and 5: horizontal lines (2 is the bottom and 5 is the top).
      //----------------------------------------------------------------------

      //@HAF: MEGET VIKTIG: KAN DENNE NUMMERERNGEN SKAPE PROBLEMER i forhold til Dune???
      if ((*simDataPtr).td_.ProblemDimension==2)
      {
	alpha1ExternalSurfaces.u_changeElement_1u(0,  (*simDataPtr).td_.BoundaryCondLeft_alpha1);
	alpha2ExternalSurfaces.u_changeElement_1u(0,  (*simDataPtr).td_.BoundaryCondLeft_alpha2);
	
	alpha1ExternalSurfaces.u_changeElement_1u(2,  (*simDataPtr).td_.BoundaryCondRight_alpha1);
	alpha2ExternalSurfaces.u_changeElement_1u(2,  (*simDataPtr).td_.BoundaryCondRight_alpha2);
	
	alpha1ExternalSurfaces.u_changeElement_1u(1,  (*simDataPtr).td_.BoundaryCondBottom_alpha1);
	alpha2ExternalSurfaces.u_changeElement_1u(1,  (*simDataPtr).td_.BoundaryCondBottom_alpha2);
	
	alpha1ExternalSurfaces.u_changeElement_1u(3,  (*simDataPtr).td_.BoundaryCondTop_alpha1);
	alpha2ExternalSurfaces.u_changeElement_1u(3,  (*simDataPtr).td_.BoundaryCondTop_alpha2);
      }
      else if ((*simDataPtr).td_.ProblemDimension==3)
      {
	cout << "No code in 3D yet!!!" << endl;
      }
    }
    else if (pressureProblemType == 1)
    {
      //Here we treat the problem for PC:
  //----------------------------------------------------------------------
  //NB: Vi maa vite hvilken av sidene nummerne som representerer hva 
  //geometrisk her. MEN DETTE BØR VÆRE presis det samme som i vanlig BC 
  //behandling (se ovenfor for psi).
  //----------------------------------------------------------------------

      if ((*simDataPtr).td_.ProblemDimension==2)
      {
	alpha1ExternalSurfaces.u_changeElement_1u(0,(*simDataPtr).IMPESd_.BoundaryCondLeft_alpha1);
	alpha2ExternalSurfaces.u_changeElement_1u(0,(*simDataPtr).IMPESd_.BoundaryCondLeft_alpha2);
   
	//*******************************************************************************
	alpha1ExternalSurfaces.u_changeElement_1u(2,(*simDataPtr).IMPESd_.BoundaryCondRight_alpha1);
	alpha2ExternalSurfaces.u_changeElement_1u(2,(*simDataPtr).IMPESd_.BoundaryCondRight_alpha2);

	//*******************************************************************************
	alpha1ExternalSurfaces.u_changeElement_1u(1,(*simDataPtr).IMPESd_.BoundaryCondBottom_alpha1);
	alpha2ExternalSurfaces.u_changeElement_1u(1,(*simDataPtr).IMPESd_.BoundaryCondBottom_alpha2);
	
	//*******************************************************************************
	alpha1ExternalSurfaces.u_changeElement_1u(3,(*simDataPtr).IMPESd_.BoundaryCondTop_alpha1);
	alpha2ExternalSurfaces.u_changeElement_1u(3,(*simDataPtr).IMPESd_.BoundaryCondTop_alpha2);
	
      }
      else if ((*simDataPtr).td_.ProblemDimension==3)
      {
	cout << "No code in 3D yet!!!" << endl;
      }
    }
    ml_createConstBoundaryConditionData(0,Igrid,alpha1ExternalSurfaces,alpha1);
    ml_createConstBoundaryConditionData(1,Igrid,alpha2ExternalSurfaces,alpha2);
  }
  
  const Array<double>& K11 = PMpp.getPermeability(0);//@HAF: OK with reference here???
  const Array<double>& K22 = PMpp.getPermeability(1);
  const Array<double>& K12 = PMpp.getPermeability(2);

  //Compute the local inverse matrix and right hand side:
  ml_generateLocalInverseMatrixAndRHS(Igrid, qhatPt, K11, K22, K12, alpha1, alpha2);
}


void SFMpfaFps::generateLocalBoundaryInverseMatrixAndRHS(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const SFPMPhysicalFieldProperties& PMpp, const SFCapPress& cPress, const int& pressureProblemType)//@HAF: OK???
{
  /*
    ------------IMPRTANT INFO.-------------------------------------------
    ---------------------------------------------------------------------
    pressureProblemType == 0 -> Pressure problem for psi.
    pressureProblemType == 1 -> Pressure problem for PC.
    ---------------------------------------------------------------------
  */
  
  //Fix BC-stuff etc.
  const PMTwoPhaseSimulationData* simDataPtr = PMpp.getSimulationDataPtr();
  Array<double> alpha1;
  Array<double> alpha2;
  Array<double> g;
  int BCLL = 2*(*simDataPtr).td_.ProblemDimension;
  Array<double> gExternalSurfaces(BCLL); 
  Array<double> alpha1ExternalSurfaces(BCLL); 
  Array<double> alpha2ExternalSurfaces(BCLL); 

  if (Igrid.boundaryIndex(qhatPt))//@HAF: Is it OK to only allocate alpha1, alpha2 and g in THIS case???
  {
    if (pressureProblemType == 0)
    {
      //Here we treat the problem for psi:
      //creates boundary data for Robin boundary conditions.

      //----------------------------------------------------------------------
      //NB: Vi maa vite hvilken av sidene nummerne som representerer hva 
      //geometrisk her.
      //Below we ASSUME the following to be correct:
      //In 2D:
      //0 and 2: verical lines (0 is left and 2 is right).
      //1 and 3: horizontal lines (1 is the bottom and 3 is the top).
      //In 3D:
      //0 and 3: verical lines (0 is left and 3 is right).
      //1 and 4: horizontal lines (1 is the front and 4 is the back).
      //2 and 5: horizontal lines (2 is the bottom and 5 is the top).
      //----------------------------------------------------------------------

      //@HAF: MEGET VIKTIG: KAN DENNE NUMMERERNGEN SKAPE PROBLEMER i forhold til Dune???
      if ((*simDataPtr).td_.ProblemDimension==2)
      {
	alpha1ExternalSurfaces.u_changeElement_1u(0,  (*simDataPtr).td_.BoundaryCondLeft_alpha1);
	alpha2ExternalSurfaces.u_changeElement_1u(0,  (*simDataPtr).td_.BoundaryCondLeft_alpha2);
	gExternalSurfaces.u_changeElement_1u(0,  (*simDataPtr).td_.BoundaryCondLeft_g);
	
	alpha1ExternalSurfaces.u_changeElement_1u(2,  (*simDataPtr).td_.BoundaryCondRight_alpha1);
	alpha2ExternalSurfaces.u_changeElement_1u(2,  (*simDataPtr).td_.BoundaryCondRight_alpha2);
	gExternalSurfaces.u_changeElement_1u(2,  (*simDataPtr).td_.BoundaryCondRight_g);
	
	alpha1ExternalSurfaces.u_changeElement_1u(1,  (*simDataPtr).td_.BoundaryCondBottom_alpha1);
	alpha2ExternalSurfaces.u_changeElement_1u(1,  (*simDataPtr).td_.BoundaryCondBottom_alpha2);
	gExternalSurfaces.u_changeElement_1u(1,  (*simDataPtr).td_.BoundaryCondBottom_g);
	
	alpha1ExternalSurfaces.u_changeElement_1u(3,  (*simDataPtr).td_.BoundaryCondTop_alpha1);
	alpha2ExternalSurfaces.u_changeElement_1u(3,  (*simDataPtr).td_.BoundaryCondTop_alpha2);
	gExternalSurfaces.u_changeElement_1u(3,  (*simDataPtr).td_.BoundaryCondTop_g);
      }
      else if ((*simDataPtr).td_.ProblemDimension==3)
      {
	cout << "No code in 3D yet!!!" << endl;
      }
    }
    else if (pressureProblemType == 1)
    {
      //Here we treat the problem for PC:
      double TOL = 0.0000001;
      double gVal = 0.0;
  //----------------------------------------------------------------------
  //NB: Vi maa vite hvilken av sidene nummerne som representerer hva 
  //geometrisk her. MEN DETTE BØR VÆRE presis det samme som i vanlig BC 
  //behandling (se ovenfor for psi).
  //----------------------------------------------------------------------

      if ((*simDataPtr).td_.ProblemDimension==2)
      {
	alpha1ExternalSurfaces.u_changeElement_1u(0,(*simDataPtr).IMPESd_.BoundaryCondLeft_alpha1);
	alpha2ExternalSurfaces.u_changeElement_1u(0,(*simDataPtr).IMPESd_.BoundaryCondLeft_alpha2);
	//Check if we have Neuman or Dirichlet BCs:
	if ((*simDataPtr).IMPESd_.BoundaryCondLeft_alpha1 < TOL)
	{//We have a Neuman BC
	  gVal = 0.0;//NB: Only treating zero flux case yet. No "CCI".
	}
	else
	{//We have a Dirichlet BC
	  double sBC = (*simDataPtr).IMPESd_.BoundaryCondLeft_g;
	  int rT = (*simDataPtr).IMPESd_.RockTypeLeft;
	  gVal = cPress.getCapPressData(0, sBC, rT);
	}
	gExternalSurfaces.u_changeElement_1u(0,gVal);
	
	//*******************************************************************************
	alpha1ExternalSurfaces.u_changeElement_1u(2,(*simDataPtr).IMPESd_.BoundaryCondRight_alpha1);
	alpha2ExternalSurfaces.u_changeElement_1u(2,(*simDataPtr).IMPESd_.BoundaryCondRight_alpha2);
	//Check if we have Neuman or Dirichlet BCs:
	if ((*simDataPtr).IMPESd_.BoundaryCondRight_alpha1 < TOL)
	{//We have a Neuman BC
	  gVal = 0.0;//NB: Only treating zero flux case yet. No "CCI".
	}
	else
	{//We have a Dirichlet BC
	  double sBC = (*simDataPtr).IMPESd_.BoundaryCondRight_g;
	  int rT = (*simDataPtr).IMPESd_.RockTypeRight;
	  gVal = cPress.getCapPressData(0, sBC, rT);
	}
	gExternalSurfaces.u_changeElement_1u(2,gVal);

	//*******************************************************************************
	alpha1ExternalSurfaces.u_changeElement_1u(1,(*simDataPtr).IMPESd_.BoundaryCondBottom_alpha1);
	alpha2ExternalSurfaces.u_changeElement_1u(1,(*simDataPtr).IMPESd_.BoundaryCondBottom_alpha2);
	//Check if we have Neuman or Dirichlet BCs:
	if ((*simDataPtr).IMPESd_.BoundaryCondBottom_alpha1 < TOL)
	{//We have a Neuman BC
	  gVal = 0.0;//NB: Only treating zero flux case yet. No "CCI".
	}
	else
	{//We have a Dirichlet BC
	  double sBC = (*simDataPtr).IMPESd_.BoundaryCondBottom_g;
	  int rT = (*simDataPtr).IMPESd_.RockTypeBottom;
	  gVal = cPress.getCapPressData(0, sBC, rT);
	}
	gExternalSurfaces.u_changeElement_1u(1,gVal);
	
	//*******************************************************************************
	alpha1ExternalSurfaces.u_changeElement_1u(3,(*simDataPtr).IMPESd_.BoundaryCondTop_alpha1);
	alpha2ExternalSurfaces.u_changeElement_1u(3,(*simDataPtr).IMPESd_.BoundaryCondTop_alpha2);
	//Check if we have Neuman or Dirichlet BCs:
	if ((*simDataPtr).IMPESd_.BoundaryCondTop_alpha1 < TOL)
	{//We have a Neuman BC
	  gVal = 0.0;//NB: Only treating zero flux case yet. No "CCI".
	}
	else
	{//We have a Dirichlet BC
	  double sBC = (*simDataPtr).IMPESd_.BoundaryCondTop_g;
	  int rT = (*simDataPtr).IMPESd_.RockTypeTop;
	  gVal = cPress.getCapPressData(0, sBC, rT);
	}
	gExternalSurfaces.u_changeElement_1u(3,gVal);
	
      }
      else if ((*simDataPtr).td_.ProblemDimension==3)
      {
	cout << "No code in 3D yet!!!" << endl;
      }
    }
    ml_createConstBoundaryConditionData(0,Igrid,alpha1ExternalSurfaces,alpha1);
    ml_createConstBoundaryConditionData(1,Igrid,alpha2ExternalSurfaces,alpha2);
    ml_createConstBoundaryConditionData(2,Igrid,gExternalSurfaces,g);

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //The following is only necessary IF the rocktype is NOT the same
    //through the whole boundary!!! See @HAF: Below
    if (pressureProblemType == 1)
    {
      //Print g
      //cout << "gSize= " << g.u_getSize_0u() << endl;
      //cout << "NumbSubBoundaries= " << Igrid.getNumberOfSubBoundaries() << endl;
      //cout << "Numb of edges at subBoundary 0= " << Igrid.getNumberOfEdgesAtSubBoundary(0) << endl;
      //cout << "Numb of edges at subBoundary 1= " << Igrid.getNumberOfEdgesAtSubBoundary(1) << endl;
      //cout << "Numb of edges at subBoundary 2 " << Igrid.getNumberOfEdgesAtSubBoundary(2) << endl;
      //cout << "Numb of edges at subBoundary 3= " << Igrid.getNumberOfEdgesAtSubBoundary(3) << endl;
      //for (int i=0; i < g.u_getSize_0u(); i++)
      // {
      //	cout << "g(i)= " << g[i] << endl; 
      //}
      //cout << "************************************************" << endl;

      //if (............)//@HAF: BURDE STYRES FRA input fila IMPESInput.dat
      //{
      ml_changeBoundaryConditionDataFor_g(PMpp,cPress,Igrid,g);
      //}
    //Print g
      //cout << "gSize= " << g.u_getSize_0u() << endl;
      //for (int i=0; i < g.u_getSize_0u(); i++)
      //{
      //cout << "g(i)= " << g[i] << endl; 
      //}
      //cout << "************************************************" << endl;
      //exit(0);
    }
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  }
  
  const Array<double>& K11 = PMpp.getPermeability(0);//@HAF: OK with reference here???
  const Array<double>& K22 = PMpp.getPermeability(1);
  const Array<double>& K12 = PMpp.getPermeability(2);

  //Compute the local inverse matrix and right hand side:
  ml_generateLocalBoundaryInverseMatrixAndRHS(Igrid, qhatPt, K11, K22, K12, alpha1, alpha2, g);
}
 


//---------------------- Private Functions ------------------------------------

double SFMpfaFps::ml_matMultMH(const int& kBar, const int& itilde) const
{
  return invMmatrix_.u_getElementFromMatMultOfNonQuadraticMatrices_0u(kBar, hmatrix_, itilde);
}


void SFMpfaFps::ml_generateLocalInverseMatrixAndRHS(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const Array<double>& K11, const Array<double>& K22, const Array<double>& K12, const Array< double >& alpha1, const Array< double >& alpha2)
{
//   int qhatMax = Igrid.cellCount(0,0);
  int cellDim = Igrid.cellDimension();
  Pair< int, int > PG(0, 0);
  Pair< int, int > PGnm1((cellDim - 1), 0);
  Pair< int, int > PGn(cellDim, 0);
  int L;
//   int khat;
  int kMax;
  int ihat;
  int itildeMax;
  int zMax;
  int nuFIndex;
  int nBournm1;
  int nBourn;
  int kBarValue;
  Boolean qhatNotAtBoundaryOrPeriodic = true;
  Boolean khatNotAtBoundary = true;
  MeshShape mShapeMSD, mShapeH;
  MeshPoint mpMSD, mpH;
  double deltaM;
  double f1Value;
  double f2Value;
//   double sum;

  //Array< int > kBarTab;
  //int kBarTab[2];//NOTE: The length of this Dune::array is bigger in 3D!!!
  std::vector<int> kBarTab;//NOTE: The length of this Dune::array is bigger in 3D!!!
  int v_shapeH[MAXDIRECTIONSIZE];

  //Find the global index from the entity pointer:
  int qhat = Igrid.getGlobalIndexFromEntityPointer(qhatPt);//@HAF: Kan vi droppe denne???
  PG . u_updateSecond_1u(qhat);
  qhatNotAtBoundaryOrPeriodic = !Igrid.boundaryIndex(qhatPt);
  nBournm1 = Igrid.getNumberOfNeighbours(qhatPt, cellDim-1);
  nBourn = Igrid.getNumberOfNeighbours(qhatPt, 0);
  (kMax = nBournm1);
  (itildeMax = nBourn);
  mShapeMSD.u_setShape_1u(2, nBournm1+1);
  mpMSD.u_setPoint_1u(mShapeMSD, 0);
  invMmatrix_.u_copy_1u(mShapeMSD, 0.0);
  v_shapeH[0] = nBournm1+1;
  v_shapeH[1] = nBourn;
  mShapeH.u_setShape_1u(2, v_shapeH);
  mpH.u_setPoint_1u(mShapeH, 0);
  hmatrix_.u_setupForNonQuadraticMatrices_1u(nBournm1+1, nBourn, 0.0);

  for (int k = 0;(k < kMax); (k ++))
  {
    mpMSD.u_changeElement_1u(0, k);
    mpH.u_changeElement_1u(0, k);
    //khat = m_C.u_gammaF_Omet_0u(ar2, mp2);
    IRISDuneGridInterface<DuneGridType>::EdgePointer khatPt = Igrid.getGlobalEntityPointerFromRelativeIndex<IRISDuneGridInterface<DuneGridType>::EdgePointer>(qhatPt, k);
    khatNotAtBoundary = !Igrid.boundaryIndex(khatPt);
    if (qhatNotAtBoundaryOrPeriodic)
    {
      (zMax = 2);
    }
    else
    {
      (zMax = (khatNotAtBoundary ? 2 : 1));
    }
    for (int z = 0;(z < zMax); (z ++))
    {
      //mp3 . u_changeElement_1u(2, z);
      //ihat = m_C.u_gammaF_Omet_0u(ar3, mp3);
      IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt = Igrid.getGlobalEntityPointerFromRelativeIndex(qhatPt, k, cellDim-1, z);
      ihat = Igrid.getGlobalIndexFromEntityPointer(ihatPt);
      PGn . u_updateSecond_1u(ihat);
      //nuFIndex = m_C.u_nuF_Omet_0u(ihat, ar3, mp3);
      nuFIndex = Igrid.getDecreasedRelativeIndexFromRelativeIndex(qhatPt, k, z);
      mpH.u_changeElement_1u(1, nuFIndex);
      //m_C.u_getSLKSetRowRightHandCSSorted_0g(qhat, ihat, kBarTab);
      Igrid.getIntersectionRelativeEdgeIndexSetRightHandCSSorted(qhatPt, ihat, kBarTab);
      //int normalVectorLabel = m_C.u_getTransformSpaceNormalVectorLabel_0u(qhat, ihat, khat);
      int normalVectorLabel = ml_getTransformSpaceNormalVectorLabel(Igrid, qhatPt, ihatPt, khatPt);
      //L = kBarTab . u_getSize_0u();
      L = 2;
      for (int kBar = 0; kBar < L; kBar++)
      {
	//kBarTab . u_getElement_0g(kBar, kBarValue);
	kBarValue = kBarTab[kBar];
	mpMSD.u_changeElement_1u(1, kBarValue);
	//*******DEBUGGING---START****************************************
	//cout << "qhat= " << qhat << " ihat= " << ihat << " ihatCheck= " << Igrid.getGlobalIndexFromEntityPointer(Igrid.getGlobalEntityPointerFromRelativeIndex<IRISDuneGridInterface<DuneGridType>::ElementPointer>(qhatPt, nuFIndex)) << " nuFIndex= " << nuFIndex << endl;
	//cout << "qhat= " << qhat << " ihat= " << ihat << " k= " << k << " kBarValue= " << kBarValue << endl;
	//cout << " qhat= " << qhat << " ihat= " << ihat << " khat= " << Igrid.getGlobalIndexFromEntityPointer(khatPt) << " normalVectorLabel= " << normalVectorLabel << endl;
	//cout << " khat= " <<  Igrid.getGlobalIndexFromEntityPointer(khatPt) << " with coords.= (" << Igrid.getCentroid(khatPt).u_getElement_0u(0) << " , " << Igrid.getCentroid(khatPt).u_getElement_0u(1) << ")" << endl;
	//cout << "Coords. for ihat= (" << Igrid.getCentroid(ihatPt).u_getElement_0u(0) << " , " << Igrid.getCentroid(ihatPt).u_getElement_0u(1) << ")" << endl;
	//*******DEBUGGING---END******************************************
	if ((qhatNotAtBoundaryOrPeriodic || khatNotAtBoundary))
	{
	  (deltaM = ((k == kBarValue) ? invMmatrix_.u_getData_0u(mpMSD) : 0.0));
	  f1Value = ml_FE1(Igrid, qhatPt, ihatPt, ihat, kBar, normalVectorLabel, K11, K22, K12);
	  double fHelp = ml_FE2(Igrid, qhatPt, ihatPt, ihat, normalVectorLabel, K11, K22, K12);
	//*******DEBUGGING---START****************************************
	  //cout << "mpMSD= (" << mpMSD[0] << "," << mpMSD[1] << ")" << " func= " << f1Value - fHelp << " MSD= " << (deltaM + (f1Value - fHelp)) << endl;
	  //cout << "mpMSD= (" << mpMSD[0] << "," << mpMSD[1] << ")" << " deltaM= " << deltaM << " f1Value= " << f1Value <<  " -fHelp= " << -fHelp << " k= " << k << " khat= " << Igrid.getGlobalIndexFromEntityPointer(khatPt) << endl;
	  //cout << "qhat= " << qhat << " ihat= " << ihat << " khat= " <<  Igrid.getGlobalIndexFromEntityPointer(khatPt) << " with coords.= (" << Igrid.getCentroid(khatPt).u_getElement_0u(0) << " , " << Igrid.getCentroid(khatPt).u_getElement_0u(1) << ")" << " kBar= " << kBar << " normalVectorLabel= " << normalVectorLabel << endl;
	//*******DEBUGGING---END******************************************
	  invMmatrix_.u_setData_1u(mpMSD, (deltaM + (f1Value - fHelp)));
	  hmatrix_.u_setData_1u(mpH, hmatrix_.u_getData_0u(mpH) + f1Value);
	}
	else
	{
	  f2Value = ml_FE1B(Igrid, qhatPt, ihatPt, ihat, kBar, khatPt, normalVectorLabel, K11, K22, K12, alpha1, alpha2);
	  double fHelp = ml_FE2B(Igrid, qhatPt, ihatPt, ihat, khatPt, normalVectorLabel, K11, K22, K12, alpha1, alpha2);//OK???
	  invMmatrix_.u_setData_1u(mpMSD, f2Value - fHelp);
	  hmatrix_.u_setData_1u(mpH, hmatrix_.u_getData_0u(mpH) + f2Value);
	  //*******DEBUGGING---START****************************************
	  //if (qhat == 0) cout << "mpH= (" << mpH[0] << "," << mpH[1] << ")" << " f2Value= " << f2Value << " hmatrix_= " << hmatrix_.u_getData_0u(mpH) << " k= " << k << " khat= " << Igrid.getGlobalIndexFromEntityPointer(khatPt) << endl;
	  //*******DEBUGGING---END******************************************
	}
      }
      if ((! (qhatNotAtBoundaryOrPeriodic || khatNotAtBoundary)))
      {
	{
	  double b_0;
	  //alpha1.u_getElement_0g(m_C.u_findBoundaryPosition_0u(PGnm1), b_0);
	  alpha1.u_getElement_0g(Igrid.findBoundaryPosition(khatPt), b_0);
	  hmatrix_.u_setData_1u(mpH, hmatrix_.u_getData_0u(mpH) - b_0);//OK ogsaa naa???Trenger vi ogsaa lignende utvidelser nedenfor???
	  //*******DEBUGGING---START****************************************
	  //if (qhat == 0) cout << "mpH= (" << mpH[0] << "," << mpH[1] << ")" << " findFCPos= " << Igrid.findBoundaryPosition(khatPt) << " b_0= " << b_0 << " hmatrix_= " << hmatrix_.u_getData_0u(mpH) << " k= " << k << " khat= " << Igrid.getGlobalIndexFromEntityPointer(khatPt) << endl;
	  //if (qhat == 0) cout << "alpha1= " << alpha1[0] << "," <<alpha1[1] << "," << alpha1[2] << "," << alpha1[3] << "," << alpha1[4] << "," << alpha1[5] << "," << alpha1[6] << "," << alpha1[7] << endl;
	  //*******DEBUGGING---END******************************************
	}
      }
      //******************************************
      //Must treat the outermost position in the matrix etc.
      mpMSD.u_changeElement_1u(1, kMax);
      double fHelp = (qhatNotAtBoundaryOrPeriodic || khatNotAtBoundary) ? ml_FE2(Igrid, qhatPt, ihatPt, ihat, normalVectorLabel, K11, K22, K12) : ml_FE2B(Igrid, qhatPt, ihatPt, ihat, khatPt, normalVectorLabel, K11, K22, K12, alpha1, alpha2);//BC-OK???
      invMmatrix_.u_setData_1u(mpMSD, invMmatrix_.u_getData_0u(mpMSD) + fHelp);
      hmatrix_.u_setData_1u(mpH, hmatrix_.u_getData_0u(mpH) - fHelp);
      //******************************************
    }
  }
  //exit(0);
  //*********************************************************
  //Here we treat the last equation over the dual grid:
  //Note the special treatment below where the dual grid 
  //is identical to the boundary: AND also the special 
  //treatment when the vertex
  //is at a boundary which has a Dirichlet boundary condition :-)
  //NOTE: This presently ONLY works for constant
  //Dirichlet boundary conditions along the boundary !!!
  //*********************************************************
  mpMSD.u_changeElement_1u(0, kMax);
  mpH.u_changeElement_1u(0, kMax);
  //CSnShapeGIndex GIB(1,0), GIE(2,0);
  bool BCk;
  int signDG = 1;
  bool DirichletBC = false;
  double alpha1Val;
  for (int i=0; i < nBourn; i++)
  {
    mpH.u_changeElement_1u(1, i);
    //mp1.u_changeElement_1u(1,i);
    //ihat = m_C.u_gammaF_Omet_0u(ar1,mp1);
    IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt = Igrid.getGlobalEntityPointerFromRelativeIndex<IRISDuneGridInterface<DuneGridType>::ElementPointer>(qhatPt, i);
    ihat = Igrid.getGlobalIndexFromEntityPointer(ihatPt);
    //m_C.u_getSLKSetRowRightHandCSSorted_0g(qhat, ihat, kBarTab);
    Igrid.getIntersectionRelativeEdgeIndexSetRightHandCSSorted(qhatPt, ihat, kBarTab);
    for (int p=0; p < 2; p++)
    {
      //mp2.u_changeElement_1u(1, kBarTab[p]);
      //khat = m_C.u_gammaF_Omet_0u(ar2,mp2);
      IRISDuneGridInterface<DuneGridType>::EdgePointer khatPt = Igrid.getGlobalEntityPointerFromRelativeIndex<IRISDuneGridInterface<DuneGridType>::EdgePointer>(qhatPt, kBarTab[p]);
      //GIB.u_updateSecond_1u(khat);
      //BCk = m_C.u_boundaryIndex_0u(GIB);
      BCk = Igrid.boundaryIndex(khatPt);
      if (BCk)
      {//Here we check for the Dirichlet boundary condition:
	//int pos = m_C.u_findBoundaryPosition_0u(GIB);
	int pos = Igrid.findBoundaryPosition(khatPt);
	double TOL = 0.000001;
	if (fabs(alpha2[pos]) < TOL) 
	{
	  alpha1Val = alpha1[pos];
	  DirichletBC = true;
	}
	if (DirichletBC) break;
      }
      //int normalVectorLabel = m_C.u_getTransformSpaceNormalVectorLabel_0u(qhat, ihat, khat);
      int normalVectorLabel = ml_getTransformSpaceNormalVectorLabel(Igrid, qhatPt, ihatPt, khatPt);
      double fHelp = signDG*ml_FE4(Igrid, qhatPt, ihatPt, ihat, normalVectorLabel, K11, K22, K12);
      double fHelp2 = (BCk) ? -ml_FE2(Igrid, qhatPt, ihatPt, ihat, normalVectorLabel, K11, K22, K12) : 0.0;//BC-OK???Fortegn OK??
      for (int j=0; j < 2; j++)
      {
	mpMSD.u_changeElement_1u(1, kBarTab[j]); //OK???
	f1Value = signDG*ml_FE3(Igrid, qhatPt, ihatPt, ihat, j, normalVectorLabel, K11, K22, K12);
	double f1Value2 = (BCk) ? -ml_FE1(Igrid, qhatPt, ihatPt, ihat, j, normalVectorLabel, K11, K22, K12) : 0.0;//BC-OK???Fortegn OK??
	invMmatrix_.u_setData_1u(mpMSD,invMmatrix_.u_getData_0u(mpMSD) + ((f1Value + f1Value2) - (fHelp + fHelp2)));
	hmatrix_.u_setData_1u(mpH, hmatrix_.u_getData_0u(mpH) + (f1Value + f1Value2));
      }
      mpMSD.u_changeElement_1u(1, kMax);
      invMmatrix_.u_setData_1u(mpMSD,invMmatrix_.u_getData_0u(mpMSD) + (fHelp + fHelp2));
      hmatrix_.u_setData_1u(mpH, hmatrix_.u_getData_0u(mpH) - (fHelp + fHelp2));
    }
    if (DirichletBC) break;
  }

  if (DirichletBC)
  {//Just to be safe (due to the above implementation) we perform
    //a "clean up" before we state the equation corresponding
    //to the Dirichlet boundary condition:
    for (int j=0; j < (kMax+1); j++)
    {
      mpMSD.u_changeElement_1u(1, j);
      invMmatrix_.u_setData_1u(mpMSD, 0.0);
    }
    for (int i=0; i < nBourn; i++)
    {
      mpH.u_changeElement_1u(1, i);
      hmatrix_.u_setData_1u(mpH, 0.0);
    }

    //We then implement the Dirichlet boundary condition:
    mpMSD.u_changeElement_1u(1, kMax);
    invMmatrix_.u_setData_1u(mpMSD, alpha1Val*1.0);
  }

  //Equation over the dual grid: FINISHED
  //******************************************************


  bool matWrite = false;
  //matWrite = true;//NB: DEBUGGING!!!
  if (matWrite)
  {
    ((cout << "------------------------------------------------------------------------------") << endl);
    ((cout << "------------------------------------------------------------------------------") << endl);
    ((cout << "generateFormFunctions: 1: Matrisen MSD foer inverteringa: qhat= ") << qhat);
    (cout << endl);
    int tlength;
    invMmatrix_.u_getShape_0u() . u_getElement_0g(0, tlength);
    for (int ti = 0;(ti < tlength); (ti ++))
    {
      mpMSD.u_changeElement_1u(0, ti);
      for (int tj = 0;(tj < tlength); (tj ++))
      {
	{
	  MeshPoint p_0 = mpMSD;
	  p_0 . u_changeElement_1u(1, tj);
	  ((cout << invMmatrix_.u_getData_0u(p_0)) << " ");
	}
      }
      (cout << endl);
    }
    ((cout << "------------------------------------------------------------------------------") << endl);
    ((cout << "------------------------------------------------------------------------------") << endl);
  }

  invMmatrix_.u_invert_1u();

  if (matWrite)
  {
    (cout << "generateFormFunctions: 1: Matrisen MSD etter inverteringa:");
    (cout << endl);
    int tlength;
    {
      invMmatrix_.u_getShape_0u() . u_getElement_0g(0, tlength);
      for (int ti = 0;(ti < tlength); (ti ++))
      {
	mpMSD.u_changeElement_1u(0, ti);
	for (int tj = 0;(tj < tlength); (tj ++))
	{
	  {
	    MeshPoint g_1 = mpMSD;
	    g_1 . u_changeElement_1u(1, tj);
	    ((cout << invMmatrix_.u_getData_0u(g_1)) << " ");
	  }
	}
	(cout << endl);
      }
    }
    ((cout << "------------------------------------------------------------------------------") << endl);
    ((cout << "------------------------------------------------------------------------------") << endl);
    (cout << "generateFormFunctions: 2: Matrisen H:");
    (cout << endl);
    hmatrix_.u_getShape_0u().u_getElement_0g(0, tlength);
    {
      for (int ti = 0;(ti < tlength); (ti ++))
      {
	mpH.u_changeElement_1u(0, ti);
	int ttlength;
	hmatrix_.u_getShape_0u().u_getElement_0g(1, ttlength);
	for (int tj = 0;(tj < ttlength); (tj ++))
	{
	  mpH.u_changeElement_1u(1, tj);
	  cout << hmatrix_.u_getData_0u(mpH) << " ";
	}
	(cout << endl);
      }
    }
    ((cout << "------------------------------------------------------------------------------") << endl);
    ((cout << "------------------------------------------------------------------------------") << endl);
  }
}




void SFMpfaFps::ml_generateLocalBoundaryInverseMatrixAndRHS(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const Array<double>& K11, const Array<double>& K22, const Array<double>& K12, const Array<double>& alpha1, const Array<double>& alpha2, const Array<double>& g)
{
//   int qhatMax = Igrid.cellCount(0,0);
  int cellDim = Igrid.cellDimension();
  Pair< int, int > PG(0, 0);
  Pair< int, int > PGnm1((cellDim - 1), 0);
  Pair< int, int > PGn(cellDim, 0);
  int L;
  int khat;
  int kMax;
  int ihat;
  int itildeMax;
  int zMax;
  int nuFIndex;
  int nBournm1;
  int nBourn;
  int kBarValue;
  bool khatNotAtBoundary = true;
  MeshShape mShapeMSD, mShapeH;
  MeshPoint mpMSD, mpH;
  double deltaM;
  double f1Value;
  double f2Value;
//   double sum;
  //Array< Array< double > > H;
  //Array< double > HInner;

  //Array< int > kBarTab;
  //int kBarTab[2];//NOTE: The length of this Dune::array is bigger in 3D!!!
  std::vector<int> kBarTab;//NOTE: The length of this Dune::array is bigger in 3D!!!
  int v_shapeH[MAXDIRECTIONSIZE];

  //Find the global index from the entity pointer:
  int qhat = Igrid.getGlobalIndexFromEntityPointer(qhatPt);//@HAF: Kan vi droppe denne???
  PG . u_updateSecond_1u(qhat);
  bool qhatAtBoundary = Igrid.boundaryIndex(qhatPt);
  if (qhatAtBoundary)
  {
    nBournm1 = Igrid.getNumberOfNeighbours(qhatPt, cellDim-1);
    nBourn = Igrid.getNumberOfNeighbours(qhatPt, 0);
    (kMax = nBournm1);
    (itildeMax = nBourn);
    mShapeMSD.u_setShape_1u(2, nBournm1+1);//HAF-Des2006
    mpMSD.u_setPoint_1u(mShapeMSD, 0);
    invMmatrix_.u_copy_1u(mShapeMSD, 0.0);
    //HInner.u_copy_1u(nBourn, 0.0);
    //H.u_copy_1u(nBournm1+1, HInner);//HAF-Des2006
    v_shapeH[0] = nBournm1+1;
    v_shapeH[1] = 1;
    mShapeH.u_setShape_1u(2, v_shapeH);//HAF-Des2006
    mpH.u_setPoint_1u(mShapeH, 0);
    hmatrix_.u_setupForNonQuadraticMatrices_1u(nBournm1+1, 1, 0.0);

    for (int k = 0;(k < kMax); (k ++))
    {
      mpMSD.u_changeElement_1u(0, k);
      mpH.u_changeElement_1u(0, k);
      //(khat = m_C . u_gammaF_Omet_0u(ar2, mp2));
      IRISDuneGridInterface<DuneGridType>::EdgePointer khatPt = Igrid.getGlobalEntityPointerFromRelativeIndex<IRISDuneGridInterface<DuneGridType>::EdgePointer>(qhatPt, k);
      PGnm1 . u_updateSecond_1u(khat);
      khatNotAtBoundary = !Igrid.boundaryIndex(khatPt);
      (zMax = (khatNotAtBoundary ? 2 : 1));
      for (int z = 0;(z < zMax); (z ++))
      {
	//mp3 . u_changeElement_1u(2, z);
	//ihat = m_C . u_gammaF_Omet_0u(ar3, mp3);
	IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt = Igrid.getGlobalEntityPointerFromRelativeIndex(qhatPt, k, cellDim-1, z);
	ihat = Igrid.getGlobalIndexFromEntityPointer(ihatPt);
	PGn . u_updateSecond_1u(ihat);
	//nuFIndex = m_C . u_nuF_Omet_0u(ihat, ar3, mp3);
	nuFIndex = Igrid.getDecreasedRelativeIndexFromRelativeIndex(qhatPt, k, z);
	//m_C.u_getSLKSetRowRightHandCSSorted_0g(qhat, ihat, kBarTab);
	Igrid.getIntersectionRelativeEdgeIndexSetRightHandCSSorted(qhatPt, ihat, kBarTab);
	//int normalVectorLabel = m_C.u_getTransformSpaceNormalVectorLabel_0u(qhat, ihat, khat);
	int normalVectorLabel = ml_getTransformSpaceNormalVectorLabel(Igrid, qhatPt, ihatPt, khatPt);
	//(L = kBarTab . u_getSize_0u());
	L = 2;
	for (int kBar = 0; kBar < L; kBar++)
	{
	  //kBarTab . u_getElement_0g(kBar, kBarValue);
	  kBarValue = kBarTab[kBar];
	  mpMSD.u_changeElement_1u(1, kBarValue);
	  if ((khatNotAtBoundary))
	  {
	    (deltaM = ((k == kBarValue) ? invMmatrix_.u_getData_0u(mpMSD) : 0.0));
	    f1Value = ml_FE1(Igrid, qhatPt, ihatPt, ihat, kBar, normalVectorLabel, K11, K22, K12);
	    double fHelp = ml_FE2(Igrid, qhatPt, ihatPt, ihat, normalVectorLabel, K11, K22, K12);
	    invMmatrix_.u_setData_1u(mpMSD, (deltaM + (f1Value - fHelp)));
	  }
	  else
	  {
	    f2Value = ml_FE1B(Igrid, qhatPt, ihatPt, ihat, kBar, khatPt, normalVectorLabel, K11, K22, K12, alpha1, alpha2);
	    double fHelp = ml_FE2B(Igrid, qhatPt, ihatPt, ihat, khatPt, normalVectorLabel, K11, K22, K12, alpha1, alpha2);//OK???
	    invMmatrix_.u_setData_1u(mpMSD, f2Value - fHelp);
	  }
	}
	if (! khatNotAtBoundary)
	{
	  {
	    double b_0;
	    //g . u_getElement_0g(m_C . u_findBoundaryPosition_0u(PGnm1), b_0);
	    g.u_getElement_0g(Igrid.findBoundaryPosition(khatPt), b_0);
	    mpH.u_changeElement_1u(1,0);//NB!!!
	    hmatrix_.u_setData_1u(mpH, b_0);
	  }
	}
	//******************************************
	//Must treat the outermost position in the matrix etc.
	mpMSD.u_changeElement_1u(1, kMax);
	double fHelp = (khatNotAtBoundary) ? ml_FE2(Igrid, qhatPt, ihatPt, ihat, normalVectorLabel, K11, K22, K12) : ml_FE2B(Igrid, qhatPt, ihatPt, ihat, khatPt, normalVectorLabel, K11, K22, K12, alpha1, alpha2);//BC-OK???
	invMmatrix_.u_setData_1u(mpMSD, invMmatrix_.u_getData_0u(mpMSD) + fHelp);
	//******************************************
      }
    }
    
    //*********************************************************
    //Here we treat the last equation over the dual grid:
    //Note the special treatment below where the dual grid is identical to the boundary: AND also the special treatment when the vertex
    //is at a boundary which has a Dirichlet boundary condition :-)
    //NOTE: This presently ONLY works for constant
    //Dirichlet boundary conditions along the boundary !!!
    //*********************************************************
    mpMSD.u_changeElement_1u(0, kMax);
    mpH.u_changeElement_1u(0, kMax);
    //CSnShapeGIndex GIB(1,0);
    bool BCk;
    bool DirichletBC = false;
    double alpha1Val, gVal;
    for (int i=0; i < nBourn; i++)
    {
      //mp1.u_changeElement_1u(1,i);
      //ihat = m_C.u_gammaF_Omet_0u(ar1,mp1);
      IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt = Igrid.getGlobalEntityPointerFromRelativeIndex<IRISDuneGridInterface<DuneGridType>::ElementPointer>(qhatPt, i);
      ihat = Igrid.getGlobalIndexFromEntityPointer(ihatPt);
      //m_C . u_getSLKSetRowRightHandCSSorted_0g(qhat, ihat, kBarTab);
      Igrid.getIntersectionRelativeEdgeIndexSetRightHandCSSorted(qhatPt, ihat, kBarTab);
      for (int p=0; p < 2; p++)
      {
	//mp2.u_changeElement_1u(1, kBarTab[p]);
	//khat = m_C.u_gammaF_Omet_0u(ar2,mp2);
	IRISDuneGridInterface<DuneGridType>::EdgePointer khatPt = Igrid.getGlobalEntityPointerFromRelativeIndex<IRISDuneGridInterface<DuneGridType>::EdgePointer>(qhatPt, kBarTab[p]);
	//GIB.u_updateSecond_1u(khat);
	//BCk = m_C.u_boundaryIndex_0u(GIB);
	BCk = Igrid.boundaryIndex(khatPt);
	if (BCk)
	{//Here we check for the Dirichlet boundary condition:
	  int pos = Igrid.findBoundaryPosition(khatPt);
	  double TOL = 0.000001;
	  if (fabs(alpha2[pos]) < TOL) 
	  {
	    alpha1Val = alpha1[pos];
	    gVal = g[pos];
#ifdef CHECK_AGAINST_ALAYTICAL_SQUARE_QUADRATIC_SOLUTION
	    //Here we must have a function that computes the
	    //appropriate value for g:
	    gVal = VertexBCValueForAnalyticalSquareQuadraticSolution(m_C, qhat);
#endif //CHECK_AGAINST_ALAYTICAL_SQUARE_QUADRATIC_SOLUTION

#ifdef CHECK_AGAINST_ALAYTICAL_MESHLOCKING_SOLUTION
	    //Here we must have a function that computes the
	    //appropriate value for g:
	    gVal = VertexBCValueForAnalyticalMeshLockingSolution(m_C, qhat);
#endif //CHECK_AGAINST_ALAYTICAL_MESHLOCKING_SOLUTION
	    DirichletBC = true;
	  }
	  if (DirichletBC) break;
	}
	//int normalVectorLabel = m_C.u_getTransformSpaceNormalVectorLabel_0u(qhat, ihat, khat);
	int normalVectorLabel = ml_getTransformSpaceNormalVectorLabel(Igrid, qhatPt, ihatPt, khatPt);
	double fHelp = ml_FE4(Igrid, qhatPt, ihatPt, ihat, normalVectorLabel, K11, K22, K12);
	double fHelp2 = (BCk) ? -ml_FE2(Igrid, qhatPt, ihatPt, ihat, normalVectorLabel, K11, K22, K12) : 0.0;//BC-OK???Fortegn OK??
	for (int j=0; j < 2; j++)
	{
	  mpMSD.u_changeElement_1u(1, kBarTab[j]); //OK???
	  f1Value = ml_FE3(Igrid, qhatPt, ihatPt, ihat, j, normalVectorLabel, K11, K22, K12);
	  double f1Value2 = (BCk) ? -ml_FE1(Igrid, qhatPt, ihatPt, ihat, j, normalVectorLabel, K11, K22, K12) : 0.0;//BC-OK???Fortegn OK??
	  invMmatrix_.u_setData_1u(mpMSD,invMmatrix_.u_getData_0u(mpMSD) + ((f1Value + f1Value2) - (fHelp + fHelp2)));
	}
	mpMSD.u_changeElement_1u(1, kMax);
	invMmatrix_.u_setData_1u(mpMSD,invMmatrix_.u_getData_0u(mpMSD) + (fHelp + fHelp2));
      }
      if (DirichletBC) break;
    }
    
    if (DirichletBC)
    {//Just to be safe (due to the above implementation) we perform
      //a "clean up" before we state the equation corresponding
      //to the Dirichlet boundary condition:
      for (int j=0; j < (kMax+1); j++)
      {
	mpMSD.u_changeElement_1u(1, j);
	invMmatrix_.u_setData_1u(mpMSD, 0.0);
      }
      
      //We then implement the Dirichlet boundary condition:
      mpMSD.u_changeElement_1u(1, kMax);
      invMmatrix_.u_setData_1u(mpMSD, alpha1Val*1.0);
      
      mpH.u_changeElement_1u(1, 0);
      hmatrix_.u_setData_1u(mpH, gVal);
    }
    
    //Equation over the dual grid: FINISHED
    //******************************************************
    

    bool matWrite = false;
    if (matWrite)
    {
      ((cout << "------------------------------------------------------------------------------") << endl);
      ((cout << "------------------------------------------------------------------------------") << endl);
      ((cout << "generateFormFunctions: 1: Matrisen MSD foer inverteringa: qhat= ") << qhat);
      (cout << endl);
      int tlength;
      invMmatrix_.u_getShape_0u() . u_getElement_0g(0, tlength);
      for (int ti = 0;(ti < tlength); (ti ++))
      {
	mpMSD.u_changeElement_1u(0, ti);
	for (int tj = 0;(tj < tlength); (tj ++))
        {
	  {
	    MeshPoint p_0 = mpMSD;
	    p_0 . u_changeElement_1u(1, tj);
	    ((cout << invMmatrix_.u_getData_0u(p_0)) << " ");
	  }
	}
	(cout << endl);
      }
      ((cout << "------------------------------------------------------------------------------") << endl);
      ((cout << "------------------------------------------------------------------------------") << endl);
    }
    
    invMmatrix_.u_invert_1u();

    if (matWrite)
    {
      (cout << "generateFormFunctions: 1: Matrisen MSD etter inverteringa:");
      (cout << endl);
      int tlength;
      {
	invMmatrix_.u_getShape_0u() . u_getElement_0g(0, tlength);
	for (int ti = 0;(ti < tlength); (ti ++))
        {
	  mpMSD.u_changeElement_1u(0, ti);
	  for (int tj = 0;(tj < tlength); (tj ++))
          {
	    {
	      MeshPoint g_1 = mpMSD;
	      g_1 . u_changeElement_1u(1, tj);
	      ((cout << invMmatrix_.u_getData_0u(g_1)) << " ");
	    }
	  }
	  (cout << endl);
	}
      }
      ((cout << "------------------------------------------------------------------------------") << endl);
      ((cout << "------------------------------------------------------------------------------") << endl);
      (cout << "generateFormFunctions: 2: Matrisen H:");
      (cout << endl);
      hmatrix_.u_getShape_0u().u_getElement_0g(0, tlength);
      {
	for (int ti = 0;(ti < tlength); (ti ++))
        {
	  mpH.u_changeElement_1u(0, ti);
	  int ttlength;
	  hmatrix_.u_getShape_0u().u_getElement_0g(1, ttlength);
	  for (int tj = 0;(tj < ttlength); (tj ++))
          {
	    mpH.u_changeElement_1u(1, tj);
	    cout << hmatrix_.u_getData_0u(mpH) << " ";
	  }
	  (cout << endl);
	}
      }
      ((cout << "------------------------------------------------------------------------------") << endl);
      ((cout << "------------------------------------------------------------------------------") << endl);
    }
  }
}


//-------------------------------------------------------------------
//-------Functions previously located in TempmodCompleteSophus------- 
//-------------------------------------------------------------------

double SFMpfaFps::ml_FE1(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt, const int &ihat, const int& j, const int& normalVector, const Array<double>& K11, const Array<double>& K22, const Array<double>& K12) 
{
  //This funtion computes analyticalintegrals of T-tensor components, used
  //in Edwards O-type discr. or "fully-continuous" versions therof.
  //Some info:
  //normalVector = 0: In this case xi is normalvector i.e. (1.0,0.0)
  //normalVector = 1: In this case eta is normalvector i.e. (0.0,1.0)
  //j can ONLY have values 0 or 1.

  int integralType;

  if (normalVector == 0)
  {
    integralType = (j==0) ? 1 : 3; 
  }
  else if (normalVector == 1)
  {
    integralType = (j==0) ? 4 : 2; 
  }

  double ret = ml_computeTransformSpaceIntegralOfTensorComponents(Igrid, false,qhatPt, ihatPt, ihat, K11, K22, K12, integralType);

  //cout << "retHanna= " << ret << endl;

  return ret;
}



double SFMpfaFps::ml_FE1B(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt, const int &ihat, const int& j, const IRISDuneGridInterface<DuneGridType>::EdgePointer khatPt, const int& normalVector, const Array<double>& K11, const Array<double>& K22, const Array<double>& K12, const Array< double > &alpha1, const Array< double > &alpha2)  
{
  //NBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNB
  //This funtion computes analyticalintegrals of T-tensor components, used
  //in Edwards O-type discr. or "fully-continuous" versions therof.
  //Some info:
  //normalVector = 0: In this case xi is normalvector i.e. (1.0,0.0)
  //normalVector = 1: In this case eta is normalvector i.e. (0.0,1.0)
  //j can ONLY have values 0 or 1.
  //NBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNB

  bool BC = false;

  int integralType; 
  double ret, integralValue;
  int pos = Igrid.findBoundaryPosition(khatPt);
  //  cout << "normalVector= " << normalVector << endl;
  if (normalVector == 0)
  {
    if (j == 0)
    {
      integralType = 1; 
      integralValue = ml_computeTransformSpaceIntegralOfTensorComponents(Igrid,BC,qhatPt, ihatPt, ihat, K11, K22, K12, integralType);
      ret = alpha1[pos]*1.0 + alpha2[pos]*integralValue*1.0;
    }
    else if (j == 1)
    {
      integralType = 3; 
      integralValue = ml_computeTransformSpaceIntegralOfTensorComponents(Igrid,BC,qhatPt, ihatPt, ihat, K11, K22, K12, integralType);
      ////      ret = alpha1[pos]*0.0 + alpha2[pos]*integralValue*0.0;
      ret = alpha1[pos]*0.0 + alpha2[pos]*integralValue*1.0;
      //ret = alpha1[pos]*(0.5) + alpha2[pos]*integralValue*1.0;
    }
  }
  else if (normalVector == 1)
  {
    if (j == 0)
    {
      integralType = 4; 
      integralValue = ml_computeTransformSpaceIntegralOfTensorComponents(Igrid,BC,qhatPt, ihatPt, ihat, K11, K22, K12, integralType);
      ////      ret = alpha1[pos]*0.0 + alpha2[pos]*integralValue*0.0;
      ret = alpha1[pos]*0.0 + alpha2[pos]*integralValue*1.0;
      //ret = alpha1[pos]*(0.5) + alpha2[pos]*integralValue*1.0;
    }
    else if (j == 1)
    {
      integralType = 2; 
      integralValue = ml_computeTransformSpaceIntegralOfTensorComponents(Igrid,BC,qhatPt, ihatPt, ihat, K11, K22, K12, integralType);
      ret = alpha1[pos]*1.0 + alpha2[pos]*integralValue*1.0;
    }
  }

  //cout << "alpha1[pos]= " << alpha1[pos] << endl;
  //cout << "alpha2[pos]= " << alpha2[pos] << endl;
  //cout << "j= " << j << endl;

  return ret;
}



double SFMpfaFps::ml_FE2(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt, const int &ihat, const int& normalVector, const Array<double>& K11, const Array<double>& K22, const Array<double>& K12) 
{
  //This funtion computes analyticalintegrals of T-tensor components, used
  //in Edwards O-type discr. or "fully-continuous" versions therof.
  //Some info:
  //normalVector = 0: In this case xi is normalvector i.e. (1.0,0.0)
  //normalVector = 1: In this case eta is normalvector i.e. (0.0,1.0)

  int integralType = (normalVector==0) ? 5 : 6; 

  double ret = ml_computeTransformSpaceIntegralOfTensorComponents(Igrid,false,qhatPt, ihatPt, ihat, K11, K22, K12, integralType);

  return ret;
}




double SFMpfaFps::ml_FE2B(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt, const int &ihat, const IRISDuneGridInterface<DuneGridType>::EdgePointer khatPt, const int& normalVector, const Array<double>& K11, const Array<double>& K22, const Array<double>& K12, const Array< double > &alpha1, const Array< double > &alpha2)  
{
  //This funtion computes analyticalintegrals of T-tensor components, used
  //in Edwards O-type discr. or "fully-continuous" versions therof.
  //Some info:
  //normalVector = 0: In this case xi is normalvector i.e. (1.0,0.0)
  //normalVector = 1: In this case eta is normalvector i.e. (0.0,1.0)

  //?????????????????????????WHAT TO DO HERE???????????????????????????????????????????
  //NB: Foreløpig er denne impl. antagligvis FEIL !!!
  //NBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNB

  bool BC = false;

  int integralType; 
  double ret, integralValue;
  int pos = Igrid.findBoundaryPosition(khatPt);

  if (normalVector == 0)
  {
    integralType = 5;
    integralValue = ml_computeTransformSpaceIntegralOfTensorComponents(Igrid,BC,qhatPt, ihatPt, ihat, K11, K22, K12, integralType);
    //ret = alpha1[pos]*0.5 + alpha2[pos]*integralValue*1.0; //OK???????????????????????????
    ret = alpha2[pos]*integralValue*1.0; //OK???????????????????????????
  }
  else if (normalVector == 1)
  {
    integralType = 6; 
    integralValue = ml_computeTransformSpaceIntegralOfTensorComponents(Igrid,BC,qhatPt, ihatPt, ihat, K11, K22, K12, integralType);
    //ret = alpha1[pos]*0.5 + alpha2[pos]*integralValue*1.0; //OK???????????????????????????
    ret = alpha2[pos]*integralValue*1.0; //OK???????????????????????????
  }

  return ret;
}



double SFMpfaFps::ml_FE3(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt, const int &ihat, const int& j, const int& normalVector, const Array<double>& K11, const Array<double>& K22, const Array<double>& K12) 
{
  //This funtion computes analyticalintegrals of T-tensor components, used
  //in Edwards O-type discr. or "fully-continuous" versions therof.
  //Some info:
  //normalVector = 0: In this case xi is normalvector i.e. (1.0,0.0)
  //normalVector = 1: In this case eta is normalvector i.e. (0.0,1.0)
  //j can ONLY have values 0 or 1.

  int integralType;

  if (normalVector == 0)
  {
    integralType = (j==0) ? 7 : 9; 
  }
  else if (normalVector == 1)
  {
    integralType = (j==0) ? 10 : 8; 
  }

  double ret = ml_computeTransformSpaceIntegralOfTensorComponents(Igrid,false,qhatPt, ihatPt, ihat, K11, K22, K12, integralType);

  //cout << "retHanna= " << ret << endl;

  return ret;
}


double SFMpfaFps::ml_FE4(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt, const int &ihat, const int& normalVector, const Array<double>& K11, const Array<double>& K22, const Array<double>& K12) 
{
  //This funtion computes analyticalintegrals of T-tensor components, used
  //in Edwards O-type discr. or "fully-continuous" versions therof.
  //Some info:
  //normalVector = 0: In this case xi is normalvector i.e. (1.0,0.0)
  //normalVector = 1: In this case eta is normalvector i.e. (0.0,1.0)

  int integralType = (normalVector==0) ? 11 : 12; 

  double ret = ml_computeTransformSpaceIntegralOfTensorComponents(Igrid, false,qhatPt, ihatPt, ihat, K11, K22, K12, integralType);

  return ret;
}


//-------------------------------------------------------------------
//------------Functions previously located in CSnShape--------------- 
//-------------------------------------------------------------------

//-----------------------------------------------------------------------------
//  ml_createConstBoundaryConditionData - 
//
//  DESCRIPTION:
//   Note first that the input parameter: BCConstNr (below) is an integer which encodes 
//the type of Robin boundary condition (BC) data such that 0 corresponds to alpha1, 
//1 corresponds to alpha2 and 2 corresponds to g. The output from this routine is a
//single Dune::array (BCData) which contains the BC data (alpha1 or alpha2 or g) for all
//boundaries (external or internal) of given configuration. Note that the data is constant
//for a given (outer or inner) boundary.  
//
//  RETURNS:
//    No return value.
//
//  NOTE:
//    
//
//  EXAMPLE:
//
//-----------------------------------------------------------------------------


void SFMpfaFps::ml_createConstBoundaryConditionData(
  const int& BCConstNR,        // integer which encodes the type of Robin BC data (see description above).
  const IRISDuneGridInterface<DuneGridType>& Igrid,   // the grid interface
  const Array<double>& valuesExternalSurfaces,   // given (constant) boundary condition data for the external boundaries.
  Array<double>& BCData) const  // The output boundary condition data.
{
  int LL = 0;

  for (int k = 0; k < Igrid.getNumberOfSubBoundaries(); k++)
  {
    LL += Igrid.getNumberOfEdgesAtSubBoundary(k);
  }

  BCData.u_copy_1u(LL,0.0);

  int boundaryStart = 0;
  //treats the external boundaries:

  for (int i = 0; i < Igrid.getNumberOfSubBoundaries(); i++)
  {
    for (int j = 0; j < Igrid.getNumberOfEdgesAtSubBoundary(i); j++)
    {
      BCData.u_changeElement_1u(j + boundaryStart, valuesExternalSurfaces[i]);
    }
    boundaryStart += Igrid.getNumberOfEdgesAtSubBoundary(i);
  }

}


void SFMpfaFps::ml_changeBoundaryConditionDataFor_g(
  const SFPMPhysicalFieldProperties& PMpp,
  const SFCapPress& cPress,
  const IRISDuneGridInterface<DuneGridType>& Igrid,   // the grid interface
  Array<double>& BCData) const  // The output boundary condition data.
{
  //The following is only relevant in the more rare cases where the rocktype is NOT
  //the same at the whole boundary.

  //The following is a very sub-optimal impl., but due to the way DUNE is organized
  //I do not know how to do it otherwise.......???????

  const PMTwoPhaseSimulationData* simDataPtr = PMpp.getSimulationDataPtr();
  for (IRISDuneGridInterface<DuneGridType>::FaceIterator khatIt = Igrid.setEntityPointerToFirst<IRISDuneGridInterface<DuneGridType>::FaceIterator>(0); khatIt!=Igrid.entityPointerIsAtEnd<IRISDuneGridInterface<DuneGridType>::FaceIterator>(0); ++khatIt)
    {
      IRISDuneGridInterface<DuneGridType>::FacePointer khatPt = khatIt;
      bool khatAtBoundary = Igrid.boundaryIndex(khatPt);
      if (khatAtBoundary)
      {
	//Find the coordinates; identify the sub-boundaries and
	//(possibly) replace the values (at the correct position) in BCData.
	RnPoint P = Igrid.getCentroid(khatPt);
	double x = P.u_getElement_0u(0);
	double y = P.u_getElement_0u(1);
	double sBC = -10000.0;
	double lengthX = (*simDataPtr).td_.X_Extent;
	double lengthY = (*simDataPtr).td_.Y_Extent;
	double TOL = 0.000001;
	if (x < TOL)
	{
	  if ((*simDataPtr).IMPESd_.BoundaryCondLeft_alpha1 > TOL)
	  {
	    sBC = (*simDataPtr).IMPESd_.BoundaryCondLeft_g;
	    int rT = ml_findRockTypeAtBoundary(PMpp, Igrid, khatPt);
	    //cout << " x= " << x << " y= " << y << " BoundaryPosLeft= " << Igrid.findBoundaryPosition(khatPt) << endl;
	    BCData.u_changeElement_1u(Igrid.findBoundaryPosition(khatPt), cPress.getCapPressData(0, sBC, rT));
	  }
	}
	else if (x > (lengthX - TOL))
	{
	  if ((*simDataPtr).IMPESd_.BoundaryCondRight_alpha1 > TOL)
	  {
	    sBC = (*simDataPtr).IMPESd_.BoundaryCondRight_g;
	    int rT = ml_findRockTypeAtBoundary(PMpp, Igrid, khatPt);
	    //cout << " x= " << x << " y= " << y << " BoundaryPosRight= " << Igrid.findBoundaryPosition(khatPt) << endl;
	    BCData.u_changeElement_1u(Igrid.findBoundaryPosition(khatPt), cPress.getCapPressData(0, sBC, rT));
	  }
	}
	else if (y < TOL)
	{
	  if ((*simDataPtr).IMPESd_.BoundaryCondBottom_alpha1 > TOL)
	  {
	    sBC = (*simDataPtr).IMPESd_.BoundaryCondBottom_g;
	    int rT = ml_findRockTypeAtBoundary(PMpp, Igrid, khatPt);
	    //cout << " x= " << x << " y= " << y << " BoundaryPosBottom= " << Igrid.findBoundaryPosition(khatPt) << endl;
	    BCData.u_changeElement_1u(Igrid.findBoundaryPosition(khatPt), cPress.getCapPressData(0, sBC, rT));
	  }
	}
	else if (y > (lengthY - TOL))
	{
	  if ((*simDataPtr).IMPESd_.BoundaryCondTop_alpha1 > TOL)
	  {
	    sBC = (*simDataPtr).IMPESd_.BoundaryCondTop_g;
	    int rT = ml_findRockTypeAtBoundary(PMpp, Igrid, khatPt);
	    //cout << " x= " << x << " y= " << y << " BoundaryPosTop= " << Igrid.findBoundaryPosition(khatPt) << endl;
	    BCData.u_changeElement_1u(Igrid.findBoundaryPosition(khatPt), cPress.getCapPressData(0, sBC, rT));
	  }
	}
      }
      
    }
}


int SFMpfaFps::ml_findRockTypeAtBoundary(const SFPMPhysicalFieldProperties& PMpp, const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::FacePointer khatPt) const
{
  //*******************************************************************
  //@HAF: NOTE:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  //This for the layered example on structured grids!!!
  //*******************************************************************
  //RnPoint P = Igrid.getCentroid(khatPt);
  //double y = P.u_getElement_0u(1);

  //const PMTwoPhaseSimulationData* simDataPtr = PMpp.getSimulationDataPtr();
  //double lengthY = (*simDataPtr).td_.Y_Extent;
  //int ret = 0;
  //if ((y < (lengthY - (1.0/3.0)*lengthY)) && (y > (1.0/3.0)*lengthY)) ret = 1;

  //return ret;
  return 1;
}


int SFMpfaFps::ml_getTransformSpaceNormalVectorLabel(const IRISDuneGridInterface<DuneGridType>& Igrid, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt, const IRISDuneGridInterface<DuneGridType>::EdgePointer khatPt) const
{
  //**********************************************************
  //NOTE that we suppose as a very important precondition for
  //using this routine that khat IS in the so called SLK set
  //defined by qhat and ihat (i.e. that khat is an edge which
  //is in the intersection set of edges belonging the
  //qhat and ihat, respectively). If not, the following implementation
  //will be absurd !!!
  //**********************************************************

  int normalVectorLabel; //NB: Can only be 0 or 1.
  //0 means that xi i.e. (1.0,0.0) is the normal vector, 
  //whereas 1 means that eta i.e. (0.0,1.0) is the normal vector.

//   int cellDim = Igrid.cellDimension();

  RnPoint x0 = Igrid.getCentroid<IRISDuneGridInterface<DuneGridType>::VertexPointer>(qhatPt);
  RnPoint x1 = Igrid.getCentroid<IRISDuneGridInterface<DuneGridType>::ElementPointer>(ihatPt);
  RnPoint x2 = Igrid.getCentroid<IRISDuneGridInterface<DuneGridType>::EdgePointer>(khatPt);

  //We now assume a right-handed xi-eta coord. system.
  normalVectorLabel = (x0.isCounterClockwise(x1,x2)) ? 0 : 1;

  return normalVectorLabel;
}


double SFMpfaFps::ml_computeTransformSpaceIntegralOfTensorComponents(const IRISDuneGridInterface<DuneGridType>& Igrid, const bool& BC, const IRISDuneGridInterface<DuneGridType>::VertexPointer qhatPt, const IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt, const int& ihat, const Array<double>& K11, const Array<double>& K22, const Array<double>& K12, const int& integralType) const
{
  //*******************************************************************
  //The parameter integralType is meant to decide which integral
  //or combination of integrals this function is supposed to deliver.
  //We have integralType = 1, 2, 3, .....,12
  //*******************************************************************

  //********************************************************
  //const double ParamADC = 0.999999;
  double ParamADC = 0.999999;
  //double ParamADC = 0.9;
  //Parameter connected to the size of the auxillary dual cell.
  //Note that we should always have 0 <= ParamADC < 1. 
  //ParamADC=0.0 gives the usual dual cell, whereas ParamADC=1
  //makes the auxillary dual cell shrink to a single point.
  //********************************************************

  //********************************************************
  //Parameters and variables etc. connected to the
  //"special elliptic scheme".
  bool SpecialEllipticScheme = true;
  //Parameter connected to the discrete flux at the grid half edges.
  //const double ParamFluxG = 1.0;//NB: Erstattet av m_ParamFluxG !!!
  //Parameter connected to the discrete flux at the auxillary dual cell edges.
  //const double ParamFluxADC = 0.9999995;
  double ParamFluxADC = 0.9999995;
  //double ParamFluxADC = 0.9;
  //Note that we should always have ParamADC <= ParamFluxADC < 1. 
  //********************************************************

  //*****************************************************************************
  //NOTE: The following is the most important flux-quadrature parameter.
  //It should in the range 0 <= m_ParamFluxG < 1
  //Should preferably have been in the input file, but a safe choice
  //is 0.0 as employed below.
  const double m_ParamFluxG = 0.0;
  //*****************************************************************************

  const double ParamFluxGNew = m_ParamFluxG;//@HAF: This is to be used particularly for the xi and eta's, which may multiply the tensor components in some of the integrals below.
  const double ParamFluxADCNew = ParamFluxADC;//@HAF: This is to be used particularly for the xi and eta's, which may multiply the tensor components in some of the integrals below.

  //----------------------NBNBNBNBNBNBNBNBNBNBNBNBNBNBNB----------------------------------------
  //Here our aim is to make an SPD T-tensor (NOTE that we in this case must have m_ParamFluxG=1.0 from the outset).
  //ParamADC = 1.0;
  //ParamFluxADC = 1.0;
  //const double ParamFluxGNew = 0.333333333*m_ParamFluxG;//@HAF: This is to be used particularly for the xi and eta's, which may multiply the tensor components in some of the integrals below. The factor 0.333333333 in this case must be >= 0 and < 1.
  //----------------------NBNBNBNBNBNBNBNBNBNBNBNBNBNBNB----------------------------------------

  int LL;
  int ipos;
  std::vector<IRISDuneGridInterface<DuneGridType>::FacePointer> nm1Neigh;

  double ret = 0.0;
 
  //std::vector<IRISDuneGridInterface<DuneGridType>::FacePointer> actualnm1Neigh(2);
  std::vector<IRISDuneGridInterface<DuneGridType>::FacePointer> actualnm1Neigh;
  actualnm1Neigh.reserve(2); //@HAF: OK???
  RnPoint x0, x1, x2, x3;
//   double x_xi, x_eta, y_xi, y_eta;

  double b, c, d, e, IT11, IT22, IT12, IT21;

  x2 = Igrid.getCentroid<IRISDuneGridInterface<DuneGridType>::VertexPointer>(qhatPt);
  
  x0 = Igrid.getCentroid<IRISDuneGridInterface<DuneGridType>::ElementPointer>(ihatPt);
  Igrid.getNeighbours(ihatPt, nm1Neigh);
  LL = nm1Neigh.size();
  ipos = 0;

  //--------------DEBUGGING----------START----------------------------------------

  //cout << "x2= (" << x2.u_getElement_0u(0) << ", " << x2.u_getElement_0u(1) << ")" << endl;
  //cout << "x0= (" << x0.u_getElement_0u(0) << ", " << x0.u_getElement_0u(1) << ")" << endl;
  //cout << "nm1Neigh[0]= (" << Igrid.getCentroid<IRISDuneGridInterface<DuneGridType>::FacePointer>(nm1Neigh[0]).u_getElement_0u(0) << ", " << Igrid.getCentroid<IRISDuneGridInterface<DuneGridType>::FacePointer>(nm1Neigh[0]).u_getElement_0u(1) << ")" << endl;
  //cout << "nm1Neigh[1]= (" << Igrid.getCentroid<IRISDuneGridInterface<DuneGridType>::FacePointer>(nm1Neigh[1]).u_getElement_0u(0) << ", " << Igrid.getCentroid<IRISDuneGridInterface<DuneGridType>::FacePointer>(nm1Neigh[1]).u_getElement_0u(1) << ")" << endl;
  //cout << "nm1Neigh[2]= (" << Igrid.getCentroid<IRISDuneGridInterface<DuneGridType>::FacePointer>(nm1Neigh[2]).u_getElement_0u(0) << ", " << Igrid.getCentroid<IRISDuneGridInterface<DuneGridType>::FacePointer>(nm1Neigh[2]).u_getElement_0u(1) << ")" << endl;
  //cout << "nm1Neigh[3]= (" << Igrid.getCentroid<IRISDuneGridInterface<DuneGridType>::FacePointer>(nm1Neigh[3]).u_getElement_0u(0) << ", " << Igrid.getCentroid<IRISDuneGridInterface<DuneGridType>::FacePointer>(nm1Neigh[3]).u_getElement_0u(1) << ")" << endl;

  //nm1Neigh.erase(nm1Neigh.begin(),nm1Neigh.end());
  //Igrid.getNeighbours(qhatPt, nm1Neigh);
  //cout << "******************************************************************" << endl;
  //cout << "nm1Neigh[0]= (" << Igrid.getCentroid<IRISDuneGridInterface<DuneGridType>::FacePointer>(nm1Neigh[0]).u_getElement_0u(0) << ", " << Igrid.getCentroid<IRISDuneGridInterface<DuneGridType>::FacePointer>(nm1Neigh[0]).u_getElement_0u(1) << ")" << endl;
  //cout << "nm1Neigh[1]= (" << Igrid.getCentroid<IRISDuneGridInterface<DuneGridType>::FacePointer>(nm1Neigh[1]).u_getElement_0u(0) << ", " << Igrid.getCentroid<IRISDuneGridInterface<DuneGridType>::FacePointer>(nm1Neigh[1]).u_getElement_0u(1) << ")" << endl;
  //cout << "nm1Neigh[2]= (" << Igrid.getCentroid<IRISDuneGridInterface<DuneGridType>::FacePointer>(nm1Neigh[2]).u_getElement_0u(0) << ", " << Igrid.getCentroid<IRISDuneGridInterface<DuneGridType>::FacePointer>(nm1Neigh[2]).u_getElement_0u(1) << ")" << endl;
  //cout << "nm1Neigh[3]= (" << Igrid.getCentroid<IRISDuneGridInterface<DuneGridType>::FacePointer>(nm1Neigh[3]).u_getElement_0u(0) << ", " << Igrid.getCentroid<IRISDuneGridInterface<DuneGridType>::FacePointer>(nm1Neigh[3]).u_getElement_0u(1) << ")" << endl;
  //exit(0);

  //std::vector<IRISDuneGridInterface<DuneGridType>::FacePointer> faceNgbs;
  //Igrid.getNeighbours(ihatPt,faceNgbs);
  //std::cout << "face-neighbours of elm:\n";
  //for (unsigned int i=0; i<faceNgbs.size(); ++i)
  // {
  //  std::cout << " " << Igrid.getGlobalIndexFromEntityPointer(faceNgbs[i]);
  //  //std::cout << "(isNeighboursOfVtx: " << Igrid.isNeighbour(qhatPt,faceNgbs[i]) << ")" << endl;
  //  std::cout << "(isNeighboursOfVtx: ";
  //  if (Igrid.isNeighbour(qhatPt,faceNgbs[i])) 
  //  {
  //    cout << " JADA)" << endl;
  //  }
  //  else
  //  {
  //    cout << " NEIDA)" << endl;
  //  }
  // }
  //std::cout << std::endl;



  //--------------DEBUGGING------END--------------------------------------------

  for (unsigned int i=0; i<nm1Neigh.size(); ++i)
  {
    //std::cout << " " << Igrid.getGlobalIndexFromEntityPointer(nm1Neigh[i]);
    //std::cout << "(isNeighboursOfVtx: ";
    if (Igrid.isNeighbour(qhatPt,nm1Neigh[i])) 
    {
      //cout << " JADA, i= " << i << endl;
      //actualnm1Neigh[ipos] = nm1Neigh[i];
      actualnm1Neigh.push_back(nm1Neigh[i]);
      ipos++;
    }
    //else
    //{
    // cout << " NEIDA)" << endl;
    //} 
  }
  //std::cout << std::endl;

  x1 = Igrid.getCentroid<IRISDuneGridInterface<DuneGridType>::FacePointer>(actualnm1Neigh[0]);
  x3 = Igrid.getCentroid<IRISDuneGridInterface<DuneGridType>::FacePointer>(actualnm1Neigh[1]);

  //**************************************************************
  //**************************************************************
  //We want a counter-clockwise rotation of (x0,x1,x2,x3), such that
  //all local coord. systems are right-handed. This can be
  //obtained in the following way:
  if (!x0.isCounterClockwise(x1,x2))
  {
    //Must change x1 and x3:
    RnPoint help = x1;
    x1 = x3;
    x3 = help;
  }


  //In the following, all integrals are computed analytically:
  int index = ihat;

  if (integralType == 1)
  {
    //We first compute the eta-integral (for xi=1.0) of tensor component T11:
    c = K11[index]*ml_y_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3)*
      ml_y_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3) +
      K22[index]*ml_x_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3)*
      ml_x_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3) -
      2.0*K12[index]*ml_x_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3)*
      ml_y_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3);

    d = ml_determinant_Forxi_eq1ANDeta_LinearPart(x0,x1,x2,x3);
    e = ml_determinant_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3);
    IT11 = (c/d)*log((e+d)/e);
    if (BC) IT11 = c/(d+e);//HAF-Edw sjekk!!!

    ret = IT11;
    if (SpecialEllipticScheme)
    {
      ret = c/((d*m_ParamFluxG)+e);
    }
  }
  else if (integralType == 2)
  {
    //We then compute the xi-integral (for eta=1.0) of tensor component T22:
    c = K11[index]*ml_y_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)*
      ml_y_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3) +
      K22[index]*ml_x_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)*
      ml_x_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3) -
      2.0*K12[index]*ml_x_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)*
      ml_y_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3);
  
    d = ml_determinant_Forxi_ANDeta_eq1_LinearPart(x0,x1,x2,x3);
    e = ml_determinant_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3);
    IT22 = (c/d)*log((e+d)/e);
    if (BC) IT22 = c/(d+e);//HAF-Edw sjekk!!!

    ret = IT22;
    if (SpecialEllipticScheme)
    {
      ret = c/((d*m_ParamFluxG)+e);
    }
  }
  else if (integralType == 3)
  {
    //We then compute the eta-integral (for xi=1.0) of tensor component T12:
    b = K12[index]*(ml_x_xi_Forxi_eq1ANDeta_LinearPart(x0,x1,x2,x3)*
		    ml_y_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3) +
		    ml_x_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3)*
		    ml_y_xi_Forxi_eq1ANDeta_LinearPart(x0,x1,x2,x3)) -
      (K11[index]*(ml_y_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3)*
		   ml_y_xi_Forxi_eq1ANDeta_LinearPart(x0,x1,x2,x3)) +
       K22[index]*(ml_x_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3)*
		   ml_x_xi_Forxi_eq1ANDeta_LinearPart(x0,x1,x2,x3)));
  
    c = K12[index]*(ml_x_xi_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3)*
		    ml_y_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3) +
		    ml_x_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3)*
		    ml_y_xi_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3)) -
      (K11[index]*(ml_y_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3)*
		   ml_y_xi_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3)) +
       K22[index]*(ml_x_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3)*
		   ml_x_xi_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3)));

    d = ml_determinant_Forxi_eq1ANDeta_LinearPart(x0,x1,x2,x3);
    e = ml_determinant_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3);
    IT12 = (b/d) + ((c*d - b*e)/(d*d))*log((e+d)/e);
    if (BC) IT12 = (b+c)/(d+e);//HAF-Edw sjekk!!!

    ret = IT12;
    if (SpecialEllipticScheme)
    {
      ret = ((b*m_ParamFluxG)+c)/((d*m_ParamFluxG)+e);
    }
  }
  else if (integralType == 4)
  {
    //We then compute the xi-integral (for eta=1.0) of tensor component T12:
    b = K12[index]*(ml_x_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)*
		    ml_y_eta_Forxi_ANDeta_eq1_LinearPart(x0,x1,x2,x3) +
		    ml_x_eta_Forxi_ANDeta_eq1_LinearPart(x0,x1,x2,x3)*
		    ml_y_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)) -
      (K11[index]*(ml_y_eta_Forxi_ANDeta_eq1_LinearPart(x0,x1,x2,x3)*
		   ml_y_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)) +
       K22[index]*(ml_x_eta_Forxi_ANDeta_eq1_LinearPart(x0,x1,x2,x3)*
		   ml_x_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)));
      
    c = K12[index]*(ml_x_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)*
		    ml_y_eta_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3) +
		    ml_x_eta_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)*
		    ml_y_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)) -
      (K11[index]*(ml_y_eta_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)*
		   ml_y_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)) +
       K22[index]*(ml_x_eta_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)*
		   ml_x_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)));
  
    d = ml_determinant_Forxi_ANDeta_eq1_LinearPart(x0,x1,x2,x3);
    e = ml_determinant_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3);
    IT21 = (b/d) + ((c*d - b*e)/(d*d))*log((e+d)/e);
    if (BC) IT21 = (b+c)/(d+e);//HAF-Edw sjekk!!!

    ret = IT21;
    if (SpecialEllipticScheme)
    {
      ret = ((b*m_ParamFluxG)+c)/((d*m_ParamFluxG)+e);
    }
  }
  else if (integralType == 5)
  {
    //We then compute the eta-integral (for xi=1.0) of T12 + T11*eta
    //We first compute the eta-integral (for xi=1.0) of tensor component T12:
    b = K12[index]*(ml_x_xi_Forxi_eq1ANDeta_LinearPart(x0,x1,x2,x3)*
		    ml_y_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3) +
		    ml_x_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3)*
		    ml_y_xi_Forxi_eq1ANDeta_LinearPart(x0,x1,x2,x3)) -
      (K11[index]*(ml_y_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3)*
		   ml_y_xi_Forxi_eq1ANDeta_LinearPart(x0,x1,x2,x3)) +
       K22[index]*(ml_x_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3)*
		   ml_x_xi_Forxi_eq1ANDeta_LinearPart(x0,x1,x2,x3)));
  
    c = K12[index]*(ml_x_xi_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3)*
		    ml_y_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3) +
		    ml_x_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3)*
		    ml_y_xi_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3)) -
      (K11[index]*(ml_y_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3)*
		   ml_y_xi_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3)) +
       K22[index]*(ml_x_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3)*
		   ml_x_xi_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3)));

    d = ml_determinant_Forxi_eq1ANDeta_LinearPart(x0,x1,x2,x3);
    e = ml_determinant_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3);
    IT12 = (b/d) + ((c*d - b*e)/(d*d))*log((e+d)/e);
    if (SpecialEllipticScheme)
    {
      IT12 = ((b*m_ParamFluxG)+c)/((d*m_ParamFluxG)+e);
    }

    //We then compute the eta-integral (for xi=1.0) of T11*eta
    b = K11[index]*ml_y_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3)*
      ml_y_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3) +
      K22[index]*ml_x_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3)*
      ml_x_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3) -
      2.0*K12[index]*ml_x_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3)*
      ml_y_eta_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3);

    d = ml_determinant_Forxi_eq1ANDeta_LinearPart(x0,x1,x2,x3);
    e = ml_determinant_Forxi_eq1ANDeta_ConstantPart(x0,x1,x2,x3);
    IT11 = (b/d) - ((b*e)/(d*d))*log((e+d)/e);
    if (SpecialEllipticScheme)
    {
      //IT11 = (b*m_ParamFluxG)/((d*m_ParamFluxG)+e);
      IT11 = (b*ParamFluxGNew)/((d*m_ParamFluxG)+e);//@HAF: 30.10.2008
    }

    ret = IT12 + IT11;
  }
  else if (integralType == 6)
  {
    //We then compute the xi-integral (for eta=1.0) of T12 + T22*xi
    //We first compute the xi-integral (for eta=1.0) of tensor component T12:
    b = K12[index]*(ml_x_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)*
		    ml_y_eta_Forxi_ANDeta_eq1_LinearPart(x0,x1,x2,x3) +
		    ml_x_eta_Forxi_ANDeta_eq1_LinearPart(x0,x1,x2,x3)*
		    ml_y_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)) -
      (K11[index]*(ml_y_eta_Forxi_ANDeta_eq1_LinearPart(x0,x1,x2,x3)*
		   ml_y_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)) +
       K22[index]*(ml_x_eta_Forxi_ANDeta_eq1_LinearPart(x0,x1,x2,x3)*
		   ml_x_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)));
      
    c = K12[index]*(ml_x_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)*
		    ml_y_eta_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3) +
		    ml_x_eta_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)*
		    ml_y_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)) -
      (K11[index]*(ml_y_eta_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)*
		   ml_y_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)) +
       K22[index]*(ml_x_eta_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)*
		   ml_x_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)));
  
    d = ml_determinant_Forxi_ANDeta_eq1_LinearPart(x0,x1,x2,x3);
    e = ml_determinant_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3);
    IT21 = (b/d) + ((c*d - b*e)/(d*d))*log((e+d)/e);
    if (SpecialEllipticScheme)
    {
      IT21 = ((b*m_ParamFluxG)+c)/((d*m_ParamFluxG)+e);
    }

    //We then compute the xi-integral (for eta=1.0) of tensor component T22:
    b = K11[index]*ml_y_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)*
      ml_y_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3) +
      K22[index]*ml_x_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)*
      ml_x_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3) -
      2.0*K12[index]*ml_x_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3)*
      ml_y_xi_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3);
  
    d = ml_determinant_Forxi_ANDeta_eq1_LinearPart(x0,x1,x2,x3);
    e = ml_determinant_Forxi_ANDeta_eq1_ConstantPart(x0,x1,x2,x3);
    IT22 = (b/d) - ((b*e)/(d*d))*log((e+d)/e);
    if (SpecialEllipticScheme)
    {
      //IT22 = (b*m_ParamFluxG)/((d*m_ParamFluxG)+e);
      IT22 = (b*ParamFluxGNew)/((d*m_ParamFluxG)+e);//@HAF: 30.10.2008
    }

    ret = IT21 + IT22;
  }
  else if (integralType == 7)
  {
    //We first compute the eta-integral (for xi=PARAM) of tensor component T11:
    c = K11[index]*ml_y_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC)*
      ml_y_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC) +
      K22[index]*ml_x_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC)*
      ml_x_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC) -
      2.0*K12[index]*ml_x_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC)*
      ml_y_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC);

    d = ml_determinant_Forxi_eqPARAMANDeta_LinearPart(x0,x1,x2,x3,ParamADC);
    e = ml_determinant_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC);
    IT11 = (c/d)*log((e+d)/(e+d*ParamADC));

    ret = IT11;
    if (SpecialEllipticScheme)
    {
      ret = c/((d*ParamFluxADC)+e);
    }
  }
  else if (integralType == 8)
  {
    //We then compute the xi-integral (for eta=PARAM) of tensor component T22:
    c = K11[index]*ml_y_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC)*
      ml_y_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC) +
      K22[index]*ml_x_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC)*
      ml_x_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC) -
      2.0*K12[index]*ml_x_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC)*
      ml_y_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC);
  
    d = ml_determinant_Forxi_ANDeta_eqPARAM_LinearPart(x0,x1,x2,x3,ParamADC);
    e = ml_determinant_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC);
    IT22 = (c/d)*log((e+d)/(e+d*ParamADC));
  
    ret = IT22;
    if (SpecialEllipticScheme)
    {
      ret = c/((d*ParamFluxADC)+e);
    }
  }
  else if (integralType == 9)
  {
    //We then compute the eta-integral (for xi=PARAM) of tensor component T12:
    b = K12[index]*(ml_x_xi_Forxi_eqPARAMANDeta_LinearPart(x0,x1,x2,x3)*
		    ml_y_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC) +
		    ml_x_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC)*
		    ml_y_xi_Forxi_eqPARAMANDeta_LinearPart(x0,x1,x2,x3)) -
      (K11[index]*(ml_y_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC)*
		   ml_y_xi_Forxi_eqPARAMANDeta_LinearPart(x0,x1,x2,x3)) +
       K22[index]*(ml_x_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC)*
		   ml_x_xi_Forxi_eqPARAMANDeta_LinearPart(x0,x1,x2,x3)));
  
    c = K12[index]*(ml_x_xi_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3)*
		    ml_y_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC) +
		    ml_x_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC)*
		    ml_y_xi_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3)) -
      (K11[index]*(ml_y_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC)*
		   ml_y_xi_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3)) +
       K22[index]*(ml_x_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC)*
		   ml_x_xi_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3)));

    d = ml_determinant_Forxi_eqPARAMANDeta_LinearPart(x0,x1,x2,x3,ParamADC);
    e = ml_determinant_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC);
    IT12 = (b/d)*(1.0-ParamADC) + ((c*d - b*e)/(d*d))*log((e+d)/(e+d*ParamADC));
  
    ret = IT12;
    if (SpecialEllipticScheme)
    {
      ret = ((b*ParamFluxADC)+c)/((d*ParamFluxADC)+e);
    }
  }
  else if (integralType == 10)
  {
    //We then compute the xi-integral (for eta=PARAM) of tensor component T12:
    b = K12[index]*(ml_x_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC)*
		    ml_y_eta_Forxi_ANDeta_eqPARAM_LinearPart(x0,x1,x2,x3) +
		    ml_x_eta_Forxi_ANDeta_eqPARAM_LinearPart(x0,x1,x2,x3)*
		    ml_y_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC)) -
      (K11[index]*(ml_y_eta_Forxi_ANDeta_eqPARAM_LinearPart(x0,x1,x2,x3)*
		   ml_y_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC)) +
       K22[index]*(ml_x_eta_Forxi_ANDeta_eqPARAM_LinearPart(x0,x1,x2,x3)*
		   ml_x_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC)));
      
    c = K12[index]*(ml_x_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC)*
		    ml_y_eta_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3) +
		    ml_x_eta_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3)*
		    ml_y_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC)) -
      (K11[index]*(ml_y_eta_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3)*
		   ml_y_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC)) +
       K22[index]*(ml_x_eta_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3)*
		   ml_x_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC)));
  
    d = ml_determinant_Forxi_ANDeta_eqPARAM_LinearPart(x0,x1,x2,x3,ParamADC);
    e = ml_determinant_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC);
    IT21 = (b/d)*(1.0-ParamADC) + ((c*d - b*e)/(d*d))*log((e+d)/(e+d*ParamADC));
  
    ret = IT21;
    if (SpecialEllipticScheme)
    {
      ret = ((b*ParamFluxADC)+c)/((d*ParamFluxADC)+e);
    }
  }
  else if (integralType == 11)
  {
    //We then compute the eta-integral (for xi=PARAM) of T12*xi + T11*eta
    //We first compute the eta-integral (for xi=PARAM) of tensor component T12*xi:
    b = K12[index]*(ml_x_xi_Forxi_eqPARAMANDeta_LinearPart(x0,x1,x2,x3)*
		    ml_y_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC) +
		    ml_x_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC)*
		    ml_y_xi_Forxi_eqPARAMANDeta_LinearPart(x0,x1,x2,x3)) -
      (K11[index]*(ml_y_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC)*
		   ml_y_xi_Forxi_eqPARAMANDeta_LinearPart(x0,x1,x2,x3)) +
       K22[index]*(ml_x_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC)*
		   ml_x_xi_Forxi_eqPARAMANDeta_LinearPart(x0,x1,x2,x3)));
  
    c = K12[index]*(ml_x_xi_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3)*
		    ml_y_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC) +
		    ml_x_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC)*
		    ml_y_xi_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3)) -
      (K11[index]*(ml_y_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC)*
		   ml_y_xi_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3)) +
       K22[index]*(ml_x_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC)*
		   ml_x_xi_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3)));

    d = ml_determinant_Forxi_eqPARAMANDeta_LinearPart(x0,x1,x2,x3,ParamADC);
    e = ml_determinant_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC);
    IT12 = (b/d)*(1.0-ParamADC) + ((c*d - b*e)/(d*d))*log((e+d)/(e+d*ParamADC));
    if (SpecialEllipticScheme)
    {
      IT12 = (((b*ParamFluxADC)+c)/((d*ParamFluxADC)+e))*ParamFluxADCNew;//@HAF: multiplication by ParamADCNew 30.10.2008
    }

    //We then compute the eta-integral (for xi=PARAM) of T11*eta
    b = K11[index]*ml_y_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC)*
      ml_y_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC) +
      K22[index]*ml_x_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC)*
      ml_x_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC) -
      2.0*K12[index]*ml_x_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC)*
      ml_y_eta_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC);

    d = ml_determinant_Forxi_eqPARAMANDeta_LinearPart(x0,x1,x2,x3,ParamADC);
    e = ml_determinant_Forxi_eqPARAMANDeta_ConstantPart(x0,x1,x2,x3,ParamADC);

    IT11 = (b/d)*(1.0-ParamADC) - ((b*e)/(d*d))*log((e+d)/(e+d*ParamADC));
    if (SpecialEllipticScheme)
    {
      //IT11 = (b*ParamFluxADC)/((d*ParamFluxADC)+e);
      IT11 = (b*ParamFluxADCNew)/((d*ParamFluxADC)+e);//@HAF:30.10.2008
    }

    ret = IT12 + IT11;//28.08.07: OK, eller kommentere ut IT12?
  }
  else if (integralType == 12)
  {
    //We then compute the xi-integral (for eta=PARAM) of T12*eta + T22*xi
    //We first compute the xi-integral (for eta=PARAM) of tensor component T12*eta:
    b = K12[index]*(ml_x_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC)*
		    ml_y_eta_Forxi_ANDeta_eqPARAM_LinearPart(x0,x1,x2,x3) +
		    ml_x_eta_Forxi_ANDeta_eqPARAM_LinearPart(x0,x1,x2,x3)*
		    ml_y_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC)) -
      (K11[index]*(ml_y_eta_Forxi_ANDeta_eqPARAM_LinearPart(x0,x1,x2,x3)*
		   ml_y_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC)) +
       K22[index]*(ml_x_eta_Forxi_ANDeta_eqPARAM_LinearPart(x0,x1,x2,x3)*
		   ml_x_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC)));
      
    c = K12[index]*(ml_x_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC)*
		    ml_y_eta_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3) +
		    ml_x_eta_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3)*
		    ml_y_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC)) -
      (K11[index]*(ml_y_eta_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3)*
		   ml_y_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC)) +
       K22[index]*(ml_x_eta_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3)*
		   ml_x_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC)));
  
    d = ml_determinant_Forxi_ANDeta_eqPARAM_LinearPart(x0,x1,x2,x3,ParamADC);
    e = ml_determinant_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC);
    IT21 = (b/d)*(1.0-ParamADC) + ((c*d - b*e)/(d*d))*log((e+d)/(e+d*ParamADC));
    if (SpecialEllipticScheme)
    {
      IT21 = (((b*ParamFluxADC)+c)/((d*ParamFluxADC)+e))*ParamFluxADCNew;//@HAF: multiplication by ParamADCNew 30.10.2008
    }

    //We then compute the xi-integral (for eta=PARAM) of tensor component T22*xi:
    b = K11[index]*ml_y_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC)*
      ml_y_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC) +
      K22[index]*ml_x_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC)*
      ml_x_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC) -
      2.0*K12[index]*ml_x_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC)*
      ml_y_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC);
  
    d = ml_determinant_Forxi_ANDeta_eqPARAM_LinearPart(x0,x1,x2,x3,ParamADC);
    e = ml_determinant_Forxi_ANDeta_eqPARAM_ConstantPart(x0,x1,x2,x3,ParamADC);

    IT22 = (b/d)*(1.0-ParamADC) - ((b*e)/(d*d))*log((e+d)/(e+d*ParamADC));
    if (SpecialEllipticScheme)
    {
      //IT22 = (b*ParamFluxADC)/((d*ParamFluxADC)+e);
      IT22 = (b*ParamFluxADCNew)/((d*ParamFluxADC)+e);//@HAF:30.10.2008
    }

    ret = IT21 + IT22;//28.08.07: OK, eller kommentere ut IT21?
  }

  return ret;
  //return 1.5*ret;//Special MGEdwards test :-)
}



//*********************************************************************************************************
//Some private functions for Edwards-type discretization:

double SFMpfaFps::ml_x_xi_Forxi_ANDeta_eq1_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x3.u_getElement_0u(0) - x4.u_getElement_0u(0));
}

double SFMpfaFps::ml_x_xi_Forxi_ANDeta_eq1_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return 0.0;
}

double SFMpfaFps::ml_x_xi_Forxi_eq1ANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x2.u_getElement_0u(0) - x1.u_getElement_0u(0));
}

double SFMpfaFps::ml_x_xi_Forxi_eq1ANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x1.u_getElement_0u(0) - x2.u_getElement_0u(0) + x3.u_getElement_0u(0) - x4.u_getElement_0u(0));
}



double SFMpfaFps::ml_y_xi_Forxi_ANDeta_eq1_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x3.u_getElement_0u(1) - x4.u_getElement_0u(1));
}

double SFMpfaFps::ml_y_xi_Forxi_ANDeta_eq1_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return 0.0;
}

double SFMpfaFps::ml_y_xi_Forxi_eq1ANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x2.u_getElement_0u(1) - x1.u_getElement_0u(1));
}

double SFMpfaFps::ml_y_xi_Forxi_eq1ANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x1.u_getElement_0u(1) - x2.u_getElement_0u(1) + x3.u_getElement_0u(1) - x4.u_getElement_0u(1));
}




double SFMpfaFps::ml_x_eta_Forxi_ANDeta_eq1_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x4.u_getElement_0u(0) - x1.u_getElement_0u(0));
}

double SFMpfaFps::ml_x_eta_Forxi_ANDeta_eq1_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x1.u_getElement_0u(0) - x2.u_getElement_0u(0) + x3.u_getElement_0u(0) - x4.u_getElement_0u(0));
}

double SFMpfaFps::ml_x_eta_Forxi_eq1ANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x3.u_getElement_0u(0) - x2.u_getElement_0u(0));
}

double SFMpfaFps::ml_x_eta_Forxi_eq1ANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return 0.0;
}




double SFMpfaFps::ml_y_eta_Forxi_ANDeta_eq1_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x4.u_getElement_0u(1) - x1.u_getElement_0u(1));
}

double SFMpfaFps::ml_y_eta_Forxi_ANDeta_eq1_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x1.u_getElement_0u(1) - x2.u_getElement_0u(1) + x3.u_getElement_0u(1) - x4.u_getElement_0u(1));
}

double SFMpfaFps::ml_y_eta_Forxi_eq1ANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x3.u_getElement_0u(1) - x2.u_getElement_0u(1));
}

double SFMpfaFps::ml_y_eta_Forxi_eq1ANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return 0.0;
}




double SFMpfaFps::ml_determinant_Forxi_ANDeta_eq1_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  double res = (x3.u_getElement_0u(0) - x4.u_getElement_0u(0))*(x4.u_getElement_0u(1) - x1.u_getElement_0u(1)) -
    (x3.u_getElement_0u(1) - x4.u_getElement_0u(1))*(x4.u_getElement_0u(0) - x1.u_getElement_0u(0));

  return res;
}

double SFMpfaFps::ml_determinant_Forxi_ANDeta_eq1_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  double res = (x3.u_getElement_0u(0) - x4.u_getElement_0u(0))*(x1.u_getElement_0u(1) - x2.u_getElement_0u(1)) -
    (x3.u_getElement_0u(1) - x4.u_getElement_0u(1))*(x1.u_getElement_0u(0) - x2.u_getElement_0u(0));

  return res;
}

double SFMpfaFps::ml_determinant_Forxi_eq1ANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  double res = (x3.u_getElement_0u(1) - x2.u_getElement_0u(1))*(x2.u_getElement_0u(0) - x1.u_getElement_0u(0)) -
    (x3.u_getElement_0u(0) - x2.u_getElement_0u(0))*(x2.u_getElement_0u(1) - x1.u_getElement_0u(1));

  return res;
}

double SFMpfaFps::ml_determinant_Forxi_eq1ANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  double res = (x3.u_getElement_0u(1) - x2.u_getElement_0u(1))*(x1.u_getElement_0u(0) - x4.u_getElement_0u(0)) -
    (x3.u_getElement_0u(0) - x2.u_getElement_0u(0))*(x1.u_getElement_0u(1) - x4.u_getElement_0u(1));

  return res;
}


//NOTE: The following is only needed for bilinear expansions:


double SFMpfaFps::ml_x_xi_Forxi_ANDeta_eq0_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x2.u_getElement_0u(0) - x1.u_getElement_0u(0));
}

double SFMpfaFps::ml_x_xi_Forxi_ANDeta_eq0_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return 0.0;
}

double SFMpfaFps::ml_x_xi_Forxi_eq0ANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x2.u_getElement_0u(0) - x1.u_getElement_0u(0));
}

double SFMpfaFps::ml_x_xi_Forxi_eq0ANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x1.u_getElement_0u(0) - x2.u_getElement_0u(0) + x3.u_getElement_0u(0) - x4.u_getElement_0u(0));
}

double SFMpfaFps::ml_y_xi_Forxi_ANDeta_eq0_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x2.u_getElement_0u(1) - x1.u_getElement_0u(1));
}

double SFMpfaFps::ml_y_xi_Forxi_ANDeta_eq0_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return 0.0;
}

double SFMpfaFps::ml_y_xi_Forxi_eq0ANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x2.u_getElement_0u(1) - x1.u_getElement_0u(1));
}

double SFMpfaFps::ml_y_xi_Forxi_eq0ANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x1.u_getElement_0u(1) - x2.u_getElement_0u(1) + x3.u_getElement_0u(1) - x4.u_getElement_0u(1));
}

double SFMpfaFps::ml_x_eta_Forxi_ANDeta_eq0_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x4.u_getElement_0u(0) - x1.u_getElement_0u(0));
}

double SFMpfaFps::ml_x_eta_Forxi_ANDeta_eq0_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
   return (x1.u_getElement_0u(0) - x2.u_getElement_0u(0) + x3.u_getElement_0u(0) - x4.u_getElement_0u(0));
}

double SFMpfaFps::ml_x_eta_Forxi_eq0ANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x4.u_getElement_0u(0) - x1.u_getElement_0u(0));
}

double SFMpfaFps::ml_x_eta_Forxi_eq0ANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return 0.0;
}

double SFMpfaFps::ml_y_eta_Forxi_ANDeta_eq0_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x4.u_getElement_0u(1) - x1.u_getElement_0u(1));
}

double SFMpfaFps::ml_y_eta_Forxi_ANDeta_eq0_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x1.u_getElement_0u(1) - x2.u_getElement_0u(1) + x3.u_getElement_0u(1) - x4.u_getElement_0u(1));
}

double SFMpfaFps::ml_y_eta_Forxi_eq0ANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x4.u_getElement_0u(1) - x1.u_getElement_0u(1));
}

double SFMpfaFps::ml_y_eta_Forxi_eq0ANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return 0.0;
}

double SFMpfaFps::ml_determinant_Forxi_ANDeta_eq0_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  double res = (x2.u_getElement_0u(0) - x1.u_getElement_0u(0))*(x4.u_getElement_0u(1) - x1.u_getElement_0u(1)) -
    (x2.u_getElement_0u(1) - x1.u_getElement_0u(1))*(x4.u_getElement_0u(0) - x1.u_getElement_0u(0));

  return res;
}

double SFMpfaFps::ml_determinant_Forxi_ANDeta_eq0_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  double res = (x2.u_getElement_0u(0) - x1.u_getElement_0u(0))*(x3.u_getElement_0u(1) - x4.u_getElement_0u(1)) -
    (x2.u_getElement_0u(1) - x1.u_getElement_0u(1))*(x3.u_getElement_0u(0) - x4.u_getElement_0u(0));

  return res;
}

double SFMpfaFps::ml_determinant_Forxi_eq0ANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  double res = (x4.u_getElement_0u(1) - x1.u_getElement_0u(1))*(x2.u_getElement_0u(0) - x1.u_getElement_0u(0)) -
    (x4.u_getElement_0u(0) - x1.u_getElement_0u(0))*(x2.u_getElement_0u(1) - x1.u_getElement_0u(1));

  return res;
}

double SFMpfaFps::ml_determinant_Forxi_eq0ANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  double res = (x4.u_getElement_0u(1) - x1.u_getElement_0u(1))*(x3.u_getElement_0u(0) - x2.u_getElement_0u(0)) -
    (x4.u_getElement_0u(0) - x1.u_getElement_0u(0))*(x3.u_getElement_0u(1) - x2.u_getElement_0u(1));

  return res;
}


//NOTE: The following is only needed for bilinear expansions in the case of 
//auxillary dual-cells

double SFMpfaFps::ml_x_xi_Forxi_ANDeta_eqPARAM_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4, const double& eta) const
{
  return (-(1.0-eta)*x1.u_getElement_0u(0) + (1.0-eta)*x2.u_getElement_0u(0) +
	  eta*x3.u_getElement_0u(0) - eta*x4.u_getElement_0u(0));
}

double SFMpfaFps::ml_x_xi_Forxi_ANDeta_eqPARAM_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return 0.0;
}

double SFMpfaFps::ml_x_xi_Forxi_eqPARAMANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x2.u_getElement_0u(0) - x1.u_getElement_0u(0));
}

double SFMpfaFps::ml_x_xi_Forxi_eqPARAMANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x1.u_getElement_0u(0) - x2.u_getElement_0u(0) + x3.u_getElement_0u(0) - x4.u_getElement_0u(0));
}

double SFMpfaFps::ml_y_xi_Forxi_ANDeta_eqPARAM_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4, const double& eta) const
{
  return (-(1.0-eta)*x1.u_getElement_0u(1) + (1.0-eta)*x2.u_getElement_0u(1) +
	  eta*x3.u_getElement_0u(1) - eta*x4.u_getElement_0u(1));
}

double SFMpfaFps::ml_y_xi_Forxi_ANDeta_eqPARAM_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return 0.0;
}

double SFMpfaFps::ml_y_xi_Forxi_eqPARAMANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x2.u_getElement_0u(1) - x1.u_getElement_0u(1));
}

double SFMpfaFps::ml_y_xi_Forxi_eqPARAMANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x1.u_getElement_0u(1) - x2.u_getElement_0u(1) + x3.u_getElement_0u(1) - x4.u_getElement_0u(1));
}

double SFMpfaFps::ml_x_eta_Forxi_ANDeta_eqPARAM_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x4.u_getElement_0u(0) - x1.u_getElement_0u(0));
}

double SFMpfaFps::ml_x_eta_Forxi_ANDeta_eqPARAM_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x1.u_getElement_0u(0) - x2.u_getElement_0u(0) + x3.u_getElement_0u(0) - x4.u_getElement_0u(0));
}

double SFMpfaFps::ml_x_eta_Forxi_eqPARAMANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4, const double& xi) const
{
  return (-(1.0-xi)*x1.u_getElement_0u(0) - xi*x2.u_getElement_0u(0) +
	  xi*x3.u_getElement_0u(0) + (1.0-xi)*x4.u_getElement_0u(0));
}

double SFMpfaFps::ml_x_eta_Forxi_eqPARAMANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return 0.0;
}

double SFMpfaFps::ml_y_eta_Forxi_ANDeta_eqPARAM_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x4.u_getElement_0u(1) - x1.u_getElement_0u(1));
}

double SFMpfaFps::ml_y_eta_Forxi_ANDeta_eqPARAM_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return (x1.u_getElement_0u(1) - x2.u_getElement_0u(1) + x3.u_getElement_0u(1) - x4.u_getElement_0u(1));
}

double SFMpfaFps::ml_y_eta_Forxi_eqPARAMANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4, const double& xi) const
{
  return (-(1.0-xi)*x1.u_getElement_0u(1) - xi*x2.u_getElement_0u(1) +
	  xi*x3.u_getElement_0u(1) + (1.0-xi)*x4.u_getElement_0u(1));
}

double SFMpfaFps::ml_y_eta_Forxi_eqPARAMANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4) const
{
  return 0.0;
}


double SFMpfaFps::ml_determinant_Forxi_ANDeta_eqPARAM_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4, const double& eta) const
{
  double a = ml_x_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x1, x2, x3, x4, eta);
  double b = ml_x_eta_Forxi_ANDeta_eqPARAM_ConstantPart(x1, x2, x3, x4);
  double d = ml_y_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x1, x2, x3, x4, eta);
  double e = ml_y_eta_Forxi_ANDeta_eqPARAM_ConstantPart(x1, x2, x3, x4);

  double res = a*e - d*b;

  return res;
}

double SFMpfaFps::ml_determinant_Forxi_ANDeta_eqPARAM_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4, const double& eta) const
{
  double a = ml_x_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x1, x2, x3, x4, eta);
  double c = ml_x_eta_Forxi_ANDeta_eqPARAM_LinearPart(x1, x2, x3, x4);
  double d = ml_y_xi_Forxi_ANDeta_eqPARAM_ConstantPart(x1, x2, x3, x4, eta);
  double f = ml_y_eta_Forxi_ANDeta_eqPARAM_LinearPart(x1, x2, x3, x4);

  double res = a*f - d*c;

  return res;
}

double SFMpfaFps::ml_determinant_Forxi_eqPARAMANDeta_ConstantPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4, const double& xi) const
{
  double a = ml_x_xi_Forxi_eqPARAMANDeta_ConstantPart(x1, x2, x3, x4);
  double c = ml_x_eta_Forxi_eqPARAMANDeta_ConstantPart(x1, x2, x3, x4, xi);
  double d = ml_y_xi_Forxi_eqPARAMANDeta_ConstantPart(x1, x2, x3, x4);
  double f = ml_y_eta_Forxi_eqPARAMANDeta_ConstantPart(x1, x2, x3, x4, xi);

  double res = a*f - d*c;

  return res;
}

double SFMpfaFps::ml_determinant_Forxi_eqPARAMANDeta_LinearPart(const RnPoint& x1, const RnPoint& x2, const RnPoint& x3, const RnPoint& x4, const double& xi) const
{
  double b = ml_x_xi_Forxi_eqPARAMANDeta_LinearPart(x1, x2, x3, x4);
  double c = ml_x_eta_Forxi_eqPARAMANDeta_ConstantPart(x1, x2, x3, x4, xi);
  double e = ml_y_xi_Forxi_eqPARAMANDeta_LinearPart(x1, x2, x3, x4);
  double f = ml_y_eta_Forxi_eqPARAMANDeta_ConstantPart(x1, x2, x3, x4, xi);

  double res = b*f - e*c;

  return res;
}

//*********************************************************************************************************
