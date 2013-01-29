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
#include "SFCentralUpwind.h"

#include <cmath>


//---------------------  Constants --------------------------------------------
#define TOLERANCE 10e-9


//---------------------- Types ------------------------------------------------








//---------------------- Public Functions -------------------------------------





// Constructors

//-----------------------------------------------------------------------------
SFCentralUpwind::SFCentralUpwind()
  : IGIPtr_(), saturationGrad_(), saturationGradNew_(), fluxIsComputed_()
//-----------------------------------------------------------------------------
{

}

//-----------------------------------------------------------------------------
SFCentralUpwind::SFCentralUpwind(const SFCentralUpwind& SFC)
  : saturationGrad_(SFC.saturationGrad_), saturationGradNew_(SFC.saturationGradNew_), fluxIsComputed_(SFC.fluxIsComputed_)
//-----------------------------------------------------------------------------
{
  IGIPtr_ = SFC.IGIPtr_;
}


// Destructors

//-----------------------------------------------------------------------------
SFCentralUpwind::~SFCentralUpwind()
//-----------------------------------------------------------------------------
{

}



// Generators

//-----------------------------------------------------------------------------
void SFCentralUpwind::setup(const IRISDuneGridInterface<DuneGridType>* IgridPtr)
//-----------------------------------------------------------------------------
{
  IGIPtr_ = IgridPtr;
  int numbOfElements = IGIPtr_->cellCount(0,0);
  RnShape Rsh(IGIPtr_->domainDimension());
  RnPoint Rp;
  Rp.u_setPoint_1u(Rsh,0.0);
  saturationGrad_.u_copy_1u(numbOfElements,Rp); 
  saturationGradNew_.u_copy_1u(numbOfElements,Rp);

  fluxIsComputed_.u_copy_1u(numbOfElements, false);
}
  

// Observers 

//----------------------------------------------------------------------------
double SFCentralUpwind::getMinGridSize() const
//----------------------------------------------------------------------------
{
  return IGIPtr_->getMinGridSize();//@HAF: This method MUST be implemented in IRISDuneGridInterface
}



//----------------------------------------------------------------------------
double SFCentralUpwind::computeTimeStep(const double& minGridSize, const double& maxVelocity) const
//----------------------------------------------------------------------------
{
  //@HAF:It is for the time being assumed that the max. velocity is 1. Eller ikke???
  //@HAF:This must later be done more general !!!
  //**************************************************************************** 
  //This routine computes a time step DT, which fulfills the CFL requirement
  //of the central-upwind scheme on triangular grids (Kurganov & Petrova (2005)).
  //It is nonadaptive and is thus computed only once, and is moreover based
  //on a maximum estimate of:
  //          totalvelcity*derivative of (water) fractional flow in the (two-phase) saturation equation.
  //i.e. the maximum of both factors.
  //Finally, it also needs the minimum altitude of the triangles in the grid.
  //****************************************************************************
  double initTStep = -1.0;

  //@HAF: Maa sjekk gridtype (e.g. structured 2D or triangular)
  IRISDuneGridInterface<DuneGridType>::GridType gridType = IGIPtr_->getGridType();//@HAF: gridType MUST of course be an enum (of type GridType see IRISDuneGridInterface)!!!

  //NOTE: Code only for 2D yet
  if (IGIPtr_->cellDimension() == 2)
  {
    if (gridType == IRISDuneGridInterface<DuneGridType>::_Triangular)
    {
      initTStep = 0.25*(minGridSize/(3.0*maxVelocity));//@HAF: MUST be checked against OLD code i.a. wrt. maxVelocity (instead of cMax...)
    }
    else if (gridType == IRISDuneGridInterface<DuneGridType>::_Structured2D)
    {
      double CFLConst = 0.75;
      //initTStep = (1.0/3.0)*(CFLConst/cMax)*deltax;
      initTStep = (1.0/3.0)*(CFLConst/maxVelocity)*minGridSize; //OK i 2D, med DY=DX. @HAF: MUST be checked against OLD code i.a. wrt. maxVelocity (instead of cMax...)
    }
  }
  else if (IGIPtr_->cellDimension() == 3)
  {
    //3D code yet to come...
  }

  return initTStep;
}


//----------------------------------------------------------------------------
void SFCentralUpwind::computeFluxOutOfAllCells(const SFRelPerm& relPerm, const SFPMPhysicalFieldProperties& PMPhysProp, SFSaturation& saturation, const SFVelocity& velocity_a, Array<double>& rhs)
//----------------------------------------------------------------------------
{
  //NOTE: Code only for 2D yet
  if (IGIPtr_->cellDimension() == 2)
  {

    //@HAF: Maa sjekk gridtype (i.e. quadrilateral or triangular)
    IRISDuneGridInterface<DuneGridType>::GridType gridType = IGIPtr_->getGridType();//@HAF: gridType MUST of course be an enum (of type GridType see IRISDuneGridInterface)!!!
    if (gridType == IRISDuneGridInterface<DuneGridType>::_Triangular)
    {
      Array<int> rockType = PMPhysProp.getRockType();//@HAF: OK???

//       double TOL = 0.000001;
//       int cellDim = IGIPtr_->cellDimension();

//       int numbOfEdges = IGIPtr_->cellCount(1,0); 
      int numbOfElements = IGIPtr_->cellCount(0,0);

      //Initialization:
      for (int i=0; i < numbOfElements; i++)
      {
	rhs.u_changeElement_1u(i, 0.0);
      }

      //Compute the saturation gradient and the limiter coefficients:

      ml_computeSaturationGradientArminjon(saturation);
      ml_MINMODLimitingOfSaturationGradients();
      
      RnPoint Mk, NVec, contrib;
      double dx_hk, uj1, uj2, veloKurg, totVelN, help, flux, flux_uj1, flux_uj2;
      int j1, j2; //, dirInfo_j1, dirInfo_j2;

      for (IRISDuneGridInterface<DuneGridType>::ElementIterator ihatIt = IGIPtr_->setEntityPointerToFirst<IRISDuneGridInterface<DuneGridType>::ElementIterator>(0); ihatIt!=IGIPtr_->entityPointerIsAtEnd<IRISDuneGridInterface<DuneGridType>::ElementIterator>(0); ++ihatIt)
      {
	IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt = ihatIt;
	j1 = IGIPtr_->getGlobalIndexFromEntityPointer(ihatPt);//@HAF: 16.11.2009, OK???
	//Fetch the element neighbours of the element (as element pointers):
	std::vector<IRISDuneGridInterface<DuneGridType>::ElementPointer> ihatPtElemNeigh;
	//ihatPtElemNeigh.reserve(3);//@HAF: 15.01.2010, OK???
	IGIPtr_->getNeighbours(ihatPt, ihatPtElemNeigh);
	int idxNeigh = 0;
	for (int k=0; k < 3; k++)
	{
	  Mk = IGIPtr_->getCentroid<1>(ihatPt, k);//@HAF: NOTE: Her er det clash med ekte RnPoint og "jukse" RnPoint som getCentroid lever...
	  //dx_hk = edgeLength[k];//@HAF: Not needed now I think...
	  //*************************************************************
	  dx_hk = 1.0;//@HAF change: Test 08.07.2008
	  //This change is very important since the velocity now is
	  //computed as a flux in contrast to the previous method.
	  //*************************************************************
	  
	  //NB: Must check for the boundary:
	  if (!IGIPtr_->boundaryIndex<1>(ihatPt, k))//@HAF: OK???
	  {
	    //The edge we are treating is NOT on the boundary !!!	    
	    j2 = IGIPtr_->getGlobalIndexFromEntityPointer(ihatPtElemNeigh[idxNeigh]);//@HAF: 16.11.2009, OK???
	    if (!fluxIsComputed_[j2])
	    {
	      //Computes uj1 and uj2:
	      contrib = Mk;
	      //contrib -= centroidOfTriangle[j1];
	      contrib -= IGIPtr_->getCentroid(ihatPt);
	      uj1 = saturation[j1] + ml_MINMOD_2ARG(contrib.u_dot_0u(saturationGradNew_[j1]), contrib.u_dot_0u(saturationGradNew_[j2]));
	    
	      contrib = Mk;
	      //contrib -= centroidOfTriangle[j2];
	      contrib -= IGIPtr_->getCentroid(ihatPtElemNeigh[idxNeigh]);
	      uj2 = saturation[j2] + ml_MINMOD_2ARG(contrib.u_dot_0u(saturationGradNew_[j1]), contrib.u_dot_0u(saturationGradNew_[j2]));
	      
	      //We compute the convective fluxes connected to uj1 and uj2:
	      //totVelN = dirInfo_j1*velocity_a[k];//Note: this is the velocity_a wrt. j1 (i.e. measured in the dir. out of j1)!
	      totVelN = velocity_a[(3*j1)+k];//Note:@HAF: OK??? this is the velocity_a wrt. j1 (i.e. measured in the dir. out of j1)!
	      
	      help = relPerm.fractionalFlowWater(uj1,rockType[j1],PMPhysProp);
	      flux_uj1 = help*totVelN;
	      
	      help = relPerm.fractionalFlowWater(uj2,rockType[j2],PMPhysProp);
	      flux_uj2 = help*totVelN;
	      
	      //Must compute the "velocity" in the central-upwind scheme:
	      //NB: MUST be changed/expanded when including types !!!
	      veloKurg = ml_computeVelocity_Kurganov(uj1,uj2,rockType[j1],rockType[j2],totVelN,relPerm,PMPhysProp);
	      
	      //Finally, we implement the total-central upwind flux for this edge from the "point og view" of j1:
	      flux = dx_hk*(-0.5*(flux_uj1 + flux_uj2) + veloKurg*(uj2 - uj1));

	      rhs.u_changeElement_1u(j1, rhs[j1] + (1.0/IGIPtr_->getVolume(ihatPt))*flux);
	      rhs.u_changeElement_1u(j2, rhs[j2] + (1.0/IGIPtr_->getVolume(ihatPtElemNeigh[idxNeigh]))*(-1.0)*flux);
	    }
	    idxNeigh++;
	  }
	  else
	  {
	    //**********************************************************************
	    //For the time being we only treat two BC's
	    //a) Closed boundaries. In which there is nothing to do here!!!
	    //b) In addition to a), we also have an imposed BC on the inlet.
	    //For the time being b) is "VERY HANDCODED". See below for _Structured2D!!!
	    //**********************************************************************
	    //The edge we are treating is ON the "real" outer boundary !!!

	  }
	}
	fluxIsComputed_.u_changeElement_1u(j1, true);
      }
      

      
      //*******************************************************************************
      //In the case of a Dirichlet BC, we must loop through the edges once more. This is due to the
      //fact that we do not want to alter the solution where we have a Dirichlet BC, 
      //and thus we must put the "total" flux to zero for the respective boundary adjacent elements:
      //*******************************************************************************
      
      //@HAF: To be continued 17.11.2009
      //Maa loope gjennom Ã¥ sjekke pÃ¥ korridnater om vi er naer nok innlet. break saa
      //ut av innerste ("k-loop") og sÃ¸org for at rhs bidraget blir null etc....!!!!!!!!!!!!!
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //Styr videre alt med en ifdef for enkelhets skyld!!!


    }
    else if (gridType == IRISDuneGridInterface<DuneGridType>::_Structured2D)
    {
      Array<int> rockType = PMPhysProp.getRockType();//@HAF: OK???

      int NX  = PMPhysProp.getSimulationDataPtr()->td_.NoGridCells_X;//@HAF: OK???
      int NY  = PMPhysProp.getSimulationDataPtr()->td_.NoGridCells_Y;//@HAF: OK???
      double DX  = PMPhysProp.getSimulationDataPtr()->td_.X_Extent/NX;//@HAF: OK???
      double DY  = PMPhysProp.getSimulationDataPtr()->td_.Y_Extent/NY;//@HAF: OK???

      //*********************************************************
      //*********************************************************

//       double TOL = 0.000001;
//       int cellDim = IGIPtr_->cellDimension();

//       int numbOfEdges = IGIPtr_->cellCount(1,0); 
      int numbOfElements = IGIPtr_->cellCount(0,0);

      //Initialization:
      for (int i=0; i < numbOfElements; i++)
      {
	rhs.u_changeElement_1u(i, 0.0);
      }

      //Compute the limiter coefficients:
      ml_generateLimitersForStructuredGrid(saturation,NX,NY,DX,DY);//@HAF: Must be more thoroughly checked???It might be something is WRONG with directions etc. pr. 02.12.2009.
      
      RnPoint Mk, NVec, contrib;
      double dx_hk, uj1, uj2, veloKurg, totVelN, help, flux, flux_uj1, flux_uj2;
      int j1, j2;//, dirInfo_j1, dirInfo_j2;

      for (IRISDuneGridInterface<DuneGridType>::ElementIterator ihatIt = IGIPtr_->setEntityPointerToFirst<IRISDuneGridInterface<DuneGridType>::ElementIterator>(0); ihatIt!=IGIPtr_->entityPointerIsAtEnd<IRISDuneGridInterface<DuneGridType>::ElementIterator>(0); ++ihatIt)
      {
	IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt = ihatIt;
	j1 = IGIPtr_->getGlobalIndexFromEntityPointer(ihatPt);//@HAF: 16.11.2009, OK???
	//Fetch the element neighbours of the element (as element pointers):
	std::vector<IRISDuneGridInterface<DuneGridType>::ElementPointer> ihatPtElemNeigh;
	//ihatPtElemNeigh.reserve(3);//@HAF: 15.01.2010, OK???
	IGIPtr_->getNeighbours(ihatPt, ihatPtElemNeigh);
	//cout << " j1= " << j1 << " ihatPtElemNeigh= " << " size= " << ihatPtElemNeigh.size() << " Elements= " << IGIPtr_->getGlobalIndexFromEntityPointer(ihatPtElemNeigh[0]) << " , " << IGIPtr_->getGlobalIndexFromEntityPointer(ihatPtElemNeigh[1]) << endl;
	//cout << "-----------------------------------------------------------" << endl;
	//exit(0);
	int idxNeigh = 0;
	for (int k=0; k < 4; k++)//@HAF: Since we are on a structured quadrilateral grid
	{
	  Mk = IGIPtr_->getCentroid<1>(ihatPt, k);
	  //dx_hk = edgeLength[k];//@HAF: Not needed now I think...
	  //*************************************************************
	  dx_hk = 1.0;//@HAF change: Test 08.07.2008
	  //This change is very important since the velocity now is
	  //computed as a flux in contrast to the previous method.
	  //************************************************************* 
	  
	  //NB: Must check for the boundary:
	  //if (j2 >= 0)
	  if (!IGIPtr_->boundaryIndex<1>(ihatPt, k))//@HAF: OK???
	  {
	    //The edge we are treating is NOT on the boundary !!!	    
	    j2 = IGIPtr_->getGlobalIndexFromEntityPointer(ihatPtElemNeigh[idxNeigh]);//@HAF: 18.02.2010, OK???
	    //*************DEBUGGING***********START****************************
	    //cout << "------------------------------------------------------" << endl;
	    //cout << "j1= " << j1 << " fluxIsComputed_[j1]= " << fluxIsComputed_[j1] << " j2= " << j2 << " fluxIsComputed_[j2]= "  << fluxIsComputed_[j2] << endl;
	    //cout << "------------------------------------------------------" << endl;
	    //*************DEBUGGING***********END******************************
	    if (!fluxIsComputed_[j2])
	    {
	      //Computes uj1 and uj2:
	      contrib = Mk;
	      //contrib -= centroidOfTriangle[j1];
	      contrib -= IGIPtr_->getCentroid(ihatPt);
	      uj1 = saturation[j1] + contrib.u_dot_0u(saturationGradNew_[j1]);
	    
	      contrib = Mk;
	      //contrib -= centroidOfTriangle[j2];
	      contrib -= IGIPtr_->getCentroid(ihatPtElemNeigh[idxNeigh]);
	      uj2 = saturation[j2] + contrib.u_dot_0u(saturationGradNew_[j2]);
	      
	      //We compute the convective fluxes connected to uj1 and uj2:
	      //totVelN = dirInfo_j1*velocity_a[k];//Note: this is the velocity_a wrt. j1 (i.e. measured in the dir. out of j1)!
	      totVelN = velocity_a[(4*j1)+k];//Note:@HAF: OK??? this is the velocity_a wrt. j1 (i.e. measured in the dir. out of j1)!//NOTE: pass paa faktoren  4 i velocity_a!!!
	      
	      help = relPerm.fractionalFlowWater(uj1,rockType[j1],PMPhysProp);
	      flux_uj1 = help*totVelN;
	      
	      help = relPerm.fractionalFlowWater(uj2,rockType[j2],PMPhysProp);
	      flux_uj2 = help*totVelN;
	      
	      //Must compute the "velocity" in the central-upwind scheme:
	      //NB: MUST be changed/expanded when including types !!!
	      veloKurg = ml_computeVelocity_Kurganov(uj1,uj2,rockType[j1],rockType[j2],totVelN,relPerm,PMPhysProp);
	      
	      //Finally, we implement the total-central upwind flux for this edge from the "point og view" of j1:
	      flux = dx_hk*(-0.5*(flux_uj1 + flux_uj2) + veloKurg*(uj2 - uj1));

	      rhs.u_changeElement_1u(j1, rhs[j1] + (1.0/IGIPtr_->getVolume(ihatPt))*flux);
	      rhs.u_changeElement_1u(j2, rhs[j2] + (1.0/IGIPtr_->getVolume(ihatPtElemNeigh[idxNeigh]))*(-1.0)*flux);
	    }
	    idxNeigh++;
	  }
	  else
	  {
	    //**********************************************************************
	    //For the time being we only treat two BC's
	    //a) Closed boundaries. In which there is nothing to do here!!!
	    //b) In addition to a), we also have an imposed BC on the inlet.
	    //For the time being b) is "VERY HANDCODED" (see below).
	    //c)Finally the fluxes must be computed at the outlet.
	    //We let the code for a) and c) (and also b) before the change below)
	    //be similar.
	    //**********************************************************************
	    //The edge we are treating is ON the "real" outer boundary !!!
	    double actualBCSaturation;
	    int BCTypeSaturation = ml_getNonPeriodicBCType(PMPhysProp,ihatPt,k);
	    if (BCTypeSaturation == 10)
	    {
	      //We have a Neumann BC for the saturation:
	      actualBCSaturation = saturation[j1];
	    }
	    else if (BCTypeSaturation == 20)
	    {
	      //We have a Dirichlet BC for the saturation:
	      actualBCSaturation = ml_getBCDirichletCond(PMPhysProp,ihatPt,k);
	    }


	    //*******************************************************************************
	    //Remark:
	    //Here the boundary is treated from the following principle:
	    //IF we have a Dirichlet type of BC, then we should NOT update the
	    //solution at the boundary. This can be performed simply by setting
	    // f[j1] = 0.0             OK ??????????????????????????????????? 
	    //*******************************************************************************

	    uj1 = actualBCSaturation;
	    //totVelN = dirInfo_j1*velocity_a[k];
	    totVelN = velocity_a[(4*j1)+k];
      
	    help = relPerm.fractionalFlowWater(uj1,rockType[j1],PMPhysProp);
	    flux_uj1 = help*totVelN;

	    //***************JUST**TO**TEST****START************************************** 
	    //@HAF: Hardcoding of a flux-BC at the inlet???
	    //if (CSaturation.u_getBoundaryNumber_0u(GInm1) == 0) flux_uj1 = velocity_a[k];//NB: We are considering the inlet.
	    //***************JUST**TO**TEST****END**************************************** 	

	    //f[j1] += dirInfo_j1*(1.0/volumes[j1])*dx_hk*(-1.0*flux_uj1);
	    rhs.u_changeElement_1u(j1, rhs[j1] + (1.0/IGIPtr_->getVolume(ihatPt))*dx_hk*(-1.0*flux_uj1));
	  }
	}
	fluxIsComputed_.u_changeElement_1u(j1, true);
      }
      

      
      //*******************************************************************************
      //In the case of a Dirichlet BC, we must loop through the edges once more. This is due to the
      //fact that we do not want to alter the solution where we have a Dirichlet BC, 
      //and thus we must put the "total" flux to zero for the respective boundary adjacent elements:
      //*******************************************************************************
      
      //#ifdef PAPER_EX2 (Dirichlet BC's at left (inlet) boundary).
      //@HAF:NOTE: This is a VERY TEMPORARY way to do it !!! *********************************
      //*********************************************************************************
      int DNUMB1 = NX;
      int DNUMB2 = NY;

      //BC for main flow in x-direction
      for (int j=0; j < DNUMB2; j++)
      {
	saturation.u_changeElement_1u(j*DNUMB1, 1.0);
	rhs.u_changeElement_1u(j*DNUMB1, 0.0);
      }
      //#endif //PAPER_EX2  
    }
  }
  else if (IGIPtr_->cellDimension() == 3)
  {
    //3D code yet to come...
  }


  //**************************************************************************
  //@HAF: Must finally reset the help-variable fluxIsComputed_:
  for (int i=0; i < IGIPtr_->cellCount(0, 0); i++)
  {
    fluxIsComputed_.u_changeElement_1u(i, false);
  }
}

 

/* ---------------------------------------------------------------------------*/
//---------------------- Private Functions -------------------------------------
/* ---------------------------------------------------------------------------*/

double SFCentralUpwind::ml_computeVelocity_Kurganov(const double& sPluss, const double& sMinus, const int& RTPluss, const int& RTMinus, const double& totalVelocityN, const SFRelPerm& relPerm, const SFPMPhysicalFieldProperties& PMPhysProp) 
{
  //This function computes the velocitiy required by the Kurganov & Tadmor algorithm
  // for the convective flux in the saturation equation.

  double maxVel;

  double fracFlowDeriv_p = relPerm.fractionalFlowWaterDerivative(sPluss,RTPluss,PMPhysProp);
  double fracFlowDeriv_m = relPerm.fractionalFlowWaterDerivative(sMinus,RTMinus,PMPhysProp);

  double fVal_p = fracFlowDeriv_p*totalVelocityN;

  double fVal_m = fracFlowDeriv_m*totalVelocityN;

  maxVel = (fabs(fVal_p) >= fabs(fVal_m)) ? fVal_p : fVal_m;

  return fabs(maxVel);
  //return maxVel; //NB: absolute value func. has been removed. MUST be thoroughly checked !!!
}



void SFCentralUpwind::ml_computeSaturationGradientArminjon(const SFSaturation& saturation)
{
  //*****************************************************************************************************
  //This function computes the saturation gradient for each triangle in an unstructured triangular grid
  //based on a least-square approach (Arminjon et al. (1997)).
  //*****************************************************************************************************

  //Array<int> elementNeigh = GridInfo.u_getDirectElementInfoFromElement_0u();
  std::vector<IRISDuneGridInterface<DuneGridType>::ElementPointer> elementNeigh;
  //elementNeigh.reserve(3);//@HAF: 15.01.2010, OK???
  //Array<RnPoint> centroidPoints = GridInfo.u_getTriangleCentroid_0u();

  //We must loop through each element. First we must find the actual number of element neighbours
  //surrounding the given element, and then implement the least square solution (Arminjon et al. (1997)).

  double D, L_YG2, L_XG2, L_XGYG, L_UTXG, L_UTYG;
  double XC, YC, XC2, YC2;
  RnPoint centroid, centroidN;
  RnPoint gradSat(saturationGrad_[0]);//@HAF: Just to initialize
  int idx, number;

  //for (int i=0; i < numbOfElements; i++)
  for (IRISDuneGridInterface<DuneGridType>::ElementIterator ihatIt = IGIPtr_->setEntityPointerToFirst<IRISDuneGridInterface<DuneGridType>::ElementIterator>(0); ihatIt!=IGIPtr_->entityPointerIsAtEnd<IRISDuneGridInterface<DuneGridType>::ElementIterator>(0); ++ihatIt)
  {
    IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt = ihatIt;
    int i = IGIPtr_->getGlobalIndexFromEntityPointer(ihatPt);
    centroid = IGIPtr_->getCentroid(ihatPt);

    //The various expressions used in the comput. must be set to zero for each new element:
    L_XG2 = 0.0; 
    L_YG2 = 0.0; 
    L_XGYG = 0.0; 
    L_UTXG = 0.0;  
    L_UTYG = 0.0; 

    IGIPtr_->getNeighbours(ihatPt, elementNeigh);
    number = 0;
    for (int j=0; j < 3; j++)
    {
      //idx = elementNeigh[(3*i)+j];
      idx = IGIPtr_->getGlobalIndexFromEntityPointer(elementNeigh[j]);
      //Must check for the boundary and skip it in the computation:
      //if (idx != -1)
      if (!IGIPtr_->boundaryIndex(elementNeigh[j]))//@HAF: FEIL???
      {
	number++;

	centroidN = IGIPtr_->getCentroid(elementNeigh[j]);
	centroidN -= centroid;
	XC = centroidN.u_getElement_0u(0);
	YC = centroidN.u_getElement_0u(1);

	L_XG2 += XC*XC;
	L_YG2 += YC*YC;
	L_XGYG += XC*YC;
	L_UTXG += (saturation[idx]  - saturation[i])*XC; 
	L_UTYG += (saturation[idx]  - saturation[i])*YC; 
      }
    }

    if (number == 3)
    {
      //We can use the least square apprach:
      //      D = D1 - D2;
      D = L_XG2*L_YG2 - L_XGYG*L_XGYG;
      //      cout << endl;
      //cout << " D= " << D;
      //cout << endl;
      gradSat.u_changeElement_1u(0,((L_YG2*L_UTXG)/D) - ((L_XGYG*L_UTYG)/D));
      gradSat.u_changeElement_1u(1,((L_XG2*L_UTYG)/D) - ((L_XGYG*L_UTXG)/D));
    }
    else if (number == 2)
    {
      //We can "solve exactly" without the least square apprach:
      int idx1, idx2;
      for (int j=0; j < 3; j++)
      {
	if (!IGIPtr_->boundaryIndex(elementNeigh[j]))
	{
	  //idx1 = elementNeigh[(3*i)+j];
	  idx1 = IGIPtr_->getGlobalIndexFromEntityPointer(elementNeigh[j]);
	  //GIN.u_updateSecond_1u(idx1);
	  centroidN = IGIPtr_->getCentroid(elementNeigh[j]);
	  centroidN -= centroid;
	  XC = centroidN.u_getElement_0u(0);
	  YC = centroidN.u_getElement_0u(1);
	  break;
	}
      }
      for (int j=0; j < 3; j++)
      {
	//if ((elementNeigh[(3*i)+j] != -1) && (elementNeigh[(3*i)+j] != idx1))
	if ((!IGIPtr_->boundaryIndex(elementNeigh[j])) && (IGIPtr_->getGlobalIndexFromEntityPointer(elementNeigh[j]) != idx1))
	{
	  //idx2 = elementNeigh[(3*i)+j];
	  idx2 = IGIPtr_->getGlobalIndexFromEntityPointer(elementNeigh[j]);
	  //GIN.u_updateSecond_1u(idx2);
	  centroidN = IGIPtr_->getCentroid(elementNeigh[j]);
	  centroidN -= centroid;
	  XC2 = centroidN.u_getElement_0u(0);
	  YC2 = centroidN.u_getElement_0u(1);
	  break;
	}
      }

      double UT = saturation[idx1] - saturation[i];
      double UT2 = saturation[idx2] - saturation[i];

      //Must now find a find a "safe" solution, i.e. avoid division by zero etc.:
      double a,b;
      double TOL = 0.000001; //Should it be less, or bigger ???
      
      if ((fabs(XC) >= TOL) && (fabs(YC) >= TOL) && (fabs(XC2) >= TOL) && (fabs(YC2) >= TOL))
      {
	//Solution will be computed in an ordinary way:
	b = (XC*UT2 - XC2*UT)/(YC2*XC - XC2*YC);
	a = (UT/XC) - (YC/XC)*b;
      }
      else if ((fabs(XC) < TOL) && (fabs(YC2) < TOL))
      {
	//Special case number 1:
	b = UT/YC;
	a = UT2/XC2;
	if ((fabs(YC) < TOL) || (fabs(XC2) < TOL))
	{
	  cout << "1: Trouble computing the (boundary) least square solution. Strange data??? Program is stopped." << endl;
	  exit(0);
	}
      }
      else if ((fabs(YC) < TOL) && (fabs(XC2) < TOL))
      {
	//Special case number 2:
	a = UT/XC;
	b = UT2/YC2;
	if ((fabs(XC) < TOL) || (fabs(YC2) < TOL))
	{
	  cout << "2: Trouble computing the (boundary) least square solution. Strange data??? Program is stopped." << endl;
	  exit(0);
	}
      }
      else if (fabs(XC) < TOL)
      {
	//Special case number 3:
	b = UT/YC;
	if (fabs(XC2) >= TOL)
	{
	  a = (1.0/XC2)*(UT2  - YC2*b);
	}
	else
	{
	  cout << "3: Trouble computing the (boundary) least square solution. Strange data??? Program is stopped." << endl;
	  exit(0);
	}
      }
      else if (fabs(YC) < TOL)
      {
	//Special case number 4:
	a = UT/XC;
	if (fabs(YC2) >= TOL)
	{
	  b = (1.0/YC2)*(UT2  - XC2*a);
	}
	else
	{
	  cout << "4: Trouble computing the (boundary) least square solution. Strange data??? Program is stopped." << endl;
	  exit(0);
	}
      }
      else if (fabs(XC2) < TOL)
      {
	//Special case number 5:
	b = UT2/YC2;
	if (fabs(XC) >= TOL)
	{
	  a = (1.0/XC)*(UT  - YC*b);
	}
	else
	{
	  cout << "5: Trouble computing the (boundary) least square solution. Strange data??? Program is stopped." << endl;
	  exit(0);
	}
      }
      else if (fabs(YC2) < TOL)
      {
	//Special case number 6:
	a = UT2/XC2;
	if (fabs(YC) >= TOL)
	{
	  b = (1.0/YC)*(UT  - XC*a);
	}
	else
	{
	  cout << "6: Trouble computing the (boundary) least square solution. Strange data??? Program is stopped." << endl;
	  exit(0);
	}
      }

      gradSat.u_changeElement_1u(0,a);
      gradSat.u_changeElement_1u(1,b);
    }
    else if (number == 1)
    {
      //We must "solve exactly" in a special way in order to obtain a solution:
      //NOTE In this case we must unfortunately have that the components of the grad
      //are equal in both the x- and y-direction.
      int idx1;
      for (int j=0; j < 3; j++)
      {
	//if (elementNeigh[(3*i)+j] != -1)
	if (!IGIPtr_->boundaryIndex(elementNeigh[j]))
	{
	  //idx1 = elementNeigh[(3*i)+j];
	  idx1 = IGIPtr_->getGlobalIndexFromEntityPointer(elementNeigh[j]);
	  //GIN.u_updateSecond_1u(idx1);
	  centroidN = IGIPtr_->getCentroid(elementNeigh[j]);
	  centroidN -= centroid;
	  XC = centroidN.u_getElement_0u(0);
	  YC = centroidN.u_getElement_0u(1);
	  break;
	}
      }

      double UT = saturation[idx1] - saturation[i];

      double a = (UT)/(XC + YC);

      gradSat.u_changeElement_1u(0,a);
      gradSat.u_changeElement_1u(1,a);
    }

    saturationGrad_.u_changeElement_1u(i,gradSat);				 
  }
}



void SFCentralUpwind::ml_MINMODLimitingOfSaturationGradients()
{
  //*****************************************************************************************************
  //This function applies the MINMOD limiter in a certain way to the components of the gradients already 
  //computed by the routine ml_computeSaturationGradientArminjon.
  //*****************************************************************************************************

  //*****************************************************************************************************
  //NOTE: The routine ml_computeSaturationGradientArminjon MUST have been called before using this routine !!!
  //*****************************************************************************************************

  //Array<int> elementNeigh = GridInfo.u_getDirectElementInfoFromElement_0u();
  std::vector<IRISDuneGridInterface<DuneGridType>::ElementPointer> elementNeigh;
  //elementNeigh.reserve(3);//@HAF: 15.01.2010, OK???

  // int numbOfElements = IGIPtr_->cellCount(0,0);
  int number;
  double arg1, arg2, arg3, arg4;

  RnShape Rsh(IGIPtr_->domainDimension());
  RnPoint gradSat;
  gradSat.u_setPoint_1u(Rsh,0.0);

  //for (int i=0; i < numbOfElements; i++)
  for (IRISDuneGridInterface<DuneGridType>::ElementIterator ihatIt = IGIPtr_->setEntityPointerToFirst<IRISDuneGridInterface<DuneGridType>::ElementIterator>(0); ihatIt!=IGIPtr_->entityPointerIsAtEnd<IRISDuneGridInterface<DuneGridType>::ElementIterator>(0); ++ihatIt)
  {
    IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt = ihatIt;
    int i = IGIPtr_->getGlobalIndexFromEntityPointer(ihatPt);
    IGIPtr_->getNeighbours(ihatPt, elementNeigh);
    number = 0;
    for (int j=0; j < 3; j++)
    {
      //idx = elementNeigh[(3*i)+j];
      //Must check for the boundary and skip it in the computation:
      //if (idx != -1) number++; 
      if (!IGIPtr_->boundaryIndex(elementNeigh[j])) number++; 
    }

    if (number == 3)
    {
      //Must apply the MINMOD limiter on each component:
      arg1 = saturationGrad_[i].u_getElement_0u(0);
      arg2 = saturationGrad_[IGIPtr_->getGlobalIndexFromEntityPointer(elementNeigh[0])].u_getElement_0u(0);
      arg3 = saturationGrad_[IGIPtr_->getGlobalIndexFromEntityPointer(elementNeigh[1])].u_getElement_0u(0);
      arg4 = saturationGrad_[IGIPtr_->getGlobalIndexFromEntityPointer(elementNeigh[2])].u_getElement_0u(0);
    
      gradSat.u_changeElement_1u(0,ml_MINMOD_4ARG(arg1,arg2,arg3,arg4));

      arg1 = saturationGrad_[i].u_getElement_0u(1);
      arg2 = saturationGrad_[IGIPtr_->getGlobalIndexFromEntityPointer(elementNeigh[0])].u_getElement_0u(1);
      arg3 = saturationGrad_[IGIPtr_->getGlobalIndexFromEntityPointer(elementNeigh[1])].u_getElement_0u(1);
      arg4 = saturationGrad_[IGIPtr_->getGlobalIndexFromEntityPointer(elementNeigh[2])].u_getElement_0u(1);

      gradSat.u_changeElement_1u(1,ml_MINMOD_4ARG(arg1,arg2,arg3,arg4));
    }
    else if (number == 2)
    {
      int idx1, idx2;
      for (int j=0; j < 3; j++)
      {
	//if (elementNeigh[(3*i)+j] != -1)
	if (!IGIPtr_->boundaryIndex(elementNeigh[j]))
	{
	  idx1 = IGIPtr_->getGlobalIndexFromEntityPointer(elementNeigh[j]);
	  break;
	}
      }
      for (int j=0; j < 3; j++)
      {
	//if ((elementNeigh[(3*i)+j] != -1) && (elementNeigh[(3*i)+j] != idx1))
	if ((!IGIPtr_->boundaryIndex(elementNeigh[j])) && (IGIPtr_->getGlobalIndexFromEntityPointer(elementNeigh[j]) != idx1))
	{
	  //idx2 = elementNeigh[(3*i)+j];
	  idx2 = IGIPtr_->getGlobalIndexFromEntityPointer(elementNeigh[j]);
	  break;
	}
      }

      //Must apply the MINMOD limiter on each component:
      arg1 = saturationGrad_[i].u_getElement_0u(0);
      arg2 = saturationGrad_[idx1].u_getElement_0u(0);
      arg3 = saturationGrad_[idx2].u_getElement_0u(0);
    
      gradSat.u_changeElement_1u(0,ml_MINMOD_3ARG(arg1,arg2,arg3));

      arg1 = saturationGrad_[i].u_getElement_0u(1);
      arg2 = saturationGrad_[idx1].u_getElement_0u(1);
      arg3 = saturationGrad_[idx2].u_getElement_0u(1);

      gradSat.u_changeElement_1u(1,ml_MINMOD_3ARG(arg1,arg2,arg3));
    }
    else if (number == 1)
    {
      int idx1;
      for (int j=0; j < 3; j++)
      {
	//if (elementNeigh[(3*i)+j] != -1)
	if (!IGIPtr_->boundaryIndex(elementNeigh[j]))
	{
	  //idx1 = elementNeigh[(3*i)+j];
	  idx1 = IGIPtr_->getGlobalIndexFromEntityPointer(elementNeigh[j]);
	  break;
	}
      }

      //Must apply the MINMOD limiter on each component:
      arg1 = saturationGrad_[i].u_getElement_0u(0);
      arg2 = saturationGrad_[idx1].u_getElement_0u(0);
    
      gradSat.u_changeElement_1u(0,ml_MINMOD_2ARG(arg1,arg2));

      arg1 = saturationGrad_[i].u_getElement_0u(1);
      arg2 = saturationGrad_[idx1].u_getElement_0u(1);

      gradSat.u_changeElement_1u(1,ml_MINMOD_2ARG(arg1,arg2));
    }

    saturationGradNew_.u_changeElement_1u(i,gradSat);				 
  }
}



void SFCentralUpwind::ml_generateLimitersForStructuredGrid(const SFSaturation& saturation, const int& NX, const int& NY, const double& DX, const double& DY)
{
  //@HAF: Must be more thoroughly checked??? Is something WRONG with directions etc. pr. 15.02.2010 due to different orderings of JMY and Dune...???

  //Must fill for saturationGradNew_

  //We first treat the limiters in the x-direction:
  int k=0;
  for (int j=0; j < NY; j++)
  {
    for (int i=1; i < NX-1; i++)
    {
      k = i + j*NX;
      //m_SLimiterX.u_changeElement_1u(k,(1.0/DX)*ul_MINMODLimiter(k,k+1,k-1));
      RnPoint gradP = saturationGradNew_[k];
      gradP.u_changeElement_1u(0, (1.0/DX)*ml_MINMODLimiterStructuredGrid(saturation,k,k+1,k-1));
      saturationGradNew_.u_changeElement_1u(k, gradP);
    }
  }


  //We then treat the limiter in the y-direction:
  for (int i=0; i < NX; i++)
  {
    for (int j=1; j < NY-1; j++)
    {
      k = i + j*NX;
      //m_SLimiterY.u_changeElement_1u(k,(1.0/DY)*ul_MINMODLimiter(k,k+NX,k-NX));
      RnPoint gradP = saturationGradNew_[k];
      //gradP.u_changeElement_1u(1,(1.0/DY)*ml_MINMODLimiterStructuredGrid(saturation,k,k+NX,k-NX));//"Sophus-JMY" solution!!!
      gradP.u_changeElement_1u(1,(1.0/DY)*ml_MINMODLimiterStructuredGrid(saturation,k,k+NX,k-NX));//@HAF: Correct Dune-solution???
      saturationGradNew_.u_changeElement_1u(k, gradP);
    }
  }

}



double SFCentralUpwind::ml_MINMODLimiterStructuredGrid(const SFSaturation& saturation, const int& ID, const int& IDdirP1, const int& IDdirM1) const
{

  //This function computes the MINMOD derivative limiter for the Kurganov-Tadmor scheme:
  // NB: The division by the characteristic grid-length is NOT performed inside this routine !!!


  //******The choice of theta *************************************
  double theta = 1.0;
  //***************************************************************

  double MM, MM_1, MM_2, MM_3;

  MM_1 = theta*(saturation[IDdirP1]-saturation[ID]);
  MM_2 = 0.5*(saturation[IDdirP1]-saturation[IDdirM1]);
  MM_3 = theta*(saturation[ID]-saturation[IDdirM1]);

  if ((MM_1 > 0.0) && (MM_2 > 0.0) &&(MM_3 > 0.0))
  {
    double help = (MM_1 <= MM_2) ? MM_1 : MM_2; 
    MM = (help <= MM_3) ? help : MM_3; 
  }
  else if ((MM_1 < 0.0) && (MM_2 < 0.0) &&(MM_3 < 0.0))
  {
    double help = (MM_1 >= MM_2) ? MM_1 : MM_2; 
    MM = (help >= MM_3) ? help : MM_3; 
  }
  else
  {
    MM = 0.0;
  }
   
  return MM;
}



double SFCentralUpwind::ml_MINMOD_2ARG(const double& arg1, const double& arg2)
{
  //This is the MINMOD limiter function for two arguments

  double ret = 0.0;

  if ((arg1 > 0.0) && (arg2 > 0.0))
  {
    ret = (arg1 <= arg2) ? arg1 : arg2; 
  }
  else if ((arg1 < 0.0) && (arg2 < 0.0))
  {
    ret = (arg1 >= arg2) ? arg1 : arg2; 
  }
  else
  {
   ret = 0.0;
  }
 
  return ret;
}




double SFCentralUpwind::ml_MINMOD_3ARG(const double& arg1, const double& arg2, const double& arg3)
{
  //*****************************************************************************************************
  //This function computes the MINMOD limiter function based on 3 arguments.
  //*****************************************************************************************************

  double ret = 0.0;

  if ((arg1 > 0.0) && (arg2 > 0.0) && (arg3 > 0.0))
  
  {    //should return the minimum:
    ret = arg1;
    if (arg2 < ret) ret = arg2; 
    if (arg3 < ret) ret = arg3; 
  }
  else if ((arg1 < 0.0) && (arg2 < 0.0) && (arg3 < 0.0))
  {
    //should return the maximum:
    ret = arg1;
    if (arg2 > ret) ret = arg2; 
    if (arg3 > ret) ret = arg3; 
  }
  else
  {
   ret = 0.0;
  }

  return ret;
}




double SFCentralUpwind::ml_MINMOD_4ARG(const double& arg1, const double& arg2, const double& arg3, const double& arg4)
{
  //*****************************************************************************************************
  //This function computes the MINMOD limiter function based on 4 arguments.
  //*****************************************************************************************************

  double ret = 0.0;

  if ((arg1 > 0.0) && (arg2 > 0.0) && (arg3 > 0.0) && (arg4 > 0.0))
  
  {    //should return the minimum:
    ret = arg1;
    if (arg2 < ret) ret = arg2; 
    if (arg3 < ret) ret = arg3; 
    if (arg4 < ret) ret = arg4; 
  }
  else if ((arg1 < 0.0) && (arg2 < 0.0) && (arg3 < 0.0) && (arg4 < 0.0))
  {
    //should return the maximum:
    ret = arg1;
    if (arg2 > ret) ret = arg2; 
    if (arg3 > ret) ret = arg3; 
    if (arg4 > ret) ret = arg4; 
  }
  else
  {
   ret = 0.0;
  }

  return ret;
}


int SFCentralUpwind::ml_getNonPeriodicBCType(const SFPMPhysicalFieldProperties& PMPhysProp, const IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt, const int& k)
{
  //**********************************************************
  //We only consider Neumann or Dirichlet BC's here !!!
  //I. e. either alpha1 or alpha2 MUST be zero.
  //NB: Neumann returns 10 and Dirichlet returns 20.
  //**********************************************************

  //*****************************************************************
  //NOTE: This code is yet ONLY for 2D !!!
  //*****************************************************************


  double TOL = 0.0000001;

  int BCType = -1;

  //We must first locate the position the position of the global edge number
  //on the boundary. I. e. on which boundary does it belong.
  RnPoint coord = IGIPtr_->getCentroid<1>(ihatPt, k);
  double XC = coord.u_getElement_0u(0);
  double YC = coord.u_getElement_0u(1);

  if ((YC < TOL) && (YC > -TOL))
  {
    //We are at the bottom boundary:
    BCType = ((PMPhysProp.getSimulationDataPtr()->IMPESd_.BoundaryCondBottom_alpha1 < TOL) && (PMPhysProp.getSimulationDataPtr()->IMPESd_.BoundaryCondBottom_alpha1 > -TOL)) ? 10 : 20; 
  }
  else if ((YC < (PMPhysProp.getSimulationDataPtr()->td_.Y_Extent+TOL)) && (YC > (PMPhysProp.getSimulationDataPtr()->td_.Y_Extent-TOL)))
  {
    //We are at the top boundary:
    BCType = ((PMPhysProp.getSimulationDataPtr()->IMPESd_.BoundaryCondTop_alpha1 < TOL) && (PMPhysProp.getSimulationDataPtr()->IMPESd_.BoundaryCondTop_alpha1 > -TOL)) ? 10 : 20; 
  }
  else if ((XC < (TOL)) && (XC > (-TOL)))
  {
    //We are at the left boundary:
    BCType = ((PMPhysProp.getSimulationDataPtr()->IMPESd_.BoundaryCondLeft_alpha1 < TOL) && (PMPhysProp.getSimulationDataPtr()->IMPESd_.BoundaryCondLeft_alpha1 > -TOL)) ? 10 : 20; 
  }
  else if ((XC < (PMPhysProp.getSimulationDataPtr()->td_.X_Extent+TOL)) && (XC > (PMPhysProp.getSimulationDataPtr()->td_.X_Extent-TOL)))
  {
    //We are at the right boundary:
    BCType = ((PMPhysProp.getSimulationDataPtr()->IMPESd_.BoundaryCondRight_alpha1 < TOL) && (PMPhysProp.getSimulationDataPtr()->IMPESd_.BoundaryCondRight_alpha1 > -TOL)) ? 10 : 20; 
  }
  
  return BCType;
}


double SFCentralUpwind::ml_getBCDirichletCond(const SFPMPhysicalFieldProperties& PMPhysProp, const IRISDuneGridInterface<DuneGridType>::ElementPointer ihatPt, const int& k)
{
  //*****************************************************************

  //NOTE: This code is yet ONLY for 2D !!!

  //*****************************************************************


  double TOL = 0.0000001;

  double cond = -77.0;

  //We must first locate the position the position of the global edge number
  //on the boundary. I. e. on which boundary does it belong.
  RnPoint coord = IGIPtr_->getCentroid<1>(ihatPt, k);
  double XC = coord.u_getElement_0u(0);
  double YC = coord.u_getElement_0u(1);

  if ((YC < TOL) && (YC > -TOL))
  {
    //We are at the bottom boundary:
    cond = PMPhysProp.getSimulationDataPtr()->IMPESd_.BoundaryCondBottom_g;
  }
  else if ((YC < (PMPhysProp.getSimulationDataPtr()->td_.Y_Extent+TOL)) && (YC > (PMPhysProp.getSimulationDataPtr()->td_.Y_Extent-TOL)))
  {
    //We are at the top boundary:
    cond = PMPhysProp.getSimulationDataPtr()->IMPESd_.BoundaryCondTop_g;
  }
  else if ((XC < (TOL)) && (XC > (-TOL)))
  {
    //We are at the left boundary:
    cond = PMPhysProp.getSimulationDataPtr()->IMPESd_.BoundaryCondLeft_g;
  }
  else if ((XC < (PMPhysProp.getSimulationDataPtr()->td_.X_Extent+TOL)) && (XC > (PMPhysProp.getSimulationDataPtr()->td_.X_Extent-TOL)))
  {
    //We are at the right boundary:
    cond = PMPhysProp.getSimulationDataPtr()->IMPESd_.BoundaryCondRight_g;
  }

  return cond;
}
