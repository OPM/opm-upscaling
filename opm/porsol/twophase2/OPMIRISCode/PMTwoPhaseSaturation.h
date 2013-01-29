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



#ifndef PMTwoPhaseSaturationH
#define PMTwoPhaseSaturationH



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

#include <opm/porsol/twophase2/OPMIRISCode/SFRelPerm.h>
#include <opm/porsol/twophase2/OPMIRISCode/SFPMPhysicalFieldProperties.h>
#include <opm/porsol/twophase2/OPMIRISCode/SFCentralUpwind.h>


//---------------------  Constants --------------------------------------------




//---------------------- Directives -------------------------------------------
using namespace std;  // This leads to that all names in the namespace std may be 
                      // used without the qualifier "std::".


//*********************************************************************************************



template <class FVMExplicit>   
class PMTwoPhaseSaturation
{
 public:
  //--------------------------------------------------------------------------
  //---------------------- Constructors --------------------------------------
  //--------------------------------------------------------------------------
  PMTwoPhaseSaturation();
  PMTwoPhaseSaturation(const FVMExplicit& FVMEMenth);

  //--------------------------------------------------------------------------
  //---------------------- Generators-----------------------------------------
  //--------------------------------------------------------------------------
  void setup(const IRISDuneGridInterface<DuneGridType>& Igrid, const SFSaturation& saturation, const PMTwoPhaseSimulationData& simulationData, const SFPMPhysicalFieldProperties& PMpp);//@HAF: OK med referanse paa SFSaturation???
  void computeForwardInTime(const double& timeInterval, const TimeIntegrationType& integrationType, const SFVelocity& velocity);//@HAF: SFVelocity er vel en typedef'et Array<double>. OK med referanse???

  //--------------------------------------------------------------------------
  //---------------------- Observers -----------------------------------------
  //--------------------------------------------------------------------------
  const SFSaturation& getSaturation() const;//@HAF: OK med aa levere en referanse for SFSaturation???

  //@HAF: Trenger vi en destructor???


  //-----------------------------------------------------------------------------------

  private : 
  //-----------------Useful local functions--------------------------------------------
  void ml_computeForwardInTimeEuler(const double& timeInterval, const SFVelocity& velocity);//@HAF. SFVelocity er vel en typedef'et Array<double>. OK med referanse???
  void ml_computeForwardInTimeSecondOrderEuler(const double& timeInterval, const SFVelocity& velocity);
  void ml_computeTimeStep(const SFVelocity& velocity);//@HAF: This routine is basically implemented within in FVMExplicit (see the implementation below).
  void ml_computePracticalTimeStep(const double& timeInterval, double& actualDeltat, int& numbOfTimeSteps);
  

  //-----------------Data structure----------------------------------------------------

  FVMExplicit FVMMeth_;
  SFSaturation saturation_;//@HAF: HUSK aat SFSaturation bare er en "typdefet" Array<double> eller lignende...
  double deltatIMPES_;//the "outer" IMPES time step size.
  double deltat_;//the "inner" numerical time step size (needed in the discr. of the saturation eq.).
  double minGridSize_;//A typical minimum grid size quantity used in time-step requirements.
  SFRelPerm rPerm_;
  SFPMPhysicalFieldProperties PMpp_;
};




/*---------------------------------------------------*/
/*              Constructors                         */
/*---------------------------------------------------*/

template <class FVMExplicit>
PMTwoPhaseSaturation<FVMExplicit>::PMTwoPhaseSaturation() 
  : FVMMeth_(), saturation_(), deltatIMPES_(0.0), deltat_(0.0), minGridSize_(0.0), rPerm_(), PMpp_()
{

}

template <class FVMExplicit>
PMTwoPhaseSaturation<FVMExplicit>::PMTwoPhaseSaturation(const FVMExplicit& FVMEMeth) 
  : FVMMeth_(FVMEMeth), saturation_(), deltatIMPES_(0.0), deltat_(0.0), minGridSize_(0.0), rPerm_(), PMpp_()
{

}


/*---------------------------------------------------*/
/*              Generators                           */
/*---------------------------------------------------*/
template <class FVMExplicit>
void PMTwoPhaseSaturation<FVMExplicit>::setup(const IRISDuneGridInterface<DuneGridType>& Igrid, const SFSaturation& saturation, const PMTwoPhaseSimulationData& simulationData, const SFPMPhysicalFieldProperties& PMpp)
{
  saturation_ = saturation;
  deltatIMPES_ = simulationData.IMPESd_.delta_T;
  FVMMeth_.setup(&Igrid);//@HAF: NOTE: The signature of this setup must be no different for a class which have central-upwind scheme, or one-point upwind etc. etc. etc.

  //ifstream relPermInStream(relPermFormula_In_file);"relPermFormula.dta";
  ifstream relPermInStream("relPermFormula.dta");//@HAF: Should be with the file name as a parameter to the setup routine???
  rPerm_.readRelPermDataInFormulaFormat(relPermInStream,simulationData.IMPESd_.NumbRockTypes);
  relPermInStream.close();
  cout << endl << "File \"" << "relPermFormula.dta" << "\" is read...." << endl;

  PMpp_ = PMpp;

  deltat_ = 0.0;//NB: Will be computed when needed in computeTimeStep
  minGridSize_ = FVMMeth_.getMinGridSize();
}


template <class FVMExplicit>
void PMTwoPhaseSaturation<FVMExplicit>::computeForwardInTime(const double& timeInterval, const TimeIntegrationType& integrationType, const SFVelocity& velocity)//@HAF: SFVelocity er vel en typedef'et Array<double>. OK med referanse???
{
  //@HAF: TimeIntegrationType should be an enum type.
  //Clearly the if-statement could be made longer with other
  //time-integration alternatives.

  if (integrationType == _Euler)
  {
    ml_computeForwardInTimeEuler(timeInterval, velocity);
  }
  else if (integrationType == _SecondOrderEuler)
  {
    ml_computeForwardInTimeSecondOrderEuler(timeInterval, velocity);
  }
}


/*---------------------------------------------------*/
/*              Observers                            */
/*---------------------------------------------------*/

template <class FVMExplicit>
const SFSaturation& PMTwoPhaseSaturation<FVMExplicit>::getSaturation() const
{
  return saturation_;
}



/*---------------------------------------------------*/
/*             Local functions                       */
/*---------------------------------------------------*/
template <class FVMExplicit>
void PMTwoPhaseSaturation<FVMExplicit>::ml_computeForwardInTimeEuler(const double& timeInterval, const SFVelocity& velocity)//@HAF. SFVelocity er vel en typedef'et Array<double>. OK med referanse???
{
  Array<double> rhs(saturation_.u_getSize_0u());//@HAF: Is this allocation becomes too costly, we should consider putting rhs in the data structure of this class. The type TimeIntegrationType should then be employed in the impl...

  double actualDeltat = -1.0;
  int numbOfTimeSteps = 1;
  ml_computeTimeStep(velocity);
  ml_computePracticalTimeStep(timeInterval,actualDeltat,numbOfTimeSteps);

  for (int i=0; i < numbOfTimeSteps; i++)
  {
    FVMMeth_.computeFluxOutOfAllCells(rPerm_,PMpp_,saturation_,velocity,rhs);
    rhs *= actualDeltat;
    saturation_ += rhs;
  }
}


template <class FVMExplicit>
void PMTwoPhaseSaturation<FVMExplicit>::ml_computeForwardInTimeSecondOrderEuler(const double& timeInterval, const SFVelocity& velocity)//@HAF. SFVelocity er vel en typedef'et Array<double>. OK med referanse???
{
  Array<double> intermediateSolution = saturation_;
  Array<double> rhs(saturation_.u_getSize_0u(), 0.0);//@HAF: Is this allocation becomes too costly, we should consider putting rhs and intermediateSolution in the data structure of this class. The type TimeIntegrationType should then be employed in the impl...

  double actualDeltat = -1.0;
  int numbOfTimeSteps = 1;
  ml_computeTimeStep(velocity);
  ml_computePracticalTimeStep(timeInterval,actualDeltat,numbOfTimeSteps);

  for (int i=0; i < numbOfTimeSteps; i++)
  {
    cout << "********************************************************" << endl;
    cout << "numbOfTimeSteps= " << numbOfTimeSteps << " i= " << i << endl;
    //1st step of algorithm:
    FVMMeth_.computeFluxOutOfAllCells(rPerm_,PMpp_,saturation_,velocity,rhs);
    cout << "rhs1= " << rhs[0] << " , " << rhs[1] << " , " << rhs[2] << " , " << rhs[3] << endl;
    cout << "-------------------------------------------------------------------" << endl;
    rhs *= actualDeltat;
    intermediateSolution += rhs;
    cout << " actualDeltat= " << actualDeltat << " intermediateSolution= " << intermediateSolution[0] << " , " << intermediateSolution[1] << " , " << intermediateSolution[2] << " , " << intermediateSolution[3] << endl;
    cout << "-------------------------------------------------------------------" << endl;

    //2nd step of algorithm:
    FVMMeth_.computeFluxOutOfAllCells(rPerm_,PMpp_,intermediateSolution,velocity,rhs);
    cout << "rhs2= " << rhs[0] << " , " << rhs[1] << " , " << rhs[2] << " , " << rhs[3] << endl;
    cout << "-------------------------------------------------------------------" << endl;
    //exit(0);
    rhs *= (0.5*actualDeltat);
    saturation_ *= 0.5;
    intermediateSolution *= 0.5;
    saturation_ += intermediateSolution;
    saturation_ += rhs;
    intermediateSolution = saturation_;//Initialization for the next time step
  }
}


template <class FVMExplicit>
void PMTwoPhaseSaturation<FVMExplicit>::ml_computeTimeStep(const SFVelocity& velocity)//@HAF. SFVelocity er vel en typedef'et Array<double>. OK med referanse???
{
  double maxV = 0.0;
  double maxVelocity = 0.0;
  for (int i=0; i < velocity.u_getSize_0u(); i++)
  {
    maxV = fabs(velocity[i]);
    if (maxV > maxVelocity) maxVelocity = maxV;
  }
  //deltat_ = FVMMeth_.computeTimeStep(minGridSize_, maxVelocity);//@HAF: OK for problem without capillary pressure
  deltat_ = FVMMeth_.computeTimeStep(minGridSize_, 2.0*1250.0);//@HAF: 50 times 75 grid 10.03.2010. "Homogeneous EX2"
  //deltat_ = FVMMeth_.computeTimeStep(minGridSize_, 2.0*2500.0);//@HAF: 50 times 75 grid 12.03.2010. "Heterogeneous EX2"
  cout << "deltat_= " << deltat_ << " deltatIMPES_= " << deltatIMPES_ << " minGridSize_= " << minGridSize_ << " maxVelocity= " << maxVelocity << endl; 

  //A final important check:
  //if (deltat_ > deltatIMPES_)
  //{
    //Here we should case an exception (but for the time-being the primitive solution is chosen)
    //cout << "Inconsistency between inner an outer time-steps in the IMPES algo. Program is stopped." << endl;
    //exit(0);
  //}
}


template <class FVMExplicit>
void PMTwoPhaseSaturation<FVMExplicit>::ml_computePracticalTimeStep(const double& timeInterval, double& actualDeltat, int& numbOfTimeSteps){
  Boolean OK = false;
    while (!OK)
   {
   if ((timeInterval/numbOfTimeSteps) <= deltat_) 
  {
    actualDeltat = timeInterval/numbOfTimeSteps;
    OK = true;
  }
  if (!OK) numbOfTimeSteps++;
  }

}


#endif //PMTwoPhaseSaturationH

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




