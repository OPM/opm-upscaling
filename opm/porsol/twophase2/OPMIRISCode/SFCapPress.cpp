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
#include "SFCapPress.h"

#include <iostream>
#include <iomanip>

//---------------------- Directives -------------------------------------------
using namespace std;  // This leads to that all names in the namespace std may be 
                      // used without the qualifier "std::".





SFCapPress::SFCapPress()  : m_variant(0), m_RowLength(0)
{

}


void SFCapPress::readCapPressDataInFormulaFormat(istream& inputStream, const int& numbRockTypes)
{
  //****************************************************************************************
  //This function reads the rel. perm. table for the case that the cap. press is represented
  //in formula format. The values are simply stored in the 1D Dune::array values.
  //****************************************************************************************

  m_variant = 0;
  m_RowLength = 14;
  m_constants = new double[numbRockTypes*m_RowLength];

  for (int j=0; j < numbRockTypes; j++)
  {
    int start=j*m_RowLength;
    int i=start;
    //    while (i < start+m_RowLength && inputStream >> m_constants[i]) i++;
    while (i < start+m_RowLength)
    {
      inputStream >> m_constants[i]; 
      i++;
    }
  }
}



double SFCapPress::getCapPressData(const int& capPressType, const double& sw, const int& rockType) const
{
  /*
    sw is the water saturation.
    capPressType=0 : oil-water capillary pressure
    capPressType=1 : oil-gas capillary pressure
  */

  //*************************************************************************************
  //This routine is basically implemented for two-phase flow, but can easily be extended
  //to three-phase flow.
  //*************************************************************************************
  
  double capp;

  if (m_variant == 0)
  {
    int start=rockType*m_RowLength;
    int step = (capPressType == 0) ? 0 : 6;
    double CL = m_constants[step+capPressType+start];
    //    double CNOL = m_constants[step+capPressType+1+start];
    //double CR = m_constants[step+capPressType+2+start];
    //double EL = m_constants[step+capPressType+3+start];
    //double ER = m_constants[step+capPressType+4+start];
    //double SL = m_constants[step+capPressType+5+start];
    //double SR = m_constants[step+capPressType+6+start];

    //NB: Note that SL and SR are left and right asymptotes of the PC-curve, respectively. 

    double TOL = 0.000001;
    //double S = sw;
    //double A = S - SL - TOL;
    //double B = SR - S - TOL;
    
    //if (A < TOL) A = TOL;
    //if (B < TOL) B = TOL;

    //capp = CL*pow(A,-EL) - CR*pow(B,-ER) + CNOL;


  //Note here we use a LOGARITHMIC cap. press profile, in order to investigate the effects of the
  //diffusion in the displacement example !!!
    //**********************************************************************
    if (sw < TOL)//@HAF: IMPORTANT!!!
    {
      capp = -CL*log(TOL);//NB: JUST A TEST !!!
    }
    else
    {
      capp = -CL*log(sw);//NB: JUST A TEST !!!
    }
    //    capp = -1000.0*CNOL*S;//NB: JUST A TEST !!!
    //    cout << "capp= " << capp << " CL= " << CL << " rockType= " << rockType << " index= " << step+capPressType+start << endl;
    //cout << "step= " << step << " start= " << start <<  " m_RowLength= " << m_RowLength << endl;
  }

  return capp;
}




double SFCapPress::getCapPressDataPartial(const int& capPressType, const double& sw, const int& rockType) const
{
  /*
    sw is the water saturation.
    capPressType=0 : oil-water capillary pressure
    capPressType=1 : oil-gas capillary pressure
  */

  //*************************************************************************************
  //This routine is basically implemented for two-phase flow, but can easily be extended
  //to three-phase flow.
  //*************************************************************************************
  
  double capp, CL;

  if (m_variant == 0)
  {
    int start=rockType*m_RowLength;
    int step = (capPressType == 0) ? 0 : 6;
    CL = m_constants[step+capPressType+start];
    //double CNOL = m_constants[step+capPressType+1+start];
    //double CR = m_constants[step+capPressType+2+start];
    //double EL = m_constants[step+capPressType+3+start];
    //double ER = m_constants[step+capPressType+4+start];
    //double SL = m_constants[step+capPressType+5+start];
    //double SR = m_constants[step+capPressType+6+start];

    //NB: Note that SL and SR are left and right asymptotes of the PC-curve, respectively. 

    double TOL = 0.000001;
    //double S = sw;
    //double A = S - SL - TOL;
    //double B = SR - S - TOL;
    
    //if (A < TOL) A = TOL;
    //if (B < TOL) B = TOL;

    //    capp = -CL*EL*pow(A,-EL-1.0) - CR*ER*pow(B,-ER-1.0) + CNOL;
    capp = -(CL/(sw+TOL));
    //    capp = -1.0; //TESTING with constant diffusion!!!
  }

  //#ifdef DisplacementProblem
  //Note here we use a LINEAR cap. press profile, just to investigate the effects of the
  //diffusion in the displacement example !!!
  //  capp = -1.0;
  //capp = -10.0;
  //capp = 0.0;

  //Note here we use a LOGARITHMIC cap. press profile, in order to investigate the effects of the
  //diffusion in the displacement example !!!
  //********************************************************************************
  //NB: Dette maa endres naar vi legger inn heterogenitet !!!!!!!!!!
  //********************************************************************************
  //double TOL = 0.000001;
  //  double constant = 2.0;
  //capp = -(constant/sw);//Note: Useless when used in saturation experiments with sw=0 in parts of the domain.
  //  capp = -(constant/(sw+TOL));
  //#endif //DisplacementProblem

  //***********************************************************
  //@HAF:
  //Special coding for so called "semi-analytical PC-test"(van Duijn & de Neef)
  cout << "Program is stopped. Do we actually need this function, i.e. the function 'getCapPressDataPartial' here ??? If so: we must understand why!!!!!" << endl;
  exit(0);
  //***********************************************************

  return capp;
}




double SFCapPress::getCapPressDataPartial2(const int& capPressType, const double& sw, const int& rockType) const
{
  /*
    sw is the water saturation.
    capPressType=0 : oil-water capillary pressure
    capPressType=1 : oil-gas capillary pressure
  */

  //*************************************************************************************
  //This routine is basically implemented for two-phase flow, but can easily be extended
  //to three-phase flow.
  //*************************************************************************************
  
  double capp, CL;
  //cout << "rockType= " << rockType << " m_RowLength= " << m_RowLength << endl;
  if (m_variant == 0)
  {
    int start=rockType*m_RowLength;
    int step = (capPressType == 0) ? 0 : 6;
    CL = m_constants[step+capPressType+start];
    //double CNOL = m_constants[step+capPressType+1+start];
    //double CR = m_constants[step+capPressType+2+start];
    //double EL = m_constants[step+capPressType+3+start];
    //double ER = m_constants[step+capPressType+4+start];
    //double SL = m_constants[step+capPressType+5+start];
    //double SR = m_constants[step+capPressType+6+start];

    //NB: Note that SL and SR are left and right asymptotes of the PC-curve, respectively. 

    double TOL = 0.000001;
    //double S = sw;
    //double A = S - SL - TOL;
    //double B = SR - S - TOL;
    
    //if (A < TOL) A = TOL;
    //if (B < TOL) B = TOL;

    //    capp = -CL*EL*pow(A,-EL-1.0) - CR*ER*pow(B,-ER-1.0) + CNOL;
    capp = CL/((sw+TOL)*(sw+TOL));
    //    capp = -1.0; //TESTING with constant diffusion!!!
  }

  //#ifdef DisplacementProblem
  //Note here we use a LINEAR cap. press profile, just to investigate the effects of the
  //diffusion in the displacement example !!!
  //  capp = -1.0;
  //capp = -10.0;
  //capp = 0.0;

  //Note here we use a LOGARITHMIC cap. press profile, in order to investigate the effects of the
  //diffusion in the displacement example !!!
  //********************************************************************************
  //NB: Dette maa endres naar vi legger inn heterogenitet !!!!!!!!!!
  //********************************************************************************
  //double TOL = 0.000001;
  //  double constant = 2.0;
  //capp = -(constant/sw);//Note: Useless when used in saturation experiments with sw=0 in parts of the domain.
  //  capp = -(constant/(sw+TOL));
  //#endif //DisplacementProblem

  return capp;
}




double SFCapPress::getInverseCapPressData(const int& capPressType, const double& pcv, const int& rockType) const
{
  /*
    sw is the water saturation.
    capPressType=0 : oil-water capillary pressure
    capPressType=1 : oil-gas capillary pressure
  */

  //*************************************************************************************
  //This routine is basically implemented for two-phase flow, but can easily be extended
  //to three-phase flow.
  //*************************************************************************************
  
  double invCapp = 0.0;

  if (m_variant == 0)
  {
    int start=rockType*m_RowLength;
    int step = (capPressType == 0) ? 0 : 6;
    double CL = m_constants[step+capPressType+start];

    double TOL = 0.000001;

  //Note here we use a LOGARITHMIC cap. press profile, in order to investigate the effects of the
  //diffusion in the displacement example !!!
    //**********************************************************************
    if (fabs(CL) > (0.01*TOL)) invCapp = exp(-pcv/CL)-TOL;//NB: JUST A TEST !!!
  }

  return invCapp;
}





double SFCapPress::getOutletSaturation(const int& capPressType, const int& rockTypeEnd) const
{
  //*************************************************************************************
  //This routine is basically implemented for two-phase flow
  //(and thus returns the water saturation)
  //,but may be extended to three-phase flow.
  //*************************************************************************************
  
  //We here solve the following nonlinear equation P_C(s_w) = 0.0.,
  // which is valid at the outlet.

  double TOL = 0.000001;//Can be discussed...
  double sw = -1.0;
  double DS = -1.0;

  //Must construct an initial guess (difficult part):
  sw = ml_constructInitialGuessForOutletSaturation(capPressType,rockTypeEnd);

  bool converged = false;
  int i = 0;
  int maxIter = 100;//Can be discussed!!!
  //The Newton iterations:
  while ((!converged) && (i < maxIter))
  {
    DS = -getCapPressData(capPressType, sw, rockTypeEnd)/getCapPressDataPartial(capPressType, sw, rockTypeEnd);
    sw += DS;
    if (fabs(getCapPressData(capPressType, sw, rockTypeEnd)) < TOL) converged= true;
    i++;
  }

  return sw;
}



SFCapPress::~SFCapPress() 
{
 
}


//-----------------------------------------------------------------
//-----------------Private---functions-----------------------------
//-----------------------------------------------------------------

double SFCapPress::ml_constructInitialGuessForOutletSaturation(const int& capPressType, const int& rockTypeEnd) const
{ 
  //We here construct a hopefully good initial guess for the
  //nonlinear equation P_C(s_w) = 0.0.,
  // which is valid at the outlet.
  //The algorithm is based on bisection:

  bool OK = false;

  //We must obviously have that 0.0 <= sw <= 1.0.
  double TOL = 0.000001;
  double sw_l = TOL;
  double sw_u = 1.0;
  double sw_m = 0.5*(sw_u + sw_l);

  double TOLS = 0.01;//Tolerance for initial guess. Can be discussed???
  //We check the start point:
  double f_sw_l = getCapPressData(capPressType, sw_l, rockTypeEnd);
  double f_sw_u = getCapPressData(capPressType, sw_u, rockTypeEnd);
  double f_sw_m;
  //  cout << "f_sw_l= " << f_sw_l << endl;
  //cout << "f_sw_u= " << f_sw_u << endl;
  //exit(0);
  //Make a check of the starting point:
  if (f_sw_l*f_sw_u > 0.0)
  {
    cout << "Something is wrong with the initial interval!!!" << endl;
    exit(0);
  }

  if ((fabs(f_sw_l)) < TOLS) 
  {
    OK = true;
    sw_m = sw_l;
  }

  if (!OK)
  {
    if ((fabs(f_sw_u)) < TOLS) 
    {
      OK = true;
      sw_m = sw_u;
    }
  }

  if (!OK)
  {
    f_sw_m = getCapPressData(capPressType, sw_m, rockTypeEnd);
    if ((fabs(f_sw_m)) < TOLS) OK = true;
  }


  //Now comes the iterative part of the algo.:
  int i = 0;
  int maxIter = 100;//Can be discussed!!!
  while ((!OK) && (i < maxIter))
  {
    if (f_sw_l*f_sw_m < 0.0)
    {
      sw_u = sw_m;
    }
    else
    {
      sw_l = sw_m;
    }
    f_sw_l = getCapPressData(capPressType, sw_l, rockTypeEnd);
    f_sw_u = getCapPressData(capPressType, sw_u, rockTypeEnd);
    sw_m = 0.5*(sw_l + sw_u);
    f_sw_m = getCapPressData(capPressType, sw_m, rockTypeEnd);
    if ((fabs(f_sw_m)) < TOLS) OK = true;
    i++;
  }

  return sw_m;
}
