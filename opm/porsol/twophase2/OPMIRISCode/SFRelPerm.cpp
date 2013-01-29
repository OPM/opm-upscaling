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
#include "SFRelPerm.h"

#include <iostream>
#include <iomanip>

//---------------------- Directives -------------------------------------------
using namespace std;  // This leads to that all names in the namespace std may be 
                      // used without the qualifier "std::".





SFRelPerm::SFRelPerm()  : m_variant(0),  m_RowLength(0)
{
}


void SFRelPerm::readRelPermDataInFormulaFormat(istream& inputStream, const int& numbRockTypes)
{
  //****************************************************************************************
  //This function reads the rel. perm. table for the case that the rel. perm is represented
  //in formula format. The values are simply stored in the 1D Dune::array m_constants.
  //****************************************************************************************

  m_variant = 0;
  m_RowLength = 6;
  m_constants = new double[numbRockTypes*m_RowLength];

  for (int j=0; j < numbRockTypes; j++)
  {
    int start=j*m_RowLength;
    int i=start;
    //    while (i < start+m_RowLength && inputStream >> m_constants[i])
    while (i < start+m_RowLength)
    {
      inputStream >> m_constants[i];
      i++;
    }
  }
}


double SFRelPerm::getRelPermData(const int& phase, const double& sw, const int& rockType) const
{
  /*
    sw is the water saturation.
    phase=0 : oil
    phase=1 : water
    phase=3 : gas
  */

  //*************************************************************************************
  //This routine is basically implemented for two-phase flow, but can easily be extended
  //to three-phase flow.
  //*************************************************************************************
  
  double relp;

  if (m_variant == 0)
  {
    int start=rockType*m_RowLength;
    double EC = m_constants[phase+start];
    double SC = m_constants[(phase+3)+start];

    double S = (phase == 1) ? sw : 1.0 - sw;
  
    double A = (S - SC)/(1.0 - SC);

    double TOL = 0.000001;

    if (A < TOL)
    {
      //      relp = TOL;
      relp = 0.0;
    }
    if (A > 1.0)
    {
      relp = 1.0;
    }
    else
    {
      relp = pow(A,EC);
    }
  }

  return relp;
}




double SFRelPerm::getRelPermDataPartial(const int& phase, const double& sw, const int& rockType) const
{
  /*
    sw is the water saturation.
    phase=0 : oil
    phase=1 : water
    phase=3 : gas
  */

  //*************************************************************************************
  //This routine is basically implemented for two-phase flow, but can easily be extended
  //to three-phase flow.
  //*************************************************************************************

  //*************************************************************************************
  //NBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNB

  //The mechanism with signDir below may NOT work for three-phase flow. It MUST be checked !!!

  //NBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNB
  //*************************************************************************************
  
  double relp;

  if (m_variant == 0)
  {
    int start=rockType*m_RowLength;
    double EC = m_constants[phase+start];
    double SC = m_constants[(phase+3)+start];

    double S = (phase == 1) ? sw : 1.0 - sw;
    double signDer = (phase == 1) ? 1.0 : -1.0;
    
    double A = (S - SC)/(1.0 - SC);

    double TOL = 0.000001;

    if (A < TOL)
    {
      relp = 0.0;
    }
    else
    {
      if (A > 1.0) A = 1.0;
      relp = signDer*((EC*pow(A,EC-1.0))/(1.0 - SC));
    }
  }

  return relp;
}



double SFRelPerm::fractionalFlowWater(const double& s, const int& rockType, const SFPMPhysicalFieldProperties& IMPESPFP) const
{
  //***********************************************************************************************
  //This function computes the water-fractional flow for a given saturation s.
  //***********************************************************************************************


  //Note: s is the water saturation.

  double mu_o = IMPESPFP.getViscosity(0);
  double mu_w = IMPESPFP.getViscosity(1);

  double lambda_w = getRelPermData(1, s, rockType)/mu_w;
  double lambda_tot = getRelPermData(0, s, rockType)/mu_o;
  lambda_tot += lambda_w;

  return lambda_w/lambda_tot;
}


double SFRelPerm::fractionalFlowWaterDerivative(const double& s, const int& rockType, const SFPMPhysicalFieldProperties& IMPESPFP) const  
{
  //***********************************************************************************************
  //This function computes the derivative of the water-fractional flow wrt. the saturation s.
  //***********************************************************************************************


  //Note: s is the water saturation.

  double mu_o = IMPESPFP.getViscosity(0);
  double mu_w = IMPESPFP.getViscosity(1);

  double lambda_w = getRelPermData(1, s, rockType)/mu_w;
  double lambda_w_deriv = getRelPermDataPartial(1, s, rockType)/mu_w;
  double lambda_o = getRelPermData(0, s, rockType)/mu_o;
  double lambda_o_deriv = getRelPermDataPartial(0, s, rockType)/mu_o;
  double lambda_tot = lambda_o + lambda_w;

  return (lambda_w_deriv*lambda_tot - lambda_w*(lambda_w_deriv + lambda_o_deriv))/(lambda_tot*lambda_tot);
}



SFRelPerm::~SFRelPerm() 
{

}
