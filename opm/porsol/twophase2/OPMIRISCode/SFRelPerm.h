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



#ifndef SFRelPerm_H
#define SFRelPerm_H

//---------------------- Include Files ----------------------------------------

#include <opm/porsol/twophase2/OPMKvasiSophusCode/GlobType.h>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/Array.h>

#include <iostream>

#include "SFPMPhysicalFieldProperties.h"

//---------------------- Directives -------------------------------------------
//using namespace std;  // This leads to that all names in the namespace std may be 
                      // used without the qualifier "std::".


//---------------------  Constants --------------------------------------------


//---------------------- Types ------------------------------------------------


class SFRelPerm
{
  public : 
  SFRelPerm() ;

  void readRelPermDataInFormulaFormat(std::istream& inputStream, const int& numbTypes);

 double getRelPermData(const int& phase, const double& sw, const int& rockType) const ;
 double getRelPermDataPartial(const int& phase, const double& sw, const int& rockType) const ;
 double fractionalFlowWater(const double& s, const int& rockType, const SFPMPhysicalFieldProperties& IMPESPFP) const;
 double fractionalFlowWaterDerivative(const double& s, const int& rockType, const SFPMPhysicalFieldProperties& IMPESPFP) const;

  ~SFRelPerm() ;


  private : 
    int m_variant;
    int m_RowLength;
  /*
    m_variant = 0 : rel. perm. in formula
    m_variant = 1 : rel. perm. in table
  */

  double *m_constants;
  /*
    These constants are to be used in conjunction with the formula representation.
    We use a 1D Dune::array repr., where the outer imaginary row has a length of m_RowLength !!
    In every "row" we store the following "Corey constants" (NB: assuming the possibility
    of 3 phase flow):
    EC(0), EC(1), EC(2), SC(0), SC(1), SC(2)
    NB: Where the indices in EC etc. just refers to the phases i.e.
    phase=0 : oil    
    phase=1 : water
    phase=3 : gas,
    BUT of course e.g. SC(0) occupies index 3 in the actual row in m_constants etc.
  */
//DECLARE_SOPHUS(1.4)
};


#endif //SFRelPerm_H



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

