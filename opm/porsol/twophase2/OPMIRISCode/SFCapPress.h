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


#ifndef SFCapPress_H
#define SFCapPress_H

/*
#ifdef DEBUG_SFCapPress_h

#define DEBUG_NOW

#endif
*/
//---------------------- Include Files ----------------------------------------

#include <opm/porsol/twophase2/OPMKvasiSophusCode/GlobType.h>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/Array.h>

#include <iostream>
#include <math.h>

//---------------------  Constants --------------------------------------------


//---------------------- Types ------------------------------------------------


class SFCapPress
{
  public : 
//Boolean DI_0u() const ;
  SFCapPress() ;

  void readCapPressDataInFormulaFormat(std::istream& inputStream, const int& numbRockTypes);

  double getCapPressData(const int& capPressType, const double& sw, const int& rockType) const ;
  double getCapPressDataPartial(const int& capPressType, const double& sw, const int& rockType) const ;
  double getCapPressDataPartial2(const int& capPressType, const double& sw, const int& rockType) const ;

  double getInverseCapPressData(const int& capPressType, const double& pcv, const int& rockType) const;

  double getOutletSaturation(const int& capPressType, const int& rockTypeEnd) const;

  ~SFCapPress() ;


  private : 
    double ml_constructInitialGuessForOutletSaturation(const int& capPressType, const int& rockTypeEnd) const;

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
    The first index in the imaginary "double Dune::array" repreresents the actual rock type.
    In every "row" we store the following "Leverett J-type constants" (NB: assuming the possibility
    of 3 phase flow):
    CL(0), CNOL(0), CR(0), EL(0), ER(0), SL(0), SR(0), CL(1), CNOL(1), CR(1), EL(1), ER(1), SL(1), SR(1)
    NB: Where the indices in EC etc. just refers to:
    capPressType=0 : oil-water capillary pressure    
    capPressType=1 : oil-gas capillary pressure    
    BUT of course e.g. SR(0) occupies index 6 in the actual row in m_constants etc.
  */
//DECLARE_SOPHUS(1.4)
};


#endif //SFCapPress_H



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

