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


#ifndef GlobConstH
#define GlobConstH




// Defines maximum spatial dimensions for the problem
//----------------------------------------------------

#define MAXDIRECTIONSIZE  4



//-------------NEW FROM HAF---START---------------------------------
//uncertain about this size; is it big enough ???
#define MAXLOCARRAYSIZE  3
//-------------NEW FROM HAF---END---------------------------------





#define MAXDIMENSIONSIZE 3
#define ARRAYMAXSIZE  8
  



// number of tensor elements being needed
#define MAXTENSORSIZE 6
#define MAXCONTSIZE MAXTENSORSIZE



#define CARRIER_SIZE 4096



#define SGML_MAXTAGLENGTH 200







#endif //GlobConstH



/*
===============================================================================
    ------------------------ CLASS DESCRIPTION --------------------------
===============================================================================

  DESCRIPTION:
  ------------------
    The file contains global constants necessary to adjust Sophus to 
    ordinary C++. Sophus is originally implemented with a set of global
    constants used to define object size etc. These constants will
    be found here.

    The intention is to include this file as one of the first include files
    and as few places as possible.
    

  USAGE/EXAMPLES:
  ---------------
    <Examples how to use the class>
    
    
  MISCELLANEOUS:
  --------------
  
  
  


===============================================================================
    --------------------- END OF CLASS DESCRIPTION ----------------------
===============================================================================
*/
