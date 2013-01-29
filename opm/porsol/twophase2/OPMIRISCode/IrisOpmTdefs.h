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


#ifndef IrisOpmTdefsH
#define IrisOpmTdefsH

/* Some dune stuff */
#include<dune/grid/sgrid.hh> // load sgrid definition
#include<dune/grid/alugrid.hh> // load alugrid definition
//#include<dune/grid/io/file/dgfparser/gridtype.hh>
//#include<dune/grid/io/file/dgfparser/dgfs.hh>// load sgrid definition through the dgf-format.
//#include<dune/grid/io/file/dgfparser/dgfalu.hh>// load alu-grid definition through the dgf-format.

//#define ALUGRID_SIMPLEX 
//#define GRIDDIM 2

//@HAF: NOTE this code must obviously be changed for another dune-grid type
typedef Dune::SGrid<2,2> DuneGridType;//This is a 2D structured grid
//typedef Dune::ALUSimplexGrid<2,2> DuneGridType;//This is a 2D triangular grid

//typedef Dune::SGrid<2,2> DuneGridTypeTest;//This is a 2D structured grid
//typedef Dune::ALUSimplexGrid<2,2> DuneGridTypeTest;//This is a 2D triangular grid


/* Some IRIS-Opm stuff */
#include <opm/porsol/twophase2/OPMKvasiSophusCode/Array.h>

typedef Array<double> SFSaturation;
typedef Array<double> SFVelocity;
typedef Array<double> SFPressure;

//An enumeration:
enum TimeIntegrationType{_Euler, _SecondOrderEuler};//@HAF: Should it be placed elsewhere???



#endif //IrisOpmTdefsH



