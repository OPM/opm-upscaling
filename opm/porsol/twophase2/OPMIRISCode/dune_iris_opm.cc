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


#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif
#include <iostream>
#include"dune/common/mpihelper.hh" // An initializer of MPI
#include"opm/porsol/blackoil/co2fluid/opm/common/exceptions.hh" // We use exceptions

//#include<dune/grid/sgrid.hh> // load sgrid definition
//#include<dune/grid/common/gridinfo.hh> // definition of gridinfo

//#include<dune/grid/common/mcmgmapper.hh>
//#include<dune/grid/common/entity.hh>
//#include<dune/grid/common/entitypointer.hh>

//#include<dune/common/fmatrix.hh>
//#include<dune/istl/bcrsmatrix.hh>
//#include <istl/solvers.hh>
//#include <istl/preconditioners.hh>

/* IRISOpm Classes to be compiled */
#include <IRISDuneGridInterface.hpp>
#include <IRISDuneISTLInterface.hpp>

#include <Array.h>
#include <MeshShape.h>
#include <MeshPoint.h>
#include <Pair.h>
#include <RnShape.h>
#include <RnPoint.h>
#include <MeshSFSD.h>

#include <PMTwoPhaseSimulationData.h>
#include <VTKHeaderWriter.h>
#include <VTKStructuredGridWriter.h>
#include <SFPMPhysicalFieldProperties.h>
#include <SFCapPress.h>
#include <SFRelPerm.h>

#include <SFCentralUpwind.h>
#include <SFMpfaFps.h>

#include <PMTwoPhaseSaturation.h>
#include <PMTwoPhaseElliptic.h>
#include <PMTwoPhaseSimulator.h>

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    std::cout << "Hello World! This is dune-iris-opm." << std::endl;
    if(Dune::MPIHelper::isFake)
    {
      std::cout<< "This is a sequential program." << std::endl;
      /*-----------------------------------------------------------------------------*/
      //                            The main program.


      //----------------------------------------------------------------
      //A) Must read necessary input stuff from files...
      //----------------------------------------------------------------
      PMTwoPhaseSimulationData simulationData;
      simulationData.readPMTwoPhaseSimulationData();
     
      //----------------------------------------------------------------
      //B) Must do necessary stuff connected to the dune-grid:
      //----------------------------------------------------------------
      //@HAF: Now tries to structure this code such that it can also work for other types of dune-grids
      //NOTE: DuneGridType is defined in IrisOpmTdefs.h
      //NOTE2: This does NOT work yet...only Sgrid can be used yet!

      //#ifdef Dune::SGridToBeUsed
      cout << endl;
      cout << simulationData.td_.NameOfGridType << endl;
      cout << simulationData.td_.NoGridCells_X << endl;
      cout << simulationData.td_.NoGridCells_Y << endl;
      DuneGridType* DGridPtr;
      //if (simulationData.td_.NameOfGridType == "Dune::SGrid2D")
      //{
	const int dim=2;                               
	Dune::FieldVector<int,dim> N(2);
	//Puts the correct number of grid-blocks in each direction:              
	N[0] = simulationData.td_.NoGridCells_X; //@HAF: OK???
	N[1] = simulationData.td_.NoGridCells_Y; //@HAF: OK???
	//The geometrical size of the domain:
	Dune::FieldVector<DuneGridType::ctype,dim> L(0.0);
	Dune::FieldVector<DuneGridType::ctype,dim> H(1.0); 
	H[0] = simulationData.td_.X_Extent; //@HAF: OK???
	H[1] = simulationData.td_.Y_Extent; //@HAF: OK???
	//Creates the Dune::SGrid:
	DuneGridType Dgrid(N,L,H);
	//DGridPtr = &Dgrid;

	//}
      //#endif //SGridToBeUsed
#ifdef AluGridToBeUsed
      if (simulationData.td_.NameOfGridType == "ALUGridSimplex2D")
      {
	Dune::GridFactory<DuneGridType> factory();
	//Must fill in the vertices and elements from Triangle:
	//insertVertices( factory, "A.1.node");	
	//insertSimplices( factory, "A.1.ele");
	//Vertices:
	Dune::FieldVector<double,2> pos;
	double x,y,dummy;
	int AnzKnoten;

	std::ifstream ifs("A.1.node");	
	
	if (ifs.is_open()) {
	  ifs >> AnzKnoten >> dummy >> dummy >> dummy;
	  for (int i = 0; i < AnzKnoten; i++) {
	    ifs >> dummy >> x >> y >> dummy >> dummy;
	    pos[0] = x; pos[1] = y;
	    factory.insertVertex(pos);
	  } 
	  ifs.close();
	}
	else
	  std::cout << "Error: Can not open file\n";

	//Elements:
	unsigned int numberOfSimplices, corner1, corner2, corner3, dummy1;
	const Dune::GeometryType type( Dune::GeometryType::simplex, 2);
	std::vector< unsigned int > cornerIDs( 3 );

	std::ifstream ifs1("A.1.ele");

	if (ifs1.is_open()) {
	  ifs1 >> numberOfSimplices >> dummy1 >> dummy1;    
	  for (unsigned int i = 0; i < numberOfSimplices; i++) {
	    ifs1 >> dummy1 >> corner1 >> corner2 >> corner3;
	    cornerIDs[0] = corner1 - 1;
	    cornerIDs[1] = corner2 - 1;
	    cornerIDs[2] = corner3 - 1;
	    factory.insertElement(type, cornerIDs);      
	  }  
	  ifs1.close();
	}
	else
	  std::cout << "Error: Can not open file\n";
	
	//-------------------------------------------------------------------------
	DGridPtr = factory.createGrid("alugrid", true);
      }
#endif //AluGridToBeUsed

      //-----------------------------------------------------------------------------------
      //C) Set up of the IRISDuneGridInterface (i.e. the IRIS interface to the Dune grid).
      //-----------------------------------------------------------------------------------
      //DuneGridType& DGrid = *DGridPtr;
      IRISDuneGridInterface<DuneGridType> IGI(Dgrid);
      IGI.setupForMPFAFPS_FVM();


	//NOTE: The dgf file stuff is NOT working properly for the time being...
        //and maybe also does NOT work together with the above non-dgf stuff.....
#ifdef dfgFileToBeUsed
      Dune::GridPtr<DuneGridType> DGridPtr(simulationData.td_.NameOfDGFFile);
      DuneGridType& DGrid = *DGridPtr; 
      IRISDuneGridInterface<DuneGridType> IGI(DGrid);
      exit(0);//Just for testing/debugging purposes
      IGI.setupForMPFAFPS_FVM();
#endif //dfgFileToBeUsed
      //----------------------------------------------------------------
      //D) Must handle and set up the physical field properties
      //----------------------------------------------------------------
      SFPMPhysicalFieldProperties PMPhysProp;
      PMPhysProp.setup(&IGI, &simulationData);
      
      //----------------------------------------------------------------
      //E) Set up the two-phase simulator
      //----------------------------------------------------------------

      //*******CHOICE OF NUMERICAL METHODS: START**************************
      typedef SFMpfaFps FVMImplicit;
      typedef SFCentralUpwind FVMExplicit;
      //*******CHOICE OF NUMERICAL METHODS: END****************************
   
      PMTwoPhaseSimulator<FVMImplicit,FVMExplicit> simulator(IGI.cellCount(0,0));
      simulator.setup(IGI, simulationData, PMPhysProp);
      
      //----------------------------------------------------------------
      //F) Perform the simulation
      //----------------------------------------------------------------
      TimeIntegrationType timeIntegrationType = _SecondOrderEuler;
      simulator.performSimulation(IGI, simulationData, timeIntegrationType);

      std::cout<< "----------------------------------------------------------" << std::endl;
      std::cout<< "Simulation completed........" << std::endl;
      std::cout<< "----------------------------------------------------------" << std::endl;
      /*-----------------------------------------------------------------------------*/
    }
    else
    {
      std::cout<<"I am rank "<<helper.rank()<<" of "<<helper.size()
        <<" processes!"<<std::endl;
    }
    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
