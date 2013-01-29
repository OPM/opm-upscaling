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


#ifndef VTKStructuredGridWriterH
#define VTKStructuredGridWriterH

#include <opm/porsol/twophase2/OPMIRISCode/IRISDuneGridInterface.hpp>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/GlobType.h>
#include <opm/porsol/twophase2/OPMIRISCode/IrisOpmTdefs.h>
#include <opm/porsol/twophase2/OPMIRISCode/VTKHeaderWriter.h>
#include <opm/porsol/twophase2/OPMKvasiSophusCode/Array.h>
#include <opm/porsol/twophase2/OPMIRISCode/PMTwoPhaseSimulationData.h>
#include <iostream>

using namespace std;


class VTKStructuredGridWriter {

  // member variables
 public:
  // first two are for z-values in the nodes, third is for z-values in 2D cells
  // planar2D: write 2D planar grid (z-values = 0)
  // non_planar2D: write non planar 2D grid (free z-values)
  // standard3D: 'normal' 3D grid
  typedef enum {planar2D, non_planar2D, standard3D} SurfaceOutput;

  // POINT_DATA or CELL_DATA
  typedef enum {on_nodes, on_cells} ValuesLocation ;

  // The list of cell types is found in the vtk library, in Common/vtkCellType.h
  typedef enum {
    // Linear cells
    VTK_VERTEX = 1,
    VTK_POLY_VERTEX = 2,
    VTK_LINE = 3,
    VTK_POLY_LINE = 4,
    VTK_TRIANGLE = 5,
    VTK_TRIANGLE_STRIP = 6,
    VTK_POLYGON = 7,
    VTK_PIXEL = 8,
    VTK_QUAD = 9,
    VTK_TETRA = 10,
    VTK_VOXEL = 11,
    VTK_HEXAHEDRON = 12,
    VTK_WEDGE = 13,
    VTK_PYRAMID = 14,
    
    // Quadratic, isoparametric cells
    VTK_QUADRATIC_EDGE = 21,
    VTK_QUADRATIC_TRIANGLE = 22,
    VTK_QUADRATIC_QUAD = 23,
    VTK_QUADRATIC_TETRA = 24,
    VTK_QUADRATIC_HEXAHEDRON = 25,
    VTK_QUADRATIC_WEDGE = 26,
    VTK_QUADRATIC_PYRAMID = 27,

    // Special class of cells formed by convex group of points
    VTK_CONVEX_POINT_SET = 41,

    // Higher order cells in parametric form
    VTK_PARAMETRIC_CURVE = 51,
    VTK_PARAMETRIC_SURFACE = 52,
    VTK_PARAMETRIC_TRI_SURFACE = 53,
    VTK_PARAMETRIC_QUAD_SURFACE = 54,
    VTK_PARAMETRIC_TETRA_REGION = 55,
    VTK_PARAMETRIC_HEX_REGION = 56
  } VTKCellType ;

 protected:
 private:
  VTKHeaderWriter                            header_;
  const IRISDuneGridInterface<DuneGridType>* grid_;
  Array<double>                        vals_;
  const PMTwoPhaseSimulationData* simulationData_;
  //char                                       *dataset_type_;
  string                                       dataset_type_;//@HAF: OK? (14.01.2010)
  SurfaceOutput                              surface_output_;
  ValuesLocation                             values_location_;
  bool                                       input_verified_ok_;


  // member functions
 public:
  VTKStructuredGridWriter(const IRISDuneGridInterface<DuneGridType>* grid,
			    const Array<double>& vals,
			    const PMTwoPhaseSimulationData* simulationData,
			    SurfaceOutput& surface_output,
			    ValuesLocation values_location,
			    string header_comment
			  //char *header_comment
			  //, const VTKHeaderWriter::OutputFormat format = 
			  //  VTKHeaderWriter::ASCII
			  );
  ~VTKStructuredGridWriter() {};

  // verify write ability of input data, mandatory to call
  bool inputVerifiedOK();
  void write(ostream& os) const ;


 protected:
  // Default constructor not allowed for user of class
  // Prevents setting up a class with illegal private variables
  // so that we possibly would have to check for legal member variables
  VTKStructuredGridWriter() { init(); };

 private:
  void writeDataSetType(ostream& os) const ;
  void writeDataSetStructure(ostream& os) const ;
  void writeDataSetAttributes(ostream& os) const ;
  void init();
  VTKCellType getVTKCellType(const int& no_nodes, const int& dim) const ;

  // prevent copying
  VTKStructuredGridWriter(const VTKStructuredGridWriter&);
  VTKStructuredGridWriter& operator=(const VTKStructuredGridWriter&);
};


#endif // VTKStructuredGridWriter
