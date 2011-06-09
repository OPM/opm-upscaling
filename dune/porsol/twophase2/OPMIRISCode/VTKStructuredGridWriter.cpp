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


//******************************************************************************
  // General note: all VTK point objects live in 3D, i.e. the geometry is always
  // 3D so that all point information found in VTK will always be a 3-tuple
  // (p_x, p_y, p_z) describing a point in 3D.
  // However, the space dimension of the more advanced topological objects
  // (vertex, line, polygon, voxel etc.) can vary between 1 and 3.
  // This is what we try to catch in getVTKCellType().


// @HAF: TODO: implement exceptions instead of cerr

#include "VTKStructuredGridWriter.h"

using namespace std;




//==============================================================================
VTKStructuredGridWriter::VTKStructuredGridWriter(const IRISDuneGridInterface<DuneGridType>* grid,
						     const Array<double>& vals,
						     const PMTwoPhaseSimulationData* simulationData,
						     SurfaceOutput& surface_output,
						     ValuesLocation values_location,
						     string header_comment
						 //char *header_comment
						 //, const VTKHeaderWriter::OutputFormat format
						 ) 
{
//==============================================================================
  init();
  header_.setHeaderComment(header_comment);
  //header_.setOutputFormat(format);
  header_.setOutputFormat(VTKHeaderWriter::ASCII);//@HAF: Change 27.01.2010.
  grid_ = grid;
  vals_ = vals;
  simulationData_ = simulationData;
  surface_output_ = surface_output;
  values_location_ = values_location;
}


//==============================================================================
void VTKStructuredGridWriter::init() {
//==============================================================================
  dataset_type_ = "DATASET STRUCTURED_POINTS";
  input_verified_ok_  = false;
}


//==============================================================================
bool VTKStructuredGridWriter::inputVerifiedOK() {
  //==============================================================================

  // grid space dimension
//  const int dim = grid.u_cellDimension_0u();

//************************************************************
//Note: This code will only work for 3D ??????????
//************************************************************
  const int dim = grid_->domainDimension();

  if (dim<1 || dim>3) {
    // VTK handles only this, can only visualise in 1D, 2D and 3D.
    cerr <<
      "VTKStructuredgridWriter.cpp: can only handle 2D and 3D grids."
	 << endl;
    return false;
  }
  
  
  if ((dim!=2 && (surface_output_==planar2D || surface_output_==non_planar2D)) ||
      (dim!=3 && surface_output_==standard3D)
      ) {
    cerr <<
      "VTKStructuredgridWriter.cpp: mismatch between space dimension of grid "
	 << "and grid_type"
	 << endl;
    return false;
  }
  
  input_verified_ok_ = true;
  return input_verified_ok_;
}




//==============================================================================
VTKStructuredGridWriter::VTKCellType
VTKStructuredGridWriter::getVTKCellType(const int& no_nodes,
					  const int& dim) const {
//==============================================================================
  // Based on number of nodes in cell and the space dimension of the input grid,
  // (typically 2D or 3D) return some cell type
  
  // See top of file for information

  // Note that all of these cases have not been tested properly.
  // Furthermore note that this function is written specifically for
  // RFSophus StructuredGrid objects.

  if (no_nodes < 1) {
    cerr <<
      "VTKStructuredgridWriter.getVTKCellType: number of nodes < 1."
	 << endl;
    return VTK_VERTEX;   // @ers: implement exception instead
  }

  switch (dim) {
  case 1:
    switch (no_nodes) {
    case 1:
      return VTK_VERTEX;
      break;
    // all cases after this might be a VTK_POLY_VERTEX...
    case 2:
      return VTK_LINE;
      break;
    default:
      return VTK_POLY_LINE;
      break;
    }
    break;
    
  case 2:
    switch (no_nodes) {
    case 1:
      return VTK_VERTEX;
      break;
    // all cases after this might be a VTK_POLY_VERTEX...
    case 2:
      return VTK_LINE;
      break;
    case 3:
      return VTK_TRIANGLE;
      break;
    case 4:
      return VTK_QUAD;
      break;
    default:
      return VTK_POLYGON;
      break;
    }
    break;

  case 3:
    switch (no_nodes) {
    case 1:
      return VTK_VERTEX;
      break;
    // all cases after this might be a VTK_POLY_VERTEX...
    case 2:
      return VTK_LINE;
      break;
    case 3:
      return VTK_TRIANGLE;
      break;
    case 4:
      return VTK_TETRA;
      break;
    case 5:
      return VTK_PYRAMID;
      break;
    case 6:
      return VTK_WEDGE;
      break;
    case 8:
      return VTK_HEXAHEDRON;
      break;
    default:
      cerr << "Undefined cell type!  Do something about it!" << endl;
      return VTK_VERTEX;   // @ers: implement exception instead
      break;
    }
    break;
    
  default:
    cerr << "Not possible with other dimensions than 1,2,3 in VTK." << endl;
    return VTK_VERTEX;   // @ers: implement exception instead
    break;
  }
}



//==============================================================================
void VTKStructuredGridWriter::write(ostream& os) const {
//==============================================================================
  if (!input_verified_ok_) {
    cerr << "VTKStructuredGridWriter.cpp: input not verified, "
	 << "cannot write to file."
	 << endl;
    return;
  }

  header_.writeHeader(os);
  writeDataSetType(os);
  writeDataSetStructure(os);
  writeDataSetAttributes(os);
}





//==============================================================================
void VTKStructuredGridWriter::writeDataSetType(ostream& os) const {
//==============================================================================
  os << dataset_type_ << endl;
}



//==============================================================================
void VTKStructuredGridWriter::writeDataSetStructure(ostream& os) const {
//==============================================================================

//************************************************************
//@HAF: NOTE: This code will only work for 3D ??????????
//************************************************************

//************************************************************
//@HAF: NOTE2: This routine assumes that the orogin always is (0.0,0.0,0.0)
//************************************************************

//  os << "DIMENSIONS " << (*simulationData_.td_.NoGridCells_X)+1 << " " << (*simulationData_.td_.NoGridCells_Y)+1 << " " << (*simulationData_.td_.NoGridCells_Z)+1 << endl;
  os << "DIMENSIONS " << (simulationData_->td_.NoGridCells_X)+1 << " " << (simulationData_->td_.NoGridCells_Y)+1 << " " << (simulationData_->td_.NoGridCells_Z)+1 << endl;

  //os << "ASPECT_RATIO " << *simulationData_.td_.X_Extent/(*simulationData_.td_.NoGridCells_X) << " " << *simulationData_.td_.Y_Extent/(*simulationData_.td_.NoGridCells_Y) << " " << *simulationData_.td_.Z_Extent/(*simulationData_.td_.NoGridCells_Z) << endl;
  os << "ASPECT_RATIO " << simulationData_->td_.X_Extent/(simulationData_->td_.NoGridCells_X) << " " << simulationData_->td_.Y_Extent/(simulationData_->td_.NoGridCells_Y) << " " << simulationData_->td_.Z_Extent/(simulationData_->td_.NoGridCells_Z) << endl;

  //os << "ORIGIN " << grid_->GetMinX() << " " << grid_->GetMinY() << " " << grid_->GetMinZ() << endl;
  os << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
  
  os << endl;

}


//==============================================================================
void VTKStructuredGridWriter::writeDataSetAttributes(ostream& os) const {
//==============================================================================

  switch (values_location_) {
  case on_nodes:
    os << "POINT_DATA ";
    break;
  case on_cells:
    os << "CELL_DATA ";
    break;
  default:
    // @ers: TODO throw exception
    break;
  }

  const int no_values = vals_.u_getSize_0u();

  os << no_values << endl;

  os << "SCALARS scalars double" << endl;//@HAF: This is best for paraview
  //os << "SCALARS scalars float" << endl;
  os << "LOOKUP_TABLE default" << endl;

  //  const int no_elements_in_row = 5;   // format output somewhat

  // assumes sorting of vals_ is identical to DuneGridType
  for (int i=0; i<no_values; ++i) {
    os << vals_[i];
     os << endl;
  }
}
