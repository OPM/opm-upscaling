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


//==============================================================================
// VTKHeaderWriter
//==============================================================================

#include "VTKHeaderWriter.h"

using namespace std;


//==============================================================================
//VTKHeaderWriter::VTKHeaderWriter(const string& header_comment,
//VTKHeaderWriter::VTKHeaderWriter(const char *header_comment,
//				 const OutputFormat& format) {
VTKHeaderWriter::VTKHeaderWriter(string header_comment,
				 const OutputFormat& format) {
//==============================================================================
  setHeaderComment(header_comment);
  setOutputFormat(format);
  init();
}


 //==============================================================================
void VTKHeaderWriter::init() {
//==============================================================================
  id_and_version_ = "# vtk DataFile Version 3.0";
}


//==============================================================================
void VTKHeaderWriter::setOutputFormat(const OutputFormat& format) {
//==============================================================================
  format_ = format;
  switch (format_) {
  case ASCII:
    // OK
    break;
  case binary:
    cerr << "VTKHeaderWriter.cpp: Does not support binary output." << endl
	 << "Currently only the ASCII format is supported." << endl
	 << "Format set to ASCII."
	 << endl;
    format_ = ASCII;
    break;
  case XML:
    cerr << "VTKHeaderWriter.cpp: Does not support XML output." << endl
	 << "Currently only the ASCII format is supported." << endl
	 << "Format set to ASCII."
	 << endl;
    format_ = ASCII;
    break;
  default:
    cerr << "VTKHeaderWriter.cpp: Unsupported format." << endl
	 << "Currently only the ASCII format is supported." << endl
	 << "Format set to ASCII."
	 << endl;
    format_ = ASCII;
    break;
  }

}


//==============================================================================
void VTKHeaderWriter::writeHeader(ostream& os) const {
//==============================================================================
  writeFileVersionAndIdentifyer(os);
  writeHeaderComment(os);
  writeOutputFormat(os);
}


//==============================================================================
void VTKHeaderWriter::writeFileVersionAndIdentifyer(ostream& os) const {
//==============================================================================
  os << id_and_version_ << endl;
}


//==============================================================================
void VTKHeaderWriter::writeHeaderComment(ostream& os) const {
//==============================================================================
  os << header_comment_ << endl;
}


//==============================================================================
void VTKHeaderWriter::writeOutputFormat(ostream& os) const {
//==============================================================================
  
  switch (format_) {
  case ASCII:
    os << "ASCII";
    break;
    // insert binary/xml handling here
  default:
    // @ers: TODO do not abort!!
    // use exceptions
    cerr << "VTKHeaderWriter.cpp: this should never happen! Aborting." << endl;
    exit(1);
    break;
  }

  os << endl;
}
