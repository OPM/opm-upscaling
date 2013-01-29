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

#ifndef VTKHeaderWriterH
#define VTKHeaderWriterH

//==============================================================================
// Class VTKHeaderWriter

// Purpose:
// Write header data to VTK stream.
//==============================================================================


#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;


class VTKHeaderWriter {

  //--------- member variables -------------------------------------------------
 public:
  typedef enum {ASCII, binary, XML} OutputFormat ;

 protected:
 private:
  string                      header_comment_;
  //const char                  *header_comment_;
  // currently only the ASCII format is supported
  OutputFormat                format_;
  string                      id_and_version_;//@HAF: OK? (14.01.2010)
  //char                        *id_and_version_;
  

  //--------- member functions -------------------------------------------------
 public:

//  VTKHeaderWriter(const string& header_comment, const OutputFormat& format); 
  //VTKHeaderWriter(const char *header_comment, const OutputFormat& format); 
  VTKHeaderWriter(string header_comment, const OutputFormat& format); 
  VTKHeaderWriter() { init(); };
  ~VTKHeaderWriter() {};

  // set functions
//  void setHeaderComment(const string& header_comment)
  //void setHeaderComment(const char *header_comment)
  void setHeaderComment(string header_comment)
    {header_comment_ = header_comment;};
  void setOutputFormat(const OutputFormat& format);

  void writeHeader(ostream& os) const;


 protected:


 private:
  void init();
  void writeFileVersionAndIdentifyer(ostream& os) const ;
  void writeHeaderComment(ostream& os) const ;
  void writeOutputFormat(ostream& os) const ;
  
};


#endif // VTKHeaderWriterH
