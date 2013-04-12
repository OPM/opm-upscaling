/*
  Copyright 2013 Statoil ASA

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

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>

#include <opm/upscaling/SinglePhaseUpscaler.hpp>

using namespace std;

typedef Opm::SinglePhaseUpscaler::permtensor_t Matrix;


// Read result file and store result as vector.
Matrix readResultFile(string filename) {

  Matrix perm(3, 3, (const double*)0);

  // Check if filename is readable
  char* file = (char*)filename.c_str();
  ifstream resultfile(file, ios::in);
  if (resultfile.fail()) {
    cerr << "Error: Filename " << filename << " not found or not readable." << endl;
    exit(1);
  }

  // Read file and put results into matrix perm
  double value;
  int row = 0;
  int col = 0;
  string line;
  while (getline(resultfile, line)) {
    if (line[0] == '#') // Skip lines starting with #
      continue;
    istringstream iss(line);
    while (iss >> value) { 
      perm(row, col) = value;
      ++col;
    }
    col = 0;
    ++row;
  }
  
  resultfile.close();

  return perm;
}


int main() {

  // Specify accepted relative tolerance
  double relTol = 1e-6;

  // Specify relative paths to files
  string gridPath = "../tests/input_data/grids/PeriodicTilted.grdecl";
  string referenceSolutionPath = "../tests/input_data/reference_solutions/upscale_perm_fixed_PeriodicTilted.txt";
  string outputPath = "../tests/temp_results.txt";

  // Run executable upscale_perm with input data. 
  string runCommand = string("./upscale_perm -output ") + outputPath + string(" ") + gridPath;
  system(runCommand.c_str());

  // Store solutions in Matrix object
  Matrix new_solution = readResultFile(outputPath);
  Matrix reference_solution = readResultFile(referenceSolutionPath);
 
  // Remove temporary result file
  string removeCommand = string("rm ") + outputPath;
  system(removeCommand.c_str());
 
  // Compare result with reference solution
  for (int r=0; r<3; ++r) {
    for (int c=0; c<3; ++c) {
      double absRelDiff = abs(reference_solution(r,c) - new_solution(r,c))/reference_solution(r,c);
      if (absRelDiff > relTol) {
	cerr << "Error: Calculated solution not equal to reference solution within an relative tolerance of "
	     << relTol << ". "  << endl
	     << "Calculated permeability tensor:" << endl
	     << new_solution << endl
	     << "Reference solution:" << endl
	     << reference_solution << endl;
	exit(1);
      }
    }
  } 

  return 0;

}
