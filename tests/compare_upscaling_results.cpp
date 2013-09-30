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

#include "config.h"

#include <string>
#include <iostream>
#include <fstream>

#include <opm/porsol/common/Matrix.hpp>

using namespace std;
using namespace Opm;

typedef FullMatrix<double, OwnData, COrdering> Matrix;


/// @brief
///    Read result file and store results in a matrix.
///
/// @param [in] filename
///    Path to result file to be read. Lines starting with '#' are ignored.
///
/// @param [in] rows
///    Number of rows in result file
///
/// @param [in] cols
///    Number of cols in result file
Matrix readResultFile(const char* filename, int rows, int cols) {

    Matrix result(rows, cols, (const double*)0);

    // Check if filename is readable
    ifstream resultfile(filename, ios::in);
    if (resultfile.fail()) {
        cerr << "Error: Filename " << filename << " not found or not readable." << endl;
        exit(1);
    }

    // Read file and put results into matrix result
    double value;
    int row = 0;
    int col = 0;
    string line;
    while (getline(resultfile, line)) {
        if (line.empty() || line[0] == '#') // Skip empty lines and lines starting with #
            continue;
        istringstream iss(line);
        while (iss >> value) {
            // Check if input file contains more rows or columns than specified
            if (row >= rows) {
                cerr << "Error: More rows in " << filename << " than specified." << endl;
                exit(1);
            }
            else if (col >= cols) {
                cerr << "Error: More columns in " << filename << " than specified." << endl;
                exit(1);
            }
            else {
                result(row, col) = value;
                ++col;
            }
        }
        col = 0;
        ++row;
    }
    resultfile.close();

    return result;
}


/// @brief
///    Test if two matrices are equal within a relative tolerance
///
/// @param [in] refSoln
///    Reference solution stored as a matrix
///
/// @param [in] newSoln
///    Solution to compare with
///
/// @param [in] tol
///     Absolute tolerance
bool matrixAlmostEqual(Matrix refSoln, Matrix newSoln, double tol) {

    assert(refSoln.numRows() == newSoln.numRows());
    assert(refSoln.numCols() == newSoln.numCols());
    // Test element by element
    for (int row=0; row<refSoln.numRows(); ++row) {
        for (int col=0; col<refSoln.numCols(); ++col) {
            double absDiff = abs(refSoln(row,col) - newSoln(row,col));
            if (absDiff > tol) return false;
        }
    }
    return true;
}


/// @brief
///    Tests if two result files are equal
///
/// Command input variables;
///    1) Path to reference solution file
///    2) Path to new solution file to compare with
///    3) Absolute tolerance
///    4) Number of rows in result files to be compared
///    5) Number of columns in result files to be compared
///
/// The results in the input files are read and compared to eachother within a given absolute tolerance.
/// Lines starting with '#' are ignored. Returns 1 if test fails, 0 otherwise. If test failes, both
/// solutions are printed to screen.
int main(int varnum, char** vararg) try {

    // Check if the correct number of variables are given
    if (varnum != 6) {
        cout << "Error: Wrong number of input variables, should be five!" << endl
             << "Usage: ./test_upscaling_results refSolnFile newSolnFile tol rows cols" << endl;
        exit(1);
    }

    // Process input
    const char* refSolnFile(vararg[1]);
    const char* newSolnFile(vararg[2]);
    double tol = atof(vararg[3]);
    int rows = atoi(vararg[4]);
    int cols = atoi(vararg[5]);

    // Read result files into Matrix objects
    Matrix refSoln = readResultFile(refSolnFile, rows, cols);
    Matrix newSoln = readResultFile(newSolnFile, rows, cols);

    // Compare results
    bool test = matrixAlmostEqual(refSoln, newSoln, tol);
    if (! test) {
        cout << endl << "Verification error: Calculated solution not equal to reference solution "
             << "within an absolute  tolerance of " << tol << "." << endl
             << "Calculated solution:" << endl
             << newSoln
             << "Reference solution:" << endl
             << refSoln << endl;
        exit(1);
    }

    return 0;
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}

