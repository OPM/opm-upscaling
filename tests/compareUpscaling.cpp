/*
  Copyright 2016 Statoil ASA.
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <opm/common/ErrorMacros.hpp>



void readData(const char* fileName, std::vector<double> &dataVec){
    std::ifstream inFile (fileName, std::ifstream::in);
    std::string line, value;
    if(inFile.fail()){
        std::string fname = fileName;
        std::cout << "Could not open file " + fname << std::endl;
        OPM_THROW(std::runtime_error, "Could not open file " + fname);
    }
    while(std::getline(inFile, line)){
	std::stringstream ss(line);
	while(ss >> value){
	    if(value.substr(0,1) == "#" ){break;}
	    dataVec.push_back(atof(value.c_str()));
	}
    }

}

void startTest(const char* testFile, const char* compareFile,  double absoluteTolerance, double relativeTolerance){
    std::vector<double> testVec, compareVec;
    readData(testFile, testVec);
    readData(compareFile, compareVec);

    if(testVec.size()!= compareVec.size()){
	OPM_THROW(std::invalid_argument, "The upscaling files do not have the same number of values.");
    }

    double absDev, relDev;
    for (size_t i = 0; i < testVec.size(); i++){

	absDev = std::abs(testVec[i]-compareVec[i]);
	relDev = absDev/double(std::max(std::abs(testVec[i]),std::abs(compareVec[i])));
	
	if (relDev > relativeTolerance && absDev > absoluteTolerance){
	    for(auto val : testVec) {
		std::cout << "For data value " << val << "- The relative deviation is "<< relDev << " greater than the limit "<< relativeTolerance << std::endl;
		OPM_THROW(std::invalid_argument, "The relative deviation exceeded the limit.");
	    }
	}

    }
}



void printHelp(){
    std::cout << "The function takes four arguments: \n 1) Path to file 1 \n 2) Path to file 2 \n 3) The maximum absolute deviation that is acceptable. \n 4) The maximum relative deviation that is acceptable." << std::endl;
}

int main(int argc, const char ** argv){
    if(argc!=5){
	printHelp();
	return EXIT_FAILURE;
    }
    const char* testFileName = argv[1];
    const char* checkFileName = argv[2];
    double absoluteTolerance = atof(argv[3]);
    double relativeTolerance = atof(argv[4]);
    try{
	startTest(testFileName, checkFileName, absoluteTolerance, relativeTolerance);
    }
    catch(std::exception& e){
	std::cout << e.what() << std::endl;
	return EXIT_FAILURE;
    }

    return 0;
}
