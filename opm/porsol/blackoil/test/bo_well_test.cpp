/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.

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

#include <opm/core/utility/have_boost_redef.hpp>

#include <iostream>
#include <opm/core/utility/Units.hpp>
#include <opm/porsol/blackoil/BlackoilWells.hpp>
#include <opm/porsol/common/Rock.hpp>
#include <iterator>

// Forward declaration
void write_fields(std::ostream& os, const Opm::EclipseGridParser& parser);

// Program for test of reading newly implemented keywords.
// Writes field data with and without unit conversion.

int main(int argc, char** argv)
{
    // Parser. 
    const std::string ecl_file = (argc == 1) ? "SPE9.DATA" : argv[1];
    std::cout << "Parsing file " << ecl_file << std::endl;

    Opm::EclipseGridParser parser1(ecl_file, false);
    
    // Stored fieldnames
    std::cout << "\nStored fieldnames.\n";
    std::vector<std::string> names = parser1.fieldNames();
    copy(names.begin(), names.end(),
	 std::ostream_iterator<std::string>(std::cout, "\n"));
    std::cout << std::endl;
    
    // Write data values for some of the keywords
    std::cout << "\nUnits are not converted to SI-units.\n";
    write_fields(std::cout, parser1);
    
    Opm::EclipseGridParser parser2(ecl_file);
    std::cout << "\nUnits are converted to SI-units.\n";
    write_fields(std::cout, parser2);
    
    return 0;
}

void write_fields(std::ostream& os, const Opm::EclipseGridParser& parser)
{
    if (parser.hasField("WELSPECS")) {
	parser.getWELSPECS().write(os);    
    }
    if (parser.hasField("COMPDAT")) {
	parser.getCOMPDAT().write(os);    
    }
    if (parser.hasField("WCONINJE")) {
	parser.getWCONINJE().write(os);    
    }
    if (parser.hasField("WCONPROD")) {
	parser.getWCONPROD().write(os);
    }
    if (parser.hasField("WELTARG")) {
	parser.getWELTARG().write(os);    
    }
    if (parser.hasField("EQUIL")) {
	parser.getEQUIL().write(os);    
    }
    if (parser.hasField("DENSITY")) {
	parser.getDENSITY().write(os);    
    }
    if (parser.hasField("PRESSURE")) {
	const std::vector<double>& pressure = 
	    parser.getFloatingPointValue("PRESSURE");
	os << "\nPRESSURE" << std::endl;
	copy(pressure.begin(), pressure.end(),
	     std::ostream_iterator<double>(os, " "));
	os << std::endl;
    }
    if (parser.hasField("SGAS")) {
	const std::vector<double>& sgas = 
	    parser.getFloatingPointValue("SGAS");
	os << "\nSGAS" << std::endl;
	copy(sgas.begin(), sgas.end(),
	     std::ostream_iterator<double>(os, " "));
	os << std::endl;
    }
    if (parser.hasField("SWAT")) {
	const std::vector<double>& swat = 
	    parser.getFloatingPointValue("SWAT");
	os << "\nSWAT" << std::endl;
	copy(swat.begin(), swat.end(),
	     std::ostream_iterator<double>(os, " "));
	os << std::endl;
    }
    if (parser.hasField("PVCDO")) {
	parser.getPVCDO().write(os);    
    }
}
