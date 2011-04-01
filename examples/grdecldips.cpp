
/*
  Copyright 2011 Statoil ASA.

  This file is part of The Open Porous Media project (OPM).

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

#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <numeric>
#include <sys/utsname.h>

#include <dune/common/EclipseGridParser.hpp>
#include <dune/common/EclipseGridInspector.hpp>
#include <dune/common/param/ParameterGroup.hpp>

#include <dune/common/CornerpointChopper.hpp>
#include <dune/upscaling/SinglePhaseUpscaler.hpp>
#include <dune/porsol/common/setupBoundaryConditions.hpp>
#include <dune/common/Units.hpp>



using namespace std;


int main(int argc, char** argv) {        
    if (argc ==  1) { // If no arguments supplied 
        cout << "Usage: grdecldips gridfilename=foo.grdecl [mincellvolume=1e-8] [skiptopbottom=true]" << endl;
        cout << "       [listallcells=false] [output=filename.txt]" << endl;
        exit(1);
    } 
    Dune::parameter::ParameterGroup param(argc, argv);
    
    std::string gridfilename = param.get<std::string>("gridfilename");
    double minCellVolume = param.getDefault("mincellvolume", 1e-8);
    bool skipTopBottom = param.getDefault("skiptopbottom", true);
    bool listallcells = param.getDefault("listallcells", false);
    std::string outputfilename = param.getDefault<std::string>("output", "");
    
    //eclParser_p = new Dune::EclipseGridParser(gridfilename);
    //Dune::EclipseGridParser& eclParser = *eclParser_p;
    
    Dune::EclipseGridParser eclParser(gridfilename);
    
    Dune::EclipseGridInspector gridinspector(eclParser);
    
    // Check that we have the information we need from the eclipse file, we will check PERM-fields later
    if (! (eclParser.hasField("SPECGRID") && eclParser.hasField("COORD") && eclParser.hasField("ZCORN"))) {  
        cerr << "Error: Did not find SPECGRID, COORD and ZCORN in Eclipse file " << gridfilename << endl;  
        exit(1);  
    }
    
    /***************************
     * Find dips for every cell.
     */
    
    vector<int>  griddims = eclParser.getSPECGRID().dimensions;
    vector<double> xdips, ydips, cellvolumes;
    vector<int> cellidxs_i, cellidxs_j, cellidxs_k;
    
    int ignoredCellCount = 0; 
    for (int k=0; k < griddims[2]; ++k) {
	for (int j=0; j < griddims[1]; ++j) {
	    for (int i=0; i < griddims[0]; ++i) { 
		double cellVolume = gridinspector.cellVolumeVerticalPillars(i, j, k);
		if (cellVolume > minCellVolume) {
		    std::pair<double,double> xydip = gridinspector.cellDips(i, j, k);
		    xdips.push_back(xydip.first);
		    ydips.push_back(xydip.second);
		    cellvolumes.push_back(cellVolume);
		    cellidxs_i.push_back(i);
		    cellidxs_j.push_back(j);
		    cellidxs_k.push_back(k);
		}
		else {
		    ignoredCellCount++;
		}
	    }
	}
    }
    
    // Average xdips and ydips 
    double xdipaverage = accumulate(xdips.begin(), xdips.end(), 0.0)/xdips.size();
    double ydipaverage = accumulate(ydips.begin(), ydips.end(), 0.0)/ydips.size();
    
    stringstream outputtmp;
    
    // Print a table of all computed values:
    outputtmp << "###############################################################################" << endl;
    outputtmp << "# Results from upscaling dips."<< endl;
    outputtmp << "#" << endl;
    time_t now = time(NULL);
    outputtmp << "# Finished: " << asctime(localtime(&now));
    utsname hostname; uname(&hostname);
    
    outputtmp << "#" << endl;
    outputtmp << "# Eclipse file: " << gridfilename << endl;
    outputtmp << "#" << endl;
    outputtmp << "# Options used:" << endl;
    outputtmp << "#    mincellvolume: " << minCellVolume << endl;
    outputtmp << "# " << endl;
    outputtmp << "# Cells included in average dip: " << xdips.size() << endl;
    outputtmp << "# Cells smaller than mincellvolume: " << ignoredCellCount << endl;
    outputtmp << "###########################################################" << endl;
    if (listallcells) {
	outputtmp << "# i j k xdip ydip  cellvolume" << endl;
	for (unsigned int i=0; i < xdips.size(); ++i) {
	    outputtmp << cellidxs_i[i] << " " << cellidxs_j[i] << " " << cellidxs_k[i] << "\t" << xdips[i] << "\t" << ydips[i] << "\t" << cellvolumes[i] << endl;
	} 
    }
    else { 
	outputtmp << "x_dip_average " << xdipaverage << endl;
	outputtmp << "y_dip_average " << ydipaverage << endl;
    } 
    
    cout << endl << outputtmp.str();
    if (outputfilename != "")  {
        cout << "Writing results to " << outputfilename << endl;
        ofstream outfile;
        outfile.open(outputfilename.c_str(), ios::out | ios::trunc);
        outfile << outputtmp.str();
        outfile.close();
    }
    return 0;
    
};
