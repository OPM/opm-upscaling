
/*
  Copyright 2010 Statoil ASA.

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

/** @file upscale_avg.C
 *  @brief Upscales using simple averages and reports statistics
 *  
 *  Processes permeability and porosity for a given Eclipse-model. 
 *  
 *  Input is Eclipse grid format specifying the corner-point
 *  grid (must be of shoebox-shape, but this condition is slightly relaxed
 *  on top and bottom surfaces).
 * 
 *  The input eclipse file must specify the permeability properties
 *  for each cell.
 * 
 *  The grid processing step from the permeability upscaling code is
 *  used as a simple code library for computing the volume of each
 *  cell, if some effort is invested, this can be circumventing by a
 *  simpler code that does nothing but interpret the coordinates and
 *  calculates the volume, which will be much faster.
 * 
 *
 */
#include <config.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <numeric> // for std::accumulate
#include <sys/utsname.h>

#include <opm/upscaling/SinglePhaseUpscaler.hpp>
#include <opm/core/io/eclipse/EclipseGridInspector.hpp>

#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
#include <dune/common/parallel/mpihelper.hh>
#else
#include <dune/common/mpihelper.hh>
#endif

using namespace Opm;
using namespace std;

void usage() {
    cout << endl <<
        "Usage: upscale_avg <options> <eclipsefile>" << endl <<
        "were the options are:" << endl <<
        "-use_actnum                  -- If used and given a number larger than 0, the cells with" << endl <<
        "                                ACTNUM 0 is removed from the calculations" << endl;
}
void usageandexit() {
    usage();
    exit(1);
}


/**
   @brief Computes simple statistics.
*/
int main(int varnum, char** vararg) try {        

    Dune::MPIHelper::instance(varnum, vararg);

    const double emptycellvolumecutoff = 1e-10;

    bool anisotropic_input = false;

    if (varnum ==  1) { // If no arguments supplied ("upscale_avg" is the first argument)
        cout << "Error: No eclipsefile provided" << endl;
        usage();
        exit(1);
    } 


   /*
     Populate options-map with default values
   */
   map<string,string> options;
   options.insert(make_pair("use_actnum",   "0"     )); //Use ACTNUM as indicator for which cells to be included
                                                        //Default is not to use actnum for this
                                                        //When use_actnum is used the cell volumes of the cells with actnum 0 is set to 0

   /* Loop over all command line options in order to look 
      for options. 

      argidx loops over all the arguments here, and updates the
      variable 'argeclindex' *if* it finds any legal options,
      'argeclindex' is so that vararg[argeclindex] = the eclipse
      filename. If options are illegal, argeclindex will be wrong, 
      
   */
   int argeclindex = 0;
   for (int argidx = 1; argidx < varnum; argidx += 2)  {
       if (string(vararg[argidx]).substr(0,1) == "-")    {
           string searchfor = string(vararg[argidx]).substr(1); // Chop off leading '-'
           /* Check if it is a match */
           if (options.count(searchfor) == 1) {
               options[searchfor] = string(vararg[argidx+1]);
               cout << "Parsed command line option: " << searchfor << " := " << vararg[argidx+1] << endl;
               argeclindex = argidx + 2;
           }
           else {
               cout << "Option -" << searchfor << " unrecognized." << endl;
               usageandexit();
           }
       }
       else { 
           // if vararg[argidx] does not start in '-', 
           // assume we have found the position of the Eclipse-file.
           argeclindex = argidx;
           break; // out of for-loop, 
       }
   }

    const char* ECLIPSEFILENAME(vararg[argeclindex]);

    // Test if filename exists and is readable
    ifstream eclipsefile(ECLIPSEFILENAME, ios::in);
    if (eclipsefile.fail()) {
        cerr << "Error: Filename " << ECLIPSEFILENAME << " not found or not readable." << endl;
        usage();
        exit(1);
    }
    eclipsefile.close();
    
    // Variables for timing/profiling
    clock_t start, finish;
    double timeused = 0;
        
    /***********************************************************************
    * Step X
    * Load geometry and data from Eclipse file
    */
    cout << "Parsing Eclipse file <" << ECLIPSEFILENAME << "> ... ";
    flush(cout);   start = clock();
    // eclParser_p is here a pointer to an object of type Opm::EclipseGridParser
    // (this pointer trick is necessary for the try-catch-clause to work)
    
    Opm::EclipseGridParser eclParser(ECLIPSEFILENAME, false);

    Opm::EclipseGridInspector eclInspector(eclParser);

    finish = clock();   timeused = (double(finish)-double(start))/CLOCKS_PER_SEC;
    cout << " (" << timeused <<" secs)" << endl;
    
    // Check that we have the information we need from the eclipse file, we will check PERM-fields later
    if (! (eclParser.hasField("SPECGRID") && eclParser.hasField("COORD") && eclParser.hasField("ZCORN"))) {  
        cerr << "Error: Did not find SPECGRID, COORD and ZCORN in Eclipse file " << ECLIPSEFILENAME << endl;  
        usage();  
        exit(1);  
    }

   SinglePhaseUpscaler upscaler;
   upscaler.init(eclParser, SinglePhaseUpscaler::Fixed,
                 Opm::unit::convert::from(1e-9, Opm::prefix::milli*Opm::unit::darcy),
                 0.0, 1e-8, 0, 1, false);


    bool use_actnum = false;
    vector<int> actnum;
    if (atoi(options["use_actnum"].c_str())>0) {
        if (eclParser.hasField("ACTNUM")) {
            use_actnum = true;
            actnum = eclParser.getIntegerValue("ACTNUM");
            cout << actnum[0] << " " << actnum[1] << endl;
        }
        else {
            cout << "Error: option use_actnum set to 1 but Gridfile " << ECLIPSEFILENAME << " does not include any field ACTNUM." << endl;
            usage();
            exit(1);            
        }
    }

    // Print header for the output
    
    cout << "Statistics for filename: " << ECLIPSEFILENAME << endl;
    cout << "-----------------------------------------------------" << endl;
    bool doporosity = false;
    if (eclParser.hasField("PORO")) {
        doporosity = true;
    }
    bool dontg = false;
    if (eclParser.hasField("NTG")) {
        // Ntg only used together with PORO
        if (eclParser.hasField("PORO")) dontg = true;
    }    

    bool doperm = false;
    if (eclParser.hasField("PERMX")) {
        doperm = true;
        if (eclParser.hasField("PERMY") && eclParser.hasField("PERMZ")) {
            anisotropic_input = true;
            cout << "Info: PERMY and PERMZ present in data file." << endl;
        }
    }
    
    

    // Global number of cells (includes inactive cells)
    vector<int>  griddims = eclParser.getSPECGRID().dimensions;
    int num_eclipse_cells = griddims[0] * griddims[1] * griddims[2];


    cout << "Active and inactive cells: " << num_eclipse_cells << " (" <<
      griddims[0] << " x " <<
      griddims[1] << " x " <<
      griddims[2] << ")" << endl;
    int pillars = (griddims[0]+1) * (griddims[1]+1);
    cout << "                  Pillars: " << pillars << " (" << griddims[0]+1 << 
        " x " << griddims[1]+1 << ")" << endl;
    
    // Find max and min in x-, y- and z-directions
    std::array<double, 6> gridlimits = eclInspector.getGridLimits();
    cout << "                 x-limits: " << gridlimits[0] << " -- " << gridlimits[1] << endl;
    cout << "                 y-limits: " << gridlimits[2] << " -- " << gridlimits[3] << endl;
    cout << "                 z-limits: " << gridlimits[4] << " -- " << gridlimits[5] << endl;

    // First do overall statistics
    vector<double> cellVolumes, cellPoreVolumes, netCellVolumes, netCellPoreVolumes;
    cellVolumes.resize(num_eclipse_cells, 0.0);
    cellPoreVolumes.resize(num_eclipse_cells, 0.0);
    netCellVolumes.resize(num_eclipse_cells, 0.0);
    netCellPoreVolumes.resize(num_eclipse_cells, 0.0);
    int active_cell_count = 0;
    vector<double> poros, ntgs;
    vector<double> permxs, permys, permzs;
    if (doporosity) {
        poros = eclParser.getFloatingPointValue("PORO");
    }
    if (dontg) {
        ntgs = eclParser.getFloatingPointValue("NTG");
    }
    if (doperm) {
        permxs = eclParser.getFloatingPointValue("PERMX");
        if (anisotropic_input) {
            permys = eclParser.getFloatingPointValue("PERMY");
            permzs = eclParser.getFloatingPointValue("PERMZ");
        }
        
    }
    const std::vector<int>& ecl_idx = upscaler.grid().globalCell();
    Dune::CpGrid::Codim<0>::LeafIterator c = upscaler.grid().leafbegin<0>();
    for (; c != upscaler.grid().leafend<0>(); ++c) {
      size_t cell_idx = ecl_idx[c->index()];
      
      //for (size_t cell_idx = 0; cell_idx < num_eclipse_cells; ++cell_idx) {
        if (!use_actnum){
          cellVolumes[cell_idx] = c->geometry().volume();
        }
        else {
          cellVolumes[cell_idx] = c->geometry().volume() * actnum[cell_idx];
        }
        if (cellVolumes[cell_idx] > emptycellvolumecutoff) {
            ++active_cell_count;
            if (doporosity) {
                cellPoreVolumes[cell_idx] = cellVolumes[cell_idx] * poros[cell_idx];
            }
            if (dontg) {
                netCellPoreVolumes[cell_idx] = cellVolumes[cell_idx] * poros[cell_idx] * ntgs[cell_idx];
                netCellVolumes[cell_idx] = cellVolumes[cell_idx] * ntgs[cell_idx];
            }
        }
    }
    
    cout << "             Active cells: " << active_cell_count << " (" << (double)active_cell_count/(double)num_eclipse_cells*100.0 << "%)" << endl;
    double volume = std::accumulate(cellVolumes.begin(),
                                    cellVolumes.end(),
                                    0.0);
    cout << "             Total volume: " << volume << endl;
    if (doporosity) {
        double poreVolume = std::accumulate(cellPoreVolumes.begin(), 
                                            cellPoreVolumes.end(),
                                            0.0);
        cout << "         Total porevolume: " << poreVolume << endl;
        cout << "        Upscaled porosity: " << poreVolume/volume << endl;
        
        int zeroporocells = 0;
        int negativeporocells = 0;
        for (size_t cell_idx = 0; cell_idx < (size_t)num_eclipse_cells; ++cell_idx) {
            if  (poros[cell_idx] == 0) {
                ++zeroporocells;
            }
            if (poros[cell_idx] < 0 ) {
                ++negativeporocells;
            }
        }
        if  (zeroporocells > 0) {
            cout << " Cells with zero porosity: " << zeroporocells << endl;
        }
        if  (negativeporocells > 0) {
            cout << "Cells with negative porosity: " << negativeporocells << endl;
        }
    }
    if (dontg) {
        double netVolume = std::accumulate(netCellVolumes.begin(),
                                               netCellVolumes.end(),
                                               0.0);
        cout << "         Total net volume: " << netVolume << endl;
        cout << "             Upscaled NTG: " << netVolume/volume << endl;        
        double netPoreVolume = std::accumulate(netCellPoreVolumes.begin(), 
                                               netCellPoreVolumes.end(),
                                               0.0);
        cout << "     Total net porevolume: " << netPoreVolume << endl;
        cout << "    Upscaled net porosity: " << netPoreVolume/netVolume << endl;        
    }
    
    double permxsum = 0.0, permysum = 0.0, permzsum = 0.0;
    double invpermxsum = 0.0, invpermysum = 0.0, invpermzsum = 0.0;
    double volpermxsum = 0.0, volpermysum = 0.0, volpermzsum = 0.0;
    double invvolpermxsum = 0.0, invvolpermysum = 0.0, invvolpermzsum = 0.0;
    double logpermxsum = 0.0, logpermysum = 0.0, logpermzsum = 0.0;
    double logvolpermxsum = 0.0, logvolpermysum = 0.0, logvolpermzsum=0.0; 
    if (doperm) {
        int zeropermcells = 0;
        int negativepermcells = 0;
        for (size_t cell_idx = 0; cell_idx < (size_t)num_eclipse_cells; ++cell_idx) {
            if (cellVolumes[cell_idx] > emptycellvolumecutoff) {
                permxsum += permxs[cell_idx];
                volpermxsum += permxs[cell_idx] * cellVolumes[cell_idx];
                invpermxsum += 1.0/permxs[cell_idx];
                invvolpermxsum += cellVolumes[cell_idx] / permxs[cell_idx];
                logpermxsum += log(permxs[cell_idx]);
                logvolpermxsum += log(permxs[cell_idx]) * cellVolumes[cell_idx];
                if (permxs[cell_idx] == 0) {
                    ++zeropermcells;
                }
                if (permxs[cell_idx] < 0) {
                    ++negativepermcells;
                }
                if (anisotropic_input) {
                    permysum += permys[cell_idx];
                    volpermysum += permys[cell_idx] * cellVolumes[cell_idx];
                    invpermysum += 1.0/permys[cell_idx];
                    invvolpermysum += cellVolumes[cell_idx] / permys[cell_idx];
                    logpermysum += log(permys[cell_idx]);
                    logvolpermysum += log(permys[cell_idx]) * cellVolumes[cell_idx];
                    
                    permzsum += permzs[cell_idx];
                    volpermzsum += permzs[cell_idx] * cellVolumes[cell_idx];
                    invpermzsum += 1.0/permzs[cell_idx];
                    invvolpermzsum += cellVolumes[cell_idx] / permzs[cell_idx];
                    logpermzsum += log(permzs[cell_idx]);
                    logvolpermzsum += log(permzs[cell_idx]) * cellVolumes[cell_idx];
                }
            }
        }
        
        cout << "Total arithmetic permeability average: " << volpermxsum/volume << endl;
        cout << "  Total harmonic permeability average: " << volume/invvolpermxsum << endl;
        cout << " Total geometric permeability average: " << exp(logvolpermxsum/volume) << endl;
        cout << "Total arithmetic permeability average: " << permxsum/((double)active_cell_count) << " (not volume-weighted)" << endl;
        cout << "  Total harmonic permeability average: " << ((double)active_cell_count)/invpermxsum << " (not volume-weighted)" << endl;
        cout << " Total geometric permeability average: " << exp(logpermxsum/(double)active_cell_count) << " (not volume-weighted)" << endl;
        
        if (anisotropic_input) {
            cout << endl;
            cout << "Total arithmetic permeability (y) average: " << volpermysum/volume << endl;
            cout << "  Total harmonic permeability (y) average: " << volume/invvolpermysum << endl;
            cout << " Total geometric permeability (y) average: " << exp(logvolpermysum/volume) << endl;
            cout << "Total arithmetic permeability (y) average: " << permysum/((double)active_cell_count) << " (not volume-weighted)" << endl;
            cout << "  Total harmonic permeability (y) average: " << ((double)active_cell_count)/invpermysum << " (not volume-weighted)" << endl;
            cout << " Total geometric permeability (y) average: " << exp(logpermysum/(double)active_cell_count) << " (not volume-weighted)" << endl;
            cout << endl;
            cout << "Total arithmetic permeability (z) average: " << volpermzsum/volume << endl;
            cout << "  Total harmonic permeability (z) average: " << volume/invvolpermzsum << endl;
            cout << " Total geometric permeability (z) average: " << exp(logvolpermzsum/volume) << endl;
            cout << "Total arithmetic permeability (z) average: " << permzsum/((double)active_cell_count) << " (not volume-weighted)" << endl;
            cout << "  Total harmonic permeability (z) average: " << ((double)active_cell_count)/invpermzsum << " (not volume-weighted)" << endl;
            cout << " Total geometric permeability (z) average: " << exp(logpermzsum/(double)active_cell_count) << " (not volume-weighted)" << endl;
            
        }

        if  (zeropermcells > 0) {
            cout << endl <<  "     Cells with zero (x) permeability: " << zeropermcells << endl;
        }
        if  (negativepermcells > 0) {
            cout << " Cells with negative (x) permeability: " << negativepermcells << endl;
        }
    }

    // Then do statistics on rocktype by rocktype basis    
    bool dosatnums = false;
    vector<int> satnums;
    if (eclParser.hasField("SATNUM")) {
        dosatnums = true;
        satnums = eclParser.getIntegerValue("SATNUM");
    } // If SATNUM was not present, maybe ROCKTYPE is there, 
    // if so, we will use it as SATNUM.
    else if (eclParser.hasField("ROCKTYPE")) {
        dosatnums = true;
        satnums = eclParser.getIntegerValue("ROCKTYPE");
    }


    if (dosatnums) {
        int maxsatnumvalue = 0;
        // Check that SATNUM are set sensibly, that is > 0 and < 1000, and find number
        // of unique satnums present ( = number of rocktypes)
        for (size_t i = 0; i < (size_t)satnums.size(); ++i) {
            if (satnums[i] > maxsatnumvalue) {
                maxsatnumvalue = satnums[i];
            }
            if (satnums[i] < 0 || satnums[i] > 1000) { 
                cerr << "satnums[" << i << "] = " << satnums[i] << ", not sane, quitting." << endl;
                exit(1);
            }
        }

        vector<double> permxsum_rocktype;
        permxsum_rocktype.resize(maxsatnumvalue+1, 0.0);
        
        vector<double> invpermxsum_rocktype;
        invpermxsum_rocktype.resize(maxsatnumvalue+1, 0.0); // for harmonic average

        vector<double> volpermxsum_rocktype;              // volume weighted
        volpermxsum_rocktype.resize(maxsatnumvalue+1, 0.0);
        
        vector<double> invvolpermxsum_rocktype;           // volume weighted    
        invvolpermxsum_rocktype.resize(maxsatnumvalue+1, 0.0); // for harmonic average

        vector<double> porevolumesum_rocktype;
        porevolumesum_rocktype.resize(maxsatnumvalue+1, 0.0);

        vector<double> porosityvariancesum_rocktype; // for estimation of porosity variance within rocktype
        porosityvariancesum_rocktype.resize(maxsatnumvalue+1, 0.0);

        vector<double> volumesum_rocktype;
        volumesum_rocktype.resize(maxsatnumvalue+1, 0.0);

        vector<int> totalcellcount_rocktype;
        totalcellcount_rocktype.resize(maxsatnumvalue+1, 0);

        vector<int> activecellcount_rocktype;
        activecellcount_rocktype.resize(maxsatnumvalue+1, 0);

        // Now loop over cells to collect statistics pr. rocktype
        for (size_t cell_idx = 0; cell_idx < (size_t)num_eclipse_cells; ++cell_idx) {
            ++totalcellcount_rocktype[satnums[cell_idx]];
            if (cellVolumes[cell_idx] > emptycellvolumecutoff) {
                ++activecellcount_rocktype[satnums[cell_idx]];
                volumesum_rocktype[satnums[cell_idx]] += cellVolumes[cell_idx];
                if (doperm) {
                    permxsum_rocktype[satnums[cell_idx]] += permxs[cell_idx];
                    invpermxsum_rocktype[satnums[cell_idx]] += 1.0/permxs[cell_idx];

                    volpermxsum_rocktype[satnums[cell_idx]] += permxs[cell_idx] * cellVolumes[cell_idx];
                    invvolpermxsum_rocktype[satnums[cell_idx]] += cellVolumes[cell_idx] / permxs[cell_idx];
                    
                }
                if (doporosity) {
                    porevolumesum_rocktype[satnums[cell_idx]] += cellVolumes[cell_idx] * poros[cell_idx];
                }
            }
        }

        // Compute the sample variance in porosity per rock type        
        /*
          phi_avg = 1/V*sum(v_i*phi_i)
          phi_var = 1/V*sum(v_i*(phi_i-phi_avg)^2)
        */
        if (doporosity) {
            for (size_t cell_idx = 0; cell_idx < (size_t)num_eclipse_cells; ++cell_idx) {
                if (cellVolumes[cell_idx] > emptycellvolumecutoff) {
                    porosityvariancesum_rocktype[satnums[cell_idx]] += cellVolumes[cell_idx] 
                        * pow((poros[cell_idx]-porevolumesum_rocktype[satnums[cell_idx]]/volumesum_rocktype[satnums[cell_idx]]),2);
                                       
                }
            }  
        }      
        
        // Now loop over rocktypes in order to print statistics

        // Does SATNUM/ROCKTYPE start with 0 or 1?
        for (int rocktype_idx = 0; rocktype_idx <= maxsatnumvalue; ++rocktype_idx) {
            if (activecellcount_rocktype[rocktype_idx] > 0) {
                cout << endl << "Statistics for rocktype " << rocktype_idx << endl;
                
                cout     << "         Volume: " << volumesum_rocktype[rocktype_idx] << " (" << volumesum_rocktype[rocktype_idx]/volume*100 << "%)" << endl;
                cout     << "    Total cells: " << totalcellcount_rocktype[rocktype_idx] << " (" << (double)totalcellcount_rocktype[rocktype_idx]/(double)num_eclipse_cells*100 << "%)" << endl;
                cout     << "   Active cells: " << activecellcount_rocktype[rocktype_idx] << " (" << 
                    (double)activecellcount_rocktype[rocktype_idx]/(double)totalcellcount_rocktype[rocktype_idx]*100.0 << "% within rocktype)" << endl;
                if (doporosity) {
                    cout << "     Porevolume: " << porevolumesum_rocktype[rocktype_idx] << endl;
                    cout << "       Porosity: " << porevolumesum_rocktype[rocktype_idx]/volumesum_rocktype[rocktype_idx] << endl;
                    cout << "   Porosity var: " << porosityvariancesum_rocktype[rocktype_idx]/volumesum_rocktype[rocktype_idx] << endl;
                }
                if (doperm) {
                    cout << "Perm arith. avg: " << volpermxsum_rocktype[rocktype_idx]/volumesum_rocktype[rocktype_idx] << endl;
                    cout << " Perm harm. avg: " << volumesum_rocktype[rocktype_idx]/invvolpermxsum_rocktype[rocktype_idx] << endl;
                    cout << "Perm arith. avg: " << permxsum_rocktype[rocktype_idx]/activecellcount_rocktype[rocktype_idx] << " (not volume-weighted)" << endl;
                    cout << " Perm harm. avg: " << activecellcount_rocktype[rocktype_idx]/invpermxsum_rocktype[rocktype_idx] << " (not volume-weighted)" << endl;
                }                
            }
        }
    }
    
    
    return 0;
    
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}

