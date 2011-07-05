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


#include <dune/common/CornerpointChopper.hpp>
#include <dune/upscaling/SinglePhaseUpscaler.hpp>
#include <dune/porsol/common/setupBoundaryConditions.hpp>
#include <dune/common/Units.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include <ios>
#include <iomanip>
#include <sys/utsname.h>
#include <ctime>
#include <sstream>
#include <fstream>
#include <iostream>


int main(int argc, char** argv)
{
    if (argc == 1) {
        std::cout << "Usage: cpchop gridfilename=filename.grdecl [subsamples=10] [ilen=5] [jlen=5] " << std::endl;
        std::cout << "       [zlen=5] [imin=] [imax=] [jmin=] [jmax=] [upscale=true] [bc=fixed]" << std::endl;
        std::cout << "       [resettoorigin=true] [seed=111] [z_tolerance=0.0] [minperm=1e-9] " << std::endl;
        std::cout << "       [dips=false] [satnumvolumes=false] [mincellvolume=1e-9]" << std::endl;
        std::cout << "       [filebase=] [resultfile=] [endpoints=false] [rock_list=]" << std::endl;
        exit(1);
    }
    Dune::parameter::ParameterGroup param(argc, argv);
    std::string gridfilename = param.get<std::string>("gridfilename");
    Dune::CornerPointChopper ch(gridfilename);

    // The cells with i coordinate in [imin, imax) are included, similar for j.
    // The z limits may be changed inside the chopper to match actual min/max z.
    const int* dims = ch.dimensions();
    int imin = param.getDefault("imin", 0);
    int imax = param.getDefault("imax", dims[0]);
    int jmin = param.getDefault("jmin", 0);
    int jmax = param.getDefault("jmax", dims[1]);
    double zmin = param.getDefault("zmin", ch.zLimits().first);
    double zmax = param.getDefault("zmax", ch.zLimits().second);
    int subsamples = param.getDefault("subsamples", 1);
    int ilen = param.getDefault("ilen", imax - imin);
    int jlen = param.getDefault("jlen", jmax - jmin);
    double zlen = param.getDefault("zlen", zmax - zmin);
    bool upscale = param.getDefault("upscale", true);
    std::string bc = param.getDefault<std::string>("bc", "fixed");
    bool resettoorigin = param.getDefault("resettoorigin", true);
    boost::mt19937::result_type userseed = param.getDefault("seed", 0);

    int outputprecision = param.getDefault("outputprecision", 8);
    std::string filebase = param.getDefault<std::string>("filebase", "");
    std::string resultfile = param.getDefault<std::string>("resultfile", "");

    double minperm = param.getDefault("minperm", 1e-9);
    double minpermSI = Dune::unit::convert::from(minperm, Dune::prefix::milli*Dune::unit::darcy);

    // Following two options are for dip upscaling (slope of cell top and bottom edges)
    bool dips = param.getDefault("dips", false);  // whether to do dip averaging
    double mincellvolume = param.getDefault("mincellvolume", 1e-9); // ignore smaller cells for dip calculations

    bool satnumvolumes = param.getDefault("satnumvolumes", false); // whether to count volumes pr. satnum

    bool endpoints = param.getDefault("endpoints", false); // whether to upscale saturation endpoints
    std::string rock_list =  param.getDefault<std::string>("rock_list", "no_list");
    //std::vector<Dune::RockAnisotropicRelperm> rock_;
    std::vector<std::vector<double> > rocksatendpoints_;
    if (endpoints) {
        if (!rock_list.compare("no_list")) {
            std::cout << "Can't do endponts without rock list (" << rock_list << ")" << std::endl;
            throw std::exception();
        }
        // Code copied from ReservoirPropertyCommon.hpp for file reading
        std::ifstream rl(rock_list.c_str());
        if (!rl) {
            THROW("Could not open file " << rock_list);
        }
        int num_rocks = -1;
        rl >> num_rocks;
        ASSERT(num_rocks >= 1);
        rocksatendpoints_.resize(num_rocks);
        
        for (int i = 0; i < num_rocks; ++i) {
            std::string spec;
            while (spec.empty()) {
                std::getline(rl, spec);
            }
            // Read the contents of the i'th rock
            std::istringstream specstream(spec);
	    std::string rockname;
	    specstream >> rockname;
            std::string rockfilename = rockname;
            std::ifstream rock_stream(rockfilename.c_str());
            if (!rock_stream) {
                THROW("Could not open file " << rockfilename);
            }
                        char c = rock_stream.peek();
            while ((c == '-' || c == '#') && rock_stream) {
                if (c == '#') {
                    std::string commentline;
                    std::getline(rock_stream, commentline);
                    c = rock_stream.peek();
                }
                else { /* Check if the initial '-' was first byte of the sequence '--' */
                    rock_stream.get(c);
                    if (c == '-') {
                        std::string commentline;
                        std::getline(rock_stream, commentline);
                        c = rock_stream.peek();
                    } else {
                        rock_stream.putback(c);
                    }
                }
            }
            if (!rock_stream) {
                THROW("Something went wrong while skipping optional comment header of rock file.");
            }
	    typedef Dune::FieldVector<double, 2> Data;
	    std::istream_iterator<Data> start(rock_stream);
	    std::istream_iterator<Data> end;
	    std::vector<Data> data(start, end);
	    if (!rock_stream.eof()) {
		THROW("Reading stopped but we're not at eof: something went wrong reading data of rock file.");
	    }
            
            double swmin = std::min(data[0][1],data[data.size()-1][1]);
            double swmax = std::max(data[0][1],data[data.size()-1][1]);
            rocksatendpoints_[i][0] = swmin;
            rocksatendpoints_[i][1] = swmax;


            // This part might be useful if one wants to implement relative permeability in cpchop
// 	    std::vector<double> pc, sw;
// 	    for (int i = 0; i < int(data.size()); ++i) {
// 		sw.push_back(data[i][1]);
// 		pc.push_back(data[i][0]);
// 	    }

//             typedef Dune::utils::NonuniformTableLinear<double> TabFunc;
//             TabFunc pc_ = TabFunc(sw, pc);
            
            //rocksatendpoints_[i] = 
            //rock_[i].read(dir, spec);
        }
    }


    double z_tolerance = param.getDefault("z_tolerance", 0.0);
    double residual_tolerance = param.getDefault("residual_tolerance", 1e-8);
    int linsolver_verbosity = param.getDefault("linsolver_verbosity", 0);
    int linsolver_type = param.getDefault("linsolver_type", 1);

    // Check that we do not have any user input
    // that goes outside the coordinates described in
    // the cornerpoint file (runtime-exception will be thrown in case of error)
    ch.verifyInscribedShoebox(imin, ilen, imax,
			      jmin, jlen, jmax,
			      zmin, zlen, zmax);

    // Random number generator from boost.
    boost::mt19937 gen;

    // Seed the random number generators with the current time, unless specified on command line
    // Warning: Current code does not allow 0 for the seed!!
    boost::mt19937::result_type autoseed = time(NULL);
    if (userseed == 0) {
        gen.seed(autoseed);
    }
    else {
        gen.seed(userseed);
    }

    Dune::SinglePhaseUpscaler::BoundaryConditionType bctype = Dune::SinglePhaseUpscaler::Fixed;
    bool isFixed, isPeriodic;
    isFixed = isPeriodic = false;
    if (upscale) {
        if (bc == "fixed") {
            isFixed = true;
            bctype = Dune::SinglePhaseUpscaler::Fixed;
        }
        else if (bc == "periodic") {
            isPeriodic = true;
            bctype = Dune::SinglePhaseUpscaler::Periodic;
        }
        else {
            std::cout << "Boundary condition type (bc=" << bc << ") not allowed." << std::endl;
            std::cout << "Only bc=fixed or bc=periodic implemented." << std::endl;
            throw std::exception();
        }
    }

    // Check for unused parameters (potential typos).
    if (param.anyUnused()) {
	std::cout << "*****     WARNING: Unused parameters:     *****\n";
	param.displayUsage();
    }
    
    // Note that end is included in interval for uniform_int.
    boost::uniform_int<> disti(imin, imax - ilen);
    boost::uniform_int<> distj(jmin, jmax - jlen);
    boost::uniform_real<> distz(zmin, std::max(zmax - zlen, zmin));
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > ri(gen, disti);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > rj(gen, distj);
    boost::variate_generator<boost::mt19937&, boost::uniform_real<> > rz(gen, distz);

    // Storage for results
    std::vector<double> porosities;
    std::vector<double> permxs;
    std::vector<double> permys;
    std::vector<double> permzs;
    std::vector<double> permyzs;
    std::vector<double> permxzs;
    std::vector<double> permxys;
    std::vector<double> minsws, maxsws;
    std::vector<double> xdips, ydips;

    // Initialize a matrix for subsample satnum volumes. 
    // Outer index is subsample index, inner index is SATNUM-value
    std::vector<std::vector<double> > rockvolumes; 
    int maxSatnum = 0; // This value is determined from the chopped cells.

    int finished_subsamples = 0; // keep explicit count of successful subsamples
    for (int sample = 1; sample <= subsamples; ++sample) {
        int istart = ri();
        int jstart = rj();
        double zstart = rz();
        ch.chop(istart, istart + ilen, jstart, jstart + jlen, zstart, zstart + zlen, resettoorigin);
        std::string subsampledgrdecl = filebase;

        // Output grdecl-data to file if a filebase is supplied.
        if (filebase != "") {
            std::ostringstream oss;
            if ((size_t) subsamples > 1) { // Only add number to filename if more than one sample is asked for
                oss << 'R' << std::setw(4) << std::setfill('0') << sample;
                subsampledgrdecl += oss.str();
            }
            subsampledgrdecl += ".grdecl";
            ch.writeGrdecl(subsampledgrdecl);
        }

        try { /* The upscaling may fail to converge on icky grids, lets just pass by those */
            if (upscale) {
                Dune::EclipseGridParser subparser = ch.subparser();
                subparser.convertToSI();
                Dune::SinglePhaseUpscaler upscaler;
                
                upscaler.init(subparser, bctype, minpermSI, z_tolerance,
                              residual_tolerance, linsolver_verbosity, linsolver_type, false);

                Dune::SinglePhaseUpscaler::permtensor_t upscaled_K = upscaler.upscaleSinglePhase();
                upscaled_K *= (1.0/(Dune::prefix::milli*Dune::unit::darcy));


                porosities.push_back(upscaler.upscalePorosity());
                permxs.push_back(upscaled_K(0,0));
                permys.push_back(upscaled_K(1,1));
                permzs.push_back(upscaled_K(2,2));
                permyzs.push_back(upscaled_K(1,2));
                permxzs.push_back(upscaled_K(0,2));
                permxys.push_back(upscaled_K(0,1));

            }

            if (endpoints) {
                // Calculate minimum and maximum water volume in each cell
                // Create single-phase upscaling object to get poro and perm values from the grid
                Dune::EclipseGridParser subparser = ch.subparser();
                subparser.convertToSI();
                Dune::SinglePhaseUpscaler upscaler;                
                upscaler.init(subparser, bctype, minpermSI, z_tolerance,
                              residual_tolerance, linsolver_verbosity, linsolver_type, false);
                std::vector<int>   satnums = subparser.getIntegerValue("SATNUM");
                std::vector<double>  poros = subparser.getFloatingPointValue("PORO");
                std::vector<double> cellVolumes, cellPoreVolumes;
                cellVolumes.resize(satnums.size(), 0.0);
                cellPoreVolumes.resize(satnums.size(), 0.0);
                int tesselatedCells = 0;
                double Swirvolume = 0, Sworvolume = 0;
                const std::vector<int>& ecl_idx = upscaler.grid().globalCell();
                Dune::CpGrid::Codim<0>::LeafIterator c = upscaler.grid().leafbegin<0>();
                for (; c != upscaler.grid().leafend<0>(); ++c) {
                    unsigned int cell_idx = ecl_idx[c->index()];
                    if (satnums[cell_idx] > 0) { // Satnum zero is "no rock"
                        cellVolumes[cell_idx] = c->geometry().volume();
                        cellPoreVolumes[cell_idx] = cellVolumes[cell_idx] * poros[cell_idx];
                        Swirvolume += rocksatendpoints_[satnums[cell_idx]-1][0] * cellPoreVolumes[cell_idx]; 
                        Sworvolume += rocksatendpoints_[satnums[cell_idx]-1][1] * cellPoreVolumes[cell_idx];
                        
                    }
                    ++ tesselatedCells;
                }
                
                double poreVolume = std::accumulate(cellPoreVolumes.begin(), 
                                                    cellPoreVolumes.end(),
                                                    0.0);
                minsws.push_back(Swirvolume/poreVolume);
                maxsws.push_back(Sworvolume/poreVolume);
            }

	    if (dips) {
		Dune::EclipseGridParser subparser = ch.subparser();
		std::vector<int>  griddims = subparser.getSPECGRID().dimensions;
		std::vector<double> xdips_subsample, ydips_subsample;

		Dune::EclipseGridInspector gridinspector(subparser);
		for (int k=0; k < griddims[2]; ++k) {
		    for (int j=0; j < griddims[1]; ++j) {
			for (int i=0; i < griddims[0]; ++i) {
			    if (gridinspector.cellVolumeVerticalPillars(i, j, k) > mincellvolume) {
				std::pair<double,double> xydip = gridinspector.cellDips(i, j, k);
				xdips_subsample.push_back(xydip.first);
				ydips_subsample.push_back(xydip.second);
			    }
			}
		    }
		}

		// Average xdips and ydips
		double xdipaverage = accumulate(xdips_subsample.begin(), xdips_subsample.end(), 0.0)/xdips_subsample.size();
		double ydipaverage = accumulate(ydips_subsample.begin(), ydips_subsample.end(), 0.0)/ydips_subsample.size();
		xdips.push_back(xdipaverage);
		ydips.push_back(ydipaverage);
	    }

	    if (satnumvolumes) {
		Dune::EclipseGridParser subparser = ch.subparser();
		Dune::EclipseGridInspector subinspector(subparser);
		std::vector<int>  griddims = subparser.getSPECGRID().dimensions;
		int number_of_subsamplecells = griddims[0] * griddims[1] * griddims[2];

		// If SATNUM is non-existent in input grid, this will fail:
		std::vector<int> satnums = subparser.getIntegerValue("SATNUM");

		std::vector<double> rockvolumessubsample;
		for (int cell_idx=0; cell_idx < number_of_subsamplecells; ++cell_idx) {
		    maxSatnum = std::max(maxSatnum, int(satnums[cell_idx]));
		    rockvolumessubsample.resize(maxSatnum); // Ensure long enough vector
		    rockvolumessubsample[int(satnums[cell_idx])-1] += subinspector.cellVolumeVerticalPillars(cell_idx);
		}

		// Normalize volumes to obtain relative volumes:
		double subsamplevolume = std::accumulate(rockvolumessubsample.begin(),
							 rockvolumessubsample.end(), 0.0);
		std::vector<double> rockvolumessubsample_normalized;
		for (size_t satnum_idx = 0; satnum_idx < rockvolumessubsample.size(); ++satnum_idx) {
		    rockvolumessubsample_normalized.push_back(rockvolumessubsample[satnum_idx]/subsamplevolume);
		}
		rockvolumes.push_back(rockvolumessubsample_normalized);
	    }

	    finished_subsamples++;
        }
        catch (...) {
            std::cerr << "Warning: Upscaling chopped subsample nr. " << sample << " failed, proceeding to next subsample\n";
        }

    }

    // Make stream of output data, to be outputted to screen and optionally to file
    std::stringstream outputtmp;

    outputtmp << "################################################################################################" << std::endl;
    outputtmp << "# Results from property analysis on subsamples" << std::endl;
    outputtmp << "#" << std::endl;
    time_t now = time(NULL);
    outputtmp << "# Finished: " << asctime(localtime(&now));

    utsname hostname;   uname(&hostname);
    outputtmp << "# Hostname: " << hostname.nodename << std::endl;
    outputtmp << "#" << std::endl;
    outputtmp << "# Options used:" << std::endl;
    outputtmp << "#     gridfilename: " << gridfilename << std::endl;
    outputtmp << "#   i; min,len,max: " << imin << " " << ilen << " " << imax << std::endl;
    outputtmp << "#   j; min,len,max: " << jmin << " " << jlen << " " << jmax << std::endl;
    outputtmp << "#   z; min,len,max: " << zmin << " " << zlen << " " << zmax << std::endl;
    outputtmp << "#       subsamples: " << subsamples << std::endl;
    if (userseed == 0) {
        outputtmp << "#      (auto) seed: " << autoseed << std::endl;
    }
    else {
        outputtmp << "#    (manual) seed: " << userseed << std::endl; 
    }        
    outputtmp << "################################################################################################" << std::endl;
    outputtmp << "# id";
    if (upscale) {
        if (isPeriodic) {
            outputtmp << "          porosity                 permx                   permy                   permz                 permyz                  permxz                  permxy";
        }
        else if (isFixed) {
            outputtmp << "          porosity                 permx                   permy                   permz";
        }
    }
    if (endpoints) {
        outputtmp << "                  Swir                    Swor";
    }
    if (dips) {
	outputtmp << "              xdip              ydip";
    }
    if (satnumvolumes) {
	for (int satnumidx = 0; satnumidx < maxSatnum; ++satnumidx) {
	    outputtmp << "         satnum_" << satnumidx+1;
	}
    }
    outputtmp << std::endl;

    const int fieldwidth = outputprecision + 8;
    for (int sample = 1; sample <= finished_subsamples; ++sample) {
	outputtmp << sample << '\t';
	if (upscale) {
	    outputtmp <<
		std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << porosities[sample-1] << '\t' <<
		std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << permxs[sample-1] << '\t' <<
		std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << permys[sample-1] << '\t' <<
		std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << permzs[sample-1] << '\t';
            if (isPeriodic) {
                outputtmp <<
                    std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << permyzs[sample-1] << '\t' <<
                    std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << permxzs[sample-1] << '\t' <<
                    std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << permxys[sample-1] << '\t';                
            }
	}
	if (endpoints) {
	    outputtmp <<
		std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << minsws[sample-1] << '\t' <<
		std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << maxsws[sample-1];
	}
	if (dips) {
	    outputtmp <<
		std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << xdips[sample-1] << '\t' <<
		std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << ydips[sample-1];
	}
	if (satnumvolumes) {
	    rockvolumes[sample-1].resize(maxSatnum, 0.0);
	    for (int satnumidx = 0; satnumidx < maxSatnum; ++satnumidx) {
		outputtmp <<
		    std::showpoint << std::setw(fieldwidth) <<
		    std::setprecision(outputprecision) << rockvolumes[sample-1][satnumidx] << '\t';
	    }
	}
	outputtmp <<  std::endl;
    }

    if (resultfile != "") {
	std::cout << "Writing results to " << resultfile << std::endl;
	std::ofstream outfile;
	outfile.open(resultfile.c_str(), std::ios::out | std::ios::trunc);
	outfile << outputtmp.str();
            outfile.close();
    }



    std::cout << outputtmp.str();
}

