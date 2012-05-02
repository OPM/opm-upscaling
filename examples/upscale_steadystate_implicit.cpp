//===========================================================================
//
// File: steadystate_test.cpp
//
// Created: Fri Aug 28 14:11:03 2009
//
// Author(s): Kari B. Skjerve <karbor@statoil.com>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

  This file is part of The Open Reservoir Simulator Project (OpenRS).

  OpenRS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OpenRS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
*/


//#define VERBOSE
#include <dune/upscaling/SteadyStateUpscalerImplicit.hpp>
#include <dune/upscaling/SteadyStateUpscalerManagerImplicit.hpp>
#include <dune/porsol/euler/EulerUpstreamImplicit.hpp>
#include <dune/porsol/common/SimulatorTraits.hpp>
#include <opm/core/utility/MonotCubicInterpolator.hpp>
#include <dune/upscaling/SinglePhaseUpscaler.hpp>
namespace Dune{
	template <class IsotropyPolicy>
    struct Implicit
    {
        template <class GridInterface, class BoundaryConditions>
        struct TransportSolver
        {
            //enum { Dimension = GridInterface::Dimension };
        	enum { Dimension = GridInterface::Dimension };
            typedef typename IsotropyPolicy::template ResProp<Dimension>::Type RP;

            typedef EulerUpstreamImplicit<GridInterface,
                                  RP,
                                  BoundaryConditions> Type;

        };
    };
	typedef SimulatorTraits<Isotropic, Implicit> UpscalingTraitsBasicImplicit;
}
using namespace Dune;

void usage()
{
    std::cout << "Usage: upscale_steadystate_implicit gridfilename=filename.grdecl  " << std::endl;
    std::cout << "       rock_list=rocklist.txt [outputWater=] [outputOil=] " << std::endl;
    std::cout << "       [bc=fixed] [min_sat] [gravity=0] [minperm=1e-9] " << std::endl;
    std::cout << "       [anisotropicrocks=false]" << std::endl;
}

void usageandexit() {
    usage();
    exit(1);
}
std::vector<std::vector<double> > getExtremeSats(std::string rock_list, bool anisorocks=false) {
    if (!rock_list.compare("no_list")) {
        std::cout << "Need rock_list to compute saturation limits (" << rock_list << ")" << std::endl;
        throw std::exception();
    }
    std::ifstream rl(rock_list.c_str());
    if (!rl) {
        THROW("Could not open file " << rock_list);
    }
    int num_rocks = -1;
    rl >> num_rocks;
    ASSERT(num_rocks >= 1);
    std::vector<std::vector<double> > rocksatendp;
    rocksatendp.resize(num_rocks);
    for (int i = 0; i < num_rocks; ++i) {
        rocksatendp[i].resize(2);
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
        
        if (! anisorocks) { //Isotropic input rocks (Sw Krw Kro J)
            MonotCubicInterpolator Jtmp;
            try {
                Jtmp = MonotCubicInterpolator(rockname, 1, 4); 
            }
            catch (const char * errormessage) {
                std::cerr << "Error: " << errormessage << std::endl;
                std::cerr << "Check filename" << std::endl;
                exit(1);
            }
            rocksatendp[i][0] = Jtmp.getMinimumX().first;
            rocksatendp[i][1] = Jtmp.getMaximumX().first;
            if (rocksatendp[i][0] < 0 || rocksatendp[i][0] > 1) {
                THROW("Minimum rock saturation (" << rocksatendp[i][0] << ") not sane for rock " 
                      << rockfilename << "." << std::endl << "Did you forget to specify anisotropicrocks=true ?");  
            }
        }
        else { //Anisotropic input rocks (Pc Sw Krxx Kryy Krzz)
            MonotCubicInterpolator Pctmp;
            try {
                Pctmp = MonotCubicInterpolator(rockname, 2, 1);
            }
            catch (const char * errormessage) {
                std::cerr << "Error: " << errormessage << std::endl;
                std::cerr << "Check filename and columns 1 and 2 (Pc and Sw)" << std::endl;
                exit(1);
            }
            rocksatendp[i][0] = Pctmp.getMinimumX().first;
            rocksatendp[i][1] = Pctmp.getMaximumX().first;
        }          
    }
    return rocksatendp;
}

template <typename T>
std::string toString(T const& value) {
    std::stringstream sstr;
    sstr << value;
    return sstr.str();
}


int main(int argc, char** argv)
{
    if (argc == 1) {
        usageandexit();
    }
    // Initialize.
    Opm::parameter::ParameterGroup param(argc, argv);
    std::string gridfilename = param.get<std::string>("gridfilename");
    Opm::EclipseGridParser eclparser(gridfilename, false);

    // Check that we have the information we need from the eclipse file:  
    if (! (eclparser.hasField("SPECGRID") && eclparser.hasField("COORD") && eclparser.hasField("ZCORN")  
           && eclparser.hasField("PORO") && eclparser.hasField("PERMX"))) {
        std::cerr << "Error: Did not find SPECGRID, COORD, ZCORN, PORO and PERMX in Eclipse file " << gridfilename << std::endl;  
        usageandexit();  
    }  

    // Set default values if not given as input
    double linsolver_tolerance = param.getDefault("residual_tolerance", 1e-8);
    int linsolver_verbosity = param.getDefault("linsolver_verbosity", 0);
    int linsolver_type = param.getDefault("linsolver_type", 1);
    std::string bc=param.getDefault<std::string>("bc","fixed");
    double gravity = param.getDefault("gravity", 0);
    double surfaceTension = param.getDefault("surfaceTension", 11);
    int num_sats = param.getDefault("num_sats", 10);
    int num_pdrops = param.getDefault("num_pdrops", 10);
    double min_pdrop = param.getDefault("min_pdrop", 1e2);
    double max_pdrop = param.getDefault("max_pdrop", 1e6);
    std::string flowdir = param.getDefault<std::string>("flowdir","x");
    if (flowdir == "x") {
        param.insertParameter("flow_direction", "0");
    }
    else if (flowdir == "y") {
        param.insertParameter("flow_direction", "1");
    }
    else if (flowdir == "z") {
        param.insertParameter("flow_direction", "2");
    }
    else {
        std::cerr << "flowdir " << flowdir << " not valid. flowdir must be either x, y or z" << std::endl;
        usageandexit();
    }
    bool start_from_cl = param.getDefault("start_from_cl", true);
    if (!param.has("outputWater")) {
        param.insertParameter("outputWater","");
    }
    if (!param.has("outputOil")) {
        param.insertParameter("outputOil","");
    }

    // Compute minimum and maximum (large scale) saturations
    std::string rock_list = param.getDefault<std::string>("rock_list", "no_list");
    std::vector<std::vector<double> > rocksatendpoints_ = getExtremeSats(rock_list);
    
    std::vector<double>  poros = eclparser.getFloatingPointValue("PORO");  
    // Anisotropic relperm not yet implemented in steadystate_implicit
    //bool anisorocks = param.getDefault("anisotropicrocks", false);
    std::vector<int> satnums(poros.size(), 1); 
    if (eclparser.hasField("SATNUM")) { 
        satnums = eclparser.getIntegerValue("SATNUM"); 
    } 
    else if (eclparser.hasField("ROCKTYPE")) { 
        satnums = eclparser.getIntegerValue("ROCKTYPE"); 
    } 
    else { 
        std::cout << "Warning: SATNUM or ROCKTYPE not found in input file, assuming only one rocktype" << std::endl; 
    } 
    // check that number of rock types in rock_list matches number of rock types in grid
    int num_rock_types_grid = int(*(max_element(satnums.begin(), satnums.end())));
    int num_rock_types_rocklist = rocksatendpoints_.size();
    if (num_rock_types_grid != num_rock_types_rocklist) {
        std::cerr << "Number of rock types in rock_list " << rock_list << "(" << num_rock_types_rocklist <<  ") does not match" << std::endl
                  << "the number of rock types in grid(" << num_rock_types_grid << ")." << std::endl;
        usageandexit();
    }
    Dune::SinglePhaseUpscaler upscaler; // needed to access porosities and cell volumes
    upscaler.init(eclparser, Dune::SinglePhaseUpscaler::Fixed,
                  0.0,0.0, linsolver_tolerance, linsolver_verbosity, linsolver_type, false);
    std::vector<double>  cellPoreVolumes; 
    cellPoreVolumes.resize(satnums.size(), 0.0);
    double swirvolume = 0.0;
    double sworvolume = 0.0;
    // cell_idx is the eclipse index.
    const std::vector<int>& ecl_idx = upscaler.grid().globalCell();
    CpGrid::Codim<0>::LeafIterator c = upscaler.grid().leafbegin<0>();
    for (; c != upscaler.grid().leafend<0>(); ++c) {
        unsigned int cell_idx = ecl_idx[c->index()];
        if (satnums[cell_idx] > 0) { // Satnum zero is "no rock"
            cellPoreVolumes[cell_idx] = c->geometry().volume() * poros[cell_idx];
            swirvolume += rocksatendpoints_[int(satnums[cell_idx])-1][0] * cellPoreVolumes[cell_idx];
            sworvolume += rocksatendpoints_[int(satnums[cell_idx])-1][1] * cellPoreVolumes[cell_idx];
        }
    }
    // Total porevolume and total volume -> upscaled porosity:
    double poreVolume = std::accumulate(cellPoreVolumes.begin(), 
                                        cellPoreVolumes.end(),
                                        0.0);
    double Swir = swirvolume/poreVolume;
    double Swor = sworvolume/poreVolume;
    // Insert computed Swir and Swor as min_sat and max_sat in param object
    param.insertParameter("min_sat",toString(Swir));
    param.insertParameter("max_sat",toString(Swor));
    
    
    
   

    
    // MPIHelper::instance(argc,argv) ;
    typedef SteadyStateUpscalerImplicit<UpscalingTraitsBasicImplicit> upscaler_t;
    SteadyStateUpscalerManagerImplicit<upscaler_t> mgr;
    mgr.upscale(param);
    
}

