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
#include <config.h>


//#define VERBOSE
#include <opm/upscaling/SteadyStateUpscalerImplicit.hpp>
#include <opm/upscaling/SteadyStateUpscalerManagerImplicit.hpp>
#include <opm/porsol/euler/EulerUpstreamImplicit.hpp>
#include <opm/porsol/common/SimulatorTraits.hpp>
#include <opm/core/utility/MonotCubicInterpolator.hpp>
#include <opm/upscaling/SinglePhaseUpscaler.hpp>
#include <sys/utsname.h>
namespace Opm{
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
using namespace Opm;

void usage()
{
    std::cout << "Usage: upscale_steadystate_implicit gridfilename=filename.grdecl  " << std::endl;
    std::cout << "       rock_list=rocklist.txt [outputWater=] [outputOil=] " << std::endl;
    std::cout << "       [bc=fixed] [num_sats=10] [num_pdrops=10] " << std::endl;
    std::cout << "       [anisotropicrocks=false]" << std::endl;
}

void usageandexit() {
    usage();
    exit(1);
}

// Assumes that permtensor_t use C ordering.
double getVoigtValue(const SinglePhaseUpscaler::permtensor_t& K, int voigt_idx)
{
    ASSERT(K.numRows() == 3 && K.numCols() == 3);
    switch (voigt_idx) {
    case 0: return K.data()[0];
    case 1: return K.data()[4];
    case 2: return K.data()[8];
    case 3: return K.data()[5];
    case 4: return K.data()[2];
    case 5: return K.data()[1];
    case 6: return K.data()[7];
    case 7: return K.data()[6];
    case 8: return K.data()[3];
    default:
        std::cout << "Voigt index out of bounds (only 0-8 allowed)" << std::endl;
        throw std::exception();
    }
}

std::vector<std::vector<double> > getExtremeSats(std::string rock_list, std::vector<std::string>& rockfilelist, bool anisorocks=false) {
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
        rockfilelist.push_back(rockfilename);
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
    int linsolver_maxit = param.getDefault("linsolver_max_iterations", 0);
    int linsolver_smooth_steps = param.getDefault("linsolver_smooth_steps", 2);
    double linsolver_prolongate_factor = param.getDefault("linsolver_prolongate_factor", 1.6);
    std::string bc=param.getDefault<std::string>("bc","fixed");
    //double gravity = param.getDefault("gravity", 0);
    double surfaceTension = param.getDefault("surfaceTension", 11);
    int num_sats = param.getDefault("num_sats", 10);
    int num_pdrops = param.getDefault("num_pdrops", 10);
    double log_min_pdrop = std::log(param.getDefault("min_pdrop", 1e2));
    double log_max_pdrop = std::log(param.getDefault("max_pdrop", 1e6));
    std::string flowdir = param.getDefault<std::string>("flowdir","x");
    int fldir;
    if (flowdir == "x") {
        //param.insertParameter("flow_direction", "0");
        fldir = 0;
    }
    else if (flowdir == "y") {
        //param.insertParameter("flow_direction", "1");
        fldir = 1;
    }
    else if (flowdir == "z") {
        //param.insertParameter("flow_direction", "2");
        fldir = 2;
    }
    else {
        std::cerr << "flowdir " << flowdir << " not valid. flowdir must be either x, y or z" << std::endl;
        usageandexit();
    }
    bool start_from_cl = param.getDefault("start_from_cl", true);
    //if (!param.has("outputWater")) {
    //    param.insertParameter("outputWater","");
    //}
    //if (!param.has("outputOil")) {
    //    param.insertParameter("outputOil","");
    //}

    int tensorElementCount; // Number of independent elements in resulting tensor
    // Set boundary condition
    if (bc == "fixed" || bc == "f") {
        param.insertParameter("boundary_condition_type","0");
        tensorElementCount = 3; // Diagonal
    }
    else if (bc == "periodic" || bc == "p") {
        param.insertParameter("boundary_condition_type","2");
        tensorElementCount = 6; // Symmetric
    }
    else if (bc == "linear" || bc == "l") {
        param.insertParameter("boundary_condition_type","1");
        tensorElementCount = 9; // Full tensor
    }
    else {
        std::cerr << "Boundary conditions (bc) must be either fixed, periodic or linear. Value is " << bc << std::endl;
        usageandexit();
    }
    
    // Compute minimum and maximum (large scale) saturations
    std::string rock_list = param.getDefault<std::string>("rock_list", "no_list");
    std::vector<std::string> rockfiles;
    std::vector<std::vector<double> > rocksatendpoints_ = getExtremeSats(rock_list,rockfiles);

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
    Opm::SinglePhaseUpscaler spupscaler; // needed to access porosities and cell volumes
    spupscaler.init(eclparser, Opm::SinglePhaseUpscaler::Fixed,
                    0.0,0.0, linsolver_tolerance, linsolver_verbosity, linsolver_type, false, linsolver_maxit,
                    linsolver_prolongate_factor, linsolver_smooth_steps);
    std::vector<double>  cellPoreVolumes; 
    cellPoreVolumes.resize(satnums.size(), 0.0);
    double swirvolume = 0.0;
    double sworvolume = 0.0;
    // cell_idx is the eclipse index.
    const std::vector<int>& ecl_idx = spupscaler.grid().globalCell();
    Dune::CpGrid::Codim<0>::LeafIterator c = spupscaler.grid().leafbegin<0>();
    for (; c != spupscaler.grid().leafend<0>(); ++c) {
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
    double min_sat = swirvolume/poreVolume;
    double max_sat = sworvolume/poreVolume;
    // Insert computed Swir and Swor as min_sat and max_sat in param object
    //param.insertParameter("min_sat",toString(min_sat));
    //param.insertParameter("max_sat",toString(max_sat));
    
    std::vector<double> saturations;
    Opm::SparseTable<double> all_pdrops;
    // Linear range of saturations
    saturations.resize(num_sats);
    for (int i = 0; i < num_sats; ++i) {
        double factor = num_sats == 1 ? 0 : double(i)/double(num_sats - 1);
        saturations[i] = (1.0 - factor)*min_sat + factor*max_sat;
    }
    // Logarithmic range of pressure drops
    std::vector<double> pdrops;
    pdrops.resize(num_pdrops);
    for (int i = 0; i < num_pdrops; ++i) {
        double factor = num_pdrops == 1 ? 0 : double(i)/double(num_pdrops - 1);
        pdrops[i] = std::exp((1.0 - factor)*log_min_pdrop + factor*log_max_pdrop);
    }
    // Assign the same pressure drops to all saturations.
    for (int i = 0; i < num_sats; ++i) {
        all_pdrops.appendRow(pdrops.begin(), pdrops.end());
    }
    typedef SteadyStateUpscalerImplicit<UpscalingTraitsBasicImplicit> Upscaler;
    typedef Upscaler::permtensor_t permtensor_t;
    Upscaler upscaler;
    upscaler.init(param);
    
    // Compute single phase permeability
    permtensor_t upscaled_K = upscaler.upscaleSinglePhase();
    permtensor_t singlephaseperm = upscaled_K;
    singlephaseperm *= (1.0/(Opm::prefix::milli*Opm::unit::darcy));
    double porosity = upscaler.upscalePorosity();
    int num_cells = upscaler.grid().size(0);
    
    // Holders for upscaled relpermvalues for each saturation-pressure-point
    std::vector<std::vector<std::vector<double> > > RelPermPhase1, RelPermPhase2;
    for (int satidx=0; satidx<num_sats; ++satidx) {
        RelPermPhase1.push_back(std::vector<std::vector<double> >(all_pdrops[satidx].size(),std::vector<double>(tensorElementCount,0.0)));
        RelPermPhase2.push_back(std::vector<std::vector<double> >(all_pdrops[satidx].size(),std::vector<double>(tensorElementCount,0.0)));
    }
    
    for (int satidx = 0; satidx < num_sats; ++satidx) {
        std::vector<double> init_sat;
        if (start_from_cl) {
            try {
                upscaler.setToCapillaryLimit(saturations[satidx], init_sat);
            }catch (...){
                init_sat.resize(num_cells, saturations[satidx]);
                std::cout << "Failed to initialize with capillary limit for s = " << saturations[satidx]
                          << ". Init with uniform distribution." << std::endl;
            }
        } {
            init_sat.resize(num_cells, saturations[satidx]);
        }
        const Opm::SparseTable<double>::row_type pdrops = all_pdrops[satidx];
        int num_pdrops = pdrops.size();
        for (int pidx = 0; pidx < num_pdrops; ++pidx) {
            upscaler.initSatLimits(init_sat);
            double pdrop = pdrops[pidx];
            bool success = false;
            std::pair<permtensor_t, permtensor_t> lambda
                = upscaler.upscaleSteadyState(fldir, init_sat, saturations[satidx], pdrop, upscaled_K, success);
            double usat = upscaler.lastSaturationUpscaled();
            std::cout << "metn: " << saturations[satidx] << " " << usat << std::endl;
            if(! success){
                std::cout << "Upscaling failed for " << usat << std::endl;
            }else{
                init_sat = upscaler.lastSaturationState();
            }
            // Store upscaled values
            for (int voigtIdx=0; voigtIdx<tensorElementCount; ++voigtIdx) {
                RelPermPhase1[satidx][pidx][voigtIdx] = getVoigtValue(lambda.first,voigtIdx);
                RelPermPhase2[satidx][pidx][voigtIdx] = getVoigtValue(lambda.second,voigtIdx);
            }
        }
    }

    // Output results to stdout and optionally to file.
    std::stringstream outputtmp;
    outputtmp << "######################################################################" << std::endl;
    outputtmp << "Results from steady state upscaling of relative permeability (impes)" << std::endl;
    outputtmp << "#" << std::endl;
    time_t now = std::time(NULL);
    outputtmp << "# Finished: " << asctime(localtime(&now));
    
    utsname hostname;   uname(&hostname);
    outputtmp << "# Hostname: " << hostname.nodename << std::endl;
    outputtmp << "#" << std::endl;
    outputtmp << "# Eclipse file: " << gridfilename << std::endl;
    outputtmp << "#        cells: " << num_cells << std::endl;
    outputtmp << "#  Pore volume: " << poreVolume << std::endl;
    outputtmp << "#     Porosity: " << porosity << std::endl;
    outputtmp << "#" << std::endl;
    outputtmp << "# Rock list: " << rock_list << std::endl;
    outputtmp << "#   with the following rock files: " << std::endl;    
    for (std::size_t ridx=0; ridx<rockfiles.size(); ++ridx) {
        outputtmp << "#   " << rockfiles[ridx] << std::endl;
    }
    outputtmp << "#" << std::endl;

    outputtmp << "# Options used:" << std::endl;
    outputtmp << "#     Boundary conditions: " << bc << std::endl;
    outputtmp << "#          flow direction: " << flowdir << std::endl;    
    outputtmp << "#       saturation points: " << num_sats << std::endl;
    outputtmp << "#    pressure drop points: " << num_pdrops << std::endl;    
    //outputtmp << "#         maxPermContrast: " << options["maxPermContrast"] << std::endl;
    //outputtmp << "#                 minPerm: " << options["minPerm"] << std::endl;
    //outputtmp << "#                 minPoro: " << options["minPoro"] << std::endl;
    outputtmp << "#          surfaceTension: " << surfaceTension << " dynes/cm" << std::endl;
    //if (includeGravity) {
    //    outputtmp << "#                 gravity: " << options["gravity"] << " m/s²" << endl;
    //    // Output values of water/gas and oil densities
    //}
    outputtmp << "#" << std::endl;
    outputtmp << "# Single phase permeability" << std::endl; 
    outputtmp << "#  |Kxx  Kxy  Kxz| = " << singlephaseperm(0,0) << "  " << singlephaseperm(0,1) << "  " << singlephaseperm(0,2) << std::endl; 
    outputtmp << "#  |Kyx  Kyy  Kyz| = " << singlephaseperm(1,0) << "  " << singlephaseperm(1,1) << "  " << singlephaseperm(1,2) << std::endl; 
    outputtmp << "#  |Kzx  Kzy  Kzz| = " << singlephaseperm(2,0) << "  " << singlephaseperm(2,1) << "  " << singlephaseperm(2,2) << std::endl; 
    outputtmp << "# " << std::endl;
    
    std::cout << outputtmp.str();
    
    // MPIHelper::instance(argc,argv) ;
    //typedef SteadyStateUpscalerImplicit<UpscalingTraitsBasicImplicit> upscaler_t;
    //SteadyStateUpscalerManagerImplicit<upscaler_t> mgr;
    //mgr.upscale(param);
    
}

// TODO
// * Include gravity?
// * Include possibility for maxPermContrast, minPerm, minPoro
// * Include possibility for interpolating relpermcurves (as functions of saturation)
