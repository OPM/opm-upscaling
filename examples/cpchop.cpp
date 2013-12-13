/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2010 Statoil ASA.

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
#include <config.h>

#include <opm/core/io/eclipse/CornerpointChopper.hpp>
#include <opm/core/io/eclipse/EclipseGridInspector.hpp>
#include <opm/upscaling/SinglePhaseUpscaler.hpp>
#include <opm/core/utility/MonotCubicInterpolator.hpp>
#include <opm/porsol/common/setupBoundaryConditions.hpp>
#include <opm/core/utility/Units.hpp>

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
#include <cfloat>
#include <cmath>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
#include <dune/common/parallel/mpihelper.hh>
#else
#include <dune/common/mpihelper.hh>
#endif

// upscaling of endpoints and capillary pressure
// Conversion factor, multiply mD numbers with this to get m² numbers
const double milliDarcyToSqMetre = 9.869233e-16;

const double saturationThreshold = 0.00001;

typedef struct {
  int imin;
  int imax;
  int jmin;
  int jmax;
  double zmin;
  double zmax;
  int ilen;
  int jlen;
  double zlen;
  bool upscale;
  std::string bc;
  Opm::SinglePhaseUpscaler::BoundaryConditionType bctype;
  bool resettoorigin;
  boost::mt19937::result_type userseed;
  std::string filebase;
  double minperm;
  double minpermSI;
  int subsamples;
  bool dips;  // whether to do dip averaging
  double azimuthdisplacement;  // posibility to add/subtract a value to/from azimuth for dip plane.
  double mincellvolume; // ignore smaller cells for dip calculations
  bool satnumvolumes; // whether to count volumes pr. satnum
  double surfaceTension; // multiply with 10^-3 to obtain SI units 
  bool endpoints; // whether to upscale saturation endpoints
  bool cappres; // whether to upscale capillary pressure
  std::string rock_list;
  bool anisorocks;
  double z_tolerance;
  double residual_tolerance;
  int linsolver_verbosity;
  int linsolver_type;
  std::vector<std::vector<double> > rocksatendpoints_;
  std::vector<std::vector<double> > jfuncendpoints_; // Used if isotropic rock input
  bool isFixed;
  bool isPeriodic;
  std::vector<Opm::MonotCubicInterpolator> InvJfunctions; // Holds the inverse of the loaded J-functions.
  // For anisotropic input rocks:
  std::vector<Opm::MonotCubicInterpolator> SwPcfunctions; // Holds Sw(Pc) for each rocktype.
  int nsatpoints;
} ChopSettings;

typedef struct {
    std::vector<double> porosities;
    std::vector<double> permxs;
    std::vector<double> permys;
    std::vector<double> permzs;
    std::vector<double> permyzs;
    std::vector<double> permxzs;
    std::vector<double> permxys;
    std::vector<double> minsws;
    std::vector<double> maxsws;
    std::vector< std::vector<double> > pcvalues;
    std::vector<double> dipangs;
    std::vector<double> azimuths;
    // Initialize a matrix for subsample satnum volumes. 
    // Outer index is subsample index, inner index is SATNUM-value
    std::vector< std::vector<double> > rockvolumes;
    std::vector<double> netporosities;
    std::vector<double> ntgs;
    std::vector<int> subsampletab;
    std::vector<int> subsample_failed;

    void resize(int size)
    {
      porosities.resize(size);
      permxs.resize(size);
      permys.resize(size);
      permzs.resize(size);
      permyzs.resize(size);
      permxzs.resize(size);
      permxys.resize(size);
      netporosities.resize(size);
      ntgs.resize(size);
      minsws.resize(size);
      maxsws.resize(size);
      rockvolumes.resize(size);
      pcvalues.resize(size);
      dipangs.resize(size);
      azimuths.resize(size);
    }
} ChopThreadContext;

void do_chop(int sample, const ChopSettings& settings,
             ChopThreadContext& tcontext,
             int ri, int rj, double rz, const Opm::CornerPointChopper& ch)
{
  Opm::CornerPointChopper::ChopContext context;
  ch.chop(ri, ri + settings.ilen, rj, rj + settings.jlen,
          rz , rz + settings.zlen, context, settings.resettoorigin);
  std::string subsampledgrdecl = settings.filebase;

  // Output grdecl-data to file if a filebase is supplied.
  if (settings.filebase != "") {
    std::ostringstream oss;
    if ((size_t) settings.subsamples > 1) { // Only add number to filename if more than one sample is asked for
      oss << 'R' << std::setw(4) << std::setfill('0') << sample;
      subsampledgrdecl += oss.str();
    }
    subsampledgrdecl += ".grdecl";
    ch.writeGrdecl(subsampledgrdecl, context);
  }


  //  Guarantee initialization
  double Pcmax = -DBL_MAX, Pcmin = DBL_MAX;

  try { /* The upscaling may fail to converge on icky grids, lets just pass by those */
    if (settings.upscale) {
      Opm::EclipseGridParser subparser = ch.subparser(context);
      subparser.convertToSI();
      Opm::SinglePhaseUpscaler upscaler;

      upscaler.init(subparser, settings.bctype, settings.minpermSI,
                    settings.z_tolerance, settings.residual_tolerance,
                    settings.linsolver_verbosity, settings.linsolver_type, false);

      Opm::SinglePhaseUpscaler::permtensor_t upscaled_K = upscaler.upscaleSinglePhase();
      upscaled_K *= (1.0/(Opm::prefix::milli*Opm::unit::darcy));

      tcontext.porosities.push_back(upscaler.upscalePorosity());
      if (!context.new_NTG_.empty()) {
        tcontext.netporosities.push_back(upscaler.upscaleNetPorosity());
        tcontext.ntgs.push_back(upscaler.upscaleNTG());
      }
      tcontext.permxs.push_back(upscaled_K(0,0));
      tcontext.permys.push_back(upscaled_K(1,1));
      tcontext.permzs.push_back(upscaled_K(2,2));
      tcontext.permyzs.push_back(upscaled_K(1,2));
      tcontext.permxzs.push_back(upscaled_K(0,2));
      tcontext.permxys.push_back(upscaled_K(0,1));
    }

    if (settings.endpoints) {
      // Calculate minimum and maximum water volume in each cell
      // Create single-phase upscaling object to get poro and perm values from the grid
      Opm::EclipseGridParser subparser = ch.subparser(context);
      std::vector<double>  perms = subparser.getFloatingPointValue("PERMX");
      subparser.convertToSI();
      Opm::SinglePhaseUpscaler upscaler;                
      upscaler.init(subparser, settings.bctype, settings.minpermSI,
                    settings.z_tolerance, settings.residual_tolerance,
                    settings.linsolver_verbosity, settings.linsolver_type, false);
      std::vector<int>   satnums = subparser.getIntegerValue("SATNUM");
      std::vector<double>  poros = subparser.getFloatingPointValue("PORO");
      std::vector<double> cellVolumes, cellPoreVolumes;
      cellVolumes.resize(satnums.size(), 0.0);
      cellPoreVolumes.resize(satnums.size(), 0.0);
      int tesselatedCells = 0;
      //double maxSinglePhasePerm = 0;
      double Swirvolume = 0;
      double Sworvolume = 0;
      const std::vector<int>& ecl_idx = upscaler.grid().globalCell();
      Dune::CpGrid::Codim<0>::LeafIterator c = upscaler.grid().leafbegin<0>();
      for (; c != upscaler.grid().leafend<0>(); ++c) {
        unsigned int cell_idx = ecl_idx[c->index()];
        if (satnums[cell_idx] > 0) { // Satnum zero is "no rock"
          cellVolumes[cell_idx] = c->geometry().volume();
          cellPoreVolumes[cell_idx] = cellVolumes[cell_idx] * poros[cell_idx];
          double Pcmincandidate = 0.0, Pcmaxcandidate = 0.0, minSw, maxSw;
          if (!settings.anisorocks) {
            if (settings.cappres) {
              Pcmincandidate = settings.jfuncendpoints_[int(satnums[cell_idx])-1][1]
                / sqrt(perms[cell_idx] * milliDarcyToSqMetre/poros[cell_idx]) * settings.surfaceTension;
              Pcmaxcandidate = settings.jfuncendpoints_[int(satnums[cell_idx])-1][0]
                / sqrt(perms[cell_idx] * milliDarcyToSqMetre/poros[cell_idx]) * settings.surfaceTension;
            }
            minSw = settings.rocksatendpoints_[int(satnums[cell_idx])-1][0];
            maxSw = settings.rocksatendpoints_[int(satnums[cell_idx])-1][1];
          }
          else { // anisotropic input, we do not to J-function scaling
            if (settings.cappres) {
              Pcmincandidate = settings.SwPcfunctions[int(satnums[cell_idx])-1].getMinimumX().first;
              Pcmaxcandidate = settings.SwPcfunctions[int(satnums[cell_idx])-1].getMaximumX().first;
            }
            minSw = settings.rocksatendpoints_[int(satnums[cell_idx])-1][0];
            maxSw = settings.rocksatendpoints_[int(satnums[cell_idx])-1][1];
          }
          if (settings.cappres) {
            Pcmin = std::min(Pcmincandidate, Pcmin);
            Pcmax = std::max(Pcmaxcandidate, Pcmax);
          }
          Swirvolume += minSw * cellPoreVolumes[cell_idx];
          Sworvolume += maxSw * cellPoreVolumes[cell_idx];
        }
        ++tesselatedCells; // keep count.
      }

      // If upscling=false, we still (may) want to have porosities together with endpoints
      if (!settings.upscale) {
        tcontext.porosities.push_back(upscaler.upscalePorosity());
      }

      // Total porevolume and total volume -> upscaled porosity:
      double poreVolume = std::accumulate(cellPoreVolumes.begin(),
                                          cellPoreVolumes.end(), 0.0);
      double Swir = Swirvolume/poreVolume;
      double Swor = Sworvolume/poreVolume;
      tcontext.minsws.push_back(Swir);
      tcontext.maxsws.push_back(Swor);
      if (settings.cappres) {
        // Upscale capillary pressure function
        Opm::MonotCubicInterpolator WaterSaturationVsCapPressure;
        double largestSaturationInterval = Swor-Swir;
        double Ptestvalue = Pcmax;
        while (largestSaturationInterval > (Swor-Swir)/double(settings.nsatpoints)) {
          if (Pcmax == Pcmin) {
            // This is a dummy situation, we go through once and then 
            // we are finished (this will be triggered by zero permeability)
            Ptestvalue = Pcmin;
            largestSaturationInterval = 0;
          }
          else if (WaterSaturationVsCapPressure.getSize() == 0) {
            /* No data values previously computed */
            Ptestvalue = Pcmax;
          }
          else if (WaterSaturationVsCapPressure.getSize() == 1) {
            /* If only one point has been computed, it was for Pcmax. So now
               do Pcmin */
            Ptestvalue = Pcmin;
          }
          else {
            /* Search for largest saturation interval in which there are no
               computed saturation points (and estimate the capillary pressure
               that will fall in the center of this saturation interval)
               */
            std::pair<double,double> SatDiff = WaterSaturationVsCapPressure.getMissingX();
            Ptestvalue = SatDiff.first;
            largestSaturationInterval = SatDiff.second;
          }
          // Check for saneness of Ptestvalue:
          if (std::isnan(Ptestvalue) || std::isinf(Ptestvalue)) {
            std::cerr << "ERROR: Ptestvalue was inf or nan" << std::endl;
            break; // Jump out of while-loop, just print out the results
            // up to now and exit the program
          }

          double waterVolume = 0.0;
          for (unsigned int i = 0; i < ecl_idx.size(); ++i) {
            unsigned int cell_idx = ecl_idx[i];
            double waterSaturationCell = 0.0;
            if (satnums[cell_idx] > 0) { // handle "no rock" cells with satnum zero
              double PtestvalueCell;

              PtestvalueCell = Ptestvalue;

              if (!settings.anisorocks) {   
                double Jvalue = sqrt(perms[cell_idx] * milliDarcyToSqMetre /poros[cell_idx]) * PtestvalueCell / settings.surfaceTension;
                waterSaturationCell 
                  = settings.InvJfunctions[int(satnums[cell_idx])-1].evaluate(Jvalue);
              }
              else { // anisotropic_input, then we do not do J-function-scaling
                waterSaturationCell = settings.SwPcfunctions[int(satnums[cell_idx])-1].evaluate(PtestvalueCell);
              }
            }
            waterVolume += waterSaturationCell  * cellPoreVolumes[cell_idx];
          }
          WaterSaturationVsCapPressure.addPair(Ptestvalue, waterVolume/poreVolume);
        }
        WaterSaturationVsCapPressure.chopFlatEndpoints(saturationThreshold);
        std::vector<double> wattest = WaterSaturationVsCapPressure.get_fVector();
        std::vector<double> cprtest = WaterSaturationVsCapPressure.get_xVector();
        Opm::MonotCubicInterpolator CapPressureVsWaterSaturation(WaterSaturationVsCapPressure.get_fVector(), 
            WaterSaturationVsCapPressure.get_xVector());
        std::vector<double> pcs;
        for (int satp=0; satp<settings.nsatpoints; ++satp) {
          pcs.push_back(CapPressureVsWaterSaturation.evaluate(Swir+(Swor-Swir)/(settings.nsatpoints-1)*satp));
        }
        tcontext.pcvalues.push_back(pcs);
      }

    }


    if (settings.dips) {
      Opm::EclipseGridParser subparser = ch.subparser(context);
      std::vector<int>  griddims = subparser.getSPECGRID().dimensions;
      std::vector<double> xdips_subsample, ydips_subsample;

      Opm::EclipseGridInspector gridinspector(subparser);
      for (int k=0; k < griddims[2]; ++k) {
        for (int j=0; j < griddims[1]; ++j) {
          for (int i=0; i < griddims[0]; ++i) {
            if (gridinspector.cellVolumeVerticalPillars(i, j, k) > settings.mincellvolume) {
              std::pair<double,double> xydip = gridinspector.cellDips(i, j, k);                       
              xdips_subsample.push_back(xydip.first);
              ydips_subsample.push_back(xydip.second);
            }
          }
        }
      }


      //  double azimuth = atan(xydip.first/xydip.second);
      //              double dip = acos(1.0/sqrt(pow(xydip.first,2.0)+pow(xydip.second,2.0)+1.0));
      //	dips_subsample.push_back( xydip.first );
      //	azims_subsample.push_back(atan(xydip.first/xydip.second));         

      // Average xdips and ydips
      double xdipaverage = accumulate(xdips_subsample.begin(), xdips_subsample.end(), 0.0)/xdips_subsample.size();
      double ydipaverage = accumulate(ydips_subsample.begin(), ydips_subsample.end(), 0.0)/ydips_subsample.size();

      // Convert to dip and azimuth
      double azimuth = atan(xdipaverage/ydipaverage)+settings.azimuthdisplacement;
      double dip = acos(1.0/sqrt(pow(xdipaverage,2.0)+pow(ydipaverage,2.0)+1.0));
      tcontext.dipangs.push_back(dip);
      tcontext.azimuths.push_back(azimuth);                    
    }

    int maxSatnum = 0; // This value is determined from the chopped cells.
    if (settings.satnumvolumes) {
      Opm::EclipseGridParser subparser = ch.subparser(context);
      Opm::EclipseGridInspector subinspector(subparser);
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
      tcontext.rockvolumes.push_back(rockvolumessubsample_normalized);
    }
  }
  catch (...) {
    std::cerr << "Warning: Upscaling chopped subsample nr. " << sample << " failed, proceeding to next subsample\n";
  }
}

int main(int argc, char** argv)
try
{
    if (argc == 1) {
        std::cout << "Usage: cpchop gridfilename=filename.grdecl [subsamples=10] [ilen=5] [jlen=5] " << std::endl;
        std::cout << "       [zlen=5] [imin=] [imax=] [jmin=] [jmax=] [upscale=true] [bc=fixed]" << std::endl;
        std::cout << "       [resettoorigin=true] [seed=111] [z_tolerance=0.0] [minperm=1e-9] " << std::endl;
        std::cout << "       [dips=false] [azimuthdisplacement=] [satnumvolumes=false] [mincellvolume=1e-9]" << std::endl;
        std::cout << "       [filebase=] [resultfile=] [endpoints=false] [cappres=false]" << std::endl;
        std::cout << "       [rock_list=] [anisotropicrocks=false]" << std::endl;
        exit(1);
    }

    Dune::MPIHelper::instance(argc, argv);

    Opm::parameter::ParameterGroup param(argc, argv);
    std::string gridfilename = param.get<std::string>("gridfilename");
    Opm::CornerPointChopper ch(gridfilename);

    // The cells with i coordinate in [imin, imax) are included, similar for j.
    // The z limits may be changed inside the chopper to match actual min/max z.
    const int* dims = ch.dimensions();
    ChopSettings settings;
    settings.imin = param.getDefault("imin", 0);
    settings.imax = param.getDefault("imax", dims[0]);
    settings.jmin = param.getDefault("jmin", 0);
    settings.jmax = param.getDefault("jmax", dims[1]);
    settings.zmin = param.getDefault("zmin", ch.zLimits().first);
    settings.zmax = param.getDefault("zmax", ch.zLimits().second);
    settings.subsamples = param.getDefault("subsamples", 1);
    settings.ilen = param.getDefault("ilen", settings.imax - settings.imin);
    settings.jlen = param.getDefault("jlen", settings.jmax - settings.jmin);
    settings.zlen = param.getDefault("zlen", settings.zmax - settings.zmin);
    settings.upscale = param.getDefault("upscale", true);
    settings.bc = param.getDefault<std::string>("bc", "fixed");
    settings.resettoorigin = param.getDefault("resettoorigin", true);
    settings.userseed = param.getDefault("seed", 0);

    int outputprecision = param.getDefault("outputprecision", 8);
    settings.filebase = param.getDefault<std::string>("filebase", "");
    std::string resultfile = param.getDefault<std::string>("resultfile", "");

    settings.minperm = param.getDefault("minperm", 1e-9);
    settings.minpermSI = Opm::unit::convert::from(settings.minperm, Opm::prefix::milli*Opm::unit::darcy);

    // Following two options are for dip upscaling (slope of cell top and bottom edges)
    settings.dips = param.getDefault("dips", false);  // whether to do dip averaging
    settings.azimuthdisplacement = param.getDefault("azimuthdisplacement", 0.0);  // posibility to add/subtract a value to/from azimuth for dip plane.
    settings.mincellvolume = param.getDefault("mincellvolume", 1e-9); // ignore smaller cells for dip calculations

    settings.satnumvolumes = param.getDefault("satnumvolumes", false); // whether to count volumes pr. satnum

    // Input for surfaceTension is dynes/cm, SI units are Joules/square metre
    settings.surfaceTension = param.getDefault("surfaceTension", 11.0) * 1e-3; // multiply with 10^-3 to obtain SI units 

    settings.endpoints = param.getDefault("endpoints", false); // whether to upscale saturation endpoints
    settings.cappres = param.getDefault("cappres", false); // whether to upscale capillary pressure
    if (settings.cappres) { settings.endpoints = true; }
    settings.rock_list = param.getDefault<std::string>("rock_list", "no_list");
    settings.anisorocks = param.getDefault("anisotropicrocks", false);
    settings.nsatpoints = 5; // number of saturation points in upscaled capillary pressure function per subsample

    // Read rock data from files specifyed in rock_list
    if (settings.endpoints) {
        if (!settings.rock_list.compare("no_list")) {
            std::cout << "Can't do endponts without rock list (" << settings.rock_list << ")" << std::endl;
            throw std::exception();
        }
        // Code copied from ReservoirPropertyCommon.hpp for file reading
        std::ifstream rl(settings.rock_list.c_str());
        if (!rl) {
            OPM_THROW(std::runtime_error, "Could not open file " << settings.rock_list);
        }
        int num_rocks = -1;
        rl >> num_rocks;
        assert(num_rocks >= 1);
        settings.rocksatendpoints_.resize(num_rocks, std::vector<double>(2, 0.0));
        settings.jfuncendpoints_.resize(num_rocks, std::vector<double>(2, 0.0));
        // Loop through rock files defined in rock_list and store the data we need
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
                OPM_THROW(std::runtime_error, "Could not open file " << rockfilename);
            }
            
            if (! settings.anisorocks) { //Isotropic input rocks (Sw Krw Kro J)
                Opm::MonotCubicInterpolator Jtmp;
                try {
                    Jtmp = Opm::MonotCubicInterpolator(rockname, 1, 4); 
                }
                catch (const char * errormessage) {
                    std::cerr << "Error: " << errormessage << std::endl;
                    std::cerr << "Check filename" << std::endl;
                    exit(1);
                }
                
                // Invert J-function, now we get saturation as a function of pressure:
                if (Jtmp.isStrictlyMonotone()) {
                    settings.InvJfunctions.push_back(Opm::MonotCubicInterpolator(Jtmp.get_fVector(), Jtmp.get_xVector()));
                }
                else {
                    std::cerr << "Error: Jfunction " << i+1 << " in rock file " << rockname << " was not invertible." << std::endl;
                    exit(1);
                }
                
                settings.jfuncendpoints_[i][0] = Jtmp.getMinimumX().second;
                settings.jfuncendpoints_[i][1] = Jtmp.getMaximumX().second;
                settings.rocksatendpoints_[i][0] = Jtmp.getMinimumX().first;
                settings.rocksatendpoints_[i][1] = Jtmp.getMaximumX().first;
                if (settings.rocksatendpoints_[i][0] < 0 || settings.rocksatendpoints_[i][0] > 1) {
                    OPM_THROW(std::runtime_error, "Minimum rock saturation (" << settings.rocksatendpoints_[i][0] << ") not sane for rock " 
                          << rockfilename << "." << std::endl << "Did you forget to specify anisotropicrocks=true ?");  
                }
            }
            else { //Anisotropic input rocks (Pc Sw Krxx Kryy Krzz)
                Opm::MonotCubicInterpolator Pctmp;
                try {
                    Pctmp = Opm::MonotCubicInterpolator(rockname, 2, 1);
                }
                catch (const char * errormessage) {
                    std::cerr << "Error: " << errormessage << std::endl;
                    std::cerr << "Check filename and columns 1 and 2 (Pc and Sw)" << std::endl;
                    exit(1);
                }
                if (settings.cappres) {
                    // Invert Pc(Sw) curve into Sw(Pc):
                    if (Pctmp.isStrictlyMonotone()) {
                        settings.SwPcfunctions.push_back(Opm::MonotCubicInterpolator(Pctmp.get_fVector(), Pctmp.get_xVector()));
                    }
                    else {
                        std::cerr << "Error: Pc(Sw) curve " << i+1 << " in rock file " << rockname << " was not invertible." << std::endl;
                        exit(1);
                    }
                }
                settings.rocksatendpoints_[i][0] = Pctmp.getMinimumX().first;
                settings.rocksatendpoints_[i][1] = Pctmp.getMaximumX().first;
            }          
        }
    }

    settings.z_tolerance = param.getDefault("z_tolerance", 0.0);
    settings.residual_tolerance = param.getDefault("residual_tolerance", 1e-8);
    settings.linsolver_verbosity = param.getDefault("linsolver_verbosity", 0);
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 3) || defined(HAS_DUNE_FAST_AMG)
    settings.linsolver_type = param.getDefault("linsolver_type", 3);
#else
    settings.linsolver_type = param.getDefault("linsolver_type", 1);
#endif

    // Check that we do not have any user input
    // that goes outside the coordinates described in
    // the cornerpoint file (runtime-exception will be thrown in case of error)
    ch.verifyInscribedShoebox(settings.imin, settings.ilen, settings.imax,
			      settings.jmin, settings.jlen, settings.jmax,
			      settings.zmin, settings.zlen, settings.zmax);

    // Random number generator from boost.
    boost::mt19937 gen;

    // Seed the random number generators with the current time, unless specified on command line
    // Warning: Current code does not allow 0 for the seed!!
    boost::mt19937::result_type autoseed = time(NULL);
    if (settings.userseed == 0) {
        gen.seed(autoseed);
    }
    else {
        gen.seed(settings.userseed);
    }

    settings.bctype = Opm::SinglePhaseUpscaler::Fixed;
    settings.isFixed = settings.isPeriodic = false;
    if (settings.upscale) {
        if (settings.bc == "fixed") {
            settings.isFixed = true;
            settings.bctype = Opm::SinglePhaseUpscaler::Fixed;
        }
        else if (settings.bc == "periodic") {
            settings.isPeriodic = true;
            settings.bctype = Opm::SinglePhaseUpscaler::Periodic;
        }
        else {
            std::cout << "Boundary condition type (bc=" << settings.bc << ") not allowed." << std::endl;
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
    boost::uniform_int<> disti(settings.imin, settings.imax - settings.ilen);
    boost::uniform_int<> distj(settings.jmin, settings.jmax - settings.jlen);
    boost::uniform_real<> distz(settings.zmin, std::max(settings.zmax - settings.zlen, settings.zmin));
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > ri(gen, disti);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > rj(gen, distj);
    boost::variate_generator<boost::mt19937&, boost::uniform_real<> > rz(gen, distz);

    // Storage for results
    int threads = 1;
#ifdef HAVE_OPENMP
    threads = omp_get_max_threads();
#endif
    std::vector<ChopThreadContext> ctx(threads);

    // draw the random numbers up front to ensure consistency with or without threads 
    std::vector<int> ris(settings.subsamples);
    std::vector<int> rjs(settings.subsamples);
    std::vector<double> rzs(settings.subsamples);
    for (int j=0;j<settings.subsamples;++j) {
      ris[j] = ri();
      rjs[j] = rj();
      rzs[j] = rz();
    }

    // first subsample has to run single threaded.
    // OPM builds up a lot of maps and tables as new entries are encountered
    // and this breaks badly with multiple threads.
    // run first subsample single threaded to get tables initialized,
    // then run the rest in parallel
    do_chop(1, settings, ctx[0], ris[0], rjs[0], rzs[0], ch);
    int finished_subsamples = 1; // keep explicit count of successful subsamples
    ctx[0].subsampletab.push_back(0);
    

#pragma omp parallel for schedule(static) reduction(+:finished_subsamples)
    for (int sample = 2; sample <= settings.subsamples; ++sample) {
        int thread = 0;
#ifdef HAVE_OPENMP
        thread = omp_get_thread_num();
#endif
        try {
            do_chop(sample, settings, ctx[thread], 
                    ris[sample-1], rjs[sample-1], rzs[sample-1], ch);
            ctx[thread].subsampletab.push_back(sample-1);
            finished_subsamples++;
        } catch (...) {
            ctx[thread].subsample_failed.push_back(sample-1);
            std::cerr << "Warning: Upscaling chopped subsample nr. " << sample << " failed, proceeding to next subsample\n";
        }
    }

    // compress data
    ChopThreadContext sc;
    int l=0;
    std::vector<int> mapping(settings.subsamples, -1);
    for (int k=0;k<settings.subsamples;++k) {
      bool append=true;
      for (int i=0;i<threads;++i) {
        if (std::find(ctx[i].subsample_failed.begin(),
                      ctx[i].subsample_failed.end(), k) != ctx[i].subsample_failed.end())
          append = false;
      }
      if (append)
        mapping[k] = l++;
    }
    sc.resize(l);

    for (int i=0;i<threads;++i) {
      for (size_t j=0;j<ctx[i].subsampletab.size();++j) {
        int idx = mapping[ctx[i].subsampletab[j]];
        if (idx == -1)
          continue;
        if (settings.upscale) {
          sc.porosities[idx] = ctx[i].porosities[j];
          sc.permxs[idx] = ctx[i].permxs[j];
          sc.permys[idx] = ctx[i].permys[j];
          sc.permzs[idx] = ctx[i].permzs[j];
          sc.permyzs[idx] = ctx[i].permyzs[j];
          sc.permxzs[idx] = ctx[i].permxzs[j];
          sc.permxys[idx] = ctx[i].permxys[j];
          if (ctx[i].netporosities.size()) {
            sc.netporosities[idx] = ctx[i].netporosities[j];
            sc.ntgs[idx] = ctx[i].ntgs[j];
          }
        }
        if (settings.endpoints) {
          sc.minsws[idx] = ctx[i].minsws[j];
          sc.maxsws[idx] = ctx[i].maxsws[j];
        }
        if (settings.cappres)
          sc.pcvalues[idx] = ctx[i].pcvalues[j];
        if (settings.dips)
          sc.azimuths[idx] = ctx[i].azimuths[j];
        if (settings.satnumvolumes)
          sc.rockvolumes[idx] = ctx[i].rockvolumes[j];
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
    outputtmp << "#   i; min,len,max: " << settings.imin << " " << settings.ilen << " " << settings.imax << std::endl;
    outputtmp << "#   j; min,len,max: " << settings.jmin << " " << settings.jlen << " " << settings.jmax << std::endl;
    outputtmp << "#   z; min,len,max: " << settings.zmin << " " << settings.zlen << " " << settings.zmax << std::endl;
    outputtmp << "#       subsamples: " << settings.subsamples << std::endl;
    if (settings.userseed == 0) {
        outputtmp << "#      (auto) seed: " << autoseed << std::endl;
    }
    else {
        outputtmp << "#    (manual) seed: " << settings.userseed << std::endl; 
    }        
    outputtmp << "################################################################################################" << std::endl;
    outputtmp << "# id";
    if (settings.upscale) {
        if (settings.isPeriodic) {
            if (!ctx[0].ntgs.empty()) {
                outputtmp << "          porosity                netporosity             ntg                      permx                   permy                   permz                   permyz                  permxz                  permxy                  netpermh";
            }
            else {
                outputtmp << "          porosity                 permx                   permy                   permz                   permyz                  permxz                  permxy";
            }
        }
        else if (settings.isFixed) {
            if (!ctx[0].ntgs.empty()) {
                outputtmp << "          porosity                netporosity             ntg                      permx                   permy                   permz                   netpermh";
            }
            else {
                outputtmp << "          porosity                 permx                   permy                   permz";
            }
        }
    }
    if (settings.endpoints) {
        if (!settings.upscale) {
            outputtmp << "          porosity";
        }
        outputtmp << "                  Swir                    Swor";
        if (settings.cappres) {
            outputtmp << "                  Pc(Swir)                Pc2                     Pc3                     Pc4                     Pc(Swor)";            
        }
    }
    if (settings.dips) {
        outputtmp << "                  dip                     azim(displacement:" << settings.azimuthdisplacement << ")";
    }
    if (settings.satnumvolumes) {
	for (size_t satnumidx = 0; satnumidx < sc.rockvolumes[0].size(); ++satnumidx) {
	    outputtmp << "               satnum_" << satnumidx+1;
	}
    }
    outputtmp << std::endl;

    const int fieldwidth = outputprecision + 8;
    for (int sample = 1; sample <= finished_subsamples; ++sample) {
	outputtmp << sample << '\t';
	if (settings.upscale) {
	    outputtmp <<
		std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << sc.porosities[sample-1] << '\t';
            if (!ctx[0].ntgs.empty()) {
                outputtmp << std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << sc.netporosities[sample-1] << '\t' <<
                    std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << sc.ntgs[sample-1] << '\t';
            }
            outputtmp <<
		std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << sc.permxs[sample-1] << '\t' <<
		std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << sc.permys[sample-1] << '\t' <<
		std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << sc.permzs[sample-1] << '\t';
            if (settings.isPeriodic) {
                outputtmp <<
                    std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << sc.permyzs[sample-1] << '\t' <<
                    std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << sc.permxzs[sample-1] << '\t' <<
                    std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << sc.permxys[sample-1] << '\t';                
            }
            if (!ctx[0].ntgs.empty()) {
                outputtmp <<
                    std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << (sc.permxs[sample-1]+sc.permys[sample-1])/(2.0*sc.ntgs[sample-1]) << '\t';
            }
	}
	if (settings.endpoints) {
            if (!settings.upscale) {
                outputtmp <<
                    std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << sc.porosities[sample-1] << '\t';
                if (!ctx[0].ntgs.empty()) {
                    outputtmp <<
                        std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << sc.netporosities[sample-1] << '\t' <<
                        std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << sc.ntgs[sample-1] << '\t';
                }

            }
	    outputtmp <<
		std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << sc.minsws[sample-1] << '\t' <<
		std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << sc.maxsws[sample-1];
            if (settings.cappres) {
                outputtmp <<
                    std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << sc.pcvalues[sample-1][0] << '\t' <<
                    std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << sc.pcvalues[sample-1][1] << '\t' <<
                    std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << sc.pcvalues[sample-1][2] << '\t' <<
                    std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << sc.pcvalues[sample-1][3] << '\t' <<
                    std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << sc.pcvalues[sample-1][4];
            }
	}
	if (settings.dips) {
	    outputtmp <<
		std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << sc.dipangs[sample-1] << '\t' <<
		std::showpoint << std::setw(fieldwidth) << std::setprecision(outputprecision) << sc.azimuths[sample-1];
	}
	if (settings.satnumvolumes) {
	    for (size_t satnumidx = 0; satnumidx < sc.rockvolumes[sample-1].size(); ++satnumidx) {
		outputtmp <<
		    std::showpoint << std::setw(fieldwidth) <<
		    std::setprecision(outputprecision) << sc.rockvolumes[sample-1][satnumidx] << '\t';
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
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}
