/*
  Copyright 2011 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_BLACKOILSIMULATOR_HEADER_INCLUDED
#define OPM_BLACKOILSIMULATOR_HEADER_INCLUDED




#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/common/Units.hpp>
#include <dune/common/EclipseGridParser.hpp>
#include <dune/common/param/ParameterGroup.hpp>
#include <dune/porsol/common/BoundaryConditions.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <numeric>


namespace Opm
{

    template<class Grid, class Rock, class Fluid, class Wells, class FlowSolver, class TransportSolver>
    class BlackoilSimulator
    {
    public:
        void init(const Dune::parameter::ParameterGroup& param);
        void simulate();

    private:
        Grid grid_;
        Rock rock_;
        Fluid fluid_;
        Wells wells_;
        FlowSolver flow_solver_;
        TransportSolver transport_solver_;

        Dune::BasicBoundaryConditions<true, false> flow_bc_;
        typename Grid::Vector gravity_;
        std::vector<double> src_;

        typedef typename Fluid::CompVec CompVec;
        typedef typename Fluid::PhaseVec PhaseVec;
        std::vector<PhaseVec> cell_pressure_;
        std::vector<PhaseVec> face_pressure_;
        std::vector<double> well_perf_pressure_;
        std::vector<double> well_perf_flux_;
        std::vector<CompVec> cell_z_;
        PhaseVec bdy_pressure_;
        CompVec bdy_z_;

        double total_time_;
        double initial_stepsize_;
        bool increase_stepsize_;
        double stepsize_increase_factor_;
        double maximum_stepsize_;
        std::vector<double> report_times_;
        bool do_impes_;
        std::string output_dir_;
        
        bool computeCellState(int iCell, 
                              int iRef,
                              double wo_contact_depth,
                              double go_contact_depth, 
                              double connate_water_saturation);


        static void output(const Grid& grid,
                           const std::vector<typename Fluid::PhaseVec>& cell_pressure,
                           const std::vector<typename Fluid::CompVec>& z,
                           const std::vector<double>& face_flux,
                           const std::vector<typename Fluid::PhaseVec>& sat, 
                           const std::vector<typename Fluid::CompVec>& mass_frac,
                           const int step,
                           const std::string& filebase);

    };




    // Method implementations below.


template<class Grid, class Rock, class Fluid, class Wells, class FlowSolver, class TransportSolver>
void
BlackoilSimulator<Grid, Rock, Fluid, Wells, FlowSolver, TransportSolver>::
init(const Dune::parameter::ParameterGroup& param)
{
    using namespace Dune;
    std::string fileformat = param.getDefault<std::string>("fileformat", "cartesian");
    if (fileformat == "eclipse") {
        Dune::EclipseGridParser parser(param.get<std::string>("filename"));
        double z_tolerance = param.getDefault<double>("z_tolerance", 0.0);
        bool periodic_extension = param.getDefault<bool>("periodic_extension", false);
        bool turn_normals = param.getDefault<bool>("turn_normals", false);
        grid_.processEclipseFormat(parser, z_tolerance, periodic_extension, turn_normals);
        double perm_threshold_md = param.getDefault("perm_threshold_md", 0.0);
        double perm_threshold = Dune::unit::convert::from(perm_threshold_md, Dune::prefix::milli*Dune::unit::darcy);
        rock_.init(parser, grid_.globalCell(), perm_threshold);
        fluid_.init(parser);
        wells_.init(parser, grid_, rock_);
    } else if (fileformat == "cartesian") {
        Dune::array<int, 3> dims = {{ param.getDefault<int>("nx", 1),
                                      param.getDefault<int>("ny", 1),
                                      param.getDefault<int>("nz", 1) }};
        Dune::array<double, 3> cellsz = {{ param.getDefault<double>("dx", 1.0),
                                           param.getDefault<double>("dy", 1.0),
                                           param.getDefault<double>("dz", 1.0) }};
        grid_.createCartesian(dims, cellsz);
        double default_poro = param.getDefault("default_poro", 1.0);
        double default_perm_md = param.getDefault("default_perm_md", 100.0);
        double default_perm = unit::convert::from(default_perm_md, prefix::milli*unit::darcy);
        MESSAGE("Warning: For generated cartesian grids, we use uniform rock properties.");
        rock_.init(grid_.size(0), default_poro, default_perm);
        EclipseGridParser parser(param.get<std::string>("filename")); // Need a parser for the fluids anyway.
        fluid_.init(parser);
        wells_.init(parser, grid_, rock_);
    } else {
        THROW("Unknown file format string: " << fileformat);
    }
    flow_solver_.init(param);
    transport_solver_.init(param);
    if (param.has("timestep_file")) {
        std::ifstream is(param.get<std::string>("timestep_file").c_str());
        std::istream_iterator<double> beg(is);
        std::istream_iterator<double> end;
        report_times_.assign(beg, end);
        // File contains deltas, we want accumulated times.
        std::partial_sum(report_times_.begin(), report_times_.end(), report_times_.begin());
        ASSERT(!report_times_.empty());
        total_time_ = report_times_.back();
        initial_stepsize_ = report_times_.front();
    } else {
        total_time_ = param.getDefault("total_time", 30*unit::day);
        initial_stepsize_ = param.getDefault("initial_stepsize", 1.0*unit::day);
        increase_stepsize_ = param.getDefault("increase_stepsize", false);
        if (increase_stepsize_) {
            stepsize_increase_factor_ = param.getDefault("stepsize_increase_factor", 1.5);
            maximum_stepsize_ = param.getDefault("maximum_stepsize", 1.0*unit::day);
        } else {
            stepsize_increase_factor_ = 1.0;
            maximum_stepsize_ = 1e100;
        }
    }
    do_impes_ = param.getDefault("do_impes", false);
    output_dir_ = param.getDefault<std::string>("output_dir", "output");

    // Boundary conditions.
    typedef Dune::FlowBC BC;
    flow_bc_.resize(7);
    bool bdy_dirichlet = param.getDefault("bdy_dirichlet", false);
    if (bdy_dirichlet) {
        flow_bc_.flowCond(1) = BC(BC::Dirichlet, param.get<double>("bdy_pressure_left"));
        flow_bc_.flowCond(2) = BC(BC::Dirichlet, param.get<double>("bdy_pressure_right"));
    }

    // Gravity.
    gravity_ = 0.0;
    if (param.has("gravity")) {
        std::string g = param.get<std::string>("gravity");
        if (g == "standard") {
            gravity_[2] = Dune::unit::gravity;
        } else {
            gravity_[2] = boost::lexical_cast<double>(g);
        }
    }

    // Flow solver setup.
    flow_solver_.setup(grid_, rock_, fluid_, wells_, gravity_, flow_bc_);

    // Transport solver setup.
    transport_solver_.setup(grid_, rock_, fluid_, wells_, flow_solver_.faceTransmissibilities(), gravity_);

    // Simple source terms.
    src_.resize(grid_.numCells(), 0.0);


    // Initial state.
    if (param.getDefault("spe9_init", false)) {
        
        double zeroDepth = param.getDefault("zero_depth", 2743.2);
        
        int nx = param.getDefault<int>("nx", 24);
        int ny = param.getDefault<int>("ny", 25);
        int nz = param.getDefault<int>("nz", 15);
        
        // double datum_depth = param.getDefault<double>("datum_depth", 2753.87) - zeroDepth;
        double datum_pressure_barsa = param.getDefault<double>("datum_pressure", 248.22);
        double datum_pressure = unit::convert::from(datum_pressure_barsa, Dune::unit::barsa);
        double wo_contact_depth = param.getDefault<double>("wo_contact_depth", 3032.76) - zeroDepth;
        double go_contact_depth = param.getDefault<double>("go_contact_depth", 2682.24) - zeroDepth;
        
        double connate_water_saturation = param.getDefault<double>("connate_water_saturation", 0.151090);
        
        double initial_mixture_gas = param.getDefault("initial_mixture_gas", 247.43);
        double initial_mixture_oil = param.getDefault("initial_mixture_oil", 1.0);
        
        // Initial fluid state
        CompVec oil_sample(0.0);
        oil_sample[Fluid::Oil] = initial_mixture_oil;
        oil_sample[Fluid::Gas] = initial_mixture_gas;
        CompVec water_sample(0.0); 
        water_sample[Fluid::Water] = 1.0;
        
        cell_z_.resize(grid_.numCells());
        cell_pressure_.resize(grid_.numCells());
        
        // Datum -cell
        // For now, assume that datum_depth corresponds the centroid of cell 0 (reasonable approx)
        cell_pressure_[0] = datum_pressure;    
        double pore_vol = grid_.cellVolume(0)*rock_.porosity(0);
        typename Fluid::FluidState state = fluid_.computeState(cell_pressure_[0],oil_sample);
        cell_z_[0] = oil_sample;
        cell_z_[0] *= (1.0-connate_water_saturation)*pore_vol/state.total_phase_volume_;
        state = fluid_.computeState(cell_pressure_[0],water_sample);
        cell_z_[0][Fluid::Water] = water_sample[Fluid::Water];
        cell_z_[0][Fluid::Water] *= connate_water_saturation*pore_vol/state.total_phase_volume_;
        
        // Rest of the cells -- NOTE: Assume uniform cell properties in y-direction
        for (int i=0; i<nx; ++i) {
            int k0=i*nz; 
            for (int k=0; k<nz; ++k) {
                int kk=k0+k;
                if (i>0 && k==0) {
                    computeCellState(kk, kk-nz, wo_contact_depth, go_contact_depth, connate_water_saturation);
                } else if (k>0) { 
                    computeCellState(kk, kk-1, wo_contact_depth, go_contact_depth, connate_water_saturation);
                }
                // Copy cell properties to y-layers
                for (int j=1; j<ny; ++j) {
                    int jj = j*nx*nz + kk;
                    cell_z_[jj] = cell_z_[kk];
                    cell_pressure_[jj] = cell_pressure_[kk];
                }
            }                
        }
        
        // Write initial state to std::cout
        /*
        for (int cell = 0; cell < grid_.numCells(); ++cell) {         
            std::cout.precision(2);
            std::cout << std::fixed << std::showpoint;
            
            std::cout << std::setw(5) << cell << std::setw(12) << grid_.cellCentroid(cell)[0]
                                              << std::setw(12) << grid_.cellCentroid(cell)[1]
                                              << std::setw(12) << grid_.cellCentroid(cell)[2] 
                                              << std::setw(20) << cell_pressure_[cell][0]
                                              << std::setw(15) << cell_z_[cell][0]
                                              << std::setw(15) << cell_z_[cell][1]
                                              << std::setw(15) << cell_z_[cell][2]
                                              << std::endl;
            if ((cell+1)%nz == 0) {
                std::cout << "------------------------------------------------------------------------------------------------------------------" << std::endl;
            }
            
        }
        */
    } else if (param.getDefault("heterogenous_initial_mix", false)) {
        CompVec init_oil(0.0);
        init_oil[Fluid::Oil] = 1.0;
        CompVec init_water(0.0);
        init_water[Fluid::Water] = 1.0;

        bdy_z_ = flow_solver_.inflowMixture();

        cell_z_.resize(grid_.numCells());
        std::fill(cell_z_.begin(), cell_z_.begin() + cell_z_.size()/2, init_oil);
        std::fill(cell_z_.begin() + cell_z_.size()/2, cell_z_.end(), init_water);
        MESSAGE("******* Assuming zero capillary pressures *******");
        PhaseVec init_p(100.0*Dune::unit::barsa);
        cell_pressure_.resize(grid_.numCells(), init_p);
//         if (gravity_.two_norm() != 0.0) {
//             double ref_gravpot = grid_.cellCentroid(0)*gravity_;
//             double rho = init_z*fluid_.surfaceDensities();  // Assuming incompressible, and constant initial z.
//             for (int cell = 1; cell < grid_.numCells(); ++cell) {
//                 double press = rho*(grid_.cellCentroid(cell)*gravity_ - ref_gravpot) + cell_pressure_[0][0];
//                 cell_pressure_[cell] = PhaseVec(press);
//             }
//         }
    } else if (param.getDefault("unstable_initial_mix", false)) {
        CompVec init_oil(0.0);
        init_oil[Fluid::Oil] = 1.0;
        init_oil[Fluid::Gas] = 0.0;
        CompVec init_water(0.0);
        init_water[Fluid::Water] = 1.0;
        CompVec init_gas(0.0);
        init_gas[Fluid::Gas] = 150.0;
        bdy_z_ = flow_solver_.inflowMixture();
        cell_z_.resize(grid_.numCells());
        std::fill(cell_z_.begin()                       , cell_z_.begin() + cell_z_.size()/3    , init_water);
        std::fill(cell_z_.begin() + cell_z_.size()/3    , cell_z_.begin() + 2*(cell_z_.size()/3), init_oil);
        std::fill(cell_z_.begin() + 2*(cell_z_.size()/3), cell_z_.end()                         , init_gas);
        MESSAGE("******* Assuming zero capillary pressures *******");
        PhaseVec init_p(100.0*Dune::unit::barsa);
        cell_pressure_.resize(grid_.numCells(), init_p);

        if (gravity_.two_norm() != 0.0) {
            
            typename Fluid::FluidState state = fluid_.computeState(cell_pressure_[0], cell_z_[0]);
            cell_z_[0] *= grid_.cellVolume(0)*rock_.porosity(0)/state.total_phase_volume_;
            for (int cell = 1; cell < grid_.numCells(); ++cell) {
                double pore_vol = grid_.cellVolume(cell)*rock_.porosity(cell);
                double fluid_vol;
                int cnt =0;    
                do {
                    double rho = 0.5*((cell_z_[cell]+cell_z_[cell-1])*fluid_.surfaceDensities());
                    double press = rho*((grid_.cellCentroid(cell) - grid_.cellCentroid(cell-1))*gravity_) + cell_pressure_[cell-1][0];
                    cell_pressure_[cell] = PhaseVec(press);
                    typename Fluid::FluidState state = fluid_.computeState(cell_pressure_[cell], cell_z_[cell]);
                    fluid_vol = state.total_phase_volume_;
                    cell_z_[cell] *= pore_vol/fluid_vol;
                    ++cnt;
                } while (std::fabs((fluid_vol-pore_vol)/pore_vol) > 1.0e-8 && cnt < 10);
                
            }  
        } else {
            std::cout << "---- Exit - BlackoilSimulator.hpp: No gravity, no fun ... ----" << std::endl;
            exit(-1);
        }       
    } else {
        CompVec init_z(0.0);
        double initial_mixture_gas = param.getDefault("initial_mixture_gas", 0.0);
        double initial_mixture_oil = param.getDefault("initial_mixture_oil", 1.0);
        double initial_mixture_water = param.getDefault("initial_mixture_water", 0.0);
        init_z[Fluid::Water] = initial_mixture_water;
        init_z[Fluid::Gas] = initial_mixture_gas;
        init_z[Fluid::Oil] = initial_mixture_oil;

        bdy_z_ = flow_solver_.inflowMixture();

        cell_z_.resize(grid_.numCells(), init_z);
        MESSAGE("******* Assuming zero capillary pressures *******");
        PhaseVec init_p(param.getDefault("initial_pressure", 100.0*Dune::unit::barsa));
        cell_pressure_.resize(grid_.numCells(), init_p);
        if (gravity_.two_norm() != 0.0) {
            double ref_gravpot = grid_.cellCentroid(0)*gravity_;
            double rho = init_z*fluid_.surfaceDensities();  // Assuming incompressible, and constant initial z.
            for (int cell = 1; cell < grid_.numCells(); ++cell) {
                double press = rho*(grid_.cellCentroid(cell)*gravity_ - ref_gravpot) + cell_pressure_[0][0];
                cell_pressure_[cell] = PhaseVec(press);
            }
        }
    }
    bdy_pressure_ = 300.0*Dune::unit::barsa;
    // PhaseVec bdy_pressure_(100.0*Dune::unit::barsa); // WELLS
    // Rescale z values so that pore volume is filled exactly
    // (to get zero initial volume discrepancy).
    for (int cell = 0; cell < grid_.numCells(); ++cell) {
        double pore_vol = grid_.cellVolume(cell)*rock_.porosity(cell);
        typename Fluid::FluidState state = fluid_.computeState(cell_pressure_[cell], cell_z_[cell]);
        double fluid_vol = state.total_phase_volume_;
        cell_z_[cell] *= pore_vol/fluid_vol;
    }
    int num_faces = grid_.numFaces();
    face_pressure_.resize(num_faces);
    for (int face = 0; face < num_faces; ++face) {
        int bid = grid_.boundaryId(face);
        if (flow_bc_.flowCond(bid).isDirichlet()) {
            face_pressure_[face] = flow_bc_.flowCond(bid).pressure();
        } else {
            int c[2] = { grid_.faceCell(face, 0), grid_.faceCell(face, 1) };
            face_pressure_[face] = 0.0;
            int num = 0;
            for (int j = 0; j < 2; ++j) {
                if (c[j] >= 0) {
                    face_pressure_[face] += cell_pressure_[c[j]];
                    ++num;
                }
            }
            face_pressure_[face] /= double(num);
        }
    }

    // Set initial well perforation pressures equal to cell pressures,
    // and perforation fluxes equal to zero.
    well_perf_pressure_.clear();
    for (int well = 0; well < wells_.numWells(); ++well) {
        int num_perf = wells_.numPerforations(well);
        for (int perf = 0; perf < num_perf; ++perf) {
            int cell = wells_.wellCell(well, perf);
            well_perf_pressure_.push_back(cell_pressure_[cell][Fluid::Liquid]);
        }
    }
    well_perf_flux_.clear();
    well_perf_flux_.resize(well_perf_pressure_.size(), 0.0);
    wells_.update(grid_.numCells(), well_perf_pressure_, well_perf_flux_);
}











template<class Grid, class Rock, class Fluid, class Wells, class FlowSolver, class TransportSolver>
void
BlackoilSimulator<Grid, Rock, Fluid, Wells, FlowSolver, TransportSolver>::
simulate()
{
    double voldisclimit = flow_solver_.volumeDiscrepancyLimit();
    double stepsize = initial_stepsize_;
    double current_time = 0.0;
    int step = 0;
    std::vector<double> face_flux;
    std::vector<double> well_perf_pressure_start;
    std::vector<double> well_perf_flux_start;
    std::vector<PhaseVec> cell_pressure_start;
    std::vector<PhaseVec> face_pressure_start;
    std::vector<CompVec> cell_z_start;
    std::string output_name = output_dir_ + "/" + "blackoil-output";
    while (current_time < total_time_) {
        cell_pressure_start = cell_pressure_;
        face_pressure_start = face_pressure_;
        well_perf_pressure_start = well_perf_pressure_;
        well_perf_flux_start = well_perf_flux_;
        cell_z_start = cell_z_;

        // Do not run past total_time_.
        if (current_time + stepsize > total_time_) {
            stepsize = total_time_ - current_time;
        }
        std::cout << "\n\n================    Simulation step number " << step
                  << "    ==============="
                  << "\n      Current time (days)     " << Dune::unit::convert::to(current_time, Dune::unit::day)
                  << "\n      Current stepsize (days) " << Dune::unit::convert::to(stepsize, Dune::unit::day)
                  << "\n      Total time (days)       " << Dune::unit::convert::to(total_time_, Dune::unit::day)
                  << "\n" << std::endl;

        // Solve flow system.
        enum FlowSolver::ReturnCode result
            = flow_solver_.solve(cell_pressure_, face_pressure_, cell_z_, face_flux,
                                 well_perf_pressure_, well_perf_flux_, src_, stepsize, do_impes_);

        // Check if the flow solver succeeded.
        if (result != FlowSolver::SolveOk) {
            THROW("Flow solver refused to run due to too large volume discrepancy.");
        }

        // Update wells with new perforation pressures and fluxes.
        wells_.update(grid_.numCells(), well_perf_pressure_, well_perf_flux_);

        // Transport and check volume discrepancy.
        bool voldisc_ok = true;
        if (!do_impes_) {
            double actual_computed_time
                = transport_solver_.transport(bdy_pressure_, bdy_z_,
                                             face_flux, cell_pressure_, face_pressure_,
                                             stepsize, voldisclimit, cell_z_);
            voldisc_ok = (actual_computed_time == stepsize);
        } else {
            voldisc_ok = flow_solver_.volumeDiscrepancyAcceptable(cell_pressure_, face_pressure_, cell_z_, stepsize);
        }

        // If discrepancy too large, redo entire pressure step.
        if (!voldisc_ok) {
            std::cout << "********* Too large volume discrepancy:  Shortening (pressure) stepsize, redoing step number " << step <<" **********" << std::endl;
            stepsize *= 0.5;
            cell_pressure_ = cell_pressure_start;
            face_pressure_ = face_pressure_start;
            well_perf_pressure_ = well_perf_pressure_start;
            well_perf_flux_ = well_perf_flux_start;
            cell_z_ = cell_z_start;
            wells_.update(grid_.numCells(), well_perf_pressure_, well_perf_flux_);
            continue;
        }

        // Compute saturations and mass fractions for output purposes.
        int num_cells = grid_.numCells();
        std::vector<typename Fluid::PhaseVec> saturation(num_cells);
        std::vector<typename Fluid::PhaseVec> mass_frac(num_cells); 
        for (int cell = 0; cell < num_cells; ++cell) {
            saturation[cell] = fluid_.computeState(cell_pressure_[cell], cell_z_[cell]).saturation_;
            double totMass = cell_z_[cell]*fluid_.surfaceDensities();
            mass_frac[cell][Fluid::Water] = cell_z_[cell][Fluid::Water]*fluid_.surfaceDensities()[Fluid::Water]/totMass;
            mass_frac[cell][Fluid::Oil] = cell_z_[cell][Fluid::Oil]*fluid_.surfaceDensities()[Fluid::Oil]/totMass;
            mass_frac[cell][Fluid::Gas] = cell_z_[cell][Fluid::Gas]*fluid_.surfaceDensities()[Fluid::Gas]/totMass;
        }

        // Adjust time.
        current_time += stepsize;
        if (voldisc_ok && increase_stepsize_ && stepsize < maximum_stepsize_) {
            stepsize *= stepsize_increase_factor_;
            stepsize = std::min(maximum_stepsize_, stepsize);
        }
        // If using given timesteps, set stepsize to match.
        if (!report_times_.empty()) {
            if (current_time >= report_times_[step]) {
                output(grid_, cell_pressure_, cell_z_, face_flux, saturation, mass_frac, step, output_name);
                ++step;
                if (step == int(report_times_.size())) {
                    break;
                }
            }
            stepsize = report_times_[step] - current_time;
        } else {
            output(grid_, cell_pressure_, cell_z_, face_flux, saturation, mass_frac, step, output_name);
            ++step;
        }
    }
}







template<class Grid, class Rock, class Fluid, class Wells, class FlowSolver, class TransportSolver>
bool
BlackoilSimulator<Grid, Rock, Fluid, Wells, FlowSolver, TransportSolver>::
computeCellState(int iCell, int iRef, double wo_contact_depth, double go_contact_depth, double connate_water_saturation)
{
   const int maxCnt = 30;
   const double eps = 1.0e-8;
   
   double pore_vol_ref = grid_.cellVolume(iRef)*rock_.porosity(iRef);
   double pore_vol = grid_.cellVolume(iCell)*rock_.porosity(iCell);    	
   cell_z_[iCell] = cell_z_[iRef];
   cell_z_[iCell] *= pore_vol/pore_vol_ref;
   bool waterOnly = false;
   if (grid_.cellCentroid(iCell)[2] > wo_contact_depth) { // Maybe a too crude??
       cell_z_[iCell][Fluid::Oil] = 0.0;
       cell_z_[iCell][Fluid::Gas] = 0.0;
       waterOnly = true;
   }
   double gZ = (grid_.cellCentroid(iCell) - grid_.cellCentroid(iRef))*gravity_;
   double fluid_vol;
   double pv_ref_inv = 1.0/pore_vol_ref;
   double pv_inv = 1.0/pore_vol;
   int cnt =0;    
   do {    
       double rho = 0.5*(pv_inv*(cell_z_[iCell]*fluid_.surfaceDensities())+pv_ref_inv*(cell_z_[iRef]*fluid_.surfaceDensities()));
       double press = rho*gZ + cell_pressure_[iRef][0];
       cell_pressure_[iCell] = PhaseVec(press);
       typename Fluid::FluidState state = fluid_.computeState(cell_pressure_[iCell], cell_z_[iCell]);
       fluid_vol = state.total_phase_volume_;
       double oil_vol = state.phase_volume_[Fluid::Liquid] + state.phase_volume_[Fluid::Vapour];
       double wat_vol = state.phase_volume_[Fluid::Aqua];
       if (waterOnly) {
           cell_z_[iCell][Fluid::Water] *= pore_vol/wat_vol;
       } else {
           cell_z_[iCell][Fluid::Oil] *= (1.0-connate_water_saturation)*pore_vol/oil_vol;
           cell_z_[iCell][Fluid::Gas] *= (1.0-connate_water_saturation)*pore_vol/oil_vol;
           cell_z_[iCell][Fluid::Water] *= connate_water_saturation*pore_vol/wat_vol;
       }       
       ++cnt;
   } while (std::fabs((fluid_vol-pore_vol)/pore_vol) > eps && cnt < maxCnt);
   
   if (cnt == maxCnt) {    
       std::cout << "z_cell_[" << iCell << "]: " << cell_z_[iCell]
                 << "  pressure: " << cell_pressure_[iCell][Fluid::Liquid]
                 <<  " cnt: " << cnt 
                 << "  eps: " << std::fabs((fluid_vol-pore_vol)/pore_vol) << std::endl;
   }
                      
   return (cnt < maxCnt);
}











template<class Grid, class Rock, class Fluid, class Wells, class FlowSolver, class TransportSolver>
void
BlackoilSimulator<Grid, Rock, Fluid, Wells, FlowSolver, TransportSolver>::
output(const Grid& grid,
       const std::vector<typename Fluid::PhaseVec>& cell_pressure,
       const std::vector<typename Fluid::CompVec>& z,
       const std::vector<double>& face_flux,
       const std::vector<typename Fluid::PhaseVec>& sat, 
       const std::vector<typename Fluid::CompVec>& mass_frac,
       const int step,
       const std::string& filebase)
{
    // Ensure directory exists.
    boost::filesystem::path fpath(filebase);
    if (fpath.has_branch_path()) {
        create_directories(fpath.branch_path());
    }

    // Output to VTK.
    std::vector<typename Grid::Vector> cell_velocity;
    estimateCellVelocitySimpleInterface(cell_velocity, grid, face_flux);
    // Dune's vtk writer wants multi-component data to be flattened.
    std::vector<double> cell_pressure_flat(&*cell_pressure.front().begin(),
                                           &*cell_pressure.back().end());
    std::vector<double> cell_velocity_flat(&*cell_velocity.front().begin(),
                                           &*cell_velocity.back().end());
    std::vector<double> z_flat(&*z.front().begin(),
                               &*z.back().end());
    std::vector<double> sat_flat(&*sat.front().begin(),
                                 &*sat.back().end());
    std::vector<double> mass_frac_flat(&*mass_frac.front().begin(),
                                       &*mass_frac.back().end());
    Dune::VTKWriter<typename Grid::LeafGridView> vtkwriter(grid.leafView());
    vtkwriter.addCellData(cell_pressure_flat, "pressure", Fluid::numPhases);
    vtkwriter.addCellData(cell_velocity_flat, "velocity", Grid::dimension);
    vtkwriter.addCellData(z_flat, "z", Fluid::numComponents);
    vtkwriter.addCellData(sat_flat, "sat", Fluid::numPhases);
    vtkwriter.addCellData(mass_frac_flat, "massFrac", Fluid::numComponents);
    vtkwriter.write(filebase + '-' + boost::lexical_cast<std::string>(step),
                    Dune::VTKOptions::ascii);

    // Dump data for Matlab.
    std::vector<double> zv[Fluid::numComponents];
    for (int comp = 0; comp < Fluid::numComponents; ++comp) {
        zv[comp].resize(grid.numCells());
        for (int cell = 0; cell < grid.numCells(); ++cell) {
            zv[comp][cell] = z[cell][comp];
        }
    }
    std::vector<double> sv[Fluid::numPhases];
    for (int phase = 0; phase < Fluid::numPhases; ++phase) {
        sv[phase].resize(grid.numCells());
        for (int cell = 0; cell < grid.numCells(); ++cell) {
            sv[phase][cell] = sat[cell][phase];
        }
    }
    std::string matlabdumpname(filebase + "-");
    matlabdumpname += boost::lexical_cast<std::string>(step);
    matlabdumpname += ".dat";
    std::ofstream dump(matlabdumpname.c_str());
    dump.precision(15);
    int num_cells = cell_pressure.size();
    std::vector<double> liq_press(num_cells);
    for (int cell = 0; cell < num_cells; ++cell) {
        liq_press[cell] = cell_pressure[cell][Fluid::Liquid];
    }
    std::copy(liq_press.begin(), liq_press.end(),
              std::ostream_iterator<double>(dump, " "));
    dump << '\n';
    for (int comp = 0; comp < Fluid::numComponents; ++comp) {
        std::copy(zv[comp].begin(), zv[comp].end(),
                  std::ostream_iterator<double>(dump, " "));
        dump << '\n';
    }
    for (int phase = 0; phase < Fluid::numPhases; ++phase) {
        std::copy(sv[phase].begin(), sv[phase].end(),
                  std::ostream_iterator<double>(dump, " "));
        dump << '\n';
    }
}





} // namespace Opm





#endif // OPM_BLACKOILSIMULATOR_HEADER_INCLUDED
