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

#ifndef OPM_BLACKOILINITIALIZATION_HEADER_INCLUDED
#define OPM_BLACKOILINITIALIZATION_HEADER_INCLUDED

#include <opm/core/utility/parameters/ParameterGroup.hpp>

namespace Opm
{

    /// Base class for initialization codes.
    template <class Simulator>
    class BlackoilInitialization
    {
    public:
        typedef typename Simulator::State State;
        typedef typename Simulator::Grid Grid;
        typedef typename Simulator::Fluid Fluid;

        virtual void init(const Opm::parameter::ParameterGroup& param,
                          const Grid& grid,
                          const Fluid& fluid,
                          typename Grid::Vector gravity,
                          State& simstate) = 0;
    };



    /// Initialize basic cases.
    template <class Simulator>
    class BasicInitialization : public BlackoilInitialization<Simulator>
    {
    public:
        typedef typename Simulator::State State;
        typedef typename Simulator::Grid Grid;
        typedef typename Simulator::Fluid Fluid;

        virtual void init(const Opm::parameter::ParameterGroup& param,
                          const Grid& grid,
                          const Fluid& fluid,
                          typename Grid::Vector gravity,
                          State& simstate)
        {
            typedef typename Fluid::CompVec CompVec;
            typedef typename Fluid::PhaseVec PhaseVec;

            if (param.getDefault("heterogenous_initial_mix", false)) {
                CompVec init_oil(0.0);
                init_oil[Fluid::Oil] = 1.0;
                CompVec init_water(0.0);
                init_water[Fluid::Water] = 1.0;
                simstate.cell_z_.resize(grid.numCells());
                std::fill(simstate.cell_z_.begin(), simstate.cell_z_.begin() + simstate.cell_z_.size()/2, init_oil);
                std::fill(simstate.cell_z_.begin() + simstate.cell_z_.size()/2, simstate.cell_z_.end(), init_water);
                MESSAGE("******* Assuming zero capillary pressures *******");
                PhaseVec init_p(100.0*Opm::unit::barsa);
                simstate.cell_pressure_.resize(grid.numCells(), init_p);
                //         if (gravity.two_norm() != 0.0) {
                //             double ref_gravpot = grid.cellCentroid(0)*gravity;
                //             double rho = init_z*fluid_.surfaceDensities();  // Assuming incompressible, and constant initial z.
                //             for (int cell = 1; cell < grid.numCells(); ++cell) {
                //                 double press = rho*(grid.cellCentroid(cell)*gravity - ref_gravpot) + simstate.cell_pressure_[0][0];
                //                 simstate.cell_pressure_[cell] = PhaseVec(press);
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
                simstate.cell_z_.resize(grid.numCells());
                std::fill(simstate.cell_z_.begin(),
                          simstate.cell_z_.begin() + simstate.cell_z_.size()/3,
                          init_water);
                std::fill(simstate.cell_z_.begin() + simstate.cell_z_.size()/3,
                          simstate.cell_z_.begin() + 2*(simstate.cell_z_.size()/3),
                          init_oil);
                std::fill(simstate.cell_z_.begin() + 2*(simstate.cell_z_.size()/3),
                          simstate.cell_z_.end(),
                          init_gas);
                MESSAGE("******* Assuming zero capillary pressures *******");
                PhaseVec init_p(100.0*Opm::unit::barsa);
                simstate.cell_pressure_.resize(grid.numCells(), init_p);

                if (gravity.two_norm() != 0.0) {
            
                    typename Fluid::FluidState state = fluid.computeState(simstate.cell_pressure_[0], simstate.cell_z_[0]);
                    simstate.cell_z_[0] *= 1.0/state.total_phase_volume_density_;
                    for (int cell = 1; cell < grid.numCells(); ++cell) {
                        double fluid_vol_dens;
                        int cnt =0;    
                        do {
                            double rho = 0.5*((simstate.cell_z_[cell]+simstate.cell_z_[cell-1])*fluid.surfaceDensities());
                            double press = rho*((grid.cellCentroid(cell) - grid.cellCentroid(cell-1))*gravity) + simstate.cell_pressure_[cell-1][0];
                            simstate.cell_pressure_[cell] = PhaseVec(press);
                            typename Fluid::FluidState state = fluid.computeState(simstate.cell_pressure_[cell], simstate.cell_z_[cell]);
                            fluid_vol_dens = state.total_phase_volume_density_;
                            simstate.cell_z_[cell] *= 1.0/fluid_vol_dens;
                            ++cnt;
                        } while (std::fabs((fluid_vol_dens-1.0)) > 1.0e-8 && cnt < 10);
                
                    }  
                } else {
                    std::cout << "---- Exit - BlackoilSimulator.hpp: No gravity, no fun ... ----" << std::endl;
                    exit(-1);
                } 
            } else if (param.getDefault("CO2-injection", false)) {
                CompVec init_water(0.0);
                // Initially water filled (use Oil-component for water in order
                // to utilise blackoil mechanisms for brine-co2 interaction)          
                init_water[Fluid::Oil] = 1.0;  
                simstate.cell_z_.resize(grid.numCells());
                std::fill(simstate.cell_z_.begin(),simstate.cell_z_.end(),init_water);

                double datum_pressure_barsa = param.getDefault<double>("datum_pressure", 200.0);
                double datum_pressure = Opm::unit::convert::from(datum_pressure_barsa, Opm::unit::barsa);
                PhaseVec init_p(datum_pressure);
                simstate.cell_pressure_.resize(grid.numCells(), init_p);

                // Simple initial condition based on "incompressibility"-assumption
                double zMin = grid.cellCentroid(0)[2];
                for (int cell = 1; cell < grid.numCells(); ++cell) {
                    if (grid.cellCentroid(cell)[2] < zMin)
                        zMin = grid.cellCentroid(cell)[2];
                }

                typename Fluid::FluidState state = fluid.computeState(init_p, init_water);
		simstate.cell_z_[0] *= 1.0/state.total_phase_volume_density_;
                double density = (init_water*fluid.surfaceDensities())/state.total_phase_volume_density_;

                for (int cell = 0; cell < grid.numCells(); ++cell) {
                    double pressure(datum_pressure + (grid.cellCentroid(cell)[2] - zMin)*gravity[2]*density);
                    simstate.cell_pressure_[cell] = PhaseVec(pressure);
                    state = fluid.computeState(simstate.cell_pressure_[cell], simstate.cell_z_[cell]);
                    simstate.cell_z_[cell] *= 1.0/state.total_phase_volume_density_;
                }       
            } else {
                CompVec init_z(0.0);
                double initial_mixture_gas = param.getDefault("initial_mixture_gas", 0.0);
                double initial_mixture_oil = param.getDefault("initial_mixture_oil", 1.0);
                double initial_mixture_water = param.getDefault("initial_mixture_water", 0.0);
                init_z[Fluid::Water] = initial_mixture_water;
                init_z[Fluid::Gas] = initial_mixture_gas;
                init_z[Fluid::Oil] = initial_mixture_oil;

                simstate.cell_z_.resize(grid.numCells(), init_z);
                MESSAGE("******* Assuming zero capillary pressures *******");
                PhaseVec init_p(param.getDefault("initial_pressure", 100.0*Opm::unit::barsa));
                simstate.cell_pressure_.resize(grid.numCells(), init_p);
                if (gravity.two_norm() != 0.0) {
                    double ref_gravpot = grid.cellCentroid(0)*gravity;
                    double rho = init_z*fluid.surfaceDensities();  // Assuming incompressible, and constant initial z.
                    for (int cell = 1; cell < grid.numCells(); ++cell) {
                        double press = rho*(grid.cellCentroid(cell)*gravity - ref_gravpot) + simstate.cell_pressure_[0][0];
                        simstate.cell_pressure_[cell] = PhaseVec(press);
                    }
                }
            }
        }
    };






    /// Initialize SPE9 type case.
    template <class Simulator>
    class SPE9Initialization : public BlackoilInitialization<Simulator>
    {
    public:
        typedef typename Simulator::State State;
        typedef typename Simulator::Grid Grid;
        typedef typename Simulator::Fluid Fluid;

        virtual void init(const Opm::parameter::ParameterGroup& param,
                          const Grid& grid,
                          const Fluid& fluid,
                          typename Grid::Vector gravity,
                          State& simstate)
        {
            typedef typename Fluid::CompVec CompVec;

            double zeroDepth = param.getDefault("zero_depth", 2743.2);
            int nx = param.getDefault<int>("nx", 24);
            int ny = param.getDefault<int>("ny", 25);
            int nz = param.getDefault<int>("nz", 15);
        
            // double datum_depth = param.getDefault<double>("datum_depth", 2753.87) - zeroDepth;
            double datum_pressure_barsa = param.getDefault<double>("datum_pressure", 248.22);
            double datum_pressure = Opm::unit::convert::from(datum_pressure_barsa, Opm::unit::barsa);
            double wo_contact_depth = param.getDefault<double>("wo_contact_depth", 3032.76) - zeroDepth;
            double go_contact_depth = param.getDefault<double>("go_contact_depth", 2682.24) - zeroDepth;
        
            double connate_water_saturation = param.getDefault<double>("connate_water_saturation", 0.151090);
            double residual_oil_saturation = param.getDefault<double>("residual_oil_saturation", 0.118510);
        
            double initial_mixture_gas = param.getDefault("initial_mixture_gas", 247.43);
            double initial_mixture_oil = param.getDefault("initial_mixture_oil", 1.0);

            // Initial fluid state
            CompVec oil_sample(0.0);
            oil_sample[Fluid::Oil] = initial_mixture_oil;
            oil_sample[Fluid::Gas] = initial_mixture_gas;
            CompVec water_sample(0.0); 
            water_sample[Fluid::Water] = 1.0;
        
            simstate.cell_z_.resize(grid.numCells());
            simstate.cell_pressure_.resize(grid.numCells());
        
            // Datum -cell
            // For now, assume that datum_depth corresponds the centroid of cell 0 (reasonable approx)
            simstate.cell_pressure_[0] = datum_pressure;    
            typename Fluid::FluidState state = fluid.computeState(simstate.cell_pressure_[0],oil_sample);
            simstate.cell_z_[0] = oil_sample;
            simstate.cell_z_[0] *= (1.0-connate_water_saturation)/state.total_phase_volume_density_;
            state = fluid.computeState(simstate.cell_pressure_[0],water_sample);
            simstate.cell_z_[0][Fluid::Water] = water_sample[Fluid::Water];
            simstate.cell_z_[0][Fluid::Water] *= connate_water_saturation/state.total_phase_volume_density_;
            // Rest of the cells -- NOTE: Assume uniform cell properties in y-direction
            for (int i=0; i<nx; ++i) {
                int k0=i*nz; 
                for (int k=0; k<nz; ++k) {
                    int kk=k0+k;
                    if (i>0 && k==0) {
                        computeCellState(grid, fluid, gravity,
                                         kk, kk-nz, wo_contact_depth, go_contact_depth, connate_water_saturation,
                                         residual_oil_saturation, simstate);
                    } else if (k>0) { 
                        computeCellState(grid, fluid, gravity,
                                         kk, kk-1, wo_contact_depth, go_contact_depth, connate_water_saturation,
                                         residual_oil_saturation, simstate);
                    }
                    // Copy cell properties to y-layers
                    for (int j=1; j<ny; ++j) {
                        int jj = j*nx*nz + kk;
                        simstate.cell_z_[jj] = simstate.cell_z_[kk];
                        simstate.cell_pressure_[jj] = simstate.cell_pressure_[kk];
                    }
                }                
            }
        }



        bool computeCellState(const Grid& grid,
                              const Fluid& fluid,
                              typename Grid::Vector gravity,
                              int iCell,
                              int iRef,
                              double wo_contact_depth,
                              double go_contact_depth,
                              double connate_water_saturation,
                              double residual_oil_saturation,
                              State& simstate)
        {
            typedef typename Fluid::PhaseVec PhaseVec;

            const int maxCnt = 30;
            const double eps = 1.0e-8;    	
            simstate.cell_z_[iCell] = simstate.cell_z_[iRef];

            bool below_wo_contact = false;
            if (grid.cellCentroid(iCell)[2] > wo_contact_depth)
                below_wo_contact = true;

            double gZ = (grid.cellCentroid(iCell) - grid.cellCentroid(iRef))*gravity;
            double fluid_vol_dens;
            int cnt =0;    
            do {
                double rho = 0.5*(simstate.cell_z_[iCell]*fluid.surfaceDensities()
                                  + simstate.cell_z_[iRef]*fluid.surfaceDensities());
                double press = rho*gZ + simstate.cell_pressure_[iRef][0];
                simstate.cell_pressure_[iCell] = PhaseVec(press);
                typename Fluid::FluidState state = fluid.computeState(simstate.cell_pressure_[iCell], simstate.cell_z_[iCell]);
                fluid_vol_dens = state.total_phase_volume_density_;
                double oil_vol_dens = state.phase_volume_density_[Fluid::Liquid]
                    + state.phase_volume_density_[Fluid::Vapour];
                double wat_vol_dens = state.phase_volume_density_[Fluid::Aqua];
                if (below_wo_contact) {
                    simstate.cell_z_[iCell][Fluid::Oil] *= residual_oil_saturation/oil_vol_dens;
                    simstate.cell_z_[iCell][Fluid::Gas] *= residual_oil_saturation/oil_vol_dens;
                    simstate.cell_z_[iCell][Fluid::Water] *= (1.0-residual_oil_saturation)/wat_vol_dens;
                } else {
                    simstate.cell_z_[iCell][Fluid::Oil] *= (1.0-connate_water_saturation)/oil_vol_dens;
                    simstate.cell_z_[iCell][Fluid::Gas] *= (1.0-connate_water_saturation)/oil_vol_dens;
                    simstate.cell_z_[iCell][Fluid::Water] *= connate_water_saturation/wat_vol_dens;
                }       
                ++cnt;
            } while (std::fabs(fluid_vol_dens-1.0) > eps && cnt < maxCnt);
   
            if (cnt == maxCnt) {    
                std::cout << "z_cell_[" << iCell << "]: " << simstate.cell_z_[iCell]
                          << "  pressure: " << simstate.cell_pressure_[iCell][Fluid::Liquid]
                          <<  " cnt: " << cnt 
                          << "  eps: " << std::fabs(fluid_vol_dens-1.0) << std::endl;
            }
                      
            return (cnt < maxCnt);
        }

    };

} // namespace Opm

#endif // OPM_BLACKOILINITIALIZATION_HEADER_INCLUDED
