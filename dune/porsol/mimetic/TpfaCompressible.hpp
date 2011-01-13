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

#ifndef OPM_TPFACOMPRESSIBLE_HEADER_INCLUDED
#define OPM_TPFACOMPRESSIBLE_HEADER_INCLUDED


#include "../opmpressure/src/TPFACompressiblePressureSolver.hpp"

#include <dune/common/ErrorMacros.hpp>
#include <dune/common/SparseTable.hpp>
#include <dune/porsol/common/LinearSolverISTL.hpp>


namespace Dune
{


    template <class GridInterface,
              class RockInterface,
              class FluidInterface,
              class BCInterface>
    class TpfaCompressible
    {
    public:
        typedef TPFACompressiblePressureSolver PressureSolver;

        /// @brief
        ///    Default constructor. Does nothing.
        TpfaCompressible()
            : pgrid_(0), prock_(0), pfluid_(0)
        {
        }


        /// @brief
        ///    Initializes run-time parameters of the solver.
        void init(const parameter::ParameterGroup& param)
        {
            // Initialize inflow mixture to a fixed, user-provided mix.
            typename FluidInterface::CompVec mix(0.0);
            const int nc = FluidInterface::numComponents;
            double inflow_mixture_gas = param.getDefault("inflow_mixture_gas", 1.0);
            double inflow_mixture_oil = param.getDefault("inflow_mixture_oil", 0.0);
            switch (nc) {
            case 2:
                mix[FluidInterface::Gas] = inflow_mixture_gas;
                mix[FluidInterface::Oil] = inflow_mixture_oil;
                break;
            case 3: {
                double inflow_mixture_water = param.getDefault("inflow_mixture_water", 0.0);
                mix[FluidInterface::Water] = inflow_mixture_water;
                mix[FluidInterface::Gas] = inflow_mixture_gas;
                mix[FluidInterface::Oil] = inflow_mixture_oil;
                break;
            }
            default:
                THROW("Unhandled number of components: " << nc);
            }
            inflow_mixture_ = mix;
            linsolver_.init(param);
            flux_rel_tol_ = param.getDefault("flux_rel_tol", 1e-5);
            max_num_iter_ = param.getDefault("max_num_iter", 15);
            max_relative_voldiscr_ = param.getDefault("max_relative_voldiscr", 0.15);
        }




        /// @brief
        ///    Accessor for the inflow mixture.
        typename FluidInterface::CompVec inflowMixture() const
        {
            return inflow_mixture_;
        }




        /// @brief
        ///    Setup routine, does grid/rock-dependent initialization.
        ///
        /// @param [in] grid
        ///    The grid.
        ///
        /// @param [in] rock
        ///    The cell-wise permeabilities and porosities.
        ///
        /// @param [in] fluid
        ///    Fluid properties.
        ///
        /// @param [in] wells
        ///    Well specifications.
        ///
        /// @param [in] grav
        ///    Gravity vector.  Its Euclidian two-norm value
        ///    represents the strength of the gravity field (in units
        ///    of m/s^2) while its direction is the direction of
        ///    gravity in the current model.
        ///
        /// @param [in] bc
        ///    Boundary conditions.
        ///
        template <class WellsInterface>
        void setup(const GridInterface&         grid,
                   const RockInterface&         rock,
                   const FluidInterface&        fluid,
                   const WellsInterface&        wells,
                   const typename GridInterface::Vector& grav,
                   const BCInterface& bc)
        {
            pgrid_ = &grid;
            prock_ = &rock;
            pfluid_ = &fluid;

            // Extract perm tensors.
            const double* perm = &(rock.permeability(0)(0,0));
            poro_.clear();
            poro_.resize(grid.numCells(), 1.0);
            for (int i = 0; i < grid.numCells(); ++i) {
                poro_[i] = rock.porosity(i);
            }
            // Initialize 
            psolver_.init(grid, wells, perm, &poro_[0], grav);

            // Build bctypes_ and bcvalues_.
            int num_faces = grid.numFaces();
            bctypes_.clear();
            bctypes_.resize(num_faces, PressureSolver::FBC_UNSET);
            bcvalues_.clear();
            bcvalues_.resize(num_faces, 0.0);
            for (int face = 0; face < num_faces; ++face) {
                int bid = pgrid_->boundaryId(face);
                if (bid == 0) {
                    bctypes_[face] = PressureSolver::FBC_UNSET;
                    continue;
                }
                FlowBC face_bc = bc.flowCond(bid);
                if (face_bc.isDirichlet()) {
                    bctypes_[face] = PressureSolver::FBC_PRESSURE;
                    bcvalues_[face] = face_bc.pressure();
                } else if (face_bc.isNeumann()) {
                    bctypes_[face] = PressureSolver::FBC_FLUX;
                    bcvalues_[face] = face_bc.outflux(); // TODO: may have to switch sign here depending on orientation.
                    if (bcvalues_[face] != 0.0) {
                        THROW("Nonzero Neumann conditions not yet properly implemented (signs must be fixed)");
                    }
                } else {
                    THROW("Unhandled boundary condition type.");
                }
            }
        }




        double volumeDiscrepancyLimit() const
        {
            return max_relative_voldiscr_;
        }




        bool volumeDiscrepancyAcceptable(const std::vector<typename FluidInterface::PhaseVec>& cell_pressure,
                                         const std::vector<typename FluidInterface::PhaseVec>& face_pressure,
                                         const std::vector<typename FluidInterface::CompVec>& cell_z,
                                         const double dt)
        {
            computeFluidProps(cell_pressure, face_pressure, cell_z, dt);
            double rel_voldiscr = *std::max_element(fp_.relvoldiscr.begin(), fp_.relvoldiscr.end());
            if (rel_voldiscr > max_relative_voldiscr_) {
                std::cout << "    Relative volume discrepancy too large: " << rel_voldiscr << std::endl;
                return false;
            } else {
                return true;
            }
        }

        enum ReturnCode { SolveOk, VolumeDiscrepancyTooLarge };


        /// @brief
        ///    Construct and solve system of linear equations for the
        ///    phase pressure values on cells and faces, also compute
        ///    total face fluxes.
        ///
        /// @param [inout] cell_pressure
        ///    Phase pressures per cell.
        ///
        /// @param [inout] face_pressure
        ///    Phase pressures per face.
        ///
        /// @param [inout] cell_z
        ///    Surface volume per cell. Only changed if the @code
        ///    transport @endcode argument is true.
        ///
        /// @param [out] face_flux
        ///    Total (summed over all phases) volume flux (signed)
        ///    across each face.
        ///
        /// @param [out] well_pressures
        ///    Pressure in each well perforation.
        ///
        /// @param [out] well_fluxes
        ///    Total (summed over all phases) volume flux (signed,
        ///    positive meaning injection) from each well perforation.
        ///
        /// @param [in] src
        ///    Explicit source terms.  One scalar value for each grid
        ///    cell representing the rate (in units of m^3/s) of fluid
        ///    being injected into (>0) or extracted from (<0) a given
        ///    grid cell.
        ///
        /// @param [in] dt
        ///    Timestep for pressure solver.
        ///
        /// @param [in] transport
        ///    If true, modify @code cell_z @endcode by IMPES scheme.
        ///
        ReturnCode solve(std::vector<typename FluidInterface::PhaseVec>& cell_pressure,
                         std::vector<typename FluidInterface::PhaseVec>& face_pressure,
                         std::vector<typename FluidInterface::CompVec>& cell_z,
                         std::vector<double>& face_flux,
                         std::vector<double>& well_pressures,
                         std::vector<double>& well_fluxes,
                         const std::vector<double>& src,
                         const double dt,
                         bool transport = false)
        {
            // Set starting pressures.
            int num_faces = pgrid_->numFaces();

            // Assemble and solve.
            int num_cells = cell_z.size();
            std::vector<double> cell_pressure_scalar(num_cells);
            // Set initial pressure to Liquid phase pressure. \TODO what is correct with capillary pressure?
            for (int cell = 0; cell < num_cells; ++cell) {
                cell_pressure_scalar[cell] = cell_pressure[cell][FluidInterface::Liquid];
            }
            std::vector<double> initial_voldiscr;
            std::vector<double> face_pressure_scalar;
            std::vector<double> start_face_flux;
            if (int(face_flux.size()) != num_faces) {
                face_flux.clear();
                face_flux.resize(num_faces, 0.0);
            }
            bool converged = false;
            for (int i = 0; i < max_num_iter_; ++i) {
                start_face_flux = face_flux;
                // (Re-)compute fluid properties.
                computeFluidProps(cell_pressure, face_pressure, cell_z, dt);
                if (i == 0) {
                    initial_voldiscr = fp_.voldiscr;
                    double rel_voldiscr = *std::max_element(fp_.relvoldiscr.begin(), fp_.relvoldiscr.end());
                    if (rel_voldiscr > max_relative_voldiscr_) {
                        std::cout << "    Relative volume discrepancy too large: " << rel_voldiscr << std::endl;
                        return VolumeDiscrepancyTooLarge;
                    }
                }

                // Assemble system matrix and rhs.
                psolver_.assemble(src, bctypes_, bcvalues_, dt,
                                  fp_.totcompr, initial_voldiscr, fp_.cellA, fp_.faceA, fp_.phasemobf,
                                  cell_pressure_scalar, &(pfluid_->surfaceDensities()[0]));
                // Solve system.
                PressureSolver::LinearSystem s;
                psolver_.linearSystem(s);
                LinearSolverISTL::LinearSolverResults res = linsolver_.solve(s.n, s.nnz, s.ia, s.ja, s.sa, s.b, s.x);
                if (!res.converged) {
                    THROW("Linear solver failed to converge in " << res.iterations << " iterations.\n"
                          << "Residual reduction achieved is " << res.reduction << '\n');
                }
                // Get pressures and face fluxes.
                psolver_.computePressuresAndFluxes(cell_pressure_scalar, face_pressure_scalar, face_flux,
                                                   well_pressures, well_fluxes);
                // Copy to phase pressures. \TODO handle capillary pressure.
                for (int cell = 0; cell < num_cells; ++cell) {
                    cell_pressure[cell] = cell_pressure_scalar[cell];
                }
                for (int face = 0; face < num_faces; ++face) {
                    face_pressure[face] = face_pressure_scalar[face];
                }

                // Test for convergence.
                double max_flux = std::max(std::fabs(*std::min_element(face_flux.begin(), face_flux.end())),
                                           std::fabs(*std::max_element(face_flux.begin(), face_flux.end())));
                double flux_change_infnorm = 0.0;
                for (int face = 0; face < num_faces; ++face) {
                    flux_change_infnorm = std::max(flux_change_infnorm, std::fabs(face_flux[face] - start_face_flux[face]));
                }
                double rel_difference = flux_change_infnorm/max_flux;
                std::cout << "    Iteration in pressure solver complete. Relative flux change: " << rel_difference << std::endl;
                if (rel_difference < flux_rel_tol_) {
                    std::cout << "    Pressure solver converged. Number of iterations: " << i + 1 << std::endl;
                    converged = true;
                    break;
                }

                // DUMP HACK
//                 std::string fname("facepress-");
//                 fname += boost::lexical_cast<std::string>(i);
//                 std::ofstream f(fname.c_str());
//                 f.precision(15);
//                 std::copy(face_pressure_scalar.begin(), face_pressure_scalar.end(),
//                           std::ostream_iterator<double>(f, "\n"));
            }

            if (!converged) {
                THROW("Pressure solver failed to converge after " << max_num_iter_ << " iterations.");
            }

            if (transport) {
                psolver_.explicitTransport(dt, &cell_pressure_scalar[0], &(cell_z[0][0]));
            }

            return SolveOk;
        }




    private:
        void computeFluidProps(const std::vector<typename FluidInterface::PhaseVec>& phase_pressure,
                               const std::vector<typename FluidInterface::PhaseVec>& phase_pressure_face,
                               const std::vector<typename FluidInterface::CompVec>& cell_z,
                               const double dt)
        {
            fp_.compute(*pgrid_, *prock_, *pfluid_, phase_pressure, phase_pressure_face, cell_z, inflow_mixture_, dt);
        }

        const GridInterface* pgrid_;
        const RockInterface* prock_;
        const FluidInterface* pfluid_;
        typename FluidInterface::FluidData fp_;
        std::vector<double> poro_;
        PressureSolver psolver_;
        LinearSolverISTL linsolver_;
        std::vector<PressureSolver::FlowBCTypes> bctypes_;
        std::vector<double> bcvalues_;

        typename FluidInterface::CompVec inflow_mixture_;
        double flux_rel_tol_;
        int max_num_iter_;
        double max_relative_voldiscr_;
    };


} // namespace Dune



#endif // OPM_TPFACOMPRESSIBLE_HEADER_INCLUDED
