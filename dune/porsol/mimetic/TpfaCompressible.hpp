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
            : pgrid_(0)
        {
        }


        /// @brief
        ///    Initializes run-time parameters of the solver.
        void init(const parameter::ParameterGroup& param)
        {
            // Initialize inflow mixture to a fixed, user-provided mix.
            typename FluidInterface::CompVec mix(0.0);
            const int nc = FluidInterface::numComponents;
            double inflow_mixture_oil = param.getDefault("inflow_mixture_oil", 0.0);
            double inflow_mixture_gas = param.getDefault("inflow_mixture_gas", nc == 3 ? 0.0 : 1.0);
            switch (nc) {
            case 2:
                mix[0] = inflow_mixture_oil;
                mix[1] = inflow_mixture_gas;
                break;
            case 3: {
                double inflow_mixture_water = param.getDefault("inflow_mixture_water", 1.0);
                mix[0] = inflow_mixture_water;
                mix[1] = inflow_mixture_oil;
                mix[0] = inflow_mixture_gas;
                break;
            }
            default:
                THROW("Unhandled number of components: " << nc);
            }
            inflow_mixture_ = mix;
            linsolver_.init(param);
            num_iter_ = param.getDefault("num_iter", 5);
        }


        /// @brief
        ///    Setup routine, does grid/rock-dependent initialization.
        /// @param [in] grid
        ///    The grid.
        ///
        /// @param [in] rock
        ///    The cell-wise permeabilities and porosities.
        ///
        /// @param [in] grav
        ///    Gravity vector.  Its Euclidian two-norm value
        ///    represents the strength of the gravity field (in units
        ///    of m/s^2) while its direction is the direction of
        ///    gravity in the current model.
        void setup(const GridInterface&         grid,
                   const RockInterface&         rock,
                   const typename GridInterface::Vector& grav)
        {
            pgrid_ = &grid;
            if (grav.two_norm() > 0.0) {
                THROW("TpfaCompressible does not handle gravity yet.");
            } 
            // Extract perm tensors.
            const double* perm = &(rock.permeability(0)(0,0));
            poro_.clear();
            poro_.resize(grid.numCells(), 1.0);
            for (int i = 0; i < grid.numCells(); ++i) {
                poro_[i] = rock.porosity(i);
            }
            // Initialize 
            psolver_.init(grid, perm, &poro_[0]);
        }





        /// @brief
        ///    Construct and solve system of linear equations for the
        ///    pressure values on each interface/contact between
        ///    neighbouring grid cells.  Recover cell pressure and
        ///    interface fluxes.  Following a call to @code solve()
        ///    @encode, you may recover the flow solution from the
        ///    @code getSolution() @endcode method.
        ///
        ///
        /// @param [in] src
        ///    Explicit source terms.  One scalar value for each grid
        ///    cell representing the rate (in units of m^3/s) of fluid
        ///    being injected into (>0) or extracted from (<0) a given
        ///    grid cell.
        ///
        /// @param [in] residual_tolerance
        ///    Control parameter for iterative linear solver software.
        ///    The iteration process is terminated when the norm of
        ///    the linear system residual is less than @code
        ///    residual_tolerance @endcode times the initial residual.
        ///
        /// @param [in] linsolver_verbosity
        ///    Control parameter for iterative linear solver software.
        ///    Verbosity level 0 prints nothing, level 1 prints
        ///    summary information, level 2 prints data for each
        ///    iteration.
        ///
        /// @param [in] linsolver_type
        ///    Control parameter for iterative linear solver software.
        ///    Type 0 selects a BiCGStab solver, type 1 selects AMG/CG.
        ///
        void solve(const FluidInterface& fluid,
                   const std::vector<typename FluidInterface::PhaseVec>& initial_phase_pressure,
                   const std::vector<typename FluidInterface::CompVec>& z,
                   const BCInterface& bc,
                   const std::vector<double>& src,
                   const double dt)
        {
            // Set starting pressures.
            int num_faces = pgrid_->numFaces();
            std::vector<typename FluidInterface::PhaseVec> phase_pressure = initial_phase_pressure;
            std::vector<typename FluidInterface::PhaseVec> phase_pressure_face(num_faces);
            for (int face = 0; face < num_faces; ++face) {
                int c[2] = { pgrid_->faceCell(face, 0), pgrid_->faceCell(face, 1) };
                phase_pressure_face[face] = 0.0;
                int num = 0;
                for (int j = 0; j < 2; ++j) {
                    if (c[j] >= 0) {
                        phase_pressure_face[face] += phase_pressure[c[j]];
                        ++num;
                    }
                }
                phase_pressure_face[face] /= double(num);
            }

            // Build bctypes and bcvalues.
            std::vector<PressureSolver::FlowBCTypes> bctypes(num_faces, PressureSolver::FBC_UNSET);
            std::vector<double> bcvalues(num_faces, 0.0);
            for (int face = 0; face < num_faces; ++face) {
                int bid = pgrid_->boundaryId(face);
                if (bid == 0) {
                    bctypes[face] = PressureSolver::FBC_UNSET;
                    continue;
                }
                FlowBC face_bc = bc.flowCond(bid);
                if (face_bc.isDirichlet()) {
                    bctypes[face] = PressureSolver::FBC_PRESSURE;
                    bcvalues[face] = face_bc.pressure();
                    phase_pressure_face[face] = face_bc.pressure();
                } else if (face_bc.isNeumann()) {
                    bctypes[face] = PressureSolver::FBC_FLUX;
                    bcvalues[face] = face_bc.outflux(); // TODO: may have to switch sign here depending on orientation.
                    if (bcvalues[face] != 0.0) {
                        THROW("Nonzero Neumann conditions not yet properly implemented (signs must be fixed)");
                    }
                } else {
                    THROW("Unhandled boundary condition type.");
                }
            }

            // Assemble and solve.
            // Set initial pressure to Liquid phase pressure. \TODO what is correct with capillary pressure?
            int num_cells = z.size();
            flow_solution_.pressure_.resize(num_cells);
            for (int cell = 0; cell < num_cells; ++cell) {
                flow_solution_.pressure_[cell] = phase_pressure[cell][FluidInterface::Liquid];
            }
            std::vector<double> face_pressure;
            std::vector<double> face_flux;
            std::vector<double> initial_voldiscr;
            for (int i = 0; i < num_iter_; ++i) {
                // (Re-)compute fluid properties.
                computeFluidProps(fluid, phase_pressure, phase_pressure_face, z);
                if (i == 0) {
                    initial_voldiscr = fp_.voldiscr;
                }

                // Assemble system matrix and rhs.
                psolver_.assemble(src, bctypes, bcvalues, dt,
                                  fp_.totcompr, initial_voldiscr, fp_.cellA, fp_.faceA, fp_.phasemobf,
                                  flow_solution_.pressure_);
                // Solve system.
                PressureSolver::LinearSystem s;
                psolver_.linearSystem(s);
                LinearSolverISTL::LinearSolverResults res = linsolver_.solve(s.n, s.nnz, s.ia, s.ja, s.sa, s.b, s.x);
                if (!res.converged) {
                    THROW("Linear solver failed to converge in " << res.iterations << " iterations.\n"
                          << "Residual reduction achieved is " << res.reduction << '\n');
                }
                // Get pressures and face fluxes.
                flow_solution_.clear();
                psolver_.computePressuresAndFluxes(flow_solution_.pressure_, face_pressure, face_flux);
                // Copy to phase pressures. \TODO handle capillary pressure.
                for (int cell = 0; cell < num_cells; ++cell) {
                    phase_pressure[cell] = flow_solution_.pressure_[cell];
                }
                for (int face = 0; face < num_faces; ++face) {
                    phase_pressure_face[face] = face_pressure[face];
                }

                // DUMP HACK
                std::string fname("facepress-");
                fname += boost::lexical_cast<std::string>(i);
                std::ofstream f(fname.c_str());
                f.precision(15);
                std::copy(face_pressure.begin(), face_pressure.end(), std::ostream_iterator<double>(f, "\n"));
            }

            // Compute and set fluxes of flow solution.
//             std::vector<double> cell_flux;
//             psolver_.faceFluxToCellFlux(face_flux, cell_flux);
//             const std::vector<int>& ncf = psolver_.numCellFaces();
            flow_solution_.faceflux_.assign(face_flux.begin(), face_flux.end());
        }

    private:
        /// A helper class for postProcessFluxes.
        class FaceFluxes
        {
        public:
            FaceFluxes(int sz)
                : fluxes_(sz, 0.0), visited_(sz, 0), max_modification_(0.0)
            {
            }
            void put(double flux, int f_ix) {
                ASSERT(visited_[f_ix] == 0 || visited_[f_ix] == 1);
                double sign = visited_[f_ix] ? -1.0 : 1.0;
                fluxes_[f_ix] += sign*flux;
                ++visited_[f_ix];
            }
            void get(double& flux, int f_ix) {
                ASSERT(visited_[f_ix] == 0 || visited_[f_ix] == 1);
                double sign = visited_[f_ix] ? -1.0 : 1.0;
                double new_flux = 0.5*sign*fluxes_[f_ix];
                double diff = std::fabs(flux - new_flux);
                max_modification_ = std::max(max_modification_, diff);
                flux = new_flux;
                ++visited_[f_ix];
            }
            void resetVisited()
            {
                std::fill(visited_.begin(), visited_.end(), 0);
            }

            double maxMod() const
            {
                return max_modification_;
            }
        private:
            std::vector<double> fluxes_;
            std::vector<int> visited_;
            double max_modification_;

        };

    public:
        /// @brief
        ///    Postprocess the solution fluxes.
        ///    This method modifies the solution object so that
        ///    out-fluxes of twin faces (that is, the two faces on a
        ///    cell-cell intersection) will be made antisymmetric.
        ///
        /// @return
        ///    The maximum modification made to the fluxes.
        double postProcessFluxes()
        {
            typedef typename GridInterface::CellIterator CI;
            typedef typename CI           ::FaceIterator FI;
            SparseTable<double>& cflux = flow_solution_.outflux_;

            FaceFluxes face_fluxes(pgrid_->numberOfFaces());
            // First pass: compute projected fluxes.
            for (CI c = pgrid_->cellbegin(); c != pgrid_->cellend(); ++c) {
                const int cell_index = c->index();
                for (FI f = c->facebegin(); f != c->faceend(); ++f) {
                    int f_ix = f->index();
                    double flux = cflux[cell_index][f->localIndex()];
                    if (f->boundary()) {
                        continue;
                    } else {
                        face_fluxes.put(flux, f_ix);
                    }
                }
            }
            face_fluxes.resetVisited();
            // Second pass: set all fluxes to the projected ones.
            for (CI c = pgrid_->cellbegin(); c != pgrid_->cellend(); ++c) {
                const int cell_index = c->index();
                for (FI f = c->facebegin(); f != c->faceend(); ++f) {
                    int f_ix = f->index();
                    double& flux = cflux[cell_index][f->localIndex()];
                    if (f->boundary()) {
                        continue;
                    } else {
                        face_fluxes.get(flux, f_ix);
                    }
                }
            }
            return face_fluxes.maxMod();
        }




        /// @brief
        ///    Type representing the solution to a given flow problem.
        class FlowSolution {
        public:
            friend class TpfaCompressible;

            /// @brief
            ///    The element type of the matrix representation of
            ///    the mimetic inner product.  Assumed to be a
            ///    floating point type, and usually, @code Scalar
            ///    @endcode is an alias for @code double @endcode.
            typedef double Scalar;

            /// @brief
            ///    Retrieve the current cell pressure in a given cell.
            ///
            /// @param [in] c
            ///    Cell for which to retrieve the current cell
            ///    pressure.
            ///
            /// @return
            ///    Current cell pressure in cell @code *c @endcode.
            Scalar cellPressure(const int cell) const
            {
                return pressure_[cell];
            }

            const std::vector<double>& cellPressure() const
            {
                return pressure_;
            }

            /// @brief
            ///    Retrieve current flux across given face in
            ///    direction of outward normal vector.
            ///
            /// @param [in] f
            ///    Face across which to retrieve the current signed
            ///    flux.
            ///
            /// @return
            ///    Current outward flux across face @code *f @endcode.
            Scalar faceFlux(const int face) const
            {
                return faceflux_[face];
            }

            const std::vector<double>& faceFlux() const
            {
                return faceflux_;
            }
        private:
            std::vector<Scalar> pressure_;
            std::vector<Scalar> faceflux_;

            void clear()
            {
                pressure_.clear();
                faceflux_.clear();
            }
        };





        /// @brief
        ///    Type representing the solution to the problem defined
        ///    by the parameters to @code solve() @endcode.  Always a
        ///    reference-to-const.  The @code SolutionType @endcode
        ///    exposes methods @code pressure(c) @endcode and @code
        ///    outflux(f) @endcode from which the cell pressure in
        ///    cell @code *c @endcode and outward-pointing flux across
        ///    interface @code *f @endcode may be recovered.
        typedef const FlowSolution& SolutionType;





        /// @brief
        ///    Recover the solution to the problem defined by the
        ///    parameters to method @code solve() @endcode.  This
        ///    solution is meaningless without a previous call to
        ///    method @code solve() @endcode.
        ///
        /// @return
        ///    The current solution.
        SolutionType getSolution()
        {
            return flow_solution_;
        }





    private:
        void computeFluidProps(const FluidInterface& fluid,
                               const std::vector<typename FluidInterface::PhaseVec>& phase_pressure,
                               const std::vector<typename FluidInterface::PhaseVec>& phase_pressure_face,
                               const std::vector<typename FluidInterface::CompVec>& z)
        {
            int num_cells = z.size();
            int num_faces = pgrid_->numFaces();
            ASSERT(num_cells == pgrid_->numCells());
            const int np = FluidInterface::numPhases;
            const int nc = FluidInterface::numComponents;
            BOOST_STATIC_ASSERT(np == nc);
            fp_.totcompr.resize(num_cells);
            fp_.voldiscr.resize(num_cells);
            fp_.cellA.resize(num_cells*nc*np);
            fp_.faceA.resize(num_faces*nc*np);
            fp_.phasemobf.resize(num_faces*np);
            fp_.phasemobc.resize(num_cells*np); // Just a helper
            typedef typename FluidInterface::PhaseVec PhaseVec;
            typedef typename FluidInterface::CompVec CompVec;
            PhaseVec mob;
            BOOST_STATIC_ASSERT(np == 3);
            for (int cell = 0; cell < num_cells; ++cell) {
                typename FluidInterface::FluidState state = fluid.computeState(phase_pressure[cell], z[cell]);
                fp_.totcompr[cell] = state.total_compressibility_;
                fp_.voldiscr[cell] = state.total_phase_volume_ - pgrid_->cellVolume(cell)*poro_[cell];
                std::copy(state.mobility_.begin(), state.mobility_.end(), fp_.phasemobc.begin() + cell*np);
                Dune::SharedFortranMatrix A(nc, np, state.phase_to_comp_);
                Dune::SharedFortranMatrix cA(nc, np, &fp_.cellA[cell*nc*np]);
                cA = A;
            }
            // Set phasemobf to average of cells' phase mobs, if pressures are equal, else use upwinding.
            // Set faceA by using average of cells' z and face pressures.
            for (int face = 0; face < num_faces; ++face) {
                int c[2] = { pgrid_->faceCell(face, 0), pgrid_->faceCell(face, 1) };
                PhaseVec phase_p[2];
                CompVec z_face(0.0);
                int num = 0;
                for (int j = 0; j < 2; ++j) {
                    if (c[j] >= 0) {
                        phase_p[j] = phase_pressure[c[j]];
                        z_face += z[c[j]];
                        ++num;
                    } else {
                        // Boundaries get essentially -inf pressure for upwinding purpose. \TODO handle BCs.
                        phase_p[j] = PhaseVec(-1e100);
                        // \TODO The two lines below are wrong for outflow faces.
                        z_face += inflow_mixture_;
                        ++num;
                    }
                }
                z_face /= double(num);
                for (int phase = 0; phase < np; ++phase) {
                    if (phase_p[0][phase] == phase_p[1][phase]) {
                        // Average mobilities.
                        double aver = 0.5*(fp_.phasemobc[np*c[0] + phase] + fp_.phasemobc[np*c[1] + phase]);
                        fp_.phasemobf[np*face + phase] = aver;
                    } else {
                        // Upwind mobilities.
                        int upwind = (phase_p[0][phase] > phase_p[1][phase]) ? 0 : 1;
                        fp_.phasemobf[np*face + phase] = fp_.phasemobc[np*c[upwind] + phase];
                    }
                }
                typename FluidInterface::FluidState face_state = fluid.computeState(phase_pressure_face[face], z_face);
                Dune::SharedFortranMatrix A(nc, np, face_state.phase_to_comp_);
                Dune::SharedFortranMatrix fA(nc, np, &fp_.faceA[face*nc*np]);
                fA = A;
            }
        }

        struct FluidProps
        {
            std::vector<double> totcompr;
            std::vector<double> voldiscr;
            std::vector<double> cellA;
            std::vector<double> faceA;
            std::vector<double> phasemobf;
            std::vector<double> phasemobc;
        };

        FluidProps fp_;
        const GridInterface* pgrid_;
        std::vector<double> poro_;
        PressureSolver psolver_;
        LinearSolverISTL linsolver_;
        FlowSolution flow_solution_;

        typename FluidInterface::CompVec inflow_mixture_;
        int num_iter_;
    };


} // namespace Dune



#endif // OPM_TPFACOMPRESSIBLE_HEADER_INCLUDED
