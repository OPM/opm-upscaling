//===========================================================================
//
// File: IfshInterface.hpp
//
// Created: Thu Sep 23 19:04:38 2010
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Bård Skaflestad     <bard.skaflestad@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2010 Statoil ASA.

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

#ifndef OPENRS_IFSHINTERFACE_HEADER
#define OPENRS_IFSHINTERFACE_HEADER

#include "../opmpressure/HybridPressureSolver.hpp"
#include <dune/common/ErrorMacros.hpp>
#include <dune/common/SparseTable.hpp>
#include <dune/porsol/common/LinearSolverISTL.hpp>


namespace Dune
{


    template <class GridInterface,
              class RockInterface,
              class BCInterface>
    class IfshInterface
    {
    public:

        /// @brief
        ///    Default constructor. Does nothing.
        IfshInterface()
            : pgrid_(0)
        {
        }


        /// @brief
        ///    All-in-one initialization routine.  Enumerates all grid
        ///    connections, allocates sufficient space, defines the
        ///    structure of the global system of linear equations for
        ///    the contact pressures, and computes the permeability
        ///    dependent inner products for all of the grid's cells.
        ///
        /// @param [in] g
        ///    The grid.
        ///
        /// @param [in] r
        ///    The reservoir properties of each grid cell.
        ///
        /// @param [in] grav
        ///    Gravity vector.  Its Euclidian two-norm value
        ///    represents the strength of the gravity field (in units
        ///    of m/s^2) while its direction is the direction of
        ///    gravity in the current model.
        ///
        /// @param [in] bc
        ///    The boundary conditions describing how the current flow
        ///    problem interacts with the outside world.  This is used
        ///    only for the purpose of introducing additional
        ///    couplings in the case of periodic boundary conditions.
        ///    The specific values of the boundary conditions are not
        ///    inspected in @code init() @endcode.
        template<class Point>
        void init(const GridInterface&      g,
                  const RockInterface&      r,
                  const Point&              grav,
                  const BCInterface&        bc)
        {
            pgrid_ = &g;
            // Extract perm tensors.
            const double* perm = &(r.permeability(0)(0,0));
            // Check that we only have noflow boundary conditions.
            for (int i = 0; i < bc.size(); ++i) {
                if (bc.flowCond(i).isPeriodic()) {
                    THROW("Periodic BCs are not handled yet by ifsh.");
                }
            }
            // Initialize ifsh
            ifsh_.init(g.grid(), perm, &grav[0]);
        }





        /// @brief
        ///    Construct and solve system of linear equations for the
        ///    pressure values on each interface/contact between
        ///    neighbouring grid cells.  Recover cell pressure and
        ///    interface fluxes.  Following a call to @code solve()
        ///    @encode, you may recover the flow solution from the
        ///    @code getSolution() @endcode method.
        ///
        /// @tparam FluidInterface
        ///    Type presenting an interface to fluid properties such
        ///    as density, mobility &c.  The type is expected to
        ///    provide methods @code phaseMobilities() @endcode and
        ///    @code phaseDensities() @endcode for phase mobility and
        ///    density in a single cell, respectively.
        ///
        /// @param [in] r
        ///    The fluid properties of each grid cell.  In method
        ///    @code solve() @endcode we query this object for the
        ///    phase mobilities (i.e., @code r.phaseMobilities()
        ///    @endcode) and the phase densities (i.e., @code
        ///    phaseDensities() @encode) of each phase.
        ///
        /// @param [in] sat
        ///    Saturation of primary phase.  One scalar value for each
        ///    grid cell.  This parameter currently limits @code
        ///    IncompFlowSolverHybrid @endcode to two-phase flow
        ///    problems.
        ///
        /// @param [in] bc
        ///    The boundary conditions describing how the current flow
        ///    problem interacts with the outside world.  Method @code
        ///    solve() @endcode inspects the actual values of the
        ///    boundary conditions whilst forming the system of linear
        ///    equations.
        ///
        ///    Specifically, the @code bc.flowCond(bid) @endcode
        ///    method is expected to yield a valid @code FlowBC
        ///    @endcode object for which the methods @code pressure()
        ///    @endcode, @code pressureDifference() @endcode, and
        ///    @code outflux() @endcode yield valid responses
        ///    depending on the type of the object.
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
        template<class FluidInterface>
        void solve(const FluidInterface&      r  ,
                   const std::vector<double>& sat,
                   const BCInterface&         bc ,
                   const std::vector<double>& src,
                   double residual_tolerance = 1e-8,
                   int linsolver_verbosity = 1,
                   int linsolver_type = 1,
                   bool same_matrix = false)
        {
            if (linsolver_type != 1) {
                MESSAGE("Requested linear solver type " << linsolver_type << ", but we only have AMG so far.");
            }
            if (same_matrix) {
                MESSAGE("Requested reuse of preconditioner, not implemented so far.");
            }

            int num_cells = sat.size();
            std::vector<double> totmob(num_cells, 1.0);
            std::vector<double> omega(num_cells, 0.0);
            for (int cell = 0; cell < num_cells; ++cell) {
                totmob[cell] = r.totalMobility(cell, sat[cell]);
                double f_w = r.fractionalFlow(cell, sat[cell]);
                omega[cell] = r.densityFirstPhase()*f_w + r.densitySecondPhase()*(1.0 - f_w);
            }
            int num_faces = pgrid_->numberOfFaces();
            std::vector<HybridPressureSolver::FlowBCTypes> bctypes(num_faces, HybridPressureSolver::FBC_UNSET);
            std::vector<double> bcvalues(num_faces, 0.0);
            for (int face = 0; face < num_faces; ++face) {
                int bid = pgrid_->grid().boundaryId(face);
                FlowBC face_bc = bc.flowCond(bid);
                if (face_bc.isDirichlet()) {
                    bctypes[face] = HybridPressureSolver::FBC_PRESSURE;
                    bcvalues[face] = face_bc.pressure();
                } else if (face_bc.isNeumann()) {
                    bctypes[face] = HybridPressureSolver::FBC_FLUX;
                    bcvalues[face] = face_bc.outflux(); // TODO: may have to switch sign here depending on orientation.
                } else {
                    THROW("Unhandled boundary condition type.");
                }
            }
//             typedef typename GridInterface::CellIterator CI;
//             typedef typename CI::FaceIterator FI;
//             for (CI c = pgrid_->cellbegin(); c != pgrid_->cellend(); ++c) {
//                 for (FI f = c->facebegin(); f != c-> faceend(); ++f) {
//                     if (f->boundary()) {
//                         int face = f->index();
//                         int bid = f->boundaryId();
//                         FlowBC face_bc = bc.flowCond(bid);
//                         if (face_bc.isDirichlet()) {
//                             bctypes[face] = HybridPressureSolver::FBC_PRESSURE;
//                             bcvalues[face] = face_bc.pressure();
//                         } else if (face_bc.isNeumann()) {
//                             bctypes[face] = HybridPressureSolver::FBC_FLUX;
//                             bcvalues[face] = face_bc.outflux(); // TODO: may have to switch sign here depending on orientation.
//                         } else {
//                             THROW("Unhandled boundary condition type.");
//                         }
//                     }
//                 }
//             }

            // Assemble system matrix and rhs.
            ifsh_.assemble(src, totmob, omega, bctypes, bcvalues);

//             static int count = 0;
//             ++count;
//             printSystem(std::string("linsys_mimetic-") + boost::lexical_cast<std::string>(count));

            HybridPressureSolver::LinearSystem s;
            ifsh_.linearSystem(s);
            LinearSolverResults res = linsolver_.solve(s.n, s.nnz, s.ia, s.ja, s.sa, s.b, s.x,
                                                       residual_tolerance, linsolver_verbosity);
            if (!res.converged) {
                THROW("Linear solver failed to converge in " << res.iterations << " iterations.\n"
                      << "Residual reduction achieved is " << res.reduction << '\n');
            }
            flow_solution_.clear();
            std::vector<double> face_flux;
            ifsh_.computePressuresAndFluxes(flow_solution_.pressure_, face_flux);
            std::vector<double> cell_flux;
            ifsh_.faceFluxToCellFlux(face_flux, cell_flux);
            const std::vector<int>& ncf = ifsh_.numCellFaces();
            flow_solution_.outflux_.assign(cell_flux.begin(), cell_flux.end(),
                                           ncf.begin(), ncf.end());
        }





        /// @brief
        ///    Type representing the solution to a given flow problem.
        class FlowSolution {
        public:
            friend class IfshInterface;

            /// @brief
            ///    The element type of the matrix representation of
            ///    the mimetic inner product.  Assumed to be a
            ///    floating point type, and usually, @code Scalar
            ///    @endcode is an alias for @code double @endcode.
            typedef typename GridInterface::Scalar       Scalar;

            /// @brief
            ///    Convenience alias for the grid interface's cell
            ///    iterator.
            typedef typename GridInterface::CellIterator CI;

            /// @brief
            ///    Convenience alias for the cell's face iterator.
            typedef typename CI           ::FaceIterator FI;

            /// @brief
            ///    Retrieve the current cell pressure in a given cell.
            ///
            /// @param [in] c
            ///    Cell for which to retrieve the current cell
            ///    pressure.
            ///
            /// @return
            ///    Current cell pressure in cell @code *c @endcode.
            Scalar pressure(const CI& c) const
            {
                return pressure_[c->index()];
            }

            /// @brief
            ///    Retrieve current flux across given face in
            ///    direction of outward normal vector.
            ///
            /// @param [in] f
            ///    Face across which to retrieve the current outward
            ///    flux.
            ///
            /// @return
            ///    Current outward flux across face @code *f @endcode.
            Scalar outflux (const FI& f) const
            {
                return outflux_[f->cellIndex()][f->localIndex()];
            }
        private:
            std::vector<Scalar> pressure_;
            SparseTable<Scalar> outflux_;

            void clear()
            {
                pressure_.clear();
                outflux_.clear();
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
        const GridInterface* pgrid_;
        HybridPressureSolver ifsh_;
        LinearSolverISTL linsolver_;
        FlowSolution flow_solution_;
    };


} // namespace Dune


#endif // OPENRS_IFSHINTERFACE_HEADER
