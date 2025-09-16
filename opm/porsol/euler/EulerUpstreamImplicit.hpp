//===========================================================================
//
// File: EulerUpstreamImplicit.hpp
//
// Created: Tue Jun 16 14:07:53 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            B�rd Skaflestad     <bard.skaflestad@sintef.no>
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

#ifndef OPENRS_EULERUPSTREAMIMPLICIT_HEADER
#define OPENRS_EULERUPSTREAMIMPLICIT_HEADER

#include <opm/common/utility/numeric/SparseVector.hpp>
#include <opm/common/utility/parameters/ParameterGroup.hpp>

#include <opm/porsol/common/ImplicitTransportDefs.hpp>
#include <opm/porsol/common/BCRSMatrixBlockAssembler.hpp>

#include <opm/grid/common/GridAdapter.hpp>

#include <opm/core/transport/implicit/ImplicitAssembly.hpp>
#include <opm/core/transport/implicit/ImplicitTransport.hpp>
#include <opm/core/transport/implicit/JacobianSystem.hpp>
#include <opm/core/transport/implicit/SinglePointUpwindTwoPhase.hpp>


namespace Opm {

    /// Class for doing simple transport by implicit Euler upstream method for general grid.
    /// @tparam
    template <class GridInterface, class ReservoirProperties, class BoundaryConditions>
    class EulerUpstreamImplicit
    {
    public:
        /// @brief
        /// @todo Doc me
        EulerUpstreamImplicit	();
        /// @brief
        /// @todo Doc me
        /// @param
        EulerUpstreamImplicit	(const GridInterface& grid,
                      const ReservoirProperties& resprop,
                      const BoundaryConditions& boundary);
        /// @brief
        /// @todo Doc me
        /// @param
        void init(const Opm::ParameterGroup& param);
        /// @brief
        /// @todo Doc me
        /// @param
        void init(const Opm::ParameterGroup& param,
                  const GridInterface& grid,
                  const ReservoirProperties& resprop,
                  const BoundaryConditions& boundary);
        /// @brief
        /// @todo Doc me
        /// @param
        void initObj(const GridInterface& grid,
                     const ReservoirProperties& resprop,
                     const BoundaryConditions& boundary);
        /// @brief
        /// @todo Doc me
        /// @param
        void display();

        /// @brief Solve transport equation.
        /// @param saturation the evolving saturation
        /// @param time Time in seconds.
        /// @param gravity Gravity
        /// @param pressure_sol Pressure solution
        /// @param injection_rates Injection ratees
        template <class PressureSolution>
        bool transportSolve(std::vector<double>& saturation,
                            const double time,
                            const typename GridInterface::Vector& gravity,
                            const PressureSolution& pressure_sol,
                            const Opm::SparseVector<double>& injection_rates) const;

    protected:
        // typedef typename GridInterface::CellIterator CIt;
        // typedef typename CIt::FaceIterator FIt;
        // typedef typename FIt::Vector Vector;

        template <class PressureSolution>
        void smallTimeStep(std::vector<double>& saturation,
                           const double time,
                           const typename GridInterface::Vector& gravity,
                           const PressureSolution& pressure_sol,
                           const Opm::SparseVector<double>& injection_rates) const;

        void checkAndPossiblyClampSat(std::vector<double>& s) const;



        // Interface to implicit solver
        typedef Opm::TwophaseFluidWrapper                          TwophaseFluid;
        typedef Opm::SinglePointUpwindTwoPhase<TwophaseFluid> TransportModel;
        typedef Dune::FieldVector<double, 1>    ScalarVectorBlockType;
        typedef Dune::FieldMatrix<double, 1, 1> ScalarMatrixBlockType;

        typedef Dune::BlockVector<ScalarVectorBlockType> ScalarBlockVector;
        typedef Dune::BCRSMatrix <ScalarMatrixBlockType> ScalarBCRSMatrix;

        typedef Opm::ImplicitTransportDefault::NewtonVectorCollection< ScalarBlockVector >          NVecColl;
        typedef Opm::ImplicitTransportDefault::JacobianSystem        < ScalarBCRSMatrix, NVecColl > JacSys;

        typedef Opm::LinearSolverBICGSTAB LinearSolver;
        typedef Opm::ImplicitTransport<TransportModel                              ,
                                       JacSys                                      ,
                                       Opm::MaxNormDune                            ,
                                       Opm::ImplicitTransportDefault::VectorNegater,
                                       Opm::ImplicitTransportDefault::VectorZero   ,
                                       Opm::ImplicitTransportDefault::MatrixZero   ,
                                       Opm::ImplicitTransportDefault::VectorAssign > TransportSolver;

        // should be initialized by param
        GridAdapter mygrid_;
        ReservoirProperties myrp_;

        bool check_sat_;
        bool clamp_sat_;
        int max_repeats_;
        std::vector<double> porevol_;

        std::vector<int> periodic_cells_;
        std::vector<int> periodic_faces_;
        std::vector<int> periodic_nbfaces_;
        std::vector<int> periodic_hfaces_;
        std::vector<int> direclet_cells_;
        std::vector<double> direclet_sat_;
        std::vector<double> direclet_hfaces_;
        std::vector<double> htrans_;
        Opm::ImplicitTransportDetails::NRControl ctrl_;
    };

} // namespace Opm

#include "EulerUpstreamImplicit_impl.hpp"

#endif // OPENRS_EULERUPSTREAM_HEADER
