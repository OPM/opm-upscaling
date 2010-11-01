//===========================================================================
//
// File: LinearSolverISTL.hpp
//
// Created: Mon Sep 27 10:07:04 2010
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

  This file is part of The Open Porous Media project (OPM).

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

#ifndef OPM_LINEARSOLVERISTL_HEADER
#define OPM_LINEARSOLVERISTL_HEADER


#include <dune/porsol/common/AbstractLinearSolver.hpp>

// TODO: clean up includes.
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/io.hh>

#include <dune/istl/overlappingschwarz.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>



namespace Dune
{

    class LinearSolverISTL : public AbstractLinearSolver
    {
    public:
        virtual ~LinearSolverISTL();
        virtual void init(const parameter::ParameterGroup& param);
        virtual LinearSolverResults solve(int size, int nonzeros,
                                          const int* ia, const int* ja, const double* sa,
                                          const double* rhs, double* solution,
                                          double relative_residual_tolerance,
                                          int verbosity_level);
    };


} // namespace Dune


#endif // OPM_LINEARSOLVERISTL_HEADER
