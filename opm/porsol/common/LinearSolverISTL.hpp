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


#include <opm/porsol/common/AbstractLinearSolver.hpp>


namespace Opm
{

    class LinearSolverISTL : public AbstractLinearSolver
    {
    public:
        LinearSolverISTL();
        virtual ~LinearSolverISTL();
        virtual void init(const Opm::parameter::ParameterGroup& param);
        virtual LinearSolverResults solve(int size, int nonzeros,
                                          const int* ia, const int* ja, const double* sa,
                                          const double* rhs, double* solution);
    private:
        double linsolver_residual_tolerance_;
        int linsolver_verbosity_;
        enum LinsolverType { CG_ILU0 = 0, CG_AMG = 1, BiCGStab_ILU0 = 2 };
        LinsolverType linsolver_type_;
        bool linsolver_save_system_;
        std::string linsolver_save_filename_;
        int linsolver_max_iterations_;
    };


} // namespace Opm


#endif // OPM_LINEARSOLVERISTL_HEADER
