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


#include <dune/common/param/ParameterGroup.hpp>

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

    struct LinearSolverResults
    {
        bool converged;
        int iterations;
        double reduction;
    };


    class LinearSolverISTL
    {
    public:
        virtual ~LinearSolverISTL()
        {
        }
        virtual void init(const parameter::ParameterGroup& param)
        {
        }
        virtual LinearSolverResults solve(int size, int nonzeros,
                                          const int* ia, const int* ja, const double* sa,
                                          const double* rhs, double* solution,
                                          double relative_residual_tolerance,
                                          int verbosity_level)
        {
            // Build ISTL structures from input.
            typedef FieldVector<double, 1   > VectorBlockType;
            typedef FieldMatrix<double, 1, 1> MatrixBlockType;
            typedef BCRSMatrix <MatrixBlockType>        Matrix;
            typedef BlockVector<VectorBlockType>        Vector;
            typedef MatrixAdapter<Matrix,Vector,Vector> Operator;
            // System matrix
            Matrix A(size, size, nonzeros, Matrix::row_wise);
            for (Matrix::CreateIterator row = A.createbegin(); row != A.createend(); ++row) {
                int ri = row.index();
                for (int i = ia[ri]; i < ia[ri + 1]; ++i) {
                    row.insert(ja[i]);
                }
            }
            for (int ri = 0; ri < size; ++ri) {
                for (int i = ia[ri]; i < ia[ri + 1]; ++i) {
                    A[ri][ja[i]] = sa[i];
                }
            }
            // System RHS
            Vector b(size);
            std::copy(rhs, rhs + size, b.begin());
            // System solution
            Vector x(size);
            x = 0.0;

            // Solve with AMG solver.
#define FIRST_DIAGONAL 1
#define SYMMETRIC 1
#define SMOOTHER_ILU 1
#define ANISOTROPIC_3D 0

#if FIRST_DIAGONAL
            typedef Amg::FirstDiagonal CouplingMetric;
#else
            typedef Amg::RowSum        CouplingMetric;
#endif

#if SYMMETRIC
            typedef Amg::SymmetricCriterion<Matrix,CouplingMetric>   CriterionBase;
#else
            typedef Amg::UnSymmetricCriterion<Matrix,CouplingMetric> CriterionBase;
#endif

#if SMOOTHER_ILU
            typedef SeqILU0<Matrix,Vector,Vector>        Smoother;
#else
            typedef SeqSSOR<Matrix,Vector,Vector>        Smoother;
#endif
            typedef Amg::CoarsenCriterion<CriterionBase> Criterion;
            typedef Amg::AMG<Operator,Vector,Smoother>   Precond;

            Operator opA(A);
            // Construct preconditioner.
            double relax = 1;
            Precond::SmootherArgs smootherArgs;
            smootherArgs.relaxationFactor = relax;

            Criterion criterion;
            criterion.setDebugLevel(verbosity_level);
#if ANISOTROPIC_3D
            criterion.setDefaultValuesAnisotropic(3, 2);
#endif
            Precond precond(opA, criterion, smootherArgs);
            CGSolver<Vector> linsolve(opA, precond, relative_residual_tolerance, size, verbosity_level);
            InverseOperatorResult result;
            linsolve.apply(x, b, result);
            std::copy(x.begin(), x.end(), solution);
            LinearSolverResults res;
            res.converged = result.converged;
            res.iterations = result.iterations;
            res.reduction = result.reduction;
            return res;
        }
    };


} // namespace Dune


#endif // OPM_LINEARSOLVERISTL_HEADER
