//===========================================================================
//
// File: blas_lapack.cpp
//
// Created: Sun Jun 21 18:56:51 2009
//
// Author(s): Bård Skaflestad     <bard.skaflestad@sintef.no>
//            Atgeirr F Rasmussen <atgeirr@sintef.no>
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

#include <opm/porsol/common/blas_lapack.hpp>

namespace Opm {
namespace BLAS_LAPACK {

template<>
void GEMV<double>(const char*   transA,
                  const int     m     , const int     n,
                  const double& a1    , const double* A, const int ldA,
                                        const double* x, const int incX,
                  const double& a2    ,       double* y, const int incY)
{
    assert((transA[0] == 'N') || (transA[0] == 'T'));

    DGEMV(F77_CHARACTER(transA[0]),
          &m, &n, &a1, A, &ldA, x, &incX, &a2, y, &incY);
}

template<>
void GEMM<double>(const char*   transA, const char*   transB,
                  const int     m     , const int     n     , const int k  ,
                  const double& a1    , const double* A     , const int ldA,
                                        const double* B     , const int ldB,
                  const double& a2    ,       double* C     , const int ldC)
{
    assert((transA[0] == 'N') || (transA[0] == 'T'));
    assert((transB[0] == 'N') || (transB[0] == 'T'));

    DGEMM(F77_CHARACTER(transA[0]), F77_CHARACTER(transB[0]),
          &m, &n, &k, &a1, A, &ldA, B, &ldB, &a2, C, &ldC);
}

template<>
void SYRK<double>(const char*   uplo, const char*   trans,
                  const int     n   , const int     k    ,
                  const double& a1  , const double* A    , const int ldA,
                  const double& a2  ,       double* C    , const int ldC)
{
    assert((uplo[0]  == 'U') || (uplo[0]  == 'L'));
    assert((trans[0] == 'N') || (trans[0] == 'T'));

    DSYRK(F77_CHARACTER(uplo[0]), F77_CHARACTER(trans[0]),
          &n, &k, &a1, A, &ldA, &a2, C, &ldC);
}

template<>
void TRMM<double>(const char*   side  , const char* uplo,
                  const char*   transA, const char* diag,
                  const int     m     , const int   n   , const double& a,
                  const double* A     , const int   ldA ,
                        double* B     , const int   ldB)
{
    assert((side[0]   == 'L') || (side[0]   == 'R'));
    assert((uplo[0]   == 'U') || (uplo[0]   == 'L'));
    assert((transA[0] == 'N') || (transA[0] == 'T'));
    assert((diag[0]   == 'N') || (diag[0]   == 'U'));

    DTRMM(F77_CHARACTER(side[0])  , F77_CHARACTER(uplo[0]),
          F77_CHARACTER(transA[0]), F77_CHARACTER(diag[0]),
          &m, &n, &a, A, &ldA, B, &ldB);
}

template<>
void GEQRF<double>(const int     m    , const int     n   ,
                         double* A    , const int     ld  ,
                         double* tau  ,       double* work,
                   const int     lwork,       int&    info)
{
    DGEQRF(&m, &n, A, &ld, tau, work, &lwork, &info);
}

template<>
void ORGQR<double>(const int     m   , const int n    , const int     k  ,
                         double* A   , const int ld   , const double* tau,
                         double* work, const int lwork,       int&    info)
{
    DORGQR(&m, &n, &k, A, &ld, tau, work, &lwork, &info);
}

template<>
void GETRF<double>(const int m, const int n , double* A,
                   const int ld, int* ipiv, int& info)
{
    DGETRF(&m, &n, A, &ld, ipiv, &info);
}

template<>
void GETRI(const int  n   , double* A   , const int ld,
           const int* ipiv, double* work, int lwork, int& info)
{
    DGETRI(&n, A, &ld, ipiv, work, &lwork, &info);
}

}
}
