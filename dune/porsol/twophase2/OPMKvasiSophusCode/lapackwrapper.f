
C  Copyright 2011 - Magne Haveraaen, Helmer André Friis and Hans Munthe-Kaas.

C  This file is part of the Open Porous Media project (OPM).

C  OPM is free software: you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation, either version 3 of the License, or
C  (at your option) any later version.

C  OPM is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  You should have received a copy of the GNU General Public License
C  along with OPM.  If not, see <http://www.gnu.org/licenses/>.



      SUBROUTINE MATINV(A,N) 
      INTEGER N,M,INFO,LDA,LWORK
      DOUBLE PRECISION A(N,N),WORK(N*2,N*2)
      INTEGER IPIV(N)
      EXTERNAL DGETRF,DGETRI
      M=N
      LDA=N
      LWORK=N
      CALL DGETRF( M, M, A, LDA, IPIV, INFO )
      CALL DGETRI( M, A, LDA, IPIV, WORK, LWORK, INFO )
      END

      SUBROUTINE EQSOLVE(A,B,N,NRHS)
      INTEGER N,NRHS,LDA,LDB,INFO
      DOUBLE PRECISION A(N,N),B(N,NRHS)
      INTEGER IPIV(N)
      EXTERNAL DGESV
      LDA = N
      LDB = N
      CALL DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO)
      END

      SUBROUTINE EQSOLVB(AB,B,N,NRHS,KL,KU,LDAB)
      INTEGER N,NRHS,LDAB,LDB,INFO,KL,KU
      DOUBLE PRECISION AB(LDAB,N),B(N,NRHS)
      INTEGER IPIV(N)
      EXTERNAL DGBSV
      LDB = N
      CALL DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO)
      END
