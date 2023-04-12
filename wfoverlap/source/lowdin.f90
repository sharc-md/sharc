!******************************************
!
!    SHARC Program Suite
!
!    Copyright (c) 2023 University of Vienna
!
!    This file is part of SHARC.
!
!    SHARC is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHARC is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    inside the SHARC manual.  If not, see <http://www.gnu.org/licenses/>.
!
!******************************************

MODULE lowdin_mod
IMPLICIT NONE

INTERFACE
  SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
    CHARACTER          JOBU, JOBVT
    INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
    DOUBLE PRECISION   A( LDA, * ), S( * ), U( LDU, * ), VT( LDVT, * ), WORK( * )
  ENDSUBROUTINE

ENDINTERFACE

 CONTAINS


! Lowdin orthogonalization of the overlap matrix
! The idea is taken from M. Persico and G. Granucci
!    see G. Granucci, M. Persico, A. Toniolo, J. Chem. Phys. 114 (2001), 10608-10615
! To do this, a singular value decomposition of the matrix is used, see
!    www.wou.edu/~beavers/Talks/LowdinJointMeetings0107.pdf
!    en.wikipedia.org/wiki/Singular_value_decomposition#Applications_of_the_SVD
!
! The orthogonalized matrix is computed as
!    S_ortho = U.V^T
! This is equivalent to
!    S_ortho = S.V.Lambda^(-1).V^T
!    used by Granucci et al., since
!    S = U.Lambda.V^T

SUBROUTINE lowdin(nstat, S, S_ortho)
    USE sysparam
    IMPLICIT NONE

    INTEGER (KIND=ilong), INTENT(IN) :: nstat
    REAL (KIND=dop), DIMENSION(nstat, nstat), INTENT(IN)    :: S
    REAL (KIND=dop), DIMENSION(nstat, nstat), INTENT(INOUT) :: S_ortho

    INTEGER (KIND=ilong):: info
    REAL (KIND=dop), DIMENSION(nstat, nstat) :: Sint
    REAL (KIND=dop), DIMENSION(nstat, nstat) :: U
    REAL (KIND=dop), DIMENSION(nstat, nstat) :: Vt
    REAL (KIND=dop), DIMENSION(nstat) :: lam
    INTEGER (KIND=ilong):: lwork
    REAL (KIND=dop), ALLOCATABLE, DIMENSION(:) :: work

    Sint = S
    
    ! do a workspace query to dgesvd
    lwork=-1
    ALLOCATE(work(1))
    CALL dgesvd('A', 'A', nstat, nstat, Sint, nstat, lam, U, nstat, Vt, nstat, work, lwork, info)
    lwork=work(1)
    DEALLOCATE(work)
!     write(*,*) 'Allocating for lapack:',lwork

    ! allocate the necessary workspace
    ALLOCATE(work(lwork))
    CALL dgesvd('A', 'A', nstat, nstat, Sint, nstat, lam, U, nstat, Vt, nstat, work, lwork, info)
    DEALLOCATE(work)
    IF (info.NE.0)THEN
        WRITE(0,201)info
        WRITE(6,201)info
!         STOP 201              ! dgesvd does not work currently with lapack/blas, only with MKL
    END IF
201 FORMAT("dgesvd failed, info = ", I12)

    IF (info.NE.0)THEN
      S_ortho=Sint
    ELSE
      CALL dgemm('N', 'N', nstat, nstat, nstat, dOne, U, nstat, Vt, nstat, dZero, S_ortho, nstat)
    ENDIF

END SUBROUTINE lowdin

ENDMODULE lowdin_mod
