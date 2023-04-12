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

MODULE outputmod
  !> module for different kinds of output
  USE sysparam
  USE global_storage
  USE calcmod
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: print_state_mat, print_renormalized_mat, print_ortho_mat, mixing_angles, print_dyson_orbs
!--------------------------------------------------------------------------------------------------------------------------------==
!--------------------------------------------------------------------------------------------------------------------------------==
CONTAINS

  !--------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE print_state_mat(A, print_debug, prupper)
    IMPLICIT NONE
    REAL(KIND=dop), DIMENSION(:,:), INTENT(IN) :: A
    LOGICAL, INTENT(IN) :: print_debug
    LOGICAL, INTENT(IN), OPTIONAL :: prupper ! print only the upper triangular part
    
    LOGICAL :: pru
    INTEGER(KIND=ilong) :: N1, N2, st1
    INTEGER(KIND=ilong) :: astate, bstate
    
    pru = .false.
    IF (PRESENT(prupper)) pru = prupper
    
    N1 = size(A,1)
    N2 = size(A,2)
    
    IF(print_debug)THEN
      DO astate=1,N1
        WRITE(6,'(999PES24.14E4)') (A(astate,:))
      END DO
    ELSE
      WRITE(6,'(16X)',ADVANCE='no')
      DO bstate=1,N2
        WRITE(6,'("|PsiB ",I2,">     ")',ADVANCE='no') bstate
      END DO
      WRITE(6,*)
      
      st1 = 1
      DO astate=1,N1
        WRITE(6,'("<PsiA ",I2,"|  ")',ADVANCE='no') astate

        IF (pru) THEN
          DO bstate = 1, astate
            WRITE(6, '(14X)', ADVANCE='no')
          END DO
          st1 = astate +1
        END IF

        WRITE(6,'(999F14.10)') (A(astate,bstate),bstate=st1,N2)
      END DO
    END IF
    
  END SUBROUTINE print_state_mat

  !--------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE print_renormalized_mat(A, L_CIcoefs, R_CIcoefs, L_nstate, R_nstate, do_sqrt, print_debug)
    IMPLICIT NONE
    REAL(KIND=dop), DIMENSION(:,:), INTENT(IN) :: A
    REAL(KIND=dop), DIMENSION(:,:), INTENT(IN) :: L_CIcoefs
    REAL(KIND=dop), DIMENSION(:,:), INTENT(IN) :: R_CIcoefs
    INTEGER(KIND=ilong), INTENT(IN) :: L_nstate
    INTEGER(KIND=ilong), INTENT(IN) :: R_nstate
    LOGICAL, INTENT(IN) :: do_sqrt
    LOGICAL, INTENT(IN) :: print_debug
    
    REAL(KIND=dop), DIMENSION(L_nstate) :: L_invCInorm
    REAL(KIND=dop) :: R_invCInorm
    REAL(KIND=dop), DIMENSION(L_nstate, R_nstate) :: A_ren
    
    INTEGER(KIND=ilong) L_state, R_state
    
    DO L_state = 1, L_nstate
      L_invCInorm(L_state) = 1.D0 / sum(L_CIcoefs(L_state,:)**2)
      IF (do_sqrt) L_invCInorm(L_state) = sqrt(L_invCInorm(L_state))
    END DO

    DO R_state = 1, R_nstate
      R_invCInorm = 1.D0 / sum(R_CIcoefs(R_state,:)**2)
      IF (do_sqrt) R_invCInorm = sqrt(R_invCInorm)
      DO L_state = 1, L_nstate
        A_ren(L_state, R_state) = A(L_state, R_state) * R_invCInorm * L_invCInorm(L_state)
      END DO
    END DO
    
    CALL print_state_mat(A_ren, print_debug)
  
  END SUBROUTINE print_renormalized_mat
  
  !--------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE print_ortho_mat(A, Aortho, print_debug)
    IMPLICIT NONE
    REAL(KIND=dop), DIMENSION(:,:), INTENT(IN) :: A
    REAL(KIND=dop), DIMENSION(:,:), INTENT(IN) :: Aortho
    LOGICAL, INTENT(IN) :: print_debug
    
    REAL(KIND=dop) :: nA, nAortho, scalar, phi
  
    CALL print_state_mat(Aortho, print_debug)
    
    nA      = SQRT(SUM(A(:,:)**2))
    nAortho = SQRT(SUM(Aortho(:,:)**2))
    scalar  = SUM(A(:,:)*Aortho(:,:))
    phi     = ACOS(scalar / (nA * nAortho))
    
    WRITE(6,'("  Angle to original matrix",F14.10)') phi
    IF (phi .GT. 0.3) THEN
      WRITE(6,*) 'WARNING: Orthogonalized matrix differs significantly from original matrix!'
      WRITE(6,*) '   There is probably mixing with external states.'
    END IF
  END SUBROUTINE print_ortho_mat
  
  !--------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE mixing_angles(A, print_debug)
  ! Compute mixing angles as the matrix logarithm of the orthogonal matrix A
    IMPLICIT NONE
    REAL(KIND=dop), DIMENSION(:,:), INTENT(IN) :: A
    LOGICAL, INTENT(IN) :: print_debug
    
    REAL(KIND=dop), DIMENSION(size(A,1),size(A,2)) :: logA
        
    WRITE(6,*)
    WRITE(6,*)'Computing mixing angles as matrix logarithm...'
    CALL matrix_log(A, logA, print_debug)
    
    WRITE(6,*)
    WRITE(6,*) 'Mixing angles'
    CALL print_state_mat(logA, print_debug, .true.)
  
  END SUBROUTINE
  
  !--------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE matrix_log(A, logA, print_debug)
  ! Compute the real matrix logarithm according to
  !   R. Shepard, S. R. Brozell, G. Gidofalvi J. Phys. Chem. A 2015, 119, 7924-7939.
    IMPLICIT NONE
    REAL(KIND=dop), DIMENSION(:,:), INTENT(IN) :: A
    REAL(KIND=dop), DIMENSION(:,:), INTENT(OUT) :: logA
    LOGICAL, INTENT(IN) :: print_debug    
    
    ! LAPACK args
    REAL(KIND=dop), DIMENSION(size(A,1),size(A,2)) :: T
    REAL(KIND=dop), DIMENSION(size(A,1),size(A,2)) :: V
    REAL(KIND=dop), DIMENSION(size(A,1),size(A,2)) :: tempmat
    REAL(KIND=dop), DIMENSION(size(A,1)) :: WR, WI
    INTEGER (KIND=ilong):: lwork, info
    REAL (KIND=dop), ALLOCATABLE, DIMENSION(:) :: work
    LOGICAL, DIMENSION(size(A,1)) :: bwork
    INTEGER (KIND=ilong) :: sdim
    
    REAL(KIND=dop), PARAMETER :: eps=1.E-10_dop
    INTEGER(KIND=ilong) :: N, i, j
    REAL(KIND=dop) :: det, temp, diag, offd
    LOGICAL :: skip
        
    N = size(A,1)
        
    det = determinant(A, N)
    WRITE(6,'("  Determinant of the original matrix:    ", F14.10)') det

    T = A
    ! Regularize the matrix to make all diagonal elements positive
    DO i = 1, N
      IF (T(i, i).lt.dZero) THEN
        DO j = 1, N
          T(j, i) = -T(j, i)
        END DO
      END IF
    END DO
    
    det = determinant(T, N)
    WRITE(6,'("  Determinant of the regularized matrix: ", F14.10)') det
    
    IF (ABS(dOne-det).le.eps) THEN
      WRITE(6,*)'  -> Using this matrix'
    ELSE IF (ABS(dOne+det).le.eps) THEN
      WRITE(6,*)'  -> Switching sign of the last column to'
      WRITE(6,*)'     create a positive definite matrix.'
      DO i = 1, N
        T(i, N) = -T(i, N)
      END DO
      IF (print_debug) THEN
        CALL print_state_mat(T, .false.)
        WRITE(6,'("  Determinant: ", F0.10)') determinant(T, N)
        
        CALL dgemm('N', 'T', N, N, N, dOne, T, N, T, N, dZero, tempmat, N)
        WRITE(6,*)'T.T^T'
        CALL print_state_mat(tempmat, .false.)
        
        CALL dgemm('T', 'N', N, N, N, dOne, T, N, T, N, dZero, tempmat, N)
        WRITE(6,*)'T^T.T'
        CALL print_state_mat(tempmat, .false.)
      END IF
    ELSE
      WRITE(6,*)'ERROR: Matrix not unitary!'
      WRITE(0,*)'ERROR: Matrix not unitary!'
      STOP 1
    END IF
    
    ! Schur decomposition
    sdim = 0
    ! workspace query
    lwork = -1
    ALLOCATE(work(1))
    CALL dgees('V', 'N', .true., N, T, N, sdim, WR, WI, V, N, work, lwork, bwork, info)
    lwork = work(1)
    DEALLOCATE(work)
    !WRITE(6,*)'lwork:', lwork
    
    ALLOCATE(work(lwork))
    CALL dgees('V', 'N', .true., N, T, N, sdim, WR, WI, V, N, work, lwork, bwork, info)
    DEALLOCATE(work)
    IF (info.NE.0)THEN
        WRITE(0,202)info
        WRITE(6,202)info
         STOP 202
    END IF
202 FORMAT("dgees failed, info = ", I12)

    ! store the diagonalized matrix log in logA
    logA = dZero
    skip = .false.
    DO i = 1, N-1
      IF (skip) THEN
        skip = .false.
        CYCLE
      END IF
      
      diag = T(i,i)
      
      ! log(1) = 0
      IF (abs(diag-dOne).le.eps) CYCLE
      
      ! Analyse a matrix of the form ((cos phi, sin phi), (-sin phi, cos phi))
      offd = T(i,i+1)
      
      IF (ABS(diag*diag + offd*offd - 1) .gt. eps) THEN
        WRITE(6,203) i, diag, offd, diag*diag + offd*offd
        WRITE(0,203) i, diag, offd, diag*diag + offd*offd
        STOP 1
203 FORMAT("Problem for matrix log: i = ", I0, ", T(i,i) = ", F14.10, ", T(i,i+1) = ", F14.10, ", squared sum:", F14.10)
      END IF
      
      temp = ATAN2(offd, diag)
      logA(i,i+1) =  temp
      logA(i+1,i) = -temp
      
      skip = .true.
    END DO

    IF (print_debug) THEN    
      WRITE(6,*)'Schur matrix'
      CALL print_state_mat(T, .false.)
    
      WRITE(6,*)'asin'
      CALL print_state_mat(logA, .false.)
      
      WRITE(6,*)'V'
      CALL print_state_mat(V, .false.)
    END IF
    
    ! tempmat = V.logA    
    CALL dgemm('N', 'N', N, N, N, dOne, V,       N, logA, N, dZero, tempmat, N)
    ! logA = tempmat.V^T
    CALL dgemm('N', 'T', N, N, N, dOne, tempmat, N,    V, N, dZero,    logA, N)
    
  END SUBROUTINE matrix_log
  !--------------------------------------------------------------------------------------------------------------------------------
  
  SUBROUTINE print_dyson_orbs(transp, R_nMO, dyson_orb_mobas)
    LOGICAL, INTENT(IN) :: transp
    INTEGER(KIND=ilong), INTENT(IN) :: R_nMO
    REAL(KIND=dop), DIMENSION(:,:,:), INTENT(IN) :: dyson_orb_mobas
    
    IF (MOprint.ge.1) CALL print_dyson_orbs_stdout(transp, R_nMO, dyson_orb_mobas)
    
    IF (MOprint.ge.2) CALL print_dyson_orbs_jmol(transp, R_nMO, dyson_orb_mobas)
  END SUBROUTINE print_dyson_orbs

  !--------------------------------------------------------------------------------------------------------------------------------
  
  SUBROUTINE print_dyson_orbs_stdout(transp, R_nMO, dyson_orb_mobas)
  ! Print the Dyson orbital coefficients to standard output
    LOGICAL, INTENT(IN) :: transp
    INTEGER(KIND=ilong), INTENT(IN) :: R_nMO
    REAL(KIND=dop), DIMENSION(:,:,:), INTENT(IN) :: dyson_orb_mobas

    INTEGER(KIND=ilong) :: MO
    INTEGER(KIND=ilong) :: astate, bstate
    
      DO astate=1,bra_nstate
        WRITE(6,'("<PsiA ",I2,"|  ")'),astate
        WRITE(6,'(6x)',ADVANCE='no')
        DO bstate=1,ket_nstate-1
          WRITE(6,'("   |PsiB ",I3,">   ")',ADVANCE='no') bstate
        END DO
        WRITE(6,'("   |PsiB ",I3,">   ")',ADVANCE='yes') ket_nstate
        DO MO=1,R_nMO
          WRITE(6,'("MO ",I5)',ADVANCE='no')MO
          IF(transp)THEN
            DO bstate=1,ket_nstate
              WRITE(6,'(x,1PE15.8)',ADVANCE='no')dyson_orb_mobas(MO,bstate,astate)
            END DO
          ELSE
            DO bstate=1,ket_nstate
              WRITE(6,'(x,1PE15.8)',ADVANCE='no')dyson_orb_mobas(MO,astate,bstate)
            END DO
          END IF
          WRITE(6,*)
        END DO
        WRITE(6,*)'-------------------------------------'
      END DO    
  
  END SUBROUTINE print_dyson_orbs_stdout

  !--------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE print_dyson_orbs_jmol(transp, R_nMO, dyson_orb_mobas)
    LOGICAL, INTENT(IN) :: transp
    INTEGER(KIND=ilong), INTENT(IN) :: R_nMO
    REAL(KIND=dop), DIMENSION(:,:,:), INTENT(IN) :: dyson_orb_mobas
    
    REAL(KIND=dop), PARAMETER :: minprt=0.005
    REAL(KIND=dop), PARAMETER :: minnorm=0.05
    
    INTEGER(KIND=ilong) :: MO
    INTEGER(KIND=ilong) :: astate, bstate
    INTEGER(KIND=ilong) :: jout, iostat
    REAL(KIND=dop), DIMENSION(R_nMO) :: cMO
    REAL(KIND=dop) :: dyson_norm
    
    jout = freeunit()
    OPEN(UNIT=jout, FILE='dyson_jmol.spt', IOSTAT=iostat)    
    
    WRITE(jout,*) 'background white'
    WRITE(jout,*) 'mo fill'
    WRITE(jout,*) 'mo cutoff 0.04'
    WRITE(jout,*) 'mo titleformat ""'
    WRITE(jout,*)
    
    DO astate=1,bra_nstate
      DO bstate=1,ket_nstate
          IF(transp)THEN
            cMO = dyson_orb_mobas(:,bstate,astate)
          ELSE
            cMO = dyson_orb_mobas(:,astate,bstate)
          END IF
          dyson_norm = sum(cMO(:)**2)          
          IF(dyson_norm.lt.minnorm) CYCLE
          
          cMO = cMO / sqrt(dyson_norm)
          WRITE(jout, '(A)', ADVANCE='no') "mo [ "          
          DO MO=1,R_nMO          
            IF(cMO(MO)*cMO(MO).ge.minprt) WRITE(jout, 401, ADVANCE='no') cMO(MO), MO
          END DO
        WRITE(jout, '(A)', ADVANCE='yes') "]"
        WRITE(jout,400) astate, bstate, dyson_norm
      END DO
    END DO
    CLOSE(jout)
    
400 FORMAT("write image png 'dyson_",I3.3,"-",I3.3,"_",F5.3,".png'")
401 FORMAT(F6.3,X,I4,X)
  END SUBROUTINE print_dyson_orbs_jmol

  !--------------------------------------------------------------------------------------------------------------------------------

END MODULE outputmod