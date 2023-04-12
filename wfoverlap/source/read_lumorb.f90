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

MODULE read_lumorb
  USE sysparam
  USE my_alloc
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_molcas_lumorb
  
!----------------------------------------------------------------------------------------------------------------------------------
! Read orbitals in MOLCAS lumorb format.
! This is a simple ASCII parser that makes no reference to any MOLCAS objects.
! Therefore, it is separate from read_molcas.f90
!----------------------------------------------------------------------------------------------------------------------------------
CONTAINS
  
  !--------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_molcas_lumorb(mofile,coeffs,nmo,nao)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: mofile
    REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: coeffs
    INTEGER(KIND=ilong), INTENT(IN) :: nmo
    INTEGER(KIND=ilong), INTENT(INOUT) :: nao
    REAL(KIND=dop) :: ver
    
    INTEGER(KIND=ilong):: nmo_read
    INTEGER:: allocstat
    INTEGER:: moflio
    INTEGER:: iost

    moflio = freeunit()
    OPEN(UNIT=moflio,FILE=trim(adjustl(mofile)),IOSTAT=iost,STATUS='old',ACTION='read')
    IF(iost .NE. 0) THEN
      WRITE(0,*) "read_lumorb_moheader: Cannot open file ",trim(adjustl(mofile))," to read"
      WRITE(6,*) "read_lumorb_moheader: Cannot open file ",trim(adjustl(mofile))," to read"
      STOP 1
    END IF
    
    CALL read_lumorb_moheader(moflio, nmo_read, nao, ver)
    IF(nmo_read.NE.nmo)THEN
      WRITE(6,'("Number of MOs in mocoef file does not agree with CIvec: ",I6," instead of ",I6)')nmo_read,nmo
      WRITE(6,*) "file: ",trim(adjustl(mofile))
      WRITE(0,'("Number of MOs in mocoef file does not agree with CIvec: ",I6," instead of ",I6)')nmo_read,nmo
      WRITE(0,*) "file: ",trim(adjustl(mofile))
      STOP 1
    END IF
    
    allocstat=myalloc(coeffs,nmo,nao,"+ mocoef")
    IF(allocstat.NE.0)THEN
      WRITE(6,'("Could not allocate mocoef; error:",I4)')allocstat
      WRITE(0,'("Could not allocate mocoef; error:",I4)')allocstat
      STOP 1
    END IF
    
    CALL read_lumorb_mos(moflio,coeffs,nmo,nao,ver)
    
    CLOSE(moflio)
  END SUBROUTINE read_molcas_lumorb
  
  !--------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_lumorb_moheader(moflio, nmo, nao, ver)
    IMPLICIT NONE
    INTEGER :: moflio
    INTEGER(KIND=ilong), INTENT(OUT) :: nmo
    INTEGER(KIND=ilong), INTENT(INOUT) :: nao
    REAL(KIND=dop) :: ver
    
    INTEGER :: n1, nirrep, n2
    CHARACTER*10 :: tmpstr
    
    READ(moflio,*) tmpstr, ver ! #INPORB x.y
    READ(moflio,*) ! #INFO
    READ(moflio,*) ! * comment
    READ(moflio,*) n1, nirrep, n2

    WRITE(6,'("INPORB version:", F5.2)') ver
    
    IF (nirrep.NE.1) THEN
      WRITE(6,'("No support for symmetry! # irreps: ", I4)')nirrep
      WRITE(0,'("No support for symmetry! # irreps: ", I4)')nirrep
      STOP 1
    END IF
    
    READ(moflio,*) nao
    READ(moflio,*) nmo
    
  END SUBROUTINE read_lumorb_moheader

  !--------------------------------------------------------------------------------------------------------------------------------
  
  SUBROUTINE read_lumorb_mos(moflio, coeffs, nmo, nao, ver)
    IMPLICIT NONE
    INTEGER :: moflio
    REAL(KIND=dop), DIMENSION(nmo,nao), INTENT(OUT) :: coeffs
    INTEGER(KIND=ilong), INTENT(IN) :: nmo
    INTEGER(KIND=ilong), INTENT(IN) :: nao
    REAL(KIND=dop) :: ver
  
    INTEGER:: iost
    INTEGER :: isym
    INTEGER :: imo
    INTEGER :: iao
    INTEGER :: lind
    CHARACTER*200:: line
    
    REAL(KIND=dop):: testvar
    
    ! find #ORB
    DO
      READ(moflio,'(A200)',IOSTAT=iost) line
      IF(iost.NE.0)THEN
        WRITE(0,*) "read_lumorb: did not find #ORB"
        WRITE(6,*) "read_lumorb: did not find #ORB"
        STOP 1
      END IF
      lind = index(line,'#ORB')
      IF (lind.NE.0) EXIT
    END DO
    
    DO imo=1,nmo
      READ(moflio,*) ! * ORBITAL
      IF ((ver>1.05).AND.(ver<1.15)) THEN
        DO iao=1,nao
          READ(moflio,'(E18.11)',IOSTAT=iost,ADVANCE='no')coeffs(imo,iao)
          IF(iost.NE.0)THEN
            WRITE(0,*) "I/O error in read_lumorb_mos. imo, iao = ", imo, iao
            WRITE(6,*) "I/O error in read_lumorb_mos. imo, iao = ", imo, iao
            STOP 1
          END IF
          IF (modulo(iao,4).EQ.0) READ(moflio,*)
        END DO
        IF (modulo(nao,4).NE.0)READ(moflio,*)
      ELSE IF ((ver>1.95).AND.(ver<2.25)) THEN
        DO iao=1,nao
          READ(moflio,'(E22.14)',IOSTAT=iost,ADVANCE='no')coeffs(imo,iao)
          IF(iost.NE.0)THEN
            WRITE(0,*) "I/O error in read_lumorb_mos. imo, iao = ", imo, iao
            WRITE(6,*) "I/O error in read_lumorb_mos. imo, iao = ", imo, iao
            STOP 1
          END IF
          IF (modulo(iao,5).EQ.0) READ(moflio,*)
        END DO
        IF (modulo(nao,5).NE.0)READ(moflio,*)
      ELSE
        WRITE(0,'("INPORB version not supported:", F5.2)') ver
        WRITE(6,'("INPORB version not supported:", F5.2)') ver
        STOP 1
      END IF
    END DO
  END SUBROUTINE read_lumorb_mos
  
  !--------------------------------------------------------------------------------------------------------------------------------

END MODULE read_lumorb