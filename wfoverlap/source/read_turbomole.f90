!******************************************
!
!    SHARC Program Suite
!
!    Copyright (c) 2019 University of Vienna
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

MODULE read_turbomole
  USE sysparam
  USE my_alloc
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_turbomole_mos
  
!----------------------------------------------------------------------------------------------------------------------------------
! Read orbitals in MOLCAS lumorb format.
! This is a simple ASCII parser that makes no reference to any MOLCAS objects.
! Therefore, it is separate from read_molcas.f90
!----------------------------------------------------------------------------------------------------------------------------------
CONTAINS
  
  !--------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_turbomole_mos(mofile,coeffs,nmo,nao)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: mofile
    REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: coeffs
    INTEGER(KIND=ilong), INTENT(IN) :: nmo
    INTEGER(KIND=ilong), INTENT(INOUT) :: nao

    INTEGER(KIND=ilong):: nmo_read
    INTEGER:: allocstat
    INTEGER:: moflio
    INTEGER:: iost

    moflio = freeunit()
    OPEN(UNIT=moflio,FILE=trim(adjustl(mofile)),IOSTAT=iost,STATUS='old',ACTION='read')
    IF(iost .NE. 0) THEN
      WRITE(0,*) "read_turbomole_moheader: Cannot open file ",trim(adjustl(mofile))," to read"
      WRITE(6,*) "read_turbomole_moheader: Cannot open file ",trim(adjustl(mofile))," to read"
      STOP 1
    END IF
    
    CALL read_turbomole_moheader(moflio, nmo_read, nao)
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
    
    CALL read_turbomole_mocoefs(moflio,coeffs,nmo,nao)
    
    CLOSE(moflio)
  END SUBROUTINE read_turbomole_mos
  
  !--------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_turbomole_moheader(moflio, nmo, nao)
    IMPLICIT NONE
    INTEGER :: moflio
    INTEGER(KIND=ilong), INTENT(OUT) :: nmo
    INTEGER(KIND=ilong), INTENT(INOUT) :: nao
    INTEGER(KIND=ilong) :: naotmp

    CHARACTER*200:: line

    naotmp=0
    nmo=0
    DO
      READ(moflio,'(A200)') line
      IF (line(16:25)=='eigenvalue') THEN
        READ(line(1:6),*) nmo
        READ(line(56:60),*) naotmp
      END IF
      IF (line(1:4)=='$end') EXIT
    END DO
    IF(naotmp.NE.nao)THEN
      WRITE(0,105) naotmp, nao
      WRITE(6,105) naotmp, nao
      STOP 1
    END IF
105 FORMAT ("read_turbomole_moheader: MOcoef file states different number of atomic orbitals &
    &than AOoverlap file. naotmp = "    , I5, ", nao = ", I5)

  END SUBROUTINE read_turbomole_moheader

  !--------------------------------------------------------------------------------------------------------------------------------
  
  SUBROUTINE read_turbomole_mocoefs(moflio, coeffs, nmo, nao)
    IMPLICIT NONE
    INTEGER :: moflio
    REAL(KIND=dop), DIMENSION(nmo,nao), INTENT(OUT) :: coeffs
    INTEGER(KIND=ilong), INTENT(IN) :: nmo
    INTEGER(KIND=ilong), INTENT(IN) :: nao

    INTEGER :: imo,iao,iost
    CHARACTER*200:: line

    REWIND(moflio)
    ! skip some lines
    READ(moflio,*)
    DO
      READ(moflio,*) line
      IF (line(1:1)/='#') EXIT
    END DO
    BACKSPACE(moflio)

    DO imo=1,nmo
      READ(moflio,*) ! * eigenvalue
      DO iao=1,nao
        READ(moflio,'(E20.14)',IOSTAT=iost,ADVANCE='no') coeffs(imo,iao)
        IF(iost.NE.0)THEN
          WRITE(0,*) "I/O error in read_turbomole_mocoefs"
          WRITE(6,*) "I/O error in read_turbomole_mocoefs"
          STOP 1
        END IF
        IF (modulo(iao,4).EQ.0) READ(moflio,*)
      END DO
      IF (modulo(nao,4).NE.0)READ(moflio,*)
    END DO
  END SUBROUTINE read_turbomole_mocoefs
  
  !--------------------------------------------------------------------------------------------------------------------------------

END MODULE read_turbomole