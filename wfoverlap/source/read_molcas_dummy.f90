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

MODULE read_molcas
  USE sysparam
  USE my_alloc
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_molcas_overlap

!--------------------------------------------------------------------------------------------------------------------------------==
!--------------------------------------------------------------------------------------------------------------------------------==
CONTAINS

  !--------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_molcas_overlap(file,ovl,nl,nr,same_aos)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: file
    REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: ovl
    INTEGER(KIND=ilong), INTENT(OUT) :: nl
    INTEGER(KIND=ilong), INTENT(OUT) :: nr
    LOGICAL, INTENT(IN) :: same_aos
    nl=-1
    nr=-1
    WRITE(6,*)'This is the dummy file read_molcas_dummy.f90.'
    WRITE(6,*)'Compilation with MOLCAS support required!'
    STOP 3
  END SUBROUTINE read_molcas_overlap

  !--------------------------------------------------------------------------------------------------------------------------------

END MODULE read_molcas
