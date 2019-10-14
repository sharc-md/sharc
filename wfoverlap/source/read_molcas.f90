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
    INTEGER(KIND=ilong):: nao
    INTEGER(KIND=ilong):: nsym
    INTEGER(KIND=ilong):: irc
    INTEGER(KIND=ilong):: nvalid
    INTEGER(KIND=ilong):: i
    INTEGER(KIND=ilong):: j
    INTEGER(KIND=ilong):: k
    INTEGER*8:: iCode=6
    INTEGER*8:: ilOne=1
    INTEGER*8:: mc1
    INTEGER(KIND=ilong), DIMENSION(:),ALLOCATABLE :: nbas
    REAL*8, DIMENSION(:),ALLOCATABLE :: temp
    INTEGER:: allocstat

    CALL link_it()
    CALL getenvinit()
    CALL fioinit()
    CALL NameRun('RUNFILE')
    CALL mcget_iscalar('nsym',nsym)
    IF(nsym.GT.iOne)THEN
      WRITE(0,'("Found ",I2," symmetry groups - this is not a C1-symmetry calculation")'),nsym
      WRITE(6,'("Found ",I2," symmetry groups - this is not a C1-symmetry calculation")'),nsym
      STOP 1
    END IF
    ALLOCATE(nbas(1:nsym))
    CALL mcget_iarray('nbas',nbas,nsym)

    nao = nbas(1)
    
    ! Determine which part of the overlap matrix should be read
    IF(same_aos)THEN
      nl = nao
      nr = nao
    ELSE
      IF ((nl.EQ.-1).AND.(nr.EQ.-1)) THEN
        ! If neither nao_a nor nao_b were specified, assume that they are the same
        nl = nao / 2
        nr = nao / 2
      ELSE
        IF(nl+nr.NE.nao)THEN
          WRITE(6,'("Inconsistent number of AOs given: ", 3I10)') nl, nr, nao
          WRITE(6,*) "file: ",trim(adjustl(file))
          WRITE(0,'("Inconsistent number of AOs given: ", 3I10)') nl, nr, nao
          WRITE(0,*) "file: ",trim(adjustl(file))
          STOP 401
        END IF
      END IF
    END IF

    allocstat=myalloc(ovl,nl,nr,'+ ovl mat')
    IF(allocstat.NE.0) STOP

    allocstat=myalloc(temp,nao*(nao+1)/2)
    IF(allocstat.NE.0) STOP

    mc1=1
    CALL rdonexx(irc,iCode,'Mltpl  0',mc1,temp,ilOne,nvalid,trim(adjustl(file)))

    k=0
    IF (same_aos)THEN
      DO i=1,nbas(1)
        DO j=1,i
          k=k+1
          ovl(i,j)=temp(k)
          ovl(j,i)=temp(k)
        END DO ! j = 1,i
      END DO ! i = 1,nbas
    ELSE
      ! Read only the mixed block
      DO i=1,nbas(1)
        DO j=1,i
          k=k+1
          IF((i > nl).AND.(j <= nl)) ovl(j,i-nl)=temp(k)
        END DO ! j = 1,i
      END DO ! i = 1,nbas
    END IF

    IF(debug)THEN
      WRITE(6,*) "Overlap matrix read from Molcas ",trim(adjustl(file))
      WRITE(6,*)nl, nr
      DO i=1,nl
        DO j=1,nr
          WRITE(6,'(x,ES24.14E3)',ADVANCE='no')ovl(i,j)
        END DO
        WRITE(6,*)
      END DO
    END IF ! debug
    
    allocstat = mydealloc(temp)
    IF(allocstat.NE.0)THEN
      WRITE(6,*) "Problem when deallocating temp matrix in read_mc_overlap - ignoring"
      WRITE(0,*) "Problem when deallocating temp matrix in read_mc_overlap - ignoring"
    END IF

  END SUBROUTINE read_molcas_overlap

  !--------------------------------------------------------------------------------------------------------------------------------

END MODULE read_molcas
