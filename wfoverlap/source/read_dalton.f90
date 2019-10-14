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

MODULE read_dalton
  USE sysparam
  USE my_alloc

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_dalton_overlap

CONTAINS

  !-------------------------------------------------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_dalton_overlap(file,ovl,nl,nr,same_aos)
    ! this is mostly copied from $COLUMBUS/source/sif/iwfmt.f and adapted to the current environment
    !  there are probably still some unnecessary variables declared.
    !
    ! Use dalton_double-mol.py of this distribution for setting up the files
    !
    ! Or use the following steps:
    !   1. Use appropriate daltcomm file (containing ".NOTWO")
    !   2. copy daltaoin from GEO1/WORK
    !   3. double the second number in the fourth line
    !   4. tail -n +5 ../GEO2/WORK/daltaoin >> daltaoin
    !   5. $COLUMBUS/dalton.x -m 1000 > hermitls

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: file
    REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: ovl
    INTEGER(KIND=ilong), INTENT(OUT) :: nl
    INTEGER(KIND=ilong), INTENT(OUT) :: nr
    LOGICAL, INTENT(IN) :: same_aos    

    INTEGER:: nao
    INTEGER:: aoints
    INTEGER:: i
    INTEGER:: j
    INTEGER:: ierr
    INTEGER:: lenbuf
    INTEGER:: nbft
    INTEGER:: nenrgy
    INTEGER:: ninfo
    INTEGER:: nmap
    INTEGER:: nmax
    INTEGER:: nsym
    INTEGER:: ntitle
    CHARACTER*80:: title(20)
    INTEGER:: nbpsy(8)
    INTEGER:: allocstat

    INTEGER, PARAMETER :: nbfmxp=511
    INTEGER, PARAMETER :: nengmx=20
    INTEGER, PARAMETER :: nmapmx=10

    INTEGER:: ietype(nengmx)
    INTEGER:: info(10)
    INTEGER:: imtype(nmapmx)
    INTEGER:: map(nbfmxp*nmapmx)
    REAL*8:: energy(nengmx)
    CHARACTER*4:: slabel(8)
    CHARACTER*8:: bfnlab(nbfmxp)

    ! Use individual allocatable arrays rather than the COLUMBUS style "core" array
    REAL*8,  DIMENSION(:),   ALLOCATABLE :: buf_alloc
    REAL*8,  DIMENSION(:),   ALLOCATABLE :: val_alloc
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ilab_alloc
    INTEGER, DIMENSION(1)                :: ibitv_st

    aoints=freeunit()
    
    OPEN(UNIT=aoints,FILE=trim(file),FORM='unformatted',STATUS='old')

    CALL sifrh1( aoints, ntitle, nsym, nbft, ninfo, nenrgy, nmap, ierr )
    IF(nsym.GT.1)THEN
      WRITE(0,*) "This is a aoints file from a symmetry run. I can handle only C1."
      WRITE(6,*) "This is a aoints file from a symmetry run. I can handle only C1."
      STOP 1
    END IF
    
    IF(debug)WRITE(6,*) 'aoints: nbft = ', nbft

    nao=nbft
    
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
    
    CALL sifrh2( aoints, ntitle, nsym, nbft, ninfo, nenrgy, nmap, &
      &  title, nbpsy, slabel, info, bfnlab, ietype, energy, imtype, map, ierr )

    lenbuf = info(2)
    nmax   = info(3)
    
    allocstat=myalloc(buf_alloc, lenbuf, '+ buf_alloc')
    IF(allocstat.NE.0) STOP
    allocstat=myalloc(ilab_alloc, 2, nmax, '+ ilab_alloc')
    IF(allocstat.NE.0) STOP
    allocstat=myalloc(val_alloc, nmax, '+ val_alloc')
    IF(allocstat.NE.0) STOP
            
    allocstat=myalloc(ovl,nl,nr,'+ ovl')
    IF(allocstat.NE.0) STOP

    CALL process1e(aoints, info, buf_alloc, ilab_alloc, val_alloc, ibitv_st, ovl, nl, same_aos)
    
    IF(debug)THEN
      WRITE(6,*) "Overlap matrix read from Dalton ",trim(adjustl(file))
      WRITE(6,*)nl, nr
      DO i=1,nl
        DO j=1,nr
          WRITE(6,'(x,ES24.14E3)',ADVANCE='no')ovl(i,j)
        END DO
        WRITE(6,*)
      END DO
    END IF ! debug     
    
    allocstat = mydealloc(buf_alloc)
    allocstat = allocstat + mydealloc(ilab_alloc)
    allocstat = allocstat + mydealloc(val_alloc)
    IF(allocstat.NE.0)THEN
      WRITE(6,109) allocstat
      WRITE(0,109) allocstat
    END IF
109 FORMAT("Problem when deallocating. allocstat = ", I8)
       
  END SUBROUTINE read_dalton_overlap
  
  !-------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE process1e(ntape, info, buf, ilab, val, ibitv, ovl, nl, same_aos)
    ! This is the analogue to prt1e of iwfmt.f

    IMPLICIT NONE
    INTEGER:: nl
    LOGICAL:: same_aos
    REAL(KIND=dop), DIMENSION(:,:), INTENT(OUT) :: ovl
    INTEGER:: nipv
    INTEGER:: msame
    INTEGER:: nmsame
    INTEGER:: nomore
    PARAMETER(nipv=2, msame=0, nmsame=1, nomore=2)

    INTEGER:: ntape
    REAL*8:: buf(*)
    REAL*8:: val(*)
    INTEGER:: info(*)
    INTEGER:: ilab(nipv,*)
    INTEGER:: ibitv(1)

    INTEGER:: ibvtyp
    INTEGER:: ierr
    INTEGER:: itypea
    INTEGER:: itypeb
    INTEGER:: last
    INTEGER:: num

    INTEGER:: iretbv
    PARAMETER(iretbv=0)
      
    INTEGER:: i

    REAL*8:: fcore

    last = msame
    ovl = dZero
    
    DO WHILE (last .NE. nomore)

      CALL sifrd1( ntape, info, nipv, iretbv, buf, num, last, itypea, &
        &    itypeb, ibvtyp, val, ilab, fcore, ibitv, ierr )

      IF ( ierr.NE.0 ) THEN
        WRITE(6,107) ierr
        WRITE(0,107) ierr
        STOP 107
      END IF
107   FORMAT("Problem reading SIFS file. ierr = ", I5)

      IF(debug)WRITE(6,*)'itypea, itypeb:',itypea,itypeb
        
      IF((itypea.EQ.0).AND.(itypeb.EQ.0))THEN
        IF (debug) WRITE(6,'("Reading SIFS overlap matrix. num = ", I6)') num
        
        IF (same_aos) THEN
          DO i = 1,num
            ovl(ilab(2,i), ilab(1,i)) = val(i)
            ovl(ilab(1,i), ilab(2,i)) = val(i)
          END DO
        ELSE
          DO i = 1,num
            IF ((ilab(1,i) > nl).AND.(ilab(2,i) <= nl)) THEN
              ovl(ilab(2,i), ilab(1,i)-nl) = val(i)
            END IF
          END DO
        END IF
        
        ! Return unless the next record also contains overlap integrals
        !   "msame" means "more of the same"
        IF (last .NE. msame) RETURN
      END IF
    END DO
      
    ! This part is only reached in the case that the overlap integrals were not on the file.
    WRITE(0,108)
    WRITE(6,108)
    STOP 108
108 FORMAT("Overlap integrals not found.")

  END SUBROUTINE process1e

  !-------------------------------------------------------------------------------------------------------------------------------

END MODULE read_dalton
