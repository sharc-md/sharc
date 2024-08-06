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

MODULE iomod
  !> drivers for reading from files
  !!
  USE sysparam
  !USE global_storage
  !USE read_ascii
  USE read_dalton
  USE read_columbus
  USE read_molcas
  USE read_lumorb
  USE read_turbomole
  USE calcmod
  USE my_alloc
#ifdef _OPENMP
  USE omp_lib
#else
#  define omp_get_num_threads() 1
#  define omp_get_max_threads() 1
#  define omp_get_thread_num() 0
#endif
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_civec, read_mo, read_ao
CONTAINS

  !-------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_CIvec(file,format,coefs,SSD_a,SSD_b,nMO,nSD,nstate,nelec,nalpha,nbeta,ncore)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: file
    INTEGER(KIND=ishort), INTENT(IN) :: format ! 0...native ascii
    REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: coefs
    INTEGER(KIND=ishort), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: SSD_a ! a=alpha, b=beta
    INTEGER(KIND=ishort), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: SSD_b ! a=alpha, b=beta
    INTEGER(KIND=ilong), INTENT(OUT) :: nMO
    INTEGER(KIND=ilong), INTENT(OUT) :: nSD
    INTEGER(KIND=ilong), INTENT(OUT) :: nstate
    INTEGER(KIND=ilong), INTENT(OUT) ::nelec
    INTEGER(KIND=ilong), INTENT(OUT) ::nalpha
    INTEGER(KIND=ilong), INTENT(OUT) ::nbeta
    INTEGER(KIND=ilong), INTENT(IN) :: ncore

    SELECT CASE(format)
      CASE(0)
        CALL read_civec_native(file,coefs,SSD_a,SSD_b,nMO,nSD,nstate,nelec,nalpha,nbeta,ncore)
    END SELECT
  END SUBROUTINE read_CIvec

  !--------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_AO(file,format,ovl,nl,nr,same_aos)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: file
    INTEGER(KIND=ishort), INTENT(IN) :: format
    REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: ovl
    INTEGER(KIND=ilong), INTENT(OUT) :: nl
    INTEGER(KIND=ilong), INTENT(OUT) :: nr
    LOGICAL, INTENT(IN) :: same_aos
    
    SELECT CASE(format)
      CASE(0) ! native ascii
        CALL read_aoovl_native(file,ovl,nl,nr)
      CASE(1) ! molcas ONEINT
        CALL read_molcas_overlap(file,ovl,nl,nr,same_aos)
      CASE(2) ! dalton
        CALL read_dalton_overlap(file,ovl,nl,nr,same_aos)
    END SELECT
    
!   TODO: In debug mode, the determinant function sometimes gives segfaults
!     IF (format.ge.0) THEN
!       IF ((debug).AND.(nl.EQ.nr)) THEN    
!        WRITE(6,'("Determinant of AO-overlap matrix:",ES24.14E3)') determinant(ovl, nl)
!       END IF
!     END IF

  END SUBROUTINE read_AO

  !--------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_MO(file,format,coefs,nMO,nAO)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: file
    INTEGER(KIND=ishort), INTENT(IN) :: format
    REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: coefs
    INTEGER(KIND=ilong), INTENT(IN) :: nMO
    INTEGER(KIND=ilong), INTENT(INOUT) :: nAO

    SELECT CASE(format)
      CASE(0) ! columbus (=native) format
        CALL get_Columbus_MO(file,coefs,nMO,nAO)
      CASE(1) ! molcas lumorb format
        CALL read_molcas_lumorb(file,coefs,nMO,nAO)
      CASE(2) ! turbomole mos format
        CALL read_turbomole_mos(file,coefs,nMO,nAO)
      CASE DEFAULT
        WRITE(6,*) "MO format not implemented"
        WRITE(0,*) "MO format not implemented"
        STOP 1
    END SELECT

  END SUBROUTINE read_MO

  !--------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_civec_native(file,coefs,SSD_a,SSD_b,nMO,nSD,nstate,nelec,nalpha,nbeta,ncore)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: file
    REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: coefs
    INTEGER(KIND=ishort), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: SSD_a
    INTEGER(KIND=ishort), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: SSD_b
    INTEGER(KIND=ilong), INTENT(OUT) :: nMO
    INTEGER(KIND=ilong), INTENT(OUT) :: nSD
    INTEGER(KIND=ilong), INTENT(OUT) :: nstate
    INTEGER(KIND=ilong), INTENT(OUT) ::nelec
    INTEGER(KIND=ilong), INTENT(OUT) ::nalpha
    INTEGER(KIND=ilong), INTENT(OUT) ::nbeta
    INTEGER(KIND=ilong), INTENT(IN) :: ncore

    INTEGER:: detflio
    LOGICAL:: test
    INTEGER:: iost
    INTEGER:: allocstat

    CHARACTER*4000:: detstring, coefstring
    INTEGER(KIND=ilong):: i
    INTEGER(KIND=ilong):: j
    INTEGER(KIND=ilong):: k
    INTEGER(KIND=ilong):: nel
    INTEGER(KIND=ilong):: na
    INTEGER(KIND=ilong):: nb
    INTEGER(KIND=ilong):: ninv

    CHARACTER*32:: formatspec


    INQUIRE(FILE=trim(adjustl(file)), EXIST=test)
    IF(.NOT.test)THEN
      WRITE(0,*) "File ",trim(adjustl(file))," not found"
      WRITE(6,*) "File ",trim(adjustl(file))," not found"
      STOP 1
    ELSE
      detflio=freeunit()
      OPEN(UNIT=detflio,FILE=trim(adjustl(file)),IOSTAT=iost,STATUS='old',ACTION='read')
      IF(iost .NE. 0) THEN
        WRITE(0,*) "Cannot open file ",trim(adjustl(file))," to read"
        WRITE(6,*) "Cannot open file ",trim(adjustl(file))," to read"
        STOP 1
      END IF
    END IF

    READ(detflio,*)nstate, nMO, nSD
    IF(debug)THEN
      WRITE(6,*) nstate, nMO, nSD
    ENDIF
    allocstat=myalloc(coefs,nstate,nSD,'+ cicoef')
    IF(allocstat.NE.0)THEN
      WRITE(6,'("Could not allocate cicoefs in read_dets; error ",I5)')allocstat
      WRITE(6,*) "file: ",trim(adjustl(file))
      WRITE(0,'("Could not allocate cicoefs in read_dets; error ",I5)')allocstat
      WRITE(0,*) "file: ",trim(adjustl(file))
      STOP 1
    END IF
    ! semi-elegant: read the first determinant string to get the number of electrons
    WRITE(formatspec,'( "(A",I0,")" )') nMO
    READ(detflio,formatspec,IOSTAT=iost)detstring
    IF(debug)THEN
      WRITE(6,*) '"'//trim(detstring)//'"'
    ENDIF
    WRITE(formatspec,'( "(A",I0,",A)" )') nMO

    na=0
    nb=0
    nel=0

    DO j=ncore+1, nMO
      SELECT CASE (detstring(j:j))
        CASE ('a')
          na = na+1
          nel = nel+1
        CASE ('b')
          nb = nb+1
          nel = nel+1
        CASE ('e')
          ! do nothing
          CYCLE
        CASE ('d')
          na = na+1
          nb = nb+1
          nel = nel+2
        CASE DEFAULT
          WRITE(0,'("error in reading determinant: ",I6," MO ",I6)')i,j
          WRITE(6,*) "file: ",trim(adjustl(file))
          WRITE(6,'("error in reading determinant: ",I6," MO ",I6)')i,j
          WRITE(0,*) "file: ",trim(adjustl(file))
          STOP 1
      END SELECT
    END DO

    nalpha=na
    nbeta=nb
    nelec=nel

    allocstat=myalloc(SSD_a,nSD,na,'+ alpha SDs')
    IF(allocstat.NE.0)THEN
      WRITE(6,'("Could not allocate alphaSDs in read_civec_native; error ",I5)')allocstat
      WRITE(6,*) "file: ",trim(adjustl(file))
      WRITE(0,'("Could not allocate alphaSDs in read_civec_native; error ",I5)')allocstat
      WRITE(0,*) "file: ",trim(adjustl(file))
      STOP 1
    END IF

    allocstat=myalloc(SSD_b,nSD,nb,'+ betadets')
    IF(allocstat.NE.0)THEN
      WRITE(6,'("Could not allocate betaSDs in read_civec_native; error ",I5)')allocstat
      WRITE(6,*) "file: ",trim(adjustl(file))
      WRITE(0,'("Could not allocate betaSDs in read_civec_native; error ",I5)')allocstat
      WRITE(0,*) "file: ",trim(adjustl(file))
      STOP 1
    END IF

    REWIND(detflio)
    READ(detflio,*) ! discard the fist line

    DO i=1,nSD
!       READ(detflio,*,IOSTAT=iost)detstring,coefs(1:nstate,i)
      READ(detflio,formatspec,IOSTAT=iost)detstring, coefstring
      READ(coefstring,*) coefs(1:nstate,i)
      IF(iost.NE.0)THEN
        WRITE(0,'("error in reading SD: ",I6,", file: ",A200)')i, trim(adjustl(file))
        PRINT*,coefs(nstate,i),detstring
        WRITE(6,'("error in reading SD: ",I6,", file: ",A200)')i, trim(adjustl(file))
        STOP 1
      END IF
      na=0
      nb=0
      nel=0
      ninv = 0

      ! Core orbitals are regarded as empty.
      ! Numerically, this corresponds to assuming a unit overlap between them.
      ! This improves the computational performance as well as the numerical stability
      !     for large displacements.

      DO j=ncore+1, nMO
        SELECT CASE (detstring(j:j))
          CASE ('a')
            na = na+1
            nel = nel+1
            SSD_a(i,na)=j
          CASE ('b')
            nb = nb+1
            nel = nel+1
            SSD_b(i,nb)=j
            ninv = ninv + nalpha - na
          CASE ('e')
            ! do nothing
            CYCLE
          CASE ('d')
            na = na+1
            nb = nb+1
            nel = nel+2
            SSD_a(i,na)=j
            SSD_b(i,nb)=j
            ninv = ninv + nalpha - na
          CASE DEFAULT
            WRITE(0,'("error in reading SD: ",I6,", file: ",A200)')i, trim(adjustl(file))
            PRINT*,coefs(nstate,i),detstring
            WRITE(6,'("error in reading SD: ",I6,", file: ",A200)')i, trim(adjustl(file))
            STOP 1
        END SELECT
      END DO ! j=ncore+1,nMO

      IF(nel.NE.nelec)THEN
        WRITE(0,'("Wrong electron count in determinant ",I6,x,I3,x,I3,", file: ",A200)')i,nel,nelec,file
        WRITE(6,'("Wrong electron count in determinant ",I6,", file: ",A200)')i, file
        STOP 1
      END IF

      coefs(:, i) = coefs(:, i) * (-1)**ninv

    END DO ! i=1,nSD

    CLOSE(detflio)

    IF(debug)THEN
      WRITE(6,*)"CI coefficients and Slater determinants read from ",trim(adjustl(file))
      DO i=1,min(nSD,maxdebug)
        WRITE(6, '(I8,x)',ADVANCE='no')i
        WRITE(6,'(10E18.8)',ADVANCE='yes')(coefs(k, i),k=1,nstate)
        WRITE(6,'(" alpha ")',ADVANCE='no')
        DO j=1,na-1
          WRITE(6,'(x,I3)',ADVANCE='no')SSD_a(i,j)
        END DO
        j=na
        WRITE(6,'(x,I3)',ADVANCE='yes')SSD_a(i,j)
        WRITE(6,'(" beta  ")',ADVANCE='no')
        DO j=1,nb-1
          WRITE(6,'(x,I3)',ADVANCE='no')SSD_b(i,j)
        END DO
        j=nb
        WRITE(6,'(x,I3)',ADVANCE='yes')SSD_b(i,j)
      END DO
    END IF

  END SUBROUTINE read_civec_native

  !--------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_aoovl_native(file,ovl,nl,nr)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: file
    REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: ovl
    INTEGER(KIND=ilong), INTENT(OUT) :: nl
    INTEGER(KIND=ilong), INTENT(OUT) :: nr


    INTEGER(KIND=ilong):: i
    INTEGER(KIND=ilong):: j
    INTEGER:: allocstat
    INTEGER:: io_status
    INTEGER:: iou_ovl

    iou_ovl=freeunit()
    OPEN(iou_ovl,FILE=file,IOSTAT=io_status)
    IF(io_status.NE.0)THEN
      WRITE(6,*)"Could not open AO overlap integral file"
      WRITE(6,*) "file: ",trim(adjustl(file))
      WRITE(0,*)"Could not open AO overlap integral file"
      WRITE(0,*) "file: ",trim(adjustl(file))
      STOP 1
    END IF
    READ(iou_ovl,*) nl, nr

    allocstat=myalloc(ovl,nl,nr,'+ ao_ovl')
    IF(allocstat.NE.0) STOP
    DO i=1,nl
      READ(iou_ovl,*)(ovl(i,j),j=1,nr)
    END DO
    IF(debug)THEN
      WRITE(6,*) "Overlap matrix read from ascii ",trim(adjustl(file))
      DO i=1,nl
        DO j=1,nr-1
          WRITE(6,'(x,ES24.14E3)',ADVANCE='no')ovl(i,j)
        END DO
        WRITE(6,'(x,ES24.14E3)',ADVANCE='yes')ovl(i,nr)
      END DO
    END IF ! debug
    CLOSE(iou_ovl)

  END SUBROUTINE read_aoovl_native

  !--------------------------------------------------------------------------------------------------------------------------------

END MODULE iomod
