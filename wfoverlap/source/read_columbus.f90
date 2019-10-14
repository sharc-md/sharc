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

MODULE read_columbus
  !> module for things to read in Columbus format
  !!   No libraries are required for running this.
  !!   mocoef files -> get_both_col_MOs
  !!
  !!  private subroutines read_col_moheader and read_col_mos

  USE sysparam
  !USE global_storage
  USE my_alloc
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: get_Columbus_MO

  !--------------------------------------------------------------------------------------------------------------------------------

CONTAINS

  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE get_Columbus_MO(file,coef,nmo,nao)

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: file
    REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: coef
    INTEGER(KIND=ilong), INTENT(IN) :: nmo
    INTEGER(KIND=ilong), INTENT(INOUT) :: nao
    INTEGER(KIND=ilong):: nmo_read

    INTEGER:: allocstat

    CALL read_col_moheader(file, nao, nmo_read)
    IF(nmo_read.NE.nmo)THEN
      WRITE(6,'("Number of MOs in mocoef file does not agree with CIvec: ",I6," instead of ",I6)')nmo_read,nmo
      WRITE(6,*) "file: ",trim(adjustl(file))
      WRITE(0,'("Number of MOs in mocoef file does not agree with CIvec: ",I6," instead of ",I6)')nmo_read,nmo
      WRITE(0,*) "file: ",trim(adjustl(file))
      STOP 1
    END IF

    allocstat=myalloc(coef,nmo,nao,"+ mocoef")
    IF(allocstat.NE.0)THEN
      WRITE(6,'("Could not allocate mocoef in get_Columbus_MO; error:",I4)')allocstat
      WRITE(6,*) "file: ",trim(adjustl(file))
      WRITE(0,'("Could not allocate mocoef in get_Columbus_MO; error:",I4)')allocstat
      WRITE(0,*) "file: ",trim(adjustl(file))
      STOP 1
    END IF

    CALL read_col_mos(file,nao,nmo,coef)

  END SUBROUTINE get_Columbus_MO

  !--------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_col_moheader(mofile, nao, nmo)
    !> subroutine to read the header of an MOCOEF-file
    !! compares the number of AOs given to the information in the header, STOPS if not equal
    !! returns the number of MOs as found in the header
    !! STOPS for symmetry!=C1
    !!
    !! Expects the AO-Overlap file to be read before and number of AOs to be known.
    !! The mofile is closed again after reading the header.
    !!
    !! variables:
    !!   input:
    !!     mofile ... character, filename to read from
    !!     nao ... integer(kind=ilong), expected number of AOs
    !!   output:
    !!     nmo ... integer(kind=ilong), number of MOs in the file
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: mofile
    INTEGER(KIND=ilong), INTENT(INOUT) :: nao
    INTEGER(KIND=ilong), INTENT(OUT) :: nmo

    LOGICAL:: test
    CHARACTER*200:: line
    CHARACTER*200:: lcline
    INTEGER:: moflio
    INTEGER:: iost
    INTEGER(KIND=ilong):: nheaderlines
    INTEGER(KIND=ilong):: sindex
    INTEGER(KIND=ilong):: i
    INTEGER(KIND=ilong):: nsym
    INTEGER(KIND=ilong):: naotmp

    ! look for file and open it
    INQUIRE(FILE=trim(adjustl(mofile)), EXIST=test)
    IF(.NOT.test)THEN
      WRITE(0,*) "File ",trim(adjustl(mofile))," not found"
      WRITE(6,*) "File ",trim(adjustl(mofile))," not found"
      STOP 1
    ELSE
      moflio=freeunit()
      OPEN(UNIT=moflio,FILE=trim(adjustl(mofile)),IOSTAT=iost,STATUS='old',ACTION='read')
      IF(iost .NE. 0) THEN
        WRITE(0,*) "read_col_moheader: Cannot open file ",trim(adjustl(mofile))," to read"
        WRITE(6,*) "read_col_moheader: Cannot open file ",trim(adjustl(mofile))," to read"
        STOP 1
      END IF
    END IF

    DO
      !read till header
      READ(moflio,'(A200)',IOSTAT=iost) line
      IF(iost.NE.0)THEN
        ! break the do-loop when file is finished
        WRITE(0,*) "read_col_moheader: MOcoef file ended before encountering header"
        WRITE(6,*) "read_col_moheader: MOcoef file ended before encountering header"
        STOP 1
      END IF
      sindex = 0
      CALL lc(line,lcline)
      sindex = index(lcline,'header')
      IF(sindex.NE.0)THEN
        EXIT
      END IF
    END DO
    READ(moflio,*) nheaderlines
    ! skip the uninteresting information
    DO i=1,nheaderlines
      READ(moflio,*)
    END DO
    ! check number of symmetries, break if > 1
    READ(moflio,*) nsym
    IF(nsym.NE.1)THEN
      WRITE(0,104) nsym
      WRITE(6,104) nsym
      STOP 1
    END IF
104 FORMAT ("read_col_moheader: This is a MOcoef file from a symmetry run. I can handle only C1. nsym = ", I4)
    ! read number of AOs and MOs. Check number of AOs
    READ(moflio,*)naotmp,nmo
    IF(nao.eq.-1)THEN
      nao = naotmp
    ENDIF
    IF(naotmp.NE.nao)THEN
      WRITE(0,105) naotmp, nao
      WRITE(6,105) naotmp, nao
      STOP 1
    END IF
105 FORMAT ("read_col_moheader: MOcoef file states different number of atomic orbitals &
    &than AOoverlap file. naotmp = "    , I5, ", nao = ", I5)

    CLOSE(moflio)

  END SUBROUTINE read_col_moheader

  !--------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_col_mos(mofile,nao,nmo,mo_coefs)
    !> subroutine to read the MO-coefficients from a MOCOEF-file
    !!
    !! Expects the header of the file to be read before (read_col_moheader) and the array for the coefficients to be allocated.
    !! The mofile is closed again after reading.
    !!
    !! variables:
    !!   input:
    !!     mofile ... character, filename to read from
    !!     nao ... integer(kind=ilong) - expected number of AOs
    !!     nmo ... integer(kind=ilong) - number of MOs in the file
    !!   output:
    !!     mo_coefs ... real(kind=ilong)(nmo,nao) - array of MO coefficients
    IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(IN) :: mofile
    INTEGER(KIND=ilong), INTENT(IN) :: nao
    INTEGER(KIND=ilong), INTENT(IN) :: nmo
    REAL(KIND=dop), DIMENSION(nmo,nao), INTENT(OUT) :: mo_coefs

    INTEGER:: moflio
    INTEGER:: iost
    INTEGER(KIND=ilong):: i
    INTEGER(KIND=ilong):: j
    INTEGER(KIND=ilong):: sindex
    CHARACTER*200:: line
    CHARACTER*200:: lcline
    CHARACTER*200:: mofmt
    LOGICAL:: test
    LOGICAL:: ldio

    ! look for file and open
    INQUIRE(FILE=trim(adjustl(mofile)), EXIST=test)
    IF(.NOT.test)THEN
      WRITE(0,*) "read_col_mos: File ",trim(adjustl(mofile))," not found"
      WRITE(6,*) "read_col_mos: File ",trim(adjustl(mofile))," not found"
      STOP 1
    ELSE
      moflio=freeunit()
      OPEN(UNIT=moflio,FILE=trim(adjustl(mofile)),IOSTAT=iost,STATUS='old',ACTION='read')
      IF(iost .NE. 0) THEN
        WRITE(0,*) "read_col_mos: Cannot open file ",trim(adjustl(mofile))," to read"
        WRITE(6,*) "read_col_mos: Cannot open file ",trim(adjustl(mofile))," to read"
        STOP 1
      END IF
    END IF

    DO
      !read till header
      READ(moflio,'(A200)',IOSTAT=iost) line
      IF(iost.NE.0)THEN
        ! break the do-loop when file is finished
        WRITE(0,*) "read_col_moheader: MOcoef file ended before encountering header"
        WRITE(6,*) "read_col_moheader: MOcoef file ended before encountering header"
        STOP 1
      END IF
      sindex = 0
      CALL lc(line,lcline)
      sindex = index(lcline,'header')
      IF(sindex.NE.0)THEN
        EXIT
      END IF
    END DO
    DO
      !read till keyword 'mocoef' -> thus skip header
      READ(moflio,'(A200)',IOSTAT=iost) line
      IF(iost.NE.0)THEN
        ! break the do-loop when file is finished
        WRITE(0,*) "read_col_mos: MOcoef file ",trim(adjustl(mofile)),"  ended before encountering mocoef-vector"
        WRITE(6,*) "read_col_mos: MOcoef file ",trim(adjustl(mofile)),"  ended before encountering mocoef-vector"
        STOP 1
      END IF
      sindex = 0
      CALL lc(line,lcline)
      sindex = index(lcline,'mocoef')
      IF(sindex.NE.0)THEN
        EXIT
      END IF
    END DO

    !now we are at the beginning of the mocoef-vector
    !read in format
    READ(moflio,*)mofmt
    IF (mofmt.EQ.'(*)') THEN
      ldio=.TRUE.
    ELSE
      ldio=.FALSE.
    ENDIF
    !read in coefficients
    DO i=1,nmo
      IF (ldio) THEN
        READ (moflio,*,IOSTAT=iost) (mo_coefs(i,j),j=1,nao)
      ELSE
        READ (moflio,FMT=trim(adjustl(mofmt)),IOSTAT=iost) (mo_coefs(i,j),j=1,nao)
      ENDIF
      IF (iost.NE.0) THEN
        WRITE(0,'("read_col_mos: Failure in reading MOcoef file ",A," , MO: ",I3)') trim(adjustl(mofile)),i
        WRITE(6,'("read_col_mos: Failure in reading MOcoef file ",A," , MO: ",I3)') trim(adjustl(mofile)),i
        STOP 1
      ENDIF
    END DO

    CLOSE(moflio)
  END SUBROUTINE read_col_mos
  !--------------------------------------------------------------------------------------------------------------------------------
END MODULE read_columbus
