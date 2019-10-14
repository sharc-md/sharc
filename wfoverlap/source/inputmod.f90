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

MODULE inputmod
  !> module for reading the input
  !! subroutine read_input opens the input file (if existent) and reads it (default filename 'cioverlap.input')
  USE sysparam
  USE global_storage
  USE my_alloc

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_input
!--------------------------------------------------------------------------------------------------------------------------------==
!--------------------------------------------------------------------------------------------------------------------------------==
CONTAINS

  !--------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_input
    IMPLICIT NONE

    CHARACTER*200, DIMENSION(:), ALLOCATABLE :: clargs ! command line arguments
    CHARACTER*200:: lcclarg   ! lower case version of a command line argument
    CHARACTER*200:: line
    CHARACTER*200:: temp
    CHARACTER*200:: varname
    CHARACTER*200:: value
    CHARACTER*200:: lcval
    CHARACTER*200:: inp_fn
    INTEGER:: nclargs
    INTEGER:: i
    INTEGER:: sindex
    INTEGER:: iost
    INTEGER:: io_inp
    LOGICAL:: readthis
    LOGICAL:: test

    WRITE(6,*),"Reading input"

    inp_fn='cioverlap.input' ! default value

    nclargs=iargc()
    ALLOCATE(clargs(1:nclargs))

    readthis=.TRUE.
    ! read command line arguments
    DO i=1,nclargs
      CALL getarg(i,clargs(i))
      clargs(i)=trim(adjustl(clargs(i)))
    END DO
    DO i=1,nclargs
      CALL lc(clargs(i),lcclarg)
      IF(.NOT.readthis) THEN
        readthis=.TRUE.
        CYCLE
      ELSE
        IF (trim(lcclarg) == '-f') THEN
          IF ((i+1).GT.nclargs) THEN
            STOP '-f must be followed by a filename'
          ELSE
            inp_fn = trim(adjustl(clargs(i+1)))
            readthis=.FALSE.
            PRINT*,"Reading input from file: ", trim(inp_fn)
          END IF
        ELSEIF (trim(lcclarg) == "--debug") THEN
          debug = .TRUE.
        ELSEIF (trim(lcclarg) == "-m") THEN
          IF ((i+1).GT.nclargs) THEN
            STOP '-m must be followed by a number'
          ELSE
            value=trim(adjustl(clargs(i+1)))
            READ(value,*)maxmem
            maxmem=maxmem*1024**2
            readthis=.FALSE.
          END IF
        ELSE
          WRITE(0,*) 'Cannot understand: ',trim(lcclarg)
          STOP 'Input error'
        END IF
      END IF
    END DO

    DEALLOCATE(clargs)

    INQUIRE(FILE=trim(adjustl(inp_fn)),EXIST=test)
    IF(.NOT. test) THEN
      WRITE(0,*) 'Input file -',trim(adjustl(inp_fn)),'- not found'
      WRITE(0,*) 'To use a different file, use the -f option'
      STOP 5
    END IF
    io_inp=freeunit()
    OPEN(UNIT=io_inp, FILE=trim(adjustl(inp_fn)), STATUS='old', IOSTAT=iost)
    IF(iost .NE. 0 )THEN
      STOP 'error reading input file'
    END IF

    mix_AOfile=''
    bra_MOfile=''
    ket_MOfile=''
    bra_CIfile=''
    ket_CIfile=''

    bra_nstate = iOne
    ket_nstate = iOne
    bra_nSD = iZero
    ket_nSD = iZero

    ncore = iZero
    ndocc = iZero
    AOformat = iZero
    force_direct_dets = .FALSE.
    force_noprecalc = .FALSE.
    do_mixing_angles = .FALSE.
    same_aos = .FALSE.

    bra_nAO = -1
    ket_nAO = -1

    bra_CIformat=0
    ket_CIformat=0
    bra_MOformat=0
    ket_MOformat=0
    MOprint = 0

    DO
      ! read and process the line
      ! split var=value lines in two strings
      line=''
      READ(io_inp,'(A195)',IOSTAT=iost) line ! Somehow this line is seen as a memory leak by Intel Inspector ...
      IF(iost .NE. 0) THEN
        EXIT
      END IF
      sindex=0
      varname=''
      value=''
      sindex=INDEX(line,'#')
      IF(sindex.NE.0)THEN
        line=line(1:sindex-1)
        sindex=0
      END IF

      ! remove empty lines
      IF(len(trim(line)).EQ.0)THEN
        CYCLE
      END IF

      sindex=INDEX(line,'=')
      IF(sindex.NE.0)THEN
        temp=trim(adjustl(line(1:sindex-1)))
        CALL lc(temp,varname)
        value=trim(adjustl(line(sindex+1:)))
        CALL lc(value,lcval)
        sindex=0
      ELSE
        temp=trim(adjustl(line))
        CALL lc(temp,varname)
        value='blank'
      END IF
      ! print*,"key: ",trim(varname),", value: ",trim(value)," | ",trim(lcval)

      IF(trim(varname) == 'mix_aoovl')THEN
        mix_AOfile = trim(value)
      ELSEIF(trim(varname) == 'a_mo')THEN
        bra_MOfile = trim(value)
      ELSEIF(trim(varname) == 'b_mo')THEN
        ket_MOfile = trim(value)
      ELSEIF(trim(varname) == 'a_det')THEN
        bra_CIfile = trim(value)
      ELSEIF(trim(varname) == 'b_det')THEN
        ket_CIfile = trim(value)
      ELSEIF(trim(varname) == 'ncore')THEN
        READ(value,*)ncore
      ELSEIF(trim(varname) == 'ndocc')THEN
        READ(value,*)ndocc
      ELSEIF(trim(varname) == 'ao_read')THEN
        READ(value,*)AOformat
      ELSEIF(trim(varname) == 'a_mo_read')THEN
        READ(value,*)bra_MOformat
      ELSEIF(trim(varname) == 'b_mo_read')THEN
        READ(value,*)ket_MOformat
      ELSEIF(trim(varname) == 'force_direct_dets')THEN
        force_direct_dets=.TRUE.
      ELSEIF(trim(varname) == 'force_noprecalc')THEN
        force_noprecalc=.TRUE.
      ELSEIF(trim(varname) == 'same_aos')THEN
        same_aos=.TRUE.
      ELSEIF(trim(varname) == 'nao_a')THEN
        READ(value,*)bra_nAO
      ELSEIF(trim(varname) == 'nao_b')THEN
        READ(value,*)ket_nAO
      ELSEIF(trim(varname) == 'moprint')THEN
        READ(value,*)MOprint
      ELSEIF(trim(varname) == 'swap')THEN
        swap=.TRUE.
      ELSEIF(trim(varname) == 'mixing_angles')THEN
        do_mixing_angles=.TRUE.
      ELSE
        WRITE(0,*)"Unknown option in input file: ",trim(varname)," = ",trim(value)
        WRITE(6,*)"Unknown option in input file: ",trim(varname)," = ",trim(value)
        STOP 1
      END IF
    END DO

    if(ndocc.lt.ncore)then
      ndocc = ncore
    end if

    IF(trim(bra_MOfile)=='') STOP "must give A mocoef filename in input"
    IF(trim(ket_MOfile)=='') STOP "must give B mocoef filename in input"
    IF(trim(bra_CIfile)=='') STOP "must give A determinants filename in input"
    IF(trim(ket_CIfile)=='') STOP "must give B determinants filename in input"

    IF(AOformat.EQ.0)THEN
      IF(mix_AOfile.EQ.'') mix_AOfile='S_mix'
    ELSEIF(AOformat.EQ.1)THEN
      IF(mix_AOfile.EQ.'') mix_AOfile='ONEINT'
    ELSEIF(AOformat.EQ.2)THEN
      IF(mix_AOfile.EQ.'') mix_AOfile='aoints'
    ELSEIF(AOformat.EQ.-1)THEN ! reconstruct according to S.C_B^T = C_B^-1
      IF(.not.same_aos)THEN
        WRITE(6,*) "MO-coefficient inversion only possible for same_aos"
        WRITE(0,*) "MO-coefficient inversion only possible for same_aos"
        STOP 1           
      ENDIF    
    ELSE
      WRITE(6,*)'Method for ao_read not implemented:',AOformat
      STOP 4
    END IF


  END SUBROUTINE read_input

  !--------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE get_next_substr(string,substring,separator,success)
    !> finds the next occurrance of the separator in the string
    !! returns as substring everything left of the separator
    !! returns a string everything right of the separator
    !! returns the string as it was if the separator was not found or if
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: separator
    CHARACTER(LEN=*), INTENT(INOUT) :: string
    CHARACTER(LEN=*), INTENT(OUT) :: substring
    LOGICAL, INTENT(OUT) :: success
    INTEGER:: sindex

    success=.FALSE.

    sindex=INDEX(string,separator)
    IF(sindex.GT.1)THEN
      success=.TRUE.
      substring=trim(adjustl(string(1:sindex-1)))
      string=trim(adjustl(string(sindex+1:)))
    END IF
  END SUBROUTINE get_next_substr

  !--------------------------------------------------------------------------------------------------------------------------------

END MODULE inputmod
