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

MODULE sysparam

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dop, sip, ilong, ishort, ivshort, debug, maxdebug, freeunit, lc, &
    & iOne, iZero, dOne, dZero, d_epsilon, pi, au2eV, au2Ang, maxmem

  INTEGER, PARAMETER :: ilong=SELECTED_INT_KIND(18) ! integer*8
  INTEGER, PARAMETER :: ishort=SELECTED_INT_KIND(9) ! INTEGER*4
  INTEGER, PARAMETER :: ivshort=SELECTED_INT_KIND(1) ! integer*1
  INTEGER, PARAMETER :: sip=SELECTED_REAL_KIND(6,31) ! "normal" real*4
  INTEGER, PARAMETER :: dop=SELECTED_REAL_KIND(15,307) ! double precision, real*8
  
  INTEGER (KIND=ilong), PARAMETER :: iOne=1_ilong
  INTEGER (KIND=ilong), PARAMETER :: iZero=0_ilong
  REAL (KIND=dop), PARAMETER :: dOne=1.0_dop
  REAL (KIND=dop), PARAMETER :: dZero=0.0_dop
  REAL (KIND=dop), PARAMETER :: pi= 2.0_dop*acos(dZero)
  REAL (KIND=dop), PARAMETER :: au2eV=27.2113961_dop
  REAL (KIND=dop), PARAMETER :: au2Ang=1.889726143_dop

  INTEGER (KIND=ilong):: maxmem=1073741824

  LOGICAL:: debug=.FALSE.
  INTEGER (KIND=ilong), PARAMETER :: maxdebug=2000

  REAL (KIND=dop), PARAMETER :: d_epsilon=epsilon(1.0_dop)


!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
CONTAINS

  !--------------------------------------------------------------------------------------------------------------------------------

  INTEGER FUNCTION freeunit()
    LOGICAL:: used
    INTEGER, PARAMETER :: iinp=5
    INTEGER, PARAMETER :: ioutp=6
    freeunit = 0
    used = .TRUE.
    DO WHILE (used)
      freeunit = freeunit + 1
      IF (freeunit.NE.iinp .AND. freeunit.NE.ioutp) THEN
        IF (freeunit .GT. 99) THEN
          STOP 'no free I/O unit'
        END IF
        INQUIRE (UNIT=freeunit,OPENED=used)
      END IF
    END DO
  END FUNCTION freeunit

  !--------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE lc(string, lcstring)
    CHARACTER (LEN=*), INTENT(IN):: string
    CHARACTER (LEN=*), INTENT(OUT) :: lcstring
    INTEGER:: capdiff
    INTEGER:: lclim
    INTEGER:: uclim
    INTEGER:: lslim
    INTEGER:: sav
    INTEGER:: i

    lcstring = ''

    lclim=IACHAR('A')
    uclim=IACHAR('Z')
    lslim=IACHAR('a')
    capdiff=lclim-lslim
    DO i=1,len(string)
      sav=IACHAR(string(i:i))
      IF((sav.GE.lclim).AND.(sav.LE.uclim))THEN
        lcstring(i:i)=ACHAR(sav-capdiff)
      ELSE
        lcstring(i:i)=achar(sav)
      END IF
    END DO
  END SUBROUTINE lc


END MODULE sysparam
