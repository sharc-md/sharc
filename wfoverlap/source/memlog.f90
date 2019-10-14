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

MODULE memlog
  USE sysparam
#ifdef _OPENMP
  USE omp_lib
#endif
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: memu, mem_used, log_memory

  INTEGER(KIND=ilong):: memu
  INTEGER(KIND=ilong):: mem_used=0

CONTAINS

  SUBROUTINE log_memory(msg,thrd)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER(KIND=ilong), INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30, SAVE:: nm=" "
    LOGICAL, SAVE :: init=.FALSE.
#ifdef _OPENMP
    REAL(KIND=dop), SAVE :: starttime
    REAL(KIND=dop):: time
#else
    INTEGER (KIND=ilong), SAVE :: starttime
    INTEGER (KIND=ilong):: time
#endif
    INTEGER (KIND=ilong), SAVE :: timeu
    INTEGER, SAVE :: tmfio

    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(present(thrd))THEN
      WRITE(nm,'(A22,x,"(",I5,")")')trim(adjustl(nm)),thrd
    END IF

    IF(.NOT. init) THEN
      tmfio=freeunit()
      OPEN(UNIT=tmfio,FILE='memlog')
      init=.TRUE.
#ifdef _OPENMP
      starttime=omp_get_wtime()
      timeu=1
#else
      CALL system_clock(starttime, timeu)
#endif
    END IF
#ifdef _OPENMP
      time=omp_get_wtime()
#else
    CALL system_clock(time)
#endif
    WRITE(tmfio,'(F9.4,X,I14,X,F12.1,X,F10.1,X,F7.1,X,A40)') &
      & real(time-starttime)/real(timeu), &
      & mem_used, real(mem_used)/1024.0, real(mem_used)/(1024.0**2), real(mem_used)/(1024.0**3), &
      & trim(adjustl(nm))
  END SUBROUTINE log_memory

END MODULE memlog
