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

MODULE my_alloc
  USE sysparam
  USE memlog
  IMPLICIT NONE
  INTERFACE myalloc
    MODULE PROCEDURE &
    &  ail1d,ail2d,ail3d,ail4d,ail5d,ail6d,ail7d &
    & ,ais1d,ais2d,ais3d,ais4d,ais5d,ais6d,ais7d &
    & ,ar1d,ar2d,ar3d,ar4d,ar5d,ar6d,ar7d &
    & ,ad1d,ad2d,ad3d,ad4d,ad5d,ad6d,ad7d &
    & ,ac1d,ac2d,ac3d,ac4d,ac5d,ac6d,ac7d &
    & ,az1d,az2d,az3d,az4d,az5d,az6d,az7d &
    & ,al1d,al2d,al3d,al4d,al5d,al6d,al7d
  END INTERFACE
  INTERFACE mydealloc
    MODULE PROCEDURE &
    &  dil1d,dil2d,dil3d,dil4d,dil5d,dil6d,dil7d &
    & ,dis1d,dis2d,dis3d,dis4d,dis5d,dis6d,dis7d &
    & ,dr1d,dr2d,dr3d,dr4d,dr5d,dr6d,dr7d &
    & ,dd1d,dd2d,dd3d,dd4d,dd5d,dd6d,dd7d &
    & ,dc1d,dc2d,dc3d,dc4d,dc5d,dc6d,dc7d &
    & ,dz1d,dz2d,dz3d,dz4d,dz5d,dz6d,dz7d &
    & ,dl1d,dl2d,dl3d,dl4d,dl5d,dl6d,dl7d
  END INTERFACE
CONTAINS

!===============================================

  INTEGER FUNCTION ail1d(arr,d1,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ilong),DIMENSION(:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1_ilong)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ail1d=-1
      WRITE(6,*)'array already allocated ail1d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ail1d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ail1d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1),STAT=ail1d)
        IF(ail1d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ail1d'
          WRITE(0,*)'error in allocation routine ail1d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ail1d

  !---------------------------------

  INTEGER FUNCTION dil1d(arr,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ilong),DIMENSION(:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dil1d=-1
      WRITE(6,*)'array not allocated dil1d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dil1d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dil1d)
      IF(dil1d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dil1d'
        WRITE(0,*)'error in deallocation routine dil1d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dil1d

  !---------------------------------

  INTEGER FUNCTION ail2d(arr,d1,d2,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ilong),DIMENSION(:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1_ilong)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ail2d=-1
      WRITE(6,*)'array already allocated ail2d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ail2d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ail2d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2),STAT=ail2d)
        IF(ail2d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ail2d'
          WRITE(0,*)'error in allocation routine ail2d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ail2d

  !---------------------------------

  INTEGER FUNCTION dil2d(arr,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ilong),DIMENSION(:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dil2d=-1
      WRITE(6,*)'array not allocated dil2d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dil2d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dil2d)
      IF(dil2d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dil2d'
        WRITE(0,*)'error in deallocation routine dil2d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dil2d

  !---------------------------------

  INTEGER FUNCTION ail3d(arr,d1,d2,d3,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ilong),DIMENSION(:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1_ilong)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ail3d=-1
      WRITE(6,*)'array already allocated ail3d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ail3d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ail3d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3),STAT=ail3d)
        IF(ail3d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ail3d'
          WRITE(0,*)'error in allocation routine ail3d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ail3d

  !---------------------------------

  INTEGER FUNCTION dil3d(arr,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ilong),DIMENSION(:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dil3d=-1
      WRITE(6,*)'array not allocated dil3d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dil3d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dil3d)
      IF(dil3d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dil3d'
        WRITE(0,*)'error in deallocation routine dil3d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dil3d

  !---------------------------------

  INTEGER FUNCTION ail4d(arr,d1,d2,d3,d4,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ilong),DIMENSION(:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1_ilong)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ail4d=-1
      WRITE(6,*)'array already allocated ail4d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ail4d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ail4d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4),STAT=ail4d)
        IF(ail4d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ail4d'
          WRITE(0,*)'error in allocation routine ail4d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ail4d

  !---------------------------------

  INTEGER FUNCTION dil4d(arr,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ilong),DIMENSION(:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dil4d=-1
      WRITE(6,*)'array not allocated dil4d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dil4d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dil4d)
      IF(dil4d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dil4d'
        WRITE(0,*)'error in deallocation routine dil4d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dil4d

  !---------------------------------

  INTEGER FUNCTION ail5d(arr,d1,d2,d3,d4,d5,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ilong),DIMENSION(:,:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    INTEGER(KIND=ilong), INTENT(IN) ::d5
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1_ilong)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ail5d=-1
      WRITE(6,*)'array already allocated ail5d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ail5d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4*d5
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ail5d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4,1:d5),STAT=ail5d)
        IF(ail5d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ail5d'
          WRITE(0,*)'error in allocation routine ail5d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ail5d

  !---------------------------------

  INTEGER FUNCTION dil5d(arr,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ilong),DIMENSION(:,:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dil5d=-1
      WRITE(6,*)'array not allocated dil5d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dil5d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dil5d)
      IF(dil5d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dil5d'
        WRITE(0,*)'error in deallocation routine dil5d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dil5d

  !---------------------------------

  INTEGER FUNCTION ail6d(arr,d1,d2,d3,d4,d5,d6,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ilong),DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    INTEGER(KIND=ilong), INTENT(IN) ::d5
    INTEGER(KIND=ilong), INTENT(IN) ::d6
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1_ilong)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ail6d=-1
      WRITE(6,*)'array already allocated ail6d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ail6d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4*d5*d6
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ail6d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4,1:d5,1:d6),STAT=ail6d)
        IF(ail6d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ail6d'
          WRITE(0,*)'error in allocation routine ail6d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ail6d

  !---------------------------------

  INTEGER FUNCTION dil6d(arr,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ilong),DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dil6d=-1
      WRITE(6,*)'array not allocated dil6d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dil6d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dil6d)
      IF(dil6d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dil6d'
        WRITE(0,*)'error in deallocation routine dil6d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dil6d

  !---------------------------------

  INTEGER FUNCTION ail7d(arr,d1,d2,d3,d4,d5,d6,d7,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ilong),DIMENSION(:,:,:,:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    INTEGER(KIND=ilong), INTENT(IN) ::d5
    INTEGER(KIND=ilong), INTENT(IN) ::d6
    INTEGER(KIND=ilong), INTENT(IN) ::d7
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1_ilong)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ail7d=-1
      WRITE(6,*)'array already allocated ail7d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ail7d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4*d5*d6*d7
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ail7d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4,1:d5,1:d6,1:d7),STAT=ail7d)
        IF(ail7d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ail7d'
          WRITE(0,*)'error in allocation routine ail7d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ail7d

  !---------------------------------

  INTEGER FUNCTION dil7d(arr,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ilong),DIMENSION(:,:,:,:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dil7d=-1
      WRITE(6,*)'array not allocated dil7d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dil7d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dil7d)
      IF(dil7d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dil7d'
        WRITE(0,*)'error in deallocation routine dil7d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dil7d

  !---------------------------------

  INTEGER FUNCTION ais1d(arr,d1,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ishort),DIMENSION(:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1_ishort)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ais1d=-1
      WRITE(6,*)'array already allocated ais1d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ais1d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ais1d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1),STAT=ais1d)
        IF(ais1d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ais1d'
          WRITE(0,*)'error in allocation routine ais1d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ais1d

  !---------------------------------

  INTEGER FUNCTION dis1d(arr,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ishort),DIMENSION(:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dis1d=-1
      WRITE(6,*)'array not allocated dis1d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dis1d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dis1d)
      IF(dis1d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dis1d'
        WRITE(0,*)'error in deallocation routine dis1d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dis1d

  !---------------------------------

  INTEGER FUNCTION ais2d(arr,d1,d2,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ishort),DIMENSION(:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1_ishort)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ais2d=-1
      WRITE(6,*)'array already allocated ais2d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ais2d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ais2d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2),STAT=ais2d)
        IF(ais2d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ais2d'
          WRITE(0,*)'error in allocation routine ais2d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ais2d

  !---------------------------------

  INTEGER FUNCTION dis2d(arr,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ishort),DIMENSION(:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dis2d=-1
      WRITE(6,*)'array not allocated dis2d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dis2d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dis2d)
      IF(dis2d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dis2d'
        WRITE(0,*)'error in deallocation routine dis2d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dis2d

  !---------------------------------

  INTEGER FUNCTION ais3d(arr,d1,d2,d3,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ishort),DIMENSION(:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1_ishort)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ais3d=-1
      WRITE(6,*)'array already allocated ais3d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ais3d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ais3d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3),STAT=ais3d)
        IF(ais3d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ais3d'
          WRITE(0,*)'error in allocation routine ais3d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ais3d

  !---------------------------------

  INTEGER FUNCTION dis3d(arr,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ishort),DIMENSION(:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dis3d=-1
      WRITE(6,*)'array not allocated dis3d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dis3d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dis3d)
      IF(dis3d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dis3d'
        WRITE(0,*)'error in deallocation routine dis3d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dis3d

  !---------------------------------

  INTEGER FUNCTION ais4d(arr,d1,d2,d3,d4,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ishort),DIMENSION(:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1_ishort)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ais4d=-1
      WRITE(6,*)'array already allocated ais4d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ais4d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ais4d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4),STAT=ais4d)
        IF(ais4d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ais4d'
          WRITE(0,*)'error in allocation routine ais4d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ais4d

  !---------------------------------

  INTEGER FUNCTION dis4d(arr,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ishort),DIMENSION(:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dis4d=-1
      WRITE(6,*)'array not allocated dis4d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dis4d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dis4d)
      IF(dis4d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dis4d'
        WRITE(0,*)'error in deallocation routine dis4d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dis4d

  !---------------------------------

  INTEGER FUNCTION ais5d(arr,d1,d2,d3,d4,d5,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ishort),DIMENSION(:,:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    INTEGER(KIND=ilong), INTENT(IN) ::d5
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1_ishort)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ais5d=-1
      WRITE(6,*)'array already allocated ais5d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ais5d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4*d5
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ais5d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4,1:d5),STAT=ais5d)
        IF(ais5d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ais5d'
          WRITE(0,*)'error in allocation routine ais5d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ais5d

  !---------------------------------

  INTEGER FUNCTION dis5d(arr,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ishort),DIMENSION(:,:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dis5d=-1
      WRITE(6,*)'array not allocated dis5d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dis5d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dis5d)
      IF(dis5d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dis5d'
        WRITE(0,*)'error in deallocation routine dis5d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dis5d

  !---------------------------------

  INTEGER FUNCTION ais6d(arr,d1,d2,d3,d4,d5,d6,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ishort),DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    INTEGER(KIND=ilong), INTENT(IN) ::d5
    INTEGER(KIND=ilong), INTENT(IN) ::d6
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1_ishort)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ais6d=-1
      WRITE(6,*)'array already allocated ais6d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ais6d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4*d5*d6
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ais6d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4,1:d5,1:d6),STAT=ais6d)
        IF(ais6d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ais6d'
          WRITE(0,*)'error in allocation routine ais6d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ais6d

  !---------------------------------

  INTEGER FUNCTION dis6d(arr,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ishort),DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dis6d=-1
      WRITE(6,*)'array not allocated dis6d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dis6d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dis6d)
      IF(dis6d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dis6d'
        WRITE(0,*)'error in deallocation routine dis6d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dis6d

  !---------------------------------

  INTEGER FUNCTION ais7d(arr,d1,d2,d3,d4,d5,d6,d7,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ishort),DIMENSION(:,:,:,:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    INTEGER(KIND=ilong), INTENT(IN) ::d5
    INTEGER(KIND=ilong), INTENT(IN) ::d6
    INTEGER(KIND=ilong), INTENT(IN) ::d7
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1_ishort)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ais7d=-1
      WRITE(6,*)'array already allocated ais7d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ais7d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4*d5*d6*d7
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ais7d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4,1:d5,1:d6,1:d7),STAT=ais7d)
        IF(ais7d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ais7d'
          WRITE(0,*)'error in allocation routine ais7d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ais7d

  !---------------------------------

  INTEGER FUNCTION dis7d(arr,msg,thrd)
    IMPLICIT NONE
    INTEGER(KIND=ishort),DIMENSION(:,:,:,:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dis7d=-1
      WRITE(6,*)'array not allocated dis7d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dis7d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dis7d)
      IF(dis7d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dis7d'
        WRITE(0,*)'error in deallocation routine dis7d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dis7d

  !---------------------------------

  INTEGER FUNCTION ar1d(arr,d1,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=sip),DIMENSION(:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1.0_sip)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ar1d=-1
      WRITE(6,*)'array already allocated ar1d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ar1d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ar1d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1),STAT=ar1d)
        IF(ar1d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ar1d'
          WRITE(0,*)'error in allocation routine ar1d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ar1d

  !---------------------------------

  INTEGER FUNCTION dr1d(arr,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=sip),DIMENSION(:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dr1d=-1
      WRITE(6,*)'array not allocated dr1d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dr1d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dr1d)
      IF(dr1d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dr1d'
        WRITE(0,*)'error in deallocation routine dr1d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dr1d

  !---------------------------------

  INTEGER FUNCTION ar2d(arr,d1,d2,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=sip),DIMENSION(:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1.0_sip)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ar2d=-1
      WRITE(6,*)'array already allocated ar2d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ar2d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ar2d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2),STAT=ar2d)
        IF(ar2d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ar2d'
          WRITE(0,*)'error in allocation routine ar2d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ar2d

  !---------------------------------

  INTEGER FUNCTION dr2d(arr,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=sip),DIMENSION(:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dr2d=-1
      WRITE(6,*)'array not allocated dr2d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dr2d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dr2d)
      IF(dr2d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dr2d'
        WRITE(0,*)'error in deallocation routine dr2d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dr2d

  !---------------------------------

  INTEGER FUNCTION ar3d(arr,d1,d2,d3,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=sip),DIMENSION(:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1.0_sip)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ar3d=-1
      WRITE(6,*)'array already allocated ar3d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ar3d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ar3d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3),STAT=ar3d)
        IF(ar3d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ar3d'
          WRITE(0,*)'error in allocation routine ar3d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ar3d

  !---------------------------------

  INTEGER FUNCTION dr3d(arr,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=sip),DIMENSION(:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dr3d=-1
      WRITE(6,*)'array not allocated dr3d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dr3d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dr3d)
      IF(dr3d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dr3d'
        WRITE(0,*)'error in deallocation routine dr3d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dr3d

  !---------------------------------

  INTEGER FUNCTION ar4d(arr,d1,d2,d3,d4,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=sip),DIMENSION(:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1.0_sip)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ar4d=-1
      WRITE(6,*)'array already allocated ar4d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ar4d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ar4d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4),STAT=ar4d)
        IF(ar4d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ar4d'
          WRITE(0,*)'error in allocation routine ar4d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ar4d

  !---------------------------------

  INTEGER FUNCTION dr4d(arr,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=sip),DIMENSION(:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dr4d=-1
      WRITE(6,*)'array not allocated dr4d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dr4d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dr4d)
      IF(dr4d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dr4d'
        WRITE(0,*)'error in deallocation routine dr4d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dr4d

  !---------------------------------

  INTEGER FUNCTION ar5d(arr,d1,d2,d3,d4,d5,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=sip),DIMENSION(:,:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    INTEGER(KIND=ilong), INTENT(IN) ::d5
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1.0_sip)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ar5d=-1
      WRITE(6,*)'array already allocated ar5d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ar5d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4*d5
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ar5d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4,1:d5),STAT=ar5d)
        IF(ar5d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ar5d'
          WRITE(0,*)'error in allocation routine ar5d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ar5d

  !---------------------------------

  INTEGER FUNCTION dr5d(arr,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=sip),DIMENSION(:,:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dr5d=-1
      WRITE(6,*)'array not allocated dr5d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dr5d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dr5d)
      IF(dr5d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dr5d'
        WRITE(0,*)'error in deallocation routine dr5d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dr5d

  !---------------------------------

  INTEGER FUNCTION ar6d(arr,d1,d2,d3,d4,d5,d6,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=sip),DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    INTEGER(KIND=ilong), INTENT(IN) ::d5
    INTEGER(KIND=ilong), INTENT(IN) ::d6
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1.0_sip)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ar6d=-1
      WRITE(6,*)'array already allocated ar6d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ar6d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4*d5*d6
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ar6d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4,1:d5,1:d6),STAT=ar6d)
        IF(ar6d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ar6d'
          WRITE(0,*)'error in allocation routine ar6d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ar6d

  !---------------------------------

  INTEGER FUNCTION dr6d(arr,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=sip),DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dr6d=-1
      WRITE(6,*)'array not allocated dr6d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dr6d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dr6d)
      IF(dr6d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dr6d'
        WRITE(0,*)'error in deallocation routine dr6d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dr6d

  !---------------------------------

  INTEGER FUNCTION ar7d(arr,d1,d2,d3,d4,d5,d6,d7,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=sip),DIMENSION(:,:,:,:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    INTEGER(KIND=ilong), INTENT(IN) ::d5
    INTEGER(KIND=ilong), INTENT(IN) ::d6
    INTEGER(KIND=ilong), INTENT(IN) ::d7
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1.0_sip)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ar7d=-1
      WRITE(6,*)'array already allocated ar7d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ar7d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4*d5*d6*d7
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ar7d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4,1:d5,1:d6,1:d7),STAT=ar7d)
        IF(ar7d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ar7d'
          WRITE(0,*)'error in allocation routine ar7d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ar7d

  !---------------------------------

  INTEGER FUNCTION dr7d(arr,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=sip),DIMENSION(:,:,:,:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dr7d=-1
      WRITE(6,*)'array not allocated dr7d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dr7d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dr7d)
      IF(dr7d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dr7d'
        WRITE(0,*)'error in deallocation routine dr7d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dr7d

  !---------------------------------

  INTEGER FUNCTION ad1d(arr,d1,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=dop),DIMENSION(:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1.0_dop)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ad1d=-1
      WRITE(6,*)'array already allocated ad1d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ad1d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ad1d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1),STAT=ad1d)
        IF(ad1d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ad1d'
          WRITE(0,*)'error in allocation routine ad1d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ad1d

  !---------------------------------

  INTEGER FUNCTION dd1d(arr,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=dop),DIMENSION(:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dd1d=-1
      WRITE(6,*)'array not allocated dd1d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dd1d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dd1d)
      IF(dd1d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dd1d'
        WRITE(0,*)'error in deallocation routine dd1d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dd1d

  !---------------------------------

  INTEGER FUNCTION ad2d(arr,d1,d2,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=dop),DIMENSION(:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1.0_dop)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ad2d=-1
      WRITE(6,*)'array already allocated ad2d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ad2d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ad2d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2),STAT=ad2d)
        IF(ad2d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ad2d'
          WRITE(0,*)'error in allocation routine ad2d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ad2d

  !---------------------------------

  INTEGER FUNCTION dd2d(arr,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=dop),DIMENSION(:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dd2d=-1
      WRITE(6,*)'array not allocated dd2d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dd2d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dd2d)
      IF(dd2d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dd2d'
        WRITE(0,*)'error in deallocation routine dd2d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dd2d

  !---------------------------------

  INTEGER FUNCTION ad3d(arr,d1,d2,d3,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=dop),DIMENSION(:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1.0_dop)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ad3d=-1
      WRITE(6,*)'array already allocated ad3d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ad3d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ad3d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3),STAT=ad3d)
        IF(ad3d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ad3d'
          WRITE(0,*)'error in allocation routine ad3d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ad3d

  !---------------------------------

  INTEGER FUNCTION dd3d(arr,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=dop),DIMENSION(:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dd3d=-1
      WRITE(6,*)'array not allocated dd3d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dd3d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dd3d)
      IF(dd3d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dd3d'
        WRITE(0,*)'error in deallocation routine dd3d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dd3d

  !---------------------------------

  INTEGER FUNCTION ad4d(arr,d1,d2,d3,d4,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=dop),DIMENSION(:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1.0_dop)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ad4d=-1
      WRITE(6,*)'array already allocated ad4d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ad4d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ad4d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4),STAT=ad4d)
        IF(ad4d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ad4d'
          WRITE(0,*)'error in allocation routine ad4d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ad4d

  !---------------------------------

  INTEGER FUNCTION dd4d(arr,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=dop),DIMENSION(:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dd4d=-1
      WRITE(6,*)'array not allocated dd4d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dd4d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dd4d)
      IF(dd4d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dd4d'
        WRITE(0,*)'error in deallocation routine dd4d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dd4d

  !---------------------------------

  INTEGER FUNCTION ad5d(arr,d1,d2,d3,d4,d5,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=dop),DIMENSION(:,:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    INTEGER(KIND=ilong), INTENT(IN) ::d5
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1.0_dop)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ad5d=-1
      WRITE(6,*)'array already allocated ad5d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ad5d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4*d5
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ad5d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4,1:d5),STAT=ad5d)
        IF(ad5d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ad5d'
          WRITE(0,*)'error in allocation routine ad5d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ad5d

  !---------------------------------

  INTEGER FUNCTION dd5d(arr,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=dop),DIMENSION(:,:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dd5d=-1
      WRITE(6,*)'array not allocated dd5d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dd5d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dd5d)
      IF(dd5d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dd5d'
        WRITE(0,*)'error in deallocation routine dd5d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dd5d

  !---------------------------------

  INTEGER FUNCTION ad6d(arr,d1,d2,d3,d4,d5,d6,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=dop),DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    INTEGER(KIND=ilong), INTENT(IN) ::d5
    INTEGER(KIND=ilong), INTENT(IN) ::d6
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1.0_dop)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ad6d=-1
      WRITE(6,*)'array already allocated ad6d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ad6d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4*d5*d6
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ad6d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4,1:d5,1:d6),STAT=ad6d)
        IF(ad6d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ad6d'
          WRITE(0,*)'error in allocation routine ad6d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ad6d

  !---------------------------------

  INTEGER FUNCTION dd6d(arr,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=dop),DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dd6d=-1
      WRITE(6,*)'array not allocated dd6d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dd6d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dd6d)
      IF(dd6d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dd6d'
        WRITE(0,*)'error in deallocation routine dd6d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dd6d

  !---------------------------------

  INTEGER FUNCTION ad7d(arr,d1,d2,d3,d4,d5,d6,d7,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=dop),DIMENSION(:,:,:,:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    INTEGER(KIND=ilong), INTENT(IN) ::d5
    INTEGER(KIND=ilong), INTENT(IN) ::d6
    INTEGER(KIND=ilong), INTENT(IN) ::d7
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(1.0_dop)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ad7d=-1
      WRITE(6,*)'array already allocated ad7d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ad7d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4*d5*d6*d7
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ad7d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4,1:d5,1:d6,1:d7),STAT=ad7d)
        IF(ad7d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ad7d'
          WRITE(0,*)'error in allocation routine ad7d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ad7d

  !---------------------------------

  INTEGER FUNCTION dd7d(arr,msg,thrd)
    IMPLICIT NONE
    REAL(KIND=dop),DIMENSION(:,:,:,:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dd7d=-1
      WRITE(6,*)'array not allocated dd7d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dd7d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dd7d)
      IF(dd7d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dd7d'
        WRITE(0,*)'error in deallocation routine dd7d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dd7d

  !---------------------------------

  INTEGER FUNCTION ac1d(arr,d1,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=sip),DIMENSION(:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=2*sizeof(1.0_sip)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ac1d=-1
      WRITE(6,*)'array already allocated ac1d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ac1d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ac1d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1),STAT=ac1d)
        IF(ac1d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ac1d'
          WRITE(0,*)'error in allocation routine ac1d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ac1d

  !---------------------------------

  INTEGER FUNCTION dc1d(arr,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=sip),DIMENSION(:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dc1d=-1
      WRITE(6,*)'array not allocated dc1d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dc1d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dc1d)
      IF(dc1d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dc1d'
        WRITE(0,*)'error in deallocation routine dc1d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dc1d

  !---------------------------------

  INTEGER FUNCTION ac2d(arr,d1,d2,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=sip),DIMENSION(:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=2*sizeof(1.0_sip)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ac2d=-1
      WRITE(6,*)'array already allocated ac2d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ac2d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ac2d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2),STAT=ac2d)
        IF(ac2d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ac2d'
          WRITE(0,*)'error in allocation routine ac2d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ac2d

  !---------------------------------

  INTEGER FUNCTION dc2d(arr,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=sip),DIMENSION(:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dc2d=-1
      WRITE(6,*)'array not allocated dc2d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dc2d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dc2d)
      IF(dc2d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dc2d'
        WRITE(0,*)'error in deallocation routine dc2d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dc2d

  !---------------------------------

  INTEGER FUNCTION ac3d(arr,d1,d2,d3,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=sip),DIMENSION(:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=2*sizeof(1.0_sip)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ac3d=-1
      WRITE(6,*)'array already allocated ac3d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ac3d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ac3d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3),STAT=ac3d)
        IF(ac3d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ac3d'
          WRITE(0,*)'error in allocation routine ac3d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ac3d

  !---------------------------------

  INTEGER FUNCTION dc3d(arr,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=sip),DIMENSION(:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dc3d=-1
      WRITE(6,*)'array not allocated dc3d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dc3d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dc3d)
      IF(dc3d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dc3d'
        WRITE(0,*)'error in deallocation routine dc3d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dc3d

  !---------------------------------

  INTEGER FUNCTION ac4d(arr,d1,d2,d3,d4,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=sip),DIMENSION(:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=2*sizeof(1.0_sip)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ac4d=-1
      WRITE(6,*)'array already allocated ac4d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ac4d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ac4d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4),STAT=ac4d)
        IF(ac4d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ac4d'
          WRITE(0,*)'error in allocation routine ac4d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ac4d

  !---------------------------------

  INTEGER FUNCTION dc4d(arr,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=sip),DIMENSION(:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dc4d=-1
      WRITE(6,*)'array not allocated dc4d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dc4d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dc4d)
      IF(dc4d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dc4d'
        WRITE(0,*)'error in deallocation routine dc4d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dc4d

  !---------------------------------

  INTEGER FUNCTION ac5d(arr,d1,d2,d3,d4,d5,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=sip),DIMENSION(:,:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    INTEGER(KIND=ilong), INTENT(IN) ::d5
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=2*sizeof(1.0_sip)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ac5d=-1
      WRITE(6,*)'array already allocated ac5d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ac5d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4*d5
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ac5d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4,1:d5),STAT=ac5d)
        IF(ac5d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ac5d'
          WRITE(0,*)'error in allocation routine ac5d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ac5d

  !---------------------------------

  INTEGER FUNCTION dc5d(arr,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=sip),DIMENSION(:,:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dc5d=-1
      WRITE(6,*)'array not allocated dc5d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dc5d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dc5d)
      IF(dc5d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dc5d'
        WRITE(0,*)'error in deallocation routine dc5d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dc5d

  !---------------------------------

  INTEGER FUNCTION ac6d(arr,d1,d2,d3,d4,d5,d6,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=sip),DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    INTEGER(KIND=ilong), INTENT(IN) ::d5
    INTEGER(KIND=ilong), INTENT(IN) ::d6
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=2*sizeof(1.0_sip)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ac6d=-1
      WRITE(6,*)'array already allocated ac6d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ac6d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4*d5*d6
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ac6d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4,1:d5,1:d6),STAT=ac6d)
        IF(ac6d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ac6d'
          WRITE(0,*)'error in allocation routine ac6d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ac6d

  !---------------------------------

  INTEGER FUNCTION dc6d(arr,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=sip),DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dc6d=-1
      WRITE(6,*)'array not allocated dc6d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dc6d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dc6d)
      IF(dc6d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dc6d'
        WRITE(0,*)'error in deallocation routine dc6d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dc6d

  !---------------------------------

  INTEGER FUNCTION ac7d(arr,d1,d2,d3,d4,d5,d6,d7,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=sip),DIMENSION(:,:,:,:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    INTEGER(KIND=ilong), INTENT(IN) ::d5
    INTEGER(KIND=ilong), INTENT(IN) ::d6
    INTEGER(KIND=ilong), INTENT(IN) ::d7
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=2*sizeof(1.0_sip)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      ac7d=-1
      WRITE(6,*)'array already allocated ac7d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated ac7d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4*d5*d6*d7
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        ac7d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4,1:d5,1:d6,1:d7),STAT=ac7d)
        IF(ac7d.NE.0)THEN
          WRITE(6,*)'error in allocation routine ac7d'
          WRITE(0,*)'error in allocation routine ac7d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION ac7d

  !---------------------------------

  INTEGER FUNCTION dc7d(arr,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=sip),DIMENSION(:,:,:,:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dc7d=-1
      WRITE(6,*)'array not allocated dc7d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dc7d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dc7d)
      IF(dc7d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dc7d'
        WRITE(0,*)'error in deallocation routine dc7d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dc7d

  !---------------------------------

  INTEGER FUNCTION az1d(arr,d1,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=dop),DIMENSION(:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=2*sizeof(1.0_dop)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      az1d=-1
      WRITE(6,*)'array already allocated az1d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated az1d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        az1d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1),STAT=az1d)
        IF(az1d.NE.0)THEN
          WRITE(6,*)'error in allocation routine az1d'
          WRITE(0,*)'error in allocation routine az1d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION az1d

  !---------------------------------

  INTEGER FUNCTION dz1d(arr,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=dop),DIMENSION(:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dz1d=-1
      WRITE(6,*)'array not allocated dz1d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dz1d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dz1d)
      IF(dz1d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dz1d'
        WRITE(0,*)'error in deallocation routine dz1d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dz1d

  !---------------------------------

  INTEGER FUNCTION az2d(arr,d1,d2,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=dop),DIMENSION(:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=2*sizeof(1.0_dop)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      az2d=-1
      WRITE(6,*)'array already allocated az2d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated az2d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        az2d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2),STAT=az2d)
        IF(az2d.NE.0)THEN
          WRITE(6,*)'error in allocation routine az2d'
          WRITE(0,*)'error in allocation routine az2d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION az2d

  !---------------------------------

  INTEGER FUNCTION dz2d(arr,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=dop),DIMENSION(:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dz2d=-1
      WRITE(6,*)'array not allocated dz2d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dz2d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dz2d)
      IF(dz2d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dz2d'
        WRITE(0,*)'error in deallocation routine dz2d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dz2d

  !---------------------------------

  INTEGER FUNCTION az3d(arr,d1,d2,d3,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=dop),DIMENSION(:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=2*sizeof(1.0_dop)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      az3d=-1
      WRITE(6,*)'array already allocated az3d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated az3d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        az3d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3),STAT=az3d)
        IF(az3d.NE.0)THEN
          WRITE(6,*)'error in allocation routine az3d'
          WRITE(0,*)'error in allocation routine az3d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION az3d

  !---------------------------------

  INTEGER FUNCTION dz3d(arr,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=dop),DIMENSION(:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dz3d=-1
      WRITE(6,*)'array not allocated dz3d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dz3d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dz3d)
      IF(dz3d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dz3d'
        WRITE(0,*)'error in deallocation routine dz3d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dz3d

  !---------------------------------

  INTEGER FUNCTION az4d(arr,d1,d2,d3,d4,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=dop),DIMENSION(:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=2*sizeof(1.0_dop)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      az4d=-1
      WRITE(6,*)'array already allocated az4d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated az4d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        az4d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4),STAT=az4d)
        IF(az4d.NE.0)THEN
          WRITE(6,*)'error in allocation routine az4d'
          WRITE(0,*)'error in allocation routine az4d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION az4d

  !---------------------------------

  INTEGER FUNCTION dz4d(arr,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=dop),DIMENSION(:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dz4d=-1
      WRITE(6,*)'array not allocated dz4d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dz4d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dz4d)
      IF(dz4d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dz4d'
        WRITE(0,*)'error in deallocation routine dz4d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dz4d

  !---------------------------------

  INTEGER FUNCTION az5d(arr,d1,d2,d3,d4,d5,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=dop),DIMENSION(:,:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    INTEGER(KIND=ilong), INTENT(IN) ::d5
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=2*sizeof(1.0_dop)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      az5d=-1
      WRITE(6,*)'array already allocated az5d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated az5d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4*d5
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        az5d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4,1:d5),STAT=az5d)
        IF(az5d.NE.0)THEN
          WRITE(6,*)'error in allocation routine az5d'
          WRITE(0,*)'error in allocation routine az5d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION az5d

  !---------------------------------

  INTEGER FUNCTION dz5d(arr,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=dop),DIMENSION(:,:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dz5d=-1
      WRITE(6,*)'array not allocated dz5d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dz5d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dz5d)
      IF(dz5d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dz5d'
        WRITE(0,*)'error in deallocation routine dz5d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dz5d

  !---------------------------------

  INTEGER FUNCTION az6d(arr,d1,d2,d3,d4,d5,d6,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=dop),DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    INTEGER(KIND=ilong), INTENT(IN) ::d5
    INTEGER(KIND=ilong), INTENT(IN) ::d6
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=2*sizeof(1.0_dop)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      az6d=-1
      WRITE(6,*)'array already allocated az6d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated az6d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4*d5*d6
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        az6d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4,1:d5,1:d6),STAT=az6d)
        IF(az6d.NE.0)THEN
          WRITE(6,*)'error in allocation routine az6d'
          WRITE(0,*)'error in allocation routine az6d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION az6d

  !---------------------------------

  INTEGER FUNCTION dz6d(arr,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=dop),DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dz6d=-1
      WRITE(6,*)'array not allocated dz6d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dz6d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dz6d)
      IF(dz6d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dz6d'
        WRITE(0,*)'error in deallocation routine dz6d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dz6d

  !---------------------------------

  INTEGER FUNCTION az7d(arr,d1,d2,d3,d4,d5,d6,d7,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=dop),DIMENSION(:,:,:,:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    INTEGER(KIND=ilong), INTENT(IN) ::d5
    INTEGER(KIND=ilong), INTENT(IN) ::d6
    INTEGER(KIND=ilong), INTENT(IN) ::d7
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=2*sizeof(1.0_dop)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      az7d=-1
      WRITE(6,*)'array already allocated az7d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated az7d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4*d5*d6*d7
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        az7d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4,1:d5,1:d6,1:d7),STAT=az7d)
        IF(az7d.NE.0)THEN
          WRITE(6,*)'error in allocation routine az7d'
          WRITE(0,*)'error in allocation routine az7d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION az7d

  !---------------------------------

  INTEGER FUNCTION dz7d(arr,msg,thrd)
    IMPLICIT NONE
    COMPLEX(KIND=dop),DIMENSION(:,:,:,:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dz7d=-1
      WRITE(6,*)'array not allocated dz7d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dz7d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dz7d)
      IF(dz7d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dz7d'
        WRITE(0,*)'error in deallocation routine dz7d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dz7d

  !---------------------------------

  INTEGER FUNCTION al1d(arr,d1,msg,thrd)
    IMPLICIT NONE
    LOGICAL,DIMENSION(:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(.TRUE.)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      al1d=-1
      WRITE(6,*)'array already allocated al1d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated al1d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        al1d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1),STAT=al1d)
        IF(al1d.NE.0)THEN
          WRITE(6,*)'error in allocation routine al1d'
          WRITE(0,*)'error in allocation routine al1d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION al1d

  !---------------------------------

  INTEGER FUNCTION dl1d(arr,msg,thrd)
    IMPLICIT NONE
    LOGICAL,DIMENSION(:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dl1d=-1
      WRITE(6,*)'array not allocated dl1d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dl1d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dl1d)
      IF(dl1d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dl1d'
        WRITE(0,*)'error in deallocation routine dl1d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dl1d

  !---------------------------------

  INTEGER FUNCTION al2d(arr,d1,d2,msg,thrd)
    IMPLICIT NONE
    LOGICAL,DIMENSION(:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(.TRUE.)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      al2d=-1
      WRITE(6,*)'array already allocated al2d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated al2d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        al2d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2),STAT=al2d)
        IF(al2d.NE.0)THEN
          WRITE(6,*)'error in allocation routine al2d'
          WRITE(0,*)'error in allocation routine al2d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION al2d

  !---------------------------------

  INTEGER FUNCTION dl2d(arr,msg,thrd)
    IMPLICIT NONE
    LOGICAL,DIMENSION(:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dl2d=-1
      WRITE(6,*)'array not allocated dl2d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dl2d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dl2d)
      IF(dl2d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dl2d'
        WRITE(0,*)'error in deallocation routine dl2d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dl2d

  !---------------------------------

  INTEGER FUNCTION al3d(arr,d1,d2,d3,msg,thrd)
    IMPLICIT NONE
    LOGICAL,DIMENSION(:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(.TRUE.)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      al3d=-1
      WRITE(6,*)'array already allocated al3d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated al3d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        al3d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3),STAT=al3d)
        IF(al3d.NE.0)THEN
          WRITE(6,*)'error in allocation routine al3d'
          WRITE(0,*)'error in allocation routine al3d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION al3d

  !---------------------------------

  INTEGER FUNCTION dl3d(arr,msg,thrd)
    IMPLICIT NONE
    LOGICAL,DIMENSION(:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dl3d=-1
      WRITE(6,*)'array not allocated dl3d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dl3d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dl3d)
      IF(dl3d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dl3d'
        WRITE(0,*)'error in deallocation routine dl3d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dl3d

  !---------------------------------

  INTEGER FUNCTION al4d(arr,d1,d2,d3,d4,msg,thrd)
    IMPLICIT NONE
    LOGICAL,DIMENSION(:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(.TRUE.)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      al4d=-1
      WRITE(6,*)'array already allocated al4d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated al4d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        al4d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4),STAT=al4d)
        IF(al4d.NE.0)THEN
          WRITE(6,*)'error in allocation routine al4d'
          WRITE(0,*)'error in allocation routine al4d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION al4d

  !---------------------------------

  INTEGER FUNCTION dl4d(arr,msg,thrd)
    IMPLICIT NONE
    LOGICAL,DIMENSION(:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dl4d=-1
      WRITE(6,*)'array not allocated dl4d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dl4d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dl4d)
      IF(dl4d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dl4d'
        WRITE(0,*)'error in deallocation routine dl4d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dl4d

  !---------------------------------

  INTEGER FUNCTION al5d(arr,d1,d2,d3,d4,d5,msg,thrd)
    IMPLICIT NONE
    LOGICAL,DIMENSION(:,:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    INTEGER(KIND=ilong), INTENT(IN) ::d5
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(.TRUE.)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      al5d=-1
      WRITE(6,*)'array already allocated al5d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated al5d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4*d5
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        al5d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4,1:d5),STAT=al5d)
        IF(al5d.NE.0)THEN
          WRITE(6,*)'error in allocation routine al5d'
          WRITE(0,*)'error in allocation routine al5d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION al5d

  !---------------------------------

  INTEGER FUNCTION dl5d(arr,msg,thrd)
    IMPLICIT NONE
    LOGICAL,DIMENSION(:,:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dl5d=-1
      WRITE(6,*)'array not allocated dl5d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dl5d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dl5d)
      IF(dl5d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dl5d'
        WRITE(0,*)'error in deallocation routine dl5d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dl5d

  !---------------------------------

  INTEGER FUNCTION al6d(arr,d1,d2,d3,d4,d5,d6,msg,thrd)
    IMPLICIT NONE
    LOGICAL,DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    INTEGER(KIND=ilong), INTENT(IN) ::d5
    INTEGER(KIND=ilong), INTENT(IN) ::d6
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(.TRUE.)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      al6d=-1
      WRITE(6,*)'array already allocated al6d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated al6d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4*d5*d6
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        al6d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4,1:d5,1:d6),STAT=al6d)
        IF(al6d.NE.0)THEN
          WRITE(6,*)'error in allocation routine al6d'
          WRITE(0,*)'error in allocation routine al6d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION al6d

  !---------------------------------

  INTEGER FUNCTION dl6d(arr,msg,thrd)
    IMPLICIT NONE
    LOGICAL,DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dl6d=-1
      WRITE(6,*)'array not allocated dl6d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dl6d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dl6d)
      IF(dl6d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dl6d'
        WRITE(0,*)'error in deallocation routine dl6d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dl6d

  !---------------------------------

  INTEGER FUNCTION al7d(arr,d1,d2,d3,d4,d5,d6,d7,msg,thrd)
    IMPLICIT NONE
    LOGICAL,DIMENSION(:,:,:,:,:,:,:),ALLOCATABLE :: arr
    INTEGER(KIND=ilong), INTENT(IN) :: &
      d1
    INTEGER(KIND=ilong), INTENT(IN) ::d2
    INTEGER(KIND=ilong), INTENT(IN) ::d3
    INTEGER(KIND=ilong), INTENT(IN) ::d4
    INTEGER(KIND=ilong), INTENT(IN) ::d5
    INTEGER(KIND=ilong), INTENT(IN) ::d6
    INTEGER(KIND=ilong), INTENT(IN) ::d7
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    INTEGER(KIND=ilong):: sizeofone
    sizeofone=sizeof(.TRUE.)
    nm='+'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(allocated(arr))THEN
      al7d=-1
      WRITE(6,*)'array already allocated al7d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array already allocated al7d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeofone*d1*d2*d3*d4*d5*d6*d7
      IF((maxmem.GT.0).AND.((mem_used+memu).GT.maxmem)) THEN
        al7d=-1
        WRITE(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)
111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")
      ELSE
        ALLOCATE(arr(1:d1,1:d2,1:d3,1:d4,1:d5,1:d6,1:d7),STAT=al7d)
        IF(al7d.NE.0)THEN
          WRITE(6,*)'error in allocation routine al7d'
          WRITE(0,*)'error in allocation routine al7d'
        ELSE
          mem_used=mem_used+memu
          IF(present(thrd))THEN
            CALL log_memory(nm,thrd)
          ELSE
            IF(present(msg))THEN
              CALL log_memory(nm)
            END IF
          END IF
        END IF
      END IF
    END IF
  END FUNCTION al7d

  !---------------------------------

  INTEGER FUNCTION dl7d(arr,msg,thrd)
    IMPLICIT NONE
    LOGICAL,DIMENSION(:,:,:,:,:,:,:),ALLOCATABLE :: arr
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg
    INTEGER, INTENT(IN), OPTIONAL :: thrd
    CHARACTER*30:: nm

    nm='-'
    IF(present(msg)) THEN
      nm=msg
    END IF
    IF(.NOT.allocated(arr))THEN
      dl7d=-1
      WRITE(6,*)'array not allocated dl7d (',trim(adjustl(nm)),')'
      WRITE(0,*)'array not allocated dl7d (',trim(adjustl(nm)),')'
    ELSE
      memu=sizeof(arr)
      DEALLOCATE(arr,STAT=dl7d)
      IF(dl7d.NE.0)THEN
        WRITE(6,*)'error in deallocation routine dl7d'
        WRITE(0,*)'error in deallocation routine dl7d'
      ELSE
        mem_used=mem_used-memu
        IF(present(thrd))THEN
          CALL log_memory(nm,thrd)
        ELSE
          IF(present(msg))THEN
            CALL log_memory(nm)
          END IF
        END IF
      END IF
    END IF
  END FUNCTION dl7d

  !---------------------------------

END MODULE my_alloc
