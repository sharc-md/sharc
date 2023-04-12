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

MODULE sortblockmod
  !> module for reading the input
  !! subroutine read_input opens the input file (if existent) and reads it (default filename 'cioverlap.input')
  !! global arrays are allocated here (or in subroutines called from here)
  USE sysparam
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
  PUBLIC ::  sort_block_these_dets, reorder_ci, debug_print_sorted

CONTAINS

  !--------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE sort_block_these_dets(SSDs,CIind,CIblkmap,blocks,blk_stind,sblk_stind,nblks,nsblks,nSD,nel,nMO)
    IMPLICIT NONE
    INTEGER(KIND=ishort), DIMENSION(:,:), INTENT(INOUT) :: SSDs
    INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: CIind
    INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: CIblkmap
    INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: blk_stind
    INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE, INTENT(OUT) ::sblk_stind
    INTEGER(KIND=ishort), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: blocks
    INTEGER(KIND=ilong), INTENT(IN) :: nSD
    INTEGER(KIND=ilong), INTENT(IN) :: nel
    INTEGER(KIND=ilong), INTENT(IN) :: nMO
    INTEGER(KIND=ilong), INTENT(OUT) :: nblks
    INTEGER(KIND=ilong), INTENT(OUT) ::nsblks

    INTEGER(KIND=ilong):: i
    INTEGER(KIND=ilong):: block
    INTEGER(KIND=ilong)::sblock
    INTEGER:: allocstat

    nblks=1
    nsblks=1

    allocstat = myalloc(Ciblkmap,nSD,'+ CI block map')
    IF (allocstat.NE.0) STOP 1

    allocstat=myalloc(CIind,nSD,'+ CI index')
    IF(allocstat.NE.0)STOP 1
    DO i=1,nSD
      CIind(i)=i
    END DO

    CALL master_sort_vec_by_excitation(SSDs,CIind,nSD,nel,iOne,nMO)

    ! count the blocks
    nblks = 1
    nsblks = 1
    DO i=2,nSD
      IF(.NOT.all(SSDs(i,:).EQ.SSDs(i-1,:)))THEN
        nblks = nblks + 1
        IF(.NOT.all(SSDs(i,1:nel-1).EQ.SSDs(i-1,1:nel-1)))THEN
          nsblks = nsblks + 1
        END IF
      END IF
    END DO

    ! find the blocks and populate the arrays
    allocstat = myalloc(blocks,nel+1_ilong,nblks,'+ blocks ')
    IF(allocstat.NE.0) STOP 1
    allocstat = myalloc(blk_stind, nblks+1_ilong, '+ block startind')
    IF(allocstat.NE.0) STOP 1
    allocstat = myalloc(sblk_stind,nsblks+1_ilong,'+ superblk startind')
    IF(allocstat.NE.0) STOP 1

    block = 1
    sblock = 1

    blk_stind(1) = 1
    sblk_stind(sblock) = 1
    DO i=2,nSD
      CIblkmap(CIind(i-1)) = block
      IF(.NOT.all(SSDs(i,:).EQ.SSDs(i-1,:)))THEN
        blocks(1,block) = i-1
        blocks(2:nel+1,block) = SSDs(i-1,:)
        block = block + 1
        blk_stind(block) = i
        IF(.NOT.all(SSDs(i,1:nel-1).EQ.SSDs(i-1,1:nel-1)))THEN
          sblock = sblock + 1
          sblk_stind(sblock) = block
        END IF
      END IF
    END DO
    CIblkmap(CIind(nSD)) = block
    blocks(1,block) = nSD
    blocks(2:nel+1, block) = SSDs(nSD,:)
    blk_stind(nblks+1) = nSD + 1
    sblk_stind(nsblks+1) = nblks + 1

  END SUBROUTINE sort_block_these_dets

  !--------------------------------------------------------------------------------------------------------------------------------

  RECURSIVE SUBROUTINE master_sort_vec_by_excitation(array,vector,nrow,ncol,col,norb)
    ! Sort the determinant strings in lexical order
    IMPLICIT NONE
    INTEGER(KIND=ishort), DIMENSION(:,:), INTENT(INOUT) :: array
    INTEGER(KIND=ishort), DIMENSION(:), INTENT(INOUT) :: vector
    INTEGER(KIND=ilong), INTENT(IN) :: nrow
    INTEGER(KIND=ilong), INTENT(IN) :: ncol
    INTEGER(KIND=ilong), INTENT(IN) :: col
    INTEGER(KIND=ilong), INTENT(IN) :: norb

    INTEGER(KIND=ilong):: irow
    INTEGER(KIND=ilong):: st_row

    INTEGER(KIND=ilong):: itmp

    IF(debug)THEN
      WRITE(6,'("master_sort_vec_by_excitation for col =", I5, " nrow =", I10)') col, nrow
      WRITE(6,'("First det.:", 100I4)') (array(1, itmp),itmp=1,ncol)
    END IF

    CALL sort_vec_by_excitation(array, vector, nrow, ncol, col, col, norb)

    IF(col.GE.ncol) RETURN ! Finished sorting all columns

    st_row = 1
    DO irow=2,nrow
      IF (array(irow, col).NE.array(irow-1, col)) THEN
        IF (irow-st_row.GE.2) THEN
          ! Sort the subarray that is constant up to col
          CALL master_sort_vec_by_excitation(array(st_row:irow-1,:), vector(st_row:irow-1), irow-st_row, ncol, col+1, norb)
        END IF
        st_row = irow
      END IF
    END DO
    IF (nrow-st_row.GE.2) THEN
      CALL master_sort_vec_by_excitation(array(st_row:nrow,:), vector(st_row:nrow), nrow-st_row+1, ncol, col+1, norb)
    END IF

  END SUBROUTINE master_sort_vec_by_excitation
  !--------------------------------------------------------------------------------------------------------------------------------

  RECURSIVE SUBROUTINE sort_vec_by_excitation(array,vector,nrow,ncol,col,low,high)
    ! Sort the determinant strings for one electron index (col)
    IMPLICIT NONE
    INTEGER(KIND=ishort), DIMENSION(:,:), INTENT(INOUT) :: array ! array may (will!) be bigger than (nrow,ncol)
    INTEGER(KIND=ishort), DIMENSION(:), INTENT(INOUT) :: vector
    INTEGER(KIND=ilong), INTENT(IN) :: nrow ! nrow and ncol are for counting how far to go in 'array'
    INTEGER(KIND=ilong), INTENT(IN) ::ncol ! nrow and ncol are for counting how far to go in 'array'
    INTEGER(KIND=ilong), INTENT(IN) ::col ! nrow and ncol are for counting how far to go in 'array'
    INTEGER(KIND=ilong), INTENT(IN) :: low ! nrow and ncol are for counting how far to go in 'array'
    INTEGER(KIND=ilong), INTENT(IN) :: high ! nrow and ncol are for counting how far to go in 'array'

    REAL(KIND=sip):: avg
    INTEGER(KIND=ilong):: i
    INTEGER(KIND=ilong):: allocstat
    INTEGER(KIND=ilong):: count_low
    INTEGER(KIND=ilong):: count_high
    INTEGER(KIND=ilong):: lowest_low
    INTEGER(KIND=ilong):: lowest_high
    INTEGER(KIND=ilong):: highest_low
    INTEGER(KIND=ilong):: highest_high
    INTEGER(KIND=ishort), DIMENSION(:,:), ALLOCATABLE :: arr_low
    INTEGER(KIND=ishort), DIMENSION(:,:), ALLOCATABLE :: arr_high
    INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE :: vec_low
    INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE :: vec_high

    IF (nrow.LE.1) RETURN

    allocstat = myalloc(arr_low,nrow,ncol)
    allocstat = allocstat + myalloc(arr_high,nrow,ncol) !,'+ arr_high')
    allocstat = allocstat + myalloc(vec_low,nrow) !,'+ vec_low')
    allocstat = allocstat + myalloc(vec_high,nrow) !,'+ vec_high')
    IF (allocstat.NE.0) THEN
      WRITE(6,*) 'Not enough memory for parallel sorting step!'
      STOP 1
    END IF

    avg=dble(high+low)/2.0_dop

    count_low=0
    count_high=0
    highest_low=low
    highest_high=int(avg+0.5)
    lowest_low=highest_high
    lowest_high=high

    DO i=1,nrow
      IF(array(i,col).LE.avg)THEN
        count_low = count_low + 1
        IF(array(i,col).LT.lowest_low)  lowest_low = array(i,col)
        IF(array(i,col).GT.highest_low) highest_low = array(i,col)

        arr_low(count_low,:) = array(i,:)
        vec_low(count_low) = vector(i)
      ELSE
        count_high = count_high + 1
        IF(array(i,col).LT.lowest_high)  lowest_high = array(i,col)
        IF(array(i,col).GT.highest_high) highest_high = array(i,col)

        arr_high(count_high,:) = array(i,:)
        vec_high(count_high) = vector(i)
      END IF
    END DO

    IF((count_low.GT.0).AND.(highest_low.GT.lowest_low))THEN
      CALL sort_vec_by_excitation(arr_low,vec_low,count_low,ncol,col,lowest_low,highest_low)
    END IF
    IF((count_high.GT.0).AND.(highest_high.GT.lowest_high))THEN
      CALL sort_vec_by_excitation(arr_high,vec_high,count_high,ncol,col,lowest_high,highest_high)
    END IF

    array(1:count_low,:) = arr_low(1:count_low,:)
    vector(1:count_low) = vec_low(1:count_low)
    array(count_low+1:count_low+count_high,:) = arr_high(1:count_high,:)
    vector(count_low+1:count_low+count_high) = vec_high(1:count_high)

    allocstat = mydealloc(arr_low)!,'- arr_low')
    allocstat = mydealloc(arr_high)!,'- arr_high')
    allocstat = mydealloc(vec_low)!,'- vec_low')
    allocstat = mydealloc(vec_high)!,'- vec_high')
  END SUBROUTINE sort_vec_by_excitation

  !--------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE debug_print_sorted(nstate, ndet, nelec, nblocks, sorted_dets, sorted_ci_index, blocks, cicoef, ci_block)
    IMPLICIT NONE

    INTEGER(KIND=ilong), INTENT(IN) :: nstate
    INTEGER(KIND=ilong), INTENT(IN) :: ndet
    INTEGER(KIND=ilong), INTENT(IN) :: nelec
    INTEGER(KIND=ilong), INTENT(IN) :: nblocks

    INTEGER(KIND=ishort), DIMENSION(ndet, nelec) :: sorted_dets
    INTEGER(KIND=ishort), DIMENSION(ndet)        :: sorted_ci_index
    INTEGER(KIND=ishort), DIMENSION(nelec+1_ilong, nblocks) :: blocks
    REAL(KIND=dop),   DIMENSION(nstate, ndet)   :: cicoef
    INTEGER(KIND=ishort), DIMENSION(ndet)        :: ci_block

    INTEGER(KIND=ilong):: i
    INTEGER(KIND=ilong)::j
    INTEGER(KIND=ilong):: block

    block = 1
    DO i=1,min(ndet,maxdebug)
      IF(i.GT.blocks(1,block))THEN
        WRITE(6,'("block ",I8," :")',ADVANCE='no')block
        WRITE(6,'("(",I9")")',ADVANCE='no')blocks(1,block)
        DO j=2,nelec+1
          WRITE(6,'(x,I3)',ADVANCE='no')blocks(j,block)
        END DO
        j=nelec+1
        WRITE(6,'(/)')
        block = block + 1
      END IF
      WRITE(6,'(I9,x,I9,x,I5,x,1PE15.8," :")',ADVANCE='no')i,sorted_ci_index(i),&
        &ci_block(sorted_ci_index(i)),cicoef(1, sorted_ci_index(i))
      DO j=1,nelec
        WRITE(6,'(x,I3)',ADVANCE='no')sorted_dets(i,j)
      END DO
      WRITE(6,*)''
    END DO
    WRITE(6,'("block ",I8," :")',ADVANCE='yes')block
    WRITE(6,'("(",I9")")',ADVANCE='no')blocks(1,block)
    DO j=2,nelec+1
      WRITE(6,'(x,I3)',ADVANCE='no')blocks(j,block)
    END DO
    WRITE(6,'(/)')
  END SUBROUTINE debug_print_sorted

  !--------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE reorder_ci(nstate, ndet, sorted_ci_index, cicoef, ciblock_alpha, ciblock_beta)
    IMPLICIT NONE
    INTEGER(KIND=ilong), INTENT(IN) :: nstate
    INTEGER(KIND=ilong), INTENT(IN) :: ndet
    INTEGER(KIND=ishort), DIMENSION(:), INTENT(IN) :: sorted_ci_index
    REAL(KIND=dop),   DIMENSION(:, :), INTENT(INOUT) :: cicoef
    INTEGER(KIND=ishort), DIMENSION(:), INTENT(INOUT) :: ciblock_alpha
    INTEGER(KIND=ishort), DIMENSION(:), INTENT(INOUT) :: ciblock_beta

    ! If ctmp is generated as a normal (rather than allocatable) array, then ifort breaks
    !    for large array sizes!
    REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE    :: ctmp
    INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE :: ciblock_alpha_tmp
    INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE :: ciblock_beta_tmp
    INTEGER(KIND=ilong):: idet
    INTEGER(KIND=ilong):: ci_ind
    INTEGER(KIND=ilong):: allocstat

    allocstat = myalloc(ciblock_alpha_tmp, ndet, '+ ciblock_alpha_tmp')
    IF (allocstat.NE.0) STOP 1
    allocstat = myalloc(ciblock_beta_tmp, ndet, '+ ciblock_beta_tmp')
    IF (allocstat.NE.0) STOP 1
    allocstat = myalloc(ctmp, nstate, ndet, '+ ctmp')
    IF (allocstat.NE.0) STOP 1

    ctmp(:,:)  = cicoef(:,:)
    ciblock_alpha_tmp(:) = ciblock_alpha(:)
    ciblock_beta_tmp(:) = ciblock_beta(:)

    DO idet=1,ndet
      ci_ind = sorted_ci_index(idet)
      
      cicoef(:, idet) = ctmp(:, ci_ind)
      ciblock_alpha(idet) = ciblock_alpha_tmp(ci_ind)
      ciblock_beta(idet)  = ciblock_beta_tmp(ci_ind)
    END DO
      
    allocstat = mydealloc(ctmp, '- ctmp')
    allocstat = mydealloc(ciblock_alpha_tmp, '- ciblock_alpha_tmp')
    allocstat = mydealloc(ciblock_beta_tmp, '- ciblock_beta_tmp')

  END SUBROUTINE reorder_ci  

  !--------------------------------------------------------------------------------------------------------------------------------
END MODULE sortblockmod
