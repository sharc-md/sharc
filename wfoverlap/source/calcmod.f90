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

MODULE calcmod
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
  PRIVATE :: laplace_vec_from_col, laplace_vec_from_lastcol
  PUBLIC :: precalc_bovl, precalc_ovl_direct, calc_cioverlap, &
    & ov_dys_mem, ov_dys_direct, determinant

CONTAINS

  !-------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE precalc_bovl(mo_ovl_mix, nelec, nsblocks_b, sblock_st_b, nblocks_a, blocks_a, blocks_b, block_ovl, maxb)
    ! Precompute the block-overlaps using an adaptive formalism that either computes determinants directly or
    !  constructs them from Laplace vectors
    IMPLICIT NONE

    REAL(KIND=dop), DIMENSION(:,:), INTENT(IN) :: mo_ovl_mix
    INTEGER(KIND=ilong), INTENT(IN) :: nelec
    INTEGER(KIND=ilong), INTENT(IN) :: nsblocks_b
    INTEGER(KIND=ishort), DIMENSION(:), INTENT(IN) :: sblock_st_b
    INTEGER(KIND=ilong), INTENT(IN) :: nblocks_a
    INTEGER(KIND=ishort), DIMENSION(:,:), INTENT(IN) :: blocks_a
    INTEGER(KIND=ishort), DIMENSION(:,:), INTENT(IN) :: blocks_b
    REAL(KIND=dop), DIMENSION(:,:), INTENT(INOUT) :: block_ovl
    INTEGER(KIND=ilong), INTENT(IN) :: maxb

    INTEGER(KIND=ilong):: ablock
    INTEGER(KIND=ilong):: bblock
    INTEGER(KIND=ilong):: i
    INTEGER(KIND=ilong):: bsblock
    INTEGER(KIND=ilong):: nblock
    INTEGER(KIND=ilong):: last_b
    REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE :: temp_ovl_mat
    REAL(KIND=dop), DIMENSION(:), ALLOCATABLE :: laplace_vec
    INTEGER:: allocstat
    INTEGER:: thrd_id

    !$omp parallel default(shared) &
    !$omp& private(bsblock, ablock, bblock, i, temp_ovl_mat, last_b, laplace_vec,thrd_id)

    thrd_id=omp_get_thread_num()+1
    allocstat=myalloc(temp_ovl_mat,nelec,nelec,'+ temp_ovl_mat',thrd_id)
    IF(allocstat.NE.0)STOP 1
    allocstat=myalloc(laplace_vec,nelec,'+ laplace_vec',thrd_id)
    IF(allocstat.NE.0)STOP 1

    !$omp do schedule(dynamic)
    DO bsblock=1,nsblocks_b
      nblock = sblock_st_b(bsblock+1)-sblock_st_b(bsblock)
      IF (nblock > nelec)THEN
        IF(debug) WRITE(6,'("Using Laplace formalism for superblock", I9, " with", I10, " blocks.")')bsblock, nblock

        DO ablock=1, nblocks_a
          ! Compute the Laplace vector, which can be used for a quick computation of the determinant
          bblock=sblock_st_b(bsblock)
          IF(bblock.GT.maxb)CYCLE ! use CYCLE rather that return to be safe for OMP
          temp_ovl_mat(:, :) = mo_ovl_mix(blocks_a(2:, ablock), blocks_b(2:, bblock))
          CALL laplace_vec_from_lastcol(nelec, temp_ovl_mat, laplace_vec)

          DO bblock=sblock_st_b(bsblock), sblock_st_b(bsblock+1)-1
            IF(bblock.GT.maxb)CYCLE
            ! Scalar product with the Laplace vector
            last_b = blocks_b(nelec+1, bblock)
            block_ovl(ablock,bblock) = sum(mo_ovl_mix(blocks_a(2:, ablock), last_b) * laplace_vec(:))
          END DO
        END DO
      ELSE
        IF(debug) WRITE(6,'("Using direct  formalism for superblock", I9, " with", I10, " blocks.")')bsblock, nblock

        DO bblock=sblock_st_b(bsblock), sblock_st_b(bsblock+1)-1
          IF(bblock.GT.maxb)CYCLE
          DO ablock=1, nblocks_a
            temp_ovl_mat(:, :) = mo_ovl_mix(blocks_a(2:, ablock), blocks_b(2:, bblock))

            block_ovl(ablock,bblock)=determinant(temp_ovl_mat,nelec)
          END DO
        END DO
      END IF
    END DO
      !$omp end do
      !$omp end parallel
  END SUBROUTINE precalc_bovl

  !-------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE precalc_ovl_direct(mo_ovl_mix, nelec, nblocks_a, nblocks_b, blocks_a, blocks_b, block_ovl,MO,maxb)
    ! Precompute the block-overlaps directly as determinants
    ! If MO is greater or equal to zero, Dyson mode is assumed
    IMPLICIT NONE

    REAL(KIND=dop), DIMENSION(:,:), INTENT(IN) :: mo_ovl_mix
    INTEGER(KIND=ilong), INTENT(IN) :: nelec
    INTEGER(KIND=ilong), INTENT(IN) :: nblocks_a
    INTEGER(KIND=ilong), INTENT(IN) :: nblocks_b
    INTEGER(KIND=ishort), DIMENSION(:,:), INTENT(IN) :: blocks_a
    INTEGER(KIND=ishort), DIMENSION(:,:), INTENT(IN) :: blocks_b
    REAL(KIND=dop), DIMENSION(:,:), INTENT(INOUT) :: block_ovl
    INTEGER(KIND=ilong), INTENT(IN) :: MO
    INTEGER(KIND=ilong), INTENT(IN) :: maxb

    INTEGER(KIND=ilong):: ablock
    INTEGER(KIND=ilong):: bblock
    INTEGER(KIND=ilong):: j
    INTEGER(KIND=ilong):: k
    REAL(KIND=dop), DIMENSION(nelec,nelec) :: temp_ovl_mat
    INTEGER(KIND=ilong):: sign_a

    !$omp parallel default(shared) &
    !$omp& private(ablock, bblock, j, k, temp_ovl_mat, sign_a)

    IF(MO.le.0)THEN ! Cioverlap case

      !$omp do schedule(dynamic)
      DO ablock=1,nblocks_a
        DO bblock=1,nblocks_b
          IF(bblock.GT.maxb)CYCLE
          temp_ovl_mat(:, :) = mo_ovl_mix(blocks_a(2:, ablock), blocks_b(2:, bblock))

          block_ovl(ablock,bblock)=determinant(temp_ovl_mat,nelec)
        END DO
      END DO
      !$omp end do

    ELSE ! Dyson case

      !$omp do schedule(dynamic)
      DO ablock=1,nblocks_a
        DO bblock=1,nblocks_b
          IF(any(blocks_b(2:,bblock).EQ.MO))THEN
            k=0
            DO j=1,nelec+1
              IF(blocks_b(j+1,bblock).NE.MO)THEN ! this is the actual annihilation
                k=k+1
                temp_ovl_mat(:, k) = mo_ovl_mix(blocks_a(2:, ablock), blocks_b(j+1, bblock))
              ELSE
                sign_a=(-1)**(nelec+j)
              END IF
            END DO
            block_ovl(ablock,bblock)=determinant(temp_ovl_mat,nelec)*sign_a
          ELSE
            block_ovl(ablock,bblock)=dZero
          END IF
        END DO
      END DO
      !$omp end do

    ENDIF
    !$omp end parallel

  END SUBROUTINE precalc_ovl_direct

  !-------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE calc_cioverlap(wf_ovl, mo_ovl_mix, L_CIcoefs, R_CIcoefs, LQ_SSDblocks, RQ_SSDblocks, LQ_CI_blockmap, LP_CI_blockmap, &
    &  RP_CI_blockmap, P_block_ovl, RQ_block_stind, RQ_superblk_stind, L_nstate, R_nstate, L_nSD, &
    & LQ_nspin, RQ_nspin, LQ_nblocks, RQ_nsuperblks, maxRP )
    ! Calculate the cioverlap in the case that only the "P" blockoverlaps were precomputed.
    ! Unfortunately the way that turned out to be the most efficient one is not very straight forward.
    ! In a first step it is decided for every superblock whether the determinants are computed using
    !   LU decomposition (direct) or Laplace vectors.
    ! Then some quantities are precomputed to limit the number of do loops as these turned out to be
    !   detrimental for the performance.

    IMPLICIT NONE
    REAL(KIND=dop), DIMENSION(:,:), INTENT(OUT) :: wf_ovl
    REAL(KIND=dop), DIMENSION(:,:), INTENT(IN) :: mo_ovl_mix
    REAL(KIND=dop), DIMENSION(:,:), INTENT(IN) :: L_CIcoefs
    REAL(KIND=dop), DIMENSION(:,:), INTENT(IN) :: R_CIcoefs
    INTEGER(KIND=ishort), DIMENSION(:,:), INTENT(IN) :: LQ_SSDblocks
    INTEGER(KIND=ishort), DIMENSION(:,:), INTENT(IN) :: RQ_SSDblocks
    INTEGER(KIND=ishort), DIMENSION(:), INTENT(IN) :: LQ_CI_blockmap
    INTEGER(KIND=ishort), DIMENSION(:), INTENT(IN) :: LP_CI_blockmap
    INTEGER(KIND=ishort), DIMENSION(:), INTENT(IN) ::  RP_CI_blockmap
    REAL(KIND=dop), DIMENSION(:,:), INTENT(IN) :: P_block_ovl
    INTEGER(KIND=ishort), DIMENSION(:), INTENT(IN) :: RQ_block_stind
    INTEGER(KIND=ishort), DIMENSION(:), INTENT(IN) :: RQ_superblk_stind
    INTEGER(KIND=ilong), INTENT(IN) :: L_nstate
    INTEGER(KIND=ilong), INTENT(IN) :: R_nstate
    INTEGER(KIND=ilong), INTENT(IN) :: L_nSD
    INTEGER(KIND=ilong), INTENT(IN) :: LQ_nspin
    INTEGER(KIND=ilong), INTENT(IN) :: RQ_nspin
    INTEGER(KIND=ilong), INTENT(IN) :: LQ_nblocks
    INTEGER(KIND=ilong), INTENT(IN) :: RQ_nsuperblks
    INTEGER(KIND=ilong), INTENT(IN) :: maxRP

    INTEGER (KIND=ilong):: i
    INTEGER (KIND=ilong):: allocstat
    INTEGER (KIND=ilong):: R_det
    INTEGER (KIND=ilong):: LQ_block
    INTEGER (KIND=ilong):: RQ_block
    INTEGER (KIND=ilong):: RP_block
    INTEGER (KIND=ilong):: R_superblk
    INTEGER (KIND=ilong):: nblock
    INTEGER (KIND=ilong):: last_b

    REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE :: wf_ovl_part ! partial contribution, computed by one thread
    REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE :: temp_ovl_mat
    REAL(KIND=dop), DIMENSION(:),   ALLOCATABLE :: Q_blockovl
    REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE :: laplace_vecs
    REAL(KIND=dop), DIMENSION(:), ALLOCATABLE :: det_prod
    REAL(KIND=dop), DIMENSION(L_nstate) :: det_prod_L
    INTEGER:: thrd_id

    IF(debug)THEN
      WRITE(6,'("sblock_st_b_alpha", 100I5)') (RQ_superblk_stind(R_superblk),R_superblk=1,RQ_nsuperblks)
    END IF

    !$omp parallel default(shared)&
    !$omp& private(allocstat, R_det, LQ_block,RQ_block,RP_block, R_superblk,nblock, last_b, det_prod, det_prod_L, wf_ovl_part, Q_blockovl, temp_ovl_mat, laplace_vecs, i, thrd_id)&
    !$omp& reduction(+:wf_ovl)

    thrd_id=omp_get_thread_num()+1

    allocstat = myalloc(wf_ovl_part,L_nstate,R_nstate,'+ wf_ovl_part',thrd_id)
    allocstat = allocstat + myalloc(temp_ovl_mat,LQ_nspin,RQ_nspin,'+ temp_ovl_mat',thrd_id)
    allocstat = allocstat + myalloc(Q_blockovl, LQ_nblocks, '+ LQ_blockovl',thrd_id)
    allocstat = allocstat + myalloc(laplace_vecs, LQ_nspin, LQ_nblocks, '+ laplace_vecs',thrd_id)
    allocstat = allocstat + myalloc(det_prod, L_nSD, '+ det_prod',thrd_id)
    IF(allocstat.NE.0)STOP 1

    ! loop over the b superblocks
    !   This can be a load balancing problem if the superblocks sizes are uneven ...

    !$omp do schedule(dynamic)
    DO R_superblk=1,RQ_nsuperblks
      wf_ovl_part=dZero

      nblock = RQ_superblk_stind(R_superblk+1)-RQ_superblk_stind(R_superblk)
      IF (nblock > RQ_nspin)THEN
        IF(debug) WRITE(6,'("Using Laplace formalism for superblock", I9, " with", I10, " blocks.")')R_superblk, nblock

        ! Precompute the Laplace vectors for all ablocks with respect to this superblock
        RQ_block=RQ_superblk_stind(R_superblk)
        DO LQ_block=1, LQ_nblocks
          temp_ovl_mat(:, :) = mo_ovl_mix(LQ_SSDblocks(2:, LQ_block), RQ_SSDblocks(2:, RQ_block))
          CALL laplace_vec_from_lastcol(LQ_nspin, temp_ovl_mat, laplace_vecs(:,LQ_block))
        END DO

        DO RQ_block=RQ_superblk_stind(R_superblk), RQ_superblk_stind(R_superblk+1)-1
          last_b = RQ_SSDblocks(RQ_nspin+1, RQ_block)
!          Q_blockovl = dZero

          ! Precompute a vector of alpha block overlaps
          DO LQ_block=1, LQ_nblocks
            Q_blockovl(LQ_block) = sum( mo_ovl_mix(LQ_SSDblocks(2:, LQ_block), last_b) * laplace_vecs(:, LQ_block) )
          END DO

          DO R_det=RQ_block_stind(RQ_block), RQ_block_stind(RQ_block+1)-1
            RP_block  = RP_CI_blockmap(R_det)

            ! Cycle in case there was not enough memory to pre-compute this block
            IF(RP_block.GT.maxRP)CYCLE

            det_prod(:) = Q_blockovl(LQ_CI_blockmap(:)) * P_block_ovl(LP_CI_blockmap(:), RP_block)
            CALL DGEMV('N', L_nstate, L_nSD, dOne, L_CIcoefs(:,:), L_nstate, det_prod, iOne, dZero, det_prod_L, iOne)
            wf_ovl_part(:, :) = wf_ovl_part(:, :) + spread(det_prod_L(:), 2, R_nstate) * spread(R_CIcoefs(:, R_det), 1, L_nstate)
          END DO
        END DO
      ELSE
        IF(debug) WRITE(6,'("Using direct  formalism for superblock", I9, " with", I10, " blocks.")')R_superblk, nblock

        DO RQ_block=RQ_superblk_stind(R_superblk), RQ_superblk_stind(R_superblk+1)-1 ! loop over the R blocks in the superblock
          ! Precompute a vector of Q block overlaps
          DO LQ_block=1, LQ_nblocks
            temp_ovl_mat(:, :) = mo_ovl_mix(LQ_SSDblocks(2:, LQ_block), RQ_SSDblocks(2:, RQ_block))
            Q_blockovl(LQ_block) = determinant(temp_ovl_mat,LQ_nspin)
          END DO

          DO R_det=RQ_block_stind(RQ_block), RQ_block_stind(RQ_block+1)-1
            RP_block  = RP_CI_blockmap(R_det)

            ! Cycle in case there was not enough memory to pre-compute this block
            IF(RP_block.GT.maxRP)CYCLE

            det_prod(:) = Q_blockovl(LQ_CI_blockmap(:)) * P_block_ovl(LP_CI_blockmap(:), RP_block)
            CALL DGEMV('N', L_nstate, L_nSD, dOne, L_CIcoefs(:,:), L_nstate, det_prod, iOne, dZero, det_prod_L, iOne)
            wf_ovl_part(:, :) = wf_ovl_part(:, :) + spread(det_prod_L(:), 2, R_nstate) * spread(R_CIcoefs(:, R_det), 1, L_nstate)
          END DO
        END DO ! b blocks
      END IF
      wf_ovl=wf_ovl+wf_ovl_part
    END DO ! b superblocks
    !$omp end do

    allocstat = mydealloc(wf_ovl_part, '- wf_ovl_part',thrd_id)
    allocstat = mydealloc(temp_ovl_mat, '- temp_ovl_mat',thrd_id)
    allocstat = mydealloc(Q_blockovl, '- Q_blockovl',thrd_id)
    allocstat = mydealloc(laplace_vecs, '- laplace_vecs',thrd_id)
    allocstat = mydealloc(det_prod, '- det_prod',thrd_id)
    !$omp end parallel
  END SUBROUTINE calc_cioverlap

  !-------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE ov_dys_mem(wf_ovl, L_CIcoefs, R_CIcoefs, LQ_CI_blockmap, LP_CI_blockmap, RQ_CI_blockmap, RP_CI_blockmap, &
    & P_block_ovl, Q_block_ovl, L_nstate, R_nstate, L_nSD, R_nSD, RQ_SSDblocks, MO)

    ! Compute the cioverlap or Dyson contribution using precomputed block overlaps
    ! The central contraction is realized as a DGEMV call
    !   Possibly, a larger part of the matrix could be constructed in memory
    !   for a DGEMM call...
    !
    ! MO specifies the eliminated MO. In the case of a normal CI-overlap run, use MO=-1
    IMPLICIT NONE
    REAL(KIND=dop), DIMENSION(:,:), INTENT(INOUT) :: wf_ovl
    REAL(KIND=dop), DIMENSION(:,:), INTENT(IN) :: L_CIcoefs
    REAL(KIND=dop), DIMENSION(:,:), INTENT(IN) :: R_CIcoefs
    INTEGER(KIND=ishort), DIMENSION(:), INTENT(IN) :: LQ_CI_blockmap
    INTEGER(KIND=ishort), DIMENSION(:), INTENT(IN) :: LP_CI_blockmap
    INTEGER(KIND=ishort), DIMENSION(:), INTENT(IN) :: RQ_CI_blockmap
    INTEGER(KIND=ishort), DIMENSION(:), INTENT(IN) :: RP_CI_blockmap
    REAL(KIND=dop), DIMENSION(:,:), INTENT(IN) :: P_block_ovl
    REAL(KIND=dop), DIMENSION(:,:), INTENT(IN) :: Q_block_ovl
    INTEGER(KIND=ilong), INTENT(IN) :: L_nstate
    INTEGER(KIND=ilong), INTENT(IN) :: R_nstate
    INTEGER(KIND=ilong), INTENT(IN) :: L_nSD
    INTEGER(KIND=ilong), INTENT(IN) :: R_nSD
    INTEGER(KIND=ishort), DIMENSION(:,:), INTENT(IN) :: RQ_SSDblocks
    INTEGER(KIND=ilong), INTENT(IN):: MO

    INTEGER (KIND=ilong):: R_det
    INTEGER (KIND=ilong):: RQ_block
    INTEGER (KIND=ilong):: RP_block
    INTEGER (KIND=ilong):: allocstat
    INTEGER (KIND=ilong):: thrd_id, nthrd
    INTEGER (KIND=ilong):: istart, iend

    REAL(KIND=dop), DIMENSION(:), ALLOCATABLE :: det_prod
    REAL(KIND=dop), DIMENSION(L_nstate) :: det_prod_L

    REAL(KIND=dop), DIMENSION(L_nstate, R_nstate) :: wf_ovl_part ! partial contribution, computed by one thread

    wf_ovl = dZero

    !$omp parallel default(shared)&
    !$omp& private(R_det, RQ_block, RP_block, det_prod, det_prod_L, thrd_id, istart, iend) &
    !$omp& reduction(+:wf_ovl)

! For parallelization divide the L CI-vector among the different threads
    nthrd = omp_get_num_threads()
    thrd_id=omp_get_thread_num()+1
    istart = ((thrd_id-1) * L_nSD) / nthrd + 1
    iend   = min((thrd_id * L_nSD) / nthrd, L_nSD)
    allocstat = myalloc(det_prod, iend-istart+1, '+ det_prod', thrd_id)
    IF (allocstat.NE.0) STOP 1


    IF (iend.GE.istart) THEN ! fix in case more threads then determinants are present
      IF(debug) WRITE(6,'("ov_dys_mem: thrd_id, istart, iend", 100I5)') thrd_id, istart, iend
      DO R_det=1, R_nSD
        IF((MO.le.0).or.any(RQ_SSDblocks(2:,RQ_CI_blockmap(R_det)).EQ.MO))THEN
          RQ_block  = RQ_CI_blockmap(R_det)
          RP_block  = RP_CI_blockmap(R_det)

          det_prod(:) = Q_block_ovl(LQ_CI_blockmap(istart:iend), RQ_block) * &
          & P_block_ovl(LP_CI_blockmap(istart:iend), RP_block)
          CALL DGEMV('N', L_nstate, iend-istart+1, dOne, L_CIcoefs(:,istart:iend), &
          & L_nstate, det_prod, iOne, dZero, det_prod_L, iOne)
          wf_ovl(:, :) = wf_ovl(:, :) + spread(det_prod_L(:), 2, R_nstate) * spread(R_CIcoefs(:, R_det), 1, L_nstate)
        END IF
      END DO
    ELSE
      IF(debug) WRITE(6,'("ov_dys_mem: thrd_id, istart, iend - skipped!", 100I5)') thrd_id, istart, iend
    END IF

    allocstat = mydealloc(det_prod, '- det_prod', thrd_id)

    !$omp end parallel
  END SUBROUTINE ov_dys_mem

  !-------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE ov_dys_direct(wf_ovl, mo_ovl_mix, L_CIcoefs, R_CIcoefs, LQ_SSDblocks, RQ_SSDblocks, LQ_CI_blockmap, &
    & LP_CI_blockmap, RP_CI_blockmap, P_block_ovl, RQ_block_stind, &
    & L_nstate, R_nstate, L_nSD, LQ_nspin, RQ_nspin, LQ_nblocks, RQ_nblocks, MO , maxRP)
    ! Compute the cioverlap using on-the-fly direct determinant computation

    IMPLICIT NONE
    REAL(KIND=dop), DIMENSION(:,:), INTENT(OUT) :: wf_ovl
    REAL(KIND=dop), DIMENSION(:,:), INTENT(IN) :: mo_ovl_mix
    REAL(KIND=dop), DIMENSION(:,:), INTENT(IN) :: L_CIcoefs
    REAL(KIND=dop), DIMENSION(:,:), INTENT(IN) :: R_CIcoefs
    INTEGER(KIND=ishort), DIMENSION(:,:), INTENT(IN) :: LQ_SSDblocks
    INTEGER(KIND=ishort), DIMENSION(:,:), INTENT(IN) :: RQ_SSDblocks
    INTEGER(KIND=ishort), DIMENSION(:), INTENT(IN) :: LQ_CI_blockmap
    INTEGER(KIND=ishort), DIMENSION(:), INTENT(IN) :: LP_CI_blockmap
    INTEGER(KIND=ishort), DIMENSION(:), INTENT(IN) :: RP_CI_blockmap
    REAL(KIND=dop), DIMENSION(:,:), INTENT(IN) :: P_block_ovl
    INTEGER(KIND=ishort), DIMENSION(:), INTENT(IN) :: RQ_block_stind
    INTEGER(KIND=ilong), INTENT(IN) :: L_nstate
    INTEGER(KIND=ilong), INTENT(IN) :: R_nstate
    INTEGER(KIND=ilong), INTENT(IN) :: L_nSD
    INTEGER(KIND=ilong), INTENT(IN) :: LQ_nspin
    INTEGER(KIND=ilong), INTENT(IN) :: RQ_nspin
    INTEGER(KIND=ilong), INTENT(IN) :: LQ_nblocks
    INTEGER(KIND=ilong), INTENT(IN) :: RQ_nblocks
    INTEGER(KIND=ilong), INTENT(IN) :: MO
    INTEGER(KIND=ilong), INTENT(IN) :: maxRP

    INTEGER (KIND=ilong):: R_det
    INTEGER (KIND=ilong):: LQ_block
    INTEGER (KIND=ilong):: RQ_block
    INTEGER (KIND=ilong):: RP_block
    INTEGER(KIND=ilong):: i
    INTEGER(KIND=ilong)::j
    INTEGER(KIND=ilong)::k
    INTEGER(KIND=ilong):: sign_a

    REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE :: wf_ovl_part ! partial contribution, computed by one thread
    REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE    :: temp_ovl_mat
    REAL(KIND=dop), DIMENSION(:),   ALLOCATABLE :: LQ_ovl
    REAL(KIND=dop), DIMENSION(:), ALLOCATABLE :: det_prod
    REAL(KIND=dop), DIMENSION(L_nstate) :: det_prod_L
    INTEGER:: allocstat
    INTEGER:: thrd_id

    !$omp parallel default(shared)&
    !$omp& private(allocstat,thrd_id,R_det,LQ_block,RQ_block,RP_block,wf_ovl_part, LQ_ovl, det_prod, det_prod_L, temp_ovl_mat, i, j, k, sign_a)&
    !$omp& reduction(+:wf_ovl)

    thrd_id=omp_get_thread_num()+1

    allocstat = myalloc(wf_ovl_part,L_nstate,R_nstate,'+ wf_ovl_part',thrd_id)
    allocstat = allocstat + myalloc(temp_ovl_mat,LQ_nspin,RQ_nspin,'+ temp_ovl_mat',thrd_id)
    allocstat = allocstat + myalloc(LQ_ovl, LQ_nblocks, '+ LQ_ovl',thrd_id)
    allocstat = allocstat + myalloc(det_prod, L_nSD, '+ det_prod',thrd_id)
    IF(allocstat.NE.0)STOP 1

    !$omp do schedule(dynamic)
    DO RQ_block=1, RQ_nblocks
      wf_ovl_part = dZero
      IF(MO.le.0)THEN ! CI-overlap case

        DO LQ_block=1, LQ_nblocks
          temp_ovl_mat(:, :) = mo_ovl_mix(LQ_SSDblocks(2:, LQ_block), RQ_SSDblocks(2:, RQ_block))
          LQ_ovl(LQ_block) = determinant(temp_ovl_mat,LQ_nspin)
        END DO

      ELSEIF(any(RQ_SSDblocks(2:,RQ_block).EQ.MO))THEN ! Dyson case

        DO LQ_block=1, LQ_nblocks
          k=0
          DO j=1,RQ_nspin
            IF(RQ_SSDblocks(j+1,RQ_block).NE.MO)THEN ! this is the actual annihilation
              k=k+1
              DO i=1,LQ_nspin
                temp_ovl_mat(i, k) = mo_ovl_mix(LQ_SSDblocks(i+1,LQ_block),RQ_SSDblocks(j+1,RQ_block))
              END DO
            ELSE
              sign_a=(-1)**(LQ_nspin+j)
            END IF
          END DO
          LQ_ovl(LQ_block) = determinant(temp_ovl_mat,LQ_nspin)*sign_a
        END DO

      ELSE ! Dyson case for inactive MO
        CYCLE
      ENDIF

        ! Precompute a vector of alpha block overlaps
        DO R_det=RQ_block_stind(RQ_block), RQ_block_stind(RQ_block+1)-1
          RP_block  = RP_CI_blockmap(R_det)

          ! Cycle in case there was not enough memory to pre-compute this block
          IF(RP_block.GT.maxRP)CYCLE

          det_prod(:) = LQ_ovl(LQ_CI_blockmap(:)) * P_block_ovl(LP_CI_blockmap(:), RP_block)
          CALL DGEMV('N', L_nstate, L_nSD, dOne, L_CIcoefs(:,:), L_nstate, det_prod, iOne, dZero, det_prod_L, iOne)
          wf_ovl_part(:, :) = wf_ovl_part(:, :) + spread(det_prod_L(:), 2, R_nstate) * spread(R_CIcoefs(:, R_det), 1, L_nstate)
        END DO
        wf_ovl=wf_ovl + wf_ovl_part
    END DO
    !$omp end do

    allocstat = mydealloc(wf_ovl_part, '- wf_ovl_part',thrd_id)
    allocstat = mydealloc(temp_ovl_mat, '- temp_ovl_mat',thrd_id)
    allocstat = mydealloc(LQ_ovl, '- LQ_ovl',thrd_id)
    allocstat = mydealloc(det_prod, '- det_prod',thrd_id)
    !$omp end parallel
  END SUBROUTINE ov_dys_direct

  !-------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE ov_nopre_direct(wf_ovl, mo_ovl_mix, L_CIcoefs, R_CIcoefs, LQ_SSDblocks, RQ_SSDblocks, LP_SSDblocks, RP_SSDblocks, &
    & LQ_CI_blockmap, LP_CI_blockmap, RQ_CI_blockmap, RP_CI_blockmap, RQ_block_stind, &
    & L_nstate, R_nstate, L_nSD, LQ_nspin, RQ_nspin, LP_nspin, RP_nspin, LQ_nblocks, LP_nblocks, RQ_nblocks,maxRP)
    ! Compute the cioverlap using on-the-fly direct determinant computation for the P and Q blocks
    ! This is slower than the other cases but requires significantly reduced memory

    IMPLICIT NONE
    REAL(KIND=dop), DIMENSION(:,:), INTENT(OUT) :: wf_ovl
    REAL(KIND=dop), DIMENSION(:,:), INTENT(IN) :: mo_ovl_mix
    REAL(KIND=dop), DIMENSION(:,:), INTENT(IN) :: L_CIcoefs
    REAL(KIND=dop), DIMENSION(:,:), INTENT(IN) :: R_CIcoefs
    INTEGER(KIND=ishort), DIMENSION(:,:), INTENT(IN) :: LQ_SSDblocks
    INTEGER(KIND=ishort), DIMENSION(:,:), INTENT(IN) :: RQ_SSDblocks
    INTEGER(KIND=ishort), DIMENSION(:,:), INTENT(IN) :: LP_SSDblocks
    INTEGER(KIND=ishort), DIMENSION(:,:), INTENT(IN) :: RP_SSDblocks
    INTEGER(KIND=ishort), DIMENSION(:), INTENT(IN) :: LQ_CI_blockmap
    INTEGER(KIND=ishort), DIMENSION(:), INTENT(IN) :: LP_CI_blockmap
    INTEGER(KIND=ishort), DIMENSION(:), INTENT(IN) :: RQ_CI_blockmap
    INTEGER(KIND=ishort), DIMENSION(:), INTENT(IN) :: RP_CI_blockmap
    INTEGER(KIND=ishort), DIMENSION(:), INTENT(IN) :: RQ_block_stind
    INTEGER(KIND=ilong), INTENT(IN) :: L_nstate
    INTEGER(KIND=ilong), INTENT(IN) :: R_nstate
    INTEGER(KIND=ilong), INTENT(IN) :: L_nSD
    INTEGER(KIND=ilong), INTENT(IN) :: LQ_nspin
    INTEGER(KIND=ilong), INTENT(IN) :: RQ_nspin
    INTEGER(KIND=ilong), INTENT(IN) :: LP_nspin
    INTEGER(KIND=ilong), INTENT(IN) :: RP_nspin
    INTEGER(KIND=ilong), INTENT(IN) :: LQ_nblocks
    INTEGER(KIND=ilong), INTENT(IN) :: LP_nblocks
    INTEGER(KIND=ilong), INTENT(IN) :: RQ_nblocks
    INTEGER(KIND=ilong), INTENT(IN) :: maxRP

    INTEGER (KIND=ilong):: R_det
    INTEGER (KIND=ilong):: LQ_block
    INTEGER (KIND=ilong):: LP_block
    INTEGER (KIND=ilong):: RQ_block
    INTEGER (KIND=ilong):: RP_block
    INTEGER (KIND=ilong):: i, prstep

    REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE :: wf_ovl_part ! partial contribution, computed by one thread
    REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE    :: Q_ovl_mat
    REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE    :: P_ovl_mat
    REAL(KIND=dop), DIMENSION(:),   ALLOCATABLE :: LQ_ovl
    REAL(KIND=dop), DIMENSION(:),   ALLOCATABLE :: LP_ovl
    REAL(KIND=dop), DIMENSION(:),   ALLOCATABLE :: det_prod
    REAL(KIND=dop), DIMENSION(L_nstate) :: det_prod_L
    INTEGER:: allocstat
    INTEGER:: thrd_id

    prstep = RQ_nblocks / 20

    !$omp parallel default(shared)&
    !$omp& private(allocstat,thrd_id,R_det,LQ_block,RQ_block,LP_block,RP_block,wf_ovl_part, LQ_ovl, LP_ovl, det_prod, det_prod_L, Q_ovl_mat, P_ovl_mat)&
    !$omp& reduction(+:wf_ovl)

    thrd_id=omp_get_thread_num()+1

    allocstat = myalloc(wf_ovl_part,L_nstate,R_nstate,'+ wf_ovl_part',thrd_id)
    allocstat = allocstat + myalloc(Q_ovl_mat,LQ_nspin,RQ_nspin,'+ temp_ovl_mat',thrd_id)
    allocstat = allocstat + myalloc(P_ovl_mat,LP_nspin,RP_nspin,'+ temp_ovl_mat',thrd_id)
    allocstat = allocstat + myalloc(LQ_ovl, LQ_nblocks, '+ LQ_ovl',thrd_id)
    allocstat = allocstat + myalloc(LP_ovl, LP_nblocks, '+ LP_ovl',thrd_id)
    allocstat = allocstat + myalloc(det_prod, L_nSD, '+ det_prod',thrd_id)
    IF(allocstat.NE.0)STOP 1

    ! Is there a problem about load balancing here?
    !  One could use a different parallelization scheme ...

    !$omp do schedule(dynamic)
    DO RQ_block=1, RQ_nblocks
        IF (mod(RQ_block,prstep).EQ.0) WRITE(6,'(I3,"%")')RQ_block/prstep*5
        wf_ovl_part = dZero

        DO LQ_block=1, LQ_nblocks
          Q_ovl_mat(:, :) = mo_ovl_mix(LQ_SSDblocks(2:, LQ_block), RQ_SSDblocks(2:, RQ_block))
          LQ_ovl(LQ_block) = determinant(Q_ovl_mat,LQ_nspin)
        END DO

        DO R_det=RQ_block_stind(RQ_block), RQ_block_stind(RQ_block+1)-1
          RP_block  = RP_CI_blockmap(R_det)

          ! Proceed only if this block has not been considered before
          IF(RP_block.LE.maxRP)CYCLE

          DO LP_block=1, LP_nblocks
            P_ovl_mat(:,:) = mo_ovl_mix(LP_SSDblocks(2:, LP_block), RP_SSDblocks(2:, RP_block))
            LP_ovl(LP_block) = determinant(P_ovl_mat,LP_nspin)
          END DO

          det_prod(:) = LQ_ovl(LQ_CI_blockmap(:)) * LP_ovl(LP_CI_blockmap(:))
          CALL DGEMV('N', L_nstate, L_nSD, dOne, L_CIcoefs(:,:), L_nstate, det_prod, iOne, dZero, det_prod_L, iOne)
          wf_ovl_part(:, :) = wf_ovl_part(:, :) + spread(det_prod_L(:), 2, R_nstate) * spread(R_CIcoefs(:, R_det), 1, L_nstate)
        END DO
        wf_ovl=wf_ovl + wf_ovl_part
    END DO
    !$omp end do

    allocstat = mydealloc(wf_ovl_part, '- wf_ovl_part',thrd_id)
    allocstat = mydealloc(Q_ovl_mat, '- temp_ovl_mat',thrd_id)
    allocstat = mydealloc(P_ovl_mat, '- temp_ovl_mat',thrd_id)
    allocstat = mydealloc(LQ_ovl, '- LQ_ovl',thrd_id)
    allocstat = mydealloc(LP_ovl, '- LP_ovl',thrd_id)
    allocstat = mydealloc(det_prod, '- det_prod',thrd_id)
    !$omp end parallel
  END SUBROUTINE ov_nopre_direct

  !-------------------------------------------------------------------------------------------------------------------------------

  !  Calculate a determinant of a square matrix
  !
  REAL(KIND=dop) FUNCTION determinant(mat,order)
    INTEGER(KIND=ilong), INTENT(IN):: order
    REAL(KIND=dop), INTENT(IN), DIMENSION(1:order,1:order) :: mat
    !
    REAL(KIND=dop), DIMENSION(1:order,1:order) :: tempm
    INTEGER(KIND=ilong):: info
    INTEGER(KIND=ilong), DIMENSION(1:order) :: ipvt
    !   REAL(KIND=dop), DIMENSION(1:order) :: work
    !   REAL(KIND=dop):: detx(2)
    INTEGER::i
    tempm = mat
    ! we now replace linpack routines by lapack/blas routines -> better performance
    !  CALL dgefa(tempm(1:order,1:order),order,order,ipvt(1:order),info)
    !  CALL dgedi(tempm(1:order,1:order),order,order,ipvt(1:order),detx,work(1:order),10_ilong)

    determinant=1.0d0
    IF (order.gt.0) THEN
      !  determinant = detx(1) * 10.0_dop**detx(2)
      CALL dgetrf(order,order,tempm(1:order,1:order),order,ipvt(1:order),info)

      DO i=1,order
      IF (ipvt(i).NE.i) THEN
        determinant=-determinant*(tempm(i,i))
      ELSE
        determinant=determinant*(tempm(i,i))
      ENDIF
      ENDDO
    ENDIF

  END FUNCTION determinant
  !-------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE laplace_vec_from_col(n,A,k,vec)
    IMPLICIT NONE
    INTEGER(KIND=ilong), INTENT(IN) :: n          ! n is number of rows, k is row to ignore
    INTEGER(KIND=ilong), INTENT(IN) ::k          ! n is number of rows, k is row to ignore
    REAL(KIND=dop), INTENT(IN) :: A(n,n)            ! input matrix
    REAL(KIND=dop), INTENT(OUT) :: vec(n)           ! output vector of cofactors

    REAL(KIND=dop):: D(n-1,n-1)                    ! temporary matrix
    INTEGER(KIND=ilong):: q                        ! loop index

    DO q=1,n
      D( 1:q-1 , 1:k-1 )=A( 1  :q-1 , 1:k-1)
      D( q:n-1 , 1:k-1 )=A( q+1:n   , 1:k-1)
      D( 1:q-1 , k:n-1 )=A( 1  :q-1 , k+1:n)
      D( q:n-1 , k:n-1 )=A( q+1:n   , k+1:n)
      vec(q)=determinant(D,n-1)*(-1)**(k+q)
    ENDDO

  END SUBROUTINE laplace_vec_from_col

  !-------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE laplace_vec_from_lastcol(n,A,vec)
    ! the input matrix is accepted as n x n matrix, but the last column is ignored, hence its content is arbitrary
    IMPLICIT NONE
    INTEGER(KIND=ilong), INTENT(IN) :: n            ! n is number of rows
    REAL(KIND=dop), INTENT(IN) :: A(n,n)            ! input matrix
    REAL(KIND=dop), INTENT(OUT) :: vec(n)           ! output vector of cofactors

    REAL(KIND=dop):: D(n-1,n-1)                    ! temporary matrix
    INTEGER(KIND=ilong):: q                        ! loop index

    DO q=1,n
      D( 1:q-1 , 1:n-1 )=A( 1  :q-1 , 1:n-1)
      D( q:n-1 , 1:n-1 )=A( q+1:n   , 1:n-1)
      vec(q)=determinant(D,n-1)*(-1)**(n+q)
    ENDDO

  END SUBROUTINE laplace_vec_from_lastcol

!-------------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE invert_matrix(A, Ainv)
    ! Compute S.C^T = C^-1
    IMPLICIT NONE
    REAL(KIND=dop), DIMENSION(:,:), INTENT(IN) :: A
    REAL(KIND=dop), DIMENSION(:,:), INTENT(OUT) :: Ainv

    INTEGER(KIND=ilong) :: N, N2, info
    REAL(KIND=dop), DIMENSION(size(A,1)) :: ipiv
    REAL(KIND=dop), DIMENSION(size(A,1)*size(A,2)) :: work

    N = size(A,1)
    N2 = size(A,2)

    IF(N.ne.N2)THEN
      WRITE(6,*)'Matrix not square, inversion not possible!'
      WRITE(0,*)'Matrix not square, inversion not possible!'
      STOP 1
    ENDIF

    Ainv = A
    CALL DGETRF(N, N, Ainv, N, ipiv, info)

    IF(info.ne.0)THEN
      WRITE(6,*)'LU factorization failed', info
      WRITE(0,*)'LU factorization failed', info
      STOP 1
    ENDIF

    CALL DGETRI(N, Ainv, N, ipiv, work, N, info)

    IF(info.ne.0)THEN
      WRITE(6,*)'Matrix inversion failed', info
      WRITE(0,*)'Matrix inversion failed', info
      STOP 1
    ENDIF
  END SUBROUTINE invert_matrix

!-------------------------------------------------------------------------------------------------------------------------------

END MODULE calcmod
