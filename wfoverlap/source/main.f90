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

PROGRAM main
  !> MyCioverlap version 1
  !! Program to compute overlaps from a COLUMBUS MRCI-calculations

  USE sysparam
  USE my_alloc
  USE memlog
  USE global_storage
  USE outputmod
  USE inputmod
  USE lowdin_mod
  USE iomod
  USE sortblockmod
  USE calcmod
#ifdef _OPENMP
  USE omp_lib
#else
#  define omp_get_num_threads() 1
#  define omp_get_max_threads() 1
#  define omp_get_thread_num() 0
#endif

  IMPLICIT NONE
  
  INTEGER (KIND=ilong):: i
  INTEGER (KIND=ilong):: j
  INTEGER (KIND=ilong):: MO
  INTEGER (KIND=ilong)::astate
  INTEGER (KIND=ilong):: bstate
  REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE :: tempm_dd
  REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE :: mo_ovl_mix
  REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE :: wf_ovl ! my result
  REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE :: wf_ovl_ortho ! after Lowdin orthogonalization
  REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE :: print_mat
  INTEGER(KIND=ilong):: starttime
  INTEGER(KIND=ilong):: lasttime
  INTEGER(KIND=ilong):: time
  INTEGER(KIND=ilong):: timeu
  INTEGER(KIND=ilong):: memleft
  INTEGER:: allocstat
  INTEGER:: mxthrd
  LOGICAL:: dyson_mode=.FALSE.
  LOGICAL:: precalc_Q

  CHARACTER:: transpAO='N'
  REAL(KIND=dop), DIMENSION(:,:), POINTER :: L_mocoefs
  REAL(KIND=dop), DIMENSION(:,:), POINTER :: R_mocoefs
  REAL(KIND=dop), DIMENSION(:,:), POINTER :: L_CIcoefs
  REAL(KIND=dop), DIMENSION(:,:), POINTER :: R_CIcoefs
  INTEGER(KIND=ishort), DIMENSION(:), POINTER :: LQ_CI_blockmap
  INTEGER(KIND=ishort), DIMENSION(:), POINTER :: RQ_CI_blockmap
  INTEGER(KIND=ishort), DIMENSION(:), POINTER :: LP_CI_blockmap
  INTEGER(KIND=ishort), DIMENSION(:), POINTER :: RP_CI_blockmap
  INTEGER(KIND=ishort), DIMENSION(:,:), POINTER :: LP_SSD
  INTEGER(KIND=ishort), DIMENSION(:,:), POINTER :: LQ_SSD
  INTEGER(KIND=ishort), DIMENSION(:,:), POINTER :: RP_SSD
  INTEGER(KIND=ishort), DIMENSION(:,:), POINTER :: RQ_SSD
  INTEGER(KIND=ishort), DIMENSION(:,:), POINTER :: LP_SSDblocks
  INTEGER(KIND=ishort), DIMENSION(:,:), POINTER :: LQ_SSDblocks
  INTEGER(KIND=ishort), DIMENSION(:,:), POINTER :: RP_SSDblocks
  INTEGER(KIND=ishort), DIMENSION(:,:), POINTER :: RQ_SSDblocks
  INTEGER(KIND=ishort), DIMENSION(:), POINTER :: RQ_superblk_stind
  INTEGER(KIND=ishort), DIMENSION(:), POINTER ::LQ_superblk_stind
  INTEGER(KIND=ishort), DIMENSION(:), POINTER :: RP_superblk_stind
  INTEGER(KIND=ishort), DIMENSION(:), POINTER ::LP_superblk_stind
  INTEGER(KIND=ishort), DIMENSION(:), POINTER :: LP_block_stind
  INTEGER(KIND=ishort), DIMENSION(:), POINTER ::LQ_block_stind
  INTEGER(KIND=ishort), DIMENSION(:), POINTER :: RP_block_stind
  INTEGER(KIND=ishort), DIMENSION(:), POINTER ::RQ_block_stind
  INTEGER(KIND=ishort), DIMENSION(:), POINTER :: LP_sorted_ci_index
  INTEGER(KIND=ishort), DIMENSION(:), POINTER :: RP_sorted_ci_index
  INTEGER(KIND=ishort), DIMENSION(:), POINTER :: LQ_sorted_ci_index
  INTEGER(KIND=ishort), DIMENSION(:), POINTER :: RQ_sorted_ci_index

  INTEGER(KIND=ilong):: L_nAO
  INTEGER(KIND=ilong):: R_nAO
  INTEGER(KIND=ilong):: L_nMO
  INTEGER(KIND=ilong):: R_nMO
  INTEGER(KIND=ilong):: L_nstate
  INTEGER(KIND=ilong):: R_nstate
  INTEGER(KIND=ilong):: L_nSD
  INTEGER(KIND=ilong):: R_nSD
  INTEGER(KIND=ilong):: L_nel
  INTEGER(KIND=ilong):: LQ_nspin
  INTEGER(KIND=ilong):: LP_nspin
  INTEGER(KIND=ilong):: R_nel
  INTEGER(KIND=ilong):: RQ_nspin
  INTEGER(KIND=ilong):: RP_nspin
  INTEGER(KIND=ilong):: LQ_nblocks
  INTEGER(KIND=ilong):: LP_nblocks
  INTEGER(KIND=ilong):: RQ_nblocks
  INTEGER(KIND=ilong):: RP_nblocks
  INTEGER(KIND=ilong):: maxRP
  INTEGER(KIND=ilong):: LQ_nsuperblks
  INTEGER(KIND=ilong):: LP_nsuperblks
  INTEGER(KIND=ilong):: RQ_nsuperblks
  INTEGER(KIND=ilong):: RP_nsuperblks

  LOGICAL:: beta_precalc
  LOGICAL:: ionization
  REAL(KIND=dop), DIMENSION(:,:,:), ALLOCATABLE :: dyson_orb_mobas
  !real(KIND=dop) :: dysnorm
  
  INTEGER(KIND=ilong):: ion_nel
  
  INTEGER(KIND=ilong):: ref_nel
  
#ifdef _OPENMP
  mxthrd = omp_get_max_threads()
  WRITE(6,'("Using OpenMP version. Maximum number of threads:", I6)') mxthrd
#else
  mxthrd = 1
  WRITE(6,*) "Using single threaded version"
#endif
  CALL system_clock(starttime,timeu)
  lasttime=starttime

  PRINT*,"-------------------------------------------------------------------------------"
  CALL read_input() ! inputmod.f90 - this populates global variables for simplicity

  !reading driver routines in iomod.f90

  ! read the CI determinants
  CALL read_CIvec(bra_CIfile,bra_CIformat,bra_CIcoefs,bra_SSD_alpha,bra_SSD_beta, &
    & bra_nMO,bra_nSD,bra_nstate,bra_nel,bra_nalpha,bra_nbeta,ncore)
  CALL read_CIvec(ket_CIfile,ket_CIformat,ket_CIcoefs,ket_SSD_alpha,ket_SSD_beta, &
    & ket_nMO,ket_nSD,ket_nstate,ket_nel,ket_nalpha,ket_nbeta,ncore)

  ! read the mixed AO overlap integrals
  CALL read_AO(mix_AOfile,AOformat,mix_AOovl,bra_nAO,ket_nAO,same_aos)

  ! read the MO coefficients
  CALL read_MO(bra_MOfile,bra_MOformat,bra_MOcoefs,bra_nMO,bra_nAO)
  CALL read_MO(ket_MOfile,ket_Moformat,ket_MOcoefs,ket_nMO,ket_nAO)

  CALL system_clock(time,timeu)
  WRITE(6, 300)(time-starttime)/(timeu*dOne), (time-lasttime)/(timeu*dOne)
  lasttime=time    
  PRINT*,"-------------------------------------------------------------------------------"

  WRITE(6,'("Number of <bra| electrons:    ",I12)') bra_nel
  WRITE(6,'("Number of <bra| alpha el:     ",I12)') bra_nalpha
  WRITE(6,'("Number of <bra| beta el:      ",I12)') bra_nbeta
  WRITE(6,'("Number of |ket> electrons:    ",I12)') ket_nel
  WRITE(6,'("Number of |ket> alpha el:     ",I12)') ket_nalpha
  WRITE(6,'("Number of |ket> beta el:      ",I12)') ket_nbeta
  WRITE(6,*)''
  WRITE(6,'("Number of <bra| MOs:          ",I12)') bra_nMO
  WRITE(6,'("Number of |ket> MOs:          ",I12)') ket_nMO
  WRITE(6,'("Number of discarded core MOs:        ",I5)') ncore
  WRITE(6,*)''
  WRITE(6,'("Number of <bra| states:       ",I12)') bra_nstate
  WRITE(6,'("Number of |ket> states:       ",I12)') ket_nstate
  WRITE(6,'("Number of <bra| determinants: ",I12)') bra_nSD
  WRITE(6,'("Number of |ket> determinants: ",I12)') ket_nSD
  WRITE(6,'("Number of det. pairs:   ",I18)')bra_nSD*ket_nSD
  CALL log_memory("after input")


  ! let's see if this is a normal cioverlap or a dyson calculation
  IF((bra_nalpha.NE.ket_nalpha).OR.(bra_nbeta.NE.ket_nbeta))THEN
    ! maximum ionization level is 1 and either number of alpha or number of beta has to be the same
    IF(((bra_nalpha.EQ.ket_nalpha).NEQV.(bra_nbeta.EQ.ket_nbeta)).AND.(abs(bra_nel-ket_nel).EQ.1))THEN
      dyson_mode=.TRUE.
      WRITE(6,*)"Joy! I calculate a dyson orbital!"
    ELSE
      WRITE(6,*)"The number of electrons in <bra| and |ket> differs by more than one!"
      WRITE(0,*)"The number of electrons in <bra| and |ket> differs by more than one!"
      STOP 1
    END IF
  ELSE
    WRITE(6,*)"Normal overlap calculation!"
  END IF

  PRINT*,"-------------------------------------------------------------------------------"
  
  !$omp parallel default(shared)

  !$omp sections

  PRINT*,"Sorting determinants and CI vectors for blocking"

  !======================================================================================
  ! A, alpha
  !$omp section
  CALL sort_block_these_dets(bra_SSD_alpha, bra_sorted_ci_index_alpha, bra_CI_blockmap_alpha, bra_SSDblocks_alpha, &
    & bra_block_stind_alpha,bra_superblk_stind_alpha,bra_nblocks_alpha,bra_nsuperblks_alpha,bra_nSD,bra_nalpha,bra_nMO)
  !======================================================================================
  ! A, beta
  !$omp section
  CALL sort_block_these_dets(bra_SSD_beta, bra_sorted_ci_index_beta, bra_CI_blockmap_beta, bra_SSDblocks_beta, &
    & bra_block_stind_beta,bra_superblk_stind_beta,bra_nblocks_beta,bra_nsuperblks_beta,bra_nSD,bra_nbeta,bra_nMO)
  !======================================================================================
  ! B, alpha
  !$omp section
  CALL sort_block_these_dets(ket_SSD_alpha, ket_sorted_ci_index_alpha, ket_CI_blockmap_alpha, ket_SSDblocks_alpha, &
    & ket_block_stind_alpha,ket_superblk_stind_alpha,ket_nblocks_alpha,ket_nsuperblks_alpha,ket_nSD,ket_nalpha,ket_nMO)
  !======================================================================================
  ! B, beta
  !$omp section
  CALL sort_block_these_dets(ket_SSD_beta, ket_sorted_ci_index_beta, ket_CI_blockmap_beta, ket_SSDblocks_beta, &
    & ket_block_stind_beta,ket_superblk_stind_beta,ket_nblocks_beta,ket_nsuperblks_beta,ket_nSD,ket_nbeta,ket_nMO)

  !$omp end sections
  !$omp end parallel
  CALL system_clock(time,timeu)
  WRITE(6, 300)(time-starttime)/(timeu*dOne), (time-lasttime)/(timeu*dOne)
  lasttime=time

  PRINT*,"-------------------------------------------------------------------------------"
  WRITE(6,'("Number of <bra| alpha blocks:       ",I12)') bra_nblocks_alpha
  WRITE(6,'("Number of |ket> alpha blocks:       ",I12)') ket_nblocks_alpha
  WRITE(6,'("Number of <bra|  beta blocks:       ",I12)') bra_nblocks_beta
  WRITE(6,'("Number of |ket>  beta blocks:       ",I12)') ket_nblocks_beta
  WRITE(6,'("Block-determinants to compute:",I18)')bra_nblocks_alpha*ket_nblocks_alpha+bra_nblocks_beta*ket_nblocks_beta
  WRITE(6,*)
  WRITE(6,'("Number of <bra| alpha superblocks:  ",I12)') bra_nsuperblks_alpha
  WRITE(6,'("Number of |ket> alpha superblocks:  ",I12)') ket_nsuperblks_alpha
  WRITE(6,'("Number of <bra|  beta superblocks:  ",I12)') bra_nsuperblks_beta
  WRITE(6,'("Number of |ket>  beta superblocks:  ",I12)') ket_nsuperblks_beta
  PRINT*,"-------------------------------------------------------------------------------"

  IF(debug)THEN
    WRITE(6,*)"sorted CI coefficients and determinants <bra| alpha "
    CALL debug_print_sorted(bra_nstate, bra_nSD, bra_nalpha, bra_nblocks_alpha, &
      &bra_SSD_alpha, bra_sorted_ci_index_alpha, bra_SSDblocks_alpha, bra_CIcoefs, bra_CI_blockmap_alpha)

    WRITE(6,*)"sorted CI coefficients and determinants <bra| beta "
    CALL debug_print_sorted(bra_nstate, bra_nSD, bra_nbeta, bra_nblocks_beta, &
      &bra_SSD_beta, bra_sorted_ci_index_beta, bra_SSDblocks_beta, bra_CIcoefs, bra_CI_blockmap_beta)

    WRITE(6,*)"sorted CI coefficients and determinants |ket> alpha"
    CALL debug_print_sorted(ket_nstate, ket_nSD, ket_nalpha, ket_nblocks_alpha, &
      &ket_SSD_alpha, ket_sorted_ci_index_alpha, ket_SSDblocks_alpha, ket_CIcoefs, ket_CI_blockmap_alpha)

    WRITE(6,*)"sorted CI coefficients and determinants |ket> beta"
    CALL debug_print_sorted(ket_nstate, ket_nSD, ket_nbeta, ket_nblocks_beta, &
      &ket_SSD_beta, ket_sorted_ci_index_beta, ket_SSDblocks_beta, ket_CIcoefs, ket_CI_blockmap_beta)
  END IF


  ! Here we decide which wavefunction is 'L'eft and which 'R'ight hand side and
  ! which spin is 'P'recomputed and which is 'Q' (you can make your own meme for 'Q' - it was just the next letter after 'P').
  !===================================================================================================================================
  IF(dyson_mode)THEN
    ! add things to be done for dyson here

    ! OK, first we look if we have ionization or electron attachment
    ! the L_cicoefs, and L_SSD are built on the fly for each annihilation
    ! This gives four possibilities: Always use the WF with more electrons as 'R' and the spin with different
    ! number of electorns as 'Q'.
    IF(bra_nel.LT.ket_nel)THEN
      !ionization
      ionization=.TRUE.
      PRINT*,"DEBUGMR: ionization!"

      ! we can already set some (spin-independet) pointers and variables
      L_CIcoefs => bra_CIcoefs
      L_mocoefs => bra_mocoefs
      L_nAO = bra_nAO
      L_nMO = bra_nMO
      L_nstate = bra_nstate
      L_nSD = bra_nSD
      L_nel = bra_nel

      R_CIcoefs => ket_CIcoefs
      R_mocoefs => ket_mocoefs
      R_nAO = ket_nAO
      R_nMO = ket_nMO
      R_nstate = ket_nstate
      R_nSD = ket_nSD
      R_nel = ket_nel

      transpAO='N'

      ion_nel = bra_nel
      ref_nel = ket_nel

      ! which spin is ionized - the one with equal number of electrons in bra and ket is always pre-computed

      IF(bra_nalpha.LT.ket_nalpha)THEN
        ! alpha is ionized, precompute beta
        beta_precalc=.TRUE.

        LP_SSD => bra_SSD_beta
        LP_SSDblocks => bra_SSDblocks_beta
        LP_block_stind => bra_block_stind_beta
        LP_superblk_stind => bra_superblk_stind_beta
        LP_CI_blockmap => bra_CI_blockmap_beta
        LP_sorted_ci_index => bra_sorted_ci_index_beta
        LP_nspin = bra_nbeta
        LP_nblocks = bra_nblocks_beta
        LP_nsuperblks = bra_nsuperblks_beta

        LQ_SSD => bra_SSD_alpha
        LQ_SSDblocks => bra_SSDblocks_alpha
        LQ_block_stind => bra_block_stind_alpha
        LQ_superblk_stind => bra_superblk_stind_alpha
        LQ_CI_blockmap => bra_CI_blockmap_alpha
        LQ_sorted_ci_index => bra_sorted_ci_index_alpha
        LQ_nspin = bra_nalpha
        LQ_nblocks = bra_nblocks_alpha
        LQ_nsuperblks = bra_nsuperblks_alpha

        RP_SSD => ket_SSD_beta
        RP_SSDblocks => ket_SSDblocks_beta
        RP_block_stind => ket_block_stind_beta
        RP_superblk_stind => ket_superblk_stind_beta
        RP_CI_blockmap => ket_CI_blockmap_beta
        RP_sorted_ci_index => ket_sorted_ci_index_beta
        RP_nspin = ket_nbeta
        RP_nblocks = ket_nblocks_beta
        RP_nsuperblks = ket_nsuperblks_beta

        RQ_SSD => ket_SSD_alpha
        RQ_SSDblocks => ket_SSDblocks_alpha
        RQ_block_stind => ket_block_stind_alpha
        RQ_superblk_stind => ket_superblk_stind_alpha
        RQ_CI_blockmap => ket_CI_blockmap_alpha
        RQ_sorted_ci_index => ket_sorted_ci_index_alpha
        RQ_nspin = ket_nalpha
        RQ_nblocks = ket_nblocks_alpha
        RQ_nsuperblks = ket_nsuperblks_alpha
      ELSE ! #<beta|.lt.#|beta>
        ! beta is ionized, precompute alpha
        beta_precalc=.FALSE.

        LP_SSD => bra_SSD_alpha
        LP_SSDblocks => bra_SSDblocks_alpha
        LP_block_stind => bra_block_stind_alpha
        LP_superblk_stind => bra_superblk_stind_alpha
        LP_CI_blockmap => bra_CI_blockmap_alpha
        LP_sorted_ci_index => bra_sorted_ci_index_alpha
        LP_nspin = bra_nalpha
        LP_nblocks = bra_nblocks_alpha
        LP_nsuperblks = bra_nsuperblks_alpha

        LQ_SSD => bra_SSD_beta
        LQ_SSDblocks => bra_SSDblocks_beta
        LQ_block_stind => bra_block_stind_beta
        LQ_superblk_stind => bra_superblk_stind_beta
        LQ_CI_blockmap => bra_CI_blockmap_beta
        LQ_sorted_ci_index => bra_sorted_ci_index_beta
        LQ_nspin = bra_nbeta
        LQ_nblocks = bra_nblocks_beta
        LQ_nsuperblks = bra_nsuperblks_beta

        RP_SSD => ket_SSD_alpha
        RP_SSDblocks => ket_SSDblocks_alpha
        RP_block_stind => ket_block_stind_alpha
        RP_superblk_stind => ket_superblk_stind_alpha
        RP_CI_blockmap => ket_CI_blockmap_alpha
        RP_sorted_ci_index => ket_sorted_ci_index_alpha
        RP_nspin = ket_nalpha
        RP_nblocks = ket_nblocks_alpha
        RP_nsuperblks = ket_nsuperblks_alpha

        RQ_SSD => ket_SSD_beta
        RQ_SSDblocks => ket_SSDblocks_beta
        RQ_block_stind => ket_block_stind_beta
        RQ_superblk_stind => ket_superblk_stind_beta
        RQ_CI_blockmap => ket_CI_blockmap_beta
        RQ_sorted_ci_index => ket_sorted_ci_index_beta
        RQ_nspin = ket_nbeta
        RQ_nblocks = ket_nblocks_beta
        RQ_nsuperblks = ket_nsuperblks_beta
      END IF ! which spin is ionized

    ELSE ! => ket_nel.lt.bra_nel

      !electron attachment - swap bra and ket for calculation
      PRINT*,"DEBUGMR: electron attachment!"

      ionization=.FALSE.

      R_CIcoefs => bra_CIcoefs
      R_mocoefs => bra_mocoefs
      R_nAO = bra_nAO
      R_nMO = bra_nMO
      R_nstate = bra_nstate
      R_nSD = bra_nSD
      R_nel = bra_nel

      L_CIcoefs => ket_CIcoefs
      L_mocoefs => ket_mocoefs
      L_nAO = ket_nAO
      L_nMO = ket_nMO
      L_nstate = ket_nstate
      L_nSD = ket_nSD
      L_nel = ket_nel

      transpAO='T'

      ion_nel = ket_nel
      ref_nel = bra_nel

      ! which spin is attached
      IF(ket_nalpha.LT.bra_nalpha)THEN
        beta_precalc=.TRUE.
        ! alpha is attached, precompute beta

        LP_SSD => ket_SSD_beta
        LP_SSDblocks => ket_SSDblocks_beta
        LP_block_stind => ket_block_stind_beta
        LP_superblk_stind => ket_superblk_stind_beta
        LP_CI_blockmap => ket_CI_blockmap_beta
        LP_sorted_ci_index => ket_sorted_ci_index_beta
        LP_nspin = ket_nbeta
        LP_nblocks = ket_nblocks_beta
        LP_nsuperblks = ket_nsuperblks_beta

        LQ_SSD => ket_SSD_alpha
        LQ_SSDblocks => ket_SSDblocks_alpha
        LQ_block_stind => ket_block_stind_alpha
        LQ_superblk_stind => ket_superblk_stind_alpha
        LQ_CI_blockmap => ket_CI_blockmap_alpha
        LQ_sorted_ci_index => ket_sorted_ci_index_alpha
        LQ_nspin = ket_nalpha
        LQ_nblocks = ket_nblocks_alpha
        LQ_nsuperblks = ket_nsuperblks_alpha

        RP_SSD => bra_SSD_beta
        RP_SSDblocks => bra_SSDblocks_beta
        RP_block_stind => bra_block_stind_beta
        RP_superblk_stind => bra_superblk_stind_beta
        RP_CI_blockmap => bra_CI_blockmap_beta
        RP_sorted_ci_index => bra_sorted_ci_index_beta
        RP_nspin = bra_nbeta
        RP_nblocks = bra_nblocks_beta
        RP_nsuperblks = bra_nsuperblks_beta

        RQ_SSD => bra_SSD_alpha
        RQ_SSDblocks => bra_SSDblocks_alpha
        RQ_block_stind => bra_block_stind_alpha
        RQ_superblk_stind => bra_superblk_stind_alpha
        RQ_CI_blockmap => bra_CI_blockmap_alpha
        RQ_sorted_ci_index => bra_sorted_ci_index_alpha
        RQ_nspin = bra_nalpha
        RQ_nblocks = bra_nblocks_alpha
        RQ_nsuperblks = bra_nsuperblks_alpha
      ELSE ! => #<beta|.gt.#|beta>
        !beta is attached, precompute alpha
        beta_precalc=.FALSE.

        LP_SSD => ket_SSD_alpha
        LP_SSDblocks => ket_SSDblocks_alpha
        LP_block_stind => ket_block_stind_alpha
        LP_superblk_stind => ket_superblk_stind_alpha
        LP_CI_blockmap => ket_CI_blockmap_alpha
        LP_sorted_ci_index => ket_sorted_ci_index_alpha
        LP_nspin = ket_nalpha
        LP_nblocks = ket_nblocks_alpha
        LP_nsuperblks = ket_nsuperblks_alpha
        LQ_SSD => ket_SSD_beta
        LQ_SSDblocks => ket_SSDblocks_beta
        LQ_block_stind => ket_block_stind_beta
        LQ_superblk_stind => ket_superblk_stind_beta
        LQ_CI_blockmap => ket_CI_blockmap_beta
        LQ_sorted_ci_index => ket_sorted_ci_index_beta
        LQ_nspin = ket_nbeta
        LQ_nblocks = ket_nblocks_beta
        LQ_nsuperblks = ket_nsuperblks_beta

        RP_SSD => bra_SSD_alpha
        RP_SSDblocks => bra_SSDblocks_alpha
        RP_block_stind => bra_block_stind_alpha
        RP_superblk_stind => bra_superblk_stind_alpha
        RP_CI_blockmap => bra_CI_blockmap_alpha
        RP_sorted_ci_index => bra_sorted_ci_index_alpha
        RP_nspin = bra_nalpha
        RP_nblocks = bra_nblocks_alpha
        RP_nsuperblks = bra_nsuperblks_alpha

        RQ_SSD => bra_SSD_beta
        RQ_SSDblocks => bra_SSDblocks_beta
        RQ_block_stind => bra_block_stind_beta
        RQ_superblk_stind => bra_superblk_stind_beta
        RQ_CI_blockmap => bra_CI_blockmap_beta
        RQ_sorted_ci_index => bra_sorted_ci_index_beta
        RQ_nspin = bra_nbeta
        RQ_nblocks = bra_nblocks_beta
        RQ_nsuperblks = bra_nsuperblks_beta
      END IF ! which spin is attached
    END IF ! ionization or electron attachment

    PRINT*,"RP_nspin:",RP_nspin
    PRINT*,"RQ_nspin:",RQ_nspin
    PRINT*,"LP_nspin:",LP_nspin
    PRINT*,"LQ_nspin:",LQ_nspin
  !===================================================================================================================================
  ELSE   ! => not dyson_mode => cioverlap

    L_CIcoefs => bra_CIcoefs
    L_mocoefs => bra_mocoefs
    L_nAO = bra_nAO
    L_nMO = bra_nMO
    L_nstate = bra_nstate
    L_nSD = bra_nSD
    L_nel = bra_nel
    R_CIcoefs => ket_CIcoefs
    R_mocoefs => ket_mocoefs
    R_nAO = ket_nAO
    R_nMO = ket_nMO
    R_nstate = ket_nstate
    R_nSD = ket_nSD
    R_nel = ket_nel

    transpAO='N'

    ! here is where the decision has to be made which spin is always pre-computed
    ! Always try the one with less blocks first for cioverlap.

    ! Do not do that yet - first eliminate all the bra_ and ket_ and _alpha and _beta entries and test.

    IF (bra_nblocks_alpha*ket_nblocks_alpha < bra_nblocks_beta*ket_nblocks_beta) THEN
      WRITE(6,*)'More beta than alpha blocks present - will pre-compute alpha'
      !      !WRITE(6,*)'       Use a different ms value to decrease memory requirements.'
      LQ_SSD => bra_SSD_beta
      LP_SSD => bra_SSD_alpha
      LQ_SSDblocks => bra_SSDblocks_beta
      LP_SSDblocks => bra_SSDblocks_alpha
      LQ_block_stind => bra_block_stind_beta
      LP_block_stind => bra_block_stind_alpha
      LQ_superblk_stind => bra_superblk_stind_beta
      LP_superblk_stind => bra_superblk_stind_alpha
      LQ_CI_blockmap => bra_CI_blockmap_beta
      LP_CI_blockmap => bra_CI_blockmap_alpha
      LQ_sorted_ci_index => bra_sorted_ci_index_beta
      LP_sorted_ci_index => bra_sorted_ci_index_alpha
      LQ_nspin = bra_nbeta
      LP_nspin = bra_nalpha
      LQ_nblocks = bra_nblocks_beta
      LP_nblocks = bra_nblocks_alpha
      LQ_nsuperblks = bra_nsuperblks_beta
      LP_nsuperblks = bra_nsuperblks_alpha

      RQ_SSD => ket_SSD_beta
      RP_SSD => ket_SSD_alpha
      RQ_SSDblocks => ket_SSDblocks_beta
      RP_SSDblocks => ket_SSDblocks_alpha
      RQ_block_stind => ket_block_stind_beta
      RP_block_stind => ket_block_stind_alpha
      RQ_superblk_stind => ket_superblk_stind_beta
      RP_superblk_stind => ket_superblk_stind_alpha
      RQ_CI_blockmap => ket_CI_blockmap_beta
      RP_CI_blockmap => ket_CI_blockmap_alpha
      RQ_sorted_ci_index => ket_sorted_ci_index_beta
      RP_sorted_ci_index => ket_sorted_ci_index_alpha
      RQ_nspin = ket_nbeta
      RP_nspin = ket_nalpha
      RQ_nblocks = ket_nblocks_beta
      RP_nblocks = ket_nblocks_alpha
      RQ_nsuperblks = ket_nsuperblks_beta
      RP_nsuperblks = ket_nsuperblks_alpha
    ELSE
      WRITE(6,*)'More alpha than beta blocks present - will pre-compute beta'
      LQ_SSD => bra_SSD_alpha
      LP_SSD => bra_SSD_beta
      LQ_SSDblocks => bra_SSDblocks_alpha
      LP_SSDblocks => bra_SSDblocks_beta
      LQ_block_stind => bra_block_stind_alpha
      LP_block_stind => bra_block_stind_beta
      LQ_superblk_stind => bra_superblk_stind_alpha
      LP_superblk_stind => bra_superblk_stind_beta
      LQ_CI_blockmap => bra_CI_blockmap_alpha
      LP_CI_blockmap => bra_CI_blockmap_beta
      LQ_sorted_ci_index => bra_sorted_ci_index_alpha
      LP_sorted_ci_index => bra_sorted_ci_index_beta
      LQ_nspin = bra_nalpha
      LP_nspin = bra_nbeta
      LQ_nblocks = bra_nblocks_alpha
      LP_nblocks = bra_nblocks_beta
      LQ_nsuperblks = bra_nsuperblks_alpha
      LP_nsuperblks = bra_nsuperblks_beta

      RQ_SSD => ket_SSD_alpha
      RP_SSD => ket_SSD_beta
      RQ_SSDblocks => ket_SSDblocks_alpha
      RP_SSDblocks => ket_SSDblocks_beta
      RQ_block_stind => ket_block_stind_alpha
      RP_block_stind => ket_block_stind_beta
      RQ_superblk_stind => ket_superblk_stind_alpha
      RP_superblk_stind => ket_superblk_stind_beta
      RQ_CI_blockmap => ket_CI_blockmap_alpha
      RP_CI_blockmap => ket_CI_blockmap_beta
      RQ_sorted_ci_index => ket_sorted_ci_index_alpha
      RP_sorted_ci_index => ket_sorted_ci_index_beta
      RQ_nspin = ket_nalpha
      RP_nspin = ket_nbeta
      RQ_nblocks = ket_nblocks_alpha
      RP_nblocks = ket_nblocks_beta
      RQ_nsuperblks = ket_nsuperblks_alpha
      RP_nsuperblks = ket_nsuperblks_beta

    ENDIF ! precomputed alpha or beta for cioverlap

  END IF ! dyson or cioverlap
  !===================================================================================================================================

  PRINT*,"-------------------------------------------------------------------------------"
  WRITE(6,*) ' Reordering CI-coefficients'

  ! Reorder the CI-coefficients according to the Q (not neccesarily precomputed) sorting order
  !    and adjust XP_ci ('P'recomputed) blockmap accordingly
  CALL reorder_ci(L_nstate, L_nSD, LQ_sorted_ci_index, L_CIcoefs, LQ_CI_blockmap, LP_CI_blockmap)
  CALL reorder_ci(R_nstate, R_nSD, RQ_sorted_ci_index, R_CIcoefs, RQ_CI_blockmap, RP_CI_blockmap)

  CALL system_clock(time,timeu)
  WRITE(6, 300)(time-starttime)/(timeu*dOne), (time-lasttime)/(timeu*dOne)
  lasttime=time
  PRINT*,"-------------------------------------------------------------------------------"
  WRITE(6,*) ' Transforming mixed overlap matrix'
  
  ! get the mixed overlap matrix in the MO basis.
  ! i.e.: construct one-particle operator matrix in cross-state molecular orbital basis.

  allocstat=myalloc(mo_ovl_mix,L_nMO,R_nMO,'+ mo_ovl_mix')
  IF(allocstat.NE.0) STOP

  ! get the mixed overlap matrix in the MO basis.
  ! i.e.: construct one-particle operator matrix in cross-state molecular orbital basis.
  allocstat=myalloc(tempm_dd,L_nAO,R_nMO)
  IF(allocstat.NE.0) STOP
  IF(AOformat.GE.0)THEN ! S-matrix read in from file
    ! multiply mixed AO-overlap with transpose of ket mocoef
    CALL DGEMM(transpAO,'T',L_nAO,R_nMO,R_nAO,dOne,mix_AOovl,L_nAO,R_MOcoefs,R_nMO,dZero,tempm_dd,L_nAO)
  ELSE ! reconstruct according to S.C_B^T = C_B^-1 (only for same_aos)
    WRITE(6,*) '   Using matrix inversion for AO overlap matrix'
    IF(R_nAO.ne.R_nMO)THEN
      WRITE(6,500) R_nAO, R_nMO
      WRITE(0,500) R_nAO, R_nMO
500 FORMAT("Inversion not possible! R_nAO = ", i0, " , R_nMO = ", i0)
      STOP 1
    ENDIF
    CALL invert_matrix(R_MOcoefs, tempm_dd)
  ENDIF
  ! multiply bra mocoef with result of above
  CALL DGEMM('N','N',L_nMO,R_nMO,L_nAO,dOne,L_MOcoefs,L_nMO,tempm_dd,L_nAO,dZero,mo_ovl_mix,L_nMO)
  allocstat=mydealloc(tempm_dd)
  IF(debug)THEN
    PRINT*,"-----------------------------------------------------------------------------"
    PRINT*,"left MO coefficients C_A"
    PRINT*,"-----------------------------------------------------------------------------"
    PRINT*, L_nAO, L_nMO
    DO i=1,L_nAO
      DO j=1,L_nMO
        WRITE(6,'(X,ES16.8e2)',ADVANCE='no')L_MOcoefs(i,j)
      END DO
      WRITE(6,*)
    END DO

    PRINT*,"-----------------------------------------------------------------------------"
    PRINT*,"right MO coefficients C_B"
    PRINT*,"-----------------------------------------------------------------------------"
    PRINT*, R_nAO, R_nMO
    DO i=1,R_nAO
      DO j=1,R_nMO
        WRITE(6,'(X,ES16.8e2)',ADVANCE='no')R_MOcoefs(i,j)
      END DO
      WRITE(6,*)
    END DO

    PRINT*,"-----------------------------------------------------------------------------"
    PRINT*,"mixed MO-overlap matrix (C_A(S_mix.C_B^T))"
    PRINT*,"-----------------------------------------------------------------------------"
    PRINT*, L_nMO, R_nMO
    DO i=1,L_nMO
      DO j=1,R_nMO-1
        WRITE(6,'(X,ES16.8e2)',ADVANCE='no')mo_ovl_mix(i,j)
      END DO
      j=R_nMO
      WRITE(6,'(X,ES16.8e2)',ADVANCE='yes')mo_ovl_mix(i,j)
    END DO
  END IF ! debug
  NULLIFY(L_MOcoefs)
  NULLIFY(R_MOcoefs)
  allocstat=mydealloc(bra_MOcoefs)
  allocstat=mydealloc(ket_MOcoefs)
  IF (debug.AND.(L_nMO.EQ.R_nMO)) THEN
    WRITE(6,'("Determinant of MO-overlap matrix:",ES24.14E3)') determinant(mo_ovl_mix, L_nMO)
  END IF
  CALL system_clock(time,timeu)
  WRITE(6, 300)(time-starttime)/(timeu*dOne), (time-lasttime)/(timeu*dOne)
  lasttime=time
  PRINT*,"-------------------------------------------------------------------------------"

  ! The P-determinants are precomputed here (or at least the part that can be stored in memory)
  
  ! Heuristic formula to get the memory needed later
  memleft = maxmem-mem_used - mxthrd*L_nSD*sizeof(1_dop) - mxthrd*LQ_nspin*LQ_nblocks*sizeof(1_dop)
  maxRP = 19*memleft / (20*sizeof(1_dop)*LP_nblocks)  
  
  IF(maxRP.LT.RP_nblocks)THEN
    WRITE(6,'(" Only ", I0, " out of ", I0, " columns of the P_ovl matrix are stored in memory (", I0," MB)!")') &
    & maxRP,RP_nblocks,memleft/1024/1024
    WRITE(6,*)'Increase the amount of memory to improve the efficiency.'
    IF(dyson_mode)THEN
      WRITE(6,*)'Low memory version not compatible with Dyson mode!'
      WRITE(0,*)'Low memory version not compatible with Dyson mode!'
      STOP 1
    END IF
  ELSE
    maxRP = RP_nblocks
  END IF
  allocstat = myalloc(P_block_ovl,LP_nblocks,maxRP,'+ P_block_ovl')
  
  IF((.NOT.force_direct_dets).and.(LP_nspin.gt.1))THEN
    WRITE(6,*),"Precomputing P determinants using adaptive formalism"
    CALL precalc_bovl(mo_ovl_mix, LP_nspin, RP_nsuperblks, RP_superblk_stind, &
      & LP_nblocks, LP_SSDblocks, RP_SSDblocks, P_block_ovl, maxRP)
  ELSE
    WRITE(6,*),"Precomputing P determinants using direct formalism"
    CALL precalc_ovl_direct(mo_ovl_mix, LP_nspin, LP_nblocks, &
      & RP_nblocks, LP_SSDblocks, RP_SSDblocks, P_block_ovl, -1, maxRP)
  END IF
  
    IF(debug)THEN
      WRITE(6,'(/,"(P)recomputed determinant block overlap matrix")')
      PRINT*,"-----------------------------------------------------------------------------"
      WRITE(6,'(I4,x,I4)')LP_nblocks,RP_nblocks
      DO i=1,LP_nblocks
        DO j=1,RP_nblocks-1
          WRITE(6,'(x,ES17.8E4)',ADVANCE='no')P_block_ovl(i,j)
        END DO
        WRITE(6,'(x,ES17.8E4)',ADVANCE='yes')P_block_ovl(i,LP_nblocks)
      END DO
    END IF  

  CALL system_clock(time,timeu)
  WRITE(6, 300)(time-starttime)/(timeu*dOne), (time-lasttime)/(timeu*dOne)
  lasttime=time

  ! now the computation of Q and actual overlap/dyson orb
  !===================================================================================================================================
  IF(dyson_mode)THEN
    ! allocations
    allocstat = myalloc(dyson_orb_mobas,R_nMO,L_nstate,R_nstate,'+ dyson_orb_mobas')
    IF(allocstat.NE.0)STOP

    dyson_orb_mobas = dZero

    ! Check if the Q determinants should also be precomputed
    IF(force_noprecalc)THEN
      precalc_Q=.FALSE.
    ELSE
      !    if (block_ovl_thresh.eq.dZero) then
      allocstat =  myalloc(Q_block_ovl,LQ_nblocks,RQ_nblocks,'+ Q_block_ovl')
      IF (allocstat.EQ.0) THEN
        precalc_Q=.TRUE.
      ELSE
        WRITE(6,*)'Allocation of Q block overlap matrix failed.'
        WRITE(6,*)' - Using on-the-fly algorithm for Q determinants.'
        precalc_Q=.FALSE.
      END IF
    END IF

    ! the big loop over the MOS
    DO MO=ndocc+1,R_nMO

      WRITE(6,*)'--------------------------------------------------------'
      WRITE(6,'("Annihilating MO ",I4," in RQ-part of wavefunction")')MO
!      WRITE(6,*)'--------------------------------------------------------'
!      WRITE(6,*)''

      ! development version - always force direct dets - adaptive formalism not implemented
      force_direct_dets=.TRUE.

      IF(precalc_Q)THEN
        IF(.NOT.force_direct_dets)THEN
          WRITE(6,*),"Precomputing Q determinants using adaptive formalism"
          WRITE(6,*)"Not implemented!"
          STOP
        ELSE
          WRITE(6,*),"Precomputing Q determinants using direct formalism"
          CALL precalc_ovl_direct(mo_ovl_mix, LQ_nspin, LQ_nblocks, &
            & RQ_nblocks, LQ_SSDblocks, RQ_SSDblocks, Q_block_ovl,MO,RQ_nblocks)
        END IF
!        PRINT*,"-------------------------------------------------------------------------------"

        CALL system_clock(time,timeu)
        WRITE(6, 300)(time-starttime)/(timeu*dOne), (time-lasttime)/(timeu*dOne)
        lasttime=time
      END IF ! precalc_Q

      dyson_orb_mobas(MO,:,:)=dZero
      IF (precalc_Q)THEN
        WRITE(6,*),"Dyson orbital part using precomputed quantities"
        CALL ov_dys_mem(dyson_orb_mobas(MO,:,:), L_CIcoefs, R_CIcoefs, LQ_CI_blockmap, LP_CI_blockmap, RQ_CI_blockmap, &
          & RP_CI_blockmap, P_block_ovl, Q_block_ovl, L_nstate, R_nstate, L_nSD, R_nSD,RQ_SSDblocks,MO)
      ELSE
        IF(.NOT.force_direct_dets)THEN
          WRITE(6,*),"Dyson orbital part using semi-direct approach and adaptive determinant formalism"
          WRITE(6,*)"Not implemented!"
          STOP
        ELSE
          WRITE(6,*),"Dyson orbital part using semi-direct approach and direct determinant formalism"
          CALL ov_dys_direct(dyson_orb_mobas(MO,:,:), mo_ovl_mix, L_CIcoefs, R_CIcoefs, &
            & LQ_SSDblocks, RQ_SSDblocks, LQ_CI_blockmap, &
            & LP_CI_blockmap, RP_CI_blockmap, P_block_ovl, RQ_block_stind, &
            & L_nstate, R_nstate, L_nSD, LQ_nspin, RQ_nspin, LQ_nblocks, RQ_nblocks, MO , maxRP)
        END IF
      END IF
      CALL system_clock(time,timeu)
      WRITE(6, 300)(time-starttime)/(timeu*dOne), (time-lasttime)/(timeu*dOne)
      lasttime=time
    END DO ! MO=ndocc+1,R_nMO
    PRINT*,"-------------------------------------------------------------------------------"

    allocstat = myalloc(print_mat, bra_nstate, ket_nstate, '+ print_mat')
    
    IF(ionization)THEN

      IF(beta_precalc)THEN
        PRINT*,"ALPHA ionization"
      ELSE
        PRINT*,"BETA ionization"
      END IF

      DO astate = 1,bra_nstate
        DO bstate = 1,ket_nstate
          print_mat(astate,bstate) = sum(dyson_orb_mobas(:,astate,bstate)**2)
        END DO
      END DO
    ELSE

      IF(beta_precalc)THEN
        PRINT*,"ALPHA attachment"
      ELSE
        PRINT*,"BETA attachment"
      END IF

      DO astate = 1,bra_nstate
        DO bstate = 1,ket_nstate
          print_mat(astate,bstate) = sum(dyson_orb_mobas(:,bstate,astate)**2)
        END DO
      END DO      
    END IF  

    WRITE(6,*)"Dyson norm matrix |<PsiA_i|PsiB_j>|^2"
    CALL print_state_mat(print_mat, debug)
      
    ! Print a renormalized result for truncated CI-vectors
    ! This is disabled for now. What is the correct scaling that should be used?
    !WRITE(6,*)
    !WRITE(6,*)"Renormalized Dyson norm matrix |<PsiA_i|PsiB_j>|^2"
    !CALL print_renormalized_mat(print_mat, L_CIcoefs, R_CIcoefs, L_nstate, R_nstate, .false., debug)
      
    IF (MOprint.ge.1) THEN
      IF (ionization) THEN  
        WRITE(6,*)" Dyson orbitals in reference |ket> MO basis:"
        CALL print_dyson_orbs(.false., R_nMO, dyson_orb_mobas)
      ELSE    
        WRITE(6,*)" Dyson orbitals in reference <bra| MO basis:"
        CALL print_dyson_orbs(.true., R_nMO, dyson_orb_mobas)
      END IF
    END IF


    !===================================================================================================================================
    ELSE ! => ci-overlap
    
      allocstat=myalloc(wf_ovl,L_nstate,R_nstate,'+ wf_ovl')
      IF(allocstat.NE.0) STOP

      ! Check if the Q determinants should be precomputed
      precalc_Q=.FALSE.

      IF(maxRP.GE.RP_nblocks)THEN
        allocstat =  myalloc(Q_block_ovl,LQ_nblocks,RQ_nblocks,'+ Q_block_ovl')
        IF (allocstat.EQ.0) THEN
          IF((.NOT.force_direct_dets).and.(LQ_nspin.gt.1))THEN
            WRITE(6,*),"Precomputing Q determinants using adaptive formalism"
            CALL precalc_bovl(mo_ovl_mix, LQ_nspin, RQ_nsuperblks, RQ_superblk_stind, &
              & LQ_nblocks, LQ_SSDblocks, RQ_SSDblocks, Q_block_ovl, RQ_nblocks)
          ELSE
            WRITE(6,*),"Precomputing Q determinants using direct formalism"
            CALL precalc_ovl_direct(mo_ovl_mix, LQ_nspin, LQ_nblocks, &
              & RQ_nblocks, LQ_SSDblocks, RQ_SSDblocks, Q_block_ovl, -1, RQ_nblocks)
          END IF
  
          CALL system_clock(time,timeu)
          WRITE(6, 300)(time-starttime)/(timeu*dOne), (time-lasttime)/(timeu*dOne)
          lasttime=time
          precalc_Q = .TRUE.
        ELSE
          WRITE(6,*)'Allocation of Q block overlap matrix failed.'
          WRITE(6,*)' - Using on-the-fly algorithm for Q determinants.'
        END IF

        PRINT*,"-------------------------------------------------------------------------------"
      END IF

      wf_ovl=dZero
      IF (precalc_Q)THEN
        WRITE(6,*),"Overlap using precomputed quantities"
        CALL ov_dys_mem(wf_ovl, L_CIcoefs, R_CIcoefs, LQ_CI_blockmap, LP_CI_blockmap, RQ_CI_blockmap, RP_CI_blockmap, &
          & P_block_ovl, Q_block_ovl, L_nstate, R_nstate, L_nSD, R_nSD, RQ_SSDblocks, -1)
      ELSE 
        IF((.NOT.force_direct_dets).and.(LQ_nspin.gt.1))THEN
          WRITE(6,*),"Overlap using semi-direct approach and adaptive determinant formalism"
          ! This subroutine does too much! Anything that needs that many arguments to work should be split up.
          CALL calc_cioverlap(wf_ovl, mo_ovl_mix, L_CIcoefs, R_CIcoefs, LQ_SSDblocks, RQ_SSDblocks, LQ_CI_blockmap, &
            & LP_CI_blockmap, RP_CI_blockmap, P_block_ovl, RQ_block_stind, RQ_superblk_stind, &
            & L_nstate, R_nstate, L_nSD, LQ_nspin, RQ_nspin, LQ_nblocks, RQ_nsuperblks, maxRP )
        ELSE
          WRITE(6,*),"Overlap using semi-direct approach and direct determinant formalism"
          CALL ov_dys_direct(wf_ovl, mo_ovl_mix, L_CIcoefs, R_CIcoefs, LQ_SSDblocks, RQ_SSDblocks, LQ_CI_blockmap, &
            & LP_CI_blockmap, RP_CI_blockmap, P_block_ovl, RQ_block_stind, &
            & L_nstate, R_nstate, L_nSD, LQ_nspin, RQ_nspin, LQ_nblocks, RQ_nblocks, -1, maxRP)
        END IF
        
      ! If there was not enough memory, compute the additional terms here
      !   Alternatively, this could be realized through a loop with several precalc calls
        IF(maxRP.LT.RP_nblocks)THEN
          CALL system_clock(time,timeu)
          WRITE(6, 300)(time-starttime)/(timeu*dOne), (time-lasttime)/(timeu*dOne)
          lasttime=time
          
          WRITE(6,*),"Computing the out-of-core terms using an on-the-fly approach"
          CALL ov_nopre_direct(wf_ovl, mo_ovl_mix, L_CIcoefs, R_CIcoefs, LQ_SSDblocks, RQ_SSDblocks, LP_SSDblocks, RP_SSDblocks, &
          & LQ_CI_blockmap, LP_CI_blockmap, RQ_CI_blockmap, RP_CI_blockmap, RQ_block_stind, &
          & L_nstate, R_nstate, L_nSD, LQ_nspin, RQ_nspin, LP_nspin, RP_nspin, LQ_nblocks, LP_nblocks, RQ_nblocks, maxRP)
        END IF
      END IF
        
      allocstat=mydealloc(bra_SSDblocks_alpha,'- blocks_a_alpha')
      allocstat=mydealloc(ket_SSDblocks_alpha,'- blocks_b_alpha')

      allocstat=mydealloc(mo_ovl_mix,'- mo_ovl_mix')
      
      allocstat=mydealloc(bra_block_stind_alpha,'- block_st_a_alpha')
      allocstat=mydealloc(ket_block_stind_alpha,'- block_st_b_alpha')
      allocstat=mydealloc(ket_superblk_stind_alpha,'- sblock_st_b_alpha')

      CALL system_clock(time,timeu)
      WRITE(6, 300)(time-starttime)/(timeu*dOne), (time-lasttime)/(timeu*dOne)
      lasttime=time
      PRINT*,"==============================================================================="
      WRITE(6,*)"Overlap matrix <PsiA_i|PsiB_j>"
      CALL print_state_mat(wf_ovl, debug)
    
      ! Print a renormalized result for truncated CI-vectors
      WRITE(6,*)
      WRITE(6,*)"Renormalized overlap matrix <PsiA_i|PsiB_j>"
      CALL print_renormalized_mat(wf_ovl, L_CIcoefs, R_CIcoefs, L_nstate, R_nstate, .true., debug)
    
      ! Do A Lowdin orthogonalization
      IF(bra_nstate.EQ.ket_nstate)THEN
        allocstat=myalloc(wf_ovl_ortho,bra_nstate,ket_nstate,'+ wf_ovl_ortho')
        IF(allocstat.NE.0) STOP 1

        WRITE(6,*)
        WRITE(6,*)'Performing Lowdin orthonormalization by SVD...'
        CALL lowdin(bra_nstate, wf_ovl, wf_ovl_ortho)
        
        WRITE(6,*)
        WRITE(6,*)"Orthonormalized overlap matrix <PsiA_i|PsiB_j>"
        CALL print_ortho_mat(wf_ovl, wf_ovl_ortho, debug)
        
        IF (do_mixing_angles) CALL mixing_angles(wf_ovl_ortho, debug)
      
        allocstat=mydealloc(wf_ovl_ortho,'- wf_ovl_ortho')
      END IF ! do Lowdin?
    !===================================================================================================================================
    
    END IF ! dyson mode or cioverlap

300 FORMAT("*** walltime: total - ",F12.3,"s, step - ",F12.3,"s ***")

  !-------------------------------------------------------------------------------------------------------------------------------

  END PROGRAM main
