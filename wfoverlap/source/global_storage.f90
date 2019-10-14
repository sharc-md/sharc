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

MODULE global_storage
  USE sysparam
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: mix_AOovl, mix_AOdip, bra_MOcoefs, ket_MOcoefs, &
    & bra_CIcoefs, ket_CIcoefs, bra_SSD_alpha, bra_SSD_beta, &
    & ket_SSD_alpha, ket_SSD_beta, bra_SSDblocks_alpha, bra_SSDblocks_beta, ket_SSDblocks_alpha, ket_SSDblocks_beta, &
    & bra_block_stind_alpha,bra_block_stind_beta, ket_block_stind_alpha,ket_block_stind_beta, &
    & bra_nblocks_alpha, bra_nblocks_beta, ket_nblocks_alpha, ket_nblocks_beta, &
    & ket_superblk_stind_beta, ket_superblk_stind_alpha,bra_superblk_stind_beta,bra_superblk_stind_alpha, &
    & bra_nsuperblks_alpha, bra_nsuperblks_beta, ket_nsuperblks_alpha, ket_nsuperblks_beta, &
    & Q_block_ovl, P_block_ovl, &
    & bra_CI_blockmap_alpha, bra_CI_blockmap_beta, ket_CI_blockmap_alpha, ket_CI_blockmap_beta, &
    & bra_nAO, ket_nAO, bra_nMO, ket_nMO, bra_nstate, ket_nstate, bra_nSD, ket_nSD, bra_nel, ket_nel, &
    & bra_nalpha, ket_nalpha, bra_nbeta, ket_nbeta, ncore, ndocc, AOformat, force_direct_dets, force_noprecalc, same_aos, &
    & bra_moformat, ket_moformat, bra_ciformat, ket_ciformat, &
    & mix_AOfile, bra_MOfile, ket_MOfile, bra_CIfile, ket_CIfile, &
    & bra_sorted_ci_index_alpha, ket_sorted_ci_index_alpha, bra_sorted_ci_index_beta, ket_sorted_ci_index_beta, &
    & swap, MOprint, do_mixing_angles

  TYPE :: CIwf
    !! for future use to make things easier
    REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE :: MOcoefs
    REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE :: CIcoefs
    INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE :: SSD_alpha
    INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE :: SSD_beta
    INTEGER(KIND=ishort), DIMENSION(:,:), ALLOCATABLE :: SSDblocks_alpha
    INTEGER(KIND=ishort), DIMENSION(:,:), ALLOCATABLE :: SSDblocks_beta
    INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE :: block_stind_alpha
    INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE :: block_stind_beta
    INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE :: superblk_stind_alpha
    INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE :: superblk_stind_beta
    INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE :: CI_blockmap_alpha
    INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE :: CI_blockmap_beta
    INTEGER(KIND=ilong):: nAO
    INTEGER(KIND=ilong):: nMO
    INTEGER(KIND=ilong):: nstate
    INTEGER(KIND=ilong):: nSD
    INTEGER(KIND=ilong):: nel
    INTEGER(KIND=ilong):: nalpha
    INTEGER(KIND=ilong):: nbeta
  END TYPE CIwf

  REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE :: mix_AOovl
  REAL(KIND=dop), DIMENSION(:,:,:), ALLOCATABLE :: mix_AOdip
  REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE, TARGET :: bra_MOcoefs
  REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE, TARGET :: ket_MOcoefs
  ! CI-coefficients. Initially they are arranged according to the input file.
  !    Then they are reordered according to Q blocks for optimized computation.
  REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE, TARGET :: bra_CIcoefs
  REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE, TARGET :: ket_CIcoefs
  INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE, TARGET :: bra_CI_blockmap_alpha
  INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE, TARGET :: bra_CI_blockmap_beta
  INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE, TARGET :: ket_CI_blockmap_alpha
  INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE, TARGET :: ket_CI_blockmap_beta
  INTEGER(KIND=ishort), DIMENSION(:,:), ALLOCATABLE, TARGET :: bra_SSD_alpha
  INTEGER(KIND=ishort), DIMENSION(:,:), ALLOCATABLE, TARGET :: bra_SSD_beta
  INTEGER(KIND=ishort), DIMENSION(:,:), ALLOCATABLE, TARGET :: ket_SSD_alpha
  INTEGER(KIND=ishort), DIMENSION(:,:), ALLOCATABLE, TARGET :: ket_SSD_beta
  INTEGER(KIND=ishort), DIMENSION(:,:), ALLOCATABLE, TARGET :: bra_SSDblocks_alpha
  INTEGER(KIND=ishort), DIMENSION(:,:), ALLOCATABLE, TARGET :: bra_SSDblocks_beta
  INTEGER(KIND=ishort), DIMENSION(:,:), ALLOCATABLE, TARGET :: ket_SSDblocks_alpha
  INTEGER(KIND=ishort), DIMENSION(:,:), ALLOCATABLE, TARGET :: ket_SSDblocks_beta
  INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE, TARGET :: ket_superblk_stind_beta
  INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE, TARGET ::bra_superblk_stind_beta
  INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE, TARGET :: ket_superblk_stind_alpha
  INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE, TARGET ::bra_superblk_stind_alpha
  INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE, TARGET :: bra_block_stind_alpha
  INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE, TARGET ::bra_block_stind_beta
  INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE, TARGET :: ket_block_stind_alpha
  INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE, TARGET ::ket_block_stind_beta

  INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE, TARGET :: bra_sorted_ci_index_alpha
  INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE, TARGET :: ket_sorted_ci_index_alpha
  INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE, TARGET :: bra_sorted_ci_index_beta
  INTEGER(KIND=ishort), DIMENSION(:), ALLOCATABLE, TARGET :: ket_sorted_ci_index_beta
  INTEGER(KIND=ilong):: bra_nblocks_alpha
  INTEGER(KIND=ilong):: bra_nblocks_beta
  INTEGER(KIND=ilong):: ket_nblocks_alpha
  INTEGER(KIND=ilong):: ket_nblocks_beta
  INTEGER(KIND=ilong):: bra_nsuperblks_alpha
  INTEGER(KIND=ilong):: bra_nsuperblks_beta
  INTEGER(KIND=ilong):: ket_nsuperblks_alpha
  INTEGER(KIND=ilong):: ket_nsuperblks_beta
  REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE :: Q_block_ovl
  REAL(KIND=dop), DIMENSION(:,:), ALLOCATABLE :: P_block_ovl
  LOGICAL:: same_aos
  INTEGER(KIND=ilong):: bra_nAO
  INTEGER(KIND=ilong):: ket_nAO
  INTEGER(KIND=ilong):: bra_nMO
  INTEGER(KIND=ilong):: ket_nMO
  INTEGER(KIND=ilong):: bra_nstate
  INTEGER(KIND=ilong):: ket_nstate
  INTEGER(KIND=ilong):: bra_nSD
  INTEGER(KIND=ilong):: ket_nSD
  INTEGER(KIND=ilong):: bra_nel
  INTEGER(KIND=ilong):: bra_nalpha
  INTEGER(KIND=ilong):: bra_nbeta
  INTEGER(KIND=ilong):: ket_nel
  INTEGER(KIND=ilong):: ket_nalpha
  INTEGER(KIND=ilong):: ket_nbeta
  INTEGER(KIND=ilong):: ncore
  INTEGER(KIND=ilong):: ndocc
  INTEGER(KIND=ishort):: AOformat
  LOGICAL:: force_direct_dets
  LOGICAL:: force_noprecalc
  LOGICAL:: do_mixing_angles
  INTEGER(KIND=ishort):: bra_MOformat
  INTEGER(KIND=ishort)::ket_MOformat
  INTEGER(KIND=ishort):: bra_CIformat
  INTEGER(KIND=ishort):: ket_CIformat
  INTEGER(KIND=ishort):: MOprint
  !  REAL(KIND=dop):: block_ovl_thresh


  CHARACTER*200:: mix_AOfile!, ion_ovl_fn
  CHARACTER*200:: bra_MOfile
  CHARACTER*200:: ket_MOfile
  CHARACTER*200:: bra_CIfile
  CHARACTER*200:: ket_CIfile

  LOGICAL:: swap=.FALSE.

CONTAINS
!--------------------------------------------------------------------------------------------------------------------------------==

!--------------------------------------------------------------------------------------------------------------------------------==

END MODULE global_storage
