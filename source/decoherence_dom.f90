!******************************************
!
!    SHARC Program Suite
!
!    Copyright (c) 2023 University of Vienna
!
!    This file is part of SHARC.
!
!    SHARC is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published
!    by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHARC is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    inside the SHARC manual.  If not, see
!    <http://www.gnu.org/licenses/>.
!
!******************************************

!> # Module  decoherence_dom
!> decay of mixing algorithm for SCP methods
!>
!> \author Yinan Shu
!> \date 04.30.2019
!>
!> This module contains subroutines that allow to perform a decoherence
!> correction according to the decay of mixing algorithm as described in
!> for NDM: M. D. Hack and D. G. Truhlar, J. Chem. Phys. 114, 9305 (2001)
!> for SCDM: C. Zhu, A. W. Jasper, and D. G. Truhlar, J. Chem. Phys. 120, 5543 (2004)
!> for CSDM: C. Zhu, S. Nangia, A. W. Jasper, and D. G. Truhlar, J. Chem. Phys. 121, 7658 (2004)

module decoherence_dom
  implicit none

  private computepvib
  private compute_svec_tau
  private compute_svec_tau_overlap
  private repropagate_coeff
  private repropagate_coeff_NPI
  private repropagate_coeff_laser
  private repropagate_coeff_NPI_laser
  private DoM_step_diag
  private DoM_step_MCH
  private DoM_coeff
  private decoherent_propagator

  public DoM_step
  public compute_svec_tau_control
  public initial_def
  public redecoherence
  public compute_decoforce

contains

! ===========================================================
!> This subroutine initialize the decoherence force
subroutine initial_def(traj, ctrl)
  use definitions
  use matrix
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

  integer :: istate, jstate, iatom, idir, ipol

  ! initialize decoherence force
  traj%decograd_ad(:,:)=0.d0

  ! compute pvib, s vector and decoherence time
  call compute_svec_tau_control(traj,ctrl)

  if (printlevel>4) then
    write(u_log,*) "the current pointer state in MCH basis is:", traj%state_MCH
    call vecwrite(ctrl%nstates, traj%decotime_MCH_s, u_log, 'decoherence time', 'F9.4')
    write(u_log,*) "the current pointer state in diag basis is:", traj%state_diag
    call vecwrite(ctrl%nstates, traj%decotime_diag_s, u_log, 'decoherence time', 'F9.4')
  endif

  ! Compute the decoherence force
  !call compute_decoforce(ctrl%natom, ctrl%nstates, traj%veloc_app_ad, traj%coeff_diag_s,&
  !  &traj%svec_diag_sad, traj%decotime_diag_s, traj%state_diag, traj%H_diag_ss, traj%decograd_ad)

  if (printlevel>4) then
    write(u_log,*) 'Decoherence force'
    do iatom=1,ctrl%natom
      write(u_log,'(3(F14.8,3X))') traj%decograd_ad(iatom,1), traj%decograd_ad(iatom,2), traj%decograd_ad(iatom,3)
    enddo
  endif

  ! Also it initializes the coherent coefficient
  traj%ccoeff_MCH_s=traj%coeff_MCH_s
  traj%ccoeff_diag_s=traj%coeff_diag_s

endsubroutine

! ===========================================================
! Control of decay of mixing

subroutine DoM_step(traj,ctrl)
  use definitions
  use matrix
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

  if (ctrl%pointer_basis==0) then 
    call DoM_step_diag(traj,ctrl)
  else if (ctrl%pointer_basis==1) then
    call DoM_step_MCH(traj,ctrl)
  endif

endsubroutine

! ===========================================================

subroutine DoM_step_diag(traj,ctrl)
  use definitions
  use matrix
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

  integer :: istate, iatom

  ! all s vector and tau are computed in diagonal basis 
  call compute_svec_tau_control(traj,ctrl)

  if (printlevel>4) then
    write(u_log,'(A,1X,I4,1X,A,1X,I4,1X,A)') 'Current pointer state:',traj%state_diag,'(diag)',traj%state_MCH,'(MCH)'
    write(u_log,'(A1,A20)') '|','decoherence time |'
    write(u_log,'(A1,A9,2X,A9,1X)') '|','diag |','MCH |'
    do istate=1,ctrl%nstates
      write(u_log,'(1X,F9.4,2X,F9.4,1X)') traj%decotime_diag_s(istate), traj%decotime_MCH_s(istate)
    enddo
    write(u_log,*) 'Coherently propagated old and new coefficients:'
    do istate=1,ctrl%nstates
      write(u_log,'(2(F14.9,1X),2X,2(F14.9,1X))') traj%coeff_diag_old_s(istate), traj%coeff_diag_s(istate)
    enddo
    write(u_log,*) 'Coherently propagated old and new densities:'
    do istate=1,ctrl%nstates
       write(u_log,'((F14.9,1X),2X,(F14.9,1X))') real(traj%coeff_diag_old_s(istate)*conjg(traj%coeff_diag_old_s(istate))), real(traj%coeff_diag_s(istate)*conjg(traj%coeff_diag_s(istate)))
    enddo
  endif

  ! perform decay of mixing of coefficients
  !call DoM_coeff(traj,ctrl)

  if (real(traj%coeff_diag_s(traj%state_diag)*conjg(traj%coeff_diag_s(traj%state_diag))) .gt. 1d-4) then
    call decoherent_propagator(&
      &ctrl%nstates, ctrl%dtstep, ctrl%nsubsteps, &
      &traj%decotime_diag_s, traj%decotime_diag_old_s, traj%state_diag, &
      &traj%coeff_diag_s, traj%Dtotal_ss)
  else 
    if (printlevel>2) then
      write(u_log,*) 'Pointer state density bellow 1d-4'
      write(u_log,*) 'Set decoherence propagator as unit'
    endif
    traj%decograd_ad=0.d0
    traj%RDtotal_ss=traj%Rtotal_ss
    traj%Dtotal_ss=dcmplx(0.d0,0.d0)
    do istate=1,ctrl%nstates
      traj%Dtotal_ss(istate,istate)=dcmplx(1.d0,0.d0)
    enddo
  endif

  ! Compute the decoherence force
  call compute_decoforce(ctrl%natom, ctrl%nstates, traj%veloc_ad, traj%coeff_diag_s,&
    &traj%svec_diag_sad, traj%decotime_diag_s, traj%state_diag, traj%H_diag_ss, traj%decograd_ad)

  if (printlevel>2) then
    write(u_log,*) 'Old and new coherent coefficients with decay of mixing in diagonal basis:'
    do istate=1,ctrl%nstates
      write(u_log,'(2(F14.9,1X),2X,2(F14.9,1X))') traj%coeff_diag_old_s(istate), traj%coeff_diag_s(istate)
    enddo
    write(u_log,*) 'Old and new densities with decay of mixing in diagonal basis:'
    do istate=1,ctrl%nstates
       write(u_log,'((F14.9,1X),2X,(F14.9,1X))') real(traj%coeff_diag_old_s(istate)*conjg(traj%coeff_diag_old_s(istate))), real(traj%coeff_diag_s(istate)*conjg(traj%coeff_diag_s(istate)))
    enddo
  endif

  if (printlevel>4) then
    write(u_log,*) 'Decoherence force'
    do iatom=1,ctrl%natom
      write(u_log,'(3(F14.8,3X))') traj%decograd_ad(iatom,1), traj%decograd_ad(iatom,2), traj%decograd_ad(iatom,3)
    enddo
  endif

  ! propagate coherent coefficients
  traj%ccoeff_diag_old_s=traj%ccoeff_diag_s
  call matvecmultiply(ctrl%nstates,traj%Rtotal_ss, traj%ccoeff_diag_old_s, traj%ccoeff_diag_s, 'n')
  if (printlevel>2) then
    write(u_log,*) 'CSDM Old and new diagonal coherent coefficients:'
    do istate=1,ctrl%nstates
      write(u_log,'(2(F14.9,1X),2X,2(F14.9,1X))') traj%ccoeff_diag_old_s(istate),traj%ccoeff_diag_s(istate)
    enddo
  endif


endsubroutine


! ===========================================================

subroutine DoM_step_MCH(traj,ctrl)
  use definitions
  use matrix
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

  integer :: istate, jstate, iatom, idir, ipol
  complex*16 :: coherent_coeff(ctrl%nstates)
  complex*16 :: Rprod(ctrl%nstates,ctrl%nstates)

!Steps of computing the decoherence term
!1. Compute pvib
!2. Compute s vector and decoherence time
!3. Repropagate the coefficients, now coeff at t+dt with decoherence part
!4. Compute the decoherence force at t+dt
   
  ! compute pvib, s vector and decoherence time
  call compute_svec_tau_control(traj,ctrl)

  if (printlevel>4) then
    write(u_log,'(A,1X,I4,1X,A,1X,I4,1X,A)') 'Current pointer state:',traj%state_diag,'(diag)',traj%state_MCH,'(MCH)'
    write(u_log,'(A1,A20)') '|','decoherence time |'
    write(u_log,'(A1,A9,2X,A9,1X)') '|','diag |','MCH |'
    do istate=1,ctrl%nstates
      write(u_log,'(1X,F9.4,2X,F9.4,1X)') traj%decotime_diag_s(istate), traj%decotime_MCH_s(istate)
    enddo
    write(u_log,*) 'Coherently propagated old and new coefficients:'
    do istate=1,ctrl%nstates 
      write(u_log,'(2(F14.9,1X),2X,2(F14.9,1X))') traj%coeff_diag_old_s(istate), traj%coeff_diag_s(istate)
    enddo
    write(u_log,*) 'Coherently propagated old and new densities:'
    do istate=1,ctrl%nstates
       coherent_coeff(istate)=traj%coeff_diag_s(istate)
       write(u_log,'((F14.9,1X),2X,(F14.9,1X))') real(traj%coeff_diag_old_s(istate)*conjg(traj%coeff_diag_old_s(istate))), real(traj%coeff_diag_s(istate)*conjg(traj%coeff_diag_s(istate)))
    enddo
  endif

  ! Repropagate the coefficients, now coeff at t+dt with decoherence part
  if (real(traj%coeff_MCH_old_s(traj%state_MCH)*conjg(traj%coeff_MCH_old_s(traj%state_MCH))) .gt. 1d-3) then
  !it is eazy to cause numerical instability when the pointer state has very low density
  
  if (ctrl%laser==0) then !laser fields

  ! call the appropriate propagator routine
  select case (ctrl%eeom)
    case (0)    ! constant interpolation, default for coupling=ddt
      call repropagate_coeff(&
        &ctrl%nstates, ctrl%dtstep, ctrl%nsubsteps, 1,&
        &traj%decotime_MCH_s, traj%decotime_MCH_old_s, traj%state_MCH,&
        &traj%coeff_diag_s, traj%coeff_diag_old_s,&
        &traj%U_ss, traj%U_old_ss,&
        &traj%H_MCH_ss, traj%H_MCH_old_ss,&
        &traj%NACdt_ss, traj%NACdt_old_ss,&
        &traj%RDtotal_ss, traj%Dtotal_ss)
    case (1)    ! linear interpolation, default for coupling=ddr,nacdr
      call repropagate_coeff(&
        &ctrl%nstates, ctrl%dtstep, ctrl%nsubsteps, 0,&
        &traj%decotime_MCH_s, traj%decotime_MCH_old_s, traj%state_MCH,&
        &traj%coeff_diag_s, traj%coeff_diag_old_s,&
        &traj%U_ss, traj%U_old_ss,&
        &traj%H_MCH_ss, traj%H_MCH_old_ss,&
        &traj%NACdt_ss, traj%NACdt_old_ss,&
        &traj%RDtotal_ss, traj%Dtotal_ss)
    case (3)    ! norm perserving interpolation 
      call repropagate_coeff_NPI(&
        &ctrl%nstates, ctrl%dtstep, ctrl%nsubsteps, ctrl%coupling,&
        &traj%decotime_MCH_s, traj%decotime_MCH_old_s, traj%state_MCH,&
        &traj%coeff_diag_s, traj%coeff_diag_old_s,&
        &traj%U_ss, traj%U_old_ss,&
        &traj%H_MCH_ss, traj%H_MCH_old_ss,&
        &traj%overlaps_ss,&
        &traj%RDtotal_ss, traj%Dtotal_ss)
  endselect

  else !laser fields
  select case (ctrl%eeom )
    case (0)    ! constant interpolation, default for coupling=ddt
      call repropagate_coeff_laser(&
        &ctrl%nstates, ctrl%dtstep, ctrl%nsubsteps, 1,&
        &traj%decotime_MCH_s, traj%decotime_MCH_old_s, traj%state_MCH,&
        &traj%coeff_diag_s, traj%coeff_diag_old_s,&
        &traj%U_ss, traj%U_old_ss,&
        &traj%H_MCH_ss, traj%H_MCH_old_ss,&
        &traj%NACdt_ss, traj%NACdt_old_ss,&
        &traj%DM_ssd,traj%DM_old_ssd,&
        &ctrl%laserfield_td( (traj%step-1)*ctrl%nsubsteps+2:traj%step*ctrl%nsubsteps+1 ,:),&
        &traj%RDtotal_ss, traj%Dtotal_ss)
    case (1)    ! linear interpolation, default for coupling=ddr,nacdr
      call repropagate_coeff_laser(&
        &ctrl%nstates, ctrl%dtstep, ctrl%nsubsteps, 0,&
        &traj%decotime_MCH_s, traj%decotime_MCH_old_s, traj%state_MCH,&
        &traj%coeff_diag_s, traj%coeff_diag_old_s,&
        &traj%U_ss, traj%U_old_ss,&
        &traj%H_MCH_ss, traj%H_MCH_old_ss,&
        &traj%NACdt_ss, traj%NACdt_old_ss,&
        &traj%DM_ssd,traj%DM_old_ssd,&
        &ctrl%laserfield_td( (traj%step-1)*ctrl%nsubsteps+2:traj%step*ctrl%nsubsteps+1 ,:),&
        &traj%RDtotal_ss, traj%Dtotal_ss)
    case (3)    ! norm perserving interpolation
      call repropagate_coeff_NPI_laser(&
        &ctrl%nstates, ctrl%dtstep, ctrl%nsubsteps, ctrl%coupling,&
        &traj%decotime_MCH_s, traj%decotime_MCH_old_s, traj%state_MCH,&
        &traj%coeff_diag_s, traj%coeff_diag_old_s,&
        &traj%U_ss, traj%U_old_ss,&
        &traj%H_MCH_ss, traj%H_MCH_old_ss,&
        &traj%overlaps_ss,&
        &traj%DM_ssd,traj%DM_old_ssd,&
        &ctrl%laserfield_td( (traj%step-1)*ctrl%nsubsteps+2:traj%step*ctrl%nsubsteps+1 ,:),&
        &traj%RDtotal_ss, traj%Dtotal_ss)
  endselect
  endif

  ! Compute the decoherence force
  call compute_decoforce(ctrl%natom, ctrl%nstates, traj%veloc_app_ad, traj%coeff_diag_s,&
    &traj%svec_diag_sad, traj%decotime_diag_s, traj%state_diag, traj%H_diag_ss, traj%decograd_ad)

  else
    if (printlevel>2) then
      write(u_log,*) 'Pointer state density bellow 1d-3'
      write(u_log,*) 'Set decoherence propagator as unit'
    endif
    traj%decograd_ad=0.d0
    traj%RDtotal_ss=traj%Rtotal_ss
    traj%Dtotal_ss=dcmplx(0.d0,0.d0)
    do istate=1,ctrl%nstates
      traj%Dtotal_ss(istate,istate)=dcmplx(1.d0,0.d0)
    enddo
  endif

  if (printlevel>2) then
    write(u_log,*) 'Old and new coherent coefficients with decay of mixing in MCH basis:'
    do istate=1,ctrl%nstates
      write(u_log,'(2(F14.9,1X),2X,2(F14.9,1X))') traj%coeff_diag_old_s(istate), traj%coeff_diag_s(istate)
    enddo
    write(u_log,*) 'Old and new densities with decay of mixing in MCH basis:'
    do istate=1,ctrl%nstates
       write(u_log,'((F14.9,1X),2X,(F14.9,1X))') real(traj%coeff_diag_old_s(istate)*conjg(traj%coeff_diag_old_s(istate))), real(traj%coeff_diag_s(istate)*conjg(traj%coeff_diag_s(istate)))
    enddo
  endif

  if (printlevel>4) then
    call matwrite(ctrl%nstates,traj%RDtotal_ss,u_log,'Total propagator','F14.9')
    call matwrite(ctrl%nstates,traj%Rtotal_ss,u_log,'Coherent propagator','F14.9')
    call matwrite(ctrl%nstates,traj%Dtotal_ss,u_log,'Decoherence propagator','F14.9')
    call matmultiply(ctrl%nstates, traj%Dtotal_ss, traj%Rtotal_ss, Rprod, 'nn')
    call matwrite(ctrl%nstates, Rprod,u_log,'C+D propagator','F14.9')
  endif

  if (printlevel>4) then
    write(u_log,*) 'Decoherence force'
    do iatom=1,ctrl%natom
      write(u_log,'(3(F14.8,3X))') traj%decograd_ad(iatom,1), traj%decograd_ad(iatom,2), traj%decograd_ad(iatom,3)
    enddo
  endif

  ! propagate coherent coefficients
  traj%ccoeff_diag_old_s=traj%ccoeff_diag_s
  call matvecmultiply(ctrl%nstates,traj%Rtotal_ss, traj%ccoeff_diag_old_s, traj%ccoeff_diag_s, 'n')
  if (printlevel>2) then
    write(u_log,*) 'CSDM Old and new diagonal coherent coefficients:'
    do istate=1,ctrl%nstates
      write(u_log,'(2(F14.9,1X),2X,2(F14.9,1X))') traj%ccoeff_diag_old_s(istate),traj%ccoeff_diag_s(istate)
    enddo
  endif
  

endsubroutine


! ===========================================================

subroutine redecoherence(traj,ctrl)
  use definitions
  use matrix
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

  integer :: istate, jstate, iatom, idir, ipol

  if (printlevel>2) then
    write(u_log,*) 'Pointer state changed, re-compute decoherence time and force'
    write(u_log,*) 'decoherence time before and after original computation in diagonal basis'
    do istate=1,ctrl%nstates
      write(u_log,'(2(F7.4,1X),4X,2(F7.4,1X))') traj%decotime_diag_old_s(istate),traj%decotime_diag_s(istate)
    enddo
  endif

  if (real(traj%coeff_diag_s(traj%state_diag)*conjg(traj%coeff_diag_s(traj%state_diag))) .gt. 1d-3) then

  ! Back propagate the decoherence time and coefficients
  do istate=1,ctrl%nstates
    traj%decotime_diag_s(istate)=traj%decotime_diag_old_s(istate)
  enddo
  do istate=1,ctrl%nstates
    traj%decotime_MCH_s(istate)=traj%decotime_MCH_old_s(istate)
  enddo

  ! recompute s vector based on the new pointer state
  call compute_svec_tau_control(traj,ctrl)

  ! Compute the decoherence force
  call compute_decoforce(ctrl%natom, ctrl%nstates, traj%veloc_app_ad, traj%coeff_diag_s,&
    &traj%svec_diag_sad, traj%decotime_diag_s, traj%state_diag, traj%H_diag_ss, traj%decograd_ad)

  if (printlevel>4) then
    write(u_log,*) 'Re-computation finished'
    write(u_log,'(A,1X,I4,1X,A,1X,I4,1X,A)') 'Current pointer state:',traj%state_diag,'(diag)',traj%state_MCH,'(MCH)'
    write(u_log,'(A1,A43)') '|','decoherence time |'
    write(u_log,'(A1,A20,2X,A20,1X)') '|','diag |','MCH |'
    write(u_log,'(A1,A9,2X,A9,2X,A9,2X,A9,1X)') '|','before |','after |','before |','after |'
    do istate=1,ctrl%nstates
      write(u_log,'(1X,3(F9.4,2X),F9.4,1X)') traj%decotime_diag_old_s(istate), traj%decotime_diag_s(istate), traj%decotime_MCH_old_s(istate), traj%decotime_MCH_s(istate)
    enddo
  endif

  else ! corresponds to if density .gt. 1d-3
    if (printlevel>2) then
      write(u_log,*) 'Pointer state density bellow 1d-3'
      write(u_log,*) 'Set decoherence force as zero, with unit decoherence propagator'
    endif
    traj%decograd_ad=0.d0

  endif


  if (printlevel>4) then
    write(u_log,*) 'Recomputed decoherence force'
    do iatom=1,ctrl%natom
      write(u_log,'(3(F14.8,3X))') traj%decograd_ad(iatom,1), traj%decograd_ad(iatom,2), traj%decograd_ad(iatom,3)
    enddo
  endif

endsubroutine

! ==========================================================
!> this subroutine computes the decoherence coefficients gradient
!> Used by Bulirsch-Stoer integrator
subroutine DoM_stepBS(traj,ctrl)
  use definitions
  use matrix
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

  integer :: istate, jstate, iatom, idir, ipol
  real*8 :: gRdcoeff_MCH_s(ctrl%nstates), gIdcoeff_MCH_s(ctrl%nstates)
  complex*16 :: den(ctrl%nstates,ctrl%nstates)

  ! This subroutine compute the following addtional terms:
  ! 1. The derivative of electronic coefficients from decoherence part in MCH basis
  ! 2. The derivative of nuclear coordinates, decoherence force in DIAG basis

  do istate=1,ctrl%nstates
    do jstate=1,ctrl%nstates
      den(istate,jstate)=traj%coeff_MCH_s(istate)*conjg(traj%coeff_MCH_s(jstate))
    enddo
  enddo

  ! compute pvib, s vector and decoherence time
  call compute_svec_tau_control(traj,ctrl)

  ! compute decoherence coefficient gradient
  gRdcoeff_MCH_s=0.d0
  gIdcoeff_MCH_s=0.d0

  do istate=1,ctrl%nstates
    if (istate.ne.traj%state_MCH) then 
      gRdcoeff_MCH_s(traj%state_MCH)=gRdcoeff_MCH_s(traj%state_MCH)+real(den(istate,istate))*traj%decotime_MCH_s(istate)
      gIdcoeff_MCH_s(traj%state_MCH)=gIdcoeff_MCH_s(traj%state_MCH)+real(den(istate,istate))*traj%decotime_MCH_s(istate)
      gRdcoeff_MCH_s(istate)=-0.5d0*traj%decotime_MCH_s(istate)*real(traj%coeff_MCH_s(istate))
      gIdcoeff_MCH_s(istate)=-0.5d0*traj%decotime_MCH_s(istate)*aimag(traj%coeff_MCH_s(istate))
    endif
  enddo
  if (abs(real(den(traj%state_MCH,traj%state_MCH))) .gt. 1.0d-3) then
    gRdcoeff_MCH_s(traj%state_MCH)=gRdcoeff_MCH_s(traj%state_MCH)*0.5d0*real(traj%coeff_MCH_s(traj%state_MCH))/real(den(traj%state_MCH,traj%state_MCH))
    gIdcoeff_MCH_s(traj%state_MCH)=gIdcoeff_MCH_s(traj%state_MCH)*0.5d0*aimag(traj%coeff_MCH_s(traj%state_MCH))/real(den(traj%state_MCH,traj%state_MCH))
  else  
    write(u_log,*) 'Pointer state density bellow 1d-3'
    write(u_log,*) 'Set decoherence force, decay of mixing as zero, with unit decoherence propagator'
    traj%decograd_ad=0.d0
    gRdcoeff_MCH_s=0.d0
    gIdcoeff_MCH_s=0.d0
  endif

  if (printlevel>4) then
    call vecwrite(ctrl%nstates, gRdcoeff_MCH_s, u_log, 'Real decoherence Re(dcoeff_mch_s) gradient','F14.9')
    call vecwrite(ctrl%nstates, gIdcoeff_MCH_s, u_log, 'Imaginary decoherence Im(dcoeff_mch_s) gradient','F14.9')
  endif

  ! add decoherence part
  traj%gRcoeff_MCH_s=traj%gRcoeff_MCH_s+gRdcoeff_MCH_s
  traj%gIcoeff_MCH_s=traj%gIcoeff_MCH_s+gIdcoeff_MCH_s

  if (printlevel>4) then
    call vecwrite(ctrl%nstates, traj%gRcoeff_MCH_s, u_log, 'total Re(coeff_mch_s) gradient','F14.9')
    call vecwrite(ctrl%nstates, traj%gIcoeff_MCH_s, u_log, 'total Im(coeff_mch_s) gradient','F14.9')
  endif

  ! Compute the decoherence force
  call compute_decoforce(ctrl%natom, ctrl%nstates, traj%veloc_app_ad, traj%coeff_diag_s,&
    &traj%svec_diag_sad, traj%decotime_diag_s, traj%state_diag, traj%H_diag_ss, traj%decograd_ad)

  if (printlevel>4) then
    write(u_log,*) 'Decoherence force'
    do iatom=1,ctrl%natom
      write(u_log,'(3(F14.8,3X))') traj%decograd_ad(iatom,1), traj%decograd_ad(iatom,2), traj%decograd_ad(iatom,3)
    enddo
  endif

endsubroutine

! ==========================================================
!> this subroutine controls how to compute the s vector
!> it returns both s vector and decoherence time in diagonal and MCH basis
!> compute_svec_control(traj,ctrl,pvib_ad)
subroutine compute_svec_tau_control(traj,ctrl)
  use definitions
  use matrix
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

  integer :: istate, iatom, idir
  real*8 :: p(ctrl%natom,3)
  real*8 :: pvib_ad(ctrl%natom,3)
  complex*16 :: cpNACdR_MCH_ssad(ctrl%nstates,ctrl%nstates,ctrl%natom,3)
  complex*16 :: cNACdR_MCH_ssad(ctrl%nstates,ctrl%nstates,ctrl%natom,3)
  complex*16 :: H_ss(ctrl%nstates,ctrl%nstates)
  complex*16 :: grad_diag_sad(ctrl%nstates,ctrl%natom,3)
  complex*16 :: cgrad_MCH_sad(ctrl%nstates,ctrl%natom,3)

  do iatom=1,ctrl%natom
    do idir=1,3
      p(iatom,idir)=traj%veloc_app_ad(iatom,idir)*traj%mass_a(iatom)
    enddo
  enddo

  grad_diag_sad=dcmplx(0.d0,0.d0)
  cgrad_MCH_sad=dcmplx(0.d0,0.d0)

  do istate=1,ctrl%nstates
     grad_diag_sad(istate,:,:)=traj%Gmatrix_ssad(istate,istate,:,:)  
     cgrad_MCH_sad(istate,:,:)=traj%grad_MCH_sad(istate,:,:)
  enddo

  cpNACdR_MCH_ssad=dcmplx(0.d0,0.d0)
  cNACdR_MCH_ssad=dcmplx(0.d0,0.d0)
  cpNACdR_MCH_ssad=traj%pNACdR_MCH_ssad
  cNACdR_MCH_ssad=traj%NACdR_ssad

  ! include laser fields
  H_ss=traj%H_MCH_ss
  if (ctrl%laser==2) then
    do idir=1,3
      H_ss=H_ss-traj%DM_ssd(:,:,idir)*real(ctrl%laserfield_td(traj%step*ctrl%nsubsteps+1,idir))
    enddo
  endif
  if (printlevel>4) then
    call matwrite(ctrl%nstates,H_ss,u_log,' H_MCH with laser field','F14.9')
  endif

! Compute pvib
    call computepvib(ctrl%natom, traj%geom_ad, traj%veloc_app_ad, traj%mass_a, pvib_ad)
  if (printlevel>4) then
    call vec3write(ctrl%natom,p,u_log,' Before removing the angular momentum, P:','F14.9')
    call vec3write(ctrl%natom,pvib_ad,u_log,' After removing the angular momentum, Pvib:','F14.9')
  endif

  if (ctrl%nac_projection==1) then ! use projected nac 
  ! Compute the s vecotr and decoherence time in diagonal basis
  select case (ctrl%neom)
    case (0)  ! ddr
      if (printlevel>4) then
        write(u_log,*) 's vector computed with projected ddr vector in diagonal basis...'
      endif
      call compute_svec_tau(ctrl%natom, ctrl%nstates, traj%veloc_app_ad, traj%mass_a, &
        &pvib_ad, traj%Etot, ctrl%decotime_method, ctrl%decoherence_alpha, &
        &ctrl%decoherence_beta, ctrl%gaussian_width, traj%pNACdR_diag_ssad, &
        &grad_diag_sad, traj%H_diag_ss, traj%state_diag, traj%svec_diag_sad, &
        &traj%decotime_diag_s, traj%decotime_diag_old_s)
    case (1)  ! effective nac
      if (printlevel>4) then
        write(u_log,*) 's vector computed with projected effective NAC in diagonal basis...'
      endif
      call compute_svec_tau_overlap(ctrl%natom, ctrl%nstates, traj%veloc_app_ad, &
        &ctrl%dtstep, traj%mass_a, pvib_ad, traj%Etot, ctrl%decotime_method, &
        &ctrl%decoherence_alpha, ctrl%decoherence_beta, ctrl%gaussian_width, &
        &traj%overlaps_ss, grad_diag_sad, traj%H_diag_ss, traj%state_diag, &
        &traj%svec_diag_sad, traj%decotime_diag_s, traj%decotime_diag_old_s, &
        &traj%pNACGV_diag_ssad)
  endselect

  ! Compute the s vecotr and decoherence time in MCH basis
  select case (ctrl%neom)
    case (0)  ! ddr
      if (printlevel>4) then
        write(u_log,*) 's vector computed with projected ddr vector in MCH basis...'
      endif
      call compute_svec_tau(ctrl%natom, ctrl%nstates, traj%veloc_app_ad, traj%mass_a, &
        &pvib_ad, traj%Etot, ctrl%decotime_method, ctrl%decoherence_alpha, &
        &ctrl%decoherence_beta, ctrl%gaussian_width, cpNACdR_MCH_ssad, &
        &cgrad_MCH_sad, H_ss, traj%state_MCH, traj%svec_MCH_sad, &
        &traj%decotime_MCH_s, traj%decotime_MCH_old_s)
    case (1)  ! effective nac
      if (printlevel>4) then
        write(u_log,*) 's vector computed with projected effective NAC in MCH basis...'
      endif
      call compute_svec_tau_overlap(ctrl%natom, ctrl%nstates, traj%veloc_app_ad, &
        &ctrl%dtstep, traj%mass_a, pvib_ad, traj%Etot, ctrl%decotime_method, &
        &ctrl%decoherence_alpha, ctrl%decoherence_beta, ctrl%gaussian_width, &
        &traj%overlaps_ss, cgrad_MCH_sad, H_ss, traj%state_MCH, &
        &traj%svec_MCH_sad, traj%decotime_MCH_s, traj%decotime_MCH_old_s, &
        &traj%pNACGV_MCH_ssad)
  endselect

  else ! use original nac
  ! Compute the s vecotr and decoherence time in diagonal basis
  select case (ctrl%neom)
    case (0)  ! ddr
      if (printlevel>4) then
        write(u_log,*) 's vector computed with ddr matrix in diagonal basis...'
      endif
      call compute_svec_tau(ctrl%natom, ctrl%nstates, traj%veloc_app_ad, traj%mass_a, &
        &pvib_ad, traj%Etot, ctrl%decotime_method, ctrl%decoherence_alpha, &
        &ctrl%decoherence_beta, ctrl%gaussian_width, traj%NACdR_diag_ssad, &
        &grad_diag_sad, traj%H_diag_ss, traj%state_diag, traj%svec_diag_sad, &
        &traj%decotime_diag_s, traj%decotime_diag_old_s)
    case (1)  ! effective nac
      if (printlevel>4) then
        write(u_log,*) 's vector computed with effective NAC in diagonal basis...'
      endif
      call compute_svec_tau_overlap(ctrl%natom, ctrl%nstates, traj%veloc_app_ad, &
        &ctrl%dtstep, traj%mass_a, pvib_ad, traj%Etot, ctrl%decotime_method, &
        &ctrl%decoherence_alpha, ctrl%decoherence_beta, ctrl%gaussian_width, &
        &traj%overlaps_ss, grad_diag_sad, traj%H_diag_ss, traj%state_diag, &
        &traj%svec_diag_sad, traj%decotime_diag_s, traj%decotime_diag_old_s, &
        &traj%NACGV_diag_ssad)
  endselect

! Compute the s vecotr and decoherence time in MCH basis
  select case (ctrl%neom)
    case (0)  ! ddr
      if (printlevel>4) then
        write(u_log,*) 's vector computed with ddr matrix in MCH basis...'
      endif
      call compute_svec_tau(ctrl%natom, ctrl%nstates, traj%veloc_app_ad, traj%mass_a, &
        &pvib_ad, traj%Etot, ctrl%decotime_method, ctrl%decoherence_alpha, &
        &ctrl%decoherence_beta, ctrl%gaussian_width, cNACdR_MCH_ssad, &
        &cgrad_MCH_sad, H_ss, traj%state_MCH, traj%svec_MCH_sad, traj%decotime_MCH_s, &
        &traj%decotime_MCH_old_s)
    case (1)  ! effective nac
      if (printlevel>4) then
        write(u_log,*) 's vector computed with effective NAC in MCH basis...'
      endif
      call compute_svec_tau_overlap(ctrl%natom, ctrl%nstates, traj%veloc_app_ad, &
        &ctrl%dtstep, traj%mass_a, pvib_ad, traj%Etot, ctrl%decotime_method, &
        &ctrl%decoherence_alpha, ctrl%decoherence_beta, ctrl%gaussian_width, &
        &traj%overlaps_ss, cgrad_MCH_sad, H_ss, traj%state_MCH, &
        &traj%svec_MCH_sad, traj%decotime_MCH_s, traj%decotime_MCH_old_s, &
        &traj%NACGV_MCH_ssad)
  endselect
  endif  ! endif of whether use projected or original nac

endsubroutine

! ==========================================================
!> this subroutine removes overall angular motion from momentum
!> the resulting momenta is placed in pvib
!> computepvib(ctrl%natom, traj%geom_ad, traj%veloc_app_ad, traj%mass_a, pvib_ad)
subroutine computepvib(n, r, v, m, pvib)
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: r(n,3), v(n,3), m(n)
  real*8, intent(out) :: pvib(n,3)

! p(n,3): momentum
! j(3): angular momentum; jmag: magnitude of angular momentum
! momi(3,3): moment of intertia
  integer :: idir, iatom
  integer :: arr, i, ik, ij, irel
  integer :: ipiv(3), infor
  real*8 :: q(n,3), p(n,3), jacob(3)
  real*8 :: jq(9), jp(9), jm(9), pvib2(9)
  real*8 :: pr(6), np(3), nq(6), Jint(6)
  real*8 :: prot(6), xdum(3,n)
  real*8 :: nj(6),nk(6),njnorm,nknorm,rr(3)
  real*8 :: summass, com(3), jmag, j(3)
  real*8 :: omega(3)
  real*8 :: momi(3,3), ivmomi(3,3)
  real*8 :: x, y, z
  real*8 :: m1, m2, qnorm, qnorm2, coschi
  real*8 :: npnorm, nqnorm, nqnorm2
  real*8 :: xk, erot, evib


  q=r

  do iatom=1,n
  do idir=1,3
    p(iatom,idir)=v(iatom,idir)*m(iatom)
  enddo
  enddo

  if (n==2) then
    pvib=p
  elseif (n==2) then ! PJsplit
    ! get Jacobis from cartesians
    ! get internals
    jacob=0.0
    do iatom=1,3
      jacob(1)=jacob(1)+(q(iatom,1)-q(iatom,2))**2
      jacob(2)=jacob(2)+(q(iatom,2)-q(iatom,3))**2
      jacob(3)=jacob(3)+(q(iatom,1)-q(iatom,3))**2
    enddo
    do iatom=1,3
      jacob(iatom)=sqrt(jacob(iatom))
    enddo
    ! pick which jacobis to use
    if (jacob(1).lt.jacob(2).and.jacob(1).lt.jacob(3)) arr=1
    if (jacob(2).lt.jacob(3).and.jacob(2).lt.jacob(3)) arr=2
    if (jacob(3).lt.jacob(1).and.jacob(3).lt.jacob(2)) arr=3
    if (arr.eq.1) then
      irel = 3
      ik = 1
      ij = 2
    endif
    if (arr.eq.2) then
      irel = 1
      ik = 2
      ij = 3
    endif 
    if (arr.eq.3) then
      irel = 2
      ik = 3
      ij = 1
    endif
    m1=m(irel)+m(ik)+m(ij)
    m2=m(ik)+m(ij)
    ! get Jac from R and P
    do i=1,3
      jq(i)=r(i,ij)-r(i,ik)
      jq(i+3)=r(i,irel)-(m(ik)*r(i,ik)+m(ij)*r(i,ij))/m2
      jq(i+6)=(m(irel)*r(i,irel)+m(ij)*r(i,ij)+m(ik)*r(i,ik))/m1
      jp(i)=-m(ij)/m2*p(i,ik)+m(ik)/m2*p(i,ij)
      jp(i+3)=m2/m1*p(i,irel)-m(irel)*(p(i,ik)+p(i,ij))/m1
      jp(i+6)=p(i,irel)+p(i,ik)+p(i,ij)
      jm(i)=m(ik)*m(ij)/m2
      jm(i+3)=m(irel)*m2/m1
      jm(i+6)=m1
    enddo
    ! transform q's and p's to mass weighted coordinates 
    do i = 1,6
      jp(i)=jp(i)/sqrt(jm(i))
      jq(i)=jq(i)*sqrt(jm(i))
    enddo
    qnorm=sqrt(jq(1)**2+jq(2)**2+jq(3)**2)
    qnorm2=sqrt(jq(4)**2+jq(5)**2 +jq(6)**2)
    coschi=(jq(1)*jq(4)+jq(2)*jq(5)+jq(3)*jq(6))/(qnorm*qnorm2)

    if (dabs(coschi) .ne. 1.0d0) then
      ! now we calculate vectors we need for the transformation
      ! np is a vector normal to the Q,q plane
      np(1)=(jq(2)*jq(6)-jq(3)*jq(5))
      np(2)=(jq(3)*jq(4)-jq(1)*jq(6))
      np(3)=(jq(1)*jq(5)-jq(2)*jq(4))
      npnorm=sqrt(np(1)**2+np(2)**2+np(3)**2)
      do i =1,3
        np(i)=np(i)/npnorm
      enddo
      ! nq(1-3) is a vector normal to q and np
      nq(1)=(np(2)*jq(3)-np(3)*jq(2))
      nq(2)=(np(3)*jq(1)-np(1)*jq(3))
      nq(3)=(np(1)*jq(2)-np(2)*jq(1))
      nqnorm=sqrt(nq(1)**2+nq(2)**2+nq(3)**2)
      do i =1,3
        nq(i) = nq(i)/nqnorm
      enddo
      ! nQ(4-6) is a vector normal to Q and np
      nq(4)=(np(2)*jq(6)-np(3)*jq(5))
      nq(5)=(np(3)*jq(4)-np(1)*jq(6))
      nq(6)=(np(1)*jq(5)-np(2)*jq(4))
      nqnorm2=sqrt(nq(4)**2+nq(5)**2+nq(6)**2)
      do i =4,6
        nq(i) = nq(i)/nqnorm2
      enddo
      ! we project p and P onto the 6 new vectors (q,nq,np,Q,nQ,np)
      pr(1) = (jp(1)*jq(1)+jp(2)*jq(2)+jp(3)*jq(3))/qnorm
      pr(2) = (jp(1)*nq(1)+jp(2)*nq(2)+jp(3)*nq(3))
      pr(3) = (jp(1)*np(1)+jp(2)*np(2)+jp(3)*np(3))
      pr(4) = (jP(4)*jQ(4)+jP(5)*jQ(5)+jP(6)*jQ(6))/Qnorm2
      pr(5) = (jP(4)*nQ(4)+jP(5)*nQ(5)+jP(6)*nQ(6))
      pr(6) = (jP(4)*np(1)+jP(5)*np(2)+jP(6)*np(3))
      ! now lets transform into internal/J coords . . 
      ! pr(1) and pr(4) are projections on r and R,
      ! pr(3) and pr(6) are out of plane motion.
      ! pr(2) and pr(5) are in plane motion . . .
      ! coords 1-4 have no J components
      Jint(1)=pr(1)
      Jint(2)=pr(4)
      Jint(3)=(Qnorm2*pr(2)-qnorm*coschi*pr(5))/sqrt(qnorm**2+Qnorm2**2)
      Jint(4)=(Qnorm2*pr(3)-qnorm*coschi*pr(6))/sqrt(qnorm**2+Qnorm2**2)
      ! coords 5-6 have no internal components
      Jint(5)=(qnorm*coschi*pr(2)+Qnorm2*pr(5))/sqrt(qnorm**2+Qnorm2**2)
      Jint(6)=(qnorm*coschi*pr(3)+Qnorm2*pr(6))/sqrt(qnorm**2+Qnorm2**2)
      ! calculatate J and energies (here's how)
      evib=0.0d0
      erot=0.0d0
      do i=1,3
        evib=evib+0.5d0*Jint(i)**2
        erot=erot+0.5d0*Jint(i+3)**2
      enddo
      xk=0.0d0
      do i=1,3
        xk = xk+(np(i)*Jint(4)*sqrt(qnorm**2+Qnorm2**2)-nq(i)*qnorm*Jint(5)-nQ(i+3)*Qnorm2*Jint(6))**2
      enddo
      xk=sqrt(xk)
      ! now transform back into prot and pvib
      ! 1st set Jint(1-3) to zero and transform 
      ! to get prot components . . .
      pr(1)=0.0d0
      pr(2)=(qnorm*Jint(4))/sqrt(qnorm**2+Qnorm2**2)
      pr(3)=Jint(5)
      pr(4)=0.0d0
      pr(5)=(Qnorm2*Jint(4))/sqrt(qnorm**2+Qnorm2**2)
      pr(6)=Jint(6)
      do i=1,3
        prot(i)=pr(1)*jq(i)/qnorm+pr(2)*nq(i)+pr(3)*np(i)
      enddo
      do i=4,6
        prot(i)=pr(4)*jq(i)/Qnorm2+pr(5)*nq(i)+pr(6)*np(i-3)
      enddo
      ! 2nd set Jint(4-6) to zero and transform 
      ! to get pvib components . . .
      pr(1)=Jint(1)
      pr(2)=(Qnorm2*Jint(3))/sqrt(qnorm**2+Qnorm2**2)
      pr(3)=0.0d0
      pr(4)=Jint(2)
      pr(5)=(-qnorm*Jint(3))/sqrt(qnorm**2+Qnorm2**2)
      pr(6)=0.0d0
      do i=1,3
        pvib2(i)=pr(1)*jq(i)/qnorm+pr(2)*nq(i)+pr(3)*np(i)
      enddo
      do i=4,6
        pvib2(i)=pr(4)*jq(i)/Qnorm2+pr(5)*nq(i)+pr(6)*np(i-3)
      enddo
      ! 3rd transform all Jint's back . . . should obtain original p's
      pr(1)=Jint(1)
      pr(2)=(qnorm*Jint(4)+Qnorm2*Jint(3))/sqrt(qnorm**2+Qnorm2**2)
      pr(3)=Jint(5)
      pr(4)=Jint(2)
      pr(5)=(Qnorm2*Jint(4)-qnorm*Jint(3))/sqrt(qnorm**2+Qnorm2**2)
      pr(6)=Jint(6)
      do i=1,3
        jp(i)=pr(1)*jq(i)/qnorm+pr(2)*nq(i)+pr(3)*np(i)
      enddo
      do i=4,6
        jp(i)=pr(4)*jq(i)/Qnorm2+pr(5)*nq(i)+pr(6)*np(i-3)
      enddo
      ! finally transform back to non mass scaled coordinates
      do i = 1,6
        jq(i)=jq(i)/sqrt(jm(i))
        jp(i)=jp(i)*sqrt(jm(i))
        prot(i)=prot(i)*sqrt(jm(i))
        pvib2(i)=pvib2(i)*sqrt(jm(i))
      enddo
    else !if (dabs(coschi) .ne. 1.0d0) then
    ! here the vectors are collinear, which means
    ! a different transformation procedure . . .
      if ((jq(2).eq.0.0d0).and.(jq(3).eq.0.0d0)) then
      ! q points along 1,0,0 . . . so
        nj(1)=0.0d0
        nj(2)=1.0d0
        nj(3)=0.0d0
        nk(1)=0.0d0
        nk(2)=0.0d0
        nk(3)=1.0d0
      else
      ! q does not point along 1,0,0.  So calculate nk = 1,0,0 x q
        nj(1)=0.0d0
        nj(2)=-jq(3)
        nj(3)=jq(2)
        njnorm=sqrt(jq(2)**2+jq(3)**2)
        do i=1,3
          nj(i)=nj(i)/njnorm
        enddo
      ! now calculate a vector nj . . .
        nk(1)=(jq(2)*nj(3)-jq(3)*nj(2))
        nk(2)=(jq(3)*nj(1)-jq(1)*nj(3))
        nk(3)=(jq(1)*nj(2)-jq(2)*nj(1))
        nknorm=sqrt(nk(1)**2+nk(2)**2+nk(3)**2)
        do i=1,3
          nk(i)=nk(i)/nknorm
        enddo
      endif
      ! we project p and P onto the 3 new vectors (q,nj,nk,Q,nj,nk)
      pr(1) = (jp(1)*jq(1)+jp(2)*jq(2)+jp(3)*jq(3))/qnorm
      pr(2) = (jp(1)*nj(1)+jp(2)*nj(2)+jp(3)*nj(3))
      pr(3) = (jp(1)*nk(1)+jp(2)*nk(2)+jp(3)*nk(3))
      pr(4) = (jP(4)*jQ(4)+jP(5)*jQ(5)+jP(6)*jQ(6))/Qnorm2
      pr(5) = (jP(4)*nj(1)+jP(5)*nj(2)+jP(6)*nj(3))
      pr(6) = (jP(4)*nk(1)+jP(5)*nk(2)+jP(6)*nk(3))
      ! now lets transform into internal/J coords . .
      ! coords 1-4 have no J components
      Jint(1)=pr(1)
      Jint(2)=pr(4)
      Jint(3)=(Qnorm2*pr(2)-qnorm*coschi*pr(5))/sqrt(qnorm**2+Qnorm2**2)
      Jint(4)=(Qnorm2*pr(3)-qnorm*coschi*pr(6))/sqrt(qnorm**2+Qnorm2**2)
      ! coords 5-6 have no internal components
      Jint(5)=(qnorm*coschi*pr(2)+Qnorm2*pr(5))/sqrt(qnorm**2+Qnorm2**2)
      Jint(6)=(qnorm*coschi*pr(3)+Qnorm2*pr(6))/sqrt(qnorm**2+Qnorm2**2)
      ! calculatate J and energies (here's how)
      evib=0.0d0
      erot=0.0d0
      do i=1,4
        evib=evib+0.5d0*Jint(i)**2
      enddo
      do i=5,6
        erot=erot+0.5d0*Jint(i)**2
      enddo
      xk=(Jint(5)**2+Jint(6)**2)*(qnorm**2+Qnorm2**2)
      xk=sqrt(xk)
      ! now transform back to prot and pvib
      ! 1st set Jint(1-4) to zero and transform
      ! to get prot components . . .
      pr(1)=0.0d0
      pr(2)=(qnorm*Jint(5)*coschi)/sqrt(qnorm**2+Qnorm2**2)
      pr(3)=(qnorm*Jint(6)*coschi)/sqrt(qnorm**2+Qnorm2**2)
      pr(4)=0.0d0
      pr(5)=(Qnorm2*Jint(5))/sqrt(qnorm**2+Qnorm2**2)
      pr(6)=(Qnorm2*Jint(6))/sqrt(qnorm**2+Qnorm2**2)
      do i=1,3
        prot(i)=pr(1)*jq(i)/qnorm+pr(2)*nj(i)+pr(3)*nk(i)
      enddo
      do i=4,6
        prot(i)=pr(4)*jq(i)/Qnorm2+pr(5)*nj(i-3)+pr(6)*nk(i-3)
      enddo
      ! 2nd set Jint(4-6) to zero and transform
      ! to get pvib components . . .
      pr(1)=Jint(1)
      pr(2)=(Qnorm2*Jint(3))/sqrt(qnorm**2+Qnorm2**2)
      pr(3)=(Qnorm2*Jint(4))/sqrt(qnorm**2+Qnorm2**2)
      pr(4)=Jint(2)
      pr(5)=(-qnorm*Jint(3)*coschi)/sqrt(qnorm**2+Qnorm2**2)
      pr(6)=(-qnorm*Jint(4)*coschi)/sqrt(qnorm**2+Qnorm2**2)
      do i=1,3
        pvib2(i)=pr(1)*jq(i)/qnorm+pr(2)*nj(i)+pr(3)*nk(i)
      enddo
      do i=4,6
        pvib2(i)=pr(4)*jq(i)/Qnorm2+pr(5)*nj(i-3)+pr(6)*nk(i-3)
      enddo
      ! 3rd transform all Jint's back . . . should obtain original p's
      pr(1)=Jint(1)
      pr(2)=(qnorm*Jint(5)*coschi+Qnorm2*Jint(3))/sqrt(qnorm**2+Qnorm2**2)
      pr(3)=(qnorm*Jint(6)*coschi+Qnorm2*Jint(4))/sqrt(qnorm**2+Qnorm2**2)
      pr(4)=Jint(2)
      pr(5)=(Qnorm2*Jint(5)-qnorm*Jint(3)*coschi)/sqrt(qnorm**2+Qnorm2**2)
      pr(6)=(Qnorm2*Jint(6)-qnorm*Jint(4)*coschi)/sqrt(qnorm**2+Qnorm2**2)
      do i=1,3
        jp(i)=pr(1)*jq(i)/qnorm+pr(2)*nj(i)+pr(3)*nk(i)
      enddo
      do i=4,6
        jp(i)=pr(4)*jq(i)/Qnorm2+pr(5)*nj(i-3)+pr(6)*nk(i-3)
      enddo
      ! finally transform back to non mass scaled coordinates
      do i = 1,6
        jq(i) = jq(i)/sqrt(jm(i))
        jp(i) = jp(i)*sqrt(jm(i))
        prot(i) = prot(i)*sqrt(jm(i))
        pvib2(i) = pvib2(i)*sqrt(jm(i))
      enddo
      do iatom=1,n
        do idir=1,3
          pvib(iatom,idir)=pvib2((iatom-1)*3+idir)
        enddo
      enddo
      endif
!######################
! for more than 3 atoms
!######################
  elseif (n>2) then
    ! center of mass
    com=0.d0
    summass=0.d0
    do iatom=1,n
      summass=summass+m(iatom)
      do idir=1,3
        com(idir)=com(idir)+q(iatom,idir)*m(iatom)
      enddo
    enddo
    com=com/summass
    ! coordinates relative to center of mass
    do iatom=1,n
      do idir=1,3
        q(iatom,idir)=q(iatom,idir)-com(idir)
      enddo
    enddo
    ! Angular momentum: J = Q X P
    ! Q is coordinate, P is momentum
    j=0.d0
    do iatom=1,n
      j(1)=j(1)+(q(iatom,2)*p(iatom,3)-q(iatom,3)*p(iatom,2))
      j(2)=j(2)-(q(iatom,1)*p(iatom,3)-q(iatom,3)*p(iatom,1))
      j(3)=j(3)+(q(iatom,1)*p(iatom,2)-q(iatom,2)*p(iatom,1))
    enddo

    ! compute the magnitude of J
    jmag=0.d0
    do idir=1,3
      jmag = jmag + j(idir)**2
    enddo
    jmag = dsqrt(jmag)

    ! compute moment of inertia
    momi(:,:)=0.d0
    do iatom=1,n
      x=q(iatom,1)
      y=q(iatom,2)
      z=q(iatom,3)
      momi(1,1)=momi(1,1)+m(iatom)*(y*y+z*z)
      momi(2,2)=momi(2,2)+m(iatom)*(x*x+z*z)
      momi(3,3)=momi(3,3)+m(iatom)*(x*x+y*y)
      momi(1,2)=momi(1,2)-m(iatom)*x*y
      momi(1,3)=momi(1,3)-m(iatom)*x*z
      momi(2,3)=momi(2,3)-m(iatom)*y*z
    enddo
    momi(2,1)=momi(1,2)
    momi(3,1)=momi(1,3)
    momi(3,2)=momi(2,3)
 
    do idir=1,3
      omega(idir)=j(idir)
    enddo

    call dgesv(3,1,momi,3,ipiv,omega,3,infor)
 
    do iatom=1,n
    do idir=1,3
      pvib(iatom,idir)=p(iatom,idir)
    enddo
    enddo

    ! pvib = p + corrections
    ! corrections = Q X Omega
    do iatom=1,n
      pvib(iatom,1)=pvib(iatom,1)+m(iatom)*(q(iatom,2)*omega(3)-q(iatom,3)*omega(2))
      pvib(iatom,2)=pvib(iatom,2)-m(iatom)*(q(iatom,1)*omega(3)-q(iatom,3)*omega(1))
      pvib(iatom,3)=pvib(iatom,3)+m(iatom)*(q(iatom,1)*omega(2)-q(iatom,2)*omega(1))
    enddo
  endif

endsubroutine

! ==========================================================
!> this subroutine computes the s vector
!> compute_svec_tau(ctrl%natom, ctrl%nstates, traj%veloc_app_ad, traj%mass_a,
!> pvib_ad, traj%Etot, ctrl%decotime_method, ctrl%decoherence_alpha, ctrl%decoherence_beta, ctrl%gaussian_width, traj%NACdR_diag_ssad, grad_diag_sad, traj%H_diag_ss, traj%state_diag, svec_sad, traj%decotime_s, traj%decotime_old_s)
subroutine compute_svec_tau(n, ns, v, m, pvib, etot, decomethod, alpha, beta, sigma, g, f, e, k, s, t, told)
  use definitions, only: u_log, printlevel
  implicit none
  integer, intent(in) :: n, ns
  integer, intent(in) :: k
  integer, intent(in) :: decomethod
  real*8, intent(in) :: etot
  real*8, intent(in) :: v(n,3), m(n), pvib(n,3)
  real*8, intent(in) :: alpha, beta, sigma
  complex*16, intent(in) :: g(ns,ns,n,3), f(ns,n,3), e(ns,ns)
  real*8, intent(out) :: s(ns,n,3)
  real*8, intent(inout) :: t(ns), told(ns)

  integer :: idir, iatom, istate, jstate
  real*8 :: dvec(ns,ns,n,3)
  complex*16 :: pdotd
  real*8 :: phase
  real*8 :: dvecmag, smag
  real*8 :: us(ns), es
  real*8 :: ekin(ns), trajekin
! decoherence time. Notice t is 1/(decoherence time)
  real*8 :: rt(ns)
  real*8 :: pi

  complex*16 :: gu(ns,ns,n,3)
  real*8 :: gv(ns), gg(ns), mgg(ns)
  real*8 :: lambda, e_available(ns)

  real*8 :: ps(ns,n,3), p(n,3)
  complex*16 :: df(n,3), dp(n,3), sp(n,3)
  real*8 :: dfdotd, dpdotd, spdotd, dpdotdf

  pi=4.d0*datan(1.d0)

! k is the current pointer state
! g is the NAC in diagonal basis, e is the Hamiltonian in diag basis
  do istate=1,ns
  do jstate=1,ns
    dvec(istate,jstate,:,:)=real(g(istate,jstate,:,:))
  enddo
  enddo

! compute the kinetic energy for each state
  do istate=1,ns
    ekin(istate)=etot-real(e(istate,istate))
  enddo

! compute the trajectory momentum
  trajekin=0.d0
  do iatom=1,n
  do idir=1,3
    p(iatom,idir)=v(iatom,idir)*m(iatom)
    trajekin=trajekin+0.5*m(iatom)*v(iatom,idir)*v(iatom,idir)
  enddo
  enddo

! compute the g direction vector gu
  gg=0.d0
  gu=dcmplx(0.d0,0.d0)
  do istate=1,ns
    if (istate .ne. k) then 
      do iatom=1,n
        do idir=1,3
          gg(istate)=gg(istate)+g(istate,k,iatom,idir)*conjg(g(istate,k,iatom,idir))
        enddo
      enddo
      if (gg(istate) .ne. 0.d0) then
        gu(istate,k,:,:)=g(istate,k,:,:)/gg(istate)
        gu(k,istate,:,:)=g(k,istate,:,:)/gg(istate)
      else
        gu(istate,k,:,:)=0.d0
        gu(k,istate,:,:)=0.d0
      endif
    endif
  enddo

! compute gudotv and mass caled gu dot gu, notice gu is complex, but only real part is considered
  gv=0.d0
  mgg=0.d0
  do istate=1,ns
    do iatom=1,n
      do idir=1,3 
        gv(istate)=gv(istate)+real(gu(istate,k,iatom,idir)*v(iatom,idir))
        mgg(istate)=mgg(istate)+real(gu(istate,k,iatom,idir)**2)/m(iatom)
      enddo
    enddo
  enddo

! compute e_available
  do istate=1,ns
    e_available(istate)=gv(istate)**2+2*mgg(istate)*(ekin(istate)-trajekin)
  enddo

! compute the momentum for each state by hopping procedure
  do istate=1,ns
    if (istate.eq.k) then ! for pointer state, momentum is rescaled on velocity direction 
      ps(istate,:,:)=p
!      if (ekin(istate).ge.0.00) then
!        ps(istate,:,:)=sqrt((ekin(istate)/trajekin))*p(:,:)
!      else
!        ps(istate,:,:)=0.d0
!        ps(istate,:,:)=-p(:,:)
!      endif
    else ! not pointer state, velocty is rescaled on NAC/effective NAC direction 
      if (e_available(istate).ge.0.d0) then
        if (gv(istate).gt.0.d0) then 
          lambda=(gv(istate)-sqrt(e_available(istate)))/gg(istate)
          ps(istate,:,:)=p(:,:)-lambda*real(gu(istate,k,:,:))
        else if (gv(istate).lt.0.d0) then
          lambda=(gv(istate)+sqrt(e_available(istate)))/gg(istate)
          ps(istate,:,:)=p(:,:)-lambda*real(gu(istate,k,:,:))
        else if (gv(istate).eq.0.d0) then ! velocty is rescaled on velocity direction, because NAC=0.d0
          ps(istate,:,:)=sqrt((ekin(istate)/trajekin))*p(:,:)
        endif
      else 
        ps(istate,:,:)=0.d0
!        ps(istate,:,:)=-p(:,:)
      endif
    endif 
  enddo

  do istate=1,ns
    told(istate)=t(istate)
  enddo

! initialize svec
  s(:,:,:)=0.d0

  do istate=1,ns
    if (istate.eq.k) then
      t(istate)=0.d0
    else
      ! first compute pdotd = p dot dvec/|devc|
      pdotd=0.d0
      dvecmag=0.d0
      do iatom=1,n
      do idir=1,3
        pdotd=pdotd + g(istate,k,iatom,idir)*pvib(iatom,idir)/m(iatom)
        dvecmag=dvecmag+(g(istate,k,iatom,idir)**2)/m(iatom)
      enddo
      enddo
      if (dvecmag.lt.1.d-8) then
        pdotd = 0.d0
      else
        dvecmag=sqrt(dvecmag)
        pdotd=pdotd/dvecmag
      endif
      ! now compute svec(istate,iatom,idir)=dvec*pdotd+ppvib
      smag=0.d0
      do iatom=1,n      
      do idir=1,3
        s(istate,iatom,idir)=real(g(istate,k,iatom,idir)*pdotd)+pvib(iatom,idir)
        smag=smag+(s(istate,iatom,idir)**2)/m(iatom)
      enddo
      enddo
      if (smag .le. 0.d0) then
        write(u_log,*)"the magnitude of s vecor is negative, calculation stop"
        stop 1
      endif
      smag=sqrt(smag)

      ! select the method to compute decoherence time
      select case (decomethod)
        case (0) ! CSDM decoherence time
          ! compute energy along s vector (es)
          us(istate)=0.d0
          do iatom=1,n
            do idir=1,3
              us(istate)=us(istate)+s(istate,iatom,idir)*v(iatom,idir)
            enddo
          enddo
          es=(us(istate)**2)/(2.d0*smag**2)
          ! compute the decoherence time, see eq.(28) in JCP 121,7658 (2004)
          ! notice this t is the inverse of tau_iK in eq.(28)
          ! in the following calculation, C=1.0d0, E0=alpha
          t(istate)=dabs(real(e(istate,istate)-e(k,k)))/(beta+alpha/es)
        case (1) ! SCDM decoherence time
          ! compute energy along s vector (es)
          us(istate)=0.d0
          do iatom=1,n
            do idir=1,3
              us(istate)=us(istate)+s(istate,iatom,idir)*v(iatom,idir)
            enddo
          enddo
          es=(us(istate)**2)/(2.d0*smag**2)
          if (dabs(real(e(istate,istate)-e(k,k)))>0.d0) then 
            rt(istate)=beta/dabs(real(e(istate,istate)-e(k,k)))+alpha/es
            t(istate)=1.d0/rt(istate)
          else if (dabs(real(e(istate,istate)-e(k,k)))==0.d0) then
            t(istate)=0.d0
          endif
        case (2) ! edc
          t(istate)=dabs(real(e(istate,istate)-e(k,k)))/(beta+alpha/trajekin)          
        case (3) ! sd
          phase=1.d0
          if (real(pdotd) .lt. 0.d0) then
            phase=-1.d0
          endif
          df(:,:)=f(istate,:,:)-f(k,:,:)
          dp(:,:)=ps(istate,:,:)-ps(k,:,:)
          sp(:,:)=ps(istate,:,:)+ps(k,:,:)
          dfdotd=0.d0
          dpdotd=0.d0
          spdotd=0.d0
          do iatom=1,n
            do idir=1,3
              dfdotd=dfdotd+phase*real(df(iatom,idir)*gu(istate,k,iatom,idir))
              dpdotd=dpdotd+real(dp(iatom,idir)*gu(istate,k,iatom,idir))/sqrt(m(iatom))
              spdotd=spdotd+real(sp(iatom,idir)*gu(istate,k,iatom,idir))
            enddo
          enddo
          rt(istate)=pi*dfdotd/spdotd+sqrt((dpdotd/(2*pi))**2*(dabs(real(e(istate,istate)-e(k,k))))+(pi*dfdotd/spdotd)**2)
          t(istate)=rt(istate)
        case(4) !fp1
          df(:,:)=f(istate,:,:)-f(k,:,:)
          dp(:,:)=ps(istate,:,:)-ps(k,:,:)
          sp(:,:)=ps(istate,:,:)+ps(k,:,:)
          rt(istate)=0.d0
          do iatom=1,n
            do idir=1,3
              rt(istate)=rt(istate)+real(sp(iatom,idir)/(2*pi*df(iatom,idir)))
            enddo
          enddo
          t(istate)=1.d0/rt(istate)
        case(5) !fp2
          df(:,:)=f(istate,:,:)-f(k,:,:)
          dp(:,:)=ps(istate,:,:)-ps(k,:,:)
          sp(:,:)=ps(istate,:,:)+ps(k,:,:)
          rt(istate)=0.d0
          do iatom=1,n
            do idir=1,3
              rt(istate)=rt(istate)+m(iatom)*real(e(istate,istate)-e(k,k))/(pi*df(iatom,idir)*dp(iatom,idir))
            enddo
          enddo
          t(istate)=1.d0/rt(istate)
      endselect

    endif
  enddo

endsubroutine

! ==========================================================
!> this subroutine computes the s vector
!> compute_svec_tau_overlap(ctrl%natom, ctrl%nstates, traj%veloc_app_ad, ctrl%dt,
!> traj%mass_a, pvib_ad, traj%Etot, ctrl%decotime_method,
!> ctrl%decoherence_alpha, ctrl%decoherence_beta, ctrl%gaussian_width, traj%overlaps_ss, grad_diag_sad, traj%H_diag_ss, traj%state_diag,
!> svec_sad, traj%decotime_s, traj%decotime_old_s, traj%NACGV_ssad)
subroutine compute_svec_tau_overlap(n, ns, v, dt, m, pvib, etot, decomethod, alpha, beta, sigma, g, f, e, k, s, t, told, gv)
  use definitions, only: u_log, printlevel
  use matrix
  implicit none
  integer, intent(in) :: n, ns
  integer, intent(in) :: k
  integer, intent(in) :: decomethod
  real*8, intent(in) :: etot
  real*8, intent(in) :: dt
  real*8, intent(in) :: v(n,3), m(n), pvib(n,3)
  real*8, intent(in) :: alpha, beta, sigma
  complex*16, intent(in) :: gv(ns,ns,n,3), f(ns,n,3)
  complex*16, intent(in) :: g(ns,ns),e(ns,ns)
  real*8, intent(out) :: s(ns,n,3)
  real*8, intent(inout) :: t(ns), told(ns)

  integer :: idir, iatom, istate, jstate
  real*8 :: dvec(ns,ns,n,3)
  complex*16 :: pdotd
  real*8 :: phase
  real*8 :: pvibmag, gvmag, smag
  real*8 :: us(ns), es
  real*8 :: ekin(ns), trajekin
! decoherence time. Notice t is 1/(decoherence time)
  real*8 :: rt(ns)
  real*8 :: pi

  complex*16 :: gu(ns,ns,n,3)
  real*8 :: gdotv(ns), gg(ns), mgg(ns)
  real*8 :: lambda, e_available(ns)

  real*8 :: ps(ns,n,3), p(n,3)
  complex*16 :: df(n,3), dp(n,3), sp(n,3)
  real*8 :: dfdotd, dpdotd, spdotd, dpdotdf

  pi=4.d0*datan(1.d0)

! k is the current pointer state
! g is the overlap matrix (time derivative * dt), e is the Hamiltonian in diag basis
  do istate=1,ns
  do jstate=1,ns
    dvec(istate,jstate,:,:)=real(gv(istate,jstate,:,:))
  enddo
  enddo

! compute the kinetic energy for each state
  do istate=1,ns
    ekin(istate)=etot-real(e(istate,istate))
  enddo

! compute the trajectory momentum 
  trajekin=0.d0
  do iatom=1,n
  do idir=1,3
    p(iatom,idir)=v(iatom,idir)*m(iatom)
    trajekin=trajekin+0.5*m(iatom)*v(iatom,idir)*v(iatom,idir)
  enddo
  enddo

! compute the gv direction vector gu
  gg=0.d0
  gu=dcmplx(0.d0,0.d0)
  do istate=1,ns
    if (istate .ne. k) then
      do iatom=1,n
        do idir=1,3
          gg(istate)=gg(istate)+gv(istate,k,iatom,idir)*conjg(gv(istate,k,iatom,idir))
        enddo
      enddo
      if (gg(istate) .ne. 0.d0) then
        gu(istate,k,:,:)=gv(istate,k,:,:)/gg(istate)
        gu(k,istate,:,:)=gv(k,istate,:,:)/gg(istate)
      else 
        gu(istate,k,:,:)=0.d0
        gu(k,istate,:,:)=0.d0
      endif
    endif
  enddo

! compute gudotv and mass caled gu dot gu, notice gu is complex, but only real
! part is considered
  gdotv=0.d0
  mgg=0.d0
  do istate=1,ns
    do iatom=1,n
      do idir=1,3
        gdotv(istate)=gdotv(istate)+real(gu(istate,k,iatom,idir)*v(iatom,idir))
        mgg(istate)=mgg(istate)+real(gu(istate,k,iatom,idir)**2)/m(iatom)
      enddo
    enddo
  enddo

! compute e_available
  do istate=1,ns
    e_available(istate)=gdotv(istate)**2+2*mgg(istate)*(ekin(istate)-trajekin)
  enddo

! compute the momentum for each state by hopping procedure
  do istate=1,ns
    if (istate.eq.k) then ! for pointer state, momentum is rescaled on velocity direction 
      ps(istate,:,:)=p
!      if (ekin(istate).ge.0.00) then
!        ps(istate,:,:)=sqrt((ekin(istate)/trajekin))*p(:,:)
!      else
!        ps(istate,:,:)=0.d0
!      endif
    else ! not pointer state, velocty is rescaled on NAC/effective NAC direction 
      if (e_available(istate).ge.0.d0) then
        if (gdotv(istate).gt.0.d0) then
          lambda=(gdotv(istate)-sqrt(e_available(istate)))/gg(istate)
          ps(istate,:,:)=p(:,:)-lambda*real(gu(istate,k,:,:))
        else if (gdotv(istate).lt.0.d0) then
          lambda=(gdotv(istate)+sqrt(e_available(istate)))/gg(istate)
          ps(istate,:,:)=p(:,:)-lambda*real(gu(istate,k,:,:))
         else if (gdotv(istate).eq.0.d0) then ! rescale along velocity direction, becasue NAC=0.d0
          ps(istate,:,:)=sqrt((ekin(istate)/trajekin))*p(:,:)
        endif
      else
        ps(istate,:,:)=0.d0
      endif
    endif
  enddo

  do istate=1,ns
    told(istate)=t(istate)
  enddo

! initialize svec
  s(:,:,:)=0.d0

  do istate=1,ns
    if (istate.eq.k) then
      t(istate)=0.d0
    else
      pdotd=0.d0
      gvmag=0.d0
      do iatom=1,n
      do idir=1,3
        pdotd=pdotd + gv(istate,k,iatom,idir)*pvib(iatom,idir)/m(iatom)
        gvmag=gvmag+(gv(istate,k,iatom,idir)**2)/m(iatom)
      enddo
      enddo
      if (gvmag.lt.1.d-8) then
        pdotd = 0.d0
      else
        gvmag=sqrt(gvmag)
        pdotd=pdotd/gvmag
      endif
      ! now compute svec(istate,iatom,idir)=gv*S(istate,jstate)
      smag=0.d0
      do iatom=1,n
      do idir=1,3
        s(istate,iatom,idir)=real(gv(istate,k,iatom,idir)*pdotd)+pvib(iatom,idir)
!        s(istate,iatom,idir)=((g(istate,k)-g(k,istate)/2/dt))*pvib(iatom,idir)/pvibmag+pvib(iatom,idir)
        smag=smag+(s(istate,iatom,idir)**2)/m(iatom)
      enddo
      enddo
      if (smag .le. 0.d0) then
        write(u_log,*)"the magnitude of s vecor is negative, calculation stop"
        stop 1
      endif
      smag=sqrt(smag)

      ! select the method to compute decoherence time
      select case (decomethod)
        case (0) ! CSDM decoherence time
          ! compute energy along s vector (es)
          us(istate)=0.d0
          do iatom=1,n
          do idir=1,3
            us(istate)=us(istate)+s(istate,iatom,idir)*v(iatom,idir)
          enddo
          enddo
          es=(us(istate)**2)/(2.d0*smag**2)
          ! compute the decoherence time, see eq.(28) in JCP 121,7658 (2004)
          ! notice this t is the inverse of tau_iK in eq.(28)
          ! in the following calculation, C=1.0d0, E0=0.1d0
          t(istate)=dabs(real(e(istate,istate)-e(k,k)))/(beta+alpha/es)
        case (1) ! SCDM decoherence time
          ! compute energy along s vector (es)
          us(istate)=0.d0
          do iatom=1,n
          do idir=1,3
            us(istate)=us(istate)+s(istate,iatom,idir)*v(iatom,idir)
          enddo
          enddo
          es=(us(istate)**2)/(2.d0*smag**2)
          if (dabs(real(e(istate,istate)-e(k,k)))>0.d0) then
            rt(istate)=beta/dabs(real(e(istate,istate)-e(k,k)))+alpha/es
            t(istate)=1.d0/rt(istate)
          else if (dabs(real(e(istate,istate)-e(k,k)))==0.d0) then 
            t(istate)=0.d0
          endif
        case (2) ! edc
          t(istate)=dabs(real(e(istate,istate)-e(k,k)))/(beta+alpha/trajekin)
        case (3) ! sd
          phase=1.d0
          if (real(pdotd) .lt. 0.d0) then 
            phase=-1.d0
          endif
          df(:,:)=f(istate,:,:)-f(k,:,:)
          dp(:,:)=ps(istate,:,:)-ps(k,:,:)
          sp(:,:)=ps(istate,:,:)+ps(k,:,:)
          dfdotd=0.d0
          dpdotd=0.d0
          spdotd=0.d0
          do iatom=1,n
          do idir=1,3
            dfdotd=dfdotd+phase*real(df(iatom,idir)*gu(istate,k,iatom,idir))
            dpdotd=dpdotd+real(dp(iatom,idir)*gu(istate,k,iatom,idir))/sqrt(m(iatom))
            spdotd=spdotd+real(sp(iatom,idir)*gu(istate,k,iatom,idir))
          enddo
          enddo
          rt(istate)=pi*dfdotd/spdotd+sqrt((dpdotd/(2*pi))**2*(dabs(real(e(istate,istate)-e(k,k))))+(pi*dfdotd/spdotd)**2)
          t(istate)=rt(istate)
        case (4) ! fp1
          df(:,:)=f(istate,:,:)-f(k,:,:)
          dp(:,:)=ps(istate,:,:)-ps(k,:,:)
          sp(:,:)=ps(istate,:,:)+ps(k,:,:)
          rt(istate)=0.d0
          do iatom=1,n
          do idir=1,3
            rt(istate)=rt(istate)+real(sp(iatom,idir)/(2*pi*df(iatom,idir)))
          enddo
          enddo
          t(istate)=1.d0/rt(istate)
        case(5) !fp2
          df(:,:)=f(istate,:,:)-f(k,:,:)
          dp(:,:)=ps(istate,:,:)-ps(k,:,:)
          sp(:,:)=ps(istate,:,:)+ps(k,:,:)
          rt(istate)=0.d0
          do iatom=1,n
          do idir=1,3
            rt(istate)=rt(istate)+m(iatom)*real(e(istate,istate)-e(k,k))/(pi*df(iatom,idir)*dp(iatom,idir))
          enddo
          enddo
          t(istate)=1.d0/rt(istate)
      endselect

    endif
  enddo

endsubroutine


! ==========================================================
!> this subroutine performs energy based decay of mixing
subroutine DoM_coeff(traj,ctrl)
  use definitions
  use matrix
  use decoherence_afssh
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl
  integer :: istate
  real*8 :: sumc
  complex*16 :: c(ctrl%nstates)

  sumc=0.d0
  do istate=1,ctrl%nstates
    if (istate/=traj%state_diag) then
      c(istate)=traj%coeff_diag_s(istate) * exp( -ctrl%dtstep * traj%decotime_diag_s(istate) )
      sumc=sumc+abs(c(istate))**2
    endif
  enddo

  if (abs(traj%coeff_diag_s(traj%state_diag))**2 > 0.d0) then
    c(traj%state_diag)=traj%coeff_diag_s(traj%state_diag) * &
    &sqrt( (1.d0-sumc) / abs(traj%coeff_diag_s(traj%state_diag))**2)
  else
    c(traj%state_diag)=dcmplx(0.d0,0.d0)
  endif

  traj%coeff_diag_s=c

endsubroutine


! =========================================================
!> decoherent_propagator(&
!>    &ctrl%nstates, ctrl%dtstep, ctrl%nsubsteps, &
!>    &traj%decotime_diag_s, traj%decotime_diag_old_s, traj%state_diag, &
!>    &traj%coeff_diag_s, traj%Dtotal_ss)
subroutine decoherent_propagator(ns, dt, nsub, tau, tauold, k, C, Dtotal)
  use definitions, only: u_log
  use matrix
  implicit none
  integer, intent(in) :: ns
  real*8, intent(in) :: dt
  integer, intent(in) :: nsub
  real*8, intent(in) :: tau(ns), tauold(ns)
  integer, intent(in) :: k
  complex*16, intent(inout) :: C(ns)
  complex*16, intent(out) :: Dtotal(ns,ns)

  integer :: istep
  real*8 :: dtsubstep
  complex*16 :: Rtotal(ns,ns)
  complex*16 :: DP(ns,ns)
  complex*16 :: Rexpd(ns,ns)
  complex*16 :: Rprod(ns,ns)
  complex*16 :: Rtmp(ns,ns)
  complex*16 :: Ctmp1(ns)
  complex*16 :: C_in(ns)
  complex*16 :: ii=dcmplx(0.d0,1.d0)
  real*8 :: taumid(ns)

  integer :: istate, jstate

  dtsubstep=dt/nsub
  Dtotal=dcmplx(0.d0,0.d0)
  do istate=1,ns
    Dtotal(istate,istate)=dcmplx(1.d0,0.d0)
  enddo
  Ctmp1(:)=C(:)
  C_in(:)=C(:)

  do istep=1,nsub
    taumid=tauold+ (tau-tauold)*istep/nsub
    ! compute decoherence propagator DP
    DP(:,:)=dcmplx(0.d0,0.d0)
    do istate=1,ns
      if (istate.ne.k) then
        DP(k,k)=DP(k,k)+(1.d0,0.d0)*tau(istate)*real(Ctmp1(istate)*conjg(Ctmp1(istate)))
        DP(istate,istate)=-0.5d0*(1.d0,0.d0)*tau(istate)
      endif
    enddo
    DP(k,k)=0.5d0*DP(k,k)/(real(Ctmp1(k)*conjg(Ctmp1(k))))

    Rexpd=dtsubstep*(DP)
    call exponentiate(ns,Rexpd,(1.d0,0.d0))
    call matmultiply(ns, Rexpd, Dtotal, Rprod, 'nn')
    Dtotal=Rprod

    ! compute coeff in substeps 
    Rtmp=Dtotal
    call matvecmultiply(ns, Rtmp, C_in, Ctmp1, 'n')
  enddo

  ! propagate the coefficients
  call matvecmultiply(ns, Dtotal, C_in, C, 'n')

endsubroutine


! ==========================================================
!> this subroutine propagates the decoherence part of the electronic coefficients 
!> repropagate_coeff(&
!    &ctrl%nstates, ctrl%dtstep, ctrl%nsubsteps, interpolation, &
!    &traj%decotime_s, traj%decotime_old_s, traj%state_diag,&
!    &traj%coeff_diag_s, traj%coeff_diag_old_s,& 
!    &traj%U_ss, traj%U_old_ss,&
!    &traj%H_MCH_ss, traj%H_MCH_old_ss,&
!    &traj%NACdt_ss, traj%NACdt_old_ss,&
!    &traj%RDtotal_ss, Dtotal_ss)
!> DC is the current decoherence coefficient, which is at time t as input
!> DC is also the output as decoherence coefficient at time t+dt
!> In this subroutine, we need to locally propagate C from t to t+0.5dt 
subroutine repropagate_coeff(ns, dt, nsub, interp, tau, tauold, k, C, Cold, U, Uold, SO, SOold, NACM, NACMold, RD, Dtotal)
!use electronic
  use definitions, only: u_log
  use matrix
  implicit none
  integer, intent(in) :: ns
  real*8, intent(in) :: dt
  integer, intent(in) :: nsub
  integer, intent(in) :: interp
  real*8, intent(in) :: tau(ns), tauold(ns)
  integer, intent(in) :: k
  complex*16, intent(in) :: Cold(ns)
  complex*16, intent(in) :: U(ns,ns), Uold(ns,ns)
  complex*16, intent(in) :: SO(ns,ns), SOold(ns,ns)
  complex*16, intent(in) :: NACM(ns,ns), NACMold(ns,ns)
  complex*16, intent(inout) :: C(ns)
  complex*16, intent(out) :: RD(ns,ns)
  complex*16, intent(out) :: Dtotal(ns,ns)

  integer :: istep
  real*8 :: dtsubstep
  complex*16 :: Rtotal(ns,ns)
  complex*16 :: H(ns,ns), T(ns,ns), DP(ns,ns)
  complex*16 :: Rexpc(ns,ns), Rexpd(ns,ns)
  complex*16 :: Rexp(ns,ns), Rprod(ns,ns)
  complex*16 :: Rtmp(ns,ns)
  complex*16 :: Ctmp1(ns)
  complex*16 :: ii=dcmplx(0.d0,1.d0)
  complex*16 :: dentmp(ns,ns),densum
  real*8 :: taumid(ns)

  integer :: istate, jstate

  !initialize Ctmp1
  call matvecmultiply(ns, Uold, Cold, Ctmp1, 'n')

  dtsubstep=dt/nsub
  Rtotal=dcmplx(0.d0,0.d0)
  do istate=1,ns
    Rtotal(istate,istate)=dcmplx(1.d0,0.d0)
  enddo
  call matmultiply(ns, Uold, Rtotal, Rprod, 'nn')
  Rtotal=Rprod

  Dtotal=dcmplx(0.d0,0.d0)
  do istate=1,ns
    Dtotal(istate,istate)=dcmplx(1.d0,0.d0)
  enddo
  call matmultiply(ns, Uold, Dtotal, Rprod, 'nn')
  Dtotal=Rprod

  do istep=1,nsub

    H=SOold+ (SO-SOold)*istep/nsub
    if (interp==1) then 
      T=NACM
    else 
      T=NACMold+ (NACM-NACMold)*istep/nsub
    endif 
    taumid=tauold+ (tau-tauold)*istep/nsub
    ! compute decoherence propagator DP
    DP(:,:)=dcmplx(0.d0,0.d0)
    do istate=1,ns
      if (istate.ne.k) then
        DP(k,k)=DP(k,k)+(1.d0,0.d0)*tau(istate)*real(Ctmp1(istate)*conjg(Ctmp1(istate)))
        DP(istate,istate)=-0.5d0*(1.d0,0.d0)*tau(istate)
      endif
    enddo
    DP(k,k)=0.5d0*DP(k,k)/(real(Ctmp1(k)*conjg(Ctmp1(k))))

    Rexpc=dtsubstep*(H-ii*T)
    Rexpd=dtsubstep*(DP)

    call exponentiate(ns,Rexpc,-ii)
    call exponentiate(ns,Rexpd,(1.d0,0.d0))

    call matmultiply(ns, Rexpc, Rexpd, Rexp, 'nn')

    call matmultiply(ns, Rexp, Rtotal, Rprod, 'nn')
    Rtotal=Rprod
 
    call matmultiply(ns, Rexpd, Dtotal, Rprod, 'nn')
    Dtotal=Rprod

    ! compute coeff in substeps in MCH representation
    Rtmp=Rtotal
    call matvecmultiply(ns, Rtmp, Cold, Ctmp1, 'n')
  
  enddo
  call matmultiply(ns, U, Rtotal, Rprod, 'tn')
  Rtotal=Rprod
  RD=Rtotal

  call matmultiply(ns, U, Dtotal, Rprod, 'tn')
  Dtotal=Rprod

  call matvecmultiply(ns, Rtotal, Cold, C, 'n')

endsubroutine

! ==========================================================
!> this subroutine propagates the decoherence part of the electronic coefficients 
!> repropagate_coeff_overlap(&
!    &ctrl%nstates, ctrl%dtstep, ctrl%nsubsteps, ctrl%coupling,&
!    &traj%decotime_s, traj%decotime_old_s, traj%state_diag,&
!    &traj%coeff_diag_s, traj%coeff_diag_old_s,& 
!    &traj%U_ss, traj%U_old_ss,&
!    &traj%H_MCH_ss, traj%H_MCH_old_ss,&
!    &traj%overlaps_ss,&
!    &traj%RDtotal_ss, Dtotal_ss)
!> DC is the current decoherence coefficient, which is at time t as input
!> DC is also the output as decoherence coefficient at time t+dt
!> In this subroutine, we need to locally propagate C from t to t+0.5dt 
subroutine repropagate_coeff_NPI(ns, dt, nsub, coup, tau, tauold, k, C, Cold, U, Uold, SO, SOold, overlap, RD, Dtotal)
!use electronic
  use definitions, only: u_log
  use matrix
  implicit none
  integer, intent(in) :: ns
  real*8, intent(in) :: dt
  integer, intent(in) :: nsub
  integer, intent(in) :: coup
  real*8, intent(in) :: tau(ns), tauold(ns)
  integer, intent(in) :: k
  complex*16, intent(in) :: Cold(ns)
  complex*16, intent(in) :: U(ns,ns), Uold(ns,ns)
  complex*16, intent(in) :: SO(ns,ns), SOold(ns,ns)
  complex*16, intent(in) :: overlap(ns,ns)
  complex*16, intent(inout) :: C(ns)
  complex*16, intent(out) :: RD(ns,ns)
  complex*16, intent(out) :: Dtotal(ns,ns)

  integer :: istep
  real*8 :: dtsubstep
  complex*16 :: Rtotal(ns,ns)
  complex*16 :: SOt(ns,ns)
  complex*16 :: H(ns,ns), T(ns,ns), DP(ns,ns)
  complex*16 :: Rexpc(ns,ns), Rexpd(ns,ns)
  complex*16 :: Rexp(ns,ns), Rprod(ns,ns)
  complex*16 :: Rtmp(ns,ns)
  complex*16 :: Ctmp1(ns)
  complex*16 :: ii=dcmplx(0.d0,1.d0)
  complex*16 :: w(ns,ns),tw(ns,ns),dw(ns,ns)
  real*8 :: taumid(ns)

  integer :: istate, jstate

  !initialize Ctmp1
  call matvecmultiply(ns, Uold, Cold, Ctmp1, 'n')

  dtsubstep=dt/nsub
  Rtotal=dcmplx(0.d0,0.d0)
  do istate=1,ns
    Rtotal(istate,istate)=dcmplx(1.d0,0.d0)
  enddo
  call matmultiply(ns, Uold, Rtotal, Rprod, 'nn')
  Rtotal=Rprod

  Dtotal=dcmplx(0.d0,0.d0)
  do istate=1,ns
    Dtotal(istate,istate)=dcmplx(1.d0,0.d0)
  enddo
  call matmultiply(ns, Uold, Dtotal, Rprod, 'nn')
  Dtotal=Rprod

  !initialize T
  T=dcmplx(0.d0,0.d0)

  do istep=1,nsub
    H=SOold+ (SO-SOold)*istep/nsub
    taumid=tauold+ (tau-tauold)*istep/nsub
    ! compute decoherence propagator DP
    DP(:,:)=dcmplx(0.d0,0.d0)
    do istate=1,ns
      if (istate.ne.k) then
        DP(k,k)=DP(k,k)+(1.d0,0.d0)*tau(istate)*real(Ctmp1(istate)*conjg(Ctmp1(istate)))
        DP(istate,istate)=-0.5d0*(1.d0,0.d0)*tau(istate)
      endif
    enddo
    DP(k,k)=0.5d0*DP(k,k)/(real(Ctmp1(k)*conjg(Ctmp1(k))))

    ! compute the NPI rotation matrix W
    do istate=1,ns
      do jstate=1,ns
        if (jstate .eq. istate) then
          w(istate,jstate)=cos(acos(overlap(istate,jstate))*istep/nsub)
          tw(jstate,istate)=cos(acos(overlap(istate,jstate))*istep/nsub)
          dw(istate,jstate)=-sin(acos(overlap(istate,jstate))*istep/nsub)*acos(overlap(istate,jstate))/dt
        else
          w(istate,jstate)=sin(asin(overlap(istate,jstate))*istep/nsub)
          tw(jstate,istate)=sin(asin(overlap(istate,jstate))*istep/nsub)
          dw(istate,jstate)=cos(asin(overlap(istate,jstate))*istep/nsub)*asin(overlap(istate,jstate))/dt
        end if
      enddo
    enddo
    !call matwrite(ns,w,u_log,'substep rotation matrix','F14.9')    
    call matmultiply(ns, tw, dw, T, 'nn')
    !call matwrite(ns,T,u_log,'substep time derivative coupling','F14.9')
    do istate=1,ns
      T(istate,istate)=dcmplx(0.d0,0.d0)
    enddo
    !call matwrite(ns,T,u_log,'substep time derivative coupling','F14.9')

    Rexpc=dtsubstep*(H-ii*T)
    Rexpd=dtsubstep*(DP)
    call exponentiate(ns,Rexpc,-ii)
    call exponentiate(ns,Rexpd,(1.d0,0.d0))
    call matmultiply(ns, Rexpc, Rexpd, Rexp, 'nn')

    call matmultiply(ns, Rexp, Rtotal, Rprod, 'nn')
    Rtotal=Rprod

    call matmultiply(ns, Rexpd, Dtotal, Rprod, 'nn')
    Dtotal=Rprod

    ! compute coeff in substeps in MCH representation
    Rtmp=Rtotal
    call matvecmultiply(ns, Rtmp, Cold, Ctmp1, 'n')
  enddo
  call matmultiply(ns, U, Rtotal, Rprod, 'tn')
  Rtotal=Rprod
  RD=Rtotal

  call matmultiply(ns, U, Dtotal, Rprod, 'tn')
  Dtotal=Rprod

  call matvecmultiply(ns, Rtotal, Cold, C, 'n')

endsubroutine

! ==========================================================
!> this subroutine propagates the decoherence part of the electronic coefficients 
!> including laser fields
!> repropagate_coeff_laser(&
!    &ctrl%nstates, ctrl%dtstep, ctrl%nsubsteps, interpolation,&
!    &traj%decotime_s, traj%decotime_old_s, traj%state_diag,&
!    &traj%coeff_diag_s, traj%coeff_diag_old_s,& 
!    &traj%U_ss, traj%U_old_ss,&
!    &traj%H_MCH_ss, traj%H_MCH_old_ss,&
!    &traj%NACdt_ss, traj%NACdt_old_ss,&
!    &traj%DM_ssd,traj%DM_old_ssd,&
!    &ctrl%laserfield_td( (traj%step-1)*ctrl%nsubsteps+2:traj%step*ctrl%nsubsteps+1 ,:),&
!    &traj%RDtotal_ss, Dtotal_ss)
!> DC is the current decoherence coefficient, which is at time t as input
!> DC is also the output as decoherence coefficient at time t+dt
!> In this subroutine, we need to locally propagate C from t to t+0.5dt 
subroutine repropagate_coeff_laser(ns, dt, nsub, interp, tau, tauold, k, C, Cold, U, Uold, SO, SOold, NACM, NACMold, DM, DMold, laserfield, RD, Dtotal)
!use electronic
  use definitions, only: u_log
  use matrix
  implicit none
  integer, intent(in) :: ns
  real*8, intent(in) :: dt
  integer, intent(in) :: nsub
  integer, intent(in) :: interp
  real*8, intent(in) :: tau(ns), tauold(ns)
  integer, intent(in) :: k
  complex*16, intent(in) :: Cold(ns)
  complex*16, intent(in) :: U(ns,ns), Uold(ns,ns)
  complex*16, intent(in) :: SO(ns,ns), SOold(ns,ns)
  complex*16, intent(in) :: NACM(ns,ns), NACMold(ns,ns)
  complex*16, intent(in) :: DM(ns,ns,3),DMold(ns,ns,3)
  complex*16, intent(in) :: laserfield(nsub,3)
  complex*16, intent(inout) :: C(ns)
  complex*16, intent(out) :: RD(ns,ns)
  complex*16, intent(out) :: Dtotal(ns,ns)

  integer :: istep
  real*8 :: dtsubstep
  complex*16 :: Rtotal(ns,ns)
  complex*16 :: H(ns,ns), T(ns,ns), DP(ns,ns)
  complex*16 :: Rexpc(ns,ns), Rexpd(ns,ns)
  complex*16 :: Rexp(ns,ns), Rprod(ns,ns)
  complex*16 :: Rtmp(ns,ns)
  complex*16 :: Ctmp1(ns)
  complex*16 :: ii=dcmplx(0.d0,1.d0)
  complex*16 :: dentmp(ns,ns),densum
  real*8 :: taumid(ns)


  integer :: istate, jstate, ixyz

  !initialize Ctmp1
  call matvecmultiply(ns, Uold, Cold, Ctmp1, 'n')

  dtsubstep=dt/nsub
  Rtotal=dcmplx(0.d0,0.d0)
  do istate=1,ns
    Rtotal(istate,istate)=dcmplx(1.d0,0.d0)
  enddo
  call matmultiply(ns, Uold, Rtotal, Rprod, 'nn')
  Rtotal=Rprod

  Dtotal=dcmplx(0.d0,0.d0)
  do istate=1,ns
    Dtotal(istate,istate)=dcmplx(1.d0,0.d0)
  enddo
  call matmultiply(ns, Uold, Dtotal, Rprod, 'nn')
  Dtotal=Rprod

  do istep=1,nsub

    H=SOold+(SO-SOold)*istep/nsub
    ! here the laser field is added to the Hamiltonian
    do ixyz=1,3
      H=H-(DMold(:,:,ixyz)+(DM(:,:,ixyz)-DMold(:,:,ixyz))*istep/nsub)*real(laserfield(istep,ixyz))
    enddo
 
    if (interp==1) then
      T=NACM
    else 
      T=NACMold+ (NACM-NACMold)*istep/nsub
    endif 

    taumid=tauold+ (tau-tauold)*istep/nsub
    ! compute decoherence propagator DP
    DP(:,:)=dcmplx(0.d0,0.d0)
    do istate=1,ns
      if (istate.ne.k) then
        DP(k,k)=DP(k,k)+(1.d0,0.d0)*tau(istate)*real(Ctmp1(istate)*conjg(Ctmp1(istate)))
        DP(istate,istate)=-0.5d0*(1.d0,0.d0)*tau(istate)
      endif
    enddo
    DP(k,k)=0.5d0*DP(k,k)/(real(Ctmp1(k)*conjg(Ctmp1(k))))

    Rexpc=dtsubstep*(H-ii*T)
    Rexpd=dtsubstep*(DP)
    call exponentiate(ns,Rexpc,-ii)
    call exponentiate(ns,Rexpd,(1.d0,0.d0))
    call matmultiply(ns, Rexpc, Rexpd, Rexp, 'nn')

    call matmultiply(ns, Rexp, Rtotal, Rprod, 'nn')
    Rtotal=Rprod

    call matmultiply(ns, Rexpd, Dtotal, Rprod, 'nn')
    Dtotal=Rprod

    ! compute coeff in substeps in MCH representation
    Rtmp=Rtotal
    call matvecmultiply(ns, Rtmp, Cold, Ctmp1, 'n')

  enddo
  call matmultiply(ns, U, Rtotal, Rprod, 'tn')
  Rtotal=Rprod
  RD=Rtotal

  call matmultiply(ns, U, Dtotal, Rprod, 'tn')
  Dtotal=Rprod

  call matvecmultiply(ns, Rtotal, Cold, C, 'n')

endsubroutine

! ==========================================================
!> this subroutine propagates the decoherence part of the electronic coefficients 
!> including laser fields
!> repropagate_coeff_overlap_laser(&
!    &ctrl%nstates, ctrl%dtstep, ctrl%nsubsteps, ctrl%coupling,&
!    &traj%decotime_s, traj%decotime_old_s, traj%state_diag,&
!    &traj%coeff_diag_s, traj%coeff_diag_old_s,& 
!    &traj%U_ss, traj%U_old_ss,&
!    &traj%H_MCH_ss, traj%H_MCH_old_ss,&
!    &traj%overlaps_ss,&
!    &traj%DM_ssd,traj%DM_old_ssd,&
!    &ctrl%laserfield_td((traj%step-1)*ctrl%nsubsteps+2:traj%step*ctrl%nsubsteps+1 ,:),&
!    &traj%RDtotal_ss, Dtotal_ss)
!> DC is the current decoherence coefficient, which is at time t as input
!> DC is also the output as decoherence coefficient at time t+dt
!> In this subroutine, we need to locally propagate C from t to t+0.5dt 
subroutine repropagate_coeff_NPI_laser(ns, dt, nsub, coup, tau, tauold, k, C, Cold, U, Uold, SO, SOold, overlap, DM, DMold, laserfield, RD, Dtotal)
!use electronic
  use definitions, only: u_log
  use matrix
  implicit none
  integer, intent(in) :: ns
  real*8, intent(in) :: dt
  integer, intent(in) :: nsub
  integer, intent(in) :: coup
  real*8, intent(in) :: tau(ns), tauold(ns)
  integer, intent(in) :: k
  complex*16, intent(in) :: Cold(ns)
  complex*16, intent(in) :: U(ns,ns), Uold(ns,ns)
  complex*16, intent(in) :: SO(ns,ns), SOold(ns,ns)
  complex*16, intent(in) :: overlap(ns,ns)
  complex*16, intent(in) :: DM(ns,ns,3),DMold(ns,ns,3)
  complex*16, intent(in) :: laserfield(nsub,3)
  complex*16, intent(inout) :: C(ns)
  complex*16, intent(out) :: RD(ns,ns)
  complex*16, intent(out) :: Dtotal(ns,ns)

  integer :: istep
  real*8 :: dtsubstep
  complex*16 :: Rtotal(ns,ns)
  complex*16 :: SOt(ns,ns)
  complex*16 :: H(ns,ns), T(ns,ns), DP(ns,ns)
  complex*16 :: Rexpc(ns,ns), Rexpd(ns,ns)
  complex*16 :: Rexp(ns,ns), Rprod(ns,ns)
  complex*16 :: Rtmp(ns,ns)
  complex*16 :: Ctmp1(ns)
  complex*16 :: ii=dcmplx(0.d0,1.d0)
  complex*16 :: w(ns,ns),tw(ns,ns),dw(ns,ns)
  real*8 :: taumid(ns)

  integer :: istate, jstate, ixyz

  !initialize Ctmp1
  call matvecmultiply(ns, Uold, Cold, Ctmp1, 'n')

  dtsubstep=dt/nsub
  Rtotal=dcmplx(0.d0,0.d0)
  do istate=1,ns
    Rtotal(istate,istate)=dcmplx(1.d0,0.d0)
  enddo
  call matmultiply(ns, Uold, Rtotal, Rprod, 'nn')
  Rtotal=Rprod

  Dtotal=dcmplx(0.d0,0.d0)
  do istate=1,ns
    Dtotal(istate,istate)=dcmplx(1.d0,0.d0)
  enddo
  call matmultiply(ns, Uold, Dtotal, Rprod, 'nn')
  Dtotal=Rprod

  !initialize T
  T=dcmplx(0.d0,0.d0)

  do istep=1,nsub

    H=SOold+ (SO-SOold)*istep/nsub
    ! here the laser field is added to the Hamiltonian
    do ixyz=1,3
      H=H-(DMold(:,:,ixyz)+(DM(:,:,ixyz)-DMold(:,:,ixyz))*istep/nsub)*real(laserfield(istep,ixyz))
    enddo

    taumid=tauold+ (tau-tauold)*istep/nsub

    ! compute decoherence propagator DP
    DP(:,:)=dcmplx(0.d0,0.d0)
    do istate=1,ns
      if (istate.ne.k) then
        DP(k,k)=DP(k,k)+(1.d0,0.d0)*tau(istate)*real(Ctmp1(istate)*conjg(Ctmp1(istate)))
        DP(istate,istate)=-0.5d0*(1.d0,0.d0)*tau(istate)
      endif
    enddo
    DP(k,k)=0.5d0*DP(k,k)/(real(Ctmp1(k)*conjg(Ctmp1(k))))

    ! compute the NPI rotation matrix W
    do istate=1,ns
      do jstate=1,ns
        if (jstate .eq. istate) then
          w(istate,jstate)=cos(acos(overlap(istate,jstate))*istep/nsub)
          tw(jstate,istate)=cos(acos(overlap(istate,jstate))*istep/nsub)
          dw(istate,jstate)=-sin(acos(overlap(istate,jstate))*istep/nsub)*acos(overlap(istate,jstate))/dt
        else
          w(istate,jstate)=sin(asin(overlap(istate,jstate))*istep/nsub)
          tw(jstate,istate)=sin(asin(overlap(istate,jstate))*istep/nsub)
          dw(istate,jstate)=cos(asin(overlap(istate,jstate))*istep/nsub)*asin(overlap(istate,jstate))/dt
        end if
      enddo
    enddo
    !call matwrite(ns,w,u_log,'substep rotation matrix','F14.9')    
    call matmultiply(ns, tw, dw, T, 'nn')
    !call matwrite(ns,T,u_log,'substep time derivative coupling','F14.9')
    do istate=1,ns
      T(istate,istate)=dcmplx(0.d0,0.d0)
    enddo
    !call matwrite(ns,T,u_log,'substep time derivative coupling','F14.9')

    Rexpc=dtsubstep*(H-ii*T)
    Rexpd=dtsubstep*(DP)
    call exponentiate(ns,Rexpc,-ii)
    call exponentiate(ns,Rexpd,(1.d0,0.d0))
    call matmultiply(ns, Rexpc, Rexpd, Rexp, 'nn')

    call matmultiply(ns, Rexp, Rtotal, Rprod, 'nn')
    Rtotal=Rprod

    call matmultiply(ns, Rexpd, Dtotal, Rprod, 'nn')
    Dtotal=Rprod

    ! compute coeff in substeps in MCH representation
    Rtmp=Rtotal
    call matvecmultiply(ns, Rtmp, Cold, Ctmp1, 'n')
  enddo
  call matmultiply(ns, U, Rtotal, Rprod, 'tn')
  Rtotal=Rprod
  RD=Rtotal

  call matmultiply(ns, U, Dtotal, Rprod, 'tn')
  Dtotal=Rprod

  call matvecmultiply(ns, Rtotal, Cold, C, 'n')

endsubroutine


! ==========================================================
!> this subroutine propagates the decoherence part of the electronic coefficients 
!> use Runge-Kutta 4 integrator. 
!> integrate_dcoeff_dt(&
!    &ctrl%nstates, ctrl%dtstep, ctrl%nsubsteps, ctrl%coupling,&
!    &traj%decotime_s, traj%decotime_old_s, tauhalf, traj%state_diag,&
!    &traj%coeff_diag_s, traj%coeff_diag_old_s,& 
!    &traj%U_ss, traj%U_old_ss,&
!    &traj%H_MCH_ss, traj%H_MCH_old_ss,&
!    &traj%NACdt_ss, traj%NACdt_old_ss, traj%overlaps_ss,&
!    &traj%dcoeff_diag_s, traj%dcoeff_diag_old_s)
!> DC is the current decoherence coefficient, which is at time t as input
!> DC is also the output as decoherence coefficient at time t+dt
!> In this subroutine, we need to locally propagate C from t to t+0.5dt 
subroutine propagate_dcoeff(ns, dt, nsub, coup, tau, tauold, tauhalf, k, C, Cold, U, Uold, SO, SOold, NACM, NACMold, Overlap, DC, DCold)
!use electronic
  use definitions, only: u_log
  use matrix
  implicit none
  integer, intent(in) :: ns
  real*8, intent(in) :: dt
  integer, intent(in) :: nsub
  integer, intent(in) :: coup
  real*8, intent(in) :: tau(ns), tauold(ns), tauhalf(ns)
  integer, intent(in) :: k
  complex*16, intent(in) :: C(ns), Cold(ns)
  complex*16, intent(in) :: U(ns,ns), Uold(ns,ns)
  complex*16, intent(in) :: SO(ns,ns), SOold(ns,ns)
  complex*16, intent(in) :: NACM(ns,ns), NACMold(ns,ns), Overlap(ns,ns)
  complex*16, intent(inout) :: DC(ns), DCold(ns)

  integer :: istep
  real*8 :: dtsubstep
  complex*16 :: H(ns,ns), T(ns,ns)
  complex*16 :: Rexp(ns,ns), Rprod(ns,ns)
  complex*16 :: ii=dcmplx(0.d0,1.d0)

  integer :: istate, jstate
  complex*16 :: Chalf(ns)
  real*8 :: dthalf
  complex*16 :: Rhalf(ns,ns)
  complex*16 :: k1(ns), k2(ns), k3(ns), k4(ns)
  complex*16 :: DCk1(ns), DCk2(ns), DCk3(ns), DCk4(ns)
  complex*16 :: DCdtk1(ns), DCdtk2(ns), DCdtk3(ns), DCdtk4(ns)

! propagate Cold(t) to Cold(t+0.5dt), call it Chalf
  dthalf=dt/2.d0

  Rhalf=dcmplx(0.d0,0.d0)
  do istate=1,ns
    Rhalf(istate,istate)=dcmplx(1.d0,0.d0)
  enddo

  dtsubstep=dthalf/nsub
  call matmultiply(ns,Uold,Rhalf,Rprod,'nn')
  Rhalf=Rprod
  do istep=1,nsub
    H=SOold+ (SO-SOold)*istep/nsub
    T=NACM
    Rexp=dtsubstep*(H-ii*T)
    call exponentiate(ns,Rexp,-ii)
    call matmultiply(ns,Rexp,Rhalf,Rprod,'nn')
    Rhalf=Rprod
  enddo
  call matmultiply(ns,U,Rhalf,Rprod,'tn')
  Rhalf=Rprod
  call matvecmultiply(ns, Rhalf, Cold, Chalf, 'n')

  do istate=1,ns
    DCold(istate)=DC(istate)
  enddo

! use RK4 to integrate DC
! first step
  call compute_dcoeff_dt(ns, Cold, DC, tauold, k, DCdtk1)
  do istate=1,ns
    DCk1(istate)=DC(istate)+0.5*dt*DCdtk1(istate)
  enddo
! second step
  call compute_dcoeff_dt(ns, Chalf, DCk1, tauhalf, k, DCdtk2)
  do istate=1,ns
    DCk2(istate)=DC(istate)+0.5*dt*DCdtk2(istate)
  enddo
! third step
  call compute_dcoeff_dt(ns, Chalf, DCk2, tauhalf, k, DCdtk3)
  do istate=1,ns
    DCk3(istate)=DC(istate)+dt*DCdtk3(istate)
    DCdtk3(istate)=DCdtk2(istate)+DCdtk3(istate)
  enddo
! fourth step
  call compute_dcoeff_dt(ns, C, DCk3, tau, k, DCdtk4)
  do istate=1,ns
    DC(istate)=DC(istate)+dt/6.d0*(DCdtk1(istate)+DCdtk4(istate)+2.d0*DCdtk3(istate))
  enddo

endsubroutine

! ==========================================================
!> this subroutine computes the time derivative of the decoherence part of the
!electronic coefficients 
!> compute_dcoeff_dt(ctrl%nstates, traj%coeff_diag_old_s, traj%dcoeff_diag_s,
!tau_s, traj%state_diag, dcoeff_diag_dt_s)
subroutine compute_dcoeff_dt(ns, c, dc, tau, k, dcdt)
  implicit none
  integer, intent(in) :: ns
  integer, intent(in) :: k
  complex*16, intent(in) :: c(ns), dc(ns)
  real*8, intent(in) :: tau(ns)
  complex*16, intent(out) :: dcdt(ns)

  integer :: istate
  complex*16 :: ctotal(ns)

  do istate=1,ns
    ctotal(istate)=c(istate)+dc(istate)
  enddo

  dcdt(k)=(0.d0,0.d0)

  do istate=1,ns
    if (istate.ne.k) then
      dcdt(k)=dcdt(k)+tau(istate)*real(ctotal(istate)*conjg(ctotal(istate)))
      dcdt(istate)=-0.5d0*tau(istate)*ctotal(istate)
    endif
  enddo

  dcdt(k)=0.5d0*dcdt(k)*ctotal(k)/(real(ctotal(k)*conjg(ctotal(k))))

endsubroutine

! ==========================================================
!> this subroutine computes the decoherence force
!> compute_decoforce(ctrl%natom, ctrl%nstates, traj%veloc_app_ad, traj%coeff_diag_s,&
!>   &svec_sad, tau_s, traj%state_diag, traj%H_diag_ss, traj%decograd_ad)
subroutine compute_decoforce(n, ns, v, c, s, t, k, h, def)
  use definitions, only: u_log
  implicit none
  integer, intent(in) :: n, ns
  integer, intent(in) :: k
  real*8, intent(in) :: v(n,3)
  complex*16, intent(in) :: c(ns)
  real*8, intent(in) :: s(ns,n,3)
  real*8, intent(in) :: t(ns)
  complex*16, intent(in) :: h(ns,ns)
  real*8, intent(out) :: def(n,3)

  integer :: istate, jstate
  integer :: iatom, idir

  real*8 :: vd(ns)
  real*8 :: df
  real*8 :: ps(ns)
  complex*16 :: den(ns,ns) 


  do istate=1,ns
    vd(istate)=0.d0
  enddo

  do istate=1,ns
    do jstate=1,ns
      den(istate,jstate)=c(istate)*conjg(c(jstate))
    enddo
  enddo 

  do istate=1,ns
      ps(istate)=0.d0
  enddo

  do istate=1,ns
    do iatom=1,n
      do idir=1,3
        ps(istate)=ps(istate)+s(istate,iatom,idir)*v(iatom,idir)
      enddo
    enddo
  enddo

  do istate=1,ns
    if (istate .ne. k) then 
      vd(istate)=vd(istate)-real(den(istate,istate)*t(istate)*h(istate,istate))
    elseif (istate .eq. k) then
      do jstate=1,ns
        if (jstate .ne. istate) then
          vd(jstate)=vd(jstate)+real(den(jstate,jstate)*t(jstate)*h(istate,istate))
        endif
      enddo
    endif
  enddo

  do idir=1,3
    do iatom=1,n
      def(iatom, idir)=0.d0
      do istate=1,ns
        if (istate .ne. k) then
          def(iatom, idir)=def(iatom,idir)-vd(istate)*s(istate,iatom,idir)/ps(istate)
        endif
      enddo
    enddo
  enddo

endsubroutine
 
endmodule
