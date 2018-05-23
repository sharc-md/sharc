!******************************************
!
!    SHARC Program Suite
!
!    Copyright (c) 2018 University of Vienna
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

!> # Module DECOHERENCE A-FSSH
!>
!> \author Felix Plasser
!> \date 10.08.2017
!>
!> This module contains subroutines that allow to perform a decoherence
!> correction according to the approximated augmented fewest-switches
!> surface hopping method as described in
!> Jain, Alguire, Subotnik JCTC 2016, 12, 5256.

module decoherence_afssh
  implicit none

  private

  public :: afssh_step, allocate_afssh

contains

! ===========================================================

subroutine allocate_afssh(traj,ctrl)
    use definitions
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl
    integer :: allocst
    integer :: natom, nstates, istate

    natom=ctrl%natom
    nstates = ctrl%nstates

    if (printlevel>1) then
        write(u_log,*) 'Initializing A-FSSH algorithm'
    endif

    if (ctrl%coupling/=2) then
      write(u_log,*) 'ERROR: A-FSSH only implemented for "coupling overlap"'
      stop 'ERROR: A-FSSH only implemented for "coupling overlap"'
    endif

    allocate(traj%auxtrajs_s(nstates), stat=allocst)
    if (allocst/=0) stop 'Could not allocate auxtrajs_s'

    do istate = 1, nstates
      allocate(traj%auxtrajs_s(istate)%mass_a(natom),stat=allocst)
      if (allocst/=0) stop 'Could not allocate mass_a'
      traj%auxtrajs_s(istate)%mass_a = traj%mass_a

      allocate(traj%auxtrajs_s(istate)%geom_ad(natom,3),stat=allocst)
      if (allocst/=0) stop 'Could not allocate geom_ad'

      allocate(traj%auxtrajs_s(istate)%veloc_ad(natom,3),stat=allocst)
      if (allocst/=0) stop 'Could not allocate veloc_ad'

      allocate(traj%auxtrajs_s(istate)%geom_tmp_ad(natom,3),stat=allocst)
      if (allocst/=0) stop 'Could not allocate geom_tmp_ad'

      allocate(traj%auxtrajs_s(istate)%veloc_tmp_ad(natom,3),stat=allocst)
      if (allocst/=0) stop 'Could not allocate veloc_tmp_ad'

      allocate(traj%auxtrajs_s(istate)%accel_ad(natom,3),stat=allocst)
      if (allocst/=0) stop 'Could not allocate accel_ad'

      allocate(traj%auxtrajs_s(istate)%grad_ad(natom,3),stat=allocst)
      if (allocst/=0) stop 'Could not allocate grad_ad'

      traj%auxtrajs_s(istate)%istate = istate
    enddo

    call reset_moments(traj,ctrl)
endsubroutine

! ===========================================================

subroutine reset_moments(traj,ctrl)
    use definitions
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl
    integer :: istate

    if (printlevel>2) then
      write(u_log, *)'Resetting all A-FSSH moments'
    endif
    do istate = 1, ctrl%nstates
      traj%auxtrajs_s(istate)%geom_ad = 0.d0
      traj%auxtrajs_s(istate)%veloc_ad = 0.d0
      traj%auxtrajs_s(istate)%accel_ad = 0.d0
      traj%auxtrajs_s(istate)%grad_ad = 0.d0
    enddo

endsubroutine

! ===========================================================

subroutine afssh_step(traj,ctrl)
  use definitions
  use matrix
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

  integer :: nstates, istate
  complex*16, dimension(ctrl%nstates, ctrl%nstates) :: trans_mat, Hnew, tmp

  nstates = ctrl%nstates

  if (printlevel>2) then
    write(u_log,*) 'Computing A-FSSH decoherence'
  endif

  if ((traj%kind_of_jump==1).or.(traj%kind_of_jump==3))then
    call reset_moments(traj,ctrl)
    return ! do not propagate moments in case of hop, correct(?)
  endif

  ! Compute the diag_old / diag transformation matrix
  call matmultiply(nstates, traj%U_old_ss, traj%overlaps_ss, tmp, 'tn')
  call matmultiply(nstates, tmp, traj%U_ss, trans_mat, 'nn')
  if (.not.isunitary(nstates, trans_mat)) then
    write(u_log,*) 'ERROR: A-FSSH transformation matrix not unitary'
    stop 'ERROR: A-FSSH transformation matrix not unitary'
  endif
  if (printlevel>2) then
    call matwrite(nstates, trans_mat, u_log, 'A-FSSH transformation matrix', 'F8.4')
  endif

  ! compute the new Hamiltonian in the old basis
  Hnew = traj%H_MCH_ss
  call transform(nstates,Hnew,traj%overlaps_ss,'uaut')
  call transform(nstates,Hnew,traj%U_old_ss, 'utau')
  if (printlevel>2) then
    call matwrite(nstates, Hnew, u_log, 'New Hamiltonian in the old basis', 'F8.4')
  endif

  ! Propagate the auxiliary trajectories
  do istate = 1, nstates
    call afssh_prop(traj, traj%auxtrajs_s(istate), ctrl, trans_mat)
  enddo

  ! Transform back to the instantaneous adiabatic basis
  do istate = 1, nstates
    call afssh_transform(traj, traj%auxtrajs_s(istate), ctrl, trans_mat)
  enddo

  ! Compute the rates and perform the collapsing / resetting steps
  do istate = 1, nstates
    call afssh_rates(traj, traj%auxtrajs_s(istate), ctrl, Hnew)
  enddo

endsubroutine

! ===========================================================

subroutine afssh_prop(traj, atraj, ctrl, trans_mat)
  use definitions
  implicit none
  type(trajectory_type) :: traj
  type(aux_trajectory_type) :: atraj
  type(ctrl_type) :: ctrl
  complex*16 :: trans_mat(ctrl%nstates, ctrl%nstates)

  integer :: istate, jstate, nstates
  real*8 :: sigma
  real*8 :: ref_grad(ctrl%natom, 3)
  real*8, dimension(ctrl%nstates) :: trans_row

  istate = atraj%istate
  nstates = ctrl%nstates

  ref_grad = traj%grad_ad !real(traj%Gmatrix_ssad(traj%state_diag,traj%state_diag,:,:))
  sigma = traj%coeff_diag_old_s(istate) ** 2
  trans_row = abs(trans_mat(istate, :)) ** 2

  ! Propagate the xstep with the density-matrix-weighted gradient from the last time step
  call VelocityVerlet_xstep_afssh(atraj, ctrl, sigma)

  ! Compute the new density-matrix-weighted gradient
  atraj%grad_ad = -ref_grad
  do jstate=1,nstates
    atraj%grad_ad = atraj%grad_ad + trans_row(jstate) * real(traj%Gmatrix_ssad(jstate,jstate,:,:))
!    print *, istate, jstate, trans_row(jstate), sum(atraj%grad_ad*atraj%grad_ad)
  enddo

  ! Propagate
  call VelocityVerlet_vstep_afssh(atraj, ctrl, sigma)

  ! store the density-matrix-weighted gradient to be used in the next time step
  atraj%grad_ad  = (real(traj%Gmatrix_ssad(istate,istate,:,:)) - ref_grad)
endsubroutine
! ===========================================================

subroutine afssh_transform(traj, atraj, ctrl, trans_mat)
  use definitions
  implicit none
  type(trajectory_type) :: traj
  type(aux_trajectory_type) :: atraj
  type(ctrl_type) :: ctrl
  complex*16 :: trans_mat(ctrl%nstates, ctrl%nstates)

  integer :: istate, jstate, nstates
  real*8, dimension(ctrl%nstates) :: trans_col

  istate = atraj%istate
  nstates = ctrl%nstates

  trans_col = abs(trans_mat(:, istate)) ** 2

  atraj%geom_ad = 0.d0
  atraj%veloc_ad = 0.d0
  do jstate = 1, nstates
    atraj%geom_ad  = atraj%geom_ad  + trans_col(jstate)*traj%auxtrajs_s(jstate)%geom_tmp_ad
    atraj%veloc_ad = atraj%veloc_ad + trans_col(jstate)*traj%auxtrajs_s(jstate)%veloc_tmp_ad
  enddo
endsubroutine

! ===========================================================

!> compute the two parts of the decoherence rate
subroutine afssh_rates(traj, atraj, ctrl, Hnew)
  use definitions
  implicit none
  type(trajectory_type) :: traj
  type(aux_trajectory_type) :: atraj
  type(ctrl_type) :: ctrl
  complex*16 :: Hnew(:,:)

  integer istate, rstate
  real*8 :: geom_disp(ctrl%natom, 3)
  complex*16 :: clam, cn

  real*8 :: tmprate2

  istate = atraj%istate
  rstate = traj%state_diag
  geom_disp = atraj%geom_ad - traj%auxtrajs_s(traj%state_diag)%geom_ad

  ! F_nn . (R_nn - R_ll) / 2 hbar
  !  Note that F_nn = -grad_ad
  atraj%rate1 = -0.5 * sum(atraj%grad_ad * geom_disp)

  tmprate2 = 2 * abs(Hnew(rstate,istate) / ctrl%dtstep * sum(geom_disp * traj%veloc_ad)/&
    &sum(traj%veloc_ad * traj%veloc_ad))
  ! rate2
  ! use the original formula using NACs if they are available
  if (ctrl%gradcorrect==1) then  !ctrl%coupling==0
    atraj%rate2 = 2 * abs( sum(real(traj%Gmatrix_ssad(rstate,istate,:,:)) * geom_disp) )
  ! consider off-diagonal Hamiltonian element
  elseif (ctrl%coupling==2) then
    atraj%rate2 = tmprate2
  else
    stop 'A-FSSH not possible without NAC vectors or overlaps!'
  endif

  if (printlevel>2) then
    write(u_log, '(A,I3,1X,A,F14.9,A,F14.9,F14.9)') 'State', istate, ', rate1 (a.u.):     ', atraj%rate1,&
    &' rate2 (a.u.):     ', atraj%rate2, tmprate2
    ! The approximate formula seems to be smaller than the full formula. But the trends are similar.
  endif

  ! decoherence
  ! All states get the same random number. This should be ok ...
  if (traj%randnum2 < ctrl%dtstep * (atraj%rate1 - atraj%rate2)) then
    if (printlevel>2) then
      write(u_log, '(A,1X,I3)') 'Collapsing amplitude of state', istate
    endif
    cn = traj%coeff_diag_s(istate)
    clam = traj%coeff_diag_s(rstate)

    traj%coeff_diag_s(istate) = 0.d0
    traj%coeff_diag_s(rstate) = clam / abs(clam) * sqrt(abs(clam)**2 + abs(cn)**2)

    atraj%geom_ad = 0.d0
    atraj%veloc_ad = 0.d0
  endif

  ! reset
  if (traj%randnum2 < -ctrl%dtstep * atraj%rate1) then
    if (printlevel>2) then
      write(u_log, '(A,1X,I3)') 'Resetting moments for state', istate
    endif
    atraj%geom_ad = 0.d0
    atraj%veloc_ad = 0.d0
  endif

endsubroutine

! ===========================================================

!> performs the geometry update of the Velocity Verlet algorithm
!> a(t)=g(t)/M
!> x(t+dt)=x(t)+v(t)*dt+0.5*a(t)*dt^2
subroutine VelocityVerlet_xstep_afssh(atraj,ctrl,sigma)
  use definitions
  use matrix
  implicit none
  type(aux_trajectory_type) :: atraj
  type(ctrl_type) :: ctrl
  real*8 :: sigma
  integer :: iatom, idir

  do iatom=1,ctrl%natom
    do idir=1,3
      atraj%accel_ad(iatom,idir)=&
      &-atraj%grad_ad(iatom,idir)/atraj%mass_a(iatom)

      atraj%geom_tmp_ad(iatom,idir)=&
      & atraj%geom_ad(iatom,idir)&
      &+atraj%veloc_ad(iatom,idir)*ctrl%dtstep&
      &+0.5d0*atraj%accel_ad(iatom,idir)*ctrl%dtstep**2*sigma
    enddo
  enddo

  if (printlevel>2) then
    write(u_log,*) '============================================================='
    write(u_log,'(5X,A,1X,I3,1X,A)') 'Velocity Verlet (A-FSSH), state', atraj%istate, '-- X-step'
    write(u_log,*) '============================================================='

    write(u_log,'(2X,A,2E12.4)') 'Geom norms (old/new):', sum(atraj%geom_ad*atraj%geom_ad),&
    &sum(atraj%geom_tmp_ad*atraj%geom_tmp_ad)
    write(u_log,'(2X,A,2E12.4)') 'grad/accel norms:', sum(atraj%grad_ad*atraj%grad_ad),&
    &sum(atraj%accel_ad*atraj%accel_ad)

    if (printlevel>3) then
      call vec3write(ctrl%natom,atraj%geom_ad,u_log,'Old geom','F12.7')
      call vec3write(ctrl%natom,atraj%geom_tmp_ad,u_log,'geom','F12.7')
    endif
  endif

endsubroutine

! ===========================================================

!> performs the velocity update of the Velocity Verlet algorithm
!> a(t+dt)=g(t+dt)/M
!> v(t+dt)=v(t)+a(t+dt)*dt
subroutine VelocityVerlet_vstep_afssh(atraj,ctrl,sigma)
  use definitions
  use matrix
  implicit none
  type(aux_trajectory_type) :: atraj
  type(ctrl_type) :: ctrl
  real*8 :: sigma
  integer :: iatom, idir


  do iatom=1,ctrl%natom
    do idir=1,3
      atraj%accel_ad(iatom,idir)=0.5d0*(atraj%accel_ad(iatom,idir)&
      &-atraj%grad_ad(iatom,idir)/atraj%mass_a(iatom) )

      atraj%veloc_tmp_ad(iatom,idir)=&
      & atraj%veloc_ad(iatom,idir)&
      &+atraj%accel_ad(iatom,idir)*ctrl%dtstep*sigma
    enddo
  enddo

  if (printlevel>2) then
    write(u_log,*) '============================================================='
    write(u_log,'(5X,A,1X,I3,1X,A)') 'Velocity Verlet (A-FSSH), state', atraj%istate, '-- V-step'
    write(u_log,*) '============================================================='

    write(u_log,'(2X,A,2E12.4)') 'Veloc norms (old/new):', sum(atraj%veloc_ad*atraj%veloc_ad),&
    &sum(atraj%veloc_tmp_ad*atraj%veloc_tmp_ad)
    write(u_log,'(2X,A,2E12.4)') 'grad/accel norms:', sum(atraj%grad_ad*atraj%grad_ad),&
    &sum(atraj%accel_ad*atraj%accel_ad)

    if (printlevel>3) then
      call vec3write(ctrl%natom,atraj%veloc_ad,u_log,'Old veloc','F12.9')
      call vec3write(ctrl%natom,atraj%veloc_tmp_ad,u_log,'veloc','F12.9')
    endif
  endif

endsubroutine

! ===========================================================

endmodule
