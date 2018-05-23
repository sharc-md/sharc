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

!> # Module ELECTRONIC_LASER
!>
!> \author Sebastian Mai
!> 13.03.2015
!>
!> This module provides modified versions of some routines in electronic.f90, which
!> are used if a laser field enters in the wavefunction propagation.
!>
!>

module electronic_laser
  contains
! ==================================================================================================
! ==================================================================================================
! ==================================================================================================

!> Calculates the propagator Rtotal from the matrices in traj 
!> (H_MCH, H_MCH_old, NACdt, NACdt_old, U, U_old) and the timestep.
!> It also updates the diagonal and MCH coefficients.
subroutine propagate_laser(traj,ctrl)
  use definitions
  use matrix
  use electronic
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl
  integer :: istate, iatom, idir

  ! initialize the propagator matrix to the unit matrix
  traj%Rtotal_ss=dcmplx(0.d0,0.d0)
  do istate=1,ctrl%nstates
    traj%Rtotal_ss(istate,istate)=dcmplx(1.d0,0.d0)
  enddo

  ! call the appropriate propagator routine
  select case (ctrl%coupling)
    case (0)    ! ddt
      ! NADdt_ss can be directly used
      ! constant interpolation
      call unitary_propagator_laser(&
        &ctrl%nstates,&
        &traj%H_MCH_ss, traj%H_MCH_old_ss,&
        &traj%NACdt_ss, traj%NACdt_old_ss,&
        &traj%U_ss,traj%U_old_ss,&
        &traj%DM_ssd,traj%DM_old_ssd,&
        &ctrl%laserfield_td( (traj%step-1)*ctrl%nsubsteps+2:traj%step*ctrl%nsubsteps+1 ,:),&
        &ctrl%dtstep, ctrl%nsubsteps, 1,&       ! 1=constant interpolation
        &traj%Rtotal_ss)
    case (1)    ! ddr
      ! NACdr_ssad has to be scalar multiplied with velocity, NACdt_old_ss already contains the old scalar products
      ! linear interpolation
      traj%NACdt_ss=dcmplx(0.d0,0.d0)
      do iatom=1,ctrl%natom
        do idir=1,3
          traj%NACdt_ss=traj%NACdt_ss+traj%NACdr_ssad(:,:,iatom,idir)*traj%veloc_ad(iatom,idir)
        enddo
      enddo
      call unitary_propagator_laser(&
        &ctrl%nstates,&
        &traj%H_MCH_ss, traj%H_MCH_old_ss,&
        &traj%NACdt_ss, traj%NACdt_old_ss,&
        &traj%U_ss,traj%U_old_ss,&
        &traj%DM_ssd,traj%DM_old_ssd,&
        &ctrl%laserfield_td( (traj%step-1)*ctrl%nsubsteps+2:traj%step*ctrl%nsubsteps+1 ,:),&
        &ctrl%dtstep, ctrl%nsubsteps, 0,&       ! 0=linear interpolation
        &traj%Rtotal_ss)
    case (2)    ! overlap
      ! overlap matrix is used
      ! use LOCAL DIABATISATION
      call LD_propagator_laser(&
        &ctrl%nstates,&
        &traj%H_MCH_ss, traj%H_MCH_old_ss,&
        &traj%U_ss,traj%U_old_ss,&
        &traj%overlaps_ss,&
        &traj%DM_ssd,traj%DM_old_ssd,&
        &ctrl%laserfield_td( (traj%step-1)*ctrl%nsubsteps+2:traj%step*ctrl%nsubsteps+1 ,:),&
        &ctrl%dtstep, ctrl%nsubsteps,&
        &traj%Rtotal_ss)
  endselect


  if (printlevel>2) then
    write(u_log,*) '============================================================='
    write(u_log,*) '            Propagating the electronic wavefunction'
    write(u_log,*) '============================================================='
    if (printlevel>3) then
      call vec3write(ctrl%nsubsteps,&
      &ctrl%laserfield_td( (traj%step-1)*ctrl%nsubsteps+2:traj%step*ctrl%nsubsteps+1 ,:),&
      &u_log,'Laser Field','F12.9')
    endif
    select case (ctrl%coupling)
      case (0)  ! ddt
        write(u_log,*) 'Propagating the coefficients using the ddt matrix...'
        if (printlevel>4) then
          call matwrite(ctrl%nstates,traj%H_MCH_old_ss,u_log,'Old H_MCH Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%H_MCH_ss,u_log,'H_MCH Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%NACdt_old_ss,u_log,'Old DDT Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%NACdt_ss,u_log,'DDT Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%U_old_ss,u_log,'U_old Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%U_ss,u_log,'U Matrix','F12.9')
        endif
      case (1)  ! ddr
        write(u_log,*) 'Propagating the coefficients using the ddr vectors...'
        if (printlevel>4) then
          call matwrite(ctrl%nstates,traj%H_MCH_old_ss,u_log,'Old H_MCH Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%H_MCH_ss,u_log,'H_MCH Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%NACdt_old_ss,u_log,'Old DDT Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%NACdt_ss,u_log,'DDT Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%U_old_ss,u_log,'U_old Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%U_ss,u_log,'U Matrix','F12.9')
        endif
      case (2)  ! overlap
        write(u_log,*) 'Propagating the coefficients using Local Diabatisation...'
        if (printlevel>4) then
          call matwrite(ctrl%nstates,traj%H_MCH_old_ss,u_log,'Old H_MCH Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%H_MCH_ss,u_log,'H_MCH Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%overlaps_ss,u_log,'Overlap Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%U_old_ss,u_log,'U_old Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%U_ss,u_log,'U Matrix','F12.9')
        endif
    endselect
    if (printlevel>3) then
      call matwrite(ctrl%nstates,traj%Rtotal_ss,u_log,'Propagator Matrix','F12.9')
    endif
  endif

  ! check for NaNs in the Propagator matrix
  if (any((real(traj%Rtotal_ss)).ne.(real(traj%Rtotal_ss))).or.any((aimag(traj%Rtotal_ss)).ne.(aimag(traj%Rtotal_ss)))) then
    write(0,*) 'The propagator matrix contains NaNs!'
    select case (ctrl%coupling)
      case (0)  ! ddt
          call matwrite(ctrl%nstates,traj%H_MCH_old_ss,0,'Old H_MCH Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%H_MCH_ss,0,'H_MCH Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%NACdt_old_ss,0,'Old DDT Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%NACdt_ss,0,'DDT Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%U_old_ss,0,'U_old Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%U_ss,0,'U Matrix','F12.9')
      case (1)  ! ddr
          call matwrite(ctrl%nstates,traj%H_MCH_old_ss,0,'Old H_MCH Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%H_MCH_ss,0,'H_MCH Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%NACdt_old_ss,0,'Old DDT Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%NACdt_ss,0,'DDT Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%U_old_ss,0,'U_old Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%U_ss,0,'U Matrix','F12.9')
      case (2)  ! overlap
          call matwrite(ctrl%nstates,traj%H_MCH_old_ss,0,'Old H_MCH Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%H_MCH_ss,0,'H_MCH Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%overlaps_ss,0,'Overlap Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%U_old_ss,0,'U_old Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%U_ss,0,'U Matrix','F12.9')
    endselect
    call matwrite(ctrl%nstates,traj%Rtotal_ss,0,'Propagator Matrix','F12.9')
  endif

  ! propagate the coefficients
  traj%coeff_diag_old_s=traj%coeff_diag_s
  call matvecmultiply(&
    &ctrl%nstates,&
    &traj%Rtotal_ss, traj%coeff_diag_old_s, traj%coeff_diag_s, &
    &'n')

  ! get coefficients in MCH picture
  call matvecmultiply(ctrl%nstates,traj%U_ss,traj%coeff_diag_s,traj%coeff_MCH_s,'n')
  traj%state_MCH=state_diag_to_MCH(ctrl%nstates,traj%state_diag,traj%U_ss)

  if (printlevel>2) then
    write(u_log,*) 'Old and new coefficients:'
    do istate=1,ctrl%nstates
      write(u_log,'(2(F7.4,1X),4X,2(F7.4,1X))') traj%coeff_diag_old_s(istate),traj%coeff_diag_s(istate)
    enddo
  endif

endsubroutine

! ==================================================================================================
! ==================================================================================================
! ==================================================================================================

!> calculates the propagator matrix in substeps, see SHARC manual for the equations.
!> \param interp 0=linear interpolation of non-adiabatic coupling matrix, 1=constant non-adiabatic coupling matrix (NACMold not used)
subroutine unitary_propagator_laser(n, SO, SOold, NACM, NACMold, U, Uold, DM, DMold, laserfield, dt, nsubsteps, interp, Rtotal)
  use definitions, only: u_log
  use matrix
! calculates the propagator matrix for a timestep
! it calculates:
!              n
!  R = U^t . PROD exp( -[iH+T]*dt ) . Uold
!             i=1
! note that Rtotal has to be initialized as a unit matrix prior to calling unitary_propagator()
!
! interp: 0 linear interpolation of T, 1 constant interpolation of T
  implicit none

  integer, intent(in) :: n, nsubsteps, interp
  complex*16, intent(in) :: U(n,n), SO(n,n), SOold(n,n), NACM(n,n), NACMold(n,n), DM(n,n,3),DMold(n,n,3)
  complex*16, intent(in) :: laserfield(nsubsteps,3)
  complex*16, intent(inout) :: Uold(n,n), Rtotal(n,n)
  real*8, intent(in) :: dt
  ! internal variables:
  integer :: istep, ixyz
  real*8 :: dtsubstep
  complex*16 :: H(n,n), T(n,n)
  complex*16 :: Rexp(n,n), Rprod(n,n)
  complex*16 :: ii=dcmplx(0.d0,1.d0)


  dtsubstep=dt/nsubsteps

  call matmultiply(n,Uold,Rtotal,Rprod,'nn')
  Rtotal=Rprod

  do istep=1,nsubsteps

    ! first ingredient, H
    H=SOold + (SO-SOold)*istep/nsubsteps
    ! here the laser field is added to the Hamiltonian
    do ixyz=1,3
      H=H - ( DMold(:,:,ixyz) + (DM(:,:,ixyz)-DMold(:,:,ixyz))*istep/nsubsteps ) * real(laserfield(istep,ixyz))
    enddo

    ! second ingredient, T
    if (interp==1) then
      T=NACM
    else
      T=NACMold+ (NACM-NACMold)*istep/nsubsteps
    endif

    ! set up the total operator: (iUHU+UTU+UdU/dt)*dtsubstep
    Rexp=dtsubstep*(H-ii*T)

    ! calculate the Operator Exponential
    call exponentiate(n,Rexp,-ii)

    ! add the propagator of the current timesubstep to the total propagator of the full timestep
    call matmultiply(n,Rexp,Rtotal,Rprod,'nn')
    Rtotal=Rprod

  enddo

  call matmultiply(n,U,Rtotal,Rprod,'tn')
  Rtotal=Rprod

  return

endsubroutine

! ==================================================================================================
! ==================================================================================================
! ==================================================================================================

!> calculates the propagator matrix in substeps, see SHARC manual for the equations.
!> Uses the local diabatization procedure
subroutine LD_propagator_laser(n, SOin, SOold, U, Uold, overlap, DMin, DMold, laserfield, dt, nsubsteps, Rtotal)
  use definitions, only: u_log
  use matrix
  ! calculates the propagator matrix for a timestep
  ! it calculates:
  !                    n
  !  R = U^t . S^t . PROD exp( -[H]*dt ) . Uold
  !                   i=1
  ! note that Rtotal has to be initialized as a unit matrix prior to calling LD_propagator()
  !
  ! H is interpolated linearly
  implicit none
  integer, intent(in) :: n, nsubsteps
  real*8, intent(in) :: dt
  complex*16, intent(in) :: U(n,n), Uold(n,n),SOin(n,n),SOold(n,n), DMin(n,n,3),DMold(n,n,3)
  complex*16, intent(in) :: laserfield(nsubsteps,3)
  complex*16, intent(inout) :: overlap(n,n)
  complex*16, intent(inout) :: Rtotal(n,n)

  integer :: i,j,k, ixyz
  complex*16 :: H(n,n), SO(n,n), DM(n,n,3), Rprod(n,n)
  real*8 :: sums, dtsubstep

  real*8,parameter :: intr_thrs=1.d-1
  complex*16,parameter :: ii=dcmplx(0.d0,1.d0)

  ! Intruder state check
  do i=1,n
    sums=0.d0
    do j=1,n
      sums=sums+abs(overlap(i,j))**2
      sums=sums+abs(overlap(j,i))**2
    enddo
    sums=sums-abs(overlap(i,i))**2

    if (sums < intr_thrs) then
      write(u_log,'(A)') '! ======== INTRUDER STATE PROBLEM ======== !'
      write(u_log,'(A,I4)') 'State: ',i
      do k=1,n
        write(u_log,'(1000(F8.5,1X))') (overlap(k,j),j=1,n)
      enddo

      overlap(i,:)=dcmplx(0.d0,0.d0)
      overlap(:,i)=dcmplx(0.d0,0.d0)
      overlap(i,i)=dcmplx(1.d0,0.d0)
    endif
  enddo

  ! Löwdin orthogonalisation
  call lowdin(n,overlap)

  ! Initialize Rtotal
  call matmultiply(n,Uold,Rtotal,Rprod,'nn')
  Rtotal=Rprod

  ! Transform SO into old basis
  SO=SOin
  call transform(n,SO,overlap,'uaut')
  DM=DMin
  do ixyz=1,3
    call transform(n,DM(:,:,ixyz),overlap,'uaut')
  enddo

  ! Evolve in the diabatic basis in substeps
  dtsubstep=dt/nsubsteps
  do k=1,nsubsteps
    H=SOold + (SO-SOold)*k/nsubsteps
    ! here the laser field is added to the Hamiltonian
    do ixyz=1,3
      H=H - ( DMold(:,:,ixyz) + (DM(:,:,ixyz)-DMold(:,:,ixyz))*k/nsubsteps ) * real(laserfield(k,ixyz))
    enddo
    H=dtsubstep*H

    call exponentiate(n,H,-ii)

    call matmultiply(n,H,Rtotal,Rprod,'nn')
    Rtotal=Rprod
  enddo

  ! Finalize Rtotal
  call matmultiply(n,overlap,Rtotal,Rprod,'tn')
  call matmultiply(n,U,Rprod,Rtotal,'tn')

endsubroutine

! ==================================================================================================
! ==================================================================================================
! ==================================================================================================

!> template for a subroutine returning the laser field for a given time
!> this is not yet fully implemented
subroutine internal_laserfield(t,field,energy)
implicit none
real*8,intent(in) :: t     ! time in atomic units
complex*16,intent(out) :: field(3)      ! x,y,z of laser field, imaginary part has to be included
complex*16,intent(out) :: energy        ! momentary energy of laser in hartree

field=dcmplx(0.d0,0.d0)
energy=dcmplx(t-t,0.d0)

return

endsubroutine




endmodule