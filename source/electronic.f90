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

!> # Module ELECTRONIC
!>
!> \author Sebastian Mai
!> \date 27.02.2015
!>
!>                   modified 11.13.2019 by Yinan Shu
!>                        add the following subroutines:
!>                        surface_switching, calculate_probabilities_CSDM,
!>                        calculate_probabilities_NDM
!>                        calculate_probabilitiesDoM - similar as calculate_probabilities
!>  
!> This module provides all routines for the propagation of the electronic wavefunction,
!> transformation between representations, decoherence and surface hopping (except for momentum adjustment)
!> 
!> The alternative module electronic_laser.f90 reimplements some of the subroutines 
!> to add a laser field

module electronic
  contains
! ==================================================================================================
! ==================================================================================================
! ==================================================================================================

!> Calculates the propagator Rtotal from the matrices in traj 
!> (H_MCH, H_MCH_old, NACdt, NACdt_old, U, U_old) and the timestep.
!> It also updates the diagonal and MCH coefficients.
subroutine propagate(traj,ctrl)
  use definitions
  use matrix
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
  select case (ctrl%eeom)
    case (0)    ! constant interpolation, default for coupling=ddt
      ! NADdt_ss can be directly used
      ! constant interpolation
      call unitary_propagator(&
        &ctrl%nstates,&
        &traj%H_MCH_ss, traj%H_MCH_old_ss,&
        &traj%NACdt_ss, traj%NACdt_old_ss,&
        &traj%U_ss,traj%U_old_ss,&
        &ctrl%dtstep, ctrl%nsubsteps, 1,&       ! 1=constant interpolation
        &traj%Rtotal_ss)
    case (1)    ! linear interpolation, default for coupling=ddr,nacdr
      ! NACdr_ssad has to be scalar multiplied with velocity, NACdt_old_ss already contains the old scalar products
      ! linear interpolation
      call unitary_propagator(&
        &ctrl%nstates,&
        &traj%H_MCH_ss, traj%H_MCH_old_ss,&
        &traj%NACdt_ss, traj%NACdt_old_ss,&
        &traj%U_ss,traj%U_old_ss,&
        &ctrl%dtstep, ctrl%nsubsteps, 0,&       ! 0=linear interpolation
        &traj%Rtotal_ss)
    case (2)    ! local diabatization, defafult for coupling=overlap
      ! overlap matrix is used
      ! use LOCAL DIABATISATION
      call LD_propagator(&
        &ctrl%nstates,&
        &traj%H_MCH_ss, traj%H_MCH_old_ss,&
        &traj%U_ss,traj%U_old_ss,&
        &traj%overlaps_ss,&
        &ctrl%dtstep, ctrl%nsubsteps,&
        &traj%Rtotal_ss)
    case (3)    ! norm perserving interporlation
      call NPI_propagator(&
        &ctrl%nstates,&
        &traj%H_MCH_ss, traj%H_MCH_old_ss,&
        &traj%U_ss,traj%U_old_ss,&
        &traj%overlaps_ss,&
        &ctrl%dtstep, ctrl%nsubsteps,&   
        &traj%Rtotal_ss)
  endselect

  if (printlevel>2) then
    write(u_log,*) '============================================================='
    write(u_log,*) '            Propagating the electronic wavefunction'
    write(u_log,*) '============================================================='
    select case (ctrl%eeom)
      case (0)  ! constant interpolation 
        write(u_log,*) 'Propagating the coefficients using the constant interpolated ddt matrix...'
        if (printlevel>4) then
          call matwrite(ctrl%nstates,traj%H_MCH_old_ss,u_log,'Old H_MCH Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%H_MCH_ss,u_log,'H_MCH Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%NACdt_old_ss,u_log,'Old DDT Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%NACdt_ss,u_log,'DDT Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%U_old_ss,u_log,'U_old Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%U_ss,u_log,'U Matrix','F12.9')
        endif
      case (1)  ! linear interpolation 
        write(u_log,*) 'Propagating the coefficients using the linearly interpolated ddt matrix...'
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
      case (3)  ! npi
        write(u_log,*) 'Propagating the coefficients using Norm Perserving Interpolation...'
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
    select case (ctrl%eeom)
      case (0)  ! constant interpolation 
          call matwrite(ctrl%nstates,traj%H_MCH_old_ss,0,'Old H_MCH Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%H_MCH_ss,0,'H_MCH Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%NACdt_old_ss,0,'Old DDT Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%NACdt_ss,0,'DDT Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%U_old_ss,0,'U_old Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%U_ss,0,'U Matrix','F12.9')
      case (1)  ! linear interpolation
          call matwrite(ctrl%nstates,traj%H_MCH_old_ss,0,'Old H_MCH Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%H_MCH_ss,0,'H_MCH Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%NACdt_old_ss,0,'Old DDT Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%NACdt_ss,0,'DDT Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%U_old_ss,0,'U_old Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%U_ss,0,'U Matrix','F12.9')
      case (2)  ! local diabatization 
          call matwrite(ctrl%nstates,traj%H_MCH_old_ss,0,'Old H_MCH Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%H_MCH_ss,0,'H_MCH Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%overlaps_ss,0,'Overlap Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%U_old_ss,0,'U_old Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%U_ss,0,'U Matrix','F12.9')
      case (3)  ! npi
          call matwrite(ctrl%nstates,traj%H_MCH_old_ss,0,'Old H_MCH Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%H_MCH_ss,0,'H_MCH Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%overlaps_ss,0,'Overlap Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%U_old_ss,0,'U_old Matrix','F12.9')
          call matwrite(ctrl%nstates,traj%U_ss,0,'U Matrix','F12.9')
    endselect
    call matwrite(ctrl%nstates,traj%Rtotal_ss,0,'Propagator Matrix','F12.9')
    stop 1
  endif

  ! save old coefficients
  traj%coeff_diag_old_s=traj%coeff_diag_s
  traj%coeff_MCH_old_s=traj%coeff_MCH_s

  ! propagate the coefficients
  call matvecmultiply(&
    &ctrl%nstates,&
    &traj%Rtotal_ss, traj%coeff_diag_old_s, traj%coeff_diag_s, &
    &'n')

  ! get coefficients in MCH picture
  call matvecmultiply(ctrl%nstates,traj%U_ss,traj%coeff_diag_s,traj%coeff_MCH_s,'n')
  traj%state_MCH=state_diag_to_MCH(ctrl%nstates,traj%state_diag,traj%U_ss)

  if (printlevel>2) then
    write(u_log,*) 'Old and new diagonal coefficients:'
    do istate=1,ctrl%nstates
      write(u_log,'(2(F7.4,1X),4X,2(F7.4,1X))') traj%coeff_diag_old_s(istate),traj%coeff_diag_s(istate)
    enddo
    write(u_log,*) 'Old and new MCH coefficients:'
    do istate=1,ctrl%nstates
      write(u_log,'(2(F7.4,1X),4X,2(F7.4,1X))') traj%coeff_mch_old_s(istate),traj%coeff_mch_s(istate)
    enddo
  endif

endsubroutine

! ==================================================================================================
! ==================================================================================================
! ==================================================================================================

!> calculates the propagator matrix in substeps, see SHARC manual for the equations.
!> \param interp 0=linear interpolation of non-adiabatic coupling matrix, 1=constant non-adiabatic coupling matrix (NACMold not used)
subroutine unitary_propagator(n, SO, SOold, NACM, NACMold, U, Uold, dt, nsubsteps, interp, Rtotal)
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
complex*16, intent(in) :: U(n,n), SO(n,n), SOold(n,n), NACM(n,n), NACMold(n,n)
complex*16, intent(inout) :: Uold(n,n), Rtotal(n,n)
real*8, intent(in) :: dt
! internal variables:
integer :: istep
real*8 :: dtsubstep
complex*16 :: H(n,n), T(n,n)
complex*16 :: Rexp(n,n), Rprod(n,n)
complex*16 :: ii=dcmplx(0.d0,1.d0)

dtsubstep=dt/nsubsteps

call matmultiply(n,Uold,Rtotal,Rprod,'nn')
Rtotal=Rprod

do istep=1,nsubsteps

  ! first ingredient, H
  H=SOold+ (SO-SOold)*istep/nsubsteps

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
subroutine LD_propagator(n, SOin, SOold, U, Uold, overlap, dt, nsubsteps, Rtotal)
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
complex*16, intent(in) :: U(n,n), Uold(n,n),SOin(n,n),SOold(n,n)
complex*16, intent(inout) :: overlap(n,n)
complex*16, intent(inout) :: Rtotal(n,n)

integer :: i,j,k
complex*16 :: H(n,n), SO(n,n), Rprod(n,n)
real*8 :: sums, dtsubstep

real*8,parameter :: intr_thrs=1.d-1
complex*16,parameter :: ii=dcmplx(0.d0,1.d0)

! ! Intruder state check
! do i=1,n
!   sums=0.d0
!   do j=1,n
!     sums=sums+abs(overlap(i,j))**2
!     sums=sums+abs(overlap(j,i))**2
!   enddo
!   sums=sums-abs(overlap(i,i))**2
! 
!   if (sums < intr_thrs) then
!     write(u_log,'(A)') '! ======== INTRUDER STATE PROBLEM ======== !'
!     write(u_log,'(A,I4)') 'State: ',i
!     do k=1,n
!       write(u_log,'(1000(F8.5,1X))') (overlap(k,j),j=1,n)
!     enddo
! 
!     overlap(i,:)=dcmplx(0.d0,0.d0)
!     overlap(:,i)=dcmplx(0.d0,0.d0)
!     overlap(i,i)=dcmplx(1.d0,0.d0)
!   endif
! enddo
! 
! ! Löwdin orthogonalisation
! call lowdin(n,overlap)

! Initialize Rtotal
call matmultiply(n,Uold,Rtotal,Rprod,'nn')
Rtotal=Rprod

! Transform SO into old basis
SO=SOin
call transform(n,SO,overlap,'uaut')

! Evolve in the diabatic basis in substeps
dtsubstep=dt/nsubsteps
do k=1,nsubsteps
  H=SOold + (SO-SOold)*k/nsubsteps
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

!> calculates the propagator matrix in substeps, see SHARC manual for the equations.
!> Uses the norm perserving interpolation. 
subroutine NPI_propagator(n, SO, SOold, U, Uold, overlap, dt, nsubsteps, Rtotal)
use definitions, only: u_log
use matrix
! calculates the propagator matrix for a timestep
! it calculates:
!                    n
!  R = U^t . S^t . PROD exp( -[iH+T]*dt ) . Uold
!                   i=1
! note that Rtotal has to be initialized as a unit matrix prior to calling NPI_propagator()
!
! T is norm perserving interpolated

implicit none
integer, intent(in) :: n, nsubsteps
real*8, intent(in) :: dt
complex*16, intent(in) :: U(n,n), Uold(n,n),SO(n,n),SOold(n,n)
complex*16, intent(inout) :: overlap(n,n)
complex*16, intent(inout) :: Rtotal(n,n)

! internal variables:
integer :: istep, istate, jstate
real*8 :: dtsubstep
complex*16 :: H(n,n), T(n,n)
complex*16 :: Rexp(n,n), Rprod(n,n)
complex*16 :: ii=dcmplx(0.d0,1.d0)
complex*16 :: w(n,n),tw(n,n),dw(n,n)

dtsubstep=dt/nsubsteps

! Initialize Rtotal
call matmultiply(n,Uold,Rtotal,Rprod,'nn')
Rtotal=Rprod

!initialize T
T=dcmplx(0.d0,0.d0)

do istep=1,nsubsteps

  H=SOold+ (SO-SOold)*istep/nsubsteps
  ! compute NPI rotation matrix W
  do istate=1,n
    do jstate=1,n 
      if (jstate .eq. istate) then
        w(istate,jstate)=cos(acos(overlap(istate,jstate))*istep/nsubsteps)
        tw(jstate,istate)=cos(acos(overlap(istate,jstate))*istep/nsubsteps)
        dw(istate,jstate)=-sin(acos(overlap(istate,jstate))*istep/nsubsteps)*acos(overlap(istate,jstate))/dt
      else
        w(istate,jstate)=sin(asin(overlap(istate,jstate))*istep/nsubsteps)
        tw(jstate,istate)=sin(asin(overlap(istate,jstate))*istep/nsubsteps)
        dw(istate,jstate)=cos(asin(overlap(istate,jstate))*istep/nsubsteps)*asin(overlap(istate,jstate))/dt
      endif
    enddo
  enddo 
  ! compute intermediate T 
  call matmultiply(n, tw, dw, T, 'nn')
  do istate=1,n
    T(istate,istate)=dcmplx(0.d0,0.d0)
  enddo

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

! ===========================================================

!> this function finds the most similar diagonal state to a given MCH state.
  integer function state_MCH_to_diag(n,state_MCH,U)
    use matrix
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: state_MCH
    complex*16, intent(in) :: U(n,n)
    complex*16 :: v(n), v2(n)
    real*8 :: maxv,currv
    integer :: i,imax

    v=dcmplx(0.d0,0.d0)
    v(state_MCH)=dcmplx(1.d0,0.d0)
    call matvecmultiply(n,U,v,v2,'t')

    imax=1
    maxv=0.d0
    do i=1,n
      currv=abs(v2(i))**2
      if (maxv<currv) then
        maxv=currv
        imax=i
      endif
    enddo

    state_MCH_to_diag=imax
    return
  endfunction

! ===========================================================

!> this function finds the most similar MCH state to a given diagonal state.
  integer function state_diag_to_MCH(n,state_diag,U)
    use matrix
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: state_diag
    complex*16, intent(in) :: U(n,n)
    complex*16 :: b(n), b2(n)
    real*8 :: maxv,currv
    integer :: i,imax

    b=dcmplx(0.d0,0.d0)
    b(state_diag)=dcmplx(1.d0,0.d0)
    call matvecmultiply(n,U,b,b2,'n')

    imax=1
    maxv=0.d0
    do i=1,n
      currv=abs(b2(i))**2
      if (maxv<currv) then
        maxv=currv
        imax=i
      endif
    enddo

    state_diag_to_MCH=imax
    return
  endfunction

! ===========================================================================
! Electronic_Gradients, used by Bulirsch_Stoer_Hack integrator
! Computes the derivatives of electronic coeffcients at current step
! dc/dt=[iH-vK]c
subroutine Electronic_gradients_MCH(traj,ctrl)
  use definitions
  use matrix
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

  complex*16 :: H_ss(ctrl%nstates,ctrl%nstates)
  complex*16 :: w(ctrl%nstates,ctrl%nstates),tw(ctrl%nstates,ctrl%nstates),dw(ctrl%nstates,ctrl%nstates)
  complex*16 :: NACT(ctrl%nstates, ctrl%nstates)
  complex*16 :: Ptotal(ctrl%nstates, ctrl%nstates)
  integer :: istate, jstate, iatom, idir
  real*8 :: pvib_ad(ctrl%natom,3)
  complex*16 :: DP(ctrl%nstates, ctrl%nstates)
  complex*16 :: gdcoeff_MCH_s(ctrl%nstates)
  real*8 :: tmp, sntmp, cstmp, tmpr, tmpi

  if (printlevel>2) then
    write(u_log,*) '============================================================='
    write(u_log,*) '           Calculating Electronic Gradients'
    write(u_log,*) '============================================================='
  endif

  ! include laser fields
  if (ctrl%laser==0) then
    H_ss=traj%H_MCH_ss
  else if (ctrl%laser==2) then
    do idir=1,3
      H_ss=traj%H_MCH_ss-traj%DM_ssd(:,:,idir)*real(ctrl%laserfield_td(traj%step*ctrl%nsubsteps+1,idir))
    enddo
  endif
  if (printlevel>4) then
    call matwrite(ctrl%nstates,H_ss,u_log,' H_MCH with laser field','F14.9')
  endif

  ! Evaluate time-derivative of electronic phase angle
  do istate=1,ctrl%nstates
    traj%gephase_s(istate)=H_ss(istate,istate)
  enddo

  ! Evaluate time-derivative coupling term vK. 
  NACT=dcmplx(0.d0,0.d0)
  if (ctrl%coupling==0) then
    NACT=traj%NACdt_ss
  elseif (ctrl%coupling==1) then  
    do iatom=1,ctrl%natom
      do idir=1,3
        NACT=NACT+traj%NACdr_ssad(:,:,iatom,idir)*traj%veloc_ad(iatom,idir)
      enddo
    enddo
  elseif (ctrl%coupling==2) then
  ! use NPI to evaluate NACT
    if (traj%step==0) then
      do iatom=1,ctrl%natom
        do idir=1,3
          NACT=NACT+traj%NACdr_ssad(:,:,iatom,idir)*traj%veloc_ad(iatom,idir)
        enddo
      enddo 
    elseif (traj%step>=1) then
       do istate=1,ctrl%nstates
         do jstate=1,ctrl%nstates
           if (jstate .eq. istate) then
             w(istate,jstate)=cos(acos(traj%overlaps_ss(istate,jstate)))
             tw(jstate,istate)=cos(acos(traj%overlaps_ss(istate,jstate)))
             dw(istate,jstate)=-sin(acos(traj%overlaps_ss(istate,jstate)))*acos(traj%overlaps_ss(istate,jstate))/ctrl%dtstep
           else
             w(istate,jstate)=sin(asin(traj%overlaps_ss(istate,jstate)))
             tw(jstate,istate)=sin(asin(traj%overlaps_ss(istate,jstate)))
             dw(istate,jstate)=cos(asin(traj%overlaps_ss(istate,jstate)))*asin(traj%overlaps_ss(istate,jstate))/ctrl%dtstep
           endif
         enddo
       enddo
       call matmultiply(ctrl%nstates, tw, dw, NACT, 'nn')
       do istate=1,ctrl%nstates
         NACT(istate,istate)=dcmplx(0.d0,0.d0)
       enddo
    endif
  elseif (ctrl%coupling==3) then
    NACT=traj%NACdt_ss
  endif
  
  if (printlevel>4) then
    call matwrite(ctrl%nstates,NACT,u_log,'time derivative coupling','F14.9')
  endif   

  traj%gRcoeff_MCH_s=0.d0
  traj%gIcoeff_MCH_s=0.d0
  traj%gRccoeff_MCH_s=0.d0
  traj%gIccoeff_MCH_s=0.d0
 
  ! compute coefficients gradient, traj%gRcoeff_MCH_s and traj%gIcoeff_MCH_s
  do istate=1,ctrl%nstates
    do jstate=1,ctrl%nstates
      if (istate.ne.jstate) then
        tmp = traj%ephase_s(jstate)-traj%ephase_s(istate)
        sntmp = sin(tmp)
        cstmp = cos(tmp)
        tmpr = -(cstmp*real(traj%coeff_MCH_s(jstate))+sntmp*aimag(traj%coeff_MCH_s(jstate)))*NACT(istate,jstate)
        tmpi = -(cstmp*aimag(traj%coeff_MCH_s(jstate))-sntmp*real(traj%coeff_MCH_s(jstate)))*NACT(istate,jstate)
        traj%gRcoeff_MCH_s(istate)=traj%gRcoeff_MCH_s(istate)+tmpr
        traj%gIcoeff_MCH_s(istate)=traj%gIcoeff_MCH_s(istate)+tmpi
      endif
    enddo
  enddo

 ! compute coherent coefficients gradient, traj%gRccoeff_MCH_s and traj%gIccoeff_MCH_s
  do istate=1,ctrl%nstates
    do jstate=1,ctrl%nstates
      if (istate.ne.jstate) then
        tmp = traj%ephase_s(jstate)-traj%ephase_s(istate)
        sntmp = sin(tmp)
        cstmp = cos(tmp)
        tmpr = -(cstmp*real(traj%ccoeff_MCH_s(jstate))+sntmp*aimag(traj%ccoeff_MCH_s(jstate)))*NACT(istate,jstate)
        tmpi = -(cstmp*aimag(traj%ccoeff_MCH_s(jstate))-sntmp*real(traj%ccoeff_MCH_s(jstate)))*NACT(istate,jstate)
        traj%gRccoeff_MCH_s(istate)=traj%gRccoeff_MCH_s(istate)+tmpr
        traj%gIccoeff_MCH_s(istate)=traj%gIccoeff_MCH_s(istate)+tmpi
      endif
    enddo
  enddo

  if (printlevel>4) then
    call vecwrite(ctrl%nstates, traj%gRcoeff_MCH_s, u_log, 'Real coefficient Re(coeff_mch_s) gradient','F14.9')
    call vecwrite(ctrl%nstates, traj%gIcoeff_MCH_s, u_log, 'Imaginary coefficient Im(coeff_mch_s) gradient','F14.9')
    call vecwrite(ctrl%nstates, traj%gRccoeff_MCH_s, u_log, 'Real coherent coefficient Re(coeff_mch_s) gradient','F14.9')
    call vecwrite(ctrl%nstates, traj%gIccoeff_MCH_s, u_log, 'Imaginary coherent coefficient Im(coeff_mch_s) gradient','F14.9')
  endif 

  ! compute pre-switching probabilities, b_ij
  if (ctrl%method==0) then
    call calculate_probabilities_preBSH(ctrl%natom, ctrl%nstates, traj%veloc_ad, traj%coeff_diag_s, traj%NACdR_diag_ssad, traj%state_diag, traj%gpreprob_s3)
  elseif (ctrl%method==1) then
    if (ctrl%decoherence==0) then
      call calculate_probabilities_preBSH(ctrl%natom, ctrl%nstates, traj%veloc_ad, traj%coeff_diag_s, traj%NACdR_diag_ssad, traj%state_diag, traj%gpreprob_s3)
    elseif (ctrl%decoherence==11) then ! decay of mixing
      if (ctrl%switching_procedure==1 .or. ctrl%switching_procedure==2) then ! CSDM and SCDM, use ccoeff
        call calculate_probabilities_preBSH(ctrl%natom, ctrl%nstates, traj%veloc_ad, traj%ccoeff_diag_s, traj%NACdR_diag_ssad, traj%state_diag, traj%gpreprob_s3)
      elseif (ctrl%switching_procedure==2) then ! NDM, use coeff
        call calculate_probabilities_preBSH(ctrl%natom, ctrl%nstates, traj%veloc_ad, traj%coeff_diag_s, traj%NACdR_diag_ssad, traj%state_diag, traj%gpreprob_s3)
      endif
    endif
  endif

  if (printlevel>4) then
    call vec3write(ctrl%nstates, traj%gpreprob_s3, u_log, 'b function, gradient of switching probability','F14.9')
  endif

endsubroutine

! ===========================================================================
! BSH integrator integrates the MCH coefficients
! Here we convert all the electronic coefficients integrated after BSH step to
! DIAG representation
subroutine Calculate_cDIAG(traj,ctrl)
  use definitions
  use matrix
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

  call matvecmultiply(ctrl%nstates,traj%U_ss,traj%coeff_MCH_s,traj%coeff_diag_s,'tn')
  call matvecmultiply(ctrl%nstates,traj%U_ss,traj%ccoeff_MCH_s,traj%ccoeff_diag_s,'tn')

  if (printlevel>4) then
    call vecwrite(ctrl%nstates, traj%coeff_MCH_s, u_log, 'coefficients in MCH basis','F14.9')
    call vecwrite(ctrl%nstates, traj%coeff_diag_s, u_log, 'coefficients in diag basis','F14.9')
    call vecwrite(ctrl%nstates, traj%ccoeff_MCH_s, u_log, 'coherent coefficients in MCH basis','F14.9')
    call vecwrite(ctrl%nstates, traj%ccoeff_diag_s, u_log, 'coherent coefficients in diag basis','F14.9')
  endif

endsubroutine

! ==================================================================================================
! ==================================================================================================
! ===                                        Surface Hopping                                    ====
! ==================================================================================================
! ==================================================================================================

!> this routine performs all steps of the surface hopping algorithm:
!> - calculate the surface hopping probabilities (update traj%hopprob_s)
!> - generate a random number
!> - find the new active state
!> - check for resonance with laser
!> - check for frustrated hops
subroutine surface_hopping(traj,ctrl)
  use definitions
  use matrix
  use nuclear
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl
  integer :: istate,jstate,kstate,iatom,idir,i,ilaser
  real*8 :: randnum, cumuprob
  real*8 :: Emax ! energy threshold for frustrated jumps
  real*8 :: Ekin_masked
  real*8 :: sum_kk, sum_vk ! temp variables for kinetic adjustment
  real*8 :: deltaE

  if (printlevel>2) then
    write(u_log,*) '============================================================='
    write(u_log,*) '                 Doing fewest switches SH'
    write(u_log,*) '============================================================='
  endif

  ! calculate the hopping probabilities
  if (ctrl%integrator==0) then
    call calculate_probabilitiesBSH(ctrl%nstates, traj%coeff_diag_s, traj%preprob_old_s3, &
    &traj%preprob_s3, traj%state_diag, traj%hopprob_s)
    traj%hopprob_s(traj%state_diag)=0.d0
  elseif (ctrl%integrator==1 .or. ctrl%integrator==2) then
    select case (ctrl%hopping_procedure)
      case (0)
        traj%hopprob_s=0.
      case (1)
        call calculate_probabilities(ctrl%nstates, traj%coeff_diag_old_s, traj%coeff_diag_s, &
        &traj%Rtotal_ss, traj%state_diag, traj%hopprob_s)
      case (2)
        call calculate_probabilities_GFSH(ctrl%nstates, traj%coeff_diag_old_s, traj%coeff_diag_s, &
        &traj%state_diag, traj%hopprob_s)
      case default
        write(0,*) 'Unknown option ',ctrl%hopping_procedure,' to hopping_procedure!'
        stop 1
    endselect
  endif

  ! forced hop to ground state if energy gap to state above is smaller than threshold
  traj%kind_of_jump=0
  if ( ctrl%force_hop_to_gs>0. ) then
    deltaE=abs( real(traj%H_diag_ss(traj%state_diag,traj%state_diag) - traj%H_diag_ss(1,1)) )
    if ( (traj%state_diag/=1).and.(deltaE<=ctrl%force_hop_to_gs) ) then
      traj%hopprob_s=0.
      traj%hopprob_s(1)=1.
      traj%kind_of_jump=4
    endif
    if (traj%state_diag==1) then
      traj%hopprob_s=0.
    endif
  endif

  if (printlevel>2) then
    write(u_log,*) 'Old and new occupancies and hopping probabilities:'
    do istate=1,ctrl%nstates
      write(u_log,'(3(F14.9,3X))') abs(traj%coeff_diag_old_s(istate))**2,&
      &abs(traj%coeff_diag_s(istate))**2, traj%hopprob_s(istate)
    enddo
  endif

  call random_number(randnum)
  traj%randnum=randnum
  cumuprob=0.d0
  if (.not. traj%kind_of_jump==4) then
    traj%kind_of_jump=0
  endif
  deltaE=0.d0

  stateloop: do istate=1,ctrl%nstates
    ! calculate cumulative probability
    cumuprob=cumuprob + traj%hopprob_s(istate)
    if (cumuprob > randnum) then

      ! check for laser resonance
      if (ctrl%laser/=0) then
        do ilaser=1,ctrl%nlasers
          deltaE=abs(real(traj%H_diag_ss(istate,istate)-traj%Epot))
          deltaE=abs(deltaE-ctrl%laserenergy_tl(traj%step*ctrl%nsubsteps+1,ilaser))
          if (deltaE<=ctrl%laser_bandwidth) then
            traj%state_diag_old=traj%state_diag
            if (ctrl%hopping_procedure/=0) then
              traj%state_diag=istate
            endif
            traj%Epot=real(traj%H_diag_ss(istate,istate))
            traj%kind_of_jump=3   ! induced hop
            exit stateloop        ! ************************************************* exit of loop
          endif
        enddo
      endif

      ! check for frustration
      select case (ctrl%ekincorrect)
        case (0)    ! no frustrated jumps: just go on
          continue

        case (1)    ! correct along v, full kinetic energy available
          Ekin_masked=Calculate_ekin_masked(ctrl%natom, traj%veloc_ad, traj%mass_a, ctrl%atommask_a)
          ! Use Epot+Ekin
!           Emax=traj%H_diag_ss(traj%state_diag,traj%state_diag) + Ekin_masked
          ! Use Etot
          Emax=traj%Etot- traj%Ekin + Ekin_masked
          if (real(traj%H_diag_ss(istate,istate)) > Emax) then
            traj%kind_of_jump=2
            traj%state_diag_frust=istate
            exit stateloop         ! ************************************************* exit of loop
          endif

        case (2)    ! correct along projected v, vibrational kinetic energy available
          call available_ekin(ctrl%natom,&
          &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag, istate,:,:),&
          &traj%mass_a, sum_kk, sum_vk)
          deltaE=4.d0*sum_kk*(traj%Etot-traj%Ekin-&
          &real(traj%H_diag_ss(istate,istate)))+sum_vk**2
          if (deltaE<0.d0) then
            traj%kind_of_jump=2
            traj%state_diag_frust=istate
            exit stateloop         ! ************************************************* exit of loop
          endif


        case (3)    ! correct along T, less energy available
          call available_ekin(ctrl%natom,&
          &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag, istate,:,:),&
          &traj%mass_a, sum_kk, sum_vk)
          deltaE=4.d0*sum_kk*(traj%Etot-traj%Ekin-&
          &real(traj%H_diag_ss(istate,istate)))+sum_vk**2
          if (deltaE<0.d0) then
            traj%kind_of_jump=2
            traj%state_diag_frust=istate
            exit stateloop         ! ************************************************* exit of loop
          endif

        case (4)    ! correct along gradient difference
          call available_ekin(ctrl%natom,&
          &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag, istate,:,:),&
          &traj%mass_a, sum_kk, sum_vk)
          deltaE=4.d0*sum_kk*(traj%Etot-traj%Ekin-&
          &real(traj%H_diag_ss(istate,istate)))+sum_vk**2
          if (deltaE<0.d0) then
            traj%kind_of_jump=2
            traj%state_diag_frust=istate
            exit stateloop         ! ************************************************* exit of loop
          endif

        case (5)    ! correct along projected NAC/M
          call available_ekin(ctrl%natom,&
          &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag, istate,:,:),&
          &traj%mass_a, sum_kk, sum_vk)
          deltaE=4.d0*sum_kk*(traj%Etot-traj%Ekin-&
          &real(traj%H_diag_ss(istate,istate)))+sum_vk**2
          if (deltaE<0.d0) then
            traj%kind_of_jump=2
            traj%state_diag_frust=istate
            exit stateloop         ! ************************************************* exit of loop
          endif

        case (6)    ! correct along projected gradient difference
          call available_ekin(ctrl%natom,&
          &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag, istate,:,:),&
          &traj%mass_a, sum_kk, sum_vk)
          deltaE=4.d0*sum_kk*(traj%Etot-traj%Ekin-&
          &real(traj%H_diag_ss(istate,istate)))+sum_vk**2
          if (deltaE<0.d0) then
            traj%kind_of_jump=2
            traj%state_diag_frust=istate
            exit stateloop         ! ************************************************* exit of loop
          endif

        case (7)    ! correct along effective NAC/M
          call available_ekin(ctrl%natom,&
          &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag, istate,:,:),&
          &traj%mass_a, sum_kk, sum_vk)
          deltaE=4.d0*sum_kk*(traj%Etot-traj%Ekin-&
          &real(traj%H_diag_ss(istate,istate)))+sum_vk**2
          if (deltaE<0.d0) then
            traj%kind_of_jump=2
            traj%state_diag_frust=istate
            exit stateloop         ! ************************************************* exit of loop
          endif
 
        case (8)    ! correct along projected effective NAC/M
          call available_ekin(ctrl%natom,&
          &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag, istate,:,:),&
          &traj%mass_a, sum_kk, sum_vk)
          deltaE=4.d0*sum_kk*(traj%Etot-traj%Ekin-&
          &real(traj%H_diag_ss(istate,istate)))+sum_vk**2
          if (deltaE<0.d0) then
            traj%kind_of_jump=2
            traj%state_diag_frust=istate
            exit stateloop         ! ************************************************* exit of loop
          endif

      endselect

      ! neither in resonance nor frustrated, we have a surface hop!
      traj%state_diag_old=traj%state_diag
      if (ctrl%hopping_procedure/=0) then
        traj%state_diag=istate
        traj%Epot=real(traj%H_diag_ss(istate,istate))
      endif
      if (.not. traj%kind_of_jump==4) then
        traj%kind_of_jump=1
      endif
      exit stateloop               ! ************************************************* exit of loop

    endif
  enddo stateloop

  if (printlevel>2) then
    write(u_log,*) 
    write(u_log,'(A,1X,F12.9)') 'Random number:',randnum
    if (ctrl%hopping_procedure==0) then
      write(u_log,*) 
      write(u_log,'(A,1X,F12.9)') '*** Hopping is forbidden: new state = old state ***'
      write(u_log,*) 
    endif
    select case (traj%kind_of_jump)
      case (0)
        write(u_log,*) 'No jump performed.'
      case (1)
        write(u_log,'(A,1X,I4,1X,A)') 'Old state:',traj%state_diag_old,'(diag)'
        write(u_log,'(A,1X,I4,1X,A)') 'New state:',traj%state_diag,'(diag)'
        write(u_log,'(A,1X,F12.9,1X,A)') 'Jump is not frustrated, new Epot=',traj%Epot*au2ev,'eV'
      case (2)
        write(u_log,'(A,1X,I4,1X,A,1X,E12.5,1X,A)') 'Jump to state',istate,'is frustrated by',deltaE*au2eV,'eV.'
        if (ctrl%laser==2) then
          write(u_log,'(A,1X,F16.9,1X,A,1X,F16.9,1X,A)') &
          &'Detuning:',deltaE*au2eV,'eV, Laser Bandwidth:',ctrl%laser_bandwidth*au2eV,'eV'
        endif
      case (3)
        write(u_log,'(A)') 'Jump is in resonance with laser.'
        write(u_log,'(A,1X,F16.9,1X,A,1X,F16.9,1X,A)') &
        &'Detuning:',deltaE*au2eV,'eV, Laser Bandwidth:',ctrl%laser_bandwidth*au2eV,'eV'
      case (4)
        write(u_log,'(A)') 'Forced hop to ground state.'
        write(u_log,'(A,1X,I4,1X,A)') 'Old state:',traj%state_diag_old,'(diag)'
        write(u_log,'(A,1X,I4,1X,A)') 'New state:',traj%state_diag,'(diag)'
    endselect
  endif

endsubroutine

! ===========================================================

!> this routine calculates the hopping probabilities based on
!> the old and new coefficients and the propagator matrix
!> negative probabilities and the probability to stay in the same state are set to zero
subroutine calculate_probabilities(n, c0, c, R, state, prob)
  implicit none
  integer, intent(in) :: n, state
  complex*16, intent(in) :: c0(n), c(n), R(n,n)
  real*8, intent(out) :: prob(n)

  complex*16 :: w,x,y,z
  integer :: i
  real*8 :: sump

  prob=0.d0

  ! this numbers are the same for all states
  w=conjg(c(state))*c(state)
  x=conjg(c0(state))*c0(state)
  y=real(x-c(state)*conjg(R(state,state))*conjg(c0(state)))

  ! if population of active state increases, all probabilities are zero
  if ( (1.d0 - real(w/x))>0.d0) then
    do i=1,n
      if (i==state) cycle
      ! this number changes for each state
      z=real(c(i)*conjg(R(i,state))*conjg(c0(state)))
      prob(i)=max(0.d0, (1.d0-real(w/x))*real(z/y) )
    enddo
  endif

  ! renormalize, if sum of probabilities is above 1
  sump=sum(prob)
  if (sump>1.d0) prob=prob/sump

endsubroutine

! ===========================================================

!> this routine calculates the hopping probabilities based on
!> the equation for GFSH of Prezhdo et al.
subroutine calculate_probabilities_GFSH(n, c0, c, state, prob)
  implicit none
  integer, intent(in) :: n, state
  complex*16, intent(in) :: c0(n), c(n)
  real*8, intent(out) :: prob(n)

  real*8 :: w,x,y,z
  real*8 :: rho0(n),rho(n),drho(n)
  integer :: i
  real*8 :: sump

  prob=0.d0

  ! calculate rhos
  do i=1,n
    rho0(i)=real(conjg(c0(i))*c0(i))
    rho(i) =real(conjg(c(i)) *c(i) )
    drho(i)=rho(i)-rho0(i)
  enddo

  ! same for all states
  w=rho(state)
  x=rho0(state)
  y=0.
  do i=1,n
    if (drho(i)<0.d0) y=y-drho(i)
  enddo

  ! if population of active state increases, all probabilities are zero
  if ( (1.d0 - w/x)>0.d0) then
    do i=1,n
      if (i==state) cycle
      if (drho(i)<0.d0) cycle
      ! this number changes for each state
      z=drho(i)
      prob(i)=max(0.d0, (1.d0-w/x)*(z/y) )
    enddo
  endif

  ! renormalize, if sum of probabilities is above 1
  sump=sum(prob)
  if (sump>1.d0) prob=prob/sump

endsubroutine

! ==================================================================================================
! ==================================================================================================
! ===                                        Surface Switching
! ==================================================================================================
! ==================================================================================================
!> - Yinan
!> this routine performs all steps of the surface switching algorithm:
!> - calculate the surface switching probabilities (update traj%switchprob_s)
!> - generate a random number
!> - find the new decoherent state
subroutine surface_switching(traj,ctrl)
  use definitions
  use matrix
  use nuclear
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

  integer :: istep


  integer :: istate, jstate, iatom, idir
  real*8 :: dsum
  complex*16 :: ccoeff_diag_old_s(ctrl%nstates)
  real*8 :: randnum, cumuprob
  complex*16 :: w(ctrl%nstates,ctrl%nstates),tw(ctrl%nstates,ctrl%nstates),dw(ctrl%nstates,ctrl%nstates)
  complex*16 :: NACT(ctrl%nstates,ctrl%nstates), pNACT(ctrl%nstates,ctrl%nstates)


  if (printlevel>2) then
    write(u_log,*) '============================================================='
    write(u_log,*) '     Doing Surface Switching, Change of The Pointer State    '
    write(u_log,*) '============================================================='
  endif

  ! calculate the switching probabilities
  select case (ctrl%switching_procedure)
    case (0)
      traj%switchprob_s=0.d0
    case (1)
      if (printlevel>2) then
        write(u_log,*) 'switching probabilities computed with CSDM'
      endif

      !==== start reinitialization
      ! compute the D, used to reinitialize
      select case (ctrl%coupling)
        case (0)  ! ddt
          ! compute magnitude of NAC in diagonal basis
          traj%dmag=0.d0
          do istate=1,ctrl%nstates
            if (istate .ne. traj%state_diag) then
              dsum=0.d0
              do iatom=1,ctrl%natom
              do idir=1,3
                dsum=dsum+traj%NACdR_diag_ssad(istate,traj%state_diag,iatom,idir)*conjg(traj%NACdR_diag_ssad(istate,traj%state_diag,iatom,idir))
              enddo
              enddo
            endif
            dsum=max(0.d0,dsum)
            dsum=sqrt(dsum)
            traj%dmag=traj%dmag+dsum
          enddo
        case (1)  ! ddr
          ! compute magnitude of NAC in diagonal
          traj%dmag=0.d0
          do istate=1,ctrl%nstates
            if (istate .ne. traj%state_diag) then
              dsum=0.d0
              do iatom=1,ctrl%natom
              do idir=1,3
                dsum=dsum+traj%NACdR_diag_ssad(istate,traj%state_diag,iatom,idir)*conjg(traj%NACdR_diag_ssad(istate,traj%state_diag,iatom,idir))
              enddo
              enddo
            endif
            dsum=max(0.d0,dsum)
            dsum=sqrt(dsum)
            traj%dmag=traj%dmag+dsum
          enddo 
        case (2)  ! overlap
          ! compute magnitude of NACdt
          traj%dmag=0.d0
          do istate=1,ctrl%nstates
            if (istate .ne. traj%state_diag) then
              dsum=abs(traj%NACdt_ss(istate,traj%state_diag))  ! CSDM-C
              !dsum=abs(traj%overlaps_ss(istate,traj%state_diag)-traj%overlaps_ss(traj%state_diag,istate))/2  ! HST-scheme
            endif
            traj%dmag=traj%dmag+dsum
          enddo
        case (3)  ! ktdc
          ! compute magnitude of NACdt
          traj%dmag=0.d0
          do istate=1,ctrl%nstates
            if (istate .ne. traj%state_diag) then
              dsum=abs(traj%NACdt_ss(istate,traj%state_diag))  ! CSDM-C
            endif
            traj%dmag=traj%dmag+dsum
          enddo

      endselect

      if (traj%step .gt. 2) then
        if (printlevel>4) then         
          write(u_log,*) 'Magnitudes of NACdr in for current and previous two steps:'
          write(u_log,'(3(F8.5,3X))') traj%dmag, traj%dmag1, traj%dmag2
        endif
      endif

      ! reinitialize the ccoeff_diag_s based on the local minimum of dmag
      if (traj%step .gt. 2) then
        if ( (traj%dmag1-traj%dmag2) .lt. 0.d0 .AND. (traj%dmag-traj%dmag1) .gt. 0.d0 ) then
          if (printlevel>2) then
            if (ctrl%coupling .eq. 2) then
              write(u_log,*) 'Using overlap coupling, reinitialization based on CSDM-C'
              write(u_log,*) 'Coherent density has been reinitialized for CSDM-C'
            else
              write(u_log,*) 'Coherent density has been reinitialized for CSDM'
            endif
          endif
          traj%ccoeff_diag_s=traj%coeff_diag_s
          traj%ccoeff_MCH_s=traj%coeff_MCH_s
          if (ctrl%integrator==0) then
            traj%gRccoeff_MCH_s=traj%gRcoeff_MCH_s
            traj%gIccoeff_MCH_s=traj%gIcoeff_MCH_s
          endif
        endif
      endif

      traj%dmag2=traj%dmag1
      traj%dmag1=traj%dmag
      !===== end reinitialization
 
      !==== compute switching probability
      if (ctrl%integrator==0) then
        call calculate_probabilitiesBSH(ctrl%nstates, traj%ccoeff_diag_s, traj%preprob_old_s3, &
        &traj%preprob_s3, traj%state_diag, traj%switchprob_s)
      elseif (ctrl%integrator==1 .or. ctrl%integrator==2) then
        call calculate_probabilitiesDoM(ctrl%nstates, traj%ccoeff_diag_old_s, traj%ccoeff_diag_s, &
        &traj%Rtotal_ss, traj%state_diag, traj%switchprob_s)
      endif

    case (2)
      if (printlevel>2) then
        write(u_log,*) 'switching probabilities computed with SCDM'
      endif
      !==== compute switching probability
      if (ctrl%integrator==0) then
        call calculate_probabilitiesBSH(ctrl%nstates, traj%ccoeff_diag_s, traj%preprob_old_s3, &
        &traj%preprob_s3, traj%state_diag, traj%switchprob_s)
      elseif (ctrl%integrator==1 .or. ctrl%integrator==2) then
        call calculate_probabilitiesDoM(ctrl%nstates, traj%ccoeff_diag_old_s, traj%ccoeff_diag_s, &
        &traj%Rtotal_ss, traj%state_diag, traj%switchprob_s)
      endif

    case (3)
      if (printlevel>2) then
        write(u_log,*) 'switching probabilities computed with NDM'
      endif
      !==== compute switching probability
      if (ctrl%integrator==0) then
        call calculate_probabilitiesBSH(ctrl%nstates, traj%coeff_diag_s, traj%preprob_old_s3, &
        &traj%preprob_s3, traj%state_diag, traj%switchprob_s)
      elseif (ctrl%integrator==1 .or. ctrl%integrator==2) then
        call calculate_probabilitiesDoM(ctrl%nstates, traj%coeff_diag_old_s, traj%coeff_diag_s, &
        &traj%RDtotal_ss, traj%state_diag, traj%switchprob_s)
      endif

    case default
      write(0,*) 'Unknown option ',ctrl%switching_procedure,' to switching_procedure!'
      stop 1
  endselect

  if (printlevel>2) then
    write(u_log,*) 'switching probabilities:'
    do istate=1,ctrl%nstates
      write(u_log,'((F14.9,3X))') traj%switchprob_s(istate)
    enddo
  endif

  call random_number(randnum)
  traj%randnum=randnum
  
  cumuprob=0.d0
  traj%kind_of_jump=0
 
  stateloop:   do istate=1,ctrl%nstates
    ! calculate cumulative probability
    cumuprob=cumuprob + traj%switchprob_s(istate)
    if (cumuprob .ge. traj%randnum) then
      !TO DO: laser resonance
      if (printlevel>4) then
        write(u_log,'(A,1X,I4,1X,A,1X,I4,1X,A)') ' previous step selected state:',traj%state_diag, '(diag)', traj%state_MCH, '(MCH)'
      endif
      traj%state_diag_old=traj%state_diag
      traj%state_MCH_old=traj%state_MCH
      if (ctrl%switching_procedure/=0) then
        traj%state_diag=istate
        traj%state_MCH=state_diag_to_MCH(ctrl%nstates,traj%state_diag,traj%U_ss)
        if (printlevel>4) then
          write(u_log,'(A,1X,I4,1X,A,1X,I4,1X,A)') ' selected state:',traj%state_diag, '(diag)', traj%state_MCH, '(MCH)'
        endif
      endif
      exit stateloop         ! ************************************************* exit of loop
    endif
  enddo stateloop

  if (traj%step .gt. 1) then
    if (traj%state_diag_old .ne. traj%state_diag) then
      traj%kind_of_jump=1
    endif
  endif

  if (printlevel>2) then
    write(u_log,*)
    write(u_log,'(A,1X,F12.9)') 'Random number:', traj%randnum
    if (ctrl%switching_procedure==0) then
      write(u_log,*)
      write(u_log,'(A,1X,F12.9)') '*** Surface Switching is turned off: new state = old state ***'
      write(u_log,'(A,1X,F12.9)') 'Propagating Ehrenfest Dynamics'
      write(u_log,*)
    endif
    select case (traj%kind_of_jump)
      case (1)
        write(u_log,*) 'Pointer state switched'
        write(u_log,'(A,1X,I4,1X,A,1X,I4,1X,A)') 'Old pointer state:',traj%state_diag_old,'(diag)',traj%state_MCH_old,'(MCH)'
        write(u_log,'(A,1X,I4,1X,A,1X,I4,1X,A)') 'New pointer state:',traj%state_diag,'(diag)',traj%state_MCH,'(MCH)'
    endselect
  endif

endsubroutine

! ===========================================================
!> this routine calculates the hopping probabilities based on
!> the old and new coefficients and the propagator matrix
!> negative probabilities are set to zero
!> the probability to stay in the same state is 1-sump
subroutine calculate_probabilitiesDoM(n, c0, c, R, state, prob)
  implicit none
  integer, intent(in) :: n, state
  complex*16, intent(in) :: c0(n), c(n), R(n,n)
  real*8, intent(out) :: prob(n)

  complex*16 :: w,x,y,z
  integer :: i
  real*8 :: sump

  prob=0.d0

  ! this numbers are the same for all states
  w=conjg(c(state))*c(state)
  x=conjg(c0(state))*c0(state)
  y=real(x-c(state)*conjg(R(state,state))*conjg(c0(state)))

  ! if population of active state increases, all probabilities are zero
  if ( (1.d0 - real(w/x))>0.d0) then
    do i=1,n
      if (i==state) cycle
      ! this number changes for each state
      z=real(c(i)*conjg(R(i,state))*conjg(c0(state)))
      prob(i)=max(0.d0, (1.d0-real(w/x))*real(z/y) )
    enddo
  endif

  ! renormalize, if sum of probabilities is above 1
  sump=sum(prob)
  if (sump>1.d0) prob=prob/sump

  prob(state)=1.0-sump

endsubroutine

! ===========================================================
!> calculate switching probability for CSDM
!> calculate_probabilities_CSDM(ctrl%nstates, ctrl%natom, ccoeff_diag_s,
!traj%veloc_ad, traj%Gmatrix_ssad,&
!> &traj%H_diag_ss, ctrl%dtstep, traj%state_diag, traj%switchprob_s)
subroutine calculate_probabilities_CSDM(ns, n, cc, v, g, e, dt, k, prob)
  use definitions, only: u_log
  implicit none
  integer, intent(in) :: ns, n, k
  complex*16, intent(in) :: cc(ns)
  real*8, intent(in) :: v(n,3)
  complex*16, intent(in) :: g(ns,ns,n,3),e(ns,ns)
  real*8, intent(in) :: dt
  real*8, intent(out) :: prob(ns)

  integer :: istate, jstate, iatom, idir
  real*8 :: dvec(ns,ns,n,3)
  complex*16 :: den(ns,ns)
  real*8 :: vdotd(ns,ns)
  real*8 :: b(ns)
  real*8 :: sump

! NAC is now in diag basis  
  do istate=1,ns
    do jstate=1,ns
      dvec(istate,jstate,:,:)=real(g(istate,jstate,:,:)/(e(jstate,jstate)-e(istate,istate)))
      den(istate,jstate)=cc(istate)*conjg(cc(jstate))
    enddo
  enddo

! compute v dot dvec
  do istate=1,ns
    do jstate=1,ns
      vdotd(istate,jstate)=0.d0
      do iatom=1,n
        do idir=1,3
          vdotd(istate,jstate)=vdotd(istate,jstate)+dvec(istate,jstate,iatom,idir)*v(iatom,idir)
        enddo
      enddo
    enddo
  enddo

! now compute b_KK'
  do istate=1,ns
    b(istate)=0.d0
    if (istate .ne. k) then
      do iatom=1,n
        do idir=1,3
          b(istate)=2.d0*real(den(istate,k)*vdotd(k,istate))
        enddo
      enddo
    endif
  enddo

! now compute switching probability 

  do istate=1,ns
    if (istate .ne. k) then
      prob(istate)=real(b(istate)*dt)/real(den(k,k))
    endif
  enddo
  prob(k)=1.d0

  do istate=1,ns
    if (istate .ne. k) then
      prob(istate)=max(prob(istate),0.d0)
      prob(istate)=min(prob(istate),1.d0)
      prob(k)=prob(k)-prob(istate)
    endif
  enddo

  ! renormalize, if sum of probabilities is above 1
  sump=sum(prob)
  if (sump>1.d0) prob=prob/sump

endsubroutine

! ===========================================================
!> calculate switching probability for NDM
subroutine calculate_probabilities_NDM(ns, n, cc, v, g, e, dt, k, prob)
  use definitions, only: u_log
  implicit none
  integer, intent(in) :: ns, n, k
  complex*16, intent(in) :: cc(ns)
  real*8, intent(in) :: v(n,3)
  complex*16, intent(in) :: g(ns,ns,n,3),e(ns,ns)
  real*8, intent(in) :: dt
  real*8, intent(out) :: prob(ns)

  integer :: istate, jstate, iatom, idir
  real*8 :: dvec(ns,ns,n,3)
  complex*16 :: den(ns,ns)
  real*8 :: vdotd(ns,ns)
  real*8 :: b(ns)
  real*8 :: sump

! NAC is now in diag basis  
  do istate=1,ns
    do jstate=1,ns
      dvec(istate,jstate,:,:)=real(g(istate,jstate,:,:)/(e(jstate,jstate)-e(istate,istate)))
    enddo
  enddo

! compute v dot dvec
  do istate=1,ns
    do jstate=1,ns
      vdotd(istate,jstate)=0.d0
      do iatom=1,n
        do idir=1,3
          vdotd(istate,jstate)=vdotd(istate,jstate)+dvec(istate,jstate,iatom,idir)*v(iatom,idir)
        enddo
      enddo
    enddo
  enddo

! now compute b_KK'
  do istate=1,ns
    b(istate)=0.d0
    if (istate .ne. k) then
      do iatom=1,n
        do idir=1,3
          b(istate)=2.d0*real(den(istate,k)*vdotd(k,istate))
        enddo
      enddo
    endif
  enddo

! now compute switching probability 
  do istate=1,ns
    if (istate .ne. k) then
      prob(istate)=real(b(istate))*dt/real(den(k,k))
    endif
  enddo
  prob(k)=1.d0

  do istate=1,ns
    if (istate .ne. k) then
      prob(istate)=max(prob(istate),0.d0)
      prob(istate)=min(prob(istate),1.d0)
      prob(k)=prob(k)-prob(istate)
    endif
  enddo

  ! renormalize, if sum of probabilities is above 1
  sump=sum(prob)
  if (sump>1.d0) prob=prob/sump

endsubroutine

!===========================================================
!> calculate_probabilities_preBSH(ctrl%natom, ctrl%nstates, traj%veloc_ad,
!>   traj%coeff, traj%NACdR_diag_ssad, traj%state_diag, traj%hopprob_s)
!> calculate pre-probabilities when use BSH integrator, these are b_ik terms
subroutine calculate_probabilities_preBSH(n, ns, v, c, dvec, k, preprob)
  implicit none
  integer, intent(in) :: n, ns
  real*8, intent(in) :: v(n,3)
  complex*16, intent(in) :: c(ns), dvec(ns,ns,n,3)
  integer, intent(in) :: k
  real*8, intent(out) :: preprob(n,3)

  integer :: istate, jstate, iatom, idir
  real*8 :: vdotd(ns,ns)
  real*8 :: sump

  preprob=0.d0

  vdotd=0.d0
  ! compute v dot dvec
  do istate=1,ns
    if (istate .ne. k) then
      do iatom=1,n
        do idir=1,3
          vdotd(istate,k)=vdotd(istate,k)+dvec(k,istate,iatom,idir)*v(iatom,idir)
        enddo
      enddo
    endif
  enddo 

  ! compute b_ik
  do istate=1,ns
    preprob(istate,1)=2*real((c(istate)*conjg(c(k)))*vdotd(istate,k))
    preprob(istate,2)=max(preprob(istate,1),0.d0)
    preprob(istate,3)=min(preprob(istate,1),0.d0)
  enddo

endsubroutine

! ===========================================================
!> calculate_probabilitiesBSH(ctrl%nstates, traj%ccoeff_diag, traj%preprob_old_s3,
!> traj%preprob_s3, traj%state_diag, traj%switchprob_s)
!> This subroutine determines the hopping/switching probabilities according to
!> prob = (integral(bjk+) + integral(bjk-))/(density(k))
!> See M.D. Hack, A.W. Jasper, Y.L. Volobuev, D.W. Schwenke, D.G. Truhlar, JPCA
!> 103, 6309 (1999).
subroutine  calculate_probabilitiesBSH(ns, c, prep_old, prep, k, prob)
  use definitions
  use matrix
  implicit none
  integer, intent(in) :: ns
  complex*16, intent(in) :: c(ns)
  real*8, intent(in) :: prep_old(ns,3), prep(ns,3) 
  integer, intent(in) :: k
  real*8, intent(out) :: prob(ns)

  integer :: istate
  real*8 :: p(ns), pmax(ns), pmin(ns) 
  real*8 :: denk
  real*8 :: sump

  denk=c(k)*conjg(c(k))
 
  p=0.d0

  do istate=1,ns
    if (istate .ne. k) then 
      p(istate)=prep(istate,1)-prep_old(istate,1)
      pmax(istate)=max((prep(istate,2)-prep_old(istate,1)), 0.0d0)
      pmin(istate)=min((prep(istate,3)-prep_old(istate,3)), 0.0d0)

      p(istate)=min(max(p(istate)/denk, 0.0d0), 1.0d0)
      pmax(istate)=min(pmax(istate)/denk, 1.0d0)

      if ((denk+pmin(istate)) .le. 1.0d0) then
        pmin(istate)=pmin(istate)/(denk+pmin(istate))
      else
        pmin(istate)=max((1.0d0-denk), 0.0d0)
      endif
    endif
  enddo

  ! renormalize, if sum of probabilities is above 1
  sump=sum(p)
  if (sump>1.d0) p=p/sump
  p(k)=1.0-sump

  sump=sum(pmax)
  if (sump>1.d0) pmax=pmax/sump
  pmax(k)=1.0-sump

  sump=sum(pmin)
  if (sump>1.d0) pmin=pmin/sump
  pmin(k)=1.0-sump

  if (printlevel>4) then
    write(u_log,*) 'probabilities in BSH'
    do istate=1,ns
      write(u_log,'(3(F14.9,1X))') p(istate), pmax(istate), pmin(istate)
    enddo
  endif

  prob=p

endsubroutine

! ===========================================================

!> applies a decoherence correction to traj%coeff_diag_s
subroutine Decoherence(traj,ctrl)
  use definitions
  use matrix
  use decoherence_afssh
  use decoherence_dom
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl
  real*8 :: randnum
  complex*16 :: cpre(ctrl%nstates)

!The decoherence methods for TSH and SCP
  if (ctrl%method==0) then !TSH
    ! draw a new random number independently of the algorithm
    if (ctrl%compat_mode==0) then
      call random_number(randnum)
      traj%randnum2=randnum
    elseif (ctrl%compat_mode==1) then
      traj%randnum2=traj%randnum
    endif
 
    cpre = traj%coeff_diag_s

    if (ctrl%decoherence==1) then
      if (printlevel>2) then
        write(u_log,*) '============================================================='
        write(u_log,*) '              Decoherence (Granucci & Persico)'
        write(u_log,*) '============================================================='
      endif
      call EDC_step(traj,ctrl)
    elseif (ctrl%decoherence==2) then
      if (printlevel>2) then
        write(u_log,*) '============================================================='
        write(u_log,*) '           Decoherence (Jain, Alguire, Subotnik)'
        write(u_log,*) '============================================================='
      endif
      call afssh_step(traj,ctrl)
    endif

    if (ctrl%decoherence>0) then
      if (printlevel>2) then
        call vecwrite(ctrl%nstates, cpre, u_log, 'Coeff before decoherence','F12.9')
        call vecwrite(ctrl%nstates, traj%coeff_diag_s, u_log, 'Coeff after decoherence','F12.9')
      endif
    endif

  else if (ctrl%method==1) then !SCP
    if (ctrl%decoherence==11) then 
      if (printlevel>2) then
        write(u_log,*) '=============================================================='
        write(u_log,*) '             Decoherence (Decay of Mixing) - SCP'
        write(u_log,*) '=============================================================='
      endif
      if (ctrl%integrator==0) then
        call DoM_stepBS(traj,ctrl)
      elseif (ctrl%integrator==1 .or. ctrl%integrator==2) then
        call DoM_step(traj,ctrl)
      endif      
    endif

    if (ctrl%decoherence>0) then
      if (printlevel>2) then
        call vecwrite(ctrl%nstates, traj%coeff_diag_s, u_log, 'Coeff after decoherence','F12.9')
      endif
    endif

  endif
 
endsubroutine

! ===========================================================

!> is based on the EDC correction by Granucci and Persico
subroutine EDC_step(traj,ctrl)
  use definitions
  use matrix
  use decoherence_afssh
  use nuclear, only: Calculate_ekin_masked
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl
  integer :: istate
  real*8 :: tau0, tau, sumc, Ekin_masked
  complex*16 :: c(ctrl%nstates)

  real*8 :: ediff

  
!  tau0=1.d0 + ctrl%decoherence_alpha/traj%Ekin
  Ekin_masked=Calculate_ekin_masked(ctrl%natom, traj%veloc_ad, traj%mass_a, ctrl%atommask_a)
  tau0=1.d0 + ctrl%decoherence_alpha/Ekin_masked

  sumc=0.d0
  do istate=1,ctrl%nstates
    if (istate/=traj%state_diag) then
      ediff=abs( real(traj%H_diag_ss(istate,istate) ) - real(traj%H_diag_ss(traj%state_diag,traj%state_diag)) )
      if (ediff>0.d0) then
        tau=tau0 / ediff
         ! c(istate)=traj%coeff_diag_s(istate) * exp( -ctrl%dtstep / tau)
         ! Equation in "Critical appraisal ..." paper is wrong, exp() should be applied to populations,
         ! not coefficients, so for coefficients we need to add a factor of 1/2
        if (ctrl%decoherence==1) then
          c(istate)=traj%coeff_diag_s(istate) * exp( -0.5d0*ctrl%dtstep / tau)
        elseif (ctrl%decoherence==-1) then
           ! old variant of EDC with wrong equation, use keyword "edc_legacy"
          c(istate)=traj%coeff_diag_s(istate) * exp( -1.0d0*ctrl%dtstep / tau)
        endif
      else if (ediff==0.d0) then
        c(istate)=traj%coeff_diag_s(istate)
      endif
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

! ===========================================================

!> updates coeff_diag_s and state_MCH
subroutine Calculate_cMCH(traj,ctrl)
  use definitions
  use matrix
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

  ! calculate the MCH coefficients
  call matvecmultiply(ctrl%nstates,traj%U_ss,traj%coeff_diag_s,traj%coeff_MCH_s,'n')
  call matvecmultiply(ctrl%nstates,traj%U_ss,traj%ccoeff_diag_s,traj%ccoeff_MCH_s,'n')
  traj%state_MCH=state_diag_to_MCH(ctrl%nstates,traj%state_diag,traj%U_ss)


endsubroutine

! ===========================================================

!> updates steps_in_gs and checks whether the trajectory should be killed
!> the trajectory is actually killed in the main routine, so that proper finalization can be conducted.
subroutine kill_after_relaxation(traj,ctrl)
  use definitions
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

  if (ctrl%killafter>0) then
    if (traj%state_diag/=1) then
      traj%steps_in_gs=0
    else
      traj%steps_in_gs=traj%steps_in_gs+1
    endif
    if (traj%steps_in_gs>ctrl%killafter) then
      if (printlevel>0) then
        write(u_log,*) '============================================================='
        write(u_log,*) '              Ground state relaxation detected'
        write(u_log,*) '============================================================='
      endif
    endif
  endif

endsubroutine

! ===================================================

!> calculates the multiplicity of a state in nmstates nomenclature
!> see "canonical ordering of states" in SHARC manual
  integer function mult_of_nmstate(n, maxmult, mults)
  ! given a state n (in nmstates nomenclature, i.e. with multiplets expanded), 
  ! this routine returns the multiplicity of state n
  !
  ! Example: nstates = (2,0,2)
  ! We have 8 states, 2 singlets and 2 triplets with 3 components each
  ! the Order of states is:
  ! mult       M_S        n
  !    1       0          1
  !    1       0          2
  !    3      -1          1
  !    3      -1          2
  !    3       0          1
  !    3       0          2
  !    3      +1          1
  !    3      +1          2
  ! Multiplicity is counted up first, then M_s, then quantum number within multiplicity
  !
  ! X = mult_of_nmstate(7, 3, (/2,0,2/) ) 
  ! would yield X=3, because state 7 is a triplet
  implicit none
  integer, intent(in) :: n, maxmult
  integer, intent(in) :: mults(maxmult)

  integer :: imult, ims, istate, i

  i=1
  do imult=1,maxmult
    do ims=1,imult
      do istate=1,mults(imult)
        if (i==n) then
          mult_of_nmstate=imult
          return
        endif
        i=i+1
      enddo
    enddo
  enddo
  write(0,*) 'Error in mult_of_nmstate'
  stop 

  endfunction

endmodule electronic
