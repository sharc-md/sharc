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

!> # Module ELECTRONIC
!>
!> \author Sebastian Mai
!> \date 27.02.2015
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
  select case (ctrl%coupling)
    case (0)    ! ddt
      ! NADdt_ss can be directly used
      ! constant interpolation
      call unitary_propagator(&
        &ctrl%nstates,&
        &traj%H_MCH_ss, traj%H_MCH_old_ss,&
        &traj%NACdt_ss, traj%NACdt_old_ss,&
        &traj%U_ss,traj%U_old_ss,&
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
      call unitary_propagator(&
        &ctrl%nstates,&
        &traj%H_MCH_ss, traj%H_MCH_old_ss,&
        &traj%NACdt_ss, traj%NACdt_old_ss,&
        &traj%U_ss,traj%U_old_ss,&
        &ctrl%dtstep, ctrl%nsubsteps, 0,&       ! 0=linear interpolation
        &traj%Rtotal_ss)
    case (2)    ! overlap
      ! overlap matrix is used
      ! use LOCAL DIABATISATION
      call LD_propagator(&
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
    stop 1
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
    write(u_log,*) 'Old and new diagonal coefficients:'
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
  integer :: istate,ilaser
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

  ! forced hop to ground state if energy gap to state above is smaller than threshold
  traj%kind_of_jump=0
  if ( ctrl%force_hop_to_gs>0. ) then
    deltaE=abs( real(traj%H_diag_ss(traj%state_diag,traj%state_diag) - traj%H_diag_ss(1,1)) )
    if ( (traj%state_diag/=1).and.(deltaE<=ctrl%force_hop_to_gs) ) then
      traj%hopprob_s=0.
      traj%hopprob_s(1)=1.
      traj%kind_of_jump=4
      ! set also populations
      traj%coeff_diag_s=dcmplx(0.d0,0.d0)
      traj%coeff_diag_s(1)=dcmplx(1.d0,0.d0)
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

        case (2)    ! correct along T, less energy available
          call available_ekin(ctrl%natom,&
          &traj%veloc_ad,real(traj%gmatrix_ssad(traj%state_diag, istate,:,:)),&
          &traj%mass_a, sum_kk, sum_vk)
          deltaE=4.d0*sum_kk*(traj%Etot-traj%Ekin-&
          &real(traj%H_diag_ss(istate,istate)))+sum_vk**2
          if (deltaE<0.d0) then
            traj%kind_of_jump=2
            traj%state_diag_frust=istate
            exit stateloop         ! ************************************************* exit of loop
          endif

        case (3)    ! correct along gradient difference
          call available_ekin(ctrl%natom,&
          &traj%veloc_ad,real(traj%gmatrix_ssad(istate, istate,:,:)-&
          &traj%gmatrix_ssad(traj%state_diag, traj%state_diag,:,:)),&
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

! ===========================================================

!> applies a decoherence correction to traj%coeff_diag_s
subroutine Decoherence(traj,ctrl)
  use definitions
  use matrix
  use decoherence_afssh
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl
  real*8 :: randnum
  complex*16 :: cpre(ctrl%nstates)

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

!  tau0=1.d0 + ctrl%decoherence_alpha/traj%Ekin
  Ekin_masked=Calculate_ekin_masked(ctrl%natom, traj%veloc_ad, traj%mass_a, ctrl%atommask_a)
  tau0=1.d0 + ctrl%decoherence_alpha/Ekin_masked

  sumc=0.d0
  do istate=1,ctrl%nstates
    if (istate/=traj%state_diag) then
      tau=tau0 / abs( real(traj%H_diag_ss(istate,istate) ) - real(traj%H_diag_ss(traj%state_diag,traj%state_diag)) )
!       c(istate)=traj%coeff_diag_s(istate) * exp( -ctrl%dtstep / tau)
      ! Equation in "Critical appraisal ..." paper is wrong, exp() should be applied to populations,
      ! not coefficients, so for coefficients we need to add a factor of 1/2
      c(istate)=traj%coeff_diag_s(istate) * exp( -0.5d0*ctrl%dtstep / tau)
      sumc=sumc+abs(c(istate))**2
    endif
  enddo

  c(traj%state_diag)=traj%coeff_diag_s(traj%state_diag) * &
  &sqrt( (1.d0-sumc) / abs(traj%coeff_diag_s(traj%state_diag))**2)

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

endmodule
