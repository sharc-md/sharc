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

!> # Module  zpepumping
!> zpepumping, a decay of mixing algorithm
!> LP-ZPE and modified LP-ZPE
!> 
!> \author Yinan Shu
!> \date 11.18.2021
!>
!> LP-ZPE and modified LP-ZPE are implemented on Sep 12, 2022
!> This module contains subroutines that allow to perform a ZPE Pumping

module zpe
  implicit none

  ! ZPE pumping
  public initial_zpe
  public ZPEcorrection
  public ZPEpumping

  private compute_local_mode
  private pumping_state_switch
  private compute_pumping_direction_time
  private propagate_pumping_coeff

  ! LP-ZPE
  public LP_ZPE

  private LP_parallel_ekin
  private compute_parallel_velocity
  private correct_bond_velocity

  private correct_scheme1
  private correct_all_bonds
  private correct_bonds
  private correct_bcbond_velocity_method1
  private correct_bcbond_velocity_method2
  private correct_bcbond_velocity_method3
 
contains

! ===========================================================
! the wrap up for ZPE correction
  subroutine initial_zpe(traj,ctrl)
    use definitions
    use matrix
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl

    integer :: istate

    if (printlevel>3) then
      write(u_log,*) '============================================================='
      write(u_log,*) '    Initialization Zero Point Energy Correction Scheme'
      write(u_log,*) '============================================================='
    endif


    if (ctrl%zpe_correction==1) then ! Do ZPE Pumping
      ! initialize all pumping state to original state, which is state 1
      traj%state_pumping_s=1
 
      ! initialize all coefficients of state 1 to original state, state 2 to zero 
      do istate=1, ctrl%nstates
        traj%coeff_zpe_s2(istate,1)=dcmplx(1.0d0,0.d0)
        traj%coeff_zpe_s2(istate,2)=dcmplx(0.d0,0.d0)
      enddo 
    else if (ctrl%zpe_correction==2) then ! Do LP-ZPE
      traj%lpzpe_cycle=0
      traj%lpzpe_ke_ah=0.d0
      traj%lpzpe_ke_bc=0.d0  
      traj%lpzpe_iter_incycle=1
      traj%in_cycle=0
    endif

  endsubroutine


! ===========================================================
! the wrap up for ZPE correction
  subroutine ZPEcorrection(traj,ctrl)
    use definitions
    use matrix
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl


    if (ctrl%zpe_correction==1) then ! Do ZPE Pumping
      call ZPEpumping(traj,ctrl)
    else if (ctrl%zpe_correction==2) then ! Do LP-ZPE 
      call LP_ZPE(traj,ctrl)
    endif 

  endsubroutine


! ===========================================================
! Using ZPE pumping method to correct each electronic state. 
  subroutine ZPEpumping(traj,ctrl)
    use definitions
    use matrix
    use nuclear 
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl

    complex*16 :: Hpump(ctrl%nstates,2,2), Gradpump(ctrl%nstates,2,ctrl%natom,3)
    real*8 :: svec_pumping_sad(ctrl%nstates,ctrl%natom,3), pumpingtime_s(ctrl%nstates)
    complex*16 :: Rpumping(ctrl%nstates,2,2)
    complex*16 :: denpump(ctrl%nstates,2,2)

    real*8 :: veloc_app_ad(ctrl%natom,3)
    real*8 :: iso_veloc_app_ad(ctrl%natom,3)
    real*8 :: iso_veloc_ad(ctrl%natom,3), iso_geom_ad(ctrl%natom,3)
    real*8 :: iso_veloc_old_ad(ctrl%natom,3)
    real*8 :: iso_accel_ad(ctrl%natom,3)

    integer :: istate, jstate, ipstate, jpstate
    integer :: iatom, idir
    character(8000) :: string

    if (printlevel>3) then
      write(u_log,*) '============================================================='
      write(u_log,*) '                Zero Point Energy Pumping'
      write(u_log,*) '============================================================='
    endif
   
    ! initialize 
    traj%Evib_local_s=dcmplx(0.d0,0.d0)
    traj%Ezpe_local_s=dcmplx(0.d0,0.d0)
    traj%grad_Ezpe_local_sad=0.d0
    
    ! compute forward velocity 
    if (printlevel>4) write(u_log,*) 'Approximate velocity with forward velocity verlet'
    call VelocityVerlet_vstep_approximate(traj%veloc_old_ad, veloc_app_ad, traj%accel_ad, ctrl%dtstep, ctrl%natom)
 
    ! Compute isoinertial coordiantes, velocities and accelerations
    ! Notice we use scaled mass mu=1.0
    do iatom=1,ctrl%natom
      do idir=1,3 
        iso_geom_ad(iatom,idir)=traj%geom_ad(iatom,idir)*sqrt(traj%mass_a(iatom))
        iso_veloc_ad(iatom,idir)=traj%veloc_ad(iatom,idir)*sqrt(traj%mass_a(iatom))
        iso_veloc_old_ad(iatom,idir)=traj%veloc_old_ad(iatom,idir)*sqrt(traj%mass_a(iatom))
        iso_veloc_app_ad(iatom,idir)=veloc_app_ad(iatom,idir)*sqrt(traj%mass_a(iatom))
        iso_accel_ad(iatom,idir)=traj%accel_ad(iatom,idir)*sqrt(traj%mass_a(iatom))
      enddo
    enddo 

    ! Compute local mode energies and gradients. 
    call compute_local_mode(ctrl%natom, ctrl%nstates, iso_veloc_ad, &
      &iso_veloc_old_ad, iso_veloc_app_ad, iso_accel_ad, traj%mass_a, traj%Evib_local_s, &
      &traj%Ezpe_local_s, traj%grad_Ezpe_local_sad, traj%grad_MCH_sad, &
      &traj%grad_MCH_old_sad, traj%grad_MCH_old2_sad, traj%grad_ad, &
      &ctrl%dtstep, ctrl%dtstep_old, traj%step)

    if (printlevel>4) then
      call vecwrite(ctrl%nstates, traj%Evib_local_s, u_log, 'vibrational energy of local mode', 'F14.9')
      call vecwrite(ctrl%nstates, traj%Ezpe_local_s, u_log, 'zero point energy of local mode', 'F14.9')
      do istate=1,ctrl%nstates
        write(string,'(A20,I3)') 'Ezpe Gradient state ',istate
        call vec3write(ctrl%natom, traj%grad_Ezpe_local_sad(istate,:,:), u_log, trim(string), 'F14.9')
      enddo
    endif

    ! Construction pumping Hamiltonian and gradient
    Hpump=dcmplx(0.d0,0.d0)
    Gradpump=dcmplx(0.d0,0.d0)
    do istate=1,ctrl%nstates
      Hpump(istate,1,1)=traj%H_MCH_ss(istate,istate)
      Hpump(istate,2,2)=traj%H_MCH_ss(istate,istate)-traj%Ezpe_local_s(istate)
      Gradpump(istate,1,:,:)=traj%grad_MCH_sad(istate,:,:)
      Gradpump(istate,2,:,:)=traj%grad_MCH_sad(istate,:,:)-traj%grad_Ezpe_local_sad(istate,:,:)
    enddo

    ! Pumping state switch
    call pumping_state_switch(ctrl%nstates, traj%Evib_local_s, traj%Ezpe_local_s, &
      &traj%state_pumping_s)

    if (printlevel>4) then
      write(u_log,*) 'pumping state:'
      write(u_log,*) (traj%state_pumping_s(istate),istate=1,ctrl%nstates)
    endif

    ! Compute pumping direction and time
    call compute_pumping_direction_time(ctrl%natom, ctrl%nstates, traj%veloc_ad,&
      &traj%mass_a, traj%Evib_local_s, traj%Ezpe_local_s, traj%state_pumping_s, &
      &svec_pumping_sad, pumpingtime_s)

    if (printlevel>4) then
      call vecwrite(ctrl%nstates, pumpingtime_s, u_log, 'pumping time', 'F14.9')
    endif

    
    ! propagate pumping population, decay of mixing
    if (printlevel>4) then
      write(u_log,*) 'old pumping coefficients:'
      do istate=1,ctrl%nstates
        write(u_log,'(2(F14.9,1X),2X,2(F14.9,1X))') traj%coeff_zpe_s2(istate,1), traj%coeff_zpe_s2(istate,2)
      enddo
      do istate=1, ctrl%nstates
        do ipstate=1, 2
          do jpstate=1, 2
            denpump(istate,ipstate,jpstate)=traj%coeff_zpe_s2(istate,ipstate)*conjg(traj%coeff_zpe_s2(istate,jpstate))
          enddo
        enddo
      enddo
      write(u_log,*) 'old diagonal pumping density:'
      do istate=1,ctrl%nstates
        write(u_log,'(2(F14.9,1X),2X,2(F14.9,1X))') denpump(istate,1,1), denpump(istate,2,2)
      enddo
    endif

    call propagate_pumping_coeff(ctrl%nstates, ctrl%dtstep, ctrl%nsubsteps, &
      &pumpingtime_s, traj%state_pumping_s, traj%coeff_zpe_s2, Rpumping) 

    ! compute pumping density matrix
    do istate=1, ctrl%nstates
      do ipstate=1, 2
        do jpstate=1, 2
          denpump(istate,ipstate,jpstate)=traj%coeff_zpe_s2(istate,ipstate)*conjg(traj%coeff_zpe_s2(istate,jpstate))
        enddo
      enddo
    enddo

    if (printlevel>4) then
      write(u_log,*) 'new pumping coefficients:'
      do istate=1,ctrl%nstates
        write(u_log,'(2(F14.9,1X),2X,2(F14.9,1X))') traj%coeff_zpe_s2(istate,1), traj%coeff_zpe_s2(istate,2)
      enddo
      do istate=1,ctrl%nstates
        write(u_log,*) 'state: ', istate
        call matwrite(2, Rpumping(istate,:,:), u_log, 'Pumping propagator','F14.9')
      enddo 
      write(u_log,*) 'new diagonal pumping density:'
      do istate=1,ctrl%nstates
        write(u_log,'(2(F14.9,1X),2X,2(F14.9,1X))') denpump(istate,1,1), denpump(istate,2,2)
      enddo
    endif

    ! adjust the true Hamiltonian and gradients by pumping Hamiltonian and gradients
    do istate=1,ctrl%nstates
      traj%H_MCH_ss(istate,istate)=(denpump(istate,1,1)+denpump(istate,2,2))*traj%H_MCH_ss(istate,istate)-&
        &denpump(istate,2,2)*traj%Ezpe_local_s(istate)
      traj%grad_MCH_sad(istate,:,:)=(denpump(istate,1,1)+denpump(istate,2,2))*traj%grad_MCH_sad(istate,:,:)-&
        &denpump(istate,2,2)*traj%grad_Ezpe_local_sad(istate,:,:)
    enddo


  endsubroutine


! ==========================================================
!> this subroutine computes the local mode vibration energy and zpe 
!> compute_local_mode(ctrl%natom, ctrl%nstates, iso_veloc_ad, &
!>      &iso_veloc_old_ad, iso_veloc_app_ad, iso_accel_ad, traj%mass_a, traj%Evib_local_s, &
!>      &traj%Ezpe_local_s, traj%grad_Ezpe_local_sad, traj%grad_MCH_sad, &
!>      &traj%grad_MCH_old_sad, traj%grad_MCH_old2_sad, traj%grad_ad, &
!>      &ctrl%dtstep, ctrl%dtstep_old, traj%step)
  subroutine compute_local_mode(n, ns, v, vold, v_app, a, m, Evib, Ezpe, gEzpe, g, gold, gold2, traj_g, dt, dtold, step)
    use definitions, only: u_log, printlevel
    use matrix
    use nuclear, only: VelocityVerlet_vstep_approximate
    implicit none
    integer, intent(in) :: n, ns
    real*8, intent(in) :: v(n,3), vold(n,3), v_app(n,3), a(n,3), m(n)
    complex*16, intent(out) :: Evib(ns), Ezpe(ns)
    real*8, intent(out) :: gEzpe(ns,n,3)
    real*8, intent(in) :: g(ns,n,3), gold(ns,n,3), gold2(ns,n,3), traj_g(n,3)
    real*8, intent(in) :: dt, dtold
    integer, intent(in) :: step

    real*8 :: gv_old(ns), gv(ns) ! dV/dt
    real*8 :: eh(ns) ! d2V/dt2
    real*8 :: Ekin_app, mass_app
    real*8 :: normv, normg
    real*8 :: vg_app_ad(n,3)
    integer :: istate, iatom, idir



    if (printlevel>3) then
      write(u_log,*) '============================================================='
      write(u_log,*) 'Calculating local vibration mode energy and elevated surface'
      write(u_log,*) '============================================================='
    endif

    ! compute local mode force constant
    if (step>=1) then 
      gv_old=0.d0
      gv=0.d0
      normv=0.d0
      normg=0.d0
 
      ! here the gradient need to be re-scaled to iso coordiantes 
      do iatom=1,n
        do idir=1,3
          gv_old=gv_old+gold(:,iatom,idir)/sqrt(m(iatom))*vold(iatom,idir)
          gv=gv+g(:,iatom,idir)/sqrt(m(iatom))*v_app(iatom,idir)
          normv=normv+v_app(iatom,idir)**2
          normg=normg+(traj_g(iatom,idir)/sqrt(m(iatom)))**2
        enddo
      enddo
      normv=sqrt(normv)
      normg=sqrt(normg)
      eh=2.0d0*(gv-gv_old)/(dt+dtold)
 
      ! compute velocity along gradient direction 
      do iatom=1,n
        do idir=1,3
          vg_app_ad(iatom,idir)=v_app(iatom,idir)*(traj_g(iatom,idir)/sqrt(m(iatom)))/normg
        enddo
      enddo

      if (printlevel>4) then
        call vec3write(n, v_app, u_log, 'original velocity', 'F14.9')
        call vec3write(n, vg_app_ad, u_log, 'velocity along g direction', 'F14.9')
      endif
          
      Ekin_app=0.d0
      do iatom=1,n
        Ekin_app=Ekin_app+0.5d0*1.d0*sum(v_app(iatom,:)**2)
      enddo
      write(u_log,*) 'original kinetic energy', Ekin_app

      ! compute local vibration mode energy
      ! notice here the scaled mass is 1
      Ekin_app=0.d0
      do iatom=1,n
        Ekin_app=Ekin_app+0.5d0*1.d0*sum(vg_app_ad(iatom,:)**2)
      enddo
      write(u_log,*) 'kinetic energy', Ekin_app

      do istate=1,ns
        if (eh(istate)>0.d0) then
          Evib(istate)=Ekin_app+0.5d0*(gv(istate)**2)/eh(istate)
          Ezpe(istate)=0.5d0*sqrt(eh(istate)/1.0)
          ! compute the derivative of local vibration mode energy 
          if (step==1) then 
            gEzpe(istate,:,:)=(g(istate,:,:)-gold(istate,:,:))/dt
          else if (step>=2) then
            gEzpe(istate,:,:)=2.0d0*((g(istate,:,:)-gold(istate,:,:))/dt-(gold(istate,:,:)-gold2(istate,:,:))/dtold)/(dt+dtold)
          endif 
        else 
          Evib(istate)=Ekin_app
          Ezpe(istate)=0.d0
          gEzpe(istate,:,:)=0.d0
      endif
      enddo ! do istate=1,ns

     endif ! if (step>=1) then

  endsubroutine

! ==========================================================
!> this subroutine switches pumping state
  subroutine pumping_state_switch(ns, Evib, Ezpe, k)
    use definitions, only: u_log, printlevel
    implicit none
    integer, intent(in) :: ns
    complex*16, intent(in) :: Evib(ns), Ezpe(ns)
    integer, intent(inout) :: k(ns) 
   
    integer :: istate

    ! the pumping state switches according to relative value of potential and Ezpe
    do istate=1, ns
      if (real(Evib(istate)) .lt. real(Ezpe(istate))) then 
        k(istate)=2
      else if (real(Evib(istate)) .ge. real(Ezpe(istate))) then 
        k(istate)=1
      endif
    enddo 

  endsubroutine

! ==========================================================
!> this subroutine computes pumping direction and time
!> compute_pumping_direction_time(ctrl%natom, ctrl%nstates, traj%veloc_ad,&
!>   &traj%mass_a, traj%Evib_local_s, traj%Ezpe_local_s, traj%state_pumping_s, &
!>   &svec_pumping_sad, pumpingtime_s)
  subroutine compute_pumping_direction_time(n, ns, v, m, Evib, Ezpe, k, s, t) 
    use definitions, only: u_log, printlevel
    implicit none
    integer, intent(in) :: n, ns
    real*8, intent(in) :: v(n,3), m(n)
    complex*16, intent(in) :: Evib(ns), Ezpe(ns)
    integer, intent(in) :: k(ns)
    real*8, intent(out) :: s(ns,n,3)
    real*8, intent(inout) :: t(ns)

    real*8 :: normv
    real*8 :: p(n,3)
    real*8 :: smag, us(ns), es(ns)
    integer :: istate, iatom, idir

    ! compute norm of v 
    normv=0.d0
    do iatom=1,n
      do idir=1,3
        normv=normv+v(iatom,idir)**2
      enddo
    enddo
    normv=sqrt(normv)
 
    smag=0.d0
    do istate=1,ns

      s(istate,:,:)=v/normv
      do iatom=1,n
        do idir=1,3
          smag=smag+(s(istate,iatom,idir)**2)/m(iatom)
        enddo
      enddo

    us(istate)=0.d0
    do iatom=1,n
      do idir=1,3
        us(istate)=us(istate)+s(istate,iatom,idir)*v(iatom,idir)
      enddo
    enddo
    es(istate)=(us(istate)**2)/(2.d0*smag**2)
    t(istate)=dabs(real(Ezpe(istate)))/(1.d0+1.d0/es(istate))

   enddo !do istate=1,ns

  endsubroutine


! ==========================================================
!> this subroutine propagates the pumping coefficients
  subroutine propagate_pumping_coeff(ns, dt, nsub, tau, k, C, R)
    use definitions, only: u_log, printlevel
    use matrix
    implicit none
    integer, intent(in) :: ns
    integer, intent(in) :: nsub
    real*8, intent(in) :: dt
    real*8, intent(in) :: tau(ns)
    integer, intent(in) :: k(ns)
    complex*16, intent(inout) :: C(ns,2)
    complex*16, intent(inout) :: R(ns,2,2)

    integer :: istep
    real*8 :: dtsubstep
    complex*16 :: Rtotal(ns,2,2), Rprod(ns,2,2)
    complex*16 :: DP(ns,2,2), Rexpd(ns,2,2)
    complex*16 :: Ctmp1(ns,2), Cold(ns,2)
    complex*16 :: ii=dcmplx(0.d0,1.d0)
    complex*16 :: Rtmp(2,2)
 
    integer :: istate, ipstate
    integer :: tracker

    !initialize Ctmp1
    Ctmp1=C
    Cold=C

    dtsubstep=dt/nsub
    Rtotal=dcmplx(0.d0,0.d0)
    do istate=1,ns
      do ipstate=1,2
        Rtotal(istate,ipstate,ipstate)=dcmplx(1.d0,0.d0)
      enddo
    enddo 
    
    ! tracker is used to avoid den(k,k) being exact zero 
    ! if this happens, we use a single state decay, and the pointer state
    ! density is computed by 1-den(i,i)
    tracker=0
    do istate=1, ns
      do ipstate=1, 2
        if (Ctmp1(istate,k(istate)) .eq. dcmplx(0.d0,0.d0)) tracker=1
      enddo
    enddo

    if (tracker==0) then ! do two state case  
      do istep=1,nsub
        DP=dcmplx(0.d0,0.d0)
        do istate=1, ns
          do ipstate=1, 2
            if (ipstate .ne. k(istate)) then
              DP(istate,k(istate),k(istate))=DP(istate,k(istate),k(istate))+(1.d0,0.d0)*tau(istate)*real(Ctmp1(istate,ipstate)*conjg(Ctmp1(istate,ipstate)))
              DP(istate,ipstate,ipstate)=-0.5d0*(1.d0,0.d0)*tau(istate)
            endif
          enddo
          DP(istate,k(istate),k(istate))=0.5d0*DP(istate,k(istate),k(istate))/(real(Ctmp1(istate,k(istate))*conjg(Ctmp1(istate,k(istate)))))
        enddo    
        Rexpd=dtsubstep*(DP)
        do istate=1, ns    
          do ipstate=1,2
            Rexpd(istate,ipstate,ipstate)=exp(Rexpd(istate,ipstate,ipstate))
          enddo
          !call exponentiate(2, Rexpd(istate,:,:), (1.d0,0.d0))
          call matmultiply(2, Rexpd(istate,:,:), Rtotal(istate,:,:), Rprod(istate,:,:), 'nn')
          Rtotal(istate,:,:)=Rprod(istate,:,:)
          call matvecmultiply(2, Rtotal(istate,:,:), Cold(istate,:), Ctmp1(istate,:), 'n')
        enddo
      enddo ! do istep=1,nsub
      do istate=1, ns
        call matvecmultiply(2, Rtotal(istate,:,:), Cold(istate,:), C(istate,:), 'n')
      enddo 
      R=Rtotal
    else ! do one state decay 
      do istep=1,nsub
        DP=dcmplx(0.d0,0.d0)
        do istate=1, ns
          do ipstate=1, 2
            if (ipstate .ne. k(istate)) then
              DP(istate,ipstate,ipstate)=-0.5d0*(1.d0,0.d0)*tau(istate)
            endif
          enddo
        enddo
        Rexpd=dtsubstep*(DP)
        do istate=1, ns
          do ipstate=1,2
            Rexpd(istate,ipstate,ipstate)=exp(Rexpd(istate,ipstate,ipstate))
          enddo
          call matmultiply(2, Rexpd(istate,:,:), Rtotal(istate,:,:), Rprod(istate,:,:), 'nn')
          Rtotal(istate,:,:)=Rprod(istate,:,:)
          call matvecmultiply(2, Rtotal(istate,:,:), Cold(istate,:), Ctmp1(istate,:), 'n')
        enddo
      enddo
    do istate=1, ns
      call matvecmultiply(2, Rtotal(istate,:,:), Cold(istate,:), C(istate,:), 'n')
      if (k(istate)==1) then
        C(istate,1)=sqrt(1.d0-C(istate,2)*conjg(C(istate,2)))
      else if (k(istate)==2) then
        C(istate,2)=sqrt(1.d0-C(istate,1)*conjg(C(istate,1)))
      endif 
    enddo
    R=Rtotal
   endif ! if (tracker==0) then

  endsubroutine


! ==========================================================
!> main driver of LP-ZPE correction 
  subroutine LP_ZPE(traj,ctrl)
    use definitions
    use matrix
    use nuclear, only: Calculate_ekin 
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl

    real*8 :: tstart, tend
 
    integer :: atom1, atom2
    real*8 :: eleak, r_m_ah
    real*8 :: factor
    integer :: scheme
    integer :: do_correction
    real*8 :: eavail
    real*8 :: Ekin_old, Ekin_new
    integer :: ibond, iatom, jatom, katom, idir
    real*8 :: bond_Ekin_AH(ctrl%lpzpe_nah), bond_Ekin_AH_old(ctrl%lpzpe_nah)
    real*8 :: bond_Ekin_BC(ctrl%lpzpe_nbc), bond_Ekin_BC_old(ctrl%lpzpe_nbc)
    real*8 :: parallel_vec_tmp(3)
    real*8 :: parallel_veloc1_tmp, parallel_veloc2_tmp

    if (printlevel>3) then
      write(u_log,*) '============================================================='
      write(u_log,*) '       Local Pair - Zero Point Energy Correction Scheme'
      write(u_log,*) '============================================================='
      write(u_log,*) 'S. Mukherjee, M. Barbatti. J. Chem. Theory Comput. 2022, 18, 4109-4116'
    endif

    tstart=traj%lpzpe_cycle*ctrl%t_cycle
    tend=traj%lpzpe_cycle*ctrl%t_cycle + ctrl%t_check

    ! within this cumulation period, computing average value of kinetic energy
    if (traj%microtime>=tstart .and. traj%microtime<tend) then

      if (traj%lpzpe_iter_incycle==1) then 
        traj%lpzpe_starttime=traj%microtime-ctrl%dtstep
        traj%in_cycle=1
        write(u_log,*) 'At the start of a cycle'
        write(u_log,*) 'Cycle start time:', tstart
        write(u_log,*) 'Cycle end time:', tend
        write(u_log,*) 'Cycle begins at:', traj%lpzpe_starttime
        write(u_log,*) 'Reset cumulative AH and BC bonds kinetic energy to zero' 
        do ibond=1,ctrl%lpzpe_nah
          traj%lpzpe_ke_ah(ibond)=0.d0
        enddo
        do ibond=1,ctrl%lpzpe_nbc
          traj%lpzpe_ke_bc(ibond)=0.d0
        enddo
      endif

      if (printlevel>4) then
        write(u_log,*) 'Within a LP-ZPE cycle, current cycle:', traj%lpzpe_cycle
        write(u_log,*) 'current ieration in cycle:', traj%lpzpe_iter_incycle
        call vecwrite(ctrl%lpzpe_nah, traj%lpzpe_ke_ah, u_log, 'Cummulative AH bonds kinetic energy of last step', 'F14.9')
        call vecwrite(ctrl%lpzpe_nbc, traj%lpzpe_ke_bc, u_log, 'Cummulative BC bonds kinetic energy of last step', 'F14.9')
      endif

      ! Compute kinetic energy parallel to AH and BC bonds
      call LP_parallel_ekin(ctrl%natom, traj%geom_ad, traj%veloc_ad, traj%mass_a, traj%element_a,&
        &ctrl%dtstep, ctrl%lpzpe_nah, ctrl%lpzpe_nbc, ctrl%lpzpe_ah, ctrl%lpzpe_bc,&
        &traj%lpzpe_ke_ah, traj%lpzpe_ke_bc)

      if (printlevel>4) then
        call vecwrite(ctrl%lpzpe_nah, traj%lpzpe_ke_ah, u_log, 'Cummulative AH bonds kinetic energy', 'F14.9')
        call vecwrite(ctrl%lpzpe_nbc, traj%lpzpe_ke_bc, u_log, 'Cummulative BC bonds kinetic energy', 'F14.9')
      endif

      traj%lpzpe_iter_incycle=traj%lpzpe_iter_incycle+1

    endif
 
    ! cummulation is done, start computing averaged ZPE, and make ZPE correction
    if (traj%microtime>=tend .and. traj%in_cycle==1) then 

      traj%lpzpe_endtime=traj%microtime
      ! compute average kinetic energy
      do ibond=1,ctrl%lpzpe_nah
        traj%lpzpe_ke_ah(ibond)=traj%lpzpe_ke_ah(ibond)/(traj%lpzpe_endtime-traj%lpzpe_starttime)
      enddo
      do ibond=1,ctrl%lpzpe_nbc
        traj%lpzpe_ke_bc(ibond)=traj%lpzpe_ke_bc(ibond)/(traj%lpzpe_endtime-traj%lpzpe_starttime)
      enddo

      if (printlevel>4) then
        write(u_log,*) 'At the end of a cycle'
        write(u_log,*) 'Cycle start time:', tstart
        write(u_log,*) 'Cycle end time:', tend
        write(u_log,*) 'Cycle begins at:', traj%lpzpe_starttime
        write(u_log,*) 'Cycle ends at:', traj%lpzpe_endtime
        call vecwrite(ctrl%lpzpe_nah, traj%lpzpe_ke_ah, u_log, 'Average AH bonds kinetic energy', 'F14.9')
        call vecwrite(ctrl%lpzpe_nbc, traj%lpzpe_ke_bc, u_log, 'Average BC bonds kinetic energy', 'F14.9')
        call vec3write(ctrl%natom, traj%veloc_ad, u_log, 'velocity before LP-ZPE correction', 'F14.9')
        write(u_log,*) 'Start LP-ZPE correction...'
      endif

      Ekin_old=Calculate_ekin(ctrl%natom, traj%veloc_ad, traj%mass_a)

      ! compute bond energy before correction 
      do ibond=1, ctrl%lpzpe_nah
        bond_Ekin_AH_old(ibond)=0.d0
        atom1=ctrl%lpzpe_ah(ibond,1)
        atom2=ctrl%lpzpe_ah(ibond,2)
        call compute_parallel_velocity(ctrl%natom, traj%geom_ad, traj%veloc_ad, traj%mass_a,&
          & atom1, atom2, parallel_vec_tmp, parallel_veloc1_tmp, parallel_veloc2_tmp, bond_Ekin_AH_old(ibond))
      enddo
      do ibond=1, ctrl%lpzpe_nbc
        bond_Ekin_BC_old(ibond)=0.d0
        atom1=ctrl%lpzpe_bc(ibond,1)
        atom2=ctrl%lpzpe_bc(ibond,2)
        call compute_parallel_velocity(ctrl%natom, traj%geom_ad, traj%veloc_ad, traj%mass_a,&
          & atom1, atom2, parallel_vec_tmp, parallel_veloc1_tmp, parallel_veloc2_tmp, bond_Ekin_BC_old(ibond))
      enddo

      if (ctrl%lpzpe_scheme==1) then ! use the atom based scheme
        call correct_scheme1(ctrl%natom, traj%geom_ad, traj%veloc_ad, traj%mass_a, &
              &traj%element_a, ctrl%lpzpe_nbc, ctrl%lpzpe_bc, ctrl%lpzpe_nah, ctrl%lpzpe_ah, &
              &ctrl%lpzpe_ke_zpe_ah, traj%lpzpe_ke_ah, traj%lpzpe_ke_bc, ctrl%ke_threshold)             
      endif 

      Ekin_new=Calculate_ekin(ctrl%natom, traj%veloc_ad, traj%mass_a)

      ! compute bond energy after correction 
      do ibond=1, ctrl%lpzpe_nah
        bond_Ekin_AH(ibond)=0.d0
        atom1=ctrl%lpzpe_ah(ibond,1)
        atom2=ctrl%lpzpe_ah(ibond,2)
        call compute_parallel_velocity(ctrl%natom, traj%geom_ad, traj%veloc_ad, traj%mass_a,&
          & atom1, atom2, parallel_vec_tmp, parallel_veloc1_tmp, parallel_veloc2_tmp, bond_Ekin_AH(ibond))
      enddo
      do ibond=1, ctrl%lpzpe_nbc
        bond_Ekin_BC(ibond)=0.d0
        atom1=ctrl%lpzpe_bc(ibond,1)
        atom2=ctrl%lpzpe_bc(ibond,2)
        call compute_parallel_velocity(ctrl%natom, traj%geom_ad, traj%veloc_ad, traj%mass_a,&
          & atom1, atom2, parallel_vec_tmp, parallel_veloc1_tmp, parallel_veloc2_tmp, bond_Ekin_BC(ibond))
      enddo

      if (printlevel>4) then
        write(u_log,*) 'AH bond summary'
        do ibond=1, ctrl%lpzpe_nah
          atom1=ctrl%lpzpe_ah(ibond,1)
          atom2=ctrl%lpzpe_ah(ibond,2)
          write(u_log,'(A2,X,A2,X,F14.9,X,F14.9,X,F14.9)') traj%element_a(atom1), traj%element_a(atom2), bond_Ekin_AH_old(ibond), bond_Ekin_AH(ibond), bond_Ekin_AH(ibond)-bond_Ekin_AH_old(ibond)
        enddo
        write(u_log,*) 'BC bond summary'
        do ibond=1, ctrl%lpzpe_nbc
          atom1=ctrl%lpzpe_bc(ibond,1)
          atom2=ctrl%lpzpe_bc(ibond,2)
          write(u_log,'(A2,X,A2,X,F14.9,X,F14.9,X,F14.9)') traj%element_a(atom1), traj%element_a(atom2), bond_Ekin_BC_old(ibond), bond_Ekin_BC(ibond), bond_Ekin_BC(ibond)-bond_Ekin_BC_old(ibond)
        enddo
      endif

      if (printlevel>4) then
        write(u_log,*) 'LP-ZPE correction finished...'
        call vec3write(ctrl%natom, traj%veloc_ad, u_log, 'velocity after LP-ZPE correction', 'F14.9')
        write(u_log,'(A33,F14.9,A19,F14.9)') 'Kinetic energy before correction:', Ekin_old, '; after correction:', Ekin_new
      endif

      ! kick trajectory out of a cycle
      traj%in_cycle=0      
      traj%lpzpe_ke_ah=0.d0
      traj%lpzpe_ke_bc=0.d0
      traj%lpzpe_iter_incycle=1
      traj%lpzpe_cycle=traj%lpzpe_cycle+1

    endif !if (traj%microtime>=tend .and. traj%in_cycle==1) then
      

  endsubroutine

! ==========================================================
!> subroutine to compute kinetic energy parallel to AH or BC bonds
  subroutine LP_parallel_ekin(n, g, v, m, element_a, dt, nah, nbc, list_ah, list_bc, ke_ah, ke_bc)
    use definitions
    use matrix

    implicit none
    integer, intent(in) :: n
    real*8, intent(in) :: g(n,3), v(n,3), m(n)
    character*2, intent(in) :: element_a(n)
    real*8, intent(in) :: dt
    integer, intent(in) :: nah, nbc
    integer, intent(in) ::list_ah(nah,2), list_bc(nbc,2)
    real*8, intent(inout) :: ke_ah(nah), ke_bc(nbc)

    integer :: atom1, atom2
    real*8 :: r_m
    real*8 :: parallel_vec(3)
    real*8 :: parallel_veloc1, parallel_veloc2
    real*8 :: pair_veloc
    real*8 :: ke_bond_ah(nah), ke_bond_bc(nbc)
    integer :: ibond, idir

    do ibond=1,nah
      atom1=list_ah(ibond,1)
      atom2=list_ah(ibond,2)
      call compute_parallel_velocity(n, g, v, m, atom1, atom2, parallel_vec, parallel_veloc1, parallel_veloc2, ke_bond_ah(ibond))
    enddo

    do ibond=1,nbc
      atom1=list_bc(ibond,1)
      atom2=list_bc(ibond,2)
      call compute_parallel_velocity(n, g, v, m, atom1, atom2, parallel_vec, parallel_veloc1, parallel_veloc2, ke_bond_bc(ibond))
    enddo
 
    ! compute cummulative kinetic energy 
    do ibond=1,nah
      ke_ah(ibond)=ke_ah(ibond)+ke_bond_ah(ibond)*dt
    enddo
    do ibond=1,nbc
      ke_bc(ibond)=ke_bc(ibond)+ke_bond_bc(ibond)*dt
    enddo

    if (printlevel>4) then
      do ibond=1,nah
        atom1=list_ah(ibond,1)
        atom2=list_ah(ibond,2)
        write(u_log, '(A8,X,A2,X,A2,X,A15,F14.9)') 'AH Bond:',element_a(atom1), element_a(atom2),' kinetic energy:', ke_bond_ah(ibond)
      enddo
      do ibond=1,nbc
        atom1=list_bc(ibond,1)
        atom2=list_bc(ibond,2)
        write(u_log, '(A8,X,A2,X,A2,X,A15,F14.9)') 'BC Bond:',element_a(atom1), element_a(atom2),' kinetic energy:', ke_bond_bc(ibond)
      enddo
    endif

  endsubroutine

! ==========================================================
!> subroutine to compute parallel velocity
  subroutine compute_parallel_velocity(n, g, v, m, atom1, atom2, parallel_vec, parallel_veloc1, parallel_veloc2, ke)
    use definitions
    use matrix

    implicit none
    integer, intent(in) :: n
    real*8, intent(in) :: g(n,3), v(n,3), m(n)
    integer, intent(in) :: atom1, atom2
    real*8, intent(inout) :: parallel_vec(3)
    real*8, intent(inout) :: parallel_veloc1, parallel_veloc2
    real*8, intent(inout) :: ke

    real*8 :: r_m
    real*8 :: norm_parallel
    real*8 :: pair_veloc

    integer :: idir

    parallel_veloc1=0.d0
    parallel_veloc2=0.d0

    r_m=m(atom1)*m(atom2)/(m(atom1)+m(atom2))

    parallel_vec(:)=g(atom1,:)-g(atom2,:)

    norm_parallel=0.d0
    do idir=1,3
      norm_parallel=norm_parallel+parallel_vec(idir)*parallel_vec(idir)
    enddo
    norm_parallel=sqrt(norm_parallel)
    parallel_vec=parallel_vec/norm_parallel

    do idir=1,3
      parallel_veloc1=parallel_veloc1+v(atom1,idir)*parallel_vec(idir)
      parallel_veloc2=parallel_veloc2+v(atom2,idir)*parallel_vec(idir)
    enddo

    pair_veloc=parallel_veloc1-parallel_veloc2

    ke=0.5*r_m*pair_veloc*pair_veloc

  endsubroutine

! ==========================================================
!> subroutine to correct bond velocity
!> correct_bond_velocity(ctrl%natom, traj%geom_ad, traj%veloc_ad, traj%mass_a, &
!>   &atom1, atom2, eleak, factor, scheme)
  subroutine correct_bond_velocity(n, g, v, m, element_a, atom1, atom2, e, factor, r_m_ah, scheme)
    use definitions
    use matrix

    implicit none
    integer, intent(in) :: n
    real*8, intent(in) :: g(n,3), m(n)
    character*2, intent(in) :: element_a(n)
    integer, intent(in) :: atom1, atom2
    real*8, intent(in) :: e, factor, r_m_ah
    integer, intent(in) :: scheme
    real*8, intent(inout) :: v(n,3)

    real*8 :: parallel_vec(3)
    real*8 :: parallel_veloc1, parallel_veloc2
    real*8 :: pair_veloc, r_m, ke_old, ke_new
    real*8 :: deltaV
    integer :: idir

    call compute_parallel_velocity(n, g, v, m, atom1, atom2, parallel_vec, parallel_veloc1, parallel_veloc2, ke_old)

    r_m=m(atom1)*m(atom2)/(m(atom1)+m(atom2))
    pair_veloc=parallel_veloc1-parallel_veloc2    

    if (scheme==1) then 
      deltaV=sqrt(pair_veloc**2+(2*e*factor/r_m_ah))-pair_veloc
    else if (scheme==2) then
      deltaV=sqrt(pair_veloc**2-(2*e*factor/r_m))-pair_veloc
    else if (scheme==3) then
      deltaV=-pair_veloc
    endif
    
   
    do idir=1,3
      v(atom1,idir)=v(atom1,idir) + r_m/m(atom1)*deltaV*(parallel_vec(idir))
      v(atom2,idir)=v(atom2,idir) - r_m/m(atom2)*deltaV*(parallel_vec(idir))
    enddo

    call compute_parallel_velocity(n, g, v, m, atom1, atom2, parallel_vec, parallel_veloc1, parallel_veloc2, ke_new)

    if (printlevel>4) then
      if (scheme==1) then
        write(u_log,'(A23,X,A2,X,A2)') '   Single AH Bond pair:',element_a(atom1), element_a(atom2)
        write(u_log,'(A24,X,F14.9,A24,X,F14.9,A26,X,F14.9)')  '   Previous bond energy:',ke_old, '; Corrected bond energy:',ke_new,'; Energy added to AH bond:',ke_new-ke_old
      else if (scheme==2) then
        write(u_log,'(A23,X,A2,X,A2)') '   Single BC Bond pair:',element_a(atom1), element_a(atom2)
        write(u_log,'(A24,X,F14.9,A24,X,F14.9,A30,X,F14.9)')  '   Previous bond energy:',ke_old, '; Corrected bond energy:',ke_new,'; Energy removed from BC bond:',ke_old-ke_new
      else if (scheme==3) then 
        write(u_log,'(A23,X,A2,X,A2)') '   Single BC Bond pair:',element_a(atom1), element_a(atom2)
        write(u_log,'(A24,X,F14.9,A24,X,F14.9,A30,X,F14.9)')  '   Previous bond energy:',ke_old, '; Corrected bond energy:',ke_new,'; Energy removed from BC bond:',ke_old-ke_new
      endif
    endif

  endsubroutine

! ==========================================================
!>correct_scheme1(ctrl%natom, traj%geom_ad, traj%veloc_ad, traj%mass_a, &
!>  &traj%element_a, ctrl%lpzpe_nbc, ctrl%lpzpe_bc, ctrl%lpzpe_nah, ctrl%lpzpe_ah, &
!>  &ctrl%lpzpe_ke_zpe_ah, traj%lpzpe_ke_ah, traj%lpzpe_ke_bc, ctrl%ke_threshold)
  subroutine correct_scheme1(n, g, v, m, element_a, nbc, list_bc, nah, list_ah, ah_zpe, ke_ah, ke_bc, ke_threshold)
    use definitions, only: u_log, printlevel
    use matrix

    implicit none
    integer, intent(in) :: n
    real*8, intent(in) :: g(n,3), m(n)
    character*2, intent(in) :: element_a(n)
    integer, intent(in) :: nbc, list_bc(nbc,2)
    integer, intent(in) :: nah, list_ah(nah,2)
    real*8, intent(in) :: ah_zpe(nah)
    real*8, intent(in) :: ke_ah(nah), ke_bc(nbc), ke_threshold
    real*8, intent(inout) :: v(n,3)


    integer :: if_leak !0=no leak, 1=leak
    real*8 :: e_leak(nah)
    real*8 :: parallel_vec(3)
    real*8 :: parallel_veloc1
    real*8 :: parallel_veloc2

    ! used for saving all energy correction 
    real*8 :: bond_Ekin_AH(nah), bond_Ekin_AH_old(nah)
    real*8 :: e_leak_sum
    integer :: total_bonds

    integer :: atom1, atom2
    integer :: ibond, iatom, jatom, katom, idir, ipair


    ! first, check ZPE leaking, if 
    if_leak=0
    e_leak(:)=0.d0
    do ibond=1,nah
      atom1=list_ah(ibond,1)
      atom2=list_ah(ibond,2)
      e_leak(ibond)=ah_zpe(ibond)-ke_ah(ibond)
      if (e_leak(ibond)>ke_threshold) then 
        if_leak=1
        if (printlevel>4) then
          write(u_log,'(A21,I3,X,A2,X,A2,X,A8)') 'ZPE leaking for bond ', ibond, element_a(atom1), element_a(atom2), 'detected'
          write(u_log,'(A15,F14.9,X,A10,F14.9)') 'leaking energy=', e_leak(ibond), 'threshold=', ke_threshold
        endif         
      else ! set e_leak(ibond) back to 0.d0
        e_leak(ibond)=0.d0
      endif 
    enddo 
    
    if (if_leak==1) then ! Start performing correction. 

      ! compute bond energy
      bond_Ekin_AH_old(:)=0.d0 
      do ibond=1, nah
        atom1=list_ah(ibond,1)
        atom2=list_ah(ibond,2)
        call compute_parallel_velocity(n, g, v, m, atom1, atom2, &
          &parallel_vec, parallel_veloc1, parallel_veloc2, bond_Ekin_AH_old(ibond))
      enddo
 
      ! now i know what the leaking energy for each AH bond is
      !save the target AH bond kinetic energy. 
      e_leak_sum=0.d0
      bond_Ekin_AH(:)=bond_Ekin_AH_old(:)
      do ibond=1,nah
        atom1=list_ah(ibond,1)
        atom2=list_ah(ibond,2) 
        bond_Ekin_AH(ibond)=bond_Ekin_AH(ibond)+e_leak(ibond)
        e_leak_sum=e_leak_sum+e_leak(ibond)
      enddo

      total_bonds=nah+nbc
      call correct_all_bonds(n, g, v, m, element_a, total_bonds, nbc, list_bc, nah, list_ah, bond_Ekin_AH_old, bond_Ekin_AH)

    else ! no leak, exit if, and no correction is performed. 
      if (printlevel>4) then 
        write(u_log,*) 'No leakage found, skip correction'
      endif 
    endif !if (if_leak==1) then

  endsubroutine

! ==========================================================
!> subroutine to correct all atoms, in an atom by atom way
  subroutine correct_all_bonds(n, g, v, m, element_a, ntotal, nbc, list_bc, nah, list_ah, AH_old, AH_target)
    use definitions, only: u_log, printlevel
    use matrix

    implicit none
    integer, intent(in) :: n
    real*8, intent(in) :: g(n,3), m(n)
    character*2, intent(in) :: element_a(n)
    integer, intent(in) :: ntotal
    integer, intent(in) :: nbc, list_bc(nbc,2)
    integer, intent(in) :: nah, list_ah(nah,2)
    real*8, intent(in) :: AH_old(nah), AH_target(nah)
    real*8, intent(inout) :: v(n,3)

    real*8 :: rm(ntotal), rm2(ntotal)
    ! used for energy distribution 
    real*8 :: factor(nbc)
    real*8 :: Ecorrection(ntotal)

    real*8 :: parallel_vec(3)
    real*8 :: parallel_veloc1, parallel_veloc2
    real*8 :: bond_Ekin_AH(nah), bond_Ekin_BC(nbc)
    real*8 :: e_leak(nah), e_leak_sum 
    real*8 :: sum_Ekin_BC


    integer :: atom1, atom2
    integer :: iatom, jatom, ibond

    ! compute reduced energy rm
    do ibond=1,nah
      atom1=list_ah(ibond,1)
      atom2=list_ah(ibond,2)
      rm(ibond)=m(atom1)*m(atom2)/(m(atom1)+m(atom2))
      rm2(ibond)=(m(atom2)-m(atom1))/(m(atom1)*m(atom2))
    enddo
    do ibond=1,nbc
      atom1=list_bc(ibond,1)
      atom2=list_bc(ibond,2)
      rm(ibond+nah)=m(atom1)*m(atom2)/(m(atom1)+m(atom2))
      rm2(ibond+nah)=(m(atom2)-m(atom1))/(m(atom1)*m(atom2))
    enddo

    e_leak_sum=0.d0
    do ibond=1,nah
      e_leak(ibond)=AH_target(ibond)-AH_old(ibond)
      e_leak_sum=e_leak_sum+e_leak(ibond)
    enddo

    ! compute Ecorrection. 
    ! notice parallel_veloc for bond AB, is defined as v_A dot u_AB
    ! v_B dot u_AB = -v_B dot u_BA
    do ibond=1, nah
      atom1=list_ah(ibond,1)
      atom2=list_ah(ibond,2)
      call compute_parallel_velocity(n, g, v, m, atom1, atom2, &
        &parallel_vec, parallel_veloc1, parallel_veloc2, bond_Ekin_AH(ibond))
    enddo
    do ibond=1, nbc
      atom1=list_bc(ibond,1)
      atom2=list_bc(ibond,2)
      call compute_parallel_velocity(n, g, v, m, atom1, atom2, &
        &parallel_vec, parallel_veloc1, parallel_veloc2, bond_Ekin_BC(ibond))
    enddo
    ! start distribute energy to BC bonds. 
    sum_Ekin_BC=0.d0
    do ibond=1,nbc
      sum_Ekin_BC=sum_Ekin_bc+bond_Ekin_BC(ibond)
    enddo
    do ibond=1,nbc
      factor(ibond)=bond_Ekin_BC(ibond)/sum_Ekin_BC
    enddo
    do ibond=1,nah
      atom1=list_ah(ibond,1)
      atom2=list_ah(ibond,2)
      Ecorrection(ibond)=e_leak(ibond)
    enddo 
    do ibond=1,nbc
      atom1=list_bc(ibond,1)
      atom2=list_bc(ibond,2)
      Ecorrection(ibond+nah)=-factor(ibond)*e_leak_sum
    enddo 

    if (printlevel>4) then
      call vecwrite(ntotal,Ecorrection,u_log,'Ecorrection','F12.9')
    endif 


    call correct_bonds(n, g, v, m, element_a, ntotal, nbc, list_bc, &
      &nah, list_ah, rm, rm2, Ecorrection) 


  endsubroutine

! ==========================================================
!> subroutine to correct single atom velocity
  subroutine correct_bonds(n, g, v, m, element_a, ntotal, nbc, list_bc, nah, list_ah, rm, rm2, Ecorrection)
    use definitions, only: u_log, printlevel
    use matrix

    implicit none
    integer, intent(in) :: n
    real*8, intent(in) :: g(n,3), m(n)
    character*2, intent(in) :: element_a(n)
    integer, intent(in) :: ntotal
    integer, intent(in) :: nbc, list_bc(nbc,2)
    integer, intent(in) :: nah, list_ah(nah,2) 
    real*8, intent(in) :: rm(ntotal), rm2(ntotal)
    real*8, intent(inout) :: Ecorrection(ntotal)
    real*8, intent(inout) :: v(n,3)

    real*8 :: parallel_vec(ntotal,3)
    real*8 :: parallel_veloc1(ntotal)
    real*8 :: parallel_veloc2(ntotal)
    real*8 :: bond_Ekin_AH(nah), bond_Ekin_BC(nbc)
    real*8 :: bond_Ekin_AH_new(nah), bond_Ekin_BC_new(nbc)

    real*8 :: dveloc(ntotal), dEvec(ntotal)
    real*8 :: bvec(ntotal), cvec(ntotal), umat(ntotal,ntotal)
    integer :: scratch_s(ntotal),info
    integer :: atom_in_bonds(n)

    integer :: atom1, atom2
    integer :: ibond, jbond, idir, iatom, jatom, katom

    ! compute bond energy 
    parallel_vec(:,:)=0.d0
    parallel_veloc1(:)=0.d0
    parallel_veloc2(:)=0.d0
    do ibond=1, nah
      atom1=list_ah(ibond,1)
      atom2=list_ah(ibond,2)
      call compute_parallel_velocity(n, g, v, m, atom1, atom2, &
        &parallel_vec(ibond,:), parallel_veloc1(ibond), &
        &parallel_veloc2(ibond), bond_Ekin_AH(ibond))
    enddo
    do ibond=1, nbc
      atom1=list_bc(ibond,1)
      atom2=list_bc(ibond,2)
      call compute_parallel_velocity(n, g, v, m, atom1, atom2, &
        &parallel_vec(ibond+nah,:), parallel_veloc1(ibond+nah), &
        &parallel_veloc2(ibond+nah), bond_Ekin_BC(ibond))
    enddo
    write(u_log,*) 'parallel_vec'
    do ibond=1,ntotal
      write(u_log,*) ibond
      write(u_log,*) parallel_vec(ibond,1),parallel_vec(ibond,2),parallel_vec(ibond,3)
      write(u_log,*) parallel_veloc1(ibond), parallel_veloc2(ibond)
    enddo

    ! now compute everything between bond pairs
    write(u_log,*) 'in ALPHA-BETA'
    do ibond=1,ntotal
      dveloc(ibond)=parallel_veloc1(ibond)-parallel_veloc2(ibond)
      dEvec(ibond)=dveloc(ibond)**2+2*Ecorrection(ibond)/rm(ibond)
      write(u_log,*) 'ibond',ibond
      write(u_log,*) 'rm2(ibond)', rm2(ibond)
      write(u_log,*) 'Ecorrection',Ecorrection(ibond)
      write(u_log,*) 'dveloc(ibond), dEvec(ibond)', dveloc(ibond), dEvec(ibond)
      if (dEvec(ibond)<0) then 
        Ecorrection(ibond)=0.5*dveloc(ibond)**2*rm(ibond)
        dEvec(ibond)=0.d0
      endif 
      bvec(ibond)=(-dveloc(ibond)+sqrt(dEvec(ibond)))/rm2(ibond)
      write(u_log,*) 'bvec(ibond)',bvec(ibond)
      do jbond=1,ntotal
        umat(ibond,jbond)=0.d0
        do idir=1,3
          umat(ibond,jbond)=umat(ibond,jbond)+parallel_vec(ibond,idir)*parallel_vec(jbond,idir)
        enddo
      enddo
    enddo
    cvec(:)=0.d0
    call matwrite(ntotal,umat,u_log,'umat','F20.14')
    call vecwrite(ntotal,bvec,u_log,'bvec','F20.14')

    ! solving cvect, umat * cvec = bvec. 
    call dgesv(ntotal,1,umat,ntotal,scratch_s,bvec,ntotal,info)
    cvec(:)=bvec(:)
   
    call vecwrite(ntotal,cvec,u_log,'cvec','F50.25')

    ! check if an atom is in AH or BC bonds
    atom_in_bonds(:)=0
    do ibond=1,nah
      atom1=list_ah(ibond,1)
      atom2=list_ah(ibond,2)
      atom_in_bonds(atom1)=1
      atom_in_bonds(atom2)=1
    enddo
    do ibond=1,nbc
      atom1=list_bc(ibond,1)
      atom2=list_bc(ibond,2)
      atom_in_bonds(atom1)=1
      atom_in_bonds(atom2)=1
    enddo

    ! update velocity 
    call vec3write(n, v, u_log, 'before!! velocity', 'F50.25')

    if (printlevel>4) then
      call vecwrite(ntotal,Ecorrection,u_log,'Ecorrection','F12.9')
    endif

    do iatom=1,n
      do ibond=1,nah
        v(iatom,:)=v(iatom,:)+cvec(ibond)*parallel_vec(ibond,:)/m(iatom)
      enddo
      do ibond=1,nbc
        v(iatom,:)=v(iatom,:)+cvec(ibond+nah)*parallel_vec(ibond+nah,:)/m(iatom)
      enddo
    enddo

    !do ibond=1,nah
    !  atom1=list_ah(ibond,1)
    !  atom2=list_ah(ibond,2)
    !  do jbond=1,ntotal
    !    v(atom1,:)=v(atom1,:)+cvec(ibond)*parallel_vec(ibond,:)/m(atom1)
    !    v(atom2,:)=v(atom2,:)+cvec(ibond)*parallel_vec(ibond,:)/m(atom2)
    !  enddo
    !enddo 
    !do ibond=1,nbc
    !  atom1=list_bc(ibond,1)
    !  atom2=list_bc(ibond,2)
    !  do jbond=1,ntotal
    !    v(atom1,:)=v(atom1,:)+cvec(ibond+nah)*parallel_vec(ibond+nah,:)/m(atom1)
    !    v(atom2,:)=v(atom2,:)+cvec(ibond+nah)*parallel_vec(ibond+nah,:)/m(atom2)
    !  enddo
    !enddo

    call vec3write(n, v, u_log, 'after!! velocity', 'F50.25') 

    ! compute bond energy 
    parallel_vec(:,:)=0.d0
    parallel_veloc1(:)=0.d0
    parallel_veloc2(:)=0.d0
    do ibond=1, nah
      atom1=list_ah(ibond,1)
      atom2=list_ah(ibond,2)
      call compute_parallel_velocity(n, g, v, m, atom1, atom2, &
        &parallel_vec(ibond,:), parallel_veloc1(ibond), &
        &parallel_veloc2(ibond), bond_Ekin_AH_new(ibond))
    enddo
    do ibond=1, nbc
      atom1=list_bc(ibond,1)
      atom2=list_bc(ibond,2)
      call compute_parallel_velocity(n, g, v, m, atom1, atom2, &
        &parallel_vec(ibond+nah,:), parallel_veloc1(ibond+nah), &
        &parallel_veloc2(ibond+nah), bond_Ekin_BC_new(ibond))
    enddo

      write(u_log,*) 'AH'
      do ibond=1, nah
      atom1=list_ah(ibond,1)
      atom2=list_ah(ibond,2)
      write(u_log,*) ibond, atom1, atom2 
      write(u_log,*) bond_Ekin_AH(ibond), bond_Ekin_AH_new(ibond), bond_Ekin_AH_new(ibond)-bond_Ekin_AH(ibond)
      enddo
      do ibond=1, nbc
      atom1=list_bc(ibond,1)
      atom2=list_bc(ibond,2)
      write(u_log,*) ibond, atom1, atom2
      write(u_log,*) bond_Ekin_BC(ibond), bond_Ekin_BC_new(ibond), bond_Ekin_BC_new(ibond)-bond_Ekin_BC(ibond)
      enddo


  endsubroutine


! ==========================================================
!> subroutine to correct BC bonds velocity - method1
  subroutine correct_bcbond_velocity_method1(n, g, v, m, element_a, r_m_ah, e, nbc, list_bc, do_correction) 
    use definitions
    use matrix

    implicit none
    integer, intent(in) :: n
    real*8, intent(in) :: g(n,3), m(n)
    character*2, intent(in) :: element_a(n)
    real*8, intent(in) :: r_m_ah, e
    integer, intent(in) :: nbc, list_bc(nbc,2)
    real*8, intent(inout) :: v(n,3)
    integer, intent(inout) :: do_correction

    real*8 :: vb(nbc), vc(nbc)
    real*8 :: parallel_vec(nbc,3), ke(nbc)
    real*8 :: sum_ke, sum_ke_new
    integer :: bc_valid(nbc)
    integer :: rescale

    real*8 :: factor(nbc)
    real*8 :: deltaE(nbc)
    integer :: scheme

    integer :: atom1, atom2
    integer :: iter, maxiter
    integer :: remaining_nbc

    integer :: ibond

    maxiter=nbc
    iter=1  
    remaining_nbc=nbc

    ! initialize all BC bonds to be valid in computing factor
    bc_valid(:)=1

    do ibond=1,nbc
      atom1=list_bc(ibond,1)
      atom2=list_bc(ibond,2)
      call compute_parallel_velocity(n, g, v, m, atom1, atom2, parallel_vec(ibond,:), vb(ibond), vc(ibond), ke(ibond))
      write(u_log,*) 'original ke', ibond, ke(ibond)
    enddo

    rescale=1 ! so that we start with the do while loop
    do while (iter.le.maxiter .and. rescale.eq.1)
      rescale=0
      sum_ke=0.d0
      factor(:)=0.d0
 
      do ibond=1,nbc
        if (bc_valid(ibond)==1) then
          sum_ke=sum_ke+ke(ibond)
        endif
      enddo
 
      do ibond=1,nbc
        if (bc_valid(ibond)==1) then
          factor(ibond)=ke(ibond)/sum_ke
          deltaE(ibond)=(vb(ibond)-vc(ibond))**2-(2*e*factor(ibond)/r_m_ah)
          if (deltaE(ibond)<0) then
            bc_valid(ibond)=0
            rescale=1
            remaining_nbc=remaining_nbc-1
          endif
        endif
      enddo

      iter=iter+1
    enddo

    if (remaining_nbc==0) then
      do_correction=0
    else if (remaining_nbc>0) then
      ! correct BC velocity
      scheme=2    
      do ibond=1,nbc
        atom1=list_bc(ibond,1)
        atom2=list_bc(ibond,2)
        if (bc_valid(ibond)==1) then
          call correct_bond_velocity(n, g, v, m, element_a, atom1, atom2, e, factor(ibond), r_m_ah, scheme)
        endif
      enddo
      ! compute new kinetic energy form BC bonds
      sum_ke_new=0.d0
      do ibond=1,nbc
        atom1=list_bc(ibond,1)
        atom2=list_bc(ibond,2)
        call compute_parallel_velocity(n, g, v, m, atom1, atom2, parallel_vec(ibond,:), vb(ibond), vc(ibond), ke(ibond))
        sum_ke_new=sum_ke_new+ke(ibond)
        write(u_log,*) 'updated ke', ibond, ke(ibond)
      enddo
    endif !if (remaining_nbc==0) then

    if (printlevel>4) then
      write(u_log,*) 'All BC bond velocity correction finished. Summary for all BC bonds'
      write(u_log,'(X,A20,X,F14.9,A28,X,F14.9,A31,X,F14.9)')  'Old BC bonds energy:',sum_ke, '; Corrected BC bonds energy:',sum_ke_new,'; Energy removed from BC bonds:', sum_ke-sum_ke_new
    endif

  endsubroutine

! ==========================================================
!> subroutine to correct BC bonds velocity - method2
  subroutine correct_bcbond_velocity_method2(n, g, v, m, element_a, r_m_ah, e, nbc, list_bc, eavail)
    use definitions
    use matrix

    implicit none
    integer, intent(in) :: n
    real*8, intent(in) :: g(n,3), m(n)
    character*2, intent(in) :: element_a(n)
    real*8, intent(in) :: r_m_ah, e
    integer, intent(in) :: nbc, list_bc(nbc,2)
    real*8, intent(inout) :: v(n,3), eavail

    real*8 :: vb(nbc), vc(nbc)
    real*8 :: parallel_vec(nbc,3), ke(nbc), ke_updated(nbc)
    real*8 :: e_updated, e_avail(nbc)
    real*8 :: sum_ke, sum_ke_new
    real*8 :: sum_e
    real*8 :: factor(nbc)
    real*8 :: deltaE(nbc)
    real*8 :: r_m
    integer :: scheme
    real*8 :: deltaE_tmp

    integer :: atom1, atom2
    integer :: ibond, jbond

    sum_ke=0.d0
    do ibond=1,nbc
      atom1=list_bc(ibond,1)
      atom2=list_bc(ibond,2)
      call compute_parallel_velocity(n, g, v, m, atom1, atom2, parallel_vec(ibond,:), vb(ibond), vc(ibond), ke(ibond))
      sum_ke=sum_ke+ke(ibond)
    enddo
    
    do ibond=1,nbc
      atom1=list_bc(ibond,1)
      atom2=list_bc(ibond,2)
      factor(ibond)=ke(ibond)/sum_ke
      r_m=m(atom1)*m(atom2)/(m(atom1)+m(atom2))
      ! update vb and vc, because vb and vc will change after the last BC bond being corrected 
      call compute_parallel_velocity(n, g, v, m, atom1, atom2, parallel_vec(ibond,:), vb(ibond), vc(ibond), ke_updated(ibond))
      ! the energy needs to be updated, because the correction on first BC bond may affect the kinetic energy of second BC bond
      e_updated=(e*factor(ibond)-(ke(ibond)-ke_updated(ibond)))/factor(ibond)
      deltaE(ibond)=(vb(ibond)-vc(ibond))**2-(2*e_updated*factor(ibond)/r_m)
      if (printlevel>4) then 
        write(u_log,*) 'Updated correction energy:', e_updated
      endif
      if (deltaE(ibond).ge.0) then
        scheme=2
        call correct_bond_velocity(n, g, v, m, element_a, atom1, atom2, e_updated, factor(ibond), r_m_ah, scheme)
      else if (deltaE(ibond).lt.0) then
        ! recompute e_avail(ibond)
        e_avail(ibond)=((vb(ibond)-vc(ibond))**2)*r_m/2/factor(ibond)
        scheme=3
        call correct_bond_velocity(n, g, v, m, element_a, atom1, atom2, e_avail(ibond), factor(ibond), r_m_ah, scheme)
      endif 
    enddo

    ! compute new kinetic energy form BC bonds
    sum_ke_new=0.d0
    do ibond=1,nbc
      atom1=list_bc(ibond,1)
      atom2=list_bc(ibond,2)
      call compute_parallel_velocity(n, g, v, m, atom1, atom2, parallel_vec(ibond,:), vb(ibond), vc(ibond), ke(ibond))
      sum_ke_new=sum_ke_new+ke(ibond)
    enddo

    eavail=sum_ke-sum_ke_new

    if (printlevel>4) then
      write(u_log,*) 'All BC bonds velocity correction finished. Summary for all BC bonds'
      write(u_log,'(A20,X,F14.9,A28,X,F14.9,A31,X,F14.9)')  'Old BC bonds energy:',sum_ke, '; Corrected BC bonds energy:',sum_ke_new,'; Energy removed from BC bonds:', sum_ke-sum_ke_new
    endif

  endsubroutine

! ==========================================================
!> subroutine to correct BC bonds velocity - method3
  subroutine correct_bcbond_velocity_method3(n, g, v, m, element_a, r_m_ah, e, nbc, list_bc, do_correction)
    use definitions
    use matrix

    implicit none
    integer, intent(in) :: n
    real*8, intent(in) :: g(n,3), m(n)
    character*2, intent(in) :: element_a(n)
    real*8, intent(in) :: r_m_ah, e
    integer, intent(in) :: nbc, list_bc(nbc,2)
    real*8, intent(inout) :: v(n,3)
    integer, intent(inout) :: do_correction

    real*8 :: vb(nbc), vc(nbc)
    real*8 :: parallel_vec(nbc,3), ke(nbc), ke_updated(nbc)
    real*8 :: sum_ke, sum_ke_old, sum_ke_new
    integer :: bc_valid(nbc)
    integer :: rescale

    real*8 :: e_updated
    real*8 :: factor(nbc)
    real*8 :: deltaE(nbc)
    integer :: scheme

    integer :: atom1, atom2
    integer :: iter, maxiter
    integer :: remaining_nbc, remaining_total
    real*8 :: factor_left, energy_left
    real*8 :: v_tmp(n,3)

    integer :: ibond, jbond

    sum_ke_old=0.d0
    do ibond=1,nbc
      atom1=list_bc(ibond,1)
      atom2=list_bc(ibond,2)
      call compute_parallel_velocity(n, g, v, m, atom1, atom2, parallel_vec(ibond,:), vb(ibond), vc(ibond), ke(ibond))
      sum_ke_old=sum_ke_old+ke(ibond)
    enddo

    v_tmp(:,:)=v(:,:)

    remaining_total=nbc
    energy_left=e
    factor_left=1.d0

    do ibond=1,nbc

      do jbond=1,nbc
        bc_valid(jbond)=1
        atom1=list_bc(jbond,1)
        atom2=list_bc(jbond,2)
        call compute_parallel_velocity(n, g, v_tmp, m, atom1, atom2, parallel_vec(jbond,:), vb(jbond), vc(jbond), ke(jbond))
      enddo

      maxiter=remaining_total
      remaining_nbc=remaining_total
      iter=1
      rescale=1 ! so that we start with the do while loop
      do while (iter.le.maxiter .and. rescale.eq.1)
        rescale=0
        sum_ke=0.d0
       
        do jbond=ibond,nbc
          factor(jbond)=0.d0
          if (bc_valid(jbond)==1) then
            sum_ke=sum_ke+ke(jbond)
          endif
        enddo

        do jbond=ibond,nbc
          if (bc_valid(jbond)==1) then
            factor(jbond)=ke(jbond)/sum_ke
            deltaE(jbond)=(vb(jbond)-vc(jbond))**2-(2*energy_left*factor(jbond)/r_m_ah)
            if (deltaE(jbond)<0) then
              bc_valid(jbond)=0
              rescale=1
              remaining_nbc=remaining_nbc-1
            endif
          endif
        enddo

        iter=iter+1
      enddo

      if (remaining_nbc==0) then
        do_correction=0
      else if (remaining_nbc>0) then ! correct the current BC bond
        scheme=2
        atom1=list_bc(ibond,1)
        atom2=list_bc(ibond,2)
        call correct_bond_velocity(n, g, v_tmp, m, element_a, atom1, atom2, energy_left, factor(ibond), r_m_ah, scheme)
      endif

      ! heading to the next BC bond
      remaining_total=remaining_total-1
      energy_left=energy_left-energy_left*factor(ibond)
      factor_left=factor_left-factor(ibond)

    enddo ! do ibond=1,nbc

    if (do_correction==1) then 
      v(:,:)=v_tmp(:,:)
    endif 

    sum_ke_new=0.d0
    do ibond=1,nbc
      atom1=list_bc(ibond,1)
      atom2=list_bc(ibond,2)
      call compute_parallel_velocity(n, g, v, m, atom1, atom2, parallel_vec(ibond,:), vb(ibond), vc(ibond), ke(ibond))
      sum_ke_new=sum_ke_new+ke(ibond)
    enddo

    if (printlevel>4) then
      write(u_log,*) 'All BC bonds velocity correction finished. Summary for all BC bonds'
      write(u_log,'(A20,X,F14.9,A28,X,F14.9,A31,X,F14.9)')  'Old BC bonds energy:',sum_ke_old, '; Corrected BC bonds energy:',sum_ke_new,'; Energy removed from BC bonds:', sum_ke_old-sum_ke_new
    endif


  endsubroutine

endmodule 
