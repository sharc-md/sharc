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

!> # Module NUCLEAR
!> 
!> \author Sebastian Mai
!> \ðate 27.02.2015
!>
!>
!>                   modified 11.13.2019 by Yinan Shu
!>                       change subroutine Calculate_etot - works for SCP now
!>                       add following subroutines used for adaptive step size:
!>                       Check_Consistency, Adaptive_Stepsize, Back_propagate
!>
!> This module defines all subroutines for the nuclear dynamics:
!> - the two velocity-verlet steps (update of geometry, update of velocity)
!> - calculation of total and kinetic energy
!> - rescaling of velocities after surface hop
!> - damping of velocities
module nuclear
 contains

! ===========================================================

!> performs the geometry update of the Velocity Verlet algorithm
!> a(t)=g(t)/M
!> x(t+dt)=x(t)+v(t)*dt+0.5*a(t)*dt^2
subroutine VelocityVerlet_xstep(traj,ctrl)
  use definitions
  use matrix
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl
  integer :: iatom, idir

  if (printlevel>2) then
    write(u_log,*) '============================================================='
    write(u_log,*) '              Velocity Verlet -- X-step'
    write(u_log,*) '============================================================='
    call vec3write(ctrl%natom,traj%accel_ad,u_log,'Old accel','F12.9')
    call vec3write(ctrl%natom,traj%geom_ad,u_log,'Old geom','F12.7')
  endif

  do iatom=1,ctrl%natom
    do idir=1,3
      traj%accel_ad(iatom,idir)=&
      &-traj%grad_ad(iatom,idir)/traj%mass_a(iatom)

      traj%geom_ad(iatom,idir)=&
      & traj%geom_ad(iatom,idir)&
      &+traj%veloc_ad(iatom,idir)*ctrl%dtstep&
      &+0.5d0*traj%accel_ad(iatom,idir)*ctrl%dtstep**2
    enddo
  enddo

  if (printlevel>2) then
    call vec3write(ctrl%natom,traj%accel_ad,u_log,'accel','F12.9')
    call vec3write(ctrl%natom,traj%geom_ad,u_log,'geom','F12.7')
  endif

endsubroutine

! ===========================================================

!> performs the velocity update of the Velocity Verlet algorithm
!> a(t+dt)=g(t+dt)/M
!> v(t+dt)=v(t)+a(t+dt)*dt
subroutine VelocityVerlet_vstep(traj,ctrl)
  use definitions
  use matrix
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl
  integer :: iatom, idir

  if (printlevel>2) then
    write(u_log,*) '============================================================='
    write(u_log,*) '              Velocity Verlet -- V-step'
    write(u_log,*) '============================================================='
    call vec3write(ctrl%natom,traj%accel_ad,u_log,'Old accel','F12.9')
    call vec3write(ctrl%natom,traj%veloc_ad,u_log,'Old veloc','F12.9')
  endif

  do iatom=1,ctrl%natom
    do idir=1,3
      traj%accel_ad(iatom,idir)=0.5d0*(traj%accel_ad(iatom,idir)&
      &-traj%grad_ad(iatom,idir)/traj%mass_a(iatom) )

      traj%veloc_ad(iatom,idir)=&
      & traj%veloc_ad(iatom,idir)&
      &+traj%accel_ad(iatom,idir)*ctrl%dtstep
    enddo
  enddo

  if (printlevel>2) then
    call vec3write(ctrl%natom,traj%accel_ad,u_log,'accel','F12.9')
    call vec3write(ctrl%natom,traj%veloc_ad,u_log,'veloc','F12.9')
  endif

endsubroutine

! ===========================================================

!> performs the velocity update of the Velocity Verlet algorithm
!> a(t+dt)=g(t+dt)/M
!> v(t+dt)=v(t)+a(t+dt)*dt
subroutine VelocityVerlet_vstep_approximate(vold, v, a, dt, n)
  use definitions
  use matrix
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: vold(n,3), a(n,3)
  real*8, intent(out) :: v(n,3)
  real*8, intent(in) :: dt
  integer :: iatom, idir


  do iatom=1,n
    do idir=1,3
      v(iatom,idir)=vold(iatom,idir)+a(iatom,idir)*dt
    enddo
  enddo

endsubroutine

! ===========================================================

!> calculates the sum of the kinetic energies of all atoms
real*8 function Calculate_ekin(n, veloc, mass) result(Ekin)
  implicit none
  integer, intent(in) :: n
  real*8,intent(in) :: veloc(n,3), mass(n)
  integer :: i

  Ekin=0.d0
  do i=1,n
    ! sum of square can be written as sum(a**2)
    Ekin=Ekin + 0.5d0*mass(i)*sum(veloc(i,:)**2)
  enddo


endfunction

! ===========================================================

!> calculates the sum of the kinetic energies of all atoms
real*8 function Calculate_ekin_masked(n, veloc, mass, mask) result(Ekin)
  implicit none
  integer, intent(in) :: n
  real*8,intent(in) :: veloc(n,3), mass(n)
  logical,intent(in) :: mask(n)
  integer :: i

  Ekin=0.d0
  do i=1,n
    ! sum of square can be written as sum(a**2)
    if (mask(i))  Ekin=Ekin + 0.5d0*mass(i)*sum(veloc(i,:)**2)
  enddo


endfunction

! ===========================================================

!> calculates the kinetic energy,
!> the potential energy (from H_diag)
!> and the total energy
subroutine Calculate_etot(traj,ctrl)
  use definitions
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

  integer :: istate, jstate
  complex*16 :: den_ss(ctrl%nstates,ctrl%nstates)

! computde density matrix
  do istate=1,ctrl%nstates
    do jstate=1,ctrl%nstates
      den_ss(istate,jstate)=traj%coeff_diag_s(istate)*conjg(traj%coeff_diag_s(jstate))
    enddo
  enddo

! save the previous step energies
  traj%Ekin_old=traj%Ekin
  traj%Epot_old=traj%Epot
  traj%Etot_old=traj%Etot

! depend on either using TSH or SCP
  if(ctrl%method==0) then !TSH
    traj%Ekin=Calculate_ekin(ctrl%natom, traj%veloc_ad, traj%mass_a)
    traj%Epot=real(traj%H_diag_ss(traj%state_diag,traj%state_diag))
    traj%Etot=traj%Ekin+traj%Epot
  else if (ctrl%method==1) then !SCP
    traj%Ekin=Calculate_ekin(ctrl%natom, traj%veloc_ad, traj%mass_a)
    traj%Epot=0.d0
    do istate=1,ctrl%nstates
      traj%Epot=traj%Epot+&
      &real(den_ss(istate,istate)*traj%H_diag_ss(istate,istate))
    enddo
    traj%Etot=traj%Ekin+traj%Epot
  endif

  if (printlevel>2) then
    write(u_log,*) '============================================================='
    write(u_log,*) '                        Energies'
    write(u_log,*) '============================================================='
    write(u_log,'(A,1X,F20.10,1X,A)') 'Ekin:',traj%Ekin*au2ev,'eV'
    write(u_log,'(A,1X,F20.10,1X,A)') 'Epot:',traj%Epot*au2ev,'eV'
    write(u_log,'(A,1X,F20.10,1X,A)') 'Etot:',traj%Etot*au2ev,'eV'
    write(u_log,*) ''
  endif

endsubroutine

! ===========================================================

!> Rescales the velocity after each timestep
!> - no rescaling if no hop occured
!> - no rescaling after field-induced hops
!> - otherwise, rescaling according to input options:
!>    * no rescaling
!>    * parallel to velocity vector
!>    * parallel to relevant non-adiabatic coupling vector
!> - reflection after frustrated hops according to input options
subroutine Rescale_velocities(traj,ctrl)
  use definitions
  use matrix
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl
  real*8 :: factor, sum_kk, sum_vk, deltaE, Ekin_masked, Ekin_new
  integer :: i

  if (printlevel>2) then
    write(u_log,*) '============================================================='
    write(u_log,*) '                     Velocity Rescaling'
    write(u_log,*) '============================================================='
  endif
  select case (traj%kind_of_jump)
    case (0)
      if (printlevel>2) write(u_log,'(A)') 'No jump occured.'
    case (1,4)
      select case (ctrl%ekincorrect)
        case (0)
          if (printlevel>2) write(u_log,*) 'Velocity is not rescaled after surface hop.'
        case (1)
          Ekin_masked=Calculate_ekin_masked(ctrl%natom, traj%veloc_ad, traj%mass_a, ctrl%atommask_a)
          ! Use Epot+Ekin
!           Ekin_new=Ekin_masked &
!           & + real(traj%H_diag_ss(traj%state_diag_old,traj%state_diag_old)) &
!           & - real(traj%H_diag_ss(traj%state_diag,traj%state_diag))
          ! Use Etot
          Ekin_new=traj%Etot&
          & - traj%Ekin + Ekin_masked &
          & - real(traj%H_diag_ss(traj%state_diag,traj%state_diag))
          factor=sqrt( Ekin_new/Ekin_masked )
          do i=1,ctrl%natom
            if (ctrl%atommask_a(i)) traj%veloc_ad(i,:)=traj%veloc_ad(i,:)*factor
          enddo
          if (printlevel>2) then
            write(u_log,'(A)') 'Velocity is rescaled along velocity vector.'
            if (.not.all(ctrl%atommask_a)) write(u_log,'(A)') 'Some velocities were not rescaled due to mask.'
            write(u_log,'(A,1X,F12.9)') 'Scaling factor is ',factor
          endif
! ! !           factor=sqrt( (traj%Etot-real(traj%H_diag_ss(traj%state_diag,traj%state_diag)))/traj%Ekin )
! ! !           traj%veloc_ad=traj%veloc_ad*factor
! ! !           if (printlevel>2) then
! ! !             write(u_log,'(A)') 'Velocity is rescaled along velocity vector.'
! ! !             write(u_log,'(A,1X,F12.9)') 'Scaling factor is ',factor
! ! !           endif
        case (2)
          call available_ekin(ctrl%natom,&
          &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag_old,traj%state_diag,:,:),&
          &traj%mass_a, sum_kk, sum_vk)
          deltaE=4.d0*sum_kk*(traj%Etot-traj%Ekin-&
          &real(traj%H_diag_ss(traj%state_diag,traj%state_diag)))+sum_vk**2
          if (sum_vk<0.d0) then
            factor=(sum_vk+sqrt(deltaE))/2.d0/sum_kk
          else
            factor=(sum_vk-sqrt(deltaE))/2.d0/sum_kk
          endif
          do i=1,3
            traj%veloc_ad(:,i)=traj%veloc_ad(:,i)-factor*&
            &traj%hopping_direction_ssad(traj%state_diag_old,traj%state_diag,:,i)/traj%mass_a(:)
          enddo
          if (printlevel>2) then
            write(u_log,'(A)') 'Velocity is rescaled along projected velocity vector.'
            write(u_log,'(A,1X,E16.8,1X,E16.8)') 'a, b: ', sum_kk, sum_vk
            write(u_log,'(A,1X,E16.8)') 'Delta is          ',deltaE
            write(u_log,'(A,1X,F12.6)') 'Scaling factor is ',factor
          endif
        case (3)
          call available_ekin(ctrl%natom,&
          &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag_old,traj%state_diag,:,:),&
          &traj%mass_a, sum_kk, sum_vk)
          deltaE=4.d0*sum_kk*(traj%Etot-traj%Ekin-&
          &real(traj%H_diag_ss(traj%state_diag,traj%state_diag)))+sum_vk**2
          if (sum_vk<0.d0) then
            factor=(sum_vk+sqrt(deltaE))/2.d0/sum_kk
          else
            factor=(sum_vk-sqrt(deltaE))/2.d0/sum_kk
          endif
          do i=1,3
            traj%veloc_ad(:,i)=traj%veloc_ad(:,i)-factor*&
            &traj%hopping_direction_ssad(traj%state_diag_old,traj%state_diag,:,i)/traj%mass_a(:)
          enddo
          if (printlevel>2) then
            write(u_log,'(A)') 'Velocity is rescaled along non-adiabatic coupling vector.'
            write(u_log,'(A,1X,E16.8,1X,E16.8)') 'a, b: ', sum_kk, sum_vk
            write(u_log,'(A,1X,E16.8)') 'Delta is          ',deltaE
            write(u_log,'(A,1X,F12.6)') 'Scaling factor is ',factor
          endif
        case (4)
          call available_ekin(ctrl%natom,&
          &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag_old,traj%state_diag,:,:),&
          &traj%mass_a, sum_kk, sum_vk)
          deltaE=4.d0*sum_kk*(traj%Etot-traj%Ekin-&
          &real(traj%H_diag_ss(traj%state_diag,traj%state_diag)))+sum_vk**2
          if (sum_vk<0.d0) then
            factor=(sum_vk+sqrt(deltaE))/2.d0/sum_kk
          else
            factor=(sum_vk-sqrt(deltaE))/2.d0/sum_kk
          endif
          do i=1,3
            traj%veloc_ad(:,i)=traj%veloc_ad(:,i)-factor*&
            &traj%hopping_direction_ssad(traj%state_diag_old,traj%state_diag,:,i)/traj%mass_a(:)
          enddo
          if (printlevel>2) then
            write(u_log,'(A)') 'Velocity is rescaled along gradient difference vector.'
            write(u_log,'(A,1X,E16.8,1X,E16.8)') 'a, b: ', sum_kk, sum_vk
            write(u_log,'(A,1X,E16.8)') 'Delta is          ',deltaE
            write(u_log,'(A,1X,F12.6)') 'Scaling factor is ',factor
          endif
        case (5)
          call available_ekin(ctrl%natom,&
          &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag_old,traj%state_diag,:,:),&
          &traj%mass_a, sum_kk, sum_vk)
          deltaE=4.d0*sum_kk*(traj%Etot-traj%Ekin-&
          &real(traj%H_diag_ss(traj%state_diag,traj%state_diag)))+sum_vk**2
          if (sum_vk<0.d0) then
            factor=(sum_vk+sqrt(deltaE))/2.d0/sum_kk
          else
            factor=(sum_vk-sqrt(deltaE))/2.d0/sum_kk
          endif
          do i=1,3
            traj%veloc_ad(:,i)=traj%veloc_ad(:,i)-factor*&
            &traj%hopping_direction_ssad(traj%state_diag_old,traj%state_diag,:,i)/traj%mass_a(:)
          enddo
          if (printlevel>2) then
            write(u_log,'(A)') 'Velocity is rescaled along projected non-adiabatic coupling vector.'
            write(u_log,'(A,1X,E16.8,1X,E16.8)') 'a, b: ', sum_kk, sum_vk
            write(u_log,'(A,1X,E16.8)') 'Delta is          ',deltaE
            write(u_log,'(A,1X,F12.6)') 'Scaling factor is ',factor
          endif
        case (6)
          call available_ekin(ctrl%natom,&
          &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag_old,traj%state_diag,:,:),&
          &traj%mass_a, sum_kk, sum_vk)
          deltaE=4.d0*sum_kk*(traj%Etot-traj%Ekin-&
          &real(traj%H_diag_ss(traj%state_diag,traj%state_diag)))+sum_vk**2
          if (sum_vk<0.d0) then
            factor=(sum_vk+sqrt(deltaE))/2.d0/sum_kk
          else
            factor=(sum_vk-sqrt(deltaE))/2.d0/sum_kk
          endif
          do i=1,3
            traj%veloc_ad(:,i)=traj%veloc_ad(:,i)-factor*&
            &traj%hopping_direction_ssad(traj%state_diag_old,traj%state_diag,:,i)/traj%mass_a(:)
          enddo
          if (printlevel>2) then
            write(u_log,'(A)') 'Velocity is rescaled along projected gradient difference vector.'
            write(u_log,'(A,1X,E16.8,1X,E16.8)') 'a, b: ', sum_kk, sum_vk
            write(u_log,'(A,1X,E16.8)') 'Delta is          ',deltaE
            write(u_log,'(A,1X,F12.6)') 'Scaling factor is ',factor
          endif
        case (7)
          call available_ekin(ctrl%natom,&
          &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag_old,traj%state_diag,:,:),&
          &traj%mass_a, sum_kk, sum_vk)
          deltaE=4.d0*sum_kk*(traj%Etot-traj%Ekin-&
          &real(traj%H_diag_ss(traj%state_diag,traj%state_diag)))+sum_vk**2
          if (sum_vk<0.d0) then
            factor=(sum_vk+sqrt(deltaE))/2.d0/sum_kk
          else
            factor=(sum_vk-sqrt(deltaE))/2.d0/sum_kk
          endif
          do i=1,3
            traj%veloc_ad(:,i)=traj%veloc_ad(:,i)-factor*&
            &traj%hopping_direction_ssad(traj%state_diag_old,traj%state_diag,:,i)/traj%mass_a(:)
          enddo
          if (printlevel>2) then
            write(u_log,'(A)') 'Velocity is rescaled along effective non-adiabatic coupling vector.'
            write(u_log,'(A,1X,E16.8,1X,E16.8)') 'a, b: ', sum_kk, sum_vk
            write(u_log,'(A,1X,E16.8)') 'Delta is          ',deltaE
            write(u_log,'(A,1X,F12.6)') 'Scaling factor is ',factor
          endif
        case (8)
          call available_ekin(ctrl%natom,&
          &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag_old,traj%state_diag,:,:),&
          &traj%mass_a, sum_kk, sum_vk)
          deltaE=4.d0*sum_kk*(traj%Etot-traj%Ekin-&
          &real(traj%H_diag_ss(traj%state_diag,traj%state_diag)))+sum_vk**2
          if (sum_vk<0.d0) then
            factor=(sum_vk+sqrt(deltaE))/2.d0/sum_kk
          else
            factor=(sum_vk-sqrt(deltaE))/2.d0/sum_kk
          endif
          do i=1,3
            traj%veloc_ad(:,i)=traj%veloc_ad(:,i)-factor*&
            &traj%hopping_direction_ssad(traj%state_diag_old,traj%state_diag,:,i)/traj%mass_a(:)
          enddo
          if (printlevel>2) then
            write(u_log,'(A)') 'Velocity is rescaled along projected effective non-adiabatic coupling vector.'
            write(u_log,'(A,1X,E16.8,1X,E16.8)') 'a, b: ', sum_kk, sum_vk
            write(u_log,'(A,1X,E16.8)') 'Delta is          ',deltaE
            write(u_log,'(A,1X,F12.6)') 'Scaling factor is ',factor
          endif
        endselect
    case (2)
      if (printlevel>2) write(u_log,'(A)') 'Frustrated jump.'
      select case (ctrl%reflect_frustrated)
        case (0)
          if (printlevel>2) write(u_log,*) 'Velocity is not reflected.'
        case (1)
          if (printlevel>2) write(u_log,*) 'Velocity is reflected completely.'
          do i=1,ctrl%natom
            if (ctrl%atommask_a(i)) traj%veloc_ad(i,:) = -traj%veloc_ad(i,:)
          enddo
        case (2)
          if (printlevel>2) write(u_log,*) 'Velocity is reflected along projected velocity vector.'
          call reflect_nac(ctrl%natom,traj%veloc_ad,traj%mass_a,&
          &traj%frustrated_hop_vec_ssad(traj%state_diag_frust, traj%state_diag,:,:),&
          &real(traj%gmatrix_ssad(traj%state_diag, traj%state_diag,:,:)),&
          &real(traj%gmatrix_ssad(traj%state_diag_frust, traj%state_diag_frust,:,:)) )
        case (3)
          if (printlevel>2) write(u_log,*) 'Velocity is reflected along nonadiabatic coupling vector'
          call reflect_nac(ctrl%natom,traj%veloc_ad,traj%mass_a,&
          &traj%frustrated_hop_vec_ssad(traj%state_diag_frust, traj%state_diag,:,:),&
          &real(traj%gmatrix_ssad(traj%state_diag, traj%state_diag,:,:)),&
          &real(traj%gmatrix_ssad(traj%state_diag_frust, traj%state_diag_frust,:,:)) )
        case (4)
          if (printlevel>2) write(u_log,*) 'Velocity is reflected along gradient difference vector'
          call reflect_nac(ctrl%natom,traj%veloc_ad,traj%mass_a,&
          &traj%frustrated_hop_vec_ssad(traj%state_diag_frust, traj%state_diag,:,:),&
          &real(traj%gmatrix_ssad(traj%state_diag, traj%state_diag,:,:)),&
          &real(traj%gmatrix_ssad(traj%state_diag_frust, traj%state_diag_frust,:,:)) )
        case (5)
          if (printlevel>2) write(u_log,*) 'Velocity is reflected along projected nonadiabatic coupling vector'
          call reflect_nac(ctrl%natom,traj%veloc_ad,traj%mass_a,&
          &traj%frustrated_hop_vec_ssad(traj%state_diag_frust, traj%state_diag,:,:),&
          &real(traj%gmatrix_ssad(traj%state_diag, traj%state_diag,:,:)),&
          &real(traj%gmatrix_ssad(traj%state_diag_frust, traj%state_diag_frust,:,:)) )
        case (6)
          if (printlevel>2) write(u_log,*) 'Velocity is reflected along projected gradient difference vector'
          call reflect_nac(ctrl%natom,traj%veloc_ad,traj%mass_a,&
          &traj%frustrated_hop_vec_ssad(traj%state_diag_frust, traj%state_diag,:,:),&
          &real(traj%gmatrix_ssad(traj%state_diag, traj%state_diag,:,:)),&
          &real(traj%gmatrix_ssad(traj%state_diag_frust, traj%state_diag_frust,:,:)) )
        case (7)
          if (printlevel>2) write(u_log,*) 'Velocity is reflected along effective nonadiabatic coupling vector'
          call reflect_nac(ctrl%natom,traj%veloc_ad,traj%mass_a,&
          &traj%frustrated_hop_vec_ssad(traj%state_diag_frust, traj%state_diag,:,:),&
          &real(traj%gmatrix_ssad(traj%state_diag, traj%state_diag,:,:)),&
          &real(traj%gmatrix_ssad(traj%state_diag_frust, traj%state_diag_frust,:,:)) )
        case (8)
          if (printlevel>2) write(u_log,*) 'Velocity is reflected along projected effective nonadiabatic coupling vector'
          call reflect_nac(ctrl%natom,traj%veloc_ad,traj%mass_a,&
          &traj%frustrated_hop_vec_ssad(traj%state_diag_frust, traj%state_diag,:,:),&
          &real(traj%gmatrix_ssad(traj%state_diag, traj%state_diag,:,:)),&
          &real(traj%gmatrix_ssad(traj%state_diag_frust, traj%state_diag_frust,:,:)) )
        case (91)
          if (printlevel>2) write(u_log,*) 'Velocity is reflected according to deltaV approach along velocity vector.'
          call reflect_deltaV(ctrl%natom,traj%veloc_ad,traj%mass_a,&
          &traj%frustrated_hop_vec_ssad(traj%state_diag_frust, traj%state_diag,:,:),&
          &real(traj%gmatrix_ssad(traj%state_diag, traj%state_diag,:,:)),&
          &real(traj%gmatrix_ssad(traj%state_diag_frust, traj%state_diag_frust,:,:)) )
        case (92)
          if (printlevel>2) write(u_log,*) 'Velocity is reflected according to deltaV approach along projected velocity vector.'
          call reflect_deltaV(ctrl%natom,traj%veloc_ad,traj%mass_a,&
          &traj%frustrated_hop_vec_ssad(traj%state_diag_frust, traj%state_diag,:,:),&
          &real(traj%gmatrix_ssad(traj%state_diag, traj%state_diag,:,:)),&
          &real(traj%gmatrix_ssad(traj%state_diag_frust, traj%state_diag_frust,:,:)) )
        case (93)
          if (printlevel>2) write(u_log,*) 'Velocity is reflected according to deltaV approach along nonadiabatic coupling vector.'
          call reflect_deltaV(ctrl%natom,traj%veloc_ad,traj%mass_a,&
          &traj%frustrated_hop_vec_ssad(traj%state_diag_frust, traj%state_diag,:,:),&
          &real(traj%gmatrix_ssad(traj%state_diag, traj%state_diag,:,:)),&
          &real(traj%gmatrix_ssad(traj%state_diag_frust, traj%state_diag_frust,:,:)) )
        case (94)
          if (printlevel>2) write(u_log,*) 'Velocity is reflected according to deltaV approach along gradient difference vector.'
          call reflect_deltaV(ctrl%natom,traj%veloc_ad,traj%mass_a,&
          &traj%frustrated_hop_vec_ssad(traj%state_diag_frust, traj%state_diag,:,:),&
          &real(traj%gmatrix_ssad(traj%state_diag, traj%state_diag,:,:)),&
          &real(traj%gmatrix_ssad(traj%state_diag_frust, traj%state_diag_frust,:,:)) )
        case (95)
          if (printlevel>2) write(u_log,*) 'Velocity is reflected according to deltaV approach along projected nonadiabatic coupling vector.'
          call reflect_deltaV(ctrl%natom,traj%veloc_ad,traj%mass_a,&
          &traj%frustrated_hop_vec_ssad(traj%state_diag_frust, traj%state_diag,:,:),&
          &real(traj%gmatrix_ssad(traj%state_diag, traj%state_diag,:,:)),&
          &real(traj%gmatrix_ssad(traj%state_diag_frust, traj%state_diag_frust,:,:)) )
        case (96)
          if (printlevel>2) write(u_log,*) 'Velocity is reflected according to deltaV approach along projected gradient difference vector.'
          call reflect_deltaV(ctrl%natom,traj%veloc_ad,traj%mass_a,&
          &traj%frustrated_hop_vec_ssad(traj%state_diag_frust, traj%state_diag,:,:),&
          &real(traj%gmatrix_ssad(traj%state_diag, traj%state_diag,:,:)),&
          &real(traj%gmatrix_ssad(traj%state_diag_frust, traj%state_diag_frust,:,:)) )
        case (97)
          if (printlevel>2) write(u_log,*) 'Velocity is reflected according to deltaV approach along effective nonadiabatic coupling vector.'
          call reflect_deltaV(ctrl%natom,traj%veloc_ad,traj%mass_a,&
          &traj%frustrated_hop_vec_ssad(traj%state_diag_frust, traj%state_diag,:,:),&
          &real(traj%gmatrix_ssad(traj%state_diag, traj%state_diag,:,:)),&
          &real(traj%gmatrix_ssad(traj%state_diag_frust, traj%state_diag_frust,:,:)) )
        case (98)
          if (printlevel>2) write(u_log,*) 'Velocity is reflected according to deltaV approach along projected effective nonadiabatic coupling vector.'
          call reflect_deltaV(ctrl%natom,traj%veloc_ad,traj%mass_a,&
          &traj%frustrated_hop_vec_ssad(traj%state_diag_frust, traj%state_diag,:,:),&
          &real(traj%gmatrix_ssad(traj%state_diag, traj%state_diag,:,:)),&
          &real(traj%gmatrix_ssad(traj%state_diag_frust, traj%state_diag_frust,:,:)) )
      endselect
    case (3)
          if (printlevel>2) write(u_log,*) 'Velocity is not rescaled after resonant surface hop.'
  endselect

endsubroutine

! ===========================================================

!> calculates the following sums:
!> sum_kk=1/2*sum_atom nac_atom*nac_atom/mass_atom
!> sum_vk=sum_atom vel_atom*nac_atom
subroutine available_ekin(natom,veloc_ad,nac_ad,mass_a, sum_kk, sum_vk)
  implicit none
  integer, intent(in) :: natom
  real*8, intent(in) :: veloc_ad(natom,3), nac_ad(natom,3), mass_a(natom)
  real*8, intent(out) :: sum_kk, sum_vk

  integer :: idir

  sum_kk=0.d0
  sum_vk=0.d0
  do idir=1,3
    sum_kk=sum_kk+0.5d0*sum( nac_ad(:,idir)*nac_ad(:,idir)/mass_a(:) )
    sum_vk=sum_vk+      sum( nac_ad(:,idir)*veloc_ad(:,idir) )
  enddo

endsubroutine

! ===========================================================

subroutine reflect_nac(natom,veloc_ad,mass_a,nac_ad,Gdiag,Gfrust)
  use definitions
  implicit none
  integer, intent(in) :: natom
  real*8, intent(inout) :: veloc_ad(natom,3)
  real*8, intent(in) :: mass_a(natom), nac_ad(natom,3), Gdiag(natom, 3), Gfrust(natom, 3)

  integer :: idir, iat
  real*8 :: mass_ad(natom,3)
  real*8 :: sum_pk, sum_Fdiagk, sum_Ffrustk
  real*8 :: factor

  do idir=1,3
    mass_ad(:,idir) = mass_a(:)
  enddo
  
  sum_Fdiagk  = sum(-Gdiag*nac_ad)
  sum_Ffrustk = sum(-Gfrust*nac_ad)
  sum_pk = sum(mass_ad*veloc_ad*nac_ad)
  if (printlevel>2)  then
    write(u_log,*)'Checking velocity reflection along the frustrated hop reflection vector.'
    write(u_log,'(A,4E14.6)') ' sum_Fdiagk, sum_Ffrustk, sum_vk, sum_kk', sum_Fdiagk, sum_Ffrustk, sum_pk,&
    &sum(nac_ad*nac_ad)
  endif
!   write(u_log,'(A,4E14.6)') 'fptmp: pk, pv, pp', sum(mass_ad*veloc_ad*nac_ad),&
!   &sum(mass_ad*veloc_ad*veloc_ad), sum(mass_ad*mass_ad*veloc_ad*veloc_ad)
    
  if ((sum_Fdiagk*sum_Ffrustk<0).and.(sum_Ffrustk*sum_pk<0))then
    if (printlevel>2) write(u_log,*) 'Conditions for reflection fulfilled.'
    ! Reverse the velocity for individual atoms
    !  -> this conserves p and Ekin
    do iat=1,natom
      factor = 2 * sum(veloc_ad(iat,:)*nac_ad(iat,:)) / sum(nac_ad(iat,:)*nac_ad(iat,:))
      veloc_ad(iat,:) = veloc_ad(iat,:) - factor * nac_ad(iat,:)
    enddo
!     write(u_log,'(A,4E14.6)') 'fptmp: pk, pv, pp', sum(mass_ad*veloc_ad*nac_ad),&
!     &sum(mass_ad*veloc_ad*veloc_ad), sum(mass_ad*mass_ad*veloc_ad*veloc_ad)
  else
    if (printlevel>2) write(u_log,*) 'Conditions for reflection not fulfilled.'
  endif

endsubroutine

! ===========================================================

subroutine reflect_deltaV(natom,veloc_ad,mass_a,nac_ad,Gdiag,Gfrust)
  use definitions
  implicit none
  integer, intent(in) :: natom
  real*8, intent(inout) :: veloc_ad(natom,3)
  real*8, intent(in) :: mass_a(natom), nac_ad(natom,3), Gdiag(natom, 3), Gfrust(natom, 3)
  integer :: idir, iat
  real*8 :: mass_ad(natom,3)
  real*8 :: sum_pk, sum_Fdiagk, sum_Ffrustk
  real*8 :: factor

  do idir=1,3
    mass_ad(:,idir) = mass_a(:)
  enddo

  ! Notice in deltaV approach, it is -Grust dot nac_ad; 
  ! the input nac_ad is nac_diagfrust_diag
  ! therefore, Gfrust*nac_ad = -Grust * nac_ad
  sum_Fdiagk  = sum(-Gdiag*nac_ad)
  sum_Ffrustk = sum(-Gfrust*nac_ad)
  sum_pk = sum(mass_ad*veloc_ad*nac_ad)
  if (printlevel>2)  then
    write(u_log,*)'Checking velocity reflection along the frustrated hop reflection vector.'
    write(u_log,'(A,4E14.6)') ' sum_Fdiagk, sum_Ffrustk, sum_vk, sum_kk', sum_Fdiagk, sum_Ffrustk, sum_pk,&
    &sum(nac_ad*nac_ad)
  endif

  if (sum_Ffrustk*sum_pk<0) then
    if (printlevel>2) write(u_log,*) 'Conditions for deltaV reflection approach fulfilled.'
    do iat=1,natom
      factor = 2 * sum(veloc_ad(iat,:)*nac_ad(iat,:)) / sum(nac_ad(iat,:)*nac_ad(iat,:))
      veloc_ad(iat,:) = veloc_ad(iat,:) - factor * nac_ad(iat,:)
    enddo
  else
    if (printlevel>2) write(u_log,*) 'Conditions for deltaV reflection not fulfilled.'
  endif

endsubroutine

! ===========================================================

!> multiplies the velocity vector by the sqrt of the damping factor
subroutine Damp_velocities(traj,ctrl)
  use definitions
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

  if (printlevel>2) then
    write(u_log,*) '============================================================='
    write(u_log,*) '                     Velocity Damping'
    write(u_log,*) '============================================================='
    write(u_log,*)
    write(u_log,*) 'Factor for the velocities is',sqrt(ctrl%dampeddyn)
  endif
  traj%veloc_ad=traj%veloc_ad*sqrt(ctrl%dampeddyn)

endsubroutine

! ===========================================================

!> Check the consistency between time steps, for adaptive time step
subroutine Check_Consistency(traj,ctrl)
  use definitions
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

  real*8 :: diff

  ! initialize consistency and discrepancy
  traj%consistency=0
  traj%discrepancy=1.d0

  ! Check energy convergence
  diff=abs(traj%Etot-traj%Etot_old)

  if (diff .gt. ctrl%convthre) then  !if energy difference is bigger than threshold
    traj%consistency=2
    traj%discrepancy=2.0
  else if (traj%step .gt. 2 .and. diff .lt. (ctrl%convthre/5.0)) then
    traj%consistency=1
    traj%discrepancy=0.5
  endif

endsubroutine

! ===========================================================

!> Change the step size
subroutine Adaptive_Stepsize(traj,ctrl)
  use definitions
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

  real*8 :: oldstep

  oldstep=ctrl%dtstep
  ctrl%dtstep=ctrl%dtstep/traj%discrepancy

  if (printlevel>2) then
    write(u_log,*) 'Time step size has been adapted'
    write(u_log,'(A,1X,F9.6)') ' Old step size in fs:', oldstep*au2fs
    write(u_log,'(A,1X,F9.6)') ' New step size in fs:', ctrl%dtstep*au2fs
  endif
  if (ctrl%dtstep.lt.ctrl%dtmin) then
    if (printlevel>2) then
      write(u_log,*) "adaptive time step increased to minimum value:",ctrl%dtmin*au2fs,"fs"
    endif
    ctrl%dtstep=ctrl%dtmin
  endif

endsubroutine

! ===========================================================

!> Read restart file and back propagate
subroutine Back_propagate(u_traj,traj,ctrl)
  use definitions
  use matrix
  use misc
  use decoherence_afssh
  implicit none
  integer :: u_ctrl,u_traj
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

  integer :: imult, iatom, i,j,k, istate,ilaser
  character(8000) :: string
  real*8 :: dummy_randnum
  integer :: time

  if (printlevel>0) then
    write(u_log,*) 'Back propagate trajectory'
    write(u_log,*) ' '
  endif

  ! These are from read_restart of the trajecotry part
  open(unit=u_traj, file='restart.traj',status='old',action='read')
  ! read everything
  read(u_traj,*) traj%RNGseed
  read(u_traj,*) traj%traj_hash
  read(u_traj,*) traj%state_MCH
  read(u_traj,*) traj%state_MCH_old
  read(u_traj,*) traj%state_diag
  read(u_traj,*) traj%state_diag_old
  read(u_traj,*) traj%state_diag_frust
  if (ctrl%zpe_correction==1) then
    read(u_traj,*) (traj%state_pumping_s(i),i=1,ctrl%nstates)
    read(u_traj,*) (traj%pumping_status_s(i),i=1,ctrl%nstates)
  else if (ctrl%zpe_correction==2) then
    read(u_traj,*) (traj%lpzpe_ke_ah(i),i=1,ctrl%lpzpe_nah)
    read(u_traj,*) (traj%lpzpe_ke_bc(i),i=1,ctrl%lpzpe_nbc)
    read(u_traj,*) traj%lpzpe_cycle
    read(u_traj,*) traj%lpzpe_iter_incycle
    read(u_traj,*) traj%lpzpe_starttime
    read(u_traj,*) traj%lpzpe_endtime
    read(u_traj,*) traj%in_cycle
  endif
  if (ctrl%time_uncertainty==1) then
    read(u_traj,*) (traj%uncertainty_time_s(i),i=1,ctrl%nstates)
    read(u_traj,*) (traj%allowed_hop_s(i),i=1,ctrl%nstates)
    read(u_traj,*) (traj%allowed_time_s(i),i=1,ctrl%nstates)
    read(u_traj,*) traj%uncertainty_time_frust
    read(u_traj,*) traj%incident_time
    read(u_traj,*) traj%in_time_uncertainty
    read(u_traj,*) traj%tu_backpropagation
    read(u_traj,*) traj%travelling_state
    read(u_traj,*) traj%target_state
  endif
  read(u_traj,*) traj%step
  read(u_traj,*) traj%microtime
  read(u_traj,*) traj%Ekin
  read(u_traj,*) traj%Epot
  read(u_traj,*) traj%Etot
  read(u_traj,*) traj%Ekin_old
  read(u_traj,*) traj%Epot_old
  read(u_traj,*) traj%Etot_old
  read(u_traj,*) traj%time_start
  read(u_traj,*) traj%time_last
  read(u_traj,*) traj%time_step
  read(u_traj,*) traj%kind_of_jump
  read(u_traj,*) traj%steps_in_gs
  read(u_traj,*) (traj%ncids(i),i=1,10)
  read(u_traj,*) traj%nc_index
  traj%nc_index=-traj%nc_index

  ! read the arrays
  read(u_traj,*) (traj%atomicnumber_a(iatom),iatom=1,ctrl%natom)
  read(u_traj,*) (traj%element_a(iatom),iatom=1,ctrl%natom)
  read(u_traj,*) (traj%mass_a(iatom),iatom=1,ctrl%natom)
  call vec3read(ctrl%natom, traj%geom_ad,  u_traj, string)
  call vec3read(ctrl%natom, traj%veloc_ad, u_traj, string)
  call vec3read(ctrl%natom, traj%veloc_old_ad, u_traj, string)
  call vec3read(ctrl%natom, traj%veloc_app_ad, u_traj, string)
  call vec3read(ctrl%natom, traj%accel_ad, u_traj, string)
  if (ctrl%nac_projection==1) then
    call matread(3*ctrl%natom, traj%trans_rot_P, u_traj, string)
  endif

  call matread(ctrl%nstates, traj%H_MCH_ss,     u_traj,   string)
  call matread(ctrl%nstates, traj%dH_MCH_ss,    u_traj,   string)
  call matread(ctrl%nstates, traj%dH_MCH_old_ss,u_traj,   string)
  call matread(ctrl%nstates, traj%H_MCH_old_ss, u_traj,   string)
  call matread(ctrl%nstates, traj%H_MCH_old2_ss,u_traj,   string)
  call matread(ctrl%nstates, traj%H_MCH_old3_ss,u_traj,   string)
  call matread(ctrl%nstates, traj%H_diag_ss,    u_traj,   string)
  call matread(ctrl%nstates, traj%H_diag_old_ss,    u_traj,   string)
  call matread(ctrl%nstates, traj%H_diag_old2_ss,    u_traj,   string)
  call matread(ctrl%nstates, traj%H_diag_old3_ss,    u_traj,   string)
  if (ctrl%zpe_correction==1) then
    call vecread(ctrl%nstates, traj%Evib_local_s, u_traj,   string)
    call vecread(ctrl%nstates, traj%Ezpe_local_s, u_traj,   string)
  endif
  call matread(ctrl%nstates, traj%U_ss,         u_traj,   string)
  call matread(ctrl%nstates, traj%U_old_ss,     u_traj,   string)
  call matread(ctrl%nstates, traj%NACdt_ss,     u_traj,   string)
  call matread(ctrl%nstates, traj%NACdt_old_ss, u_traj,   string)
  call matread(ctrl%nstates, traj%overlaps_ss,  u_traj,   string)
  call matread(ctrl%nstates, traj%DM_ssd(:,:,1),  u_traj, string)
  call matread(ctrl%nstates, traj%DM_ssd(:,:,2),  u_traj, string)
  call matread(ctrl%nstates, traj%DM_ssd(:,:,3),  u_traj, string)
  call matread(ctrl%nstates, traj%DM_old_ssd(:,:,1),  u_traj, string)
  call matread(ctrl%nstates, traj%DM_old_ssd(:,:,2),  u_traj, string)
  call matread(ctrl%nstates, traj%DM_old_ssd(:,:,3),  u_traj, string)
  call matread(ctrl%nstates, traj%DM_print_ssd(:,:,1),  u_traj, string)
  call matread(ctrl%nstates, traj%DM_print_ssd(:,:,2),  u_traj, string)
  call matread(ctrl%nstates, traj%DM_print_ssd(:,:,3),  u_traj, string)
!     call matread(ctrl%nstates, traj%Property_ss,  u_traj,   string)
  call matread(ctrl%nstates, traj%Rtotal_ss,    u_traj,   string)
  call matread(ctrl%nstates, traj%RDtotal_ss,    u_traj,   string)
  call matread(ctrl%nstates, traj%Dtotal_ss,    u_traj,   string)
  call matread(ctrl%nstates, traj%dendt_MCH_ss,    u_traj,   string)
  call matread(ctrl%nstates, traj%dendt_diag_ss,    u_traj,   string)
  call matread(ctrl%nstates, traj%Kmatrix_MCH_ss,    u_traj,   string)
  call matread(ctrl%nstates, traj%Kmatrix_diag_ss,    u_traj,   string)
  call vecread(ctrl%nstates, traj%phases_s, u_traj,       string)
  call vecread(ctrl%nstates, traj%phases_old_s, u_traj,   string)
  call vecread(ctrl%nstates, traj%hopprob_s, u_traj,      string)
  call vecread(ctrl%nstates, traj%switchprob_s, u_traj,      string)
  call vec3read(ctrl%nstates, traj%gpreprob_s3,  u_traj, string)
  call vec3read(ctrl%nstates, traj%preprob_s3,  u_traj, string)
  call vec3read(ctrl%nstates, traj%preprob_old_s3,  u_traj, string)
  call vecread(ctrl%nstates, traj%decotime_diag_s, u_traj,       string)
  call vecread(ctrl%nstates, traj%decotime_diag_old_s, u_traj,   string)
  call vecread(ctrl%nstates, traj%decotime_MCH_s, u_traj,       string)
  call vecread(ctrl%nstates, traj%decotime_MCH_old_s, u_traj,   string)
  do i=1,ctrl%nstates
    call vec3read(ctrl%natom,traj%svec_MCH_sad(i,:,:),u_traj,string)
  enddo
  do i=1,ctrl%nstates
    call vec3read(ctrl%natom,traj%svec_diag_sad(i,:,:),u_traj,string)
  enddo
  do i=1,ctrl%nstates
    call vec3read(ctrl%natom,traj%psvec_MCH_sad(i,:,:),u_traj,string)
  enddo
  do i=1,ctrl%nstates
    call vec3read(ctrl%natom,traj%psvec_diag_sad(i,:,:),u_traj,string)
  enddo

  read(u_traj,*) traj%randnum
  read(u_traj,*) traj%randnum2
  read(u_traj,*) traj%dmag
  read(u_traj,*) traj%dmag1
  read(u_traj,*) traj%dmag2

  if (ctrl%calc_dipolegrad>-1) then
    do i=1,ctrl%nstates
      do j=1,ctrl%nstates
        do k=1,3
          call vec3read(ctrl%natom,traj%DMgrad_ssdad(i,j,k,:,:),u_traj,string)
        enddo
      enddo
    enddo
  endif
  if (ctrl%calc_nacdr>-1) then
    do i=1,ctrl%nstates
      do j=1,ctrl%nstates
        call vec3read(ctrl%natom,traj%NACdr_ssad(i,j,:,:),u_traj,string)
      enddo
    enddo
    do i=1,ctrl%nstates
      do j=1,ctrl%nstates
        call vec3read(ctrl%natom,traj%NACdr_old_ssad(i,j,:,:),u_traj,string)
      enddo
    enddo
    do i=1,ctrl%nstates
      do j=1,ctrl%nstates
        call vec3read(ctrl%natom,traj%NACdr_diag_ssad(i,j,:,:),u_traj,string)
      enddo
    enddo
    if (ctrl%nac_projection==1) then
      do i=1,ctrl%nstates
         do j=1,ctrl%nstates
          call vec3read(ctrl%natom,traj%pNACdr_MCH_ssad(i,j,:,:),u_traj,string)
         enddo
      enddo
      do i=1,ctrl%nstates
        do j=1,ctrl%nstates
          call vec3read(ctrl%natom,traj%pNACdr_diag_ssad(i,j,:,:),u_traj,string)
        enddo
      enddo
    endif
  endif
  do i=1,ctrl%nstates
    call vec3read(ctrl%natom,traj%grad_MCH_sad(i,:,:),u_traj,string)
  enddo
  do i=1,ctrl%nstates
    call vec3read(ctrl%natom,traj%grad_MCH_old_sad(i,:,:),u_traj,string)
  enddo
  if (ctrl%zpe_correction==1) then
    do i=1,ctrl%nstates
      call vec3read(ctrl%natom,traj%grad_MCH_old2_sad(i,:,:),u_traj,string)
    enddo
    do i=1,ctrl%nstates
      call vec3read(ctrl%natom,traj%grad_Ezpe_local_sad(i,:,:),u_traj,string)
    enddo
  endif
  do i=1,ctrl%nstates
    do j=1,ctrl%nstates
      call vec3read(ctrl%natom,traj%hopping_direction_ssad(i,j,:,:),u_traj,string)
    enddo
  enddo
  do i=1,ctrl%nstates
    do j=1,ctrl%nstates
      call vec3read(ctrl%natom,traj%frustrated_hop_vec_ssad(i,j,:,:),u_traj,string)
    enddo
  enddo
  do i=1,ctrl%nstates
    do j=1,ctrl%nstates
      call vec3read(ctrl%natom,traj%Gmatrix_ssad(i,j,:,:),u_traj,string)
    enddo
  enddo
  if (ctrl%calc_effectivenac==1) then
    do i=1,ctrl%nstates
      do j=1,ctrl%nstates
        call vec3read(ctrl%natom,traj%NACGV_MCH_ssad(i,j,:,:),u_traj,string)
      enddo
    enddo
    do i=1,ctrl%nstates
      do j=1,ctrl%nstates
        call vec3read(ctrl%natom,traj%NACGV_diag_ssad(i,j,:,:),u_traj,string)
      enddo
    enddo
    if (ctrl%nac_projection==1) then
      do i=1,ctrl%nstates
        do j=1,ctrl%nstates
          call vec3read(ctrl%natom,traj%pNACGV_MCH_ssad(i,j,:,:),u_traj,string)
        enddo
      enddo
      do i=1,ctrl%nstates
        do j=1,ctrl%nstates
          call vec3read(ctrl%natom,traj%pNACGV_diag_ssad(i,j,:,:),u_traj,string)
        enddo
      enddo
    endif
  endif
  call vec3read(ctrl%natom,traj%grad_ad(:,:),u_traj,string)
  call vec3read(ctrl%natom,traj%decograd_ad(:,:),u_traj,string)

  call vecread(ctrl%nstates, traj%coeff_diag_s, u_traj, string)
  call vecread(ctrl%nstates, traj%coeff_diag_old_s, u_traj, string)
  call vecread(ctrl%nstates, traj%coeff_mch_s, u_traj, string)
  call vecread(ctrl%nstates, traj%coeff_mch_old_s, u_traj, string)
  call vecread(ctrl%nstates, traj%ccoeff_diag_s, u_traj, string)
  call vecread(ctrl%nstates, traj%ccoeff_diag_old_s, u_traj, string)
  call vecread(ctrl%nstates, traj%ccoeff_mch_s, u_traj, string)
  if (ctrl%zpe_correction==1) then
    call vec2read(ctrl%nstates, traj%coeff_zpe_s2, u_traj, string)
  endif
  if (ctrl%integrator==0) then 
    call vecread(ctrl%nstates, traj%gRcoeff_mch_s, u_traj, string)
    call vecread(ctrl%nstates, traj%gIcoeff_mch_s, u_traj, string)
    call vecread(ctrl%nstates, traj%gRccoeff_mch_s, u_traj, string)
    call vecread(ctrl%nstates, traj%gIccoeff_mch_s, u_traj, string)
    call vecread(ctrl%nstates, traj%ephase_s, u_traj, string)
    call vecread(ctrl%nstates, traj%gephase_s, u_traj, string)
  endif

  read(u_traj,*) (traj%selg_s(i),i=1,ctrl%nstates)
  do i=1,ctrl%nstates
    read(u_traj,*) (traj%selt_ss(i,j),j=1,ctrl%nstates)
  enddo
  if (ctrl%calc_dipolegrad>-1) then
    do i=1,ctrl%nstates
      read(u_traj,*) (traj%seldm_ss(i,j),j=1,ctrl%nstates)
    enddo
  endif
  read(u_traj,*) traj%phases_found

  call vecread(ctrl%n_property1d, traj%Property1d_labels_y, u_traj, string)
  call vecread(ctrl%n_property2d, traj%Property2d_labels_x, u_traj, string)
  do i=1,ctrl%n_property1d
    call vecread(ctrl%nstates, traj%Property1d_ys(i,:), u_traj, string)
  enddo
  do i=1,ctrl%n_property2d
    call matread(ctrl%nstates, traj%Property2d_xss(i,:,:), u_traj, string)
  enddo

  ! read restart info for the auxilliary trajectories
  if (ctrl%decoherence==2) then
    call allocate_afssh(traj, ctrl)
    do i=1,ctrl%nstates
      read(u_traj,*) traj%auxtrajs_s(i)%istate
      read(u_traj,*) traj%auxtrajs_s(i)%rate1
      read(u_traj,*) traj%auxtrajs_s(i)%rate2
      call vec3read(ctrl%natom, traj%auxtrajs_s(i)%geom_ad,   u_traj, string)
      call vec3read(ctrl%natom, traj%auxtrajs_s(i)%veloc_ad,  u_traj, string)
      call vec3read(ctrl%natom, traj%auxtrajs_s(i)%accel_ad,  u_traj, string)
      call vec3read(ctrl%natom, traj%auxtrajs_s(i)%grad_ad,   u_traj, string)
      call vec3read(ctrl%natom, traj%auxtrajs_s(i)%geom_tmp_ad,   u_traj, string)
      call vec3read(ctrl%natom, traj%auxtrajs_s(i)%veloc_tmp_ad,  u_traj, string)
    enddo
  endif

  if (ctrl%army_ants==1) then
    read(u_traj,*) traj%army_ants_weight
    read(u_traj,*) traj%randnum_branching
    read(u_traj,*) traj%branching_likehood
  endif

  if (ctrl%pointer_basis==2) then
    read(u_traj,*) traj%entropy
    read(u_traj,*) traj%entropy_old
    read(u_traj,*) traj%linear_entropy
    read(u_traj,*) traj%linear_entropy_old
    read(u_traj,*) traj%entropy_grad
    read(u_traj,*) traj%linear_entropy_grad
    call matread(ctrl%nstates, traj%U_pointer_ss,  u_traj,   string)
    call matread(ctrl%nstates, traj%U_pointer_old_ss,  u_traj,   string)
  endif

  ! Trajectory consistency
  read(u_traj,*) traj%discrepancy

  close(u_traj)
  ! Now everything about trajectory has been read from the restart file

endsubroutine

endmodule
