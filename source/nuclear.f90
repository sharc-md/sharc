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

!> # Module NUCLEAR
!> 
!> \author Sebastian Mai
!> \ðate 27.02.2015
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

  ! variables for thermostat
  real*8 :: b

  ! variables for constraints
  real*8 :: initdistvec(ctrl%n_constraints,3)
  logical :: check_constraints(ctrl%n_constraints)
  integer :: iconstr, iA, iB, iiter
  real*8 :: relpos(3)
  real*8 :: d2t, coeff


  if (printlevel>2) then
    write(u_log,*) '============================================================='
    write(u_log,*) '              Velocity Verlet -- X-step'
    write(u_log,*) '============================================================='
    call vec3write(ctrl%natom,traj%accel_ad,u_log,'Old accel','F12.9')
    call vec3write(ctrl%natom,traj%geom_ad,u_log,'Old geom','F12.7')
  endif

  ! prepare if RATTLE
  if (ctrl%do_constraints==1) then
    do iconstr = 1,ctrl%n_constraints
      iA=ctrl%constraints_ca(iconstr,1)
      iB=ctrl%constraints_ca(iconstr,2)
      initdistvec(iconstr,:) = traj%geom_ad(iA,:)-traj%geom_ad(iB,:)
    enddo
  endif

  ! propagate positions
  select case (ctrl%thermostat)
    case (0) ! no thermostat
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
    case (1) ! Langevin thermostat
      traj%thermostat_random=gaussian_random(ctrl%natom,real(0,8),ctrl%temperature) !ctrl%temperature is variance here
      !write(u_log,*) traj%thermostat_random
      do iatom=1,ctrl%natom
        b=1/(1+ctrl%thermostat_const(1)*ctrl%dtstep/(2*traj%mass_a(iatom)))
        do idir=1,3                 ! propagate positions according to Langevin equation
          traj%accel_ad(iatom,idir)=&
          &-traj%grad_ad(iatom,idir)/traj%mass_a(iatom)

          traj%geom_ad(iatom,idir)=&
          & traj%geom_ad(iatom,idir)&
          &+b*traj%veloc_ad(iatom,idir)*ctrl%dtstep&
          &+0.5d0*b*traj%accel_ad(iatom,idir)*ctrl%dtstep**2&
          &+0.5d0*b*traj%thermostat_random(3*(iatom-1)+idir)*ctrl%dtstep/traj%mass_a(iatom)
        enddo
      enddo
  end select

  ! carry out RATTLE
  if (ctrl%do_constraints==1) then
    if (printlevel>2) then
      write(u_log,'(A)') 'RATTLE iterations'
      write(u_log,'(A3,X,A3,X,A5,X,A5,X,A12,X,A12,X,A2)') 'It','Con','AtomA','AtomB','Resid','coeff','OK' 
    endif
    do iiter=1,1000
      ! initialize logical variable that controls whether constraints are enforced
      check_constraints(:) = .TRUE.
        ! loop over the constrained bonds
        do iconstr = 1, ctrl%n_constraints
          ! define the atomic id of the atoms that have fixed distance
          iA=ctrl%constraints_ca(iconstr,1)
          iB=ctrl%constraints_ca(iconstr,2)
          ! compute the relative position of the two constrained atoms, and the relative norm squared
          relpos = traj%geom_ad(iA,:) - traj%geom_ad(iB,:)
          D2t = DOT_PRODUCT(relpos, relpos)
          ! when the difference with the fixed distance is significantly different from zero, do RATTLE
          coeff=0.d0
          if ( abs(D2t - ctrl%constraints_dist_c(iconstr) ) > ctrl%constraints_tol ) then
            ! compute RATTLE coefficient
            coeff = (D2t - ctrl%constraints_dist_c(iconstr)) / &
            & (2.D0 * ctrl%dtstep * DOT_PRODUCT(initdistvec(iconstr,:), relpos) * &
            & (1.D0 / traj%mass_a(iA) + 1.D0 / traj%mass_a(iB)))
            ! correct positions
            traj%geom_ad(iA,:) = traj%geom_ad(iA,:) - (coeff * ctrl%dtstep * initdistvec(iconstr,:) / traj%mass_a(iA) )
            traj%geom_ad(iB,:) = traj%geom_ad(iB,:) + (coeff * ctrl%dtstep * initdistvec(iconstr,:) / traj%mass_a(iB) )
            ! correct velocities
            traj%veloc_ad(iA,:) = traj%veloc_ad(iA,:) - (coeff * initdistvec(iconstr,:) / traj%mass_a(iA) )
            traj%veloc_ad(iB,:) = traj%veloc_ad(iB,:) + (coeff * initdistvec(iconstr,:) / traj%mass_a(iB) )
            ! this constraint was not ok
            check_constraints(iconstr) = .FALSE.
          endif
          ! print
          if (printlevel>2) then
            write(u_log,'(I3,X,I3,X,I5,X,I5,X,F12.7,X,F12.7,X,L1)') iiter,iconstr,iA,iB,&
            &abs(D2t - ctrl%constraints_dist_c(iconstr)),coeff,check_constraints(iconstr)
          endif
        end do ! end of the loop over the constraints
      ! break the loop when all constraints are satisfied
      if ( all(check_constraints) ) exit
    end do
  endif
  if (.not. all(check_constraints) ) then
    write(0,*) 'Could not satisfy RATTLE constraints in 1000 iterations (xstep)!'
    stop 1
  endif



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

  ! variables for thermostat
  real*8 :: a, b

  ! variables for constraints
  logical :: check_constraints(ctrl%n_constraints)
  integer :: iconstr, iA, iB, iiter
  real*8 :: relpos(3), relvel(3)
  real*8 :: d2t, coeff

  if (printlevel>2) then
    write(u_log,*) '============================================================='
    write(u_log,*) '              Velocity Verlet -- V-step'
    write(u_log,*) '============================================================='
    call vec3write(ctrl%natom,traj%accel_ad,u_log,'Old accel','F12.9')
    call vec3write(ctrl%natom,traj%veloc_ad,u_log,'Old veloc','F12.9')
  endif

  ! propagate velocities
  select case (ctrl%thermostat)
    case (0) ! no thermostat
      do iatom=1,ctrl%natom
        do idir=1,3
          traj%accel_ad(iatom,idir)=0.5d0*(traj%accel_ad(iatom,idir)&
          &-traj%grad_ad(iatom,idir)/traj%mass_a(iatom) )
    
          traj%veloc_ad(iatom,idir)=&
          & traj%veloc_ad(iatom,idir)&
          &+traj%accel_ad(iatom,idir)*ctrl%dtstep
        enddo
      enddo
    case (1) ! Langevin thermostat
      do iatom=1,ctrl%natom
        a=ctrl%thermostat_const(1)*ctrl%dtstep/(2*traj%mass_a(iatom))
        b=1/(1+a)
        a=(1-a)*b
        do idir=1,3
          traj%accel_ad(iatom,idir)=0.5d0*(a*traj%accel_ad(iatom,idir)&
          &-traj%grad_ad(iatom,idir)/traj%mass_a(iatom) )
    
          traj%veloc_ad(iatom,idir)=&
          & a*traj%veloc_ad(iatom,idir)&
          &+traj%accel_ad(iatom,idir)*ctrl%dtstep&
          &+b*traj%thermostat_random(3*(iatom-1)+idir)/traj%mass_a(iatom)
        enddo
      enddo
  endselect    

  ! carry out RATTLE
  if (ctrl%do_constraints==1) then
    if (printlevel>2) then
      write(u_log,'(A)') 'RATTLE iterations'
      write(u_log,'(A3,X,A3,X,A5,X,A5,X,A12,X,A12,X,A2)') 'it','Con','AtomA','AtomB','Resid','coeff','OK'
    endif
    do iiter=1,1000
      ! initialize logical variable that controls whether constraints are enforced
      check_constraints(:) = .TRUE.
        ! loop over the constrained bonds
        do iconstr = 1, ctrl%n_constraints
          ! define the atomic id of the atoms that have fixed distance
          iA = ctrl%constraints_ca(iconstr,1)
          iB = ctrl%constraints_ca(iconstr,2)
          ! compute projection of the relative velocity with respect to the distance vector of the bond
          relvel = traj%veloc_ad(iA,:) - traj%veloc_ad(iB,:)
          relpos = traj%geom_ad(iA,:) - traj%geom_ad(iB,:)
          d2t = DOT_PRODUCT(relpos, relvel)
          ! when this projection is significantly different from zero, do RATTLE
          coeff=0.d0
          if ( abs(d2t) > ctrl%constraints_tol ) then
            ! compute RATTLE coefficient
            coeff = d2t / ((1.D0 / traj%mass_a(iA) + 1.D0 / traj%mass_a(iB)) * ctrl%constraints_dist_c(iconstr))
            ! correct velocities
            traj%veloc_ad(iA,:) = traj%veloc_ad(iA,:) - (coeff * relpos(:) / traj%mass_a(iA) )
            traj%veloc_ad(iB,:) = traj%veloc_ad(iB,:) + (coeff * relpos(:) / traj%mass_a(iB) )
            ! this constraint was not ok
            check_constraints(iconstr) = .FALSE.
          endif
          ! print
          if (printlevel>2) then
            write(u_log,'(I3,X,I3,X,I5,X,I5,X,F12.7,X,F12.7,X,L1)') iiter,iconstr,iA,iB,&
            &abs(D2t),coeff,check_constraints(iconstr)
          endif
        end do ! end of the loop over the constraints
      ! break the loop when all constraints are satisfied
      if ( all(check_constraints) ) exit
    end do
  endif
  if (.not. all(check_constraints) ) then
    write(0,*) 'Could not satisfy RATTLE constraints in 1000 iterations (vstep)!'
    stop 1
  endif
  


  if (printlevel>2) then
    call vec3write(ctrl%natom,traj%accel_ad,u_log,'accel','F12.9')
    call vec3write(ctrl%natom,traj%veloc_ad,u_log,'veloc','F12.9')
  endif

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

  traj%Ekin=Calculate_ekin(ctrl%natom, traj%veloc_ad, traj%mass_a)
  traj%Epot=real(traj%H_diag_ss(traj%state_diag,traj%state_diag))
  traj%Etot=traj%Ekin+traj%Epot

  if (printlevel>2) then
    write(u_log,*) '============================================================='
    write(u_log,*) '                        Energies'
    write(u_log,*) '============================================================='
    write(u_log,'(A,1X,F14.9,1X,A)') 'Ekin:',traj%Ekin*au2ev,'eV'
    write(u_log,'(A,1X,F14.9,1X,A)') 'Epot:',traj%Epot*au2ev,'eV'
    write(u_log,'(A,1X,F14.9,1X,A)') 'Etot:',traj%Etot*au2ev,'eV'
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
          &traj%veloc_ad,real(traj%gmatrix_ssad(traj%state_diag_old, traj%state_diag,:,:)),&
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
            &real(traj%gmatrix_ssad(traj%state_diag_old, traj%state_diag,:,i))/traj%mass_a(:)
          enddo
          if (printlevel>2) then
            write(u_log,'(A)') 'Velocity is rescaled along non-adiabatic coupling vector.'
            write(u_log,'(A,1X,E16.8,1X,E16.8)') 'a, b: ', sum_kk, sum_vk
            write(u_log,'(A,1X,E16.8)') 'Delta is          ',deltaE
            write(u_log,'(A,1X,F12.6)') 'Scaling factor is ',factor
          endif
        case (3)
          call available_ekin(ctrl%natom,&
          &traj%veloc_ad,real(traj%gmatrix_ssad(traj%state_diag, traj%state_diag,:,:)-&
          &traj%gmatrix_ssad(traj%state_diag_old, traj%state_diag_old,:,:)),&
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
            &real(traj%gmatrix_ssad(traj%state_diag, traj%state_diag,:,i)-&
            &traj%gmatrix_ssad(traj%state_diag_old, traj%state_diag_old,:,i))/traj%mass_a(:)
          enddo
          if (printlevel>2) then
            write(u_log,'(A)') 'Velocity is rescaled along gradient difference vector.'
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
          call reflect_nac(ctrl%natom,traj%veloc_ad,traj%mass_a,&
          &real(traj%gmatrix_ssad(traj%state_diag_frust, traj%state_diag,:,:)),&
          &real(traj%gmatrix_ssad(traj%state_diag, traj%state_diag,:,:)),&
          &real(traj%gmatrix_ssad(traj%state_diag_frust, traj%state_diag_frust,:,:)) )
        case (3)
          call reflect_nac(ctrl%natom,traj%veloc_ad,traj%mass_a,&
          &real(traj%gmatrix_ssad(traj%state_diag_frust, traj%state_diag_frust,:,:)-&
          &traj%gmatrix_ssad(traj%state_diag, traj%state_diag,:,:)),&
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
    write(u_log,*)'Checking velocity reflection along the nonadiabatic coupling vector.'
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

!> returns 3*natoms gaussian distributed random numbers with mean=mu and variance=var
function gaussian_random(natoms,mu,var)
  use definitions
  implicit none
  integer,intent(in) :: natoms
  integer :: i
  real*8,intent(in) :: mu, var
  real*8 :: gaussian_random(2*((natoms*3+1)/2))   !note:fortran integer division -> truncation
  real*8 :: theta,r
  
  do i=1,(3*natoms+1)/2
    theta=rand()
    r=rand()
    theta=2*pi*rand()
    r=sqrt(-2*log(rand())*var)
    gaussian_random(2*i-1)=r*dcos(theta) +mu
    gaussian_random(2*i)=r*dsin(theta) +mu
  end do

  return
endfunction


endmodule
