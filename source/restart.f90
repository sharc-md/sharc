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

!> # Module RESTART
!>
!> \author Sebastian Mai
!> \date 27.02.2015
!>
!>                   modified 11.13.2019 by Yinan Shu
!>                       add new variables (see definitions.F90)
!>
!> This module provides all subroutines necessary for the restart feature
!> These are:
!> - creating a restart/ directory where the interfaces put the wavefunction files, etc.
!> - writing static restart information to restart.ctrl
!> - writing timestep-dependent information to restart.traj
!> - reading restart.ctrl and restart.traj to initialize all arrays
module restart
 contains

!> creates a directory $CWD/restart/
!> Via QM.in the interfaces are given the path to this directory
!> where they have to save all files necessary for restart
  subroutine mkdir_restart(ctrl)
    use definitions
    implicit none
    type(ctrl_type) :: ctrl
    character*1023 :: filename
    logical :: exists

    filename=trim(ctrl%cwd)//'/restart'
    inquire(file=filename, exist=exists)
    if (.not.exists) then
      call system('mkdir '//filename)
    endif

  endsubroutine

!> writes the static ctrl variable to restart.ctrl
  subroutine write_restart_ctrl(u,ctrl)
    use definitions
    use matrix
    implicit none
    integer :: u
    type(ctrl_type) :: ctrl

    integer :: imult, istate, ilaser, iatom, iconstr, ipair

    ! the ctrl restart file is only written once at the beginning to avoid writing the laser field
    ! each timestep
    open(unit=u, file='restart.ctrl',status='replace',action='write')

    ! write the ctrl compound, including all entries (see definition)
    ! even the restart keyword is written for completeness
    ! however, the read_restart routine must set ctrl%restart afterwards to 1

    ! printlevel is not part of ctrl, but write it anyways
    write(u,*) printlevel

    ! write ctrl
    write(u,'(A)') trim(ctrl%cwd)
    write(u,*) ctrl%output_version
    write(u,*) ctrl%compat_mode
    write(u,*) ctrl%natom, '! natom'
    write(u,*) ctrl%maxmult
    write(u,*) (ctrl%nstates_m(imult),imult=1,ctrl%maxmult)
    write(u,*) ctrl%nstates
    write(u,*) ctrl%nsteps
    write(u,*) ctrl%nsubsteps
    write(u,*) ctrl%tmax
    write(u,*) ctrl%dtstep_min
    write(u,*) ctrl%dtstep_max
    write(u,*) ctrl%dtstep
    write(u,*) ctrl%dtstep_old
    write(u,*) ctrl%ezero
    write(u,*) ctrl%convthre
    write(u,*) ctrl%scalingfactor
    write(u,*) ctrl%soc_scaling
    write(u,*) ctrl%eselect_grad
    write(u,*) ctrl%eselect_nac
    write(u,*) ctrl%eselect_dmgrad
    write(u,*) ctrl%dampeddyn
    write(u,*) ctrl%decoherence_alpha
    write(u,*) ctrl%force_hop_to_gs
    write(u,*) (ctrl%actstates_s(istate),istate=1,ctrl%nstates)
    write(u,*) (ctrl%output_steps_stride(istate),istate=1,3)
    write(u,*) (ctrl%output_steps_limits(istate),istate=1,3)
    write(u,*) ctrl%restart
    write(u,*) ctrl%restart_rerun_last_qm_step
    write(u,*) ctrl%method
    write(u,*) ctrl%integrator
    write(u,*) ctrl%staterep
    write(u,*) ctrl%initcoeff
    write(u,*) ctrl%laser, '! laser'
    write(u,*) ctrl%coupling
    write(u,*) ctrl%ktdc_method
    write(u,*) ctrl%kmatrix_method
    write(u,*) ctrl%eeom
    write(u,*) ctrl%neom
    write(u,*) ctrl%neom_rep
    write(u,*) ctrl%surf
    write(u,*) ctrl%decoherence
    write(u,*) ctrl%ekincorrect
    write(u,*) ctrl%reflect_frustrated
    write(u,*) ctrl%time_uncertainty
    write(u,*) ctrl%gradcorrect
    write(u,*) ctrl%dipolegrad, '! dipolegrad'
    write(u,*) ctrl%nac_projection
    write(u,*) ctrl%zpe_correction
    write(u,*) ctrl%lpzpe_scheme
    write(u,*) ctrl%lpzpe_nah
    write(u,*) ctrl%lpzpe_nbc
    write(u,*) (ctrl%lpzpe_ah(ipair,1),ipair=1,ctrl%lpzpe_nah)
    write(u,*) (ctrl%lpzpe_ah(ipair,2),ipair=1,ctrl%lpzpe_nah)
    write(u,*) (ctrl%lpzpe_bc(ipair,1),ipair=1,ctrl%lpzpe_nbc)
    write(u,*) (ctrl%lpzpe_bc(ipair,2),ipair=1,ctrl%lpzpe_nbc)
    write(u,*) (ctrl%lpzpe_ke_zpe_ah(ipair),ipair=1,ctrl%lpzpe_nah)
    write(u,*) (ctrl%lpzpe_ke_zpe_bc(ipair),ipair=1,ctrl%lpzpe_nbc)
    write(u,*) ctrl%ke_threshold
    write(u,*) ctrl%t_cycle
    write(u,*) ctrl%t_check
    write(u,*) ctrl%pointer_basis
    write(u,*) ctrl%pointer_maxiter
    write(u,*) ctrl%calc_soc
    write(u,*) ctrl%calc_grad
    write(u,*) ctrl%calc_overlap
    write(u,*) ctrl%calc_nacdt
    write(u,*) ctrl%calc_nacdr
    write(u,*) ctrl%calc_effectivenac
    write(u,*) ctrl%calc_dipolegrad, '!calc_dipolegrad'
    write(u,*) ctrl%calc_second
    write(u,*) ctrl%calc_phases
    write(u,*) ctrl%killafter
    write(u,*) ctrl%ionization
    write(u,*) ctrl%theodore
    write(u,*) ctrl%track_phase
    write(u,*) ctrl%track_phase_at_zero
    write(u,*) ctrl%hopping_procedure
    write(u,*) ctrl%switching_procedure
    write(u,*) ctrl%army_ants
    write(u,*) ctrl%output_format, '! output_format'

    ! thermostat
    write(u,*) ctrl%thermostat
    if (ctrl%thermostat/=0) then
      write(u,*) ctrl%temperature
      if (ctrl%thermostat==1) then   !Langevin: only 1 Thermostat constant
        write(u,*) ctrl%thermostat_const(1)
      endif
      write(u,*) ctrl%restart_thermostat_random
    endif

    ! constraints
    write(u,*) ctrl%do_constraints
    write(u,*) ctrl%constraints_tol
    write(u,*) ctrl%n_constraints
    if (ctrl%do_constraints==1) then
      do iconstr=1,ctrl%n_constraints
        write(u,*) ctrl%constraints_ca(iconstr,1),ctrl%constraints_ca(iconstr,2)
      enddo
      do iconstr=1,ctrl%n_constraints
        write(u,*) ctrl%constraints_dist_c(iconstr)
      enddo
    endif

    ! write the laser field
    if (ctrl%laser==2) then
      write(u,*) ctrl%laser_bandwidth
      write(u,*) ctrl%nlasers
      call vec3write(ctrl%nsteps*ctrl%nsubsteps+1, ctrl%laserfield_td, u, 'Laser field','ES24.16E3')
      do ilaser=1,ctrl%nlasers
        call vecwrite(ctrl%nsteps*ctrl%nsubsteps+1, ctrl%laserenergy_tl(:,ilaser), u, 'Laser Energy','ES24.16E3')
      enddo
    endif
    
    write(u,*) ctrl%write_soc
    write(u,*) ctrl%write_overlap
    write(u,*) ctrl%write_grad
    write(u,*) ctrl%write_nacdr
    write(u,*) ctrl%write_property1d
    write(u,*) ctrl%write_property2d
    write(u,*) ctrl%n_property1d
    write(u,*) ctrl%n_property2d
    
    do iatom=1,ctrl%natom
      write(u,*) ctrl%atommask_a(iatom)
    enddo
    
    close(u)

  endsubroutine

! =========================================================== !

!> writes the traj variable to restart.traj
!> while restart.ctrl is written only once at initialization,
!> restart.traj is rewritten after each timestep
  subroutine write_restart_traj(u,ctrl,traj)
    use definitions
    use matrix
    implicit none
    integer :: u
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl

    integer :: iatom, i,j,k
    character(8000) :: string

    ! the restart file has to be opened each time, since its contents have to be replaced
    open(unit=u, file='restart.traj',status='replace',action='write')
    ! TODO perhaps replace this with a rewind command

    ! write everything
    write(u,*) traj%RNGseed
    write(u,*) traj%RNGseed_thermostat
    write(u,*) traj%traj_hash
    write(u,*) traj%state_MCH
    write(u,*) traj%state_MCH_old
    write(u,*) traj%state_diag
    write(u,*) traj%state_diag_old
    write(u,*) traj%state_diag_frust
    if (ctrl%zpe_correction==1) then
      write(u,*) (traj%state_pumping_s(i),i=1,ctrl%nstates)
      write(u,*) (traj%pumping_status_s(i),i=1,ctrl%nstates)
    else if (ctrl%zpe_correction==2) then 
      write(u,*) (traj%lpzpe_ke_ah(i),i=1,ctrl%lpzpe_nah)
      write(u,*) (traj%lpzpe_ke_bc(i),i=1,ctrl%lpzpe_nbc)
      write(u,*) traj%lpzpe_cycle
      write(u,*) traj%lpzpe_iter_incycle
      write(u,*) traj%lpzpe_starttime
      write(u,*) traj%lpzpe_endtime
      write(u,*) traj%in_cycle
    endif
    if (ctrl%time_uncertainty==1) then
      write(u,*) (traj%uncertainty_time_s(i),i=1,ctrl%nstates)
      write(u,*) (traj%allowed_hop_s(i),i=1,ctrl%nstates)
      write(u,*) (traj%allowed_time_s(i),i=1,ctrl%nstates)
      write(u,*) traj%uncertainty_time_frust
      write(u,*) traj%incident_time
      write(u,*) traj%in_time_uncertainty
      write(u,*) traj%tu_backpropagation
      write(u,*) traj%travelling_state
      write(u,*) traj%target_state
    endif
    write(u,*) traj%step
    write(u,*) traj%microtime
    write(u,*) traj%Ekin
    write(u,*) traj%Epot
    write(u,*) traj%Etot
    write(u,*) traj%Ekin_old
    write(u,*) traj%Epot_old
    write(u,*) traj%Etot_old
    write(u,*) traj%time_start
    write(u,*) traj%time_last
    write(u,*) traj%time_step
    write(u,*) traj%kind_of_jump
    write(u,*) traj%steps_in_gs
    write(u,*) (traj%ncids(i),i=1,10)
    write(u,*) traj%nc_index

    ! write the arrays
    write(u,*) (traj%atomicnumber_a(iatom),iatom=1,ctrl%natom)
    write(u,'(99999(A3,1X))') (traj%element_a(iatom),iatom=1,ctrl%natom)
    write(u,*) (traj%mass_a(iatom),iatom=1,ctrl%natom)
    call vec3write(ctrl%natom, traj%geom_ad,  u, 'Geometry','ES24.16E3')
    call vec3write(ctrl%natom, traj%veloc_ad, u, 'Velocity','ES24.16E3')
    call vec3write(ctrl%natom, traj%veloc_old_ad, u, 'Velocity Old','ES24.16E3')
    call vec3write(ctrl%natom, traj%veloc_app_ad, u, 'Velocity Old','ES24.16E3')
    call vec3write(ctrl%natom, traj%accel_ad, u, 'Acceleration','ES24.16E3')
    if ((ctrl%method==1 .and. ctrl%nac_projection==1) .or. &
      &(ctrl%method==0 .and. (ctrl%ekincorrect==2 .or. ctrl%ekincorrect==5 .or. ctrl%ekincorrect==6 .or. ctrl%ekincorrect==8)) .or. &
      &(ctrl%method==0 .and. (ctrl%reflect_frustrated==2 .or. ctrl%reflect_frustrated==5 .or. ctrl%reflect_frustrated==6 .or. &
      &ctrl%reflect_frustrated==8 .or. ctrl%reflect_frustrated==92 .or. ctrl%reflect_frustrated==95 .or. &
      &ctrl%reflect_frustrated==96 .or. ctrl%reflect_frustrated==98))) then 
      call matwrite(3*ctrl%natom, traj%trans_rot_P, u, 'trans_rot_P','ES24.16E3')
    endif 

    call matwrite(ctrl%nstates, traj%H_MCH_ss,     u, 'H_MCH_ss','ES24.16E3')
    call matwrite(ctrl%nstates, traj%dH_MCH_ss,    u, 'dH_MCH_ss','ES24.16E3')
    call matwrite(ctrl%nstates, traj%dH_MCH_old_ss,u, 'dH_MCH_old_ss','ES24.16E3')
    call matwrite(ctrl%nstates, traj%H_MCH_old_ss, u, 'H_MCH_old_ss','ES24.16E3')
    call matwrite(ctrl%nstates, traj%H_MCH_old2_ss, u, 'H_MCH_old2_ss','ES24.16E3')
    call matwrite(ctrl%nstates, traj%H_MCH_old3_ss, u, 'H_MCH_old3_ss','ES24.16E3')
    call matwrite(ctrl%nstates, traj%H_diag_ss,    u, 'H_diag_ss','ES24.16E3')
    call matwrite(ctrl%nstates, traj%H_diag_old_ss,    u, 'H_diag_old_ss','ES24.16E3')
    call matwrite(ctrl%nstates, traj%H_diag_old2_ss,    u, 'H_diag_old2_ss','ES24.16E3')
    call matwrite(ctrl%nstates, traj%H_diag_old3_ss,    u, 'H_diag_old3_ss','ES24.16E3')
    if (ctrl%zpe_correction==1) then
      call vecwrite(ctrl%nstates, traj%Evib_local_s, u, 'Evib_local_s','ES24.16E3')
      call vecwrite(ctrl%nstates, traj%Ezpe_local_s, u, 'Ezpe_local_s','ES24.16E3')
    endif
    call matwrite(ctrl%nstates, traj%U_ss,         u, 'U_ss','ES24.16E3')
    call matwrite(ctrl%nstates, traj%U_old_ss,     u, 'U_old_ss','ES24.16E3')
    call matwrite(ctrl%nstates, traj%NACdt_ss,     u, 'NACdt_ss','ES24.16E3')
    call matwrite(ctrl%nstates, traj%NACdt_old_ss, u, 'NACdt_old_ss','ES24.16E3')
    call matwrite(ctrl%nstates, traj%overlaps_ss,  u, 'overlaps_ss','ES24.16E3')
    call matwrite(ctrl%nstates, traj%DM_ssd(:,:,1),  u, 'DM_ssd(x)','ES24.16E3')
    call matwrite(ctrl%nstates, traj%DM_ssd(:,:,2),  u, 'DM_ssd(y)','ES24.16E3')
    call matwrite(ctrl%nstates, traj%DM_ssd(:,:,3),  u, 'DM_ssd(z)','ES24.16E3')
    call matwrite(ctrl%nstates, traj%DM_old_ssd(:,:,1),  u, 'DM_old_ssd(x)','ES24.16E3')
    call matwrite(ctrl%nstates, traj%DM_old_ssd(:,:,2),  u, 'DM_old_ssd(y)','ES24.16E3')
    call matwrite(ctrl%nstates, traj%DM_old_ssd(:,:,3),  u, 'DM_old_ssd(z)','ES24.16E3')
    call matwrite(ctrl%nstates, traj%DM_print_ssd(:,:,1),  u, 'DM_print_ssd(x)','ES24.16E3')
    call matwrite(ctrl%nstates, traj%DM_print_ssd(:,:,2),  u, 'DM_print_ssd(y)','ES24.16E3')
    call matwrite(ctrl%nstates, traj%DM_print_ssd(:,:,3),  u, 'DM_print_ssd(z)','ES24.16E3')
!     call matwrite(ctrl%nstates, traj%Property_ss,  u, 'Property_ss','ES24.16E3')
    call matwrite(ctrl%nstates, traj%Rtotal_ss,    u, 'Rtotal_ss','ES24.16E3')
    call matwrite(ctrl%nstates, traj%RDtotal_ss,    u, 'RDtotal_ss','ES24.16E3')
    call matwrite(ctrl%nstates, traj%Dtotal_ss,    u, 'Dtotal_ss','ES24.16E3')
    call matwrite(ctrl%nstates, traj%dendt_MCH_ss,    u, 'dendt_MCH_ss','ES24.16E3')
    call matwrite(ctrl%nstates, traj%dendt_diag_ss,    u, 'dendt_diag_ss','ES24.16E3')
    call matwrite(ctrl%nstates, traj%Kmatrix_MCH_ss,    u, 'Kmatrix_MCH_ss','ES24.16E3')
    call matwrite(ctrl%nstates, traj%Kmatrix_diag_ss,    u,  'Kmatrix_diag_ss','ES24.16E3')
    call vecwrite(ctrl%nstates, traj%phases_s, u, 'phases_s','ES24.16E3')
    call vecwrite(ctrl%nstates, traj%phases_old_s, u, 'phases_old_s','ES24.16E3')
    call vecwrite(ctrl%nstates, traj%hopprob_s, u, 'hopprob_s','ES24.16E3')
    call vecwrite(ctrl%nstates, traj%switchprob_s, u, 'switchingprob_s','ES24.16E3')
    call vec3write(ctrl%nstates, traj%gpreprob_s3,  u, 'gprepro_s3','ES24.16E3')
    call vec3write(ctrl%nstates, traj%preprob_s3,  u, 'prepro_s3','ES24.16E3')
    call vec3write(ctrl%nstates, traj%preprob_old_s3,  u, 'prepro_old_s3','ES24.16E3')
    call vecwrite(ctrl%nstates, traj%decotime_diag_s, u, 'decotime_diag_s','ES24.16E3')
    call vecwrite(ctrl%nstates, traj%decotime_diag_old_s, u, 'decotime_diag_old_s','ES24.16E3')
    call vecwrite(ctrl%nstates, traj%decotime_MCH_s, u, 'decotime_MCH_s','ES24.16E3')
    call vecwrite(ctrl%nstates, traj%decotime_MCH_old_s, u, 'decotime_MCH_old_s','ES24.16E3')
    do i=1,ctrl%nstates
      write(string,'(A45,I3,1X,I3)') 'svec_MCH_sad',i
      call vec3write(ctrl%natom,traj%svec_MCH_sad(i,:,:),u,trim(string),'ES24.16E3')
    enddo
    do i=1,ctrl%nstates
      write(string,'(A45,I3,1X,I3)') 'svec_diag_sad',i
      call vec3write(ctrl%natom,traj%svec_diag_sad(i,:,:),u,trim(string),'ES24.16E3')
    enddo
    do i=1,ctrl%nstates
      write(string,'(A45,I3,1X,I3)') 'psvec_MCH_sad',i
      call vec3write(ctrl%natom,traj%psvec_MCH_sad(i,:,:),u,trim(string),'ES24.16E3')
    enddo
    do i=1,ctrl%nstates
      write(string,'(A45,I3,1X,I3)') 'psvec_diag_sad',i
      call vec3write(ctrl%natom,traj%psvec_diag_sad(i,:,:),u,trim(string),'ES24.16E3')
    enddo

    write(u,*) traj%randnum
    write(u,*) traj%randnum2
    write(u,*) traj%dmag
    write(u,*) traj%dmag1
    write(u,*) traj%dmag2

    if (ctrl%calc_dipolegrad>-1) then
      do i=1,ctrl%nstates
        do j=1,ctrl%nstates
          do k=1,3
            write(string,'(A45,I3,1X,I3,1X,I3)') 'DMgrad_ssdad',i,j,k
            call vec3write(ctrl%natom,traj%DMgrad_ssdad(i,j,k,:,:),u,trim(string),'ES24.16E3')
          enddo
        enddo
      enddo
    endif
    if (ctrl%calc_nacdr>-1) then
      do i=1,ctrl%nstates
        do j=1,ctrl%nstates
          write(string,'(A45,I3,1X,I3)') 'naddr_ssad',i,j
          call vec3write(ctrl%natom,traj%NACdr_ssad(i,j,:,:),u,trim(string),'ES24.16E3')
        enddo
      enddo
      do i=1,ctrl%nstates
        do j=1,ctrl%nstates
          write(string,'(A45,I3,1X,I3)') 'naddr_old_ssad',i,j
          call vec3write(ctrl%natom,traj%NACdr_old_ssad(i,j,:,:),u,trim(string),'ES24.16E3')
        enddo
      enddo
      do i=1,ctrl%nstates
        do j=1,ctrl%nstates
          write(string,'(A45,I3,1X,I3)') 'naddr_diag_ssad',i,j
          call vec3write(ctrl%natom,traj%NACdr_diag_ssad(i,j,:,:),u,trim(string),'ES24.16E3')
        enddo
      enddo
      if (ctrl%nac_projection==1) then
        do i=1,ctrl%nstates
          do j=1,ctrl%nstates
            write(string,'(A45,I3,1X,I3)') 'pnaddr_MCH_ssad',i,j
            call vec3write(ctrl%natom,traj%pNACdr_MCH_ssad(i,j,:,:),u,trim(string),'ES24.16E3')
          enddo
        enddo
        do i=1,ctrl%nstates
          do j=1,ctrl%nstates
            write(string,'(A45,I3,1X,I3)') 'pnaddr_diag_ssad',i,j
            call vec3write(ctrl%natom,traj%pNACdr_diag_ssad(i,j,:,:),u,trim(string),'ES24.16E3')
          enddo
        enddo
      endif
    endif
    do i=1,ctrl%nstates
      write(string,'(A45,I3,1X,I3)') 'grad_MCH_sad',i
      call vec3write(ctrl%natom,traj%grad_MCH_sad(i,:,:),u,trim(string),'ES24.16E3')
    enddo
    do i=1,ctrl%nstates
      write(string,'(A45,I3,1X,I3)') 'grad_MCH_old_sad',i
      call vec3write(ctrl%natom,traj%grad_MCH_old_sad(i,:,:),u,trim(string),'ES24.16E3')
    enddo
    if (ctrl%zpe_correction==1) then
      do i=1,ctrl%nstates
        write(string,'(A45,I3,1X,I3)') 'grad_MCH_old2_sad',i
        call vec3write(ctrl%natom,traj%grad_MCH_old2_sad(i,:,:),u,trim(string),'ES24.16E3')
      enddo
      do i=1,ctrl%nstates
        write(string,'(A45,I3,1X,I3)') 'grad_Ezpe_local_sad',i
        call vec3write(ctrl%natom,traj%grad_Ezpe_local_sad(i,:,:),u,trim(string),'ES24.16E3')
      enddo
    endif   
    do i=1,ctrl%nstates
      do j=1,ctrl%nstates
        write(string,'(A45,I3,1X,I3)') 'hopping_direction_ssad',i,j
        call vec3write(ctrl%natom,traj%hopping_direction_ssad(i,j,:,:),u,trim(string),'ES24.16E3')
      enddo
    enddo
    do i=1,ctrl%nstates
      do j=1,ctrl%nstates
        write(string,'(A45,I3,1X,I3)') 'frustrated_hop_vec_ssad',i,j
        call vec3write(ctrl%natom,traj%frustrated_hop_vec_ssad(i,j,:,:),u,trim(string),'ES24.16E3')
      enddo
    enddo
    do i=1,ctrl%nstates
      do j=1,ctrl%nstates
        write(string,'(A45,I3,1X,I3)') 'Gmatrix_ssad',i,j
        call vec3write(ctrl%natom,traj%Gmatrix_ssad(i,j,:,:),u,trim(string),'ES24.16E3')
      enddo
    enddo
    if (ctrl%calc_effectivenac==1) then
      do i=1,ctrl%nstates
        do j=1,ctrl%nstates
          write(string,'(A45,I3,1X,I3)') 'NACGV_MCH_ssad',i,j
          call vec3write(ctrl%natom,traj%NACGV_MCH_ssad(i,j,:,:),u,trim(string),'ES24.16E3')
        enddo
      enddo
      do i=1,ctrl%nstates
        do j=1,ctrl%nstates
          write(string,'(A45,I3,1X,I3)') 'NACGV_diag_ssad',i,j
          call vec3write(ctrl%natom,traj%NACGV_diag_ssad(i,j,:,:),u,trim(string),'ES24.16E3')
        enddo
      enddo
      if (ctrl%nac_projection==1) then
      do i=1,ctrl%nstates
        do j=1,ctrl%nstates
          write(string,'(A45,I3,1X,I3)') 'pNACGV_MCH_ssad',i,j
          call vec3write(ctrl%natom,traj%pNACGV_MCH_ssad(i,j,:,:),u,trim(string),'ES24.16E3')
        enddo
      enddo
      do i=1,ctrl%nstates
        do j=1,ctrl%nstates
          write(string,'(A45,I3,1X,I3)') 'pNACGV_diag_ssad',i,j
          call vec3write(ctrl%natom,traj%pNACGV_diag_ssad(i,j,:,:),u,trim(string),'ES24.16E3')
        enddo
      enddo
      endif
    endif

    call vec3write(ctrl%natom,traj%grad_ad(:,:),u,'grad_ad','ES24.16E3')
    call vec3write(ctrl%natom,traj%decograd_ad(:,:),u,'decograd_ad','ES24.16E3')

    call vecwrite(ctrl%nstates, traj%coeff_diag_s, u, 'coeff_diag_s','ES24.16E3')
    call vecwrite(ctrl%nstates, traj%coeff_diag_old_s, u, 'coeff_diag_old_s','ES24.16E3')
    call vecwrite(ctrl%nstates, traj%coeff_mch_s, u, 'coeff_mch_s','ES24.16E3')
    call vecwrite(ctrl%nstates, traj%coeff_mch_old_s, u, 'coeff_mch_old_s','ES24.16E3')
    call vecwrite(ctrl%nstates, traj%ccoeff_diag_s, u, 'ccoeff_diag_s','ES24.16E3')
    call vecwrite(ctrl%nstates, traj%ccoeff_diag_old_s, u, 'ccoeff_diag_old_s','ES24.16E3')
    call vecwrite(ctrl%nstates, traj%ccoeff_mch_s, u, 'ccoeff_mch_s','ES24.16E3')
    if (ctrl%zpe_correction==1) then
      call vec2write(ctrl%nstates, traj%coeff_zpe_s2, u, 'coeff_zpe_s2','ES24.16E3')
    endif
    if (ctrl%integrator==0) then
      call vecwrite(ctrl%nstates, traj%gRcoeff_mch_s, u, 'gRcoeff_mch_s','ES24.16E3')
      call vecwrite(ctrl%nstates, traj%gIcoeff_mch_s, u, 'gIcoeff_mch_s','ES24.16E3')
      call vecwrite(ctrl%nstates, traj%gRcoeff_mch_s, u, 'gRccoeff_mch_s','ES24.16E3')
      call vecwrite(ctrl%nstates, traj%gIcoeff_mch_s, u, 'gIccoeff_mch_s','ES24.16E3')
      call vecwrite(ctrl%nstates, traj%ephase_s, u, 'ephase_s','ES24.16E3')
      call vecwrite(ctrl%nstates, traj%gephase_s, u, 'gephase_s','ES24.16E3')
    endif

    write(u,*) (traj%selg_s(i),i=1,ctrl%nstates)
    do i=1,ctrl%nstates
      write(u,*) (traj%selt_ss(i,j),j=1,ctrl%nstates)
    enddo
    if (ctrl%calc_dipolegrad>-1) then
      do i=1,ctrl%nstates
        write(u,*) (traj%seldm_ss(i,j),j=1,ctrl%nstates)
      enddo
    endif
    write(u,*) traj%phases_found

    call vecwrite(ctrl%n_property1d, traj%Property1d_labels_y, u, 'Property1d_labels_y','A40')
    call vecwrite(ctrl%n_property2d, traj%Property2d_labels_x, u, 'Property2d_labels_x','A40')
    do i=1,ctrl%n_property1d
      write(string,'(A45,I3,1X,I3)') 'Property1d_ys',i
      call vecwrite(ctrl%nstates, traj%Property1d_ys(i,:), u, string,'ES24.16E3')
    enddo
    do i=1,ctrl%n_property2d
      write(string,'(A45,I3,1X,I3)') 'Property2d_xss',i
      call matwrite(ctrl%nstates, traj%Property2d_xss(i,:,:), u, string,'ES24.16E3')
    enddo

    ! save restart info for the auxilliary trajectories
    if (ctrl%decoherence==2) then
      do i=1,ctrl%nstates
        write(u,*) traj%auxtrajs_s(i)%istate
        write(u,*) traj%auxtrajs_s(i)%rate1
        write(u,*) traj%auxtrajs_s(i)%rate2
        call vec3write(ctrl%natom, traj%auxtrajs_s(i)%geom_ad,   u, 'AuxGeometry','ES24.16E3')
        call vec3write(ctrl%natom, traj%auxtrajs_s(i)%veloc_ad,  u, 'AuxVeloc','ES24.16E3')
        call vec3write(ctrl%natom, traj%auxtrajs_s(i)%accel_ad,  u, 'AuxAccel','ES24.16E3')
        call vec3write(ctrl%natom, traj%auxtrajs_s(i)%grad_ad,   u, 'AuxGrad','ES24.16E3')
        call vec3write(ctrl%natom, traj%auxtrajs_s(i)%geom_tmp_ad,   u, 'AuxGeometry2','ES24.16E3')
        call vec3write(ctrl%natom, traj%auxtrajs_s(i)%veloc_tmp_ad,  u, 'AuxVeloc2','ES24.16E3')
      enddo
    endif

    ! army ants sampling
    if (ctrl%army_ants==1) then 
      write(u,*) traj%army_ants_weight
      write(u,*) traj%randnum_branching
      write(u,*) traj%branching_likehood
    endif

    ! pointer basis optimization 
    if (ctrl%pointer_basis==2) then
      write(u,*) traj%entropy
      write(u,*) traj%entropy_old
      write(u,*) traj%linear_entropy
      write(u,*) traj%linear_entropy_old
      write(u,*) traj%entropy_grad
      write(u,*) traj%linear_entropy_grad
      call matwrite(ctrl%nstates, traj%U_pointer_ss,  u, 'U_ss','ES24.16E3')
      call matwrite(ctrl%nstates, traj%U_pointer_old_ss,  u, 'U_old_ss','ES24.16E3')
    endif

    ! Trajectory consistency
    write(u,*) traj%discrepancy

    close(u)

  endsubroutine

! =========================================================== !

!> this routine reads all restart information from restart.ctrl and restart.traj
!> and initializes the ctrl and traj compounds
!> it also does:
!> - allocation of all arrays
!> - setting ctrl%restart to .true.
!> - fast-forwards the random number generator so that the 
!>     restarted trajectory uses the same random number sequence as if it was not restarted
!> - sets steps_in_gs correctly
!> - sets the wallclock timing 
  subroutine read_restart(u_ctrl,u_traj,ctrl,traj)
    use definitions
    use matrix
    use misc
    use decoherence_afssh
    implicit none
    integer :: u_ctrl,u_traj
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl

    ! allocation of trajectory is in read_restart_ctrl
    call read_restart_ctrl(u_ctrl,u_traj,ctrl,traj)
    call read_restart_traj(u_ctrl,u_traj,ctrl,traj)
    call read_restart_miscllaneous(u_ctrl,u_traj,ctrl,traj)

  endsubroutine

! =========================================================== !

!> this routine reads all restart information from restart.ctrl
!> and initializes the ctrl and traj compounds
!> it also does:
!> - allocation of all arrays including traj 
!> - setting ctrl%restart to .true.
  subroutine read_restart_ctrl(u_ctrl,u_traj,ctrl,traj)
    use definitions
    use matrix
    use misc
    use decoherence_afssh
    implicit none
    integer :: iconstr
    integer :: u_ctrl,u_traj
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl

    integer :: imult, iatom, i,j,k, istate, ilaser, ipair
    character(8000) :: string

    open(unit=u_ctrl, file='restart.ctrl',status='old',action='read')

    ! read the ctrl compound, including all entries (see definition)
    ! even the restart keyword is written for completeness
    ! however, the read_restart routine must set ctrl%restart afterwards to 1

    ! printlevel is not part of ctrl, but write it anyways
    read(u_ctrl,*) printlevel

    if (printlevel>0) then
      write(u_log,'(A)')      '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<============================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(u_log,'(A)')      '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    Initializing Restart    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(u_log,'(A)')      '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<============================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(u_log,*)
    endif

    if (printlevel>0) then
      write(u_log,*) '============================================================='
      write(u_log,*) '                          Allocation'
      write(u_log,*) '============================================================='
      write(u_log,*)
    endif

    ! read ctrl
    read(u_ctrl,'(A)') ctrl%cwd
    call getcwd(ctrl%cwd)
    read(u_ctrl,*) ctrl%output_version
    read(u_ctrl,*) ctrl%compat_mode
    read(u_ctrl,*) ctrl%natom
    read(u_ctrl,*) ctrl%maxmult
    allocate( ctrl%nstates_m(ctrl%maxmult) )
    read(u_ctrl,*) (ctrl%nstates_m(imult),imult=1,ctrl%maxmult)
    read(u_ctrl,*) ctrl%nstates
    read(u_ctrl,*) ctrl%nsteps
    read(u_ctrl,*) ctrl%nsubsteps
    read(u_ctrl,*) ctrl%tmax
    read(u_ctrl,*) ctrl%dtstep_min
    read(u_ctrl,*) ctrl%dtstep_max
    read(u_ctrl,*) ctrl%dtstep
    read(u_ctrl,*) ctrl%dtstep_old
    read(u_ctrl,*) ctrl%ezero
    read(u_ctrl,*) ctrl%convthre
    read(u_ctrl,*) ctrl%scalingfactor
    read(u_ctrl,*) ctrl%soc_scaling
    read(u_ctrl,*) ctrl%eselect_grad
    read(u_ctrl,*) ctrl%eselect_nac
    read(u_ctrl,*) ctrl%eselect_dmgrad
    read(u_ctrl,*) ctrl%dampeddyn
    read(u_ctrl,*) ctrl%decoherence_alpha
    read(u_ctrl,*) ctrl%force_hop_to_gs
    allocate( ctrl%actstates_s(ctrl%nstates) )
    read(u_ctrl,*) (ctrl%actstates_s(istate),istate=1,ctrl%nstates)
    read(u_ctrl,*) (ctrl%output_steps_stride(istate),istate=1,3)
    read(u_ctrl,*) (ctrl%output_steps_limits(istate),istate=1,3)
    read(u_ctrl,*) ctrl%restart
    read(u_ctrl,*) ctrl%restart_rerun_last_qm_step
    read(u_ctrl,*) ctrl%method
    read(u_ctrl,*) ctrl%integrator
    read(u_ctrl,*) ctrl%staterep
    read(u_ctrl,*) ctrl%initcoeff
    read(u_ctrl,*) ctrl%laser
    read(u_ctrl,*) ctrl%coupling
    read(u_ctrl,*) ctrl%ktdc_method
    read(u_ctrl,*) ctrl%kmatrix_method
    read(u_ctrl,*) ctrl%eeom
    read(u_ctrl,*) ctrl%neom
    read(u_ctrl,*) ctrl%neom_rep
    read(u_ctrl,*) ctrl%surf
    read(u_ctrl,*) ctrl%decoherence
    read(u_ctrl,*) ctrl%ekincorrect
    read(u_ctrl,*) ctrl%reflect_frustrated
    read(u_ctrl,*) ctrl%time_uncertainty
    read(u_ctrl,*) ctrl%gradcorrect
    read(u_ctrl,*) ctrl%dipolegrad
    read(u_ctrl,*) ctrl%nac_projection
    read(u_ctrl,*) ctrl%zpe_correction
    read(u_ctrl,*) ctrl%lpzpe_scheme
    read(u_ctrl,*) ctrl%lpzpe_nah
    read(u_ctrl,*) ctrl%lpzpe_nbc
    allocate( ctrl%lpzpe_ah(ctrl%lpzpe_nah,2) )
    read(u_ctrl,*) (ctrl%lpzpe_ah(ipair,1),ipair=1,ctrl%lpzpe_nah)
    read(u_ctrl,*) (ctrl%lpzpe_ah(ipair,2),ipair=1,ctrl%lpzpe_nah)
    allocate( ctrl%lpzpe_bc(ctrl%lpzpe_nbc,2) )
    read(u_ctrl,*) (ctrl%lpzpe_bc(ipair,1),ipair=1,ctrl%lpzpe_nbc)
    read(u_ctrl,*) (ctrl%lpzpe_bc(ipair,2),ipair=1,ctrl%lpzpe_nbc)
    allocate( ctrl%lpzpe_ke_zpe_ah(ctrl%lpzpe_nah) )
    read(u_ctrl,*) (ctrl%lpzpe_ke_zpe_ah(ipair),ipair=1,ctrl%lpzpe_nah)
    allocate( ctrl%lpzpe_ke_zpe_bc(ctrl%lpzpe_nbc) )
    read(u_ctrl,*) (ctrl%lpzpe_ke_zpe_bc(ipair),ipair=1,ctrl%lpzpe_nbc)
    read(u_ctrl,*) ctrl%ke_threshold
    read(u_ctrl,*) ctrl%t_cycle
    read(u_ctrl,*) ctrl%t_check
    read(u_ctrl,*) ctrl%pointer_basis
    read(u_ctrl,*) ctrl%pointer_maxiter
    read(u_ctrl,*) ctrl%calc_soc
    read(u_ctrl,*) ctrl%calc_grad
    read(u_ctrl,*) ctrl%calc_overlap
    read(u_ctrl,*) ctrl%calc_nacdt
    read(u_ctrl,*) ctrl%calc_nacdr
    read(u_ctrl,*) ctrl%calc_effectivenac
    read(u_ctrl,*) ctrl%calc_dipolegrad
    read(u_ctrl,*) ctrl%calc_second
    read(u_ctrl,*) ctrl%calc_phases
    read(u_ctrl,*) ctrl%killafter
    read(u_ctrl,*) ctrl%ionization
    read(u_ctrl,*) ctrl%theodore
    read(u_ctrl,*) ctrl%track_phase
    read(u_ctrl,*) ctrl%track_phase_at_zero
    read(u_ctrl,*) ctrl%hopping_procedure
    read(u_ctrl,*) ctrl%switching_procedure
    read(u_ctrl,*) ctrl%army_ants
    read(u_ctrl,*) ctrl%output_format

    ! thermostat
    read(u_ctrl,*) ctrl%thermostat
    if (ctrl%thermostat/=0) then
      read(u_ctrl,*) ctrl%temperature
      if (ctrl%thermostat==1) then   !Langevin: only 1 Thermostat constant
        allocate(ctrl%thermostat_const(1))
        read(u_ctrl,*) ctrl%thermostat_const(1)
      endif
      read(u_ctrl,*) ctrl%restart_thermostat_random
    endif


    ! constraints
    read(u_ctrl,*) ctrl%do_constraints
    read(u_ctrl,*) ctrl%constraints_tol
    read(u_ctrl,*) ctrl%n_constraints
    allocate( ctrl%constraints_ca(ctrl%n_constraints,2) )
    allocate( ctrl%constraints_dist_c(ctrl%n_constraints) )
    if (ctrl%do_constraints==1) then
      do iconstr=1,ctrl%n_constraints
        read(u_ctrl,*) ctrl%constraints_ca(iconstr,1),ctrl%constraints_ca(iconstr,2)
      enddo
      do iconstr=1,ctrl%n_constraints
        read(u_ctrl,*) ctrl%constraints_dist_c(iconstr)
      enddo
    endif


    ! read the laser field
    ! with an external laser, increasing the simulation time necessitates that
    ! the laserfield in
    ! the control file is enlarged 
    if (ctrl%laser==2) then
      read(u_ctrl,*) ctrl%laser_bandwidth
      read(u_ctrl,*) ctrl%nlasers
      allocate( ctrl%laserfield_td(ctrl%nsteps*ctrl%nsubsteps+1,3) )
      allocate( ctrl%laserenergy_tl(ctrl%nsteps*ctrl%nsubsteps+1,ctrl%nlasers) )
      call vec3read(ctrl%nsteps*ctrl%nsubsteps+1, ctrl%laserfield_td, u_ctrl, string)
      do ilaser=1,ctrl%nlasers
        call vecread(ctrl%nsteps*ctrl%nsubsteps+1, ctrl%laserenergy_tl(:,ilaser), u_ctrl, string)
      enddo
    endif

    read(u_ctrl,*) ctrl%write_soc
    read(u_ctrl,*) ctrl%write_overlap
    read(u_ctrl,*) ctrl%write_grad
    read(u_ctrl,*) ctrl%write_nacdr
    read(u_ctrl,*) ctrl%write_property1d
    read(u_ctrl,*) ctrl%write_property2d
    read(u_ctrl,*) ctrl%n_property1d
    read(u_ctrl,*) ctrl%n_property2d

    allocate( ctrl%atommask_a(ctrl%natom))
    do iatom=1,ctrl%natom
      read(u_ctrl,*) ctrl%atommask_a(iatom)
    enddo

    close(u_ctrl)

    ! -------------------------

    ctrl%restart=.true.

    call allocate_traj(traj,ctrl)
    call additional_allocate_traj(traj,ctrl)

    if (printlevel>1) then
      write(u_log,'(a,1x,i4,1x,a,1x,i4,1x,a)') 'Allocation with nstates=',ctrl%nstates,'and natom=',ctrl%natom,'successful.'
      write(u_log,*)
    endif

    call flush(u_log)
  endsubroutine

! =========================================================== !

!> this routine reads all restart information from restart.traj
  subroutine read_restart_traj(u_ctrl,u_traj,ctrl,traj)
    use definitions
    use matrix
    use misc
    use decoherence_afssh
    implicit none
    integer :: u_ctrl,u_traj
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl

    integer :: imult, iatom, i,j,k, istate, ilaser, ipair
    character(8000) :: string

    open(unit=u_traj, file='restart.traj',status='old',action='read')
  ! read everything
  read(u_traj,*) traj%RNGseed
  read(u_traj,*) traj%RNGseed_thermostat
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
  if ((ctrl%method==1 .and. ctrl%nac_projection==1) .or. &
    &(ctrl%method==0 .and. (ctrl%ekincorrect==2 .or. ctrl%ekincorrect==5 .or. ctrl%ekincorrect==6 .or. ctrl%ekincorrect==8)) .or. &
    &(ctrl%method==0 .and. (ctrl%reflect_frustrated==2 .or. ctrl%reflect_frustrated==5 .or. ctrl%reflect_frustrated==6 .or. &
    &ctrl%reflect_frustrated==8 .or. ctrl%reflect_frustrated==92 .or. ctrl%reflect_frustrated==95 .or. &
    &ctrl%reflect_frustrated==96 .or. ctrl%reflect_frustrated==98))) then
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

! =========================================================== !
!> - fast-forwards the random number generator so that the 
!>     restarted trajectory uses the same random number sequence as if it was not restarted
!> - sets steps_in_gs correctly
!> - sets the wallclock timing 
  subroutine read_restart_miscllaneous(u_ctrl,u_traj,ctrl,traj)
    use definitions
    use matrix
    use misc
    use decoherence_afssh
    implicit none
    integer :: u_ctrl,u_traj
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl

    integer :: imult, iatom, i,j,k, istate, ilaser, ipair
    character(8000) :: string
    real*8 :: dummy_randnum
    integer :: time

    ! call the random number generator until it is in the same status as before
    ! the restart
    call init_random_seed(traj%RNGseed)
    do i=1,2*traj%step
      call random_number(dummy_randnum)
    enddo

    ! since the relaxation check is done after writing the restart file,
    ! add one to the relaxation counter
    traj%steps_in_gs=traj%steps_in_gs+1

    ! set time so that timing is correct
    traj%time_start=time()
    traj%time_last=time()

  endsubroutine

endmodule restart
