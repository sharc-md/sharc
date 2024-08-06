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

!> # Module  tsh_tu
!> tsh_tu, trajectory surface hopping with time uncertainty
!>
!> \author Yinan Shu
!> \date 2.5.2022
!>
!> This module contains subroutines that allow to perform a trajectory surface
!> hopping with time uncertainty. Original publication:
!> A. W. Jasper, S. N. Stechmann, D. G. Truhlar, J. Chem. Phys. 2002, 116,
!> 5424-5431

module tsh_tu
  implicit none

  public time_uncertainty_initialize
  public time_uncertainty_surface_hopping
  public tshtu_time_travelling
  public record_time_travelling_point

  private check_frustration
  private save_time_travelling_point
  private available_e
  private write_tshtu_restart_traj
  private time_uncertainty_hopping
  private read_tshtu_restart_traj

contains

! ===========================================================
! time uncertainty hopping process initialization
  subroutine time_uncertainty_initialize(traj,ctrl)
    use definitions
    use matrix
    use electronic
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl

    ! in regular time uncertainty process
    traj%in_time_uncertainty=0
    ! not do time travelling to a record point
    traj%tu_backpropagation=0
    ! frustrated hop has not happened yet 
    traj%incident_time=0.d0

  endsubroutine

! ===========================================================
! time uncertainty hopping process wrap up
  subroutine time_uncertainty_surface_hopping(traj,ctrl)
    use definitions
    use matrix
    use electronic
    use army_ants
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl

  if (printlevel>2) then
    write(u_log,*) '============================================================='
    write(u_log,*) '       Doing fewest switches with time uncertainty'
    write(u_log,*) '============================================================='
  endif

    ! check frustration for each state
    call check_frustration(traj,ctrl)

    if (traj%in_time_uncertainty==0) then
      ! saving is only done when outside the time uncertainty region 
      call save_time_travelling_point(traj,ctrl)
      ! record time travelling point, and check for hop
      !call record_time_travelling_point(traj,ctrl)
      if (ctrl%army_ants==0) then
        call surface_hopping(traj,ctrl)
      else if (ctrl%army_ants==1) then
        call army_ants_surface_hopping(traj,ctrl)
      endif
      if (traj%kind_of_jump==2) then
         ! frustrated, turn on time uncertainty, record uncertainty time
         traj%in_time_uncertainty=1
         traj%uncertainty_time_frust=traj%uncertainty_time_s(traj%state_diag_frust)
         traj%target_state=traj%state_diag_frust
         traj%incident_time=traj%microtime
      endif 
    endif

    if (traj%in_time_uncertainty==1) then
      ! in time uncertainty process
      call time_uncertainty_hopping(traj,ctrl)
    endif

    if (traj%in_time_uncertainty==2) then
      ! right after a time travelling, do forced hop, or frustrated hop
      if (traj%travelling_state==traj%state_diag) then 
        ! travel back to where frustrated hop happened, and do frustrated hop.
        traj%state_diag_frust=traj%target_state
        traj%kind_of_jump=2
        ! then we go out of time uncertainty process
        traj%in_time_uncertainty=0
        write(u_log,*) "frustrated hop finished, going out of time uncertainty process"
      else if (traj%travelling_state==traj%state_diag_frust) then 
        ! travel back to where a point where frustrated hop is energetically available
        traj%state_diag_old=traj%state_diag
        traj%state_diag=traj%target_state
        traj%kind_of_jump=1
        ! then we go out of time uncertainty process
        traj%in_time_uncertainty=0
        write(u_log,*) "forced hop finished, going out of time uncertainty process"
      endif
    endif


  endsubroutine

! ===========================================================
! check_frustration at every time step
  subroutine check_frustration(traj,ctrl)
    use definitions
    use matrix
    use nuclear
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl

    real*8 :: Emax ! energy threshold for frustrated jumps
    real*8 :: Ekin_masked
    real*8 :: sum_kk, sum_vk ! temp variables for kinetic adjustment
    real*8 :: deltaE_s(ctrl%nstates)

    integer :: istate, ilaser

    ! initialize traj%allowed_hop_s, and deltaE_s
    traj%allowed_hop_s=1
    deltaE_s=0.d0

    ! check each state, whether it's energetically allowed or frustrated
    stateloop: do istate=1,ctrl%nstates

      if (istate.ne.traj%state_diag) then
     
      ! check for laser resonance
      ! if it is resonance, this state is allowed to hop 
      if (ctrl%laser/=0) then
        do ilaser=1,ctrl%nlasers
          deltaE_s(istate)=abs(real(traj%H_diag_ss(istate,istate)-traj%Epot))
          deltaE_s(istate)=abs(deltaE_s(istate)-ctrl%laserenergy_tl(traj%step*ctrl%nsubsteps+1,ilaser))
          if (deltaE_s(istate)<=ctrl%laser_bandwidth) then
            traj%allowed_hop_s(istate)=1
          endif
        enddo
      endif

      ! check for frustration
      select case (ctrl%ekincorrect)
        case (0)    ! no frustrated jumps: all states are allowed
          traj%allowed_hop_s(istate)=1

        case (1)    ! correct along v, full kinetic energy available
          Ekin_masked=Calculate_ekin_masked(ctrl%natom, traj%veloc_ad, traj%mass_a, ctrl%atommask_a)
          Emax=traj%Etot- traj%Ekin + Ekin_masked
          deltaE_s(istate)=Emax-real(traj%H_diag_ss(istate,istate))
          if (deltaE_s(istate)<0.d0) then
            traj%allowed_hop_s(istate)=0
          else
            traj%allowed_hop_s(istate)=1
          endif

        case (2)    ! correct along projected v, vibrational kinetic energy available
          call available_ekin(ctrl%natom,&
          &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag, istate,:,:),&
          &traj%mass_a, sum_kk, sum_vk)
          deltaE_s(istate)=sum_vk**2/(4.d0*sum_kk)+(traj%Etot-traj%Ekin-&
          &real(traj%H_diag_ss(istate,istate)))
          if (deltaE_s(istate)<0.d0) then
            traj%allowed_hop_s(istate)=0
          else
            traj%allowed_hop_s(istate)=1
          endif

        case (3)    ! correct along T, less energy available
          call available_ekin(ctrl%natom,&
          &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag, istate,:,:),&
          &traj%mass_a, sum_kk, sum_vk)
          deltaE_s(istate)=sum_vk**2/(4.d0*sum_kk)+(traj%Etot-traj%Ekin-&
          &real(traj%H_diag_ss(istate,istate)))
          if (deltaE_s(istate)<0.d0) then
            traj%allowed_hop_s(istate)=0
          else
            traj%allowed_hop_s(istate)=1
          endif

        case (4)    ! correct along gradient difference
          call available_ekin(ctrl%natom,&
          &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag, istate,:,:),&
          &traj%mass_a, sum_kk, sum_vk)
          deltaE_s(istate)=sum_vk**2/(4.d0*sum_kk)+(traj%Etot-traj%Ekin-&
          &real(traj%H_diag_ss(istate,istate)))
          if (deltaE_s(istate)<0.d0) then
            traj%allowed_hop_s(istate)=0
          else
            traj%allowed_hop_s(istate)=1
          endif

        case (5)    ! correct along projected NAC/M
          call available_ekin(ctrl%natom,&
          &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag, istate,:,:),&
          &traj%mass_a, sum_kk, sum_vk)
          deltaE_s(istate)=sum_vk**2/(4.d0*sum_kk)+(traj%Etot-traj%Ekin-&
          &real(traj%H_diag_ss(istate,istate)))
          if (deltaE_s(istate)<0.d0) then
            traj%allowed_hop_s(istate)=0
          else
            traj%allowed_hop_s(istate)=1
          endif

        case (6)    ! correct along projected gradient difference
          call available_ekin(ctrl%natom,&
          &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag, istate,:,:),&
          &traj%mass_a, sum_kk, sum_vk)
          deltaE_s(istate)=sum_vk**2/(4.d0*sum_kk)+(traj%Etot-traj%Ekin-&
          &real(traj%H_diag_ss(istate,istate)))
          if (deltaE_s(istate)<0.d0) then
            traj%allowed_hop_s(istate)=0
          else
            traj%allowed_hop_s(istate)=1
          endif

        case (7)    ! correct along effective NAC
          call available_ekin(ctrl%natom,&
          &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag, istate,:,:),&
          &traj%mass_a, sum_kk, sum_vk)
          deltaE_s(istate)=sum_vk**2/(4.d0*sum_kk)+(traj%Etot-traj%Ekin-&
          &real(traj%H_diag_ss(istate,istate)))
          if (deltaE_s(istate)<0.d0) then
            traj%allowed_hop_s(istate)=0
          else
            traj%allowed_hop_s(istate)=1
          endif

        case (8)    ! correct along projected effective NAC
          call available_ekin(ctrl%natom,&
          &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag, istate,:,:),&
          &traj%mass_a, sum_kk, sum_vk)
          deltaE_s(istate)=sum_vk**2/(4.d0*sum_kk)+(traj%Etot-traj%Ekin-&
          &real(traj%H_diag_ss(istate,istate)))
          if (deltaE_s(istate)<0.d0) then
            traj%allowed_hop_s(istate)=0
          else
            traj%allowed_hop_s(istate)=1
          endif

      endselect
    
      if (deltaE_s(istate) .lt. 0.d0) then
        traj%uncertainty_time_s(istate)=0.5d0/abs(deltaE_s(istate))
      else
        traj%uncertainty_time_s(istate)=0.d0
      endif

      else ! if (istate.ne.traj%state_diag) then
        traj%allowed_hop_s(istate)=1 
      endif ! if (istate.ne.traj%state_diag) then

    enddo stateloop

    if (printlevel>2) then
      do istate=1,ctrl%nstates
        write(u_log,*) "state:", istate, "; allowed(1)/frustrated(0):", traj%allowed_hop_s(istate)
      enddo
      call vecwrite(ctrl%nstates, deltaE_s, u_log, 'available energy', 'F24.12')
      call vecwrite(ctrl%nstates, traj%uncertainty_time_s, u_log, 'uncertainty time', 'F24.12')
    endif

  endsubroutine

! ===========================================================
! compute available energy along hopping direction
! notice this is different compared to subroutine available_ekin
! in terms of units
  subroutine available_e(natom,veloc_ad,nac_ad,mass_a, sum_kk, sum_vk)
    implicit none
    integer, intent(in) :: natom
    real*8, intent(in) :: veloc_ad(natom,3), nac_ad(natom,3), mass_a(natom)
    real*8, intent(out) :: sum_kk, sum_vk

    integer :: iatom, idir

    sum_kk=0.d0
    sum_vk=0.d0
    do iatom=1,natom
      do idir=1,3
        sum_kk=sum_kk+2.0d0*sum( nac_ad(:,idir)*nac_ad(:,idir) )
        sum_vk=sum_vk+      sum( nac_ad(:,idir)*veloc_ad(:,idir)*sqrt(mass_a(:)) )
      enddo
    enddo

  endsubroutine

! ===========================================================
! wrap up subroutine for writing restart file for tshtu if hopping is allowed
  subroutine record_time_travelling_point(traj,ctrl)
    use definitions
    use matrix
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl

    integer :: istate 

    do istate=1,ctrl%nstates
      ! write restart file
      call write_tshtu_restart_traj(u_rest,ctrl,traj,istate)
    enddo

  endsubroutine

! ===========================================================
! wrap up subroutine for writing restart file for tshtu if hopping is allowed
  subroutine save_time_travelling_point(traj,ctrl)
    use definitions
    use matrix
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl

    integer :: istate


    do istate=1,ctrl%nstates
 
      if (traj%allowed_hop_s(istate)==1) then ! this state is hopping allowed 
        ! record the time: notice the records are for last time step. 
        traj%allowed_time_s(istate)=max(traj%microtime,0.d0)
        ! write restart file
        call copy_tshtu_restart_traj(ctrl,traj,istate)
      endif

    enddo

    if (printlevel>2) then
      call vecwrite(ctrl%nstates, traj%allowed_time_s*au2fs, u_log, 'current record point', 'F14.9')
    endif

  endsubroutine

! ===========================================================
! copy the restart file
  subroutine copy_tshtu_restart_traj(ctrl,traj,l)
    use definitions
    use matrix
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl
    integer, intent (in) :: l

    character*10 :: state_id
    character*100 :: file_name1, file_name2

    write(state_id, '(i0)') l
    file_name1='restart_save_tshtu_state'//trim(adjustl(state_id))//'.traj'
    file_name2='restart_tshtu_state'//trim(adjustl(state_id))//'.traj'
    call system("cp "//trim(file_name1)//" "//trim(file_name2))

  endsubroutine


! ===========================================================
! write the trajectory restart file for each state 
  subroutine write_tshtu_restart_traj(u,ctrl,traj,l)
    use definitions
    use matrix
    implicit none
    integer :: u
    integer, intent (in) :: l
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl

    character*10 :: state_id
    character*100 :: file_name
    integer :: iatom, i,j,k
    character(8000) :: string

    write(state_id, '(i0)') l
    file_name='restart_save_tshtu_state'//trim(adjustl(state_id))//'.traj'

    open(unit=u,file=trim(file_name),status='replace',action='write')
 
    ! notice the following information is not recorded
    ! because when time traveling happens, the following information should not
    ! be reverted.
    ! traj%state_diag_frust
    ! traj%uncertainty_time_frust
    ! traj%incident_time
    ! traj%in_time_uncertainty
    ! traj%tu_backpropagation
    ! traj%travelling_state
    ! traj%target_state


    ! write everything
    write(u,*) traj%RNGseed
    write(u,*) traj%traj_hash
    write(u,*) traj%state_MCH
    write(u,*) traj%state_MCH_old
    write(u,*) traj%state_diag
    write(u,*) traj%state_diag_old
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
    call matwrite(ctrl%nstates, traj%Kmatrix_diag_ss,    u, 'Kmatrix_diag_ss','ES24.16E3')
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

    if (ctrl%army_ants==1) then
      write(u,*) traj%army_ants_weight
      write(u,*) traj%randnum_branching
      write(u,*) traj%branching_likehood
    endif

    if (ctrl%pointer_basis==2) then
      write(u,*) traj%entropy
      write(u,*) traj%entropy_old
      write(u,*) traj%linear_entropy
      write(u,*) traj%linear_entropy_old
      write(u,*) traj%entropy_grad
      write(u,*) traj%linear_entropy_grad
      call matwrite(ctrl%nstates, traj%U_pointer_ss,   u, 'U_ss','ES24.16E3')
      call matwrite(ctrl%nstates, traj%U_pointer_old_ss,   u, 'U_old_ss','ES24.16E3')
    endif

    ! Trajectory consistency
    write(u,*) traj%discrepancy

    close(u)

  endsubroutine

! ===========================================================
! if a frustrated hopping happens, do time uncertainty hopping
  subroutine time_uncertainty_hopping(traj,ctrl)
    use definitions
    use matrix
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl

    real*8 :: time_back, time_forward

    ! compute the uncertainty time for the frustrated state 
    ! i.e., the maximum of allowed travelling time for a hop

    time_back=traj%incident_time-traj%allowed_time_s(traj%state_diag_frust)
    time_forward=traj%microtime-traj%incident_time

    if (printlevel>4) then
      write(u_log,*) "active state:", traj%state_diag
      write(u_log,*) "frustrated state:", traj%state_diag_frust
      write(u_log,*) "target state:", traj%target_state
      write(u_log,*) "incident time:", traj%incident_time*au2fs, "fs"
      write(u_log,*) "uncertainty time for frustrated state:", traj%uncertainty_time_frust*au2fs, "fs"
      !write(u_log,*) "time travelling records for active state", traj%allowed_time_s(traj%state_diag)*au2fs, "fs"
      write(u_log,*) "time travelling records for frustrated state:", traj%allowed_time_s(traj%state_diag_frust)*au2fs, "fs"
      write(u_log,*) "current time:", traj%microtime*au2fs, "fs"
      write(u_log,*) "backward travelling time:", time_back*au2fs, "fs"
      write(u_log,*) "forward travelling time:", time_forward*au2fs, "fs"
    endif
   
    if (traj%uncertainty_time_frust.ge.0.d0 .and. time_back.lt.time_forward .and. time_back.le.traj%uncertainty_time_frust) then ! do backward hop
      traj%tu_backpropagation=1 
      traj%travelling_state=traj%state_diag_frust
      if (printlevel>2) write(u_log,*) "Backward time travelling to an energetically available record point"
      if (printlevel>2) write(u_log,*) "travelling state:", traj%travelling_state

    else if (time_forward.le.traj%uncertainty_time_frust) then ! do forward hop 
      if (printlevel>2) write(u_log,*) "In forward time travelling"
      if (traj%allowed_hop_s(traj%state_diag_frust)==1) then ! the frustrated state is now allowed
        if (printlevel>2) write(u_log,*) "Forward time travelling successfully"
        if (printlevel>2) write(u_log,*) traj%state_diag, "->", traj%state_diag_frust
        traj%state_diag_old=traj%state_diag
        traj%state_diag=traj%state_diag_frust
        traj%kind_of_jump=1
        traj%in_time_uncertainty=0        
      endif  

    else if (time_forward.gt.traj%uncertainty_time_frust) then ! do frustrated hop at original position
      ! travel back to original position, and turn off time uncertainty process
      traj%tu_backpropagation=1
      traj%travelling_state=traj%state_diag
      if (printlevel>2) write(u_log,*) "Backward time travelling to original place"
      if (printlevel>2) write(u_log,*) "travelling state:", traj%travelling_state
    endif

  endsubroutine

! ===========================================================
! tsh time uncertainty back travelling
  subroutine tshtu_time_travelling(u_traj,traj,ctrl)
    use definitions
    use matrix
    use misc
    implicit none
    integer :: u_traj
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl

    call read_tshtu_restart_traj(u_traj,traj,ctrl) 

    ! once time travelling is done, set up status
    traj%in_time_uncertainty=2 
    traj%tu_backpropagation=0

    if (printlevel>2) write(u_log,*) "backward travelling finished" 

  endsubroutine

! ===========================================================
! tsh time uncertainty back travelling
  subroutine read_tshtu_restart_traj(u_traj,traj,ctrl)
    use definitions
    use matrix
    use misc
    use decoherence_afssh
    implicit none
    integer :: u_traj
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl
 
    integer :: imult, iatom, i,j,k, istate,ilaser
    character(8000) :: string
    real*8 :: dummy_randnum
    integer :: time

    character*10 :: state_id
    character*50 :: file_name

    ! notice we are travelling back to a state either it is state_diag or
    ! state_diag_frust
    write(state_id, '(i0)') traj%travelling_state
    file_name='restart_tshtu_state'//trim(adjustl(state_id))//'.traj'

    open(unit=u_traj,file=trim(file_name),status='old',action='read')

    ! read everything
    read(u_traj,*) traj%RNGseed
    read(u_traj,*) traj%traj_hash
    read(u_traj,*) traj%state_MCH
    read(u_traj,*) traj%state_MCH_old
    read(u_traj,*) traj%state_diag
    read(u_traj,*) traj%state_diag_old
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
      call matread(ctrl%nstates, traj%U_pointer_ss,   u_traj,   string)
      call matread(ctrl%nstates, traj%U_pointer_old_ss,   u_traj,   string)
    endif

    ! Trajectory consistency
    read(u_traj,*) traj%discrepancy

    close(u_traj)
    ! Now everything about trajectory has been read from the restart file

  endsubroutine



endmodule
