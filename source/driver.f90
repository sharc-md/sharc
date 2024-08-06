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

!> # Module driver
!> 
!> \author Yinan Shu
!> \date 1.1.2020
!>
!> This module contains the integrators that control the flow of dynamics propagation
!> includes:
!> Bulirsch-Stoer
!> adaptive_velocity_verlet
!> fixed_velocity_verlet

module driver
  contains

!====================================================================
!The main driver for adaptive velocity verlet integrator
subroutine adaptive_velocity_verlet(traj,ctrl)
  use zpe
  use pointer_basis
  use decoherence_afssh
  use decoherence_dom
  use definitions
  use electronic
  use electronic_laser
  use tsh_tu
  use army_ants
  use input
  use matrix
  use misc
  use nuclear
  use qm
  use restart
  use output

  implicit none

  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

do while (traj%microtime.lt.ctrl%tmax)
!do i_step=traj%step+1,ctrl%nsteps
  traj%step=traj%step+1
  ctrl%nsteps=traj%step
  if ( (traj%microtime+ctrl%dtstep) .le. ctrl%tmax) then
    traj%microtime=traj%microtime+ctrl%dtstep
  else
    ctrl%dtstep=ctrl%tmax-traj%microtime
    traj%microtime=ctrl%tmax
  endif

  call write_logtimestep(u_log,traj%step,traj%microtime)

  ! Velocity Verlet x
  call VelocityVerlet_xstep(traj,ctrl)
  ! QM Calculation
  call do_qm_calculations(traj,ctrl)
  ! QM Processing
  call QM_processing(traj,ctrl)
  ! Adjust Phases
  call Adjust_phases(traj,ctrl)
  ! Compute NAC in diagonal basis
  call NAC_processing(traj, ctrl)

  ! Optimize pointer basis
  if (ctrl%pointer_basis==2) then 
    call pointer_basis_opt(traj, ctrl)
  endif

  if (ctrl%method==0) then !TSH
    ! Mix Gradients
    call Mix_gradients(traj,ctrl)
    ! Velocity Verlet v    (before SH)
    call VelocityVerlet_vstep(traj,ctrl)
    if (ctrl%dampeddyn/=1.d0) call Damp_Velocities(traj,ctrl)
    traj%Ekin=Calculate_ekin(ctrl%natom, traj%veloc_ad, traj%mass_a)
    !   call Calculate_etot(traj,ctrl)
    ! Propagation
    if (ctrl%laser==0) then
      call propagate(traj,ctrl)
    else
      call propagate_laser(traj,ctrl)
    endif
    if (ctrl%time_uncertainty==0) then ! SH
      if (ctrl%army_ants==0) then 
        call surface_hopping(traj,ctrl)
      else if (ctrl%army_ants==1) then 
        call army_ants_surface_hopping(traj,ctrl)
      endif
    else if (ctrl%time_uncertainty==1) then ! time uncertainty SH
      call time_uncertainty_surface_hopping(traj,ctrl) 
    endif
    ! Rescale v
    call Rescale_Velocities(traj,ctrl)
    ! Do ZPE correction
    if (ctrl%zpe_correction .ne. 0) then
      call ZPEcorrection(traj,ctrl)
    endif
    call Calculate_etot(traj,ctrl)
    ! Decoherence
    call Decoherence(traj,ctrl)
    ! obtain the correct gradient
    call Calculate_cMCH(traj,ctrl)
    if (ctrl%calc_grad>=1) call redo_qm_gradients(traj,ctrl)
    if (traj%kind_of_jump/=0) call Mix_gradients(traj,ctrl)

  else if (ctrl%method==1) then !SCP
    ! Propagation coherent coefficients
    if (ctrl%laser==0) then
      call propagate(traj,ctrl)
    else
      call propagate_laser(traj,ctrl) !This needs to be done for SCP
    endif
    ! Decoherence, decay of mixing
    call Decoherence(traj,ctrl)
    call Calculate_cMCH(traj,ctrl)
    ! Switching of the pointer state
    if (ctrl%army_ants==0) then 
      call surface_switching(traj,ctrl)
    else if (ctrl%army_ants==1) then
      call army_ants_surface_switching(traj,ctrl)
    endif
    !if (traj%kind_of_jump/=0) then
    !  call redecoherence(traj,ctrl)
    !endif
    ! Mix Gradients
    call Mix_gradients(traj,ctrl)
    ! Velocity Verlet v
    call VelocityVerlet_vstep(traj,ctrl)
    ! Do ZPE correction
    if (ctrl%zpe_correction .ne. 0) then
      call ZPEcorrection(traj,ctrl)
    endif
    if (ctrl%dampeddyn/=1.d0) call Damp_Velocities(traj,ctrl)
    traj%Ekin=Calculate_ekin(ctrl%natom, traj%veloc_ad, traj%mass_a)
    !endif
    call Calculate_etot(traj,ctrl)
  endif

  if (ctrl%dtstep.gt.ctrl%dtstep_min) then
    call Check_Consistency(traj,ctrl) ! check total energy conservation
  else if (ctrl%dtstep.eq.ctrl%dtstep_min) then
    call Check_Consistency(traj,ctrl)
    traj%consistency=1
  else
    traj%consistency=0
    traj%discrepancy=1.d0
  endif

  if (traj%consistency/=2) then ! write output, heading to the next step
    ! Finalization: Variable update, Output, Restart File, Consistency Checks
    call Update_old(traj, ctrl)
    call set_time(traj)
    call write_list_line(u_lis,traj,ctrl)
    call write_dat(u_dat, traj, ctrl)
    call write_geom(u_geo, traj, ctrl)
    ! write_restart_traj must be the last command
    call write_restart_traj(u_rest,ctrl,traj)
    call allflush()
    ! kill trajectory 
    call kill_after_relaxation(traj,ctrl)
    if ((ctrl%killafter>=0).and.(traj%steps_in_gs>ctrl%killafter)) exit
    if (check_stop(ctrl%cwd)) exit
    ctrl%restart=.false.
    if (traj%consistency==1) then
      call Adaptive_Stepsize(traj,ctrl) ! increase the stepsize for next step
    endif
  else ! else, energy conserves badly, back propagate and reduce step size
    traj%step=traj%step-1
    call Adaptive_Stepsize(traj,ctrl) ! decrease the stepsize for next step
    call read_restart_traj(u_resc,u_rest,ctrl,traj) ! back propagate to previous step
  endif

  if (ctrl%time_uncertainty==1 .and. traj%in_time_uncertainty==0) then
    call record_time_travelling_point(traj,ctrl)
  endif

  if (ctrl%time_uncertainty==1 .and. traj%in_time_uncertainty==1 .and. traj%tu_backpropagation==1) then 
    ! back propagate to a point where hopping is energetically allowed, and
    ! forced hop to that state
    call tshtu_time_travelling(u_rest,traj,ctrl)
  endif

enddo

endsubroutine

!====================================================================
!The main driver for fixed stepsize velocity verlet integrator
subroutine fixed_velocity_verlet(traj,ctrl)
  use zpe
  use pointer_basis
  use decoherence_afssh
  use decoherence_dom
  use definitions
  use electronic
  use electronic_laser
  use tsh_tu
  use army_ants
  use input
  use matrix
  use misc
  use nuclear
  use qm
  use restart
  use output

  implicit none

  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

  integer :: itotal_step

  itotal_step=ctrl%tmax/ctrl%dtstep

do while (traj%step.lt.itotal_step)
!do while (traj%microtime.lt.ctrl%tmax)
!do i_step=traj%step+1,ctrl%nsteps
  traj%step=traj%step+1
  ctrl%nsteps=traj%step
  if ( (traj%microtime+ctrl%dtstep) .le. ctrl%tmax) then
    traj%microtime=traj%microtime+ctrl%dtstep
  else
    ctrl%dtstep=ctrl%tmax-traj%microtime
    traj%microtime=ctrl%tmax
  endif

  call write_logtimestep(u_log,traj%step,traj%microtime)

  ! Velocity Verlet x
  call VelocityVerlet_xstep(traj,ctrl)
  ! QM Calculation
  call do_qm_calculations(traj,ctrl)
  ! QM Processing
  call QM_processing(traj,ctrl)
  ! Adjust Phases
  call Adjust_phases(traj,ctrl)
  ! Compute NAC in diagonal basis
  call NAC_processing(traj, ctrl)

  ! Optimize pointer basis
  if (ctrl%pointer_basis==2) then
    call pointer_basis_opt(traj, ctrl)
  endif

  if (ctrl%method==0) then !TSH
    ! Mix Gradients
    call Mix_gradients(traj,ctrl)
    ! Velocity Verlet v    (before SH)
    call VelocityVerlet_vstep(traj,ctrl)
    if (ctrl%dampeddyn/=1.d0) call Damp_Velocities(traj,ctrl)
    traj%Ekin=Calculate_ekin(ctrl%natom, traj%veloc_ad, traj%mass_a)
    !   call Calculate_etot(traj,ctrl)
    ! Propagation
    if (ctrl%laser==0) then
      call propagate(traj,ctrl)
    else
      call propagate_laser(traj,ctrl)
    endif
    if (ctrl%time_uncertainty==0) then ! SH
      if (ctrl%army_ants==0) then
        call surface_hopping(traj,ctrl)
      else if (ctrl%army_ants==1) then
        call army_ants_surface_hopping(traj,ctrl)
      endif
    else if (ctrl%time_uncertainty==1) then ! time uncertainty SH
      call time_uncertainty_surface_hopping(traj,ctrl)
    endif
    ! Rescale v
    call Rescale_Velocities(traj,ctrl)
    ! Do ZPE correction
    if (ctrl%zpe_correction .ne. 0) then
      call ZPEcorrection(traj,ctrl)
    endif
    call Calculate_etot(traj,ctrl)
    ! Decoherence
    call Decoherence(traj,ctrl)
    ! obtain the correct gradient
    call Calculate_cMCH(traj,ctrl)
    if (ctrl%calc_grad>=1) then
      call redo_qm_gradients(traj,ctrl)
      call NAC_processing(traj, ctrl)
    endif
    if (traj%kind_of_jump/=0) then
      call Mix_gradients(traj,ctrl)
    endif

  else if (ctrl%method==1) then !SCP
    ! Propagation coherent coefficients
    if (ctrl%laser==0) then
      call propagate(traj,ctrl)
    else
      call propagate_laser(traj,ctrl) !This needs to be done for SCP
    endif
    ! Decoherence, decay of mixing
    call Decoherence(traj,ctrl)
    call Calculate_cMCH(traj,ctrl)
    ! Switching of the pointer state
    if (ctrl%army_ants==0) then
      call surface_switching(traj,ctrl)
    else if (ctrl%army_ants==1) then
      call army_ants_surface_switching(traj,ctrl)
    endif
    !if (traj%kind_of_jump/=0) then
    !  call redecoherence(traj,ctrl)
    !endif
    ! Mix Gradients
    call Mix_gradients(traj,ctrl)
    ! Velocity Verlet v
    call VelocityVerlet_vstep(traj,ctrl)
    ! Do ZPE correction
    if (ctrl%zpe_correction .ne. 0) then
      call ZPEcorrection(traj,ctrl)
    endif
    if (ctrl%dampeddyn/=1.d0) call Damp_Velocities(traj,ctrl)
    traj%Ekin=Calculate_ekin(ctrl%natom, traj%veloc_ad, traj%mass_a)
    !endif
    call Calculate_etot(traj,ctrl)
  endif

  ! Finalization: Variable update, Output, Restart File, Consistency Checks
  call Update_old(traj, ctrl)
  call set_time(traj)
  call write_list_line(u_lis,traj,ctrl)
  call write_dat(u_dat, traj, ctrl)
  call write_geom(u_geo, traj, ctrl)
  ! write_restart_traj must be the last command
  call write_restart_traj(u_rest,ctrl,traj)
  call allflush()
  ! kill trajectory 
  call kill_after_relaxation(traj,ctrl)
  if ((ctrl%killafter>=0).and.(traj%steps_in_gs>ctrl%killafter)) exit
  if (check_stop(ctrl%cwd)) exit
  ctrl%restart=.false.

  if (ctrl%time_uncertainty==1 .and. traj%in_time_uncertainty==0) then
    call record_time_travelling_point(traj,ctrl)
  endif

  if (ctrl%time_uncertainty==1 .and. traj%in_time_uncertainty==1 .and. traj%tu_backpropagation==1) then
    ! back propagate to a point where hopping is energetically allowed, and
    ! forced hop to that state
    call tshtu_time_travelling(u_rest,traj,ctrl)
  endif

enddo

endsubroutine

!====================================================================
!The main driver for Bulirsch-Stoer integrator
subroutine Bulirsch_Stoer_Hack(traj,ctrl)
  use zpe
  use pointer_basis
  use decoherence_afssh
  use decoherence_dom
  use definitions
  use electronic
  use electronic_laser
  use tsh_tu
  use army_ants
  use input
  use matrix
  use misc
  use nuclear
  use qm
  use bsh
  use restart
  use output

  implicit none

  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

  integer :: bsv

  bsv=6*ctrl%natom+8*ctrl%nstates

do while (traj%microtime.lt.ctrl%tmax)

  ! traj%microtime is updated after BSH_BSH_propagation
  traj%step=traj%step+1
  ctrl%nsteps=traj%step
  if ( (traj%microtime+ctrl%dtstep) .gt. ctrl%tmax) then
    ctrl%dtstep=ctrl%tmax-traj%microtime
    traj%microtime=ctrl%tmax
  endif

  call write_logtimestep(u_log,traj%step,traj%microtime)

  ! in BSH integrator:
  ! electronic coefficients are propagated in MCH basis
  ! nuclear coordinates are propagated on PES in DIAG basis
  ! Later on, we convert MCH coefficients to DIAG coefficients
  call BSH_propagation(traj,ctrl,bsv)

  ! geometry updated, need a new calculation.
  ! QM Calculation
  call do_qm_calculations(traj,ctrl)
  ! QM Processing
  call QM_processing(traj,ctrl)
  ! Adjust Phases
  call Adjust_phases(traj,ctrl)
  ! Compute NAC in diagonal basis
  call NAC_processing(traj, ctrl)

  ! Optimize pointer basis
  if (ctrl%pointer_basis==2) then 
    call pointer_basis_opt(traj, ctrl)
  endif

  ! Have a new U_ss, we can convert coefficients from MCH to DIAG
  call Calculate_cDIAG(traj,ctrl)
 
  if (ctrl%method==0) then !TSH
    if (ctrl%time_uncertainty==0) then ! SH
      if (ctrl%army_ants==0) then 
        call surface_hopping(traj,ctrl)
      else if (ctrl%army_ants==1) then
        call army_ants_surface_hopping(traj,ctrl)
      endif
    else if (ctrl%time_uncertainty==1) then ! time uncertainty SH
      call time_uncertainty_surface_hopping(traj,ctrl)
    endif
    ! Rescale v
    call Rescale_Velocities(traj,ctrl)
  elseif (ctrl%method==1) then !SCP
    ! Switching of the pointer state
    if (ctrl%army_ants==0) then
      call surface_switching(traj,ctrl) 
    else if (ctrl%army_ants==1) then
      call army_ants_surface_switching(traj,ctrl)
    endif
    ! Notice the decay of mixing decoherence part is in Electronic_gradient
  endif
 
  if (ctrl%dampeddyn/=1.d0) call Damp_Velocities(traj,ctrl)
  ! Do ZPE correction
  if (ctrl%zpe_correction .ne. 0) then
    call ZPEcorrection(traj,ctrl)
  endif
  traj%Ekin=Calculate_ekin(ctrl%natom, traj%veloc_ad, traj%mass_a)

  ! Compute electronic and nuclear gradients
  call Electronic_gradients_MCH(traj,ctrl)
  ! Now all the gradients are done, apply decoherence
  ! For TSH, it changes the coefficents
  ! For SCP, decay of mixing, it changes the electronic coefficients gradient
  call Decoherence(traj,ctrl)
  call Mix_gradients(traj,ctrl)

  ! Now done with everything
  if (ctrl%dampeddyn/=1.d0) call Damp_Velocities(traj,ctrl)
  traj%Ekin=Calculate_ekin(ctrl%natom, traj%veloc_ad, traj%mass_a)
  call Calculate_etot(traj,ctrl)

  ! Finalization: Variable update, Output, Restart File
  call Update_old(traj, ctrl)
  call set_time(traj)
  call write_list_line(u_lis,traj,ctrl)
  call write_dat(u_dat, traj, ctrl)
  call write_geom(u_geo, traj, ctrl)
  ! write_restart_traj must be the last command
  call write_restart_traj(u_rest,ctrl,traj)
  call allflush()
  ! kill trajectory 
  call kill_after_relaxation(traj,ctrl)
  if ((ctrl%killafter>=0).and.(traj%steps_in_gs>ctrl%killafter)) exit
  if (check_stop(ctrl%cwd)) exit
  ctrl%restart=.false.

  if (ctrl%time_uncertainty==1 .and. traj%in_time_uncertainty==0) then
    call record_time_travelling_point(traj,ctrl)
  endif

  if (ctrl%time_uncertainty==1 .and. traj%in_time_uncertainty==1 .and. traj%tu_backpropagation==1) then
    ! back propagate to a point where hopping is energetically allowed, and
    ! forced hop to that state
    call tshtu_time_travelling(u_rest,traj,ctrl)
  endif

enddo

endsubroutine

endmodule driver
