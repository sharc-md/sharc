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

!> # Main program
!> \author Sebastian Mai
!> \date 10.07.2014
!> 
!> small modifications using the new defined function in interface.f90
!> \author Maximilian F.S.J. Menger
!> \date 19.04.2018
!>
!> add Ehrenfest and CSDM by using self-consistent potential (SCP)
!> trajectory loop over time instead of steps
!> capability of adapative step size
!> \author Yinan Shu
!> \date 13.11.2019 
!>
!> This is the main code of SHARC
!> It contains the following tasks:
!> - Reading of the input files (in read_input)
!> - Performing all steps of the initialization:
!>   - allocation
!>   - calling the interface for the initial quantum chemistry calculation
!>   - digesting the output of the interface (mixing gradients)
!>   - initial write-out
!> - Main loop of the dynamics:
!>   - Velocity Verlet (nuclear coordinates)
!>   - calling the interface
!>   - digesting the interface output (mixing gradients, phase adjustment)
!>   - electronic wavefunction propagation
!>   - surface hopping
!>   - Calculation of energy (kinetic, potential, total)
!>   - writing of output
!> - Timing
program sharc


#ifndef __PYSHARC__
use zpe
use pointer_basis
use decoherence_afssh
use decoherence_dom
use definitions
use driver
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

!> \param traj Contains all data which would be private to each trajectory in an ensemble
type(trajectory_type) :: traj
!> \param ctrl Contains all data which would be shared in an ensemble
type(ctrl_type) :: ctrl
!> \param i_step Loop variable for the dynamics loop
integer :: i_step
!> \param time Define the integer function time()
integer :: time

! open(0,file='output.err',status='replace',action='write')

traj%time_start=time()
traj%time_last=traj%time_start

call read_input(traj,ctrl)
call allocate_lapack(ctrl%nstates)

if (.not.ctrl%restart) then
  if (ctrl%decoherence==2) call allocate_afssh(traj, ctrl)
  call write_list_header(u_lis)
  call do_initial_qm(traj,ctrl)
  call QM_processing(traj,ctrl)
  if (ctrl%zpe_correction.ne.0) call initial_zpe(traj, ctrl)
  if (ctrl%army_ants==1) call army_ants_initialize(traj, ctrl)
  if (ctrl%pointer_basis==2) call pointer_basis_initialize(traj, ctrl)
  call NAC_processing(traj, ctrl)
  call Calculate_etot(traj,ctrl)
  if (ctrl%time_uncertainty==1) call time_uncertainty_initialize(traj,ctrl)
  if (ctrl%decoherence==11) call initial_def(traj, ctrl)
  call Mix_gradients(traj,ctrl)
  if (ctrl%integrator==0) call Electronic_gradients_MCH(traj,ctrl)
  call Update_old(traj,ctrl)
!  call Calculate_etot(traj,ctrl)
  call set_time(traj)
  call write_dat(u_dat, traj, ctrl)
  call write_list_line(u_lis,traj,ctrl)
  call write_geom(u_geo, traj, ctrl)
  call write_restart_ctrl(u_resc,ctrl)
  call write_restart_traj(u_rest,ctrl,traj)
  call mkdir_restart(ctrl)
endif

if (ctrl%integrator==0) then
  call Bulirsch_Stoer_Hack(traj,ctrl)
elseif (ctrl%integrator==1) then
  call adaptive_velocity_verlet(traj,ctrl)
elseif (ctrl%integrator==2) then
  call fixed_velocity_verlet(traj,ctrl)
endif

call write_final(traj)
    
#else
    use memory_module, only: traj, ctrl
    use qm, only: do_initial_qm, do_qm_calculations, redo_qm_gradients
    implicit none

    !> \param i_step Loop variable for the dynamics loop
    integer :: i_step
    !> \param time Define the integer function time()
    integer :: IRestart, IExit, IRedo
    character*255 :: filename
    integer :: narg
    integer, parameter :: iskip = 1 ! do not skip any steps for writeout
    
    ! get the input filename from the command line argument
    ! get number of command line arguments
    narg=iargc()
    ! input filename must be present as argument
    if (narg==0) then
        write(0,*) 'Usage: sharc <inputfile>'
        stop 1
    endif
    call getarg(1,filename)
    
    ! setup sharc
    call setup_sharc(filename, IRestart)
    ! do initial qm
    if (IRestart .eq. 0) then
      call do_initial_qm(traj,ctrl)
      call initial_step(IRestart)
    endif
    ! everything is set up for the loop
    do i_step=traj%step+1,ctrl%nsteps
      call Verlet_xstep(i_step)
      ! QM Calculation
      call do_qm_calculations(traj,ctrl)
      call Verlet_vstep(IRedo)
      if (IRedo .eq. 1) call redo_qm_gradients(traj,ctrl)
      call Verlet_finalize(IExit, iskip)
      if (IExit .eq. 1) exit
    enddo
    ! finalize, should also deallocate memory for traj, ctrl!
    call finalize_sharc(traj)
#endif

endprogram
