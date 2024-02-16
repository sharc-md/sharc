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

!> # Module QM
!> 
!> \author Sebastian mai
!> \date 27.02.2015
!>                   modified 11.17.2022 by Yinan Shu
!>                       added new QM_processing to post processing quantum chemistry data
!>                   modified 06.22.2022 by Yinan Shu
!>                       NAC_processing subroutine can now deal with new method - kmatrix_method
!>                       to conserve energy in diagonal basis.
!> 
!>                   modified 09.28.2021 by Yinan Shu
!>                       added NAC_processing subroutine, which deals with all NAC related calculations
!> 
!>                   modified 11.13.2019 by Yinan Shu
!>                       changed subroutines select_grad and select_nacdr - for SCP, all states are selected
!>                       changed subroutine Mix_gradients - add the SCP gradient
!> 
!> This module implements the SHARC-QM interface.
!> It writes the QM.in file, calls the interfaces, 
!> and reads the QM.out file to update the electronic matrices.
!> Also performs post-processing of the matrices:
!> - applies frozen-state filters
!> - tracks U matrix phase
!> - corrects state phases
!> - updates _old matrices after a timestep
!> - gradient, non-adiabatic coupling vector selection
!> - calculation of diagonal gradients
!> 
!> In summary, this module updates all relevant quantities in traj at a new timestep,
!> so that the propagation, surface hopping and nuclear dynamics can be performed.

module qm
  contains

  !> Calls the QM calculation for the zero-th timestep and
  !> then performs the remaining initialization:
  !> calculation of initial coefficients and state in diagonal basis
  !> calculation of initial NACdt matrix
  subroutine do_initial_qm(traj,ctrl)
    use definitions
    use electronic
    use matrix
    use output
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl
    integer :: i, iatom, idir

    if (printlevel>1) then
      call write_logtimestep(u_log,traj%step,traj%microtime)
    endif

    ! initial QM calculation
    call do_qm_calculations(traj,ctrl)

    ! correct phases if required
    if (ctrl%track_phase_at_zero==1) then
      call Adjust_phases(traj,ctrl)
    endif

    ! now finalize the initial state, if this was not done in the read_input routine
    if (printlevel>1) then
      write(u_log,*) '============================================================='
      write(u_log,*) '             Initializing states and coefficients'
      write(u_log,*) '============================================================='
    endif

    ! initialize the magnitude of NAC as zero
    traj%dmag=0.d0
    traj%dmag1=0.d0
    traj%dmag2=0.d0

    ! using the U matrix, calculate the remaining coefficient vectors
    select case (ctrl%staterep)
      case (0)      ! coeff is in diag, transform to MCH for printing
        call matvecmultiply(ctrl%nstates,traj%U_ss,traj%coeff_diag_s,traj%coeff_MCH_s,'n')
        traj%state_MCH=state_diag_to_MCH(ctrl%nstates,traj%state_diag,traj%U_ss)
      case (1)      ! coeff is in MCH, transform to diag
        call matvecmultiply(ctrl%nstates,traj%U_ss,traj%coeff_MCH_s,traj%coeff_diag_s,'t')
        traj%state_diag=state_MCH_to_diag(ctrl%nstates,traj%state_MCH,traj%U_ss)
    endselect
 
    ! initialize electronic phase used in BSH integrator
    if (ctrl%integrator==0) then
      traj%ephase_s=0.d0
    endif

    ! check whether the initial state is active
    ! (active states are defined in MCH basis, initial state might be given in diagonal basis
    if (ctrl%actstates_s(traj%state_MCH).eqv..false.) then
      write(0,*) 'Initial state is not active!'
      stop 1
    endif

    if (printlevel>1) then
      write(u_log,'(a,1x,i3,1x,a)') 'Initial state is ',traj%state_mch,'in the MCH basis. '
      write(u_log,'(a,1x,i3,1x,a)') 'Initial state is ',traj%state_diag,'in the DIAG basis. '
      write(u_log,*) 'Coefficients (MCH):'
      write(u_log,'(a3,1x,A12,1X,A12)') '#','Real(c)','Imag(c)'
      do i=1,ctrl%nstates
        write(u_log,'(i3,1x,F12.9,1X,F12.9)') i,traj%coeff_MCH_s(i)
      enddo
      write(u_log,*) 'Coefficients (diag):'
      write(u_log,'(a3,1x,A12,1X,A12)') '#','Real(c)','Imag(c)'
      do i=1,ctrl%nstates
        write(u_log,'(i3,1x,F12.9,1X,F12.9)') i,traj%coeff_diag_s(i)
      enddo
      write(u_log,*)
    endif
    if (abs(traj%coeff_diag_s(traj%state_diag))<1.d-9) then
      write(0,*) 'Initial state has zero population!'
      stop 1
    endif

    ! we have to set up the initial NACdt_ss matrix here -- moved to NAC_processing

    ! initialize preprobability
    traj%preprob_s3=0.d0
    traj%preprob_old_s3=0.d0

  endsubroutine

! ===========================================================

!> This routine performs all steps of the QM calculation, including post-processing 
!> (more post-processing steps) are performed in main.
  subroutine do_qm_calculations(traj,ctrl)
    use definitions
    use electronic
    use matrix
    use qm_out
    use restart
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl
    integer :: stat,i,j,imult,istate,jstate,iatom,idir

    if (printlevel>3) then
      write(u_log,*) '============================================================='
      write(u_log,*) '                       QM calculation'
      write(u_log,*) '============================================================='
      write(u_log,*) 'QMin file="QM/QM.in"'
    endif

    ! open QM.in and write geometry + some keywords (task keywords are written by write_tasks_*)
    open(u_qm_qmin,file='QM/QM.in',status='replace',action='write')
    call write_infos(traj,ctrl)

    ! if necessary, select the quantities for calculation
    if (ctrl%calc_grad==1) call select_grad(traj,ctrl)
    if (ctrl%calc_nacdr==1) call select_nacdr(traj,ctrl)
    if (ctrl%calc_dipolegrad==1) call select_dipolegrad(traj,ctrl)

    ! write tasks for first QM call
    call write_tasks_first(traj,ctrl)
    close(u_qm_qmin)

    ! run QM interface
    if (printlevel>3) write(u_log,*) 'Running file="QM/runQM.sh"'
    call call_runqm(traj)

    ! open QM.out
    if (printlevel>3) write(u_log,*) 'QMout file="QM/QM.out"'
    call open_qmout(u_qm_qmout, 'QM/QM.out')

    ! get Hamiltonian
    call get_hamiltonian(ctrl%nstates, traj%H_MCH_ss)
    ! apply reference energy shift
    do i=1,ctrl%nstates
      traj%H_MCH_ss(i,i)=traj%H_MCH_ss(i,i)-ctrl%ezero
    enddo
    ! apply scaling factor
    if (ctrl%scalingfactor/=1.d0) then
      traj%H_MCH_ss=traj%H_MCH_ss*ctrl%scalingfactor
    endif
    ! apply SOC scaling factor
    if (ctrl%soc_scaling/=1.d0) then
      do istate=1, ctrl%nstates
        do jstate=1, ctrl%nstates
          if (istate.ne.jstate) then
            traj%H_MCH_ss(istate,jstate)=traj%H_MCH_ss(istate,jstate)*ctrl%soc_scaling
          endif
        enddo
      enddo
    endif
    ! apply frozen-state mask
    do i=1,ctrl%nstates
      do j=1,ctrl%nstates
        if (ctrl%actstates_s(i).neqv.ctrl%actstates_s(j)) traj%H_MCH_ss(i,j)=dcmplx(0.d0,0.d0)
        if ((ctrl%calc_soc/=1).and.(i/=j)) traj%H_MCH_ss(i,j)=dcmplx(0.d0,0.d0)
      enddo
    enddo
    if (printlevel>3) write(u_log,'(A31,A2)') 'Hamiltonian:                   ','OK'

    ! get Dipole moments
    if (ctrl%calc_dipole==1) then
      call get_dipoles(ctrl%nstates, traj%DM_ssd)
      if (printlevel>3) write(u_log,'(A31,A2)') 'Dipole Moments:                ','OK'
      traj%DM_print_ssd=traj%DM_ssd
      ! apply frozen-state mask 
      do i=1,ctrl%nstates
        do j=1,ctrl%nstates
          if (ctrl%actstates_s(i).neqv.ctrl%actstates_s(j)) traj%DM_ssd(i,j,:)=dcmplx(0.d0,0.d0)
        enddo
      enddo
    endif

!     ! get Property matrix, if necessary
!     if ((ctrl%ionization>0).and.(mod(traj%step,ctrl%ionization)==0)) then
!       call get_property(ctrl%nstates, traj%Property_ss,stat)
!     else
!       traj%Property_ss=dcmplx(0.d0,0.d0)
!     endif

    ! get all available Properties
    call get_properties_new(ctrl, traj)



    ! if gradients were calculated in first call, get them
    if (ctrl%calc_grad<=1) then
      call get_gradients(ctrl%nstates, ctrl%natom, traj%grad_MCH_sad)
      ! apply scaling factor to gradients
      if (ctrl%scalingfactor/=1.d0) then
        traj%grad_MCH_sad=traj%grad_MCH_sad*ctrl%scalingfactor
      endif
      if (printlevel>3) write(u_log,'(A31,A2)') 'Gradients:                     ','OK'
    endif

    ! if this is not the initial QM calculation
    if (traj%step>=1) then
      ! get NACdt matrix
      if (ctrl%calc_nacdt==1) then
        call get_nonadiabatic_ddt(ctrl%nstates, traj%NACdt_ss)
        ! apply frozen-state mask
        do i=1,ctrl%nstates
          do j=1,ctrl%nstates
            if (ctrl%actstates_s(i).neqv.ctrl%actstates_s(j)) traj%NACdt_ss(i,j)=dcmplx(0.d0,0.d0)
          enddo
        enddo
        if (printlevel>3) write(u_log,'(A31,A2)') 'Non-adiabatic couplings (DDT): ','OK'
      endif

      ! get overlap matrix
      if (ctrl%calc_overlap==1) then
        call get_overlap(ctrl%nstates, traj%overlaps_ss)
        ! apply frozen-state mask
        do i=1,ctrl%nstates
          do j=1,ctrl%nstates
            if (ctrl%actstates_s(i).neqv.ctrl%actstates_s(j)) traj%overlaps_ss(i,j)=dcmplx(0.d0,0.d0)
          enddo
        enddo
        if (printlevel>3) write(u_log,'(A31,A2)') 'Overlap matrix:                ','OK'
      endif

      ! get wavefunction phases
      call get_phases(ctrl%nstates,traj%phases_s,stat)
      if (stat==0) then
        traj%phases_found=.true.
        if (printlevel>3) write(u_log,'(A31,A2)') 'Phases:                        ','OK'
      else
        traj%phases_found=.false.
        if (printlevel>3) write(u_log,'(A31,A9)') 'Phases:                        ','NOT FOUND'
      endif
    endif
    if (traj%step==0) then
      traj%phases_s=dcmplx(1.d0,0.d0)
      if (ctrl%track_phase_at_zero==1) then
        call get_phases(ctrl%nstates,traj%phases_s,stat)
        if (stat==0) then
          traj%phases_found=.true.
          if (printlevel>3) write(u_log,'(A31,A2)') 'Phases:                        ','OK'
        else
          traj%phases_found=.false.
          if (printlevel>3) write(u_log,'(A31,A9)') 'Phases:                        ','NOT FOUND'
        endif
      endif
    endif

    ! get non-adiabatic couplings
    if ( (ctrl%calc_nacdr==0).or.(ctrl%calc_nacdr==1) ) then
      call get_nonadiabatic_ddr(ctrl%nstates, ctrl%natom, traj%NACdr_ssad)
      if (printlevel>3) write(u_log,'(A31,A2)') 'Non-adiabatic couplings (DDR): ','OK'
    endif

    ! get dipole moment derivatives
    if ( (ctrl%calc_dipolegrad==0).or.(ctrl%calc_dipolegrad==1) ) then
      call get_dipolegrad(ctrl%nstates, ctrl%natom, traj%DMgrad_ssdad)
      if (printlevel>3) write(u_log,'(A31,A2)') 'Dipole moment gradients:       ','OK'
    endif
    call close_qmout
    if (printlevel>3) write(u_log,*) ''

    ! ===============================
    ! all quantities read, post-processing
    ! ===============================

    ! here the Hamiltonian is diagonalized, but without phase adjustment
    ! phase adjusted diagonalization is carried out later
    ! here we need to diagonalize only for selection of gradients/couplings/...
    traj%H_diag_ss=traj%H_MCH_ss
    ! if laser field, add it here, without imaginary part
    if (ctrl%laser==2) then
      do i=1,3
        traj%H_diag_ss=traj%H_diag_ss - traj%DM_ssd(:,:,i)*real(ctrl%laserfield_td(traj%step*ctrl%nsubsteps+1,i))
      enddo
    endif

    ! diagonalize, if SHARC dynamics
    if (ctrl%surf==0) then
      call diagonalize(ctrl%nstates,traj%H_diag_ss,traj%U_ss)
    elseif (ctrl%surf==1) then
      traj%U_ss=dcmplx(0.d0,0.d0)
      do i=1,ctrl%nstates
        traj%U_ss(i,i)=dcmplx(1.d0,0.d0)
      enddo
    endif

!     call check_allocation(u_log,ctrl,traj)

    ! get state in all representations
    if ((traj%step==0).and.(ctrl%staterep==1)) then
      traj%state_diag=state_MCH_to_diag(ctrl%nstates,traj%state_MCH,traj%U_ss)
    endif
    traj%state_MCH=state_diag_to_MCH(ctrl%nstates,traj%state_diag,traj%U_ss)

    ! ===============================
    ! now perform a second QM call, where selected gradients/couplings/... are calculated
    ! ===============================

    if (ctrl%calc_second==1) then
      if (printlevel>3) write(u_log,*) 'Doing a second calculation...'
      if (printlevel>3) write(u_log,*) ''

      ! select quantities
      if (ctrl%calc_grad==2) call select_grad(traj,ctrl)
      if (ctrl%calc_nacdr==2) call select_nacdr(traj,ctrl)
      if (ctrl%calc_dipolegrad==2) call select_dipolegrad(traj,ctrl)

      ! write a new QM.in file
      if (printlevel>3) write(u_log,*) ''
      if (printlevel>3) write(u_log,*) 'QMin file="QM/QM.in"'
      open(u_qm_qmin,file='QM/QM.in',status='replace',action='write')
      call write_infos(traj,ctrl)

      ! write tasks
      call write_tasks_second(traj,ctrl)
      close(u_qm_qmin)

      ! call QM interface
      if (printlevel>3) write(u_log,*) 'Running file="QM/runQM.sh"'
      call call_runqm(traj)

      ! open QM.out
      if (printlevel>3) write(u_log,*) 'QMout file="QM/QM.out"'
      call open_qmout(u_qm_qmout, 'QM/QM.out')

      ! read properties: gradients, non-adiabatic couplings, dipole moment derivatives
      if (ctrl%calc_grad==2) then
        call get_gradients(ctrl%nstates, ctrl%natom, traj%grad_MCH_sad)
        ! apply scaling factor to gradient
        if (ctrl%scalingfactor/=1.d0) then
          traj%grad_MCH_sad=traj%grad_MCH_sad*ctrl%scalingfactor
        endif
        if (printlevel>3) write(u_log,'(A31,A2)') 'Gradients:                     ','OK'
      endif
      if (ctrl%calc_nacdr==2) then
        call get_nonadiabatic_ddr(ctrl%nstates, ctrl%natom, traj%NACdr_ssad)
        if (printlevel>3) write(u_log,'(A31,A2)') 'Non-adiabatic couplings (DDR): ','OK'
      endif
      if (ctrl%calc_dipolegrad==2) then
        call get_dipolegrad(ctrl%nstates, ctrl%natom, traj%DMgrad_ssdad)
        if (printlevel>3) write(u_log,'(A31,A2)') 'Dipole moment gradients:       ','OK'
      endif
      call close_qmout
    endif

    if (printlevel>3) write(u_log,*)

    if (printlevel>4) call print_qm(u_log,traj,ctrl)

  endsubroutine

! ===========================================================

!> this routine checks whether after a surface hops there is a need
!> to calculate additional gradients.
!> If necessary, performs the QM call and reads the additional gradients
  subroutine redo_qm_gradients(traj,ctrl)
    use definitions
    use electronic
    use matrix
    use qm_out
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl
    integer :: i

    logical :: old_selg_s(ctrl%nstates)
    real*8 :: old_grad_MCH_sad(ctrl%nstates,ctrl%natom,3)

    if (printlevel>2) then
      write(u_log,*) '============================================================='
      write(u_log,*) '           Checking for additional gradient calculation'
      write(u_log,*) '============================================================='
    endif

    ! old selection mask
    old_selg_s=traj%selg_s
    ! previously calculated gradients
    old_grad_MCH_sad=traj%grad_MCH_sad

    ! make new selection mask
    call select_grad(traj,ctrl)

    if (printlevel>2) then
      write(u_log,*)
      write(u_log,*) 'Previously calculated gradients'
      write(u_log,*) old_selg_s
      write(u_log,*) 'Necessary gradients'
      write(u_log,*) traj%selg_s
    endif
    ! compare old and new selection masks
    traj%selg_s=.not.old_selg_s.and.traj%selg_s
    if (printlevel>2) then
      write(u_log,*) 'Missing gradients'
      write(u_log,*) traj%selg_s
      write(u_log,*)
    endif

    ! if any element of old_selg_s now is not false, the corresponding gradient has to be calculated
    if (any(traj%selg_s)) then
      if (printlevel>3) then
        write(u_log,*) 'Additional calculation necessary'
        write(u_log,*) 
        write(u_log,*) '============================================================='
        write(u_log,*) '                       QM calculation'
        write(u_log,*) '============================================================='
        write(u_log,*) 'Doing a third calculation...'
        write(u_log,*) ''
        write(u_log,*) 'QMin file="QM/QM.in"'
      endif
      open(u_qm_qmin,file='QM/QM.in',status='replace',action='write')
      call write_infos(traj,ctrl)
      call write_tasks_third(traj,ctrl)
      close(u_qm_qmin)
      if (printlevel>3) write(u_log,*) 'Running file="QM/runQM.sh"'
      call call_runqm(traj)
      if (printlevel>3) write(u_log,*) 'QMout file="QM/QM.out"'
      call open_qmout(u_qm_qmout, 'QM/QM.out')
      call get_gradients(ctrl%nstates, ctrl%natom, traj%grad_MCH_sad)
        if (ctrl%scalingfactor/=1.d0) then
          traj%grad_MCH_sad=traj%grad_MCH_sad*ctrl%scalingfactor
        endif
      call close_qmout

      ! insert the previously calculated gradients
      do i=1,ctrl%nstates
        if (.not.traj%selg_s(i)) then
          traj%grad_MCH_sad(i,:,:)=old_grad_MCH_sad(i,:,:)
        endif
      enddo
      if (printlevel>3) write(u_log,'(A31,A2)') 'Gradients:                     ','OK'
      if (printlevel>3) write(u_log,*)
      if (printlevel>4) call print_qm(u_log,traj,ctrl)
    else
      if (printlevel>3) then
        write(u_log,*) 'No further calculation necessary'
      endif
    endif

  endsubroutine

! ===========================================================

!> performs the system call to the QM interface
!> calls QM/runQM.sh
!> In QM/runQM.sh, the correct interface has to be called.
!> sharc.x does not know which interface is employed.
  subroutine call_runqm(traj)
    use definitions
    use output
    implicit none
    type(trajectory_type) :: traj
    integer(KIND=4):: status        ! TODO: check integer/integer(KIND=2)
    character(255) :: command
    integer(KIND=4) :: system       ! TODO: check integer/integer(KIND=2)

    call flush(u_log)
    command='sh QM/runQM.sh'
    status=system(command)          ! TODO: /2**8 factor?

    if (status/=0) then
      write(0,*) 
      write(0,*) '#===================================================#'
      write(0,*) 'QM call was not successful, aborting the run.'
      write(0,*) 'Error code: ',status
      write(0,*) '#===================================================#'
      write(0,*) 
      if (printlevel>0) then
        write(u_log,*) '============================================================='
        write(u_log,*) 'QM call was not successful, aborting the run.'
        write(u_log,*) 'Error code: ',status
        write(u_log,*) '============================================================='
        call write_final(traj)
      endif
      stop 1
    endif

  endsubroutine

! ===========================================================

!> writes ctrl parameters to QM.in, but not task keywords
  subroutine write_infos(traj,ctrl)
    use definitions
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl
    integer :: i,j
    character*1023 :: cwd

    call getcwd(cwd)

    write(u_qm_qmin,'(I6)') ctrl%natom
    write(u_qm_qmin,*) traj%traj_hash
    do i=1,ctrl%natom
      write(u_qm_qmin,'(A2,3(F12.7,1X),3X,3(F12.7,1X))') &
      &traj%element_a(i),(au2a*traj%geom_ad(i,j),j=1,3),(traj%veloc_ad(i,j),j=1,3)
    enddo
    write(u_qm_qmin,'(A)') 'unit angstrom'
    write(u_qm_qmin,'(A)', advance='no') 'states '
    do i=1,ctrl%maxmult
      write(u_qm_qmin,'(I3)', advance='no') ctrl%nstates_m(i)
    enddo
    write(u_qm_qmin,*) 
    write(u_qm_qmin,'(a,1x,F12.6)') 'dt',ctrl%dtstep
    write(u_qm_qmin,'(a,1x,I7)') 'step',traj%step
    write(u_qm_qmin,'(a,1x,a)') 'savedir',trim(cwd)//'/restart'

  endsubroutine

! ===========================================================

!> writes task keywords for all non-selected quantities
!> or for selected quantities which need to be calculated in the first QM call 
!> (select_directly keyword)
  subroutine write_tasks_first(traj,ctrl)
    use definitions
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl
    integer :: i,j

    if ((traj%step==0).and..not.(ctrl%track_phase_at_zero==1)) then
      write(u_qm_qmin,'(A)') 'init'
    endif
    if (ctrl%restart_rerun_last_qm_step) then
      write(u_qm_qmin,'(A)') 'restart'
      ctrl%restart_rerun_last_qm_step=.false.
    endif
    if (ctrl%calc_soc==1) then
      write(u_qm_qmin,'(A)') 'SOC'
    else
      write(u_qm_qmin,'(A)') 'H'
    endif
    if (ctrl%calc_dipole==1) then
      write(u_qm_qmin,'(A)') 'DM'
    endif

    select case (ctrl%calc_grad)
      case (0)
        write(u_qm_qmin,'(A)') 'GRAD all'
      case (1)
        write(u_qm_qmin,'(A)',advance='no') 'GRAD'
        do i=1,ctrl%nstates
          if (traj%selg_s(i)) write(u_qm_qmin,'(1X,I3)',advance='no') i
        enddo
        write(u_qm_qmin,'(1X)')
      case (2)
        write(u_qm_qmin,*)
    endselect

    if ((traj%step==0).and.(ctrl%track_phase_at_zero==1)) then
      write(u_qm_qmin,'(A)') 'PHASES'
    endif
    if (traj%step>=1) then
      if (ctrl%calc_nacdt==1) write(u_qm_qmin,'(A)') 'NACDT'
      if (ctrl%calc_overlap==1) write(u_qm_qmin,'(A)') 'OVERLAP'
      if (ctrl%calc_phases==1) write(u_qm_qmin,'(A)') 'PHASES'
    endif

    select case (ctrl%calc_nacdr)
      case (-1)
        write(u_qm_qmin,*)
      case (0)
        write(u_qm_qmin,'(A)') 'NACDR'
      case (1)
        write(u_qm_qmin,'(A)') 'NACDR SELECT'
        do i=1,ctrl%nstates
          do j=1,ctrl%nstates
            if (traj%selt_ss(i,j)) write(u_qm_qmin,'(I3,1X,I3)') i,j
          enddo
        enddo
        write(u_qm_qmin,'(A)') 'END'
      case (2)
        write(u_qm_qmin,*)
    endselect

    select case (ctrl%calc_dipolegrad)
      case (-1)
        write(u_qm_qmin,*)
      case (0)
        write(u_qm_qmin,'(A)') 'DMDR'
      case (1)
        write(u_qm_qmin,'(A)') 'DMDR SELECT'
        do i=1,ctrl%nstates
          do j=1,ctrl%nstates
            if (traj%seldm_ss(i,j)) write(u_qm_qmin,'(I3,1X,I3)') i,j
          enddo
        enddo
        write(u_qm_qmin,'(A)') 'END'
      case (2)
        write(u_qm_qmin,*)
    endselect

    if (ctrl%ionization>0) then
      if (mod(traj%step,ctrl%ionization)==0) then
        write(u_qm_qmin,'(A)') 'ION'
      endif
    endif

    if (ctrl%theodore>0) then
      if (mod(traj%step,ctrl%theodore)==0) then
        write(u_qm_qmin,'(A)') 'THEODORE'
      endif
    endif

  endsubroutine

! ===========================================================

!> writes task keywords for all selected quantities (grad, nacdr, dmdr)
  subroutine write_tasks_second(traj,ctrl)
    use definitions
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl
    integer :: i,j

    write(u_qm_qmin,'(A)') 'samestep'

    if (ctrl%calc_grad==2) then
      write(u_qm_qmin,'(A)',advance='no') 'GRAD'
      do i=1,ctrl%nstates
        if (traj%selg_s(i)) write(u_qm_qmin,'(1X,I3)',advance='no') i
      enddo
      write(u_qm_qmin,'(1X)')
    endif

    if (ctrl%calc_nacdr==2) then
      write(u_qm_qmin,'(A)') 'NACDR SELECT'
      do i=1,ctrl%nstates
        do j=1,ctrl%nstates
          if (traj%selt_ss(i,j)) write(u_qm_qmin,'(I3,1X,I3)') i,j
        enddo
      enddo
      write(u_qm_qmin,'(A)') 'END'
    endif

    if (ctrl%calc_dipolegrad==2) then
      write(u_qm_qmin,'(A)') 'DMDR SELECT'
      do i=1,ctrl%nstates
        do j=1,ctrl%nstates
          if (traj%seldm_ss(i,j)) write(u_qm_qmin,'(I3,1X,I3)') i,j
        enddo
      enddo
      write(u_qm_qmin,'(A)') 'END'
    endif

  endsubroutine

! ===========================================================

!> writes task keywords for additional gradient calculation
  subroutine write_tasks_third(traj,ctrl)
    use definitions
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl
    integer :: i

    write(u_qm_qmin,'(A)') 'samestep'

!     if (ctrl%calc_grad==2) then
      write(u_qm_qmin,'(A)',advance='no') 'GRAD'
      do i=1,ctrl%nstates
        if (traj%selg_s(i)) write(u_qm_qmin,'(1X,I3)',advance='no') i
      enddo
      write(u_qm_qmin,'(1X)')
!     endif

  endsubroutine

! ===========================================================

!> selects gradients to be calculated by the QM interface.
!> The selection is based on an energy criterion: Only gradients of states
!> closer than <eselect_grad> to the active state are calculated.
  subroutine select_grad(traj,ctrl)
    use definitions
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl
    integer :: i
    real*8 :: E=0.d0

! depend on either use TSH or SCP (Ehrenfest)
    if (ctrl%method==0) then !TSH
      ! in first timestep always calculate all gradients
      if (traj%step==0) then
        traj%selg_s=.true.
        if (ctrl%surf==1) then    ! except if we have non-diagonal dynamics
          traj%selg_s=.false.
          traj%selg_s(traj%state_MCH)=.true.
        endif
      else
        traj%selg_s=.false.
        ! energy-based selection for SHARC
        select case (ctrl%surf)
          case (0) ! SHARC
            E=real(traj%H_diag_ss(traj%state_diag,traj%state_diag))
            do i=1,ctrl%nstates
              if ( abs( traj%H_MCH_ss(i,i) - E )< ctrl%eselect_grad ) traj%selg_s(i)=.true.
            enddo
            traj%selg_s(traj%state_MCH)=.true.

          ! only active state for non-diagonal dynamics
          case (1) ! non-SHARC
            traj%selg_s(traj%state_MCH)=.true. 
        endselect
      endif
      if (ctrl%coupling==3 .and. ctrl%ktdc_method==0) then ! kTDC computed by gradients
        traj%selg_s(:)=.true.
      endif
      if (ctrl%gradcorrect==2 .and. ctrl%kmatrix_method==0) then ! Kmatrix computed by gradients
        traj%selg_s(:)=.true.
      endif
      ! never calculate gradients of frozen states
      traj%selg_s=traj%selg_s.and.ctrl%actstates_s
    else if (ctrl%method==1) then !SCP
      traj%selg_s=.true.
      ! never calculate gradients of frozen states
      traj%selg_s=traj%selg_s.and.ctrl%actstates_s
    endif

    if (printlevel>3) then
      write(u_log,*) '-------------------- Gradient selection ---------------------'
      if (ctrl%method==0) then !TSH
        if (traj%step==0) then
          write(u_log,*) 'Selecting all states in first timestep.'
        else
          write(u_log,*) 'Select gradients:'
          write(u_log,*) 'State(diag)= ',traj%state_diag,'State(MCH)=',traj%state_MCH
          write(u_log,*) 'Selected States(MCH)=',(traj%selg_s(i),i=1,ctrl%nstates)
          if (printlevel>4) then
            write(u_log,'(A,1X,F16.9,1X,F16.9)') 'Energy window:',E-ctrl%eselect_grad,E+ctrl%eselect_grad
            do i=1,ctrl%nstates
              write(u_log,'(I4,1X,F16.9)') i,real( traj%H_MCH_ss(i,i))
            enddo
          endif
        endif
      else if (ctrl%method==1) then !SCP
        write(u_log,*) 'Selecting all states for Ehrenfest dynamics'
      endif
    endif

  endsubroutine

! ===========================================================

!> selects non-adiabatic couplings to be calculated by the QM interface.
!> The selection is based on an energy criterion: Only gradients of states
!> closer than <eselect_nac> to the active state are calculated.
!> 
!> The interface has to take care to only calculate non-zero NACs 
!> (e.g., not between singlets and triplets)
  subroutine select_nacdr(traj,ctrl)
    use definitions
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl
    integer :: i,j
    real*8 :: E=0.d0

! depend on either use TSH or SCP (Ehrenfest)
    if (ctrl%method==0) then !TSH
      ! calculate all non-adiabatic couplings in first timestep
      if (traj%step==0) then
        traj%selt_ss=.true.
      else
        traj%selt_ss=.false.
        ! energy-based selection for SHARC and non-SHARC
        if (ctrl%surf==0) then
          E=real(traj%H_diag_ss(traj%state_diag,traj%state_diag))
        elseif (ctrl%surf==1) then
          E=real(traj%H_MCH_ss(traj%state_MCH,traj%state_MCH))
        endif
        do i=1,ctrl%nstates
          do j=1,ctrl%nstates
            if ( ( abs( traj%H_MCH_ss(i,i) - E )< ctrl%eselect_nac ) .and.&
            &( abs( traj%H_MCH_ss(j,j) - E )< ctrl%eselect_nac ) ) traj%selt_ss(i,j)=.true.
          enddo
        enddo
      endif
  
      ! apply frozen-state mask
      do i=1,ctrl%nstates
        do j=1,ctrl%nstates
          traj%selt_ss(i,j)=traj%selt_ss(i,j).and.ctrl%actstates_s(i).and.ctrl%actstates_s(j)
        enddo
      enddo
    else if (ctrl%method==1) then !SCP
      if (traj%step==0) then
        traj%selt_ss=.true.
      else
        select case (ctrl%coupling)
          case (0)
            traj%selt_ss=.true.
          case (1)
            traj%selt_ss=.true.
          case (2)
            traj%selt_ss=.false.
          case (3) ! approximated time derivative coupling.
            traj%selt_ss=.false.
        endselect
        ! if uses NAC in nuclear eom, select all NACs
        select case (ctrl%neom)
          case(0)
            traj%selt_ss=.true.
        endselect
       endif
   endif

    if (printlevel>3) then
      write(u_log,*) '------------- Non-adiabatic coupling selection --------------'
      if (ctrl%method==0) then !TSH
        if (traj%step==0) then
          write(u_log,*) 'Selecting all states in first timestep.'
        else
          write(u_log,*) 'Select nacs:'
          write(u_log,*) 'State(diag)= ',traj%state_diag,'State(MCH)=',traj%state_MCH
          write(u_log,*) 'Selected Coupled States(MCH)='
          do i=1,ctrl%nstates
            write(u_log,*) (traj%selt_ss(i,j),j=1,ctrl%nstates)
          enddo
        endif
      else if (ctrl%method==1) then !SCP
        if (ctrl%neom==0) then
          write(u_log,*) 'Nuclear EOM using NAC, Selecting all states for Ehrenfest dynamics'
        else
          if (ctrl%coupling==1) then
            write(u_log,*) 'Electronic coupling requires NAC, Selecting all states for Ehrenfest dynamics'
          elseif (ctrl%coupling==2) then
            write(u_log,*) 'Electronic coupling requires overlap, no selection'
          elseif (ctrl%coupling==3) then
            write(u_log,*) 'Electronic coupling using approximated time derivative coupling, no selection'
          endif
        endif
      endif
    endif

  endsubroutine

! ===========================================================

!> selects dipole moment derivatives to be calculated by the QM interface.
!> The selection is based on an energy criterion: Only gradients of states
!> closer than <eselect_dmgrad> to the active state are calculated.
!> 
!> The interface has to take care to only calculate non-zero NACs 
!> (e.g., not between singlets and triplets)
!> 
!> TODO: The threshold eselect_dmgrad should be dependent on the laser energy
  subroutine select_dipolegrad(traj,ctrl)
    use definitions
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl
    integer :: i,j
    real*8 :: E=0.d0

! depend on either use TSH or SCP (Ehrenfest)
    if (ctrl%method==0) then !TSH       
      if (traj%step==0) then
        traj%seldm_ss=.true.
      else
        traj%seldm_ss=.false.
        if (ctrl%surf==0) then
          E=real(traj%H_diag_ss(traj%state_diag,traj%state_diag))
        elseif (ctrl%surf==1) then
          E=real(traj%H_MCH_ss(traj%state_MCH,traj%state_MCH))
        endif
        do i=1,ctrl%nstates
          do j=1,ctrl%nstates
            if ( ( abs( traj%H_MCH_ss(i,i) - E )< ctrl%eselect_dmgrad ) .and.&
            &( abs( traj%H_MCH_ss(j,j) - E )< ctrl%eselect_dmgrad ) ) traj%seldm_ss(i,j)=.true.
          enddo
        enddo
      endif

      ! apply frozen-state mask
      do i=1,ctrl%nstates
        do j=1,ctrl%nstates
          traj%seldm_ss(i,j)=traj%seldm_ss(i,j).and.ctrl%actstates_s(i).and.ctrl%actstates_s(j)
        enddo
      enddo
    else if (ctrl%method==1) then !SCP
      if (traj%step==0) then
        traj%seldm_ss=.true.
      else
        traj%seldm_ss=.true.
      endif
    endif

    if (printlevel>3) then
      write(u_log,*) '------------- Dipole moment gradient selection --------------'
      if (ctrl%method==0) then !TSH
        if (traj%step==0) then
          write(u_log,*) 'Selecting all states in first timestep.'
        else
          write(u_log,*) 'Select dipole gradients:'
          write(u_log,*) 'State(diag)= ',traj%state_diag,'State(MCH)=',traj%state_MCH
          write(u_log,*) 'Selected Coupled States(MCH)='
          do i=1,ctrl%nstates
            write(u_log,*) (traj%seldm_ss(i,j),j=1,ctrl%nstates)
          enddo
        endif
      else if (ctrl%method==1) then !SCP
        write(u_log,*) 'Selecting all states for Ehrenfest dynamics'
      endif
    endif

  endsubroutine

! ===========================================================

!> This routine prints all electronic quantities after a QM call
  subroutine print_qm(u,traj,ctrl)
    use definitions
    use matrix
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl
    integer :: u,i,j,k
    character,dimension(3) :: xyz=(/'x','y','z'/)
    character(255) :: string

    write(u,*) '============================================================='
    write(u,*) '                         QM results'
    write(u,*) '============================================================='
    write(u,*)
    call matwrite(ctrl%nstates,traj%H_MCH_ss,u,'Hamiltonian (MCH basis)','F9.4')
    call matwrite(ctrl%nstates,traj%H_diag_ss,u,'Hamiltonian (diag basis)','F9.4')
    call matwrite(ctrl%nstates,traj%U_ss,u,'U matrix','F9.4')
    write(u,*)
    do i=1,3
      call matwrite(ctrl%nstates,traj%DM_ssd(:,:,i),u,'Dipole matrix (MCH basis) '//xyz(i)//' direction','F9.4')
    enddo
    write(u,*)
    do i=1,ctrl%nstates
      write(string,'(A27,I3)') 'Gradient (MCH basis) state ',i
      call vec3write(ctrl%natom,traj%grad_MCH_sad(i,:,:),u,trim(string),'F9.4')
    enddo
    write(u,*)
    if (traj%step>=1) then
      if (ctrl%calc_nacdt.gt.0) then
        call matwrite(ctrl%nstates,traj%NACdt_ss,u,'Non-adiabatic couplings DDT (MCH basis)','F9.4')
      write(u,*)
      endif
    endif
    if (ctrl%calc_nacdr.ge.0) then
      do i=1,ctrl%nstates
        do j=1,ctrl%nstates
          write(string,'(A45,I3,1X,I3)') 'Non-adiabatic coupling DDR (MCH basis) state ',i,j
          call vec3write(ctrl%natom,traj%NACdr_ssad(i,j,:,:),u,trim(string),'F9.4')
        enddo
      enddo
      write(u,*)
    endif
    if (ctrl%calc_dipolegrad.ge.0) then
      do i=1,ctrl%nstates
        do j=1,ctrl%nstates
          do k=1,3
            write(string,'(A45,I3,1X,I3,1X,A,1X,1A)') 'Dipole moment gradient (MCH basis) state ',i,j,'polarization ',xyz(k)
            call vec3write(ctrl%natom,traj%DMgrad_ssdad(i,j,k,:,:),u,trim(string),'F9.4')
          enddo
        enddo
      enddo
      write(u,*)
    endif
    if (traj%step>=1) then
      if (ctrl%calc_overlap.gt.0) then
        call matwrite(ctrl%nstates,traj%overlaps_ss,u,'Overlap matrix (MCH basis)','F12.9')
      write(u,*)
      endif
    endif
    if (traj%phases_found) then
      call vecwrite(ctrl%nstates,traj%phases_s,u,'Wavefunction phases (MCH basis)','F9.4')
    endif

  endsubroutine

! =========================================================== !

!> Updates all *_old matrices when advancing to a new timestep
  subroutine Update_old(traj, ctrl)
    use definitions
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl

    if (printlevel>3) then
      write(u_log,*) '============================================================='
      write(u_log,*) '                   Advancing to next step'
      write(u_log,*) '============================================================='
    endif
    ! initialize old variables from current ones

    traj%dH_MCH_old_ss=traj%dH_MCH_ss
    traj%dH_MCH_ss=traj%H_MCH_ss-traj%H_MCH_old_ss

    traj%H_MCH_old3_ss=traj%H_MCH_old2_ss
    traj%H_MCH_old2_ss=traj%H_MCH_old_ss
    traj%H_MCH_old_ss=traj%H_MCH_ss
   
    traj%H_diag_old3_ss=traj%H_diag_old2_ss
    traj%H_diag_old2_ss=traj%H_diag_old_ss
    traj%H_diag_old_ss=traj%H_diag_ss
 
    traj%DM_old_ssd=traj%DM_ssd
    traj%U_old_ss=traj%U_ss
    traj%NACdt_old_ss=traj%NACdt_ss
    traj%NACdr_old_ssad=traj%NACdr_ssad
    traj%phases_old_s=traj%phases_s
    traj%grad_MCH_old2_sad=traj%grad_MCH_old_sad
    traj%grad_MCH_old_sad=traj%grad_MCH_sad
    traj%Gmatrix_old_ssad=traj%Gmatrix_ssad

    !nuclear information
    traj%veloc_old_ad=traj%veloc_ad
    ctrl%dtstep_old=ctrl%dtstep

  endsubroutine

! ===========================================================
!> This subroutine computes the NAC in diagonal basis
!> Not use anymore, Yinan merged it into Nac_processing
  subroutine Nac_diagonal(traj, ctrl)
    use definitions
    use matrix
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl
    integer :: i, iatom, idir, istate, jstate, ipol
    complex*16 :: U_temp(ctrl%nstates,ctrl%nstates), H_temp(ctrl%nstates,ctrl%nstates)
    complex*16 :: NACdR_diag_ss(ctrl%nstates,ctrl%nstates), pNACdR_diag_ss(ctrl%nstates,ctrl%nstates)
    character(255) :: string

    if (printlevel>3) then
      write(u_log,*) '============================================================='
      write(u_log,*) '             Current Rotational Matrix U'
      write(u_log,*) '============================================================='
    endif

    ! get correct U matrix (or laser fields etc.)
    if (ctrl%surf==0) then
      if (ctrl%laser==0) then
        U_temp=traj%U_ss
      elseif (ctrl%laser==2) then
        H_temp=traj%H_MCH_ss
        do idir=1,3
          H_temp=H_temp - traj%DM_ssd(:,:,idir)*real(ctrl%laserfield_td(traj%step*ctrl%nsubsteps+1,idir))
        enddo
        call diagonalize(ctrl%nstates,H_temp,U_temp)
      endif
    elseif (ctrl%surf==1) then
      U_temp=traj%U_ss
    endif

    if (printlevel>4) call matwrite(ctrl%nstates,U_temp,u_log,'U_ss','F12.9')

    if (printlevel>3) then
      write(u_log,*) '============================================================='
      write(u_log,*) '       Calculating NAC and pNAC in diagonal basis'
      write(u_log,*) '============================================================='
    endif

    ! compute NACdR in diag basis
    traj%NACdR_diag_ssad=dcmplx(0.d0,0.d0)
    do iatom=1,ctrl%natom
      do idir=1,3
        ! initialize NACdR_diag_ss
        NACdR_diag_ss=dcmplx(0.d0,0.d0)
        ! save NACdR to NACdR_diag_ss
        do istate=1,ctrl%nstates
          do jstate=1,ctrl%nstates
            NACdR_diag_ss(istate,jstate)=traj%NACdR_ssad(istate,jstate,iatom,idir)
          enddo
        enddo
        ! if available, add dipole moment derivatives
        if (ctrl%dipolegrad==1) then
          do istate=1,ctrl%nstates
            do jstate=1,ctrl%nstates
              do ipol=1,3
                NACdR_diag_ss(istate,jstate)=NACdR_diag_ss(istate,jstate)-&
                &traj%DMgrad_ssdad(istate,jstate,ipol,iatom,idir)*ctrl%laserfield_td(traj%step*ctrl%nsubsteps+1,ipol)
              enddo
            enddo
          enddo
        endif
        ! transform NACdR_ss to diagonal basis
        call transform(ctrl%nstates,NACdR_diag_ss,U_temp,'utau')
        ! save full NACdR_diag_ss matrix
        do istate=1,ctrl%nstates
          do jstate=1,ctrl%nstates
            if (istate .ne. jstate) then
              traj%NACdR_diag_ssad(istate,jstate,iatom,idir)=NACdR_diag_ss(istate,jstate)
            else
              traj%NACdR_diag_ssad(istate,jstate,iatom,idir)=dcmplx(0.d0,0.d0)
            endif
          enddo
        enddo
      enddo
    enddo
    if (printlevel>4) then
      do istate=1,ctrl%nstates
        do jstate=1,ctrl%nstates
          write(string,'(A46,I3,1X,I3)') 'Non-adiabatic coupling DDR (diag basis) state ',istate,jstate
          call vec3write(ctrl%natom,traj%NACdR_diag_ssad(istate,jstate,:,:),u_log,trim(string),'F12.9')
        enddo
      enddo
    endif

    ! compute pNACdr in diag basis 
    traj%pNACdR_diag_ssad=dcmplx(0.d0,0.d0)
    do iatom=1,ctrl%natom
      do idir=1,3
        ! initialize pNACdR_diag_ss
        pNACdR_diag_ss=dcmplx(0.d0,0.d0)
        ! save pNACdR to pNACdR_diag_ss
        do istate=1,ctrl%nstates
          do jstate=1,ctrl%nstates
            pNACdR_diag_ss(istate,jstate)=traj%pNACdR_MCH_ssad(istate,jstate,iatom,idir)
          enddo
        enddo
        ! if available, add dipole moment derivatives
        if (ctrl%dipolegrad==1) then
          do istate=1,ctrl%nstates
            do jstate=1,ctrl%nstates
              do ipol=1,3
                pNACdR_diag_ss(istate,jstate)=pNACdR_diag_ss(istate,jstate)-&
                &traj%DMgrad_ssdad(istate,jstate,ipol,iatom,idir)*ctrl%laserfield_td(traj%step*ctrl%nsubsteps+1,ipol)
              enddo
            enddo
          enddo
        endif
        ! transform pNACdR_diag_ss to diagonal basis
        call transform(ctrl%nstates,pNACdR_diag_ss,U_temp,'utau')
        ! save full pNACdR_diag_ss matrix
        do istate=1,ctrl%nstates
          do jstate=1,ctrl%nstates
            if (istate .ne. jstate) then
              traj%pNACdR_diag_ssad(istate,jstate,iatom,idir)=pNACdR_diag_ss(istate,jstate)
            else
              traj%pNACdR_diag_ssad(istate,jstate,iatom,idir)=dcmplx(0.d0,0.d0)
            endif
          enddo
        enddo
      enddo
    enddo
    if (printlevel>4) then
      do istate=1,ctrl%nstates
        do jstate=1,ctrl%nstates
          write(string,'(A56,I3,1X,I3)') 'Projected Non-adiabatic coupling DDR (diag basis) state ',istate,jstate
          call vec3write(ctrl%natom,traj%pNACdR_diag_ssad(istate,jstate,:,:),u_log,trim(string),'F12.9')
        enddo
      enddo
    endif



  endsubroutine


! ===========================================================

!> This routine calculates the diagonal gradients by mixing the MCH gradients
!> and non-adiabatic couplings.
!> See IJQC paper.
  subroutine Mix_gradients(traj,ctrl)
    use definitions
    use matrix
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl
    integer :: i, iatom, idir, istate, jstate, lstate, ipol
! The following variables are used by dynamics uses self-consistent potential 
    complex*16 :: den_diag_ss(ctrl%nstates,ctrl%nstates)
    complex*16 :: den_MCH_ss(ctrl%nstates,ctrl%nstates)
    complex*16 :: Gvector_diag_ssad(ctrl%nstates,ctrl%nstates,ctrl%natom,3)
    complex*16 :: GT_diag(ctrl%nstates,ctrl%nstates)
    complex*16 :: Gvector_MCH_ssad(ctrl%nstates,ctrl%nstates,ctrl%natom,3)
    complex*16 :: GT_MCH(ctrl%nstates,ctrl%nstates)
    complex*16 :: NACT(ctrl%nstates,ctrl%nstates), pNACT(ctrl%nstates,ctrl%nstates)
    real*8 :: grad1_MCH(ctrl%natom,3), grad1_diag(ctrl%natom,3)

    real*8 :: vdotg(ctrl%nstates,ctrl%nstates)
    character(255) :: string

    ! Computation of Gmatrix is moved to NAC_processing.

    if (printlevel>3) then
      write(u_log,*) '============================================================='
      write(u_log,*) '                     Nuclear Gradient'
      write(u_log,*) '============================================================='
    endif

    ! pick the classical gradient for the nuclei
    ! =============================
    ! TSH gradient 
    ! =============================
    if (ctrl%method==0) then !TSH
      traj%grad_ad(:,:)=real(traj%Gmatrix_ssad(traj%state_diag,traj%state_diag,:,:))
      if (printlevel>3) then
        write(u_log,*) ''
        write(u_log,*) 'Gradient of diagonal state',traj%state_diag,'picked.'
        write(u_log,*) ''
        ! TODO: print gradient for printlevel>4
      endif
      if (printlevel>4) then
        write(u_log,*) 'current TSH gradient'
        do iatom=1,ctrl%natom
          write(u_log,'(3(F8.5,3X))') traj%grad_ad(iatom,1),traj%grad_ad(iatom,2),traj%grad_ad(iatom,3)
        enddo
      endif

    ! =============================
    ! SCP gradient 
    ! =============================
    else if (ctrl%method==1) then !SCP
      ! initialize traj%grad_ad
      traj%grad_ad(:,:)=0.d0

      ! compute density matrix
      do istate=1,ctrl%nstates
        do jstate=1,ctrl%nstates
          den_diag_ss(istate,jstate)=traj%coeff_diag_s(istate)*conjg(traj%coeff_diag_s(jstate))
          den_MCH_ss(istate,jstate)=traj%coeff_MCH_s(istate)*conjg(traj%coeff_MCH_s(jstate))
        enddo
      enddo

      ! compute time derivative of density matrix
      do istate=1,ctrl%nstates
        do jstate=1,ctrl%nstates
          traj%dendt_MCH_ss(istate,jstate)=dcmplx(0.d0,0.d0)
          traj%dendt_diag_ss(istate,jstate)=dcmplx(0.d0,0.d0)
          do lstate=1,ctrl%nstates
            traj%dendt_MCH_ss(istate,jstate)=traj%dendt_MCH_ss(istate,jstate)+&
              &(-ii*traj%H_MCH_ss(istate,lstate)-traj%NACdt_ss(istate,lstate))*den_MCH_ss(lstate,jstate)+&
              &(ii*conjg(traj%H_MCH_ss(jstate,lstate))-conjg(traj%NACdt_ss(jstate,lstate)))*den_MCH_ss(istate,lstate)
            traj%dendt_diag_ss(istate,jstate)=traj%dendt_diag_ss(istate,jstate)+&
              &(-ii*traj%H_diag_ss(istate,lstate)-traj%NACdt_diag_ss(istate,lstate))*den_diag_ss(lstate,jstate)+&
              &(ii*conjg(traj%H_diag_ss(jstate,lstate))-conjg(traj%NACdt_diag_ss(jstate,lstate)))*den_diag_ss(istate,lstate)
          enddo 
        enddo
      enddo
 
      ! select the direction in generalized coherent nuclear EOM 
      Gvector_diag_ssad(:,:,:,:)=dcmplx(0.d0,0.d0)
      GT_diag(:,:)=dcmplx(0.d0,0.d0)
      Gvector_MCH_ssad(:,:,:,:)=dcmplx(0.d0,0.d0)
      GT_MCH(:,:)=dcmplx(0.d0,0.d0)
      select case (ctrl%neom)
        case (0)  ! use original propagator, based on NAC
          if (ctrl%nac_projection==1) then
            if (printlevel>4) then
              write(u_log,*) 'Nonadiabatic component of coherent SCP gradient from projected NAC'
            endif
            Gvector_diag_ssad(:,:,:,:)=traj%pNACdR_diag_ssad(:,:,:,:)
            Gvector_MCH_ssad(:,:,:,:)=traj%pNACdR_MCH_ssad(:,:,:,:)
          else if (ctrl%nac_projection==0) then
            if (printlevel>4) then
              write(u_log,*) 'Nonadiabatic component of coherent SCP gradient from NAC'
            endif
            Gvector_diag_ssad(:,:,:,:)=traj%NACdR_diag_ssad(:,:,:,:)
            Gvector_MCH_ssad(:,:,:,:)=traj%NACdR_ssad(:,:,:,:)
          endif
        case (1)  ! use effective NAC
          if (ctrl%nac_projection==1) then
            if (printlevel>4) then
              write(u_log,*) 'Nonadiabatic component of coherent SCP gradient from projected effective NAC'
            endif
            Gvector_diag_ssad(:,:,:,:)=traj%pNACGV_diag_ssad(:,:,:,:)
            Gvector_MCH_ssad(:,:,:,:)=traj%pNACGV_MCH_ssad(:,:,:,:)
          else if (ctrl%nac_projection==0) then
            if (printlevel>4) then
              write(u_log,*) 'Nonadiabatic component of coherent SCP gradient from effective NAC'
            endif
            Gvector_diag_ssad(:,:,:,:)=traj%NACGV_diag_ssad(:,:,:,:)
            Gvector_MCH_ssad(:,:,:,:)=traj%NACGV_MCH_ssad(:,:,:,:)
          endif
      endselect

      ! compute Gvector dot velocity
      do istate=1,ctrl%nstates
        do jstate=1,ctrl%nstates
          do iatom=1,ctrl%natom
            do idir=1,3
              GT_diag(istate,jstate)=GT_diag(istate,jstate)+traj%veloc_ad(iatom,idir)*Gvector_diag_ssad(istate,jstate,iatom,idir)
              GT_MCH(istate,jstate)=GT_MCH(istate,jstate)+traj%veloc_ad(iatom,idir)*Gvector_MCH_ssad(istate,jstate,iatom,idir)
            enddo
          enddo
        enddo
      enddo

      if (ctrl%neom_rep==0) then ! compute SCP force in diagonal representation 

        ! compute adiabatic component 
        do istate=1,ctrl%nstates
          traj%grad_ad(:,:)=traj%grad_ad(:,:)+&
          &real(den_diag_ss(istate,istate)*traj%Gmatrix_ssad(istate,istate,:,:))
        enddo   
        if (printlevel>4) then
          write(u_log,*) 'SCP force is computed with diagonal representation'
          write(u_log,*) 'adiabatic component of coherent SCP gradient'
          do iatom=1,ctrl%natom
            write(u_log,'(3(F8.5,3X))') traj%grad_ad(iatom,1),traj%grad_ad(iatom,2),traj%grad_ad(iatom,3)
          enddo
        endif

        ! compute nonadiabatic component 
        do istate=1,ctrl%nstates
          do jstate=1,ctrl%nstates
            if (istate .ne. jstate) then
              if (GT_diag(istate,jstate) .ne. 0.d0) then
                traj%grad_ad(:,:)=traj%grad_ad(:,:)-&
                  &real(den_diag_ss(jstate,istate)*(traj%H_diag_ss(istate,istate)-traj%H_diag_ss(jstate,jstate))*traj%NACdt_diag_ss(istate,jstate)/&
                  &GT_diag(istate,jstate)*Gvector_diag_ssad(istate,jstate,:,:))
              endif
            endif
          enddo
        enddo
        if (printlevel>4) then
          write(u_log,*) 'total (adiabatic and nonadiabatic) coherent SCP gradient'
          do iatom=1,ctrl%natom
            write(u_log,'(3(F8.5,3X))') traj%grad_ad(iatom,1),traj%grad_ad(iatom,2),traj%grad_ad(iatom,3)
          enddo
        endif

      else if (ctrl%neom_rep==1) then ! compute SCP force in MCH representation 

        ! compute adiabatic component 
        do istate=1,ctrl%nstates
          traj%grad_ad(:,:)=traj%grad_ad(:,:)+&
            &real(den_MCH_ss(istate,istate)*traj%grad_MCH_sad(istate,:,:))
        enddo
        if (printlevel>4) then
          write(u_log,*) 'SCP force is computed with MCH representation'
          write(u_log,*) 'adiabatic component of coherent SCP gradient'
          do iatom=1,ctrl%natom
            write(u_log,'(3(F8.5,3X))') traj%grad_ad(iatom,1),traj%grad_ad(iatom,2),traj%grad_ad(iatom,3)
          enddo
        endif

        do istate=1,ctrl%nstates
          do jstate=1,ctrl%nstates
            do lstate=1,ctrl%nstates
              if (GT_MCH(istate,lstate) .ne. 0.d0) then
                traj%grad_ad(:,:)=traj%grad_ad(:,:)-&
                  &real(traj%NACdt_ss(istate,lstate)*den_MCH_ss(lstate,jstate)*traj%H_MCH_ss(jstate,istate)/&
                  &GT_MCH(istate,lstate)*Gvector_MCH_ssad(istate,lstate,:,:))
              endif
              if (GT_MCH(lstate,jstate) .ne. 0.d0) then
                traj%grad_ad(:,:)=traj%grad_ad(:,:)-&
                  &real(traj%NACdt_ss(lstate,jstate)*den_MCH_ss(istate,lstate)*traj%H_MCH_ss(jstate,istate)/&
                  &GT_MCH(lstate,jstate)*Gvector_MCH_ssad(lstate,jstate,:,:))
              endif
            enddo
          enddo
        enddo
        if (printlevel>4) then
          write(u_log,*) 'total (adiabatic and nonadiabatic) coherent SCP gradient'
          do iatom=1,ctrl%natom
            write(u_log,'(3(F8.5,3X))') traj%grad_ad(iatom,1),traj%grad_ad(iatom,2),traj%grad_ad(iatom,3)
          enddo
        endif

      endif ! if (ctrl%neom_rep==0) then

      ! include decoherence force
      if (ctrl%decoherence>0) then
        do iatom=1,ctrl%natom
          do idir=1,3
            traj%grad_ad(iatom,idir)= traj%grad_ad(iatom,idir)-traj%decograd_ad(iatom,idir)
          enddo
        enddo
      endif
      if (printlevel>4) then
        write(u_log,*) 'SCP gradient with decoehrence force'
        do iatom=1,ctrl%natom
          write(u_log,'(3(F8.5,3X))') traj%grad_ad(iatom,1),traj%grad_ad(iatom,2),traj%grad_ad(iatom,3)
        enddo
      endif
      ! print info
      if (printlevel>3) then
        write(u_log,*) ''
        write(u_log,*) 'Self-consistent potential energy gradient done'
        write(u_log,*) ''
      endif

    endif ! The big if for ctrl%method


  endsubroutine

! ===========================================================

  subroutine Adjust_phases(traj,ctrl)
    use definitions
    use matrix
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl
    integer :: istate, jstate, ixyz
    complex*16:: scalarProd(ctrl%nstates,ctrl%nstates)
    complex*16 :: Utemp(ctrl%nstates,ctrl%nstates), Htemp(ctrl%nstates,ctrl%nstates)

    ! if phases were not found in the QM output, try to obtain it
    if (traj%phases_found.eqv..false.) then

      ! from overlap matrix diagonal
      if (ctrl%calc_overlap==1) then

        if (printlevel>4) then 
          write(u_log,*) 'phase correction based on overlaps'
        endif 
        traj%phases_s=traj%phases_old_s
        do istate=1,ctrl%nstates
          if (real(traj%overlaps_ss(istate,istate))<0.d0) then
            traj%phases_s(istate)=traj%phases_s(istate)*(-1.d0)
          endif
        enddo

      ! from scalar products of old and new NAC vectors
      else if (ctrl%calc_overlap==0 .and. ctrl%calc_nacdr>=0) then

        if (printlevel>4) then
          write(u_log,*) 'phase correction based on NACs'
        endif

        scalarProd=dcmplx(0.d0,0.d0)

        do istate=1,ctrl%nstates
          do jstate=1,ctrl%nstates
            scalarProd(istate,jstate)=phase_from_NAC(ctrl%natom, &
            &traj%nacdr_ssad(istate,jstate,:,:),traj%nacdr_old_ssad(istate,jstate,:,:) )
          enddo
        enddo

        ! case of ddt couplings
        ! TODO

        ! phase from SOC
        if (traj%step>1) then
          do istate=1,ctrl%nstates
            do jstate=1,ctrl%nstates
              if (istate==jstate) cycle
              if (scalarProd(istate,jstate)==dcmplx(0.d0,0.d0) ) then

                scalarProd(istate,jstate)=phase_from_SOC(&
                &traj%H_MCH_ss(istate,jstate),traj%H_MCH_old_ss(istate,jstate), traj%dH_MCH_ss(istate,jstate))

              endif
            enddo
          enddo
        endif

        call fill_phase_matrix(ctrl%nstates,scalarProd)
        if (printlevel>4) call matwrite(ctrl%nstates,scalarProd,u_log,'scalarProd matrix','F4.1')
        traj%phases_s=scalarProd(:,1)

      ! from scalar products of old and new TDCs
      else if (ctrl%calc_overlap==0 .and. ctrl%calc_nacdr<0) then

        if (printlevel>4) then
          write(u_log,*) 'phase correction based on TDCs'
        endif

        scalarProd=dcmplx(0.d0,0.d0)

        do istate=1,ctrl%nstates
          do jstate=1,ctrl%nstates
            scalarProd(istate,jstate)=phase_from_TDC(ctrl%natom, &
            &traj%NACdt_ss(istate,jstate),traj%NACdt_old_ss(istate,jstate))
          enddo
        enddo

        ! phase from SOC
        if (traj%step>1) then
          do istate=1,ctrl%nstates
            do jstate=1,ctrl%nstates
              if (istate==jstate) cycle
              if (scalarProd(istate,jstate)==dcmplx(0.d0,0.d0) ) then

                scalarProd(istate,jstate)=phase_from_SOC(&
                &traj%H_MCH_ss(istate,jstate),traj%H_MCH_old_ss(istate,jstate), traj%dH_MCH_ss(istate,jstate))

              endif
            enddo
          enddo
        endif

        call fill_phase_matrix(ctrl%nstates,scalarProd)
        if (printlevel>4) call matwrite(ctrl%nstates,scalarProd,u_log,'scalarProd matrix','F4.1')
        traj%phases_s=scalarProd(:,1)
      
      endif ! if (ctrl%calc_overlap==1) then
    endif

    ! Patch phases for Hamiltonian, DM matrix ,NACs, Overlap
    ! Bra
    do istate=1,ctrl%nstates
      traj%H_MCH_ss(istate,:)=traj%H_MCH_ss(istate,:)*traj%phases_s(istate)
      traj%DM_ssd(istate,:,:)=traj%DM_ssd(istate,:,:)*traj%phases_s(istate)
      traj%DM_print_ssd(istate,:,:)=traj%DM_print_ssd(istate,:,:)*traj%phases_s(istate)
      !this if is taken off because we add QM processing subroutine
      !if (ctrl%calc_nacdt==1) then
        traj%NACdt_ss(istate,:)=traj%NACdt_ss(istate,:)*traj%phases_s(istate)
      !endif
      if (ctrl%calc_nacdr>=0) then      ! calc_nacdr=0 computes all nacs
        traj%NACdr_ssad(istate,:,:,:)=traj%NACdr_ssad(istate,:,:,:)*real(traj%phases_s(istate))
      endif
      !this if is taken off because we add QM processing subroutine
      if (ctrl%calc_overlap==1) then
        traj%overlaps_ss(istate,:)=traj%overlaps_ss(istate,:)*traj%phases_old_s(istate)
      endif
    enddo
    ! Ket
    do istate=1,ctrl%nstates
      traj%H_MCH_ss(:,istate)=traj%H_MCH_ss(:,istate)*traj%phases_s(istate)
      traj%DM_ssd(:,istate,:)=traj%DM_ssd(:,istate,:)*traj%phases_s(istate)
      traj%DM_print_ssd(:,istate,:)=traj%DM_print_ssd(:,istate,:)*traj%phases_s(istate)
      !if (ctrl%calc_nacdt==1) then
        traj%NACdt_ss(:,istate)=traj%NACdt_ss(:,istate)*traj%phases_s(istate)
      !endif
      if (ctrl%calc_nacdr>=0) then
        traj%NACdr_ssad(:,istate,:,:)=traj%NACdr_ssad(:,istate,:,:)*real(traj%phases_s(istate))
      endif
      !if (ctrl%calc_overlap==1) then
        traj%overlaps_ss(:,istate)=traj%overlaps_ss(:,istate)*traj%phases_s(istate)
      !endif
    enddo

    ! electronic structure phase patching finished
    if (traj%step>0) then
      ! U matrix phase patching follows

      ! ================== S matrix Löwdin

      if (ctrl%calc_overlap==1) then
        call intruder(ctrl%nstates, traj%overlaps_ss)
        call lowdin(ctrl%nstates, traj%overlaps_ss)
      endif

      ! ==================

      traj%H_diag_ss=traj%H_MCH_ss
      if (ctrl%laser==2) then
        do ixyz=1,3
          traj%H_diag_ss=traj%H_diag_ss - traj%DM_ssd(:,:,ixyz)*real(ctrl%laserfield_td(traj%step*ctrl%nsubsteps+1,ixyz))
        enddo
      endif
      if (ctrl%surf==0) then
        ! obtain the diagonal Hamiltonian
        if (printlevel>4) then
          write(u_log,*) '============================================================='
          write(u_log,*) '             Adjusting phase of U matrix'
          write(u_log,*) '============================================================='
          call matwrite(ctrl%nstates,traj%H_diag_ss,u_log,'H_MCH + Field','F12.9')
        endif
        call diagonalize(ctrl%nstates,traj%H_diag_ss,traj%U_ss)
        if (printlevel>4) call matwrite(ctrl%nstates,traj%U_old_ss,u_log,'Old U','F12.9')
        if (printlevel>4) call matwrite(ctrl%nstates,traj%U_ss,u_log,'U before adjustment','F12.9')
        ! obtain the U matrix with the correct phase
        if (ctrl%track_phase/=0) then
          if (ctrl%laser==0) then
            if (ctrl%calc_overlap==1) then
              Utemp=traj%U_old_ss
              call matmultiply(ctrl%nstates,traj%overlaps_ss,traj%U_old_ss,Utemp,'tn')
              Htemp=traj%H_MCH_old_ss
              call transform(ctrl%nstates,Htemp,traj%overlaps_ss,'utau')
              if (printlevel>4) call matwrite(ctrl%nstates,Utemp,u_log,'Old U transformed','F12.9')
              if (printlevel>4) call matwrite(ctrl%nstates,Htemp,u_log,'Old H transformed','F12.9')
              if (printlevel>4) call matwrite(ctrl%nstates,traj%H_MCH_ss,u_log,'New H','F12.9')
              call project_recursive(ctrl%nstates, Htemp, traj%H_MCH_ss, Utemp, traj%U_ss,&
              &ctrl%dtstep, printlevel, u_log)
            else
              call project_recursive(ctrl%nstates, traj%H_MCH_old_ss, traj%H_MCH_ss, traj%U_old_ss, traj%U_ss,&
              &ctrl%dtstep, printlevel, u_log)
            endif
          else
            if (printlevel>4) then
              write(u_log,*) 'Phase tracking turned off with laser fields.'
            endif
          endif
        else
          if (printlevel>4) then
            write(u_log,*) 'Phase tracking turned off.'
          endif
        endif
        if (printlevel>4) then
          call matwrite(ctrl%nstates,traj%U_ss,u_log,'U after adjustment','F12.9')
          Htemp=traj%H_MCH_ss
          call transform(ctrl%nstates,Htemp,traj%U_ss,'utau')
          call matwrite(ctrl%nstates,Htemp,u_log,'H_MCH transformed','F12.9')
          call matwrite(ctrl%nstates,traj%H_diag_ss,u_log,'H_diag','F12.9')
        endif
  !       call diagonalize_and_project(ctrl%nstates,traj%H_diag_ss,traj%U_ss,traj%U_old_ss)
      elseif (ctrl%surf==1) then
        traj%U_ss=dcmplx(0.d0,0.d0)
        do istate=1,ctrl%nstates
          traj%U_ss(istate,istate)=dcmplx(1.d0,0.d0)
        enddo
      endif
    endif

  endsubroutine

! ===========================================================

  complex*16 function phase_from_NAC(natom,nac1,nac2) result(phase)
  implicit none
  integer :: natom
  real*8 :: nac1(natom,3),nac2(natom,3)

  integer :: iatom,idir
  real*8 :: prod
  real*8, parameter :: threshold=0.1d-6

  prod=0.d0
  do iatom=1,natom
    do idir=1,3
      prod=prod+nac1(iatom,idir)*nac2(iatom,idir)
    enddo
  enddo
  if (abs(prod)<threshold) then
    phase=dcmplx(0.d0,0.d0)
    return
  endif
  if (prod<0.d0) then
    phase=dcmplx(-1.d0,0.d0)
  else
    phase=dcmplx(1.d0,0.d0)
  endif

  endfunction

! ===========================================================

  complex*16 function phase_from_TDC(natom,tdc1,tdc2) result(phase)
  implicit none
  integer :: natom
  complex*16 :: tdc1,tdc2

  integer :: iatom,idir
  real*8 :: prod
  real*8, parameter :: threshold=0.1d-6

  prod=real(tdc1)*real(tdc2)

  if (abs(prod)<threshold) then
    phase=dcmplx(0.d0,0.d0)
    return
  endif
  if (prod<0.d0) then
    phase=dcmplx(-1.d0,0.d0)
  else
    phase=dcmplx(1.d0,0.d0)
  endif

  endfunction

! ===========================================================

  complex*16 function phase_from_SOC(soc1,soc2,dsoc) result(phase)
  implicit none
  complex*16 :: soc1,soc2,dsoc       ! old soc, new soc

  complex*16 :: prod, diff
  real*8,parameter :: threshold=(1.d0/219474.d0)**(2)

  prod=conjg(soc1)*soc2

  if (abs(prod)<threshold) then
    phase=dcmplx(0.d0,0.d0)
    return
  endif
  if (real(prod)<0.d0) then
    diff=soc2-soc1
    if (abs(diff)>2.d0*abs(dsoc)) then
      phase=dcmplx(-1.d0,0.d0)
    else
      phase=dcmplx(1.d0,0.d0)
    endif
  else
    phase=dcmplx(1.d0,0.d0)
  endif


  endfunction

! ===========================================================

  subroutine fill_phase_matrix(n,A)
  ! fills up an incomplete phase matrix, e.g.
!  1  0  0  0 
!  0  1 -1  1
!  0 -1  0  0
!  1  1  0  0
  ! is completed to
!  1  1 -1  1
!  1  1 -1  1
! -1 -1 -1 -1
!  1  1 -1  1
  ! actually, it is sufficient to complete the first row, since it contains the necessary phase information
  implicit none
  integer :: n 
  complex*16 :: A(n,n)

  integer :: i, j, k

  do i=1,n
    if (A(i,1)==dcmplx(0.d0,0.d0)) then
      do j=2,n
        if (j==i) cycle
        do k=1,n
          if (k==i) cycle
          if (k==j) cycle
          if ( (A(i,j)/=dcmplx(0.d0,0.d0)).and.(A(k,1)/=dcmplx(0.d0,0.d0)).and.(A(k,j)/=dcmplx(0.d0,0.d0)) ) then
            A(i,1)=A(i,j)*A(k,1)*A(k,j)
          endif
        enddo
      enddo
    endif
  enddo

  do i=1,n
    if (A(i,1)==dcmplx(0.d0,0.d0)) A(i,1)=dcmplx(1.d0,0.d0)
  enddo

  endsubroutine

! ===========================================================
!> This subroutine performs all the post-processing of QM calculations
  subroutine QM_processing(traj,ctrl)
    use definitions
    use matrix
    use nuclear
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl
    integer :: stat,i,j,istate,jstate,iatom,idir,ipol
    integer :: imult

    real*8 :: overlap_sum
    ! npi variables
    complex*16 :: w(ctrl%nstates,ctrl%nstates),tw(ctrl%nstates,ctrl%nstates),dw(ctrl%nstates,ctrl%nstates)
    ! npi new variables
    complex*16 :: wv1,wv2,wv3,wv4,wa,wb,wc,wd,we
    complex*16 :: w_lj,w_lk
    ! ktdc variables
    real*8 :: gv_old(ctrl%nstates), gv(ctrl%nstates)
    complex*16 :: eh(ctrl%nstates)
    complex*16 :: de(ctrl%nstates,ctrl%nstates)
    complex*16 :: de_old(ctrl%nstates,ctrl%nstates)
    complex*16 :: de_old2(ctrl%nstates,ctrl%nstates)
    complex*16 :: de_old3(ctrl%nstates,ctrl%nstates)
    real*8 :: fmag
    integer :: sum_state, first_state, final_state
    complex*16:: scalarProd(ctrl%nstates,ctrl%nstates)
    ! effective NAC variables
    complex*16 :: gvec(ctrl%nstates,ctrl%nstates,ctrl%natom,3) ! difference gradient vector
    real*8 :: normg(ctrl%nstates,ctrl%nstates), normv_old, normv
    real*8 :: vdotg(ctrl%nstates,ctrl%nstates) ! v dot difference gradient vector
    real*8 :: fac1(ctrl%nstates,ctrl%nstates),fac2(ctrl%nstates,ctrl%nstates)
    character(255) :: string

    if (printlevel>3) then
      write(u_log,*) '============================================================='
      write(u_log,*) '                        QM processing'
      write(u_log,*) '============================================================='
      write(u_log,*) 'In QM processing, we will:'
      write(u_log,*) '[1].Compute time derivative couplings in MCH basis'
      write(u_log,*) '    Norm-Perserving Interpolation: G. A. Meek, B. G. Levine,'
      write(u_log,*) '    J. Phys. Chem. Lett. 2014, 5, 2351-2356'
      write(u_log,*) '    Curvature TDC: Y. Shu, L. Zhang, X. Chen, S. Sun, Y. Huang,'
      write(u_log,*) '    D. G. Truhlar, J. Chem. Theory Comput. 2022, 18, 1320-1328'
      write(u_log,*) '[2].If overlaps is not computed, approximate overlaps from TDCs'
      write(u_log,*) '[3].Compute effective NAC in MCH basis if needed'
      write(u_log,*) '    Y. Shu, L. Zhang, S. Sun, D. G. Truhlar,'
      write(u_log,*) '    J. Chem. Theory Comput. 2020, 11, 1135-1140'
      write(u_log,*) 'The following table provides the suggested approches for each coupling type' 
      write(u_log,*) 'Couplings types:       NACdr        |     Overlap      |   Curvature'
      write(u_log,*) 'NAC:                 ab initio      |   effective NAC  | effective NAC' 
      write(u_log,*) 'TDC:                 NAC*velocity   |       NPI        |   curvature'
      write(u_log,*) 'Overlap:         from TDC/ab initio |    ab initio     |   from TDC'
      write(u_log,*) 'Note: One can always chose to compute NAC or Overlap ab initially'
      write(u_log,*) 'TDC and overlap computed in this section will be phase corrected'
    endif

    ! ===============================
    ! 0. Save Hamiltonian, and compute approximated velocity 
    ! ===============================

    ! Approximate velocity by forward verlet
    if (traj%step>=1) then
      if (printlevel>4) write(u_log,*) 'Approximate velocity with forward velocity verlet'
      call VelocityVerlet_vstep_approximate(traj%veloc_old_ad, traj%veloc_app_ad, traj%accel_ad, ctrl%dtstep, ctrl%natom)
    else
      traj%veloc_app_ad=traj%veloc_ad
    endif
    if (printlevel>2) then
      call vec3write(ctrl%natom,traj%veloc_app_ad,u_log,'approximated veloc','F12.9')
    endif


    ! ===============================
    ! 1. Computing time deriative coupling
    ! ===============================

    if (printlevel>3) then
      write(u_log,*) '============================================================='
      write(u_log,*) '         [1].Calculating time derivative coupling'
      write(u_log,*) '============================================================='
    endif

    select case (ctrl%coupling)
      case(1) ! coupling using ddr, NACdt is computed by direct product between NAC and velocity
        if (printlevel>3) write(u_log,*) 'Computing time derivative coupling by direct product'
        traj%NACdt_ss=dcmplx(0.d0,0.d0)
        do iatom=1,ctrl%natom
          do idir=1,3
            traj%NACdt_ss=traj%NACdt_ss+traj%NACdr_ssad(:,:,iatom,idir)*traj%veloc_app_ad(iatom,idir)
          enddo
        enddo
      case(2) ! coupling using overlap, NACdt is computed by NPI
        if (printlevel>3) write(u_log,*) 'Computing time derivative coupling by Norm-Perserving Interpolation'
        traj%NACdt_ss=dcmplx(0.d0,0.d0)
        if (traj%step>=1) then
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
          call matmultiply(ctrl%nstates, tw, dw, traj%NACdt_ss, 'nn')
        endif
        !if (printlevel>3)  call matwrite(ctrl%nstates,traj%NACdt_ss,u_log,'NPI Time Derivative Coupling in MCH basis...','F12.9')
        !traj%NACdt_ss=dcmplx(0.d0,0.d0)
        !if (traj%step>=1) then
        !  do istate=1,ctrl%nstates
        !    do jstate=1,ctrl%nstates
        !      wv1=acos(traj%overlaps_ss(jstate,jstate))-asin(traj%overlaps_ss(jstate,istate))
        !      wa=-sin(wv1)/wv1
        !      wv2=acos(traj%overlaps_ss(jstate,jstate))+asin(traj%overlaps_ss(jstate,istate))
        !      wb=sin(wv2)/wv2
        !      wv3=acos(traj%overlaps_ss(istate,istate))-asin(traj%overlaps_ss(istate,jstate))
        !      wc=sin(wv3)/wv3
        !      wv4=acos(traj%overlaps_ss(istate,istate))+asin(traj%overlaps_ss(istate,jstate))
        !      wd=sin(wv4)/wv4
        !      w_lj=sqrt(1-traj%overlaps_ss(jstate,jstate)**2-traj%overlaps_ss(istate,jstate)**2)
        !      w_lk=(-traj%overlaps_ss(jstate,istate)*traj%overlaps_ss(jstate,jstate)-traj%overlaps_ss(istate,istate)*traj%overlaps_ss(istate,jstate))&
        !           &/w_lj
        !      we=2*asin(w_lj)*(w_lj*w_lk*asin(w_lj)+(sqrt((1-w_lj**2)*(1-w_lk**2))-1)*asin(w_lk))&
        !         &/((asin(w_lj))**2-(asin(w_lk))**2)
        !      traj%NACdt_ss(istate,jstate)=(acos(traj%overlaps_ss(jstate,jstate))*(wa+wb)+(asin(traj%overlaps_ss(istate,jstate))*(wc+wd))+we)/2/ctrl%dtstep
        !    enddo
        !  enddo
        !endif
        !if (printlevel>3)  call matwrite(ctrl%nstates,traj%NACdt_ss,u_log,'NPI (at t+0.5dt) Time Derivative Coupling in MCH basis...','F12.9')
      case(3) ! coupling using approximated TDC, NACdt is computed by approximation.
        if (printlevel>3) write(u_log,*) 'Computing time derivative coupling by Curvature Approximation'
        traj%NACdt_ss=dcmplx(0.d0,0.d0)
        if (ctrl%ktdc_method==0) then
          if (printlevel>4) write(u_log,*) 'Curvature TDC is computed by first order difference of gradients'
          ! Compute d(deltaV)/dt = d(deltaV)/dR * velocity for old and current step
          if (traj%step>=1) then
            gv_old(:)=0.d0
            gv(:)=0.d0
            do iatom=1,ctrl%natom
              do idir=1,3
                gv_old(:)=gv_old(:)+traj%grad_MCH_old_sad(:,iatom,idir)*traj%veloc_old_ad(iatom,idir)
                gv(:)=gv(:)+traj%grad_MCH_sad(:,iatom,idir)*traj%veloc_app_ad(iatom,idir)
              enddo
            enddo
            eh(:)=(gv(:)-gv_old(:))/(ctrl%dtstep)
            ! only the states with same spin will be non-zero 
            sum_state=0
            do imult=1, ctrl%maxmult
              first_state=1+sum_state
              sum_state=sum_state+ctrl%nstates_m(imult)
              final_state=sum_state
              do istate=first_state, final_state, 1
                do jstate=istate+1, final_state, 1
                  if (((traj%H_MCH_ss(jstate,jstate)-traj%H_MCH_ss(istate,istate))) .eq. 0.d0) then
                    fmag=(eh(jstate)-eh(istate))/1.d-8
                  else
                    fmag=(eh(jstate)-eh(istate))/(traj%H_MCH_ss(jstate,jstate)-traj%H_MCH_ss(istate,istate))
                  endif
                  if (fmag>0.d0) then
                    traj%NACdt_ss(istate,jstate)=0.5*sqrt(fmag)
                  else
                    traj%NACdt_ss(istate,jstate)=dcmplx(0.d0,0.d0)
                  endif
                enddo
              enddo
            enddo
            do istate=1, ctrl%nstates
              do jstate=1, ctrl%nstates
                if (jstate<istate) then
                  traj%NACdt_ss(istate,jstate)=-traj%NACdt_ss(jstate,istate)
                endif
              enddo
            enddo
          endif !if (traj%step>=1) then
        else if (ctrl%ktdc_method==1) then
          if (printlevel>4) write(u_log,*) 'Curvature TDC is computed by second order difference of energies'
          do istate=1, ctrl%nstates
            do jstate=1, ctrl%nstates
              if (istate.ne.jstate) then 
                de(istate,jstate)=traj%H_MCH_ss(istate,istate)-traj%H_MCH_ss(jstate,jstate)
                de_old(istate,jstate)=traj%H_MCH_old_ss(istate,istate)-traj%H_MCH_old_ss(jstate,jstate)
                de_old2(istate,jstate)=traj%H_MCH_old2_ss(istate,istate)-traj%H_MCH_old2_ss(jstate,jstate)
                de_old3(istate,jstate)=traj%H_MCH_old3_ss(istate,istate)-traj%H_MCH_old3_ss(jstate,jstate)
                if (de(istate,jstate).eq.0.d0) then
                  if (printlevel>4)  write(u_log,*) "WARNING, two MCH states ", istate, jstate, "reach exact degeneracy, set de=1.d-8"
                  de(istate,jstate)=1.d-8
                endif
              else
                de(istate,jstate)=0.d0
                de_old(istate,jstate)=0.d0
                de_old2(istate,jstate)=0.d0
                de_old3(istate,jstate)=0.d0
              endif
            enddo
          enddo
          if (ctrl%integrator==2) then ! fixed time step.
            if (traj%step==2) then
            ! only the states with same spin will be non-zero
              sum_state=0
              do imult=1, ctrl%maxmult
                first_state=1+sum_state
                sum_state=sum_state+ctrl%nstates_m(imult)
                final_state=sum_state
                do istate=first_state, final_state, 1
                  do jstate=istate+1, final_state, 1
                    fmag=(de(istate,jstate)-2*de_old(istate,jstate)+de_old2(istate,jstate))/(ctrl%dtstep**2)/de(istate,jstate)
                    if (fmag>0.d0) then
                      traj%NACdt_ss(istate,jstate)=0.5*sqrt(fmag)
                    else
                      traj%NACdt_ss(istate,jstate)=dcmplx(0.d0,0.d0)
                    endif
                  enddo
                enddo
              enddo
              do istate=1, ctrl%nstates
                do jstate=1, ctrl%nstates
                  if (jstate<istate) then
                    traj%NACdt_ss(istate,jstate)=-traj%NACdt_ss(jstate,istate)
                  endif
                enddo
              enddo
            else if (traj%step>=3) then
            ! only the states with same spin will be non-zero
              sum_state=0
              do imult=1, ctrl%maxmult
                first_state=1+sum_state
                sum_state=sum_state+ctrl%nstates_m(imult)
                final_state=sum_state
                do istate=first_state, final_state, 1
                  do jstate=istate+1, final_state, 1
                    fmag=(2*de(istate,jstate)-5*de_old(istate,jstate)+4*de_old2(istate,jstate)-de_old3(istate,jstate))/(ctrl%dtstep**2)/de(istate,jstate)
                    if (fmag>0.d0) then
                      traj%NACdt_ss(istate,jstate)=0.5*sqrt(fmag)
                    else
                      traj%NACdt_ss(istate,jstate)=dcmplx(0.d0,0.d0)
                    endif
                  enddo
                enddo
              enddo
              do istate=1, ctrl%nstates
                do jstate=1, ctrl%nstates
                  if (jstate<istate) then
                    traj%NACdt_ss(istate,jstate)=-traj%NACdt_ss(jstate,istate)
                  endif
                enddo
              enddo
            endif ! if (traj%step==2) then
          else if (ctrl%integrator==1) then ! adaptive time step, only use first order finite difference
            sum_state=0
            do imult=1, ctrl%maxmult
              first_state=1+sum_state
              sum_state=sum_state+ctrl%nstates_m(imult)
              final_state=sum_state
              do istate=first_state, final_state, 1
                do jstate=istate+1, final_state, 1
                  fmag=(de(istate,jstate)-2*de_old(istate,jstate)+de_old2(istate,jstate))/(ctrl%dtstep**2)/de(istate,jstate)
                  if (fmag>0.d0) then
                    traj%NACdt_ss(istate,jstate)=0.5*sqrt(fmag)
                  else
                    traj%NACdt_ss(istate,jstate)=dcmplx(0.d0,0.d0)
                  endif
                enddo
              enddo
            enddo
            do istate=1, ctrl%nstates
              do jstate=1, ctrl%nstates
                if (jstate<istate) then
                  traj%NACdt_ss(istate,jstate)=-traj%NACdt_ss(jstate,istate)
                endif
              enddo
            enddo
          endif ! if (ctrl%integrator==2) then 
        endif ! if (ctrl%ktdc_method==0) then
    endselect
 
    if (printlevel>3)  call matwrite(ctrl%nstates,traj%NACdt_ss,u_log,'Time Derivative Coupling in MCH basis...','F12.9')

    ! do phase correction on curvature TDC therefore it is consistent with spin-orbit coupling
    if (ctrl%coupling==3) then ! for curvature drivent TDC

      scalarProd=dcmplx(0.d0,0.d0)
      if (traj%step>1) then
        do istate=1,ctrl%nstates
          do jstate=1,ctrl%nstates
            if (istate==jstate) cycle
            if (scalarProd(istate,jstate)==dcmplx(0.d0,0.d0) ) then
              scalarProd(istate,jstate)=phase_from_SOC(&
              &traj%H_MCH_ss(istate,jstate),traj%H_MCH_old_ss(istate,jstate), traj%dH_MCH_ss(istate,jstate))
            endif 
          enddo
        enddo
      endif
      call fill_phase_matrix(ctrl%nstates,scalarProd)
      do istate=1,ctrl%nstates
        traj%NACdt_ss(istate,:)=traj%NACdt_ss(istate,:)*traj%phases_s(istate)
      enddo
      do istate=1,ctrl%nstates
        traj%NACdt_ss(:,istate)=traj%NACdt_ss(:,istate)*traj%phases_s(istate)
      enddo

      if (printlevel>3)  call matwrite(ctrl%nstates,traj%NACdt_ss,u_log,'Time Derivative Coupling in MCH basis after phase consistency...','F12.9')

    endif 


    ! ===============================
    ! 2. Approximate overlaps from TDC
    ! ===============================

    if (ctrl%calc_overlap==0) then

      if (printlevel>3) then
        write(u_log,*) '============================================================='
        write(u_log,*) '           [2].Calculating overlaps from TDC'
        write(u_log,*) '============================================================='
      endif

      ! approximate overlaps from curvature approximated TDC 
      traj%overlaps_ss=dcmplx(0.d0,0.d0)
      do istate=1,ctrl%nstates
        overlap_sum=0.d0
        do jstate=1,ctrl%nstates
          if (istate.ne.jstate) then
            traj%overlaps_ss(istate,jstate)=traj%NACdt_old_ss(istate,jstate)*ctrl%dtstep
            overlap_sum=overlap_sum+real(traj%overlaps_ss(istate,jstate)**2)
          endif
        enddo
        ! rescale in case traj%overlaps_ss(istate,jstate) exceeds unity. 
        if (overlap_sum>1.d0) then
          do jstate=1,ctrl%nstates
            if (istate.ne.jstate) then
              traj%overlaps_ss(istate,jstate)=traj%overlaps_ss(istate,jstate)/overlap_sum
            endif
          enddo
        endif
        traj%overlaps_ss(istate,istate)=sqrt(1.d0-overlap_sum)
      enddo

      if (printlevel>3) then
        call matwrite(ctrl%nstates,traj%overlaps_ss,u_log,'Approximated overlap matrix from TDC (MCH basis)','F12.9')
      endif

    endif ! if (ctrl%calc_overlap==0) then


    ! ===============================
    ! 3. Compute effective NAC in MCH basis
    ! ===============================

    if (ctrl%calc_effectivenac==1) then !one employs effective NAC in nuclear EOM. 

      if (printlevel>3) then
        write(u_log,*) '============================================================='
        write(u_log,*) '        [3].Calculating effective NAC in MCH basis'
        write(u_log,*) '============================================================='
      endif

      ! initialize
      traj%NACGV_MCH_ssad=0.d0

      ! compute norm of v
      normv=0.d0
      do iatom=1,ctrl%natom
        do idir=1,3
          normv=normv+traj%veloc_ad(iatom,idir)*traj%veloc_ad(iatom,idir)
        enddo
      enddo
      normv=sqrt(normv)

      !first in MCH basis...
      if (printlevel>4) write(u_log,*) 'Computing effective NAC in MCH basis...'
      gvec=dcmplx(0.d0,0.d0)
      ! compute guess g vector
      vdotg=0.d0
      normg=0.d0
      do istate=1,ctrl%nstates
        do jstate=1,ctrl%nstates
          gvec(istate,jstate,:,:)=traj%grad_MCH_sad(istate,:,:)-traj%grad_MCH_sad(jstate,:,:)
          do iatom=1,ctrl%natom
            do idir=1,3
              normg(istate,jstate)=normg(istate,jstate)+gvec(istate,jstate,iatom,idir)*gvec(istate,jstate,iatom,idir)
              vdotg(istate,jstate)=vdotg(istate,jstate)+traj%veloc_ad(iatom,idir)*real(gvec(istate,jstate,iatom,idir))
            enddo
          enddo
          normg(istate,jstate)=sqrt(normg(istate,jstate))
        enddo
      enddo
      ! rotate g vector, such that vdotg(istate,jstate)=traj%NACdt_ss(istate,jstate)
      fac1=0.d0
      fac2=1.d0
      do istate=1,ctrl%nstates
        do jstate=istate+1,ctrl%nstates
          if (real(traj%NACdt_ss(istate,jstate)) .ne. 0.000) then
            fac1(istate,jstate)=(real(traj%NACdt_ss(istate,jstate))-vdotg(istate,jstate))/normv/normv
            fac2(istate,jstate)=1.d0
            gvec(istate,jstate,:,:)=fac1(istate,jstate)*traj%veloc_ad(:,:)+fac2(istate,jstate)*gvec(istate,jstate,:,:)
            fac1(jstate,istate)=-fac1(istate,jstate)
            fac2(jstate,istate)=-fac2(istate,jstate)
            gvec(jstate,istate,:,:)=-gvec(istate,jstate,:,:)
          else ! if traj%NACdt_ss is exactly zero, set effective NAC as zero.
            fac1(istate,jstate)=-123.d0
            fac2(istate,jstate)=-123.d0
            fac1(jstate,istate)=fac1(istate,jstate)
            fac2(jstate,istate)=fac2(istate,jstate)
            gvec(istate,jstate,:,:)=0.d0
            gvec(jstate,istate,:,:)=0.d0
            endif
        enddo
      enddo
      if (printlevel>5) then
        call matwrite(ctrl%nstates,fac1,u_log,'Mixing coefficients of velocity vector','F14.9')
        call matwrite(ctrl%nstates,fac2,u_log,'Mixing coefficients of difference gradient vector','F14.9')
      endif
      traj%NACGV_MCH_ssad=gvec

      if (printlevel>4) then
        do i=1,ctrl%nstates
          do j=1,ctrl%nstates
            write(string,'(A55,I3,1X,I3)') 'Effective Non-adiabatic coupling DDR (MCH basis) state ',i,j
            call vec3write(ctrl%natom,traj%NACGV_MCH_ssad(i,j,:,:),u_log,trim(string),'F9.4')
          enddo
        enddo
      endif

    endif ! if (ctrl%calc_effectivenac==1) then


  endsubroutine

! ===========================================================
!> This subroutine performs all the post-processing of NACs
  subroutine NAC_processing(traj,ctrl)
    use definitions
    use matrix
    use electronic, only:state_diag_to_MCH
    use nuclear
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl
    integer :: stat,i,j,istate,jstate,kstate,iatom,idir,ipol
    integer :: imult

    ! 0. compute rotation matrix variables.
    complex*16 :: U_temp(ctrl%nstates,ctrl%nstates), H_temp(ctrl%nstates,ctrl%nstates)
    ! 1. tdh variables
    complex*16 :: K1matrix_MCH_ss(ctrl%nstates,ctrl%nstates)
    complex*16 :: K1matrix_diag_ss(ctrl%nstates,ctrl%nstates)
    real*8 :: gv_old(ctrl%nstates), gv(ctrl%nstates)
    ! 2. ngh, Gmatrices variables 
    complex*16 :: G1matrix_ssad(ctrl%nstates,ctrl%nstates,ctrl%natom,3)
    complex*16 :: G1matrix_ss(ctrl%nstates,ctrl%nstates)
    complex*16 :: G2matrix_ssad(ctrl%nstates,ctrl%nstates,ctrl%natom,3)
    complex*16 :: G2matrix_ss(ctrl%nstates,ctrl%nstates)
    ! 3. projection variables
    real*8 :: NACtmp_MCH(3*ctrl%natom), pNACtmp_MCH(3*ctrl%natom)
    complex*16 :: NACtmp_diag(3*ctrl%natom), pNACtmp_diag(3*ctrl%natom)
    complex*16 :: ctrans_rot_P(3*ctrl%natom,3*ctrl%natom)
    ! 5. Patch gmatrix
    complex*16 :: Gmatrix_ss(ctrl%nstates,ctrl%nstates)
    ! 6. hopping direction and frustared hop velocity reflection vector variables
    real*8 :: hopping_tmp(3*ctrl%natom), phopping_tmp(3*ctrl%natom)

    character(255) :: string

    if (printlevel>3) then
      write(u_log,*) '============================================================='
      write(u_log,*) '                NAC and Gradient processing'
      write(u_log,*) '============================================================='
      write(u_log,*) 'In NAC and Gradient processing, we will:'
      write(u_log,*) '[0].Compute rotational matrix'
      write(u_log,*) '[1].Compute time derivative Hamiltonian matrix and TDC in diagonal' 
      write(u_log,*) '    basis:'
      write(u_log,*) '    Y. Shu, L. Zhang, D. Wu, X. Chen, S. Sun, D. G. Truhlar'
      write(u_log,*) '    Submitted to J. Chem. Theory Comput. 2022'
      write(u_log,*) '[2].Compute nulcear gradient Hamiltonian matrix and NAC and effective'
      write(u_log,*) '    NAC in diagonal basis'
      write(u_log,*) '[3].Compute projected NAC and effective NAC in MCH and diagonal basis'
      write(u_log,*) '    Y. Shu, L. Zhang, Z. Varga, K. A. Parker, S. Kanchanakungwankul,'
      write(u_log,*) '    S. Sun, D. G. Truhlar, J. Chem. Theory Comput. 2020, 11, 1135-1140'
      write(u_log,*) '[4].Patch Kmatrix - used for TDH gradient correction scheme'
      write(u_log,*) '    Y. Shu, L. Zhang, D. Wu, X. Chen, S. Sun, D. G. Truhlar'
      write(u_log,*) '    Submitted to J. Chem. Theory Comput. 2022' 
      write(u_log,*) '[5].Patch Gmatrix'
      write(u_log,*) '[6].Compute hopping direction and frustrated hop velocity reflection vector'
    endif


    ! ===============================
    ! 0. Compute Rotational Matrix U
    ! ===============================

    if (printlevel>3) then
      write(u_log,*) '============================================================='
      write(u_log,*) '            [0].Current Rotational Matrix U'
      write(u_log,*) '============================================================='
    endif

    ! get correct U matrix (or laser fields etc.)
    if (ctrl%surf==0) then
      if (ctrl%laser==0) then
        U_temp=traj%U_ss
      elseif (ctrl%laser==2) then
        H_temp=traj%H_MCH_ss
        do idir=1,3
          H_temp=H_temp - traj%DM_ssd(:,:,idir)*real(ctrl%laserfield_td(traj%step*ctrl%nsubsteps+1,idir))
        enddo
        call diagonalize(ctrl%nstates,H_temp,U_temp)
      endif
    elseif (ctrl%surf==1) then
      U_temp=traj%U_ss
    endif
    if (printlevel>4) call matwrite(ctrl%nstates,U_temp,u_log,'U_ss','F12.9')


    ! ===============================
    ! 1. Compute time derivative Hamiltonian matrix (Kmatrix, TDH matrix)
    ! ===============================

    if (printlevel>3) then
      write(u_log,*) '============================================================='
      write(u_log,*) '    [1].Time derivative Hamiltonian matrix - Kmatrix'
      write(u_log,*) '============================================================='
    endif

    ! Compute TDH matrix
    gv(:)=0.d0
    do iatom=1,ctrl%natom
      do idir=1,3
        gv(:)=gv(:)+traj%grad_MCH_sad(:,iatom,idir)*traj%veloc_app_ad(iatom,idir)
      enddo
    enddo
    ! diagonal elements are time derivative of potential energies
    do istate=1, ctrl%nstates
      K1matrix_MCH_ss(istate,istate)=gv(istate)
    enddo
    ! off-diagonal elements are TDCs
    do istate=1, ctrl%nstates
      do jstate=1, ctrl%nstates
        if (istate .ne. jstate) then
          K1matrix_MCH_ss(istate,jstate)=&
          &-(traj%H_MCH_ss(istate,istate)-traj%H_MCH_ss(jstate,jstate))&
          &*traj%NACdt_ss(istate,jstate)
        endif
      enddo
    enddo
    ! transform K matrix to diagonal basis
    K1matrix_diag_ss=K1matrix_MCH_ss
    call transform(ctrl%nstates,K1matrix_diag_ss,U_temp,'utau')

    if (printlevel>4) then
      call matwrite(ctrl%nstates,K1matrix_MCH_ss,u_log,'MCH time derivative Hamiltonian matrix','F12.9')
      call matwrite(ctrl%nstates,K1matrix_diag_ss,u_log,'diagonal time derivative Hamiltonian matrix','F12.9')
    endif


    ! ===============================
    ! 1.1 Now we can assign TDC in diagonal basis 
    ! ===============================
    traj%NACdt_diag_ss=dcmplx(0.d0,0.d0)
    do istate=1, ctrl%nstates
      do jstate=1, ctrl%nstates
        if (istate .ne. jstate) then
          if (traj%H_diag_ss(istate,istate)-traj%H_diag_ss(jstate,jstate) .ne. 0.d0) then 
            traj%NACdt_diag_ss(istate,jstate)=&
            &K1matrix_diag_ss(istate,jstate)/(-(traj%H_diag_ss(istate,istate)-traj%H_diag_ss(jstate,jstate))) 
          else 
            traj%NACdt_diag_ss(istate,jstate)=&
            &K1matrix_diag_ss(istate,jstate)/1.d-8
          endif
        endif
      enddo
    enddo
 
    if (printlevel>3)  call matwrite(ctrl%nstates,traj%NACdt_diag_ss,u_log,'Time Derivative Coupling in diag basis...','F12.9')

    ! ===============================
    ! 2. Compute nuclear gradient Hamiltonian matrix (Gmatrix, NGH matrix)
    ! ===============================

    if (printlevel>3) then
      write(u_log,*) '============================================================='
      write(u_log,*) '    [2].Nuclear gradient Hamiltonian matrix - Gmatrix'
      write(u_log,*) '============================================================='
    endif

    ! Notice that there are two such matrices. 
    ! 1. diagonal are gradients, off-diagonal are NACs  -> G1matrix
    ! 2. diagonal are gradients, off-diagonal are effective NACs  -> G2matrix

    G1matrix_ssad=dcmplx(0.d0,0.d0)
    G2matrix_ssad=dcmplx(0.d0,0.d0)

    ! loop over all atom displacements
    do iatom=1,ctrl%natom
      do idir=1,3

        ! initialize G matrix, MCH gradients on diagonal
        G1matrix_ss=dcmplx(0.d0,0.d0)
        G2matrix_ss=dcmplx(0.d0,0.d0)
        do istate=1,ctrl%nstates
          G1matrix_ss(istate,istate)=traj%grad_MCH_sad(istate,iatom,idir)
          G2matrix_ss(istate,istate)=traj%grad_MCH_sad(istate,iatom,idir)
        enddo
 
        ! if NACs are available, add non-adiabatic coupling terms
        if ( (ctrl%calc_nacdr==0).or.(ctrl%calc_nacdr==1) ) then
          do istate=1,ctrl%nstates
            do jstate=istate+1,ctrl%nstates
              G1matrix_ss(istate,jstate)=&
               &-(traj%H_MCH_ss(istate,istate)-traj%H_MCH_ss(jstate,jstate))&
               &*traj%NACdr_ssad(istate,jstate,iatom,idir)
              G1matrix_ss(jstate,istate)=G1matrix_ss(istate,jstate)
            enddo
          enddo
        endif
        ! if effective NACs are available, add effective NAC terms
        if (ctrl%calc_effectivenac==1) then
          do istate=1,ctrl%nstates
            do jstate=istate+1,ctrl%nstates
              G2matrix_ss(istate,jstate)=&
               &-(traj%H_MCH_ss(istate,istate)-traj%H_MCH_ss(jstate,jstate))&
               &*traj%NACGV_MCH_ssad(istate,jstate,iatom,idir)
              G2matrix_ss(jstate,istate)=G2matrix_ss(istate,jstate)
            enddo
          enddo
        endif

        ! if available, add dipole moment derivatives
        if (ctrl%dipolegrad==1) then
          do istate=1,ctrl%nstates
            do jstate=1,ctrl%nstates
              do ipol=1,3
                G1matrix_ss(istate,jstate)=G1matrix_ss(istate,jstate)-&
                &traj%DMgrad_ssdad(istate,jstate,ipol,iatom,idir)*ctrl%laserfield_td(traj%step*ctrl%nsubsteps+1,ipol)
                G2matrix_ss(istate,jstate)=G2matrix_ss(istate,jstate)-&
                &traj%DMgrad_ssdad(istate,jstate,ipol,iatom,idir)*ctrl%laserfield_td(traj%step*ctrl%nsubsteps+1,ipol)
              enddo
            enddo
          enddo
        endif
  
        if (printlevel>6) then
          write(u_log,'(A,1X,I4,1X,A,1X,I4)') 'Gmatrix calculation... iatom=',iatom,'idir=',idir
          call matwrite(ctrl%nstates,G1matrix_ss,u_log,'G1matrix MCH','F12.9')
          call matwrite(ctrl%nstates,G2matrix_ss,u_log,'G2matrix MCH','F12.9')
          write(u_log,*)
        endif

        ! transform G matrix to diagonal basis
        call transform(ctrl%nstates,G1matrix_ss,U_temp,'utau')
        call transform(ctrl%nstates,G2matrix_ss,U_temp,'utau')

        ! save full G matrix
        G1matrix_ssad(:,:,iatom,idir)=G1matrix_ss
        G2matrix_ssad(:,:,iatom,idir)=G2matrix_ss

        if (printlevel>6) then
          call matwrite(ctrl%nstates,G1matrix_ss,u_log,'G1matrix diag','F12.9')
          call matwrite(ctrl%nstates,G2matrix_ss,u_log,'G2matrix diag','F12.9')
          write(u_log,*)
        endif

      enddo !do idir=1,3
    enddo !do iatom=1,ctrl%natom

    ! ===============================
    ! 2.1 Now we can assign NAC and effective NAC in diagonal basis 
    ! ===============================
    if ( (ctrl%calc_nacdr==0).or.(ctrl%calc_nacdr==1) ) then
      traj%NACdR_diag_ssad=dcmplx(0.d0,0.d0)
      do istate=1, ctrl%nstates
        do jstate=1, ctrl%nstates
          if (istate .ne. jstate) then
            if (traj%H_diag_ss(istate,istate)-traj%H_diag_ss(jstate,jstate) .ne. 0.d0) then
              traj%NACdR_diag_ssad(istate,jstate,:,:)=&
              &G1matrix_ssad(istate,jstate,:,:)/(-(traj%H_diag_ss(istate,istate)-traj%H_diag_ss(jstate,jstate)))
            else
              traj%NACdR_diag_ssad(istate,jstate,:,:)=&
              &G1matrix_ssad(istate,jstate,:,:)/1.d-8
            endif
          endif
        enddo
      enddo
      if (printlevel>4) then
        do istate=1,ctrl%nstates
          do jstate=1,ctrl%nstates
            write(string,'(A46,I3,1X,I3)') 'Non-adiabatic coupling DDR (diag basis) state ',istate,jstate
            call vec3write(ctrl%natom,traj%NACdR_diag_ssad(istate,jstate,:,:),u_log,trim(string),'F12.9')
          enddo
        enddo
      endif
    endif ! if ( (ctrl%calc_nacdr==0).or.(ctrl%calc_nacdr==1) ) then

    if (ctrl%calc_effectivenac==1) then
      traj%NACGV_diag_ssad=dcmplx(0.d0,0.d0)
      do istate=1, ctrl%nstates
        do jstate=1, ctrl%nstates
          if (istate .ne. jstate) then
            if (traj%H_diag_ss(istate,istate)-traj%H_diag_ss(jstate,jstate) .ne. 0.d0) then
              traj%NACGV_diag_ssad(istate,jstate,:,:)=&
              &G2matrix_ssad(istate,jstate,:,:)/(-(traj%H_diag_ss(istate,istate)-traj%H_diag_ss(jstate,jstate)))
            else
              traj%NACGV_diag_ssad(istate,jstate,:,:)=&
              &G2matrix_ssad(istate,jstate,:,:)/1.d-8
            endif
          endif
        enddo
      enddo
      if (printlevel>4) then
        do istate=1,ctrl%nstates
          do jstate=1,ctrl%nstates
            write(string,'(A56,I3,1X,I3)') 'Effective non-adiabatic coupling DDR (diag basis) state ',istate,jstate
            call vec3write(ctrl%natom,traj%NACGV_diag_ssad(istate,jstate,:,:),u_log,trim(string),'F12.9')
          enddo
        enddo
      endif
    endif !if (ctrl%calc_effectivenac==1) then


    ! ===============================
    ! 3. perform a NAC projection, project out rotational and translational motions from NAC
    ! ===============================
 
    if (ctrl%nac_projection==1) then
      if (printlevel>3) then
        write(u_log,*) '============================================================='
        write(u_log,*) '               [3].Performing NAC projection'
        write(u_log,*) '============================================================='
      endif

      ! Compute projection operator.
      traj%trans_rot_P=0.d0
      ctrans_rot_P=dcmplx(0.d0,0.d0)
      call compute_projection(traj%geom_ad, ctrl%natom, traj%trans_rot_P)
      if (printlevel>5) then
        call matwrite(3*ctrl%natom,traj%trans_rot_P,u_log,'Translational and Rotational Project Operator','F14.9')
      endif
      ctrans_rot_P=traj%trans_rot_P 

      ! initialize
      traj%pNACdR_MCH_ssad=0.d0
      traj%pNACdR_diag_ssad=dcmplx(0.d0,0.d0)
      traj%pNACGV_MCH_ssad=0.d0
      traj%pNACGV_diag_ssad=dcmplx(0.d0,0.d0)
  
      if ( (ctrl%calc_nacdr==0).or.(ctrl%calc_nacdr==1) ) then
        ! Project NACdR_MCH and NACdR_diag
        NACtmp_MCH=0.d0
        pNACtmp_MCH=0.d0
        NACtmp_diag=dcmplx(0.d0,0.d0)
        pNACtmp_diag=dcmplx(0.d0,0.d0)
        do istate=1,ctrl%nstates
          do jstate=1,ctrl%nstates
            do iatom=1,ctrl%natom
              do idir=1,3
                i=3*(iatom-1)+idir
                NACtmp_MCH(i)=traj%NACdr_ssad(istate,jstate,iatom,idir)
                NACtmp_diag(i)=traj%NACdR_diag_ssad(istate,jstate,iatom,idir)
              enddo
            enddo
            call matvecmultiply(3*ctrl%natom, traj%trans_rot_P, NACtmp_MCH, pNACtmp_MCH, 'n')
            call matvecmultiply(3*ctrl%natom, ctrans_rot_P, NACtmp_diag, pNACtmp_diag, 'n')
            do iatom=1,ctrl%natom
              do idir=1,3
                i=3*(iatom-1)+idir
                traj%pNACdR_MCH_ssad(istate,jstate,iatom,idir)=pNACtmp_MCH(i)
                traj%pNACdR_diag_ssad(istate,jstate,iatom,idir)=pNACtmp_diag(i)
              enddo
            enddo
          enddo
        enddo

        if (printlevel>4) then
          do i=1,ctrl%nstates
            do j=1,ctrl%nstates
              write(string,'(A55,I3,1X,I3)') 'Projected Non-adiabatic coupling DDR (MCH basis) state ',i,j
              call vec3write(ctrl%natom,traj%pNACdr_MCH_ssad(i,j,:,:),u_log,trim(string),'F9.4')
            enddo
          enddo
          write(u_log,*)
          do i=1,ctrl%nstates
            do j=1,ctrl%nstates
              write(string,'(A56,I3,1X,I3)') 'Projected Non-adiabatic coupling DDR (diag basis) state ',i,j
              call vec3write(ctrl%natom,traj%pNACdr_diag_ssad(i,j,:,:),u_log,trim(string),'F9.4')
            enddo
          enddo
          write(u_log,*)
        endif
      endif !if ( (ctrl%calc_nacdr==0).or.(ctrl%calc_nacdr==1) ) then
 
      if (ctrl%calc_effectivenac==1) then !one employs effective NAC in nuclear EOM. 
        NACtmp_MCH=0.d0
        pNACtmp_MCH=0.d0
        NACtmp_diag=dcmplx(0.d0,0.d0)
        pNACtmp_diag=dcmplx(0.d0,0.d0)
        do istate=1,ctrl%nstates
          do jstate=1,ctrl%nstates
            do iatom=1,ctrl%natom
              do idir=1,3
                i=3*(iatom-1)+idir
                NACtmp_MCH(i)=traj%NACGV_MCH_ssad(istate,jstate,iatom,idir)
                NACtmp_diag(i)=traj%NACGV_diag_ssad(istate,jstate,iatom,idir)
              enddo
            enddo
            call matvecmultiply(3*ctrl%natom, traj%trans_rot_P, NACtmp_MCH, pNACtmp_MCH, 'n')
            call matvecmultiply(3*ctrl%natom, ctrans_rot_P, NACtmp_diag, pNACtmp_diag, 'n')
            do iatom=1,ctrl%natom
              do idir=1,3
                i=3*(iatom-1)+idir
                traj%pNACGV_MCH_ssad(istate,jstate,iatom,idir)=pNACtmp_MCH(i)
                traj%pNACGV_diag_ssad(istate,jstate,iatom,idir)=pNACtmp_diag(i)
              enddo
            enddo
          enddo
        enddo

        if (printlevel>4) then
          do i=1,ctrl%nstates
            do j=1,ctrl%nstates
              write(string,'(A65,I3,1X,I3)') 'Projected effective Non-adiabatic coupling DDR (MCH basis) state ',i,j
              call vec3write(ctrl%natom,traj%pNACGV_MCH_ssad(i,j,:,:),u_log,trim(string),'F9.4')
            enddo
          enddo
          write(u_log,*)
          do i=1,ctrl%nstates
            do j=1,ctrl%nstates
              write(string,'(A66,I3,1X,I3)') 'Projected effective Non-adiabatic coupling DDR (diag basis) state ',i,j
              call vec3write(ctrl%natom,traj%pNACGV_diag_ssad(i,j,:,:),u_log,trim(string),'F9.4')
            enddo
          enddo
          write(u_log,*)
        endif
 
      endif !if (ctrl%calc_effectivenac==1) then

    endif !if (ctrl%nac_projection==1) then

 

    ! ===============================
    ! 4. Patch Kmatrix. 
    ! ===============================

    if (printlevel>3) then
      write(u_log,*) '============================================================='
      write(u_log,*) '                [4].Patch Kmatrix'
      write(u_log,*) '============================================================='
    endif

    ! Compute Kmatrix in MCH representation first.
    ! There are two ways to compute dH/dt, gradient or energy difference
    traj%Kmatrix_MCH_ss=dcmplx(0.d0,0.d0)
    if (ctrl%kmatrix_method==0) then

      if (printlevel>4) write(u_log,*) 'Time derivative matrix (Kmatrix) is computed by gradient dot velocity'
      traj%Kmatrix_MCH_ss=K1matrix_MCH_ss
      traj%Kmatrix_diag_ss=K1matrix_diag_ss

    else if (ctrl%kmatrix_method==1) then
      if (printlevel>4) write(u_log,*) 'Time derivative matrix (Kmatrix) is computed by energy difference'
      traj%Kmatrix_MCH_ss=K1matrix_MCH_ss
      traj%Kmatrix_diag_ss=K1matrix_diag_ss

      if (ctrl%integrator==2) then ! fixed time step 
        if (traj%step==0) then 
          do istate=1, ctrl%nstates
            traj%Kmatrix_MCH_ss(istate,istate)=0.d0
            traj%Kmatrix_diag_ss(istate,istate)=0.d0
          enddo
        else if (traj%step==1) then 
          do istate=1, ctrl%nstates
            traj%Kmatrix_MCH_ss(istate,istate)=(traj%H_MCH_ss(istate,istate)-traj%H_MCH_old_ss(istate,istate))/ctrl%dtstep
            traj%Kmatrix_diag_ss(istate,istate)=(traj%H_diag_ss(istate,istate)-traj%H_diag_old_ss(istate,istate))/ctrl%dtstep
          enddo
        else if (traj%step>=2) then 
          do istate=1, ctrl%nstates
            traj%Kmatrix_MCH_ss(istate,istate)=(3*traj%H_MCH_ss(istate,istate)-4*traj%H_MCH_old_ss(istate,istate)+&
              traj%H_MCH_old2_ss(istate,istate))/(2*ctrl%dtstep)
            traj%Kmatrix_diag_ss(istate,istate)=(3*traj%H_diag_ss(istate,istate)-4*traj%H_diag_old_ss(istate,istate)+&
              traj%H_diag_old2_ss(istate,istate))/(2*ctrl%dtstep)
          enddo
        !else if (traj%step>=3) then
        !  do istate=1, ctrl%nstates
        !    traj%Kmatrix_MCH_ss(istate,istate)=(29*traj%H_MCH_ss(istate,istate)-54*traj%H_MCH_old_ss(istate,istate)+&
        !      27*traj%H_MCH_old2_ss(istate,istate)-2*traj%H_MCH_old3_ss(istate,istate))/(6*ctrl%dtstep)
        !    traj%Kmatrix_diag_ss(istate,istate)=(29*traj%H_diag_ss(istate,istate)-54*traj%H_diag_old_ss(istate,istate)+&
        !      27*traj%H_diag_old2_ss(istate,istate)-2*traj%H_diag_old3_ss(istate,istate))/(6*ctrl%dtstep)
        !  enddo
        endif
      else if (ctrl%integrator==1) then ! adaptive time step, only use first order finite difference
        if (traj%step==0) then 
          do istate=1, ctrl%nstates
            traj%Kmatrix_MCH_ss(istate,istate)=1.d0
            traj%Kmatrix_diag_ss(istate,istate)=1.d0
          enddo
        else if (traj%step>=1) then 
          do istate=1, ctrl%nstates
            traj%Kmatrix_MCH_ss(istate,istate)=(traj%H_MCH_ss(istate,istate)-traj%H_MCH_old_ss(istate,istate))/ctrl%dtstep
            traj%Kmatrix_diag_ss(istate,istate)=(traj%H_diag_ss(istate,istate)-traj%H_diag_old_ss(istate,istate))/ctrl%dtstep
          enddo
        endif 
      endif
    endif ! (ctrl%kmatrix_method==0) then
 
    if (printlevel>4) then
      call matwrite(ctrl%nstates,traj%Kmatrix_MCH_ss,u_log,'Final patched Kmatrix MCH','F12.9')
      call matwrite(ctrl%nstates,traj%Kmatrix_diag_ss,u_log,'Final patched Kmatrix diag','F12.9')
    endif


    ! ===============================
    ! 5. Patch Gmatrix. 
    ! ===============================

    if (printlevel>3) then
      write(u_log,*) '============================================================='
      write(u_log,*) '          [5].Patch mixed gradient matrix'
      write(u_log,*) '============================================================='
    endif

    if (ctrl%gradcorrect==0 .or. ctrl%gradcorrect==1) then 
      ! loop over all atom displacements
      do iatom=1,ctrl%natom
        do idir=1,3

          ! initialize G matrix, MCH gradients on diagonal
          Gmatrix_ss=dcmplx(0.d0,0.d0)
          do istate=1,ctrl%nstates
            Gmatrix_ss(istate,istate)=traj%grad_MCH_sad(istate,iatom,idir)
          enddo

          ! if available, add non-adiabatic coupling terms
          if (ctrl%gradcorrect==1) then
            do istate=1,ctrl%nstates
              do jstate=istate+1,ctrl%nstates
                Gmatrix_ss(istate,jstate)=&
                 &-(traj%H_MCH_ss(istate,istate)-traj%H_MCH_ss(jstate,jstate))&
                 &*traj%NACdr_ssad(istate,jstate,iatom,idir)
                Gmatrix_ss(jstate,istate)=Gmatrix_ss(istate,jstate)
              enddo
            enddo
          endif

          ! if available, add dipole moment derivatives
          if (ctrl%dipolegrad==1) then
            do istate=1,ctrl%nstates
              do jstate=1,ctrl%nstates
                do ipol=1,3
                  Gmatrix_ss(istate,jstate)=Gmatrix_ss(istate,jstate)-&
                  &traj%DMgrad_ssdad(istate,jstate,ipol,iatom,idir)*ctrl%laserfield_td(traj%step*ctrl%nsubsteps+1,ipol)
                enddo
              enddo
            enddo
          endif

          if (printlevel>4) then
            write(u_log,'(A,1X,I4,1X,A,1X,I4)') 'Gmatrix calculation... iatom=',iatom,'idir=',idir
            call matwrite(ctrl%nstates,Gmatrix_ss,u_log,'Gmatrix MCH','F12.9')
          endif

          ! transform G matrix to diagonal basis
          call transform(ctrl%nstates,Gmatrix_ss,U_temp,'utau')
          ! save full G matrix in traj
          traj%Gmatrix_ssad(:,:,iatom,idir)=Gmatrix_ss

          if (printlevel>4) then
            call matwrite(ctrl%nstates,Gmatrix_ss,u_log,'Gmatrix diag','F12.9')
            write(u_log,*)
          endif

        enddo !do idir=1,3
      enddo !do iatom=1,ctrl%natom

    elseif (ctrl%gradcorrect==2) then
      ! assign diagonal gradient
      do istate=1,ctrl%nstates
        jstate=state_diag_to_MCH(ctrl%nstates,istate,U_temp)
        if (printlevel>4) then
          write(u_log,'(A,1X,I4,1X,A,1X,I4)') ' gradient of diagonal state',istate, 'is approximated by gradient of MCH state', jstate
        endif
        if (ctrl%kmatrix_method==1 .and. traj%step==0) then 
          traj%Gmatrix_ssad(istate,istate,:,:)=traj%grad_MCH_sad(jstate,:,:)
        else 
          traj%Gmatrix_ssad(istate,istate,:,:)=traj%Kmatrix_diag_ss(istate,istate)/traj%Kmatrix_MCH_ss(jstate,jstate)*traj%grad_MCH_sad(jstate,:,:)
        endif
      enddo
      ! assign off diagonal elements, which are effective NACs
      do istate=1,ctrl%nstates
        do jstate=1,ctrl%nstates
          if (istate .ne. jstate) then
            traj%Gmatrix_ssad(istate,jstate,:,:)=traj%NACGV_diag_ssad(istate,jstate,:,:)
          endif
        enddo
      enddo

      ! finally, add diagonal contribution from dipole gradients. Notice this has to be diagonal contribution
      if (ctrl%dipolegrad==1) then
        ! loop over all atom displacements
        do iatom=1,ctrl%natom
          do idir=1,3
            Gmatrix_ss=dcmplx(0.d0,0.d0)
            do istate=1,ctrl%nstates
              do jstate=1,ctrl%nstates
                do ipol=1,3
                  Gmatrix_ss(istate,jstate)=Gmatrix_ss(istate,jstate)-&
                  &traj%DMgrad_ssdad(istate,jstate,ipol,iatom,idir)*ctrl%laserfield_td(traj%step*ctrl%nsubsteps+1,ipol)
                enddo
              enddo
            enddo
            ! transform dipole gradient matrix to diagonal basis
            call transform(ctrl%nstates,Gmatrix_ss,U_temp,'utau')
            ! add dipole gradient part to final Gmatrix
            traj%Gmatrix_ssad(:,:,iatom,idir)=traj%Gmatrix_ssad(:,:,iatom,idir)+Gmatrix_ss(:,:)
          enddo
        enddo
      endif
      ! print Gmatrix
      if (printlevel>4) then
        do iatom=1,ctrl%natom
          do idir=1,3
            write(u_log,'(A,1X,I4,1X,A,1X,I4)') 'Gmatrix calculation... iatom=',iatom,'idir=',idir
            call matwrite(ctrl%nstates,traj%Gmatrix_ssad(:,:,iatom,idir),u_log,'Gmatrix diag','F12.9')
          enddo
        enddo
      endif

    endif ! if ((ctrl%gradcorrect==0 .or. ctrl%gradcorrect==1) then


    ! ===============================
    ! 6. Compute Hopping and Frustratede hop direction 
    ! ===============================

    if (ctrl%method==0) then ! TSH method

      if (printlevel>3) then
        write(u_log,*) '============================================================='
        write(u_log,*) '    [6].Caculating Hopping and Frustrated Hop  Direction '
        write(u_log,*) '============================================================='
      endif

      !====================
      ! Hopping direction 
      !====================

      select case (ctrl%ekincorrect)
        case (0)    ! no frustrated jumps: just go on
          if (printlevel>4) write(u_log,*) 'No hopping direction needed'
          continue

        case (1)    ! correct along v, full kinetic energy available
          if (printlevel>4) write(u_log,*) 'Hopping direction set to velocity vector'
          do istate=1,ctrl%nstates
            do jstate=1,ctrl%nstates
              traj%hopping_direction_ssad(istate,jstate,:,:)=traj%veloc_ad(:,:)
            enddo
          enddo

        case (2)    ! correct along projected v, vibrational kinetic energy available
          if (printlevel>4) write(u_log,*) 'Hopping direction set to projected velocity vector'
          do istate=1,ctrl%nstates
            do jstate=1,ctrl%nstates
              traj%hopping_direction_ssad(istate,jstate,:,:)=traj%veloc_ad(:,:)
            enddo
          enddo

        case (3)    ! correct along NAC
          if (printlevel>4) write(u_log,*) 'Hopping direction set to NAC'
          do istate=1,ctrl%nstates
            do jstate=1,ctrl%nstates
              traj%hopping_direction_ssad(istate,jstate,:,:)=real(traj%Gmatrix_ssad(istate,jstate,:,:))
            enddo
          enddo

        case (4)    ! correct along gradient difference
          if (printlevel>4) write(u_log,*) 'Hopping direction set to difference gradient'
          do istate=1,ctrl%nstates
            do jstate=1,ctrl%nstates
              traj%hopping_direction_ssad(istate,jstate,:,:)=real(traj%Gmatrix_ssad(istate,istate,:,:)-&
                &traj%Gmatrix_ssad(jstate,jstate,:,:))
            enddo
          enddo
 
        case (5)    ! correct along projected NAC/M
          if (printlevel>4) write(u_log,*) 'Hopping direction set to projected NAC'
          do istate=1,ctrl%nstates
            do jstate=1,ctrl%nstates
              traj%hopping_direction_ssad(istate,jstate,:,:)=real(traj%Gmatrix_ssad(istate,jstate,:,:))
            enddo
          enddo

        case (6)    ! correct along projected gradient difference
          if (printlevel>4) write(u_log,*) 'Hopping direction set to projected difference gradient'
          do istate=1,ctrl%nstates
            do jstate=1,ctrl%nstates
              traj%hopping_direction_ssad(istate,jstate,:,:)=real(traj%Gmatrix_ssad(istate,istate,:,:)-&
                &traj%Gmatrix_ssad(jstate,jstate,:,:))
            enddo
          enddo

        case (7)    ! correct along effective NAC/M
          if (printlevel>4) write(u_log,*) 'Hopping direction set to effective NAC'
          do istate=1,ctrl%nstates
            do jstate=1,ctrl%nstates
              traj%hopping_direction_ssad(istate,jstate,:,:)=real(traj%NACGV_diag_ssad(istate,jstate,:,:))
            enddo
          enddo

        case (8)    ! correct along projected effective NAC/M
          if (printlevel>4) write(u_log,*) 'Hopping direction set to projected effective NAC'
          do istate=1,ctrl%nstates
            do jstate=1,ctrl%nstates
              traj%hopping_direction_ssad(istate,jstate,:,:)=real(traj%NACGV_diag_ssad(istate,jstate,:,:))
            enddo
          enddo

      endselect

      ! for projected ones, we should do projection 
      if (ctrl%ekincorrect==2 .or. ctrl%ekincorrect==5 .or. ctrl%ekincorrect==6 .or. ctrl%ekincorrect==8) then
        if (printlevel>4) write(u_log,*) 'hopping direction projected'

        hopping_tmp=0.d0
        do istate=1,ctrl%nstates
          do jstate=1,ctrl%nstates
            do iatom=1,ctrl%natom
              do idir=1,3
                i=3*(iatom-1)+idir
                hopping_tmp(i)=traj%hopping_direction_ssad(istate,jstate,iatom,idir)
              enddo
            enddo
            call matvecmultiply(3*ctrl%natom, traj%trans_rot_P, hopping_tmp, phopping_tmp, 'n')
            do iatom=1,ctrl%natom
              do idir=1,3
                i=3*(iatom-1)+idir
                traj%hopping_direction_ssad(istate,jstate,iatom,idir)=phopping_tmp(i)
              enddo
            enddo
          enddo
        enddo

      endif 

      !====================
      ! Frustrated hop direction 
      !====================

      if (ctrl%reflect_frustrated==0) then 
        if (printlevel>4) write(u_log,*) 'Frustrated hops are ignored'
        continue
      else if (ctrl%reflect_frustrated==1 .or. ctrl%reflect_frustrated==91) then    
        ! reverse along v or delV approach with reverse along v
        if (printlevel>4) write(u_log,*) 'Frustrated hop direction set to reverse along velocity vector'
        do istate=1,ctrl%nstates
          do jstate=1,ctrl%nstates
            traj%frustrated_hop_vec_ssad(istate,jstate,:,:)=traj%veloc_ad(:,:)
          enddo
        enddo
      else if (ctrl%reflect_frustrated==2 .or. ctrl%reflect_frustrated==92) then
        ! reverse along projected v or delV approach with reverse along projected v
        if (printlevel>4) write(u_log,*) 'Frustrated hop direction set to reverse along projected velocity vector'
        do istate=1,ctrl%nstates
          do jstate=1,ctrl%nstates
            traj%frustrated_hop_vec_ssad(istate,jstate,:,:)=traj%veloc_ad(:,:)
          enddo
        enddo
      else if (ctrl%reflect_frustrated==3 .or. ctrl%reflect_frustrated==93) then
        ! reverse along NAC or delV approach with reverse along NAC
        if (printlevel>4) write(u_log,*) 'Frustrated hop direction set to reverse along NAC'
        do istate=1,ctrl%nstates
          do jstate=1,ctrl%nstates
            traj%frustrated_hop_vec_ssad(istate,jstate,:,:)=real(traj%Gmatrix_ssad(istate,jstate,:,:))
          enddo
        enddo
      else if (ctrl%reflect_frustrated==4 .or. ctrl%reflect_frustrated==94) then
        ! reverse along NAC or delV approach with reverse along gradient difference
        if (printlevel>4) write(u_log,*) 'Frustrated hop direction set to reverse along difference gradient'
        do istate=1,ctrl%nstates
          do jstate=1,ctrl%nstates
            traj%frustrated_hop_vec_ssad(istate,jstate,:,:)=real(traj%Gmatrix_ssad(istate,istate,:,:)-&
              &traj%Gmatrix_ssad(jstate,jstate,:,:))
          enddo
        enddo
      else if (ctrl%reflect_frustrated==5 .or. ctrl%reflect_frustrated==95) then
        ! reverse along NAC or delV approach with reverse along projected NAC
        if (printlevel>4) write(u_log,*) 'Frustrated hop direction set to reverse along projected NAC'
        do istate=1,ctrl%nstates
          do jstate=1,ctrl%nstates
            traj%frustrated_hop_vec_ssad(istate,jstate,:,:)=real(traj%Gmatrix_ssad(istate,jstate,:,:))
          enddo
        enddo
      else if (ctrl%reflect_frustrated==6 .or. ctrl%reflect_frustrated==96) then
        ! reverse along NAC or delV approach with reverse along projected gradient difference
        if (printlevel>4) write(u_log,*) 'Frustrated hop direction set to reverse along projected difference gradient'
        do istate=1,ctrl%nstates
          do jstate=1,ctrl%nstates
            traj%frustrated_hop_vec_ssad(istate,jstate,:,:)=real(traj%Gmatrix_ssad(istate,istate,:,:)-&
              &traj%Gmatrix_ssad(jstate,jstate,:,:))
          enddo
        enddo
      else if (ctrl%reflect_frustrated==7 .or. ctrl%reflect_frustrated==97) then
        ! reverse along NAC or delV approach with reverse along effective NAC
        if (printlevel>4) write(u_log,*) 'Frustrated hop direction set to reverse along effective NAC'
        do istate=1,ctrl%nstates
          do jstate=1,ctrl%nstates
            traj%frustrated_hop_vec_ssad(istate,jstate,:,:)=real(traj%NACGV_diag_ssad(istate,jstate,:,:))
          enddo
        enddo
      else if (ctrl%reflect_frustrated==8 .or. ctrl%reflect_frustrated==98) then
        ! reverse along NAC or delV approach with reverse along projected effective NAC
        if (printlevel>4) write(u_log,*) 'Frustrated hop direction set to reverse along projected effective NAC'
        do istate=1,ctrl%nstates
          do jstate=1,ctrl%nstates
            traj%frustrated_hop_vec_ssad(istate,jstate,:,:)=real(traj%NACGV_diag_ssad(istate,jstate,:,:))
          enddo
        enddo
      endif

      ! for projected ones, we should do projection 
      if (ctrl%reflect_frustrated==2 .or. ctrl%reflect_frustrated==5 .or. ctrl%reflect_frustrated==6 .or. ctrl%reflect_frustrated==8 .or. &
         &ctrl%reflect_frustrated==92 .or. ctrl%reflect_frustrated==95 .or. ctrl%reflect_frustrated==96 .or. ctrl%reflect_frustrated==98) then
        if (printlevel>4) write(u_log,*) 'Frustrated hop direction projected'
        hopping_tmp=0.d0
        do istate=1,ctrl%nstates
          do jstate=1,ctrl%nstates
            do iatom=1,ctrl%natom
              do idir=1,3
                i=3*(iatom-1)+idir
                hopping_tmp(i)=traj%frustrated_hop_vec_ssad(istate,jstate,iatom,idir)
              enddo
            enddo
            call matvecmultiply(3*ctrl%natom, traj%trans_rot_P, hopping_tmp, phopping_tmp, 'n')
            do iatom=1,ctrl%natom
              do idir=1,3
                i=3*(iatom-1)+idir
                traj%frustrated_hop_vec_ssad(istate,jstate,iatom,idir)=phopping_tmp(i)
              enddo
            enddo
          enddo
        enddo

      endif

    endif ! if (ctrl%method==0)

  endsubroutine

! ===========================================================


! ===========================================================

! Subroutine compute_projection, used to compute the projection matrix
! The projection matrix projects out the translational and rotational part
  subroutine compute_projection(geom, n, P)
  use definitions, only: u_log
  use matrix
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: geom(n,3)
  real*8, intent(out) :: P(3*n,3*n)

  integer :: iatom, jatom, idir, jdir
  integer :: ia, ib, ic, ja, jb, jc, jend
  integer :: i, j, ii, jj
  real*8 :: sym(3,3,3)
  real*8 :: mgeom(n,3),x,y,z
  real*8 :: centroid(3), momi(3,3), invmomi(3,3), detinv
  real*8 :: prodmomi,det,tmp
  real*8 :: sumall

! initialize Levi-Civita 3X3X3 pseudo tensor
  sym=0.d0
  sym(1,2,3)=1.0d0
  sym(1,3,2)=-1.0d0
  sym(2,1,3)=-1.0d0
  sym(2,3,1)=1.0d0
  sym(3,1,2)=1.0d0
  sym(3,2,1)=-1.0d0

! compute the centroid
  centroid=0.d0
  do iatom=1,n
    do idir=1,3
      centroid(idir)=centroid(idir)+geom(iatom,idir)
    enddo
  enddo
  centroid=centroid/n

! compute deltaR
  do iatom=1,n
    do idir=1,3
      mgeom(iatom,idir)=geom(iatom,idir)-centroid(idir)
    enddo
  enddo

! compute moment of inertia
  momi=0.d0
  do iatom=1,n
    x=mgeom(iatom,1)
    y=mgeom(iatom,2)
    z=mgeom(iatom,3)
    momi(1,1)=momi(1,1)+y*y+z*z
    momi(2,2)=momi(2,2)+x*x+z*z
    momi(3,3)=momi(3,3)+x*x+y*y
    momi(1,2)=momi(1,2)-x*y
    momi(1,3)=momi(1,3)-x*z
    momi(2,3)=momi(2,3)-y*z
  enddo
  momi(2,1)=momi(1,2)
  momi(3,1)=momi(1,3)
  momi(3,2)=momi(2,3)

! Invert the moment of inertia matrix
! check moment of inertia with zeros
  prodmomi=abs(momi(1,1)*momi(2,2)*momi(3,3))
  if (prodmomi .le. 1.d-8) then
    ! if X=0, Y=0, Z=0
    if (abs(momi(1,1)) .le. 1.d-8 .and. abs(momi(2,2)) .le. 1.d-8 .and. abs(momi(3,3)) .le. 1.d-8) then 
      write(u_log, *) ' Warning: all digaonal elements of moment of inertia equals zero!'
      write(u_log, '(3(F20.10))')  momi(1,1), momi(2,2), momi(3,3)
      stop 1
    ! if X=0, Y.NE.0, Z.NE.0
    else if (abs(momi(1,1)) .le. 1.d-8 .and. abs(momi(2,2)) .gt. 1.d-8 .and. abs(momi(3,3)) .gt. 1.d-8) then
      det=momi(3,3)*momi(2,2)-momi(3,2)*momi(2,3)
      tmp=momi(3,3)
      momi(3,3)=momi(2,2)/det
      momi(2,2)=tmp/det
      momi(3,2)=-momi(3,2)/det
      momi(2,3)=-momi(2,3)/det
    ! if X.NE.0, Y=0, Z.NE.0
    else if (abs(momi(1,1)) .gt. 1.d-8 .and. abs(momi(2,2)) .le. 1.d-8 .and. abs(momi(3,3)) .gt. 1.d-8) then
      det=momi(1,1)*momi(3,3)-momi(1,3)*momi(3,1)
      tmp=momi(1,1)
      momi(1,1)=momi(3,3)/det
      momi(3,3)=tmp/det
      momi(1,3)=-momi(1,3)/det
      momi(3,1)=-momi(3,1)/det
    ! if X.NE.0, Y.NE.0, Z=0
    else if (abs(momi(1,1)) .gt. 1.d-8 .and. abs(momi(2,2)) .gt. 1.d-8 .and. abs(momi(3,3)) .le. 1.d-8) then
      det=momi(1,1)*momi(2,2)-momi(1,2)*momi(2,1)
      tmp=momi(1,1)
      momi(1,1)=momi(2,2)/det
      momi(2,2)=tmp/det
      momi(1,2)=-momi(1,2)/det
      momi(2,1)=-momi(2,1)/det
    ! if X=0, Y=0, Z.NE.0
    else if (abs(momi(1,1)) .le. 1.d-8 .and. abs(momi(2,2)) .le. 1.d-8 .and. abs(momi(3,3)) .gt. 1.d-8) then
      momi(3,3)=1.d0/momi(3,3)
    ! if X=0, Y.NE.0, Z=0 
    else if (abs(momi(1,1)) .le. 1.d-8 .and. abs(momi(2,2)) .gt. 1.d-8 .and. abs(momi(3,3)) .le. 1.d-8) then
      momi(2,2)=1.d0/momi(2,2)
    ! if X.NE.0, Y=0, Z=0 
    else if (abs(momi(1,1)) .gt. 1.d-8 .and. abs(momi(2,2)) .le. 1.d-8 .and. abs(momi(3,3)) .le. 1.d-8) then
      momi(1,1)=1.d0/momi(1,1)
    endif
    invmomi=0.d0
    invmomi=momi
  else !if (prodmomi .le. 1.d-8) then
    invmomi=0.d0
    call mat3invert(momi, invmomi)
  endif


! now compute P matrix
!    ---------------- 
  do iatom=1,n
    do jatom=1,iatom
      i=3*(iatom-1)
      j=3*(jatom-1)
      do ic=1,3
        jend=3
        if (jatom.eq.iatom) jend=ic
        do jc=1,jend
          sumall=0.d0
          do ia=1,3
          do ib=1,3
            if (sym(ia,ib,ic).ne.0) then 
              do ja=1,3
              do jb=1,3
                if (sym(ja,jb,jc).ne.0) then 
                  sumall=sumall+sym(ia,ib,ic)*sym(ja,jb,jc)*invmomi(ia,ja)* &
                    mgeom(iatom,ib)*mgeom(jatom,jb)
                endif   
              enddo    
              enddo
            endif
          enddo
          enddo
          ii=i+ic
          jj=j+jc
          P(ii,jj)=sumall
          if (ic.eq.jc) P(ii,jj)=P(ii,jj)+1.d0/n
        enddo
      enddo
    enddo
  enddo
     
! compute I-P
  do i=1,3*n
    do j=1,i
      P(i,j)=-P(i,j)
      if (i.eq.j) P(i,j)=1+P(i,j)
    enddo
  enddo

! remove small values less than 1.0d-08 and assign the lower triangle of P
  do i=1,3*n
    do j=1,i
      if (abs(P(i,j)) .lt. 1.d-8) P(i,j)=0.d0
      P(j,i)=P(i,j)
    enddo
  enddo

  endsubroutine

endmodule

