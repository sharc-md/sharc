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

!> # Module QM
!> 
!> \author Sebastian mai
!> \date 27.02.2015
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
      call write_logtimestep(u_log,traj%step)
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

    ! using the U matrix, calculate the remaining coefficient vectors
    select case (ctrl%staterep)
      case (0)      ! coeff is in diag, transform to MCH for printing
        call matvecmultiply(ctrl%nstates,traj%U_ss,traj%coeff_diag_s,traj%coeff_MCH_s,'n')
        traj%state_MCH=state_diag_to_MCH(ctrl%nstates,traj%state_diag,traj%U_ss)
      case (1)      ! coeff is in MCH, transform to diag
        call matvecmultiply(ctrl%nstates,traj%U_ss,traj%coeff_MCH_s,traj%coeff_diag_s,'t')
        traj%state_diag=state_MCH_to_diag(ctrl%nstates,traj%state_MCH,traj%U_ss)
    endselect

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

    ! we have to set up the initial NACdt_ss matrix here
    if (ctrl%coupling==1) then
      traj%NACdt_ss=dcmplx(0.d0,0.d0)
      do iatom=1,ctrl%natom
        do idir=1,3
          traj%NACdt_ss=traj%NACdt_ss+traj%NACdr_ssad(:,:,iatom,idir)*traj%veloc_ad(iatom,idir)
        enddo
      enddo
!       call matwrite(ctrl%nstates,traj%NACdt_ss,u_log,'Old DDT Matrix','F12.9')
    endif

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
    integer :: stat,i,j

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
    ! apply frozen-state mask
    do i=1,ctrl%nstates
      do j=1,ctrl%nstates
        if (ctrl%actstates_s(i).neqv.ctrl%actstates_s(j)) traj%H_MCH_ss(i,j)=dcmplx(0.d0,0.d0)
        if ((ctrl%calc_soc/=1).and.(i/=j)) traj%H_MCH_ss(i,j)=dcmplx(0.d0,0.d0)
      enddo
    enddo
    if (printlevel>3) write(u_log,'(A31,A2)') 'Hamiltonian:                   ','OK'

    ! get Dipole moments
    call get_dipoles(ctrl%nstates, traj%DM_ssd)
    if (printlevel>3) write(u_log,'(A31,A2)') 'Dipole Moments:                ','OK'
    traj%DM_print_ssd=traj%DM_ssd
    ! apply frozen-state mask 
    do i=1,ctrl%nstates
      do j=1,ctrl%nstates
        if (ctrl%actstates_s(i).neqv.ctrl%actstates_s(j)) traj%DM_ssd(i,j,:)=dcmplx(0.d0,0.d0)
      enddo
    enddo

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
    write(u_qm_qmin,'(A)') 'DM'

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
    ! never calculate gradients of frozen states
    traj%selg_s=traj%selg_s.and.ctrl%actstates_s

    if (printlevel>3) then
      write(u_log,*) '-------------------- Gradient selection ---------------------'
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

    if (printlevel>3) then
      write(u_log,*) '------------- Non-adiabatic coupling selection --------------'
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

    do i=1,ctrl%nstates
      do j=1,ctrl%nstates
        traj%seldm_ss(i,j)=traj%seldm_ss(i,j).and.ctrl%actstates_s(i).and.ctrl%actstates_s(j)
      enddo
    enddo

    if (printlevel>3) then
      write(u_log,*) '------------- Dipole moment gradient selection --------------'
      if (traj%step==0) then
        write(u_log,*) 'Selecting all states in first timestep.'
      else
        write(u_log,*) 'Select nacs:'
        write(u_log,*) 'State(diag)= ',traj%state_diag,'State(MCH)=',traj%state_MCH
        write(u_log,*) 'Selected Coupled States(MCH)='
        do i=1,ctrl%nstates
          write(u_log,*) (traj%seldm_ss(i,j),j=1,ctrl%nstates)
        enddo
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
        call matwrite(ctrl%nstates,traj%overlaps_ss,u,'Overlap matrix (MCH basis)','F9.4')
      write(u,*)
      endif
    endif
    if (traj%phases_found) then
      call vecwrite(ctrl%nstates,traj%phases_s,u,'Wavefunction phases (MCH basis)','F9.4')
    endif

  endsubroutine

! =========================================================== !

!> Updates all *_old matrices when advancing to a new timestep
  subroutine Update_old(traj)
    use definitions
    implicit none
    type(trajectory_type) :: traj

    if (printlevel>3) then
      write(u_log,*) '============================================================='
      write(u_log,*) '                   Advancing to next step'
      write(u_log,*) '============================================================='
    endif
    ! initialize old variables from current ones
    traj%dH_MCH_ss=traj%H_MCH_ss-traj%H_MCH_old_ss
    traj%H_MCH_old_ss=traj%H_MCH_ss
    traj%DM_old_ssd=traj%DM_ssd
    traj%U_old_ss=traj%U_ss
    traj%NACdt_old_ss=traj%NACdt_ss
    traj%NACdr_old_ssad=traj%NACdr_ssad
    traj%phases_old_s=traj%phases_s

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
    integer :: iatom, idir, istate, jstate, ipol
    complex*16 :: Gmatrix_ss(ctrl%nstates,ctrl%nstates)
    complex*16 :: U_temp(ctrl%nstates,ctrl%nstates), H_temp(ctrl%nstates,ctrl%nstates)

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

    if (printlevel>3) then
      write(u_log,*) '============================================================='
      write(u_log,*) '             Calculating mixed gradient matrix'
      write(u_log,*) '============================================================='
    endif
    if (printlevel>4) call matwrite(ctrl%nstates,U_temp,u_log,'U_ss','F12.9')
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


      enddo
    enddo

    ! pick the classical gradient for the nuclei
    traj%grad_ad(:,:)=real(traj%Gmatrix_ssad(traj%state_diag,traj%state_diag,:,:))
    if (printlevel>3) then
      write(u_log,*) ''
      write(u_log,*) 'Gradient of diagonal state',traj%state_diag,'picked.'
      write(u_log,*) ''
      ! TODO: print gradient for printlevel>4
    endif

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
        traj%phases_s=traj%phases_old_s
        do istate=1,ctrl%nstates
          if (real(traj%overlaps_ss(istate,istate))<0.d0) then
            traj%phases_s(istate)=traj%phases_s(istate)*(-1.d0)
          endif
        enddo
      ! from scalar products of old and new NAC vectors
      else

        scalarProd=dcmplx(0.d0,0.d0)

        if (ctrl%calc_nacdr>=0) then
          do istate=1,ctrl%nstates
            do jstate=1,ctrl%nstates
              scalarProd(istate,jstate)=phase_from_NAC(ctrl%natom, &
              &traj%nacdr_ssad(istate,jstate,:,:),traj%nacdr_old_ssad(istate,jstate,:,:) )
            enddo
          enddo
        endif

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
!         traj%phases_s=scalarProd(:,1)
      endif
    endif

    ! Patch phases for Hamiltonian, DM matrix ,NACs, Overlap
    ! Bra
    do istate=1,ctrl%nstates
      traj%H_MCH_ss(istate,:)=traj%H_MCH_ss(istate,:)*traj%phases_s(istate)
      traj%DM_ssd(istate,:,:)=traj%DM_ssd(istate,:,:)*traj%phases_s(istate)
      traj%DM_print_ssd(istate,:,:)=traj%DM_print_ssd(istate,:,:)*traj%phases_s(istate)
      if (ctrl%calc_nacdt==1) then
        traj%NACdt_ss(istate,:)=traj%NACdt_ss(istate,:)*traj%phases_s(istate)
      endif
      if (ctrl%calc_nacdr>=0) then      ! calc_nacdr=0 computes all nacs
        traj%NACdr_ssad(istate,:,:,:)=traj%NACdr_ssad(istate,:,:,:)*real(traj%phases_s(istate))
      endif
      if (ctrl%calc_overlap==1) then
        traj%overlaps_ss(istate,:)=traj%overlaps_ss(istate,:)*traj%phases_old_s(istate)
      endif
    enddo
    ! Ket
    do istate=1,ctrl%nstates
      traj%H_MCH_ss(:,istate)=traj%H_MCH_ss(:,istate)*traj%phases_s(istate)
      traj%DM_ssd(:,istate,:)=traj%DM_ssd(:,istate,:)*traj%phases_s(istate)
      traj%DM_print_ssd(:,istate,:)=traj%DM_print_ssd(:,istate,:)*traj%phases_s(istate)
      if (ctrl%calc_nacdt==1) then
        traj%NACdt_ss(:,istate)=traj%NACdt_ss(:,istate)*traj%phases_s(istate)
      endif
      if (ctrl%calc_nacdr>=0) then
        traj%NACdr_ssad(:,istate,:,:)=traj%NACdr_ssad(:,istate,:,:)*real(traj%phases_s(istate))
      endif
      if (ctrl%calc_overlap==1) then
        traj%overlaps_ss(:,istate)=traj%overlaps_ss(:,istate)*traj%phases_s(istate)
      endif
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

endmodule
