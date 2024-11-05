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

!> # Module INPUT
!>
!> \author Sebastian Mai
!> \date 27.02.2015, modified 27.02.2017 by Philipp Marquetand
!>
!>                   modified 11.13.2019 by Yinan Shu
!>                       add keyword "method" 
!>                       put "tmax" into global variable
!>                       modified keyword "surf", option ehrenfest is deleted
!>                       modified keyword "decoherence", add option decay of mixng
!>                       add keyword "switching_procedure"
!>                       add keyword "integrator"
!>                       add keyword "convthre"
!>                       add keyword "dtmin"
!>                       add keyword "nac_projection"
!>                       
!> This module provides the central input parsing routine and some 
!> auxilliary routines for input parsing
!> 
!> 
module input
 contains

! =================================================================== !

  !> initialises traj and ctrl by reading the input, geom, veloc, laser and coeff files
  !> or alternatively by calling the restart reader
  !> excluded from initialitzation:
  !> - coefficients will be set after the first QM calculation if generated automatically
  !> - all matrices are only initialized after the initial QM calculation (or the restart)
#ifdef __PYSHARC__
  subroutine read_input(nchars,input_name, traj, ctrl)
#else
  subroutine read_input(traj, ctrl)
#endif
  use definitions
  use input_list
  use misc
  use output
  use restart
  use string
  implicit none
#ifdef __PYSHARC__
  integer :: nchars
  character(len=nchars) :: input_name
#endif
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl
  character*255 :: filename
  character*8000 :: geomfilename, line, rattlefilename
  character*8000, allocatable :: values(:)
  integer :: narg, io, nlines, selg, selt
  integer :: i,j,k,n
  integer :: min_order, max_order
  integer :: imult,ims
  real*8 :: a,b,tmax2
  character*24 :: ctime, date
  integer :: idate,time
  character*8000 :: string1
  character*8000, allocatable :: string2(:)
  logical :: selectdirectly_bool

  
#ifndef __PYSHARC__
  ! get the input filename from the command line argument

  ! get number of command line arguments
  narg=iargc()
  ! input filename must be present as argument
  if (narg==0) then
    write(0,*) 'Usage: sharc <inputfile>'
    stop 1
  endif
  call getarg(1,filename)
#endif

  ! =====================================================
#ifdef __PYSHARC__
  
  if ( (trim(input_name)=='-v').or.(trim(input_name)=='--version').or.(trim(input_name)=='--info') ) then
    call write_license(0,version)
    stop
  endif
#else
  if ( (trim(filename)=='-v').or.(trim(filename)=='--version').or.(trim(filename)=='--info') ) then
    call write_license(0,version)
    stop
  endif
#endif

  ! =====================================================
  ! open the input file
#ifdef __PYSHARC__
  open(u_i_input,file=input_name, status='old', action='read', iostat=io)
#else
  open(u_i_input,file=filename, status='old', action='read', iostat=io)
#endif
  if (io/=0) then
    write(0,*) 'Could not find input file "',trim(filename),'"!'
    stop 1
  endif

  call read_input_list_from_file(u_i_input)
  
  ! input file is now in input list
  close(u_i_input)
  
  ! =====================================================

    ! default is no restart
    ctrl%restart=.false.
    ctrl%restart_rerun_last_qm_step=.false.
    ! look for restart keyword
    line=get_value_from_key('restart',io)
    if (io==0) then
      ctrl%restart=.true.
    endif
    ! look for norestart keyword
    line=get_value_from_key('norestart',io)
    if (io==0) then
      ctrl%restart=.false.
    endif
    ! if both keywords are present, norestart will take precendence

  ! =====================================================

  ! open the four main output files
  ! log file, list file, data file and geometry file
    if (ctrl%restart) then
    ! append to old files if restart
      open(unit=u_log, file='output.log', status='old', position='append', iostat=io)
      if (io/=0) then
        write(0,*) 'Restart: Could not find output file "output.log"'
        stop 
      endif
      open(unit=u_lis, file='output.lis', status='old', position='append', iostat=io)
      if (io/=0) then
        write(0,*) 'Restart: Could not find output file "output.lis"'
        stop 
      endif
      open(unit=u_dat, file='output.dat', status='old', position='append', iostat=io)
      if (io/=0) then
        write(0,*) 'Restart: Could not find output file "output.dat"'
        stop 
      endif
      open(unit=u_geo, file='output.xyz', status='old', position='append', iostat=io)
      if (io/=0) then
        write(0,*) 'Restart: Could not find output file "output.xyz"'
        stop 
      endif
    else
    ! make new files if no restart
      open(unit=u_log, file='output.log', status='replace', position='rewind', iostat=io)
      if (io/=0) then
        write(0,*) 'Could not create output file "output.log"'
        stop 
      endif
      open(unit=u_lis, file='output.lis', status='replace', position='rewind', iostat=io)
      if (io/=0) then
        write(0,*) 'Could not create output file "output.lis"'
        stop 
      endif
      open(unit=u_dat, file='output.dat', status='replace', position='rewind', iostat=io)
      if (io/=0) then
        write(0,*) 'Could not create output file "output.dat"'
        stop 
      endif
      open(unit=u_geo, file='output.xyz', status='replace', position='rewind', iostat=io)
      if (io/=0) then
        write(0,*) 'Could not create output file "output.xyz"'
        stop 
      endif
    endif

  ! =====================================================

  call write_logheader(u_log,version)

  ! =====================================================

    ! determine print level
    line=get_value_from_key('printlevel',io)
    if (io==0) then
      read(line,*) printlevel
    else
      ! printlevel=2 is default
      printlevel=2
    endif

    write(u_log,*) 'Print level: ',printlevel
    write(u_log,*)

  ! =====================================================

    ! print complete input (without comments and blank lines)
    if (printlevel>0) then
      write(u_log,*) '============================================================='
      write(u_log,*) '                         Input File'
      write(u_log,*) '============================================================='
      ! get_ncurr() gives the length of the input list
      n=get_ncurr()
      do i=1,n
        write(u_log,*) trim(get_key(i,io)),' ',trim(get_value(i,io))
      enddo
      write(u_log,*)
    endif

  ! =====================================================

    if (printlevel>0) then
      write(u_log,*) '============================================================='
      write(u_log,*) '                            Restart'
      write(u_log,*) '============================================================='
      if (ctrl%restart) then
        write(u_log,*) 'RESTART requested. Reading the dump file...'
      else
        write(u_log,*) 'NO RESTART requested. Setting up the initial data from input files...'
      endif
      write(u_log,*)
    endif
    call flush(u_log)

    if (ctrl%restart) then
!       close(u_i_input)
      call read_restart(u_resc,u_rest,ctrl,traj)

      if (ctrl%integrator==2) then
        ! if explicit laser field is used, the simulation time cannot be changed
        ! otherwise, the simulation time is read from the input file instead of
        ! the restart file
        if (ctrl%laser/=2) then
          line=get_value_from_key('nsteps',io)
          if (io==0) then
            read(line,*) ctrl%nsteps
          else
            line=get_value_from_key('tmax',io)
            if (io==0) then
              read(line,*) ctrl%tmax
              ctrl%nsteps=int(ctrl%tmax/ctrl%dtstep/au2fs+0.01d0)
            endif
          endif
          if (printlevel>0) then
            write(u_log,*) '============================================================='
            write(u_log,*) '                       Simulation Time'
            write(u_log,*) '============================================================='
            write(u_log,*) 'Using Fixed stepsize Velocity-Verlet integrator'
            write(u_log,'(a,1x,i6,1x,a,1x,f6.3,1x,a)') 'Found nsteps=',ctrl%nsteps,'and stepsize=',ctrl%dtstep*au2fs,'fs.'
            if (printlevel>1) then
              write(u_log,'(a,1x,f9.3,1x,a)') 'This makes a total simulation time of ',ctrl%dtstep*ctrl%nsteps*au2fs,'fs.'
              write(u_log,'(a,1x,f7.4,1x,a)') 'The electronic wavefunction will be propagated using a ',&
              &ctrl%dtstep/ctrl%nsubsteps*au2fs, 'fs step.'
            endif
            write(u_log,*)
          endif
        endif
      elseif (ctrl%integrator==0 .or. ctrl%integrator==1) then
        line=get_value_from_key('tmax',io)
        if (io==0) then
          read(line,*) ctrl%tmax
        endif
        if (printlevel>0) then
          write(u_log,*) '============================================================='
          write(u_log,*) '                       Simulation Time'
          write(u_log,*) '============================================================='
          write(u_log,'(a,1x,f9.3,1x,a)') 'Found trial stepsize: ',ctrl%dtstep*au2fs,'fs.'
          write(u_log,'(a,1x,f9.3,1x,a)') 'Total simulation time: ',ctrl%tmax,'fs.'
          if (ctrl%integrator .eq. 0) then
            write(u_log,'(a,1x,f20.10,1x,a)') 'Using Bulirsch-Stoer integrator, threshold:',ctrl%convthre, 'a.u.'
          else if (ctrl%integrator .eq. 1) then
            write(u_log,'(a,1x,f20.10,1x,a)') 'Using adaptive Velocity-Verlet integrator, threshold:', ctrl%convthre*au2eV, 'eV'
          endif
          write(u_log,*)
        endif
      endif
      line=get_value_from_key('restart_rerun_last_qm_step',io)
      if (io==0) then
        ctrl%restart_rerun_last_qm_step=.true.
        write(u_log,'(a)') 'Assuming that restart/ directory corresponds to upcoming time step)'
        write(u_log,'(a)') '(e.g., after a crash or job kill).'
      else
        line=get_value_from_key('restart_goto_new_qm_step',io)
        if (io==0) then
          ctrl%restart_rerun_last_qm_step=.false.
          write(u_log,'(a)') 'Assuming that restart/ directory corresponds to time step in restart files'
          write(u_log,'(a)') '(e.g., after using STOP file, killafter mechanism, or reaching time step limit).'
        else
          write(0,*) 'Please add "restart_rerun_last_qm_step" or "restart_goto_new_qm_step" keyword!'
          stop
        endif
      endif

      ! convert tmax to atomic unit 
      ctrl%tmax=ctrl%tmax/au2fs

      ! rest of this routine is skipped
      return
    else
      traj%step=0
    endif

    idate=time()
    date=ctime(idate)
    if (printlevel>0) then
      write(u_log,'(A)')      '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<============================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(u_log,'(A)')      '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   Initializing Dynamics    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(u_log,'(A)')      '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<============================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(u_log,'(68X,A)')      trim(date)
      write(u_log,*)
    endif

  ! =====================================================
  ! Here starts parsing of the input file step by step
  ! Parsing is not line by line, but rather option by option
  ! =====================================================

  ! next find compatibility mode
  line=get_value_from_key('compatibility',io)
  if (io==0) then
    call split(line,' ',values,n)
    read(values(1),*) ctrl%compat_mode
    deallocate(values)
    write(u_log,*)
    write(u_log,'(A)')      '  #############################################'
    write(u_log,'(A,I3,A)') '  # Compatibility Mode active! Value is = ',ctrl%compat_mode,' #'
    write(u_log,'(A)')      '  #############################################'
    write(u_log,*)
  else
    ! no keyword => no special actions
    ctrl%compat_mode=0
  endif



  ! =====================================================

#ifdef __PYSHARC__
  ctrl%output_format = 1
#else
  ctrl%output_format = 0
#endif
  line=get_value_from_key('output_format', io)
  if (io==0) then
    select case (trim(line))
      case ('ascii') 
        ctrl%output_format=0
      case ('netcdf') 
#ifdef __PYSHARC__
        ctrl%output_format=1
#else
        write(u_log, '(A)') 'Error: Cannot write NetCDF format. Rebuild pysharc with NetCDF support.'
        stop 1
#endif
      case default
        write(0,*) 'Unknown keyword ',trim(line),' to "output_format"!'
        stop 1
    endselect
  endif
  if (printlevel>1) then
    select case (ctrl%output_format)
      case (0)
        write(u_log, '(A)') 'Saving output data in ASCII format (output.dat)'
        write(u_log, '(A)') 'Use data_extractor.x'
      case (1)
        write(u_log, '(A)') 'Saving output data in NetCDF format (output.dat [header] + output.dat.nc)'
        write(u_log, '(A)') 'Use data_extractor_NetCDF.x'
    endselect
    write(u_log, *) 
  endif

  ! =====================================================


  ! next find the number of states

    ! look up nstates keyword
    line=get_value_from_key('nstates',io)
    if (io==0) then
      ! value needs to be split into values (each one is a string)
      call split(line,' ',values,n)
      ! n is number of multiplicities = maxmult
      ctrl%maxmult=n
      allocate(ctrl%nstates_m(n))
      ctrl%nstates=0
      ! read number of states per multiplicity
      do i=1,n
        read(values(i),*) ctrl%nstates_m(i)
        ctrl%nstates=ctrl%nstates+ctrl%nstates_m(i)*i
      enddo
      ! values is not needed anymore
      deallocate(values)
    else
      ! no nstates keywords => adiabatic dynamics, 1 singlet state
      ctrl%maxmult=1
      allocate(ctrl%nstates_m(1))
      ctrl%nstates_m(1)=1
      ctrl%nstates=1
    endif

    if (printlevel>0) then
      write(u_log,*) '============================================================='
      write(u_log,*) '                 Number of States and Atoms'
      write(u_log,*) '============================================================='
      if (printlevel>1) then
        if (io==0) then
          write(u_log,*) 'Keyword NSTATES found.'
          do i=1,ctrl%maxmult
            if (i<=8) then
              write(u_log,'(i4,1x,a)') ctrl%nstates_m(i), trim(multnames(i))
            else
              write(u_log,'(i4,1x,i3,a)') ctrl%nstates_m(i), i,'-tet'
            endif
          enddo
          write(u_log,'(a,1x,i4)') 'Total number of states: ',ctrl%nstates
        else
          write(u_log,*) 'No keyword NSTATES found.'
          write(u_log,*) 'Doing dynamics on a single singlet surface.'
        endif
      endif
      write(u_log,*)
    endif

    ! look up actstates keyword
    allocate( ctrl%actstates_s(ctrl%nstates))
    line=get_value_from_key('actstates',io)
    if (io==0) then
      call split(line,' ',values,n)
      if (n<ctrl%maxmult) then
        write(0,*) 'The length of the actstates_m and nstates_m lists must match!'
        stop 1
      endif
      k=0
      ctrl%actstates_s=.false.
      do imult=1,ctrl%maxmult
        read(values(imult),*) j
        if (j>ctrl%nstates_m(imult)) then
          write(0,*) 'actstates cannot be larger than nstates for any multiplicity!'
          stop 1
        endif
        do ims=1,imult
          ctrl%actstates_s(k+1:k+j)=.true.
          k=k+ctrl%nstates_m(imult)
        enddo
      enddo
      deallocate(values)
      if ((printlevel>1)) then!.and.(.not.all(ctrl%actstates_s) )) then
        write(u_log,*) 'Dynamics is constrained to a subset of active states.'
        write(u_log,*) 'States (T=active, F=frozen):',ctrl%actstates_s
        write(u_log,*)
      endif
    else
      ctrl%actstates_s=.true.
    endif

  ! =====================================================

    ! look up geom filename
    line=get_value_from_key('geomfile',io)
    if (io==0) then
      ! extract string within quotes (can contain spaces)
      call get_quoted(line,geomfilename)
      filename=trim(geomfilename)
    else
      ! default geom filename
      filename='geom'
    endif
  ! open the geometry file
    open(u_i_geom,file=filename, status='old', action='read', iostat=io)
    if (io/=0) then
      write(0,*) 'Could not find geometry file "',trim(filename),'"!'
      stop 1
    endif

  ! =====================================================

  ! next find the number of atoms from file geom (in COLUMBUS geom file format)
  ! count number of consecutive lines with 6 entries

    nlines=0
    do 
      read(u_i_geom,'(A)',iostat=io) line
      if (io/=0) exit
      call split(line,' ', values,n)
      deallocate(values)
      if (n/=6) exit
      nlines=nlines+1
    enddo
    ctrl%natom=nlines
    if (ctrl%natom==0) then
      write(0,*) 'No atoms found in ',geomfilename
      stop 1
    endif

    if (printlevel>1) then
      write(u_log,'(3a,1x,i4)') 'Total number of atoms from "',trim(geomfilename),'":',ctrl%natom
      write(u_log,*)
    endif
    call flush(u_log)

  ! =====================================================

  ! find number of properties

    line=get_value_from_key('n_property1d',io)
    if (io==0) then
      read(line,*) ctrl%n_property1d
    else
      ctrl%n_property1d=1
    endif
    line=get_value_from_key('n_property2d',io)
    if (io==0) then
      read(line,*) ctrl%n_property2d
    else
      ctrl%n_property2d=1
    endif

  ! =====================================================

  ! then call the allocator for the trajectory
    call allocate_traj(traj,ctrl)

    if (printlevel>0) then
      write(u_log,*) '============================================================='
      write(u_log,*) '                          Allocation'
      write(u_log,*) '============================================================='
      if (printlevel>1) then
        write(u_log,'(a,1x,i4,1x,a,1x,i4,1x,a)') 'Allocation with nstates=',ctrl%nstates,'and natom=',ctrl%natom,'successful.'
!         n=sizeof(traj)
!         write(u_log,'(a,1x,i10,1x,a)') 'Using',n,'bytes for trajectory data.'
      endif
      write(u_log,*)
    endif
    call flush(u_log)

  ! =====================================================

  ! setting up simulation time and time step

    ! total simulation time
    line=get_value_from_key('tmax',io)
    if (io==0) then
      read(line,*) ctrl%tmax
    else
      ctrl%tmax=5.0
    endif

    ! stepsize is dt, default is 1 fs
    line=get_value_from_key('stepsize',io)
    if (io==0) then
      read(line,*) ctrl%dtstep
    else
      ctrl%dtstep=1.d0
    endif

    ! minimum timestep allowed in adapative
    line=get_value_from_key('stepsize_min',io)
    if (io==0) then
      read(line,*) ctrl%dtstep_min
    else
      ctrl%dtstep_min=ctrl%dtstep/16
    endif

    ! maximum timestep allowed in adapative
    line=get_value_from_key('stepsize_max',io)
    if (io==0) then
      read(line,*) ctrl%dtstep_max
    else
      ctrl%dtstep_max=ctrl%dtstep*2
    endif
 
    ! Alternative step for minimum timestep allowed in adapative
    line=get_value_from_key('stepsize_min_exp',io)
    if (io==0) then
      read(line,*) min_order
    else
      min_order=-4
    endif
    ctrl%dtstep_min=ctrl%dtstep*2**(min_order)

    ! Alternative step for maximum timestep allowed in adapative
    line=get_value_from_key('stepsize_max_exp',io)
    if (io==0) then
      read(line,*) max_order
    else
      max_order=1
    endif
    ctrl%dtstep_max=ctrl%dtstep*2**(max_order)


    ! nsteps is computed by simulation time/stepsize
    ! nsteps is only used in fixed stepsize Velocity-Verlet integrator
    ctrl%nsteps=int(ctrl%tmax/ctrl%dtstep)

    ! number of substeps for electronic interpolation
    line=get_value_from_key('nsubsteps',io)
    if (io==0) then
      read(line,*) ctrl%nsubsteps
    else
      ctrl%nsubsteps=25
    endif

    ! integrator
    line=get_value_from_key('integrator',io)
    if (io==0) then
      select case (trim(line))
        case ('bsh')
          ctrl%integrator=0
        case ('avv')
          ctrl%integrator=1
        case ('fvv')
          ctrl%integrator=2
        case default
          write(0,*) 'Unknown keyword "',trim(line),'" to "integrator"!'
          stop 1
      endselect
    else
      ctrl%integrator=2
    endif


    if (printlevel>0) then
      if (ctrl%integrator .eq. 2) then
        write(u_log,*) '============================================================='
        write(u_log,*) '                       Simulation Time'
        write(u_log,*) '============================================================='
        write(u_log,'(a)') 'Using fixed stepsize Velocity-Verlet integrator'
        write(u_log,'(a,1x,i6,1x,a,1x,f9.3,1x,a)') 'Found nsteps=',ctrl%nsteps,'and stepsize=',ctrl%dtstep,'fs.'
        if (printlevel>1) then
          write(u_log,'(a,1x,f9.3,1x,a)') 'This makes a total simulation time of ',ctrl%dtstep*ctrl%nsteps,'fs.'
          write(u_log,'(a,1x,f7.4,1x,a)') 'The electronic wavefunction will be propagated using a ',&
          &ctrl%dtstep/ctrl%nsubsteps, 'fs step.'
        endif
        write(u_log,*)
      elseif (ctrl%integrator==0 .or. ctrl%integrator==1) then
        write(u_log,*) '============================================================='
        write(u_log,*) '                       Simulation Time'
        write(u_log,*) '============================================================='
        write(u_log,'(a,1x,f9.3,1x,a)') 'Found trial stepsize: ',ctrl%dtstep,'fs.'
        write(u_log,'(a,1x,f9.3,1x,a)') 'Total simulation time: ',ctrl%tmax,'fs.'
        if (ctrl%integrator .eq. 0) then
          write(u_log,'(a)') 'Using Bulirsch-Stoer integrator'
        elseif (ctrl%integrator .eq. 1) then
          write(u_log,'(a)') 'Using adaptive Velocity-Verlet integrator'
        endif
      endif
      write(u_log,*)
    endif

    ! kill after n timesteps in gs ============================================

    line=get_value_from_key('killafter',io)
    if (io==0) then
      read(line,*) tmax2
      if (tmax2<=0) then
        ctrl%killafter=-1
        traj%steps_in_gs=-123
      else
        ctrl%killafter=int(anint(tmax2/ctrl%dtstep))
        traj%steps_in_gs=0
      endif
    else
      ctrl%killafter=-1
      traj%steps_in_gs=-123
    endif

    if (ctrl%killafter>=1) then
      if (printlevel>1) then
        write(u_log,'(A)') 'Trajectories will be killed if staying in the lowest state'
        write(u_log,'(A,1X,I6,1X,A)') 'for longer than',ctrl%killafter,'steps.'
        write(u_log,*)
      endif
    endif

  ! ionization keyword

  line=get_value_from_key('ionization',io)
  if (io==0) then
    ctrl%ionization=1
    line=get_value_from_key('ionization_step',io)
    if (io==0) then
      read(line,*) ctrl%ionization
    endif
    if (printlevel>1) then
      write(u_log,'(A,1X,I6,1X,A)') 'Calculating ionizations every',ctrl%ionization,'steps.'
      write(u_log,*)
    endif
  else
    ctrl%ionization=-1
  endif
  line=get_value_from_key('noionization',io)
  if (io==0) then
    ctrl%ionization=-1
  endif

  ! theodore keyword

  line=get_value_from_key('theodore',io)
  if (io==0) then
    ctrl%theodore=1
    line=get_value_from_key('theodore_step',io)
    if (io==0) then
      read(line,*) ctrl%theodore
    endif
    if (printlevel>1) then
      write(u_log,'(A,1X,I6,1X,A)') 'Calculating TheoDORE properties every',ctrl%theodore,'steps.'
      write(u_log,*)
    endif
  else
    ctrl%theodore=-1
  endif
  line=get_value_from_key('notheodore',io)
  if (io==0) then
    ctrl%theodore=-1
  endif

  ! =====================================================

  ! setting up simulation method
  ! either trajectory surface hopping or self-consistent potential
  ! notice self-consistent potential is same as semi-classical Ehrenfest 

   ! Method type
    line=get_value_from_key('method',io)
    if (io==0) then
      select case (trim(line))
        case ('tsh')
          ctrl%method=0
        case ('scp')
          ctrl%method=1
        case ('ehrenfest')
          ctrl%method=1
        case ('csdm')
          ctrl%method=1
        case default
          write(0,*) 'Unknown keyword "',trim(line),'" to "method"!'
          stop 1
      endselect
    else ! set the default as surface hopping
      ctrl%method=0 
    endif

  ! =====================================================

  ! other keywords, including consistency checks

    ! Dynamics type
    line=get_value_from_key('surf',io)
    if (io==0) then
      select case (trim(line))
        case ('sharc') 
          ctrl%surf=0
        case ('diagonal') 
          ctrl%surf=0
        case ('fish')
          ctrl%surf=1
        case ('mch')
          ctrl%surf=1
        case ('ehrenfest')
          write(0,*) 'Ehrenfest dynamics is not implemented in this SHARC version.'
          stop 1
        case default
          write(0,*) 'Unknown keyword "',trim(line),'" to "surf"!'
          stop 1
      endselect
    else
      ctrl%surf=0
    endif

    ! coupling quantity
    line=get_value_from_key('coupling',io)
    if (io==0) then
      select case (trim(line))
        case ('ddt') 
          ctrl%coupling=0
        case ('nacdt') 
          ctrl%coupling=0
        case ('ddr')
          ctrl%coupling=1
        case ('nacdr')
          ctrl%coupling=1
        case ('overlap')
          ctrl%coupling=2
        case ('ktdc') ! kappa tdc
          ctrl%coupling=3
        case default
          write(0,*) 'Unknown keyword ',trim(line),' to "coupling"!'
          stop 1
      endselect
    else
      ctrl%coupling=2
    endif

    ! turn program off if adaptive is used for overlap
    if (ctrl%integrator==0 .or. ctrl%integrator==1) then
      if (ctrl%coupling==2) then 
        write(0,*) 'Adaptive integrator is not yet working for coupling overlap'
        stop 1
      endif
    endif

    ! method to compute kappa TDC
    if (ctrl%coupling==3) then
      line=get_value_from_key('ktdc_method',io)
      if (io==0) then
        select case (trim(line))
          case ('gradient')
            ctrl%ktdc_method=0
          case ('energy')
            ctrl%ktdc_method=1
          case default
            write(0,*) 'Unknown keyword ',trim(line),' to "ktdc_method"!'
            stop 1
        endselect
      else  ! set default
        if (ctrl%method==0) then ! for TSH, use energy based
          ctrl%ktdc_method=1
        else if (ctrl%method==1) then ! for SCP, use gradient based 
          ctrl%ktdc_method=0
        endif 
      endif
    endif

    ! method to compute kmatrix
    line=get_value_from_key('kmatrix_method',io)
    if (io==0) then
      select case (trim(line))
        case ('gradient')
          ctrl%kmatrix_method=0
        case ('energy')
          ctrl%kmatrix_method=1
        case default
          write(0,*) 'Unknown keyword ',trim(line),' to "kmatrix_method"!'
          stop 1
      endselect
    else ! set the default
      if (ctrl%coupling==1 .or. ctrl%coupling==2) then 
        ctrl%kmatrix_method=0
      else if (ctrl%coupling==3) then 
        ctrl%kmatrix_method=1
      endif 
    endif

    ! method to compute kmatrix, another keyword
    line=get_value_from_key('tdm_method',io)
    if (io==0) then
      select case (trim(line))
        case ('gradient')
          ctrl%kmatrix_method=0
        case ('energy')
          ctrl%kmatrix_method=1
        case default
          write(0,*) 'Unknown keyword ',trim(line),' to "kmatrix_method"!'
          stop 1
      endselect
    else ! set the default
      if (ctrl%coupling==1 .or. ctrl%coupling==2) then
        ctrl%kmatrix_method=0
      else if (ctrl%coupling==3) then
        ctrl%kmatrix_method=1
      endif
    endif

    ! electronic eom method for TSH and SCP
    line=get_value_from_key('eeom',io)
    if (io==0) then
      select case (trim(line))
        case ('ci')
          ctrl%eeom=0
        case ('li')
          ctrl%eeom=1
        case ('ld')
          ctrl%eeom=2
        case ('npi')
          ctrl%eeom=3
        case default
          write(0,*) 'Unknown keyword ',trim(line),' to "eeom"!'
          stop 1
      endselect
    else ! set the default nuclear propagators
      if (ctrl%coupling==0) then
        ctrl%eeom=0
      elseif (ctrl%coupling==1 .or. ctrl%coupling==3) then
        ctrl%eeom=1
      elseif (ctrl%method==0 .and. ctrl%coupling==2) then
        ctrl%eeom=2
      elseif (ctrl%method==1 .and. ctrl%coupling==2) then
        ctrl%eeom=3
      endif
    endif

    ! nuclear eom method for scp
    line=get_value_from_key('neom',io)
    if (io==0) then
      select case (trim(line))
        case ('ddr')
          ctrl%neom=0
        case ('nacdr')
          ctrl%neom=0
        case ('gdiff')
          ctrl%neom=1
        case default
          write(0,*) 'Unknown keyword ',trim(line),' to "neom"!'
          stop 1
      endselect
    else ! set the default nuclear propagators
      if (ctrl%coupling==0) then
        ctrl%neom=0
      elseif (ctrl%coupling==1) then
        ctrl%neom=0
      elseif (ctrl%coupling==2) then
        ctrl%neom=1
      elseif (ctrl%coupling==3) then
        ctrl%neom=1
      endif
    endif

    line=get_value_from_key('neom_rep',io)
    if (io==0) then
      select case (trim(line))
        case ('diag')
          ctrl%neom_rep=0
        case ('mch')
          ctrl%neom_rep=1
        case default
          write(0,*) 'Unknown keyword ',trim(line),' to "neom"!'
          stop 1
      endselect
    else ! set the default 
      ctrl%neom_rep=0
    endif

    ! default is do projection for system with more than 4 atoms
    if (ctrl%method==1 .and. ctrl%natom.ge.4) then
      ctrl%nac_projection=1
    else
      ctrl%nac_projection=0
    endif
    line=get_value_from_key('nac_projection',io)
    if (io==0) then
      ctrl%nac_projection=1
    endif
    line=get_value_from_key('nonac_projection',io)
    if (io==0) then
      ctrl%nac_projection=0
    endif

    ! method to maintain ZPE
    ctrl%zpe_correction=0
    line=get_value_from_key('zpe_correction',io)
    if (io==0) then
      select case (trim(line))
        case ('pumping')
          ctrl%zpe_correction=1
        case ('lp')
          ctrl%zpe_correction=2
        case default
          write(0,*) 'Unknown keyword ',trim(line),' to "zpe_correction"!'
          stop 1
      endselect
    endif

    ! LP-ZPE correction scheme input
 
    ! LP-ZPE correction scheme
    line=get_value_from_key('lpzpe_scheme',io)
    if (io==0) then
      read(line,*) ctrl%lpzpe_scheme
    else
      ! default: using original scheme
      ctrl%lpzpe_scheme=0
    endif

    ! number of AH bonds
    line=get_value_from_key('number_ah',io)
    if (io==0) then
      read(line,*) ctrl%lpzpe_nah
    endif

    ! number of BC bonds
    line=get_value_from_key('number_bc',io)
    if (io==0) then
      read(line,*) ctrl%lpzpe_nbc
    endif

    ! list of AH bonds
    line=get_value_from_key('ah_list',io)
    if (io==0) then
      ! value needs to be split into values (each one is a string)
      call split(line,' ',values,n)
      ! n is twice the number of AH bonds
      ctrl%lpzpe_nah=n/2
      allocate(ctrl%lpzpe_ah(ctrl%lpzpe_nah,2)) 
      ! read AH bonds 
      do i=1,ctrl%lpzpe_nah
        do j=1,2
          k=(i-1)*2+j
          read(values(k),*) ctrl%lpzpe_ah(i,j)
        enddo
      enddo
      ! values is not needed anymore
      deallocate(values)
    endif

    ! list of BC bonds
    line=get_value_from_key('bc_list',io)
    if (io==0) then
      ! value needs to be split into values (each one is a string)
      call split(line,' ',values,n)
      ! n is twice the number of AH bonds
      ctrl%lpzpe_nbc=n/2
      allocate(ctrl%lpzpe_bc(ctrl%lpzpe_nbc,2))
      ! read AH bonds 
      do i=1,ctrl%lpzpe_nbc
        do j=1,2
          k=(i-1)*2+j
          read(values(k),*) ctrl%lpzpe_bc(i,j)
        enddo
      enddo
      ! values is not needed anymore
      deallocate(values)
    endif

    ! zpe of AH bonds
    line=get_value_from_key('ah_zpe',io)
    if (io==0) then
      ! value needs to be split into values (each one is a string)
      call split(line,' ',values,n)
      allocate(ctrl%lpzpe_ke_zpe_ah(n))
      ! read AH bonds 
      do i=1,n
        read(values(i),*) ctrl%lpzpe_ke_zpe_ah(i)
      enddo
      ! values is not needed anymore
      deallocate(values)
    endif

    ! zpe of BC bonds
    line=get_value_from_key('bc_zpe',io)
    if (io==0) then
      ! value needs to be split into values (each one is a string)
      call split(line,' ',values,n)
      allocate(ctrl%lpzpe_ke_zpe_bc(n))
      ! read AH bonds 
      do i=1,n
        read(values(i),*) ctrl%lpzpe_ke_zpe_bc(i)
      enddo
      ! values is not needed anymore
      deallocate(values)
    else 
      allocate(ctrl%lpzpe_ke_zpe_bc(ctrl%lpzpe_nbc))
      do i=1,ctrl%lpzpe_nbc
        ctrl%lpzpe_ke_zpe_bc(i)=0.0
      enddo
    endif

    ! read kinetic energy threshold
    line=get_value_from_key('ke_threshold',io)
    if (io==0) then
      read(line,*) ctrl%ke_threshold
    endif

    ! read time cycle of zpe checking
    line=get_value_from_key('t_cycle',io)
    if (io==0) then
      read(line,*) ctrl%t_cycle
    endif
    ctrl%t_cycle=ctrl%t_cycle/au2fs

    ! read period of time for kinetic energy averaging
    line=get_value_from_key('t_check',io)
    if (io==0) then
      read(line,*) ctrl%t_check
    endif
    ctrl%t_check=ctrl%t_check/au2fs

    ! Done LP-ZPE correction scheme input
   
    ! pointer basis selection
    line=get_value_from_key('pointer_basis',io)
    if (io==0) then
      select case (trim(line))
        case ('diag')
          ctrl%pointer_basis=0
        case ('mch')
          ctrl%pointer_basis=1
        case ('opt')
          ctrl%pointer_basis=2
        case default
          write(0,*) 'Unknown keyword ',trim(line),' to "pointer_basis"!'
          stop 1
      endselect
    else
      ctrl%pointer_basis=0
    endif

    ! maximum iteration of pointer basis optimization
    line=get_value_from_key('pointer_maxiter',io)
    if (io==0) then
      read(line,*) ctrl%pointer_maxiter
    else
      ctrl%pointer_maxiter=1000
    endif

    ! spin-orbit couplings
    ctrl%calc_soc=1
    line=get_value_from_key('nospinorbit',io)
    if (io==0) then
      ctrl%calc_soc=0
    endif
    line=get_value_from_key('spinorbit',io)
    if (io==0) then
      ctrl%calc_soc=1
    endif

    ! request phase corrections from interface
    if (ctrl%coupling.ne.3) then 
      ctrl%calc_phases=1
    elseif (ctrl%coupling==3) then
      ctrl%calc_phases=0
    endif
    line=get_value_from_key('nophases_from_interface',io)
    if (io==0) then
      ctrl%calc_phases=0
    endif
    line=get_value_from_key('phases_from_interface',io)
    if (io==0) then
      ctrl%calc_phases=1
    endif

    ! request phase corrections from interface at time step zero (only works if something is in savedir)
    ctrl%track_phase_at_zero=0
    line=get_value_from_key('phases_at_zero',io)
    if (io==0) then
      ctrl%track_phase_at_zero=1
    endif

    ! non-adiabatic couplings for gradients
    line=get_value_from_key('gradcorrect',io)
    if (io==0) then
      select case (trim(line))
        case ('')
          ctrl%gradcorrect=1
        case ('nac')
          ctrl%gradcorrect=1
        case ('ngt')
          ctrl%gradcorrect=1
        case ('kmatrix')
          ctrl%gradcorrect=2
        case ('tdm')
          ctrl%gradcorrect=2
        case ('enac')
          ctrl%gradcorrect=3
        case default
          write(0,*) 'Unknown keyword ',trim(line),' to "gradcorrect"!'
          stop 1
      endselect
    else
      ctrl%gradcorrect=0
    endif
    line=get_value_from_key('nogradcorrect',io)
    if (io==0) then
      ctrl%gradcorrect=0
    endif

    ! adjustment of kinetic energy
    line=get_value_from_key('ekincorrect',io)
    if (io==0) then
      select case (trim(line))
        case ('none')
          ctrl%ekincorrect=0
        case ('parallel_vel')
          ctrl%ekincorrect=1
        case ('parallel_pvel')
          ctrl%ekincorrect=2
        case ('parallel_nac')
          ctrl%ekincorrect=3
        case ('parallel_diff')
          ctrl%ekincorrect=4
        case ('parallel_pnac')
          ctrl%ekincorrect=5
        case ('parallel_pdiff')
          ctrl%ekincorrect=6
        case ('parallel_enac') ! effective NAC
          ctrl%ekincorrect=7
        case ('parallel_penac') ! projected effective NAC
          ctrl%ekincorrect=8
        case default
          write(0,*) 'Unknown keyword ',trim(line),' to "ekincorrect"!'
          stop 1
      endselect
    else
      ctrl%ekincorrect=1
    endif

    ! reflection after frustrated hops
    line=get_value_from_key('reflect_frustrated',io)
    if (io==0) then
      select case (trim(line))
        case ('none') 
          ctrl%reflect_frustrated=0
        case ('parallel_vel') 
          ctrl%reflect_frustrated=1
        case ('parallel_pvel')
          ctrl%reflect_frustrated=2
        case ('parallel_nac') 
          ctrl%reflect_frustrated=3
        case ('parallel_diff') 
          ctrl%reflect_frustrated=4
        case ('parallel_pnac')
          ctrl%reflect_frustrated=5
        case ('parallel_pdiff')
          ctrl%reflect_frustrated=6
        case ('parallel_enac') ! effective NAC
          ctrl%reflect_frustrated=7
        case ('parallel_penac') ! projected effective NAC
          ctrl%reflect_frustrated=8
        case ('delV_vel')
          ctrl%reflect_frustrated=91
        case ('delV_pvel')
          ctrl%reflect_frustrated=92
        case ('delV_nac')
          ctrl%reflect_frustrated=93
        case ('delV_diff')
          ctrl%reflect_frustrated=94
        case ('delV_pnac')
          ctrl%reflect_frustrated=95
        case ('delV_pdiff')
          ctrl%reflect_frustrated=96
        case ('delV_enac')
          ctrl%reflect_frustrated=97
        case ('delV_penac')
          ctrl%reflect_frustrated=98
        case default
          write(0,*) 'Unknown keyword ',trim(line),' to "reflect_frustrated"!'
          stop 1
      endselect
    else
      ctrl%reflect_frustrated=0
    endif

    ! initialize time uncertainty
    ctrl%time_uncertainty=0
  
    ! do time uncertainty or not
    line=get_value_from_key('time_uncertainty',io)
    if (io==0) then
      ctrl%time_uncertainty=1
    endif
    line=get_value_from_key('notime_uncertainty',io)
    if (io==0) then
      ctrl%time_uncertainty=0
    endif

    ! initialize time uncertainty process
    traj%in_time_uncertainty=0

    ! selection of gradients/non-adiabatic couplings
    selg=0
    line=get_value_from_key('grad_select',io)
    if (io==0) then
      selg=1
    endif
    line=get_value_from_key('grad_all',io)
    if (io==0) then
      selg=0
    endif
    line=get_value_from_key('nograd_select',io)
    if (io==0) then
      selg=0
    endif
    selt=0
    line=get_value_from_key('nac_select',io)
    if (io==0) then
      selt=1
    endif
    line=get_value_from_key('nac_all',io)
    if (io==0) then
      selt=0
    endif
    line=get_value_from_key('nonac_select',io)
    if (io==0) then
      selt=0
    endif



    ! ===========================================
    ! process which quantities need to be calculated by the interface
    ! initialize to "do not calculate"
    ctrl%calc_grad=0
    ctrl%calc_nacdt=0
    ctrl%calc_nacdr=-1
    ctrl%calc_overlap=0
    ctrl%calc_second=0
    ctrl%calc_effectivenac=0

    ! for debug purposes, direct set up of keywords 
    line=get_value_from_key('calc_overlap',io)
    if (io==0) then
      ctrl%calc_overlap=1
    endif
    line=get_value_from_key('calc_effectivenac',io)
    if (io==0) then
      ctrl%calc_effectivenac=1
    endif

    ! calculate the active coupling quantity (nacdr, nacdt, overlaps)
    select case (ctrl%coupling)
      case (0)  ! ddt
        ctrl%calc_nacdt=1
      case (1)  ! ddr
        ctrl%calc_nacdr=0
      case (2)  ! overlap
        ctrl%calc_overlap=1
      case (3)  ! ktdc
        ctrl%calc_effectivenac=1
      case default
        write(0,*) 'Internal error 1!'
        stop 1
    endselect

    if (ctrl%method==1) then 
      select case (ctrl%neom)
        case (0)  ! ddr
          ctrl%calc_nacdr=0
        case (1)  ! effective nac
          ctrl%calc_effectivenac=1
        case default
          write(0,*) 'Internal error 1!'
          stop 1
      endselect
    endif 

!     if (ctrl%coupling==1) then   ! having ddr couplings
!       ctrl%gradcorrect=1         ! they can be used for gradient correction
!     endif

    if (ctrl%surf/=0) then   ! not doing SHARC dynamics
      if (ctrl%gradcorrect==1) then
        write(0,*) 'Info: keyword "gradcorrect" has no effect in dynamics on MCH potentials.'
      endif
      ctrl%gradcorrect=0     ! no need for transformed gradients
    endif

    if (ctrl%gradcorrect==1) then ! for gradcorrect
      ctrl%calc_nacdr=0           ! we need nacdr
    endif
    if (ctrl%ekincorrect==3 .or. ctrl%ekincorrect==5) then ! for kinetic energy correction parallel to nac or projected nac, we need nacdr
      ctrl%calc_nacdr=0
      ctrl%gradcorrect=1        ! NACs must be transformed
    endif
    if (ctrl%ekincorrect==7 .or. ctrl%ekincorrect==8) then ! for kinetic energy correction parallel to effective nac or projected effective nac, we need effective nac
      ctrl%calc_effectivenac=1
    endif 
    if (ctrl%reflect_frustrated==3 .or. ctrl%reflect_frustrated==5 .or. ctrl%reflect_frustrated==93 .or. ctrl%reflect_frustrated==95) then ! for reflection parallel to nac or projected nac, we need nacdr
      ctrl%calc_nacdr=0
      ctrl%gradcorrect=1        ! NACs must be transformed
    endif
    if (ctrl%reflect_frustrated==7 .or. ctrl%reflect_frustrated==8 .or. ctrl%reflect_frustrated==97 .or. ctrl%reflect_frustrated==98) then ! for reflection parallel to effective nac or projected effective nac, we need effective nac
      ctrl%calc_effectivenac=1
    endif
! !     if (ctrl%decoherence==2) then ! for A-FSSH we need nacdr
! !       ctrl%calc_nacdr=0
! !     endif

    if (ctrl%surf==1) then   ! doing MCH dynamics
      ctrl%calc_grad=1       ! always select grad
    endif

    if (selg==1) then           ! select gradients
      if (ctrl%surf==0) then
        ctrl%calc_grad=2        ! for SHARC in second calculation
      else
        ctrl%calc_grad=1        ! else in first calculation
      endif
    endif

    if (selt==1.and.ctrl%calc_nacdr==0) then    ! if selection of non-adiabatic couplings
      ctrl%calc_nacdr=2                         ! calculate them in second calculation
    endif

    if (ctrl%calc_grad==2.or.ctrl%calc_nacdr==2) then
      ctrl%calc_second=1                                ! do a second interface call
    endif

    selectdirectly_bool=.true.
    line=get_value_from_key('select_directly',io)       ! do not do a second interface call
    if (io==0) then
      selectdirectly_bool=.true.
    endif
    line=get_value_from_key('noselect_directly',io)       ! do a second interface call
    if (io==0) then
      selectdirectly_bool=.false.
    endif
    if (selectdirectly_bool) then
      if (ctrl%calc_grad==2) ctrl%calc_grad=1
      if (ctrl%calc_nacdr==2) ctrl%calc_nacdr=1
      ctrl%calc_second=0
    endif

    ! obtain the energy threshold for selection, default 0.5 eV
    line=get_value_from_key('eselect',io)
    if (io==0) then
      read(line,*) ctrl%eselect_grad
      ctrl%eselect_nac=ctrl%eselect_grad
    else
      ctrl%eselect_grad=0.5d0
      ctrl%eselect_nac=0.5d0
    endif

    ! if we have frozen states, setup selection
    if ( (ctrl%calc_nacdr==0).and.(.not.all(ctrl%actstates_s) ) ) then
      ctrl%calc_nacdr=1
      ctrl%eselect_nac=99999.9d0
    endif

    ! if we have frozen states, setup selection
    if ( (ctrl%calc_grad==0).and.(.not.all(ctrl%actstates_s) ) ) then
      ctrl%calc_grad=1
      ctrl%eselect_grad=99999.9d0
    endif




    ! Flags for switching on/off writing data to output.dat
    !
    ctrl%write_soc=1                     !< write SOC to output.dat or not \n 0=no soc, 1=write soc )
!     line=get_value_from_key('write_soc',io)
!     if (io==0) then
!       ctrl%write_soc=1
!     endif
!     line=get_value_from_key('nowrite_soc',io)
!     if (io==0) then
!       ctrl%write_soc=0
!     endif
!     if (ctrl%calc_soc==0 .and. ctrl%write_soc==1) then
!       write(u_log,*) 'Warning! Requested writing SOCs but no SOC calculated. Writing of SOC disabled.'
!       ctrl%write_soc = 0
!     endif

    ctrl%write_grad=0                     !< write gradients:   \n        0=no gradients, 1=write gradients
    line=get_value_from_key('write_grad',io)
    if (io==0) then
      ctrl%write_grad=1
    endif
    line=get_value_from_key('nowrite_grad',io)
    if (io==0) then
      ctrl%write_grad=0
    endif

    ctrl%write_overlap=1                  !< write overlap matrix:   \n        0=no overlap, 1=write overlap
    line=get_value_from_key('write_overlap',io)
    if (io==0) then
      ctrl%write_overlap=1
    endif
    line=get_value_from_key('nowrite_overlap',io)
    if (io==0) then
      ctrl%write_overlap=0
    endif
    if (ctrl%calc_overlap==0 .and. ctrl%write_overlap==1) then
      write(u_log,*) 'Warning! Requested writing overlaps but no overlap calculated. Writing of overlap disabled.'
      ctrl%write_overlap = 0
    endif

    ctrl%write_NACdr=0                   !< write nonadiabatic couplings:   \n        0=no NACs, 1=write NACs
    line=get_value_from_key('write_nacdr',io)
    if (io==0) then
      ctrl%write_NACdr=1
    endif
    line=get_value_from_key('nowrite_nacdr',io)
    if (io==0) then
      ctrl%write_NACdr=0
    endif
    if (ctrl%calc_nacdr==-1 .and. ctrl%write_NACdr==1) then
      write(u_log,*) 'Warning! Requested writing NonAdiabatic Coupling (NACs) but no NACs calculated. Writing of NACs disabled.'
      ctrl%write_NACdr = 0
    endif

    if (ctrl%theodore==1) then
      ctrl%write_property1d=1
    else
      ctrl%write_property1d=0                !< write property vectors:   \n        0=no property, 1=write property
    endif
    line=get_value_from_key('write_property1d',io)
    if (io==0) then
      ctrl%write_property1d=1
    endif
    line=get_value_from_key('nowrite_property1d',io)
    if (io==0) then
      ctrl%write_property1d=0
    endif

    if (ctrl%ionization==1) then
      ctrl%write_property2d=1
    else
      ctrl%write_property2d=0                !< write property vectors:   \n        0=no property, 1=write property
    endif
    ctrl%write_property2d=0                !< write property matrices:   \n        0=no property, 1=write property
    line=get_value_from_key('write_property2d',io)
    if (io==0) then
      ctrl%write_property2d=1
    endif
    line=get_value_from_key('nowrite_property2d',io)
    if (io==0) then
      ctrl%write_property2d=0
    endif


    ! how often output.dat is written
    ctrl%output_steps_limits=0
    ctrl%output_steps_stride=1
    line=get_value_from_key('output_dat_steps',io)
    if (io==0) then
         if (ctrl%integrator/=2) then
             write(0,*) 'Bulirsch-Stoer-Hack and adaptive velocity Verlet integrators not yet compatible with variable output stride!'
             stop 1
         endif 
      call split(line,' ',values,n)
      if (n>=1) then
        read(values(1),*) i
        ctrl%output_steps_stride=max(i,1)
      endif
      if (n>=3) then
        read(values(2),*) i
        ctrl%output_steps_limits(2)=max(i,0)
        ctrl%output_steps_limits(3)=max(i,0)
        read(values(3),*) i
        ctrl%output_steps_stride(2)=max(i,1)
        ctrl%output_steps_stride(3)=max(i,1)
      endif
      if (n>=5) then
        read(values(4),*) i
        ctrl%output_steps_limits(3)=max(i,0)
        read(values(5),*) i
        ctrl%output_steps_stride(3)=max(i,1)
      endif
    endif



    line=get_value_from_key('output_version',io)
    if (io==0) then
      read(line,*) ctrl%output_version
    else
      string1=version
      call split(string1,' ',string2,n)
      read(string2(1),*) ctrl%output_version
      deallocate(string2)
    endif



    if (printlevel>0) then
      write(u_log,*) '============================================================='
      write(u_log,*) '                       Dynamics options'
      write(u_log,*) '============================================================='
      if (printlevel>1) then
        select case (ctrl%method)
          case (0)
            write(u_log,'(a)') 'Doing Trajectory Surface Hopping Dynamics.'
          case (1)
            write(u_log,'(a)') 'Doing semi-classical Ehrenfest Dynamics using self-consistent potentials'
        endselect
        if (ctrl%time_uncertainty==1) then
            write(u_log,'(a)') 'Doing Trajectory Surface Hopping with Time Uncertainty'
        endif
        select case (ctrl%surf)
          case (0)
            write(u_log,'(a)') 'Doing SHARC dynamics (on diagonal surfaces).'
          case (1)
            write(u_log,'(a)') 'Doing dynamics on MCH surfaces.'
        endselect
        select case (ctrl%zpe_correction)
          case (0)
            write(u_log,'(a)') 'No Zero Point Energy correction scheme employed'
          case (1)
            write(u_log,'(a)') 'Using Zero Point Energy pumping method tomaintain ZPE'
          case (2)
            write(u_log,'(a)') 'Using Local Pair Zero Point Energy to maintain ZPE'
            if (ctrl%lpzpe_scheme==0) then 
              write(u_log,'(a)') 'Using original LP-ZPE scheme: skip correction if BC kinetic energy is not enough'
            else if (ctrl%lpzpe_scheme==1) then
              write(u_log,'(a)') 'Using new LP-ZPE scheme: adjust the correction energy based on BC kinetic energy'
            endif
        endselect
        select case (ctrl%coupling)
          case (0)
            write(u_log,'(a)') 'Using TIME DERIVATIVES <a|d/dt|b> as wavefunction coupling.'
          case (1)
            write(u_log,'(a)') 'Using SPATIAL DERIVATIVES <a|d/dR|b> as wavefunction coupling.'
          case (2)
            write(u_log,'(a)') 'Using OVERLAPS as wavefunction coupling.'
          case (3)
            write(u_log,'(a)') 'Using curvature TDC approximation as wavefunction coupling.'
        endselect
        select case (ctrl%eeom)
          case (0)
            write(u_log,'(a)') 'Using constant interpolation of TDC for wavefunction propagation.'
          case (1)
            write(u_log,'(a)') 'Using linear interpolation of TDC for wavefunction propagation.'
          case (2)
            write(u_log,'(a)') 'Doing LOCAL DIABATISATION propagation.'
          case (3)
            write(u_log,'(a)') 'Using norm perserving interpolation of TDC for wavefunction propagation.'
        endselect
        if (ctrl%method==1) then 
          select case (ctrl%neom)   ! this should this only be printed if SCP is used
            case (0)
              write(u_log,'(a)') 'Using NAC for nuclear equation of motion.'
            case (1)
              write(u_log,'(a)') 'Using effective NAC for nuclear equation of motion.'
          endselect
        endif
        if (ctrl%gradcorrect==1) then
          write(u_log,'(a)') 'Including non-adiabatic coupling vectors in the gradient transformation.'
        elseif (ctrl%gradcorrect==2) then
          write(u_log,'(a)') 'Using Kmatrix in the gradient transformation.'
        endif
        select case (ctrl%ekincorrect)
          case (0)
            write(u_log,'(a)') 'No correction to the kinetic energy after surface hop.'
          case (1)
            write(u_log,'(a)') 'Correction to the kinetic energy after surface hop parallel to velocity vector.'
          case (2)
            write(u_log,'(a)') 'Correction to the kinetic energy after surface hop parallel to projected velocity vector.'
          case (3)
            write(u_log,'(a)') 'Correction to the kinetic energy after surface hop parallel to non-adiabatic coupling vector.'
          case (4)
            write(u_log,'(a)') 'Correction to the kinetic energy after surface hop parallel to gradient difference vector.'
          case (5)
            write(u_log,'(a)') 'Correction to the kinetic energy after surface hop parallel to projected non-adiabatic coupling vector.'
          case (6)
            write(u_log,'(a)') 'Correction to the kinetic energy after surface hop parallel to projected gradient difference vector.'
          case (7)
            write(u_log,'(a)') 'Correction to the kinetic energy after surface hop parallel to effective non-adiabatic coupling vector.'
          case (8)
            write(u_log,'(a)') 'Correction to the kinetic energy after surface hop parallel to projected effective non-adiabatic coupling vector.'
        endselect
        write(u_log,*)
        select case (ctrl%calc_grad)
          case (0)
            write(u_log,'(a)') 'Including all gradients in the dynamics.'
          case (1)
            if (ctrl%surf==1) then
              write(u_log,'(a)') 'Including only gradient of classical state.'
            else
              write(u_log,'(a,1x,f8.3,1x,a)') 'Gradients are included for Delta E < ',ctrl%eselect_grad,'eV.'
              write(u_log,'(a)') 'Gradients are calculated in the first QM calculation.'
              if (.not.all(ctrl%actstates_s) ) then
                write(u_log,'(a,1x,f6.3)') 'Including only gradients of the active states in the dynamics.'
              endif
            endif
          case (2)
            write(u_log,'(a,1x,f8.3,1x,a)') 'Gradients are included for Delta E < ',ctrl%eselect_grad,'eV.'
            write(u_log,'(a)') 'Gradients are calculated in the second QM calculation.'
        endselect
        if (ctrl%calc_nacdt==1) then
          write(u_log,'(a)') 'Calculating time derivatives.'
        endif
        if (ctrl%calc_overlap==1) then
          write(u_log,'(a)') 'Calculating wavefunction overlaps.'
        endif
        if (ctrl%calc_phases==1) then
          write(u_log,'(a)') 'Calculating wavefunction phases.'
        endif
        select case (ctrl%calc_nacdr)
          case (0)
            write(u_log,'(a)') 'Including all non-adiabatic coupling vectors in the dynamics.'
          case (1)
            write(u_log,'(a,1x,f8.3,1x,a)') 'Non-adiabatic coupling vectors are included for Delta E < ',ctrl%eselect_nac,'eV.'
            write(u_log,'(a)') 'Non-adiabatic coupling vectors are calculated in the first QM calculation.'
            if (.not.all(ctrl%actstates_s) ) then
              write(u_log,'(a)') 'Including only non-adiabatic coupling vectors of the active states.'
            endif
          case (2)
            write(u_log,'(a,1x,f8.3,1x,a)') 'Non-adiabatic coupling vectors are included for Delta E < ',ctrl%eselect_nac,'eV.'
            write(u_log,'(a)') 'Non-adiabatic coupling vectors are calculated in the second QM calculation.'
        endselect
        if (ctrl%calc_second==1) then
          write(u_log,'(a)') 'Doing two QM calculations per step for selection.'
        endif
        if (ctrl%calc_soc==1) then
          write(u_log,'(a)') 'Calculating spin-orbit couplings.'
        else
          write(u_log,'(a)') 'Not calculating spin-orbit couplings.'
        endif
      endif

      if (ctrl%calc_soc/=1) then
        n=0
        do i=1,ctrl%maxmult
          if (ctrl%nstates_m(i)>0) n=n+1
        enddo
        if (n>1) then
          write(u_log,'(a)') 'Warning: More than one multiplicity, but spin-orbit couplings are disabled.'
        endif
      endif
      write(u_log,*)







    if (printlevel>1) then
      if (ctrl%write_overlap==0) then
        write(u_log,'(a)') 'Not writing overlap matrix.'
      else
        write(u_log,'(a)') 'Writing overlap matrix.'
      endif
      ! ---------------------
      if (ctrl%write_grad==0) then
        write(u_log,'(a)') 'Not writing gradients.'
      else
        write(u_log,'(a)') 'Writing gradients.'
        if (ctrl%output_format==1) then
          write(u_log,'(a)') 'Error: Currently, NetCDF output is not compatible with write_grad'
          stop 1
        endif
      endif
      ! ---------------------
      if (ctrl%write_NACdr==0) then
        write(u_log,'(a)') 'Not writing nonadiabatic couplings.'
      else
        write(u_log,'(a)') 'Writing nonadiabatic couplings.'
        if (ctrl%output_format==1) then
          write(u_log,'(a)') 'Error: Currently, NetCDF output is not compatible with write_NACdr'
          stop 1
        endif
      endif
      ! ---------------------
      if (ctrl%write_property1d==0) then
        write(u_log,'(a)') 'Not writing property vectors.'
      else
        if (ctrl%theodore==1) then
          write(u_log,'(a)') 'Writing property vectors (due to TheoDORE being active).'
        else
          write(u_log,'(a)') 'Writing property vectors.'
        endif
        if (ctrl%output_format==1) then
          write(u_log,'(a)') 'Error: Currently, NetCDF output is not compatible with write_property1d'
          stop 1
        endif
      endif
      ! ---------------------
      if (ctrl%write_property2d==0) then
        write(u_log,'(a)') 'Not writing property matrices.'
      else
        if (ctrl%ionization==1) then
          write(u_log,'(a)') 'Writing property matrices (due to ionization being active).'
        else
          write(u_log,'(a)') 'Writing property matrices.'
        endif
        if (ctrl%output_format==1) then
          write(u_log,'(a)') 'Error: Currently, NetCDF output is not compatible with write_property2d'
          stop 1
        endif
      endif
      ! ---------------------
      write(u_log,*)
        write(u_log,'(a,i6,a,i6)') 'First,   writing to output.dat every ',ctrl%output_steps_stride(1),&
        &' steps if step is >= ',ctrl%output_steps_limits(1)
        write(u_log,'(a,i6,a,i6)') 'Then,    writing to output.dat every ',ctrl%output_steps_stride(2),&
        &' steps if step is >= ',ctrl%output_steps_limits(2)
        write(u_log,'(a,i6,a,i6)') 'Finally, writing to output.dat every ',ctrl%output_steps_stride(3),&
        &' steps if step is >= ',ctrl%output_steps_limits(3)
      write(u_log,*)
    endif

    endif

  ! =====================================================

  ! fill up the other ctrl parameters:
  ! cwd, ezero, scalingfactor, printlevel, RNGseed, dampeddyn, decoherence, decoherence_scheme, decoherence_alpha

    ! current directory ============================================

    call getcwd(ctrl%cwd)

    if (printlevel>1) then
      write(u_log,'(A)') 'Current working directory is'
      write(u_log,'(a,a)') '        ',trim(ctrl%cwd)
      write(u_log,*)
    endif

    ! energy shift ============================================

    line=get_value_from_key('ezero',io)
    if (io==0) then
      read(line,*) ctrl%ezero
    else
      ctrl%ezero=0.d0
    endif

    if (printlevel>1) write(u_log,'(a,1x,f15.9)') 'Shift to the diagonal energies is',ctrl%ezero

    ! convergence threshold ============================================

    line=get_value_from_key('convthre',io)
    if (io==0) then
      read(line,*) ctrl%convthre
    else
      ctrl%convthre=0.0001
    endif

    if (printlevel>1) then
      if (ctrl%integrator==0) then
        write(u_log,'(a,1x,f15.9,1x,a)') 'Bulirsch-Stoer convergence threshold:',ctrl%convthre, 'a.u.'
      else if (ctrl%integrator==1) then
        write(u_log,'(a,1x,f15.9,1x,a)') 'Adaptive Velocity Verlet convergence threshold:',ctrl%convthre,'eV'
        write(u_log,'(a,1x,f15.9,1x,a)') 'Minimum time step:',ctrl%dtstep_min,"fs"
        write(u_log,'(a,1x,f15.9,1x,a)') 'Maximum time step:',ctrl%dtstep_max,"fs"
      endif
    endif

    ! scaling ============================================

    line=get_value_from_key('scaling',io)
    if (io==0) then
      read(line,*) ctrl%scalingfactor
    else
      ctrl%scalingfactor=1.d0
    endif
    if (ctrl%scalingfactor<=0.d0) then
      write(0,*) 'Scaling must not be smaller than or equal to 0.0!'
      stop 1
    endif

    ! spin orbit coupling scaling factor==============================
    line=get_value_from_key('soc_scaling',io)
    if (io==0) then
      read(line,*) ctrl%soc_scaling
    else
      ctrl%soc_scaling=1.d0
    endif
    if (ctrl%soc_scaling<=0.d0) then
      write(0,*) 'SOC scaling must not be smaller than or equal to 0.0!'
      stop 1
    endif

    if (printlevel>1) then
      write(u_log,'(a,1x,f6.3)') 'Scaling factor to the Hamiltonian and gradients is',ctrl%scalingfactor
      write(u_log,'(a,1x,f6.3)') 'SOC scaling factor is',ctrl%soc_scaling 
      write(u_log,*)
    endif

    ! random number seed ============================================

    line=get_value_from_key('rngseed',io)
    if (io==0) then
      read(line,*) traj%rngseed
    else
      traj%rngseed=1099279      ! some prime number
    endif
    call init_random_seed(traj%rngseed)
    traj%randnum=0.d0

    if (printlevel>1) then
      write(u_log,'(a,1x,i10)') 'Random number seed is',traj%rngseed
      write(u_log,*)
    endif

    ! damped dynamics ============================================

    line=get_value_from_key('dampeddyn',io)
    if (io==0) then
      read(line,*)ctrl%dampeddyn
    else
      ctrl%dampeddyn=1.d0
    endif
    if (ctrl%dampeddyn>1.d0) then
      write(0,*) 'Dampeddyn must not be larger than 1.0!'
      stop 1
    endif
    if (ctrl%dampeddyn<0.d0) then
      write(0,*) 'Dampeddyn must not be smaller than 0.0!'
      stop 1
    endif

    if (printlevel>1) then
      write(u_log,'(a,1x,f6.3)') 'Damping factor for kinetic energy is',ctrl%dampeddyn
      write(u_log,*)
    endif

    ! decoherence ============================================
    ! alternatively one can use the decoherence/nodecoherence or decoherence_scheme keywords

    line=get_value_from_key('decoherence',io)
    if (io==0) then
      if (ctrl%method==0) then
        ctrl%decoherence=1
        ctrl%decoherence_alpha=0.1d0
      else if (ctrl%method==1) then
        ctrl%decoherence=11
        ctrl%decoherence_alpha=0.1d0
        ctrl%decoherence_beta=1.0d0
      endif
    else
      ctrl%decoherence=0
    endif
    line=get_value_from_key('nodecoherence',io)
    if (io==0) then
      ctrl%decoherence=0
    endif

    line=get_value_from_key('decoherence_scheme',io)
    if (io==0) then
      select case (trim(line))
        case ('none')
          ctrl%decoherence=0
        case ('edc')
          ctrl%decoherence=1
        case ('afssh') 
          ctrl%decoherence=2
        case ('dom')
          ctrl%decoherence=11
          ctrl%decoherence_alpha=0.1d0
          ctrl%decoherence_beta=1.0d0
        case('edc_legacy')
          ctrl%decoherence=-1
        case default
          write(0,*) 'Unknown keyword ',trim(line),' to "decoherence_scheme"!'
          stop 1
      endselect
    endif   

    if (ctrl%decoherence==0) write(0,*) 'Warning: Decoherence correction turned off!'
    if (ctrl%decoherence==1) then
      line=get_value_from_key('decoherence_param',io)
      if (io==0) read(line,*) ctrl%decoherence_alpha
      if (ctrl%decoherence_alpha<0.d0) then
        write(0,*) 'Decoherence parameter must not be smaller than 0.0!'
        stop 1
      endif
    endif
    if (ctrl%decoherence==11) then
      line=get_value_from_key('decoherence_param_alpha',io)
      if (io==0) read(line,*) ctrl%decoherence_alpha
      if (ctrl%decoherence_alpha<0.d0) then
        write(0,*) 'Decoherence parameter must not be smaller than 0.0!'
        stop 1
      endif
      line=get_value_from_key('decoherence_param_beta',io)
      if (io==0) read(line,*) ctrl%decoherence_beta
    endif

    if (printlevel>1) then
      if (ctrl%decoherence==-1) then
        write(u_log,'(a)') 'Decoherence is 1 (EDC by Granucci, Persico, Zoccante)'
        write(u_log,'(a)') '(legacy version with equation from original paper)'
        write(u_log,'(a,1x,f6.3,1x,a)') 'Decoherence constant is',ctrl%decoherence_alpha,'Hartree'
        write(u_log,*)
      elseif (ctrl%decoherence==1) then
        write(u_log,'(a)') 'Decoherence is 1 (EDC by Granucci, Persico, Zoccante)'
        write(u_log,'(a,1x,f6.3,1x,a)') 'Decoherence constant is',ctrl%decoherence_alpha,'Hartree'
        write(u_log,*)
      elseif (ctrl%decoherence==2) then
        write(u_log,'(a)') 'Decoherence is 2 (A-FSSH by Jain, Alguire, Subotnik)'
        write(u_log,*)
        if (ctrl%compat_mode==1) then
          write(u_log,'(a)') 'A-FSSH decoherence scheme cannot be used in compatibility mode = 1 !'
          write(u_log,*)
          stop 1
        endif
      elseif (ctrl%decoherence==11) then
        write(u_log,'(a)') 'Decoherence is 11 (Decay of Mixing by Zhu, Nangia, Jasper, Truhlar)'
        write(u_log,'(a,1x,f6.3,1x,a)') 'Alpha decoherence constant is',ctrl%decoherence_alpha,'Hartree'
        write(u_log,'(a,1x,f6.3,1x,a)') 'Beta decoherence constant is',ctrl%decoherence_beta,'Hartree'
        write(u_log,*)
      else
        write(u_log,'(a)') 'Decoherence is OFF'
        write(u_log,*)
      endif
    endif

    line=get_value_from_key('decotime_method',io)
    if (io==0) then
      select case (trim(line))
        case ('csdm')
          ctrl%decotime_method=0
        case ('scdm')
          ctrl%decotime_method=1
        case ('edc')
          ctrl%decotime_method=2
        case ('sd')
          ctrl%decotime_method=3
        case ('fp1')
          ctrl%decotime_method=4
        case ('fp2')
          ctrl%decotime_method=5
        case default
          write(0,*) 'Unknown keyword ',trim(line),' to "decotime_method"!'
          stop 1
      endselect
    else
      ctrl%decotime_method=0
    endif

    if (ctrl%decotime_method==5) then !for fp2 method, we need to get gaussian width parameter
      line=get_value_from_key('gaussian_width',io)
      if (io==0) read(line,*) ctrl%gaussian_width
      if (ctrl%gaussian_width<0.d0) then
        write(0,*) 'Gaussian width must not be smaller than 0.0!'
        stop 1
      endif
    endif

    if ( (printlevel>1) .and. (ctrl%method==1)) then   ! should this only be printed if SCP?
      if (ctrl%decotime_method==0) then
        write(u_log,'(a)') 'Decoherence time is computed with CSDM method'
      elseif (ctrl%decotime_method==1) then
        write(u_log,'(a)') 'Decoherence time is computed with SCDM method'
      elseif (ctrl%decotime_method==2) then
        write(u_log,'(a)') 'Decoherence time is computed with EDC (energy based decoherence) method'
      elseif (ctrl%decotime_method==3) then
        write(u_log,'(a)') 'Decoherence time is computed with SD (stochastic decoherence) method'
      elseif (ctrl%decotime_method==4) then
        write(u_log,'(a)') 'Decoherence time is computed with FP1 (force momentum 1) method'
      elseif (ctrl%decotime_method==5) then
        write(u_log,'(a)') 'Decoherence time is computed with FP2 (force momentum 2) method'
      endif
    endif

    ! phase tracking algorithm, actually not necessary in release SHARC version
    ! (see IJQC paper)
    line=get_value_from_key('notrack_phase',io)
    if (io==0) then
      ctrl%track_phase=0
      if (printlevel>1) then
        write(u_log,'(a)') 'Phase tracking of MCH-diagonal transformation matrix (U) is OFF'
        write(u_log,*)
      endif
    else
      ctrl%track_phase=1
    endif
    line=get_value_from_key('track_phase',io)
    if (io==0) then
      ctrl%track_phase=1
    endif

    ! surface hopping procedure
    ctrl%hopping_procedure=1
    line=get_value_from_key('hopping_procedure',io)
    if (io==0) then
      select case (trim(line))
        case ('standard') 
          ctrl%hopping_procedure=1
          if (printlevel>1) then
            write(u_log,'(a)') 'Surface Hopping procedure is SHARC (Mai, Marquetand, Gonzalez)'
          endif
        case ('sharc') 
          ctrl%hopping_procedure=1
          if (printlevel>1) then
            write(u_log,'(a)') 'Surface Hopping procedure is SHARC (Mai, Marquetand, Gonzalez)'
          endif
        case ('gfsh') 
          ctrl%hopping_procedure=2
          if (printlevel>1) then
            write(u_log,'(a)') 'Surface Hopping procedure is GFSH (Wang, Trivedi, Prezhdo)'
          endif
        case ('off') 
          ctrl%hopping_procedure=0
          ctrl%ekincorrect=0
          if (printlevel>1) then
            write(u_log,'(a)') 'Surface Hopping is OFF (will stay in initial diagonal state)'
          endif
        case default
          write(0,*) 'Unknown keyword ',trim(line),' to "hopping_procedure"!'
          stop 1
      endselect
    endif

    ! force hops to ground state
    ctrl%force_hop_to_gs=-1.
    line=get_value_from_key('force_hop_to_gs',io)
    if (io==0) then
      read(line,*)ctrl%force_hop_to_gs
      ctrl%force_hop_to_gs=ctrl%force_hop_to_gs/au2eV
      if (printlevel>1) then
        write(u_log,'(a,F6.3,a)') 'Forcing hops to lowest state if active-lowest gap is <='&
        &,ctrl%force_hop_to_gs*au2eV,'eV'
        write(u_log,*)
      endif
    endif
    write(u_log,*)

  ! =====================================================

  ! find constraints
    ctrl%constraints_tol=1e-7
    ctrl%do_constraints=0
    line=get_value_from_key('rattle',io)
    if (io==0) then
      ctrl%do_constraints=1
    endif
    line=get_value_from_key('norattle',io)
    if (io==0) then
      ctrl%do_constraints=0
    endif

     line=get_value_from_key('rattletolerance',io)
     if (io==0) then
       read(line,*) ctrl%constraints_tol
     endif

    if (ctrl%do_constraints==1) then

      ! look up rattle filename
      line=get_value_from_key('rattlefile',io)
      if (io==0) then
        ! extract string within quotes (can contain spaces)
        call get_quoted(line,rattlefilename)
        filename=trim(rattlefilename)
      else
        ! default rattle filename
        filename='rattle'
        rattlefilename=filename
      endif
    ! open the rattle file
      open(u_i_rattle,file=filename, status='old', action='read', iostat=io)
      if (io/=0) then
        write(0,*) 'Could not find rattle file "',trim(filename),'"!'
        stop 1
      endif

    ! find number of constraints
      nlines=0
      do
        read(u_i_rattle,'(A)',iostat=io) line
        if (io/=0) exit
        call split(line,' ', values,n)
        deallocate(values)
        if (n<2) exit
        nlines=nlines+1
      enddo
      ctrl%n_constraints=nlines
      if (ctrl%n_constraints==0) then
        write(0,*) 'No constraints found in ',rattlefilename
        stop 1
      endif

      allocate( ctrl%constraints_ca(ctrl%n_constraints,2) )
      allocate( ctrl%constraints_dist_c(ctrl%n_constraints) )

      rewind(u_i_rattle)
      do i=1,ctrl%n_constraints
        read(u_i_rattle,'(A)') line
        call split(line,' ',values,n)
        if (n<2) then
          write(0,*) 'Problem reading the rattle file!'
          stop 1
        endif
        read(values(1),*) ctrl%constraints_ca(i,1)
        read(values(2),*) ctrl%constraints_ca(i,2)
        if ( (ctrl%constraints_ca(i,1)>ctrl%natom).or.&
            &(ctrl%constraints_ca(i,2)>ctrl%natom).or.&
            &(ctrl%constraints_ca(i,1)<1).or.&
            &(ctrl%constraints_ca(i,2)<1) ) then
          write(0,*) 'Atom index in constraint wrong:'
          write(0,*) trim(line)
          stop 1
        endif
        if (n>=3) then
          read(values(3),*) ctrl%constraints_dist_c(i)
          ctrl%constraints_dist_c(i)=ctrl%constraints_dist_c(i)**2
        else
          ctrl%constraints_dist_c(i)=-1.d0
        endif
      enddo

    else   ! if no constraints
      ctrl%n_constraints=0
    endif

    ! surface switching procedure
    ctrl%switching_procedure=1
    line=get_value_from_key('switching_procedure',io)
    if (io==0) then
      select case (trim(line))
        case ('csdm')
          ctrl%switching_procedure=1
          if (printlevel>1) then
            write(u_log,'(a)') 'Surface Switching procedure is CSDM (Shu, Zhang, Mai, Sun, Truhlar, Gonzalez)'
            write(u_log,*)
          endif
        case ('scdm')
          ctrl%switching_procedure=2
          if (printlevel>1) then
            write(u_log,'(a)') 'Surface Switching procedure is SCDM (Shu, Zhang, Mai, Sun, Truhlar, Gonzalez)'
            write(u_log,*)
          endif
        case ('ndm')
          ctrl%switching_procedure=3
          if (printlevel>1) then
            write(u_log,'(a)') 'Surface Switching procedure is NDM (Shu, Zhang, Mai, Sun, Truhlar, Gonzalez)'
            write(u_log,*)
          endif
        case ('off')
          ctrl%switching_procedure=0
          ctrl%ekincorrect=0
          if (printlevel>1) then
            write(u_log,'(a)') 'Surface Switching if OFF (will stay in initial diagonal state)'
            write(u_log,*)
          endif
        case default
          write(0,*) 'Unknown keyword ',trim(line),' to "switching_procedure"!'
          stop 1
      endselect
    endif
  




  ! =====================================================

  ! Reading the molecule

    ! the geomfile was already opened before to determine natom
    rewind(u_i_geom)
    ! read in line by line
    ! reads element, atomic number, position(3), mass
    do i=1,ctrl%natom
      read(u_i_geom,'(A)') line
      call split(line,' ',values,n)
      if (n<6) then
        write(0,*) 'Problem reading the geometry file!'
        stop 1
      endif
      traj%element_a(i)=values(1)(1:2)
      read(values(2),*) traj%atomicnumber_a(i)
      do j=1,3
        read(values(2+j),*) traj%geom_ad(i,j)
      enddo
      read(values(6),*) traj%mass_a(i)
      deallocate(values)
    enddo
    close(u_i_geom)

    if (printlevel>0) then
      write(u_log,*) '============================================================='
      write(u_log,*) '                            Geometry'
      write(u_log,*) '============================================================='
      if (printlevel>1) then
        write(u_log,'(3a)') 'Reading from geometry file: "',trim(geomfilename),'"'
        write(u_log,*) 'Geometry (Bohr):'
        write(u_log,'(A2,1X,3(A9,1X),3X,A3,1X,A12)') 'El','x','y','z','#','mass'
        do i=1,ctrl%natom
          write(u_log,'(A2,1X,3(F9.6,1X),3X,F4.0,1X,F12.6)') traj%element_a(i),&
          &(traj%geom_ad(i,j),j=1,3),traj%atomicnumber_a(i),traj%mass_a(i)
        enddo
      endif
      write(u_log,*)
    endif

  ! =====================================================

    if (ctrl%do_constraints==1) then

      do i=1,ctrl%n_constraints
        if (ctrl%constraints_dist_c(i)<0.d0) then
          a=0.
          k=ctrl%constraints_ca(i,1)
          n=ctrl%constraints_ca(i,2)
          do j=1,3
            a=a+(traj%geom_ad(k,j)-traj%geom_ad(n,j))**2
          enddo
          ctrl%constraints_dist_c(i)=a
        endif
      enddo


      if (printlevel>1) then
        write(u_log,'(3a)') 'Constraints from file: "',trim(rattlefilename),'"'
        write(u_log,*) 'List of constraints:'
        write(u_log,'(3X,A6,3X,A6,3X,A15)') 'Atom A','Atom B','Distance (Bohr)'
        do i=1,ctrl%n_constraints
          write(u_log,'(3X,I6,3X,I6,3X,F15.6)') ctrl%constraints_ca(i,1),&
          &ctrl%constraints_ca(i,2),dsqrt(ctrl%constraints_dist_c(i))
        enddo
      endif
      write(u_log,*)

    endif










  ! =====================================================

  ! Reading the velocities

    line=get_value_from_key('veloc',io)
    if (io==0) then
      call split(line,' ',values,n)
      select case (trim(values(1)))
        ! initialize velocities to zero
        case ('zero')
          if (printlevel>1) write(u_log,'(A)') 'Velocities are assumed to be zero!'
          traj%veloc_ad=0.
        ! initialize velocities randomly
        case ('random')
          if (size(values)<2) then
            write(0,*) 'Give kinetic energy per atom (in eV) for "veloc random"!'
            stop 1
          endif
          read(values(2),*) a   ! a is kinetic energy per atom in hartree
          a=dabs(a)
          if (printlevel>1) write(u_log,'(A,1X,F7.4,1X,A)') 'Picking velocities at random with E_atom=',a,'eV'
          a=a/au2eV
          call set_random_velocities(traj,ctrl,a)
          call init_random_seed(traj%rngseed)
        ! initialize velocities from file
        case ('external')

          line=get_value_from_key('velocfile',io)
          if (io==0) then
            call get_quoted(line,geomfilename)
            filename=trim(geomfilename)
          else
            filename='veloc'
          endif
          open(u_i_veloc,file=filename, status='old', action='read', iostat=io)

          if (printlevel>1) write(u_log,'(3A)') 'Reading velocities from file "',trim(filename),'"'
          if (io/=0) then
            write(0,*) 'Could not find velocity file!'
            stop 1
          endif
          do i=1,ctrl%natom
            read(u_i_veloc,'(A)') line
            call split(line,' ',values,n)
            if (n<3) then
              write(0,*) 'Problem reading the velocity file!'
              stop 1
            endif
            do j=1,3
              read(values(j),*) traj%veloc_ad(i,j)
            enddo
          enddo
          close(u_i_veloc)
        case default
          write(0,*) 'Options for keyword veloc are "zero", "random" and "external"!'
          stop 1
      endselect
      deallocate(values)
    else
      if (printlevel>1) write(u_log,'(A)') 'No keyword "veloc", assuming zero initial velocities!'
      traj%veloc_ad=0.d0
    endif

    if (printlevel>1) then
      write(u_log,*) 'Velocities (Bohr/atu):'
      write(u_log,'(A2,1X,3(A9,1X))') 'El','dx/dt','dy/dt','dz/dt'
      do i=1,ctrl%natom
        write(u_log,'(A2,1X,3(F9.6,1X))') traj%element_a(i),&
        &(traj%veloc_ad(i,j),j=1,3)
      enddo
      write(u_log,*)
    endif

  ! =====================================================

  ! process the atom mask
  allocate( ctrl%atommask_a(ctrl%natom))
  ctrl%atommask_a=.true.
  line=get_value_from_key('atommask',io)
  if (io==0) then
    ! ------- some combinations have no effect
    if (ctrl%ekincorrect==3) then
      write(u_log,*) 'HINT: "ekincorrect parallel_nac" ignores "atommask".'
    endif
    if (ctrl%ekincorrect==4) then
      write(u_log,*) 'HINT: "ekincorrect parallel_diff" ignores "atommask".'
    endif
    if (ctrl%ekincorrect==5) then
      write(u_log,*) 'HINT: "ekincorrect parallel_pnac" ignores "atommask".'
    endif
    if (ctrl%ekincorrect==6) then
      write(u_log,*) 'HINT: "ekincorrect parallel_pdiff" ignores "atommask".'
    endif
    if (ctrl%ekincorrect==7) then
      write(u_log,*) 'HINT: "ekincorrect parallel_enac" ignores "atommask".'
    endif
    if (ctrl%ekincorrect==8) then
      write(u_log,*) 'HINT: "ekincorrect parallel_penac" ignores "atommask".'
    endif
    if (ctrl%decoherence==2) then
      write(u_log,*) 'HINT: "decoherence_scheme afssh" ignores "atommask".'
    endif
    if (ctrl%reflect_frustrated==2) then
      write(u_log,*) 'HINT: "reflect_frustrated parallel_nac" ignores "atommask".'
    endif
    if (ctrl%reflect_frustrated==3) then
      write(u_log,*) 'HINT: "reflect_frustrated parallel_nac" ignores "atommask".'
    endif
    ! TODO: add more as needed
    ! -------
    call split(line,' ',values,n)
    select case (trim(values(1)))
      ! initialize atommask from file
      case ('external')
          line=get_value_from_key('atommaskfile',io)
          if (io==0) then
            call get_quoted(line,geomfilename)
            filename=trim(geomfilename)
          else
            filename='atommask'
          endif
          open(u_i_atommask,file=filename, status='old', action='read', iostat=io)

          if (printlevel>1) write(u_log,'(3A)') 'Reading atom mask from file "',trim(filename),'"'
          if (io/=0) then
            write(0,*) 'Could not find atom mask file!'
            stop 1
          endif
          do i=1,ctrl%natom
            read(u_i_atommask,*) ctrl%atommask_a(i)
          enddo
          close(u_i_atommask)
      case ('none')
          continue
      case default
        write(0,*) 'Unknown option for keyword atommask!'
        stop 1
    endselect
    deallocate(values)
  endif
  if (printlevel>1) then
    if (printlevel>2) then
      write(u_log,*) 'Atom mask:'
      do i=1,ctrl%natom
        write(u_log,*) i,traj%element_a(i),ctrl%atommask_a(i)
      enddo
    else
      write(u_log,*) 'Atom mask (only active):'
      do i=1,ctrl%natom
        if (ctrl%atommask_a(i)) write(u_log,*) i,traj%element_a(i),ctrl%atommask_a(i)
      enddo
    endif
    if (ctrl%ekincorrect==1) then
      write(u_log,*) 'Atom mask will be used for kinetic energy adjustment.'
    endif
    if (ctrl%decoherence==1) then
      write(u_log,*) 'Atom mask will be used for kinetic energy in EDC decoherence.'
    endif
    if (ctrl%reflect_frustrated==1) then
      write(u_log,*) 'Atom mask will be used for reflection after frustrated hops.'
    endif
  endif

  ! =====================================================

  ! Reading the coefficients

    if (printlevel>0) then
      write(u_log,*) '============================================================='
      write(u_log,*) '               Initial State and Coefficients'
      write(u_log,*) '============================================================='
      write(u_log,*) 
    endif

    select case (ctrl%surf)
      ! =====================================================
      case (0)  !SHARC
        if (printlevel>1) write(u_log,'(A)') 'Setting state and coefficients for SHARC dynamics.'
        if (printlevel>1) write(u_log,*)
        line=get_value_from_key('state',io)
        if (io/=0) then
          write(0,*) 'You have to give an initial state!'
          stop 1
        endif
        call split(line,' ',values,n)
        if (n<1) then
          write(0,*) 'You have to give an initial state!'
          stop 1
        endif
        read(values(1),*) j
        if (j>ctrl%nstates) then
          write(0,*) 'Initial state can''t be larger than nstates=',ctrl%nstates,'!'
          stop 1
        endif
        if (j<=0) then
           write(0,*) 'Initial state can''t be zero or negative!'
           stop 1
        endif
        if (n==1) then
          write(0,*) 'Please specify representation of initial state!'
          stop 1
        endif
        select case (trim(values(2)))
          case ('mch')
            traj%state_MCH=j
            ctrl%staterep=1
            if (printlevel>1) write(u_log,'(a,1x,i3,1x,a)') 'Initial state is ',traj%state_mch,'in the MCH basis. '
          case ('diag')
            traj%state_diag=j
            ctrl%staterep=0
            if (printlevel>1) write(u_log,'(a,1x,i3,1x,a)') 'Initial state is ',traj%state_diag,'in the diag basis.'
          case default
            write(0,*) 'Unknown keyword for "state"!'
            stop 1
        endselect
        deallocate(values)

        line=get_value_from_key('coeff',io)
        if (io/=0) then
          write(0,*) 'Please specify origin of coefficients ("auto" or "external")!'
          stop 1
        endif
        select case (trim(line))
          case ('auto')
            ctrl%initcoeff=2+ctrl%staterep
            if (printlevel>1) write(u_log,'(a)') 'Initial coefficients will be set automatically.'
          case ('external')
            ctrl%initcoeff=ctrl%staterep
            if (printlevel>1) write(u_log,'(a)') 'Initial coefficients will be read from file.'
          case default
            write(0,*) 'Unknown keyword for "coeff"!'
            stop 1
        endselect
      ! =====================================================
      case (1)  !non-SHARC
        if (printlevel>1) write(u_log,'(A)') 'Setting state and coefficients for dynamics on MCH surfaces.'
        if (printlevel>1) write(u_log,*)
        line=get_value_from_key('state',io)
        if (io/=0) then
          write(0,*) 'You have to give an initial state!'
          stop 1
        endif
        read(line,*) j
        if (j>ctrl%nstates) then
          write(0,*) 'Initial state can''t be larger than nstates=',ctrl%nstates,'!'
          stop 1
        endif
        if (j<=0) then
           write(0,*) 'Initial state can''t be zero or negative!'
           stop 1
        endif
        traj%state_MCH=j
        ctrl%staterep=1
        if (printlevel>1) write(u_log,'(a,1x,i3)') 'Initial state is ',traj%state_MCH

        line=get_value_from_key('coeff',io)
        if (io/=0) then
          write(0,*) 'Please specify origin of coefficients ("auto" or "external")!'
          stop 1
        endif
        select case (trim(line))
          case ('auto')
            ctrl%initcoeff=3
            if (printlevel>1) write(u_log,'(a)') 'Initial coefficients will be set automatically.'
          case ('external')
            ctrl%initcoeff=ctrl%staterep
            if (printlevel>1) write(u_log,'(a)') 'Initial coefficients will be read from file.'
          case default
            write(0,*) 'Unknown keyword for "coeff"!'
            stop 1
        endselect
      ! =====================================================
      case default
    endselect

    ! =====================================================
    if (ctrl%initcoeff<2) then

      ! get filename
      line=get_value_from_key('coefffile',io)
      if (io==0) then
        call get_quoted(line,geomfilename)
        filename=trim(geomfilename)
      else
        filename='coeff'
      endif
      open(u_i_coeff,file=filename, status='old', action='read', iostat=io)

      ! read the coefficients
      if (printlevel>1) write(u_log,'(3a)') 'Reading from coefffile "',trim(filename),'"'
      if (io/=0) then
        write(0,*) 'Could not find coefficients file!'
        stop 1
      endif
      do i=1,ctrl%nstates
        read(u_i_coeff,'(A)') line
        call split(line,' ',values,n)
        if (n<2) then
          write(0,*) 'Problem reading the coefficients file!'
          stop 1
        endif
        read(line,*) a,b
        traj%coeff_MCH_s(i)=dcmplx(a,b)
      enddo
      close(u_i_coeff)

      ! Normalize the coefficients
      a=0.
      do i=1,ctrl%nstates
        a=a+real(traj%coeff_MCH_s(i)*conjg(traj%coeff_MCH_s(i)))
      enddo
      if (a<1.d-9) then
        write(0,*) 'Sum of initial coefficients is very small!'
        stop
      endif
      a=1.d0/sqrt(a)
      traj%coeff_MCH_s=dcmplx(a,0.d0)*traj%coeff_MCH_s
      if (ctrl%initcoeff==0) then
        traj%coeff_diag_s=traj%coeff_MCH_s
      endif
    elseif (ctrl%initcoeff==2) then
      traj%coeff_diag_s=dcmplx(0.d0,0.d0)
      traj%coeff_diag_s(traj%state_diag)=dcmplx(1.d0,0.d0)
    elseif (ctrl%initcoeff==3) then
      traj%coeff_MCH_s=dcmplx(0.d0,0.d0)
      traj%coeff_MCH_s(traj%state_MCH)=dcmplx(1.d0,0.d0)
    endif

    if (printlevel>1) then
      write(u_log,*)
      ! Writing the coefficients
      select case (ctrl%initcoeff)
        case (0)
          write(u_log,*) 'Coefficients (diag):'
          write(u_log,'(a3,1x,A12,1X,A12)') '#','Real(c)','Imag(c)'
          do i=1,ctrl%nstates
            write(u_log,'(i3,1x,F12.9,1X,F12.9)') i,traj%coeff_diag_s(i)
          enddo
          write(u_log,*)
        case (1)
          write(u_log,*) 'Coefficients (MCH):'
          write(u_log,'(a3,1x,A12,1X,A12)') '#','Real(c)','Imag(c)'
          do i=1,ctrl%nstates
            write(u_log,'(i3,1x,F12.9,1X,F12.9)') i,traj%coeff_MCH_s(i)
          enddo
          write(u_log,*)
          if (ctrl%surf==0) then
            write(u_log,*) 'Coefficients will be transformed after the initial energy calculation.'
            write(u_log,*)
          endif
        case (2)
          write(u_log,*) 'Coefficients (diag):'
          write(u_log,'(a3,1x,A12,1X,A12)') '#','Real(c)','Imag(c)'
          do i=1,ctrl%nstates
            write(u_log,'(i3,1x,F12.9,1X,F12.9)') i,traj%coeff_diag_s(i)
          enddo
          write(u_log,*)
        case (3)
          write(u_log,*) 'Coefficients (MCH):'
          write(u_log,'(a3,1x,A12,1X,A12)') '#','Real(c)','Imag(c)'
          do i=1,ctrl%nstates
            write(u_log,'(i3,1x,F12.9,1X,F12.9)') i,traj%coeff_MCH_s(i)
          enddo
          write(u_log,*)
          if (ctrl%surf==0) then
            write(u_log,*) 'Coefficients will be transformed after the initial energy calculation.'
          endif
          write(u_log,*)
      endselect
    endif

  ! =====================================================

  ! check for laser keyword

    line=get_value_from_key('laser',io)
    if (io==0) then
      select case (trim(line))
        case ('none') 
          ctrl%laser=0
        case ('internal')
          ctrl%laser=1
        case ('external')
          ctrl%laser=2
        case default
          ctrl%laser=0
      endselect
    else
      ctrl%laser=0
    endif

    if (printlevel>0) then
      write(u_log,*) '============================================================='
      write(u_log,*) '                       Laser Field'
      write(u_log,*) '============================================================='
      if (printlevel>1) then
        select case (ctrl%laser)
          case (0)
            write(u_log,'(a)') 'No laser field will be applied.'
          case (1)
            write(u_log,'(a)') 'Laser field is calculated from function internal_laserfield().'
          case (2)
            write(u_log,'(a)') 'Laser field is read from file.'
        endselect
      endif
      write(u_log,*)
    endif

    if (ctrl%laser/=0) then
      line=get_value_from_key('laserwidth',io)
      if (io==0) then
        read(line,*) ctrl%laser_bandwidth
        ctrl%laser_bandwidth=ctrl%laser_bandwidth/au2eV
      else
        ctrl%laser_bandwidth=1./au2eV
      endif
    endif

    ctrl%calc_dipole=1

    if (ctrl%laser/=0) then
      ctrl%calc_dipole=1
    endif

    ctrl%dipolegrad=0
    ctrl%calc_dipolegrad=-1

    if (ctrl%laser/=0) then
      line=get_value_from_key('dipole_gradient',io)
      if (io==0) then
        ctrl%dipolegrad=1
        ctrl%calc_dipolegrad=0
        if ((ctrl%calc_grad>=1).or.(ctrl%calc_nacdr>=1)) then
          ctrl%calc_dipolegrad=1
          ctrl%eselect_dmgrad=ctrl%eselect_nac
          if (ctrl%calc_second==1) then
            ctrl%calc_dipolegrad=2
          endif
        endif
      endif
      line=get_value_from_key('nodipole_gradient',io)
      if (io==0) then
        ctrl%dipolegrad=0
        ctrl%calc_dipolegrad=-1
      endif
    endif


    if (ctrl%laser==2) then

      ! get filename
      line=get_value_from_key('laserfile',io)
      if (io==0) then
        call get_quoted(line,geomfilename)
        filename=trim(geomfilename)
      else
        filename='laser'
      endif
      open(u_i_laser,file=filename, status='old', action='read', iostat=io)

      if (printlevel>1) write(u_log,'(3a)') 'Reading from laser file "',trim(filename),'"'
      if (io/=0) then
        write(0,*) 'Could not find laser file!'
        stop 1
      endif

      ! find number of lasers
      read(u_i_laser,'(A)',iostat=io) line
      if (io/=0) then
        write(0,*) 'EOF encountered during read of laser file!'
        stop 1
      endif
      call split(line,' ',values,n)
      ctrl%nlasers=n-7
      if (ctrl%nlasers<1) then
        write(0,*) 'No central energies for lasers found in ',filename
        stop 1
      endif
      rewind(u_i_laser)

      allocate(ctrl%laserfield_td(ctrl%nsteps*ctrl%nsubsteps+1,3))
      allocate(ctrl%laserenergy_tl(ctrl%nsteps*ctrl%nsubsteps+1,ctrl%nlasers))

      ! read laser field line by line, checking the time with the substeps given above
      do i=1,ctrl%nsteps*ctrl%nsubsteps+1
        read(u_i_laser,'(A)',iostat=io) line
        if (io/=0) then
          write(0,*) 'EOF encountered during read of laser file!'
          stop 1
        endif
        call split(line,' ',values,n)
        if (n<8) then
          write(0,*) 'Laser file malformatted! Line=',i
          stop 1
        endif
        read(values(1),*) a
        if (i==1) then
          if (dabs(a)>0.001d0) then
            write(0,*) 'Laser field must start at t=0 fs!'
            stop 1
          endif
        endif
        b=ctrl%dtstep/ctrl%nsubsteps
        if ( dabs(a-b*(i-1))>0.001d0) then 
          write(0,*) 'Laser field spacing does not match substep spacing!'
          stop 1
        endif
        do j=1,3
          read(values(2*j),*) a
          read(values(2*j+1),*) b
          ctrl%laserfield_td(i,j)=dcmplx(a,b)
        enddo
        do j=1,ctrl%nlasers
          read(values(7+j),*) a
          ctrl%laserenergy_tl(i,j)=dcmplx(a,0.d0)
        enddo
      enddo
      close(u_i_laser)

      if (printlevel>1) then
        write(u_log,'(a,1x,i8,1x,a)') 'Laser field with',(ctrl%nsteps*ctrl%nsubsteps+1), 'steps has been read successfully.'
!         n=sizeof(ctrl%laserfield_td)
!         write(u_log,'(a,1x,i10,1x,a)') 'Using',n,'bytes for laser data.'
        write(u_log,'(a)') 'Step size has been checked.'
        write(u_log,'(a,1x,i2,1x,a)') 'Laser central frequencies for',ctrl%nlasers,'lasers read.'
        write(u_log,*)
        if (ctrl%dipolegrad==1) then
          write(u_log,'(a)') 'Will include the cartesian gradient of the dipole moments in the gradient transformation.'
          select case (ctrl%calc_dipolegrad)
            case (-1) 
              write(0,'(a)') 'Internal error 4'
              stop
            case (0)
              write(u_log,'(a)') 'Calculating all dipole moment gradients.'
            case (1)
              write(u_log,'(a)') 'Selecting dipole moment gradients in first quantum chemistry calculation.'
            case (2)
              write(u_log,'(a)') 'Selecting dipole moment gradients in second quantum chemistry calculation.'
            case default
              write(0,'(a)') 'Internal error 4'
              stop
          endselect
        else
          write(u_log,'(a)') 'Gradients of dipole moments will be neglected.'
        endif
      endif
    endif

  ! =====================================================

  ! check for thermostat
  line=get_value_from_key('thermostat',io)
    if (io==0) then
      select case (trim(line))
        case ('none')
          ctrl%thermostat=0
        case ('langevin')
          ctrl%thermostat=1
        case default
          ctrl%thermostat=0
      endselect
    else
      ctrl%thermostat=0
    endif

   if (printlevel>0) then
      write(u_log,*) '============================================================='
      write(u_log,*) '                       Thermostat'
      write(u_log,*) '============================================================='
      if (printlevel>1) then
        select case (ctrl%thermostat)
          case (0)
            write(u_log,'(a)') 'No thermostat will be applied.'
          case (1)
            write(u_log,'(a)') 'Langevin thermostat will be applied.'
            write(u_log,'(a)') 'Temperature (in K) and friction coeffitient (in m_e*fs^-1): '
        endselect
      endif
    endif


    ! set up values needed for thermostat
    if (ctrl%thermostat/=0) then

      ! random number seed for thermostat
      line=get_value_from_key('rngseed',io)
      if (io==0) then
        read(line,*) traj%rngseed_thermostat !for now: use same rngseed for thermostat as given for initial velocities. Maybe change later.
      else
        traj%rngseed_thermostat=1099279      ! some prime number
      endif
      call init_random_seed_thermostat(traj%rngseed_thermostat)
     ! call srand(traj%rngseed_thermostat) alternatively
      if (ctrl%thermostat==1) then
        allocate (traj%thermostat_random(2*((3*ctrl%natom+1)/2))) ! allocate randomnes for all atoms in all directions
      endif

      ! restart with same random number sequence?
      ! default is restarting with same random number sequence
      ctrl%restart_thermostat_random=.true.
      ! look for norestart_thermostat_random keyword
      line=get_value_from_key('norestart_thermostat_random',io)
      if (io==0) then
        ctrl%restart_thermostat_random=.false.
      endif

      ! get temperature
      line=get_value_from_key('temperature',io)
      if (io==0) then
        read(line,*) ctrl%temperature
      else
        ctrl%temperature=293.15 !default temperature
      endif
      write(u_log,'(1x,F11.4)') ctrl%temperature


      ! get constants needed for thermostat
      if (ctrl%thermostat==1) allocate(ctrl%thermostat_const(1)) !allocate right amount of thermostat constants
      line=get_value_from_key('thermostat_const',io) !provide in fs^-1
      if (io==0) then
        call split(line,' ',values,n)
        if (n/=size(ctrl%thermostat_const)) then
          write(0,*) 'Wrong number of thermostat constants!'
          stop 1
        else
          do i=1,n
            read(values(i),*) a
            ctrl%thermostat_const = a  ! set the thermostat constants
            if (printlevel>1) then
              write(u_log,'(1x,ES11.4)') ctrl%thermostat_const(i)
            endif
          enddo
        endif
        deallocate(values)
      else
        write(0,*) 'No thermostat constants given!'
        stop 1
      endif

      if (ctrl%thermostat==1) then !save variance in ctrl%temperature
        ctrl%temperature=2*ctrl%thermostat_const(1)*1.38064852e-23*ctrl%temperature*2.2937126583579e+17*ctrl%dtstep ! var=2*alpha*k_bT*dt (1 J = 2.29...e+17 a.u.)
      endif
    endif


#ifdef __PYSHARC__

    ! catch all options that are not yet compatible with pysharc

    if (ctrl%zpe_correction/=0) then
        write(0,*) 'ZPE correction not yet compatible with pysharc!'
        stop 1
    endif

    if (ctrl%army_ants/=0) then
        write(0,*) 'Army Ants not yet compatible with pysharc!'
        stop 1
    endif

    if (ctrl%time_uncertainty/=0) then
        write(0,*) 'Fewest switches with time uncertainty not yet compatible with pysharc!'
        stop 1
    endif

    if (ctrl%decoherence==11) then
        write(0,*) 'Decay of mixing decoherence not yet compatible with pysharc!'
        stop 1
    endif

    if (ctrl%integrator/=2) then
        write(0,*) 'Bulirsch-Stoer-Hack and adaptive velocity Verlet integrators not yet compatible with pysharc!'
        stop 1
    endif

    if (ctrl%method/=0) then
        write(0,*) 'Self-consistent potential/Ehrenfest not yet compatible with pysharc!'
        stop 1
    endif



#endif



    if (printlevel>0) then
      write(u_log,*)
    endif

  ! =====================================================
  ! This part is added to deal with memory-heavy allocations 
  ! The variables were allocated in subroutine allocate_traj, 
  ! now moved to subroutine additional_allocate_traj(traj,ctrl)
   call additional_allocate_traj(traj,ctrl)


  ! =====================================================

!   ! check for floquet keyword
!     if (printlevel>1) then
!       write(u_log,*) '============================================================='
!       write(u_log,*) '                     Floquet settings'
!       write(u_log,*) '============================================================='
!       write(u_log,'(A)') 'Not yet implemented.'
!     endif

  ! =====================================================
  ! make trajectory hash
    if (printlevel>0) then
      write(u_log,*) '============================================================='
      write(u_log,*) '                     Finalizing Input'
      write(u_log,*) '============================================================='
    endif

    traj%traj_hash=hash_input(ctrl,traj)

    if (printlevel>1) then
      write(u_log,'(a,I12)') 'Trajectory input checksum is',traj%traj_hash
      write(u_log,*) 
    endif

  ! =====================================================

    ! initialize hopping probabilities and acceleration
    traj%hopprob_s=0.d0
    traj%accel_ad=0.d0

  ! =====================================================

  ! convert all numbers to atomic units
    ctrl%tmax=ctrl%tmax/au2fs
    ctrl%dtstep_min=ctrl%dtstep_min/au2fs
    ctrl%dtstep_max=ctrl%dtstep_max/au2fs
    ctrl%dtstep=ctrl%dtstep/au2fs
    ctrl%eselect_grad=ctrl%eselect_grad/au2eV
    ctrl%eselect_nac=ctrl%eselect_nac/au2eV
    traj%mass_a=traj%mass_a/au2u
    if (ctrl%integrator==1) then
      ctrl%convthre=ctrl%convthre/au2eV
    endif
    ! geometry is read in bohrs, as specified in the COLUMBUS geom format
    ! velocity is read in as bohrs/atu
    ! laser field must be in a.u.
    ! langevin friction coefficient in a.u.^-1
    if (ctrl%thermostat==1) ctrl%thermostat_const(1)=ctrl%thermostat_const(1)*au2fs

  ! =====================================================
    ! write some basic information into the data file
    call write_dat_initial(u_dat, ctrl, traj)

  endsubroutine


! =================================================================== !

!> calculates random velocities for each atom with given kinetic energy
!> random velocity is uniformly distributed over the sphere
!> \param a kinetic energy per atom
  subroutine set_random_velocities(traj,ctrl,a)
    ! uniform distribution on a sphere with radius sqrt(2a/m)
    use definitions
    implicit none
    type(trajectory_type),intent(inout) :: traj
    type(ctrl_type),intent(in) :: ctrl
    real*8, intent(in) :: a
    integer :: iatom
    real*8 :: theta,phi,v

    ! per atom
    do iatom=1,ctrl%natom
      call random_number(theta)
      call random_number(phi)
      theta=2*pi*theta
      phi=dacos(2*phi-1)
      v=sqrt(2*a/traj%mass_a(iatom)*au2u)

      traj%veloc_ad(iatom,1)=v*dcos(theta)*dsin(phi)
      traj%veloc_ad(iatom,2)=v*dsin(theta)*dsin(phi)
      traj%veloc_ad(iatom,3)=v*dcos(phi)
    enddo

  endsubroutine

! ===================================================

!> concatenates all relevant input files to a long string and generates a hash from the string
!> the following infos are hashed:
!> - input file (cleaned from comments, newlines, duplicates, ...)
!> - geometry
!> - velocities
!> - coefficients
!> - up to 40 steps from laser
  integer*8 function hash_input(ctrl,traj)
    use definitions
    use input_list
    use misc
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl
    character*4096 :: string
    character*512 :: key
    integer :: i,j,n,io
    integer*8 :: temp

    string=''
    n=get_ncurr()
    do i=1,n
      key=get_key(i,io)
      string=trim(string)//trim(key)
      key=get_value(i,io)
      string=trim(string)//trim(key)
    enddo
    do i=1,ctrl%natom
      write(key,'(a2,6(F9.6))') traj%element_a(i),(traj%geom_ad(i,j),j=1,3),(traj%veloc_ad(i,j),j=1,3)
      string=trim(string)//trim(key)
    enddo
    do i=1,ctrl%nstates
      write(key,'(F9.6,F9.6)') traj%coeff_MCH_s(i)
      string=trim(string)//trim(key)
    enddo
    if (ctrl%laser==2) then
      do i=1,min(40,ctrl%nsteps*ctrl%nsubsteps+1)
        write(key,'(6(F9.6))') (ctrl%laserfield_td(i,j),j=1,3)
        string=trim(string)//trim(key)
      enddo
    endif
!     write(*,*) trim(string)
    temp=djb_hash(trim(string))
!     write(*,*) n
    hash_input=temp

  endfunction

! =================================================================== !

endmodule input
