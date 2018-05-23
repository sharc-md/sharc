!******************************************
!
!    SHARC Program Suite
!
!    Copyright (c) 2018 University of Vienna
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
  subroutine read_input(traj, ctrl)
  use definitions
  use input_list
  use misc
  use output
  use restart
  use string
  implicit none

  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl
  character*255 :: filename
  character*8000 :: geomfilename, line
  character*8000, allocatable :: values(:)
  integer :: narg, io, nlines, selg, selt
  integer :: i,j,k,n
  integer :: imult,ims
  real*8 :: a,b,tmax
  character*24 :: ctime, date
  integer :: idate,time
  character*8000 :: string1
  character*8000, allocatable :: string2(:)

  ! get the input filename from the command line argument

    ! get number of command line arguments
    narg=iargc()
    ! input filename must be present as argument
    if (narg==0) then
      write(0,*) 'Usage: sharc <inputfile>'
      stop 1
    endif
    call getarg(1,filename)

  ! =====================================================
  if ( (trim(filename)=='-v').or.(trim(filename)=='--version').or.(trim(filename)=='--info') ) then
    call write_license(0,version)
    stop
  endif

  ! =====================================================
  ! open the input file
  open(u_i_input,file=filename, status='old', action='read', iostat=io)
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
      close(u_i_input)
      call read_restart(u_resc,u_rest,ctrl,traj)

      ! if explicit laser field is used, the simulation time cannot be changed
      ! otherwise, the simulation time is read from the input file instead of the restart file
      if (ctrl%laser/=2) then
        line=get_value_from_key('nsteps',io)
        if (io==0) then
          read(line,*) ctrl%nsteps
        else
          line=get_value_from_key('tmax',io)
          if (io==0) then
            read(line,*) tmax
            ctrl%nsteps=int(tmax/ctrl%dtstep/au2fs+0.01d0)
          endif
        endif
        if (printlevel>0) then
          write(u_log,*) '============================================================='
          write(u_log,*) '                       Simulation Time'
          write(u_log,*) '============================================================='
          write(u_log,'(a,1x,i6,1x,a,1x,f6.3,1x,a)') 'Found nsteps=',ctrl%nsteps,'and stepsize=',ctrl%dtstep*au2fs,'fs.'
          if (printlevel>1) then
            write(u_log,'(a,1x,f9.3,1x,a)') 'This makes a total simulation time of ',ctrl%dtstep*ctrl%nsteps*au2fs,'fs.'
            write(u_log,'(a,1x,f7.4,1x,a)') 'The electronic wavefunction will be propagated using a ',&
            &ctrl%dtstep/ctrl%nsubsteps*au2fs, 'fs step.'
          endif
          write(u_log,*)
        endif
      endif

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

    ! stepsize is dt, default is 1 fs
    line=get_value_from_key('stepsize',io)
    if (io==0) then
      read(line,*) ctrl%dtstep
    else
      ctrl%dtstep=1.d0
    endif

    ! number of substeps for electronic interpolation
    line=get_value_from_key('nsubsteps',io)
    if (io==0) then
      read(line,*) ctrl%nsubsteps
    else
      ctrl%nsubsteps=25
    endif

    ! nsteps is either taken from input 
    ! or calculated based on maximum simulation time
    ! detault is 3 steps
    line=get_value_from_key('nsteps',io)
    if (io==0) then
      read(line,*) ctrl%nsteps
    else
      line=get_value_from_key('tmax',io)
      if (io==0) then
        read(line,*) tmax
        ctrl%nsteps=int(tmax/ctrl%dtstep)
      else
        ctrl%nsteps=3
      endif
    endif

    if (printlevel>0) then
      write(u_log,*) '============================================================='
      write(u_log,*) '                       Simulation Time'
      write(u_log,*) '============================================================='
      write(u_log,'(a,1x,i6,1x,a,1x,f6.3,1x,a)') 'Found nsteps=',ctrl%nsteps,'and stepsize=',ctrl%dtstep,'fs.'
      if (printlevel>1) then
        write(u_log,'(a,1x,f9.3,1x,a)') 'This makes a total simulation time of ',ctrl%dtstep*ctrl%nsteps,'fs.'
        write(u_log,'(a,1x,f7.4,1x,a)') 'The electronic wavefunction will be propagated using a ',&
        &ctrl%dtstep/ctrl%nsubsteps, 'fs step.'
      endif
      write(u_log,*)
    endif

    ! kill after n timesteps in gs ============================================

    line=get_value_from_key('killafter',io)
    if (io==0) then
      read(line,*) tmax
      if (tmax<=0) then
        ctrl%killafter=-1
        traj%steps_in_gs=-123
      else
        ctrl%killafter=int(anint(tmax/ctrl%dtstep))
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
        case default
          write(0,*) 'Unknown keyword ',trim(line),' to "coupling"!'
          stop 1
      endselect
    else
      ctrl%coupling=2
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
    ctrl%calc_phases=1
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
      ctrl%gradcorrect=1
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
        case ('parallel_nac')
          ctrl%ekincorrect=2
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
        case ('parallel_nac') 
          ctrl%reflect_frustrated=2
        case default
          write(0,*) 'Unknown keyword ',trim(line),' to "reflect_frustrated"!'
          stop 1
      endselect
    else
      ctrl%reflect_frustrated=0
    endif

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

    ! calculate the active coupling quantity (nacdr, nacdt, overlaps)
    select case (ctrl%coupling)
      case (0)  ! ddt
        ctrl%calc_nacdt=1
      case (1)  ! ddr
        ctrl%calc_nacdr=0
      case (2)  ! overlap
        ctrl%calc_overlap=1
      case default
        write(0,*) 'Internal error 1!'
        stop 1
    endselect

!     if (ctrl%coupling==1) then   ! having ddr couplings
!       ctrl%gradcorrect=1         ! they can be used for gradient correction
!     endif

    if (ctrl%surf/=0) then   ! not doing SHARC dynamics
      if (ctrl%gradcorrect==1) then
        write(0,*) 'Info: keyword "gradcorrect" has no effect in FISH dynamics.'
      endif
      ctrl%gradcorrect=0     ! no need for transformed gradients
    endif

    if (ctrl%gradcorrect==1) then ! for gradcorrect
      ctrl%calc_nacdr=0           ! we need nacdr
    endif
    if (ctrl%ekincorrect==2) then ! for kinetic energy correction parallel to nac, we need nacdr
      ctrl%calc_nacdr=0
!       ctrl%gradcorrect=1
    endif
    if (ctrl%reflect_frustrated==2) then ! for reflection parallel to nac, we need nacdr
      ctrl%calc_nacdr=0
    endif
! !     if (ctrl%decoherence==2) then ! for A-FSSH we need nacdr
! !       ctrl%calc_nacdr=0
! !     endif

    if (ctrl%surf==1) then   ! doing FISH
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

    line=get_value_from_key('select_directly',io)       ! do not do a second interface call
    if (io==0) then
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
        select case (ctrl%surf)
          case (0)
            write(u_log,'(a)') 'Doing SHARC dynamics.'
          case (1)
            write(u_log,'(a)') 'Doing FISH dynamics.'
        endselect
        select case (ctrl%coupling)
          case (0)
            write(u_log,'(a)') 'Using TIME DERIVATIVES <a|d/dt|b> for wavefunction propagation.'
          case (1)
            write(u_log,'(a)') 'Using SPATIAL DERIVATIVES <a|d/dR|b> for wavefunction propagation.'
          case (2)
            write(u_log,'(a)') 'Doing LOCAL DIABATISATION propagation.'
        endselect
        if (ctrl%gradcorrect==1) then
          write(u_log,'(a)') 'Including non-adiabatic coupling vectors in the gradient transformation.'
        endif
        select case (ctrl%ekincorrect)
          case (0)
            write(u_log,'(a)') 'No correction to the kinetic energy after surface hop.'
          case (1)
            write(u_log,'(a)') 'Correction to the kinetic energy after surface hop parallel to velocity vector.'
          case (2)
            write(u_log,'(a)') 'Correction to the kinetic energy after surface hop parallel to non-adiabatic coupling vector.'
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
          write(u_log,'(a)') 'Calculating Spin-Orbit couplings.'
        else
          write(u_log,'(a)') 'Not calculating Spin-Orbit couplings.'
        endif
      endif

      if (ctrl%calc_soc/=1) then
        n=0
        do i=1,ctrl%maxmult
          if (ctrl%nstates_m(i)>0) n=n+1
        enddo
        if (n>1) then
          write(u_log,'(a)') 'Warning: More than one multiplicity, but Spin-Orbit couplings are disabled.'
        endif
      endif
      write(u_log,*)

!       if (ctrl%write_soc==0) then
!         write(u_log,'(a)') 'Not writing Spin-Orbit couplings.'
!       else
!         write(u_log,'(a)') 'Writing Spin-Orbit couplings.'
!       endif
    if (printlevel>1) then
      if (ctrl%write_overlap==0) then
        write(u_log,'(a)') 'Not writing overlap matrix.'
      else
        write(u_log,'(a)') 'Writing overlap matrix.'
      endif
      if (ctrl%write_grad==0) then
        write(u_log,'(a)') 'Not writing gradients.'
      else
        write(u_log,'(a)') 'Writing gradients.'
      endif
      if (ctrl%write_NACdr==0) then
        write(u_log,'(a)') 'Not writing nonadiabatic couplings.'
      else
        write(u_log,'(a)') 'Writing nonadiabatic couplings.'
      endif
      if (ctrl%write_property1d==0) then
        write(u_log,'(a)') 'Not writing property vectors.'
      else
        write(u_log,'(a)') 'Writing property vectors.'
      endif
      if (ctrl%write_property2d==0) then
        write(u_log,'(a)') 'Not writing property matrices.'
      else
        write(u_log,'(a)') 'Writing property matrices.'
      endif
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

    if (printlevel>1) then
      write(u_log,'(a,1x,f6.3)') 'Scaling factor to the hamiltonian and gradients is',ctrl%scalingfactor
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
      ctrl%decoherence=1
      ctrl%decoherence_alpha=0.1d0
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

    if (printlevel>1) then
      if (ctrl%decoherence==1) then
        write(u_log,'(a)') 'Decoherence is 1 (EDC by Granucci, Persico, Zoccante)'
        write(u_log,'(a,1x,f6.3,1x,a)') 'Decoherence constant is',ctrl%decoherence_alpha,'Hartree'
        write(u_log,*)
      elseif (ctrl%decoherence==2) then
        write(u_log,'(a)') 'Decoherence is 2 (A-FSSH by Jain, Alguire, Subotnik)'
        write(u_log,*)
      else
        write(u_log,'(a)') 'Decoherence is OFF'
        write(u_log,*)
      endif
    endif

    ! phase tracking algorithm, actually not necessary in release SHARC version
    ! (see IJQC paper)
    line=get_value_from_key('notrack_phase',io)
    if (io==0) then
      ctrl%track_phase=0
      if (printlevel>1) then
        write(u_log,'(a)') 'Phase tracking is OFF'
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
            write(u_log,*)
          endif
        case ('sharc') 
          ctrl%hopping_procedure=1
          if (printlevel>1) then
            write(u_log,'(a)') 'Surface Hopping procedure is SHARC (Mai, Marquetand, Gonzalez)'
            write(u_log,*)
          endif
        case ('gfsh') 
          ctrl%hopping_procedure=2
          if (printlevel>1) then
            write(u_log,'(a)') 'Surface Hopping procedure is GFSH (Wang, Trivedi, Prezhdo)'
            write(u_log,*)
          endif
        case ('off') 
          ctrl%hopping_procedure=0
          ctrl%ekincorrect=0
          if (printlevel>1) then
            write(u_log,'(a)') 'Surface Hopping is OFF (will stay in initial diagonal state)'
            write(u_log,*)
          endif
        case default
          write(0,*) 'Unknown keyword ',trim(line),' to "hopping_procedure"!'
          stop 1
      endselect
    endif

    ! TODO: could delete this keyword
    line=get_value_from_key('no_hops',io)
    if (io==0) then
      ctrl%hopping_procedure=0  ! negate hopping
      ctrl%ekincorrect=0        ! negate kinetic energy adjustment at a negated hopping
      if (printlevel>1) then
        write(u_log,'(a)') 'Surface Hopping is OFF (will stay in initial diagonal state)'
        write(u_log,*)
      endif
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
    if (ctrl%ekincorrect==2) then
      write(u_log,*) 'HINT: "ekincorrect parallel_nac" ignores "atommask".'
    endif
    if (ctrl%decoherence==2) then
      write(u_log,*) 'HINT: "decoherence_scheme afssh" ignores "atommask".'
    endif
    if (ctrl%reflect_frustrated==2) then
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
  if (printlevel>2) then
    write(u_log,*) 'Atom mask:'
    do i=1,ctrl%natom
      write(u_log,*) i,traj%element_a(i),ctrl%atommask_a(i)
    enddo
    if (ctrl%ekincorrect==1) then
      write(u_log,*) 'Atom mask will be used for kinetic energy adjustment.'
    endif
    if (ctrl%decoherence==1) then
      write(u_log,*) 'Atom mask will be used for EDC decoherence.'
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
      case (1)  !FISH
        if (printlevel>1) write(u_log,'(A)') 'Setting state and coefficients for FISH dynamics.'
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
    ctrl%dtstep=ctrl%dtstep/au2fs
    ctrl%eselect_grad=ctrl%eselect_grad/au2eV
    ctrl%eselect_nac=ctrl%eselect_nac/au2eV
    traj%mass_a=traj%mass_a/au2u
    ! geometry is read in bohrs, as specified in the COLUMBUS geom format
    ! velocity is read in as bohrs/atu
    ! laser field must be in a.u.

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

endmodule
