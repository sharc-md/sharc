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

!> # Program DATA_EXTRACTOR.X
!> \authors Sebastian Mai, Philipp Marquetand
!> \date 22.10.2018
!>
!> This program reads the output.dat file of a trajectory and calculates various
!> properties per timestep, which are printed in plottable tables.
!>
!> Usage: '$SHARC/data_extractor.x <output.dat>'
!>
!> If a file "Reference/QM.out" exists, the overlap matrix from this file is read 
!> and used as reference overlap for along-trajectory-diabatization.
!>
!> Output files:
!> - energy.out 
!> - fosc.out
!> - fosc_act.out
!> - coeff_diab.out
!> - coeff_MCH.out
!> - coeff_diag.out
!> - spin.out
!> - prob.out
!> - expec.out
!> - expec_MCH.out
!> - dip_mom_diag.out
!> - dip_mom_proj.out
!> 
!> Additionally, the output file <input.file>.ext contains build infos of the data_extractor program
program data_extractor
  use matrix
  use definitions, only: au2a, au2fs, au2u, au2rcm, au2eV, au2debye
  use qm_out
  use string
  use input_list
  implicit none

  !> # Parameters: unit numbers for all files
  integer, parameter :: u_dat=11          !< unit for output.dat file
  integer, parameter :: u_ener=21         !< energy.out
  integer, parameter :: u_dm=22           !< fosc.out
  integer, parameter :: u_spin=23         !< spin.out
  integer, parameter :: u_coefd=24        !< coeff_diag.out
  integer, parameter :: u_coefm=25        !< coeff_MCH.out
  integer, parameter :: u_prob=26         !< prob.out
  integer, parameter :: u_expec=27        !< expec.out
  integer, parameter :: u_coefdiab=28     !< coeff_diab.out
  integer, parameter :: u_expec_mch=29    !< expec_MCH.out
  integer, parameter :: u_fosc_act=30     !< fosc_act.out
  integer, parameter :: u_ref=31          !< Reference/QM.out
  integer, parameter :: u_info=42         !< output.dat.ext
  integer, parameter :: u_ion_diag=51     !< ion_diag.out
  integer, parameter :: u_ion_mch=52      !< ion_mch.out
  integer, parameter :: u_dm_diag=53      !< dip_mom_diag.out
  integer, parameter :: u_dm_proj=55      !< dip_mom_proj.out
  integer, parameter :: u_proj=56         !< projections.inp

  integer, parameter :: u_classd=61        !< JCP 139, 211101 (2013), Method 1 (diag)
  integer, parameter :: u_classm=62        !< JCP 139, 211101 (2013), Method 1 (MCH)
  integer, parameter :: u_classdiab=63     !< JCP 139, 211101 (2013), Method 1 (diabatic)
  integer, parameter :: u_cmixd=64         !< JCP 139, 211101 (2013), Method 3 (diag)
  integer, parameter :: u_cmixm=65         !< JCP 139, 211101 (2013), Method 3 (MCH)
  integer, parameter :: u_cmixdiab=66      !< JCP 139, 211101 (2013), Method 3 (diabatic)


  !> # Information which is constant throughout all timesteps
  integer :: nstates                      !< total number of states
  integer :: nargs                        !< number of command line arguments
  integer :: maxmult                      !< maximum multiplicity
  integer :: natom                        !< number of atoms
  integer, allocatable :: nstates_m(:)    !< number of states per multiplicity
  real*8 :: dtstep                        !< nuclear timestep
  real*8 :: ezero                         !< reference energy
  integer :: have_overlap                 !< whether overlap matrices are in the dat file (0=no, 1=yes)
  integer :: have_grad                    !< whether gradients are in the dat file (0=no, 1=yes)
  integer :: have_NAC                     !< whether nonadiabatic couplings are in the dat file (0=no, 1=yes)
  integer :: have_property1d                !< whether property vectors are in the dat file (0=no, 1=yes)
  integer :: have_property2d                !< whether property matrices are in the dat file (0=no, 1=yes)
  integer :: n_property1d                !< 
  integer :: n_property2d                !< 
  integer :: laser                        !< whether a laser field is in the dat file (0=no, 1=, 2=yes)
  integer :: nsteps                       !< number of timesteps from dat file (needed to read the laser field)
  integer :: nsubsteps                    !< number of substeps (needed to read the laser field)
  integer :: nprojections                    !< number of vectors to project dipole moment onto

  !> # Information which is updated per time step
  !> Most of these are equivalent to their definition in definitions.f90
  integer :: step
  complex*16, allocatable :: H_MCH_ss(:,:),U_ss(:,:),DM_mch_ssd(:,:,:),DM_diag_ssd(:,:,:)
  complex*16, allocatable :: Prop2d_xss(:,:,:)
  real*8, allocatable     :: Prop1d_ys(:,:)
  complex*16, allocatable :: coeff_diag_s(:),overlaps_ss(:,:), ref_ovl_ss(:,:)
  real*8, allocatable :: geom_ad(:,:), veloc_ad(:,:)
  real*8 , allocatable :: hopprob_s(:)
  real*8 :: Ekin, Epot, randnum
  integer :: state_diag, state_MCH, runtime

  !> # Information which is calculated at each time step
  complex*16, allocatable :: H_diag_ss(:,:)       !< diagonal Hamiltonian (including laser field)
  complex*16, allocatable :: A_ss(:,:)            !< temporary matrix
  complex*16, allocatable :: coeff_MCH_s(:)       !< MCH coefficient vector
  complex*16, allocatable :: laser_td(:,:)        !< laser field for all timesteps
  complex*16, allocatable :: coeff_diab_s(:)      !< diabatic coefficient vector
  real*8,allocatable :: expec_pop(:)                !< spin expectation value per state
  real*8,allocatable :: expec_s(:)                !< spin expectation value per state
  real*8,allocatable :: expec_dm(:)               !< oscillator strength per state
  real*8,allocatable :: expec_dm_mch(:)           !< oscillator strength per state in MCH basis
  real*8,allocatable :: expec_ion_diag(:)           !< oscillator strength per state in MCH basis
  real*8,allocatable :: expec_ion_mch(:)           !< oscillator strength per state in MCH basis
  real*8,allocatable :: expec_dm_act(:)           !< oscillator strength per state with active state as source state
  real*8,allocatable :: dm_proj_sp(:,:)           !< dipole moment per state projected onto vector yet to be defined
  real*8,allocatable :: proj_vec_pd(:,:)           !< vectors to project dipole moment onto
  real*8,allocatable :: norm_proj_vec_p(:)          !< norm of vectors to project dipole moment onto
  real*8,allocatable :: proj_point(:,:)          !< norm of vectors to project dipole moment onto
  real*8,allocatable :: spin0_s(:)                !< spin value per MCH state (initialized in the beginning)
  real*8,allocatable :: grad_mch_sad(:,:,:)       !< gradient per MCH state per atom per xyz
  real*8,allocatable :: NAC_ssad(:,:,:,:)         !< nonadiabatic coupling per element (MCH state, MCH state) per atom per xyz
  real*8,allocatable :: distance_scaling(:,:)         !< scaling factor to define point between atom positions
  real*8 :: sumc                                  !< sum of coefficients
  integer,allocatable :: proj_atomindex(:,:,:)         !< atomindex to be read from projections.inp

  ! helper
  character*8000 :: filename, string1, string3, line
  character*8000, allocatable :: args(:)
  character*8000, allocatable :: values(:)
  character*21 :: string2
  integer :: i, io, idir,istate,jstate,kstate,imult,ims,j,n,iproj
  logical :: exists
  logical :: is_integer
  logical :: write_energy
  logical :: write_dip
  logical :: write_spin
  logical :: write_coeffdiag
  logical :: write_coeffmch
  logical :: write_prob
  logical :: write_expec
  logical :: write_expecmch
  logical :: write_coeffdiab
  logical :: write_dipact
  logical :: write_dm_diag
  logical :: write_dm_proj
  logical :: write_iondiag
  logical :: write_ionmch
  logical :: skip_geom_vel_grad_nac
  logical :: anyoptions
  integer :: skipthese
  integer :: ncolumns

  ! build_info.inc is written by the Makefile and contains the 
  ! date and host, when/where SHARC was built
  include 'build_info.inc'

  ! =============================================================================================
  !                                Command line argument processing
  ! =============================================================================================

  ! print usage if no arguments are given
  nargs=iargc()
  if (nargs<1) then
    call print_usage(0)
    stop 
  endif

  ! get the command line arguments
  allocate(args(1:nargs))
  do i=1,nargs
    call getarg(i,args(i))
    args(i)=trim(adjustl(args(i)))
  end do

  ! defaults for writing options
  write_energy    = .false.
  write_dip       = .false.
  write_spin      = .false.
  write_coeffdiag = .false.
  write_coeffmch  = .false.
  write_prob      = .false.
  write_expec     = .false.
  write_expecmch  = .false.
  write_coeffdiab = .false.
  write_dipact    = .false.
  write_dm_diag   = .false.
  write_dm_proj   = .false.
  write_iondiag   = .false.
  write_ionmch    = .false.
  skip_geom_vel_grad_nac = .false.

  ! read command line arguments
  skipthese=0
  anyoptions=.false.
  filename=''

  do i=1,nargs
    ! skip arguments to options
    if (skipthese>0) then
      skipthese=skipthese-1
      cycle
    endif
    ! make option flags lowercase
    if (args(i)(1:1)=='-') then
      call lowercase(args(i))
      anyoptions=.true.
    else
      ! we have an argument which is not an option, hence it is the filename
      filename = trim(adjustl(args(i)))
      cycle
    endif

    ! process the options
    if (trim(args(i)) == '-f') then
      if ((i+1).GT.nargs) stop '-f must be followed by a filename'
      filename = trim(adjustl(args(i+1)))
      skipthese=1
    elseif (trim(args(i)) == "-e") then
      write_energy = .true.
    elseif (trim(args(i)) == "-d") then
      write_dip = .true.
    elseif (trim(args(i)) == "-sp") then
      write_spin = .true.
    elseif (trim(args(i)) == "-cd") then
      write_coeffdiag = .true.
    elseif (trim(args(i)) == "-cm") then
      write_coeffmch = .true.
    elseif (trim(args(i)) == "-p") then
      write_prob = .true.
    elseif (trim(args(i)) == "-x") then
      write_expec = .true.
    elseif (trim(args(i)) == "-xm") then
      write_expecmch = .true.
    elseif (trim(args(i)) == "-cb") then
      write_coeffdiab = .true.
    elseif (trim(args(i)) == "-da") then
      write_dipact = .true.
    elseif (trim(args(i)) == "-dd") then
      write_dm_diag = .true.
    elseif (trim(args(i)) == "-dp") then
      write_dm_proj = .true.
    elseif (trim(args(i)) == "-id") then
      write_iondiag = .true.
    elseif (trim(args(i)) == "-im") then
      write_ionmch = .true.
    elseif (trim(args(i)) == "-sk") then
      skip_geom_vel_grad_nac = .true.
    ! all flags true
    elseif (trim(args(i)) == "-xl") then
      write_energy = .true.
      write_dip = .true.
      write_spin = .true.
      write_coeffdiag = .true.
      write_coeffmch = .true.
      write_prob = .true.
      write_expec = .true.
      write_expecmch = .true.
      write_coeffdiab = .true.
      write_dipact = .true.
      write_dm_diag = .true.
      write_dm_proj = .true.
      write_iondiag = .true.
      write_ionmch = .true.
    ! large set of flags true
    elseif (trim(args(i)) == "-l") then
      write_energy = .true.
      write_dip = .true.
      write_spin = .true.
      write_coeffdiag = .true.
      write_coeffmch = .true.
      write_prob = .true.
      write_expec = .true.
      write_expecmch = .true.
      write_coeffdiab = .true.
      write_dipact = .true.
      write_iondiag = .true.
      write_ionmch = .true.
      skip_geom_vel_grad_nac = .true.
    ! small set of flags true
    elseif (trim(args(i)) == "-s") then
      write_energy = .true.
      write_dip = .true.
      write_spin = .true.
      write_coeffdiag = .true.
      write_coeffmch = .true.
      write_prob = .true.
      write_expec = .true.
      write_expecmch = .true.
      write_coeffdiab = .true.
      write_dipact = .true.
      skip_geom_vel_grad_nac = .true.
    ! very small set of flags true
    elseif (trim(args(i)) == "-xs") then
      write_energy = .true.
      write_dip = .true.
      write_coeffdiag = .true.
      write_coeffmch = .true.
      write_prob = .true.
      write_expec = .true.
      skip_geom_vel_grad_nac = .true.
    elseif (trim(args(i)) == "-xyz") then
!       write_options%write_geometry = .true.
      write(0,*) 'Ignoring -xyz option in ASCII mode.'
    elseif (trim(args(i)) == "-h") then
      call print_usage(0)
      stop 
    else
      write(0,*) 'Cannot understand argument: ',trim(args(i))
      !stop 'Input error'
    endif
  enddo

  if (.not.anyoptions) then
    ! defaults for writing options
    write_energy    = .true.
    write_dip       = .true.
    write_spin      = .true.
    write_coeffdiag = .true.
    write_coeffmch  = .true.
    write_prob      = .true.
    write_expec     = .true.
    write_expecmch  = .true.
    write_coeffdiab = .true.
    write_dipact    = .true.
    write_dm_diag   = .false.
    write_dm_proj   = .false.
    write_iondiag   = .false.
    write_ionmch    = .false.
    skip_geom_vel_grad_nac = .true.
  endif

  if (skip_geom_vel_grad_nac .and. write_dm_proj) then
    write_dm_proj   = .false.
    write(*,*) 'Skip reading of geometries, writing of dipole moment projections switched of'
  endif
  
  deallocate(args)
  
  ! =============================================================================================
  !                                Open projections file 
  ! =============================================================================================

  
  ! open the projections file
  if (write_dm_proj) then
    open(unit=u_proj, file='projections.inp', status='old', action='read', iostat=io)
    if (io/=0) then
      write(*,*) 'Warning!'
      write(*,*) 'File projections.inp not found'
      write(*,*) 'Switching off writing of dipole moment projections'
      write(*,*) 
      write_dm_proj   = .false.
    else
      read(u_proj,*,iostat=io) nprojections
      if (io<0) then
        write(*,*) 'First entry of projections.inp should be number of projections to be carried out.'
        write(*,*)
        stop
      endif
      
      allocate(proj_atomindex(nprojections,2,2))
      allocate(proj_point(2,2))
      allocate(distance_scaling(nprojections,2))
      proj_atomindex=-1
      distance_scaling=1.
      do iproj=1,nprojections
        read(u_proj,*) ! empty line / comment
        ! 1st row
        read(u_proj,'(A)', iostat=io) line
        if (io<0) then
          write(*,*) 'Error in format of projections.inp'
          write(*,*) 'Format example for 1 projection (first line)'
          write(*,*) 'on vector from coordinates of atom with index 2'
          write(*,*) 'to point with coordinates half way between atoms with indices 3 and 4:'
          write(*,*) '1'
          write(*,*) 'comment'
          write(*,*) '2'
          write(*,*) '3 4 0.5'
          stop
        endif
        ncolumns = count([len_trim(line) > 0,(line(i:i)/=' '.and.line(i:i)/=','&
         &  .and.line(i+1:i+1)==' '.or.line(i+1:i+1)==',', i=1,len_trim(line)-1)])
        if (ncolumns == 1) then
          read(line,*) proj_atomindex(iproj,1,1)
        elseif (ncolumns == 3) then
          read(line,*) proj_atomindex(iproj,1,1), proj_atomindex(iproj,1,2), distance_scaling(iproj,1)
        else
          write(*,*) 'projections.inp malformatted'
        endif
        ! 2nd row
        read(u_proj,'(A)', iostat=io) line
        if (io<0) then
          write(*,*) 'Error in format of projections.inp'
          write(*,*) 'Format example for 1 projection (first line)'
          write(*,*) 'on vector from coordinates of atom with index 2'
          write(*,*) 'to point with coordinates half way between atoms with indices 3 and 4:'
          write(*,*) '1'
          write(*,*) 'comment'
          write(*,*) '2'
          write(*,*) '3 4 0.5'
          stop
        endif
        ncolumns = count([len_trim(line) > 0,(line(i:i)/=' '.and.line(i:i)/=','&
         &  .and.line(i+1:i+1)==' '.or.line(i+1:i+1)==',', i=1,len_trim(line)-1)])
        if (ncolumns == 1) then
          read(line,*) proj_atomindex(iproj,2,1)
        elseif (ncolumns == 3) then
          read(line,*) proj_atomindex(iproj,2,1), proj_atomindex(iproj,2,2), distance_scaling(iproj,2)
        else
          write(*,*) 'projections.inp malformatted'
          write(*,*) 'Format example for 1 projection (first line)'
          write(*,*) 'on vector from coordinates of atom with index 2'
          write(*,*) 'to point with coordinates half way between atoms with indices 3 and 4:'
          write(*,*) '1'
          write(*,*) 'comment'
          write(*,*) '2'
          write(*,*) '3 4 0.5'
          stop
        endif
      enddo
    endif ! io/=0
  endif
  
  ! =============================================================================================
  !                                Open dat file and write build info
  ! =============================================================================================

  ! open the dat file
  if (trim(filename)=='') stop 'No filename given!'
  open(unit=u_dat, file=filename, status='old', action='read', iostat=io)
  if (io/=0) then
    write(*,*) 'File ',trim(filename),' not found'
    stop
  endif

  ! write build infos to file, filename is automatically generated from dat file + '.ext'
  string1=trim(filename)//'.ext'
  open(u_info,file=string1,status='replace',action='write')
  write(u_info,*) 'BUILD INFORMATION:'
  write(u_info,*) 'Build date: ',trim(build_date)
  write(u_info,*) 'Build host: ',trim(build_host)
  write(u_info,*) 'Build directory: ',trim(build_dir)
  write(u_info,*) 'Compiler: ',trim(build_compiler)
  close(u_info)

  ! =============================================================================================
  !                                Read dat file header
  ! =============================================================================================

  ! Check whether first entry is an integer, which indicates that we have the old format starting with "<integer> ! maxmult"
  read(u_dat,'(I8,A)',iostat=io) maxmult, string3
  if (io==0) then
    is_integer = .true.
  else
    is_integer = .false.
  endif
  rewind u_dat


  ! set default switches for different properties
  have_NAC=0
  have_grad=0
  have_overlap=0
  have_property1d=0
  have_property2d=0


  if (is_integer) then
    write(*,*) 'Found SHARC v1.0 format'
    ! old format has property matrix by default
    have_property2d=1
    n_property1d=1
    n_property2d=1
    ! determine number of states per mult and total number of states
    read(u_dat,*) maxmult
    allocate( nstates_m(maxmult) )
    read(u_dat,*) (nstates_m(i),i=1,maxmult)
    nstates=0
    do i=1,maxmult
      nstates=nstates+i*nstates_m(i)
    enddo
    write(*,*) 'Found nstates=',nstates
    read(u_dat,*) natom

    ! allocate everything
    allocate( H_MCH_ss(nstates,nstates), H_diag_ss(nstates,nstates) )
    allocate( U_ss(nstates,nstates) )
    allocate( Prop2d_xss(n_property2d,nstates,nstates),Prop1d_ys(n_property1d,nstates) )
    allocate( overlaps_ss(nstates,nstates), ref_ovl_ss(nstates,nstates) )
    allocate( DM_mch_ssd(nstates,nstates,3) )
    allocate( DM_diag_ssd(nstates,nstates,3) )
    allocate( dm_proj_sp(nstates,nprojections) )
    allocate( proj_vec_pd(nprojections,3) )
    allocate( norm_proj_vec_p(nprojections) )
    allocate( coeff_diag_s(nstates), coeff_MCH_s(nstates), coeff_diab_s(nstates) )
    allocate( hopprob_s(nstates) )
    allocate( A_ss(nstates,nstates) )
    allocate( expec_s(nstates),expec_dm(nstates),expec_dm_mch(nstates),expec_dm_act(nstates) )
    allocate( expec_ion_diag(nstates),expec_ion_mch(nstates),expec_pop(nstates) )
    allocate( spin0_s(nstates) )
    if (.not. skip_geom_vel_grad_nac) then
      allocate( geom_ad(natom,3), veloc_ad(natom,3) )
    endif
    call allocate_lapack(nstates)
    overlaps_ss=dcmplx(0.d0,0.d0)

    ! obtain the timestep, ezero and have_overlap 
    ! most of the information from dat file header is needed in order to interpret the matrices per timestep
    read(u_dat,*) dtstep
    dtstep=dtstep*au2fs
    read(u_dat,*) ezero
    read(u_dat,*) have_overlap
    read(u_dat,*) laser
    read(u_dat,*) nsteps
    read(u_dat,*) nsubsteps
  
  else
    write(*,*) 'Found SHARC v2.0 format'
    ! else we have the format of SHARC 2.0, which is list-based
    ! =====================================================
    call read_input_list_from_file(u_dat)

    ! look up nstates keyword
    line=get_value_from_key('nstates_m',io)
    if (io==0) then
      ! value needs to be split into values (each one is a string)
      call split(line,' ',values,n)
      ! n is number of multiplicities = maxmult
      maxmult=n
      allocate(nstates_m(n))
      nstates=0
      ! read number of states per multiplicity
      do i=1,n
        read(values(i),*) nstates_m(i)
        nstates=nstates+nstates_m(i)*i
      enddo
      write(*,*) 'Found nstates=',nstates
      ! values is not needed anymore
      deallocate(values)
    else
      ! no nstates keywords => adiabatic dynamics, 1 singlet state
      maxmult=1
      allocate(nstates_m(1))
      nstates_m(1)=1
      nstates=1
    endif

    ! look up nstates keyword
    line=get_value_from_key('natom',io)
    if (io==0) then
      read(line,*) natom
    else
      ! currently natom is needed
      stop 'Error! Number of atoms (keyword: natom) is required!'
    endif

    ! look up properties
    line=get_value_from_key('n_property1d',io)
    if (io==0) then
      read(line,*) n_property1d
    else
      n_property1d=1
    endif
    line=get_value_from_key('n_property2d',io)
    if (io==0) then
      read(line,*) n_property2d
    else
      n_property2d=1
    endif

    ! allocate everything
    allocate( H_MCH_ss(nstates,nstates), H_diag_ss(nstates,nstates) )
    allocate( U_ss(nstates,nstates) )
    allocate( Prop2d_xss(n_property2d,nstates,nstates), Prop1d_ys(n_property1d,nstates) )
    allocate( overlaps_ss(nstates,nstates), ref_ovl_ss(nstates,nstates) )
    allocate( DM_mch_ssd(nstates,nstates,3) )
    allocate( DM_diag_ssd(nstates,nstates,3) )
    allocate( dm_proj_sp(nstates,nprojections) )
    allocate( proj_vec_pd(nprojections,3) )
    allocate( norm_proj_vec_p(nprojections) )
    allocate( coeff_diag_s(nstates), coeff_MCH_s(nstates), coeff_diab_s(nstates) )
    allocate( hopprob_s(nstates) )
    allocate( A_ss(nstates,nstates) )
    allocate( expec_s(nstates),expec_dm(nstates),expec_dm_mch(nstates),expec_dm_act(nstates) )
    allocate( expec_ion_diag(nstates),expec_ion_mch(nstates),expec_pop(nstates) )
    allocate( spin0_s(nstates) )
    if (.not. skip_geom_vel_grad_nac) then
      allocate( geom_ad(natom,3), veloc_ad(natom,3) )
    endif
    call allocate_lapack(nstates)
    overlaps_ss=dcmplx(0.d0,0.d0)

    ! look up dtstep keyword
    line=get_value_from_key('dtstep',io)
    if (io==0) then
      read(line,*) dtstep
      dtstep=dtstep*au2fs
    else
      ! dtstep is needed
      stop 'Error! Time step (keyword: dtstep) is required!'
    endif

    ! look up have_overlap keyword
    line=get_value_from_key('write_overlap',io)
    if (io==0) then
      read(line,*) have_overlap
    endif
    if (have_overlap == 0 .and. write_coeffdiab) then
      write(*,*)  'Warning! Writing diabatic coefficients is impossible without overlaps present.'
      write(*,*)  'Unsetting flag -cb'
      write_coeffdiab = .false.
    endif

    ! look up have_grad keyword
    line=get_value_from_key('write_grad',io)
    if (io==0) then
      read(line,*) have_grad
    endif

    ! look up have_NAC keyword
    line=get_value_from_key('write_nacdr',io)
    if (io==0) then
      read(line,*) have_NAC
    endif

    ! look up have_property keyword
    line=get_value_from_key('write_property1d',io)
    if (io==0) then
      read(line,*) have_property1d
    endif
    line=get_value_from_key('write_property2d',io)
    if (io==0) then
      read(line,*) have_property2d
    elseif (io==-1) then
      line=get_value_from_key('write_property',io)      ! backwards compatibility
      if (io==0) read(line,*) have_property2d
      n_property2d=1
    endif
    if (have_property2d == 0 .and. (write_iondiag .or. write_ionmch)) then
      write(*,*)  'Warning! Writing ionization probabilities does not make sense if property matrix is not present.'
      write(*,*)  'Unsetting flags -im and -id'
      write_ionmch = .false.
      write_iondiag = .false.
    endif

    ! look up laser keyword
    line=get_value_from_key('laser',io)
    if (io==0) then
      read(line,*) laser
      if (laser==2) then
        line=get_value_from_key('nsteps',io)
        if (io==0) then
          read(line,*) nsteps
        else
          stop 'Error! If laser is present, keyword nsteps must be given.'
        endif
        line=get_value_from_key('nsubsteps',io)
        if (io==0) then
          read(line,*) nsubsteps
        else
          stop 'Error! If laser is present, keyword nsubsteps must be given.'
        endif
      endif
    else
      laser=0
    endif

    if (have_grad == 1) then
      if (.not. skip_geom_vel_grad_nac) then
        allocate( grad_mch_sad(nstates,natom,3) )
      endif
    endif

    if (have_NAC == 1) then
      if (.not. skip_geom_vel_grad_nac) then
        allocate( NAC_ssad(nstates,nstates,natom,3) )
      endif
    endif

    ! Now we skip over the header array data (atomic numbers, elements, masses)
    do i=1,3*(1+natom)
      read(u_dat,*) string1
    enddo

    ! if an explicit laser file is in the dat file, read it now
    ! laser field comes before the time step data
    if (laser==2) then
      allocate( laser_td(nsteps*nsubsteps+1,3) )
      call vec3read(nsteps*nsubsteps+1,laser_td,u_dat,string1)
    endif

    ! skip the "End of header array data" separator line
    read(u_dat,*)

  endif

  ! =============================================================================================
  !                                Create output directory and files
  ! =============================================================================================

  ! create output directory "output_data"
  ! inquire will not work with Intel compiler, so mkdir is always attempted
  inquire(file="output_data", exist=exists)
  if (.not.exists) then
    write(*,'(A)') ' Creating directory "output_data"'
    call system('mkdir output_data')
  else
    write(*,'(A)') ' Writing to directory "output_data"'
  endif



  ! open output files
  if (write_energy)    open(unit=u_ener, file='output_data/energy.out', status='replace', action='write')           ! -e
  if (write_dip)       open(unit=u_dm, file='output_data/fosc.out', status='replace', action='write')               ! -d
  if (write_spin)      open(unit=u_spin, file='output_data/spin.out', status='replace', action='write')             ! -s

  if (write_coeffdiag) open(unit=u_coefd, file='output_data/coeff_diag.out', status='replace', action='write')      ! -cd
  if (write_coeffmch)  open(unit=u_coefm, file='output_data/coeff_MCH.out', status='replace', action='write')       ! -cm
  if (write_coeffdiab) open(unit=u_coefdiab, file='output_data/coeff_diab.out', status='replace', action='write')   ! -cb

  if (write_coeffdiag) open(unit=u_classd, file='output_data/coeff_class_diag.out', status='replace', action='write')      ! -cd
  if (write_coeffmch)  open(unit=u_classm, file='output_data/coeff_class_MCH.out', status='replace', action='write')       ! -cm
  if (write_coeffdiab) open(unit=u_classdiab, file='output_data/coeff_class_diab.out', status='replace', action='write')   ! -cb

  if (write_coeffdiag) open(unit=u_cmixd, file='output_data/coeff_mixed_diag.out', status='replace', action='write')      ! -cd
  if (write_coeffmch)  open(unit=u_cmixm, file='output_data/coeff_mixed_MCH.out', status='replace', action='write')       ! -cm
  if (write_coeffdiab) open(unit=u_cmixdiab, file='output_data/coeff_mixed_diab.out', status='replace', action='write')   ! -cb

  if (write_prob)      open(unit=u_prob, file='output_data/prob.out', status='replace', action='write')             ! -p
  if (write_expec)     open(unit=u_expec, file='output_data/expec.out', status='replace', action='write')           ! -x
  if (write_expecmch)  open(unit=u_expec_mch, file='output_data/expec_MCH.out', status='replace', action='write')   ! -xm
  if (write_dipact)    open(unit=u_fosc_act, file='output_data/fosc_act.out', status='replace', action='write')     ! -da
  if (write_iondiag)   open(unit=u_ion_diag, file='output_data/ion_diag.out', status='replace', action='write')     ! -id
  if (write_ionmch)    open(unit=u_ion_mch, file='output_data/ion_mch.out', status='replace', action='write')       ! -im
  if (write_dm_diag)   open(unit=u_dm_diag, file='output_data/dip_mom_diag.out', status='replace', action='write')   ! -dd
  if (write_dm_proj)   open(unit=u_dm_proj, file='output_data/dip_mom_proj.out', status='replace', action='write')   ! -dd
                                                                                                                    ! -a



  ! write output file headers
  if (write_energy) then
    write(u_ener,'(A1,1X,1000(I20,1X))') '#',(i,i=1,nstates+4)
    write(u_ener,'(A1,1X,5(A20,1X))') '#','Time |','Ekin |','Epot |','Etot |','=== Energy ===>'
    write(u_ener,'(A1,1X,5(A20,1X))') '#','[fs] |','[eV] |','[eV] |','[eV] |','[eV] |'
  endif

  if (write_dip) then
    write(u_dm,'(A1,1X,1000(I20,1X))') '#',(i,i=1,nstates+2)
    write(u_dm,'(A1,1X,3(A20,1X))') '#','Time |','f_osc (state) |','=== f_osc ===>'
    write(u_dm,'(A1,1X,3(A20,1X))') '#','[fs] |','[] |','[] |'
  endif

  if (write_dipact) then
    write(u_fosc_act,'(A1,1X,1000(I20,1X))') '#',(i,i=1,2*nstates+1)
    write(string1, '(A1,1X,1(A20,1X))') '#','Time |'
    do i=1,nstates
      write(string2,'(1X,A8,I10,A2)') 'dE ',i,' |'
      string1=trim(string1)//string2
    enddo
    do i=1,nstates
      write(string2,'(1X,A6,I12,A2)') 'f_osc ',i,' |'
      string1=trim(string1)//string2
    enddo
    write(u_fosc_act,'(A)') trim(string1)
    write(string1, '(A1,1X,1(A20,1X))') '#','[fs] |'
    do i=1,nstates
      write(string2,'(1X,A20)') '[eV] |'
      string1=trim(string1)//string2
    enddo
    do i=1,nstates
      write(string2,'(1X,A20)') '[] |'
      string1=trim(string1)//string2
    enddo
    write(u_fosc_act,'(A)') trim(string1)
  endif
  
  if (write_iondiag) then
    write(u_ion_diag,'(A1,1X,1000(I20,1X))') '#',(i,i=1,2*nstates+2)
    write(string1, '(A1,1X,2(A20,1X))') '#','Time |','State (diag) |'
    do i=1,nstates
      write(string2,'(1X,A8,I10,A2)') 'dE ',i,' |'
      string1=trim(string1)//string2
    enddo
    do i=1,nstates
      write(string2,'(1X,A6,I12,A2)') 'f_osc ',i,' |'
      string1=trim(string1)//string2
    enddo
    write(u_ion_diag,'(A)') trim(string1)
    write(string1, '(A1,1X,2(A20,1X))') '#','[fs] |','|'
    do i=1,nstates
      write(string2,'(1X,A20)') '[eV] |'
      string1=trim(string1)//string2
    enddo
    do i=1,nstates
      write(string2,'(1X,A20)') '[] |'
      string1=trim(string1)//string2
    enddo
    write(u_ion_diag,'(A)') trim(string1)
  endif
  
  if (write_ionmch) then 
    write(u_ion_mch,'(A1,1X,1000(I20,1X))') '#',(i,i=1,2*nstates+2)
    write(string1, '(A1,1X,2(A20,1X))') '#','Time |','State (diag) |'
    do i=1,nstates
      write(string2,'(1X,A8,I10,A2)') 'dE ',i,' |'
      string1=trim(string1)//string2
    enddo
    do i=1,nstates
      write(string2,'(1X,A6,I12,A2)') 'f_osc ',i,' |'
      string1=trim(string1)//string2
    enddo
    write(u_ion_mch,'(A)') trim(string1)
    write(string1, '(A1,1X,2(A20,1X))') '#','[fs] |','|'
    do i=1,nstates
      write(string2,'(1X,A20)') '[eV] |'
      string1=trim(string1)//string2
    enddo
    do i=1,nstates
      write(string2,'(1X,A20)') '[] |'
      string1=trim(string1)//string2
    enddo
    write(u_ion_mch,'(A)') trim(string1)
  endif
  
  if (write_spin) then
    write(u_spin,'(A1,1X,1000(I20,1X))') '#',(i,i=1,nstates+2)
    write(u_spin,'(A1,1X,3(A20,1X))') '#','Time |','Spin_s |','=== Spins ===>'
    write(u_spin,'(A1,1X,3(A20,1X))') '#','[fs] |','[] |','[] |'
  endif
  
  if (write_coeffdiag) then
    write(u_coefd,'(A1,1X,1000(I20,1X))') '#',(i,i=1,2*nstates+2)
    write(u_coefd,'(A1,1X,3(A20,1X))') '#','Time |','Sum c**2 |','=== coeff_diag ===>'
    write(u_coefd,'(A1,1X,3(A20,1X))') '#','[fs] |','[] |','[] |'
    write(u_classd,'(A1,1X,1000(I20,1X))') '#',(i,i=1,2*nstates+2)
    write(u_classd,'(A1,1X,3(A20,1X))') '#','Time |','Sum c**2 |','=== coeff_diag ===>'
    write(u_classd,'(A1,1X,3(A20,1X))') '#','[fs] |','[] |','[] |'
    write(u_cmixd,'(A1,1X,1000(I20,1X))') '#',(i,i=1,2*nstates+2)
    write(u_cmixd,'(A1,1X,3(A20,1X))') '#','Time |','Sum c**2 |','=== coeff_diag ===>'
    write(u_cmixd,'(A1,1X,3(A20,1X))') '#','[fs] |','[] |','[] |'
  endif
  
  if (write_coeffmch) then
    write(u_coefm,'(A1,1X,1000(I20,1X))') '#',(i,i=1,2*nstates+2)
    write(u_coefm,'(A1,1X,3(A20,1X))') '#','Time |','Sum c**2 |','=== coeff_MCH ===>'
    write(u_coefm,'(A1,1X,3(A20,1X))') '#','[fs] |','[] |','[] |'
    write(u_classm,'(A1,1X,1000(I20,1X))') '#',(i,i=1,2*nstates+2)
    write(u_classm,'(A1,1X,3(A20,1X))') '#','Time |','Sum c**2 |','=== coeff_MCH ===>'
    write(u_classm,'(A1,1X,3(A20,1X))') '#','[fs] |','[] |','[] |'
    write(u_cmixm,'(A1,1X,1000(I20,1X))') '#',(i,i=1,2*nstates+2)
    write(u_cmixm,'(A1,1X,3(A20,1X))') '#','Time |','Sum c**2 |','=== coeff_MCH ===>'
    write(u_cmixm,'(A1,1X,3(A20,1X))') '#','[fs] |','[] |','[] |'
  endif
  
  if (write_coeffdiab) then
    write(u_coefdiab,'(A1,1X,1000(I20,1X))') '#',(i,i=1,2*nstates+2)
    write(u_coefdiab,'(A1,1X,3(A20,1X))') '#','Time |','Sum c**2 |','=== coeff_diab ===>'
    write(u_coefdiab,'(A1,1X,3(A20,1X))') '#','[fs] |','[] |','[] |'
    write(u_classdiab,'(A1,1X,1000(I20,1X))') '#',(i,i=1,2*nstates+2)
    write(u_classdiab,'(A1,1X,3(A20,1X))') '#','Time |','Sum c**2 |','=== coeff_diab ===>'
    write(u_classdiab,'(A1,1X,3(A20,1X))') '#','[fs] |','[] |','[] |'
    write(u_cmixdiab,'(A1,1X,1000(I20,1X))') '#',(i,i=1,2*nstates+2)
    write(u_cmixdiab,'(A1,1X,3(A20,1X))') '#','Time |','Sum c**2 |','=== coeff_diab ===>'
    write(u_cmixdiab,'(A1,1X,3(A20,1X))') '#','[fs] |','[] |','[] |'
  endif
  
  if (write_prob) then
    write(u_prob,'(A1,1X,1000(I20,1X))') '#',(i,i=1,nstates+2)
    write(u_prob,'(A1,1X,3(A20,1X))') '#','Time |','Random Number |','=== cumu Prob ===>'
    write(u_prob,'(A1,1X,3(A20,1X))') '#','[fs] |','[] |','[] |'
  endif
  
  if (write_expec .or. write_expecmch) then
    ! Strings for expec and expec_mch
    write(string1, '(A1,1X,4(A20,1X))') '#','Time |','Ekin |','Epot |','Etot |'
    do i=1,nstates
      write(string2,'(1X,A8,I10,A2)') 'Energy ',i,' |'
      string1=trim(string1)//string2
    enddo
    do i=1,nstates
      write(string2,'(1X,A5,I13,A2)') 'Spin ',i,' |'
      string1=trim(string1)//string2
    enddo
    do i=1,nstates
      write(string2,'(1X,A6,I12,A2)') 'f_osc ',i,' |'
      string1=trim(string1)//string2
    enddo 

    write(string3, '(A1,1X,4(A20,1X))') '#','[fs] |','[eV] |','[eV] |','[eV] |'
    do i=1,nstates
      write(string2,'(1X,A20)') '[eV] |'
      string3=trim(string3)//string2
    enddo
    do i=1,nstates
      write(string2,'(1X,A20)') '[] |'
      string3=trim(string3)//string2
    enddo
    do i=1,nstates
      write(string2,'(1X,A20)') '[] |'
      string3=trim(string3)//string2
    enddo
  endif    
     
  if (write_expec) then
    write(u_expec,'(A1,1X,1000(I20,1X))') '#',(i,i=1,3*nstates+4)
    write(u_expec,'(A)') trim(string1)
    write(u_expec,'(A)') trim(string3)
  endif
    
  if (write_expecmch) then
    write(u_expec_mch,'(A1,1X,1000(I20,1X))') '#',(i,i=1,3*nstates+4)
    write(u_expec_mch,'(A)') trim(string1)
    write(u_expec_mch,'(A)') trim(string3)
  endif


  ! =============================================================================================
  !                                Initialize data
  ! =============================================================================================


  ! spin values in MCH basis
  ! spin values in diagonal basis are calculated from these 
  spin0_s=0.d0
  i=0
  do imult=1,maxmult
    do ims=1,imult
      do istate=1,nstates_m(imult)
        i=i+1
        spin0_s(i)=real(imult-1)
      enddo
    enddo
  enddo

  ! Initial overlap for diabatic populations
  if (write_coeffdiab) then
    ! reference overlap
    ! by default, the reference overlap is the unit matrix
    ref_ovl_ss=dcmplx(0.d0,0.d0)
    do i=1,nstates
      ref_ovl_ss(i,i)=dcmplx(1.d0,0.d0)
    enddo

    ! if "Reference/QM.out" exists, reference overlap is read from there and
    ! Löwdin orthogonalized
    filename='Reference/QM.out'
    inquire(file=filename,exist=exists)
    if (exists) then
      call open_qmout(u_ref,filename)
      call get_overlap(nstates,ref_ovl_ss)
      call close_qmout
      call lowdin(nstates,ref_ovl_ss)
    else
      write(6,*) 'WARNING: Reference overlap not available! Data in coeff_diab.out will be incompatible with other trajectories.'
    endif
  endif

  ! =============================================================================================
  !                                Main loop
  ! =============================================================================================

  write(6,*) 
  write(6,*) 'Running...'
  do
    ! read everything: H, U, DM, overlap, coeff_diag, hopprob, Ekin, active states
    ! random number, runtime for the timestep, geometry, velocity, property matrix
    read(u_dat,*,iostat=io) string1
    if (io/=0) exit
    read(u_dat,*) step
    call matread(nstates,H_MCH_ss,u_dat,string1)
    call matread(nstates,U_ss,u_dat,string1)
    do idir=1,3
      call matread(nstates,DM_mch_ssd(:,:,idir),u_dat,string1)
    enddo
    if (have_overlap==1) then
      call matread(nstates,overlaps_ss,u_dat,string1)
    endif
    call vecread(nstates,coeff_diag_s,u_dat,string1)
    call vecread(nstates,hopprob_s,u_dat,string1)
    read(u_dat,*)
    read(u_dat,*) Ekin
    read(u_dat,*)
    read(u_dat,*) state_diag,state_MCH
    read(u_dat,*)
    read(u_dat,*) randnum
    read(u_dat,*)
    read(u_dat,*) runtime
    if (skip_geom_vel_grad_nac) then
      do i=1,2*natom+2
        read(u_dat,*)
      enddo
    else
      call vec3read(natom,geom_ad,u_dat,string1)
      call vec3read(natom,veloc_ad,u_dat,string1)
    endif
    if (have_property2d==1) then
!       if (.not.is_integer) read(u_dat,*)
      do i=1,n_property2d
        call matread(nstates,Prop2d_xss(i,:,:),u_dat,string1)
      enddo
    endif
    if (have_property1d==1) then
!       read(u_dat,*)
      do i=1,n_property1d
        call vecread(nstates,Prop1d_ys(i,:),u_dat,string1)
      enddo
    endif
    if (have_grad == 1) then
!       read(u_dat,*) 
      if (skip_geom_vel_grad_nac) then
        do i=1,nstates*natom+nstates
          read(u_dat,*)
        enddo
      else
        do i=1,nstates
          call vec3read(natom,grad_mch_sad(i,:,:),u_dat,string1)
        enddo
      endif
    endif
    if (have_NAC == 1) then
!       read(u_dat,*) 
      if (skip_geom_vel_grad_nac) then
        do i=1,nstates*nstates*natom+nstates*nstates
          read(u_dat,*)
        enddo
      else
        do i=1,nstates
          do j=1,nstates
            call vec3read(natom,NAC_ssad(i,j,:,:),u_dat,string1)
          enddo
        enddo
      endif
    endif
    ! ========== Reading is done for this time step =============


    ! calculate Hamiltonian including laser field
    H_diag_ss=H_MCH_ss
    if (laser==2) then
      do idir=1,3
        H_diag_ss=H_diag_ss - DM_mch_ssd(:,:,idir)*real(laser_td(step*nsubsteps+1,idir))
      enddo
    endif
!     call matwrite(nstates,H_diag_ss,0,'Hbefore','F12.9')
!     call matwrite(nstates,U_ss,0,'U','F12.9')
    call transform(nstates,H_diag_ss,U_ss,'utau')
!     call matwrite(nstates,H_diag_ss,0,'Hafter','F12.9')

    ! calculate MCH coefficients and potential energy
    call matvecmultiply(nstates,U_ss,coeff_diag_s,coeff_MCH_s,'n')
    Epot=real(H_diag_ss(state_diag,state_diag))
    
    if (write_coeffdiab) then
      ! calculate diabatic coefficients
      if (step>0) then
        call matmultiply(nstates,ref_ovl_ss,overlaps_ss,A_ss,'nn')
        ref_ovl_ss=A_ss
      endif
      call matvecmultiply(nstates,ref_ovl_ss,coeff_MCH_s,coeff_diab_s,'n')
    endif

    ! calculate oscillator strengths
    if (write_dip .or. write_dipact .or. write_expec .or. write_expecmch .or. write_dm_diag) then
      expec_dm=0.d0
      expec_dm_mch=0.d0
      expec_dm_act=0.d0
      do idir=1,3
        A_ss=DM_mch_ssd(:,:,idir)
        expec_dm_mch=expec_dm_mch+real(A_ss(:,1)*A_ss(1,:))
        call transform(nstates,A_ss,U_ss,'utau')
        DM_diag_ssd(:,:,idir)=A_ss
        expec_dm=expec_dm+real(A_ss(:,1)*A_ss(1,:))
        expec_dm_act=expec_dm_act+real(A_ss(:,state_diag)*A_ss(state_diag,:))
      enddo
      expec_dm=expec_dm*2./3.
      expec_dm_mch=expec_dm_mch*2./3.
      expec_dm_act=expec_dm_act*2./3.
      do i=1,nstates
        expec_dm(i)=expec_dm(i)*real(H_diag_ss(i,i)-H_diag_ss(1,1))
        expec_dm_mch(i)=expec_dm_mch(i)*real(H_MCH_ss(i,i)-H_MCH_ss(1,1))
        expec_dm_act(i)=expec_dm_act(i)*real(H_diag_ss(i,i)-H_diag_ss(state_diag,state_diag))
      enddo
    endif
    
    ! calculate length of projection of dipole moment onto vector
    ! projection length is scalar product of dip mom with unit vector of desired direction (proj_vec/norm_proj_vec)
    if (write_dm_proj) then
      ! calculate desired vectors to project onto
      do iproj=1,nprojections
        do idir=1,3
	  do i=1,2
	    do j=1,2
	      if (proj_atomindex(iproj,i,j) > 0) then
    	        proj_point(i,j) = geom_ad(proj_atomindex(iproj,i,j),idir)
	      else
                ! in principle, two atom indeces can be given per line in projections.inp
                ! if a single atom position is given in a line of projections.inp, then the atom index of the second (undefined) atom is -1
	        proj_point(i,j) = 0.
	      endif
	    enddo
	  enddo
    	  proj_vec_pd(iproj,idir) = (proj_point(1,2) + (proj_point(1,1) - proj_point(1,2)) * distance_scaling(iproj,1)) &
    	  &			  - (proj_point(2,2) + (proj_point(2,1) - proj_point(2,2)) * distance_scaling(iproj,2)) 
	enddo
      enddo
      ! calculate norm
      do iproj=1,nprojections
        norm_proj_vec_p(iproj) = (proj_vec_pd(iproj,1)**2+proj_vec_pd(iproj,2)**2+proj_vec_pd(iproj,3)**2)
      enddo
      ! do projection statewise
      dm_proj_sp = 0.
      do i=1,nstates
        ! scalar product
        do iproj=1,nprojections
          do idir=1,3
            dm_proj_sp(i,iproj) = dm_proj_sp(i,iproj) + abs(DM_diag_ssd(1,i,idir)) * proj_vec_pd(iproj,idir)/norm_proj_vec_p(iproj)
          enddo
        enddo
      enddo
    endif

    ! calculate ionization
    if (write_iondiag .or. write_ionmch) then
      expec_ion_diag=0.d0
      expec_ion_mch=0.d0
      A_ss=Prop2d_xss(1,:,:)
      expec_ion_mch=expec_ion_mch + real(A_ss(:,state_mch)*A_ss(state_mch,:))
      call transform(nstates,A_ss,U_ss,'utau')
      expec_ion_diag=expec_ion_diag + real(A_ss(:,state_diag)*A_ss(state_diag,:))
    endif


    ! calculate spin expectation value
    if (write_spin .or. write_expec .or. write_expecmch) then
      expec_s=0.d0
      do istate=1,nstates
        do jstate=1,nstates
          expec_s(istate)=expec_s(istate) + spin0_s(jstate) * real(U_ss(jstate,istate) * conjg(U_ss(jstate,istate)))
        enddo
      enddo
    endif
    ! ========== Calculating is done for this time step =============


    if (write_energy) then
      ! write to energy.out
      write(u_ener,'(2X,1000(ES20.12E3,1X))') &
      &step*dtstep, Ekin*au2eV, Epot*au2eV, (Epot+Ekin)*au2eV,&
      (real(H_diag_ss(istate,istate)*au2eV),istate=1,nstates)
    endif


    if (write_dip) then
      ! write to fosc.out
      write(u_dm,'(2X,1000(ES20.12E3,1X))') &
      &step*dtstep, expec_dm(state_diag),&
      (expec_dm(istate),istate=1,nstates)
    endif
    if (write_dipact)  then
      ! write to fosc_act.out
      write(u_fosc_act,'(2X,1000(ES20.12E3,1X))') &
      &step*dtstep,(abs(real(H_diag_ss(istate,istate)-H_diag_ss(state_diag,state_diag)))*au2eV,istate=1,nstates),&
      (expec_dm_act(istate),istate=1,nstates)
!       write(u_fosc_act,'(2X,1000(ES20.12E3,1X))') &
!       &step*dtstep,(real(H_diag_ss(istate,istate)-H_diag_ss(state_diag,state_diag))*au2eV,istate=1,nstates),&
!       (expec_dm_act(istate),istate=1,nstates)
    endif


    if (write_iondiag) then
      ! write to ion_diag.out
      write(u_ion_diag,'(2X,ES20.12E3,1X,I20,1X,1000(ES20.12E3,1X))') &
      &step*dtstep,state_diag,(real(H_diag_ss(istate,istate)-H_diag_ss(state_diag,state_diag))*au2eV,istate=1,nstates),&
      (expec_ion_diag(istate),istate=1,nstates)
    endif
    if (write_ionmch) then
      ! write to ion_mch.out
      write(u_ion_mch,'(2X,ES20.12E3,1X,I20,1X,1000(ES20.12E3,1X))') &
      &step*dtstep,state_diag,(real(H_mch_ss(istate,istate)-H_mch_ss(state_mch,state_mch))*au2eV,istate=1,nstates),&
      (expec_ion_mch(istate),istate=1,nstates)
    endif


    if (write_spin) then
      ! write to spin.out
      write(u_spin,'(2X,1000(ES20.12E3,1X))') &
      &step*dtstep, expec_s(state_diag),&
      (expec_s(istate),istate=1,nstates)
    endif


    if (write_coeffdiag) then
      ! calculate sumsq of diagonal coefficients
      sumc=0.d0
      do istate=1,nstates
        sumc=sumc + real(conjg(coeff_diag_s(istate))*coeff_diag_s(istate))
      enddo
      ! write to coeff_diag.out
      write(u_coefd,'(2X,1000(ES20.12E3,1X))') &
      &step*dtstep, sumc,&
      (coeff_diag_s(istate),istate=1,nstates)
      ! write to coeff_class_diag.out
      write(u_classd,'(2X,1000(ES20.12E3,1X))') &
      &step*dtstep, 1.d0,&
      (delta(istate,state_diag),istate=1,nstates)
      ! write to coeff_mixed_diag.out
      write(u_cmixd,'(2X,1000(ES20.12E3,1X))') &
      &step*dtstep, 1.d0,&
      (delta(istate,state_diag),istate=1,nstates)
    endif


    if (write_coeffmch) then
      ! calculate sumsq of MCH coefficients
      sumc=0.d0
      do istate=1,nstates
        sumc=sumc + real(conjg(coeff_MCH_s(istate))*coeff_MCH_s(istate))
      enddo
      ! write to coeff_MCH.out
      write(u_coefm,'(2X,1000(ES20.12E3,1X))') &
      &step*dtstep, sumc,&
      (coeff_MCH_s(istate),istate=1,nstates)
      ! write to coeff_class_MCH.out
      write(u_classm,'(2X,1000(ES20.12E3,1X))') &
      &step*dtstep, 1.d0,&
      (real(abs(U_ss(istate,state_diag))**2),istate=1,nstates)
      ! write to coeff_mixed_MCH.out
      expec_pop=0.d0
      do istate=1,nstates
        expec_pop(istate)=real(abs(U_ss(istate,state_diag))**2)
      enddo
      do kstate=1,nstates
        do istate=1,nstates
          do jstate=1,nstates
            if (istate<jstate) then
              expec_pop(kstate)=expec_pop(kstate)&
              &+2.d0*real( U_ss(kstate,istate)*conjg(U_ss(kstate,jstate))&
              &*coeff_diag_s(istate)*conjg(coeff_diag_s(jstate)))
            endif
          enddo
        enddo
      enddo
      write(u_cmixm,'(2X,1000(ES20.12E3,1X))') &
      &step*dtstep, 1.d0,&
      (expec_pop(istate),istate=1,nstates)
    endif


    if (write_coeffdiab) then
      ! calculate sumsq of diabatic coefficients
      sumc=0.d0
      do istate=1,nstates
        sumc=sumc + dconjg(coeff_diab_s(istate))*coeff_diab_s(istate)
      enddo
      ! write to coeff_diab.out
      write(u_coefdiab,'(2X,1000(ES20.12E3,X))') &
      &step*dtstep, sumc,&
      (coeff_diab_s(istate),istate=1,nstates)
      ! write to coeff_class_diab.out
      call matmultiply(nstates,ref_ovl_ss,U_ss,A_ss,'nn')
      write(u_classdiab,'(2X,1000(ES20.12E3,1X))') &
      &step*dtstep, 1.d0,&
      (real(abs(A_ss(istate,state_diag))**2),istate=1,nstates)
      ! write to coeff_mixed_diab.out
      expec_pop=0.d0
      do istate=1,nstates
        expec_pop(istate)=real(abs(A_ss(istate,state_diag))**2)
      enddo
      do kstate=1,nstates
        do istate=1,nstates
          do jstate=1,nstates
            if (istate<jstate) then
              expec_pop(kstate)=expec_pop(kstate)&
              &+2.d0*real( A_ss(kstate,istate)*conjg(A_ss(kstate,jstate))&
              &*coeff_diag_s(istate)*conjg(coeff_diag_s(jstate)))
            endif
          enddo
        enddo
      enddo
      write(u_cmixdiab,'(2X,1000(ES20.12E3,1X))') &
      &step*dtstep, 1.d0,&
      (expec_pop(istate),istate=1,nstates)
    endif


    if (write_prob) then
      ! calculate cumulative hopping probabilities
      do istate=2,nstates
        hopprob_s(istate)=hopprob_s(istate)+hopprob_s(istate-1)
      enddo
      ! write to prob.out
      write(u_prob,'(2X,1000(ES20.12E3,1X))') &
      &step*dtstep, randnum,&
      (hopprob_s(istate),istate=1,nstates)
    endif


    if (write_expec) then
      ! write to expec.out
      ! this infos are also in energy.out, spin.out and fosc.out
      ! but in order to plot them together they are also written in one file
      write(u_expec,'(2X,1000(ES20.12E3,1X))') &
      &step*dtstep, Ekin*au2eV, Epot*au2eV, (Epot+Ekin)*au2eV,&
      &(real(H_diag_ss(istate,istate)*au2eV),istate=1,nstates),&
      &(expec_s(istate),istate=1,nstates),&
      &(expec_dm(istate),istate=1,nstates)
    endif


    if (write_expecmch) then
      write(u_expec_mch,'(2X,1000(ES20.12E3,1X))') &
      &step*dtstep, Ekin*au2eV, Epot*au2eV, (Epot+Ekin)*au2eV,&
      &(real(H_MCH_ss(istate,istate)*au2eV),istate=1,nstates),&
      &(spin0_s(istate),istate=1,nstates),&
      &(expec_dm_mch(istate),istate=1,nstates)
    endif
    
    if (write_DM_diag) then
      write(u_dm_diag,'(A)') '! 0 Step'
      write(u_dm_diag,'(I12)') step
      call matwrite(nstates, DM_diag_ssd(:,:,1), u_dm_diag, '! 1 Dipole moments X (diag) in a.u.', 'F12.8')
      call matwrite(nstates, DM_diag_ssd(:,:,2), u_dm_diag, '! 1 Dipole moments Y (diag) in a.u.', 'F12.8')
      call matwrite(nstates, DM_diag_ssd(:,:,3), u_dm_diag, '! 1 Dipole moments Z (diag) in a.u.', 'F12.8')
    endif
    
    if (write_dm_proj) then
      write(u_dm_proj,'(A,I12)') '# ',step
      do i=1,nstates
        write(u_dm_proj,'(2X,1000(F12.8,1X))') (dm_proj_sp(i,iproj), iproj=1,nprojections)
      enddo
    endif

    
    ! ========== Writing is done for this time step =============




    ! write progress to screen
    write(*,'(A,A,F9.2,A)',advance='no') achar(13), 't=',step*dtstep,' fs'
  enddo
  write(*,*)



! subroutines for data extractor
  contains

  subroutine print_usage(u)
    implicit none
    integer :: u
    write(u,*) 'Usage: ./data_extractor <flags> <data-file>'
    write(u,*) '       -xl : extralarge = write all output files'
    write(u,*) '       -l  : large = write all output files except diagonal dipole and projection'
    write(u,*) '       -s  : small = write all output files except ionization data, diagonal dipole/projection'
    write(u,*) '       -xs : extrasmall = energy (-e), coeffdiag (-cd), coeffmch (-cm), prob (-p), expec (-x), dip (-d), skip (-sk)'
    write(u,*) '       -e  : write energy file              (output_data/energy.out)'
    write(u,*) '       -d  : write dipole file              (output_data/fosc.out)'
    write(u,*) '       -sp : write spin expec file          (output_data/spin.out)'
    write(u,*) '       -cd : write diag coefficient file    (output_data/coeff_diag.out,'
    write(u,*) '                                             output_data/coeff_class_diag.out,'
    write(u,*) '                                             output_data/coeff_mixed_diag.out)'
    write(u,*) '       -cm : write MCH coefficient file     (output_data/coeff_MCH.out,'
    write(u,*) '                                             output_data/coeff_class_MCH.out,'
    write(u,*) '                                             output_data/coeff_mixed_MCH.out)'
    write(u,*) '       -cb : write diab coefficient file    (output_data/coeff_diab.out,'
    write(u,*) '                                             output_data/coeff_class_diab.out,'
    write(u,*) '                                             output_data/coeff_mixed_diab.out)'
    write(u,*) '       -p  : write hop probability file     (output_data/prob.out)'
    write(u,*) '       -x  : write expec (E,S^2,mu) file    (output_data/expec.out)'
    write(u,*) '       -xm : write MCH expec file           (output_data/expec_MCH.out)'
    write(u,*) '       -da : write dip of active state file (output_data/fosc_act.out)'
    write(u,*) '       -dd : write dip in diag. represent.  (output_data/dip_mom_diag.out)'
!     write(u,*) '       -dp : write projection of dip diag   (output_data/dip_mom_proj.out)'
    write(u,*) '       -id : write diag ion file            (output_data/ion_diag.out)'
    write(u,*) '       -im : write MCH ion file             (output_data/ion_mch.out)'
    write(u,*) '       -sk : skip reading geometries, velocities, gradients, NACs'
  endsubroutine



  real*8 function delta(i,j)
    implicit none
    integer :: i,j
    if (i==j) then
      delta=1.d0
    else
      delta=0.d0
    endif
    return
  endfunction


endprogram


