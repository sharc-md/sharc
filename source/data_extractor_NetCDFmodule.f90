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


module data_extractor_NetCDFmodule
  implicit none

  public
  save

  type Tncdata
      integer :: id
      integer :: energy_id
  endtype

  type Tncxyz
      integer :: id
      integer :: ian_id
      integer :: crd_id
  endtype

  type Twrite_options
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
    logical :: write_iondiag
    logical :: write_ionmch
    logical :: write_geometry
  end type



  type Tprop_info
    integer :: have_overlap                 !< whether overlap matrices are in the dat file (0=no, 1=yes)
    integer :: have_grad                    !< whether gradients are in the dat file (0=no, 1=yes)
    integer :: have_NAC                     !< whether nonadiabatic couplings are in the dat file (0=no, 1=yes)
    integer :: have_property1d                !< whether property vectors are in the dat file (0=no, 1=yes)
    integer :: have_property2d                !< whether property matrices are in the dat file (0=no, 1=yes)
  end type

  type Tshdata
    integer              :: state_diag, state_MCH, time_step
    real*8               :: Ekin, Epot, randnum, Etot
    integer, allocatable :: nstates_m(:)    !< number of states per multiplicity
    complex*16, allocatable :: H_MCH_ss(:,:),U_ss(:,:),DM_ssd(:,:,:)
    complex*16, allocatable :: Prop2d_xss(:,:,:)
    real*8, allocatable     :: Prop1d_ys(:,:)
    complex*16, allocatable :: coeff_diag_s(:),overlaps_ss(:,:), ref_ovl_ss(:,:)
    real*8, allocatable :: geom_ad(:,:), veloc_ad(:,:)
    integer, allocatable :: ian(:)
    real*8 , allocatable :: hopprob_s(:)
    character*2,allocatable :: element_a(:)
  !> # Information which is calculated at each time step
    complex*16, allocatable :: H_diag_ss(:,:)       !< diagonal Hamiltonian (including laser field)
    complex*16, allocatable :: A_ss(:,:)            !< temporary matrix
    complex*16, allocatable :: coeff_MCH_s(:)       !< MCH coefficient vector
    complex*16, allocatable :: laser_td(:,:)        !< laser field for all timesteps
    complex*16, allocatable :: coeff_diab_s(:)      !< diabatic coefficient vector
    real*8,allocatable :: expec_s(:)                !< spin expectation value per state
    real*8,allocatable :: expec_dm(:)               !< oscillator strength per state
    real*8,allocatable :: expec_dm_mch(:)           !< oscillator strength per state in MCH basis
    real*8,allocatable :: expec_ion_diag(:)           !< oscillator strength per state in MCH basis
    real*8,allocatable :: expec_ion_mch(:)           !< oscillator strength per state in MCH basis
    real*8,allocatable :: expec_dm_act(:)           !< oscillator strength per state with active state as source state
    real*8,allocatable :: spin0_s(:)                !< spin value per MCH state (initialized in the beginning)
    real*8,allocatable :: grad_mch_sad(:,:,:)       !< gradient per MCH state per atom per xyz
    real*8,allocatable :: NAC_ssad(:,:,:,:)         !< nonadiabatic coupling per element (MCH state, MCH state) per atom per xyz
  end type

  type Tgeneral_infos
    integer :: maxmult                      !< maximum multiplicity
    integer :: natom                        !< number of atoms
    integer :: n_property1d                !< 
    integer :: n_property2d                !< 
    integer :: laser                        !< whether a laser field is in the dat file (0=no, 1=, 2=yes)
    integer :: nsteps                       !< number of timesteps from dat file (needed to read the laser field)
    integer :: nsubsteps                    !< number of substeps (needed to read the laser field)

    real*8 :: dtstep                        !< nuclear timestep
    real*8 :: ezero                         !< reference energy
  end type 
contains

! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------------

!   subroutine print_usage(u)
!     implicit none
!     integer :: u
!     write(u,*) 'Usage: ./data_extractor <flags> -f <data-file>'
!     write(u,*) '        -a  : write all output files'
!     write(u,*) '        -s  : standard = write all output files except ionization data'
!     write(u,*) '        -e  : write energy file              (output_data/energy.out)'
!     write(u,*) '        -d  : write dipole file              (output_data/fosc.out)'
!     write(u,*) '        -sp : write spin expec file          (output_data/spin.out)'
! !     write(u,*) '        -cd : write diag coefficient file    (output_data/coeff_diag.out)'
! !     write(u,*) '        -cm : write MCH coefficient file     (output_data/coeff_MCH.out)'
! !     write(u,*) '        -cb : write diab coefficient file    (output_data/coeff_diab.out)'
!     write(u,*) '        -cd : write diag coefficient file    (output_data/coeff_diag.out, output_data/class_diag.out, output_data/cmix_diag)'
!     write(u,*) '        -cm : write MCH coefficient file     (output_data/coeff_MCH.out, output_data/class_MCH.out, output_data/cmix_MCH.out)'
!     write(u,*) '        -cb : write diab coefficient file    (output_data/coeff_diab.out, output_data/class_diab.out, output_data/cmix_diab.out)'
!     write(u,*) '        -p  : write hop probability file     (output_data/prob.out)'
!     write(u,*) '        -x  : write expec (E,S^2,mu) file    (output_data/expec.out)'
!     write(u,*) '        -xm : write MCH expec file           (output_data/expec_MCH.out)'
!     write(u,*) '        -da : write dip of active state file (output_data/fosc_act.out)'
! !     write(u,*) '        -id : write diag ion file            (output_data/ion_diag.out)'
! !     write(u,*) '        -im : write MCH ion file             (output_data/ion_mch.out)'
!     write(u,*) '        -xyz  : write XYZ geometry file (output.xyz)'
!   endsubroutine
  subroutine print_usage(u)
    implicit none
    integer :: u
    write(u,*) 'Usage: ./data_extractor_NetCDF <flags> <data-file>'
    write(u,*) '       -xl : extralarge = write all output files'
    write(u,*) '       -l  : large = write all output files except diagonal dipole and projection'
    write(u,*) '       -s  : small = write all output files except ionization data, diagonal dipole/projection'
    write(u,*) '       -xs : extrasmall = energy (-e), coeffdiag (-cd), coeffmch (-cm), prob (-p), expec (-x), dip (-d), skip (-sk)'
    write(u,*) '       -e  : write energy file              (output_data/energy.out)'
    write(u,*) '       -d  : write dipole file              (output_data/fosc.out)'
    write(u,*) '       -sp : write spin expec file          (output_data/spin.out)'
    write(u,*) '       -cd : write diag coefficient file    (output_data/coeff_diag.out, output_data/class_diag.out, output_data/cmix_diag.out)'
    write(u,*) '       -cm : write MCH coefficient file     (output_data/coeff_MCH.out, output_data/class_MCH.out, output_data/cmix_MCH.out)'
    write(u,*) '       -cb : write diab coefficient file    (output_data/coeff_diab.out, output_data/class_diab.out, output_data/cmix_diab.out)'
    write(u,*) '       -p  : write hop probability file     (output_data/prob.out)'
    write(u,*) '       -x  : write expec (E,S^2,mu) file    (output_data/expec.out)'
    write(u,*) '       -xm : write MCH expec file           (output_data/expec_MCH.out)'
    write(u,*) '       -da : write dip of active state file (output_data/fosc_act.out)'
!     write(u,*) '       -dd : write dip in diag. represent.  (output_data/dip_mom_diag.out)'
!     write(u,*) '       -dp : write projection of dip diag   (output_data/dip_mom_proj.out)'
!     write(u,*) '       -id : write diag ion file            (output_data/ion_diag.out)'
!     write(u,*) '       -im : write MCH ion file             (output_data/ion_mch.out)'
!     write(u,*) '       -sk : skip reading geometries, velocities, gradients, NACs'
    write(u,*) '        -xyz  : write XYZ geometry file (output.xyz)'
  endsubroutine

! -----------------------------------------------------------------------------

  subroutine write_build_info(filename, u_info, u_dat)
      implicit none
      include 'build_info.inc'
      integer, intent(in)              :: u_info, u_dat
      character*8000, intent(in)       :: filename
      character*8000                   :: string1
      integer                          :: io

      ! open the dat file
      if (trim(filename)=='') stop 'No filename given!'
      open(unit=u_dat, file=filename, status='old', action='read', iostat=io)
      if (io/=0) then
        write(*,*) 'File ',trim(filename),' not found'
        stop 1
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

  end subroutine

! -----------------------------------------------------------------------------

  subroutine get_commandline_input(filename, write_options)
      use string
      use input_list
      implicit none
      type(Twrite_options), intent(out) :: write_options
      character*8000, intent(out)       :: filename

      character*8000, allocatable       :: args(:)

      logical                           :: anyoptions
      integer                           :: skipthese
      integer                           :: nargs
      integer                           :: i
  ! print usage if no arguments are given
  nargs=iargc()
  if (nargs<1) then
    call print_usage(0)
    stop 1
  endif

  ! get the command line arguments
  allocate(args(1:nargs))
  do i=1,nargs
    call getarg(i,args(i))
    args(i)=trim(adjustl(args(i)))
  end do

  ! defaults for writing options
  write_options%write_energy    = .false.
  write_options%write_dip       = .false.
  write_options%write_spin      = .false.
  write_options%write_coeffdiag = .false.
  write_options%write_coeffmch  = .false.
  write_options%write_prob      = .false.
  write_options%write_expec     = .false.
  write_options%write_expecmch  = .false.
  write_options%write_coeffdiab = .false.
  write_options%write_dipact    = .false.
  write_options%write_iondiag   = .false.
  write_options%write_ionmch    = .false.
  write_options%write_geometry  = .false.

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
      if ((i+1).GT.nargs) then 
          stop '-f must be followed by a filename'
      endif
      filename = trim(adjustl(args(i+1)))
      skipthese=1
    elseif (trim(args(i)) == "-e") then
      write_options%write_energy = .true.
    elseif (trim(args(i)) == "-d") then
      write_options%write_dip = .true.
    elseif (trim(args(i)) == "-sp") then
      write_options%write_spin = .true.
    elseif (trim(args(i)) == "-cd") then
      write_options%write_coeffdiag = .true.
    elseif (trim(args(i)) == "-cm") then
      write_options%write_coeffmch = .true.
    elseif (trim(args(i)) == "-p") then
      write_options%write_prob = .true.
    elseif (trim(args(i)) == "-x") then
      write_options%write_expec = .true.
    elseif (trim(args(i)) == "-xm") then
      write_options%write_expecmch = .true.
    elseif (trim(args(i)) == "-cb") then
      write_options%write_coeffdiab = .true.
    elseif (trim(args(i)) == "-da") then
      write_options%write_dipact = .true.
    elseif (trim(args(i)) == "-dd") then
      write(0,*) 'Ignoring -dd option in NetCDF mode.'
!       write_dm_diag = .true.
    elseif (trim(args(i)) == "-dp") then
      write(0,*) 'Ignoring -dp option in NetCDF mode.'
!       write_dm_proj = .true.
    elseif (trim(args(i)) == "-id") then
!       write_options%write_iondiag = .true.
      write(0,*) 'Ignoring -id option in NetCDF mode.'
    elseif (trim(args(i)) == "-im") then
!       write_options%write_ionmch = .true.
      write(0,*) 'Ignoring -im option in NetCDF mode.'
    elseif (trim(args(i)) == "-sk") then
!       skip_geom_vel_grad_nac = .true.
      write(0,*) 'Ignoring -sk option in NetCDF mode.'
!     elseif (trim(args(i)) == "-a") then
!       write_options%write_energy = .true.
!       write_options%write_dip = .true.
!       write_options%write_spin = .true.
!       write_options%write_coeffdiag = .true.
!       write_options%write_coeffmch = .true.
!       write_options%write_prob = .true.
!       write_options%write_expec = .true.
!       write_options%write_expecmch = .true.
!       write_options%write_coeffdiab = .true.
!       write_options%write_dipact = .true.
! !       write_options%write_iondiag = .true.
! !       write_options%write_ionmch = .true.
!     elseif (trim(args(i)) == "-s") then
!       write_options%write_energy = .true.
!       write_options%write_dip = .true.
!       write_options%write_spin = .true.
!       write_options%write_coeffdiag = .true.
!       write_options%write_coeffmch = .true.
!       write_options%write_prob = .true.
!       write_options%write_expec = .true.
!       write_options%write_expecmch = .true.
!       write_options%write_coeffdiab = .true.
!       write_options%write_dipact = .true.
!     elseif (trim(args(i)) == "-z") then
!       write_options%write_energy = .true.
!       write_options%write_coeffdiag = .true.
!       write_options%write_coeffmch = .true.
!       write_options%write_prob = .true.
!       write_options%write_expec = .true.
    ! all flags true
    elseif (trim(args(i)) == "-xl") then
      write_options%write_energy = .true.
      write_options%write_dip = .true.
      write_options%write_spin = .true.
      write_options%write_coeffdiag = .true.
      write_options%write_coeffmch = .true.
      write_options%write_prob = .true.
      write_options%write_expec = .true.
      write_options%write_expecmch = .true.
      write_options%write_coeffdiab = .true.
      write_options%write_dipact = .true.
!       write_dm_diag = .true.
!       write_dm_proj = .true.
!       write_iondiag = .true.
!       write_ionmch = .true.
    ! large set of flags true
    elseif (trim(args(i)) == "-l") then
      write_options%write_energy = .true.
      write_options%write_dip = .true.
      write_options%write_spin = .true.
      write_options%write_coeffdiag = .true.
      write_options%write_coeffmch = .true.
      write_options%write_prob = .true.
      write_options%write_expec = .true.
      write_options%write_expecmch = .true.
      write_options%write_coeffdiab = .true.
      write_options%write_dipact = .true.
!       write_iondiag = .true.
!       write_ionmch = .true.
    ! small set of flags true
    elseif (trim(args(i)) == "-s") then
      write_options%write_energy = .true.
      write_options%write_dip = .true.
      write_options%write_spin = .true.
      write_options%write_coeffdiag = .true.
      write_options%write_coeffmch = .true.
      write_options%write_prob = .true.
      write_options%write_expec = .true.
      write_options%write_expecmch = .true.
     write_options% write_coeffdiab = .true.
      write_options%write_dipact = .true.
    ! very small set of flags true
    elseif (trim(args(i)) == "-xs") then
      write_options%write_energy = .true.
      write_options%write_dip = .true.
      write_options%write_coeffdiag = .true.
      write_options%write_coeffmch = .true.
      write_options%write_prob = .true.
      write_options%write_expec = .true.
! ----------------------
    elseif (trim(args(i)) == "-xyz") then
      write_options%write_geometry = .true.
    elseif (trim(args(i)) == "-h") then
      call print_usage(0)
      stop 1
    else
      write(0,*) 'Cannot understand argument: ',trim(args(i))
      !stop 'Input error'
    endif
  enddo

  if (.not.anyoptions) then
    ! defaults for writing options
    write_options%write_energy    = .true.
    write_options%write_dip       = .true.
    write_options%write_spin      = .true.
    write_options%write_coeffdiag = .true.
    write_options%write_coeffmch  = .true.
    write_options%write_prob      = .true.
    write_options%write_expec     = .true.
    write_options%write_expecmch  = .true.
    write_options%write_coeffdiab = .true.
    write_options%write_dipact    = .true.
!     write_options%write_iondiag   = .false.
!     write_options%write_ionmch    = .false.
  endif

  deallocate(args)
  end subroutine get_commandline_input

! -----------------------------------------------------------------------------

  subroutine read_data_file(prop_info, write_options, shdata, general_infos, nstates, u_dat)
  use matrix
  use definitions, only: au2a, au2fs, au2u, au2rcm, au2eV, au2debye
  use qm_out
  use string
  use input_list
  implicit none

  integer, intent(in)           :: u_dat
  type(Twrite_options), intent(inout) :: write_options

  integer, intent(out)          :: nstates
  type(Tprop_info), intent(out) :: prop_info
  type(Tshdata),    intent(out) :: shdata
  type(Tgeneral_infos),    intent(out) :: general_infos 

  ! local variables
  character*8000                :: line, string1, string3
  character*8000, allocatable       :: values(:)
  real*8                        :: fhelper
  logical                       :: is_integer
  integer                       :: i, io, idir,istate,jstate,imult,ims,j,n


  ! Check whether first entry is an integer, which indicates that we have the old format starting with "<integer> ! maxmult"
  read(u_dat,'(I8,A)',iostat=io) general_infos%maxmult, string3
  if (io==0) then
    is_integer = .true.
  else
    is_integer = .false.
  endif
  rewind u_dat


  ! set default switches for different properties
  prop_info%have_NAC=0
  prop_info%have_grad=0
  prop_info%have_overlap=0
  prop_info%have_property1d=0
  prop_info%have_property2d=0


  if (is_integer) then
    write(*,*) 'Found SHARC v1.0 format'
    ! old format has property matrix by default
    prop_info%have_property2d=1
    general_infos%n_property1d=1
    general_infos%n_property2d=1
    ! determine number of states per mult and total number of states
    read(u_dat,*) general_infos%maxmult
    allocate( shdata%nstates_m(general_infos%maxmult) )
    read(u_dat,*) (shdata%nstates_m(i),i=1,general_infos%maxmult)
    nstates=0
    do i=1,general_infos%maxmult
      nstates=nstates+i*shdata%nstates_m(i)
    enddo
    write(*,*) 'Found nstates=',nstates
    read(u_dat,*) general_infos%natom

    ! allocate everything
    allocate( shdata%H_MCH_ss(nstates,nstates), shdata%H_diag_ss(nstates,nstates) )
    allocate( shdata%U_ss(nstates,nstates) )
    allocate( shdata%Prop2d_xss(general_infos%n_property2d,nstates,nstates),shdata%Prop1d_ys(general_infos%n_property1d,nstates) )
    allocate( shdata%overlaps_ss(nstates,nstates), shdata%ref_ovl_ss(nstates,nstates) )
    allocate( shdata%DM_ssd(nstates,nstates,3) )
    allocate( shdata%coeff_diag_s(nstates), shdata%coeff_MCH_s(nstates), shdata%coeff_diab_s(nstates) )
    allocate( shdata%hopprob_s(nstates) )
    allocate( shdata%A_ss(nstates,nstates) )
    allocate( shdata%expec_s(nstates),shdata%expec_dm(nstates),shdata%expec_dm_mch(nstates),shdata%expec_dm_act(nstates) )
    allocate( shdata%expec_ion_diag(nstates),shdata%expec_ion_mch(nstates) )
    allocate( shdata%spin0_s(nstates) )
    allocate( shdata%geom_ad(general_infos%natom,3), shdata%veloc_ad(general_infos%natom,3) )
    allocate( shdata%ian(general_infos%natom) )
    allocate( shdata%element_a(general_infos%natom) )
    call allocate_lapack(nstates)
    shdata%overlaps_ss=dcmplx(0.d0,0.d0)

    ! obtain the timestep, ezero and have_overlap 
    ! most of the information from dat file header is needed in order to interpret the matrices per timestep
    read(u_dat,*) general_infos%dtstep
    general_infos%dtstep=general_infos%dtstep*au2fs
    read(u_dat,*) general_infos%ezero
    read(u_dat,*) prop_info%have_overlap
    read(u_dat,*) general_infos%laser
    read(u_dat,*) general_infos%nsteps
    read(u_dat,*) general_infos%nsubsteps
  
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
      general_infos%maxmult=n
      allocate(shdata%nstates_m(n))
      nstates=0
      ! read number of states per multiplicity
      do i=1,n
        read(values(i),*) shdata%nstates_m(i)
        nstates=nstates+shdata%nstates_m(i)*i
      enddo
      write(*,*) 'Found nstates=',nstates
      ! values is not needed anymore
      deallocate(values)
    else
      ! no nstates keywords => adiabatic dynamics, 1 singlet state
      general_infos%maxmult=1
      allocate(shdata%nstates_m(1))
      shdata%nstates_m(1)=1
      nstates=1
    endif

    ! look up nstates keyword
    line=get_value_from_key('natom',io)
    if (io==0) then
      read(line,*) general_infos%natom
    else
      ! currently natom is needed
      stop 'Error! Number of atoms (keyword: natom) is required!'
    endif

    ! look up properties
    line=get_value_from_key('n_property1d',io)
    if (io==0) then
      read(line,*) general_infos%n_property1d
    else
      general_infos%n_property1d=1
    endif
    line=get_value_from_key('n_property2d',io)
    if (io==0) then
      read(line,*) general_infos%n_property2d
    else
      general_infos%n_property2d=1
    endif

    ! allocate everything
    allocate( shdata%H_MCH_ss(nstates,nstates), shdata%H_diag_ss(nstates,nstates) )
    allocate( shdata%U_ss(nstates,nstates) )
    allocate( shdata%Prop2d_xss(general_infos%n_property2d,nstates,nstates), shdata%Prop1d_ys(general_infos%n_property1d,nstates) )
    allocate( shdata%overlaps_ss(nstates,nstates), shdata%ref_ovl_ss(nstates,nstates) )
    allocate( shdata%DM_ssd(nstates,nstates,3) )
    allocate( shdata%coeff_diag_s(nstates), shdata%coeff_MCH_s(nstates), shdata%coeff_diab_s(nstates) )
    allocate( shdata%hopprob_s(nstates) )
    allocate( shdata%A_ss(nstates,nstates) )
    allocate( shdata%expec_s(nstates),shdata%expec_dm(nstates),shdata%expec_dm_mch(nstates),shdata%expec_dm_act(nstates) )
    allocate( shdata%expec_ion_diag(nstates),shdata%expec_ion_mch(nstates) )
    allocate( shdata%spin0_s(nstates) )
    allocate( shdata%geom_ad(general_infos%natom,3), shdata%veloc_ad(general_infos%natom,3) )
    allocate( shdata%ian(general_infos%natom) )
    allocate( shdata%element_a(general_infos%natom) )
    call allocate_lapack(nstates)
    shdata%overlaps_ss=dcmplx(0.d0,0.d0)

    ! look up dtstep keyword
    line=get_value_from_key('dtstep',io)
    if (io==0) then
      read(line,*) general_infos%dtstep
      general_infos%dtstep=general_infos%dtstep*au2fs
    else
      ! dtstep is needed
      stop 'Error! Time step (keyword: dtstep) is required!'
    endif

    ! look up have_overlap keyword
    line=get_value_from_key('write_overlap',io)
    if (io==0) then
      read(line,*) prop_info%have_overlap
    endif
    if (prop_info%have_overlap == 0 .and. write_options%write_coeffdiab) then
      write(*,*)  'Warning! Writing diabatic coefficients is impossible without overlaps present.'
      write(*,*)  'Unsetting flag -cb'
      write_options%write_coeffdiab = .false.
    endif

    ! look up have_grad keyword
    line=get_value_from_key('write_grad',io)
    if (io==0) then
      read(line,*) prop_info%have_grad
    endif

    ! look up have_NAC keyword
    line=get_value_from_key('write_nacdr',io)
    if (io==0) then
      read(line,*) prop_info%have_NAC
    endif

    ! look up have_property keyword
    line=get_value_from_key('write_property1d',io)
    if (io==0) then
      read(line,*) prop_info%have_property1d
    endif
    line=get_value_from_key('write_property2d',io)
    if (io==0) then
      read(line,*) prop_info%have_property2d
    elseif (io==-1) then
      line=get_value_from_key('write_property',io)      ! backwards compatibility
      if (io==0) read(line,*) prop_info%have_property2d
      general_infos%n_property2d=1
    endif
    if (prop_info%have_property2d == 0 .and. (write_options%write_iondiag .or. write_options%write_ionmch)) then
      write(*,*)  'Warning! Writing ionization probabilities does not make sense if property matrix is not present.'
      write(*,*)  'Unsetting flags -im and -id'
      write_options%write_ionmch = .false.
      write_options%write_iondiag = .false.
    endif

    ! look up laser keyword
    line=get_value_from_key('laser',io)
    if (io==0) then
      read(line,*) general_infos%laser
      if (general_infos%laser==2) then
        line=get_value_from_key('nsteps',io)
        if (io==0) then
          read(line,*) general_infos%nsteps
        else
          stop 'Error! If laser is present, keyword nsteps must be given.'
        endif
        line=get_value_from_key('nsubsteps',io)
        if (io==0) then
          read(line,*) general_infos%nsubsteps
        else
          stop 'Error! If laser is present, keyword nsubsteps must be given.'
        endif
      endif
    else
      general_infos%laser=0
    endif

    if (prop_info%have_grad == 1) then
      allocate( shdata%grad_mch_sad(nstates,general_infos%natom,3) )
    endif

    if (prop_info%have_NAC == 1) then
      allocate( shdata%NAC_ssad(nstates,nstates,general_infos%natom,3) )
    endif

    ! Now we skip over the header array data (atomic numbers, elements, masses)
    read(u_dat,*) string1
    do i=1,general_infos%natom
      read(u_dat,*) string1
      read(string1,*) fhelper 
      shdata%ian(i) = int(fhelper)
    enddo
    read(u_dat,*) string1
    do i=1,general_infos%natom
      read(u_dat,*) shdata%element_a(i)
    enddo
    do i=1,(1+general_infos%natom)
      read(u_dat,*) string1
    enddo

    ! if an explicit laser file is in the dat file, read it now
    ! laser field comes before the time step data
    if (general_infos%laser==2) then
      allocate( shdata%laser_td(general_infos%nsteps*general_infos%nsubsteps+1,3) )
      call vec3read(general_infos%nsteps*general_infos%nsubsteps+1,shdata%laser_td,u_dat,string1)
    endif

    ! skip the "End of header array data" separator line
    read(u_dat,*)

  endif
  end subroutine read_data_file

! -----------------------------------------------------------------------------

  subroutine mk_output_folder()
      implicit none
      logical        :: exists

  ! inquire will not work with Intel compiler, so mkdir is always attempted
      inquire(file="output_data", exist=exists)
      if (.not.exists) then
        write(*,'(A)') ' Creating directory "output_data"'
        call system('mkdir -p output_data')
      else
        write(*,'(A)') ' Writing to directory "output_data"'
      endif
  end subroutine mk_output_folder

! -----------------------------------------------------------------------------

  ! write output file headers
  subroutine write_output_file_headers(nstates, write_options)
      implicit none

      integer, intent(in)                :: nstates
      type(Twrite_options), intent(in)   :: write_options
      character*8000                     :: string1, string3
      character*21                       :: string2
  include 'u_parameter_data_extractor.inc'
      integer                            :: i
  ! write output file headers
  if (write_options%write_energy) then
    write(u_ener,'(A1,1X,1000(I20,1X))') '#',(i,i=1,nstates+4)
    write(u_ener,'(A1,1X,5(A20,1X))') '#','Time |','Ekin |','Epot |','Etot |','=== Energy ===>'
    write(u_ener,'(A1,1X,5(A20,1X))') '#','[fs] |','[eV] |','[eV] |','[eV] |','[eV] |'
  endif

  if (write_options%write_dip) then
    write(u_dm,'(A1,1X,1000(I20,1X))') '#',(i,i=1,nstates+2)
    write(u_dm,'(A1,1X,3(A20,1X))') '#','Time |','f_osc (state) |','=== f_osc ===>'
    write(u_dm,'(A1,1X,3(A20,1X))') '#','[fs] |','[] |','[] |'
  endif

  if (write_options%write_dipact) then
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
  
  if (write_options%write_iondiag) then
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
  
  if (write_options%write_ionmch) then 
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
  
  if (write_options%write_spin) then
    write(u_spin,'(A1,1X,1000(I20,1X))') '#',(i,i=1,nstates+2)
    write(u_spin,'(A1,1X,3(A20,1X))') '#','Time |','Spin_s |','=== Spins ===>'
    write(u_spin,'(A1,1X,3(A20,1X))') '#','[fs] |','[] |','[] |'
  endif
  
  if (write_options%write_coeffdiag) then
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
  
  if (write_options%write_coeffmch) then
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
  
  if (write_options%write_coeffdiab) then
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
  
  if (write_options%write_prob) then
    write(u_prob,'(A1,1X,1000(I20,1X))') '#',(i,i=1,nstates+2)
    write(u_prob,'(A1,1X,3(A20,1X))') '#','Time |','Random Number |','=== cumu Prob ===>'
    write(u_prob,'(A1,1X,3(A20,1X))') '#','[fs] |','[] |','[] |'
  endif
  
  if (write_options%write_expec .or. write_options%write_expecmch) then
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
     
  if (write_options%write_expec) then
    write(u_expec,'(A1,1X,1000(I20,1X))') '#',(i,i=1,3*nstates+4)
    write(u_expec,'(A)') trim(string1)
    write(u_expec,'(A)') trim(string3)
  endif
    
  if (write_options%write_expecmch) then
    write(u_expec_mch,'(A1,1X,1000(I20,1X))') '#',(i,i=1,3*nstates+4)
    write(u_expec_mch,'(A)') trim(string1)
    write(u_expec_mch,'(A)') trim(string3)
  endif



  end subroutine write_output_file_headers

! -----------------------------------------------------------------------------

  subroutine open_output_files(write_options)
      implicit none

      type(Twrite_options) :: write_options

      include 'u_parameter_data_extractor.inc'

      if (write_options%write_energy)    open(unit=u_ener, file='output_data/energy.out', status='replace', action='write')           ! -e
      if (write_options%write_dip)       open(unit=u_dm, file='output_data/fosc.out', status='replace', action='write')               ! -d
      if (write_options%write_spin)      open(unit=u_spin, file='output_data/spin.out', status='replace', action='write')             ! -s
      if (write_options%write_prob)      open(unit=u_prob, file='output_data/prob.out', status='replace', action='write')             ! -p
      if (write_options%write_expec)     open(unit=u_expec, file='output_data/expec.out', status='replace', action='write')           ! -x
      if (write_options%write_expecmch)  open(unit=u_expec_mch, file='output_data/expec_MCH.out', status='replace', action='write')   ! -xm
      if (write_options%write_dipact)    open(unit=u_fosc_act, file='output_data/fosc_act.out', status='replace', action='write')     ! -da
      if (write_options%write_iondiag)   open(unit=u_ion_diag, file='output_data/ion_diag.out', status='replace', action='write')     ! -id
      if (write_options%write_ionmch)    open(unit=u_ion_mch, file='output_data/ion_mch.out', status='replace', action='write')       ! -im


      if (write_options%write_coeffdiag) open(unit=u_coefd, file='output_data/coeff_diag.out', status='replace', action='write')      ! -cd           
      if (write_options%write_coeffmch)  open(unit=u_coefm, file='output_data/coeff_MCH.out', status='replace', action='write')       ! -cm           
      if (write_options%write_coeffdiab) open(unit=u_coefdiab, file='output_data/coeff_diab.out', status='replace', action='write')   ! -cb           

      if (write_options%write_coeffdiag) open(unit=u_classd, file='output_data/coeff_class_diag.out', status='replace', action='write')      ! -cd    
      if (write_options%write_coeffmch)  open(unit=u_classm, file='output_data/coeff_class_MCH.out', status='replace', action='write')       ! -cm    
      if (write_options%write_coeffdiab) open(unit=u_classdiab, file='output_data/coeff_class_diab.out', status='replace', action='write')   ! -cb    

      if (write_options%write_coeffdiag) open(unit=u_cmixd, file='output_data/coeff_mixed_diag.out', status='replace', action='write')      ! -cd
      if (write_options%write_coeffmch)  open(unit=u_cmixm, file='output_data/coeff_mixed_MCH.out', status='replace', action='write')       ! -cm
      if (write_options%write_coeffdiab) open(unit=u_cmixdiab, file='output_data/coeff_mixed_diab.out', status='replace', action='write')   ! -cb

      if (write_options%write_geometry) open(unit=u_xyz, file='output.xyz', status='replace', action='write')   ! -cb


  end subroutine open_output_files

! -----------------------------------------------------------------------------

  subroutine initialize_data(nstates, general_infos, shdata, write_options)
      use qm_out, only: open_qmout, get_overlap, close_qmout
      use matrix, only: lowdin
      implicit none

      integer, intent(in)               :: nstates
      type(Tgeneral_infos), intent(in)  :: general_infos
      type(Twrite_options), intent(in)  :: write_options
      type(Tshdata), intent(inout)      :: shdata
      include 'u_parameter_data_extractor.inc'
      character*8000 :: filename
      logical        :: exists
      integer        :: i, imult, ims, istate
  ! spin values in MCH basis
  ! spin values in diagonal basis are calculated from these 
  shdata%spin0_s=0.d0
  i=0
  do imult=1,general_infos%maxmult
    do ims=1,imult
      do istate=1,shdata%nstates_m(imult)
        i=i+1
        shdata%spin0_s(i)=real(imult-1)
      enddo
    enddo
  enddo

  ! Initial overlap for diabatic populations
  if (write_options%write_coeffdiab) then
    ! reference overlap
    ! by default, the reference overlap is the unit matrix
    shdata%ref_ovl_ss=dcmplx(0.d0,0.d0)
    do i=1,nstates
      shdata%ref_ovl_ss(i,i)=dcmplx(1.d0,0.d0)
    enddo

    ! if "Reference/QM.out" exists, reference overlap is read from there and
    ! Löwdin orthogonalized
    filename='Reference/QM.out'
    inquire(file=filename,exist=exists)
    if (exists) then
      call open_qmout(u_ref,filename)
      call get_overlap(nstates,shdata%ref_ovl_ss)
      call close_qmout
      call lowdin(nstates,shdata%ref_ovl_ss)
    else
      write(6,*) 'WARNING: Reference overlap not available! Data in coeff_diab.out will be incompatible with other trajectories.'
    endif
  endif
  end subroutine initialize_data

! -----------------------------------------------------------------------------

  subroutine read_properties_from_output(nstates, step, u_dat, general_infos, prop_info, shdata, ierr)
      use matrix
      implicit none
      integer, intent(in)          :: u_dat, nstates
      integer, intent(out)         :: step
      integer, intent(out)         :: ierr

      type(Tshdata), intent(inout) :: shdata
      type(Tprop_info), intent(in) :: prop_info
      type(Tgeneral_infos), intent(in) :: general_infos

    character*8000 :: string1 
    integer        :: i, j, idir
    ! read everything: H, U, DM, overlap, coeff_diag, hopprob, Ekin, active states
    ! random number, time_step for the timestep, geometry, velocity, property matrix
    read(u_dat,*,iostat=ierr) string1
    if (ierr/=0) return
    read(u_dat,*) step
    call matread(nstates,shdata%H_MCH_ss,u_dat,string1)
    call matread(nstates,shdata%U_ss,u_dat,string1)
    do idir=1,3
      call matread(nstates,shdata%DM_ssd(:,:,idir),u_dat,string1)
    enddo
    if (prop_info%have_overlap==1) then
      call matread(nstates,shdata%overlaps_ss,u_dat,string1)
    endif
    call vecread(nstates,shdata%coeff_diag_s,u_dat,string1)
    call vecread(nstates,shdata%hopprob_s,u_dat,string1)
    read(u_dat,*)
    read(u_dat,*) shdata%Ekin
    read(u_dat,*)
    read(u_dat,*) shdata%state_diag, shdata%state_MCH
    read(u_dat,*)
    read(u_dat,*) shdata%randnum
    read(u_dat,*)
    read(u_dat,*) shdata%time_step
    shdata%time_step=step
    call vec3read(general_infos%natom,shdata%geom_ad,u_dat,string1)
    call vec3read(general_infos%natom,shdata%veloc_ad,u_dat,string1)
    if (prop_info%have_property2d==1) then
!       if (.not.is_integer) read(u_dat,*)
      do i=1,general_infos%n_property2d
        call matread(nstates,shdata%Prop2d_xss(i,:,:),u_dat,string1)
      enddo
    endif
    if (prop_info%have_property1d==1) then
!       read(u_dat,*)
      do i=1,general_infos%n_property1d
        call vecread(nstates,shdata%Prop1d_ys(i,:),u_dat,string1)
      enddo
    endif
    if (prop_info%have_grad == 1) then
!       read(u_dat,*) 
      do i=1,nstates
        call vec3read(general_infos%natom,shdata%grad_mch_sad(i,:,:),u_dat,string1)
      enddo
    endif
    if (prop_info%have_NAC == 1) then
!       read(u_dat,*) 
      do i=1,nstates
        do j=1,nstates
          call vec3read(general_infos%natom,shdata%NAC_ssad(i,j,:,:),u_dat,string1)
        enddo
      enddo
    endif
  end subroutine read_properties_from_output

! -----------------------------------------------------------------------------

  subroutine write_data_to_file(nstates, step, general_infos, write_options, shdata)
      use definitions, only: au2a, au2fs, au2u, au2rcm, au2eV, au2debye
      use matrix, only: matmultiply
      implicit none

      integer, intent(in)              :: nstates, step
      
      type(Twrite_options), intent(in) :: write_options
      type(Tgeneral_infos), intent(in) :: general_infos
      type(Tshdata), intent(inout)     :: shdata

      real*8 :: expec_pop(nstates)
      complex*16 :: A_ss(nstates,nstates)
      integer :: kstate, jstate


      include 'u_parameter_data_extractor.inc'

      integer :: istate, time_step, iatom, idir
      real*8 :: sumc                                  !< sum of coefficients


    time_step=shdata%time_step


    if (write_options%write_geometry) then
      write(u_xyz,'(I12)') general_infos%natom
      write(u_xyz,'(A5, 1X, F14.5, 1X, I4, 1X, i4)') 't= ',time_step*general_infos%dtstep, shdata%state_diag, shdata%state_MCH
      do iatom=1,general_infos%natom
        write(u_xyz,'(A2,3(1X,F16.9))') shdata%element_a(iatom), (shdata%geom_ad(iatom,idir)*au2a,idir=1,3)
      enddo
    endif




    if (write_options%write_energy) then
      ! write to energy.out
      write(u_ener,'(2X,1000(ES20.12E3,1X))') &
      &time_step*general_infos%dtstep, shdata%Ekin*au2eV, shdata%Epot*au2eV, (shdata%Epot+shdata%Ekin)*au2eV,&
      (real(shdata%H_diag_ss(istate,istate)*au2eV),istate=1,nstates)
    endif


    if (write_options%write_dip) then
      ! write to fosc.out
      write(u_dm,'(2X,1000(ES20.12E3,1X))') &
      &time_step*general_infos%dtstep, shdata%expec_dm(shdata%state_diag),&
      (shdata%expec_dm(istate),istate=1,nstates)
    endif
    if (write_options%write_dipact)  then
      ! write to fosc_act.out
      write(u_fosc_act,'(2X,1000(ES20.12E3,1X))') &
      & time_step*general_infos%dtstep,(abs(real(shdata%H_diag_ss(istate,istate)&
      & -shdata%H_diag_ss(shdata%state_diag,shdata%state_diag)))*au2eV,istate=1,nstates),&
      & (shdata%expec_dm_act(istate),istate=1,nstates)
!       write(u_fosc_act,'(2X,1000(ES20.12E3,1X))') &
!       &time_step*dtstep,(real(H_diag_ss(istate,istate)-H_diag_ss(state_diag,state_diag))*au2eV,istate=1,nstates),&
!       (expec_dm_act(istate),istate=1,nstates)
    endif


    if (write_options%write_iondiag) then
      ! write to ion_diag.out
      write(u_ion_diag,'(2X,ES20.12E3,1X,I20,1X,1000(ES20.12E3,1X))') &
      & time_step*general_infos%dtstep,shdata%state_diag,(real(shdata%H_diag_ss(istate,istate)&
      & -shdata%H_diag_ss(shdata%state_diag,shdata%state_diag))*au2eV,istate=1,nstates),&
      & (shdata%expec_ion_diag(istate),istate=1,nstates)
    endif
    if (write_options%write_ionmch) then
      ! write to ion_mch.out
      write(u_ion_mch,'(2X,ES20.12E3,1X,I20,1X,1000(ES20.12E3,1X))') &
      & time_step*general_infos%dtstep,shdata%state_diag,(real(shdata%H_mch_ss(istate,istate) &
      & -shdata%H_mch_ss(shdata%state_mch,shdata%state_mch))*au2eV,istate=1,nstates),&
      & (shdata%expec_ion_mch(istate),istate=1,nstates)
    endif


    if (write_options%write_spin) then
      ! write to spin.out
      write(u_spin,'(2X,1000(ES20.12E3,1X))') &
      &time_step*general_infos%dtstep, shdata%expec_s(shdata%state_diag),&
      (shdata%expec_s(istate),istate=1,nstates)
    endif


    if (write_options%write_coeffdiag) then
      ! calculate sumsq of diagonal coefficients
      sumc=0.d0
      do istate=1,nstates
        sumc=sumc + real(conjg(shdata%coeff_diag_s(istate))*shdata%coeff_diag_s(istate))
      enddo
      ! write to coeff_diag.out
      write(u_coefd,'(2X,1000(ES20.12E3,1X))') &
      &time_step*general_infos%dtstep, sumc,&
      (shdata%coeff_diag_s(istate),istate=1,nstates)
      ! write to coeff_class_diag.out
      write(u_classd,'(2X,1000(ES20.12E3,1X))') &
      &time_step*general_infos%dtstep, 1.d0,&
      (delta(istate,shdata%state_diag),istate=1,nstates)
      ! write to coeff_mixed_diag.out
      write(u_cmixd,'(2X,1000(ES20.12E3,1X))') &
      &time_step*general_infos%dtstep, 1.d0,&
      (delta(istate,shdata%state_diag),istate=1,nstates)
    endif


    if (write_options%write_coeffmch) then
      ! calculate sumsq of MCH coefficients
      sumc=0.d0
      do istate=1,nstates
        sumc=sumc + real(conjg(shdata%coeff_MCH_s(istate))*shdata%coeff_MCH_s(istate))
      enddo
      ! write to coeff_MCH.out
      write(u_coefm,'(2X,1000(ES20.12E3,1X))') &
      &time_step*general_infos%dtstep, sumc,&
      (shdata%coeff_MCH_s(istate),istate=1,nstates)
      ! write to coeff_class_MCH.out
      write(u_classm,'(2X,1000(ES20.12E3,1X))') &
      &time_step*general_infos%dtstep, 1.d0,&
      (real(abs(shdata%U_ss(istate,shdata%state_diag))**2),istate=1,nstates)
      ! write to coeff_mixed_MCH.out
      expec_pop=0.d0
      do istate=1,nstates
        expec_pop(istate)=real(abs(shdata%U_ss(istate,shdata%state_diag))**2)
      enddo
      do kstate=1,nstates
        do istate=1,nstates
          do jstate=1,nstates
            if (istate<jstate) then
              expec_pop(kstate)=expec_pop(kstate)&
              &+2.d0*real( shdata%U_ss(kstate,istate)*conjg(shdata%U_ss(kstate,jstate))&
              &*shdata%coeff_diag_s(istate)*conjg(shdata%coeff_diag_s(jstate)))
            endif
          enddo
        enddo
      enddo
      write(u_cmixm,'(2X,1000(ES20.12E3,1X))') &
      &time_step*general_infos%dtstep, 1.d0,&
      (expec_pop(istate),istate=1,nstates)
    endif


    if (write_options%write_coeffdiab) then
      ! calculate sumsq of diabatic coefficients
      sumc=0.d0
      do istate=1,nstates
        sumc=sumc + dconjg(shdata%coeff_diab_s(istate))*shdata%coeff_diab_s(istate)
      enddo
      ! write to coeff_diab.out
      write(u_coefdiab,'(2X,1000(ES20.12E3,X))') &
      &time_step*general_infos%dtstep, sumc,&
      (shdata%coeff_diab_s(istate),istate=1,nstates)
      ! write to coeff_class_diab.out
      call matmultiply(nstates,shdata%ref_ovl_ss,shdata%U_ss,A_ss,'nn')
      write(u_classdiab,'(2X,1000(ES20.12E3,1X))') &
      &time_step*general_infos%dtstep, 1.d0,&
      (real(abs(A_ss(istate,shdata%state_diag))**2),istate=1,nstates)
      ! write to coeff_mixed_diab.out
      expec_pop=0.d0
      do istate=1,nstates
        expec_pop(istate)=real(abs(A_ss(istate,shdata%state_diag))**2)
      enddo
      do kstate=1,nstates
        do istate=1,nstates
          do jstate=1,nstates
            if (istate<jstate) then
              expec_pop(kstate)=expec_pop(kstate)&
              &+2.d0*real( A_ss(kstate,istate)*conjg(A_ss(kstate,jstate))&
              &*shdata%coeff_diag_s(istate)*conjg(shdata%coeff_diag_s(jstate)))
            endif
          enddo
        enddo
      enddo
      write(u_cmixdiab,'(2X,1000(ES20.12E3,1X))') &
      &time_step*general_infos%dtstep, 1.d0,&
      (expec_pop(istate),istate=1,nstates)
    endif


    if (write_options%write_prob) then
      ! calculate cumulative hopping probabilities
      do istate=2,nstates
        shdata%hopprob_s(istate)=shdata%hopprob_s(istate)+shdata%hopprob_s(istate-1)
      enddo
      ! write to prob.out
      write(u_prob,'(2X,1000(ES20.12E3,1X))') &
      &time_step*general_infos%dtstep, shdata%randnum,&
      (shdata%hopprob_s(istate),istate=1,nstates)
    endif


    if (write_options%write_expec) then
      ! write to expec.out
      ! this infos are also in energy.out, spin.out and fosc.out
      ! but in order to plot them together they are also written in one file
      write(u_expec,'(2X,1000(ES20.12E3,1X))') &
      &time_step*general_infos%dtstep, shdata%Ekin*au2eV, shdata%Epot*au2eV, (shdata%Epot+shdata%Ekin)*au2eV,&
      &(real(shdata%H_diag_ss(istate,istate)*au2eV),istate=1,nstates),&
      &(shdata%expec_s(istate),istate=1,nstates),&
      &(shdata%expec_dm(istate),istate=1,nstates)
    endif


    if (write_options%write_expecmch) then
      write(u_expec_mch,'(2X,1000(ES20.12E3,1X))') &
      &time_step*general_infos%dtstep, shdata%Ekin*au2eV, shdata%Epot*au2eV, (shdata%Epot+shdata%Ekin)*au2eV,&
      &(real(shdata%H_MCH_ss(istate,istate)*au2eV),istate=1,nstates),&
      &(shdata%spin0_s(istate),istate=1,nstates),&
      &(shdata%expec_dm_mch(istate),istate=1,nstates)
    endif
    ! ========== Writing is done for this time step =============

    end subroutine write_data_to_file

! -----------------------------------------------------------------------------

    subroutine process_data(nstates, step, general_infos, write_options, shdata)
        use matrix, only: transform, matvecmultiply, matmultiply
        implicit none
  integer, intent(in) :: nstates, step
  type(Tgeneral_infos), intent(in) :: general_infos
  type(Twrite_options), intent(in) :: write_options
  type(Tshdata), intent(inout)     :: shdata
  integer :: i, idir, istate, jstate, time_step
  
    time_step=shdata%time_step
    ! calculate Hamiltonian including laser field
    shdata%H_diag_ss=shdata%H_MCH_ss
    if (general_infos%laser==2) then
      do idir=1,3
        shdata%H_diag_ss=shdata%H_diag_ss - shdata%DM_ssd(:,:,idir)*real(shdata%laser_td(time_step*general_infos%nsubsteps+1,idir))
      enddo
    endif
!     call matwrite(nstates,H_diag_ss,0,'Hbefore','F12.9')
!     call matwrite(nstates,U_ss,0,'U','F12.9')
    call transform(nstates,shdata%H_diag_ss,shdata%U_ss,'utau')
!     call matwrite(nstates,H_diag_ss,0,'Hafter','F12.9')

    ! calculate MCH coefficients and potential energy
    call matvecmultiply(nstates,shdata%U_ss,shdata%coeff_diag_s,shdata%coeff_MCH_s,'n')
    shdata%Epot=real(shdata%H_diag_ss(shdata%state_diag,shdata%state_diag))
    
    if (write_options%write_coeffdiab) then
      ! calculate diabatic coefficients
      if (time_step>0) then
        call matmultiply(nstates,shdata%ref_ovl_ss,shdata%overlaps_ss,shdata%A_ss,'nn')
        shdata%ref_ovl_ss=shdata%A_ss
      endif
      call matvecmultiply(nstates,shdata%ref_ovl_ss,shdata%coeff_MCH_s,shdata%coeff_diab_s,'n')
    endif

    ! calculate oscillator strengths
    if (write_options%write_dip .or. write_options%write_dipact .or.&
       &write_options%write_expec .or. write_options%write_expecmch) then
      shdata%expec_dm=0.d0
      shdata%expec_dm_mch=0.d0
      shdata%expec_dm_act=0.d0
      do idir=1,3
        shdata%A_ss=shdata%DM_ssd(:,:,idir)
        shdata%expec_dm_mch=shdata%expec_dm_mch+real(shdata%A_ss(:,1)*shdata%A_ss(1,:))
        call transform(nstates,shdata%A_ss,shdata%U_ss,'utau')
        shdata%expec_dm=shdata%expec_dm+real(shdata%A_ss(:,1)*shdata%A_ss(1,:))
        shdata%expec_dm_act=shdata%expec_dm_act+real(shdata%A_ss(:,shdata%state_diag)*shdata%A_ss(shdata%state_diag,:))
      enddo
      shdata%expec_dm=shdata%expec_dm*2./3.
      shdata%expec_dm_mch=shdata%expec_dm_mch*2./3.
      shdata%expec_dm_act=shdata%expec_dm_act*2./3.
      do i=1,nstates
        shdata%expec_dm(i)=shdata%expec_dm(i)*real(shdata%H_diag_ss(i,i)-shdata%H_diag_ss(1,1))
        shdata%expec_dm_mch(i)=shdata%expec_dm_mch(i)*real(shdata%H_MCH_ss(i,i)-shdata%H_MCH_ss(1,1))
        shdata%expec_dm_act(i)=shdata%expec_dm_act(i)*real(shdata%H_diag_ss(i,i) &
        &-shdata%H_diag_ss(shdata%state_diag,shdata%state_diag))
      enddo
    endif

    ! calculate ionization
    if (write_options%write_iondiag .or. write_options%write_ionmch) then
      shdata%expec_ion_diag=0.d0
      shdata%expec_ion_mch=0.d0
      shdata%A_ss=shdata%Prop2d_xss(1,:,:)
      shdata%expec_ion_mch=shdata%expec_ion_mch + real(shdata%A_ss(:,shdata%state_mch)*shdata%A_ss(shdata%state_mch,:))
      call transform(nstates,shdata%A_ss,shdata%U_ss,'utau')
      shdata%expec_ion_diag=shdata%expec_ion_diag + real(shdata%A_ss(:,shdata%state_diag)*shdata%A_ss(shdata%state_diag,:))
    endif


    ! calculate spin expectation value
    if (write_options%write_spin .or. write_options%write_expec .or. write_options%write_expecmch) then
      shdata%expec_s=0.d0
      do istate=1,nstates
        do jstate=1,nstates
          shdata%expec_s(istate)=shdata%expec_s(istate) &
              &+ shdata%spin0_s(jstate) * real(shdata%U_ss(jstate,istate) * conjg(shdata%U_ss(jstate,istate)))
        enddo
      enddo
    endif
    end subroutine process_data

! -----------------------------------------------------------------------------

    subroutine xyz_from_sharc_to_normal(NAtoms, geom_ad, geom_da)
        implicit none

        integer, intent(in)           :: NAtoms
        real*8, intent(in)            :: geom_ad(NAtoms, 3)
        real*8, intent(out)            :: geom_da(3, NAtoms)

        integer                       :: i, j

        do i=1, NAtoms
           geom_da(:, i) = geom_ad(i, :)
           !do j=1, 3
           !   geom_da(j, i) = geom_ad(i, j)
           !enddo
        enddo

    end subroutine xyz_from_sharc_to_normal
 
! -----------------------------------------------------------------------------

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


end module data_extractor_NetCDFmodule
