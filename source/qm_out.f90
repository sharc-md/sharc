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

!> # Module QM_OUT
!>
!> \author Sebastian Mai
!> \date 27.02.2015
!>
!> This module realizes the QM -> SHARC half of the QM interface.
!> It provides routines to open and close the QM.out file,
!> and to extract the different quantities found in this file.
!>
module qm_out
implicit none

! =================================================================== !

! private global variables
private

integer, save :: qmout_unit=-1


! =================================================================== !

! private routines
! private 


! =================================================================== !

! public routines
public open_qmout
public close_qmout
public get_hamiltonian
public get_dipoles
public get_gradients
public get_phases
public get_nonadiabatic_ddt
public get_nonadiabatic_ddr
public get_overlap
public get_property
public get_properties_new
public get_dipolegrad
public get_QMruntime


! =================================================================== !

! interfaces


! =================================================================== !

! subroutines

 contains

! =================================================================== !

!> opens a file with filename
!> and assigns unit nunit to it
!> for subsequent extraction of QMout data
subroutine open_qmout(nunit, filename)
  implicit none

  integer, intent(in) :: nunit
  character(len=*), intent(in) :: filename

  integer :: io

  qmout_unit=nunit

  open(nunit,file=filename,status='old',action='read',iostat=io)
  if (io/=0) then
    ! unit 0 is standard error
    write(0,*) 'Failed to open QMout file "',trim(filename),'" !'
    stop 1
  endif

  return

endsubroutine

! =================================================================== !

!> closes a QMout file 
!> and assigns -1 to qmout_unit
subroutine close_qmout
!
  implicit none

  integer :: io

  close(qmout_unit,iostat=io)
  if (io/=0) then
    ! unit 0 is standard error
    write(0,*) 'Failed to close QMout file, unit=',qmout_unit,'!'
    stop 1
  endif

  qmout_unit=-1

  return

endsubroutine

! =================================================================== !

!> check whether a QM.out file is currently open
subroutine check_qmout_unit(routine)
  implicit none
  character(len=*) :: routine

  if (qmout_unit==-1) then
    write(0,*) 'Tried to read from QMout file before opening!'
    write(0,*) 'Routine=',trim(routine)
    stop 1
  endif

endsubroutine

! =================================================================== !

!> rewinds and searches the QM.out file for the given flag
!> each flag marks a quantity (e.g., Hamiltonian)
!> \param flag1 requested flag
subroutine goto_flag(flag1,routine)
  implicit none
  character(len=*) :: routine
!   character :: marker
  integer :: flag1, flag
  integer :: io
  character(len=20) :: string

  rewind(qmout_unit)
  do
!     read(qmout_unit,*,iostat=io) marker,flag
    read(qmout_unit,'(a)',iostat=io) string     !marker,flag
    if (io==-1) then
      write(0,*) 'Quantity not found in QMout file, unit=',qmout_unit
      write(0,*) 'Routine=',trim(routine)
      stop 1
    endif
!     if ( (marker=='!').and.(flag==flag1) ) exit
    if ( string(1:1)=='!' ) then
      read(string(2:20),*,iostat=io) flag
      if ( (io==0).and.(flag==flag1) ) then
!         stat=0
        exit
      endif
    endif
  enddo

endsubroutine

! =================================================================== !

!> rewinds and searches the QM.out file for the given flag
!> each flag marks a quantity (e.g., Hamiltonian)
!> In opposite to goto_flag, this routine does not give an error
!> if the requested flag is not found.
!> \param flag1 requested flag
subroutine goto_flag_nostop(flag1,stat)
  implicit none
!   character :: marker
  integer :: flag1, flag, stat
  integer :: io
  character(len=20) :: string

  rewind(qmout_unit)
  do
    read(qmout_unit,'(a)',iostat=io) string     !marker,flag
    if (io==-1) then
      stat=-1
      return
    endif
    if ( string(1:1)=='!' ) then
      read(string(2:20),*,iostat=io) flag
      if ( (io==0).and.(flag==flag1) ) then
        stat=0
        return
      endif
    endif
  enddo

endsubroutine

! =================================================================== !

!> reads the Hamiltonian matrix from the already opened QMout file
!> the reference energy shift is not applied here
subroutine get_hamiltonian(n, H_ss)
  use matrix
  implicit none
  integer,intent(in) :: n       ! size of the matrix
  complex*16,intent(out) :: H_ss(n,n)
  integer :: icol,irow
  character(len=8000) title

  call check_qmout_unit('get_hamiltonian')

  call goto_flag(1,'get_hamiltonian')

  call matread(n, H_ss, qmout_unit, title)
  read(title,*) irow,icol
  if ( (irow==n).and.(icol==n) ) then
    continue
  else
    write(0,*) 'Hamiltonian matrix has wrong format! nrow=',irow,'ncol=',icol
    stop 1
  endif

  ! if degenerate diagonals, add tiny amount to avoid numerical problems
  do icol=2,n
    do irow=1,icol-1
    if (H_ss(icol,icol)==H_ss(irow,irow)) then
      H_ss(icol,icol)=H_ss(icol,icol)+1.d-12
    endif
    enddo
  enddo

endsubroutine

! =================================================================== !

!> reads the Dipole moment matrix from the already opened QMout file
subroutine get_dipoles(n, DM_ssd)
  use matrix
  implicit none
  integer,intent(in) :: n       ! size of the matrix
  complex*16,intent(out) :: DM_ssd(n,n,3)
  integer :: icol,irow,idir
  character(len=8000) title

  call check_qmout_unit('get_dipoles')

  call goto_flag(2,'get_dipoles')

  do idir=1,3
    call matread(n, DM_ssd(:,:,idir), qmout_unit, title)
    read(title,*) irow,icol
    if ( (irow==n).and.(icol==n) ) then
      continue
    else
      write(0,*) 'Dipole matrix has wrong format! nrow=',irow,'ncol=',icol
      stop 1
    endif
  enddo

endsubroutine

! =================================================================== !

!> reads the gradients from the already opened QMout file
subroutine get_gradients(nstates, natom, grad_sad)
  use matrix
  implicit none
  integer,intent(in) :: nstates,natom
  real*8,intent(out) :: grad_sad(nstates,natom,3)
  integer :: istate,iatom,idir
  character(len=8000) title

  call check_qmout_unit('get_gradients')

  call goto_flag(3,'get_gradients')

  do istate=1,nstates
    call vec3read(natom, grad_sad(istate,:,:), qmout_unit, title)
    read(title,*) iatom,idir
    if ( (iatom==natom).and.(idir==3) ) then
      continue
    else
      write(0,*) 'Gradient has wrong format! natom=',iatom,'ndir=',idir
      stop 1
    endif
  enddo

endsubroutine

! =================================================================== !

!> reads the phases from QMout file
!> if phases are present, stat=0
!> if no phases are present, phase_s=1.0 and stat=-1
subroutine get_phases(n,phase_s,stat)
  use matrix
  implicit none

  integer,intent(in) :: n
  complex*16,intent(out) :: phase_s(n)
  complex*16 :: phase_tmp(n)
  integer,intent(out) :: stat
  integer :: io
  character(len=8000) title

  call check_qmout_unit('get_phases')

  call goto_flag_nostop(7,io)
  if (io==-1) then
    phase_s=dcmplx(1.d0,0.d0)
    stat=-1
    return
  endif

  call vecread(n,phase_tmp,qmout_unit, title)
  phase_s=phase_s*phase_tmp
  read(title,*) io
  if ( io==n ) then
    stat=0
  else
    write(0,*) 'Phase has wrong format! nstates=',io
    stop 1
  endif

endsubroutine

! =================================================================== !

!> reads the NACdt matrix (nstates x nstates) from the already opened QMout file
!> does not read the vectorial couplings, which are read by get_nonadiabatic_ddr
subroutine get_nonadiabatic_ddt(n, T_ss)
  use matrix
  implicit none
  integer,intent(in) :: n       ! size of the matrix
  complex*16,intent(out) :: T_ss(n,n)
  integer :: icol,irow
  character(len=8000) title

  call check_qmout_unit('get_nonadiabatic_ddt')

  call goto_flag(4,'get_nonadiabatic_ddt')

  call matread(n, T_ss, qmout_unit, title)
  read(title,*) irow,icol
  if ( (irow==n).and.(icol==n) ) then
    continue
  else
    write(0,*) 'NAC matrix has wrong format! nrow=',irow,'ncol=',icol
    stop 1
  endif

endsubroutine

! =================================================================== !

!> reads the NACdR vectors (nstates x nstates vectors) from the already opened QMout file
!> does not read the NACdt matrix, which is read by get_nonadiabatic_ddt
subroutine get_nonadiabatic_ddr(nstates, natom, T_ssad)
  use matrix
  implicit none
  integer,intent(in) :: nstates,natom
  real*8,intent(out) :: T_ssad(nstates,nstates,natom,3)
  integer :: icol,irow,iatom,idir
  character(len=8000) title

  call check_qmout_unit('get_nonadiabatic_ddr')

  call goto_flag(5,'get_nonadiabatic_ddr')

  do icol=1,nstates
    do irow=1,nstates
      call vec3read(natom, T_ssad(icol,irow,:,:), qmout_unit, title)
      read(title,*) iatom,idir
      if ( (iatom==natom).and.(idir==3) ) then
        continue
      else
        write(0,*) 'NAC has wrong format! natom=',iatom,'ndir=',idir
        stop 1
      endif
    enddo
  enddo


endsubroutine

! =================================================================== !

!> reads the overlap matrix (nstates x nstates) from the already opened QMout 
subroutine get_overlap(n, S_ss)
  use matrix
  implicit none
  integer,intent(in) :: n       ! size of the matrix
  complex*16,intent(out) :: S_ss(n,n)
  integer :: icol,irow
  character(len=8000) title

  call check_qmout_unit('get_overlap')

  call goto_flag(6,'get_overlap')

  call matread(n, S_ss, qmout_unit, title)
  read(title,*) irow,icol
  if ( (irow==n).and.(icol==n) ) then
    continue
  else
    write(0,*) 'Overlap matrix has wrong format! nrow=',irow,'ncol=',icol
    stop 1
  endif

endsubroutine

! =================================================================== !

!> reads the runtime of the quantum mechanics call
subroutine get_QMruntime(runtime)
  implicit none
  real*8 :: runtime
  integer :: io

  call check_qmout_unit('get_QMruntime')

  call goto_flag_nostop(8,io)
  if (io==-1) then
    runtime=0.
    return
  endif

  read(qmout_unit,*) runtime
  return

endsubroutine

! =================================================================== !

! reads the property matrix from QMout file
! if property matrix present, stat=0
! if no property matrix: property_ss=(0,-123) and stat=-1
subroutine get_property(n,property_ss,stat)
  use matrix
  implicit none

  integer,intent(in) :: n
  complex*16,intent(out) :: property_ss(n,n)
  integer,intent(out) :: stat
  integer :: io, icol,irow
  character(len=8000) title

  call check_qmout_unit('get_property')

  call goto_flag_nostop(11,io)
  if (io==-1) then
    property_ss=dcmplx(0.d0,-123.d0)
    stat=-1
    return
  endif

  call matread(n, property_ss, qmout_unit, title)
  read(title,*) irow,icol
  if ( (irow==n).and.(icol==n) ) then
    stat=0
  else
    write(0,*) 'Property matrix has wrong format! nrow=',irow,'ncol=',icol
    property_ss=dcmplx(0.d0,-123.d0)
    stat=-1
    return
  endif

endsubroutine

! =================================================================== !

! read 1d and 2d properties from QMout file, reading all available data sets
! all allocated parts which are not found in the QMout file are given the value (0,-123)
subroutine get_properties_new(ctrl,traj)
  use matrix
  use definitions
  implicit none

  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

  integer :: i,j,io,nprop,nread
  character(len=8000) :: string

  call check_qmout_unit('get_properties_new')

  ! first, the 1d properties will be processed
  traj%Property1d_ys=dcmplx(0.d0,-123.d0)
  traj%Property1d_labels_y='N/A'

  call goto_flag_nostop(21,io)
  if (io/=-1) then
    read(qmout_unit,*) nprop
    nread=nprop
    if (nprop>ctrl%n_property1d) nread=ctrl%n_property1d
    call vecread(nread, traj%Property1d_labels_y(:nread), qmout_unit, string)
    do j=1,1+nprop-nread
      read(qmout_unit,*) 
    enddo
    do i=1,nread
      call vecread(ctrl%nstates, traj%Property1d_ys(i,:), qmout_unit, string)
    enddo
  endif


  ! second, the 2d properties will be processed
  traj%Property2d_xss=dcmplx(0.d0,-123.d0)
  traj%Property2d_labels_x='N/A'

  ! get property matrix in old way
  call get_property(ctrl%nstates,traj%Property2d_xss(1,:,:),i)

  ! new way
  call goto_flag_nostop(20,io)
  if (io/=-1) then
    read(qmout_unit,*) nprop
    nread=nprop
    if (nprop>ctrl%n_property1d) nread=ctrl%n_property1d
    call vecread(nread, traj%Property2d_labels_x(:nread), qmout_unit, string)
    do j=1,1+nprop-nread
      read(qmout_unit,*) 
    enddo
    do i=1,nread
      call matread(ctrl%nstates, traj%Property2d_xss(i,:,:), qmout_unit, string)
    enddo
  endif

endsubroutine

! =================================================================== !

!> reads the dipole moment derivatives from QMout
subroutine get_dipolegrad(nstates, natom, DMDR_ssdad)
  use matrix
  implicit none
  integer,intent(in) :: nstates,natom
  real*8,intent(out) :: DMDR_ssdad(nstates,nstates,3,natom,3)
  integer :: icol,irow,iatom,idir,ipol
  character(len=8000) title

  call check_qmout_unit('get_dipolegrad')

  call goto_flag(12,'get_dipolegrad')

  do icol=1,nstates
    do irow=1,nstates
      do ipol=1,3
        call vec3read(natom, DMDR_ssdad(icol,irow,ipol,:,:), qmout_unit, title)
        read(title,*) iatom,idir
        if ( (iatom==natom).and.(idir==3) ) then
          continue
        else
          write(0,*) 'NAC has wrong format! natom=',iatom,'ndir=',idir
          stop 1
        endif
      enddo
    enddo
  enddo

endsubroutine

! =================================================================== !

endmodule