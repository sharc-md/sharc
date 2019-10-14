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

module LASER_input
  integer :: Nlasers
  integer :: Nt
  integer,allocatable :: type_envelope(:)
  real(kind=8) :: t0
  real(kind=8) :: dt
  real(kind=8),allocatable :: polarization(:,:)
  real(kind=8),allocatable :: field_strength(:)
  real(kind=8),allocatable :: fwhm(:)
  real(kind=8),allocatable :: pulse_begin(:)
  real(kind=8),allocatable :: pulse_center(:)
  real(kind=8),allocatable :: pulse_center2(:)
  real(kind=8),allocatable :: pulse_end(:)
  real(kind=8),allocatable :: omega_0(:)
  real(kind=8),allocatable :: phase(:)
  real(kind=8),allocatable :: b_1(:)
  real(kind=8),allocatable :: b_2(:)
  real(kind=8),allocatable :: b_3(:)
  real(kind=8),allocatable :: b_4(:)
  real(kind=8),allocatable :: threshold(:)
  logical :: realvalued
  
  contains
subroutine read_params
  
  use LASER_definitions
  
  implicit none
  
  integer :: ixyz
  integer :: ilasers
  integer :: unit_case
  
  real(kind=8) :: tEnd
  real(kind=8) :: polarization_norm

  character*255 :: line

  open(42, file='KEYSTROKES.temp',status='replace',action='write')
  
  write(6,*) 'Number of lasers:'
  read(5,*) line
  call removecomment(line)
  write(42,'(A,A)') trim(line),' ! Number of lasers'
  read(line,*)  Nlasers
  write(6,*) Nlasers

  write(6,*) 'Real-valued field (T) or not (F):'
  read(5,*) line
  call removecomment(line)
  write(42,'(A,A)') trim(line),' ! Real-valued field (T) or not (F)'
  read(line,*) realvalued
  write(6,*) realvalued

  write(6,*) 'Set starting time, end of time and number of time steps (t0[fs],tEnd[fs],Nt) :'
  read(5,'(A)') line
  call removecomment(line)
  write(42,'(A,A)') trim(line),' ! Set starting time, end of time and number of time steps (t0[fs],tEnd[fs],Nt)'
  read(line,*) t0,tEnd,Nt
  write(6,*) t0,tEnd,Nt
  t0   = t0 / au2fs
  tEnd = tEnd / au2fs
  dt   = (tEnd - t0) / (Nt-1)
  write(6,*) 'consequently, we have a step size of ',dt*au2fs
  
  write(6,*) 'Write additional files for debugging (T) or not (F):'
  read(5,*) line
  call removecomment(line)
  write(42,'(A,A)') trim(line),' ! Write additional files for debugging (T) or not (F)'
  read(line,*) debug
  write(6,*) debug

  
  allocate (polarization(3,Nlasers))
  allocate (type_envelope(Nlasers))
  allocate (field_strength(Nlasers))
  allocate (fwhm(Nlasers))
  allocate (pulse_begin(Nlasers))
  allocate (pulse_center(Nlasers))
  allocate (pulse_center2(Nlasers))
  allocate (pulse_end(Nlasers))
  allocate (omega_0(Nlasers))
  allocate (phase(Nlasers))
  allocate (b_1(Nlasers))
  allocate (b_2(Nlasers))
  allocate (b_3(Nlasers))
  allocate (b_4(Nlasers))
  allocate (threshold(Nlasers))
  
  do ilasers = 1, Nlasers
    write(6,*) '             (Empty line to increase readability. Press Enter.)'
    read(5,*) 
    write(42,'(255A)')
    write(6,*)
    
    write(6,*) 'Choose polarization vector (e.g. 2.,0.,0. will be normalized):'
    read(5,'(A)') line
    call removecomment(line)
    write(42,'(255A)') trim(line),' ! Choose polarization vector (e.g. 2.,0.,0. will be normalized)'
    read(line,*)  (polarization(ixyz,ilasers),ixyz=1,3)
    polarization_norm = 0.
    do ixyz = 1,3
      polarization_norm = polarization_norm + (polarization(ixyz,ilasers))**2
    enddo
    polarization(:,ilasers) = polarization(:,ilasers) / sqrt(polarization_norm)
    write(6,*) (polarization(ixyz,ilasers),ixyz=1,3)
     
    write(6,*) 'Choose type of envelope (1=Gaussian,2=Sinusoidal):'
    read(5,'(A)') line
    call removecomment(line)
    write(42,'(255A)') trim(line),' ! Choose type of envelope (1=Gaussian,2=Sinusoidal)'
    read(line,*)  type_envelope(ilasers)
    write(6,*) type_envelope(ilasers)
           
    write(6,*) 'Choose field strength in (1) [GV/m] (2) [TW/cm^2] (3) [a.u.]:'
    read(5,'(A)') line
    call removecomment(line)
    write(42,'(255A)') trim(line),' ! Choose field strength in (1) [GV/m] (2) [TW/cm^2] (3) [a.u.]'
    read(line,*) unit_case
    write(6,*) 'Enter field strength:'
    read(5,'(A)') line
    call removecomment(line)
    write(42,'(255A)') trim(line),' ! Enter field strength'
    read(line,*) field_strength(ilasers)
    select case (unit_case)
      case (1) ! [GV/m]
        write(6,*) field_strength(ilasers),' GV/m'
        field_strength(ilasers) = field_strength(ilasers) / au2GV_m
        write(6,*) '= ',field_strength(ilasers)**2 * au2I/1.d12,' TW/cm^2'
        write(6,*) '= ',field_strength(ilasers),' a.u.'
      case (2) ! [TW/cm^2]
        write(6,*) field_strength(ilasers),' TW/cm^2'
        field_strength(ilasers) = sqrt(field_strength(ilasers)/au2I*1.d12)
        write(6,*) '= ',field_strength(ilasers)*au2GV_m,' GV/m'
        write(6,*) '= ',field_strength(ilasers),' a.u.'
      case (3) ! [a.u.]
        write(6,*) field_strength(ilasers),' a.u.'
        write(6,*) '= ',field_strength(ilasers)**2 * au2I/1.d12,' TW/cm^2'
        write(6,*) '= ',field_strength(ilasers)*au2GV_m,' GV/m'
      case default
        print*, 'Error! Choose 1,2 or 3 (field strength in GV/m,TW/cm^2,a.u.)'
        stop
    end select
    
    write(6,*) 'Set FWHM,begin,center,center2 and end time of pulse (1 + 3 affect Gaussian, 2-5 affect Sinus) [fs]:'
    read(5,'(A)') line
    call removecomment(line)
    write(42,'(255A)') trim(line),' ! Set FWHM,begin,center,center2 and end time of pulse ', &
                       '(1 + 3 affect Gaussian, 2-5 affect Sinus) [fs]'
    read(line,*)  fwhm(ilasers),pulse_begin(ilasers),pulse_center(ilasers),pulse_center2(ilasers),pulse_end(ilasers)
    write(6,*) fwhm(ilasers),pulse_begin(ilasers),pulse_center(ilasers),pulse_center2(ilasers),pulse_end(ilasers)
    fwhm(ilasers)           = fwhm(ilasers) / au2fs
    pulse_begin(ilasers)   = pulse_begin(ilasers) / au2fs
    pulse_center(ilasers)  = pulse_center(ilasers) / au2fs
    pulse_center2(ilasers) = pulse_center2(ilasers) / au2fs
    pulse_end(ilasers)           = pulse_end(ilasers) / au2fs
        
    write(6,*) 'Choose central frequency in (1) [nm] (2) [eV] (3) [a.u.]:'
    read(5,'(A)') line
    call removecomment(line)
    write(42,'(255A)') trim(line),' ! Choose central frequency in (1) [nm] (2) [eV] (3) [a.u.]'
    read(line,*) unit_case
    write(6,*) 'Enter central frequency:'
    read(5,'(A)') line
    call removecomment(line)
    write(42,'(255A)') trim(line),' ! Enter central frequency'
    read(line,*) omega_0(ilasers)
    select case (unit_case)
      case (1) ! [nm]
        write(6,*) omega_0(ilasers),' nm'
        omega_0(ilasers) = (1.d0 / (omega_0(ilasers) * 1.d-7)) * cm2au
        write(6,*) '= ',omega_0(ilasers) * au2eV,' eV'
        write(6,*) '= ',omega_0(ilasers),' a.u.'
      case (2) ! [eV]
        write(6,*) omega_0(ilasers),' eV'
        omega_0(ilasers) = omega_0(ilasers) / au2eV
        write(6,*) '= ',(1.d0 / ((omega_0(ilasers)/cm2au) * 1.d-7)),' nm'
        write(6,*) '= ',omega_0(ilasers),' a.u.'
      case (3) ! [a.u.]
        write(6,*) omega_0(ilasers),' a.u.'
        write(6,*) '= ',(1.d0 / ((omega_0(ilasers)/cm2au) * 1.d-7)),' nm'
        write(6,*) '= ',omega_0(ilasers) * au2eV,' eV'
      case default
        print*, 'Error! Choose 1,2 or 3 (central frequency in nm,eV,a.u.)'
        stop
    end select
      
    write(6,*) 'Choose phase (as multiple of pi):'
    read(5,'(A)') line
    call removecomment(line)
    write(42,'(255A)') trim(line),' ! Choose phase (as multiple of pi)'
    read(line,*)  phase(ilasers)
    write(6,*) phase(ilasers)    
    phase(ilasers)   = phase(ilasers) * pi

    write(6,*) 'Choose colored double pulse (b_1[fs]), which is multiplied with abs(omega-omega_0):'
    read(5,'(A)') line
    call removecomment(line)
    write(42,'(255A)') trim(line),' ! Choose colored double pulse (b_1[fs]), which is multiplied with abs(omega-omega_0)'
    read(line,*)  b_1(ilasers)
    write(6,*) b_1(ilasers)
    b_1(ilasers)     = b_1(ilasers) / au2fs

    write(6,*) 'Choose linear chirp (b_2[fs^2]):'
    read(5,'(A)') line
    call removecomment(line)
    write(42,'(255A)') trim(line),' ! Choose linear chirp (b_2[fs^2])'
    read(line,*)  b_2(ilasers)
    write(6,*) b_2(ilasers)
    b_2(ilasers)     = b_2(ilasers) / au2fs**2

    write(6,*) 'Choose third-order chirp (b_3[fs^3]):'
    read(5,'(A)') line
    call removecomment(line)
    write(42,'(255A)') trim(line),' ! Choose third-order chirp (b_3[fs^3])'
    read(line,*)  b_3(ilasers)
    write(6,*) b_3(ilasers)
    b_3(ilasers)     = b_3(ilasers) / au2fs**3

    write(6,*) 'Choose fourth-order chirp (b_4[fs^4]):'
    read(5,'(A)') line
    call removecomment(line)
    write(42,'(255A)') trim(line),' ! Choose fourth-order chirp (b_4[fs^4])'
    read(line,*)  b_4(ilasers)
    write(6,*) b_4(ilasers)
    b_4(ilasers)     = b_4(ilasers) / au2fs**4

!     write(6,*) 'Choose threshold to set laser field to zero [a.u.]):'
!     read(5,'(A)') line
!     call removecomment(line)
!     write(42,'(255A)') trim(line),' ! Choose threshold to set laser field to zero [a.u.])'
!     read(line,*)  threshold(ilasers)
!     write(6,*) threshold(ilasers)
    threshold(ilasers)=0.d0
    
  enddo
  
  write(6,*) 
  write(6,*) 'Done with input.'

  call flush(6)
  close(42)
  call system('mv KEYSTROKES.temp KEYSTROKES.laser')
    
  return

end subroutine read_params

subroutine removecomment(string)
  implicit none
  character*255, intent(inout) :: string
  integer :: i
  character*255 :: tempstring

  tempstring=' '
  do i=1,len_trim(string)
    if (string(i:i)/='!') then
      tempstring(i:i)=string(i:i)
    else
      exit
    endif
  enddo
  string=trim(tempstring)

  return

  endsubroutine removecomment

end module LASER_input
