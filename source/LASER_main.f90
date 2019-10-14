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

program create_laser

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Simulate arbitrary laser pulses
!! written by Philipp Marquetand
!! www.marquetand.net
!! this version is part of the SHARC suite of programs
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  use LASER_definitions
  use LASER_input
  use LASER_calc
  
  implicit none
  
  integer :: ilasers
  integer :: ixyz
  integer :: it
!   integer :: NE
!   integer :: max_polarizaton_index(1)
  
  real(kind=8) :: t
!   real(kind=8) :: bandwidth_factor
  real(kind=8), allocatable :: envelope(:)
  real(kind=8), allocatable :: env(:,:)
  real(kind=8), allocatable :: momentary_frequency(:)
  real(kind=8), allocatable :: mom_freq(:,:)
  
  complex(kind=8),allocatable :: laser(:,:)
  complex(kind=8),allocatable :: laser_t(:)
!   complex(kind=8),allocatable :: Wig(:,:)
  
  call read_params
  
  allocate (laser(Nt,3), STAT=allocatestatus)
  if (allocatestatus /= 0) stop "*** Not enough memory 1 ***"
  allocate (laser_t(Nt), STAT=allocatestatus)
  if (allocatestatus /= 0) stop "*** Not enough memory 2 ***"
  allocate (envelope(Nt), STAT=allocatestatus)
  if (allocatestatus /= 0) stop "*** Not enough memory 3 ***"
  allocate (env(Nt,Nlasers), STAT=allocatestatus)
  if (allocatestatus /= 0) stop "*** Not enough memory 4 ***"
  allocate (momentary_frequency(Nt), STAT=allocatestatus)
  if (allocatestatus /= 0) stop "*** Not enough memory 5 ***"
  allocate (mom_freq(Nt,Nlasers), STAT=allocatestatus)
  if (allocatestatus /= 0) stop "*** Not enough memory 6 ***"
  laser = 0.
  
  do ilasers = 1, Nlasers
    call field_transform(laser_t,type_envelope(ilasers),field_strength(ilasers),fwhm(ilasers), &
                         pulse_begin(ilasers),pulse_center(ilasers),pulse_center2(ilasers), &
                         pulse_end(ilasers),omega_0(ilasers),phase(ilasers),b_1(ilasers), &
                         b_2(ilasers),b_3(ilasers),b_4(ilasers),dt,t0,Nt,envelope,momentary_frequency) 
    env(:,ilasers) = sqrt(dble(laser_t(:))**2+aimag(laser_t(:))**2)
    mom_freq(:,ilasers) = momentary_frequency(:)
    do ixyz = 1,3
      do it = 1,Nt
        if (abs(laser_t(it)) < threshold(ilasers)) then
          laser_t(it) = 0.
        endif
        laser(it,ixyz) = laser(it,ixyz) + polarization(ixyz,ilasers)*laser_t(it)
      enddo
    enddo
  enddo
  
  write(6,*) 'Writing out laser field'
  
  open (10,file='laser')
  do it = 1,Nt
    t = t0 + (it-1) * dt
    if (realvalued) then
      write(10,'(107(e16.8))') t*au2fs, &
                             dble(laser(it,1)), 0.d0, &
                             dble(laser(it,2)), 0.d0, &
                             dble(laser(it,3)), 0.d0, &
                             (mom_freq(it,ilasers),ilasers=1,Nlasers) 
    else
      write(10,'(107(e16.8))') t*au2fs, &
                             dble(laser(it,1)), aimag(laser(it,1)), &
                             dble(laser(it,2)), aimag(laser(it,2)), &
                             dble(laser(it,3)), aimag(laser(it,3)), &
                             (mom_freq(it,ilasers),ilasers=1,Nlasers)
    endif 
  enddo
  close (10)
  

  deallocate (polarization)
  deallocate (type_envelope)
  deallocate (field_strength)
  deallocate (fwhm)
  deallocate (pulse_begin)
  deallocate (pulse_center)
  deallocate (pulse_center2)
  deallocate (pulse_end)
  deallocate (omega_0)
  deallocate (phase)
  deallocate (b_2)
  deallocate (b_3)
  deallocate (b_4)
  deallocate (laser)
  deallocate (laser_t)
  
end


