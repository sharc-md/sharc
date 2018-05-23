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

module LASER_definitions
  
  implicit none
  save
  ! conversion constants
        
  double precision, parameter :: au2A = 0.5291772083d0        ! transforms length in a.u. to Angstrom
  double precision, parameter :: cm2au = 4.5554927d-6         ! transforms energy in wavenumbers to a.u.
  double precision, parameter :: au2fs = 2.4188843265d-2      ! transforms time in a.u. to femtoseconds
  double precision, parameter :: ram2au = 1822.888516331d0    ! transforms mass in relative atomic mass units to a.u.
  double precision, parameter :: J2eV = 6.242d18              ! transforms energy in J to eV
  double precision, parameter :: D2au = 0.3934302014076827d0  ! transforms dipole moment in debye to a.u.
  double precision, parameter :: au2V_m = 5.14220624d11       ! transforms electric field in a.u. to V/m
  double precision, parameter :: au2J = 4.35974381d-18        ! transforms energy in a.u. to J
  double precision, parameter :: au2eV = 27.2113835095688d0   ! transforms energy in a.u. to eV
  double precision, parameter :: D2Cm = 3.336d-30             ! transforms dipole moment in d to C * m
  double precision, parameter :: au2I = 3.509d16              ! transforms (intensity in a.u.)^2 to W/cm^2, attention: watch the square!!!
  double precision, parameter :: au2GV_m = 514.502            ! transforms (field strength in a.u.) to GV/m
  
  double precision, parameter :: pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164d0 ! that's just pi
  
  logical :: debug
  integer :: allocatestatus
  
end module LASER_definitions
