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


module definitions_NetCDF
    implicit none
    save
    public

type Tsharc_ncoutput
    sequence
    integer :: id
    integer :: H_MCH_id
    integer :: U_id
    integer :: DM_id
    integer :: overlaps_id
    integer :: coeff_diag_id
    integer :: e_id
    integer :: hopprob_id
    integer :: crd_id
    integer :: veloc_id
    integer :: randnum_id
    integer :: state_diag_id
    integer :: state_MCH_id
    integer :: time_step_id
end type

contains

end module definitions_NetCDF
