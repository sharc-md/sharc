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


!C
!C @author: Maximilian F.S.J. Menger
!C @date: 18.04.2018
!C @version: 0.1.1
!C 
!C wrapper functions for the SHARC python module
!C 
!C All functions etc. are modified version from the original sharc library
!C
!C memory_module is used to store the trajectory and ctrl data in
!C 
!C

#include "sharc_fortran.inc"

module memory_module
!C
!C  stores the traj und ctrl types for the dynamics as "global variables"
!C  same for ncdata
!C
    use definitions
    use definitions_NetCDF
    implicit none
    save
    public

!> \param traj Contains all data which would be private to each trajectory in an ensemble
    type(trajectory_type) :: traj
!> \param ctrl Contains all data which would be shared in an ensemble
    type(ctrl_type) :: ctrl
!> netcdf stuff
    type(Tsharc_ncoutput) :: ncdat

end module memory_module

!C ****************************************************************************
!C
!C  SYSTEM INFORMATION
!C
!C ****************************************************************************

subroutine set_qmin_pointers(Crd_ptr) bind(C, name='setQMinPointers')
    use, intrinsic :: iso_c_binding
    use memory_module, only: traj, ctrl

    implicit none

    type(c_ptr), intent(inout) :: Crd_ptr

    Crd_ptr = C_NULL_PTR

    if (associated(traj%geom_ad)) then
        Crd_ptr = c_loc(traj%geom_ad(1,1))
    endif

    return 
end subroutine set_qmin_pointers

! ------------------------------------------------------

subroutine set_pointers(H, dm, overlap, grad, nac) bind(C, name='setPointers')
    use, intrinsic :: iso_c_binding
    use memory_module, only: traj, ctrl

    implicit none

    type(c_ptr), intent(inout) :: H, dm, overlap, grad
    type(c_ptr), intent(inout) :: nac

    H = C_NULL_PTR
    dm = C_NULL_PTR
    overlap = C_NULL_PTR
    grad = C_NULL_PTR
    nac = C_NULL_PTR

    if (associated(traj%H_MCH_ss)) then
        H = c_loc(traj%H_MCH_ss(1,1))
    endif

    if (associated(traj%DM_ssd)) then
        dm = c_loc(traj%DM_ssd(1,1,1))
    endif

    if (associated(traj%overlaps_ss)) then
        overlap = c_loc(traj%overlaps_ss(1,1))
    endif

    if (associated(traj%overlaps_ss)) then
        overlap = c_loc(traj%overlaps_ss(1,1))
    endif

    if (associated(traj%grad_MCH_sad)) then
        grad = c_loc(traj%grad_MCH_sad(1,1,1))
    endif

    if (associated(traj%NACdr_ssad)) then
        nac = c_loc(traj%NACdr_ssad(1,1,1,1))
    endif
    return 
end subroutine set_pointers

! ------------------------------------------------------

subroutine get_NAtoms(NAtoms)
    use memory_module
    implicit none
    __INT__, intent(out) :: NAtoms

    NAtoms = ctrl%natom

    return
end subroutine get_NAtoms

! ------------------------------------------------------

subroutine get_NSteps(NSteps)
    use memory_module
    implicit none
    __INT__, intent(out) :: NSteps 

    NSteps = ctrl%nsteps

    return
end subroutine get_NSteps

! ------------------------------------------------------

subroutine get_trajstep(IStep)
    use memory_module
    implicit none
    __INT__, intent(out) :: IStep

    IStep = traj%step

    return
end subroutine get_trajstep

! ------------------------------------------------------

subroutine get_element_name(NAtoms, ATyp)
    use iso_c_binding
    use memory_module
    implicit none
    __INT__, intent(in) :: NAtoms
    character(kind=c_char, len=3), dimension(NAtoms) :: ATyp
    __INT__ :: i

    do i=1,NAtoms
        ATyp(i) = traj%element_a(i) // CHAR(0)
    end do

    return
end subroutine get_element_name

! ------------------------------------------------------

subroutine get_IAn(NAtoms, IAn)
    use iso_c_binding
    use memory_module
    implicit none
    __INT__, intent(in) :: NAtoms
    __INT__, dimension(NAtoms) :: IAn
    __INT__ :: i

    do i=1,NAtoms
        IAn(i) = int(traj%atomicnumber_a(i))
    end do

    return
end subroutine get_IAn

! ------------------------------------------------------

subroutine get_current_coordinates(NAtoms, Crd, Ang)
!C
!C 
!C
    use memory_module
    use definitions, only: au2a
    implicit none
    __INT__, intent(in) :: NAtoms
    __REAL__, dimension(3, NAtoms), intent(out) :: Crd
    __INT__, intent(in) :: Ang
    __INT__ :: i,j

    if (Ang .eq. 1) then
        do i=1,NAtoms 
            do j=1,3
                Crd(j, i) = au2a*traj%geom_ad(i,j)
            end do
        end do
    else
        do i=1,NAtoms 
            do j=1,3
                Crd(j, i) = traj%geom_ad(i,j)
            end do
        end do
    end if

    return
end subroutine get_current_coordinates

! ------------------------------------------------------

subroutine get_IPrint(IPrint)
    use definitions, only: printlevel
    implicit none

    __INT__, intent(out) :: IPrint

    IPrint = printlevel

    return
end subroutine get_IPrint

! ------------------------------------------------------

subroutine get_Constants(consts)
    use definitions, only: au2a, au2fs, au2u, au2rcm, &
        au2eV, au2debye
    implicit none

    __REAL__, dimension(6), intent(out) :: consts

    consts(1) = au2a       !< length
    consts(2) = au2fs      !< time
    consts(3) = au2u       !< mass
    consts(4) = au2rcm     !< energy
    consts(5) = au2eV      !< energy
    consts(6) = au2debye   !< energy
    return

end subroutine

!C ****************************************************************************
!C
!C  SHARC funtions to getQMin INFOS 
!C
!C ****************************************************************************

!C ****************************************************
!C  Get static information, does not change over time!
!C ****************************************************

subroutine get_States(string)
    use memory_module, only: traj, ctrl
    implicit none
    integer :: i
    __C_OUT_STRING_S_ :: string


    string = ""
    do i=1,ctrl%maxmult
      write(string,'(A, I3)') trim(string) // ' ', ctrl%nstates_m(i)
    enddo
    write(string,'(A)') trim(string) // CHAR(0)
    
   
    return 
endsubroutine

! ------------------------------------------------------

subroutine get_dt(string)
    use memory_module, only: traj, ctrl
    implicit none
    __C_OUT_STRING_S_ :: string

    string = ""
    write(string,'(a,1x,I7,a)') 'step',traj%step, CHAR(0)
    
endsubroutine

! ------------------------------------------------------

subroutine get_Savedir(string)
!C
!C 
!C
    implicit none
    __C_OUT_STRING_S_ :: string
    character(len=256) :: cwd

    call getcwd(cwd)

    string = ""
    write(string,'(A)') trim(cwd)//'/restart' // CHAR(0)
    
    return
endsubroutine

! ------------------------------------------------------

subroutine get_scalingfactor(scl,soc_scl)
    use memory_module, only: ctrl
    implicit none
    __REAL__, intent(out) :: scl
    __REAL__, intent(out) :: soc_scl

    scl = ctrl%scalingfactor
    soc_scl = ctrl%soc_scaling

    return
endsubroutine


!C ****************************************************
!C  Get static information, does not change over time!
!C ****************************************************

subroutine get_tasks(string, ICALL)
!C
!C Get all single key tasks as a single string, that can be split with
!C string.split() to get all keys
!C
    use memory_module, only: traj, ctrl
    use definitions, only: printlevel, u_log
    use qm, only: select_grad, select_nacdr, select_dipolegrad
    implicit none
    __C_OUT_STRING_S_ :: string
    __INT__, intent(in) :: ICALL
    ! only needed for ICALL .eq. 3
    logical :: old_selg_s(ctrl%nstates)


    if (ICALL .eq. 1) then
        ! if necessary, select the quantities for calculation
        if (ctrl%calc_grad==1) call select_grad(traj,ctrl)
        if (ctrl%calc_nacdr==1) call select_nacdr(traj,ctrl)
        if (ctrl%calc_dipolegrad==1) call select_dipolegrad(traj,ctrl)


        string = ''
        if ((traj%step==0).and..not.(ctrl%track_phase_at_zero==1)) then
          write(string,'(A)') trim(string)  // ' init'
        endif
        if (ctrl%restart) then
          write(string,'(A)') trim(string)  // ' restart'
        endif
        if (ctrl%calc_soc==1) then
          write(string,'(A)') trim(string)  //  ' SOC'
        else
          write(string,'(A)')  trim(string) // ' H'
        endif
        write(string,'(A)')  trim(string)   // ' DM'

        if ((traj%step==0).and.(ctrl%track_phase_at_zero==1)) then
          write(string,'(A)') trim(string)  // ' PHASES'
        endif

        if (traj%step>=1) then
          if (ctrl%calc_nacdt==1) write(string,'(A)')   trim(string) // ' NACDT'
          if (ctrl%calc_overlap==1) write(string,'(A)') trim(string) // ' OVERLAP'
          if (ctrl%calc_phases==1) write(string,'(A)')  trim(string) // ' PHASES'
        endif

        if (ctrl%ionization>0) then
          if (mod(traj%step,ctrl%ionization)==0) then
            write(string,'(A)') trim(string) // ' ION'
          endif
        endif

        if (ctrl%theodore>0) then
          if (mod(traj%step,ctrl%theodore)==0) then
            write(string,'(A)') trim(string) // ' THEODORE'
          endif
        endif

        write(string, '(A)') trim(string) // CHAR(0)

    else if (ICALL .eq. 2) then
        ! select quantities
        if (ctrl%calc_grad==2) call select_grad(traj,ctrl)
        if (ctrl%calc_nacdr==2) call select_nacdr(traj,ctrl)
        if (ctrl%calc_dipolegrad==2) call select_dipolegrad(traj,ctrl)
        !
        string = 'samestep' // CHAR(0)
    else if (ICALL .eq. 3) then
        old_selg_s = traj%selg_s
        call select_grad(traj,ctrl)
        ! compare old and new selection masks
        traj%selg_s=.not.old_selg_s.and.traj%selg_s
        if (printlevel>2) then
            write(u_log,*) 'Missing gradients'
            write(u_log,*) traj%selg_s
            write(u_log,*)
        endif
        string = 'samestep' // CHAR(0)
    else
        write(*,*) "tasks can only be called with icall 0 < icall =< 3"
        call Exit(100)
    endif

    return
endsubroutine

! ------------------------------------------------------

subroutine get_grad(string, ICALL)
    use memory_module, only: traj, ctrl
    implicit none
    __C_OUT_STRING_L_ :: string
    __INT__, intent(in) :: ICALL
    integer :: i

    string = ''
    if (ICALL .eq. 1) then
        select case (ctrl%calc_grad)
          case (0)
            write(string,'(A)') 'all'
          case (1)
            do i=1,ctrl%nstates
              if (traj%selg_s(i)) write(string,'(A,1X,I3)')  trim(string), i
            enddo
          case (2)
              string = ''
        endselect

    else if (ICALL .eq. 2) then
        if (ctrl%calc_grad == 2) then
          do i=1,ctrl%nstates
            if (traj%selg_s(i)) write(string,'(A,1X,I3)') trim(string), i
          enddo
        endif
    else if (ICALL .eq. 3) then
        do i=1,ctrl%nstates
            if (traj%selg_s(i)) write(string,'(A,1X,I3)') trim(string), i
        enddo
    else
        write(*,*) "tasks can only be called with icall 0 < icall =< 3"
        call Exit(100)
    endif
    write(string,'(A)') trim(string) // CHAR(0)
    return
endsubroutine

! ------------------------------------------------------

subroutine get_nacdr(string, ICALL)
    use iso_c_binding
    use memory_module, only: traj, ctrl
    implicit none
    ! __C_OUT_STRING_XL_ :: string
    __C_OUT_STRING_L_ :: string
    __INT__, intent(in) :: ICALL
    integer :: i,j

    string = ''
    if (ICALL .eq. 1) then
        select case (ctrl%calc_nacdr)
          case (-1)
!             write(*,*) 'nonac'
          case (0)
            write(string,'(A)')  'NACDR'
          case (1)
            write(string,'(A)') 'NACDR SELECT' 
            do i=1,ctrl%nstates
              do j=1,ctrl%nstates
                if (traj%selt_ss(j,i)) write(string,'(A,1X,I3,1X,I3)') trim(string) // C_NEW_LINE , i,j
                ! if (traj%selt_ss(j,i)) write(string,'(A,1X,I3,1X,I3)') trim(string) , i,j
              enddo
            enddo
            write(string,'(A)') trim(string) // C_NEW_LINE // 'END'
          case (2)
            write(*,*)
        endselect
    else if (ICALL .eq. 2) then
        if (ctrl%calc_nacdr == 2) then
          do i=1,ctrl%nstates
            do j=1,ctrl%nstates
              if (traj%selt_ss(j,i)) write(string,'(A,1X,I3,1X,I3)') trim(string) // C_NEW_LINE , i,j
            enddo
          enddo
        endif
    else if (ICALL .eq. 3) then
        string = '' 
    else
        write(*,*) "tasks can only be called with icall 0 < icall =< 3"
        call Exit(100)
    endif

    write(string, '(A)') trim(string) // CHAR(0)

    return
endsubroutine

! ------------------------------------------------------

subroutine get_dipolegrad(string, ICALL)
    use iso_c_binding
    use memory_module, only: traj, ctrl
    implicit none
    __C_OUT_STRING_L_ :: string
    __INT__, intent(in) :: ICALL
    integer :: i,j

    string = ''
    select case (ctrl%calc_dipolegrad)
      case (-1)
        string = ''
      case (0)
        write(string,'(A)') 'DMDR'
      case (1)
        write(string,'(A)') 'DMDR SELECT'
        do i=1,ctrl%nstates
          do j=1,ctrl%nstates
            if (traj%seldm_ss(i,j)) write(string,'(A,1X,I3,1X,I3)') trim(string) // C_NEW_LINE ,i,j
          enddo
        enddo
        write(string,'(A)')trim(string) // C_NEW_LINE // 'END'
      case (2)
          string = ''
    endselect

    return
endsubroutine


!C ****************************************************************************
!C
!C  SHARC funtions to setQMout 
!C
!C ****************************************************************************

!if pointer are use, also call the scaling etc.!

subroutine postprocess_qmout_data(IH, IDM, IGrad, IOverlap, INac)
    use memory_module, only: traj, ctrl
    use definitions, only: printlevel, u_log
!C
!C  if pointers are used, do still the postprocessing!
!C
    implicit none
    __INT__, intent(inout) :: IH, IDM, IGrad, IOverlap, INac
    integer :: i,j,istate,jstate

!    write(*,*) "Postprocess setting data", IH, IDM, IGrad, IOverlap

    if (IH .eq. 1) then
        do i=1,ctrl%nstates
            traj%H_MCH_ss(i,i)=traj%H_MCH_ss(i,i)-ctrl%ezero
        enddo

        ! apply scaling factor
        if (ctrl%scalingfactor/=1.d0) then
          traj%H_MCH_ss=traj%H_MCH_ss*ctrl%scalingfactor
        endif

        if (ctrl%soc_scaling/=1.d0) then 
          do istate=1,ctrl%nstates
            do jstate=1,ctrl%nstates
              if (istate.ne.jstate) then 
                traj%H_MCH_ss(istate,jstate)=traj%H_MCH_ss(istate,jstate)*ctrl%soc_scaling
              endif
            enddo
          enddo
        endif

        ! apply frozen-state mask
        do i=1,ctrl%nstates
          do j=1,ctrl%nstates
            if (ctrl%actstates_s(i).neqv.ctrl%actstates_s(j)) traj%H_MCH_ss(j,i)=dcmplx(0.d0,0.d0)
            if ((ctrl%calc_soc/=1).and.(i/=j)) traj%H_MCH_ss(j,i)=dcmplx(0.d0,0.d0)
          enddo
        enddo
        if (printlevel>3) write(u_log,'(A31,A2)') 'Hamiltonian:                   ','OK'
        IH = 0
    end if

    if (IDM .eq. 1) then
        if (printlevel>3) write(u_log,'(A31,A2)') 'Dipole Moments:                ','OK'
        traj%DM_print_ssd=traj%DM_ssd
        ! apply frozen-state mask 
        do i=1,ctrl%nstates
          do j=1,ctrl%nstates
            if (ctrl%actstates_s(i).neqv.ctrl%actstates_s(j)) traj%DM_ssd(j,i,:)=dcmplx(0.d0,0.d0)
          end do
        end do
        IDM = 0
    end if


    if (IGrad .eq. 1) then
        if (printlevel>3) write(u_log,'(A31,A2)') 'Gradients:                     ','OK'
        if (ctrl%scalingfactor/=1.d0) then
            traj%grad_MCH_sad=traj%grad_MCH_sad*ctrl%scalingfactor
        endif
        IGrad = 0
    endif

    if (IOverlap .eq. 1) then
        if (ctrl%calc_overlap == 1) then
            do i=1,ctrl%nstates
              do j=1,ctrl%nstates
                if (ctrl%actstates_s(i).neqv.ctrl%actstates_s(j)) traj%overlaps_ss(j,i)=dcmplx(0.d0,0.d0)
              enddo
            enddo
            if (printlevel>3) write(u_log,'(A31,A2)') 'Overlap matrix:                ','OK'
        endif
        IOverlap = 0
    endif

    ! if it was set, reset to 0, else stay
    if (INac .eq. 1) then
      if (printlevel>3) write(u_log,'(A31,A2)') 'Non-adiabatic couplings (DDR): ','OK'
    endif
    INac = 0

end subroutine postprocess_qmout_data

! ------------------------------------------------------

subroutine set_hamiltonian(N, H_MCH_ss)
    use memory_module, only: traj, ctrl
    implicit none
    integer, intent(in)    :: N
    __COMPLEX__, intent(in) :: H_MCH_ss(N, N) 
    
    integer :: i,j,istate,jstate

    if ( ctrl%nstates .ne. N) then
        write(*,*) "Hamiltonian has wrong dimension!"
        call Exit(1)
    end if

    do i=1,N
        do j=1,N
            ! apply reference energy shift
            if (i .eq. j) then
                traj%H_MCH_ss(j,i) =  H_MCH_ss(j, i)-ctrl%ezero
            else
                traj%H_MCH_ss(j,i) =  H_MCH_ss(j, i)
            endif
        end do
    end do
    ! apply scaling factor
    if (ctrl%scalingfactor/=1.d0) then
      traj%H_MCH_ss=traj%H_MCH_ss*ctrl%scalingfactor
    endif

    if (ctrl%soc_scaling/=1.d0) then
      do istate=1,ctrl%nstates
        do jstate=1,ctrl%nstates
          if (istate.ne.jstate) then
            traj%H_MCH_ss(istate,jstate)=traj%H_MCH_ss(istate,jstate)*ctrl%soc_scaling
          endif
        enddo
      enddo
    endif

    ! apply frozen-state mask
    do i=1,ctrl%nstates
      do j=1,ctrl%nstates
        if (ctrl%actstates_s(i).neqv.ctrl%actstates_s(j)) traj%H_MCH_ss(j,i)=dcmplx(0.d0,0.d0)
        if ((ctrl%calc_soc/=1).and.(i/=j)) traj%H_MCH_ss(j,i)=dcmplx(0.d0,0.d0)
      enddo
    enddo


endsubroutine

! ------------------------------------------------------

subroutine set_phases()
    use memory_module, only: traj, ctrl
!C
!C  Currently phases key word not implemented
!C
    implicit none
    __INT__ :: Istart

    traj%phases_s=dcmplx(1.d0,0.d0)
    traj%phases_found = .false.

    return 
endsubroutine

! ------------------------------------------------------

subroutine set_dipolemoments(N, DM_ssd)
!C
!C  
!C
    use memory_module, only: traj, ctrl
    use definitions, only: printlevel, u_log
    implicit none
    integer, intent(in)    :: N
    __COMPLEX__, intent(in) :: DM_ssd(N,N,3) 
    
    integer :: i,j,k

    if ( ctrl%nstates .ne. N) then
        write(*,*) "Dipole Matrix has wrong dimension!"
        call Exit(1)
    end if

    do k = 1,3
        do i=1,N
            do j=1,N
                traj%DM_ssd(j, i, k) =  DM_ssd(j, i, k)
            end do
        end do
    end do
    if (printlevel>3) write(u_log,'(A31,A2)') 'Dipole Moments:                ','OK'
    traj%DM_print_ssd=traj%DM_ssd
    ! apply frozen-state mask 
    do i=1,ctrl%nstates
      do j=1,ctrl%nstates
        if (ctrl%actstates_s(i).neqv.ctrl%actstates_s(j)) traj%DM_ssd(j,i,:)=dcmplx(0.d0,0.d0)
      end do
    end do

endsubroutine
!C
! ------------------------------------------------------
!C
subroutine set_properties(N)
!C
!C  properties are not loaded! 
!C
    use memory_module, only: traj, ctrl
    implicit none
    integer, intent(in)    :: N
    
    integer :: i,j,k

    write(*,*) "Property not supported, yet"
    call Exit(-1)

    if ( ctrl%nstates .ne. N) then
        write(*,*) ""
        call Exit(1)
    end if


endsubroutine

! ------------------------------------------------------

subroutine set_overlap(N, overlap)
    use memory_module, only: traj, ctrl
    use definitions, only: printlevel, u_log
    implicit none
    integer, intent(in)    :: N
    complex*16, intent(in) :: overlap(N, N) 

    integer :: i,j

    if ( ctrl%nstates .ne. N) then
        write(*,*) "Overlap is of wrong dimension!"
        call Exit(1)
    end if
    if (ctrl%calc_overlap == 1) then
        do i=1,N
            do j=1,N
                traj%overlaps_ss(j, i) = overlap(j, i)
            end do
        end do

        do i=1,ctrl%nstates
          do j=1,ctrl%nstates
            if (ctrl%actstates_s(i).neqv.ctrl%actstates_s(j)) traj%overlaps_ss(j,i)=dcmplx(0.d0,0.d0)
          enddo
        enddo
        if (printlevel>3) write(u_log,'(A31,A2)') 'Overlap matrix:                ','OK'
    endif

endsubroutine

! ------------------------------------------------------

subroutine set_gradients(N, NAtoms, grad_sad) 
    use memory_module, only: traj, ctrl
    use definitions, only: printlevel, u_log
    implicit none
    integer, intent(in)    :: N, NAtoms
    __REAL__, intent(in) :: grad_sad(3, NAtoms, N) 

    integer :: i,j,k

    do i=1,N
        do j=1,NAtoms
            do k=1,3
                traj%grad_MCH_sad(i, j, k) = grad_sad(k, j, i)
            end do
        end do
    end do
    if (ctrl%scalingfactor/=1.d0) then
        traj%grad_MCH_sad=traj%grad_MCH_sad*ctrl%scalingfactor
    endif
endsubroutine

! ------------------------------------------------------

subroutine set_nacs(NStates, NAtoms, nacs)
    use memory_module, only: traj, ctrl
    use definitions, only: printlevel, u_log
    implicit none
    __INT__, intent(in)  :: NStates, NAtoms
    __REAL__, intent(in) :: nacs(3, NAtoms, NStates, NStates)


    __INT__ :: i,j,k,l

    if ( ctrl%nstates .ne. NStates) then
        write(*,*) "Overlap is of wrong dimension!"
        call Exit(1)
    end if

    do i=1, NStates
        do j=1, NStates
            do k=1, NAtoms
                do l=1, 3
                    traj%NACdr_ssad(i,j,k,l) = nacs(l, k, j, i)
                end do
            end do
        end do
    end do

    if (printlevel>3) write(u_log,'(A31,A2)') 'Non-adiabatic couplings (DDR): ','OK'

end subroutine set_nacs

! ------------------------------------------------------

subroutine initial_qm_pre()
    use memory_module, only: traj, ctrl
    use definitions, only: u_log, printlevel
    use output, only: write_logtimestep
    implicit none


    if (printlevel>1) then
      call write_logtimestep(u_log,traj%step,traj%microtime)
    endif

endsubroutine

! ------------------------------------------------------

subroutine initial_qm_post()
    use memory_module, only: traj, ctrl
    use qm
    use definitions
    use electronic
    use matrix
    use output
    implicit none
    integer :: i, iatom, idir

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

! ------------------------------------------------------

subroutine post_process_data(ISecond)
    use memory_module, only: traj, ctrl
    use electronic, only: state_MCH_to_diag, state_diag_to_MCH
    use matrix, only: diagonalize, matwrite
    use definitions, only: printlevel, u_log
    use qm, only: print_qm
    implicit none
    integer :: i
    __INT__, intent(out) :: ISecond
    ! ===============================
    ! all quantities read, post-processing
    ! ===============================
    
    ! here the Hamiltonian is diagonalized, but without phase adjustment
    ! phase adjusted diagonalization is carried out later
    ! here we need to diagonalize only for selection of gradients/couplings/...

    ISecond = 0


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

!    call check_allocation(u_log,ctrl,traj)

    ! get state in all representations
    if ((traj%step==0).and.(ctrl%staterep==1)) then
      traj%state_diag=state_MCH_to_diag(ctrl%nstates,traj%state_MCH,traj%U_ss)
    endif
    traj%state_MCH=state_diag_to_MCH(ctrl%nstates,traj%state_diag,traj%U_ss)

    if (printlevel>4) call print_qm(u_log,traj,ctrl)

    if (ctrl%calc_second==1) then
        ISecond = 1
    endif

endsubroutine

! ------------------------------------------------------

subroutine c_to_f_string(s, str, nchars) 
    use iso_c_binding
    character(kind=c_char, len=1), intent(in) :: s(*)
    character(len=256), intent(out) :: str
    integer i, nchars
    str = ''
    i = 1
    do 
        if (s(i) == c_null_char) exit
        i = i + 1
    end do
    nchars = i - 1 ! Exclude null character form Fortran string

    str = transfer(s(1:nchars), str)

end subroutine c_to_f_string

! ------------------------------------------------------

subroutine write_dat_new(u, traj, ctrl)
  use definitions
  use output, only: write_dat
  use matrix
  use memory_module, only: ncdat
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl
  integer :: istep, u
  select case (ctrl%output_format)
    case (0)
      call write_dat(u, traj, ctrl)
    case (1)
      call write_data_netcdf()
  endselect
  
endsubroutine

subroutine write_data_netcdf()
  use definitions
  use memory_module, only: ncdat, traj, ctrl

  implicit none

  real*8 :: E(3)
  integer :: stride

  E(1) = traj%Etot
  E(2) = traj%Epot
  E(3) = traj%Ekin

  ! check if writing
  stride=ctrl%output_steps_stride(1)
  if (traj%step>=ctrl%output_steps_limits(2)) then
    stride=ctrl%output_steps_stride(2)
  endif
  if (traj%step>=ctrl%output_steps_limits(3)) then
    stride=ctrl%output_steps_stride(3)
  endif
  
  ! TODO: striding is deactivated currently
!   stride=1
  if (modulo(traj%step,stride)==0) then

    call write_sharc_ncoutputdat_istep(&
        & traj%nc_index, &
        & ctrl%natom, &
        & ctrl%nstates, &
          !
        & traj%H_MCH_ss, &
        & traj%U_ss, &
        & traj%DM_print_ssd, &
        & traj%overlaps_ss, &
        & traj%coeff_diag_s, &
        & E, &
        & traj%hopprob_s, &
        & traj%geom_ad, &
        & traj%veloc_ad, &
        & traj%randnum, &
        & traj%state_diag, &
        & traj%state_MCH, &
        & traj%step, &
        & ncdat)
    if (traj%nc_index<0) traj%nc_index=-traj%nc_index
    traj%nc_index=traj%nc_index+1
  endif

end subroutine write_data_netcdf

! ------------------------------------------------------

subroutine close_files()
    use memory_module, only: ncdat, ctrl
    implicit none
  select case (ctrl%output_format)
    case (0)
      continue
    case (1)
      call close_ncfile(ncdat%id)
  endselect
end subroutine close_files

!C ****************************************************************************
!C
!C  SHARC MAIN ROUTINES, without qm_calls, these 
!C
!C ****************************************************************************

subroutine setup_sharc(input_name, IRestart)
!C
!C  setup memory allocation for sharc!
!C  if IRestart = 0, it is not a Restart calculation
!C  if IRestart = 1, it is a Restart calculation!
!C
    use iso_c_binding
    use memory_module, only: traj, ctrl
    use definitions
    use decoherence_afssh,  only: allocate_afssh
    use matrix,  only: allocate_lapack
    use output, only: write_list_header, write_logtimestep
    use input, only: read_input

    implicit none
    ! input file name
    character(len=256, kind=c_char),  intent(in) :: input_name
    character(len=256) :: str
    ! IRestart
    integer, intent(out)   :: IRestart
    !> \param time Define the integer function time()
    integer :: time, nchars

    IRestart = 1
    traj%nc_index = 0

    call c_to_f_string(input_name, str, nchars) 

    traj%time_start=time()
    traj%time_last=traj%time_start

    call read_input(nchars, str(1:nchars), traj, ctrl)
    call allocate_lapack(ctrl%nstates)

    if (.not.ctrl%restart) then
        if (ctrl%decoherence==2) call allocate_afssh(traj, ctrl)
        call write_list_header(u_lis)
        IRestart = 0
    endif

    return
end subroutine setup_sharc

! ------------------------------------------------------

subroutine initial_step(IRestart)
    use memory_module, only: traj, ctrl
    use definitions
    use misc, only: set_time
    use nuclear,  only: Calculate_etot
    use qm, only: Mix_gradients, Update_old, do_initial_qm, QM_processing, NAC_processing
    use restart, only: mkdir_restart, write_restart_ctrl!, write_restart_traj 
    use output, only: write_list_header, write_dat, &
                      write_list_line, write_geom
    implicit none
    __INT__, intent(in)   :: IRestart

    if ( IRestart .eq. 0 ) then
        call QM_processing(traj,ctrl)
        call NAC_processing(traj, ctrl)
        call Calculate_etot(traj,ctrl)   ! not sure if necessary here...
        call Mix_gradients(traj,ctrl)
        call Update_old(traj,ctrl)
        call Calculate_etot(traj,ctrl)
        call set_time(traj)
        traj%microtime=ctrl%dtstep*traj%step
        call write_dat_new(u_dat, traj, ctrl)    
        if (ctrl%output_format==0) then
          call write_geom(u_geo,traj,ctrl)
        endif
        call write_list_line(u_lis,traj,ctrl)
        call mkdir_restart(ctrl)
        call write_restart_ctrl(u_resc,ctrl)
    end if

    return
end subroutine initial_step

! ------------------------------------------------------

subroutine Verlet_xstep(i_step)
    use memory_module, only: traj, ctrl
    use definitions
    use nuclear, only: VelocityVerlet_xstep
    use output, only: write_logtimestep
    implicit none
    __INT__, intent(in) :: i_step

    traj%step=i_step
    call write_logtimestep(u_log, i_step, traj%microtime)
    ! Velocity Verlet x
    call VelocityVerlet_xstep(traj, ctrl)
    return
end subroutine Verlet_xstep

! ------------------------------------------------------

subroutine Verlet_vstep(IRedo)
    use memory_module, only: traj, ctrl
    use definitions
    use qm, only: Adjust_phases, Mix_gradients, QM_processing, NAC_processing
    use electronic, only: propagate, surface_hopping, decoherence, &
                          Calculate_cMCH
    use electronic_laser, only: propagate_laser
    use nuclear, only: VelocityVerlet_vstep, Damp_Velocities, Calculate_ekin, &
                       Calculate_etot, Rescale_Velocities 
    use matrix, only: matwrite

    implicit none
    __INT__, intent(out) :: IRedo

    IRedo = 0
  ! QM Processing
  call QM_processing(traj,ctrl)
  ! Adjust Phases
  call Adjust_phases(traj,ctrl)
  ! Compute NAC in diagonal basis
  call NAC_processing(traj, ctrl)
    ! Mix Gradients
    call Mix_gradients(traj,ctrl)

    ! Velocity Verlet v    (before SH)
    call VelocityVerlet_vstep(traj,ctrl)
    if (ctrl%dampeddyn/=1.d0) call Damp_Velocities(traj,ctrl)
    traj%Ekin=Calculate_ekin(ctrl%natom, traj%veloc_ad, traj%mass_a)
    ! Propagation
    if (ctrl%laser==0) then
      call propagate(traj,ctrl)
    else
      call propagate_laser(traj,ctrl)
    endif
    ! SH
    call surface_hopping(traj,ctrl)
    ! Rescale v
    call Rescale_Velocities(traj,ctrl)
    call Calculate_etot(traj,ctrl)
    ! Decoherence
    call Decoherence(traj,ctrl)
    ! obtain the correct gradient
    call Calculate_cMCH(traj,ctrl)
    traj%microtime=ctrl%dtstep*traj%step
    if (ctrl%calc_grad>=1) then
        IRedo = 1
    endif

    return
end subroutine Verlet_vstep

! ------------------------------------------------------

subroutine Verlet_finalize(IExit, iskip)
    use memory_module, only: traj, ctrl
    use misc
    use definitions
    use qm, only: Update_old, Mix_gradients, NAC_processing
    use electronic, only: kill_after_relaxation
    use output, only: allflush, write_dat, write_list_line, write_geom
   use restart, only: write_restart_traj
    implicit none

    __INT__, intent(out) :: IExit ! if IExit = 0 end loop, else continue
    __INT__, intent(in)  :: iskip ! if IExit = 0 end loop, else continue

    if (traj%kind_of_jump/=0) then
        call NAC_processing(traj, ctrl)
        call Mix_gradients(traj, ctrl)
    endif
    ! Finalization: Variable update, Output, Restart File, Consistency Checks
    call Update_old(traj,ctrl)
    call set_time(traj)
    call write_list_line(u_lis, traj, ctrl)
    call write_dat_new(u_dat, traj, ctrl)


    ! write_restart_traj must be the last command
    if (ctrl%output_format==0) then
      call write_restart_traj(u_rest,ctrl,traj)
      call write_geom(u_geo,traj,ctrl)
    endif

    call allflush()
    ! kill trajectory 
    call kill_after_relaxation(traj, ctrl)
    if ((ctrl%killafter >= 0).and.(traj%steps_in_gs > ctrl%killafter)) then
        IExit = 1 
    endif
    if (check_stop(ctrl%cwd)) then
        IExit = 1
    endif

    if (IExit .eq. 0) then
        ctrl%restart=.false.
    endif

    return
end subroutine Verlet_finalize

! ------------------------------------------------------

subroutine finalize_sharc()
    use memory_module, only: traj, ctrl
    use output, only: write_final
    use definitions, only: deallocate_traj, deallocate_ctrl
    implicit none
    call write_final(traj)
    call write_restart()
    ! add further functions to deallocate sharc main!
    call deallocate_traj(traj)
    call deallocate_ctrl(ctrl)
    ! 
    call close_files()
end subroutine finalize_sharc

! ------------------------------------------------------

subroutine error_finalize_sharc()
    use memory_module, only: traj, ctrl
    use output, only: write_final
    use definitions, only:  u_log, deallocate_traj, deallocate_ctrl
    implicit none
    write(u_log, *) "Error called in finalize SHARC, ending Calculation"
    write(u_log, *) "Restart files etc. saved!"
    call write_restart()
    ! add further functions to deallocate sharc main!
    call deallocate_traj(traj)
    call deallocate_ctrl(ctrl)
    ! 
    call close_files()
end subroutine error_finalize_sharc

! ------------------------------------------------------

subroutine write_restart()

    use memory_module, only: traj, ctrl
    use definitions, only: u_rest, u_resc
    use restart, only: write_restart_ctrl, write_restart_traj 
    implicit none

    call write_restart_ctrl(u_resc, ctrl)
    call write_restart_traj(u_rest, ctrl,traj)

end subroutine

! ------------------------------------------------------
