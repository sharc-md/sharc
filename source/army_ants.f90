!******************************************
!
!    SHARC Program Suite
!
!    Copyright (c) 2023 University of Vienna
!
!    This file is part of SHARC.
!
!    SHARC is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published
!    by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHARC is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    inside the SHARC manual.  If not, see
!    <http://www.gnu.org/licenses/>.
!
!******************************************

!> # Module  army_ants
!> army_ants, army ants hopping/switching procedure for rare event sampling
!>
!> \author Yinan Shu
!> \date June 2, 2022
!>
!> This module contains subroutines that allow to perform army ants algorithm
!> that is capable doing rare event sampling. Original publication:
!> S. Nangia, A. W. Jasper, T. F. Miller III, D. G. Truhlar, J. Chem. Phys.
!> 2004, 120, 3586.

module army_ants
  implicit none

  public army_ants_initialize
  public army_ants_surface_hopping
  public army_ants_surface_switching


contains

! ===========================================================
! army ants weight initialization
  subroutine army_ants_initialize(traj, ctrl)
    use definitions
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl

    traj%army_ants_weight=1.d0


  endsubroutine

! ===========================================================
! army ants surface hopping algorithm
  subroutine army_ants_surface_hopping(traj,ctrl)
    use definitions
    use matrix
    use nuclear
    use electronic
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl

    integer :: istate,jstate,kstate,iatom,idir,i,ilaser
    real*8 :: randnum, cumuprob
    real*8 :: Emax ! energy threshold for frustrated jumps
    real*8 :: Ekin_masked
    real*8 :: sum_kk, sum_vk ! temp variables for kinetic adjustment
    real*8 :: deltaE

    real*8 :: braching_probability(ctrl%nstates)
    real*8 :: weighting(ctrl%nstates)
    integer :: allowed_branching(ctrl%nstates)
    real*8 :: branching_randum(ctrl%nstates)

    if (printlevel>2) then
      write(u_log,*) '============================================================='
      write(u_log,*) '             Doing army ants fewest switches SH'
      write(u_log,*) '============================================================='
    endif

    ! calculate the hopping probabilities
    if (ctrl%integrator==0) then
      call calculate_probabilitiesBSH(ctrl%nstates, traj%coeff_diag_s, traj%preprob_old_s3, &
      &traj%preprob_s3, traj%state_diag, traj%hopprob_s)
      traj%hopprob_s(traj%state_diag)=0.d0
    elseif (ctrl%integrator==1 .or. ctrl%integrator==2) then
      select case (ctrl%hopping_procedure)
        case (0)
          traj%hopprob_s=0.
          braching_probability=0.
        case (1)
          call calculate_probabilities(ctrl%nstates, traj%coeff_diag_old_s, traj%coeff_diag_s, &
          &traj%Rtotal_ss, traj%state_diag, traj%hopprob_s)
          ! Calculate branching probabilities
          do istate=1,ctrl%nstates
            braching_probability(istate)=max(traj%hopprob_s(istate), traj%branching_likehood)
          enddo
        case (2)
          call calculate_probabilities_GFSH(ctrl%nstates, traj%coeff_diag_old_s, traj%coeff_diag_s, &
          &traj%state_diag, traj%hopprob_s)
          ! Calculate branching probabilities
          do istate=1,ctrl%nstates
            braching_probability(istate)=max(traj%hopprob_s(istate), traj%branching_likehood)
          enddo
        case default
          write(0,*) 'Unknown option ',ctrl%hopping_procedure,' to hopping_procedure!'
        stop 1
      endselect
    endif

    if (printlevel>2) then
      write(u_log,*) 'Old and new occupancies and hopping probabilities:'
      do istate=1,ctrl%nstates
        write(u_log,'(3(F14.9,3X))') abs(traj%coeff_diag_old_s(istate))**2,&
        &abs(traj%coeff_diag_s(istate))**2, traj%hopprob_s(istate)
      enddo
    endif

    ! Compute current step contribution of weighting
    weighting(traj%state_diag)=1.d0
    do istate=1,ctrl%nstates
      if (istate.ne.traj%state_diag) then 
        if (braching_probability(istate)>0.d0) then
          weighting(istate)=2.d0*traj%hopprob_s(istate)/braching_probability(istate)
          weighting(traj%state_diag)=weighting(traj%state_diag)-weighting(istate)/2.d0
        else
          weighting(istate)=0.d0
          weighting(traj%state_diag)=weighting(traj%state_diag)-weighting(istate)/2.d0
        endif
      endif
    enddo
    weighting(traj%state_diag)=2.d0*weighting(traj%state_diag)

    if (printlevel>2) then
      call vecwrite(ctrl%nstates, braching_probability,  u_log, 'branching probability', 'F14.9')
      call vecwrite(ctrl%nstates, weighting,  u_log, 'weighting', 'F14.9')
    endif

    ! forced hop to ground state if energy gap to state above is smaller than threshold
    traj%kind_of_jump=0
    if ( ctrl%force_hop_to_gs>0. ) then
      deltaE=abs( real(traj%H_diag_ss(traj%state_diag,traj%state_diag) - traj%H_diag_ss(1,1)) )
      if ( (traj%state_diag/=1).and.(deltaE<=ctrl%force_hop_to_gs) ) then
        traj%hopprob_s=0.
        traj%hopprob_s(1)=1.
        braching_probability=0.
        braching_probability(1)=1.
        traj%kind_of_jump=4
        traj%state_diag=1
      endif
      if (traj%state_diag==1) then
        traj%hopprob_s=0.
        braching_probability=0.
      endif

    else ! not forced hop to ground state, and check branching.

      ! Step 1.
      ! Check whether we will do a branching
      ! initialize allowed branching to 0, i.e. all states are not allowed
      allowed_branching=0
      do istate=1,ctrl%nstates
        if (istate.ne.traj%state_diag) then
          call random_number(randnum)
          branching_randum(istate)=randnum  
          if (branching_randum(istate)<braching_probability(istate)) then 
            allowed_branching(istate)=1
          else
            allowed_branching(istate)=0
          endif
        else 
          branching_randum(istate)=1
        endif
      enddo

      ! Step2. 
      ! Check whether we will do a hopping.
      cumuprob=0.d0
      if (.not. traj%kind_of_jump==4) then
        traj%kind_of_jump=0
      endif
      deltaE=0.d0
      call random_number(randnum)
      traj%randnum=randnum
      stateloop: do istate=1,ctrl%nstates
        ! calculate cumulative probability
        cumuprob=cumuprob + 1.d0/ctrl%nstates
        if (cumuprob > randnum .and. allowed_branching(istate)==1) then

          ! check for laser resonance
          if (ctrl%laser/=0) then
            do ilaser=1,ctrl%nlasers
              deltaE=abs(real(traj%H_diag_ss(istate,istate)-traj%Epot))
              deltaE=abs(deltaE-ctrl%laserenergy_tl(traj%step*ctrl%nsubsteps+1,ilaser))
              if (deltaE<=ctrl%laser_bandwidth) then
                traj%state_diag_old=traj%state_diag
                if (ctrl%hopping_procedure/=0) then
                  traj%state_diag=istate
                endif
                traj%Epot=real(traj%H_diag_ss(istate,istate))
                traj%kind_of_jump=3   ! induced hop
                exit stateloop        ! ************************************************* exit of loop
              endif
            enddo
          endif

          ! check for frustration
          select case (ctrl%ekincorrect)
            case (0)    ! no frustrated jumps: just go on
              continue

            case (1)    ! correct along v, full kinetic energy available
              Ekin_masked=Calculate_ekin_masked(ctrl%natom, traj%veloc_ad, traj%mass_a, ctrl%atommask_a)
              Emax=traj%Etot- traj%Ekin + Ekin_masked
              if (real(traj%H_diag_ss(istate,istate)) > Emax) then
                traj%kind_of_jump=2
                traj%state_diag_frust=istate
                exit stateloop         ! ************************************************* exit of loop
              endif

            case (2)    ! correct along projected v, vibrational kinetic energy available
              call available_ekin(ctrl%natom,&
              &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag,istate,:,:),&
              &traj%mass_a, sum_kk, sum_vk)
              deltaE=4.d0*sum_kk*(traj%Etot-traj%Ekin-&
              &real(traj%H_diag_ss(istate,istate)))+sum_vk**2
              if (deltaE<0.d0) then
                traj%kind_of_jump=2
                traj%state_diag_frust=istate
                exit stateloop         ! ************************************************* exit of loop
              endif

            case (3)    ! correct along T, less energy available
              call available_ekin(ctrl%natom,&
              &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag,istate,:,:),&
              &traj%mass_a, sum_kk, sum_vk)
              deltaE=4.d0*sum_kk*(traj%Etot-traj%Ekin-&
              &real(traj%H_diag_ss(istate,istate)))+sum_vk**2
              if (deltaE<0.d0) then
                traj%kind_of_jump=2
                traj%state_diag_frust=istate
                exit stateloop         ! ************************************************* exit of loop
              endif

            case (4)    ! correct along gradient difference
              call available_ekin(ctrl%natom,&
              &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag,istate,:,:),&
              &traj%mass_a, sum_kk, sum_vk)
              deltaE=4.d0*sum_kk*(traj%Etot-traj%Ekin-&
              &real(traj%H_diag_ss(istate,istate)))+sum_vk**2
              if (deltaE<0.d0) then
                traj%kind_of_jump=2
                traj%state_diag_frust=istate
                exit stateloop         ! ************************************************* exit of loop
              endif

            case (5)    ! correct along projected NAC/M
              call available_ekin(ctrl%natom,&
              &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag,istate,:,:),&
              &traj%mass_a, sum_kk, sum_vk)
              deltaE=4.d0*sum_kk*(traj%Etot-traj%Ekin-&
              &real(traj%H_diag_ss(istate,istate)))+sum_vk**2
              if (deltaE<0.d0) then
                traj%kind_of_jump=2
                traj%state_diag_frust=istate
                exit stateloop         ! ************************************************* exit of loop
              endif

            case (6)    ! correct along projected gradient difference
              call available_ekin(ctrl%natom,&
              &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag,istate,:,:),&
              &traj%mass_a, sum_kk, sum_vk)
              deltaE=4.d0*sum_kk*(traj%Etot-traj%Ekin-&
              &real(traj%H_diag_ss(istate,istate)))+sum_vk**2
              if (deltaE<0.d0) then
                traj%kind_of_jump=2
                traj%state_diag_frust=istate
                exit stateloop         ! ************************************************* exit of loop
              endif

            case (7)    ! correct along effective NAC/M
              call available_ekin(ctrl%natom,&
              &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag,istate,:,:),&
              &traj%mass_a, sum_kk, sum_vk)
              deltaE=4.d0*sum_kk*(traj%Etot-traj%Ekin-&
              &real(traj%H_diag_ss(istate,istate)))+sum_vk**2
              if (deltaE<0.d0) then
                traj%kind_of_jump=2
                traj%state_diag_frust=istate
                exit stateloop         ! ************************************************* exit of loop
              endif

            case (8)    ! correct along projected effective NAC/M
              call available_ekin(ctrl%natom,&
              &traj%veloc_ad,traj%hopping_direction_ssad(traj%state_diag,istate,:,:),&
              &traj%mass_a, sum_kk, sum_vk)
              deltaE=4.d0*sum_kk*(traj%Etot-traj%Ekin-&
              &real(traj%H_diag_ss(istate,istate)))+sum_vk**2
              if (deltaE<0.d0) then
                traj%kind_of_jump=2
                traj%state_diag_frust=istate
                exit stateloop         ! ************************************************* exit of loop
              endif

          endselect

          ! neither in resonance nor frustrated, we have a surface hop!
          traj%state_diag_old=traj%state_diag
          traj%state_MCH_old=traj%state_MCH
          if (ctrl%hopping_procedure/=0) then
            traj%state_diag=istate
            traj%Epot=real(traj%H_diag_ss(istate,istate))
          endif
          if (.not. traj%kind_of_jump==4) then
            traj%kind_of_jump=1
          endif
          exit stateloop               ! ************************************************* exit of loop

        endif ! if (cumuprob > randnum .and. allowed_branching(istate)==1) then
      enddo stateloop

    endif ! if ( ctrl%force_hop_to_gs>0. ) then

    ! finalize trajectory weight
    traj%army_ants_weight=traj%army_ants_weight*weighting(traj%state_diag)

    if (printlevel>2) then
      write(u_log,*)
      write(u_log,'(A,1X,F12.9)') 'Random number:',randnum
      if (ctrl%hopping_procedure==0) then
        write(u_log,*)
        write(u_log,'(A,1X,F12.9)') '*** Hopping is forbidden: new state = old state ***'
        write(u_log,*)
      endif
      select case (traj%kind_of_jump)
        case (0)
          write(u_log,*) 'No jump performed.'
        case (1)
          write(u_log,'(A,1X,I4,1X,A)') 'Old state:',traj%state_diag_old,'(diag)'
          write(u_log,'(A,1X,I4,1X,A)') 'New state:',traj%state_diag,'(diag)'
          write(u_log,'(A,1X,F12.9,1X,A)') 'Jump is not frustrated, new Epot=',traj%Epot*au2ev,'eV'
        case (2)
          write(u_log,'(A,1X,I4,1X,A,1X,E12.5,1X,A)') 'Jump to state',istate,'is frustrated by',deltaE*au2eV,'eV.'
          if (ctrl%laser==2) then
            write(u_log,'(A,1X,F16.9,1X,A,1X,F16.9,1X,A)') &
            &'Detuning:',deltaE*au2eV,'eV, Laser Bandwidth:',ctrl%laser_bandwidth*au2eV,'eV'
          endif
        case (3)
          write(u_log,'(A)') 'Jump is in resonance with laser.'
          write(u_log,'(A,1X,F16.9,1X,A,1X,F16.9,1X,A)') &
          &'Detuning:',deltaE*au2eV,'eV, Laser Bandwidth:',ctrl%laser_bandwidth*au2eV,'eV'
        case (4)
          write(u_log,'(A)') 'Forced hop to ground state.'
          write(u_log,'(A,1X,I4,1X,A)') 'Old state:',traj%state_diag_old,'(diag)'
          write(u_log,'(A,1X,I4,1X,A)') 'New state:',traj%state_diag,'(diag)'
      endselect
      write(u_log,*)
      write(u_log,'(A,1X,F20.12)') 'trajectory weight',traj%army_ants_weight
    endif


  endsubroutine

! ===========================================================
! army ants surface switching algorithm 
  subroutine army_ants_surface_switching(traj,ctrl)
    use definitions
    use electronic
    use matrix
    use nuclear
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl
 
    integer :: istep

    integer :: istate, jstate, iatom, idir, ilaser
    real*8 :: dsum
    complex*16 :: ccoeff_diag_old_s(ctrl%nstates)
    real*8 :: randnum, cumuprob
    complex*16 :: w(ctrl%nstates,ctrl%nstates),tw(ctrl%nstates,ctrl%nstates),dw(ctrl%nstates,ctrl%nstates)
    complex*16 :: NACT(ctrl%nstates,ctrl%nstates), pNACT(ctrl%nstates,ctrl%nstates)
    real*8 :: deltaE

    real*8 :: braching_probability(ctrl%nstates)
    real*8 :: weighting(ctrl%nstates)
    integer :: allowed_branching(ctrl%nstates)
    real*8 :: branching_randum(ctrl%nstates)


    if (printlevel>2) then
      write(u_log,*) '========================================================================='
      write(u_log,*) '     Doing army ants surface switching, change of the pointer state '
      write(u_log,*) '========================================================================='
    endif

    ! calculate the switching probabilities
    select case (ctrl%switching_procedure)
      case (0)
        traj%switchprob_s=0.d0
        braching_probability=0.d0
      case (1)
        if (printlevel>2) then
          write(u_log,*) 'switching probabilities computed with CSDM'
        endif

        !==== start reinitialization
        ! compute the D, used to reinitialize
        select case (ctrl%coupling)
          case (0)  ! ddt
            ! compute magnitude of NAC in diagonal basis
            traj%dmag=0.d0
            do istate=1,ctrl%nstates
              if (istate .ne. traj%state_diag) then
                dsum=0.d0
                do iatom=1,ctrl%natom
                do idir=1,3
                  dsum=dsum+traj%NACdR_diag_ssad(istate,traj%state_diag,iatom,idir)*conjg(traj%NACdR_diag_ssad(istate,traj%state_diag,iatom,idir))
                enddo
                enddo
              endif
              dsum=max(0.d0,dsum)
              dsum=sqrt(dsum)
              traj%dmag=traj%dmag+dsum
            enddo
          case (1)  ! ddr
            ! compute magnitude of NAC in diagonal
            traj%dmag=0.d0
            do istate=1,ctrl%nstates
              if (istate .ne. traj%state_diag) then
                dsum=0.d0
                do iatom=1,ctrl%natom
                do idir=1,3
                  dsum=dsum+traj%NACdR_diag_ssad(istate,traj%state_diag,iatom,idir)*conjg(traj%NACdR_diag_ssad(istate,traj%state_diag,iatom,idir))
                enddo
                enddo
              endif
              dsum=max(0.d0,dsum)
              dsum=sqrt(dsum)
              traj%dmag=traj%dmag+dsum
            enddo
          case (2)  ! overlap
            ! compute magnitude of NACdt
            traj%dmag=0.d0
            do istate=1,ctrl%nstates
              if (istate .ne. traj%state_diag) then
                dsum=abs(traj%NACdt_ss(istate,traj%state_diag))  ! CSDM-C
                !dsum=abs(traj%overlaps_ss(istate,traj%state_diag)-traj%overlaps_ss(traj%state_diag,istate))/2
                !! HST-scheme
              endif
              traj%dmag=traj%dmag+dsum
            enddo
          case (3)  ! ktdc
            ! compute magnitude of NACdt
            traj%dmag=0.d0
            do istate=1,ctrl%nstates
              if (istate .ne. traj%state_diag) then
                dsum=abs(traj%NACdt_ss(istate,traj%state_diag))  ! CSDM-C
              endif
              traj%dmag=traj%dmag+dsum
            enddo

        endselect

        if (traj%step .gt. 2) then
          if (printlevel>4) then
            write(u_log,*) 'Magnitudes of NACdr in for current and previous two steps:'
            write(u_log,'(3(F8.5,3X))') traj%dmag, traj%dmag1, traj%dmag2
          endif
        endif

        ! reinitialize the ccoeff_diag_s based on the local minimum of dmag
        if (traj%step .gt. 2) then
          if ( (traj%dmag1-traj%dmag2) .lt. 0.d0 .AND. (traj%dmag-traj%dmag1) .gt. 0.d0 ) then
            if (printlevel>2) then
              if (ctrl%coupling .eq. 2) then
                write(u_log,*) 'Using overlap coupling, reinitialization based on CSDM-C'
                write(u_log,*) 'Coherent density has been reinitialized for CSDM-C'
              else
                write(u_log,*) 'Coherent density has been reinitialized for CSDM'
              endif
            endif
            traj%ccoeff_diag_s=traj%coeff_diag_s
            traj%ccoeff_MCH_s=traj%coeff_MCH_s
            if (ctrl%integrator==0) then
              traj%gRccoeff_MCH_s=traj%gRcoeff_MCH_s
              traj%gIccoeff_MCH_s=traj%gIcoeff_MCH_s
            endif
          endif
        endif
        traj%dmag2=traj%dmag1
        traj%dmag1=traj%dmag
        !===== end reinitialization

        !==== compute switching probability
        if (ctrl%integrator==0) then
          call calculate_probabilitiesBSH(ctrl%nstates, traj%ccoeff_diag_s, traj%preprob_old_s3, &
          &traj%preprob_s3, traj%state_diag, traj%switchprob_s)
        elseif (ctrl%integrator==1 .or. ctrl%integrator==2) then
          call calculate_probabilitiesDoM(ctrl%nstates, traj%ccoeff_diag_old_s, traj%ccoeff_diag_s, &
          &traj%Rtotal_ss, traj%state_diag, traj%switchprob_s)
        endif
        ! Calculate branching probabilities
        do istate=1,ctrl%nstates
          braching_probability(istate)=max(traj%switchprob_s(istate), traj%branching_likehood)
        enddo

      case (2)
        if (printlevel>2) then
          write(u_log,*) 'switching probabilities computed with SCDM'
        endif
        !==== compute switching probability
        if (ctrl%integrator==0) then
          call calculate_probabilitiesBSH(ctrl%nstates, traj%ccoeff_diag_s, traj%preprob_old_s3, &
          &traj%preprob_s3, traj%state_diag, traj%switchprob_s)
        elseif (ctrl%integrator==1 .or. ctrl%integrator==2) then
          call calculate_probabilitiesDoM(ctrl%nstates, traj%ccoeff_diag_old_s, traj%ccoeff_diag_s, &
          &traj%Rtotal_ss, traj%state_diag, traj%switchprob_s)
        endif
        ! Calculate branching probabilities
        do istate=1,ctrl%nstates
          braching_probability(istate)=max(traj%switchprob_s(istate), traj%branching_likehood)
        enddo

      case (3)
        if (printlevel>2) then
          write(u_log,*) 'switching probabilities computed with NDM'
        endif
        !==== compute switching probability
        if (ctrl%integrator==0) then
          call calculate_probabilitiesBSH(ctrl%nstates, traj%coeff_diag_s, traj%preprob_old_s3, &
          &traj%preprob_s3, traj%state_diag, traj%switchprob_s)
        elseif (ctrl%integrator==1 .or. ctrl%integrator==2) then
          call calculate_probabilitiesDoM(ctrl%nstates, traj%coeff_diag_old_s, traj%coeff_diag_s, & 
          &traj%RDtotal_ss, traj%state_diag, traj%switchprob_s)
        endif
        ! Calculate branching probabilities
        do istate=1,ctrl%nstates
          braching_probability(istate)=max(traj%switchprob_s(istate), traj%branching_likehood)
        enddo

      case default
        write(0,*) 'Unknown option ',ctrl%switching_procedure,' to switching_procedure!'
      stop 1
    endselect

    if (ctrl%switching_procedure==1 .or. ctrl%switching_procedure==2) then
      if (printlevel>2) then
        write(u_log,*) 'Old and new coherent occupancies and switching probabilities:'
        do istate=1,ctrl%nstates
          write(u_log,'(3(F14.9,3X))') abs(traj%ccoeff_diag_old_s(istate))**2,&
          &abs(traj%ccoeff_diag_s(istate))**2, traj%switchprob_s(istate)
        enddo
      endif
    else if (ctrl%switching_procedure==3) then 
      if (printlevel>2) then
        write(u_log,*) 'Old and new occupancies and switching probabilities:'
        do istate=1,ctrl%nstates
          write(u_log,'(3(F14.9,3X))') abs(traj%coeff_diag_old_s(istate))**2,&
          &abs(traj%coeff_diag_s(istate))**2, traj%switchprob_s(istate)
        enddo
      endif
    endif

    ! Compute current step contribution of weighting
    weighting(traj%state_diag)=1.d0
    do istate=1,ctrl%nstates
      if (istate.ne.traj%state_diag) then
        if (braching_probability(istate)>0.d0) then
          weighting(istate)=2.d0*traj%switchprob_s(istate)/braching_probability(istate)
          weighting(traj%state_diag)=weighting(traj%state_diag)-weighting(istate)/2.d0
        else
          weighting(istate)=0.d0
          weighting(traj%state_diag)=weighting(traj%state_diag)-weighting(istate)/2.d0
        endif
      endif
    enddo
    weighting(traj%state_diag)=2.d0*weighting(traj%state_diag)

    if (printlevel>2) then
      call vecwrite(ctrl%nstates, braching_probability,  u_log, 'branching probability', 'F14.9')
      call vecwrite(ctrl%nstates, weighting,  u_log, 'weighting', 'F14.9')
    endif

    ! forced switch to ground state if energy gap to state above is smaller than
    ! threshold
    traj%kind_of_jump=0
    if ( ctrl%force_hop_to_gs>0. ) then
      deltaE=abs( real(traj%H_diag_ss(traj%state_diag,traj%state_diag) - traj%H_diag_ss(1,1)) )
      if ( (traj%state_diag/=1).and.(deltaE<=ctrl%force_hop_to_gs) ) then
        traj%switchprob_s=0.
        traj%switchprob_s(1)=1.
        braching_probability=0.
        braching_probability(1)=1.
        traj%kind_of_jump=4
        traj%state_diag=1
      endif
      if (traj%state_diag==1) then
        traj%switchprob_s=0.
        braching_probability=0.
      endif

    else ! not forced hop to ground state, and check branching.

      ! Step 1.
      ! Check whether we will do a branching
      ! initialize allowed branching to 0, i.e. all states are not allowed
      allowed_branching=0
      do istate=1,ctrl%nstates
        if (istate.ne.traj%state_diag) then
          call random_number(randnum)
          branching_randum(istate)=randnum
          if (branching_randum(istate)<braching_probability(istate)) then
            allowed_branching(istate)=1
          else
            allowed_branching(istate)=0
          endif
        else
          branching_randum(istate)=1
        endif
      enddo

      ! Step2. 
      ! Check whether we will do a hopping.
      cumuprob=0.d0
      if (.not. traj%kind_of_jump==4) then
        traj%kind_of_jump=0
      endif
      deltaE=0.d0
      call random_number(randnum)
      traj%randnum=randnum
      stateloop: do istate=1,ctrl%nstates
        ! calculate cumulative probability
        cumuprob=cumuprob + 1.d0/ctrl%nstates
        if (cumuprob .ge. randnum .and. allowed_branching(istate)==1) then

        ! check for laser resonance
          if (ctrl%laser/=0) then
            do ilaser=1,ctrl%nlasers
              deltaE=abs(real(traj%H_diag_ss(istate,istate)-traj%Epot))
              deltaE=abs(deltaE-ctrl%laserenergy_tl(traj%step*ctrl%nsubsteps+1,ilaser))
              if (deltaE<=ctrl%laser_bandwidth) then
                traj%state_diag_old=traj%state_diag
                if (ctrl%switching_procedure/=0) then
                  traj%state_diag=istate
                endif
                traj%Epot=real(traj%H_diag_ss(istate,istate))
                traj%kind_of_jump=3   ! induced hop
                exit stateloop        ! ************************************************* exit of loop
              endif
            enddo
          endif

          if (printlevel>4) then
            write(u_log,'(A,1X,I4,1X,A,1X,I4,1X,A)') ' previous step selected state:',traj%state_diag, '(diag)', traj%state_MCH, '(MCH)'
          endif
          traj%state_diag_old=traj%state_diag
          traj%state_MCH_old=traj%state_MCH
          if (ctrl%switching_procedure/=0) then
            traj%state_diag=istate
            traj%state_MCH=state_diag_to_MCH(ctrl%nstates,traj%state_diag,traj%U_ss)
            if (printlevel>4) then
              write(u_log,'(A,1X,I4,1X,A,1X,I4,1X,A)') ' selected state:',traj%state_diag, '(diag)', traj%state_MCH, '(MCH)'
            endif
          endif
          if (.not. traj%kind_of_jump==4) then
            traj%kind_of_jump=1
          endif

          exit stateloop         ! ************************************************* exit of loop

        endif ! if (cumuprob .ge. randnum .and. allowed_branching(istate)==1) then
      enddo stateloop
    
    endif ! if ( ctrl%force_hop_to_gs>0. ) then

    ! finalize trajectory weight
    traj%army_ants_weight=traj%army_ants_weight*weighting(traj%state_diag)

    if (printlevel>2) then
      write(u_log,*)
      write(u_log,'(A,1X,F12.9)') 'Random number:', traj%randnum
      if (ctrl%switching_procedure==0) then
        write(u_log,*)
        write(u_log,'(A,1X,F12.9)') '*** Surface Switching is turned off: new state = old state ***'
        write(u_log,*)
      endif
      select case (traj%kind_of_jump)
        case (0)
          write(u_log,*) 'No jump performed.'
        case (1)
          write(u_log,*) 'Pointer state switched'
          write(u_log,'(A,1X,I4,1X,A,1X,I4,1X,A)') 'Old pointer state:',traj%state_diag_old,'(diag)',traj%state_MCH_old,'(MCH)'
          write(u_log,'(A,1X,I4,1X,A,1X,I4,1X,A)') 'New pointer state:',traj%state_diag,'(diag)',traj%state_MCH,'(MCH)'
        case (3)
          write(u_log,'(A)') 'Jump is in resonance with laser.'
          write(u_log,'(A,1X,F16.9,1X,A,1X,F16.9,1X,A)') &
          &'Detuning:',deltaE*au2eV,'eV, Laser Bandwidth:',ctrl%laser_bandwidth*au2eV,'eV'
        case (4)
          write(u_log,'(A)') 'Forced switch to ground state.'
          write(u_log,'(A,1X,I4,1X,A)') 'Old pointer state:',traj%state_diag_old,'(diag)'
          write(u_log,'(A,1X,I4,1X,A)') 'New pointer state:',traj%state_diag,'(diag)'
      endselect
      write(u_log,*)
      write(u_log,'(A,1X,F20.12)') 'trajectory weight',traj%army_ants_weight
    endif


  endsubroutine





endmodule
