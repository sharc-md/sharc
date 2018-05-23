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

!> # Module DEFINITIONS
!> \author Sebastian Mai
!> \date 10.07.2014, modified 27.02.2017 by Philipp Marquetand
!> 
!> This module defines the trajectory and control types.
!> 
!> All arrays defined here have their order and meaning of the
!> indices in the last part of the name, e.g.:
!> mass_a is the array containing the masses of all atoms,
!> where the index runs over the atoms.
!>
!> Indices:
!> - _a Atoms
!> - _d Cartesian direction (1=x,2=y,3=z)
!> - _s Electronic states
!> - _m Multiplicities (1=singlets, 2=doublets, triplets, ...)
!> - _t Timesteps
!>
module definitions
implicit none
save

!> # Auxiliary trajectory type:
!> This is a type with the data relevant for auxiliary trajectories
!> needed for the decoherence correction in A-FSSH.
!> geom, veloc, and accel are given relative to the main trajectory
type aux_trajectory_type
  sequence             ! to store the constituents of the type contiguously

  ! nuclear information
  real*8,allocatable :: mass_a(:)                        !< atomic mass in a.u. (1 a.u. = rest mass of electron m_e)
  real*8,allocatable :: geom_ad(:,:)                     !< Cartesian displacement from SH trajectory in a.u. (bohr)
  real*8,allocatable :: veloc_ad(:,:)                    !< Cartesian displacment velocity in a.u. (bohr/atu)
  real*8,allocatable :: accel_ad(:,:)                    !< Cartesian displacment acceleration in a.u. (bohr/atu/atu)
  real*8,allocatable :: grad_ad(:,:)                     !< Difference gradient
  real*8,allocatable :: geom_tmp_ad(:,:)                    !< Temporary Cartesian displacement during A-FSSH step
  real*8,allocatable :: veloc_tmp_ad(:,:)                   !< Temporary Cartesian displacment velocity during A-FSSH step

  ! other quantities
  integer :: istate                                       !< Which state
  real*8 :: rate1, rate2                                  !< Term occuring in the decoherence rate
endtype

! =========================================================== !

public
!> # Trajectory type:
!> This is a type with all data which would be private to a trajectory,
!> if several trajectories would be handled at the same time.
!> This data includes the following:
!> - General information (timestep, RNGseed, traj_hash)
!> - The current and old state in MCH and diag representation
!> - Energies
!> - Information about the atoms (mass, coordinate, velocity, acceleration)
!> - Information about the electronic states (for current and last step):
!>   - Hamiltonian
!>   - Transformation matrix
!>   - non-adiabatic couplings/overlaps
!>   - (Transition) Dipole moments
!>   - Property matrices
!>   - Propagator
!>   - Hopping probabilities
!> - Gradients and non-adiabatic coupling vectors
!> - Selection masks
type trajectory_type
  sequence             ! to store the constituents of the type contiguously

  ! general information
  integer :: RNGseed                                     !< seed for random nuber generator
  integer :: step                                        !< current timestep (step=0 => t=0)
  integer*8 :: traj_hash                                 !< Trajectory ID based on hashing the input files

  integer :: state_MCH                                   !< currently occupied state in MCH basis
  integer :: state_diag                                  !< currently occupied state in diag basis
  integer :: state_diag_old                              !< diag state occupied in the last timestep
  integer :: state_diag_frust                            !< diag state of a frustrated hop

  real*8 :: Ekin                                         !< kinetic energy
  real*8 :: Epot                                         !< potential energy (diag energy of state_diag)
  real*8 :: Etot                                         !< total energy, Ekin+Epot

  integer :: time_start                                  !< system time when sharc is started
  integer :: time_last                                   !< system time after completion of the previous timestep
  integer :: time_step                                   !< system time(after timestep) - system time(before timestep)
  integer :: kind_of_jump                                !< 0=no jump, 1=regular, 2=frustrated, 3=resonant
  integer :: steps_in_gs                                 !< counter for the number of timesteps in the lowest state
  logical :: phases_found                                !< whether wavefunction phases were found in QM.out

  ! nuclear information
  real*8,allocatable :: atomicnumber_a(:)                !< atomic number
  character*2,allocatable :: element_a(:)                !< element descriptor
  real*8,allocatable :: mass_a(:)                        !< atomic mass in a.u. (1 a.u. = rest mass of electron m_e)
  real*8,allocatable :: geom_ad(:,:)                     !< Cartesian coordinates of atom in a.u. (bohr)
  real*8,allocatable :: veloc_ad(:,:)                    !< Cartesian velocity in a.u. (bohr/atu)
  real*8,allocatable :: accel_ad(:,:)                    !< Cartesian acceleration in a.u. (bohr/atu/atu)

  ! electronic information
  complex*16,allocatable :: H_MCH_ss(:,:)                !< MCH Hamiltonian as read from QM.out (no laser)
                                                         !< Laser interaction is added during propagation
  complex*16,allocatable :: dH_MCH_ss(:,:)               !< time derivative of MCH Hamiltonian
  complex*16,allocatable :: H_MCH_old_ss(:,:)            !< MCH Hamiltonian of last timestep (no laser)
  complex*16,allocatable :: H_diag_ss(:,:)               !< diag Hamiltonian
  complex*16,allocatable :: U_ss(:,:)                    !< transformation matrix
  complex*16,allocatable :: U_old_ss(:,:)                !< transformation matrix of last timestep
  complex*16,allocatable :: NACdt_ss(:,:)                !< time-derivatives of wavefunctions
  complex*16,allocatable :: NACdt_old_ss(:,:)            !< time-derivatives of wavefunctions of last timestep
  complex*16,allocatable :: overlaps_ss(:,:)             !< overlaps for LD propagation
  complex*16,allocatable :: DM_ssd(:,:,:)                !< (transition) dipole moment matrix 
                                                         !< transition dipoles between active and inactive states are zero.
  complex*16,allocatable :: DM_old_ssd(:,:,:)            !< old dipole moment matrix
  complex*16,allocatable :: DM_print_ssd(:,:,:)          !< dipole moment matrix used for the output routines
                                                         !< transition dipoles between active and inactive states are not zero.
  complex*16,allocatable :: Rtotal_ss(:,:)               !< total propagator for the current timestep
  complex*16,allocatable :: phases_s(:)                  !< electronic state phases of the current step
  complex*16,allocatable :: phases_old_s(:)              !< electronic state phases of the last step
  real*8, allocatable :: hopprob_s(:)                    !< hopping probabilities
  real*8 :: randnum, randnum2                            !< random number for surface hopping and A-FSSH

  ! arbitrary properties
  complex*16,allocatable :: Property2d_xss(:,:,:)        !< list of matrices containing arbitrary data (not used in propagation)
  real*8,    allocatable :: Property1d_ys(:,:)           !< list of vectors containing arbitrary data (not used in propagation)
  character*40,allocatable :: Property2d_labels_x(:)     !< list of descriptions for 2d properties
  character*40,allocatable :: Property1d_labels_y(:)     !< list of descriptions for 1d properties

  ! vector information
  real*8,allocatable :: DMgrad_ssdad(:,:,:,:,:)          !< Cartesian gradient of the dipole moments (bra, ket, polarization, atom, cartesian component of atom displacement)
  real*8,allocatable :: NACdR_ssad(:,:,:,:)              !< vectorial non-adiabatic couplings in a.u.
  real*8,allocatable :: NACdR_old_ssad(:,:,:,:)          !< vectorial non-adiabatic couplings of last timestep
  real*8,allocatable :: grad_MCH_sad(:,:,:)              !< Cartesian gradient in a.u (hartree/bohr) of all states
  complex*16,allocatable :: Gmatrix_ssad(:,:,:,:)        !< Cartesian gradient in a.u (hartree/bohr) of the current state in diag basis
  real*8, allocatable :: grad_ad(:,:)                    !< final gradient used in velocity verlet

  ! coefficient information
  complex*16,allocatable :: coeff_diag_s(:)              !< coefficients of electronic wavefunction in diag representation
  complex*16,allocatable :: coeff_diag_old_s(:)          !< coefficients of electronic wavefunction (of previous timestep) in diag representation
  complex*16,allocatable :: coeff_MCH_s(:)               !< coefficients of electronic wavefunction in MCH representation

  ! gradient and nac selection information
  logical,allocatable :: selG_s(:)                       !< selection mask for gradients
  logical,allocatable :: selT_ss(:,:)                    !< selection mask for non-adiabatic coupling vectors
  logical,allocatable :: selDM_ss(:,:)                   !< selection mask for dipole moment gradients

  ! Auxiliary trajectories for A-FSSH
  type(aux_trajectory_type),allocatable :: auxtrajs_s(:)

endtype

! =========================================================== !

!> # Control type:
!> This is a type containing data which are shared between trajectories of the ensemble.
!>This data includes the following:
!> - Number of atoms, states, multiplicities, active states
!> - number of timesteps, substeps, length of timesteps
!> - Energy shift, scaling factor, damping factor
!> - energy-based-selection thresholds, decoherence parameter
!> - method switches (SHARC/regular, coupling type, laser, quantities to calculate, ...)
type ctrl_type
  sequence

  character*1023 :: cwd                     !< working directory for SHARC
  real*8 :: output_version                         !< version as float for checks during writing output

! numerical constants
  integer :: natom                          !< number of atoms
  logical,allocatable :: atommask_a(:)      !< atoms which are considered for decoherence, rescaling, ...
  integer :: maxmult                        !< highest spin quantum number (determines length of nstates_m)
  integer,allocatable :: nstates_m(:)       !< numer of states considered in each multiplicy
  integer :: nstates                        !< total number of states
  integer :: nsteps                         !< total number of simulation steps 
  integer :: nsubsteps                      !< number of steps for the electron propagation
  real*8 :: dtstep                          !< length of timestep in a.u (atu)
  real*8 :: ezero                           !< energy offset in a.u. (e.g. ground state equilibrium energy)
  real*8 :: scalingfactor                   !< scales the Hamiltonian and gradients
  real*8 :: eselect_grad                    !< energy difference for neglecting gradients
  real*8 :: eselect_nac                     !< energy difference for neglecting na-couplings
  real*8 :: eselect_dmgrad                  !< energy difference for neglecting dipole gradients
  real*8 :: dampeddyn                       !< damping factor for kinetic energy
  real*8 :: decoherence_alpha               !< decoherence parameter (a.u.) for energy-based decoherence
  logical,allocatable :: actstates_s(:)     !< mask of the active states

! methods and switches
  logical :: restart                        !< restart yes or no
  integer :: staterep                       !< 0=initial state is given in diag representation, 1=in MCH representation
  integer :: initcoeff                      !< 0=initial coefficients are diag, 1=initial coefficients are MCH, 2=auto diag, 3=auto MCH
  integer :: laser                          !< 0=none, 1=internal, 2=external
  integer :: coupling                       !< 0=ddt, 1=ddr, 2=overlap
  integer :: surf                           !< 0=propagation in diag surfaces (SHARC), 1=on MCH surfaces (regular SH)
  integer :: decoherence                    !< 0=off, 1=activate energy-based decoherence correction
  integer :: ekincorrect                    !< 0=none, 1=adjust momentum along velocity, 2=adjust momentum along nac vector
  integer :: reflect_frustrated             !< 0=none, 1=reflect along velocity, 2=reflect along nac vector
  integer :: gradcorrect                    !< 0=no, 1=include nac vectors in gradient transformation
  integer :: dipolegrad                     !< 0=no, 1=include dipole gradients in gradient transformation

  integer :: calc_soc                       !< request SOC, otherwise only the diagonal elements of H (plus any laser interactions) are taken into account\n 0=no soc, 1=soc enabled
  integer :: calc_grad                      !< request gradients:   \n        0=all in step 1, 1=select in step 1, 2=select in step 2
  integer :: calc_overlap                   !< 0=no, 1=request overlap matrices
  integer :: calc_nacdt                     !< 0=no, 1=request time derivatives
  integer :: calc_nacdr                     !< request nac vectors: \n -1=no, 0=all in step 1, 1=select in step 1, 2=select in step 2
  integer :: calc_dipolegrad                !< request dipole gradient vectors: \n -1=no, 0=all in step 1, 1=select in step 1, 2=select in step 2
  integer :: calc_second                    !< 0=no, 1=do two interface calls per timestep
  integer :: calc_phases                    !< 0=nom 1=yes

  integer :: write_soc                      !< write SOC to output.dat or not \n 0=no soc, 1=write soc )
  integer :: write_grad                     !< write gradients:   \n        0=no gradients, 1=write gradients
  integer :: write_overlap                  !< write overlap matrix:   \n        0=no overlap, 1=write overlap
  integer :: write_NACdr                    !< write nac vectors:   \n        0=no vectors, 1=write vectors

  integer :: write_property2d               !< write property matrices:   \n        0=no property, 1=write property
  integer :: write_property1d               !< write property vectors:   \n        0=no property, 1=write property
  integer :: n_property2d                   !< number of property matrices
  integer :: n_property1d                   !< number of property vectors

  integer :: killafter                      !< -1=no, >1=kill after that many steps in the ground state
  integer :: ionization                     !<  -1=no, n=ionization every n steps
  integer :: theodore                       !<  -1=no, n=theodore every n steps
  integer :: track_phase                    !< 0=no, 1=track phase of U matrix through the propagation (turn off only for debugging purposes)
  integer :: track_phase_at_zero            !< 0=nothing, 1=at time zero, get phases from whatever is in the savedir
  integer :: hopping_procedure              !< 0=no hops, 1=hops (standard formula), 2=GFSH

! thresholds
!   real*8 :: propag_sharc_UdUdiags=1.d-2           ! Threshold for the size of diagonal elements in UdU (needed for dynamic substeps)        in hartree
!   real*8 :: min_dynamic_substep=1.d-5             ! In dynamic substepping, the shortest substep allowed                                        in atomic time units
!   real*8 :: diagonalize_degeneracy_diff=1.d-9     ! Energy difference threshold for treating states as degenerate                                in hartree

  ! only array in ctrl
  real*8 :: laser_bandwidth                       !< for detecting induced hops (in a.u.)
  integer :: nlasers
  complex*16, allocatable :: laserfield_td(:,:)   !< complex valued laser field
  complex*16, allocatable :: laserenergy_tl(:,:)     !< momentary central energy of laser (for detecting induced hops)

endtype

! =========================================================== !
integer :: printlevel
!< verbosity of the log file
!< -0=build and execution info (hostname, date, cwd, compiler)
!< -1=+ internal steps
!< -2=+ input parsing infos
!< -3 and higher=+ print various numerical values per timestep
! =========================================================== !

real*8,parameter:: au2a=0.529177211d0             !< length
real*8,parameter:: au2fs=0.024188843d0            !< time 
real*8,parameter:: au2u=5.4857990943d-4           !< mass
real*8,parameter:: au2rcm=219474.631370d0         !< energy
real*8,parameter:: au2eV=27.21138386d0            !< energy
real*8,parameter:: au2debye=2.5417469d0           !< dipole moment

complex*16,parameter:: ii=dcmplx(0.d0,1.d0)       !< imaginary unit
real*8,parameter:: pi=4.d0*datan(1.d0)            !< pi

character*20,parameter :: multnames(8)=(/'Singlet','Doublet','Triplet','Quartet','Quintet',' Sextet',' Septet','  Octet'/)
!< strings used to represent the multiplicities
! =========================================================== !

character*255, parameter :: version='2.0 (March 13, 2017)'    !< string holding the version number

integer, parameter :: u_log=1                !< long output file
integer, parameter :: u_lis=2                !< short output file
integer, parameter :: u_dat=3                !< compressed data output file
integer, parameter :: u_geo=4                !< geometry output file
integer, parameter :: u_resc=7               !< restart file ctrl
integer, parameter :: u_rest=8               !< restart file traj
! 
integer, parameter :: u_i_input=12           !< trajectory input (control variables, initial state, ...)
integer, parameter :: u_i_geom=13            !< initial geometry
integer, parameter :: u_i_veloc=14           !< initial velocity
integer, parameter :: u_i_coeff=15           !< initial coefficients
integer, parameter :: u_i_laser=16           !< numerical laser field
integer, parameter :: u_i_atommask=17        !< which atoms are active for rescaling/decoherence/...

integer, parameter :: u_qm_QMin=41           !< here SHARC writes information for the QM interface (like geometry, number of states, what kind of data is requested)
integer, parameter :: u_qm_QMout=42          !< here SHARC retrieves the results of the QM run (Hamiltonian, gradients, couplings, etc.)

! =========================================================== !

  contains

! =========================================================== !

    subroutine allocate_traj(traj,ctrl)
      !< Allocates all arrays in traj
      !< Does not allocate arrays in ctrl (laser-related, actstates_s, nstates_m)
      !< Reads natom and nstates from ctrl
      !< Initializes all elements of all arrays to -123, 'Q' or .true.
      implicit none
      type(ctrl_type), intent(inout) :: ctrl
      type(trajectory_type), intent(inout) :: traj
      integer :: status
      integer :: natom,nstates,n1d,n2d

      natom=ctrl%natom
      nstates=ctrl%nstates
      n1d=ctrl%n_property1d
      n2d=ctrl%n_property2d

      allocate(traj%atomicnumber_a(natom),stat=status)
      if (status/=0) stop 'Could not allocate atomicnumber_a'
      traj%atomicnumber_a=-123.d0

      allocate(traj%element_a(natom),stat=status)
      if (status/=0) stop 'Could not allocate element_a'
      traj%element_a='Q'

      allocate(traj%mass_a(natom),stat=status)
      if (status/=0) stop 'Could not allocate mass_a'
      traj%mass_a=-123.d0


      allocate(traj%geom_ad(natom,3),stat=status)
      if (status/=0) stop 'Could not allocate geom_ad'
      traj%geom_ad=-123.d0

      allocate(traj%veloc_ad(natom,3),stat=status)
      if (status/=0) stop 'Could not allocate veloc_ad'
      traj%veloc_ad=-123.d0

      allocate(traj%accel_ad(natom,3),stat=status)
      if (status/=0) stop 'Could not allocate accel_ad'
      traj%accel_ad=-123.d0


      allocate(traj%H_MCH_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate H_MCH_ss'
      traj%H_MCH_ss=-123.d0

      allocate(traj%dH_MCH_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate dH_MCH_ss'
      traj%dH_MCH_ss=-123.d0

      allocate(traj%H_MCH_old_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate H_MCH_old_ss'
      traj%H_MCH_old_ss=-123.d0

      allocate(traj%H_diag_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate H_diag_ss'
      traj%H_diag_ss=-123.d0

      allocate(traj%U_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate U_ss'
      traj%U_ss=-123.d0

      allocate(traj%U_old_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate U_old_ss'
      traj%U_old_ss=-123.d0

      allocate(traj%NACdt_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate NACdt_ss'
      traj%NACdt_ss=-123.d0

      allocate(traj%NACdt_old_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate NACdt_old_ss'
      traj%NACdt_old_ss=-123.d0

      allocate(traj%overlaps_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate overlaps_ss'
      traj%overlaps_ss=-123.d0

      allocate(traj%DM_ssd(nstates,nstates,3),stat=status)
      if (status/=0) stop 'Could not allocate DM_ssd'
      traj%DM_ssd=-123.d0

      allocate(traj%DM_old_ssd(nstates,nstates,3),stat=status)
      if (status/=0) stop 'Could not allocate DM_old_ssd'
      traj%DM_old_ssd=-123.d0

      allocate(traj%DM_print_ssd(nstates,nstates,3),stat=status)
      if (status/=0) stop 'Could not allocate DM_print_ssd'
      traj%DM_print_ssd=-123.d0

      allocate(traj%Rtotal_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate Rtotal_ss'
      traj%Rtotal_ss=-123.d0

      allocate(traj%hopprob_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate hopprob_s'
      traj%hopprob_s=-123.d0


      allocate(traj%Property2d_xss(n2d,nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate Property2d_xss'
      traj%Property2d_xss=-123.d0

      allocate(traj%Property1d_ys(n1d,nstates),stat=status)
      if (status/=0) stop 'Could not allocate Property1d_ys'
      traj%Property1d_ys=-123.d0

      allocate(traj%Property2d_labels_x(n2d),stat=status)
      if (status/=0) stop 'Could not allocate Property2d_labels_x'
      traj%Property2d_labels_x='N/A'

      allocate(traj%Property1d_labels_y(n1d),stat=status)
      if (status/=0) stop 'Could not allocate Property1d_labels_y'
      traj%Property1d_labels_y='N/A'


      allocate(traj%phases_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate phases_s'
      traj%phases_s=-123.d0

      allocate(traj%phases_old_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate phases_old_s'
      traj%phases_old_s=-123.d0

      allocate(traj%coeff_diag_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate coeff_diag_s'
      traj%coeff_diag_s=-123.d0

      allocate(traj%coeff_diag_old_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate coeff_old_diag_s'
      traj%coeff_diag_old_s=-123.d0

      allocate(traj%coeff_MCH_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate coeff_MCH_s'
      traj%coeff_MCH_s=-123.d0

      allocate(traj%selG_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate selG_s'
      traj%selG_s=.true.

      allocate(traj%selT_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate selT_ss'
      traj%selT_ss=.true.

      allocate(traj%selDM_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate selDM_ss'
      traj%selDM_ss=.true.


      allocate(traj%DMgrad_ssdad(nstates,nstates,3,natom,3),stat=status)
      if (status/=0) stop 'Could not allocate DMgrad_ssdad'
      traj%DMgrad_ssdad=-123.d0

      allocate(traj%NACdR_ssad(nstates,nstates,natom,3),stat=status)
      if (status/=0) stop 'Could not allocate NACdR_ssad'
      traj%NACdR_ssad=-123.d0

      allocate(traj%NACdR_old_ssad(nstates,nstates,natom,3),stat=status)
      if (status/=0) stop 'Could not allocate NACdR_old_ssad'
      traj%NACdR_old_ssad=-123.d0

      allocate(traj%grad_MCH_sad(nstates,natom,3),stat=status)
      if (status/=0) stop 'Could not allocate gradMCH_sad'
      traj%grad_MCH_sad=-123.d0

      allocate(traj%grad_ad(natom,3),stat=status)
      if (status/=0) stop 'Could not allocate grad_ad'
      traj%grad_ad=-123.d0

      allocate(traj%Gmatrix_ssad(nstates,nstates,natom,3),stat=status)
      if (status/=0) stop 'Could not allocate Gmatrix_ssad'
      traj%Gmatrix_ssad=-123.d0

    endsubroutine

! =========================================================== !
  !> checks whether all members of ctrl and traj are allocated
  !> also checks for NaNs
  !> for debugging purposes, currently not used
  !> \todo add printlevel
  !> \param u Unit to which the report should be printed
  subroutine check_allocation(u,ctrl,traj)

    implicit none
    integer,intent(in) :: u
    type(ctrl_type), intent(inout) :: ctrl
    type(trajectory_type), intent(inout) :: traj

    write(u,*) '________________ CHECKING ISALLOCATED ___________________'

    write(u,'(A20,1X,L1)') 'atomicnumber_a',  allocated(traj%atomicnumber_a  )
    write(u,'(A20,1X,L1)') 'element_a',       allocated(traj%element_a       )
    write(u,'(A20,1X,L1)') 'mass_a',          allocated(traj%mass_a          )
    write(u,'(A20,1X,L1)') 'geom_ad',         allocated(traj%geom_ad         )
    write(u,'(A20,1X,L1)') 'veloc_ad',        allocated(traj%veloc_ad        )
    write(u,'(A20,1X,L1)') 'accel_ad',        allocated(traj%accel_ad        )
    write(u,'(A20,1X,L1)') 'H_MCH_ss',        allocated(traj%H_MCH_ss        )
    write(u,'(A20,1X,L1)') 'dH_MCH_ss',       allocated(traj%dH_MCH_ss       )
    write(u,'(A20,1X,L1)') 'H_MCH_old_ss',    allocated(traj%H_MCH_old_ss    )
    write(u,'(A20,1X,L1)') 'H_diag_ss',       allocated(traj%H_diag_ss       )
    write(u,'(A20,1X,L1)') 'U_ss',            allocated(traj%U_ss            )
    write(u,'(A20,1X,L1)') 'U_old_ss',        allocated(traj%U_old_ss        )
    write(u,'(A20,1X,L1)') 'NACdt_ss',        allocated(traj%NACdt_ss        )
    write(u,'(A20,1X,L1)') 'NACdt_old_ss',    allocated(traj%NACdt_old_ss    )
    write(u,'(A20,1X,L1)') 'NACdR_ssad',      allocated(traj%NACdR_ssad      )
    write(u,'(A20,1X,L1)') 'NACdR_old_ssad',  allocated(traj%NACdR_old_ssad  )
    write(u,'(A20,1X,L1)') 'overlaps_ss',     allocated(traj%overlaps_ss     )
    write(u,'(A20,1X,L1)') 'DM_ssd',          allocated(traj%DM_ssd          )
    write(u,'(A20,1X,L1)') 'DM_old_ssd',      allocated(traj%DM_old_ssd      )
    write(u,'(A20,1X,L1)') 'DM_print_ssd',    allocated(traj%DM_print_ssd    )
    write(u,'(A20,1X,L1)') 'DM_ssd',          allocated(traj%DM_ssd          )
    write(u,'(A20,1X,L1)') 'Rtotal_ss',       allocated(traj%Rtotal_ss       )
    write(u,'(A20,1X,L1)') 'phases_s',        allocated(traj%phases_s        )
    write(u,'(A20,1X,L1)') 'phases_old_s',    allocated(traj%phases_old_s    )
    write(u,'(A20,1X,L1)') 'hopprob_s',       allocated(traj%hopprob_s       )
    write(u,'(A20,1X,L1)') 'Property2d_xss',  allocated(traj%Property2d_xss     )
    write(u,'(A20,1X,L1)') 'Property1d_ys',   allocated(traj%Property1d_ys      )
    write(u,'(A20,1X,L1)') 'Property2d_labels_x',     allocated(traj%Property2d_labels_x     )
    write(u,'(A20,1X,L1)') 'Property1d_labels_y',     allocated(traj%Property1d_labels_y     )
    write(u,'(A20,1X,L1)') 'grad_MCH_sad',    allocated(traj%grad_MCH_sad    )
    write(u,'(A20,1X,L1)') 'Gmatrix_ssad',    allocated(traj%Gmatrix_ssad    )
    write(u,'(A20,1X,L1)') 'grad_ad',         allocated(traj%grad_ad         )
    write(u,'(A20,1X,L1)') 'coeff_diag_s',    allocated(traj%coeff_diag_s    )
    write(u,'(A20,1X,L1)') 'coeff_diag_old_s',allocated(traj%coeff_diag_old_s)
    write(u,'(A20,1X,L1)') 'coeff_MCH_s',     allocated(traj%coeff_MCH_s     )
    write(u,'(A20,1X,L1)') 'selG_s',          allocated(traj%selG_s          )
    write(u,'(A20,1X,L1)') 'selT_ss',         allocated(traj%selT_ss         )
    write(u,'(A20,1X,L1)') 'nstates_m',       allocated(ctrl%nstates_m       )
    write(u,'(A20,1X,L1)') 'actstates_s',     allocated(ctrl%actstates_s     )

    write(u,*) '_______________________ CHECKING NaNs _______________________'

    write(u,'(A20,1X,L1)') 'mass_a',          any((traj%mass_a          ).ne.(traj%mass_a          ))
    write(u,'(A20,1X,L1)') 'geom_ad',         any((traj%geom_ad         ).ne.(traj%geom_ad         ))
    write(u,'(A20,1X,L1)') 'veloc_ad',        any((traj%veloc_ad        ).ne.(traj%veloc_ad        ))
    write(u,'(A20,1X,L1)') 'accel_ad',        any((traj%accel_ad        ).ne.(traj%accel_ad        ))
    write(u,'(A20,1X,L1)') 'NACdR_ssad',      any((traj%NACdR_ssad      ).ne.(traj%NACdR_ssad      ))
    write(u,'(A20,1X,L1)') 'NACdR_old_ssad',  any((traj%NACdR_old_ssad  ).ne.(traj%NACdR_old_ssad  ))
    write(u,'(A20,1X,L1)') 'hopprob_s',       any((traj%hopprob_s       ).ne.(traj%hopprob_s       ))
    write(u,'(A20,1X,L1)') 'grad_MCH_sad',    any((traj%grad_MCH_sad    ).ne.(traj%grad_MCH_sad    ))
    write(u,'(A20,1X,L1)') 'grad_ad',         any((traj%grad_ad         ).ne.(traj%grad_ad         ))
    write(u,*) 'Real parts:'
    write(u,'(A20,1X,L1)') 'H_MCH_ss',        any((real(traj%H_MCH_ss        )).ne.(real(traj%H_MCH_ss        )))
    write(u,'(A20,1X,L1)') 'dH_MCH_ss',       any((real(traj%dH_MCH_ss       )).ne.(real(traj%dH_MCH_ss       )))
    write(u,'(A20,1X,L1)') 'H_MCH_old_ss',    any((real(traj%H_MCH_old_ss    )).ne.(real(traj%H_MCH_old_ss    )))
    write(u,'(A20,1X,L1)') 'H_diag_ss',       any((real(traj%H_diag_ss       )).ne.(real(traj%H_diag_ss       )))
    write(u,'(A20,1X,L1)') 'U_ss',            any((real(traj%U_ss            )).ne.(real(traj%U_ss            )))
    write(u,'(A20,1X,L1)') 'U_old_ss',        any((real(traj%U_old_ss        )).ne.(real(traj%U_old_ss        )))
    write(u,'(A20,1X,L1)') 'NACdt_ss',        any((real(traj%NACdt_ss        )).ne.(real(traj%NACdt_ss        )))
    write(u,'(A20,1X,L1)') 'NACdt_old_ss',    any((real(traj%NACdt_old_ss    )).ne.(real(traj%NACdt_old_ss    )))
    write(u,'(A20,1X,L1)') 'overlaps_ss',     any((real(traj%overlaps_ss     )).ne.(real(traj%overlaps_ss     )))
    write(u,'(A20,1X,L1)') 'DM_ssd',          any((real(traj%DM_ssd          )).ne.(real(traj%DM_ssd          )))
    write(u,'(A20,1X,L1)') 'DM_old_ssd',      any((real(traj%DM_old_ssd      )).ne.(real(traj%DM_old_ssd      )))
    write(u,'(A20,1X,L1)') 'Property2d_xss',  any((real(traj%Property2d_xss  )).ne.(real(traj%Property2d_xss  )))
    write(u,'(A20,1X,L1)') 'Property1d_ys',   any((real(traj%Property1d_ys   )).ne.(real(traj%Property1d_ys   )))
    write(u,'(A20,1X,L1)') 'DM_print_ssd',    any((real(traj%DM_print_ssd    )).ne.(real(traj%DM_print_ssd    )))
    write(u,'(A20,1X,L1)') 'Rtotal_ss',       any((real(traj%Rtotal_ss       )).ne.(real(traj%Rtotal_ss       )))
    write(u,'(A20,1X,L1)') 'phases_s',        any((real(traj%phases_s        )).ne.(real(traj%phases_s        )))
    write(u,'(A20,1X,L1)') 'phases_old_s',    any((real(traj%phases_old_s    )).ne.(real(traj%phases_old_s    )))
    write(u,'(A20,1X,L1)') 'Gmatrix_ssad',    any((real(traj%Gmatrix_ssad    )).ne.(real(traj%Gmatrix_ssad    )))
    write(u,'(A20,1X,L1)') 'coeff_diag_s',    any((real(traj%coeff_diag_s    )).ne.(real(traj%coeff_diag_s    )))
    write(u,'(A20,1X,L1)') 'coeff_diag_old_s',any((real(traj%coeff_diag_old_s)).ne.(real(traj%coeff_diag_old_s)))
    write(u,'(A20,1X,L1)') 'coeff_MCH_s',     any((real(traj%coeff_MCH_s     )).ne.(real(traj%coeff_MCH_s     )))
    write(u,*) 'Imag parts:'
    write(u,'(A20,1X,L1)') 'H_MCH_ss',        any((aimag(traj%H_MCH_ss        )).ne.(aimag(traj%H_MCH_ss        )))
    write(u,'(A20,1X,L1)') 'dH_MCH_ss',       any((aimag(traj%dH_MCH_ss       )).ne.(aimag(traj%dH_MCH_ss       )))
    write(u,'(A20,1X,L1)') 'H_MCH_old_ss',    any((aimag(traj%H_MCH_old_ss    )).ne.(aimag(traj%H_MCH_old_ss    )))
    write(u,'(A20,1X,L1)') 'H_diag_ss',       any((aimag(traj%H_diag_ss       )).ne.(aimag(traj%H_diag_ss       )))
    write(u,'(A20,1X,L1)') 'U_ss',            any((aimag(traj%U_ss            )).ne.(aimag(traj%U_ss            )))
    write(u,'(A20,1X,L1)') 'U_old_ss',        any((aimag(traj%U_old_ss        )).ne.(aimag(traj%U_old_ss        )))
    write(u,'(A20,1X,L1)') 'NACdt_ss',        any((aimag(traj%NACdt_ss        )).ne.(aimag(traj%NACdt_ss        )))
    write(u,'(A20,1X,L1)') 'NACdt_old_ss',    any((aimag(traj%NACdt_old_ss    )).ne.(aimag(traj%NACdt_old_ss    )))
    write(u,'(A20,1X,L1)') 'overlaps_ss',     any((aimag(traj%overlaps_ss     )).ne.(aimag(traj%overlaps_ss     )))
    write(u,'(A20,1X,L1)') 'DM_ssd',          any((aimag(traj%DM_ssd          )).ne.(aimag(traj%DM_ssd          )))
    write(u,'(A20,1X,L1)') 'DM_old_ssd',      any((aimag(traj%DM_old_ssd      )).ne.(aimag(traj%DM_old_ssd      )))
    write(u,'(A20,1X,L1)') 'DM_print_ssd',    any((aimag(traj%DM_print_ssd    )).ne.(aimag(traj%DM_print_ssd    )))
    write(u,'(A20,1X,L1)') 'Property2d_xss',  any((aimag(traj%Property2d_xss  )).ne.(aimag(traj%Property2d_xss  )))
!     write(u,'(A20,1X,L1)') 'Property1d_ys',   any((aimag(traj%Property1d_ys   )).ne.(aimag(traj%Property1d_ys   )))
    write(u,'(A20,1X,L1)') 'Rtotal_ss',       any((aimag(traj%Rtotal_ss       )).ne.(aimag(traj%Rtotal_ss       )))
    write(u,'(A20,1X,L1)') 'phases_s',        any((aimag(traj%phases_s        )).ne.(aimag(traj%phases_s        )))
    write(u,'(A20,1X,L1)') 'phases_old_s',    any((aimag(traj%phases_old_s    )).ne.(aimag(traj%phases_old_s    )))
    write(u,'(A20,1X,L1)') 'Gmatrix_ssad',    any((aimag(traj%Gmatrix_ssad    )).ne.(aimag(traj%Gmatrix_ssad    )))
    write(u,'(A20,1X,L1)') 'coeff_diag_s',    any((aimag(traj%coeff_diag_s    )).ne.(aimag(traj%coeff_diag_s    )))
    write(u,'(A20,1X,L1)') 'coeff_diag_old_s',any((aimag(traj%coeff_diag_old_s)).ne.(aimag(traj%coeff_diag_old_s)))
    write(u,'(A20,1X,L1)') 'coeff_MCH_s',     any((aimag(traj%coeff_MCH_s     )).ne.(aimag(traj%coeff_MCH_s     )))

    write(u,*) '____________________________________________________________'

  endsubroutine

! =========================================================== !

endmodule
















