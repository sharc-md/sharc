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

!> # Module DEFINITIONS
!> \author Sebastian Mai
!> \date 10.07.2014, modified 27.02.2017 by Philipp Marquetand
!>                   modified 19.04.2018 by Maximilian F.S.J. Menger:
!>                        traj%H_MCH_ss, traj%DM_ssd
!>                        traj%overlaps_ss, traj%grad_MCH_sad
!>                   are now as well pointer tragets
!>
!>                   modified 13.11.2019 by Yinan Shu
!>                        add the following new varibales:
!>                        traj%microtime - record of trajecotry time
!>                        traj%Ekin_old, traj%Epot_old, traj%Etot_old - used by adaptive step size
!>                        traj%RDtotal_ss, traj%Dtotal_ss - CSDM propagators
!>                        traj%switchprob_s - CSDM switching probabilities
!>                        traj%gpreprob_s3, traj%preprob_s3, preprob_old_s3 - CSDM or SE pre switching probabilitirs
!>                                                                            only used in BSH integrator
!>                        traj%decotime_diag_s, traj%decotime_diag_old_s  - CSDM decoherent time 
!>                        traj%decotime_MCH_s, traj%decotime_MCH_old_s - CSDM decoherent time 
!>                        traj%svec_diag_sad, traj%svec_MCH_sad - CSDM decoherent direction vector
!>                        traj%dmag, traj%dmag1, traj%dmag2 - CSDM record, magnitude of NAC
!>                        traj%NACdR_diag_ssad - NAC in diagonal representation
!>                        traj%trans_rot_P - projection operator, used in projected NAC
!>                        traj%pNACdR_MCH_ssad, traj%pNACdR_diag_ssad - projected NAC
!>                        traj%NACGV_ssad - CSDM TDC approximated NAC
!>                        traj%decograd_ad - CSDM decoherent force
!>                        traj%ccoeff_MCH_s, traj%ccoeff_diag_s, traj%ccoeff_diag_old_s - CSDM coherent coefficient
!>                        traj%gRcoeff_MCH_s, traj%gRccoeff_MCH_s, traj%gIcoeff_MCH_s,
!>                        traj%gIccoeff_MCH_s, traj%ephase_s, traj%gephase_s - gradients of coefficients and phase angles
!>                                                                             only used in BSH integrator
!>                        traj%consistency, traj%discrepancy - used by adaptive step size
!>                        ctrl%tmax - maximum simulation time
!>                        ctrl%dtmin - minimum allowed timestep in adapative integrators - used by adaptive step size
!>                        ctrl%convthre - convergence threshold for adaptive integrators
!>                        ctrl%method - control of method, TSH or SCP
!>                        ctrl%integrator - control of the integrator, fvv, avv, or bsh 
!>                        ctrl%nac_projection - control of using projected NAC 
!>                        ctrl%decoherence - add decay of mixing for Ehrenfest dynamics -> CSDM, SCDM, or NDM
!>                        ctrl%switching_procedure - switching procedure used by CSDM
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
  real*8,allocatable :: geom_tmp_ad(:,:)                 !< Temporary Cartesian displacement during A-FSSH step
  real*8,allocatable :: veloc_tmp_ad(:,:)                !< Temporary Cartesian displacment velocity during A-FSSH step

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
  integer :: RNGseed_thermostat                          !< seed for random number generator used in thermostat
  integer :: step                                        !< current timestep (step=0 => t=0)
  real*8 :: microtime                                    !< current trajecotry time in a.u.
  integer*8 :: traj_hash                                 !< Trajectory ID based on hashing the input files

  integer :: state_MCH                                   !< currently occupied state in MCH basis
  integer :: state_MCH_old                               !< MCH state occupied in the last timestep
  integer :: state_diag                                  !< currently occupied state in diag basis
  integer :: state_diag_old                              !< diag state occupied in the last timestep
  integer :: state_diag_frust                            !< diag state of a frustrated hop

  real*8 :: Ekin                                         !< kinetic energy
  real*8 :: Epot                                         !< potential energy (diag energy of state_diag)
  real*8 :: Etot                                         !< total energy, Ekin+Epot
  real*8 :: Ekin_old                                     !< kinetic energy of last step
  real*8 :: Epot_old                                     !< potential energy (diag energy of state_diag) of last step
  real*8 :: Etot_old                                     !< total energy, Ekin+Epot of last step

  integer :: time_start                                  !< system time when sharc is started
  integer :: time_last                                   !< system time after completion of the previous timestep
  integer :: time_step                                   !< system time(after timestep) - system time(before timestep)
  integer :: kind_of_jump                                !< 0=no jump, 1=regular, 2=frustrated, 3=resonant, 4=forced
  integer :: steps_in_gs                                 !< counter for the number of timesteps in the lowest state
  integer :: ncids(10)                                   !< NetCDF indices
  integer :: nc_index                                    !< number of steps written to NetCDF

  logical :: phases_found                                !< whether wavefunction phases were found in QM.out

  ! nuclear information
  real*8,allocatable :: atomicnumber_a(:)                !< atomic number
  character*2,allocatable :: element_a(:)                !< element descriptor
  real*8,allocatable :: mass_a(:)                        !< atomic mass in a.u. (1 a.u. = rest mass of electron m_e)
#ifdef __PYSHARC__
  real*8,pointer :: geom_ad(:,:)                         !< Cartesian coordinates of atom in a.u. (bohr)
#else
  real*8,allocatable :: geom_ad(:,:)                     !< Cartesian coordinates of atom in a.u. (bohr)
#endif
  real*8,allocatable :: veloc_ad(:,:)                    !< Cartesian velocity in a.u. (bohr/atu)
  real*8,allocatable :: veloc_old_ad(:,:)                !< Cartesian velocity in a.u. (bohr/atu) of last timestep
  real*8,allocatable :: veloc_app_ad(:,:)                !< Forward verlet approximated Cartesian velocity in a.u. (bohr/atu)

  real*8,allocatable :: accel_ad(:,:)                    !< Cartesian acceleration in a.u. (bohr/atu/atu)

  ! electronic information
#ifdef __PYSHARC__
  complex*16, pointer :: H_MCH_ss(:,:)                   !< MCH Hamiltonian as read from QM.out (no laser)
#else
  complex*16, allocatable :: H_MCH_ss(:,:)               !< MCH Hamiltonian as read from QM.out (no laser)
#endif
                                                         !< Laser interaction is added during propagation
  complex*16,allocatable :: dH_MCH_ss(:,:)               !< time derivative of MCH Hamiltonian
  complex*16,allocatable :: dH_MCH_old_ss(:,:)           !< time derivative of MCH Hamiltonian of last timestep
  complex*16,allocatable :: H_MCH_old_ss(:,:)            !< MCH Hamiltonian of last timestep (no laser)
  complex*16,allocatable :: H_MCH_old2_ss(:,:)           !< MCH Hamiltonian of second last timestep (no laser)
  complex*16,allocatable :: H_MCH_old3_ss(:,:)           !< MCH Hamiltonian of third last timestep (no laser)
  complex*16,allocatable :: H_diag_ss(:,:)               !< diag Hamiltonian
  complex*16,allocatable :: H_diag_old_ss(:,:)           !< diag Hamiltonian of last timestep (no laser)
  complex*16,allocatable :: H_diag_old2_ss(:,:)          !< diag Hamiltonian of second last timestep (no laser)
  complex*16,allocatable :: H_diag_old3_ss(:,:)          !< diag Hamiltonian of third last timestep (no laser)
  complex*16,allocatable :: U_ss(:,:)                    !< transformation matrix
  complex*16,allocatable :: U_old_ss(:,:)                !< transformation matrix of last timestep
  complex*16,allocatable :: NACdt_ss(:,:)                !< time-derivatives of wavefunctions
  complex*16,allocatable :: NACdt_old_ss(:,:)            !< time-derivatives of wavefunctions of last timestep
  complex*16,allocatable :: NACdt_diag_ss(:,:)           !< time-derivatives of wavefunctions in diagonal basis
  complex*16,allocatable :: NACdt_diag_old_ss(:,:)       !< time-derivatives of wavefunctions in diagonal basis of last timestep
#ifdef __PYSHARC__
  complex*16, pointer :: overlaps_ss(:,:)                !< overlaps for LD propagation
  complex*16, pointer :: DM_ssd(:,:,:)                   !< (transition) dipole moment matrix
                                                         !< transition dipoles between active and inactive states are zero.
#else
  complex*16, allocatable :: overlaps_ss(:,:)            !< overlaps for LD propagation
  complex*16, allocatable :: DM_ssd(:,:,:)               !< (transition) dipole moment matrix
                                                         !< transition dipoles between active and inactive states are zero.
#endif
  complex*16,allocatable :: DM_old_ssd(:,:,:)            !< old dipole moment matrix
  complex*16,allocatable :: DM_print_ssd(:,:,:)          !< dipole moment matrix used for the output routines
                                                         !< transition dipoles between active and inactive states are not zero.
  complex*16,allocatable :: Rtotal_ss(:,:)               !< total propagator for the current timestep
  complex*16,allocatable :: RDtotal_ss(:,:)              !< total propagator with decay of mixing for the current timestep
  complex*16,allocatable :: Dtotal_ss(:,:)               !< total decoherent propagator for the current timestep
  complex*16,allocatable :: phases_s(:)                  !< electronic state phases of the current step
  complex*16,allocatable :: phases_old_s(:)              !< electronic state phases of the last step
  real*8, allocatable :: hopprob_s(:)                    !< hopping probabilities
  real*8, allocatable :: switchprob_s(:)                 !< switching probabilities
  real*8, allocatable :: gpreprob_s3(:,:)                !< gradient of pre-probabilities used by BSH integrator, has dimension of nstates*3, 3 corresponds to b_ij, b_ij(+), and b_ij(-)
  real*8, allocatable :: preprob_s3(:,:)                 !< pre-probabilities used by BSH integrator, has dimension of nstates*3, 3 corresponds to b_ij, b_ij(+), and b_ij(-)
  real*8, allocatable :: preprob_old_s3(:,:)             !< pre-probabilities used by BSH integrator of the last step
  real*8, allocatable :: decotime_diag_s(:)              !< decoherence time in diag basis
  real*8, allocatable :: decotime_diag_old_s(:)          !< old decoherence time in diag basis
  real*8, allocatable :: decotime_MCH_s(:)               !< decoherence time in MCH basis
  real*8, allocatable :: decotime_MCH_old_s(:)           !< old decoherence time in MCH basis
  real*8, allocatable :: svec_diag_sad(:,:,:)            !< s vector in diag basis
  real*8, allocatable :: svec_MCH_sad(:,:,:)             !< s vector in MCH basis
  real*8, allocatable :: psvec_diag_sad(:,:,:)           !< projected s vector in diag basis
  real*8, allocatable :: psvec_MCH_sad(:,:,:)            !< projected s vector in MCH basis
  real*8 :: randnum, randnum2                            !< random number for surface hopping and A-FSSH
  real*8 :: dmag, dmag1, dmag2                           !< magnitude of NACdR for the current and previous two steps
 
  ! arbitrary properties
  complex*16,allocatable :: Property2d_xss(:,:,:)        !< list of matrices containing arbitrary data (not used in propagation)
  real*8,    allocatable :: Property1d_ys(:,:)           !< list of vectors containing arbitrary data (not used in propagation)
  character*40,allocatable :: Property2d_labels_x(:)     !< list of descriptions for 2d properties
  character*40,allocatable :: Property1d_labels_y(:)     !< list of descriptions for 1d properties

  ! vector information
  real*8,allocatable :: DMgrad_ssdad(:,:,:,:,:)          !< Cartesian gradient of the dipole moments (bra, ket, polarization, atom, cartesian component of atom displacement)
#ifdef __PYSHARC__
  real*8, pointer :: NACdR_ssad(:,:,:,:)                 !< vectorial non-adiabatic couplings in a.u.
  real*8, pointer :: grad_MCH_sad(:,:,:)                 !< Cartesian gradient in a.u (hartree/bohr) of all states
#else
  real*8, allocatable :: NACdR_ssad(:,:,:,:)             !< vectorial non-adiabatic couplings in a.u.
  real*8, allocatable :: grad_MCH_sad(:,:,:)             !< Cartesian gradient in a.u (hartree/bohr) of all states
#endif
  real*8, allocatable :: grad_MCH_old_sad(:,:,:)         !< Cartesian gradient in a.u (hartree/bohr) of all states of last timestep
  real*8, allocatable :: grad_MCH_old2_sad(:,:,:)        !< Cartesian gradient in a.u (hartree/bohr) of all states of second last timestep
  real*8,allocatable :: NACdR_old_ssad(:,:,:,:)          !< vectorial non-adiabatic couplings of last timestep
  complex*16, allocatable :: NACdR_diag_ssad(:,:,:,:)    !< vectorial non-adiabatic couplings in diagonal basis in a.u.
  real*8, allocatable :: trans_rot_P(:,:)                !< project operator that removes translational and rotational motions
  real*8, allocatable :: pNACdR_MCH_ssad(:,:,:,:)        !< projected vectorial non-adiabatic couplings in MCH in a.u.
  complex*16, allocatable :: pNACdR_diag_ssad(:,:,:,: )  !< projected vectorial non-adiabatic couplings in diag basis in a.u.
  complex*16, allocatable :: Gmatrix_ssad(:,:,:,:)       !< Cartesian gradient in a.u (hartree/bohr) of the current state in diag basis
  complex*16, allocatable :: Gmatrix_old_ssad(:,:,:,:)   !< Cartesian gradient in a.u (hartree/bohr) of the current state in diag basis of last timestep
  complex*16, allocatable :: Kmatrix_MCH_ss(:,:)         !< Time gradient in a.u (hartree/bohr) of the current state in MCH basis
  complex*16, allocatable :: Kmatrix_diag_ss(:,:)        !< Time gradient in a.u (hartree/bohr) of the current state in diag basis
  complex*16, allocatable :: NACGV_MCH_ssad(:,:,:,:)     !< effective NAC in MCH basis in a.u.
  complex*16, allocatable :: NACGV_diag_ssad(:,:,:,:)    !< effective NAC in diagonal basis in a.u.
  complex*16, allocatable :: pNACGV_MCH_ssad(:,:,:,:)    !< projected effective NAC in MCH basis in a.u.
  complex*16, allocatable :: pNACGV_diag_ssad(:,:,:,:)   !< projected effective NAC in diagonal basis in a.u.
  complex*16, allocatable :: dendt_MCH_ss(:,:)           !< time gradient of density matrix in MCH basis
  complex*16, allocatable :: dendt_diag_ss(:,:)          !< time gradient of density matrix in diagonal basis
  real*8, allocatable :: hopping_direction_ssad(:,:,:,:) !< hopping direction 
  real*8, allocatable :: frustrated_hop_vec_ssad(:,:,:,:)!< the vector to perform frustrated hop 

  real*8, allocatable :: grad_ad(:,:)                    !< final gradient used in velocity verlet
  real*8, allocatable :: decograd_ad(:,:)                !< decoherent gradient

  ! coefficient information
  complex*16, allocatable :: coeff_diag_s(:)             !< coefficients of electronic wavefunction in diag representation
  complex*16, allocatable :: coeff_diag_old_s(:)         !< coefficients of electronic wavefunction (of previous timestep) in diag representation
  complex*16, allocatable :: coeff_MCH_s(:)              !< coefficients of electronic wavefunction in MCH representation
  complex*16, allocatable :: coeff_MCH_old_s(:)          !< coefficients of electronic wavefunction (of previous timestep) in MCH representation

  ! coherent coefficient information
  complex*16, allocatable :: ccoeff_MCH_s(:)             !< coherent coefficients of electronic wavefunction in MCH representation, used in CSDM 
  complex*16, allocatable :: ccoeff_diag_s(:)            !< coherent coefficients of electronic wavefunction in diag representation, used in CSDM
  complex*16, allocatable :: ccoeff_diag_old_s(:)        !< coherent coefficients of electronic wavefunction (of previous timestep) in diag representation
 
  ! coeficient gradient-used by Bulirsch_Stoer_Hack integrator
  real*8, allocatable :: gRcoeff_MCH_s(:)                !< Real part coefficients gradient of electronic wavefunction in MCH representation
  real*8, allocatable :: gRccoeff_MCH_s(:)               !< Real part coherent coefficients gradiet of electronic wavefunction in diag representation, used in CSDM
  real*8, allocatable :: gIcoeff_MCH_s(:)                !< Imaginary part coefficients gradient of electronic wavefunction in MCH representation
  real*8, allocatable :: gIccoeff_MCH_s(:)               !< Imaginary part coherent coefficients gradiet of electronic wavefunction in diag representation, used in CSDM
  real*8, allocatable :: ephase_s(:)                     !< phase angle of electronic coordinates 
  real*8, allocatable :: gephase_s(:)                    !< phase angle derivative of electronic coordinates

  ! Hamiltonian, energy and gradients used in ZPE pumping
  complex*16, allocatable :: Evib_local_s(:)             !< local mode vibrational energy, used in ZPE pumping
  complex*16, allocatable :: Ezpe_local_s(:)             !< local mode zero point energy, used in ZPE pumping
  real*8, allocatable :: grad_Ezpe_local_sad(:,:,:)      !< Cartesian gradient in a.u (hartree/bohr) of all states of Ezpe
  complex*16, allocatable :: coeff_zpe_s2(:,:)           !< coeffcient of ZPE pumping states
  integer, allocatable :: state_pumping_s(:)             !< the pointer state in pumping Hamiltonian, one for each state, either 1 or 2
  integer, allocatable :: pumping_status_s(:)            !< whether it's in a pumping process or not 0=no, 1=yes

  ! LP-ZPE correction scheme
  real*8, allocatable :: lpzpe_ke_ah(:)                  !< kinetic energy of A-H mode
  real*8, allocatable :: lpzpe_ke_bc(:)                  !< kinetic energy of B-C mode
  integer :: lpzpe_cycle                                 !< cycles of vibrations
  integer :: lpzpe_iter_incycle                          !< iterations within a cycle
  real*8 :: lpzpe_starttime                              !< the start time of cummulation
  real*8 :: lpzpe_endtime                                !< the end time of cummulation
  integer :: in_cycle                                    !< indicator of whether in cycles, 0=out of cycles, 1=in cycles 
 
  ! gradient and nac selection information
  logical,allocatable :: selG_s(:)                       !< selection mask for gradients
  logical,allocatable :: selT_ss(:,:)                    !< selection mask for non-adiabatic coupling vectors
  logical,allocatable :: selDM_ss(:,:)                   !< selection mask for dipole moment gradients

  ! Auxiliary trajectories for A-FSSH
  type(aux_trajectory_type),allocatable :: auxtrajs_s(:)

  ! TSH with time uncertainty
  real*8, allocatable :: uncertainty_time_s(:)           !< maximum allowed travelling time in tsh time uncertainty
  real*8 :: uncertainty_time_frust                       !< maximum allowed travelling time for the frustrated state at position where frustrated hop happened
  real*8 :: incident_time                                !< the time when frustrated hop happens
  integer, allocatable :: allowed_hop_s(:)               !< whether a state is an energetically allowed hop or not, 0=frustrated, 1=allowed
  real*8, allocatable :: allowed_time_s(:)               !< the nearest time when a hop is allowed
  integer :: in_time_uncertainty                         !< an indicator of whether the trajectory is in time uncertainty process or not, 0=not, 1=in, 2=right after time travelling
  integer :: tu_backpropagation                          !< an indicator of backpropagation in time uncertainty
  integer :: travelling_state                            !< the target state when trajectory travels back in time uncertainty process
  integer :: target_state                                !< recorded frustrated state

  ! Army ants weight
  real*8 :: army_ants_weight                             !< trajectory weight for army ants rare event sampling
  real*8 :: randnum_branching                            !< the random number to decide a branching event
  real*8 :: branching_likehood                           !< the gamma parameter, used to decide the branching probability

  ! Trajectory entropy
  complex*16 :: entropy                                  !< trajectory entropy 
  complex*16 :: entropy_old                              !< trajectory entropy of last timestep 
  complex*16 :: linear_entropy                           !< trajectory linearized entropy 
  complex*16 :: linear_entropy_old                       !< trajectory linearized entropy of last timestep 
  complex*16 :: entropy_grad                             !< derivative of trajectory entropy w.r.t. time
  complex*16 :: linear_entropy_grad                      !< derivative of trajectory inearized entropy w.r.t. time
  complex*16,allocatable :: U_pointer_ss(:,:)            !< transformation matrix for pointer state optimization 
  complex*16,allocatable :: U_pointer_old_ss(:,:)        !< transformation matrix for pointer state optimization of last time step 

  ! Trajectory consistency
  integer :: consistency                                 !< consistency between current and previous step, 0=consistent, 1=increase time step, 2=decrease time step
  real*8 :: discrepancy                                  !< the difference between current and previous step

  ! Thermostat randomness
  real*8,allocatable :: thermostat_random(:)
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
  integer :: compat_mode                    !< compatibility mode with older versions
                                            ! currently: 
                                            ! 0 : no compatibility mode
                                            ! 1 : in EDC, do not draw a second random number

! numerical constants
  integer :: natom                          !< number of atoms
  logical,allocatable :: atommask_a(:)      !< atoms which are considered for decoherence, rescaling, ...
  integer :: maxmult                        !< highest spin quantum number (determines length of nstates_m)
  integer,allocatable :: nstates_m(:)       !< numer of states considered in each multiplicy
  integer :: nstates                        !< total number of states
  integer :: nsteps                         !< total number of simulation steps
  integer :: nsubsteps                      !< number of steps for the electron propagation
  real*8 :: tmax                            !< maximum time step in a.u (atu)
  real*8 :: dtstep_min                      !< minimum time step allowed in adapative algorithm (atu)
  real*8 :: dtstep_max                      !< maximum time step allowed in adapative algorithm (atu)
  real*8 :: dtstep                          !< length of timestep in a.u (atu)
  real*8 :: dtstep_old                      !< length of timestep in a.u (atu) of last timestep
  real*8 :: ezero                           !< energy offset in a.u. (e.g. ground state equilibrium energy)
  real*8 :: convthre                        !< convergence threshold for Bulisrch-Stoer integrator
  real*8 :: scalingfactor                   !< scales the Hamiltonian and gradients
  real*8 :: soc_scaling                     !< scales spin orbit coupling only
  real*8 :: eselect_grad                    !< energy difference for neglecting gradients
  real*8 :: eselect_nac                     !< energy difference for neglecting na-couplings
  real*8 :: eselect_dmgrad                  !< energy difference for neglecting dipole gradients
  real*8 :: dampeddyn                       !< damping factor for kinetic energy
  real*8 :: decoherence_alpha               !< decoherence parameter (a.u.) for energy-based decoherence
  real*8 :: decoherence_beta                !< decoherence parameter (a.u.) for energy-based decoherence
  real*8 :: gaussian_width                  !< Gaussian width for force-momentum2 method to compute decoherence time
  real*8 :: force_hop_to_gs                 !< if positive, trajectories automatically jump to lowest state if active-lowest gap is below value
  logical,allocatable :: actstates_s(:)     !< mask of the active states
  integer :: output_steps_stride(3)         !< how often output.dat is written
  integer :: output_steps_limits(3)         !< switches stride for output.dat writing

! methods and switches
  logical :: restart                        !< restart yes or no
  logical :: restart_rerun_last_qm_step     !< if true, then qm.f90 will write "restart" instruction
  integer :: method                         !< 0=trajectory surface hopping(tsh), 1=self-consistent potential(scp)
  integer :: integrator                     !< integrator used, 0=Bulirsch-Stoer, 1=adaptive Velocity Verlet, 2=fixed stepzie Velocity Verlet
  integer :: staterep                       !< 0=initial state is given in diag representation, 1=in MCH representation
  integer :: initcoeff                      !< 0=initial coefficients are diag, 1=initial coefficients are MCH, 2=auto diag, 3=auto MCH
  integer :: laser                          !< 0=none, 1=internal, 2=external
  integer :: coupling                       !< 0=ddt, 1=ddr, 2=overlap, 3=ktdc
  integer :: ktdc_method                    !< 0=gradient based approximation, 1=energy based approximation
  integer :: kmatrix_method                 !< 0=gradient based approximation, 1=energy based approximation
  integer :: eeom                           !< method to control the electronic eom for both tsh and scp, 0=constant interpolation, 1=linear interpolation, 2=local diabatization, 3=norm perserving interpolation
  integer :: neom                           !< method to control the nuclear eom for scp, 0=ddr, 1=gradient difference
  integer :: neom_rep                       !< representation to compute SCP force
  integer :: surf                           !< 0=propagation in diag surfaces (SHARC), 1=on MCH surfaces (regular SH)
  integer :: decoherence                    !< 0=off, 1=EDC, 2=AFSSH, 11=DoM
  integer :: decotime_method                !< method to control computation of decoherence time, 0=originalCSDM, 1=masslessCSDM, 2=edc, 3=scdm, 4=FP, 5=FPA
  integer :: ekincorrect                    !< 0=none, 1=adjust momentum along velocity, 2=adjust momentum along projected velocity,
                                            !< 3=adjust momentum along nac vector, 4=adjust momentum along gradient difference, 
                                            !< 5=adjust momentum along projected nac vector, 6=adjust momentum along projected gradient difference,
                                            !< 7=adjust momentum along effective nac vector, 8=adjust momentum along projected effective nac vector
  integer :: reflect_frustrated             !< 0=none, 1=reflect along velocity, 2=reflect along projected velocity, 
                                            !< 3=reflect along nac vector, 4=reflect momentum along gradient difference, 
                                            !< 5=reflect momentum along projected nac vector, 6=reflect momentum along projected gradient difference,
                                            !< 7=reflect momentum along effective nac vector, 8=reflect momentum along projected effective nac vector
  integer :: time_uncertainty               !< 0=no time uncertainty, 1=do TSH with time uncertainty, a method to improve frustrated hops
  integer :: gradcorrect                    !< 0=no, 1=include nac vectors in gradient transformation, 2=using Kmatrix to correct
  integer :: dipolegrad                     !< 0=no, 1=include dipole gradients in gradient transformation
  integer :: nac_projection                 !< 0=no, 1=do nac projection, remove translational and rotational motions.
  integer :: zpe_correction                 !< method to control ZPE correction in nonadiabatic dynamics, 0=no correction, 1=pumping, 2=lp
  integer :: pointer_basis                  !< compute pointer basis, 0=adiabatic/diagonal basis, 1=mch basis, 2=optimized pointer basis
  integer :: pointer_maxiter                !< maximum iterations allowed in pointer state optimization
  
  integer :: thermostat                     !< 0=none, 1=Langevin thermostat
  logical :: restart_thermostat_random      !< F=no, T=yes (default) to use same random number sequence if restarted

! lp-zpe
  integer :: lpzpe_scheme                   !< correction_scheme=0 skip correction if BC bonds do not have enough kinetic energy; 1=forced correction by moving all kinetic energies of BC bonds to AH bonds
  integer :: lpzpe_nah                      !< number of A-H bond pairs
  integer :: lpzpe_nbc                      !< number of B-C bond pairs
  integer, allocatable :: lpzpe_ah(:,:)     !< atom indicies of A-H bond pairs
  integer, allocatable :: lpzpe_bc(:,:)     !< atom indicies of B-C bond pairs
  real*8, allocatable :: lpzpe_ke_zpe_ah(:) !< kinetic energy of A-H mode
  real*8, allocatable :: lpzpe_ke_zpe_bc(:) !< kinetic energy of B-C mode
  real*8 :: ke_threshold                    !< threshold for zpe leaking
  real*8 :: t_cycle                         !< time cycle of zpe checking
  real*8 :: t_check                         !< period of time for kinetic energy averaging
 
  integer :: calc_soc                       !< request SOC, otherwise only the diagonal elements of H (plus any laser interactions) are taken into account\n 0=no soc, 1=soc enabled
  integer :: calc_grad                      !< request gradients:   \n        0=all in step 1, 1=select in step 1, 2=select in step 2
  integer :: calc_overlap                   !< 0=no, 1=request overlap matrices
  integer :: calc_nacdt                     !< 0=no, 1=request time derivatives
  integer :: calc_nacdr                     !< request nac vectors: \n -1=no, 0=all in step 1, 1=select in step 1, 2=select in step 2
  integer :: calc_effectivenac              !< 0=no, 1=request effective nac vectors
  integer :: calc_dipole                    !< 0=no, 1=request dipolemoment 
  integer :: calc_dipolegrad                !< request dipole gradient vectors: \n -1=no, 0=all in step 1, 1=select in step 1, 2=select in step 2
  integer :: calc_second                    !< 0=no, 1=do two interface calls per timestep
  integer :: calc_phases                    !< 0=no, 1=yes

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
  integer :: switching_procedure            !< 0=no switches, 1=CSDM, 2=SCDM, 3=NDM
  integer :: army_ants                      !< 0=no army ants,i.e. anteater algorithm, 1=army ants algorithm
  integer :: output_format                  !< 0 ASCII, 1 NetCDF

! thresholds
!   real*8 :: propag_sharc_UdUdiags=1.d-2           ! Threshold for the size of diagonal elements in UdU (needed for dynamic substeps)        in hartree
!   real*8 :: min_dynamic_substep=1.d-5             ! In dynamic substepping, the shortest substep allowed                                        in atomic time units
!   real*8 :: diagonalize_degeneracy_diff=1.d-9     ! Energy difference threshold for treating states as degenerate                                in hartree

  ! RATTLE
  integer :: do_constraints                    !< 0=none, 1=rattle
  real*8 :: constraints_tol                    !< tolerance for RATTLE
  integer :: n_constraints                     !< number of constraints
  integer, allocatable :: constraints_ca(:,:)  !< atom pairs, first index is constraint, second is atom (1 or 2)
  real*8, allocatable :: constraints_dist_c(:) !< squared distance for each constraint pair

  ! laser
  real*8 :: laser_bandwidth                       !< for detecting induced hops (in a.u.)
  integer :: nlasers
  complex*16, allocatable :: laserfield_td(:,:)   !< complex valued laser field
  complex*16, allocatable :: laserenergy_tl(:,:)  !< momentary central energy of laser (for detecting induced hops)

  ! thermostat
  real*8 :: temperature                     !< temperature used for thermostat
  real*8,allocatable :: thermostat_const(:) !< constants needed for thermostat. Langevin: friction coeffitient

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

character*255, parameter :: version='3.0 (January 26, 2023)'    !< string holding the version number

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
integer, parameter :: u_i_rattle=18          !< atoms for constraints

integer, parameter :: u_qm_QMin=41           !< here SHARC writes information for the QM interface (like geometry, number of states, what kind of data is requested)
integer, parameter :: u_qm_QMout=42          !< here SHARC retrieves the results of the QM run (Hamiltonian, gradients, couplings, etc.)

! =========================================================== !

  contains

! =========================================================== !

    subroutine allocate_traj(traj,ctrl)
      !< Allocates almost all arrays in traj
      !< Memory-heavy allocations are done in additional_allocate_traj
      !< Does not allocate arrays in ctrl (laser-related, actstates_s,
      !  nstates_m, lpzpe_ah, lpzpe_bc, lpzpe_ke_zpe_ah, lpzpe_ke_zpe_bc)
      !< Reads natom and nstates from ctrl
      !< Initializes all elements of all arrays to -123, 'Q' or .true.
      implicit none
      type(ctrl_type), intent(inout) :: ctrl
      type(trajectory_type), intent(inout) :: traj
      integer :: status
      integer :: natom,nstates,n1d,n2d,nah,nbc

      natom=ctrl%natom
      nstates=ctrl%nstates
      n1d=ctrl%n_property1d
      n2d=ctrl%n_property2d

      nah=ctrl%lpzpe_nah
      nbc=ctrl%lpzpe_nbc

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

      allocate(traj%veloc_old_ad(natom,3),stat=status)
      if (status/=0) stop 'Could not allocate veloc_old_ad'
      traj%veloc_old_ad=-123.d0

      allocate(traj%veloc_app_ad(natom,3),stat=status)
      if (status/=0) stop 'Could not allocate veloc_app_ad'
      traj%veloc_app_ad=-123.d0

      allocate(traj%accel_ad(natom,3),stat=status)
      if (status/=0) stop 'Could not allocate accel_ad'
      traj%accel_ad=-123.d0


      allocate(traj%H_MCH_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate H_MCH_ss'
      traj%H_MCH_ss=-123.d0

      allocate(traj%dH_MCH_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate dH_MCH_ss'
      traj%dH_MCH_ss=-123.d0

      allocate(traj%dH_MCH_old_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate dH_MCH_old_ss'
      traj%dH_MCH_old_ss=-123.d0

      allocate(traj%H_MCH_old_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate H_MCH_old_ss'
      traj%H_MCH_old_ss=-123.d0

      allocate(traj%H_MCH_old2_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate H_MCH_old2_ss'
      traj%H_MCH_old2_ss=-123.d0

      allocate(traj%H_MCH_old3_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate H_MCH_old3_ss'
      traj%H_MCH_old3_ss=-123.d0

      allocate(traj%H_diag_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate H_diag_ss'
      traj%H_diag_ss=-123.d0

      allocate(traj%H_diag_old_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate H_diag_old_ss'
      traj%H_diag_old_ss=-123.d0

      allocate(traj%H_diag_old2_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate H_diag_old2_ss'
      traj%H_diag_old2_ss=-123.d0

      allocate(traj%H_diag_old3_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate H_diag_old3_ss'
      traj%H_diag_old3_ss=-123.d0

      allocate(traj%Evib_local_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate Evib_local_s'
      traj%Evib_local_s=-123.d0

      allocate(traj%Ezpe_local_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate Ezpe_local_s'
      traj%Ezpe_local_s=-123.d0

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

      allocate(traj%NACdt_diag_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate NACdt_diag_ss'
      traj%NACdt_diag_ss=-123.d0

      allocate(traj%NACdt_diag_old_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate NACdt_diag_old_ss'
      traj%NACdt_diag_old_ss=-123.d0

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

      allocate(traj%RDtotal_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate RDtotal_ss'
      traj%RDtotal_ss=-123.d0

      allocate(traj%Dtotal_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate Dtotal_ss'
      traj%Dtotal_ss=-123.d0

      allocate(traj%hopprob_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate hopprob_s'
      traj%hopprob_s=-123.d0

      allocate(traj%switchprob_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate switchprob_s'
      traj%switchprob_s=-123.d0

      allocate(traj%gpreprob_s3(nstates,3),stat=status)
      if (status/=0) stop 'Could not allocate gpreprob_s3'
      traj%gpreprob_s3=-123.d0

      allocate(traj%preprob_s3(nstates,3),stat=status)
      if (status/=0) stop 'Could not allocate preprob_s3'
      traj%preprob_s3=-123.d0

      allocate(traj%preprob_old_s3(nstates,3),stat=status)
      if (status/=0) stop 'Could not allocate preprob_old_s3'
      traj%preprob_old_s3=-123.d0

      allocate(traj%decotime_diag_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate decotime_diag_s'
      traj%decotime_diag_s=-123.d0

      allocate(traj%decotime_diag_old_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate decotime_diag_old_s'
      traj%decotime_diag_old_s=-123.d0

      allocate(traj%decotime_MCH_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate decotime_MCH_s'
      traj%decotime_MCH_s=-123.d0

      allocate(traj%decotime_MCH_old_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate decotime_MCH_old_s'
      traj%decotime_MCH_old_s=-123.d0

      allocate(traj%svec_diag_sad(nstates,natom,3),stat=status)
      if (status/=0) stop 'Could not allocate svec_diag_sad'
      traj%svec_diag_sad=-123.d0

      allocate(traj%svec_MCH_sad(nstates,natom,3),stat=status)
      if (status/=0) stop 'Could not allocate svec_MCH_sad'
      traj%svec_MCH_sad=-123.d0

      allocate(traj%psvec_diag_sad(nstates,natom,3),stat=status)
      if (status/=0) stop 'Could not allocate psvec_diag_sad'
      traj%psvec_diag_sad=-123.d0

      allocate(traj%psvec_MCH_sad(nstates,natom,3),stat=status)
      if (status/=0) stop 'Could not allocate psvec_MCH_sad'
      traj%psvec_MCH_sad=-123.d0

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

      allocate(traj%coeff_MCH_old_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate coeff_MCH_old_s'
      traj%coeff_MCH_old_s=-123.d0

      allocate(traj%ccoeff_MCH_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate ccoeff_MCH_s'
      traj%ccoeff_MCH_s=-123.d0

      allocate(traj%ccoeff_diag_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate ccoeff_diag_s'
      traj%ccoeff_diag_s=-123.d0

      allocate(traj%ccoeff_diag_old_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate ccoeff_diag_old_s'
      traj%ccoeff_diag_old_s=-123.d0

      allocate(traj%coeff_zpe_s2(nstates,2),stat=status)
      if (status/=0) stop 'Could not allocate coeff_zpe_s2'
      traj%coeff_zpe_s2=-123.d0

      allocate(traj%gRcoeff_MCH_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate gRcoeff_MCH_s'
      traj%gRcoeff_MCH_s=-123.d0

      allocate(traj%gRccoeff_MCH_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate gRccoeff_MCH_s'
      traj%gRccoeff_MCH_s=-123.d0

      allocate(traj%gIcoeff_MCH_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate gIcoeff_MCH_s'
      traj%gIcoeff_MCH_s=-123.d0

      allocate(traj%gIccoeff_MCH_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate gIccoeff_MCH_s'
      traj%gIccoeff_MCH_s=-123.d0

      allocate(traj%ephase_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate ephase_s'
      traj%ephase_s=-123.d0

      allocate(traj%gephase_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate gephase_s'
      traj%gephase_s=-123.d0

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

      allocate(traj%NACdR_diag_ssad(nstates,nstates,natom,3),stat=status)
      if (status/=0) stop 'Could not allocate NACdR_diag_ssad'
      traj%NACdR_diag_ssad=-123.d0

      allocate(traj%dendt_MCH_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate dendt_MCH_ss'
      traj%dendt_MCH_ss=-123.d0

      allocate(traj%dendt_diag_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate dendt_diag_ss'
      traj%dendt_diag_ss=-123.d0

      allocate(traj%grad_MCH_sad(nstates,natom,3),stat=status)
      if (status/=0) stop 'Could not allocate gradMCH_sad'
      traj%grad_MCH_sad=-123.d0

      allocate(traj%grad_MCH_old_sad(nstates,natom,3),stat=status)
      if (status/=0) stop 'Could not allocate gradMCH_old_sad'
      traj%grad_MCH_old_sad=-123.d0

      allocate(traj%grad_MCH_old2_sad(nstates,natom,3),stat=status)
      if (status/=0) stop 'Could not allocate gradMCH_old2_sad'
      traj%grad_MCH_old2_sad=-123.d0

      allocate(traj%grad_Ezpe_local_sad(nstates,natom,3),stat=status)
      if (status/=0) stop 'Could not allocate grad_Ezpe_local_sad'
      traj%grad_Ezpe_local_sad=-123.d0

      allocate(traj%grad_ad(natom,3),stat=status)
      if (status/=0) stop 'Could not allocate grad_ad'
      traj%grad_ad=-123.d0

      allocate(traj%decograd_ad(natom,3),stat=status)
      if (status/=0) stop 'Could not allocate decograd_ad'
      traj%decograd_ad=-123.d0

      allocate(traj%hopping_direction_ssad(nstates,nstates,natom,3),stat=status)
      if (status/=0) stop 'Could not allocate hopping_direction_ssad'
      traj%hopping_direction_ssad=-123.d0

      allocate(traj%frustrated_hop_vec_ssad(nstates,nstates,natom,3),stat=status)
      if (status/=0) stop 'Could not allocate frustrated_hop_vec_ssad'
      traj%frustrated_hop_vec_ssad=-123.d0

      allocate(traj%Gmatrix_ssad(nstates,nstates,natom,3),stat=status)
      if (status/=0) stop 'Could not allocate Gmatrix_ssad'
      traj%Gmatrix_ssad=-123.d0

      allocate(traj%Gmatrix_old_ssad(nstates,nstates,natom,3),stat=status)
      if (status/=0) stop 'Could not allocate Gmatrix_ssad'
      traj%Gmatrix_ssad=-123.d0

      allocate(traj%Kmatrix_MCH_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate Kmatrix_MCH_ss'
      traj%Kmatrix_MCH_ss=-123.d0

      allocate(traj%Kmatrix_diag_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate Kmatrix_diag_ss'
      traj%Kmatrix_diag_ss=-123.d0

      allocate(traj%state_pumping_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate state_pumping_s'
      traj%state_pumping_s=0

      allocate(traj%pumping_status_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate pumping_status_s'
      traj%pumping_status_s=0

      allocate(traj%uncertainty_time_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate uncertainty_time_s'
      traj%uncertainty_time_s=-123.d0

      allocate(traj%allowed_hop_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate allowed_hop_s'
      traj%allowed_hop_s=0

      allocate(traj%allowed_time_s(nstates),stat=status)
      if (status/=0) stop 'Could not allocate allowed_time_s'
      traj%allowed_time_s=-123.d0

      allocate(traj%U_pointer_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate U_pointer_ss'
      traj%U_pointer_ss=-123.d0

      allocate(traj%U_pointer_old_ss(nstates,nstates),stat=status)
      if (status/=0) stop 'Could not allocate U_pointer_old_ss'
      traj%U_pointer_old_ss=-123.d0

      allocate(traj%lpzpe_ke_ah(nah),stat=status)
      if (status/=0) stop 'Could not allocate lpzpe_ke_ah'
      traj%lpzpe_ke_ah=-123.d0

      allocate(traj%lpzpe_ke_bc(nbc),stat=status)
      if (status/=0) stop 'Could not allocate lpzpe_ke_bc'
      traj%lpzpe_ke_bc=-123.d0

    endsubroutine


! ===========================================================
! used to manage memory-heavy allocations
! perform conditional allocation 
    subroutine additional_allocate_traj(traj,ctrl)
      implicit none
      type(ctrl_type), intent(inout) :: ctrl
      type(trajectory_type), intent(inout) :: traj
      integer :: status
      integer :: natom,nstates

      ! GB denotes a general basis, meaning it can be either MCH or diag
      integer :: allocate_trans_rot_P
      integer :: allocate_pNACdR_GB_ssad
      integer :: allocate_NACGV_GB_ssad
      integer :: allocate_pNACGV_GB_ssad

      natom=ctrl%natom
      nstates=ctrl%nstates

      if (printlevel>0) then
        write(u_log,*) '============================================================='
        write(u_log,*) '                  Memory-heavy Allocations'
        write(u_log,*) '============================================================='
      endif

      ! projection operator trans_rot_P is only computed when:
      ! 1. in SCP, use projected nonadiabatic force direction; 
      ! 2. in TSH, use projected hopping direction or velocity reflection vector. 
      allocate_trans_rot_P=0
      if (ctrl%method==1 .and. ctrl%nac_projection==1) then
        allocate_trans_rot_P=1
      else if (ctrl%method==0 .and. &
        &(ctrl%ekincorrect==2 .or. ctrl%ekincorrect==5 .or. ctrl%ekincorrect==6 .or. ctrl%ekincorrect==8)) then 
        allocate_trans_rot_P=1
      else if (ctrl%method==0 .and. &
        &(ctrl%reflect_frustrated==2 .or. ctrl%reflect_frustrated==5 .or. ctrl%reflect_frustrated==6 .or. ctrl%reflect_frustrated==8 .or. &
        &ctrl%reflect_frustrated==92 .or. ctrl%reflect_frustrated==95 .or. ctrl%reflect_frustrated==96 .or. ctrl%reflect_frustrated==98)) then 
        allocate_trans_rot_P=1
      endif 

      if (allocate_trans_rot_P==1) then
        write(u_log,*) "allocating trans_rot_P"
        allocate(traj%trans_rot_P(3*natom,3*natom),stat=status)
        if (status/=0) stop 'Could not allocate trans_rot_P'
        traj%trans_rot_P=-123.d0
      endif

      ! projected NAC is only used in SCP for nonadiabatic force direction 
      ! in TSH, if one is using projected NAC, the vector is given in
      ! traj%hopping_direction_ssad and traj%frustrated_hop_vec_ssad
      allocate_pNACdR_GB_ssad=0
      if (ctrl%method==1 .and. ctrl%nac_projection==1) then
        allocate_pNACdR_GB_ssad=1
      endif

      if (allocate_pNACdR_GB_ssad==1) then 
        write(u_log,*) "allocating projected NAC"
        allocate(traj%pNACdR_MCH_ssad(nstates,nstates,natom,3),stat=status)
        if (status/=0) stop 'Could not allocate pNACdR_MCH_ssad'
        traj%pNACdR_MCH_ssad=-123.d0
        allocate(traj%pNACdR_diag_ssad(nstates,nstates,natom,3),stat=status)
        if (status/=0) stop 'Could not allocate pNACdR_diag_ssad'
        traj%pNACdR_diag_ssad=-123.d0
      endif

      ! effective NAC is only used when calc_effectivenac is set to 1
      allocate_NACGV_GB_ssad=0
      if (ctrl%calc_effectivenac==1) then
        allocate_NACGV_GB_ssad=1
      endif
      if (allocate_NACGV_GB_ssad==1) then
        write(u_log,*) "allocating effective NAC"
        allocate(traj%NACGV_MCH_ssad(nstates,nstates,natom,3),stat=status)
        if (status/=0) stop 'Could not allocate NACGV_MCH_ssad'
        traj%NACGV_MCH_ssad=-123.d0
        allocate(traj%NACGV_diag_ssad(nstates,nstates,natom,3),stat=status)
        if (status/=0) stop 'Could not allocate NACGV_diag_ssad'
        traj%NACGV_diag_ssad=-123.d0
      endif
        
      ! projected effective NAC is only used when calc_effectivenac is set to 1 
      allocate_pNACGV_GB_ssad=0
      if (ctrl%calc_effectivenac==1 .and. ctrl%nac_projection==1) then
        allocate_pNACGV_GB_ssad=1
      endif
      if (allocate_pNACGV_GB_ssad==1) then 
        write(u_log,*) "allocating projected effective NAC"
        allocate(traj%pNACGV_MCH_ssad(nstates,nstates,natom,3),stat=status)
        if (status/=0) stop 'Could not allocate pNACGV_MCH_ssad'
        traj%pNACGV_MCH_ssad=-123.d0
        allocate(traj%pNACGV_diag_ssad(nstates,nstates,natom,3),stat=status)
        if (status/=0) stop 'Could not allocate pNACGV_diag_ssad'
        traj%pNACGV_diag_ssad=-123.d0
      endif 

    endsubroutine


    subroutine deallocate_ctrl(ctrl)
      implicit none
      type(ctrl_type), intent(inout) :: ctrl

      if (allocated(ctrl%nstates_m))                  deallocate(ctrl%nstates_m)
      if (allocated(ctrl%actstates_s))                deallocate(ctrl%actstates_s)
      if (allocated(ctrl%atommask_a))                 deallocate(ctrl%atommask_a)
      if (allocated(ctrl%laserfield_td))              deallocate(ctrl%laserfield_td)
      if (allocated(ctrl%laserenergy_tl))             deallocate(ctrl%laserenergy_tl)
      if (allocated(ctrl%lpzpe_ah))                   deallocate(ctrl%lpzpe_ah)
      if (allocated(ctrl%lpzpe_bc))                   deallocate(ctrl%lpzpe_bc)
      if (allocated(ctrl%lpzpe_ke_zpe_ah))            deallocate(ctrl%lpzpe_ke_zpe_ah)
      if (allocated(ctrl%lpzpe_ke_zpe_bc))            deallocate(ctrl%lpzpe_ke_zpe_bc)
      if (allocated(ctrl%thermostat_const))           deallocate(ctrl%thermostat_const)

    endsubroutine

    subroutine deallocate_traj(traj)

      implicit none
      type(trajectory_type), intent(inout) :: traj

#ifdef __PYSHARC__
    ! Pointer routines
    if (associated(traj%H_MCH_ss))                  deallocate(traj%H_MCH_ss)
    if (associated(traj%DM_ssd))                    deallocate(traj%DM_ssd)
    if (associated(traj%overlaps_ss))               deallocate(traj%overlaps_ss)
    if (associated(traj%grad_MCH_sad))              deallocate(traj%grad_MCH_sad)
    if (associated(traj%NACdR_ssad))                deallocate(traj%NACdR_ssad)
    if (associated(traj%geom_ad))                   deallocate(traj%geom_ad)
#else
    if (allocated(traj%H_MCH_ss))                   deallocate(traj%H_MCH_ss)
    if (allocated(traj%DM_ssd))                     deallocate(traj%DM_ssd)
    if (allocated(traj%overlaps_ss))                deallocate(traj%overlaps_ss)
    if (allocated(traj%grad_MCH_sad))               deallocate(traj%grad_MCH_sad)
    if (allocated(traj%NACdR_ssad))                 deallocate(traj%NACdR_ssad)
    if (allocated(traj%geom_ad))                    deallocate(traj%geom_ad)
#endif

    if (allocated(traj%atomicnumber_a))             deallocate(traj%atomicnumber_a)
    if (allocated(traj%element_a))                  deallocate(traj%element_a)
    if (allocated(traj%mass_a))                     deallocate(traj%mass_a)
    if (allocated(traj%veloc_ad))                   deallocate(traj%veloc_ad)
    if (allocated(traj%veloc_old_ad))               deallocate(traj%veloc_old_ad)
    if (allocated(traj%veloc_app_ad))               deallocate(traj%veloc_app_ad)
    if (allocated(traj%accel_ad))                   deallocate(traj%accel_ad)
    if (allocated(traj%dH_MCH_ss))                  deallocate(traj%dH_MCH_ss)
    if (allocated(traj%dH_MCH_old_ss))              deallocate(traj%dH_MCH_old_ss)
    if (allocated(traj%H_MCH_old_ss))               deallocate(traj%H_MCH_old_ss)
    if (allocated(traj%H_MCH_old2_ss))              deallocate(traj%H_MCH_old2_ss)
    if (allocated(traj%H_MCH_old3_ss))              deallocate(traj%H_MCH_old3_ss)
    if (allocated(traj%H_diag_ss))                  deallocate(traj%H_diag_ss)
    if (allocated(traj%H_diag_old_ss))              deallocate(traj%H_diag_old_ss)
    if (allocated(traj%H_diag_old2_ss))             deallocate(traj%H_diag_old2_ss)
    if (allocated(traj%H_diag_old3_ss))             deallocate(traj%H_diag_old3_ss)
    if (allocated(traj%Evib_local_s))               deallocate(traj%Evib_local_s)
    if (allocated(traj%Ezpe_local_s))               deallocate(traj%Ezpe_local_s)
    if (allocated(traj%U_ss))                       deallocate(traj%U_ss)
    if (allocated(traj%U_old_ss))                   deallocate(traj%U_old_ss)
    if (allocated(traj%NACdt_ss))                   deallocate(traj%NACdt_ss)
    if (allocated(traj%NACdt_old_ss))               deallocate(traj%NACdt_old_ss)
    if (allocated(traj%NACdt_diag_ss))              deallocate(traj%NACdt_diag_ss)
    if (allocated(traj%NACdt_diag_old_ss))          deallocate(traj%NACdt_diag_old_ss)
    if (allocated(traj%grad_MCH_old_sad))           deallocate(traj%grad_MCH_old_sad)
    if (allocated(traj%grad_MCH_old2_sad))          deallocate(traj%grad_MCH_old2_sad)
    if (allocated(traj%grad_Ezpe_local_sad))        deallocate(traj%grad_Ezpe_local_sad)
    if (allocated(traj%NACdR_old_ssad))             deallocate(traj%NACdR_old_ssad)
    if (allocated(traj%NACdR_diag_ssad))            deallocate(traj%NACdR_diag_ssad)
    if (allocated(traj%trans_rot_P))                deallocate(traj%trans_rot_P)
    if (allocated(traj%pNACdR_MCH_ssad))            deallocate(traj%pNACdR_MCH_ssad)
    if (allocated(traj%pNACdR_diag_ssad))           deallocate(traj%pNACdR_diag_ssad)
    if (allocated(traj%DM_old_ssd))                 deallocate(traj%DM_old_ssd)
    if (allocated(traj%DM_print_ssd))               deallocate(traj%DM_print_ssd)
    if (allocated(traj%Rtotal_ss))                  deallocate(traj%Rtotal_ss)
    if (allocated(traj%RDtotal_ss))                 deallocate(traj%RDtotal_ss)
    if (allocated(traj%Dtotal_ss))                  deallocate(traj%Dtotal_ss)
    if (allocated(traj%Kmatrix_MCH_ss))             deallocate(traj%Kmatrix_MCH_ss)
    if (allocated(traj%Kmatrix_diag_ss))            deallocate(traj%Kmatrix_diag_ss)
    if (allocated(traj%phases_s))                   deallocate(traj%phases_s)
    if (allocated(traj%phases_old_s))               deallocate(traj%phases_old_s)
    if (allocated(traj%hopprob_s))                  deallocate(traj%hopprob_s)
    if (allocated(traj%switchprob_s))               deallocate(traj%switchprob_s)
    if (allocated(traj%gpreprob_s3))                deallocate(traj%gpreprob_s3)
    if (allocated(traj%preprob_s3))                 deallocate(traj%preprob_s3)
    if (allocated(traj%preprob_old_s3))             deallocate(traj%preprob_old_s3)
    if (allocated(traj%decotime_diag_s))            deallocate(traj%decotime_diag_s)
    if (allocated(traj%decotime_diag_old_s))        deallocate(traj%decotime_diag_old_s)
    if (allocated(traj%decotime_MCH_s))             deallocate(traj%decotime_MCH_s)
    if (allocated(traj%decotime_MCH_old_s))         deallocate(traj%decotime_MCH_old_s)
    if (allocated(traj%svec_diag_sad))              deallocate(traj%svec_diag_sad)
    if (allocated(traj%svec_MCH_sad))               deallocate(traj%svec_MCH_sad)
    if (allocated(traj%psvec_diag_sad))             deallocate(traj%psvec_diag_sad)
    if (allocated(traj%psvec_MCH_sad))              deallocate(traj%psvec_MCH_sad)
    if (allocated(traj%Property2d_xss))             deallocate(traj%Property2d_xss)
    if (allocated(traj%Property1d_ys))              deallocate(traj%Property1d_ys)
    if (allocated(traj%Property2d_labels_x))        deallocate(traj%Property2d_labels_x)
    if (allocated(traj%Property1d_labels_y))        deallocate(traj%Property1d_labels_y)
    if (allocated(traj%hopping_direction_ssad))     deallocate(traj%hopping_direction_ssad)
    if (allocated(traj%frustrated_hop_vec_ssad))    deallocate(traj%frustrated_hop_vec_ssad)
    if (allocated(traj%Gmatrix_ssad))               deallocate(traj%Gmatrix_ssad)
    if (allocated(traj%Gmatrix_old_ssad))           deallocate(traj%Gmatrix_old_ssad)
    if (allocated(traj%NACGV_MCH_ssad))             deallocate(traj%NACGV_MCH_ssad)
    if (allocated(traj%NACGV_diag_ssad))            deallocate(traj%NACGV_diag_ssad)
    if (allocated(traj%pNACGV_MCH_ssad))            deallocate(traj%pNACGV_MCH_ssad)
    if (allocated(traj%pNACGV_diag_ssad))           deallocate(traj%pNACGV_diag_ssad)
    if (allocated(traj%dendt_MCH_ss))               deallocate(traj%dendt_MCH_ss)
    if (allocated(traj%dendt_diag_ss))              deallocate(traj%dendt_diag_ss)
    if (allocated(traj%grad_ad))                    deallocate(traj%grad_ad)
    if (allocated(traj%decograd_ad))                deallocate(traj%decograd_ad)
    if (allocated(traj%coeff_diag_s))               deallocate(traj%coeff_diag_s)
    if (allocated(traj%coeff_diag_old_s))           deallocate(traj%coeff_diag_old_s)
    if (allocated(traj%coeff_MCH_s))                deallocate(traj%coeff_MCH_s)
    if (allocated(traj%coeff_MCH_old_s))            deallocate(traj%coeff_MCH_old_s)
    if (allocated(traj%ccoeff_MCH_s))               deallocate(traj%ccoeff_MCH_s)
    if (allocated(traj%ccoeff_diag_s))              deallocate(traj%ccoeff_diag_s)
    if (allocated(traj%ccoeff_diag_old_s))          deallocate(traj%ccoeff_diag_old_s)
    if (allocated(traj%coeff_zpe_s2))               deallocate(traj%coeff_zpe_s2)
    if (allocated(traj%gRcoeff_MCH_s))              deallocate(traj%gRcoeff_MCH_s)
    if (allocated(traj%gRccoeff_MCH_s))             deallocate(traj%gRccoeff_MCH_s)
    if (allocated(traj%gIcoeff_MCH_s))              deallocate(traj%gIcoeff_MCH_s)
    if (allocated(traj%gIccoeff_MCH_s))             deallocate(traj%gIccoeff_MCH_s)
    if (allocated(traj%ephase_s))                   deallocate(traj%ephase_s)
    if (allocated(traj%gephase_s))                  deallocate(traj%gephase_s)
    if (allocated(traj%selG_s))                     deallocate(traj%selG_s)
    if (allocated(traj%selT_ss))                    deallocate(traj%selT_ss)
    if (allocated(traj%state_pumping_s))            deallocate(traj%state_pumping_s)
    if (allocated(traj%pumping_status_s))           deallocate(traj%pumping_status_s)
    if (allocated(traj%uncertainty_time_s))         deallocate(traj%uncertainty_time_s)
    if (allocated(traj%allowed_hop_s))              deallocate(traj%allowed_hop_s)
    if (allocated(traj%allowed_time_s))             deallocate(traj%allowed_time_s)
    if (allocated(traj%U_pointer_ss))               deallocate(traj%U_pointer_ss)
    if (allocated(traj%U_pointer_old_ss))           deallocate(traj%U_pointer_old_ss)
    if (allocated(traj%lpzpe_ke_ah))                deallocate(traj%lpzpe_ke_ah)
    if (allocated(traj%lpzpe_ke_bc))                deallocate(traj%lpzpe_ke_bc)
    if (allocated(traj%thermostat_random))          deallocate(traj%thermostat_random)
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
#ifdef __PYSHARC__
    ! Pointer routines
    write(u,'(A20,1X,L1)') 'H_MCH_ss',        associated(traj%H_MCH_ss        )
    write(u,'(A20,1X,L1)') 'DM_ssd',          associated(traj%DM_ssd          )
    write(u,'(A20,1X,L1)') 'overlaps_ss',     associated(traj%overlaps_ss     )
    write(u,'(A20,1X,L1)') 'grad_MCH_sad',    associated(traj%grad_MCH_sad    )
    write(u,'(A20,1X,L1)') 'NACdR_ssad',      associated(traj%NACdR_ssad      )
    write(u,'(A20,1X,L1)') 'geom_ad',         associated(traj%geom_ad         )
#else
    write(u,'(A20,1X,L1)') 'H_MCH_ss',        allocated(traj%H_MCH_ss        )
    write(u,'(A20,1X,L1)') 'DM_ssd',          allocated(traj%DM_ssd          )
    write(u,'(A20,1X,L1)') 'overlaps_ss',     allocated(traj%overlaps_ss     )
    write(u,'(A20,1X,L1)') 'grad_MCH_sad',    allocated(traj%grad_MCH_sad    )
    write(u,'(A20,1X,L1)') 'NACdR_ssad',      allocated(traj%NACdR_ssad      )
    write(u,'(A20,1X,L1)') 'geom_ad',         allocated(traj%geom_ad         )
#endif


    write(u,'(A20,1X,L1)') 'atomicnumber_a',  allocated(traj%atomicnumber_a  )
    write(u,'(A20,1X,L1)') 'element_a',       allocated(traj%element_a       )
    write(u,'(A20,1X,L1)') 'mass_a',          allocated(traj%mass_a          )
    write(u,'(A20,1X,L1)') 'veloc_ad',        allocated(traj%veloc_ad        )
    write(u,'(A20,1X,L1)') 'veloc_old_ad',    allocated(traj%veloc_old_ad    )
    write(u,'(A20,1X,L1)') 'veloc_app_ad',    allocated(traj%veloc_app_ad    )
    write(u,'(A20,1X,L1)') 'accel_ad',        allocated(traj%accel_ad        )
    write(u,'(A20,1X,L1)') 'dH_MCH_ss',       allocated(traj%dH_MCH_ss       )
    write(u,'(A20,1X,L1)') 'dH_MCH_old_ss',   allocated(traj%dH_MCH_old_ss   )
    write(u,'(A20,1X,L1)') 'H_MCH_old_ss',    allocated(traj%H_MCH_old_ss    )
    write(u,'(A20,1X,L1)') 'H_MCH_old2_ss',   allocated(traj%H_MCH_old2_ss   )
    write(u,'(A20,1X,L1)') 'H_MCH_old3_ss',   allocated(traj%H_MCH_old3_ss   )
    write(u,'(A20,1X,L1)') 'H_diag_ss',       allocated(traj%H_diag_ss       )
    write(u,'(A20,1X,L1)') 'H_diag_old_ss',   allocated(traj%H_diag_old_ss   )
    write(u,'(A20,1X,L1)') 'H_diag_old2_ss',  allocated(traj%H_diag_old2_ss  )
    write(u,'(A20,1X,L1)') 'H_diag_old3_ss',  allocated(traj%H_diag_old3_ss  )
    write(u,'(A20,1X,L1)') 'Evib_local_s',    allocated(traj%Evib_local_s    )
    write(u,'(A20,1X,L1)') 'Ezpe_local_s',    allocated(traj%Ezpe_local_s    )
    write(u,'(A20,1X,L1)') 'U_ss',            allocated(traj%U_ss            )
    write(u,'(A20,1X,L1)') 'U_old_ss',        allocated(traj%U_old_ss        )
    write(u,'(A20,1X,L1)') 'NACdt_ss',        allocated(traj%NACdt_ss        )
    write(u,'(A20,1X,L1)') 'NACdt_old_ss',    allocated(traj%NACdt_old_ss    )
    write(u,'(A20,1X,L1)') 'NACdt_diag_ss',        allocated(traj%NACdt_diag_ss      )
    write(u,'(A20,1X,L1)') 'NACdt_diag_old_ss',    allocated(traj%NACdt_diag_old_ss  )
    write(u,'(A20,1X,L1)') 'grad_MCH_old_sad',     allocated(traj%grad_MCH_old_sad   )
    write(u,'(A20,1X,L1)') 'grad_MCH_old2_sad',    allocated(traj%grad_MCH_old2_sad  )
    write(u,'(A20,1X,L1)') 'grad_Ezpe_local_sad',  allocated(traj%grad_Ezpe_local_sad)
    write(u,'(A20,1X,L1)') 'NACdR_old_ssad',  allocated(traj%NACdR_old_ssad  )
    write(u,'(A20,1X,L1)') 'NACdR_diag_ssad', allocated(traj%NACdR_diag_ssad )
    write(u,'(A20,1X,L1)') 'trans_rot_P',     allocated(traj%trans_rot_P     )
    write(u,'(A20,1X,L1)') 'pNACdR_MCH_ssad', allocated(traj%pNACdR_MCH_ssad )
    write(u,'(A20,1X,L1)') 'pNACdR_diag_ssad',allocated(traj%pNACdR_diag_ssad)
    write(u,'(A20,1X,L1)') 'NACGV_MCH_ssad',  allocated(traj%NACGV_MCH_ssad  )
    write(u,'(A20,1X,L1)') 'NACGV_diag_ssad', allocated(traj%NACGV_diag_ssad )
    write(u,'(A20,1X,L1)') 'pNACGV_MCH_ssad', allocated(traj%pNACGV_MCH_ssad )
    write(u,'(A20,1X,L1)') 'pNACGV_diag_ssad',allocated(traj%pNACGV_diag_ssad)
    write(u,'(A20,1X,L1)') 'dendt_MCH_ss',    allocated(traj%dendt_MCH_ss    )
    write(u,'(A20,1X,L1)') 'dendt_diag_ss',   allocated(traj%dendt_diag_ss   )
    write(u,'(A20,1X,L1)') 'DM_old_ssd',      allocated(traj%DM_old_ssd      )
    write(u,'(A20,1X,L1)') 'DM_print_ssd',    allocated(traj%DM_print_ssd    )
    write(u,'(A20,1X,L1)') 'Rtotal_ss',       allocated(traj%Rtotal_ss       )
    write(u,'(A20,1X,L1)') 'RDtotal_ss',      allocated(traj%RDtotal_ss      )
    write(u,'(A20,1X,L1)') 'Dtotal_ss',       allocated(traj%Dtotal_ss       )
    write(u,'(A20,1X,L1)') 'Kmatrix_MCH_ss',  allocated(traj%Kmatrix_MCH_ss  )
    write(u,'(A20,1X,L1)') 'Kmatrix_diag_ss', allocated(traj%Kmatrix_diag_ss )
    write(u,'(A20,1X,L1)') 'phases_s',        allocated(traj%phases_s        )
    write(u,'(A20,1X,L1)') 'phases_old_s',    allocated(traj%phases_old_s    )
    write(u,'(A20,1X,L1)') 'hopprob_s',       allocated(traj%hopprob_s       )
    write(u,'(A20,1X,L1)') 'switchprob_s',    allocated(traj%switchprob_s    )
    write(u,'(A20,1X,L1)') 'gpreprob_s3',     allocated(traj%gpreprob_s3     )
    write(u,'(A20,1X,L1)') 'preprob_s3',      allocated(traj%preprob_s3      )
    write(u,'(A20,1X,L1)') 'preprob_old_s3',  allocated(traj%preprob_old_s3  )
    write(u,'(A20,1X,L1)') 'decotime_diag_s',      allocated(traj%decotime_diag_s      )
    write(u,'(A20,1X,L1)') 'decotime_diag_old_s',  allocated(traj%decotime_diag_old_s  )
    write(u,'(A20,1X,L1)') 'decotime_MCH_s',       allocated(traj%decotime_MCH_s       )
    write(u,'(A20,1X,L1)') 'decotime_MCH_old_s',   allocated(traj%decotime_MCH_old_s   )
    write(u,'(A20,1X,L1)') 'svec_diag_sad',        allocated(traj%svec_diag_sad        )
    write(u,'(A20,1X,L1)') 'svec_MCH_sad',         allocated(traj%svec_MCH_sad         )
    write(u,'(A20,1X,L1)') 'psvec_diag_sad',       allocated(traj%psvec_diag_sad       )
    write(u,'(A20,1X,L1)') 'psvec_MCH_sad',        allocated(traj%psvec_MCH_sad        )
    write(u,'(A20,1X,L1)') 'Property2d_xss',  allocated(traj%Property2d_xss     )
    write(u,'(A20,1X,L1)') 'Property1d_ys',   allocated(traj%Property1d_ys      )
    write(u,'(A20,1X,L1)') 'Property2d_labels_x',     allocated(traj%Property2d_labels_x     )
    write(u,'(A20,1X,L1)') 'Property1d_labels_y',     allocated(traj%Property1d_labels_y     )
    write(u,'(A20,1X,L1)') 'hopping_direction_ssad',  allocated(traj%hopping_direction_ssad  )
    write(u,'(A20,1X,L1)') 'frustrated_hop_vec_ssad', allocated(traj%frustrated_hop_vec_ssad )
    write(u,'(A20,1X,L1)') 'Gmatrix_ssad',     allocated(traj%Gmatrix_ssad     )
    write(u,'(A20,1X,L1)') 'Gmatrix_old__ssad',allocated(traj%Gmatrix_old_ssad )
    write(u,'(A20,1X,L1)') 'grad_ad',          allocated(traj%grad_ad          )
    write(u,'(A20,1X,L1)') 'decograd_ad',      allocated(traj%decograd_ad      )
    write(u,'(A20,1X,L1)') 'coeff_diag_s',     allocated(traj%coeff_diag_s     )
    write(u,'(A20,1X,L1)') 'coeff_diag_old_s', allocated(traj%coeff_diag_old_s )
    write(u,'(A20,1X,L1)') 'coeff_MCH_s',      allocated(traj%coeff_MCH_s      )
    write(u,'(A20,1X,L1)') 'coeff_MCH_old_s',  allocated(traj%coeff_MCH_old_s  )
    write(u,'(A20,1X,L1)') 'ccoeff_MCH_s',     allocated(traj%ccoeff_MCH_s     )
    write(u,'(A20,1X,L1)') 'ccoeff_diag_s',    allocated(traj%ccoeff_diag_s    )
    write(u,'(A20,1X,L1)') 'ccoeff_diag_old_s',allocated(traj%ccoeff_diag_old_s)
    write(u,'(A20,1X,L1)') 'coeff_zpe_s2',     allocated(traj%coeff_zpe_s2     )
    write(u,'(A20,1X,L1)') 'gRcoeff_MCH_s',    allocated(traj%gRcoeff_MCH_s    )
    write(u,'(A20,1X,L1)') 'gRccoeff_MCH_s',   allocated(traj%gRccoeff_MCH_s   )
    write(u,'(A20,1X,L1)') 'gIcoeff_MCH_s',    allocated(traj%gIcoeff_MCH_s    )
    write(u,'(A20,1X,L1)') 'gIccoeff_MCH_s',   allocated(traj%gIccoeff_MCH_s   )
    write(u,'(A20,1X,L1)') 'ephase_s',         allocated(traj%ephase_s         )
    write(u,'(A20,1X,L1)') 'gephase_s',        allocated(traj%gephase_s        )
    write(u,'(A20,1X,L1)') 'selG_s',           allocated(traj%selG_s           )
    write(u,'(A20,1X,L1)') 'selT_ss',          allocated(traj%selT_ss          )
    write(u,'(A20,1X,L1)') 'state_pumping_s',  allocated(traj%state_pumping_s  )
    write(u,'(A20,1X,L1)') 'pumping_status_s', allocated(traj%pumping_status_s )
    write(u,'(A20,1X,L1)') 'uncertainty_time_s',      allocated(traj%uncertainty_time_s     )
    write(u,'(A20,1X,L1)') 'allowed_hop_s',    allocated(traj%allowed_hop_s    )
    write(u,'(A20,1X,L1)') 'allowed_time_s',   allocated(traj%allowed_time_s   )
    write(u,'(A20,1X,L1)') 'U_pointer_ss',     allocated(traj%U_pointer_ss      )
    write(u,'(A20,1X,L1)') 'U_pointer_old_ss', allocated(traj%U_pointer_old_ss  )
    write(u,'(A20,1X,L1)') 'lpzpe_ah',         allocated(ctrl%lpzpe_ah         )
    write(u,'(A20,1X,L1)') 'lpzpe_bc',         allocated(ctrl%lpzpe_bc         )
    write(u,'(A20,1X,L1)') 'lpzpe_ke_zpe_ah',  allocated(ctrl%lpzpe_ke_zpe_ah  )
    write(u,'(A20,1X,L1)') 'lpzpe_ke_zpe_bc',  allocated(ctrl%lpzpe_ke_zpe_bc  )
    write(u,'(A20,1X,L1)') 'lpzpe_ke_ah',      allocated(traj%lpzpe_ke_ah      )
    write(u,'(A20,1X,L1)') 'lpzpe_ke_bc',      allocated(traj%lpzpe_ke_bc      )
    write(u,'(A20,1X,L1)') 'nstates_m',        allocated(ctrl%nstates_m        )
    write(u,'(A20,1X,L1)') 'actstates_s',      allocated(ctrl%actstates_s      )

    write(u,*) '_______________________ CHECKING NaNs _______________________'

    write(u,'(A20,1X,L1)') 'mass_a',          any((traj%mass_a          ).ne.(traj%mass_a          ))
    write(u,'(A20,1X,L1)') 'geom_ad',         any((traj%geom_ad         ).ne.(traj%geom_ad         ))
    write(u,'(A20,1X,L1)') 'veloc_ad',        any((traj%veloc_ad        ).ne.(traj%veloc_ad        ))
    write(u,'(A20,1X,L1)') 'veloc_old_ad',    any((traj%veloc_old_ad    ).ne.(traj%veloc_old_ad    ))
    write(u,'(A20,1X,L1)') 'veloc_app_ad',    any((traj%veloc_app_ad    ).ne.(traj%veloc_app_ad    ))
    write(u,'(A20,1X,L1)') 'accel_ad',        any((traj%accel_ad        ).ne.(traj%accel_ad        ))
    write(u,'(A20,1X,L1)') 'NACdR_ssad',      any((traj%NACdR_ssad      ).ne.(traj%NACdR_ssad      ))
    write(u,'(A20,1X,L1)') 'NACdR_old_ssad',  any((traj%NACdR_old_ssad  ).ne.(traj%NACdR_old_ssad  ))
    write(u,'(A20,1X,L1)') 'trans_rot_P',     any((traj%trans_rot_P     ).ne.(traj%trans_rot_P     ))
    write(u,'(A20,1X,L1)') 'pNACdR_MCH_ssad', any((traj%pNACdR_MCH_ssad ).ne.(traj%pNACdR_MCH_ssad ))
    write(u,'(A20,1X,L1)') 'hopprob_s',       any((traj%hopprob_s       ).ne.(traj%hopprob_s       ))
    write(u,'(A20,1X,L1)') 'switchprob_s',    any((traj%switchprob_s    ).ne.(traj%switchprob_s    ))
    write(u,'(A20,1X,L1)') 'gpreprob_s3',     any((traj%gpreprob_s3     ).ne.(traj%gpreprob_s3     ))
    write(u,'(A20,1X,L1)') 'preprob_s3',      any((traj%preprob_s3      ).ne.(traj%preprob_s3      ))
    write(u,'(A20,1X,L1)') 'preprob_old_s3',  any((traj%preprob_old_s3  ).ne.(traj%preprob_old_s3  ))
    write(u,'(A20,1X,L1)') 'decotime_diag_s',      any((traj%decotime_diag_s      ).ne.(traj%decotime_diag_s      ))
    write(u,'(A20,1X,L1)') 'decotime_diag_old_s',  any((traj%decotime_diag_old_s  ).ne.(traj%decotime_diag_old_s  ))
    write(u,'(A20,1X,L1)') 'decotime_MCH_s',       any((traj%decotime_MCH_s       ).ne.(traj%decotime_MCH_s       ))
    write(u,'(A20,1X,L1)') 'decotime_MCH_old_s',   any((traj%decotime_MCH_old_s   ).ne.(traj%decotime_MCH_old_s   ))
    write(u,'(A20,1X,L1)') 'svec_diag_sad',        any((traj%svec_diag_sad        ).ne.(traj%svec_diag_sad        ))
    write(u,'(A20,1X,L1)') 'svec_MCH_sad',         any((traj%svec_MCH_sad         ).ne.(traj%svec_MCH_sad         ))
    write(u,'(A20,1X,L1)') 'psvec_diag_sad',       any((traj%psvec_diag_sad       ).ne.(traj%psvec_diag_sad       ))
    write(u,'(A20,1X,L1)') 'psvec_MCH_sad',        any((traj%psvec_MCH_sad        ).ne.(traj%psvec_MCH_sad        ))
    write(u,'(A20,1X,L1)') 'hopping_direction_ssad', any((traj%hopping_direction_ssad  ).ne.(traj%hopping_direction_ssad  ))
    write(u,'(A20,1X,L1)') 'frustrated_hop_vec_ssad',any((traj%frustrated_hop_vec_ssad ).ne.(traj%frustrated_hop_vec_ssad ))
    write(u,'(A20,1X,L1)') 'grad_MCH_sad',         any((traj%grad_MCH_sad         ).ne.(traj%grad_MCH_sad         ))
    write(u,'(A20,1X,L1)') 'grad_MCH_old_sad',     any((traj%grad_MCH_old_sad     ).ne.(traj%grad_MCH_old_sad     ))
    write(u,'(A20,1X,L1)') 'grad_MCH_old2_sad',    any((traj%grad_MCH_old2_sad    ).ne.(traj%grad_MCH_old2_sad    ))
    write(u,'(A20,1X,L1)') 'grad_Ezpe_local_sad',  any((traj%grad_Ezpe_local_sad  ).ne.(traj%grad_Ezpe_local_sad  ))
    write(u,'(A20,1X,L1)') 'grad_ad',         any((traj%grad_ad         ).ne.(traj%grad_ad         ))
    write(u,'(A20,1X,L1)') 'decograd_ad',     any((traj%decograd_ad     ).ne.(traj%decograd_ad     ))
    write(u,'(A20,1X,L1)') 'gRcoeff_MCH_s',   any((traj%gRcoeff_MCH_s   ).ne.(traj%gRcoeff_MCH_s   ))
    write(u,'(A20,1X,L1)') 'gRccoeff_MCH_s',  any((traj%gRccoeff_MCH_s  ).ne.(traj%gRccoeff_MCH_s  ))
    write(u,'(A20,1X,L1)') 'gIcoeff_MCH_s',   any((traj%gIcoeff_MCH_s   ).ne.(traj%gIcoeff_MCH_s   ))
    write(u,'(A20,1X,L1)') 'gIccoeff_MCH_s',  any((traj%gIccoeff_MCH_s  ).ne.(traj%gIccoeff_MCH_s  ))
    write(u,'(A20,1X,L1)') 'ephase_s',        any((traj%ephase_s        ).ne.(traj%ephase_s        ))
    write(u,'(A20,1X,L1)') 'gephase_s',       any((traj%gephase_s       ).ne.(traj%gephase_s       ))
    write(u,'(A20,1X,L1)') 'state_pumping_s', any((traj%state_pumping_s ).ne.(traj%state_pumping_s ))
    write(u,'(A20,1X,L1)') 'pumping_status_s',any((traj%pumping_status_s).ne.(traj%pumping_status_s))
    write(u,'(A20,1X,L1)') 'uncertainty_time_s',   any((traj%uncertainty_time_s   ).ne.(traj%uncertainty_time_s   ))
    write(u,'(A20,1X,L1)') 'allowed_hop_s',        any((traj%allowed_hop_s        ).ne.(traj%allowed_hop_s        ))
    write(u,'(A20,1X,L1)') 'allowed_time_s',       any((traj%allowed_time_s       ).ne.(traj%allowed_time_s       ))
    write(u,'(A20,1X,L1)') 'lpzpe_ke_ah',     any((traj%lpzpe_ke_ah     ).ne.(traj%lpzpe_ke_ah     ))
    write(u,'(A20,1X,L1)') 'lpzpe_ke_bc',     any((traj%lpzpe_ke_bc     ).ne.(traj%lpzpe_ke_bc     ))

    write(u,*) 'Real parts:'
    write(u,'(A20,1X,L1)') 'H_MCH_ss',        any((real(traj%H_MCH_ss        )).ne.(real(traj%H_MCH_ss        )))
    write(u,'(A20,1X,L1)') 'dH_MCH_ss',       any((real(traj%dH_MCH_ss       )).ne.(real(traj%dH_MCH_ss       )))
    write(u,'(A20,1X,L1)') 'dH_MCH_old_ss',   any((real(traj%dH_MCH_old_ss   )).ne.(real(traj%dH_MCH_old_ss   )))
    write(u,'(A20,1X,L1)') 'H_MCH_old_ss',    any((real(traj%H_MCH_old_ss    )).ne.(real(traj%H_MCH_old_ss    )))
    write(u,'(A20,1X,L1)') 'H_MCH_old2_ss',    any((real(traj%H_MCH_old2_ss  )).ne.(real(traj%H_MCH_old2_ss   )))
    write(u,'(A20,1X,L1)') 'H_MCH_old3_ss',   any((real(traj%H_MCH_old3_ss   )).ne.(real(traj%H_MCH_old3_ss   )))
    write(u,'(A20,1X,L1)') 'H_diag_ss',       any((real(traj%H_diag_ss       )).ne.(real(traj%H_diag_ss       )))
    write(u,'(A20,1X,L1)') 'H_diag_old_ss',   any((real(traj%H_diag_old_ss   )).ne.(real(traj%H_diag_old_ss   )))
    write(u,'(A20,1X,L1)') 'H_diag_old2_ss',  any((real(traj%H_diag_old2_ss  )).ne.(real(traj%H_diag_old2_ss  )))
    write(u,'(A20,1X,L1)') 'H_diag_old3_ss',  any((real(traj%H_diag_old3_ss  )).ne.(real(traj%H_diag_old3_ss  )))
    write(u,'(A20,1X,L1)') 'Evib_local_s',    any((real(traj%Evib_local_s    )).ne.(real(traj%Evib_local_s    )))
    write(u,'(A20,1X,L1)') 'Ezpe_local_s',    any((real(traj%Ezpe_local_s    )).ne.(real(traj%Ezpe_local_s    )))
    write(u,'(A20,1X,L1)') 'U_ss',            any((real(traj%U_ss            )).ne.(real(traj%U_ss            )))
    write(u,'(A20,1X,L1)') 'U_old_ss',        any((real(traj%U_old_ss        )).ne.(real(traj%U_old_ss        )))
    write(u,'(A20,1X,L1)') 'NACdt_ss',        any((real(traj%NACdt_ss        )).ne.(real(traj%NACdt_ss        )))
    write(u,'(A20,1X,L1)') 'NACdt_old_ss',    any((real(traj%NACdt_old_ss    )).ne.(real(traj%NACdt_old_ss    )))
    write(u,'(A20,1X,L1)') 'NACdt_diag_ss',        any((real(traj%NACdt_diag_ss     )).ne.(real(traj%NACdt_diag_ss     )))
    write(u,'(A20,1X,L1)') 'NACdt_diag_old_ss',    any((real(traj%NACdt_diag_old_ss )).ne.(real(traj%NACdt_diag_old_ss )))
    write(u,'(A20,1X,L1)') 'NACdR_diag_ssad', any((real(traj%NACdR_diag_ssad )).ne.(real(traj%NACdR_diag_ssad )))
    write(u,'(A20,1X,L1)') 'overlaps_ss',     any((real(traj%overlaps_ss     )).ne.(real(traj%overlaps_ss     )))
    write(u,'(A20,1X,L1)') 'DM_ssd',          any((real(traj%DM_ssd          )).ne.(real(traj%DM_ssd          )))
    write(u,'(A20,1X,L1)') 'DM_old_ssd',      any((real(traj%DM_old_ssd      )).ne.(real(traj%DM_old_ssd      )))
    write(u,'(A20,1X,L1)') 'Property2d_xss',  any((real(traj%Property2d_xss  )).ne.(real(traj%Property2d_xss  )))
    write(u,'(A20,1X,L1)') 'Property1d_ys',   any((real(traj%Property1d_ys   )).ne.(real(traj%Property1d_ys   )))
    write(u,'(A20,1X,L1)') 'DM_print_ssd',    any((real(traj%DM_print_ssd    )).ne.(real(traj%DM_print_ssd    )))
    write(u,'(A20,1X,L1)') 'Rtotal_ss',       any((real(traj%Rtotal_ss       )).ne.(real(traj%Rtotal_ss       )))
    write(u,'(A20,1X,L1)') 'RDtotal_ss',      any((real(traj%RDtotal_ss      )).ne.(real(traj%RDtotal_ss      )))
    write(u,'(A20,1X,L1)') 'Dtotal_ss',       any((real(traj%Dtotal_ss       )).ne.(real(traj%Dtotal_ss       )))
    write(u,'(A20,1X,L1)') 'Kmatrix_MCH_ss',  any((real(traj%Kmatrix_MCH_ss  )).ne.(real(traj%Kmatrix_MCH_ss  )))
    write(u,'(A20,1X,L1)') 'Kmatrix_diag_ss', any((real(traj%Kmatrix_diag_ss )).ne.(real(traj%Kmatrix_diag_ss )))
    write(u,'(A20,1X,L1)') 'phases_s',        any((real(traj%phases_s        )).ne.(real(traj%phases_s        )))
    write(u,'(A20,1X,L1)') 'phases_old_s',    any((real(traj%phases_old_s    )).ne.(real(traj%phases_old_s    )))
    write(u,'(A20,1X,L1)') 'pNACdR_diag_ssad',any((real(traj%pNACdR_diag_ssad)).ne.(real(traj%pNACdR_diag_ssad)))
    write(u,'(A20,1X,L1)') 'Gmatrix_ssad',    any((real(traj%Gmatrix_ssad    )).ne.(real(traj%Gmatrix_ssad    )))
    write(u,'(A20,1X,L1)') 'Gmatrix_old_ssad',any((real(traj%Gmatrix_old_ssad)).ne.(real(traj%Gmatrix_old_ssad)))
    write(u,'(A20,1X,L1)') 'NACGV_MCH_ssad',  any((real(traj%NACGV_MCH_ssad  )).ne.(real(traj%NACGV_MCH_ssad  )))
    write(u,'(A20,1X,L1)') 'NACGV_diag_ssad', any((real(traj%NACGV_diag_ssad )).ne.(real(traj%NACGV_diag_ssad )))
    write(u,'(A20,1X,L1)') 'pNACGV_MCH_ssad', any((real(traj%pNACGV_MCH_ssad )).ne.(real(traj%pNACGV_MCH_ssad )))
    write(u,'(A20,1X,L1)') 'pNACGV_diag_ssad',any((real(traj%pNACGV_diag_ssad)).ne.(real(traj%pNACGV_diag_ssad)))
    write(u,'(A20,1X,L1)') 'dendt_MCH_ss',    any((real(traj%dendt_MCH_ss    )).ne.(real(traj%dendt_MCH_ss    )))
    write(u,'(A20,1X,L1)') 'dendt_diag_ss',   any((real(traj%dendt_diag_ss   )).ne.(real(traj%dendt_diag_ss   )))
    write(u,'(A20,1X,L1)') 'coeff_diag_s',    any((real(traj%coeff_diag_s    )).ne.(real(traj%coeff_diag_s    )))
    write(u,'(A20,1X,L1)') 'coeff_diag_old_s',any((real(traj%coeff_diag_old_s)).ne.(real(traj%coeff_diag_old_s)))
    write(u,'(A20,1X,L1)') 'coeff_MCH_s',     any((real(traj%coeff_MCH_s     )).ne.(real(traj%coeff_MCH_s     )))
    write(u,'(A20,1X,L1)') 'coeff_MCH_old_s', any((real(traj%coeff_MCH_old_s )).ne.(real(traj%coeff_MCH_old_s )))
    write(u,'(A20,1X,L1)') 'ccoeff_MCH_s',    any((real(traj%ccoeff_MCH_s    )).ne.(real(traj%ccoeff_MCH_s    )))
    write(u,'(A20,1X,L1)') 'ccoeff_diag_s',   any((real(traj%ccoeff_diag_s   )).ne.(real(traj%ccoeff_diag_s   )))
    write(u,'(A20,1X,L1)') 'ccoeff_diag_old_s',any((real(traj%ccoeff_diag_old_s)).ne.(real(traj%ccoeff_diag_old_s)))
    write(u,'(A20,1X,L1)') 'coeff_zpe_s2',    any((real(traj%coeff_zpe_s2    )).ne.(real(traj%coeff_zpe_s2    )))
    write(u,'(A20,1X,L1)') 'U_pointer_ss',    any((real(traj%U_pointer_ss    )).ne.(real(traj%U_pointer_ss    )))
    write(u,'(A20,1X,L1)') 'U_pointer_oldss', any((real(traj%U_pointer_old_ss)).ne.(real(traj%U_pointer_old_ss)))
    write(u,*) 'Imag parts:'
    write(u,'(A20,1X,L1)') 'H_MCH_ss',        any((aimag(traj%H_MCH_ss        )).ne.(aimag(traj%H_MCH_ss        )))
    write(u,'(A20,1X,L1)') 'dH_MCH_ss',       any((aimag(traj%dH_MCH_ss       )).ne.(aimag(traj%dH_MCH_ss       )))
    write(u,'(A20,1X,L1)') 'dH_MCH_old_ss',   any((aimag(traj%dH_MCH_old_ss   )).ne.(aimag(traj%dH_MCH_old_ss   )))
    write(u,'(A20,1X,L1)') 'H_MCH_old_ss',    any((aimag(traj%H_MCH_old_ss    )).ne.(aimag(traj%H_MCH_old_ss    )))
    write(u,'(A20,1X,L1)') 'H_MCH_old2_ss',   any((aimag(traj%H_MCH_old2_ss   )).ne.(aimag(traj%H_MCH_old2_ss   )))
    write(u,'(A20,1X,L1)') 'H_MCH_old3_ss',   any((aimag(traj%H_MCH_old3_ss   )).ne.(aimag(traj%H_MCH_old3_ss   )))
    write(u,'(A20,1X,L1)') 'H_diag_ss',       any((aimag(traj%H_diag_ss       )).ne.(aimag(traj%H_diag_ss       )))
    write(u,'(A20,1X,L1)') 'H_diag_old_ss',   any((aimag(traj%H_diag_old_ss   )).ne.(aimag(traj%H_diag_old_ss   )))
    write(u,'(A20,1X,L1)') 'H_diag_old2_ss',  any((aimag(traj%H_diag_old2_ss  )).ne.(aimag(traj%H_diag_old2_ss  )))
    write(u,'(A20,1X,L1)') 'H_diag_old3_ss',  any((aimag(traj%H_diag_old3_ss  )).ne.(aimag(traj%H_diag_old3_ss  )))
    write(u,'(A20,1X,L1)') 'Evib_local_s',    any((aimag(traj%Evib_local_s    )).ne.(aimag(traj%Evib_local_s    )))
    write(u,'(A20,1X,L1)') 'Ezpe_local_s',    any((aimag(traj%Ezpe_local_s    )).ne.(aimag(traj%Ezpe_local_s    )))
    write(u,'(A20,1X,L1)') 'U_ss',            any((aimag(traj%U_ss            )).ne.(aimag(traj%U_ss            )))
    write(u,'(A20,1X,L1)') 'U_old_ss',        any((aimag(traj%U_old_ss        )).ne.(aimag(traj%U_old_ss        )))
    write(u,'(A20,1X,L1)') 'NACdt_ss',        any((aimag(traj%NACdt_ss        )).ne.(aimag(traj%NACdt_ss        )))
    write(u,'(A20,1X,L1)') 'NACdt_old_ss',    any((aimag(traj%NACdt_old_ss    )).ne.(aimag(traj%NACdt_old_ss    )))
    write(u,'(A20,1X,L1)') 'NACdt_diag_ss',        any((aimag(traj%NACdt_diag_ss      )).ne.(aimag(traj%NACdt_diag_ss     )))
    write(u,'(A20,1X,L1)') 'NACdt_diag_old_ss',    any((aimag(traj%NACdt_diag_old_ss  )).ne.(aimag(traj%NACdt_diag_old_ss )))
    write(u,'(A20,1X,L1)') 'NACdR_diag_ssad', any((aimag(traj%NACdR_diag_ssad )).ne.(aimag(traj%NACdR_diag_ssad )))
    write(u,'(A20,1X,L1)') 'overlaps_ss',     any((aimag(traj%overlaps_ss     )).ne.(aimag(traj%overlaps_ss     )))
    write(u,'(A20,1X,L1)') 'DM_ssd',          any((aimag(traj%DM_ssd          )).ne.(aimag(traj%DM_ssd          )))
    write(u,'(A20,1X,L1)') 'DM_old_ssd',      any((aimag(traj%DM_old_ssd      )).ne.(aimag(traj%DM_old_ssd      )))
    write(u,'(A20,1X,L1)') 'DM_print_ssd',    any((aimag(traj%DM_print_ssd    )).ne.(aimag(traj%DM_print_ssd    )))
    write(u,'(A20,1X,L1)') 'Property2d_xss',  any((aimag(traj%Property2d_xss  )).ne.(aimag(traj%Property2d_xss  )))
!     write(u,'(A20,1X,L1)') 'Property1d_ys',   any((aimag(traj%Property1d_ys   )).ne.(aimag(traj%Property1d_ys   )))
    write(u,'(A20,1X,L1)') 'Rtotal_ss',       any((aimag(traj%Rtotal_ss       )).ne.(aimag(traj%Rtotal_ss       )))
    write(u,'(A20,1X,L1)') 'RDtotal_ss',      any((aimag(traj%RDtotal_ss      )).ne.(aimag(traj%RDtotal_ss      )))
    write(u,'(A20,1X,L1)') 'Dtotal_ss',       any((aimag(traj%Dtotal_ss       )).ne.(aimag(traj%Dtotal_ss       )))
    write(u,'(A20,1X,L1)') 'Kmatrix_MCH_ss',  any((aimag(traj%Kmatrix_MCH_ss  )).ne.(aimag(traj%Kmatrix_MCH_ss  )))
    write(u,'(A20,1X,L1)') 'Kmatrix_diag_ss', any((aimag(traj%Kmatrix_diag_ss )).ne.(aimag(traj%Kmatrix_diag_ss )))
    write(u,'(A20,1X,L1)') 'phases_s',        any((aimag(traj%phases_s        )).ne.(aimag(traj%phases_s        )))
    write(u,'(A20,1X,L1)') 'phases_old_s',    any((aimag(traj%phases_old_s    )).ne.(aimag(traj%phases_old_s    )))
    write(u,'(A20,1X,L1)') 'pNACdR_diag_ssad',any((aimag(traj%pNACdR_diag_ssad)).ne.(aimag(traj%pNACdR_diag_ssad)))
    write(u,'(A20,1X,L1)') 'Gmatrix_ssad',    any((aimag(traj%Gmatrix_ssad    )).ne.(aimag(traj%Gmatrix_ssad    )))
    write(u,'(A20,1X,L1)') 'Gmatrix_old_ssad',any((aimag(traj%Gmatrix_old_ssad)).ne.(aimag(traj%Gmatrix_old_ssad)))
    write(u,'(A20,1X,L1)') 'NACGV_MCH_ssad',  any((aimag(traj%NACGV_MCH_ssad  )).ne.(aimag(traj%NACGV_MCH_ssad  )))
    write(u,'(A20,1X,L1)') 'NACGV_diag_ssad', any((aimag(traj%NACGV_diag_ssad )).ne.(aimag(traj%NACGV_diag_ssad )))
    write(u,'(A20,1X,L1)') 'pNACGV_MCH_ssad', any((aimag(traj%pNACGV_MCH_ssad )).ne.(aimag(traj%pNACGV_MCH_ssad )))
    write(u,'(A20,1X,L1)') 'pNACGV_diag_ssad',any((aimag(traj%pNACGV_diag_ssad)).ne.(aimag(traj%pNACGV_diag_ssad)))
    write(u,'(A20,1X,L1)') 'dendt_MCH_ss',    any((aimag(traj%dendt_MCH_ss    )).ne.(aimag(traj%dendt_MCH_ss    )))
    write(u,'(A20,1X,L1)') 'dendt_diag_ss',   any((aimag(traj%dendt_diag_ss   )).ne.(aimag(traj%dendt_diag_ss   )))
    write(u,'(A20,1X,L1)') 'coeff_diag_s',    any((aimag(traj%coeff_diag_s    )).ne.(aimag(traj%coeff_diag_s    )))
    write(u,'(A20,1X,L1)') 'coeff_diag_old_s',any((aimag(traj%coeff_diag_old_s)).ne.(aimag(traj%coeff_diag_old_s)))
    write(u,'(A20,1X,L1)') 'coeff_MCH_s',     any((aimag(traj%coeff_MCH_s     )).ne.(aimag(traj%coeff_MCH_s     )))
    write(u,'(A20,1X,L1)') 'coeff_MCH_old_s', any((aimag(traj%coeff_MCH_old_s )).ne.(aimag(traj%coeff_MCH_old_s )))
    write(u,'(A20,1X,L1)') 'ccoeff_MCH_s',    any((aimag(traj%ccoeff_MCH_s    )).ne.(aimag(traj%ccoeff_MCH_s    )))
    write(u,'(A20,1X,L1)') 'ccoeff_diag_s',   any((aimag(traj%ccoeff_diag_s   )).ne.(aimag(traj%ccoeff_diag_s   )))
    write(u,'(A20,1X,L1)') 'ccoeff_diag_old_s',any((aimag(traj%ccoeff_diag_old_s)).ne.(aimag(traj%ccoeff_diag_old_s)))
    write(u,'(A20,1X,L1)') 'coeff_zpe_s2',    any((aimag(traj%coeff_zpe_s2    )).ne.(aimag(traj%coeff_zpe_s2    )))
    write(u,'(A20,1X,L1)') 'U_pointer_ss',    any((aimag(traj%U_pointer_ss    )).ne.(aimag(traj%U_pointer_ss    )))
    write(u,'(A20,1X,L1)') 'U_pointer_old_ss',any((aimag(traj%U_pointer_old_ss)).ne.(aimag(traj%U_pointer_old_ss)))

    write(u,*) '____________________________________________________________'

  endsubroutine

! =========================================================== !

endmodule
















