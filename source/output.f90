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

!> # Module OUTPUT
!>
!> \author Sebastian Mai
!> \date 27.02.2015, modified 27.02.2017 by Philipp Marquetand
!>
!>                   modified 11.13.2019 by Yinan Shu
!>                       modified subroutine write_logtimestep so it tells the trajecotry time
!>                       modified output.lis, so it prints total density
!>
!> This module defines a number of subroutines which print various
!> information.
!> 
!> The following files are written here:
!> - parts of output.log
!> - output.lis
!> - output.dat
!> - output.xyz
!>
!> Note that the print routines for the restart files are in restart.f90
!> Also note that output.log is written to by many routines from all modules, depending on printlevel.
module output

 contains

!> prints the header of the SHARC log file, containing:
!> - PWD, hostname, start time
!> - Logo
!> - Build/verson information
!> 
subroutine write_logheader(u,version)
  implicit none
  character*1023 :: hostname, cwd
  character*24 :: ctime, date
  integer :: idate, u
  character(len=*) :: version
  integer :: time

  ! build_info.inc is written by the Makefile and contains the 
  ! date and host, when/where SHARC was built
  include 'build_info.inc'

  call getcwd(cwd)
  call hostnm(hostname)
  idate=time()
  date=ctime(idate)

  call write_sharc_l(u)
  write(u,*)
  write(u,*) 'EXECUTION INFORMATION:'
  write(u,*) 'Start date: ',trim(date)
  write(u,*) 'Run host: ',trim(hostname)
  write(u,*) 'Run directory: ',trim(cwd)
  write(u,*)
  write(u,*) 'BUILD INFORMATION:'
  write(u,*) 'Build date: ',trim(build_date)
  write(u,*) 'Build host: ',trim(build_host)
  write(u,*) 'Build directory: ',trim(build_dir)
  write(u,*) 'Compiler: ',trim(build_compiler)
  write(u,*) 'PYSHARC activated: ',trim(use_pysharc)
  write(u,*)

  write(u,'(A)')      '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<============================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
  write(u,'(A)')      '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   Program SHARC started    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
  write(u,'(A)')      '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<============================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
  write(u,*)
  write(u,*) 'Version: ', trim(version)
  write(u,*) ''

endsubroutine

! ===================================================

!> writes the headline for a new timestep and the date where the step was entered
subroutine write_logtimestep(u,step,trajtime)
  use definitions
  implicit none
  character*24 :: ctime, date
  integer :: idate,time
  integer :: u, step
  real*8 :: trajtime

  idate=time()
  date=ctime(idate)
  if (printlevel>0) then
    write(u,'(A)')      '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<============================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
    write(u,'(A,I6,A)') '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  Entering timestep ',step,'  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
    write(u,'(A,F12.5,A)') '<<<<<<<<<<<<<<<<<<<<<<<<<<<<  Trajectory time in fs',trajtime*au2fs,'  >>>>>>>>>>>>>>>>>>>>>>>>>>>>'
    write(u,'(A)')      '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<============================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
    write(u,'(56X,A12,A)')     'Start time: ', trim(date)
    write(u,*)
  endif

endsubroutine

! ===================================================

!> writes the old SHARC logo from the pre-release versions
  subroutine write_sharc_r(u)
  implicit none
  integer :: u

    write(u,'(A)') ''
    write(u,'(A)') '                                                      .ggp                                   '
    write(u,'(A)') '                                                      :$WV#ap                                '
    write(u,'(A)') '                                                       ]$[-?4yw                              '
    write(u,'(A)') '                                                       ]3k  -!Uc,                            '
    write(u,'(A)') '                                                        ]A.    46                            '
    write(u,'(A)') '                                                        )d|    ))g,                          '
    write(u,'(A)') '                                                     _._)ds,,    !kp                         '
    write(u,'(A)') '                                              __aawmZUVTTTTV#mgag_3m_                        '
    write(u,'(A)') '                                         _qaamZT!"~^         -""?VWBqg,                      '
    write(u,'(A)') '                                      _qwUT!~`                    -"!T$gg.                   '
    write(u,'(A)') '                                   _gw2Y~`                            ""Y$ga_                '
    write(u,'(A)') '                             ..qa_M"!^                                    ]"?saap,           '
    write(u,'(A)') '     ____gpw_. .           _jwUF^.                                            -"4$gg,        '
    write(u,'(A)') '  _qwm#VY??YTTAgwaw_,,. ,.gwZ"`                                                  -)Tmw_      '
    write(u,'(A)') '  ]?TTA#waaw_, -!"?T9X##ZV?^`     __._w,ggagagg_,                                   ""#gp,   '
    write(u,'(A)') '        -""TTmwg,,         __gawwmmQWQVTY??!??9T$mwwag_,            _              Sa,^^$L_  '
    write(u,'(A)') '             -7?$6_:    _qwBD?!^~`]QV"`           "^??V$mmwaagg,_.. ww     :    ,   ~`  -]Wp,'
    write(u,'(A)') '                -4Wa,  jjE~`      ]^-                   -"^"!?T$QWqqWq   ]je   :mmga,,    TGc'
    write(u,'(A)') '                  "Hc, jEf                                    _jmgw#QZ   j#[    3Q@TQma.   ]#'
    write(u,'(A)') '                   |GcaF`                                     "Y?!~5QN _wD$qa,, ]3Qg94WmwawZf'
    write(u,'(A)') '                    3mm`                                           ]S vdY` ""$gg/)3Qc/!VV!"` '
    write(u,'(A)') '                    ]]V                                           jdiwF`     -9T$mWQ(        '
    write(u,'(A)') '                     -.                                          .yQU!           ^"^`        '
    write(u,'(A)') '                                                                _wWF-                        '
    write(u,'(A)') '                                                               =]T`                          '
    write(u,'(A)') '                                                                 :                           '
    write(u,'(A)') ''

  endsubroutine

! ===================================================

! writes the official SHARC logo in ASCII art
  subroutine write_sharc_l(u)
  implicit none
  integer :: u

    write(u,'(A)') ''
    write(u,'(A)') '                                ,._gz,                                        ,._\'
    write(u,'(A)') '                               .g@@@p                                       ._Q@$+'
    write(u,'(A)') '                            ,.Q@@@@@f                                    ._Q@@@+  '
    write(u,'(A)') '                           .g@@@@@@@I                                 ._g@@@@F!   '
    write(u,'(A)') '                       ,,_zQ@@@@@@@@Q_,                            ,_zQ@@@@$+     '
    write(u,'(A)') '                .__ggQ@@@@@@@@@@@@@@@@@L_.              ,         _Q@@@@@@v       '
    write(u,'(A)') '       , __zg@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@QQz_z_zzQ@@L,    .zQ@@@@@@F`        '
    write(u,'(A)') '   .__gQ@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@Q\__zQ@@@@@@$+          '
    write(u,'(A)') '  G@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@D            '
    write(u,'(A)') '    =4A@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@A@$@@@@@L,           '
    write(u,'(A)') '        =vVAA@$@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@p=VX@@@B*   <V@$@@@Q_,         '
    write(u,'(A)') '                   @@@@@@@@@@P^^^^^^--`     -q@@$=`     =4V`       =X@@@@_        '
    write(u,'(A)') '                   1@@@@@@@$!                  ^`                       -`        '
    write(u,'(A)') '                   -q@@@@@@]                                                      '
    write(u,'(A)') '                    ($@@@@$`                                                      '
    write(u,'(A)') '                     4@@@@p                                                       '
    write(u,'(A)') '                      d@@@p                                                       '
    write(u,'(A)') '                       d@@b                                                       '
    write(u,'(A)') '                        =4!                                                       '
    write(u,'(A)') ''

  endsubroutine

  ! =====================================================

!> prints the version, authors, adress, citations, license and build information
!> can be accessed by executing $SHARC/sharc.x -v | --version | --info
  subroutine write_license(u,version)
  implicit none
  integer :: u
  character(len=*) :: version
  include 'build_info.inc'

    write(u,'(A)') '**********************************************************'
    write(u,'(A)') '* SHARC -- Surface Hopping including ARbitrary Couplings *'
    write(u,'(A)') '*        Ab initio non-adiabatic dynamics program        *'
    write(u,'(A)') '**********************************************************'
    write(u,'(A)') 
    write(u,'(A,A)') 'Version: ', trim(version)
    write(u,'(A)') 
    write(u,'(A)') 
    write(u,'(A)') 'Authors: Sebastian Mai, Yinan Shu, Martin Richter, '
    write(u,'(A)') 'Matthias Ruckenbauer,Markus Oppel, Philipp Marquetand,' 
    write(u,'(A)') 'Donald Truhlar, and Leticia González'
    write(u,'(A)')
    write(u,'(A)') 'Institute of Theoretical Chemistry, University of Vienna'
    write(u,'(A)') 'Währinger Straße 17'
    write(u,'(A)') '1090 Vienna, Austria'
    write(u,'(A)')
    write(u,'(A)') 'Please cite the following:'
    write(u,'(A)') 
    write(u,'(A)') '* M. Richter, P. Marquetand, J. González-Vázquez, '
    write(u,'(A)') '  I. Sola and L. González:'
    write(u,'(A)') '  J. Chem. Theory Comput., 7, 1253--1258 (2011).'
    write(u,'(A)') 
    write(u,'(A)') '* S. Mai, M. Richter, M. Ruckenbauer,'
    write(u,'(A)') '  M. Oppel, P. Marquetand, and L. Gonzalez:'
    write(u,'(A)') '  "SHARC: Surface Hopping Including Arbitrary Couplings --'
    write(u,'(A)') '  Program Package for Non-Adiabatic Dynamics" (2014),'
    write(u,'(A)') '  sharc-md.org'
    write(u,'(A)')
    write(u,'(A)') 'Please report bugs and feature suggestions to:'
    write(u,'(A)') 'sharc@univie.ac.at'
    write(u,'(A)') 'philipp.marquetand@univie.ac.at'
    write(u,'(A)')
    write(u,'(A)') 'Copyright (c) 2023, University of Vienna'
    write(u,'(A)') 
    write(u,'(A)') 'Permission is hereby granted, free of charge, to any person obtaining a copy '
    write(u,'(A)') 'of this software and associated documentation files (the "Software"), to deal'
    write(u,'(A)') 'in the Software without restriction, including without limitation the rights '
    write(u,'(A)') 'to use, copy, modify, merge, publish, distribute, sublicense, and/or sell    '
    write(u,'(A)') 'copies of the Software, and to permit persons to whom the Software is        '
    write(u,'(A)') 'furnished to do so, subject to the following conditions:                     '
    write(u,'(A)') '                                                                             '
    write(u,'(A)') 'The above copyright notice and this permission notice shall be included in   '
    write(u,'(A)') 'all copies or substantial portions of the Software.                          '
    write(u,'(A)') '                                                                             '
    write(u,'(A)') 'The Software is provided "AS IS", without warranty of any kind, express or   '
    write(u,'(A)') 'implied, including but not limited to the warranties of merchantability,     '
    write(u,'(A)') 'fitness for a particular purpose and noninfringement. In no event shall the  '
    write(u,'(A)') 'authors or copyright holders be liable for any claim, damages or other       '
    write(u,'(A)') 'liability, whether in an action of contract, tort or otherwise, arising from,'
    write(u,'(A)') 'out of or in connection with the software or the use or other dealings in    '
    write(u,'(A)') 'the Software.                                                                '
    write(u,'(A)')
    write(u,'(A)')
    write(u,'(A)') 'BUILD INFORMATION:'
    write(u,'(A,A)') 'Build date: ',trim(build_date)
    write(u,'(A,A)') 'Build host: ',trim(build_host)
    write(u,'(A,A)') 'Build directory: ',trim(build_dir)
    write(u,'(A,A)') 'Compiler: ',trim(build_compiler)
  endsubroutine

  ! =====================================================

!> writes a formatted header for the listing file
subroutine write_list_header(u)
  implicit none
  integer :: u

  write(u,'(a1,a)') '#',repeat('=',175)
  write(u,'(a1,A11,1X,A14,1X,A15,1X,A44,1X,3(A14,1X),A29,1X,A12)') '#',&
  &'Step |','Time |','State |','Energy |','Angular |','Gradient |','Density |','Expectation Value |','Runtime |'
  write(u,'(a1,A11,1X,A14,1X,2(A7,1X),8(A14,1X),A12)') '#',&
  &'|','|','diag |','MCH |','kin |','pot |','tot |','Momentum |','RMS |','Total |','DM |','S |','|'
  write(u,'(a1,A11,1X,A14,1X,2(A7,1X),8(A14,1X),A12)') '#',&
  &'|','[fs] |','|','|','[eV] |','[eV] |','[eV] |','[hbar] |','[eV/Ang] |','|','[Debye] |','|','[sec] |'
  write(u,'(a1,a)') '#',repeat('=',175)

endsubroutine

  ! =====================================================

!> calculates and prints properties for the listing file
!> 
subroutine write_list_line(u, traj, ctrl)
  use definitions
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl
  integer :: u, imult,ims,istate,jstate,i,iatom,idir
  real*8 :: expec_dm, expec_s, grad_length, temp_dm, den
  real*8 :: p(ctrl%natom,3), r(ctrl%natom,3), summass, com(3), jmag, j(3)

  ! calculate properties
  ! gradient
  grad_length=sqrt(sum(traj%grad_ad(:,:)**2)/3/ctrl%natom)*au2eV/au2a
  ! expec_dm
  expec_dm=0.d0
  do idir=1,3
    temp_dm=0.d0
    do istate=1,ctrl%nstates
      do jstate=1,ctrl%nstates
        temp_dm=temp_dm&
        &+real( conjg(traj%U_ss(istate,traj%state_diag))&
        &*traj%DM_print_ssd(istate,jstate,idir)&
        &*traj%U_ss(jstate,traj%state_diag) )
      enddo
    enddo
    expec_dm=expec_dm+temp_dm**2
  enddo
  expec_dm=dsqrt(expec_dm)*au2debye

  ! expec_s
  expec_s=0.d0
  i=1
  do imult=1,ctrl%maxmult
    do ims=1,imult
      do istate=1,ctrl%nstates_m(imult)
        expec_s=expec_s + &
        &real(conjg(traj%U_ss(i,traj%state_diag))*traj%U_ss(i,traj%state_diag))*(imult-1)
        i=i+1
      enddo
    enddo
  enddo

  ! total density
  den=0.d0
  do istate=1,ctrl%nstates
    den=den+real(traj%coeff_diag_s(istate)*conjg(traj%coeff_diag_s(istate)))
  enddo

  ! angular momentum 
  !-momentum
  j=0.d0
  jmag=0.d0
  do iatom=1,ctrl%natom
    do idir=1,3
      p(iatom,idir)=traj%veloc_ad(iatom,idir)*traj%mass_a(iatom)
    enddo
  enddo
  !-center of mass
  com=0.d0
  summass=0.d0
  do iatom=1,ctrl%natom
    summass=summass+traj%mass_a(iatom)
    do idir=1,3
      com(idir)=com(idir)+traj%geom_ad(iatom,idir)*traj%mass_a(iatom)
    enddo
  enddo
  com=com/summass
  !-coordinates relative to center of mass
  do iatom=1,ctrl%natom
    do idir=1,3
      r(iatom,idir)=traj%geom_ad(iatom,idir)-com(idir)
    enddo
  enddo
  !-angular momentum J
  do iatom=1,ctrl%natom
    j(1)=j(1)+(r(iatom,2)*p(iatom,3)-r(iatom,3)*p(iatom,2))
    j(2)=j(2)-(r(iatom,1)*p(iatom,3)-r(iatom,3)*p(iatom,1))
    j(3)=j(3)+(r(iatom,1)*p(iatom,2)-r(iatom,2)*p(iatom,1))
  enddo
  !-norm of J
  do idir=1,3
    jmag = jmag + j(idir)**2
  enddo
  jmag = dsqrt(jmag)

  if ( (ctrl%time_uncertainty==1 .and. traj%in_time_uncertainty==0) .or. ctrl%time_uncertainty==0 ) then 
  select case (traj%kind_of_jump)
    case (0)
      continue
    case (1)
      write(u,'(A,1X,A,1X,I4,1X,A,1X,I4,1X,A,1X,F12.9)') &
      &'#','Surface Hop/Pointer State Switch: new state=',traj%state_diag,'old state=',traj%state_diag_old,'randnum=',traj%randnum
    case (2)
      write(u,'(A)') '# Jump frustrated.'
    case (3)
      write(u,'(A,1X,A,1X,I4,1X,A,1X,I4,1X,A,1X,F12.9)') &
      &'#','Surface Hop: new state=',traj%state_diag,'old state=',traj%state_diag_old,'randnum=',traj%randnum
      write(u,'(A)') '# Transition resonant to laser.'
    case (4)
      write(u,'(A)') '# Forced jump to ground state.'
      write(u,'(A,1X,A,1X,I4,1X,A,1X,I4,1X,A,1X,F12.9)') &
      &'#','Surface Hop: new state=',traj%state_diag,'old state=',traj%state_diag_old,'randnum=',traj%randnum
  endselect
  endif

  write(u,'(1X,I9,3X,F12.5,3X,2(I5,3X),8(F12.6,3X),I10)') &
  &traj%step, au2fs*traj%microtime, &
  &traj%state_diag, traj%state_MCH, &
  &traj%Ekin*au2eV, traj%Epot*au2eV, traj%Etot*au2eV, &
  &jmag, grad_length, den, expec_dm, expec_s, &
  &traj%time_step

endsubroutine

  ! =====================================================

!> writes the new geometry in xyz format (for output.xyz)
subroutine write_geom(u,traj,ctrl)
  use definitions
  use matrix
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl
  integer :: u
  integer :: iatom, idir

  write(u,'(I12)') ctrl%natom
  write(u,'(A5, 1X, F14.5, 1X, I4, 1X, i4)') 't= ',au2fs*traj%microtime, traj%state_diag, traj%state_MCH
  do iatom=1,ctrl%natom
    write(u,'(A2,3(1X,F16.9))') traj%element_a(iatom), (traj%geom_ad(iatom,idir)*au2a,idir=1,3)
  enddo

endsubroutine

  ! =====================================================

!> writes the header of the dat file. 
subroutine write_dat_initial(u, ctrl, traj)
  use definitions
  use matrix
  use string
  implicit none
  type(ctrl_type) :: ctrl
  type(trajectory_type) :: traj
  integer :: u, n
  character*8000 :: string1
  character*8000, allocatable :: string2(:)

  string1=version
  call split(string1,' ',string2,n)
  read(string2(1),*) ctrl%output_version
  deallocate(string2)
  
  if (ctrl%output_version == 1.0) then
    ! old header for SHARC v1.0
    write(u,*) ctrl%maxmult, '! maxmult'
    write(u,*) ctrl%nstates_m, '! nstates_m'
    write(u,*) ctrl%natom, '! natom'
    write(u,*) ctrl%dtstep, '! dtstep'
    write(u,*) ctrl%ezero, '! ezero'
    write(u,*) ctrl%calc_overlap, '! calc_overlap'
    write(u,*) ctrl%laser, '! laser'
    write(u,*) ctrl%nsteps,'! nsteps'
    write(u,*) ctrl%nsubsteps,'! nsubsteps'
    if (ctrl%laser==2) call vec3write(ctrl%nsteps*ctrl%nsubsteps+1, ctrl%laserfield_td, u, '! Laser field','E21.13e3')    
  elseif   (ctrl%output_version >= 2.0) then
    ! header for SHARC v2.0
    write(u,'(a14,f5.1)') 'SHARC_version ',  ctrl%output_version
    write(u,*) 'method',           ctrl%method
    write(u,*) 'integrator',       ctrl%integrator
    write(u,*) 'maxmult',          ctrl%maxmult
    write(u,'(1X,A9,64(1X,I4))') 'nstates_m',        ctrl%nstates_m
    write(u,*) 'natom',            ctrl%natom
    write(u,*) 'dtstep',           ctrl%dtstep
    write(u,*) 'nsteps',           ctrl%nsteps
    write(u,*) 'nsubsteps',        ctrl%nsubsteps
    write(u,*) 'ezero',            ctrl%ezero
    write(u,*) 'write_overlap',    ctrl%write_overlap
    write(u,*) 'write_grad',       ctrl%write_grad
    write(u,*) 'write_nacdr',      ctrl%write_NACdr
    write(u,*) 'write_property1d', ctrl%write_property1d
    write(u,*) 'write_property2d', ctrl%write_property2d
    write(u,*) 'n_property1d',     ctrl%n_property1d
    write(u,*) 'n_property2d',     ctrl%n_property2d
    write(u,*) 'laser',            ctrl%laser
    write(u,'(a)') '************************************* End of settings *************************************'
    call vecwrite(ctrl%natom,traj%atomicnumber_a,u,'! Atomic numbers','E21.13e3')
    call vecwrite(ctrl%natom,traj%element_a,     u,'! Elements',      'A3'  )
    call vecwrite(ctrl%natom,traj%mass_a,        u,'! Atomic masses', 'E21.13e3')
    if (ctrl%laser==2) call vec3write(ctrl%nsteps*ctrl%nsubsteps+1, ctrl%laserfield_td, u, '! Laser field','E21.13e3')    
    write(u,'(a)') '********************************* End of header array data ********************************'
  endif

endsubroutine

  ! =====================================================

!> writes all information for a timestep to the dat file.
subroutine write_dat(u, traj, ctrl)
  use definitions
  use matrix
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl
  integer :: u, i, j
  character(8000) :: string

  integer :: nstates, natom, stride

  ! check if writing
  stride=ctrl%output_steps_stride(1)
  if (traj%step>=ctrl%output_steps_limits(2)) then
    stride=ctrl%output_steps_stride(2)
  endif
  if (traj%step>=ctrl%output_steps_limits(3)) then
    stride=ctrl%output_steps_stride(3)
  endif
  if (modulo(traj%step,stride)==0) then

    nstates=ctrl%nstates
    natom=ctrl%natom

    if (ctrl%integrator==2) then 
      write(u,'(A)') '! 0 Step'
      write(u,'(I12)') traj%step
    else if (ctrl%integrator==1 .or. ctrl%integrator==0) then 
      write(u,'(A)') '! 0 Step'
      write(u,'(I12,F12.6,F12.6)') traj%step, traj%microtime*au2fs, ctrl%dtstep*au2fs
    endif
    call matwrite(nstates, traj%H_MCH_ss, u, '! 1 Hamiltonian (MCH) in a.u.', 'E21.13e3')
    call matwrite(nstates, traj%U_ss, u, '! 2 U matrix', 'E21.13e3')
    call matwrite(nstates, traj%DM_print_ssd(:,:,1), u, '! 3 Dipole moments X (MCH) in a.u.', 'E21.13e3')
    call matwrite(nstates, traj%DM_print_ssd(:,:,2), u, '! 3 Dipole moments Y (MCH) in a.u.', 'E21.13e3')
    call matwrite(nstates, traj%DM_print_ssd(:,:,3), u, '! 3 Dipole moments Z (MCH) in a.u.', 'E21.13e3')

    if (ctrl%write_overlap==1) then
      call matwrite(nstates, traj%overlaps_ss, u, '! 4 Overlap matrix (MCH)', 'E21.13e3')
    endif

    call vecwrite(nstates, traj%coeff_diag_s, u, '! 5 Coefficients (diag)','E21.13e3')
    if (ctrl%method==0) then !TSH
      call vecwrite(nstates, traj%hopprob_s, u, '! 6 Hopping Probabilities (diag)','E21.13e3')
    else !SCP
      if (ctrl%decoherence==11) then !DoM 
        call vecwrite(nstates, traj%switchprob_s, u, '! 6 Hopping Probabilities (diag)','E21.13e3')
      endif
    endif

    write(u,'(A)') '! 7 Ekin (a.u.)'
    write(u,'(E21.13e3)') traj%Ekin
    write(u,'(A)') '! 8 states (diag, MCH)'
    write(u,'(I12,1X,I12)') traj%state_diag, traj%state_MCH
    write(u,'(A)') '! 9 Random number'
    write(u,'(E21.13e3)') traj%randnum
    write(u,'(A)') '! 10 Runtime (sec)'
    write(u,'(I12)') traj%time_step

    call vec3write(natom, traj%geom_ad, u, '! 11 Geometry in a.u.','E21.13e3')
    call vec3write(natom, traj%veloc_ad, u, '! 12 Velocities in a.u.','E21.13e3')

    if (ctrl%output_version <= 1.0) then
      call matwrite(nstates, traj%Property2d_xss(1,:,:), u, '! 13 Property matrix (MCH)', 'E21.13e3')
    elseif (ctrl%output_version >= 2.0) then
      if (ctrl%write_property2d==1) then
  !       write(u,'(A)') '! 17 Property matrices (MCH)'
        do i=1,ctrl%n_property2d
          write(string,'(A26,I3,A3,A)') '! 13 Property matrix (MCH)',i,' : ',traj%Property2d_labels_x(i)
          call matwrite(nstates, traj%Property2d_xss(i,:,:), u, string, 'E21.13e3')
        enddo
      endif
      if (ctrl%write_property1d==1) then
  !       write(u,'(A)') '! 18 Property vectors (MCH)'
        do i=1,ctrl%n_property1d
          write(string,'(A26,I3,A3,A)') '! 14 Property vector (MCH)',i,' : ',traj%Property1d_labels_y(i)
          call vecwrite(nstates, traj%Property1d_ys(i,:), u, string, 'E21.13e3')
        enddo
      endif
    endif


    if (ctrl%write_grad == 1) then
  !     write(u,'(A)') '! 19 Gradient matrix (MCH) as x,y,z (per line) for each atom (per newline)'
      do i=1,ctrl%nstates
        write(string,'(A26,I3)') '! 15 Gradients (MCH) State',i
        call vec3write(ctrl%natom,traj%grad_mch_sad(i,:,:),u,trim(string),'E21.13e3')
      enddo
    endif


    if (ctrl%write_NACdr == 1) then
  !     write(u,'(A)') '! 20 NAC matrix (MCH) as x,y,z (per line) for each atom (per newline)'
      do i=1,ctrl%nstates
        do j=1,ctrl%nstates
          write(string,'(A31,I3,1X,I3)') '! 16 NACdr matrix element (MCH)',i,j
          call vec3write(ctrl%natom,traj%NACdr_ssad(i,j,:,:),u,trim(string),'E21.13e3')
        enddo
      enddo
    endif

  endif  ! <-- end of stride check

endsubroutine

  ! =====================================================

!> writes the total runtime of SHARC at the end and appends a SHARC fun fact.
subroutine write_final(traj)
  use definitions
  use misc
  implicit none
  type(trajectory_type) :: traj

  integer :: runtime, days, hours, minutes, seconds
  integer :: time

  runtime=time()-traj%time_start
  days=runtime/86400
  hours=mod(runtime,86400)/3600
  minutes=mod(runtime,3600)/60
  seconds=mod(runtime,60)

  write(u_log,'(A,4(1X,I4,1X,A))') 'Total wallclock time:',days,'days',hours,'h',minutes,'min',seconds,'sec'

  call write_sharcfact(u_log)

endsubroutine

  ! =====================================================

!> flushes all output files
subroutine allflush()
use definitions
implicit none

call flush(u_log)
call flush(u_lis)
call flush(u_dat)
call flush(u_geo)
call flush(u_resc)
call flush(u_rest)
call flush(u_qm_QMin)

end subroutine










endmodule
