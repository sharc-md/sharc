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

!> # Module BSH
!> 
!> \author Yinan Shu
!> \date 1.1.2020
!>
!> BSH module contains subroutines used by Bulirsch-Stoer-Hack integrator
!>
module bsh
 contains

! ===========================================================
!> the main driver of Bulirsch-Stoer-Hack integrator
subroutine BSH_propagation(traj,ctrl,bsv)
  use definitions
  use matrix
  implicit none

  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

  integer, intent(in) :: bsv

  real*8 :: y(bsv),dydt(bsv)
  real*8 :: dtstep_next
  integer :: printlevel_save

  integer :: i

  traj%preprob_old_s3=traj%preprob_s3

  if (printlevel>1) then
    write(u_log,*) '============================================================='
    write(u_log,*) '  Start Bulirsch-Stoer-Hack integrator intermediate steps'
    write(u_log,*) '============================================================='
    write(u_log,*) ' '
    write(u_log,'(A,1X,F14.9,1X,A)') 'target step:', ctrl%dtstep*au2fs, 'fs'
  endif

  printlevel_save=printlevel 
  printlevel=0

! Create an BSH array
  call BSHarray_create(traj,ctrl,y,dydt,bsv)

  call BSHstep(traj,ctrl,y,dydt,bsv,traj%microtime,ctrl%dtstep,ctrl%convthre,dtstep_next)
  ctrl%dtstep=dtstep_next

! convert the BSH array back to variables
  call BSHarray_convertback(traj,ctrl,y,bsv)

  printlevel=printlevel_save

  if (printlevel>1) then
    write(u_log,*) '============================================================='
    write(u_log,*) '  End Bulirsch-Stoer-Hack integrator intermediate steps'
    write(u_log,*) '============================================================='
  endif

  if (printlevel>4) then
    write(u_log,*) '============================================================='
    write(u_log,*) '              Current geometry and velocity'
    write(u_log,*) '============================================================='
    call vec3write(ctrl%natom,traj%veloc_ad,u_log,'velocity','F12.9')
    call vec3write(ctrl%natom,traj%geom_ad,u_log,'geometry','F12.9')
  endif


endsubroutine

!==========================================================
!> The Bulirsch-Stoer-Hack integrator integrate a single array, called y,
!> and useing its derivatives, called dydt, to integrate
!> subroutine BSHarray_create creates the array, y and dydt
subroutine BSHarray_create(traj,ctrl,y,dydt,bsv)
  use definitions
  use matrix
  implicit none

  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

  integer, intent(in) :: bsv
  real*8, intent(out) :: y(bsv), dydt(bsv)

  integer :: istate, jstate, iatom, idir
  integer :: i, track

! The y array contains elements:
! 1 to ctrl%natom                                            -> X coordinates of all atoms traj%geom_ad(iatom,1)
! ctrl%natom+1 to 2*ctrl%natom                               -> Y coordinates of all atoms traj%geom_ad(iatom,2)
! 2*ctrl%natom+1 to 3*ctrl%natom                             -> Z coordinates of all atoms traj%geom_ad(iatom,3)
! 3*ctrl%natom+1 to 4*ctrl%natom                             -> X momenta of all atoms traj%mass_a(iatom)*traj%veloc_ad(iatom,1)
! 4*ctrl%natom+1 to 5*ctrl%natom                             -> Y momenta of all atoms traj%mass_a(iatom)*traj%veloc_ad(iatom,2)
! 5*ctrl%natom+1 to 6*ctrl%natom                             -> Z momenta of all atoms traj%mass_a(iatom)*traj%veloc_ad(iatom,3) 
! 6*ctrl%natom+1 to 6*trl%natom+ctrl%nstates                 -> Real part of electronic coefficients real(traj%coeff_MCH_s(istate))
! 6*trl%natom+ctrl%nstates+1 to 6*trl%natom+2*ctrl%nstates   -> Imaginary part of electronic coefficients aimag(traj%coeff_MCH_s(istate)) 
! 6*trl%natom+2*ctrl%nstates+1 to 6*trl%natom+3*ctrl%nstates -> Real part of coherent electronic coefficients real(traj%ccoeff_MCH_s(istate)) - used by CSDM
! 6*trl%natom+3*ctrl%nstates+1 to 6*trl%natom+4*ctrl%nstates -> Imaginary part of coherent electronic coefficients aimag(traj%ccoeff_MCH_s(istate)) - used by CSDM
! 6*trl%natom+4*ctrl%nstates+1 to 6*trl%natom+5*ctrl%nstates -> electronic phase
! 6*trl%natom+5*ctrl%nstates+1 to 6*trl%natom+6*ctrl%nstates -> b_ik
! 6*trl%natom+6*ctrl%nstates+1 to 6*trl%natom+7*ctrl%nstates -> b_ik(+)
! 6*trl%natom+7*ctrl%nstates+1 to 6*trl%natom+8*ctrl%nstates -> b_ik(-)

  ! coordinates
  do iatom=1,ctrl%natom
    do idir=1,3
      i=ctrl%natom*(idir-1)+iatom
      y(i)=traj%geom_ad(iatom,idir)
      dydt(i)=traj%veloc_ad(iatom,idir)
    enddo
  enddo
  track=3*ctrl%natom

  ! momenta
  do iatom=1,ctrl%natom
    do idir=1,3
      i=track+ctrl%natom*(idir-1)+iatom
      y(i)=traj%mass_a(iatom)*traj%veloc_ad(iatom,idir)
      dydt(i)=-traj%grad_ad(iatom,idir)
    enddo
  enddo
  track=track+3*ctrl%natom

  ! Real traj%coeff_MCH_s(istate)
  do istate=1,ctrl%nstates
    i=track+istate
    y(i)=real(traj%coeff_MCH_s(istate))
    dydt(i)=traj%gRcoeff_MCH_s(istate)
  enddo
  track=track+ctrl%nstates

  ! Imaginary traj%coeff_MCH_s(istate)
  do istate=1,ctrl%nstates
    i=track+istate
    y(i)=aimag(traj%coeff_MCH_s(istate))
    dydt(i)=traj%gIcoeff_MCH_s(istate)
  enddo
  track=track+ctrl%nstates

  ! Real traj%ccoeff_MCH_s(istate)
  do istate=1,ctrl%nstates
    i=track+istate
    y(i)=real(traj%ccoeff_MCH_s(istate))
    dydt(i)=traj%gRccoeff_MCH_s(istate)
  enddo
  track=track+ctrl%nstates

  ! Imaginary traj%ccoeff_MCH_s(istate)
  do istate=1,ctrl%nstates
    i=track+istate
    y(i)=aimag(traj%ccoeff_MCH_s(istate))
    dydt(i)=traj%gIccoeff_MCH_s(istate)
  enddo
  track=track+ctrl%nstates

  ! Electronic phase
  do istate=1,ctrl%nstates
    i=track+istate
    y(i)=traj%ephase_s(istate)
    dydt(i)=traj%gephase_s(istate)
  enddo
  track=track+ctrl%nstates

  ! b_ik, b_ik(+), b_ik(-)
  do istate=1,ctrl%nstates
    i=track+istate
    y(i)=0
    dydt(i)=traj%gpreprob_s3(istate,1)
  enddo
  track=track+ctrl%nstates

  do istate=1,ctrl%nstates
    i=track+istate
    y(i)=0
    dydt(i)=traj%gpreprob_s3(istate,2)
  enddo
  track=track+ctrl%nstates

  do istate=1,ctrl%nstates
    i=track+istate
    y(i)=0
    dydt(i)=traj%gpreprob_s3(istate,3)
  enddo
  track=track+ctrl%nstates


endsubroutine


!==========================================================
!> The Bulirsch-Stoer-Hack integrator integrate a single array, called y,
!> and useing its derivatives, called dydt, to integrate
!> subroutine BSHarray_convertback reads y and put these values to traj, ctrl
subroutine BSHarray_convertback(traj,ctrl,y,bsv)
  use definitions
  use matrix
  implicit none

  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

  integer, intent(in) :: bsv
  real*8, intent(in) :: y(bsv)

  integer :: istate, jstate, iatom, idir
  integer :: i, track

  ! coordinates
  do iatom=1,ctrl%natom
    do idir=1,3
      i=ctrl%natom*(idir-1)+iatom
      traj%geom_ad(iatom,idir)=y(i)
    enddo
  enddo
  track=3*ctrl%natom

  ! momenta
  do iatom=1,ctrl%natom
    do idir=1,3
      i=track+ctrl%natom*(idir-1)+iatom
      traj%veloc_ad(iatom,idir)=y(i)/traj%mass_a(iatom)
    enddo
  enddo
  track=track+3*ctrl%natom

  ! Real traj%coeff_MCH_s(istate)
  do istate=1,ctrl%nstates
    i=track+istate
    traj%coeff_MCH_s(istate)=y(i)
  enddo
  track=track+ctrl%nstates

  ! Imaginary traj%coeff_MCH_s(istate)
  do istate=1,ctrl%nstates
    i=track+istate
    traj%coeff_MCH_s(istate)=traj%coeff_MCH_s(istate)+ii*(y(i))
  enddo
  track=track+ctrl%nstates

  ! Real traj%ccoeff_MCH_s(istate)
  do istate=1,ctrl%nstates
    i=track+istate
    traj%ccoeff_MCH_s(istate)=y(i)
  enddo
  track=track+ctrl%nstates

  ! Imaginary traj%ccoeff_MCH_s(istate)
  do istate=1,ctrl%nstates
    i=track+istate
    traj%ccoeff_MCH_s(istate)=traj%ccoeff_MCH_s(istate)+ii*(y(i))
  enddo
  track=track+ctrl%nstates

  ! Electronic phase
  do istate=1,ctrl%nstates
    i=track+istate
    traj%ephase_s(istate)=y(i)
  enddo
  track=track+ctrl%nstates

  ! b_ik, b_ik(+), b_ik(-)
  do istate=1,ctrl%nstates
    i=track+istate
    traj%preprob_s3(istate,1)=y(i)
  enddo
  track=track+ctrl%nstates

  do istate=1,ctrl%nstates
    i=track+istate
    traj%preprob_s3(istate,2)=y(i)
  enddo
  track=track+ctrl%nstates

  do istate=1,ctrl%nstates
    i=track+istate
    traj%preprob_s3(istate,3)=y(i)
  enddo
  track=track+ctrl%nstates

endsubroutine

!==========================================================
!> This is the Bulirsch-Stoer integrator taken from Numerical Recipies
!> BSHstep(traj,ctrl,y,dydt,bsv,traj%microtime,ctrl%dtstep,ctrl%convthre,dtstep_next)
subroutine BSHstep(traj,ctrl,y,dydx,nv,x,htry,eps,hnext)
  use definitions

  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

  real*8, intent(inout) :: y(nv), dydx(nv)
  integer, intent(in) :: nv
  real*8, intent(inout) :: x
  real*8, intent(in) :: htry, eps
  real*8, intent(out) :: hnext

  integer, parameter :: kmaxx=8, imax=kmaxx+1
  real*8, parameter :: minstep=0.001/au2fs
  real*8, parameter :: safe1=0.25, safe2=0.7, redmax=1.d-5, redmin=0.7,&
                       tiny=1.d-30, scalmx=0.1

  integer :: i,iq,k,kk,km,kmax,kopt,nseq(imax),nsurf
  real*8 :: eps1, epsold, errmax, fact, h, red, scale, work, wrkmin, &
            xest, xnew, a(imax), alf(kmaxx,kmaxx), err(kmaxx), hdid, &
            yscal(nv)
  real*8 :: yerr(nv), ysav(nv), yseq(nv)
  logical :: first, reduct
  save a,alf,epsold,first,kmax,kopt,nseq,xnew
  data first/.true./,epsold/-1.d0/
  data nseq /2,4,6,8,10,12,14,16,18/ 

  yscal=1.d0

      if(eps.ne.epsold)then
        hnext=-1.d29
        xnew=-1.d29
        eps1=safe1*eps
        a(1)=nseq(1)+1
        do 11 k=1,kmaxx
          a(k+1)=a(k)+nseq(k+1)
11      continue
        do 13 iq=2,kmaxx
          do 12 k=1,iq-1
            alf(k,iq)=eps1**((a(k+1)-a(iq+1))/((a(iq+1)-a(1)+1.)*(2*k+1)))
12        continue
13      continue
        epsold=eps
        do 14 kopt=2,kmaxx-1
          if(a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt))goto 1
14      continue
1       kmax=kopt
      endif

      h=htry
      do 15 i=1,nv
        ysav(i)=y(i)
15    continue
      if(h.ne.hnext.or.x.ne.xnew)then
        first=.true.
        kopt=kmax
      endif
      reduct=.false.
2     do 17 k=1,kmax
        xnew=x+h
        if(xnew.eq.x) stop 'step size underflow in bsstep'

        write(u_log,'(A,1X,I2)') "BSH_midpoint:",nseq(k) 
        call BSH_midpoint(traj,ctrl,ysav,dydx,nv,x,h,nseq(k),yseq)
        xest=(h/nseq(k))**2
        call BSH_extrapolation(k,xest,yseq,y,yerr,nv)
        if(k.ne.1)then
          errmax=TINY
          do 16 i=1,nv
            errmax=max(errmax,abs(yerr(i)/yscal(i)))
16        continue
          errmax=errmax/eps
          km=k-1
          err(km)=(errmax/SAFE1)**(1./(2*km+1))
        endif
        if(k.ne.1.and.(k.ge.kopt-1.or.first))then
          if(errmax.lt.1.)goto 4
          if(k.eq.kmax.or.k.eq.kopt+1)then
            red=SAFE2/err(km)
            goto 3
          else if(k.eq.kopt)then
            if(alf(kopt-1,kopt).lt.err(km))then
              red=1./err(km)
              goto 3
            endif
          else if(kopt.eq.kmax)then
            if(alf(km,kmax-1).lt.err(km))then
              red=alf(km,kmax-1)*SAFE2/err(km)
              goto 3
            endif
          else if(alf(km,kopt).lt.err(km))then
            red=alf(km,kopt-1)/err(km)
            goto 3
          endif
        endif
17    continue
3     red=min(red,REDMIN)
      red=max(red,REDMAX)
      h=h*red
      reduct=.true.
      goto 2
4     x=xnew
      hdid=h
      first=.false.
      wrkmin=1.e35
      do 18 kk=1,km
        fact=max(err(kk),SCALMX)
        work=fact*a(kk+1)
        if(work.lt.wrkmin)then
          scale=fact
          wrkmin=work
          kopt=kk+1
        endif
18    continue
      hnext=h/scale
      if(kopt.ge.k.and.kopt.ne.kmax.and..not.reduct)then
        fact=max(scale/alf(kopt-1,kopt),SCALMX)
        if(a(kopt+1)*fact.le.wrkmin)then
          hnext=h/fact
          kopt=kopt+1
        endif
      endif 

endsubroutine

!==========================================================
!> BSH_midpoint is midpoint method used by Bulirsch-Stoer-Hack integrator
!> call BSH_midpoint(traj,ctrl,ysav,dydx,nv,x,h,nseq(k),yseq)
subroutine BSH_midpoint(traj,ctrl,y,dydt,nv,time,htot,nstep,yout)
  use definitions
  implicit none
  
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

  real*8, intent(in) :: y(nv),dydt(nv)
  integer, intent(in) :: nv, nstep
  real*8, intent(in) :: time, htot
  real*8, intent(out) :: yout(nv)

  integer :: i, j, k
  real*8 :: h,time2,h2
  real*8 :: y1(nv), y2(nv), swap
  real*8 :: dy2dt(nv)

  integer :: ix

  h=htot/dble(nstep)

  do i=1,nv
    y1(i)=y(i)
    y2(i)=y(i)+h*dydt(i)
  enddo 

  time2=time+h

  call Get_all_gradients(traj,ctrl,y2,time2,dy2dt,nv)

  h2=2.d0*h
  do j=2,nstep
    do i=1,nv
      swap=y1(i)+h2*dy2dt(i)
      y1(i)=y2(i)
      y2(i)=swap
    enddo
    time2=time2+h
    call Get_all_gradients(traj,ctrl,y2,time2,dy2dt,nv)
  enddo

  do i=1,nv
    yout(i)=0.5*(y1(i)+y2(i)+h*dy2dt(i))
  enddo
 
endsubroutine

!==========================================================
subroutine BSH_extrapolation(iest,xest,yest,yz,dy,nv)
  implicit none
  integer :: iest,nv
  real*8 :: xest,dy(nv),yest(nv),yz(nv)
  integer, parameter :: imax=13,nmax=10000
  integer :: j, k1
  real*8 :: delta, f1,f2,q,d(nv),qcol(nmax,imax),x(imax) 
  save qcol, x

      if (NMAX.lt.nv) then
          print *,"increase NMAX to at least ",nv," in PZEXTR"
          stop
      endif

      x(iest)=xest
      do 11 j=1,nv
        dy(j)=yest(j)
        yz(j)=yest(j)
11    continue
      if(iest.eq.1) then
        do 12 j=1,nv
          qcol(j,1)=yest(j)
12      continue
      else
        do 13 j=1,nv
          d(j)=yest(j)
13      continue
        do 15 k1=1,iest-1
          delta=1./(x(iest-k1)-xest)
          f1=xest*delta
          f2=x(iest-k1)*delta
          do 14 j=1,nv
            q=qcol(j,k1)
            qcol(j,k1)=dy(j)
            delta=d(j)-q
            dy(j)=f1*delta
            d(j)=f2*delta
            yz(j)=yz(j)+dy(j)
14        continue
15      continue
        do 16 j=1,nv
          qcol(j,iest)=dy(j)
16      continue
      endif

endsubroutine

!==========================================================
!> This subroutine computes all gradients
subroutine Get_all_gradients(traj,ctrl,y,time,dydt,nv)
  use decoherence_dom
  use definitions
  use electronic
  use qm
  use matrix

  implicit none

  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl

  real*8, intent(inout) :: y(nv)
  integer, intent(in) :: nv
  real*8, intent(in) :: time
  real*8, intent(out) :: dydt(nv)

  call BSHarray_convertback(traj,ctrl,y,nv)
 
  if (printlevel>4) then
    call vecwrite(ctrl%nstates, traj%coeff_MCH_s, u_log, 'Coefficient coeff_mch_s','F14.9')
    call vecwrite(ctrl%nstates, traj%ccoeff_MCH_s, u_log, 'Coherent coefficient ccoeff_mch_s','F14.9')
    call vec3write(ctrl%nstates, traj%preprob_s3, u_log, 'electronic phase','F14.9')
  endif
  
  ! QM Calculation
  call do_qm_calculations(traj,ctrl)
  call Adjust_phases(traj,ctrl)
  call Nac_diagonal(traj, ctrl)
  call Calculate_cDIAG(traj,ctrl)
  ! Compute electronic and nuclear gradients
  call Electronic_gradients_MCH(traj,ctrl)
  if (ctrl%method==1 .and. ctrl%decoherence==11) then
    call DoM_stepBS(traj,ctrl)
  endif
  call Mix_gradients(traj,ctrl)

  call BSHarray_create(traj,ctrl,y,dydt,nv)

endsubroutine

endmodule
