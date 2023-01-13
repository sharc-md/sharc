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

!> # Module  pointer_basis
!> optimization of pointer basis
!>
!> \author Yinan Shu
!> \date July 3, 2022
!>
!> This module contains subroutines that allow to one to optimize 
!> pointer basis

module pointer_basis
  implicit none

  public pointer_basis_initialize
  public pointer_basis_opt

  private entropy_and_gradient_MCH

contains

! ===========================================================
! initialize  
  subroutine pointer_basis_initialize(traj, ctrl)
    use definitions
    use matrix
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl
  
    traj%entropy=dcmplx(0.d0,0.d0)
    traj%entropy_old=dcmplx(0.d0,0.d0)
    traj%linear_entropy=dcmplx(0.d0,0.d0)
    traj%linear_entropy_old=dcmplx(0.d0,0.d0)

    traj%U_pointer_old_ss=dcmplx(0.d0,0.d0)
    traj%U_pointer_ss=dcmplx(0.d0,0.d0)

  endsubroutine

! ===========================================================
! the main driver to optimize pointer basis 
  subroutine pointer_basis_opt(traj, ctrl)
    use definitions
    use matrix
    implicit none
    type(trajectory_type) :: traj
    type(ctrl_type) :: ctrl

    complex*16 :: A(ctrl%nstates,ctrl%nstates)
    complex*16 :: B(ctrl%nstates,ctrl%nstates)
    integer :: istate 


    if (printlevel>2) then
      write(u_log,*) '============================================================='
      write(u_log,*) '               Pointer Basis Optimization'
      write(u_log,*) '============================================================='
    endif

    ! update old
    traj%entropy_old=traj%entropy
    traj%linear_entropy_old=traj%linear_entropy    

    ! compute entropy and time derivative of entropy
    call entropy_and_gradient_MCH(ctrl%nstates, traj%coeff_MCH_s, traj%coeff_MCH_old_s,&
      &traj%state_MCH, traj%decotime_MCH_s, traj%decotime_MCH_old_s, traj%H_MCH_ss, traj%NACdt_ss,&
      &traj%entropy, traj%entropy_old, traj%entropy_grad,&
      &traj%linear_entropy, traj%linear_entropy_old, traj%linear_entropy_grad,&
      &A, B)
    
    if (printlevel>4) then
      write(u_log,'(A,2(F14.9,1X))') "von Neumann entropy:", traj%entropy
      write(u_log,'(A,2(F14.9,1X))') "von Neumann entropy of last timestep:", traj%entropy_old
      write(u_log,'(A,2(F14.9,1X))') "time derivative of von Neumann entropy:", traj%entropy_grad
      write(u_log,'(A,2(F14.9,1X))') "linear von Neumann entropy:", traj%linear_entropy
      write(u_log,'(A,2(F14.9,1X))') "linear von Neumann entropy of last timestep:", traj%linear_entropy_old
      write(u_log,'(A,2(F14.9,1X))') "time derivative of linear von Neumann entropy:", traj%linear_entropy_grad
    endif

    ! optimize the U based on minimizing time derivative of linear von Neumann entropy 
    traj%U_pointer_old_ss=traj%U_pointer_ss
    traj%U_pointer_ss=dcmplx(1.d0,0.d0) 
    call optimization_Upointer(ctrl%nstates, traj%coeff_MCH_s, ctrl%pointer_maxiter,&
      &A, B, traj%linear_entropy, traj%linear_entropy_grad,&
      &traj%U_pointer_old_ss, traj%U_pointer_ss)

    if (printlevel>4) then
      write(u_log,*) "pointer basis optimization finished"
      call matwrite(ctrl%nstates, traj%U_pointer_old_ss, u_log, 'rotation matrix of last time step', 'F14.9') 
      call matwrite(ctrl%nstates, traj%U_pointer_ss, u_log, 'optimized rotation matrix', 'F14.9')
    endif


  endsubroutine

! ===========================================================
! compute entropy and time derivative of entropy in MCH basis
! entropy_and_gradient_MCH(ctrl%nstates, traj%coeff_MCH_s, traj%coeff_MCH_old_s,&
!  &traj%state_MCH, traj%decotime_s, traj%decotime_old_s, traj%H_MCH_ss, traj%NACdt_ss,&
!  &traj%entropy, traj%entropy_old, traj%entropy_grad,&
!  &traj%linear_entropy, traj%linear_entropy_old, traj%linear_entropy_grad,&
!  &A, B)
! A=dendt_den_ss, B=den_dendt_ss
  subroutine entropy_and_gradient_MCH(ns, c, c_old, k, tau, tau_old, h, NACT, s, s_old, s_grad, ls, ls_old, ls_grad, A, B)
    use definitions, only: u_log, printlevel
    use matrix
    implicit none
    integer, intent(in) :: ns
    complex*16, intent(in) :: c(ns), c_old(ns)
    integer, intent(in) :: k
    real*8, intent(in) :: tau(ns), tau_old(ns)
    complex*16, intent(in) :: h(ns,ns)
    complex*16, intent(in) :: NACT(ns,ns)
    complex*16, intent(in) :: s_old, ls_old
    complex*16, intent(inout) :: s, s_grad
    complex*16, intent(inout) :: ls, ls_grad
    complex*16, intent(inout) :: A(ns,ns), B(ns,ns)

    complex*16 :: den_ss(ns,ns)
    real*8 :: den_norm_ss(ns,ns)
    complex*16 :: log_den_ss(ns,ns)
    complex*16 :: entropy_matrix(ns,ns)

    complex*16 :: DP(ns,ns)
    complex*16 :: DP_old(ns,ns)
    complex*16 :: propagator_matrix_ss(ns,ns)
    complex*16 :: coeff_diag_grad_s(ns)
    complex*16 :: ii=dcmplx(0.d0,1.d0)

    complex*16 :: den_square_ss(ns,ns)
    complex*16 :: den_square_old_ss(ns,ns)
    complex*16 :: den_grad_ss(ns,ns)

    complex*16 :: dendt_den_ss(ns,ns)
    complex*16 :: den_dendt_ss(ns,ns)

    complex*16 :: entropy, linear_entropy

    integer :: istate, jstate


    ! density matrix 
    do istate=1,ns
      do jstate=1,ns
        den_ss(istate,jstate)=c(istate)*conjg(c(jstate))
        log_den_ss(istate,jstate)=log(den_ss(istate,jstate))
        den_norm_ss(istate,jstate)=sqrt(den_ss(istate,jstate)*conjg(den_ss(istate,jstate)))
      enddo
    enddo

    call matmultiply(ns,den_ss,log_den_ss,entropy_matrix,'nn')
    call matmultiply(ns,den_ss,den_ss,den_square_ss,'nn')
 
    call matwrite(ns, entropy_matrix, u_log, 'entropy matrix', 'F14.9')
    call matwrite(ns, log_den_ss, u_log, 'log_den', 'F14.9')

    s=dcmplx(0.d0,0.d0)
    ls=dcmplx(0.d0,0.d0)

    do istate=1,ns
      s=s-entropy_matrix(istate,istate)
      ls=ls+den_square_ss(istate,istate)
    enddo
    ls=1.d0-ls

    ! compute decoherent operator (decay of mixing operator)
    DP=dcmplx(0.d0,0.d0)
    do istate=1,ns
      if (istate.ne.k) then
        DP(k,k)=DP(k,k)+(1.d0,0.d0)*tau(istate)*real(c(istate)*conjg(c(istate)))
        DP(istate,istate)=-0.5d0*(1.d0,0.d0)*tau(istate)
      endif
    enddo
    DP(k,k)=0.5d0*DP(k,k)/(real(c(k)*conjg(c(k))))

    ! compute totoal propagator 
    propagator_matrix_ss=-ii*h-NACT+DP
    ! compute time derivative of coefficients
    call matvecmultiply(ns, propagator_matrix_ss, c, coeff_diag_grad_s, 'n')

    ! compute time detivative of density matrix 
    do istate=1,ns
      do jstate=1,ns
        den_grad_ss(istate,jstate)=coeff_diag_grad_s(istate)*conjg(c(jstate))+c(istate)*conjg(coeff_diag_grad_s(jstate))
      enddo
    enddo 

    ! time derivative of von Neumann entropy 
    s_grad=dcmplx(0.d0,0.d0)
    do istate=1,ns
      do jstate=1,ns
        if (den_norm_ss(jstate,istate).ge.1.d-8) then
          s_grad=ls_grad-den_grad_ss(istate,jstate)*log(den_ss(jstate,istate))-&
            &den_ss(istate,jstate)*den_grad_ss(jstate,istate)/den_ss(jstate,istate)
        endif
      enddo
    enddo

    ! time derivative of linear von Neumann entropy 
    call matmultiply(ns,den_grad_ss,den_ss,dendt_den_ss,'nn')
    call matmultiply(ns,den_ss,den_grad_ss,den_dendt_ss,'nn')

    ls_grad=dcmplx(0.d0,0.d0)
    do istate=1,ns
      ls_grad=-dendt_den_ss(istate,istate)-den_dendt_ss(istate,istate)
    enddo

    A=dendt_den_ss
    B=den_dendt_ss

  endsubroutine

! ===========================================================
! optimization of pointer basis
! optimization_Upointer(ctrl%nstates, traj%coeff_MCH_s, traj%state_MCH,&
!   &traj%decotime_s, traj%H_MCH_ss, traj%NACdt_ss, traj%pointer_maxiter,&
!   &A, B, traj%linear_entropy, traj%linear_entropy_grad,&
!   &traj%U_pointer_old_ss, traj%U_pointer_ss)
  subroutine optimization_Upointer(ns, c, maxiter, A, B, ls, ls_grad, U_old, U)
    use definitions
    use matrix
    implicit none
    integer, intent(in) :: ns
    complex*16, intent(in) :: c(ns)
    integer, intent(in) :: maxiter
    complex*16, intent(in) :: A(ns,ns), B(ns,ns)
    complex*16, intent(in) :: ls, ls_grad
    complex*16, intent(in) :: U_old(ns,ns)
    complex*16, intent(inout) :: U(ns,ns)

    ! matrix elements of a skew-Hermitian matrix
    complex*16 :: Udiag(ns),Uoffdiag(ns*(ns-1)/2)
    complex*16 :: U_tmp(ns,ns)
    complex*16 :: ls_tmp
    complex*16 :: ls_tmp_old
    complex*16 :: ls_grad_tmp
    complex*16 :: ls_grad_tmp_old
    ! converged: 0=not converged, 1=converged
    integer :: converged
    complex*16 :: c_tmp(ns)
    complex*16 :: f
    ! gradient
    complex*16 :: df(ns*(ns+1)/2) 
    ! hessian
    complex*16 :: ddf(ns*(ns+1)/2, ns*(ns+1)/2)
    complex*16 :: U_hess(ns*(ns+1)/2,ns*(ns+1)/2)
    complex*16 :: trans_df(ns*(ns+1)/2)
    complex*16 :: x(ns*(ns+1)/2)
    complex*16 :: trans_x(ns*(ns+1)/2)
    ! rotation matrix 
    complex*16 :: Uexp(ns,ns)

    integer :: istate, jstate, pstate
    integer :: iter
    integer :: ipair


    ! initialize diagonal and off-diagonal matrix elements. 
    do istate=1,ns
      Udiag(istate)=U_old(istate,istate)
      do jstate=istate+1,ns
        pstate=(istate-1)*(2*ns-istate)/2
        Uoffdiag(pstate+jstate)=U_old(istate,jstate)
      enddo
    enddo

    ! construction of rotation matrix
    do istate=1,ns
      U_tmp(istate,istate)=Udiag(istate)
      do jstate=istate+1,ns
        pstate=(istate-1)*(2*ns-istate)/2
        U_tmp(istate,jstate)=Uoffdiag(pstate+jstate)
        U_tmp(jstate,istate)=-conjg(U_tmp(istate,jstate))
      enddo
    enddo

    call linear_entropy_and_grad(ns, U_tmp, c, A, B, ls_tmp, ls_grad_tmp, f)

    converged=0
    iter=0
    

    ! optimization
    do while (converged==0)

      iter=iter+1
      if (printlevel>4) write(u_log,*) "iteration:", iter    
      ls_tmp_old=ls_tmp
      ls_grad_tmp_old=ls_grad_tmp

      ! compute liner entropy, time derivative of linear entropy, and 
      ! derivative of time derivative of linear entropy w.r.t. rotation matrix
      call compute_numerical_gradients(ns, U_tmp, c, A, B, f, df)
      write(u_log,*) "numerical_gradients done"
      call compute_numerical_hessian(ns, U_tmp, c, A, B, f, df, ddf)
      write(u_log,*) "numerical_hessian done"

      Uexp=U_tmp
      call exponentiate(ns,Uexp,(1.d0,0.d0))

      if (printlevel>4) then
        write(u_log,'(A,2(F14.9,1X))') "linear entropy:", ls_tmp
        write(u_log,'(A,2(F14.9,1X))') "time derivative of linear entropy (tdle):", ls_grad_tmp
        call matwrite(ns, U_tmp, u_log, 'skew-hermitian matrix', 'F14.9')
        call matwrite(ns, Uexp, u_log, 'rotation matrix', 'F14.9')
        call vecwrite(ns, df, u_log, 'gradient of tdle w.r.t. rotation matrix', 'F14.9')
        call matwrite(ns*(ns+1)/2, ddf, u_log, 'hessian of tdle w.r.t. rotation matrix', 'F14.9')
      endif
 
      call diagonalize(ns*(ns+1)/2,ddf,U_hess)
      call matvecmultiply(ns, U_hess, df, trans_df, 't')
      do ipair=1, ns*(ns+1)/2
        trans_x(ipair)=trans_df(ipair)/ddf(ipair,ipair)
      enddo
      call matvecmultiply(ns, U_hess, trans_x, x, 'n')

      call matwrite(ns*(ns+1)/2, U_hess, u_log, 'U_hess', 'F14.9')

      ! assign x to U_tmp
      do istate=1,ns
        U_tmp(istate,istate)=x(istate)
      enddo
      do istate=1,ns
        do jstate=istate+1, ns
          pstate=(istate-1)*(2*ns-istate)/2
          U_tmp(istate,jstate)=x(ns+pstate+jstate)
          U_tmp(jstate,istate)=-conjg(U_tmp(istate,jstate))
        enddo
      enddo   

      call linear_entropy_and_grad(ns, U_tmp, c, A, B, ls_tmp, ls_grad_tmp, f)
      Uexp=U_tmp
      call exponentiate(ns,Uexp,(1.d0,0.d0))

      if (printlevel>4) then
        write(u_log,'(A,2(F14.9,1X))') "updated linear entropy:", ls_tmp
        call matwrite(ns, Uexp, u_log, 'updated rotation matrix', 'F14.9')
      endif

      if (abs(ls_grad_tmp-ls_grad_tmp_old).lt.1.d-10 .and. iter.le.maxiter) then 
        converged=1 
        U=U_tmp
      else if (iter.gt.maxiter) then 
        write(u_log,*) "optimization reaches maximum iteration"
        converged=1
      else 
        converged=0
      endif

    enddo

  endsubroutine

! ===========================================================
! compute liner entropy, time derivative of linear entropy 
  subroutine linear_entropy_and_grad(ns, U_tmp, c, A, B, ls_tmp, ls_grad_tmp, f)
    use definitions
    use matrix
    implicit none
    integer, intent(in) :: ns
    complex*16, intent(in) :: U_tmp(ns,ns)
    complex*16, intent(in) :: c(ns)
    complex*16, intent(in) :: A(ns,ns), B(ns,ns)
    complex*16, intent(inout) :: ls_tmp, ls_grad_tmp
    complex*16, intent(inout) :: f

    complex*16 :: Uexp(ns,ns)
    complex*16 :: c_tmp(ns)

    complex*16 :: den_ss(ns,ns)
    complex*16 :: den_square_ss(ns,ns)
    complex*16 :: A_tmp(ns,ns), A_tmp1(ns,ns)
    complex*16 :: B_tmp(ns,ns), B_tmp1(ns,ns)

    complex*16 :: dfdUexp(ns,ns)
    complex*16 :: transA_tmp(ns,ns), transB_tmp(ns,ns)
    complex*16 :: m1(ns,ns),m2(ns,ns),m3(ns,ns),m4(ns,ns)

    integer :: istate, jstate

    ! obtain rotation matrix
    Uexp=U_tmp
    call exponentiate(ns,Uexp,(1.d0,0.d0))

    ! compute coefficients
    call matvecmultiply(ns, Uexp, c, c_tmp, 'n')

    ! density matrix 
    do istate=1,ns
      do jstate=1,ns
        den_ss(istate,jstate)=c_tmp(istate)*conjg(c_tmp(jstate))
      enddo
    enddo
    call matmultiply(ns,den_ss,den_ss,den_square_ss,'nn')

    ! compute linear entropy 
    ls_tmp=dcmplx(0.d0,0.d0)
    do istate=1,ns
      ls_tmp=ls_tmp+den_square_ss(istate,istate)
    enddo
    ls_tmp=1.d0-ls_tmp
 
    ! compute time derivative of linear entropy 
    A_tmp=A
    B_tmp=B
    call transform(ns, A_tmp, Uexp, 'atau')
    call transform(ns, B_tmp, Uexp, 'atau')

    ls_grad_tmp=dcmplx(0.d0,0.d0)
    do istate=1,ns
      ls_grad_tmp=-A_tmp(istate,istate)-B_tmp(istate,istate)
    enddo

    f=ls_grad_tmp

  endsubroutine

! ===========================================================
! compute numerical gradient of time derivative of linear entropy 
! w.r.t. the rotation matrix 
! compute_numerical_gradients(ns, U_tmp, c, A, B, df)
  subroutine compute_numerical_gradients(ns, U_tmp, c, A, B, f, df)
    use definitions
    use matrix
    implicit none
    integer, intent(in) :: ns
    complex*16, intent(in) :: U_tmp(ns,ns)
    complex*16, intent(in) :: c(ns)
    complex*16, intent(in) :: A(ns,ns), B(ns,ns)
    complex*16, intent(in) :: f
    complex*16, intent(inout) :: df(ns*(ns+1)/2)

    complex*16 :: U_tmp_tmp(ns,ns)
    complex*16 :: f_forward_real, f_forward_imaginary
    complex*16 :: f_backward_real, f_backward_imaginary
    complex*16 :: ls_tmp
    complex*16 :: ls_grad_tmp

    integer :: istate, jstate, pstate

    ! compute numerical gradient of diagonal elements
    ! notice diagonal elements for skew-Hermitian matrix is always pure imarginary
    do istate=1,ns
      U_tmp_tmp=U_tmp
      U_tmp_tmp(istate,istate)=U_tmp(istate,istate)+(0.d0, 1.d-4)
      call linear_entropy_and_grad(ns, U_tmp_tmp, c, A, B, ls_tmp, ls_grad_tmp, f_forward_imaginary)
      U_tmp_tmp=U_tmp
      U_tmp_tmp(istate,istate)=U_tmp(istate,istate)-(0.d0, 1.d-4)
      call linear_entropy_and_grad(ns, U_tmp_tmp, c, A, B, ls_tmp, ls_grad_tmp, f_backward_imaginary)
      df(istate)=(((f_forward_imaginary-f)/(1.d-4))+(((f_backward_imaginary-f)/(1.d-4))))/2*(0.d0, 1.d0)
    enddo

    ! compute numerical gradient of off-diagonal elements
    do istate=1,ns
      do jstate=istate+1,ns
        U_tmp_tmp=U_tmp
        U_tmp_tmp(istate,jstate)=U_tmp(istate,jstate)+(1.d-4, 0.d0)
        U_tmp_tmp(jstate,istate)=-conjg(U_tmp_tmp(istate,jstate))
        call linear_entropy_and_grad(ns, U_tmp_tmp, c, A, B, ls_tmp, ls_grad_tmp, f_forward_real)
        U_tmp_tmp=U_tmp
        U_tmp_tmp(istate,jstate)=U_tmp(istate,jstate)-(1.d-4, 0.d0)
        U_tmp_tmp(jstate,istate)=-conjg(U_tmp_tmp(istate,jstate))
        call linear_entropy_and_grad(ns, U_tmp_tmp, c, A, B, ls_tmp, ls_grad_tmp, f_backward_real)
        U_tmp_tmp=U_tmp
        U_tmp_tmp(istate,jstate)=U_tmp(istate,jstate)+(0.d0, 1.d-4)
        U_tmp_tmp(jstate,istate)=-conjg(U_tmp_tmp(istate,jstate))
        call linear_entropy_and_grad(ns, U_tmp_tmp, c, A, B, ls_tmp, ls_grad_tmp, f_forward_imaginary)
        U_tmp_tmp=U_tmp
        U_tmp_tmp(istate,jstate)=U_tmp(istate,jstate)-(0.d0, 1.d-4)
        U_tmp_tmp(jstate,istate)=-conjg(U_tmp_tmp(istate,jstate))
        call linear_entropy_and_grad(ns, U_tmp_tmp, c, A, B, ls_tmp, ls_grad_tmp, f_backward_imaginary)
        pstate=(istate-1)*(2*ns-istate)/2
        df(ns+pstate+jstate)=&
         &(((f_forward_real-f)/(1.d-4))+(((f_backward_real-f)/(1.d-4))))/2*(1.d0,0.d0) +&
         &(((f_forward_imaginary-f)/(1.d-4))+(((f_backward_imaginary-f)/(1.d-4))))/2*(0.d0, 1.d0)
      enddo
    enddo

  endsubroutine

! ===========================================================
! compute numerical hessian of time derivative of linear entropy 
! w.r.t. the rotation matrix 
! compute_numerical_hessian(ns, U_tmp, c, A, B, f, df, ddf)
  subroutine compute_numerical_hessian(ns, U_tmp, c, A, B, f, df, ddf)
    use definitions
    use matrix
    implicit none
    integer, intent(in) :: ns
    complex*16, intent(in) :: U_tmp(ns,ns)
    complex*16, intent(in) :: c(ns)
    complex*16, intent(in) :: A(ns,ns), B(ns,ns)
    complex*16, intent(in) :: f
    complex*16, intent(in) :: df(ns*(ns+1)/2)
    complex*16, intent(inout) :: ddf(ns*(ns+1)/2, ns*(ns+1)/2)

    complex*16 :: U_tmp_tmp(ns,ns)
    complex*16 :: df_forward_real(ns*(ns+1)/2), df_forward_imaginary(ns*(ns+1)/2)
    complex*16 :: df_backward_real(ns*(ns+1)/2), df_backward_imaginary(ns*(ns+1)/2)
    complex*16 :: ls_tmp
    complex*16 :: ls_grad_tmp

    integer :: istate, jstate, pstate

    write(u_log,*) "in nh"

    ! compute numerical hessian of diagonal elements
    do istate=1,ns
      U_tmp_tmp=U_tmp
      U_tmp_tmp(istate,istate)=U_tmp(istate,istate)+(0.d0, 1.d-4)
      call compute_numerical_gradients(ns, U_tmp_tmp, c, A, B, f, df_forward_imaginary)
      U_tmp_tmp=U_tmp
      U_tmp_tmp(istate,istate)=U_tmp(istate,istate)-(0.d0, 1.d-4)
      call compute_numerical_gradients(ns, U_tmp_tmp, c, A, B, f, df_backward_imaginary)
      ddf(istate,:)=(((df_forward_imaginary-df)/(1.d-4))+(((df_backward_imaginary-df)/(1.d-4))))/2*(0.d0, 1.d0)
    enddo

    ! compute numerical hessian of off-diagonal elemenets
    do istate=1,ns
      do jstate=istate+1,ns
        U_tmp_tmp=U_tmp
        U_tmp_tmp(istate,jstate)=U_tmp(istate,jstate)+(1.d-4, 0.d0)
        U_tmp_tmp(jstate,istate)=-conjg(U_tmp_tmp(istate,jstate))
        call compute_numerical_gradients(ns, U_tmp_tmp, c, A, B, f, df_forward_real)
        U_tmp_tmp=U_tmp
        U_tmp_tmp(istate,jstate)=U_tmp(istate,jstate)-(1.d-4, 0.d0)
        U_tmp_tmp(jstate,istate)=-conjg(U_tmp_tmp(istate,jstate))
        call compute_numerical_gradients(ns, U_tmp_tmp, c, A, B, f, df_backward_real)
        U_tmp_tmp=U_tmp
        U_tmp_tmp(istate,jstate)=U_tmp(istate,jstate)+(0.d0, 1.d-4)
        U_tmp_tmp(jstate,istate)=-conjg(U_tmp_tmp(istate,jstate))
        call compute_numerical_gradients(ns, U_tmp_tmp, c, A, B, f, df_forward_imaginary)
        U_tmp_tmp=U_tmp
        U_tmp_tmp(istate,jstate)=U_tmp(istate,jstate)-(0.d0, 1.d-4)
        U_tmp_tmp(jstate,istate)=-conjg(U_tmp_tmp(istate,jstate))
        call compute_numerical_gradients(ns, U_tmp_tmp, c, A, B, f, df_backward_imaginary)
        pstate=(istate-1)*(2*ns-istate)/2
        ddf(ns+pstate+jstate,:)=&
         &(((df_forward_real-df)/(1.d-4))+(((df_backward_real-df)/(1.d-4))))/2*(1.d0,0.d0) +&
         &(((df_forward_imaginary-df)/(1.d-4))+(((df_backward_imaginary-df)/(1.d-4))))/2*(0.d0, 1.d0)
      enddo
    enddo

  endsubroutine


endmodule
