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

!> # Module MATRIX
!>
!> \author Sebastian Mai
!> \date 10.07.2014
!>
!> Contains matrix routines used in SHARC.
!> Matrix multiplication and diagonalisation use Lapack routines
!> 
!> Contains:
!> - lapack workspace allocation
!> - overloaded (for real*8 and complex*16) routines for:
!>   - lapack matrix multiplication (A^c.B^c and U^c.A.U^c)
!>   - lapack matrix diagonalisation
!>   - Loewdin orthogonalisation
!>   - matrix normalisation
!>   - diagonalisation with projection of U matrix onto Uold
!>   - read and write routines for matrices
!>   - check functions for symmetry/hermiticity and anti-symmetry/anti-hermiticity
!>   - check functions for orthogonality/unitarity
!> 
!> the lapack workspaces are private arrays, global to this module
!> they can only be changed by calling allocatelapack()
module matrix

implicit none

! =================================================================== !

private

! lapack workspace variables
complex*16, save, allocatable :: lapack_work_z(:)               !< WORK in zheev
real*8,     save, allocatable :: lapack_rwork_z(:)              !< RWORK in zheev
integer,    save              :: lapack_lwork_z                 !< LWORK in zheev
real*8,     save, allocatable :: lapack_work_d(:)               !< WORK in dsyev
integer,    save              :: lapack_lwork_d                 !< LWORK in dsyev
logical,    save              :: lapack_isalloc=.false.         !< whether lapack was already allocated

!> this threshold is used in the diagonalisation with subsequent projection on the Uold matrix
!> can be changed by calling set_project_diffthr()
real*8, save :: diagonalize_degeneracy_diff=1.d-12

! =================================================================== !

! this routines can't be called from outside

private dnormalize,znormalize
private dlowdin,zlowdin
private dtransform,ztransform
private dmultiply,zmultiply
private ddiagonalize,zdiagonalize
private dexponential,zexponential
private ddiagonalize_and_project,zdiagonalize_and_project
private dwrite,zwrite,iwrite
private dvecmultiply,zvecmultiply
private dvecwrite,zvecwrite
private d3vecwrite,z3vecwrite
private dishermitian,zishermitian
private disantihermitian,zisantihermitian
private disunitary,zisunitary
private z_project_recursive
private d3project_a_on_b, z3project_a_on_b
private zintruder

! this routines can be called from outside

public normalize
public lowdin
public transform
public matmultiply
public matvecmultiply
public diagonalize
public exponentiate
public diagonalize_and_project
public allocate_lapack
public set_project_diffthr
public matwrite
public vecwrite
public vec3write
public matread
public vecread
public vec3read
public ishermitian
public isantihermitian
public isunitary
public project_recursive
public project_a_on_b
public intruder

! =================================================================== !

! Overloading all procedures for use with real*8 and complex*16 matrices

interface normalize
  module procedure dnormalize,znormalize
endinterface

interface lowdin
  module procedure dlowdin,zlowdin
endinterface

interface transform
  module procedure dtransform,ztransform
endinterface

interface matmultiply
  module procedure dmultiply,zmultiply
endinterface

interface matvecmultiply
  module procedure dvecmultiply,zvecmultiply
endinterface

interface diagonalize
  module procedure ddiagonalize,zdiagonalize
endinterface

interface exponentiate
  module procedure dexponential,zexponential
endinterface

interface diagonalize_and_project
  module procedure ddiagonalize_and_project,zdiagonalize_and_project
endinterface

interface project_recursive
  module procedure z_project_recursive
endinterface

interface matwrite
  module procedure dwrite,zwrite,iwrite
endinterface

interface vecwrite
  module procedure dvecwrite,zvecwrite,svecwrite
endinterface

interface vec3write
  module procedure d3vecwrite,z3vecwrite
endinterface

interface matread
  module procedure dread,zread,iread
endinterface

interface vecread
  module procedure dvecread,zvecread,svecread
endinterface

interface vec3read
  module procedure d3vecread,z3vecread
endinterface

interface ishermitian
  module procedure dishermitian,zishermitian
endinterface

interface isantihermitian
  module procedure disantihermitian,zisantihermitian
endinterface

interface isunitary
  module procedure disunitary,zisunitary
endinterface

interface project_a_on_b
  module procedure d3project_a_on_b, z3project_a_on_b
endinterface

interface intruder
  module procedure zintruder
endinterface

! =================================================================== !
! interfaces to lapack routines

interface

  subroutine DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
    CHARACTER          JOBZ, UPLO
    INTEGER            INFO, LDA, LWORK, N
    DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
  endsubroutine

  subroutine ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO )
    CHARACTER          JOBZ, UPLO
    INTEGER            INFO, LDA, LWORK, N
    DOUBLE PRECISION   RWORK( * ), W( * )
    COMPLEX*16         A( LDA, * ), WORK( * )
  endsubroutine

  subroutine DGEMM( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
    CHARACTER*1        TRANSA, TRANSB
    INTEGER            M, N, K, LDA, LDB, LDC
    DOUBLE PRECISION   ALPHA, BETA
    DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
  endsubroutine

  subroutine ZGEMM( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
    CHARACTER*1        TRANSA, TRANSB
    INTEGER            M, N, K, LDA, LDB, LDC
    COMPLEX*16         ALPHA, BETA
    COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * )
  endsubroutine

endinterface


! =================================================================== !

 contains

! =================================================================== !

!> Is called with the size of the matrices handled in SHARC (i.e. nstates)
!> and allocates the workspace arrays for the lapack routines.
!> \param n the size of the matrices for which the workspace should be allocated
subroutine allocate_lapack(n)
  implicit none

  integer,intent(in) :: n
  complex*16 :: Uz(n,n)  !< complex dummy matrix for zheev
  real*8 :: Ud(n,n)      !< real dummy matrix for dsyev
  real*8 :: w(n)         !< dummy eigenvalue array
  integer :: ierr

  Ud=0.d0
  lapack_lwork_d=-1
  allocate(lapack_work_d(1))
  call dsyev('V','L',n,Ud,n,w,lapack_work_d,lapack_lwork_d,ierr)
  if (ierr/=0) then
    write(0,*) 'LAPACK could not allocate workspace, ierr=',ierr,'lwork=',lapack_work_d(1)
    stop
  endif
  lapack_lwork_d=int(lapack_work_d(1))
  deallocate(lapack_work_d)
  allocate(lapack_work_d(lapack_lwork_d))

  Uz=dcmplx(0.d0,0.d0)
  lapack_lwork_z=-1
  allocate(lapack_work_z(1))
  allocate(lapack_rwork_z(max(1,3*n-2)))
  call zheev('V','L',n,Uz,n,w,lapack_work_z,lapack_lwork_z,lapack_rwork_z,ierr)
  if (ierr/=0) then
    write(0,*) 'LAPACK could not allocate workspace, ierr=',ierr,'lwork=',lapack_work_z(1)
    stop
  endif
  lapack_lwork_z=int(lapack_work_z(1))
  deallocate(lapack_work_z)
  allocate(lapack_work_z(lapack_lwork_z))

  lapack_isalloc=.true.

endsubroutine

! =================================================================== !

!> deallocates the lapack workspaces (needed if the workspaces need to be reallocated for larger matrices)
subroutine deallocate_lapack
  implicit none

  deallocate( lapack_work_d )
  deallocate( lapack_work_z, lapack_rwork_z )

endsubroutine

! =================================================================== !
! =================================================================== !
!                             Normalization                           !
! =================================================================== !
! =================================================================== !

!> normalises the columns of matrix A_ss
!>
!> A_ss is normalised in place
subroutine dnormalize(n,A_ss)
  implicit none
! parameters
  integer, intent(in) :: n                      !< size of matrix A
  real*8, intent(inout) :: A_ss(n,n)            !< the matrix to normalize
! internals
  real*8 :: sum_column  !< sum for normalisation
  integer :: i,j

  do j=1,n
    sum_column=0.d0
    do i=1,n
      sum_column=sum_column+A_ss(i,j)**2
    enddo
    sum_column=1/sqrt(sum_column)
    A_ss(:,j)=sum_column*A_ss(:,j)
  enddo

  return

endsubroutine

! =================================================================== !

!> normalises the columns of matrix A_ss
!> A_ss is normalised in place
subroutine znormalize(n,A_ss)
  implicit none
! parameters
  integer, intent(in) :: n                  !< size of matrix A
  complex*16, intent(inout) :: A_ss(n,n)    !< the matrix to normalize
! internals
  real*8 :: sum_column  !< sum for normalisation
  integer :: i,j

  do j=1,n
    sum_column=0.d0
    do i=1,n
      sum_column=real(sum_column+A_ss(i,j)*conjg(A_ss(i,j)))
    enddo
    sum_column=1.d0/sqrt(sum_column)
    A_ss(:,j)=dcmplx(sum_column,0.d0)*A_ss(:,j)
  enddo

  return

endsubroutine

! =================================================================== !
! =================================================================== !
!                          Lowdin orthogonalisation                   !
! =================================================================== !
! =================================================================== !

subroutine zintruder(n,A_ss)
  use definitions, only: u_log
  implicit none
  integer, intent(in) :: n
  complex*16, intent(inout) :: A_ss(n,n)
  real*8 :: sums
  integer :: i,j,k
  real*8,parameter :: intr_thrs=1.d-1

    ! Intruder state check
    do i=1,n
      sums=0.d0
      do j=1,n
        sums=sums+abs(A_ss(i,j))**2
        sums=sums+abs(A_ss(j,i))**2
      enddo
      sums=sums-abs(A_ss(i,i))**2

      if (sums < intr_thrs) then
        write(u_log,'(A)') '! ======== INTRUDER STATE PROBLEM ======== !'
        write(u_log,'(A,I4)') 'State: ',i
        do k=1,n
          write(u_log,'(1000(F8.5,1X))') (A_ss(k,j),j=1,n)
        enddo

        A_ss(i,:)=dcmplx(0.d0,0.d0)
        A_ss(:,i)=dcmplx(0.d0,0.d0)
        A_ss(i,i)=dcmplx(1.d0,0.d0)
      endif
    enddo

endsubroutine

! =================================================================== !

!> takes a real matrix A and orthonormalises it
!> using Löwdin's symmetric orthogonalisation
subroutine dlowdin(n,A_ss)
  implicit none
! parameters
  integer,intent(in) :: n                          !< size of matrix A
  real*8,intent(inout) :: A_ss(n,n)                !< matrix A
! internal variables
  integer :: i                                     !< loop counters
  real*8 :: S_ss(n,n),U_ss(n,n),Scratch_ss(n,n)    !< three matrices for matrix-multiplication

! build real symmetric matrix S=A^T.A
  call dgemm('T','N',n,n,n,1.d0,A_ss,n,A_ss,n,0.d0,S_ss,n)

! calculate eigenvalues and eigenvectors of S via dsyev wrapper
  call ddiagonalize(n,S_ss,U_ss)

! calculate S^(-1/2)
! first get inverse square root of diagonal matrix
  do i=1,n
    S_ss(i,i)=1.d0/sqrt(S_ss(i,i))
  enddo
! second transform back 
  call dtransform(n,S_ss,U_ss,'uaut')

! calculate the orthonormalised matrix A=S.A
  call dgemm('N','N',n,n,n,1.d0,A_ss,n,S_ss,n,0.d0,Scratch_ss,n)
  A_ss=Scratch_ss

! normalisation of the columns
  call dnormalize(n,A_ss)

  return

endsubroutine

! =================================================================== !

!> takes a complex matrix A and orthonormalises it
!> uses Löwdin's symmetric orthogonalisation
subroutine zlowdin(n,A_ss)
  implicit none
! parameters
  integer,intent(in) :: n                               !< size of matrix A
  complex*16,intent(inout) :: A_ss(n,n)                 !< matrix A
! internal variables
  integer :: i                                          !< loop counters
  complex*16 :: S_ss(n,n),U_ss(n,n),Scratch_ss(n,n)     !< three matrices for matrix-multiplication

! build hermitian matrix S=A^dagger.A
  call zgemm('C','N',n,n,n,dcmplx(1.d0,0.d0),A_ss,n,A_ss,n,dcmplx(0.d0,0.d0),S_ss,n)

! calculate eigenvalues and eigenvectors of S via zheev wrapper
  call zdiagonalize(n,S_ss,U_ss)

! calculate S^(-1/2)
! first get inverse square root of diagonal matrix
  do i=1,n
    S_ss(i,i)=dcmplx(1.d0,0.d0)/sqrt(S_ss(i,i))
  enddo
! second transform back 
  call ztransform(n,S_ss,U_ss,'uaut')

! calculate the orthonormalised matrix A=S.A
  call zgemm('N','N',n,n,n,dcmplx(1.d0,0.d0),A_ss,n,S_ss,n,dcmplx(0.d0,0.d0),Scratch_ss,n)
  A_ss=Scratch_ss

! normalisation of the columns
  call znormalize(n,A_ss)

  return

endsubroutine

! =================================================================== !
! =================================================================== !
!                             Transformation                          !
! =================================================================== !
! =================================================================== !

!> calculates the (orthogonal) transformations U^T.A.U or U.A.U^T
!>
!> \param mode can be 'utau' or 'uaut' to calculate U^T.A.U or U.A.U^T, respectively
subroutine dtransform(n,A_ss,U_ss,mode)
  implicit none
! parameters
  integer, intent(in) :: n
  real*8, intent(inout) :: A_ss(n,n)
  real*8, intent(in) :: U_ss(n,n)
  character*4, intent(in):: mode
! internal variables
  real*8 :: Scratch_ss(n,n)

  if (mode=='utau') then
    call dgemm('N','N',n,n,n,1.d0,A_ss,n,U_ss,n,0.d0,Scratch_ss,n)
    call dgemm('T','N',n,n,n,1.d0,U_ss,n,Scratch_ss,n,0.d0,A_ss,n)
  elseif (mode=='uaut') then
    call dgemm('N','T',n,n,n,1.d0,A_ss,n,U_ss,n,0.d0,Scratch_ss,n)
    call dgemm('N','N',n,n,n,1.d0,U_ss,n,Scratch_ss,n,0.d0,A_ss,n)
  else
    write(0,*) 'Unknown transformation mode in dtransform'
  endif

  return

endsubroutine

! =================================================================== !

!> calculates the (unitary) transformations U^dagger.A.U or U.A.U^dagger
!>
!> \param mode can be 'utau' or 'uaut' to calculate U^dagger.A.U or U.A.U^dagger, respectively
subroutine ztransform(n,A_ss,U_ss,mode)
  implicit none
! parameters
  integer, intent(in) :: n
  complex*16, intent(inout) :: A_ss(n,n)
  complex*16, intent(in) :: U_ss(n,n)
  character*4, intent(in):: mode
! internal variables
  complex*16 :: Scratch_ss(n,n)

  if (mode=='utau') then
    call zgemm('N','N',n,n,n,dcmplx(1.d0,0.d0),A_ss,n,U_ss,n,dcmplx(0.d0,0.d0),Scratch_ss,n)
    call zgemm('C','N',n,n,n,dcmplx(1.d0,0.d0),U_ss,n,Scratch_ss,n,dcmplx(0.d0,0.d0),A_ss,n)
  elseif (mode=='uaut') then
    call zgemm('N','C',n,n,n,dcmplx(1.d0,0.d0),A_ss,n,U_ss,n,dcmplx(0.d0,0.d0),Scratch_ss,n)
    call zgemm('N','N',n,n,n,dcmplx(1.d0,0.d0),U_ss,n,Scratch_ss,n,dcmplx(0.d0,0.d0),A_ss,n)
  else
    write(0,*) 'Unknown transformation mode in ztransform'
  endif

  return

endsubroutine

! =================================================================== !
! =================================================================== !
!                             Diagonalization                         !
! =================================================================== !
! =================================================================== !

!> decomposes a symmetric matrix A into U.A_diag.U^T
!>
!> is a wrapper around dsyev
!> \param A_ss on entry the matrix to be diagonalized, on exit the diagonal matrix with the eigenvalues on the diagonal (with increasing values)
!> \param U_ss on exit the matrix of eigenvectors
subroutine ddiagonalize(n,A_ss,U_ss)
implicit none
! parameters
  integer, intent(in) :: n
  real*8, intent(inout) :: A_ss(n,n)
  real*8, intent(out) :: U_ss(n,n)
! internal variables
  real*8 :: EV_s(n)
  integer :: io,i

  U_ss=A_ss
  call dsyev('V','L',n,U_ss,n,EV_s,lapack_work_d,lapack_lwork_d,io)
! U_ss already holds the transformation matrix
! now  building the diagonal matrix A_ss
  A_ss=0.d0
  do i=1,n
    A_ss(i,i)=EV_s(i)
  enddo

  return

endsubroutine

! =================================================================== !

!> decomposes a hermitian matrix A into U.A_diag.U^dagger
!>
!> is a wrapper around zheev
!> \param A_ss on entry the matrix to be diagonalized, on exit the diagonal matrix with the eigenvalues on the diagonal (with increasing values)
!> \param U_ss on exit the matrix of eigenvectors
subroutine zdiagonalize(n,A_ss,U_ss)
  implicit none
! parameters
  integer, intent(in) :: n
  complex*16, intent(inout) :: A_ss(n,n)
  complex*16, intent(out) :: U_ss(n,n)
! internal variables
  real*8 :: EV_s(n)     ! holds eigenvalues after dsyev
  integer :: io,i

  U_ss=A_ss
  call zheev('V','L',n,U_ss,n,EV_s,lapack_work_z,lapack_lwork_z,lapack_rwork_z,io)      ! actual diagonalisation
! U_ss already holds the transformation matrix
! now  building the diagonal matrix A_ss
  A_ss=dcmplx(0.d0,0.d0)
  do i=1,n
    A_ss(i,i)=dcmplx(EV_s(i),0.d0)
  enddo

  return

endsubroutine

! =================================================================== !
! =================================================================== !
! =================================================================== !

!> sets the variable diagonalize_degeneracy_diff
!> which is used in zdiagonalize_and_project() and ddiagonalize_and_project()
subroutine set_project_diffthr(value)
  implicit none
  real*8 :: value

  diagonalize_degeneracy_diff=value

  return
endsubroutine

! =================================================================== !

!> Diagonalizes matrix H to obtain the eigenvalues and eigenvectors in U
!> 
!> then transforms U to make it as similar as possible to U_old
!> by adjusting the complex phases of the columns of U
!> while maintaining U^dagger.H.U=diag
!>
subroutine zdiagonalize_and_project(n,H,U,Uold)
  implicit none
  ! parameters
  integer, intent(in) :: n
  complex*16, intent(inout) :: H(n,n)
  complex*16, intent(out) :: U(n,n)
  complex*16, intent(in) :: Uold(n,n)
  ! internal variables
  complex*16 :: S(n,n), P(n,n), Pafter(n,n), tempM(n,n), temp(n), eigv(n), temph, Hinput(n,n)
  logical :: mask(n,n)
  real*8 :: sums
  integer :: i,j,k

  ! diagonalise H
!           Hinput=H
!           call matwrite(n,H,6,'H','ES24.16E2')
  call zdiagonalize(n,H,U)
!           call transform(n,Hinput,U,'utau')
!           call matwrite(n,Hinput,6,'U^tHU','F12.9')
  do i=1,n
    eigv(i)=H(i,i)
  enddo

!           call matwrite(n,Uold,6,'Uold','ES24.16E2')
!           call matwrite(n,U,6,'U','ES24.16E2')

  ! calculate overlap of U and Uold
  call zmultiply(n,U,Uold,S,'tn')

!           call matwrite(n,S,6,'S before','ES24.16E2')

  ! scan rows for values abs(s) close to 1
  ! reorder these on the diagonal
  ! to, e.g., give vector 1 the phase of the old vector 2, if their overlap is very large
  do i=1,n
    do j=1,n
      if (abs(S(i,j))**2>0.5d0) then
!         write(6,*) i,j
        ! do not reorder if degenerate
        if (abs(eigv(i)-eigv(j))<diagonalize_degeneracy_diff) exit
        if (i==j) then
          exit
        elseif (i<j) then
          ! sort to the right
          temp=S(:,i)
          S(:,i)=S(:,j)
          S(:,i+2:j)=S(:,i+1:j-1)
          S(:,i+1)=temp

          temph=eigv(i)
          eigv(i)=eigv(j)
          eigv(i+2:j)=eigv(i+1:j-1)
          eigv(i+1)=temph
        elseif (i>j) then
          ! sort to the left
          temp=S(:,i)
          S(:,i)=S(:,j)
          S(:,j:i-2)=S(:,j+1:i-1)
          S(:,i-1)=temp

          temph=eigv(i)
          eigv(i)=eigv(j)
          eigv(j:i-2)=eigv(j+1:i-1)
          eigv(i-1)=temph
        endif
!         call matwrite(n,S,6,'S during','ES24.16E2')
        exit
      endif
    enddo
  enddo
  tempM=S

!           call matwrite(n,S,6,'S after','ES24.16E2')
  ! initialise projection matrix
  P=dcmplx(0.d0,0.d0)

  mask=.false.
  ! construct the projection matrix, depending on the degeneracy of H
  do i=1,n
    do j=1,n
      if (abs(eigv(i)-eigv(j))<diagonalize_degeneracy_diff) then
        mask(i,j)=.true.
        P(i,j)=S(i,j)
      endif
    enddo
  enddo

!           call matwrite(n,P,6,'P before','ES24.16E2')

  ! If a vector in P is zero, replace it with a (0..1..0) vector
  ! should not occur often after the resort above
  do i=1,n
    sums=0.d0
    do j=1,n
      sums=sums+abs(P(j,i))**2
    enddo

    if (sums < 1d-9) then
      P(:,i)=dcmplx(0.d0,0.d0)

      do j=1,n
        if (mask(j,i)) then
          sums=0.d0
          do k=1,n
            sums=sums+abs(P(j,k))**2
          enddo
          if (sums < 1d-9) then
            P(j,i)=dcmplx(1.d0,0.d0)
            exit
          endif
        endif
      enddo
    endif
  enddo

!           call matwrite(n,P,6,'P after','ES24.16E2')

  ! orthonormalize P
  Pafter=P
  call zlowdin(n,P)
!           call matwrite(n,P,6,'P matrix','ES24.16E2')

  ! check for NaNs in the P matrix
  if ( any( (real(P)).ne.(real(P)) ).or.any( (aimag(P)).ne.(aimag(P)) ) ) then
    write(0,*) 'NaN during zdiagonalize_and_project.'
    write(0,*) 'This happens if the phase-tracking algorithm fails.'
    call matwrite(n,Uold,0,'Uold','ES24.16E2')
    call matwrite(n,U,0,'U','ES24.16E2')
    call matwrite(n,H,0,'H','ES24.16E2')
    write(0,*) 'mask'
    do i=1,n
      write(0,*) (mask(i,j),j=1,n)
    enddo
    call matwrite(n,S,0,'S','ES24.16E2')
    call matwrite(n,tempM,0,'S after reordering','ES24.16E2')
    call matwrite(n,Pafter,0,'P before Loewdin step','ES24.16E2')
    call matwrite(n,P,0,'P after Loewdin step','ES24.16E2')
    write(0,*) 
    write(0,*) '#=============================================#'
    write(0,*) '      Restarting U matrix phase tracking'
    write(0,*) '#=============================================#'
    write(0,*) 
    P=dcmplx(0.d0,0.d0)
    do i=1,n
      P(i,i)=dcmplx(1.d0,0.d0)
    enddo
  endif

  ! the continuous U matrix is obtained from U=U.P
  ! use S as temporary memory, since it is not needed anymore
  call zmultiply(n,U,P,S,'nn')

!           call transform(n,Hinput,U,'utau')
!           call matwrite(n,Hinput,6,'U^tHU after all','ES24.16E2')
  U=S
!           call matwrite(n,Uold,6,'Uold after all','ES24.16E2')
!           call matwrite(n,U,6,'U after all','ES24.16E2')
!           call transform(n,Hinput2,U,'utau')
!           call matwrite(n,Hinput2,6,'U^tHU after all','ES24.16E2')

  return

endsubroutine

! =================================================================== !
!> Diagonalizes matrix H to obtain the eigenvalues and eigenvectors in U
!> 
!> then transforms U to make it as similar as possible to U_old
!> by adjusting the complex phases of the columns of U
!> while maintaining U^dagger.H.U=diag
!>
!> \todo some code snippets from zdiagonalize_and_project() are missing, since ddiagonalize_and_project() is not used in SHARC
subroutine ddiagonalize_and_project(n,H,U,Uold)
  implicit none
  ! parameters
  integer, intent(in) :: n
  real*8, intent(inout) :: H(n,n)
  real*8, intent(out) :: U(n,n)
  real*8, intent(in) :: Uold(n,n)
  ! internal variables
  real*8 :: S(n,n), P(n,n)
  real*8 :: sums
  integer :: i,j

  ! diagonalise H
  call ddiagonalize(n,H,U)
  ! should make this inline

  ! calculate overlap of U and Uold
  call dmultiply(n,U,Uold,S,'tn')
  ! initialise projection matrix
  P=0.d0

  ! construct the projection matrix, depending on the degeneracy of H
  do i=1,n
    do j=1,n
      if (dabs(H(i,i)-H(j,j))<diagonalize_degeneracy_diff) P(i,j)=S(i,j)
    enddo
  enddo

  ! TODO: this may be solved more elegantly
  ! If a vector in P is zero, replace it with a (0..1..0) vector
  do i=1,n
    sums=0.d0
    do j=1,n
      sums=sums+P(j,i)**2
    enddo

    if (sums < 1d-9) then
      P(:,i)=0.d0
      P(i,i)=1.d0
    endif
  enddo

  ! orthonormalize P
  call dlowdin(n,P)

  ! the continuous U matrix is obtained from U=U.P
  ! use S as temporary memory, since it is not needed anymore
  call dmultiply(n,U,P,S,'nn')

  U=S

  return

endsubroutine

! =================================================================== !

!> linearly interpolates between SOold and SO,
!> diagonalizing the interpolated matrices and 
!> in each step adjusts the phases of the U matrix to the previous one
!> if the U matrix changes too fast, call this subroutine recursively with a smaller interpolation step
recursive subroutine z_project_recursive(n, SOold, SO, Uold, U, dt, printlevel, u_print)
  implicit none

  integer, intent(in) :: n                              !< size of the matrices
  complex*16, intent(in) :: SOold(n,n)                  !< old matrix with known U matrix
  complex*16, intent(in) :: SO(n,n)                     !< new matrix to diagonalize
  complex*16, intent(in) :: Uold(n,n)                   !< known old U matrix
  complex*16, intent(out) :: U(n,n)                     !< the goal is to calculate this as the matrix which diagonalizes SO, while being as similar to Uold as possible
  real*8, intent(in) :: dt                              !< time in a.u. between SOold and SO, needed for limiting the recursion depth
  integer, intent(in) :: printlevel, u_print

  integer :: i,istep
  complex*16 :: H(n,n), Hold(n,n), dU(n,n), UdU(n,n), Uold_dummy(n,n)
  integer, parameter :: nsubsteps=10
  real*8 :: dtsubstep
  logical :: sub

  real*8, parameter :: diagthrs=1.d-3           !< threshold for the diagonal elements of U^dagger.dU
                                                !< controls whether to enter the next recursion level
  real*8, parameter :: mindtsubstep=1.d-3       !< minimum dt
                                                !< controls how deep the recursion progresses
  real*8, parameter :: dt_reduce=0.1d0          !< factor by which dt is reduced in substeps
                                                !< reciproke value of nsubsteps \todo calculate dt_reduce from nsubsteps

dtsubstep=dt/nsubsteps
Uold_dummy=Uold

do istep=1,nsubsteps

!   write(0,*) k, dtsubstep

  H=SOold + (SO-SOold)*istep/nsubsteps
  call zdiagonalize_and_project(n, H, U, Uold_dummy)

  dU=(U-Uold_dummy)/dtsubstep
  call zmultiply(n,U,dU,UdU,'tn')
  sub=.false.
  do i=1,n
    if (abs(UdU(i,i))>diagthrs) then
      sub=.true.
      exit
    endif
  enddo

  if ( (sub) .and. (dtsubstep*dt_reduce>mindtsubstep) ) then

    H=SOold + (SO-SOold)*istep/nsubsteps
    Hold=SOold + (SO-SOold)*(istep-1)/nsubsteps
    if (printlevel>6) then
      write(u_print,*) 'dtsubstep=',dtsubstep
    endif
    call z_project_recursive(n, Hold, H, Uold_dummy, U, dtsubstep, printlevel, u_print)

  else

    Uold_dummy=U

  endif

enddo

endsubroutine

! =================================================================== !
! =================================================================== !
!                             Matrix functionals                      !
! =================================================================== !
! =================================================================== !

!> calculates the matrix exponential e^A
!> by e^A = U.e^(factor*U^T.A.U).U^T
!>
!> is a wrapper around dsyev
subroutine dexponential(n,A_ss,factor)
implicit none
! parameters
  integer, intent(in) :: n
  real*8, intent(inout) :: A_ss(n,n)
  real*8, intent(in) :: factor
! internal variables
  real*8 :: EV_s(n),Scratch_ss(n,n),U_ss(n,n)
  integer :: io,i

  U_ss=A_ss
  call dsyev('V','L',n,U_ss,n,EV_s,lapack_work_d,lapack_lwork_d,io)
! U_ss already holds the transformation matrix
! now  building the diagonal matrix A_ss
  A_ss=0.d0
  do i=1,n
    A_ss(i,i)=exp(factor*EV_s(i))
  enddo
  call dgemm('N','T',n,n,n,1.d0,A_ss,n,U_ss,n,0.d0,Scratch_ss,n)
  call dgemm('N','N',n,n,n,1.d0,U_ss,n,Scratch_ss,n,0.d0,A_ss,n)

  return

endsubroutine

! =================================================================== !

!> calculates the matrix exponential e^A
!> by e^A = U.e^(factor*U^T.A.U).U^T
!>
!> is a wrapper around zheev
!> \param factor allows to calculate matrix exponentials of non-hermitian matrices
!> for an anti-hermitian matrix A, call dexponential(n,ii*A,-ii), then ii*A is hermitian, but the result is e^A
subroutine zexponential(n,A_ss,factor)
  implicit none
! parameters
  integer, intent(in) :: n
  complex*16, intent(inout) :: A_ss(n,n)
  complex*16, intent(in) :: factor
! internal variables
  real*8 :: EV_s(n)
  integer :: io,i
  complex*16 :: U_ss(n,n),Scratch_ss(n,n)

  U_ss=A_ss
  call zheev('V','L',n,U_ss,n,EV_s,lapack_work_z,lapack_lwork_z,lapack_rwork_z,io)
! U_ss already holds the transformation matrix
! now  building the diagonal matrix A_ss
  A_ss=dcmplx(0.d0,0.d0)
  do i=1,n
    A_ss(i,i)=exp(factor*dcmplx(EV_s(i),0.d0))
  enddo
  call zgemm('N','C',n,n,n,dcmplx(1.d0,0.d0),A_ss,n,U_ss,n,dcmplx(0.d0,0.d0),Scratch_ss,n)
  call zgemm('N','N',n,n,n,dcmplx(1.d0,0.d0),U_ss,n,Scratch_ss,n,dcmplx(0.d0,0.d0),A_ss,n)

  return

endsubroutine

! =================================================================== !
! =================================================================== !
!                             Matrix multiplication                   !
! =================================================================== !
! =================================================================== !

!> calculates A.B, A^dagger.B, A.B^dagger or A^dagger.B^dagger
!>
!> wrapper around zgemm
!> \param mode meaning: 'nn' for A.B, 'tn' for A^dagger.B, etc.
subroutine zmultiply(n,A_ss,B_ss,R_ss,mode)
  implicit none
! parameters
  integer, intent(in) :: n
  complex*16, intent(in) :: A_ss(n,n),B_ss(n,n)
  complex*16, intent(out) :: R_ss(n,n)
  character*2, intent(in) :: mode
! internal variables
  character :: Amode, Bmode

  select case (mode)
    case ('nn')
      Amode='N'
      Bmode='N'
    case ('tn')
      Amode='C'
      Bmode='N'
    case ('nt')
      Amode='N'
      Bmode='C'
    case ('tt')
      Amode='C'
      Bmode='C'
    case default
      stop 'Unknown multiplication mode in zmultiply'
  endselect

  call zgemm(Amode,Bmode,n,n,n,dcmplx(1.d0,0.d0),A_ss,n,B_ss,n,dcmplx(0.d0,0.d0),R_ss,n)

  return

endsubroutine

! =================================================================== !

! calculates A.B, A^T.B, A.B^T or A^T.B^T
!
! wrapper around dgemm
!> \param mode meaning: 'nn' for A.B, 'tn' for A^T.B, etc.
subroutine dmultiply(n,A_ss,B_ss,R_ss,mode)
  implicit none
! parameters
  integer, intent(in) :: n
  real*8, intent(in) :: A_ss(n,n),B_ss(n,n)
  real*8, intent(out) :: R_ss(n,n)
  character*2, intent(in) :: mode
! internal variables
  character :: Amode, Bmode

  select case (mode)
    case ('nn')
      Amode='N'
      Bmode='N'
    case ('tn')
      Amode='T'
      Bmode='N'
    case ('nt')
      Amode='N'
      Bmode='T'
    case ('tt')
      Amode='T'
      Bmode='T'
    case default
      stop 'Unknown multiplication mode in dmultiply'
  endselect

  call dgemm(Amode,Bmode,n,n,n,1.d0,A_ss,n,B_ss,n,0.d0,R_ss,n)

  return

endsubroutine

! =================================================================== !

! calculates matrix-vector products U^T.c and U.c
!
! wrapper around dgemm
!> \param mode meaning: 'n' for A.c, 't' for A^T.c
subroutine dvecmultiply(n,A_ss,c_s,cresult_s,mode)
  implicit none
! parameters
  integer, intent(in) :: n
  real*8, intent(in) :: A_ss(n,n),c_s(n)
  real*8, intent(out) :: cresult_s(n)
  character, intent(in) :: mode
! internal variables
  character :: Amode

  select case (mode)
    case ('n')
      Amode='N'
    case ('t')
      Amode='T'
    case default
      stop 'Unknown multiplication mode in dvecmultiply'
  endselect

  call dgemm(Amode,'N',n,1,n,1.d0,A_ss,n,c_s,n,0.d0,cresult_s,n)

  return

endsubroutine

! =================================================================== !

! calculates U^T.c and U.c
!
! wrapper around zgemm
!> \param mode meaning: 'n' for A.c, 't' for A^dagger.c
subroutine zvecmultiply(n,A_ss,c_s,cresult_s,mode)
  implicit none
! parameters
  integer, intent(in) :: n
  complex*16, intent(in) :: A_ss(n,n),c_s(n)
  complex*16, intent(out) :: cresult_s(n)
  character, intent(in) :: mode
! internal variables
  character :: Amode

  select case (mode)
    case ('n')
      Amode='N'
    case ('t')
      Amode='C'
    case default
      stop 'Unknown multiplication mode in zvecmultiply'
  endselect

  call zgemm(Amode,'N',n,1,n,dcmplx(1.d0,0.d0),A_ss,n,c_s,n,dcmplx(0.d0,0.d0),cresult_s,n)

  return

endsubroutine

! =================================================================== !
! =================================================================== !
!                                    Write                            !
! =================================================================== !
! =================================================================== !

!> writes matrix A to unit wrunit
!> writes the title before the matrix
!> uses precstring as format string
subroutine dwrite(n,A,wrunit,title,precstring)
  implicit none
  ! parameters
  integer, intent(in) :: n
  real*8,intent(in) :: A(n,n)
  integer, intent(in) :: wrunit
  character(len=*), intent(in) :: title
  character(len=*), intent(in) :: precstring
  ! internal variables
  integer :: i,j
  character*255 :: fmtstring

  write(wrunit,'(A)') trim(title)

  write(fmtstring,'(I10)') n
  fmtstring='('//trim(adjustl(fmtstring))//'('//trim(adjustl(precstring))//',1X))'
  do i=1,n
    write(wrunit,fmtstring) (A(i,j),j=1,n)
  enddo

endsubroutine

! =================================================================== !

!> writes matrix A to unit wrunit
!> writes the title before the matrix
!> uses precstring as format string
subroutine zwrite(n,A,wrunit,title,precstring)
  implicit none
  ! parameters
  integer, intent(in) :: n
  complex*16,intent(in) :: A(n,n)
  integer, intent(in) :: wrunit
  character(len=*), intent(in) :: title
  character(len=*), intent(in) :: precstring
  ! internal variables
  character*100 :: fmtstring
  integer :: i,j

  write(wrunit,'(A)') trim(title)

  write(fmtstring,'(I10)') n
  fmtstring='('//trim(adjustl(fmtstring))//'('//trim(adjustl(precstring))//',1X,'//trim(adjustl(precstring))//',4X))'
  do i=1,n
    write(wrunit,fmtstring) (A(i,j),j=1,n)
  enddo

endsubroutine

! =================================================================== !

!> writes matrix A to unit wrunit
!> writes the title before the matrix
!> uses precstring as format string
subroutine iwrite(n,A,wrunit,title,precstring)
  implicit none
  ! parameters
  integer, intent(in) :: n
  integer, intent(in) :: A(n,n)
  integer, intent(in) :: wrunit
  character(len=*), intent(in) :: title
  character(len=*), intent(in) :: precstring
  ! internal variables
  integer :: i,j
  character*255 :: fmtstring

  write(wrunit,'(A)') trim(title)

  write(fmtstring,'(I10)') n
  fmtstring='('//trim(adjustl(fmtstring))//'('//trim(adjustl(precstring))//',1X))'
  do i=1,n
    write(wrunit,fmtstring) (A(i,j),j=1,n)
  enddo

endsubroutine

! =================================================================== !

!> writes vector c to unit wrunit
!> writes the title before the vector
!> uses precstring as format string
subroutine dvecwrite(n,c,wrunit,title,precstring)
  implicit none
  ! parameters
  integer, intent(in) :: n
  real*8,  intent(in) :: c(n)
  integer, intent(in) :: wrunit
  character(len=*), intent(in) :: title
  character(len=*), intent(in) :: precstring
  ! internal variables
  integer :: i
  character*255 :: fmtstring

  write(wrunit,'(A)') trim(title)

  write(fmtstring,'(I10)') n
  fmtstring='('//trim(adjustl(fmtstring))//'('//trim(adjustl(precstring))//',1X))'
  do i=1,n
    write(wrunit,fmtstring) c(i)
  enddo

endsubroutine

! =================================================================== !

!> writes vector c to unit wrunit
!> writes the title before the vector
!> uses precstring as format string
subroutine zvecwrite(n,c,wrunit,title,precstring)
  implicit none
  ! parameters
  integer, intent(in) :: n
  complex*16,intent(in) :: c(n)
  integer, intent(in) :: wrunit
  character(len=*), intent(in) :: title
  character(len=*), intent(in) :: precstring
  ! internal variables
  integer :: i
  character*255 :: fmtstring

  write(wrunit,'(A)') trim(title)

  write(fmtstring,'(I10)') n
  fmtstring='('//trim(adjustl(fmtstring))//'('//trim(adjustl(precstring))//',1X,'//trim(adjustl(precstring))//',4X))'
  do i=1,n
    write(wrunit,fmtstring) c(i)
  enddo

endsubroutine

! =================================================================== !

!> writes vector c to unit wrunit
!> writes the title before the vector
!> uses precstring as format string
subroutine svecwrite(n,c,wrunit,title,precstring)
  implicit none
  ! parameters
  integer, intent(in) :: n
  character(len=*),  intent(in) :: c(n)
  integer, intent(in) :: wrunit
  character(len=*), intent(in) :: title
  character(len=*), intent(in) :: precstring
  ! internal variables
  integer :: i
  character*255 :: fmtstring

  write(wrunit,'(A)') trim(title)

  write(fmtstring,'(I10)') n
  fmtstring='('//trim(adjustl(fmtstring))//'('//trim(adjustl(precstring))//',1X))'
  do i=1,n
    write(wrunit,fmtstring) c(i)
  enddo

endsubroutine

! =================================================================== !

!> writes vector c of dimension 3xn to unit wrunit
!> writes the title before the vector
!> uses precstring as format string
subroutine d3vecwrite(n,c,wrunit,title,precstring)
  implicit none
  ! parameters
  integer, intent(in) :: n
  real*8,intent(in) :: c(n,3)
  integer, intent(in) :: wrunit
  character(len=*), intent(in) :: title
  character(len=*), intent(in) :: precstring
  ! internal variables
  integer :: i,j
  character*255 :: fmtstring

  write(wrunit,'(A)') trim(title)

  write(fmtstring,'(I10)') 3
  fmtstring='('//trim(adjustl(fmtstring))//'('//trim(adjustl(precstring))//',1X))'
  do i=1,n
    write(wrunit,fmtstring) (c(i,j),j=1,3)
  enddo

endsubroutine

! =================================================================== !

!> writes vector c of dimension 3xn to unit wrunit
!> writes the title before the vector
!> uses precstring as format string
subroutine z3vecwrite(n,c,wrunit,title,precstring)
  implicit none
  ! parameters
  integer, intent(in) :: n
  complex*16,intent(in) :: c(n,3)
  integer, intent(in) :: wrunit
  character(len=*), intent(in) :: title
  character(len=*), intent(in) :: precstring
  ! internal variables
  integer :: i,j
  character*255 :: fmtstring

  write(wrunit,'(A)') trim(title)

  write(fmtstring,'(I10)') 3
  fmtstring='('//trim(adjustl(fmtstring))//'('//trim(adjustl(precstring))//',1X,'//trim(adjustl(precstring))//',4X))'
  do i=1,n
    write(wrunit,fmtstring) (c(i,j),j=1,3)
  enddo

endsubroutine

! =================================================================== !
! =================================================================== !
!                                    Read                             !
! =================================================================== !
! =================================================================== !

!> reads matrix A from unit runit
!> also reads the title, THIS MEANS THAT the current line of runit has to be the title line
subroutine dread(n,A,runit,title)
  implicit none
  ! parameters
  integer, intent(in) :: n
  integer, intent(in) :: runit
  real*8,intent(out) :: A(n,n)
  character(len=8000), intent(out) :: title
  ! internal variables
  integer :: i,j, io

  read(runit,'(A)', iostat=io) title

  do i=1,n
    read(runit,*, iostat=io) (A(i,j),j=1,n)
    if (io/=0) then
      write(*,*) 'Could not read matrix'
      write(*,*) 'routine=dread(), n=',n,', unit=',runit
      write(*,*) 'title=',trim(title)
    endif
  enddo

endsubroutine

! =================================================================== !

!> reads matrix A from unit runit
!> also reads the title, THIS MEANS THAT the current line of runit has to be the title line
subroutine zread(n,A,runit,title)
  implicit none
  ! parameters
  integer, intent(in) :: n
  integer, intent(in) :: runit
  complex*16,intent(out) :: A(n,n)
  character(len=8000), intent(out) :: title
  ! internal variables
  integer :: i,j, io
  real*8 :: line(2*n)
  ! variables for line splitting
!   character*8000, allocatable :: values(:)
!   integer :: nvalues
!   character*8000 :: line
!   real*8 :: re, im

  read(runit,'(A)', iostat=io) title

  do i=1,n
    read(runit,*,iostat=io) (line(j),j=1,2*n)
    if (io/=0) then
      write(*,*) 'Could not read matrix'
      write(*,*) 'routine=zread(), n=',n,', unit=',runit
      write(*,*) 'title=',trim(title)
      stop 1
    endif
    do j=1,n
      A(i,j)=dcmplx(line(2*j-1),line(2*j))
    enddo
  enddo

endsubroutine

! =================================================================== !

!> reads matrix A from unit runit
!> also reads the title, THIS MEANS THAT the current line of runit has to be the title line
subroutine iread(n,A,runit,title)
  implicit none
  ! parameters
  integer, intent(in) :: n
  integer, intent(in) :: runit
  integer,intent(out) :: A(n,n)
  character(len=8000), intent(out) :: title
  ! internal variables
  integer :: i,j, io

  read(runit,'(A)', iostat=io) title

  do i=1,n
    read(runit,*, iostat=io) (A(i,j),j=1,n)
    if (io/=0) then
      write(*,*) 'Could not read matrix'
      write(*,*) 'routine=iread(), n=',n,', unit=',runit
      write(*,*) 'title=',trim(title)
    endif
  enddo

endsubroutine

! =================================================================== !

!> reads vector c from unit runit
!> also reads the title, THIS MEANS THAT the current line of runit has to be the title line
subroutine dvecread(n,c,runit,title)
  implicit none
  ! parameters
  integer, intent(in) :: n
  real*8,intent(out) :: c(n)
  integer, intent(in) :: runit
  character(len=8000), intent(out) :: title
  ! internal variables
  integer :: i, io

  read(runit,'(A)') title

  do i=1,n
    read(runit,*, iostat=io) c(i)
    if (io/=0) then
      write(*,*) 'Could not read vector'
      write(*,*) 'routine=dvecread(), n=',n,', unit=',runit
      write(*,*) 'title=',trim(title)
    endif
  enddo

endsubroutine

! =================================================================== !

!> reads vector c from unit runit
!> also reads the title, THIS MEANS THAT the current line of runit has to be the title line
subroutine zvecread(n,c,runit,title)
  implicit none
  ! parameters
  integer, intent(in) :: n
  complex*16,intent(out) :: c(n)
  integer, intent(in) :: runit
  character(len=8000), intent(out) :: title
  ! internal variables
  integer :: i,io
  real*8 :: re,im

  read(runit,'(A)') title

  do i=1,n
    read(runit,*, iostat=io) re,im
    if (io/=0) then
      write(*,*) 'Could not read vector'
      write(*,*) 'routine=zvecread(), n=',n,', unit=',runit
      write(*,*) 'title=',trim(title)
    endif
    c(i)=dcmplx(re,im)
  enddo

endsubroutine

! =================================================================== !

!> reads vector c from unit runit
!> also reads the title, THIS MEANS THAT the current line of runit has to be the title line
subroutine svecread(n,c,runit,title)
  implicit none
  ! parameters
  integer, intent(in) :: n
  character(len=*),intent(out) :: c(n)
  integer, intent(in) :: runit
  character(len=8000), intent(out) :: title
  ! internal variables
  integer :: i, io

  read(runit,'(A)') title

  do i=1,n
    read(runit,'(A)', iostat=io) c(i)
    if (io/=0) then
      write(*,*) 'Could not read vector'
      write(*,*) 'routine=svecread(), n=',n,', unit=',runit
      write(*,*) 'title=',trim(title)
    endif
  enddo

endsubroutine

! =================================================================== !

!> reads a 3-vector c from unit runit
!> also reads the title, THIS MEANS THAT the current line of runit has to be the title line
subroutine d3vecread(n,c,runit,title)
  implicit none
  ! parameters
  integer, intent(in) :: n
  real*8,intent(out) :: c(n,3)
  integer, intent(in) :: runit
  character(len=8000), intent(out) :: title
  ! internal variables
  integer :: i,j, io

  read(runit,'(A)', iostat=io) title

  do i=1,n
    read(runit,*, iostat=io) (c(i,j),j=1,3)
    if (io/=0) then
      write(*,*) 'Could not read 3-vector'
      write(*,*) 'routine=d3vecread(), n=',n,', unit=',runit
      write(*,*) 'title=',trim(title)
    endif
  enddo

endsubroutine

! =================================================================== !

!> reads a 3-vector c from unit runit
!> also reads the title, THIS MEANS THAT the current line of runit has to be the title line
subroutine z3vecread(n,c,runit,title)
  implicit none
  ! parameters
  integer, intent(in) :: n
  complex*16,intent(out) :: c(n,3)
  integer, intent(in) :: runit
  character(len=8000), intent(out) :: title
  ! internal variables
  integer :: i,j, io
  real*8 :: re_im(6)

  read(runit,'(A)', iostat=io) title

  do i=1,n
    read(runit,*) (re_im(j),j=1,6)
    if (io/=0) then
      write(*,*) 'Could not read 3-vector'
      write(*,*) 'routine=d3vecread(), n=',n,', unit=',runit
      write(*,*) 'title=',trim(title)
    endif
    do j=1,3
      c(i,j)=dcmplx(re_im(2*j-1),re_im(2*j))
    enddo
  enddo

endsubroutine

! =================================================================== !
! =================================================================== !
!                                Check functions                      !
! =================================================================== !
! =================================================================== !

logical function dishermitian(n,H) result(res)
! determines whether matrix H is symmetric
!
!
  implicit none

  integer,intent(in) :: n
  real*8,intent(in) :: H(n,n)

  integer :: i,j

  res=.true.
  do i=1,n
    do j=1,i
      if (.not.(H(i,j)==H(j,i)) ) res=.false.
    enddo
  enddo

endfunction

! =================================================================== !

logical function zishermitian(n,H) result(res)
! determines whether matrix H is hermitian
!
!
  implicit none

  integer,intent(in) :: n
  complex*16,intent(in) :: H(n,n)

  integer :: i,j

  res=.true.
  do i=1,n
    do j=1,i
      if (.not.(H(i,j)==conjg(H(j,i))) ) res=.false.
    enddo
  enddo

endfunction

! =================================================================== !

logical function disantihermitian(n,H) result(res)
! determines whether matrix H is antisymmetric
!
!
  implicit none

  integer,intent(in) :: n
  real*8,intent(in) :: H(n,n)

  integer :: i,j

  res=.true.
  do i=1,n
    do j=1,i
      if (.not.(H(i,j)==-H(j,i)) ) res=.false.
    enddo
  enddo

endfunction

! =================================================================== !

logical function zisantihermitian(n,H) result(res)
! determines whether matrix H is antihermitian
!
!
  implicit none

  integer,intent(in) :: n
  complex*16,intent(in) :: H(n,n)

  integer :: i,j

  res=.true.
  do i=1,n
    do j=1,i
      if (.not.(H(i,j)==-conjg(H(j,i))) ) res=.false.
    enddo
  enddo

endfunction

! =================================================================== !

logical function disunitary(n,H) result(res)
! determines whether matrix H is orthogonal
! calculates H^TH
! diagonal elements must be >1-1e-12, off-diagonal elements <1e-12
  implicit none

  integer,intent(in) :: n
  real*8,intent(in) :: H(n,n)

  integer :: i,j
  real*8 :: element_thres=1.d-12
  real*8 :: S(n,n)

  call dmultiply(n,H,H,S,'tn')

  res=.true.
  do i=1,n
    do j=1,i-1
      if ( dabs(S(i,j))>element_thres) res=.false.
    enddo
    if ( 1.d0-dabs(S(i,i))>element_thres) res=.false.
  enddo

endfunction

! =================================================================== !

logical function zisunitary(n,H) result(res)
! determines whether matrix H is orthogonal
! calculates H^TH
! diagonal elements must be >1-1e-12, off-diagonal elements <1e-12
  implicit none

  integer,intent(in) :: n
  complex*16,intent(in) :: H(n,n)

  integer :: i,j
  real*8 :: element_thres=1.d-12
  complex*16 :: S(n,n)

  call zmultiply(n,H,H,S,'tn')

  res=.true.
  do i=1,n
    do j=1,i-1
      if ( abs(S(i,j))>element_thres) res=.false.
    enddo
    if ( 1.d0-abs(S(i,i))>element_thres) res=.false.
  enddo

endfunction

! =================================================================== !

subroutine d3project_a_on_b(n, a, b, new_a)
  ! projects a on the direction of b and returns the new vector
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: a(n,3), b(n,3)
  real*8, intent(out) :: new_a(n,3)
  real*8 :: f

  f=sum(a*b)/sum(b*b)
  new_a=f*b

endsubroutine

! =================================================================== !

subroutine z3project_a_on_b(n, a, b, new_a)
  ! projects a on the direction of b and returns the new vector
  implicit none
  integer, intent(in) :: n
  complex*16, intent(in) :: a(n,3), b(n,3)
  complex*16, intent(out) :: new_a(n,3)
  complex*16 :: f

  f=sum(a*conjg(b))/sum(b*conjg(b))
  new_a=f*b

endsubroutine

endmodule