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

!> # Module INPUT_LIST
!>
!> \authors Sebastian Mai, Philipp Marquetand
!> \date 03.03.2017
!>
!> This module defines a variable-size list of pairs of strings.
!> Features:
!> - a (key,value) pair can be added with add_key
!> - pairs cannot be deleted
!> - the list is automatically allocated when the first pair is added
!> - the module takes care of reallocating the list if too many keys are added
!> - key and value can be obtained based on index
!> - values can be obtained based on keys
!> - it can be checked whether a given key is in the list
module input_list

private alloc_list, resize_list

  character*255,allocatable, private :: list(:,:)       !< the list
  integer, private :: ncurr=-1, nalloc=-1               !< position counters
  integer,parameter, private :: nstep=10                !< always allocate a multiple of nstep elements

  contains

! =========================================================== !

!> used to initialize the list (allocate with minimum number of elements)
!> \param n allocate at least n elements
  subroutine read_input_list_from_file(nunit)
  
    use string
    
    implicit none
    
    integer, intent(in) :: nunit
    character*8000 :: keyword, value
    integer :: io
    
    ! read all keywords into the variable size list and copy to output
    do 
      io=0 !no idea why something has to be here but without it the next line is not executed correctly
      ! reads the next non-comment line from input file
      call next_keyword(nunit, keyword, value, io)
      if (io/=0) exit
      if (keyword(1:4) == '****') then
!         write(*,*) 'Found ****, finished reading keywords.' 
        exit
      endif
      ! keyword,value pair is added to input list (in input_list)
      call add_key(keyword,value)
    enddo


  endsubroutine
  
! =================================================================== !

!> this routine returns the next valid (key,value) pair from unit nunit
!> removes comments, skips empty/comment-only lines, makes lowercase, adjusts to left
!> and splits the line into key and value
  subroutine next_keyword(nunit, keyword, values, stat)
    use string
    implicit none
    integer, intent(in) :: nunit
    character*8000, intent(out) :: values
    character*8000, intent(out) :: keyword
    integer, intent(out) :: stat
    character*8000 :: line, line_precomment, line_postcomment
    integer :: io

    do
      ! read a line
      read(nunit,'(A)', iostat=io) line
      if (io<0) then
        stat=-1
        return
      endif

      ! remove comment, spaces, make lowercase and split the string at spaces
      call cut(line, '#', line_precomment, line_postcomment)
      if (trim(line_precomment)=='') cycle
      call lowercase(line_precomment)
      call compact(line_precomment,' ')
      call cut(line_precomment, ' ', line, line_postcomment)
      line_postcomment=adjustl(line_postcomment)

      keyword=line
      values=line_postcomment
      return
    enddo

  endsubroutine


! =========================================================== !

!> used to initialize the list (allocate with minimum number of elements)
!> \param n allocate at least n elements
  subroutine alloc_list(n)
    implicit none
    integer, intent(in) :: n

    nalloc=(n/nstep+1)*nstep
!       write(0,*) 'allocate',nalloc
    allocate(list(nalloc,2))
    ncurr=n
  endsubroutine

! ============================================= !

!> copies the list to a temporary array, deallocates and allocates the list
!> with a larger number of elements and copies back the list
!> \param n allocate at least so many elements
  subroutine resize_list(n)
    implicit none
    integer, intent(in) :: n
    character*255 :: temp(min(n,ncurr),2)

    if (n<=nalloc) then
      ncurr=n
    else
      nalloc=(n/nstep+1)*nstep
!         write(0,*) 'allocate',nalloc
      temp=list
      deallocate(list)
      allocate(list(nalloc,2))
      list(1:min(n,ncurr),:)=temp
      ncurr=n
    endif
  endsubroutine

! ============================================= !

!> add a (key,value) pair to the list at position ncurr
  subroutine add_key(key,value)
    implicit none
    character(*),intent(in) :: key,value
    integer :: nold, n=1

    if (ncurr==-1) then
      call alloc_list(n)
      ncurr=n
    else
      nold=ncurr
      call resize_list(ncurr+n)
    endif
    list(ncurr,1)=key
    list(ncurr,2)=value
  endsubroutine

! ============================================= !

!> return the key of element number i
!> \param io is -1 if there is no element i, and 0 otherwise
  character*255 function get_key(i,io)
    implicit none
    integer,intent(in) :: i
    integer,intent(out) :: io

    if (i>ncurr) then
      io=-1
      get_key=''
    else
      io=0
      get_key=list(i,1)
    endif
  endfunction

! ============================================= !

!> return the value of element number i
!> \param io is -1 if there is no element i, and 0 otherwise
  character*255 function get_value(i,io)
    implicit none
    integer,intent(in) :: i
    integer,intent(out) :: io

    if (i>ncurr) then
      io=-1
      get_value=''
    else
      io=0
      get_value=list(i,2)
    endif
  endfunction

! ============================================= !

!> return the value of the first element with the given key
!> returns an empty string if key was not found
!> \param io is -1 if the key was not found, 0 otherwise
  character*255 function get_value_from_key(key,io)
    implicit none
    character(len=*),intent(in) :: key
    integer,intent(out) :: io
    integer :: i

    do i=1,ncurr
      if (trim(key)==trim(list(i,1))) then
        io=0
        get_value_from_key=list(i,2)
        return
      endif
    enddo
    io=-1
    get_value_from_key=''
  endfunction

! ============================================= !

!> is .true. if the key is in the list
!> .false. otherwise
  logical function key_in_list(key)
    implicit none
    character*255,intent(in) :: key
    integer :: i

    key_in_list=.false.
    do i=1,ncurr
      if (key==list(i,1)) then
        key_in_list=.true.
        return
      endif
    enddo
  endfunction

! ============================================= !

!> get the current number of elements
  integer function get_ncurr()
    implicit none

    get_ncurr=ncurr
  endfunction

endmodule
