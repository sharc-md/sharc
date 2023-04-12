#!/usr/bin/perl -w
#
#******************************************
#
#    SHARC Program Suite
#
#    Copyright (c) 2023 University of Vienna
#
#    This file is part of SHARC.
#
#    SHARC is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    SHARC is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    inside the SHARC manual.  If not, see <http://www.gnu.org/licenses/>.
#
#******************************************

use strict;

my($type,$ndim,$d,$maxdim);

$maxdim=7;

print STDOUT "module my_alloc\n";
print STDOUT "  use sysparam\n";
print STDOUT "  use memlog\n";
print STDOUT "  implicit none\n";
print STDOUT "  interface myalloc\n";
print STDOUT "    module procedure &\n    & ";
$type='il';
print STDOUT " a${type}1d";
for($ndim=2;$ndim<=$maxdim;$ndim++)
{
  print STDOUT ",a${type}${ndim}d";
  if((${ndim}%20 eq 0) and (${ndim} lt ${maxdim})){print STDOUT " &\n    & ";}
}
print STDOUT " &\n    & ";
foreach $type qw(is r d c z l)
{
  for($ndim=1;$ndim<=$maxdim;$ndim++)
  {
    print STDOUT ",a${type}${ndim}d";
    if((${ndim}%20 eq 0) and (${ndim} lt ${maxdim})){print STDOUT " &\n    & ";}
  }
  if($type ne "l"){ print STDOUT " &\n    & ";}
}
print STDOUT "\n  end interface\n";
print STDOUT "  interface mydealloc\n";
print STDOUT "    module procedure &\n    & ";
$type='il';
print STDOUT " d${type}1d";
for($ndim=2;$ndim<=$maxdim;$ndim++)
{
  print STDOUT ",d${type}${ndim}d";
  if((${ndim}%20 eq 0) and (${ndim} lt ${maxdim})){print STDOUT " &\n    & ";}
}
print STDOUT " &\n    & ";
foreach $type qw(is r d c z l)
{
  for($ndim=1;$ndim<=$maxdim;$ndim++)
  {
    print STDOUT ",d${type}${ndim}d";
    if((${ndim}%20 eq 0) and (${ndim} lt ${maxdim})){print STDOUT " &\n    & ";}
  }
  if($type ne "l"){ print STDOUT " &\n    & ";}
}
print STDOUT "\n  end interface\n";
print STDOUT "contains\n\n!===============================================\n\n";
foreach $type qw(il is r d c z l)
{
  for($ndim=1;$ndim<=$maxdim;$ndim++)
  {
    print STDOUT "  integer function a${type}${ndim}d(arr";
    for($d=1;$d<=$ndim;$d++)
    {
      print STDOUT ",d${d}";
      if((${d}%15 eq 0) and (${d} lt ${ndim})){print STDOUT " &\n      & ";}
    }
    print STDOUT ",msg,thrd)\n";
    print STDOUT "    implicit none\n";
    if($type eq "il") {print STDOUT "    integer(kind=ilong)";}
    elsif($type eq "is") {print STDOUT "    integer(kind=ishort)";}
    elsif($type eq "r") {print STDOUT "    real(kind=sip)";}
    elsif($type eq "d") {print STDOUT "    real(kind=dop)";}
    elsif($type eq "c") {print STDOUT "    complex(kind=sip)";}
    elsif($type eq "z") {print STDOUT "    complex(kind=dop)";}
    elsif($type eq "l") {print STDOUT "    logical";}
    print STDOUT ",dimension(:";
    for($d=2;$d<=$ndim;$d++)
    {
      print STDOUT ",:";
      if((${d}%20 eq 0) and (${d} lt ${ndim})){print STDOUT " &\n      & ";}
    }
    print STDOUT "),allocatable :: arr\n";
    print STDOUT "    integer(kind=ilong), intent(in) :: &\n      &  d1"; 
    for($d=2;$d<=$ndim;$d++)
    {
      print STDOUT ",d${d}";
      if((${d}%10 eq 0) and (${d} lt ${ndim})){print STDOUT " &\n      & ";}
    }
    print STDOUT "\n";
    print STDOUT "    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg\n";
    print STDOUT "    INTEGER, INTENT(IN), OPTIONAL :: thrd\n";
    #print STDOUT "    INTEGER :: a${type}${ndim}d\n"; # MBR
    print STDOUT "    CHARACTER*30 :: nm\n\n";
    print STDOUT "    integer(kind=ilong) :: sizeofone\n";
    if($type eq "il") {print STDOUT "    sizeofone=sizeof(1_ilong)\n";}
    elsif($type eq "is") {print STDOUT "    sizeofone=sizeof(1_ishort)\n";}
    elsif($type eq "r") {print STDOUT "    sizeofone=sizeof(1.0_sip)\n";}
    elsif($type eq "d") {print STDOUT "    sizeofone=sizeof(1.0_dop)\n";}
    elsif($type eq "c") {print STDOUT "    sizeofone=2*sizeof(1.0_sip)\n";}
    elsif($type eq "z") {print STDOUT "    sizeofone=2*sizeof(1.0_dop)\n";}
    elsif($type eq "l") {print STDOUT "    sizeofone=sizeof(.true.)\n";}
    print STDOUT "    nm=\'+\'\n";
    print STDOUT "    IF(present(msg)) THEN\n      nm=msg\n    END IF\n";
    print STDOUT "    if(allocated(arr))then\n";
    print STDOUT "      a${type}${ndim}d=-1\n";
    print STDOUT "      write(6,*)\'array already allocated a${type}${ndim}d (\',trim(adjustl(nm)),\')\'\n";
    print STDOUT "      write(0,*)\'array already allocated a${type}${ndim}d (\',trim(adjustl(nm)),\')\'\n";
    print STDOUT "    else\n";
    print STDOUT "      memu=sizeofone";
    for($d=1;$d<=$ndim;$d++){print STDOUT "*d${d}";}
    print STDOUT "\n";
    print STDOUT "      if((maxmem.gt.0).and.((mem_used+memu).gt.maxmem)) then\n";
    print STDOUT "        a${type}${ndim}d=-1\n";
    print STDOUT "        write(6,111)maxmem/(1024*1024), mem_used/(1024*1024), memu/(1024*1024)\n";
#    print STDOUT "        write(0,111)maxmem, mem_used, memu\n";
    print STDOUT '111 FORMAT("Exceeding memory limit. (maxmem, mem_used, memu) =", "(", I10, ",", I10, ",", I10, ") MB")'."\n";
    print STDOUT "      else\n";
    print STDOUT "        allocate(arr(1:d1"; 
    for($d=2;$d<=$ndim;$d++)
    {
      print STDOUT ",1:d${d}";
      if((${d}%10 eq 0) and (${d} lt ${ndim})){print STDOUT " &\n        & ";}
    }
    print STDOUT "),STAT=a${type}${ndim}d)\n";
    print STDOUT "        if(a${type}${ndim}d.ne.0)then\n";
    print STDOUT "          write(6,*)\'error in allocation routine a${type}${ndim}d\'\n";
    print STDOUT "          write(0,*)\'error in allocation routine a${type}${ndim}d\'\n";
    print STDOUT "        else\n";
    # print STDOUT "          memu=sizeof(arr)\n";
    print STDOUT "          mem_used=mem_used+memu\n";
    print STDOUT "          if(present(thrd))then\n";
    print STDOUT "            call log_memory(nm,thrd)\n";
    print STDOUT "          else\n";
    print STDOUT "            if(present(msg))then\n";
    print STDOUT "              call log_memory(nm)\n";
    print STDOUT "            end if\n";
    print STDOUT "          end if\n";
    print STDOUT "        end if\n";
    print STDOUT "      end if\n";
    print STDOUT "    end if\n";
    print STDOUT "  end function a${type}${ndim}d\n";
    print STDOUT "\n";
    print STDOUT "  !---------------------------------\n\n";
    print STDOUT "  integer function d${type}${ndim}d(arr,msg,thrd)\n";
    print STDOUT "    implicit none\n";
    if($type eq "il") {print STDOUT "    integer(kind=ilong)";}
    elsif($type eq "is") {print STDOUT "    integer(kind=ishort)";}
    elsif($type eq "r") {print STDOUT "    real(kind=sip)";}
    elsif($type eq "d") {print STDOUT "    real(kind=dop)";}
    elsif($type eq "c") {print STDOUT "    complex(kind=sip)";}
    elsif($type eq "z") {print STDOUT "    complex(kind=dop)";}
    elsif($type eq "l") {print STDOUT "    logical";}
    print STDOUT ",dimension(:";
    for($d=2;$d<=$ndim;$d++)
    {
      print STDOUT ",:";
      if((${d}%20 eq 0) and (${d} lt ${ndim})){print STDOUT " &\n      & ";}
    }
    print STDOUT "),allocatable :: arr\n";
    print STDOUT "    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg\n";
    print STDOUT "    INTEGER, INTENT(IN), OPTIONAL :: thrd\n";
    #print STDOUT "    INTEGER :: d${type}${ndim}d\n"; # MBR
    print STDOUT "    CHARACTER*30 :: nm\n\n";
    print STDOUT "    nm=\'-\'\n";
    print STDOUT "    IF(present(msg)) THEN\n      nm=msg\n    END IF\n";
    print STDOUT "    if(.not.allocated(arr))then\n";
    print STDOUT "      d${type}${ndim}d=-1\n";
    print STDOUT "      write(6,*)\'array not allocated d${type}${ndim}d (\',trim(adjustl(nm)),\')\'\n";
    print STDOUT "      write(0,*)\'array not allocated d${type}${ndim}d (\',trim(adjustl(nm)),\')\'\n";
    print STDOUT "    else\n";
    print STDOUT "      memu=sizeof(arr)\n";
    print STDOUT "      deallocate(arr,STAT=d${type}${ndim}d)\n";
    print STDOUT "      if(d${type}${ndim}d.ne.0)then\n";
    print STDOUT "        write(6,*)\'error in deallocation routine d${type}${ndim}d\'\n";
    print STDOUT "        write(0,*)\'error in deallocation routine d${type}${ndim}d\'\n";
    print STDOUT "      else\n";
    print STDOUT "        mem_used=mem_used-memu\n";
    print STDOUT "        if(present(thrd))then\n";
    print STDOUT "          call log_memory(nm,thrd)\n";
    print STDOUT "        else\n";
    print STDOUT "          if(present(msg))then\n";
    print STDOUT "            call log_memory(nm)\n";
    print STDOUT "          end if\n";
    print STDOUT "        end if\n";
    print STDOUT "      end if\n";
    print STDOUT "    end if\n";
    print STDOUT "  end function d${type}${ndim}d\n";
    print STDOUT "\n";
    print STDOUT "  !---------------------------------\n\n";
  }
}
print STDOUT "end module\n";
