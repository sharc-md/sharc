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
# Makefile for WFoverlap 
# SHARC Version 2.0
# Feb-01-2018
#

############## Compiler definition ##############



# **** ifort ****

# #STATIC = -static-intel -qopenmp-link=static
# PROFILE = # -pg
# OMP = -qopenmp
# FC = ifort
# DEBUG = #-g #-warn all # -traceback -check bounds
# OPT = -O3 -ipo
# FCFLAGS = $(OPT) $(OMP) $(PROFILE) $(DEBUG) -fpp -i8 -DEXTBLAS   -I"${MKLROOT}/include" -z muldefs
# LINKFLAGS = $(STATIC) $(PROFILE)
# 
# # openmp (multithreaded) compilation
# #LALIB =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm -ldl
# #LALIB  = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lm
# LALIB  =  -Wl,--start-group ${MKLROOT}/lib/libmkl_intel_ilp64.a ${MKLROOT}/lib/libmkl_intel_thread.a ${MKLROOT}/lib/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl



# **** gfortran ****

FC = gfortran
FCFLAGS = -O0 -cpp -g3 -fdefault-integer-8 -Wall -fbacktrace -DEXTBLAS
LINKFLAGS = 
LALIB = -lblas -llapack -fopenmp

############## Main objects and libraries #######

MAINOBS=sysparam.o memlog.o my_alloc.o global_storage.o read_columbus.o read_lumorb.o read_turbomole.o dgetrf.o lowdin.o sortblockmod.o calcmod.o inputmod.o outputmod.o
DUMMYOBS=read_dalton_dummy.o read_molcas_dummy.o


############## Dalton and SEWARD ################
# no Dalton and SEWARD support
#LIBS = $(LALIB)
#OPTOBS = $(DUMMYOBS)

# activate direct reading of Dalton and Seward files
LIBS = $(COLUMBUS)/libmolcas_col.a $(COLUMBUS)/colib.a $(COLUMBUS)/blaswrapper.a  $(LALIB)
OPTOBS = read_dalton.o read_molcas.o

############## Compilation routines #############

.PHONY: all clean test alloc

%.o: %.f90
	$(FC) $(FCFLAGS) $(PROFILE) -c $<
%.o: %.f
	$(FC) $(FCFLAGS) $(PROFILE) -c $<

wfoverlap.x : main.f90 $(MAINOBS) $(OPTOBS) iomod.o
	$(FC) $(FCFLAGS) $(LINKFLAGS) $(PROFILE) $^ -o $@ $(LIBS)
	cp $@ ../../bin

# executable that is linked to dummy I/O interfaces and that reads only ASCII files
wfoverlap_ascii.x : main.f90 $(MAINOBS) $(DUMMYOBS) iomod.o
	$(FC) $(FCFLAGS) $(LINKFLAGS) $(PROFILE) $^ -o $@ $(LALIB)
	cp $@ ../../bin
	ln -fs $@ ../../bin/wfoverlap.x

alloc:
	./write_allocmod.pl > my_alloc.f90

all: wfoverlap_ascii.x wfoverlap.x

clean:
	rm -f *.o *.mod *.x *_genmod.f90

test: 
	cd ../; ./scripts/ovl_test.bash $(shell pwd)/.. | tee test.out

