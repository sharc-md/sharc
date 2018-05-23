#******************************************
#
#    SHARC Program Suite
#
#    Copyright (c) 2018 University of Vienna
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

#******************************************
# Makefile for Sharc 
# Version 2.0
# Feb-01-2018
#
# uncomment the following line and comment the 2nd next line to 
# use the intel fortran compiler.
F90 = ifort 
# F90 = gfortran


# uncomment the following line and comment the 2nd next line to use intel's MKL
# Note: you probably need to adjust the search path if you haven't set
# $(MKLROOT) 
LD =  -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
# LD = -lblas -llapack -lfftw3


# works for GFORTRAN and IFORT
F90FLAGS     = -O2
DEBUGFLAGS   = 
#DEBUGFLAGS   = -Wall -Wextra -Winline -Wno-unused-parameter -Wimplicit-interface -fbounds-check -fimplicit-none -fbacktrace -g 


# will be included in sharc source code so that compiler, host, etc. can be printed at runtime
DATE  = $(shell date)
HOST  = $(shell hostname)
COMPILEDIR = $(shell pwd)
COMPILER   = $(shell which ${F90})

# where the binaries should be put
EXEDIR     = ../bin

# sources for sharc.x, data_extractor.x and diagonalizer.x
TOOLS = definitions.o \
        matrix.o \
 	string.o \
	input_list.o

SHARC  =  main.o \
	  misc.o \
	  output.o \
	  restart.o \
   	  qm_out.o \
       	  nuclear.o \
          electronic.o \
          qm.o \
          electronic_laser.o \
	  input.o \
          decoherence_afssh.o 



DATA_EXTRACTOR =      qm_out.o \
  	              data_extractor.o

LASER=  LASER_definitions.o \
	LASER_calc_fftw.o \
        LASER_input.o \
        LASER_main.o 


DIAGONALIZER=string.o \
             matrix.o \
             diagonalizer.o


all: sharc data_extractor diagonalizer laser

sharc: build_info $(TOOLS) $(SHARC)  
	$(F90)  $(SHARC) $(TOOLS) -o  $@.x $(LD)

data_extractor: build_info $(DATA_EXTRACTOR) 
	$(F90)  $(DATA_EXTRACTOR) $(TOOLS) -o $@.x $(LD)

diagonalizer: $(DIAGONALIZER) 
	$(F90)  $(DIAGONALIZER) -o $@.x $(LD)


laser:  $(LASER) 
	$(F90)   $(LASER) -o $@.x $(LD)

clean:
	rm -rvf *.o
	rm -rvf *.x
	rm -rvf *.mod
	rm -rvf build_info.inc


install: 
	cp -f *.x $(EXEDIR)	

# the following lines define the dependencies of the Fortran90 module files
main.o: output.o qm.o restart.o misc.o electronic_laser.o input.o decoherence_afssh.o
misc.o:	definitions.o
output.o: $(TOOLS) build_info misc.o
nuclear.o: $(TOOLS)
electronic.o: nuclear.o $(TOOLS) decoherence_afssh.o
electronic_laser.o: electronic.o $(TOOLS)
input.o: $(TOOLS) output.o restart.o misc.o
qm.o: $(TOOLS) electronic.o qm_out.o restart.o
qm_out.o: $(TOOLS)
restart.o: $(TOOLS) misc.o
decoherence_afssh.o: $(TOOLS) definitions.o 

LASER_input.o: LASER_definitions.o 
LASER_calc_fftw.o: LASER_definitions.o



# do not call this target build_info.inc, because the file should always be updated when compiling
build_info: 
	@echo "Updating \"build_info.inc\""
	@echo "character*8000, parameter :: build_date=\"$(DATE)\"" > build_info.inc
	@echo "character*8000, parameter :: build_host=\"$(HOST)\"" >> build_info.inc
	@echo "character*8000, parameter :: build_dir=\"$(COMPILEDIR)\"" >> build_info.inc
	@echo "character*8000, parameter :: build_compiler=\"$(COMPILER)\"" >> build_info.inc

%.o: %.f90 
	$(F90) $(DEBUGFLAGS) $(F90FLAGS) -c $<

