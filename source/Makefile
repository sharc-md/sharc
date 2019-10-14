#******************************************
#
#SHARC Program Suite
#
#    Copyright (c) 2019 University of Vienna
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
# Makefile for Sharc 
# Version 2.1
# 1.9.2019


# -------------------------------------- configure ---------------------------------------

# =======================================
#true, false     # also implies NETCDF
USE_PYSHARC := false

#intel, gnu
USE_COMPILER := intel

#mkl,gnu
USE_LIBS := mkl

#needed for PYSHARC
ANACONDA := /user/maximilian/anaconda2/
# =======================================




# -------------------------------------- settings ---------------------------------------

# path
SHARC_SOURCE := $(shell pwd)



# -------------------------------------- conditionals ---------------------------------------

# 
ifeq ($(USE_PYSHARC),true)
  FPREPROCESSOR := -D__PYSHARC__
  NETCDF_LIB = -L$(ANACONDA)/lib \
           -lhdf5 \
           -lhdf5_hl \
           -lnetcdf \
           -lmfhdf \
           -ldf \
           -ljpeg \
           -L../pysharc/lib \
           -lsharcnc
else
  FPREPROCESSOR := 
  NETCDF_LIB = 
endif


ifeq ($(USE_COMPILER), )
  USE_COMPILER:=gnu
endif
ifeq ($(USE_COMPILER), gnu)
  CC :=gcc
  F90 :=gfortran 
  EXTRAFLAGS := -ffree-line-length-none 
else ifeq ($(USE_COMPILER), intel)
  CC:=icc
  F90:=ifort 
  EXTRAFLAGS :=
endif


# Note: you probably need to adjust the search path if you haven't set $(MKLROOT) 
ifeq ($(USE_LIBS), mkl)
  LD= -L$(MKLROOT)/lib/intel64 -lmkl_rt -lpthread -lm -lgfortran $(NETCDF_LIB)
#   LD= -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -lgfortran $(NETCDF_LIB)
else ifeq ($(USE_LIBS), gnu)
  LD= -lblas -llapack -lfftw3 -lgfortran $(NETCDF_LIB)
endif


# -------------------------------------- settings ---------------------------------------



# include files
# FINCLUDE= -I$(SHARC_SOURCE)/
FINCLUDE= 
CINCLUDE= -I$(ANACONDA)/include

# flags
F90FLAGS = -O3 -fPIC -g  $(EXTRAFLAGS)

# debug
DEBUGFLAGS  = #-Wall -Wextra -Winline -Wno-unused-parameter -Wimplicit-interface -fbounds-check -fimplicit-none -fbacktrace -g 

# C flags
CLIB= -lgfortran
CPREPROCESSOR=
CSHARE= -fPIC
CDEBUG=
CFLAGS= -Wall -std=c99 -g -pedantic



# will be included in sharc source code so that compiler, host, etc. can be printed at runtime
DATE  = $(shell date)
HOST  = $(shell hostname)
COMPILEDIR = $(shell pwd)
COMPILER   = $(shell which ${F90})








# -------------------------------------- lists ---------------------------------------

# where the binaries should be put
EXEDIR     = ../bin

# sources for sharc.x, data_extractor.x and diagonalizer.x
TOOLS = definitions.o \
	matrix.o \
	string.o \
	input_list.o

SHARC  =  misc.o \
	  output.o \
	  restart.o \
	  qm_out.o \
	  nuclear.o \
	  electronic.o \
	  qm.o \
	  electronic_laser.o \
	  input.o \
	  decoherence_afssh.o
ifeq ($(USE_PYSHARC),true)
SHARC += interface.o
endif

ifeq ($(USE_PYSHARC),true)
SHARCTARGET := sharc2.x
else
SHARCTARGET := sharc1.x
endif


DATA_EXTRACTOR =  qm_out.o \
		  data_extractor.o

LASER=  LASER_definitions.o \
	LASER_calc_fftw.o \
	LASER_input.o \
	LASER_main.o 


DIAGONALIZER = string.o \
	      matrix.o \
	      diagonalizer.o


DATA_EXTRACTOR_NETCDF = data_extractor_NetCDF.o data_extractor_NetCDFmodule.o definitions_NetCDF.o qm_out.o 


DATA_CONVERTER = data_converter.o data_extractor_NetCDFmodule.o definitions_NetCDF.o qm_out.o 



# -------------------------------------- rules ---------------------------------------

ifeq ($(USE_PYSHARC),true)
all: sharc data_extractor diagonalizer laser sharcvars libsharc data_extractor_NetCDF data_converter
else
all: sharc data_extractor diagonalizer laser sharcvars
endif


sharc: build_info $(TOOLS) $(SHARC) main.o
	$(F90)    $(TOOLS) $(SHARC) main.o  -o  $(SHARCTARGET) $(LD)
	ln -sf $(SHARCTARGET) sharc.x


libsharc: build_info   $(TOOLS) $(SHARC)
	$(F90) -shared $(TOOLS) $(SHARC)  -o ../pysharc/lib/libsharc.so $(LD)


data_extractor: build_info $(DATA_EXTRACTOR)  $(TOOLS)
	$(F90)  $(DATA_EXTRACTOR) $(TOOLS) -o $@.x $(LD)


data_extractor_NetCDF: $(DATA_EXTRACTOR_NETCDF) $(TOOLS)
	$(F90)  $(DATA_EXTRACTOR_NETCDF) $(TOOLS) -o $@.x $(LD)


data_converter: $(DATA_CONVERTER) $(TOOLS)
	$(F90)  $(DATA_CONVERTER) $(TOOLS) -o $@.x $(LD)


diagonalizer: $(DIAGONALIZER) 
	$(F90)  $(DIAGONALIZER) -o $@.x $(LD)


laser:  $(LASER) 
	$(F90)   $(LASER) -o $@.x $(LD)


install: all
	cp -f *.x $(EXEDIR)
	cp -f sharcvars.sh $(EXEDIR)
	cp -f sharcvars.csh $(EXEDIR)


# write sharcvars.sh file
ifeq ($(USE_PYSHARC),true)
sharcvars:
	@echo 'export SHARC=$(SHARC_SOURCE)/../bin' > sharcvars.sh
	@echo 'export PYSHARC=$(SHARC_SOURCE)/../pysharc' >> sharcvars.sh
	@echo 'export SHARCLIB=$(SHARC_SOURCE)/../lib' >> sharcvars.sh
	@echo 'export ANACONDA=$(ANACONDA)' >> sharcvars.sh
	@echo 'export PYTHONPATH=$$SHARCLIB:$$PYSHARC:$$PYTHONPATH' >> sharcvars.sh
	@echo 'export LD_LIBRARY_PATH=$$SHARCLIB:$$ANACONDA/lib:$$LD_LIBRARY_PATH' >> sharcvars.sh
	@echo 'setenv SHARC=$(SHARC_SOURCE)/../bin' > sharcvars.csh
	@echo 'setenv PYSHARC=$(SHARC_SOURCE)/../pysharc' >> sharcvars.csh
	@echo 'setenv SHARCLIB=$(SHARC_SOURCE)/../lib' >> sharcvars.csh
	@echo 'setenv ANACONDA=$(ANACONDA)' >> sharcvars.csh
	@echo 'setenv PYTHONPATH=$$SHARCLIB:$$PYSHARC:$$PYTHONPATH' >> sharcvars.csh
	@echo 'setenv LD_LIBRARY_PATH=$$SHARCLIB:$$ANACONDA/lib:$$LD_LIBRARY_PATH' >> sharcvars.csh
else
sharcvars:
	@echo 'export SHARC=$(SHARC_SOURCE)/../bin' > sharcvars.sh
	@echo 'export SHARCLIB=$(SHARC_SOURCE)/../lib' >> sharcvars.sh
	@echo 'export PYTHONPATH=$$SHARCLIB:$$PYTHONPATH' >> sharcvars.sh
	@echo 'export LD_LIBRARY_PATH=$$SHARCLIB:$$LD_LIBRARY_PATH' >> sharcvars.sh
	@echo 'setenv SHARC=$(SHARC_SOURCE)/../bin' > sharcvars.csh
	@echo 'setenv SHARCLIB=$(SHARC_SOURCE)/../lib' >> sharcvars.csh
	@echo 'setenv PYTHONPATH=$$SHARCLIB:$$PYTHONPATH' >> sharcvars.csh
	@echo 'setenv LD_LIBRARY_PATH=$$SHARCLIB:$$LD_LIBRARY_PATH' >> sharcvars.csh
endif


# do not call this target build_info.inc
build_info: 
	@echo "Updating \"build_info.inc\""
	@echo "character(len=500), parameter :: build_date=\"$(DATE)\"" > build_info.inc
	@echo "character(len=500), parameter :: build_host=\"$(HOST)\"" >> build_info.inc
	@echo "character(len=500), parameter :: build_dir=\"$(COMPILEDIR)\"" >> build_info.inc
	@echo "character(len=500), parameter :: build_compiler=\"$(COMPILER)\"" >> build_info.inc
	@echo "character(len=500), parameter :: use_pysharc=\"$(USE_PYSHARC)\"" >> build_info.inc


# clean function
clean:
	rm -rvf *.o
	rm -rvf *.x
	rm -rvf *.mod
	rm -rvf build_info.inc
	rm -rvf *.so
	rm -rvf ../pysharc/lib/libsharc.so
	rm -rvf sharcvars.sh

test_path:
	echo $(MAKEFILEPATH)

%.o: %.f90 
	$(F90) $(FINCLUDE) $(DEBUGFLAGS) $(F90FLAGS) -c $<

%.o: %.F90 
	$(F90) $(FINCLUDE) $(FPREPROCESSOR) $(DEBUGFLAGS) $(F90FLAGS) -c $<

%.o: %.c
	$(CC) $(CINCLUDE) $(CPREPROCESSOR) $(CFLAGS) $(CSHARE) $(CDEBUG) -c $<



# -------------------------------------- dependencies ---------------------------------------


# the following lines define the dependencies of the Fortran90 module files
ifeq ($(USE_PYSHARC),true)
main.o: output.o qm.o restart.o misc.o electronic_laser.o input.o decoherence_afssh.o interface.o
else
main.o: output.o qm.o restart.o misc.o electronic_laser.o input.o decoherence_afssh.o
endif

interface.o: electronic.o definitions.o definitions_NetCDF.o output.o qm.o restart.o misc.o electronic_laser.o input.o decoherence_afssh.o
misc.o: definitions.o
output.o: $(TOOLS) build_info misc.o
nuclear.o: $(TOOLS)
electronic.o: nuclear.o $(TOOLS) decoherence_afssh.o
electronic_laser.o: electronic.o $(TOOLS)
input.o: $(TOOLS) output.o restart.o misc.o
qm.o: $(TOOLS) electronic.o qm_out.o restart.o
qm_out.o: $(TOOLS)
restart.o: $(TOOLS) misc.o decoherence_afssh.o
decoherence_afssh.o: $(TOOLS) definitions.o 

LASER_input.o: LASER_definitions.o 
LASER_calc_fftw.o: LASER_definitions.o

data_extractor_NetCDFmodule.o: $(TOOLS) qm_out.o
data_extractor_NetCDF.o: data_extractor_NetCDFmodule.o definitions.o matrix.o definitions_NetCDF.o
data_converter.o:        data_extractor_NetCDFmodule.o definitions.o matrix.o definitions_NetCDF.o

