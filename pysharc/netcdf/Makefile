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







# -------------------------------------- configure ---------------------------------------

# =======================================
ANACONDA := ${CONDA_PREFIX}
# =======================================


INC= -I.  -I$(ANACONDA)/include
OPT= -O3 -Wall
SHARED= -fPIC

CC=gcc


OBJ_LIB=nc_basic.o nc_sharc.o
library=libsharcnc.so




all: lib

lib: $(OBJ_LIB)
	[ -d ../lib ] || mkdir ../lib
	ar -cvq ../lib/$(library) $(OBJ_LIB)

clean:
	rm ../lib/$(library)
	rm *.o

%.o: %.c
	$(CC) $(SHARED) $(INC) $(OPT) -c $<

