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





all: build_ext




build_ext: lib/libsharc.so
	python sharc_setup build_ext --build-lib .

install: lib/libsharc.so lib/libsharcnc.so build_ext
	cp lib/*.so ../lib
	cp bin/*.py ../bin

clean:
	python sharc_setup clean --all
	rm lib/libsharc.so lib/libsharcnc.so
	rm sharc/sharc*.so
	rm netcdf/*.o

lib/libsharc.so: lib/libsharcnc.so
	cd ../source; make libsharc

lib/libsharcnc.so:
	cd netcdf; make lib



