//******************************************
//
//    SHARC Program Suite
//
//    Copyright (c) 2019 University of Vienna
//
//    This file is part of SHARC.
//
//    SHARC is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    SHARC is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    inside the SHARC manual.  If not, see <http://www.gnu.org/licenses/>.
//
//******************************************


#ifndef NETCDEF_ERROR_H_
#define NETCDEF_ERROR_H_

#include <netcdf.h>

// check_nccall macro
#define check_nccall(iret, ncall) \
    if (( iret = ncall)) \
        prt_error(iret);


#ifdef __cplusplus
extern"C" {
#endif

void prt_error(int istatus);

int open_ncfile(const char *str, int imode);
int create_ncfile(const char *str, int imode);
void close_ncfile_(const int* ncid);
int define_dimension(const int ncid, const char *name, const size_t N);


#ifdef __cplusplus
}
#endif

// end code
#endif
