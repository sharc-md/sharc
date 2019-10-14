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



#include <stdlib.h>
#include <stdio.h>
#include "nc_basic.h"

inline
void
prt_error(int istatus)
{
    if (istatus != NC_NOERR) {
             fprintf(stderr, "%s\n", nc_strerror(istatus)); 
             exit(-1);
    }

};

inline
int
create_ncfile(const char *str, int mode)
{
    int iretval;
    int ncid = 0;
    if ((iretval = nc_create(str, mode, &ncid)))
        prt_error(iretval);
    return ncid;
};

inline
int
open_ncfile(const char *str, int mode)
{
    int iretval;
    int ncid = 0;
    if ((iretval = nc_open(str, mode, &ncid)))
        prt_error(iretval);
    return ncid;
};

inline
void 
close_ncfile_(const int* ncid)
{
    int iretval;
    if ((iretval = nc_close(*ncid)))
        prt_error(iretval);
};


inline
int
define_dimension(const int ncid, const char *name,const size_t N)
{

    int iretval;
    int id = 0;
    if ((iretval = nc_def_dim(ncid, name, N, &id)))
        prt_error(iretval);
    return id;
};
