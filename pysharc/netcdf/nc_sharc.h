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


#ifndef NC_SHARC_H_
#define NC_SHARC_H_

#include "nc_basic.h"

#ifdef __cplusplus
extern "C" {
#endif

struct sharc_ncxyz {
    // File ID
    int id;
    // Variable IDs 
    int ian_id;
    int crd_id;
//    int veloc_id;
}; 

struct sharc_ncdat {
    // File ID
    int id;
    // Variable IDs 
    int energy_id;
}; 

struct sharc_ncoutput {
    // File ID
    int id;
    // Variable IDs 
    int H_MCH_id;
    int U_id;
    int DM_id;
    int overlaps_id;
    int coeff_diag_id;
    // other
    int e_id;
    int hopprop_id;
    int crd_id;
    int veloc_id;
    // one dimensional
    int randnum_id;
    int state_diag_id;
    int state_MCH_id;
    int time_step_id;
}; 


void write_sharc_ncxyz_traj_(const int* istep, const int* NAtoms, double* Crd, struct sharc_ncxyz* ncxyz);
void write_sharc_ncxyz_init_(const int* NAtoms, int* IAn, double* Crd, struct sharc_ncxyz* ncxyz);

void write_sharc_ncdat_init_(const double* E, struct sharc_ncdat* ncdat);
void write_sharc_ncdat_traj_(const int* istep, const double* E, struct sharc_ncdat* ncdat);


void write_sharc_ncoutputdat_istep_(
        // 
        const int* istep,
        const int* natoms,
        const int* nstates,
        // Multidimensional
        const double* H_MCH_ss,           // complex, (frame, 2*nstates, nstates)
        const double* U_ss,               // complex, (frame, 2*nstates, nstates)
        const double* DM_print_ssd,       // complex, (frame, 2*nstates, nstates, 3)
        const double* overlaps_ss,        // complex, (frame, 2*nstates, nstates)
        const double* coeff_diag_s,       // complex, (frame, 2*nstates)
        const double* E,                  // real, contains Etot, Epot and Ekin, (frame, 3)
        const double* hopprob_s,          // real, (frame, nstates)
        const double* geom_ad,            // real, (frame, natoms, 3)
        const double* veloc_ad,           // real, (frame, natoms, 3)
        // One dimensional opjects
        const double* randnum,            // real, (frame, 1)
        const int* state_diag,            // int, (frame, 1)
        const int* state_MCH,             // int, (frame, 1)
        const int* time_step,             // int, (frame, 1)
        //
        struct sharc_ncoutput* ncdat
);

#ifdef __cplusplus
}
#endif

#endif
