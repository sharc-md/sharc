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


#include "nc_sharc.h"
#include "stdio.h"

#define INATOMS 1
#define ISPATIAL 2

void
write_sharc_ncoutputdat_init_()
{
};

// --------------------------------------------------------------------------

void 
setup_ncoutputdat(int natoms, int nstates, struct sharc_ncoutput* ncdat)
{
    /*
     * defines all variables used in ncoutput
     */
    
    // error handler
    int iret = 0;
    // open sharc file
    ncdat->id      = create_ncfile("output.dat.nc", NC_CLOBBER);
    // define dimensions
    int ONE_ID = define_dimension(ncdat->id, "1", 1);
    int THREE_ID = define_dimension(ncdat->id, "3", 3);
    int NATOMS_ID = define_dimension(ncdat->id, "natoms", natoms);
    int NSTATES_ID = define_dimension(ncdat->id, "nstates", nstates);
    int CNSTATES_ID = define_dimension(ncdat->id, "cnstates", 2*nstates);
    int FRAME_ID = define_dimension(ncdat->id, "frame", NC_UNLIMITED);
    // define dimensions
    int dimids[4];
    dimids[0] = FRAME_ID;
    dimids[1] = CNSTATES_ID;
    dimids[2] = NSTATES_ID;
    dimids[3] = THREE_ID;
    // define nc variable ids
    check_nccall(iret, 
         nc_def_var(ncdat->id, "H_MCH", NC_DOUBLE, 3, &dimids[0], &ncdat->H_MCH_id)
    );
    check_nccall(iret, 
         nc_def_var(ncdat->id, "U", NC_DOUBLE, 3, &dimids[0], &ncdat->U_id)
    );
    check_nccall(iret, 
         nc_def_var(ncdat->id, "Ovlap", NC_DOUBLE, 3, &dimids[0], &ncdat->overlaps_id)
    );
    check_nccall(iret, 
         nc_def_var(ncdat->id, "DM", NC_DOUBLE, 4, &dimids[0], &ncdat->DM_id)
    );
    check_nccall(iret, 
         nc_def_var(ncdat->id, "coeff_diag", NC_DOUBLE, 2, &dimids[0], &ncdat->coeff_diag_id)
    );
    // change in dimensions  as hop prob is real
    dimids[1] = NSTATES_ID;
    check_nccall(iret, 
         nc_def_var(ncdat->id, "hopprob", NC_DOUBLE, 2, &dimids[0], &ncdat->hopprop_id)
    );
    // change in dimensions 
    dimids[1] = THREE_ID;
    check_nccall(iret, 
         nc_def_var(ncdat->id, "Energy", NC_DOUBLE, 2, &dimids[0], &ncdat->e_id)
    );
    // change in dimensions for geom, veloc (frames, natoms, 3) 
    dimids[1] = NATOMS_ID;
    dimids[2] = THREE_ID;
    check_nccall(iret, 
         nc_def_var(ncdat->id, "geom", NC_DOUBLE, 3, &dimids[0], &ncdat->crd_id)
    );
    check_nccall(iret, 
         nc_def_var(ncdat->id, "veloc", NC_DOUBLE, 3, &dimids[0], &ncdat->veloc_id)
    );
    // change in dimensions  for one dimensional objects
    dimids[1] = ONE_ID;
    check_nccall(iret, 
         nc_def_var(ncdat->id, "randnum", NC_DOUBLE, 2, &dimids[0], &ncdat->randnum_id)
    );
    check_nccall(iret, 
         nc_def_var(ncdat->id, "state_diag", NC_INT, 2, &dimids[0], &ncdat->state_diag_id)
    );
    check_nccall(iret, 
         nc_def_var(ncdat->id, "state_MCH", NC_INT, 2, &dimids[0], &ncdat->state_MCH_id)
    );
    check_nccall(iret, 
         nc_def_var(ncdat->id, "time_step", NC_INT, 2, &dimids[0], &ncdat->time_step_id)
    );
    // end definition section
    check_nccall(iret, nc_enddef(ncdat->id));
};

// --------------------------------------------------------------------------

void 
reopen_ncoutputdat(int natoms, int nstates, struct sharc_ncoutput* ncdat)
{
    /*
     * defines all variables used in ncoutput
     */
    
    // error handler
    int iret = 0;
    // open sharc file
//     printf("REOPENING!\n");
    ncdat->id      = open_ncfile("output.dat.nc", NC_WRITE);
//     printf("REOPENED!\n");

    // init nsteps
    int nsteps = 0;

    int unlim_id = 0;

    check_nccall(iret, 
          nc_inq_unlimdim(ncdat->id, &unlim_id)
    );

    check_nccall(iret,
            nc_inq_dimlen(ncdat->id, unlim_id, nsteps)
    );

//     printf("found %d steps\n", nsteps);
    
    check_nccall(iret, 
            nc_inq_varid(ncdat->id, "H_MCH", &ncdat->H_MCH_id)
    );
    check_nccall(iret, nc_inq_varid(ncdat->id, "U", &ncdat->U_id));
    check_nccall(iret, nc_inq_varid(ncdat->id, "Ovlap", &ncdat->overlaps_id));
    check_nccall(iret, nc_inq_varid(ncdat->id, "DM", &ncdat->DM_id));
    check_nccall(iret, nc_inq_varid(ncdat->id, "coeff_diag", &ncdat->coeff_diag_id));
    check_nccall(iret, nc_inq_varid(ncdat->id, "hopprob", &ncdat->hopprop_id));
    check_nccall(iret, nc_inq_varid(ncdat->id, "Energy", &ncdat->e_id));
    check_nccall(iret, nc_inq_varid(ncdat->id, "geom", &ncdat->crd_id));
    check_nccall(iret, nc_inq_varid(ncdat->id, "veloc", &ncdat->veloc_id));
    check_nccall(iret, nc_inq_varid(ncdat->id, "randnum", &ncdat->randnum_id));
    check_nccall(iret, nc_inq_varid(ncdat->id, "state_diag", &ncdat->state_diag_id));
    check_nccall(iret, nc_inq_varid(ncdat->id, "state_MCH", &ncdat->state_MCH_id));
    check_nccall(iret, nc_inq_varid(ncdat->id, "time_step", &ncdat->time_step_id));

};

// --------------------------------------------------------------------------

void
write_sharc_ncoutputdat_istep_(
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
        )
{
    int iret = 0;
//     printf("Entering write_sharc_ncoutputdat_istep_ %d\n",*istep);
    if (*istep == 0) {
        setup_ncoutputdat(*natoms, *nstates, ncdat);
    } else if(*istep < 0) {
        reopen_ncoutputdat(*natoms, *nstates, ncdat);
    }

    // counter does not change 
    size_t count[4] = {1, 2* *nstates, *nstates, 3};

    size_t start[4] = {*istep, 0, 0, 0};
    if(*istep < 0) {
      start[0] *= -1;
    }
    check_nccall(iret, 
            nc_put_vara_double(ncdat->id, ncdat->H_MCH_id, start, count, H_MCH_ss)
    );
    check_nccall(iret, 
            nc_put_vara_double(ncdat->id, ncdat->U_id, start, count, U_ss)
    );
    check_nccall(iret, 
            nc_put_vara_double(ncdat->id, ncdat->DM_id, start, count, DM_print_ssd)
    );
    check_nccall(iret, 
            nc_put_vara_double(ncdat->id, ncdat->overlaps_id, start, count, overlaps_ss)
    );
    check_nccall(iret, 
            nc_put_vara_double(ncdat->id, ncdat->coeff_diag_id, start, count, coeff_diag_s)
    );
    // change dimension
    count[1] = 3;
    check_nccall(iret, 
            nc_put_vara_double(ncdat->id, ncdat->e_id, start, count, E)
    );
    // change dimension
    count[1] = *nstates;
    check_nccall(iret, 
            nc_put_vara_double(ncdat->id, ncdat->hopprop_id, start, count, hopprob_s)
    );
    // change dimensions
    count[1] = *natoms;
    count[2] = 3;
    check_nccall(iret, 
            nc_put_vara_double(ncdat->id, ncdat->crd_id, start, count, geom_ad)
    );
    check_nccall(iret, 
            nc_put_vara_double(ncdat->id, ncdat->veloc_id, start, count, veloc_ad)
    );
    // change dimensions
    count[1] = 1;
    check_nccall(iret, 
            nc_put_vara_double(ncdat->id, ncdat->randnum_id, start, count, randnum)
    );
    check_nccall(iret, 
            nc_put_vara_int(ncdat->id, ncdat->state_diag_id, start, count, state_diag)
    );
    check_nccall(iret, 
            nc_put_vara_int(ncdat->id, ncdat->state_MCH_id, start, count, state_MCH)
    );
    check_nccall(iret, 
            nc_put_vara_int(ncdat->id, ncdat->time_step_id, start, count, time_step)
    );
    

};


// --------------------------------------------------------------------------

void
read_sharc_ncoutputdat_istep_(
        // 
        int* nsteps,
        const int* istep,
        const int* natoms,
        const int* nstates,
        // Multidimensional
        double* H_MCH_ss,           // complex, (frame, 2*nstates, nstates)
        double* U_ss,               // complex, (frame, 2*nstates, nstates)
        double* DM_print_ssd,       // complex, (frame, 2*nstates, nstates, 3)
        double* overlaps_ss,        // complex, (frame, 2*nstates, nstates)
        double* coeff_diag_s,       // complex, (frame, 2*nstates)
        double* E,                  // real, contains Etot, Epot and Ekin, (frame, 3)
        double* hopprob_s,          // real, (frame, nstates)
        double* geom_ad,            // real, (frame, natoms, 3)
        double* veloc_ad,           // real, (frame, natoms, 3)
        // One dimensional opjects
        double* randnum,            // real, (frame, 1)
        int* state_diag,            // int, (frame, 1)
        int* state_MCH,             // int, (frame, 1)
        int* time_step,             // int, (frame, 1)
        //
        struct sharc_ncoutput* ncdat
)
{
   int iret = 0;

   if (*istep == 0) {
        ncdat->id = open_ncfile("output.dat.nc", NC_NOWRITE);

        // init nsteps
        *nsteps = 0;

        int unlim_id = 0;

        check_nccall(iret, 
             nc_inq_unlimdim(ncdat->id, &unlim_id)
        );

        check_nccall(iret,
                nc_inq_dimlen(ncdat->id, unlim_id, nsteps)
        );

        printf("found %d steps\n", *nsteps);
        
        check_nccall(iret, 
                nc_inq_varid(ncdat->id, "H_MCH", &ncdat->H_MCH_id)
        );
        check_nccall(iret, nc_inq_varid(ncdat->id, "U", &ncdat->U_id));
        check_nccall(iret, nc_inq_varid(ncdat->id, "Ovlap", &ncdat->overlaps_id));
        check_nccall(iret, nc_inq_varid(ncdat->id, "DM", &ncdat->DM_id));
        check_nccall(iret, nc_inq_varid(ncdat->id, "coeff_diag", &ncdat->coeff_diag_id));
        check_nccall(iret, nc_inq_varid(ncdat->id, "hopprob", &ncdat->hopprop_id));
        check_nccall(iret, nc_inq_varid(ncdat->id, "Energy", &ncdat->e_id));
        check_nccall(iret, nc_inq_varid(ncdat->id, "geom", &ncdat->crd_id));
        check_nccall(iret, nc_inq_varid(ncdat->id, "veloc", &ncdat->veloc_id));
        check_nccall(iret, nc_inq_varid(ncdat->id, "randnum", &ncdat->randnum_id));
        check_nccall(iret, nc_inq_varid(ncdat->id, "state_diag", &ncdat->state_diag_id));
        check_nccall(iret, nc_inq_varid(ncdat->id, "state_MCH", &ncdat->state_MCH_id));
        check_nccall(iret, nc_inq_varid(ncdat->id, "time_step", &ncdat->time_step_id));
    }

   printf("processing step = %d\n", *istep);
   size_t start[4] = {*istep, 0, 0, 0};
   size_t count[4] = {1, *nstates*2, *nstates, 3};

   check_nccall(iret, 
            nc_get_vara_double(ncdat->id, 
                               ncdat->H_MCH_id, 
                               start, 
                               count, 
                               H_MCH_ss)
   );

   check_nccall(iret, 
            nc_get_vara_double(ncdat->id, 
                               ncdat->U_id, 
                               start, 
                               count, 
                               U_ss)
   );

   check_nccall(iret, 
            nc_get_vara_double(ncdat->id, 
                               ncdat->overlaps_id, 
                               start, 
                               count, 
                               overlaps_ss)
   );

   check_nccall(iret, 
            nc_get_vara_double(ncdat->id, 
                               ncdat->DM_id, 
                               start, 
                               count, 
                               DM_print_ssd)
   );

   check_nccall(iret, 
            nc_get_vara_double(ncdat->id, 
                               ncdat->coeff_diag_id, 
                               start, 
                               count, 
                               coeff_diag_s)
   );

   count[1] = 3;

   check_nccall(iret, 
            nc_get_vara_double(ncdat->id, 
                               ncdat->e_id, 
                               start, 
                               count, 
                               E)
   );

   count[1] = *nstates;

   check_nccall(iret, 
            nc_get_vara_double(ncdat->id, 
                               ncdat->hopprop_id, 
                               start, 
                               count, 
                               hopprob_s)
   );

   count[1] = *natoms;
   count[2] = 3;

   check_nccall(iret, 
            nc_get_vara_double(ncdat->id, 
                               ncdat->crd_id, 
                               start, 
                               count, 
                               geom_ad)
   );

   check_nccall(iret, 
            nc_get_vara_double(ncdat->id, 
                               ncdat->veloc_id, 
                               start, 
                               count, 
                               veloc_ad)
   );

   count[1] = 1;

   check_nccall(iret, 
            nc_get_vara_double(ncdat->id, 
                               ncdat->randnum_id, 
                               start, 
                               count, 
                               randnum)
   );

   check_nccall(iret, 
            nc_get_vara_int(ncdat->id, 
                            ncdat->state_diag_id, 
                            start, 
                            count, 
                            state_diag)
   );

   check_nccall(iret, 
            nc_get_vara_int(ncdat->id, 
                            ncdat->state_MCH_id, 
                            start, 
                            count, 
                            state_MCH)
   );

   check_nccall(iret, 
            nc_get_vara_int(ncdat->id, 
                            ncdat->time_step_id, 
                            start, 
                            count, 
                            time_step)
   );
}

// --------------------------------------------------------------------------


void
write_sharc_ncdat_init_(const double* E, struct sharc_ncdat* ncdat)
{
    // error handler
    int iret = 0;
    // open sharc file
    ncdat->id      = create_ncfile("output.dat.nc", NC_CLOBBER);
    // define dimensions
    int energy_id = define_dimension(ncdat->id, "one", 1);
    int frames_id = define_dimension(ncdat->id, "frame", NC_UNLIMITED);
    // define dimensions
    int dimids[2];
    dimids[0] = frames_id;
    dimids[1] = energy_id;
    // define variable IAN
    check_nccall(iret, 
         nc_def_var(ncdat->id, "Energy", NC_DOUBLE, 2, &dimids[0], &ncdat->energy_id)
    );
    // end definition section
    check_nccall(iret, nc_enddef(ncdat->id));
    // counter per unit
    size_t count[2];
    count[0] = 1;      // Frame
    count[1] = 1;      // 
    // everything starts at 0
    size_t start[2] = {0, 0};
    // write data to file
    // Energy
    check_nccall(iret, 
            nc_put_vara_double(ncdat->id, ncdat->energy_id, start, count, E)
    );
    //
};

// --------------------------------------------------------------------------

void
write_sharc_ncdat_traj_(const int* istep, const double* E, struct sharc_ncdat* ncdat)
{
    // error handler
    int iret = 0;
    // counter per unit
    size_t count[2] = { 1, 1} ; // Frame, 
    // everything starts at 0
    size_t start[2] = {*istep-1, 0};     // istep > 2!
    // write data to file
    // Crd
    check_nccall(iret, 
            nc_put_vara_double(ncdat->id, ncdat->energy_id, start, count, E)
    );
    //
//    close_ncfile_(ncxyz->id);
};


// --------------------------------------------------------------------------

void
write_sharc_ncxyz_init_(const int* NAtoms, int* IAn, double* Crd, struct sharc_ncxyz* ncxyz)
{
    // error handler
    int iret = 0;
    // open sharc file
    ncxyz->id      = create_ncfile("sharc_traj.nc", NC_CLOBBER);
    // define dimensions
    int xyz_id    = define_dimension(ncxyz->id, "spatial", 3);
    int natoms_id = define_dimension(ncxyz->id, "atom", *NAtoms);
    int frames_id = define_dimension(ncxyz->id, "frame", NC_UNLIMITED);
    // define dimensions
    int dimids[3];
    dimids[0] = frames_id;
    dimids[INATOMS] = natoms_id;
    dimids[ISPATIAL] = xyz_id;
    // define variable IAN
    check_nccall(iret, 
         nc_def_var(ncxyz->id, "IAn", NC_INT, 1, &natoms_id, &ncxyz->ian_id)
    );
    // define variable Crd 
    check_nccall(iret, 
         nc_def_var(ncxyz->id, "Crd", NC_DOUBLE, 3, dimids, &ncxyz->crd_id)
    );
    // end definition section
    check_nccall(iret, nc_enddef(ncxyz->id));
    // counter per unit
    size_t count[3];
    count[0] = 1;      // Frame
    count[INATOMS] = *NAtoms; // NAtoms
    count[ISPATIAL] = 3;      // Spatial
    // everything starts at 0
    size_t start[3] = {0, 0, 0};
    // write data to file
    // IAn
    check_nccall(iret, 
            nc_put_var_int(ncxyz->id, ncxyz->ian_id, IAn)
    );
    // Crd
    check_nccall(iret, 
            nc_put_vara_double(ncxyz->id, ncxyz->crd_id, start, count, Crd)
    );
    //
};

// --------------------------------------------------------------------------

void
write_sharc_ncxyz_traj_(const int* istep, const int* NAtoms, double* Crd, struct sharc_ncxyz* ncxyz)
{
    // error handler
    int iret = 0;
    // counter per unit
    size_t count[3];
    count[0] = 1;      // Frame
    count[INATOMS] = *NAtoms; // NAtoms
    count[ISPATIAL] = 3;      // Spatial
    // everything starts at 0
    size_t start[3] = {*istep-1, 0, 0};     // istep > 2!
    // write data to file
    // Crd
    check_nccall(iret, 
            nc_put_vara_double(ncxyz->id, ncxyz->crd_id, start, count, Crd)
    );
    //
//    close_ncfile_(ncxyz->id);
};
