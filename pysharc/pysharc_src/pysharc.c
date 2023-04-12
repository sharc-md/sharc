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


/*
 * @author: Maximilian F.S.J. Menger
 * @date: 18.04.2018
 * @version: 0.1.1
 *
 * Python Wrapper for the SHARC LIBRARY
 *
 * Main routine to setup the sharc driver.
 * uses interface.f90 to call sharc.
 */

#include <Python.h>
#include "structmember.h"
#ifdef __NUMPY__ALLOWED__
    #include <numpy/arrayobject.h>
#endif
/* DEFINITIONS USED TO COMMUNICATE */
#include "data.inc"
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
// basic tools
#include "pysharc_tools.h"
#include "libsharc.h"
#include "python_2_3.h"

// pysharc.h
PyObject * get_atomid(void);
//PyObject * get_atom_names(void);

/*********************** GET INFO ********************************************/

/* get constants */
static char get_constants_docstring[] =
    "get_contstans()\n\
    :return: dict";

static PyObject * get_constants(PyObject * self)
{

    PyObject * dct;

    int NConsts = 6;
    char * const_names[] =
        {"au2a", "au2fs", "au2u", "au2rcm", "au2eV", "au2debye"};

    double * consts;
    consts = (double *)malloc(NConsts * sizeof(double));

    get_constants_(consts);

    dct = PyDict_New();

    for (int i=0; i < NConsts; i++) {
        PyObject * pyfloat = PyFloat_FromDouble(consts[i]);
        if (pyfloat == NULL) {
            return NULL;
        }
        PyDict_SetItemString(dct, const_names[i], pyfloat);
    }
    free(consts);

    return dct;
}

/* get current atomid */
PyObject * get_atomid(void)
{

    PyObject * lst;

    int NAtoms = 0;
    get_natoms_(&NAtoms);

    int * IAn;
    IAn = (int *)malloc(NAtoms*sizeof(int));

    get_ian_(&NAtoms, IAn);

    lst = PyList_New(NAtoms);
    for (int i=0; i<NAtoms; i++){
        PyObject * pyfloat = PyInt_FromLong( *(IAn+i) );
        PyList_SetItem(lst, i, pyfloat);
    }
    free(IAn);

    return lst;
}
/* get atom names */
/* Function does not work! */
/*
PyObject * get_atom_names(void)
{

    PyObject * lst;

    int NAtoms = 0;
    get_natoms_(&NAtoms);

    char * atom_names;

    atom_names = (char *)malloc(4*NAtoms*sizeof(char));

    get_element_name_(&NAtoms, atom_names);

    lst = PyList_New(NAtoms);
    for (int i=0; i<NAtoms; i++){
        PyObject * pystring = PyString_FromString( *(atom_names+i*3) );
        PyList_SetItem(lst, i, pystring);
    }
    free(atom_names);

    return lst;

}

*/



/* get current coordinates */
static char get_current_coordinates_docstring[] =
    "get_current_coordinates()\n\
    :return: lst";

static PyObject * get_current_coordinates(PyObject * self, PyObject * args)
{

    PyObject * lst;

    int ang;

    if (!PyArg_ParseTuple(args, "i", &ang))
        return NULL;


    int NAtoms = 0;
    get_natoms_(&NAtoms);

    double * Crd;
    Crd = (double *)malloc(3*NAtoms*sizeof(double));

    get_current_coordinates_(&NAtoms, Crd, &ang);

    lst = PyList_New(NAtoms);
    for (int i=0; i<NAtoms; i++){
        PyObject * pylst = PyList_New(3);
        for (int j =0; j<3; j++){
            PyList_SetItem(pylst, j, PyFloat_FromDouble( *(Crd+i*3+j)));
        }
        PyList_SetItem(lst, i, pylst);
    }
    free(Crd);

    return lst;
}


/* get_basic_info */
static char get_basic_info_docstring[] =
    "setup_sharc(fileName)\n\
    :return: int ";

static PyObject * get_basic_info(PyObject * self)
{
    int N_func_str = 3;
    char * info_names_str [] =
    { "states", "dt", "savedir"  } ;
    void (*get_info_str []) (char *) =
    { get_states_, get_dt_, get_savedir_  };

    int N_func_int = 3;
    char * info_names_int [] =
    { "NAtoms", "NSteps", "istep" } ;
    void (*get_info_int []) (int *) =
    { get_natoms_, get_nsteps_, get_trajstep_ };


    int N_func_pyobj = 1;
    char * info_names_pyobj [] =
        {"IAn"//, "AtNames"
        };
    PyObject * (*get_info_pyobj []) (void) =
        { get_atomid//,   get_atom_names 
        };

    char * string;

    PyObject * dct;

    dct = PyDict_New();
//  Strings
    string = (char *)malloc(STRING_SIZE_S_*sizeof(char));
    for (int i=0; i < N_func_str; i++) {
        get_info_str[i](string);
        PyObject * pystring = PyString_FromString(string);
        if (pystring == NULL) {
            goto fail_string;
        }
        PyDict_SetItemString(dct, info_names_str[i], pystring);
    }
    free(string);
// Integers
    int ivalue = 0;
    for (int i=0; i < N_func_int; i++) {
        get_info_int[i](&ivalue);
        PyObject * pyint = PyInt_FromLong(ivalue);
        if (pyint == NULL) {
            goto fail_int;
        }
        PyDict_SetItemString(dct, info_names_int[i], pyint);
    }
// PyObject
    for (int i=0; i < N_func_pyobj; i++) {
        PyObject * pyobj = get_info_pyobj[i]();
        if (pyobj == NULL) {
            goto fail_int;
        }
        PyDict_SetItemString(dct, info_names_pyobj[i], pyobj);
    }
    

    return dct;
    fail_string:
        free(string);
        Py_XDECREF(dct);
        return NULL;
    fail_int:
        Py_XDECREF(dct);
        return NULL;
}


/* get_all_tasks */
static char get_all_tasks_docstring[] =
    "setup_sharc(fileName)\n\
    :return: int ";

static PyObject * get_all_tasks(PyObject * self, PyObject * args)
{

    int N_func = 3;
    char * task_names_str [] =
    { "tasks", "grad", "nacdr" } ;
    void (*get_task_str []) (char *, int *) =
    { get_tasks_, get_grad_, get_nacdr_ };

    char * string;

    int icall = 0;

    PyObject * dct;


    if (!PyArg_ParseTuple(args, "i", &icall))
        return NULL;

    dct = PyDict_New();

    string = (char *)malloc(STRING_SIZE_L_*sizeof(char));
    for (int i=0; i < N_func; i++) {
        get_task_str[i](string, &icall);
        PyObject * pystring = PyString_FromString(string);
        if (pystring == NULL) {
            goto fail_string;
        }
        PyDict_SetItemString(dct, task_names_str[i], pystring);
    }
    free(string);
    return dct;

    fail_string:
        free(string);
        Py_DECREF(dct);
        return NULL;
}


/* get_tasks */
static char get_tasks_docstring[] =
    "setup_sharc(fileName)\n\
    :return: int ";

static PyObject * get_tasks(PyObject * self, PyObject * args)
{

    char * string;
    int icall = 0;

    if (!PyArg_ParseTuple(args, "i", &icall))
        return NULL;
    // get tasks
    string = (char *)malloc(STRING_SIZE_L_*sizeof(char));
    get_tasks_(string, &icall);
    // create pystring from c string
    PyObject * pystring = PyString_FromString(string);
    if (pystring == NULL) {
        goto fail_string;
    }
    // free memory
    free(string);
    // return pystring
    return pystring;

    fail_string:
        free(string);
        Py_DECREF(pystring);
        return NULL;
}

/* include sharc main */
#include "pysharc_main.c"
/* include everything related to qmout */
#include "pysharc_QMout.c"
/* include everything related to qmin */
#include "pysharc_QMin.c"

// -----------------------------------------------------------------

/* set QMout */
static char set_qmout_docstring[] =
    "setup_sharc(qmout)\n\
    :qmout qmout: needs to be of type qmout ! \n\
    :return: None";

static PyObject * set_qmout(PyObject * self, PyObject * args)
{
    QMout * qmout;
    int icall = 0;
    if (!PyArg_ParseTuple(args, "Oi", &qmout, &icall))
        return NULL;

    if (!PyObject_TypeCheck(qmout, &QMoutType)){
        PyErr_SetString(PyExc_TypeError, "arg #1 needs to be of type QMout! ");
        return NULL;
    }

    const int iset_g = qmout->iset_g;
    const int iset_nacdr = qmout->iset_nacdr;

#ifdef __OWN_SPACE_QMout__
    // qmout are not pointers to the traj object anymore
    // but are seperate memory
    /* Hamiltonian */
    if (qmout->iset_h == 1){
        printf("sharc.c: set h\n");
        set_hamiltonian_(&qmout->NStates, qmout->hamiltonian);
        qmout->iset_h = 0;
    }
    /* DIPOLE MOMENT */
    if (qmout->iset_d == 1){
        printf("sharc.c: set dm\n");
        set_dipolemoments_(&qmout->NStates, qmout->dipole_mom);
        qmout->iset_d = 0;
    }
    /* Gradient */
    if (qmout->iset_g == 1){
        printf("sharc.c: set g\n");
        set_gradients_(&qmout->NStates, &qmout->NAtoms, qmout->gradient);
        qmout->iset_g = 0;
    }
    /* OVERLAP */
    if (qmout->iset_o == 1){
        printf("sharc.c: set o\n");
        set_overlap_(&qmout->NStates, qmout->overlap);
        qmout->iset_o = 0;
    }
    /* Non-adiabatic couplings */
    if (qmout->iset_nacdr == 1) {
        printf("sharc.c: set nacs\n");
        set_nacs_(&qmout->NStates, &qmout->NAtoms, qmout->nac);
        qmout->iset_nacdr = 0;
    }
#else
    // only properties that need to be changed, are done here
    postprocess_qmout_data_(&qmout->iset_h,
                              &qmout->iset_d,
                              &qmout->iset_g,
                              &qmout->iset_o,
                              &qmout->iset_nacdr
                            );
#endif
    /* set phases */
    //set_phases_();
    // Post process data after setting it
    int ISecond = 0;
    if (icall == 1) {
        post_process_data_(&ISecond);
    }
    /*    if nacdr/gradients were not in icall 1
     *    but iscond is true they need to be cleared!
     */
    if (ISecond == 1) {
        // clear memory!
        if (iset_g == 0) {
            clear_double(qmout->NStates * qmout->NAtoms * 3, qmout->gradient);
        }
        if (iset_nacdr == 0) {
            clear_double(qmout->NStates * qmout->NStates * qmout->NAtoms * 3, 
                    qmout->nacdr);
        }
    }
    // Return ISecond
    return Py_BuildValue("i", ISecond);
}

// -----------------------------------------------------------------

/* SHARC METHODS */
static PyMethodDef SHARC_METHODS[] = {
    /* QMout */
    {"set_qmout", (PyCFunction)set_qmout, METH_VARARGS, set_qmout_docstring},
    /* GET INFO  */
    {"get_constants", (PyCFunction)get_constants, METH_NOARGS, get_constants_docstring},
    {"get_basic_info", (PyCFunction)get_basic_info, METH_NOARGS, get_basic_info_docstring},
    {"get_tasks", (PyCFunction)get_tasks, METH_VARARGS, get_tasks_docstring},
    {"get_all_tasks", (PyCFunction)get_all_tasks, METH_VARARGS, get_all_tasks_docstring},
    {"get_crd", (PyCFunction)get_current_coordinates, METH_VARARGS, get_current_coordinates_docstring},
    /* sharc initial qm */
    {"initial_qm_pre", (PyCFunction)initial_qm_pre, METH_NOARGS, initial_qm_pre_docstring},
    {"initial_qm_post", (PyCFunction)initial_qm_post, METH_NOARGS, initial_qm_post_docstring},
    /* SHARC MAIN ROUTINES*/
    {"setup_sharc", (PyCFunction)setup_sharc, METH_VARARGS, setup_sharc_docstring},
    {"initial_step", (PyCFunction)initial_step, METH_VARARGS, initial_step_docstring},
    {"verlet_xstep", (PyCFunction)verlet_xstep, METH_VARARGS, verlet_xstep_docstring},
    {"verlet_vstep", (PyCFunction)verlet_vstep, METH_NOARGS, verlet_vstep_docstring},
    {"verlet_finalize", (PyCFunction)verlet_finalize, METH_VARARGS, verlet_finalize_docstring},
    {"finalize_sharc", (PyCFunction)finalize_sharc, METH_NOARGS, finalize_sharc_docstring},
    {"error_finalize_sharc", (PyCFunction)error_finalize_sharc, METH_NOARGS, error_finalize_sharc_docstring},
    /* SENTINEL */
    {NULL, NULL, 0, NULL}
};

// define sharc_module_init
static PyObject *
sharc_module_init(void)
{
    // check if Python Type is ready!
    if (PyType_Ready(&QMoutType) < 0 )
        return NULL;
    // check if Python Type is ready!
    if (PyType_Ready(&QMinType) < 0 )
        return NULL;
    // Define Module
    PyObject * mod;
    MOD_DEF(mod, sharc, "sharc",
      "Python API for the SHARC MD code",
      SHARC_METHODS,
      NULL, NULL, NULL, NULL)

    if (mod == NULL)
        return NULL;
    // set  QMout module
    Py_INCREF(&QMoutType);
    PyModule_AddObject(mod, "QMout",
            (PyObject *)&QMoutType);
    // set  QMin module
    Py_INCREF(&QMinType);
    PyModule_AddObject(mod, "QMin",
            (PyObject *)&QMinType);
    /* Load `numpy` */
#ifdef __NUMPY__ALLOWED__
        import_array();
#endif
    return mod;
}


/* DEFINE NEW MODULE SHARC */
MOD_INIT(sharc)
{
#if PY_MAJOR_VERSION >= 3
    return sharc_module_init();
#else
    PyObject * mod=sharc_module_init();
    if (mod == NULL)
        return ;
#endif
}
