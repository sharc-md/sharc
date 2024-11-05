//******************************************
//
//    SHARC Program Suite
//
//    Copyright (c) 2023 University of Vienna
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
 * DEFINES the QMin object of the SHARC LIBRARY
 */

typedef struct {
    PyObject_HEAD
    int NAtoms;
    double * Crd;
} QMin;

static void
QMin_dealloc(QMin * self)
{
#ifdef __OWN_SPACE_QMout__
    free(self->Crd);
#endif
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject *
QMin_new(PyTypeObject * type, PyObject *args, PyObject *kwds)
{
    QMin * self;

    self = (QMin *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->NAtoms = 0;
    }
    return (PyObject *)self;
}

static int
QMin_init(QMin *self, PyObject *args)
{

    if (! PyArg_ParseTuple(args, "i", &self->NAtoms))
        goto fail;

#ifdef __OWN_SPACE_QMout__
    /* Allocate memory for the properties ! */
    self->Crd = (double *)malloc(self->NAtoms * 3 * sizeof(double));
    /* if fail goto fail */
    if ( (self->Crd == NULL) ) 
        goto fail;
#else
    double ** Crd_ptr = &self->Crd;
    setQMinPointers( (double **)Crd_ptr );
#endif
    return 0;
  fail:
    //self->tp_dealloc(type, 0);
    return -1;
}


static PyMemberDef QMin_members[] = {
    // do not expose anything else to the python api
    {NULL}  /* Sentinel */
};

static PyObject *
QMin_updateCrd(QMin * self)
{
#ifdef __OWN_SPACE_QMout__
    static int iang=0;
    get_current_coordinates_(&self->NAtoms, self->Crd, &iang);
#endif
    Py_RETURN_NONE;
}


static PyMethodDef QMin_methods[] = {
    {"update_geometry", (PyCFunction)QMin_updateCrd, METH_NOARGS,
     "update the current geometries"},
    {NULL}  /* Sentinel */
};


static PyTypeObject QMinType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "sharc.QMin",             /* tp_name */
    sizeof(QMin),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)QMin_dealloc, /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_compare */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE,   /* tp_flags */
    "QMin object for sharc",  /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    QMin_methods,             /* tp_methods */
    QMin_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)QMin_init,      /* tp_init */
    0,                         /* tp_alloc */
    QMin_new,                 /* tp_new */
};

