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
 * Preprocessor flags to ensure compatible with 
 * python 2 and python 3
 *
 */
// MOD_INIT
#if PY_MAJOR_VERSION >= 3
    #define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
#else
    #define MOD_INIT(name) PyMODINIT_FUNC init##name(void)
#endif
// MOD_DEF
#if PY_MAJOR_VERSION >= 3
    #define MOD_DEF(ob,modname, name, doc, methods, reload, traverse, clear, free) \
        static struct PyModuleDef MODULE_DEF_##modname = {                    \
                    PyModuleDef_HEAD_INIT,                                 \
                    name,                                                \
                    doc,                                                   \
                    -1,                                                    \
                    methods,                                               \
                    reload,                                                \
                    traverse,                                              \
                    clear,                                                 \
                    free,                                                  \
        };                                                                 \
        ob = PyModule_Create(&MODULE_DEF_##modname);
#else
    #define MOD_DEF(ob, modname, name, doc, methods, reload, traverse, clear, free) \
        ob = Py_InitModule3(name, methods, doc);
#endif
// INIT_ERROR
#if PY_MAJOR_VERSION >= 3
    #define INIT_ERROR return NULL
#else
    #define INIT_ERROR return
#endif
// others
#if PY_MAJOR_VERSION >= 3
    #define PyInt_FromLong PyLong_FromLong
    #define PyInt_AsLong PyLong_AsLong
// Strings...
    #define PyString_FromString PyUnicode_FromString
#endif


