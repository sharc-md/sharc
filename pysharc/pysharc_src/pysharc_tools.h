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
 */

#ifndef __TOOLS_H_
#define __TOOLS_H_
#ifdef __cplusplus
extern "C" {
#endif
void clear_double(int N, double * vec);
void clear_complex_double(int N, complex double * vec);
void set_gradient(double * gradient, int NAtoms, int IState, double * state_gradient, double scale);
void set_gradient_in_sharc_order(double * gradient, int NAtoms, int NStates, int IState, double * state_gradient, double scale);
void set_nacdr(double * nac, int NAtoms, int NStates,
        int IState, int JState, double * nac_i_j);
void set_nacdr_in_sharc_order(double * nac, 
        int NAtoms, int NStates, int IState, int JState, double * nac_i_j);
#ifdef __cplusplus
}
#endif
#endif
