#!/usr/bin/env python3

# ******************************************
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
# ******************************************

import datetime
import os
import sys
import json
import itertools
import numpy as np
from optparse import OptionParser



def json_load_byteified(file_handle):
    return _byteify(
        json.load(file_handle, object_hook=_byteify),
        ignore_dicts=True
    )


def json_loads_byteified(json_text):
    return _byteify(
        json.loads(json_text, object_hook=_byteify),
        ignore_dicts=True
    )


def _byteify(data, ignore_dicts=False):
    # if this is a list of values, return list of byteified values
    if isinstance(data, list):
        return [_byteify(item, ignore_dicts=True) for item in data]
    # if this is a dictionary, return dictionary of byteified keys and values
    # but only if we haven't already byteified it
    if isinstance(data, dict) and not ignore_dicts:
        return {
            _byteify(key, ignore_dicts=True): _byteify(value, ignore_dicts=True)
            for key, value in data.items()
        }
    # if it's anything else, return it in its original form
    return data


if sys.version_info[0] != 3:
    print('This is a script for Python 3!')
    sys.exit(0)

version = '3.0'
versionneeded = [0.2, 1.0, 2.0, 2.1, float(version)]
versiondate = datetime.date(2023, 4, 1)


# ======================================================================= #


U_TO_AMU = 1. / 5.4857990943e-4   # conversion from g/mol to amu

pthresh = 1.e-5**2

IToMult = {
    1: 'Singlet',
    2: 'Doublet',
    3: 'Triplet',
    4: 'Quartet',
    5: 'Quintet',
    6: 'Sextet',
    7: 'Septet',
    8: 'Octet',
    'Singlet': 1,
    'Doublet': 2,
    'Triplet': 3,
    'Quartet': 4,
    'Quintet': 5,
    'Sextet': 6,
    'Septet': 7,
    'Octet': 8
}
# ======================================================================= #


def displaywelcome():
    print('Script for setup of displacements started...\n')
    string = '\n'
    string += '  ' + '=' * 80 + '\n'
    string += '||' + '{:^80}'.format('') + '||\n'
    string += '||' + '{:^80}'.format('Compute LVC parameters') + '||\n'
    string += '||' + '{:^80}'.format('') + '||\n'
    string += '||' + '{:^80}'.format('Author: Simon Kropf, Sebastian Mai, Severin Polonius') + '||\n'
    string += '||' + '{:^80}'.format('') + '||\n'
    string += '||' + '{:^80}'.format('Version:' + version) + '||\n'
    string += '||' + '{:^80}'.format(versiondate.strftime("%d.%m.%y")) + '||\n'
    string += '||' + '{:^80}'.format('') + '||\n'
    string += '  ' + '=' * 80 + '\n\n'
    string += 'This script automatizes the setup of excited-state calculations for displacements\nfor SHARC dynamics.'
    print(string)


# ======================================================================= #
def itnmstates(states):
    '''Takes an array of the number of states in each multiplicity and generates an iterator over all states specified. Iterates also over all MS values of all states.

    Example:
    [3,0,3] yields 12 iterations with
    1,1,0
    1,2,0
    1,3,0
    3,1,-1
    3,2,-1
    3,3,-1
    3,1,0
    3,2,0
    3,3,0
    3,1,1
    3,2,1
    3,3,1

    Arguments:
    1 list of integers: States specification

    Returns:
    1 integer: multiplicity
    2 integer: state
    3 float: MS value'''

    for i in range(len(states)):
        if states[i] < 1:
            continue
        for k in range(i + 1):
            for j in range(states[i]):
                yield i + 1, j + 1, k - i / 2.
    return


# ======================================================================= #
def readfile(filename):
    try:
        f = open(filename)
        out = f.readlines()
        f.close()
    except IOError:
        print('File %s does not exist!' % (filename))
        sys.exit(12)
    return out
# ======================================================================= #


def read_QMout(path, nstates, natom, request):
    targets = {'h': {'flag': 1,
                     'type': complex,
                     'dim': (nstates, nstates)},
               'dm': {'flag': 2,
                      'type': complex,
                      'dim': (3, nstates, nstates)},
               'grad': {'flag': 3,
                        'type': float,
                        'dim': (nstates, natom, 3)},
               'nacdr': {'flag': 5,
                         'type': float,
                         'dim': (nstates, nstates, natom, 3)},
               'overlap': {'flag': 6,
                           'type': complex,
                           'dim': (nstates, nstates)}
               }

    # read QM.out
    lines = readfile(path)

    # obtain all targets
    QMout = {}
    for t in targets:
        if t in request:
            iline = -1
            while True:
                iline += 1
                if iline >= len(lines):
                    print('Could not find "%s" (flag "%i") in file %s!' % (t, targets[t]['flag'], path))
                    sys.exit(11)
                line = lines[iline]
                if '! %i' % (targets[t]['flag']) in line:
                    break
            values = []
            # =========== single matrix
            if len(targets[t]['dim']) == 2:
                iline += 1
                for irow in range(targets[t]['dim'][0]):
                    iline += 1
                    line = lines[iline].split()
                    if targets[t]['type'] == complex:
                        row = [complex(float(line[2 * i]), float(line[2 * i + 1])) for i in range(targets[t]['dim'][1])]
                    elif targets[t]['type'] == float:
                        row = [float(line[i]) for i in range(targets[t]['dim'][1])]
                    values.append(row)
            # =========== list of matrices
            elif len(targets[t]['dim']) == 3:
                for iblocks in range(targets[t]['dim'][0]):
                    iline += 1
                    block = []
                    for irow in range(targets[t]['dim'][1]):
                        iline += 1
                        line = lines[iline].split()
                        if targets[t]['type'] == complex:
                            row = [complex(float(line[2 * i]), float(line[2 * i + 1])) for i in range(targets[t]['dim'][2])]
                        elif targets[t]['type'] == float:
                            row = [float(line[i]) for i in range(targets[t]['dim'][2])]
                        block.append(row)
                    values.append(block)
            # =========== matrix of matrices
            elif len(targets[t]['dim']) == 4:
                for iblocks in range(targets[t]['dim'][0]):
                    sblock = []
                    for jblocks in range(targets[t]['dim'][1]):
                        iline += 1
                        block = []
                        for irow in range(targets[t]['dim'][2]):
                            iline += 1
                            line = lines[iline].split()
                            if targets[t]['type'] == complex:
                                row = [complex(float(line[2 * i]), float(line[2 * i + 1])) for i in range(targets[t]['dim'][3])]
                            elif targets[t]['type'] == float:
                                row = [float(line[i]) for i in range(targets[t]['dim'][3])]
                            block.append(row)
                        sblock.append(block)
                    values.append(sblock)
            QMout[t] = values

    # pprint.pprint(QMout)
    return QMout


# ======================================================================= #

def LVC_complex_mat(header, mat, deldiag=False, oformat=' % .7e'):
    rnonzero = False
    inonzero = False

    rstr = header + ' R\n'
    istr = header + ' I\n'
    for i in range(len(mat)):
        for j in range(len(mat)):
            val = mat[i][j].real
            if deldiag and i == j:
                val = 0.
            rstr += oformat % val
            if val * val > pthresh:
                rnonzero = True

            val = mat[i][j].imag
            if deldiag and i == j:
                val = 0.
            istr += oformat % val
            if val * val > pthresh:
                inonzero = True

        rstr += '\n'
        istr += '\n'

    retstr = ''
    if rnonzero:
        retstr += rstr
    if inonzero:
        retstr += istr

    return retstr

# ======================================================================= #


def loewdin_orthonormalization(A, debug=False):
    '''
    returns loewdin orthonormalized matrix
    '''

    # S = A^T * A
    S = np.dot(A.T, A)
    if debug:
        print(S)

    # S^d = U^T * S * U
    S_diag_only, U = np.linalg.eigh(S)
    if debug:
        print('test')
        print(U)
        print(U.shape)
        print(U[0].shape)

    # calculate the inverse sqrt of the diagonal matrix
    S_diag_only_inverse_sqrt = [1. / (float(d) ** 0.5) for d in S_diag_only]
    S_diag_inverse_sqrt = np.diag(S_diag_only_inverse_sqrt)

    # calculate inverse sqrt of S
    if debug:
        printmatrix(S_diag_inverse_sqrt)
        printmatrix(U.T,'bla')
    S_inverse_sqrt = np.dot(np.dot(U, S_diag_inverse_sqrt), U.T)

    # calculate loewdin orthonormalized matrix
    A_lo = np.dot(A, S_inverse_sqrt)

    # normalize A_lo
    A_lo = A_lo.T
    length = len(A_lo)
    A_lon = np.zeros((length, length), dtype=complex)

    if debug:
        printmatrix(A_lo,'A_lo before normalization')
    for i in range(length):
        norm_of_col = np.linalg.norm(A_lo[i])
        A_lon[i] = [e / (norm_of_col ** 0.5) for e in A_lo[i]][0]

    return A_lon.T

# ======================================================================= #


def partition_matrix(matrix, multiplicity, states):
    '''
    return the first partitioned matrix of the given multiplicity

    e. g.: (3 0 2) states

      [111, 121, 131,   0,   0,   0,   0,   0,   0]       returns for multiplicity of 1:
      [112, 122, 132,   0,   0,   0,   0,   0,   0]             [111, 121, 131]
      [113, 123, 133,   0,   0,   0,   0,   0,   0]             [112, 122, 132]
      [  0,   0,   0, 311, 321,   0,   0,   0,   0]             [113, 123, 133]
      [  0,   0,   0, 312, 322,   0,   0,   0,   0] ====>
      [  0,   0,   0,   0,   0, 311, 321,   0,   0]       returns for multiplicity of 3:
      [  0,   0,   0,   0,   0, 312, 322,   0,   0]               [311, 321]
      [  0,   0,   0,   0,   0,   0,   0, 311, 321]               [312, 322]
      [  0,   0,   0,   0,   0,   0,   0, 312, 322]

      123 ^= 1...multiplicity
             2...istate
             3...jstate
    '''
    # get start index based on given multiplicity
    start_index = 0
    for i, state in enumerate(states):
        if (i + 1) == multiplicity:
            break
        else:
            start_index += state

    # size of the partition ^= state for given multiplicity
    size = states[multiplicity - 1]

    # create empty partition
    partition = np.zeros((size, size), dtype=complex)

    # get the partition out of the matrix
    for i in range(start_index, start_index + size):
        for j in range(start_index, start_index + size):
            partition[i - start_index][j - start_index] = matrix[i][j]

    return partition

# ======================================================================= #

def phase_correction(matrix, debug=False):
    U = np.array(matrix).real
    #print(U)
    det = np.linalg.det(U)
    #print(det)
    if det < 0:
        U[:,0]*=-1.   # this row/column convention is correct
    #printmatrix(U)

    # sweeps
    l = len(U)
    while True:
        done = True
        for j in range(l):
            for k in range(j+1,l):
                delta  = 3.*(U[j,j]**2+U[k,k]**2)
                delta += 6.*U[k,j]*U[j,k]
                delta += 8.*(U[k,k]+U[j,j])
                for i in range(l):
                    delta -= 3.*(U[j,i]*U[i,j]+U[k,i]*U[i,k])
                #print(j,k,delta)
                if delta<0:
                    U[:,j]*=-1.   # this row/column convention is correct
                    U[:,k]*=-1.   # this row/column convention is correct
                    #printmatrix(U)
                    #print(np.linalg.det(U))
                    done = False
        if done:
            break

    return U.tolist()


# ======================================================================= #
def phase_correction_old(matrix):
    length = len(matrix)
    phase_corrected_matrix = [[.0 for x in range(length)] for x in range(length)]

    for i in range(length):
        diag = matrix[i][i].real

        # look if diag is significant and negative & switch phase
        if diag ** 2 > 0.5 and diag < 0:
            for j in range(length):
                phase_corrected_matrix[j][i] = matrix[j][i] * -1
        # otherwise leave values as is
        else:
            for j in range(length):
                phase_corrected_matrix[j][i] = matrix[j][i]

    return phase_corrected_matrix

# ======================================================================= #


def check_overlap_diagonal(matrix, states, normal_mode, displacement, ignore_problematic_states):
    '''
    Checks for problematic states (diagonals**2 of overlap matrix smaller than 0.5)
    '''
    problematic_states = {}

    for imult in range(len(states)):
        part_matrix = partition_matrix(matrix, imult + 1, states)

        for state in range(len(part_matrix)):
            sum_column = sum([part_matrix[j][state] ** 2 for j in range(len(part_matrix))])
            if sum_column < 0.5:
                print('* Problematic state %i in %i%s: %s' % (state + 1, int(normal_mode), displacement, IToMult[imult + 1]))
                problematic_states[str(normal_mode) + displacement] = imult + 1

    return problematic_states

# ======================================================================= #


def calculate_W_dQi(H, S, e_ref, normal_mode, displ, debug=False):
    '''
    Calculates the displacement matrix
    '''

    # get diagonalised hamiltonian
    H = np.diag([e - e_ref for e in np.diag(H)])

    # do phase correction if necessary
    # if any([x for x in np.diag(S) if x < 0]):
    if debug:
        print(normal_mode, displ)
        printmatrix(S,'before phase 1')
    S = phase_correction(S)
    if debug:
        printmatrix(S,'after phase 1')

    # do loewdin orthonorm. on overlap matrix
    if debug:
        printmatrix(S,'before Loewdin')
    U = loewdin_orthonormalization(np.matrix(S), debug=debug)
    if debug:
        printmatrix(U,'after Loewdin')
    U = np.array(phase_correction(U))
    if debug:
        printmatrix(U,'after phase 2')

    return np.dot(np.dot(U, H), U.T)   # TODO: or transposed the other way around?
    # return np.dot(np.dot(U.T, H), U)

# ======================================================================= #

def printmatrix(M,title=''):
    s='%s\n' % title
    for i in M:
      print(i)
      for j in i:
        print(j)
        s+='%6.3f ' % j.real
      s+='\n'
    print(s)
    


def write_LVC_template(INFOS):
    lvc_template_content = '%s\n' % (INFOS['v0f'])
    lvc_template_content += str(INFOS['states'])[1:-1].replace(',', '') + '\n'

    # print INFOS

    # print some infos
    print('\nData extraction started ...')
    print('Number of states:', INFOS['nstates'])
    print('Number of atoms:', len(INFOS['atoms']))
    print('Kappas:', ['numerical', 'analytical'][INFOS['ana_grad']])
    print('Lambdas:', ['numerical', 'analytical'][INFOS['ana_nac']])
    print
    print('Reading files ...')
    print

    # extract data from central point
    requests = ['h', 'dm']
    if INFOS['ana_grad']:
        requests.append('grad')
    if INFOS['ana_nac']:
        requests.append('nacdr')
    path = os.path.join(INFOS['paths']['0eq'], 'QM.out')
    print(path, requests)
    QMout_eq = read_QMout(path, INFOS['nstates'], len(INFOS['atoms']), requests)

    # ------------------ epsilon ----------------------
    epsilon_str_list = []

    i = 0
    e_ref = QMout_eq['h'][0][0]

    # run through all multiplicities
    for imult in range(len(INFOS['states'])):
        # partition matrix for every multiplicity
        partition = partition_matrix(QMout_eq['h'], imult + 1, INFOS['states'])

        # run over diagonal and get epsilon values
        for istate in range(len(partition)):
            epsilon_str_list.append('%3i %3i % .10f\n' % (imult + 1, istate + 1, (partition[istate][istate] - e_ref).real))

    # add results to template string
    lvc_template_content += 'epsilon\n'
    lvc_template_content += '%i\n' % (len(epsilon_str_list))
    lvc_template_content += ''.join(sorted(epsilon_str_list))

    # ------------------- kappa -----------------------
    nkappa = 0
    kappa_str_list = []
    r3N = [i for i in range(3 * len(INFOS['atoms']))]

    # run through all possible states
    if INFOS['ana_grad']:
        for i, sti in enumerate(itnmstates(INFOS['states'])):
            imult, istate, ims = sti

            if ims == (imult - 1) / 2.:
                # puts the gradient matrix into a list, has form: [ x, y, z, x, y, z, x, y, z]
                gradient = list(itertools.chain(*QMout_eq['grad'][i]))

                # runs through normal modes
                for normal_mode in INFOS['fmw_normal_modes'].keys():

                    # calculates kappa from normal modes and grad
                    kappa = sum([INFOS['fmw_normal_modes'][normal_mode][ixyz] * gradient[ixyz] for ixyz in r3N])

                    # writes kappa to result string
                    if kappa ** 2 > pthresh:
                        kappa_str_list.append('%3i %3i %5i % .5e\n' % (imult, istate, int(normal_mode), kappa))
                        nkappa += 1

    # ------------------------ lambda --------------------------
    lam = 0
    nlambda = 0
    lambda_str_list = []

    if INFOS['ana_nac']:

        for i, sti in enumerate(itnmstates(INFOS['states'])):
            imult, istate, ims = sti

            if ims != (imult - 1) / 2.:
                continue

            for j, stj in enumerate(itnmstates(INFOS['states'])):
                jmult, jstate, jms = stj

                if jms != (jmult - 1) / 2.:
                    continue

                if i >= j:
                    continue

                if imult != jmult:
                    continue

                if ims != jms:
                    continue

                nacvector = list(itertools.chain(*QMout_eq['nacdr'][i][j]))

                # runs through normal modes
                for normal_mode in INFOS['fmw_normal_modes'].keys():

                    # calculates lambd from normal modes and grad
                    dE = (QMout_eq['h'][j][j] - QMout_eq['h'][i][i]).real
                    lambd = sum([INFOS['fmw_normal_modes'][normal_mode][ixyz] * nacvector[ixyz] for ixyz in r3N]) * dE

                    # writes lambd to result string
                    if lambd ** 2 > pthresh:
                        # lambda_str_list.append('%3i %3i %5i % .5e\n' % (imult, istate, int(normal_mode), lambd))
                        lambda_str_list.append('%3i %3i %3i %3i % .5e\n' % (imult, istate, jstate, int(normal_mode), lambd))
                        nlambda += 1


    # ------------------------ numerical kappas and lambdas --------------------------

    if not (INFOS['ana_nac'] and INFOS['ana_grad']):
        if 'displacements' not in INFOS:
            print('No displacement info found in "displacements.json"!')
            sys.exit(1)

        if not INFOS['ana_nac'] and not INFOS['ana_grad']:
            whatstring = 'kappas and lambdas'
        elif not INFOS['ana_grad']:
            whatstring = 'kappas'
        elif not INFOS['ana_nac']:
            whatstring = 'lambdas'

        # running through all normal modes
        for normal_mode, v in INFOS['normal_modes'].items():

            twosided = False

            # get pos displacement
            pos_displ_mag = INFOS['displacement_magnitudes'][normal_mode]

            # get hamiltonian & overlap matrix from QM.out
            path = os.path.join(INFOS['paths'][str(normal_mode) + 'p'], 'QM.out')
            requests = ['h', 'overlap']
            print(path, requests)
            pos_H, pos_S = read_QMout(path, INFOS['nstates'], len(INFOS['atoms']), requests).values()

            # check diagonal of S & print warning
            INFOS['problematic_mults'] = check_overlap_diagonal(pos_S, INFOS['states'], normal_mode, 'p', INFOS['ignore_problematic_states'])

            # calculate displacement matrix
            pos_W_dQi = calculate_W_dQi(pos_H, pos_S, e_ref, normal_mode, 'p', debug=INFOS['debug'])
            print('Mode %s positive' % normal_mode)
            # printmatrix(pos_H, 'Hpos')
            # printmatrix(pos_S, 'Spos')
            # print(np.linalg.det(pos_S))
            #printmatrix(pos_W_dQi, 'Wpos')


            # Check for two-sided differentiation
            if str(normal_mode) + 'n' in INFOS['displacements']:
                twosided = True
                # get neg displacement
                neg_displ_mag = INFOS['displacement_magnitudes'][normal_mode]

                # get hamiltonian & overlap matrix from QM.out
                path = os.path.join(INFOS['paths'][str(normal_mode) + 'n'], 'QM.out')
                requests = ['h', 'overlap']
                print(path, requests)
                neg_H, neg_S = read_QMout(path, INFOS['nstates'], len(INFOS['atoms']), requests).values()

                # check diagonal of S & print warning if wanted
                INFOS['problematic_mults'].update(check_overlap_diagonal(neg_S, INFOS['states'], normal_mode, 'n', INFOS['ignore_problematic_states']))

                # calculate displacement matrix
                neg_W_dQi = calculate_W_dQi(neg_H, neg_S, e_ref, normal_mode, 'n', debug=INFOS['debug'])
                print('Mode %s negative' % normal_mode)
                # printmatrix(neg_H, 'Hneg')
                # printmatrix(neg_S, 'Sneg')
                # print(np.linalg.det(neg_S))
                #printmatrix(neg_W_dQi, 'Wneg')
            #if twosided:
                #printmatrix(pos_W_dQi-neg_W_dQi, 'Difference')



            # Loop over multiplicities to get kappas and lambdas
            for imult in range(len(INFOS['states'])):

                # checking problematic states
                if INFOS['ignore_problematic_states']:
                    if str(normal_mode) + 'p' in INFOS['problematic_mults']:
                        if INFOS['problematic_mults'][str(normal_mode) + 'p'] == imult + 1:
                            print('Not producing %s for normal mode: %s' % (whatstring, normal_mode))
                            continue
                    if str(normal_mode) + 'n' in INFOS['problematic_mults']:
                        if twosided and INFOS['problematic_mults'][str(normal_mode) + 'n'] == imult + 1:
                            print('! Not producing %s for multiplicity %i for normal mode: %s' % (whatstring, imult + 1, normal_mode))
                            continue

                # partition matrices
                pos_partition = partition_matrix(pos_W_dQi, imult + 1, INFOS['states'])
                if twosided:
                    neg_partition = partition_matrix(neg_W_dQi, imult + 1, INFOS['states'])
                partition_length = len(pos_partition)

                # get lambdas and kappas
                for i in range(partition_length):
                    if not INFOS['ana_grad']:
                        if not twosided:
                            kappa = pos_partition[i][i].real / pos_displ_mag
                        else:
                            kappa = (pos_partition[i][i] - neg_partition[i][i]).real / (pos_displ_mag + neg_displ_mag)
                        if kappa ** 2 > pthresh:
                            kappa_str_list.append('%3i %3i %5i % .5e\n' % (imult + 1, i + 1, int(normal_mode), kappa))
                            nkappa += 1

                    if not INFOS['ana_nac']:
                        for j in range(partition_length):
                            if i >= j:
                                continue
                            if not twosided:
                                lam = pos_partition[i][j].real / pos_displ_mag
                            else:
                                lam = (pos_partition[i][j] - neg_partition[i][j]).real / (pos_displ_mag + neg_displ_mag)
                            if lam ** 2 > pthresh:
                                lambda_str_list.append('%3i %3i %3i %3i % .5e\n' % (imult + 1, i + 1, j + 1, int(normal_mode), lam))
                                nlambda += 1






    # add results to template string
    lvc_template_content += 'kappa\n'
    lvc_template_content += '%i\n' % (nkappa)
    lvc_template_content += ''.join(sorted(kappa_str_list))

    lvc_template_content += 'lambda\n'
    lvc_template_content += '%i\n' % (nlambda)
    lvc_template_content += ''.join(sorted(lambda_str_list))


    # ----------------------- matrices ------------------------------
    lvc_template_content += LVC_complex_mat('SOC', QMout_eq['h'], deldiag=True)
    lvc_template_content += LVC_complex_mat('DMX', QMout_eq['dm'][0])
    lvc_template_content += LVC_complex_mat('DMY', QMout_eq['dm'][1])
    lvc_template_content += LVC_complex_mat('DMZ', QMout_eq['dm'][2])


    # -------------------- write to file ----------------------------
    print('\nFinished!\nLVC parameters written to file: LVC.template\n')
    lvc_template = open('LVC.template', 'w')
    lvc_template.write(lvc_template_content)
    lvc_template.close()


# ======================================================================= #
# ======================================================================= #
# ======================================================================= #

def main():
    '''Main routine'''
    script_name = sys.argv[0].split('/')[-1]

    usage = '''python %s''' % (script_name)

    parser = OptionParser(usage=usage, description='')
    parser.add_option("-d", "--debug",
                  action="store_true", dest="debug", default=False,
                  help="Print all involved matrices and quantities for debugging")
    (options, args) = parser.parse_args()

    displaywelcome()

    # load INFOS object from file
    displacement_info_filename = 'displacements.json'
    try:
        with open(displacement_info_filename, 'r') as displacement_info:
            INFOS = json_load_byteified(displacement_info)
            displacement_info.close()
    except IOError:
        print('IOError during opening readable %s - file. Quitting.' % (displacement_info_filename))
        quit(1)

    # set manually for old calcs
    # INFOS['ignore_problematic_states'] = True
    INFOS['debug'] = options.debug

    # write LVC.template
    write_LVC_template(INFOS)


# ======================================================================= #


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print('\nCtrl+C occured. Exiting.\n')
        sys.exit()
