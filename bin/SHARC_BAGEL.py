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

import traceback
from socket import gethostname
import time
from multiprocessing import Pool
from copy import deepcopy
import datetime
import math
import pprint
import re
import sys
import subprocess as sp
import shutil
import os
global DEBUG
global PRINT
DEBUG = False
PRINT = True


# Modules:
# Operating system, isfile and related routines, move files, create directories
# External Calls to MOLCAS
# Command line arguments
# Regular expressions
# debug print(for dicts and arrays)
# sqrt and other math
# runtime measurement
# copy of arrays of arrays
# parallel calculations
# hostname
# write debug traces when in pool threads
# parse Python literals from input


# =========================================================0

version = '1.1'
versiondate = datetime.date(2021, 10, 25)
changelogstring = '''
25.10.2021: 
- removed PyQuante dependency in get_smat_from_Molden
- rewrote get_smat_from_Molden with pyscf

'''
# ======================================================================= #
# holds the system time when the script was started
starttime = datetime.datetime.now()

# global variables for printing (PRINT gives formatted output, DEBUG gives raw output)
DEBUG = False
PRINT = True

# hash table for conversion of multiplicity to the keywords used in MOLCAS
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

# hash table for conversion of polarisations to the keywords used in MOLCAS
IToPol = {
    0: 'X',
    1: 'Y',
    2: 'Z',
    'X': 0,
    'Y': 1,
    'Z': 2
}

# Number of frozen core orbitals
FROZENS = {'H': 0, 'He': 0,
           'Li': 1, 'Be': 1, 'B': 1, 'C': 1, 'N': 1, 'O': 1, 'F': 1, 'Ne': 1,
           'Na': 1, 'Mg': 1, 'Al': 5, 'Si': 5, 'P': 5, 'S': 5, 'Cl': 5, 'Ar': 5,
           'K': 5, 'Ca': 5,
           'Sc': 5, 'Ti': 5, 'V': 5, 'Cr': 5, 'Mn': 5, 'Fe': 5, 'Co': 5, 'Ni': 5, 'Cu': 5, 'Zn': 5,
           'Ga': 9, 'Ge': 9, 'As': 9, 'Se': 9, 'Br': 9, 'Kr': 9,
           'Rb': 9, 'Sr': 9,
           'Y': 14, 'Zr': 14, 'Nb': 14, 'Mo': 14, 'Tc': 14, 'Ru': 14, 'Rh': 14, 'Pd': 14, 'Ag': 14, 'Cd': 14,
           'In': 18, 'Sn': 18, 'Sb': 18, 'Te': 18, 'I': 18, 'Xe': 18,
           'Cs': 18, 'Ba': 18,
           'La': 23,
           'Ce': 23, 'Pr': 23, 'Nd': 23, 'Pm': 23, 'Sm': 23, 'Eu': 23, 'Gd': 23, 'Tb': 23, 'Dy': 23, 'Ho': 23, 'Er': 23, 'Tm': 23, 'Yb': 23, 'Lu': 23,
           'Hf': 23, 'Ta': 23, 'W': 23, 'Re': 23, 'Os': 23, 'Ir': 23, 'Pt': 23, 'Au': 23, 'Hg': 23,
           'Tl': 23, 'Pb': 23, 'Bi': 23, 'Po': 23, 'At': 23, 'Rn': 23,
           'Fr': 30, 'Ra': 30,
           'Ac': 30,
           'Th': 30, 'Pa': 30, 'U': 30, 'Np': 30, 'Pu': 30, 'Am': 30, 'Cm': 30, 'Bk': 30, 'Cf': 30, 'Es': 30, 'Fm': 30, 'Md': 30, 'No': 30, 'Lr': 30,
           'Rf': 30, 'Db': 30, 'Sg': 30, 'Bh': 30, 'Hs': 30, 'Mt': 30, 'Ds': 30, 'Rg': 30, 'Cn': 30,
           'Nh': 39, 'Fl': 39, 'Mc': 39, 'Lv': 39, 'Ts': 39, 'Og': 39
           }


ATOMCHARGE = {'H': 1, 'He': 2,
              'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
              'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18,
              'K': 19, 'Ca': 20,
              'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30,
              'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36,
              'Rb': 37, 'Sr': 38,
              'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48,
              'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54,
              'Cs': 55, 'Ba': 56,
              'La': 57,
              'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71,
              'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80,
              'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86,
              'Fr': 87, 'Ra': 88,
              'Ac': 89,
              'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103,
              'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112,
              'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118
              }

# conversion factors
au2a = 0.529177211
rcm_to_Eh = 4.556335e-6
D2au = 0.393430307






# =============================================================================================== #
# =============================================================================================== #
# =========================================== general routines ================================== #
# =============================================================================================== #
# =============================================================================================== #

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


def writefile(filename, content):
    # content can be either a string or a list of strings
    try:
        f = open(filename, 'w')
        if isinstance(content, list):
            for line in content:
                f.write(line)
        elif isinstance(content, str):
            f.write(content)
        else:
            print('Content %s cannot be written to file!' % (content))
            sys.exit(13)
        f.close()
    except IOError:
        print('Could not write to file %s!' % (filename))
        sys.exit(14)


# ======================================================================= #
def eformat(f, prec, exp_digits):
    '''Formats a float f into scientific notation with prec number of decimals and exp_digits number of exponent digits.

    String looks like:
    [ -][0-9]\\.[0-9]*E[+-][0-9]*

    Arguments:
    1 float: Number to format
    2 integer: Number of decimals
    3 integer: Number of exponent digits

    Returns:
    1 string: formatted number'''

    s = "% .*e" % (prec, f)
    mantissa, exp = s.split('e')
    return "%sE%+0*d" % (mantissa, exp_digits + 1, int(exp))

# ======================================================================= #


def measuretime():
    '''Calculates the time difference between global variable starttime and the time of the call of measuretime.

    Prints the Runtime, if PRINT or DEBUG are enabled.

    Arguments:
    none

    Returns:
    1 float: runtime in seconds'''

    endtime = datetime.datetime.now()
    runtime = endtime - starttime
    hours = runtime.seconds // 3600
    minutes = runtime.seconds // 60 - hours * 60
    seconds = runtime.seconds % 60
    print('==> Runtime:\n%i Days\t%i Hours\t%i Minutes\t%i Seconds\n\n' % (runtime.days, hours, minutes, seconds))
    total_seconds = runtime.days * 24 * 3600 + runtime.seconds + runtime.microseconds // 1.e6
    return total_seconds

# ======================================================================= #


def removekey(d, key):
    '''Removes an entry from a dictionary and returns the dictionary.

    Arguments:
    1 dictionary
    2 anything which can be a dictionary keyword

    Returns:
    1 dictionary'''

    if key in d:
        r = dict(d)
        del r[key]
        return r
    return d

# ======================================================================= #         OK


def containsstring(string, line):
    '''Takes a string (regular expression) and another string. Returns True if the first string is contained in the second string.

    Arguments:
    1 string: Look for this string
    2 string: within this string

    Returns:
    1 boolean'''

    a = re.search(string, line)
    if a:
        return True
    else:
        return False



# =============================================================================================== #
# =============================================================================================== #
# ============================= iterator routines  ============================================== #
# =============================================================================================== #
# =============================================================================================== #

# ======================================================================= #
def itmult(states):

    for i in range(len(states)):
        if states[i] < 1:
            continue
        yield i + 1
    return

# ======================================================================= #


def itnmstates(states):

    for i in range(len(states)):
        if states[i] < 1:
            continue
        for k in range(i + 1):
            for j in range(states[i]):
                yield i + 1, j + 1, k - i / 2.
    return

# =============================================================================================== #
# =============================================================================================== #
# =========================================== print(routines ==================================== #)
# =============================================================================================== #
# =============================================================================================== #

# ======================================================================= #


def printheader():
    '''Prints the formatted header of the log file. Prints version number and version date

    Takes nothing, returns nothing.'''

    print(starttime, gethostname(), os.getcwd())
    if not PRINT:
        return
    string = '\n'
    string += '  ' + '=' * 80 + '\n'
    string += '||' + ' ' * 80 + '||\n'
    string += '||' + ' ' * 29 + 'SHARC - BAGEL - Interface' + ' ' * 28 + '||\n'
    string += '||' + ' ' * 80 + '||\n'
    string += '||' + ' ' * 20 + 'Authors: Moritz Heindl and Sebastian Mai' + ' ' * 20 + '||\n'
    string += '||' + ' ' * 80 + '||\n'
    string += '||' + ' ' * (36 - (len(version) + 1) // 2) + 'Version: %s' % (version) + ' ' * (35 - (len(version)) // 2) + '||\n'
    lens = len(versiondate.strftime("%d.%m.%y"))
    string += '||' + ' ' * (37 - lens // 2) + 'Date: %s' % (versiondate.strftime("%d.%m.%y")) + ' ' * (37 - (lens + 1) // 2) + '||\n'
    string += '||' + ' ' * 80 + '||\n'
    string += '  ' + '=' * 80 + '\n\n'
    print(string)
    if DEBUG:
        print(changelogstring)

# ======================================================================= #


def printQMin(QMin):

    if not PRINT:
        return
    print('==> QMin Job description for:\n%s' % (QMin['comment']))

    string = 'Mode:   '
    if 'init' in QMin:
        string += '\tINIT'
    if 'restart' in QMin:
        string += '\tRESTART'
    if 'samestep' in QMin:
        string += '\tSAMESTEP'
    if 'newstep' in QMin:
        string += '\tNEWSTEP'

    string += '\nTasks:  '
    if 'h' in QMin:
        string += '\tH'
    if 'soc' in QMin:
        string += '\tSOC'
    if 'dm' in QMin:
        string += '\tDM'
    if 'grad' in QMin:
        string += '\tGrad'
    if 'nacdr' in QMin:
        string += '\tNac(ddr)'
    if 'nacdt' in QMin:
        string += '\tNac(ddt)'
    if 'overlap' in QMin:
        string += '\tOverlap'
    if 'angular' in QMin:
        string += '\tAngular'
    if 'ion' in QMin:
        string += '\tDyson'
    if 'dmdr' in QMin:
        string += '\tDM-Grad'
    if 'socdr' in QMin:
        string += '\tSOC-Grad'
    if 'theodore' in QMin:
        string += '\tTheoDORE'
    if 'phases' in QMin:
        string += '\tPhases'
    print(string)

    string = 'States:        '
    for i in itmult(QMin['states']):
        string += '% 2i %7s  ' % (QMin['states'][i - 1], IToMult[i])
    print(string)

    string = 'Charges:       '
    for i in itmult(QMin['states']):
        string += '%+2i %7s  ' % (QMin['chargemap'][i], '')
    print(string)


    string = 'Method: \t'
    string += 'SA(%i' % (QMin['template']['nstate'][0])
    for i in QMin['template']['nstate'][1:]:
        string += '|%i' % (i)
    string += ')-'
    if QMin['template']['method'] == 'caspt2':
        if QMin['template']['xms'] == "true":
            string += 'XMS-'
        elif QMin['template']['ms'] == "true":
            string += 'MS-'
    string += '%s' % (QMin['template']['method'].upper())
    string += '(%i' % (QMin['Atomcharge'] - QMin['template']['nclosed'] * 2 - QMin['template']['charge'][0])
    for i in range(len(QMin['template']['nstate'][1:])):
        string += '|%i' % (QMin['Atomcharge'] - QMin['template']['nclosed'] * 2 - QMin['template']['charge'][i + 1])
    string += '/%i' % QMin['template']['nact']
    string += ')-'
    string += '/%s' % (os.path.basename(QMin['template']['basis']))
    if QMin['template']['dkh'] == "true":
        string += '\t(Douglas-Kroll)'

    print(string)

    string = 'Found Geo'
    if 'veloc' in QMin:
        string += ' and Veloc! '
    else:
        string += '! '
    string += 'NAtom is %i.\n' % (QMin['natom'])
    print(string)

    string = 'Geometry in Bohrs (%i atoms):\n' % QMin['natom']
    if DEBUG:
        for i in range(QMin['natom']):
            string += '%2s ' % (QMin['geo'][i][0])
            for j in range(3):
                string += '% 7.4f ' % (QMin['geo'][i][j + 1])
            string += '\n'
    else:
        for i in range(min(QMin['natom'], 5)):
            string += '%2s ' % (QMin['geo'][i][0])
            for j in range(3):
                string += '% 7.4f ' % (QMin['geo'][i][j + 1])
            string += '\n'
        if QMin['natom'] > 5:
            string += '..     ...     ...     ...\n'
            string += '%2s ' % (QMin['geo'][-1][0])
            for j in range(3):
                string += '% 7.4f ' % (QMin['geo'][-1][j + 1])
            string += '\n'
    print(string)

    if 'veloc' in QMin and DEBUG:
        string = ''
        for i in range(QMin['natom']):
            string += '%s ' % (QMin['geo'][i][0])
            for j in range(3):
                string += '% 7.4f ' % (QMin['veloc'][i][j])
            string += '\n'
        print(string)

    if 'grad' in QMin:
        string = 'Gradients requested:   '
        for i in range(1, QMin['nmstates'] + 1):
            if i in QMin['grad']:
                string += 'X '
            else:
                string += '. '
        string += '\n'
        print(string)

    if 'nacdr' in QMin:
        string = 'Non-adiabatic couplings:\n'
        for i in range(1, QMin['nmstates'] + 1):
            for j in range(1, QMin['nmstates'] + 1):
                if [i, j] in QMin['nacdr'] or [j, i] in QMin['nacdr']:
                    string += 'X '
                else:
                    string += '. '
            string += '\n'
        print(string)

    if 'overlap' in QMin:
        string = 'Overlaps:\n'
        for i in range(1, QMin['nmstates'] + 1):
            for j in range(1, QMin['nmstates'] + 1):
                if [i, j] in QMin['overlap'] or [j, i] in QMin['overlap']:
                    string += 'X '
                else:
                    string += '. '
        string += '\n'
        print(string)

    print('State map:')
    pprint.pprint(QMin['statemap'])
    print

    for i in sorted(QMin):
        if not any([i == j for j in ['h', 'dm', 'soc', 'dmdr', 'socdr', 'theodore', 'geo', 'veloc', 'states', 'comment', 'grad', 'nacdr', 'ion', 'overlap', 'template', 'statemap']]):
            if not any([i == j for j in ['ionlist']]) or DEBUG:
                string = i + ': '
                string += str(QMin[i])
                print(string)
    print('\n')
    sys.stdout.flush()


# ======================================================================= #
def printcomplexmatrix(matrix, states):
    '''Prints a formatted matrix. Zero elements are not printed, blocks of different mult and MS are delimited by dashes. Also prints a matrix with the imaginary parts, of any one element has non-zero imaginary part.

    Arguments:
    1 list of list of complex: the matrix
    2 list of integers: states specs'''

    nmstates = 0
    for i in range(len(states)):
        nmstates += states[i] * (i + 1)
    string = 'Real Part:\n'
    string += '-' * (11 * nmstates + nmstates // 3)
    string += '\n'
    istate = 0
    for imult, i, ms in itnmstates(states):
        jstate = 0
        string += '|'
        for jmult, j, ms2 in itnmstates(states):
            if matrix[istate][jstate].real == 0.:
                string += ' ' * 11
            else:
                string += '% .3e ' % (matrix[istate][jstate].real)
            if j == states[jmult - 1]:
                string += '|'
            jstate += 1
        string += '\n'
        if i == states[imult - 1]:
            string += '-' * (11 * nmstates + nmstates // 3)
            string += '\n'
        istate += 1
    print(string)
    imag = False
    string = 'Imaginary Part:\n'
    string += '-' * (11 * nmstates + nmstates // 3)
    string += '\n'
    istate = 0
    for imult, i, ms in itnmstates(states):
        jstate = 0
        string += '|'
        for jmult, j, ms2 in itnmstates(states):
            if matrix[istate][jstate].imag == 0.:
                string += ' ' * 11
            else:
                imag = True
                string += '% .3e ' % (matrix[istate][jstate].imag)
            if j == states[jmult - 1]:
                string += '|'
            jstate += 1
        string += '\n'
        if i == states[imult - 1]:
            string += '-' * (11 * nmstates + nmstates // 3)
            string += '\n'
        istate += 1
    string += '\n'
    if imag:
        print(string)

# ======================================================================= #


def printgrad(grad, natom, geo):
    '''Prints a gradient or nac vector. Also prints the atom elements. If the gradient is identical zero, just prints one line.

    Arguments:
    1 list of list of float: gradient
    2 integer: natom
    3 list of list: geometry specs'''

    string = ''
    iszero = True
    for atom in range(natom):
        if not DEBUG:
            if atom == 5:
                string += '...\t...\t     ...\t     ...\t     ...\n'
            if 5 <= atom < natom - 1:
                continue
        string += '%i\t%s\t' % (atom + 1, geo[atom][0])
        for xyz in range(3):
            if grad[atom][xyz] != 0:
                iszero = False
            g = grad[atom][xyz]
            if isinstance(g, float):
                string += '% .5f\t' % (g)
            elif isinstance(g, complex):
                string += '% .5f\t% .5f\t\t' % (g.real, g.imag)
        string += '\n'
    if iszero:
        print('\t\t...is identical zero...\n')
    else:
        print(string)

# ======================================================================= #


def printtheodore(matrix, QMin):
    string = '%6s ' % 'State'
    for i in QMin['template']['theodore_prop']:
        string += '%6s ' % i
    for i in range(len(QMin['template']['theodore_fragment'])):
        for j in range(len(QMin['template']['theodore_fragment'])):
            string += '  Om%1i%1i ' % (i + 1, j + 1)
    string += '\n' + '-------' * (1 + QMin['template']['theodore_n']) + '\n'
    istate = 0
    for imult, i, ms in itnmstates(QMin['states']):
        istate += 1
        string += '%6i ' % istate
        for i in matrix[istate - 1]:
            string += '%6.4f ' % i.real
        string += '\n'
    print(string)

# ======================================================================= #


def printQMout(QMin, QMout):
    '''If PRINT, prints a summary of all requested QM output values. Matrices are formatted using printcomplexmatrix, vectors using printgrad.

    Arguments:
    1 dictionary: QMin
    2 dictionary: QMout'''

    # if DEBUG:
    # pprint.pprint(QMout)
    if not PRINT:
        return
    states = QMin['states']
    nstates = QMin['nstates']
    nmstates = QMin['nmstates']
    natom = QMin['natom']
    print('\n\n>>>>>>>>>>>>> Results\n')
    # Hamiltonian matrix, real or complex
    if 'h' in QMin or 'soc' in QMin:
        eshift = math.ceil(QMout['h'][0][0].real)
        print('=> Hamiltonian Matrix:\nDiagonal Shift: %9.2f' % (eshift))
        matrix = deepcopy(QMout['h'])
        for i in range(nmstates):
            matrix[i][i] -= eshift
        printcomplexmatrix(matrix, states)
    # Dipole moment matrices
    if 'dm' in QMin:
        print('=> Dipole Moment Matrices:\n')
        for xyz in range(3):
            print('Polarisation %s:' % (IToPol[xyz]))
            matrix = QMout['dm'][xyz]
            printcomplexmatrix(matrix, states)
    # Gradients
    if 'grad' in QMin:
        print('=> Gradient Vectors:\n')
        istate = 0
        for imult, i, ms in itnmstates(states):
            print('%s\t%i\tMs= % .1f:' % (IToMult[imult], i, ms))
            printgrad(QMout['grad'][istate], natom, QMin['geo'])
            istate += 1
    # Nonadiabatic couplings
    if 'nacdr' in QMin:
        print('=> Analytical Non-adiabatic coupling vectors:\n')
        istate = 0
        for imult, i, msi in itnmstates(states):
            jstate = 0
            for jmult, j, msj in itnmstates(states):
                if imult == jmult and msi == msj:
                    print('%s\tStates %i - %i\tMs= % .1f:' % (IToMult[imult], i, j, msi))
                    printgrad(QMout['nacdr'][istate][jstate], natom, QMin['geo'])
                jstate += 1
            istate += 1
    # Overlaps
    if 'overlap' in QMin:
        print('=> Overlap matrix:\n')
        matrix = QMout['overlap']
        printcomplexmatrix(matrix, states)
        if 'phases' in QMout:
            print('=> Wavefunction Phases:\n')
            for i in range(nmstates):
                print('% 3.1f % 3.1f' % (QMout['phases'][i].real, QMout['phases'][i].imag))
            print('\n')
    # Spin-orbit coupling derivatives
    if 'socdr' in QMin:
        print('=> Spin-Orbit Gradient Vectors:\n')
        istate = 0
        for imult, i, ims in itnmstates(states):
            jstate = 0
            for jmult, j, jms in itnmstates(states):
                print('%s\t%i\tMs= % .1f -- %s\t%i\tMs= % .1f:' % (IToMult[imult], i, ims, IToMult[jmult], j, jms))
                printgrad(QMout['socdr'][istate][jstate], natom, QMin['geo'])
                jstate += 1
            istate += 1
    # Dipole moment derivatives
    if 'dmdr' in QMin:
        print('=> Dipole moment derivative vectors:\n')
        istate = 0
        for imult, i, msi in itnmstates(states):
            jstate = 0
            for jmult, j, msj in itnmstates(states):
                if imult == jmult and msi == msj:
                    for ipol in range(3):
                        print('%s\tStates %i - %i\tMs= % .1f\tPolarization %s:' % (IToMult[imult], i, j, msi, IToPol[ipol]))
                        printgrad(QMout['dmdr'][ipol][istate][jstate], natom, QMin['geo'])
                jstate += 1
            istate += 1
    # Property matrix (dyson norms)
    if 'ion' in QMin and 'prop' in QMout:
        print('=> Property matrix:\n')
        matrix = QMout['prop']
        printcomplexmatrix(matrix, states)
    # TheoDORE
    if 'theodore' in QMin:
        print('=> TheoDORE results:\n')
        matrix = QMout['theodore']
        printtheodore(matrix, QMin)
    sys.stdout.flush()


# =============================================================================================== #
# =============================================================================================== #
# ======================================= Matrix initialization ================================= #
# =============================================================================================== #
# =============================================================================================== #

# ======================================================================= #         OK
def makecmatrix(a, b):
    '''Initialises a complex axb matrix.

    Arguments:
    1 integer: first dimension
    2 integer: second dimension

    Returns;
    1 list of list of complex'''

    mat = [[complex(0., 0.) for i in range(a)] for j in range(b)]
    return mat

# ======================================================================= #         OK


def makermatrix(a, b):
    '''Initialises a real axb matrix.

    Arguments:
    1 integer: first dimension
    2 integer: second dimension

    Returns;
    1 list of list of real'''

    mat = [[0. for i in range(a)] for j in range(b)]
    return mat


# =============================================================================================== #
# =============================================================================================== #
# =========================================== QMout writing ===================================== #
# =============================================================================================== #
# =============================================================================================== #


# ======================================================================= #
def writeQMout(QMin, QMout, QMinfilename):
    '''Writes the requested quantities to the file which SHARC reads in. The filename is QMinfilename with everything after the first dot replaced by "out".

    Arguments:
    1 dictionary: QMin
    2 dictionary: QMout
    3 string: QMinfilename'''

    k = QMinfilename.find('.')
    if k == -1:
        outfilename = QMinfilename + '.out'
    else:
        outfilename = QMinfilename[:k] + '.out'
    if PRINT:
        print('===> Writing output to file %s in SHARC Format\n' % (outfilename))
    string = ''
    if 'h' in QMin or 'soc' in QMin:
        string += writeQMoutsoc(QMin, QMout)
    if 'dm' in QMin:
        string += writeQMoutdm(QMin, QMout)
    if 'grad' in QMin:
        string += writeQMoutgrad(QMin, QMout)
    if 'overlap' in QMin:
        string += writeQMoutnacsmat(QMin, QMout)
    if 'nacdr' in QMin:
        string += writeQMoutnacana(QMin, QMout)
    if 'socdr' in QMin:
        string += writeQMoutsocdr(QMin, QMout)
    if 'dmdr' in QMin:
        string += writeQMoutdmdr(QMin, QMout)
    if 'ion' in QMin:
        string += writeQMoutprop(QMin, QMout)
    if 'theodore' in QMin:  # or QMin['template']['qmmm']:
        string += writeQMoutTHEODORE(QMin, QMout)
    if 'phases' in QMin:
        string += writeQmoutPhases(QMin, QMout)
    string += writeQMouttime(QMin, QMout)
    outfile = os.path.join(QMin['pwd'], outfilename)
    writefile(outfile, string)
    return

# ======================================================================= #


def writeQMoutsoc(QMin, QMout):
    '''Generates a string with the Spin-Orbit Hamiltonian in SHARC format.

    The string starts with a ! followed by a flag specifying the type of data. In the next line, the dimensions of the matrix are given, followed by nmstates blocks of nmstates elements. Blocks are separated by a blank line.

    Arguments:
    1 dictionary: QMin
    2 dictionary: QMout

    Returns:
    1 string: multiline string with the SOC matrix'''

    states = QMin['states']
    nstates = QMin['nstates']
    nmstates = QMin['nmstates']
    natom = QMin['natom']
    string = ''
    string += '! %i Hamiltonian Matrix (%ix%i, complex)\n' % (1, nmstates, nmstates)
    string += '%i %i\n' % (nmstates, nmstates)
    for i in range(nmstates):
        for j in range(nmstates):
            string += '%s %s ' % (eformat(QMout['h'][i][j].real, 9, 3), eformat(QMout['h'][i][j].imag, 9, 3))
        string += '\n'
    string += '\n'
    return string

# ======================================================================= #


def writeQMoutdm(QMin, QMout):
    '''Generates a string with the Dipole moment matrices in SHARC format.

    The string starts with a ! followed by a flag specifying the type of data. In the next line, the dimensions of the matrix are given, followed by nmstates blocks of nmstates elements. Blocks are separated by a blank line. The string contains three such matrices.

    Arguments:
    1 dictionary: QMin
    2 dictionary: QMout

    Returns:
    1 string: multiline string with the DM matrices'''

    states = QMin['states']
    nstates = QMin['nstates']
    nmstates = QMin['nmstates']
    natom = QMin['natom']
    string = ''
    string += '! %i Dipole Moment Matrices (3x%ix%i, complex)\n' % (2, nmstates, nmstates)
    for xyz in range(3):
        string += '%i %i\n' % (nmstates, nmstates)
        for i in range(nmstates):
            for j in range(nmstates):
                string += '%s %s ' % (eformat(QMout['dm'][xyz][i][j].real, 9, 3), eformat(QMout['dm'][xyz][i][j].imag, 9, 3))
            string += '\n'
        # string+='\n'
    string += '\n'
    return string

# ======================================================================= #


def writeQMoutgrad(QMin, QMout):
    '''Generates a string with the Gradient vectors in SHARC format.

    The string starts with a ! followed by a flag specifying the type of data. On the next line, natom and 3 are written, followed by the gradient, with one line per atom and a blank line at the end. Each MS component shows up (nmstates gradients are written).

    Arguments:
    1 dictionary: QMin
    2 dictionary: QMout

    Returns:
    1 string: multiline string with the Gradient vectors'''

    states = QMin['states']
    nstates = QMin['nstates']
    nmstates = QMin['nmstates']
    natom = QMin['natom']
    string = ''
    string += '! %i Gradient Vectors (%ix%ix3, real)\n' % (3, nmstates, natom)
    i = 0
    for imult, istate, ims in itnmstates(states):
        string += '%i %i ! %i %i %i\n' % (natom, 3, imult, istate, ims)
        for atom in range(natom):
            for xyz in range(3):
                string += '%s ' % (eformat(QMout['grad'][i][atom][xyz], 9, 3))
            string += '\n'
        # string+='\n'
        i += 1
    string += '\n'
    return string

# ======================================================================= #


def writeQMoutnacsmat(QMin, QMout):
    '''Generates a string with the adiabatic-diabatic transformation matrix in SHARC format.

    The string starts with a ! followed by a flag specifying the type of data. In the next line, the dimensions of the matrix are given, followed by nmstates blocks of nmstates elements. Blocks are separated by a blank line.

    Arguments:
    1 dictionary: QMin
    2 dictionary: QMout

    Returns:
    1 string: multiline string with the transformation matrix'''

    states = QMin['states']
    nstates = QMin['nstates']
    nmstates = QMin['nmstates']
    natom = QMin['natom']
    string = ''
    string += '! %i Overlap matrix (%ix%i, complex)\n' % (6, nmstates, nmstates)
    string += '%i %i\n' % (nmstates, nmstates)
    for j in range(nmstates):
        for i in range(nmstates):
            string += '%s %s ' % (eformat(QMout['overlap'][j][i].real, 9, 3), eformat(QMout['overlap'][j][i].imag, 9, 3))
        string += '\n'
    string += '\n'
    return string

# ======================================================================= #


def writeQMoutnacana(QMin, QMout):
    '''Generates a string with the NAC vectors in SHARC format.

    The string starts with a ! followed by a flag specifying the type of data. On the next line, natom and 3 are written, followed by the gradient, with one line per atom and a blank line at the end. Each MS component shows up (nmstates x nmstates vectors are written).

    Arguments:
    1 dictionary: QMin
    2 dictionary: QMout

    Returns:
    1 string: multiline string with the NAC vectors'''

    states = QMin['states']
    nstates = QMin['nstates']
    nmstates = QMin['nmstates']
    natom = QMin['natom']
    string = ''
    string += '! %i Non-adiabatic couplings (ddr) (%ix%ix%ix3, real)\n' % (5, nmstates, nmstates, natom)
    i = 0
    for imult, istate, ims in itnmstates(states):
        j = 0
        for jmult, jstate, jms in itnmstates(states):
            string += '%i %i ! %i %i %i %i %i %i\n' % (natom, 3, imult, istate, ims, jmult, jstate, jms)
            for atom in range(natom):
                for xyz in range(3):
                    string += '%s ' % (eformat(QMout['nacdr'][i][j][atom][xyz], 12, 3))
                string += '\n'
            string += ''
            j += 1
        i += 1
    return string
# ======================================================================= #


def writeQMoutdmdr(QMin, QMout):

    states = QMin['states']
    nmstates = QMin['nmstates']
    natom = QMin['natom']
    string = ''
    string += '! %i Dipole moment derivatives (%ix%ix3x%ix3, real)\n' % (12, nmstates, nmstates, natom)
    i = 0
    for imult, istate, ims in itnmstates(states):
        j = 0
        for jmult, jstate, jms in itnmstates(states):
            for ipol in range(3):
                string += '%i %i ! m1 %i s1 %i ms1 %i   m2 %i s2 %i ms2 %i   pol %i\n' % (natom, 3, imult, istate, ims, jmult, jstate, jms, ipol)
                for atom in range(natom):
                    for xyz in range(3):
                        string += '%s ' % (eformat(QMout['dmdr'][ipol][i][j][atom][xyz], 12, 3))
                    string += '\n'
                string += ''
            j += 1
        i += 1
    string += '\n'
    return string

# ======================================================================= #


def writeQMoutsocdr(QMin, QMout):

    states = QMin['states']
    nmstates = QMin['nmstates']
    natom = QMin['natom']
    string = ''
    string += '! %i Spin-Orbit coupling derivatives (%ix%ix3x%ix3, complex)\n' % (13, nmstates, nmstates, natom)
    i = 0
    for imult, istate, ims in itnmstates(states):
        j = 0
        for jmult, jstate, jms in itnmstates(states):
            string += '%i %i ! m1 %i s1 %i ms1 %i   m2 %i s2 %i ms2 %i\n' % (natom, 3, imult, istate, ims, jmult, jstate, jms)
            for atom in range(natom):
                for xyz in range(3):
                    string += '%s %s ' % (eformat(QMout['socdr'][i][j][atom][xyz].real, 12, 3), eformat(QMout['socdr'][i][j][atom][xyz].imag, 12, 3))
            string += '\n'
            string += ''
            j += 1
        i += 1
    string += '\n'
    return string

# ======================================================================= #


def writeQMoutprop(QMin, QMout):

    nmstates = QMin['nmstates']

    # print(property matrix (flag 11) for backwards compatibility)
    string = ''
    string += '! %i Property Matrix (%ix%i, complex)\n' % (11, nmstates, nmstates)
    string += '%i %i\n' % (nmstates, nmstates)
    for i in range(nmstates):
        for j in range(nmstates):
            string += '%s %s ' % (eformat(QMout['prop'][i][j].real, 12, 3), eformat(QMout['prop'][i][j].imag, 12, 3))
        string += '\n'
    string += '\n'

    # print(property matrices (flag 20) in new format)
    string += '! %i Property Matrices\n' % (20)
    string += '%i    ! number of property matrices\n' % (1)

    string += '! Property Matrix Labels (%i strings)\n' % (1)
    string += 'Dyson norms\n'

    string += '! Property Matrices (%ix%ix%i, complex)\n' % (1, nmstates, nmstates)
    string += '%i %i   ! Dyson norms\n' % (nmstates, nmstates)
    for i in range(nmstates):
        for j in range(nmstates):
            string += '%s %s ' % (eformat(QMout['prop'][i][j].real, 12, 3), eformat(QMout['prop'][i][j].imag, 12, 3))
        string += '\n'
    string += '\n'

    return string

# ======================================================================= #


def writeQMoutTHEODORE(QMin, QMout):

    nmstates = QMin['nmstates']
    nprop = QMin['template']['theodore_n']
    # if QMin['template']['qmmm']:
    # nprop+=7
    if nprop <= 0:
        return '\n'

    string = ''

    string += '! %i Property Vectors\n' % (21)
    string += '%i    ! number of property vectors\n' % (nprop)

    string += '! Property Vector Labels (%i strings)\n' % (nprop)
    descriptors = []
    if 'theodore' in QMin:
        for i in QMin['template']['theodore_prop']:
            descriptors.append('%s' % i)
            string += descriptors[-1] + '\n'
        for i in range(len(QMin['template']['theodore_fragment'])):
            for j in range(len(QMin['template']['theodore_fragment'])):
                descriptors.append('Om_{%i,%i}' % (i + 1, j + 1))
                string += descriptors[-1] + '\n'
    # if QMin['template']['qmmm']:
        # for label in QMout['qmmm_energies']:
            # descriptors.append(label)
            # string+=label+'\n'

    string += '! Property Vectors (%ix%i, real)\n' % (nprop, nmstates)
    if 'theodore' in QMin:
        for i in range(QMin['template']['theodore_n']):
            string += '! TheoDORE descriptor %i (%s)\n' % (i + 1, descriptors[i])
            for j in range(nmstates):
                string += '%s\n' % (eformat(QMout['theodore'][j][i].real, 12, 3))
    # if QMin['template']['qmmm']:
        # for label in QMout['qmmm_energies']:
            #string+='! QM/MM energy contribution (%s)\n' % (label)
            # for j in range(nmstates):
                #string+='%s\n' % (eformat(QMout['qmmm_energies'][label],12,3))
    string += '\n'

    return string

# ======================================================================= #


def writeQmoutPhases(QMin, QMout):

    string = '! 7 Phases\n%i ! for all nmstates\n' % (QMin['nmstates'])
    for i in range(QMin['nmstates']):
        string += '%s %s\n' % (eformat(QMout['phases'][i].real, 9, 3), eformat(QMout['phases'][i].imag, 9, 3))
    return string

# ======================================================================= #


def writeQMouttime(QMin, QMout):
    '''Generates a string with the quantum mechanics total runtime in SHARC format.

    The string starts with a ! followed by a flag specifying the type of data. In the next line, the runtime is given

    Arguments:
    1 dictionary: QMin
    2 dictionary: QMout

    Returns:
    1 string: multiline string with the runtime'''

    string = '! 8 Runtime\n%s\n' % (eformat(QMout['runtime'], 9, 3))
    return string


# =============================================================================================== #
# =============================================================================================== #
# =========================================== SUBROUTINES TO readQMin =========================== #
# =============================================================================================== #
# =============================================================================================== #

# ======================================================================= #
def checkscratch(SCRATCHDIR):
    '''Checks whether SCRATCHDIR is a file or directory. If a file, it quits with exit code 1, if its a directory, it passes. If SCRATCHDIR does not exist, tries to create it.

    Arguments:
    1 string: path to SCRATCHDIR'''

    exist = os.path.exists(SCRATCHDIR)
    if exist:
        isfile = os.path.isfile(SCRATCHDIR)
        if isfile:
            print('$SCRATCHDIR=%s exists and is a file!' % (SCRATCHDIR))
            sys.exit(15)
    else:
        try:
            os.makedirs(SCRATCHDIR)
        except OSError:
            print('Can not create SCRATCHDIR=%s\n' % (SCRATCHDIR))
            sys.exit(16)

# ======================================================================= #


def removequotes(string):
    if string.startswith("'") and string.endswith("'"):
        return string[1:-1]
    elif string.startswith('"') and string.endswith('"'):
        return string[1:-1]
    else:
        return string

# ======================================================================= #


def getBAGELversion(BAGELBIN):
    data = readfile(os.path.join(BAGELBIN, 'version'))
    for line in data:
        if 'release=' in line:
            s = line.split('=')
            x = s[-1].split('.')
            return (int(x[0]), int(x[1]))
    else:
        print('WARNING: Could not detect BAGEL version!')
        return (1900, 0)

# ======================================================================= #


def getsh2BAGELkey(sh2BAGEL, key):
    i = -1
    while True:
        i += 1
        try:
            line = re.sub('#.*$', '', sh2BAGEL[i])
        except IndexError:
            break
        line = line.split(None, 1)
        if line == []:
            continue
        if key.lower() in line[0].lower():
            return line
    return ['', '']

# ======================================================================= #


def get_sh2BAGEL_environ(sh2BAGEL, key, environ=True, crucial=True):
    line = getsh2BAGELkey(sh2BAGEL, key)
    if line[0]:
        LINE = line[1]
        LINE = removequotes(LINE).strip()
    else:
        if environ:
            LINE = os.getenv(key.upper())
            if not LINE:
                if crucial:
                    print('Either set $%s or give path to %s in BAGEL.resources!' % (key.upper(), key.upper()))
                    sys.exit(17)
                else:
                    return None
        else:
            if crucial:
                print('Give path to %s in BAGEL.resources!' % (key.upper()))
                sys.exit(18)
            else:
                return None
    LINE = os.path.expandvars(LINE)
    LINE = os.path.expanduser(LINE)
    if containsstring(';', LINE):
        print("$%s contains a semicolon. Do you probably want to execute another command after %s? I can't do that for you..." % (key.upper(), key.upper()))
        sys.exit(19)
    return LINE

# ======================================================================= #


def get_pairs(QMinlines, i):
    nacpairs = []
    while True:
        i += 1
        try:
            line = QMinlines[i].lower()
        except IndexError:
            print('"keyword select" has to be completed with an "end" on another line!')
            sys.exit(20)
        if 'end' in line:
            break
        fields = line.split()
        try:
            nacpairs.append([int(fields[0]), int(fields[1])])
        except ValueError:
            print('"nacdr select" is followed by pairs of state indices, each pair on a new line!')
            sys.exit(21)
    return nacpairs, i

# ======================================================================= #         OK


def readQMin(QMinfilename):
    '''Reads the time-step dependent information from QMinfilename.

    Arguments:
    1 string: name of the QMin file

    Returns:
    1 dictionary: QMin'''


# --------------------------------------------- QM.in ----------------------------------

    QMinlines = readfile(QMinfilename)
    QMin = {}

    # Get natom
    try:
        natom = int(QMinlines[0])
    except ValueError:
        print('first line must contain the number of atoms!')
        sys.exit(22)
    QMin['natom'] = natom
    if len(QMinlines) < natom + 4:
        print('Input file must contain at least:\nnatom\ncomment\ngeometry\nkeyword "states"\nat least one task')
        sys.exit(23)

    # Save Comment line
    QMin['comment'] = QMinlines[1]

    # Get geometry and possibly velocity (for backup-analytical non-adiabatic couplings)
    QMin['geo'] = []
    QMin['veloc'] = []
    hasveloc = True
    QMin['frozcore'] = 0
    QMin['Atomcharge'] = 0
    for i in range(2, natom + 2):
        if not containsstring('[a-zA-Z][a-zA-Z]?[0-9]*.*[-]?[0-9]+[.][0-9]*.*[-]?[0-9]+[.][0-9]*.*[-]?[0-9]+[.][0-9]*', QMinlines[i]):
            print('Input file does not comply to xyz file format! Maybe natom is just wrong.')
            sys.exit(24)
        fields = QMinlines[i].split()
        fields[0] = fields[0].title()
        symb = fields[0]
        QMin['frozcore'] += FROZENS[symb]
        QMin['Atomcharge'] += ATOMCHARGE[symb]
        for j in range(1, 4):
            fields[j] = float(fields[j])
        QMin['geo'].append(fields[0:4])
        if len(fields) >= 7:
            for j in range(4, 7):
                fields[j] = float(fields[j])
            QMin['veloc'].append(fields[4:7])
        else:
            hasveloc = False
    if not hasveloc:
        QMin = removekey(QMin, 'veloc')


    # Parse remaining file
    i = natom + 1
    while i + 1 < len(QMinlines):
        i += 1
        line = QMinlines[i]
        line = re.sub('#.*$', '', line)
        if len(line.split()) == 0:
            continue
        key = line.lower().split()[0]
        if 'savedir' in key:
            args = line.split()[1:]
        else:
            args = line.lower().split()[1:]
        if key in QMin:
            print('Repeated keyword %s in line %i in input file! Check your input!' % (key, i + 1))
            continue  # only first instance of key in QM.in takes effect
        if len(args) >= 1 and args[0] == 'select':
            pairs, i = get_pairs(QMinlines, i)
            QMin[key] = pairs
        else:
            QMin[key] = args

    if 'unit' in QMin:
        if QMin['unit'][0] == 'angstrom':
            factor = 1. / au2a
        elif QMin['unit'][0] == 'bohr':
            factor = 1.
        else:
            print('Dont know input unit %s!' % (QMin['unit'][0]))
            sys.exit(25)
    else:
        factor = 1. / au2a

    for iatom in range(len(QMin['geo'])):
        for ixyz in range(3):
            QMin['geo'][iatom][ixyz + 1] *= factor



    if 'states' not in QMin:
        print('Number of states not given in QM input file %s!' % (QMinfilename))
        sys.exit(26)
    # Calculate states, nstates, nmstates
    for i in range(len(QMin['states'])):
        QMin['states'][i] = int(QMin['states'][i])
    reduc = 0
    for i in reversed(QMin['states']):
        if i == 0:
            reduc += 1
        else:
            break
    for i in range(reduc):
        del QMin['states'][-1]
    nstates = 0
    nmstates = 0
    for i in range(len(QMin['states'])):
        nstates += QMin['states'][i]
        nmstates += QMin['states'][i] * (i + 1)
    QMin['nstates'] = nstates
    QMin['nmstates'] = nmstates


    # Various logical checks
    possibletasks = ['h', 'soc', 'dm', 'grad', 'overlap', 'dmdr', 'socdr', 'ion', 'theodore', 'phases']
    if not any([i in QMin for i in possibletasks]):
        print('No tasks found! Tasks are %s.' % possibletasks)
        sys.exit(27)

    if 'h' not in QMin:  # and not 'soc' in QMin:
        QMin['h'] = []

#    if 'h' in QMin and 'soc' in QMin:
#        QMin=removekey(QMin,'h')

    if 'soc' in QMin:
        print('Spin-orbit couplings ("soc") are not (yet) implemented')
        sys.exit(28)

#    if 'soc' in QMin and (len(QMin['states'])<3 or QMin['states'][2]<=0):
#        QMin=removekey(QMin,'soc')
#        QMin['h']=[]
#        print('HINT: No triplet states requested, turning off SOC request.')

    if 'samestep' in QMin and 'init' in QMin:
        print('"Init" and "Samestep" cannot be both present in QM.in!')
        sys.exit(29)

    if 'restart' in QMin and 'init' in QMin:
        print('"Init" and "Restart" cannot be both present in QM.in!')
        sys.exit(30)

    if 'phases' in QMin:
        QMin['overlap'] = []

    if 'overlap' in QMin and 'init' in QMin:
        print('"overlap" and "phases" cannot be calculated in the first timestep! Delete either "overlap" or "init"')
        sys.exit(31)

    if 'init' not in QMin and 'samestep' not in QMin and 'restart' not in QMin:
        QMin['newstep'] = []

    if not any([i in QMin for i in ['h', 'soc', 'dm', 'grad']]) and 'overlap' in QMin:
        QMin['h'] = []

    if 'nacdt' in QMin:
        print('Within the SHARC-BAGEL interface, "nacdt" couplings are not supported.')
        sys.exit(32)

    if 'dmdr' in QMin:
        print('Dipole derivatives ("dmdr") not currently supported')
        sys.exit(33)

    if 'socdr' in QMin:
        print('Spin-orbit coupling derivatives ("socdr") are not implemented')
        sys.exit(34)


    # Check for correct gradient list
    if 'grad' in QMin:
        if len(QMin['grad']) == 0 or QMin['grad'][0] == 'all':
            QMin['grad'] = [i + 1 for i in range(nmstates)]
            # pass
        else:
            for i in range(len(QMin['grad'])):
                try:
                    QMin['grad'][i] = int(QMin['grad'][i])
                except ValueError:
                    print('Arguments to keyword "grad" must be "all" or a list of integers!')
                    sys.exit(35)
                if QMin['grad'][i] > nmstates:
                    print('State for requested gradient does not correspond to any state in QM input file state list!')
                    sys.exit(36)

    # Process the non-adiabatic coupling requests
    # type conversion has already been done
    if 'nacdr' in QMin:
        if len(QMin['nacdr']) >= 1:
            nacpairs = QMin['nacdr']
            for i in range(len(nacpairs)):
                if nacpairs[i][0] > nmstates or nacpairs[i][1] > nmstates:
                    print('State for requested non-adiabatic couplings does not correspond to any state in QM input file state list!')
                    sys.exit(37)
        else:
            QMin['nacdr'] = [[j + 1, i + 1] for i in range(nmstates) for j in range(i)]







# --------------------------------------------- BAGEL.resources ----------------------------------

    QMin['pwd'] = os.getcwd()

    # open BAGEL.resources
    filename = 'BAGEL.resources'
    if os.path.isfile(filename):
        sh2BAGEL = readfile(filename)
    else:
        print('No BAGEL.resources file found!')
        sys.exit(38)


    # Set up scratchdir
    line = get_sh2BAGEL_environ(sh2BAGEL, 'scratchdir', False, False)
    if line is None:
        line = QMin['pwd'] + '/SCRATCHDIR/'
    line = os.path.expandvars(line)
    line = os.path.expanduser(line)
    line = os.path.abspath(line)
    # checkscratch(line)
    QMin['scratchdir'] = line
    link(QMin['scratchdir'], os.path.join(QMin['pwd'], 'SCRATCH'), False, False)



    # Set up savedir
    if 'savedir' in QMin:
        # savedir may be read from QM.in file
        line = QMin['savedir'][0]
    else:
        line = get_sh2BAGEL_environ(sh2BAGEL, 'savedir', False, False)
        if line is None:
            line = QMin['pwd'] + '/SAVEDIR/'
    line = os.path.expandvars(line)
    line = os.path.expanduser(line)
    line = os.path.abspath(line)
    if 'init' in QMin:
        checkscratch(line)
    QMin['savedir'] = line
    link(QMin['savedir'], os.path.join(QMin['pwd'], 'SAVE'), False, False)




    # setup environment for BAGEL
    QMin['BAGEL'] = get_sh2BAGEL_environ(sh2BAGEL, 'bagel')
    if not os.path.isfile(os.path.join(QMin['BAGEL'], 'bin', 'BAGEL')):
        print('BAGEL executable at "%s" not found!' % os.path.join(QMin['BAGEL'], 'bin', 'BAGEL'))
        sys.exit(39)
    os.environ['LD_LIBRARY_PATH'] = os.environ['LD_LIBRARY_PATH'] + ':%s' % (QMin['BAGEL'] + '/lib')

    # set up boost library
    QMin['LD_LIB'] = get_sh2BAGEL_environ(sh2BAGEL, 'ld_lib', False, False)
    if QMin['LD_LIB']:
        os.environ['LD_LIBRARY_PATH'] = os.environ['LD_LIBRARY_PATH'] + ':%s' % QMin['LD_LIB']


    # debug option
    line = getsh2BAGELkey(sh2BAGEL, 'debug')
    if line[0]:
        if len(line) <= 1 or 'true' in line[1].lower():
            global DEBUG
            DEBUG = True


    # debug option
    line = getsh2BAGELkey(sh2BAGEL, 'no_print')
    if line[0]:
        if len(line) <= 1 or 'true' in line[1].lower():
            global PRINT
            PRINT = False



    # resources
    QMin['ncpu'] = 1
    line = getsh2BAGELkey(sh2BAGEL, 'ncpu')
    if line[0]:
        try:
            QMin['ncpu'] = int(line[1])
        except ValueError:
            print('Number of CPUs does not evaluate to numerical value!')
            sys.exit(40)
    if os.environ.get('NSLOTS') is not None:
        QMin['ncpu'] = int(os.environ.get('NSLOTS'))
        print('Detected $NSLOTS variable. Will use ncpu=%i' % (QMin['ncpu']))
    elif os.environ.get('SLURM_NTASKS_PER_NODE') is not None:
        QMin['ncpu'] = int(os.environ.get('SLURM_NTASKS_PER_NODE'))
        print('Detected $SLURM_NTASKS_PER_NODE variable. Will use ncpu=%i' % (QMin['ncpu']))
    QMin['ncpu'] = max(1, QMin['ncpu'])

    QMin['mpi'] = False
    line = getsh2BAGELkey(sh2BAGEL, 'mpi_parallel')
    if line[0]:
        if len(line) <= 1 or 'true' in line[1].lower():
            QMin['mpi'] = True

    QMin['delay'] = 0.0
    line = getsh2BAGELkey(sh2BAGEL, 'delay')
    if line[0]:
        try:
            QMin['delay'] = float(line[1])
        except ValueError:
            print('Submit delay does not evaluate to numerical value!')
            sys.exit(41)


    QMin['schedule_scaling'] = 0.9
    line = getsh2BAGELkey(sh2BAGEL, 'schedule_scaling')
    if line[0]:
        try:
            x = float(line[1])
            if 0 < x <= 2.:
                QMin['schedule_scaling'] = x
        except ValueError:
            print('"schedule_scaling" does not evaluate to numerical value!')
            sys.exit(42)


    # initial MO guess settings
    # if neither keyword is present, the interface will reuse MOs from savedir, or let BAGEL generate a guess
    line = getsh2BAGELkey(sh2BAGEL, 'always_orb_init')
    if line[0]:
        QMin['always_orb_init'] = []
    line = getsh2BAGELkey(sh2BAGEL, 'always_guess')
    if line[0]:
        QMin['always_guess'] = []
    if 'always_orb_init' in QMin and 'always_guess' in QMin:
        print('Keywords "always_orb_init" and "always_guess" cannot be used together!')
        sys.exit(43)


    # wfoverlap settings
    if 'overlap' in QMin or 'ion' in QMin:
        # wfoverlap
        QMin['wfoverlap'] = get_sh2BAGEL_environ(sh2BAGEL, 'wfoverlap', False, False)
        if QMin['wfoverlap'] is None:
            ciopath = os.path.join(os.path.expandvars(os.path.expanduser('$SHARC')), 'wfoverlap.x')
            if os.path.isfile(ciopath):
                QMin['wfoverlap'] = ciopath
            else:
                print('Give path to wfoverlap.x in BAGEL.resources!')
                sys.exit(44)

        # PyQuante
        # QMin['pyquante'] = get_sh2BAGEL_environ(sh2BAGEL, 'pyquante', False, False)
        # if QMin['pyquante'] is None or not os.path.isdir(QMin['pyquante']):
        #     print('Give path to the PyQuante installation directory in BAGEL.resources!')
        #     sys.exit(45)
        # if 'PYTHONPATH' in os.environ:
        #     # os.environ['PYTHONPATH']=os.path.join(QMin['pyquante']) + os.pathsep + os.environ['PYTHONPATH']

        #     # os.pathsep + QMin['pyquante']
        #     #    print("already pypath")
        #     #    print(os.environ['PYTHONPATH'])
        #     # else:
        #     #    os.environ['PYTHONPATH']=QMin['pyquante']
        #     #    print("no pypath")
        #     sys.path.append(QMin['pyquante'])

    # memory cannot set memory! #TODO
    QMin['memory'] = 100
    line = getsh2BAGELkey(sh2BAGEL, 'memory')
    if line[0]:
        QMin['memory'] = float(line[1])

    # dipolelevel
    QMin['dipolelevel'] = 0
    line = getsh2BAGELkey(sh2BAGEL, 'dipolelevel')
    if line[0]:
        QMin['dipolelevel'] = int(line[1])

    # truncation threshold
    QMin['wfthres'] = 0.99
    line = getsh2BAGELkey(sh2BAGEL, 'wfthres')
    if line[0]:
        QMin['wfthres'] = float(line[1])

    # get the nooverlap keyword: no dets will be extracted if present
    line = getsh2BAGELkey(sh2BAGEL, 'nooverlap')
    if line[0]:
        QMin['nooverlap'] = []


    # neglected gradients
    QMin['neglected_gradient'] = 'zero'
    if 'grad' in QMin:
        line = getsh2BAGELkey(sh2BAGEL, 'neglected_gradient')
        if line[0]:
            if line[1].lower().strip() == 'zero':
                QMin['neglected_gradient'] = 'zero'
            elif line[1].lower().strip() == 'gs':
                QMin['neglected_gradient'] = 'gs'
            elif line[1].lower().strip() == 'closest':
                QMin['neglected_gradient'] = 'closest'
            else:
                print('Unknown argument to "neglected_gradient"!')
                sys.exit(46)



# --------------------------------------------- BAGEL.template ----------------------------------

    # define classes and defaults
    bools = {'angstrom': "false",
             'dkh': "false",
             'ms': "false",
             'xms': "false",
             'msmr': "false",
             'dipole': "true",
             'shift_imag': "false",
             'orthogonal_basis': "false",
             'numerical': "false"
             }
    strings = {'basis': 'svp',
               'df_basis': 'svp-jkfit',
               'method': 'casscf'
               }
    integers = {'maxiter': 500,
                'maxziter': 100,
                'nact': 0,
                'nclosed': 0,
                'frozen': -1
                }
    floats = {'shift': 0.0
              }
    special = {'basis_per_element': {},
               'nstate': QMin['states'],
               'charge': [i % 2 for i in range(len(QMin['states']))]
               }

    # create QMin subdictionary
    QMin['template'] = {}
    for i in bools:
        QMin['template'][i] = bools[i]
    for i in strings:
        QMin['template'][i] = strings[i]
    for i in integers:
        QMin['template'][i] = integers[i]
    for i in floats:
        QMin['template'][i] = floats[i]
    for i in special:
        QMin['template'][i] = special[i]

    # open template
    template = readfile('BAGEL.template')

    # go through template
    for line in template:
        orig = re.sub('#.*$', '', line).strip()
        line = orig.lower().split()
        if len(line) == 0:
            continue
        elif line[0] in bools:
            QMin['template'][line[0]] = "true"
        elif line[0] in strings:
            QMin['template'][line[0]] = orig.split(None, 1)[1]
        elif line[0] in integers:
            QMin['template'][line[0]] = int(float(line[1]))
        elif line[0] in floats:
            QMin['template'][line[0]] = float(line[1])
        elif line[0] in special:

            # basis_per_element can occur several times
            if line[0] == 'basis_per_element':
                line2 = orig.split(None, 2)
                QMin['template']['basis_per_element'][line2[1]] = os.path.expandvars(os.path.expanduser(line2[2]))


            # charge needs to be autoexpanded, checked and assigned to the multiplicities
            elif line[0] == 'charge':
                if len(line) == 2:
                    charge = int(float(line[1]))
                    if (QMin['Atomcharge'] + charge) % 2 == 1 and len(QMin['states']) > 1:
                        print('HINT: Charge shifted by -1 to be compatible with multiplicities.')
                        charge -= 1
                    QMin['template']['charge'] = [i % 2 + charge for i in range(len(QMin['states']))]
                    print('HINT: total charge per multiplicity automatically assigned, please check (%s).' % QMin['template']['charge'])
                    print('You can set the charge in the template manually for each multiplicity ("charge 0 +1 0 ...")\n')
                elif len(line) - 1 >= len(QMin['states']):
                    QMin['template']['charge'] = [int(float(line[1 + i])) for i in range(len(QMin['states']))]
                    compatible = True
                    for imult, cha in enumerate(QMin['template']['charge']):
                        if not (QMin['Atomcharge'] + cha + imult) % 2 == 0:
                            compatible = False
                    if not compatible:
                        print('WARNING: Charges from template not compatible with multiplicities!  (this is probably OK if you use QM/MM)')
                        sys.exit(47)
#                else:
#                    print('Length of "charge" does not match length of "states"!')
#                    sys.exit(48)
            elif line[0] == 'nstate':
                rootstates = [int(x) for x in line[1:]]
                while rootstates[-1] == 0:
                    rootstates.pop()
                if len(QMin['states']) != len(rootstates):
                    print('WARNING: Different set of multiplicities specified in input and BAGEL.template, exiting')
                    sys.exit(48)
                for i in range(len(QMin['states'])):
                    if QMin['states'][i] > rootstates[i]:
                        print('WARNING: Requesting more states in input than are calculated in the QM call!')
                        sys.exit(49)
                QMin['template']['nstate'] = rootstates



    if QMin['template']['shift_imag'] == "true" and QMin['template']['orthogonal_basis'] == "false":
        print('Use of the imaginary shift is only possible in the orthogonal basis. The corresponding keyword has been set')
        QMin['template']['orthogonal_basis'] = "true"


    line = getsh2BAGELkey(sh2BAGEL, 'numfrozcore')
    if line[0]:
        numfroz = int(line[1])
        if numfroz == 0:
            QMin['frozcore'] = 0
        elif numfroz > 0:
            QMin['frozcore'] = numfroz
        elif numfroz < 0:
            pass        # here we take frozcore from above
    # overrides the old frozcore keyword in the resources file
    if QMin['template']['frozen'] != -1:
        QMin['frozcore'] = QMin['template']['frozen']



    # number of doubly occupied orbitals for Dyson
    line = getsh2BAGELkey(sh2BAGEL, 'numocc')
    if line[0]:
        numfroz = int(line[1])
        if numfroz <= 0:
            QMin['ndocc'] = 0
        elif numfroz > 0:
            QMin['ndocc'] = max(0, numfroz - QMin['frozcore'])
    else:
        QMin['ndocc'] = 0

# --------------------------------------------- Logic ----------------------------------

    # check what type of caspt2 was requested
    QMin['template']['method'] = QMin['template']['method'].lower()
    if QMin['template']['method'] == 'caspt2' and any([True if x > 1 else False for x in QMin['states']]):
        print('CASPT2 nuclear gradients are not implemented in BAGEL for SS-CASPT2 with multiple reference states')
        sys.exit(50)
    elif QMin['template']['method'] == 'ms-caspt2':
        QMin['template']['method'] = 'caspt2'
        QMin['template']['xms'] = 'false'
        QMin['template']['ms'] = 'true'
        print('MS-CASPT2 gradients do not seem to be stable at the moment, violating total energy conversion throughout the dynamics in coupling regions. If static properties are of interest, it is sufficient to comment out the exit statement after this comment.')
        # sys.exit()
    elif QMin['template']['method'] == 'xms-caspt2':
        QMin['template']['method'] = 'caspt2'
        QMin['template']['xms'] = 'true'
        QMin['template']['ms'] = 'true'
    if QMin['template']['msmr'] == 'true':
        QMin['template']['sssr'] = 'false'
    else:
        QMin['template']['sssr'] = 'true'

    # obtain the statemap
    statemap = {}
    i = 1
    for imult, istate, ims in itnmstates(QMin['states']):
        statemap[i] = [imult, istate, ims]
        i += 1
    QMin['statemap'] = statemap

    # obtain the states to actually compute
    QMin['states_to_do'] = deepcopy(QMin['states'])


    # make the jobs
    jobs = {}
    if QMin['states_to_do'][0] > 0:
        jobs[1] = {'mults': [1], 'restr': True}
    if len(QMin['states_to_do']) >= 2 and QMin['states_to_do'][1] > 0:
        jobs[2] = {'mults': [2], 'restr': False}
    if len(QMin['states_to_do']) >= 3 and QMin['states_to_do'][2] > 0:
        jobs[3] = {'mults': [3], 'restr': False}
    if len(QMin['states_to_do']) >= 4:
        for imult, nstate in enumerate(QMin['states_to_do'][3:]):
            if nstate > 0:
                jobs[len(jobs) + 1] = {'mults': [imult + 4], 'restr': False}
    QMin['jobs'] = jobs


    # make the multmap (mapping between multiplicity and job)
    multmap = {}
    for ijob in jobs:
        job = jobs[ijob]
        for imult in job['mults']:
            multmap[imult] = ijob
        multmap[-(ijob)] = job['mults']
    QMin['multmap'] = multmap

    # get the joblist
    joblist = set()
    for i in jobs:
        joblist.add(i)
    joblist = sorted(joblist)
    QMin['joblist'] = joblist
    njobs = len(joblist)
    QMin['njobs'] = njobs


    # make the gsmap
    gsmap = {}
    for i in range(QMin['nmstates']):
        m1, s1, ms1 = tuple(QMin['statemap'][i + 1])
        gs = (m1, 1, ms1)
        job = QMin['multmap'][m1]
        if m1 == 3 and QMin['jobs'][job]['restr']:
            gs = (1, 1, 0.0)
        for j in range(QMin['nmstates']):
            m2, s2, ms2 = tuple(QMin['statemap'][j + 1])
            if (m2, s2, ms2) == gs:
                break
        gsmap[i + 1] = j + 1
    QMin['gsmap'] = gsmap

    # need NAC calculation for transition dipole moments  and grads for static
    if QMin['dipolelevel'] > 0:
        if 'grad' not in QMin and QMin['dipolelevel'] == 2:
            QMin['grad'] = []
        if 'nacdr' not in QMin:
            QMin['nacdr'] = []

    # get the set of states for which gradients actually need to be calculated
    gradmap = set()
    if 'grad' in QMin:
        if QMin['dipolelevel'] == 2:
            for i in range(nmstates):
                QMin['grad'].append(i + 1)
        for i in QMin['grad']:
            gradmap.add(tuple(statemap[i][0:2]))
            #gradmap.add(      (statemap[i][0],1) )
    gradmap = sorted(gradmap)
    QMin['gradmap'] = gradmap


    # get the list of statepairs for NACdr calculation
    nacmap = set()
    if 'nacdr' in QMin:
        if QMin['dipolelevel'] == 1:
            for i in range(2, QMin['states'][0] + 1):
                QMin['nacdr'].append([1, i])
                #print [1,i]
        if QMin['dipolelevel'] == 2:
            for i in range(nmstates):
                for j in range(nmstates):
                    if i != j:
                        QMin['nacdr'].append([i + 1, j + 1])
        for i in QMin['nacdr']:
            s1 = statemap[i[0]][0:2]
            s2 = statemap[i[1]][0:2]
            if s1[0] != s2[0] or s1 == s2:
                continue
            if s1[1] > s2[1]:
                continue
            nacmap.add(tuple(s1 + s2))
    nacmap = list(nacmap)
    nacmap.sort()
    QMin['nacmap'] = nacmap


    # make the chargemap
    QMin['chargemap'] = {}
    for i, c in enumerate(QMin['template']['charge']):
        QMin['chargemap'][i + 1] = c

    # make the ionmap
    if 'ion' in QMin:
        ionmap = []
        for m1 in itmult(QMin['states']):
            job1 = QMin['multmap'][m1]
            el1 = QMin['chargemap'][m1]
            for m2 in itmult(QMin['states']):
                if m1 >= m2:
                    continue
                job2 = QMin['multmap'][m2]
                el2 = QMin['chargemap'][m2]
                if abs(m1 - m2) == 1 and abs(el1 - el2) == 1:
                    ionmap.append((m1, job1, m2, job2))
        QMin['ionmap'] = ionmap

    # number of properties/entries calculated by TheoDORE
#    if 'theodore' in QMin:
#        QMin['template']['theodore_n']=len(QMin['template']['theodore_prop']) + len(QMin['template']['theodore_fragment'])**2
#    else:
#        QMin['template']['theodore_n']=0

# --------------------------------------------- File setup ----------------------------------

    # check for initial orbitals
    initorbs = {}
    if 'always_guess' in QMin:
        QMin['initorbs'] = {}
    elif 'init' in QMin or 'always_orb_init' in QMin:
        for job in QMin['joblist']:
            filename = os.path.join(QMin['pwd'], 'archive.%i.init' % (job))
            if os.path.isfile(filename):
                initorbs[job] = filename
        if 'always_orb_init' in QMin and len(initorbs) < njobs:
            print('Initial orbitals missing for some jobs!')
            sys.exit(52)
        QMin['initorbs'] = initorbs
    elif 'newstep' in QMin:
        for job in QMin['joblist']:
            filename = os.path.join(QMin['savedir'], 'archive.%i' % (job))
            if os.path.isfile(filename):
                initorbs[job] = filename + '.old.archive'     # file will be moved to .old
            else:
                print('File %s missing in savedir!' % (filename))
                sys.exit(53)
        QMin['initorbs'] = initorbs
    elif 'samestep' in QMin:
        for job in QMin['joblist']:
            filename = os.path.join(QMin['savedir'], 'archive.%i' % (job))
            if os.path.isfile(filename):
                initorbs[job] = filename
            else:
                print('File %s missing in savedir!' % (filename))
                sys.exit(54)
        QMin['initorbs'] = initorbs
    elif 'restart' in QMin:
        for job in QMin['joblist']:
            filename = os.path.join(QMin['savedir'], 'archive.%i.old.archive' % (job))
            if os.path.isfile(filename):
                initorbs[job] = filename
            else:
                print('File %s missing in savedir!' % (filename))
                sys.exit(55)
        QMin['initorbs'] = initorbs

    # check for atomic fragment files
#    frags_there=(not 'init' in QMin)
#    for i in QMin['atom_frags']:
#        filename=os.path.join(QMin['savedir'],'frag.t21.%s' % (i))
#        if not os.path.isfile(filename):
#            frags_there=False
#    QMin['frags_there']=frags_there

    # make name for backup directory
    if 'backup' in QMin:
        backupdir = QMin['savedir'] + '/backup'
        backupdir1 = backupdir
        i = 0
        while os.path.isdir(backupdir1):
            i += 1
            if 'step' in QMin:
                backupdir1 = backupdir + '/step%s_%i' % (QMin['step'][0], i)
            else:
                backupdir1 = backupdir + '/calc_%i' % (i)
        QMin['backup'] = backupdir

    if DEBUG:
        print('======= DEBUG print(for QMin =======')
        pprint.pprint(QMin)
        print('====================================')
    return QMin

# =============================================================================================== #
# =============================================================================================== #
# =========================================== Job Scheduling ==================================== #
# =============================================================================================== #
# =============================================================================================== #


def parallel_speedup(N, scaling):
    # computes the parallel speedup from Amdahls law
    # with scaling being the fraction of parallelizable work and (1-scaling) being the serial part
    return 1. / ((1 - scaling) + scaling / N)


def divide_slots(ncpu, ntasks, scaling):
    # this routine figures out the optimal distribution of the tasks over the CPU cores
    #   returns the number of rounds (how many jobs each CPU core will contribute to),
    #   the number of slots which should be set in the Pool,
    #   and the number of cores for each job.
    minpar = 1
    ntasks_per_round = ncpu // minpar
    if ncpu == 1:
        ntasks_per_round = 1
    ntasks_per_round = min(ntasks_per_round, ntasks)
    optimal = {}
    for i in range(1, 1 + ntasks_per_round):
        nrounds = int(math.ceil(float(ntasks) // i))
        ncores = ncpu // i
        optimal[i] = nrounds / parallel_speedup(ncores, scaling)
    # print(optimal)
    best = min(optimal, key=optimal.get)
    nrounds = int(math.ceil(float(ntasks) // best))
    ncores = ncpu // best

    cpu_per_run = [0 for i in range(ntasks)]
    if nrounds == 1:
        itask = 0
        for icpu in range(ncpu):
            cpu_per_run[itask] += 1
            itask += 1
            if itask >= ntasks:
                itask = 0
        nslots = ntasks
    else:
        for itask in range(ntasks):
            cpu_per_run[itask] = ncores
        nslots = ncpu // ncores
    # print(nrounds,nslots,cpu_per_run)
    return nrounds, nslots, cpu_per_run

# =============================================================================================== #


def generate_joblist(QMin):

    # pprint.pprint(QMin)

    # sort the gradients into the different jobs
    gradjob = {}
    for ijob in QMin['joblist']:
        gradjob['master_%i' % ijob] = {}
    for grad in QMin['gradmap']:
        gradjob['grad_%i_%i' % grad] = {}
        gradjob['grad_%i_%i' % grad][grad] = {'gs': False}
        n = 0
        for gradx in gradjob['master_%i' % ijob]:
            n += 1
        if n > 0:
            gradjob['grad_%i_%i' % grad] = {}
            gradjob['grad_%i_%i' % grad][grad] = {'gs': False}

    # make map for states onto gradjobs
    jobgrad = {}
    for job in gradjob:
        for state in gradjob[job]:
            jobgrad[state] = (job, gradjob[job][state]['gs'])
    QMin['jobgrad'] = jobgrad
    schedule = []
    QMin['nslots_pool'] = []

    # add the master calculations
    ntasks = 0
    for i in gradjob:
        if 'master' in i:
            ntasks += 1
    nrounds, nslots, cpu_per_run = divide_slots(QMin['ncpu'], ntasks, QMin['schedule_scaling'])
    QMin['nslots_pool'].append(nslots)
    schedule.append({})
    icount = 0
    for i in sorted(gradjob):
        if 'master' in i:
            QMin1 = deepcopy(QMin)
            QMin1['master'] = True
            QMin1['IJOB'] = int(i.split('_')[1])
            remove = ['ncpu']
            for r in remove:
                QMin1 = removekey(QMin1, r)
            # QMin1['gradmap']=list(gradjob[i])
            QMin1['ncpu'] = cpu_per_run[icount]
            icount += 1
            schedule[-1][i] = QMin1





    return QMin, schedule


# =============================================================================================== #
# =============================================================================================== #
# ======================================= BAGEL Job Execution ===================================== #
# =============================================================================================== #
# =============================================================================================== #

def runjobs(schedule, QMin):

    if 'newstep' in QMin:
        moveOldFiles(QMin)

    print('>>>>>>>>>>>>> Starting the BAGEL job execution')
    errorcodes = {}
    for ijobset, jobset in enumerate(schedule):
        if not jobset:
            continue
        pool = Pool(processes=QMin['nslots_pool'][ijobset])
        for job in jobset:
            QMin1 = jobset[job]
            WORKDIR = os.path.join(QMin['scratchdir'], job)
            errorcodes[job] = pool.apply_async(run_calc, [WORKDIR, QMin1])
            time.sleep(QMin['delay'])
        pool.close()
        pool.join()


    for i in errorcodes:
        errorcodes[i] = errorcodes[i].get()
    j = 0
    string = 'Error Codes:\n'
    for i in errorcodes:
        string += '\t%s\t%i' % (i + ' ' * (10 - len(i)), errorcodes[i])
        j += 1
        if j == 4:
            j = 0
            string += '\n'
    if any((i != 0 for i in errorcodes.values())):
        print('Some subprocesses did not finish successfully!')
        print('See %s:%s for error messages in BAGEL output.' % (gethostname(), QMin['scratchdir']))
        sys.exit(56)
    print
    if PRINT:
        print('>>>>>>>>>>>>> Saving files')
        starttime = datetime.datetime.now()

    for ijobset, jobset in enumerate(schedule):
        if not jobset:
            continue
        for job in jobset:
            if 'master' in job:
                WORKDIR = os.path.join(QMin['scratchdir'], job)
                if 'samestep' not in QMin:
                    saveFiles(WORKDIR, jobset[job])
                if 'ion' in QMin and ijobset == 0:
                    saveAOmatrix(WORKDIR, QMin, job)
    if PRINT:
        endtime = datetime.datetime.now()
        print('Saving Runtime: %s' % (endtime - starttime))
    print

    return errorcodes

# ======================================================================= #


def run_calc(WORKDIR, QMin):
    try:
        setupWORKDIR(WORKDIR, QMin)
        strip = False
        err = runBAGEL(WORKDIR, QMin['BAGEL'], QMin['ncpu'], QMin['mpi'], QMin['savedir'], strip)
    except Exception as problem:
        print('*' * 50 + '\nException in run_calc(%s)!' % (WORKDIR))
        traceback.print_exc()
        print('*' * 50 + '\n')
        raise problem

    return err

# ======================================================================= #


def setupWORKDIR(WORKDIR, QMin):
    # mkdir the WORKDIR, or clean it if it exists, then copy all necessary files from pwd and savedir
    # then put the BAGEL.input file

    # setup the directory
    # print("\n new working directory", WORKDIR, "\n")
    mkdir(WORKDIR)
    # write BAGEL.input
    inputstring = writeBAGELinput(QMin)
    # for string in inputstrings:
    #  print(string)

    filename = os.path.join(WORKDIR, 'BAGEL.run')
    writefile(filename, inputstring)

    if DEBUG:
        print('================== DEBUG input file for WORKDIR %s =================' % (shorten_DIR(WORKDIR)))
        print(inputstring)
        print('BAGEL input written to: %s' % (filename))
        print('====================================================================')

    # wf file copying
    if 'master' in QMin:
        job = QMin['IJOB']
        # print QMin['initorbs'][job]
        if job in QMin['initorbs']:
            fromfile = QMin['initorbs'][job]
            tofile = os.path.join(WORKDIR, 'archive.archive')
            print('copying from file ', fromfile, ' to file ', tofile)
            shutil.copy(fromfile, tofile)
    elif 'grad' in QMin:
        job = QMin['IJOB']
        fromfile = os.path.join(QMin['scratchdir'], 'master_%i' % job, 'archive.archive')
        tofile = os.path.join(WORKDIR, 'archive.archive')
        shutil.copy(fromfile, tofile)

    # copy basis set files
    if not os.path.isabs(QMin['template']['basis']):
        fromfile = os.path.join(QMin['pwd'], QMin['template']['basis'])
        if os.path.isfile(fromfile):
            tofile = os.path.join(WORKDIR, os.path.basename(QMin['template']['basis']))
            shutil.copy(fromfile, tofile)
    if not os.path.isabs(QMin['template']['df_basis']):
        fromfile = os.path.join(QMin['pwd'], QMin['template']['df_basis'])
        if os.path.isfile(fromfile):
            tofile = os.path.join(WORKDIR, os.path.basename(QMin['template']['df_basis']))
            shutil.copy(fromfile, tofile)

    return

# ======================================================================= #


def writeBAGELinput(QMin):

    # pprint.pprint(QMin)

    # general setup
    job = QMin['IJOB']
    gsmult = QMin['multmap'][-job][0]
    restr = QMin['jobs'][job]['restr']
    charge = QMin['chargemap'][gsmult]

    # excited states to calculate
    states_to_do = QMin['states_to_do']
    for imult in range(len(states_to_do)):
        if not imult + 1 in QMin['multmap'][-job]:
            states_to_do[imult] = 0

    # do minimum number of states for gradient jobs
    if 'gradonly' in QMin:
        gradmult = QMin['gradmap'][0][0]
        gradstat = QMin['gradmap'][0][1]
        for imult in range(len(states_to_do)):
            if imult + 1 == gradmult:
                states_to_do[imult] = gradstat - (gradmult == gsmult)
            else:
                states_to_do[imult] = 0

    # number of states to calculate
    if restr:
        ncalc = max(states_to_do)
        sing = states_to_do[0] > 0
        trip = (len(states_to_do) >= 3 and states_to_do[2] > 0)
        if sing and not trip:
            onlysing = True
            onlytrip = False
        elif not sing and trip:
            onlysing = False
            onlytrip = True
        elif sing and trip:
            onlysing = False
            onlytrip = False
    else:
        ncalc = max(states_to_do)
        onlysing = False
        onlytrip = False

    # gradients
    grads = []
    for i in QMin['gradmap']:
        if i[0] == gsmult:
            grads.append(i)

    if len(grads) != 0:
        dograd = True
    else:
        dograd = False

    # non-adiabatic couplings
    pairs = []
    if 'nacdr' in QMin:
        nacmap = QMin['nacmap']
        for nacpair in nacmap:
            if nacpair[0] == job:
                pairs.append((nacpair[1] - 1, nacpair[3] - 1))

    if len(pairs) != 0:
        donac = True
    else:
        donac = False

    # construct the input string
    string = ''

    string += '{ "bagel" : [\n\n'

    # molecule
    string += '  { \n    "title" : "molecule",\n'
    string += '    "%s" : "%s",\n' % ('dkh', QMin['template']['dkh'])
    string += '    "cartesian": "true",\n'

    # basis
    if not os.path.isabs(QMin['template']['basis']):
        basis = os.path.basename(QMin['template']['basis'])
    else:
        basis = QMin['template']['basis']
    string += '    "%s" : "%s",\n' % ('basis', basis)

    # df basis
    if not os.path.isabs(QMin['template']['df_basis']):
        basis = os.path.basename(QMin['template']['df_basis'])
    else:
        basis = QMin['template']['df_basis']
    string += '    "%s" : "%s",\n' % ('df_basis', basis)

    #keys=['basis', 'df_basis','dkh']
    # for key in keys:
    #string+='    "%s" : "%s",\n' % (key,QMin['template'][key])
    string += '    "geometry" : [\n'
    for atom in QMin['geo']:
        string += '    { "atom" : "%s", "xyz" : [ % 7.5f,  % 7.5f,  % 7.5f] },\n' % (atom[0], atom[1], atom[2], atom[3])
    string = string[:-2] + '\n    ]\n'
    string += '  },\n\n'

    # load_ref to restart orbitals from previous calculation
    if job in QMin['initorbs']:
        string += '  { \n    "title" : "load_ref",\n'
        string += '    "continue_geom" : "false",\n'
        string += '    "file" : "./archive"\n  },\n\n'
    else:
        string += '  { \n    "title" : "rohf",\n    "charge" : %i\n  },\n\n' % (charge)

    casscfkeys = ['nact', 'nclosed', 'maxiter']
    caspt2keys = ['ms', 'xms', 'sssr', 'shift', 'shift_imag', 'orthogonal_basis', 'maxiter']

    # we must do a force-less calculation first, otherwise we get garbage Molden files
    string += '  {\n    "title" : "casscf",\n'
    for key in casscfkeys:
        string += '    "%s" : %s,\n' % (key, QMin['template'][key])
    string += '    "nspin" : %i,\n' % (gsmult - 1)
    string += '    "print_thresh" : 1e-10,\n'
    string += '    "charge" : %i,\n' % (charge)
    string += '    "nstate" : %i,\n' % (QMin['template']['nstate'][gsmult - 1])
    string = string[:-2] + '\n  },\n\n'
    if QMin['template']['method'] == 'caspt2':
        string += '  { \n  '
        string += '  "title" : "smith",\n'
        string += '    "method" : "caspt2",\n'
        for key in caspt2keys:
            string += '    "%s" : "%s",\n' % (key, QMin['template'][key])
        string = string[:-2] + '\n   },\n\n'

    # save_ref
    string += '  { \n    "title" : "save_ref",\n'
    string += '    "file" : "%s"\n  },\n\n' % (os.path.join('./archive'))

    # write molden file
    string += '  { \n    "title" : "print",\n'
    string += '    "file" : "%s",\n' % (os.path.join('./orbitals.molden.%i' % QMin['IJOB']))
    string += '    "orbitals" : "true"\n  },\n\n'

    # forces
    if dograd or donac:
        string += '  { \n    "title" : "forces",\n'
        string += '    "grads" : [\n'
        for grad in grads:
            string += '      { "title" : "force", "target" : %i , "maxziter" : %i, "numerical" : %s },\n' % (grad[1] - 1, QMin['template']['maxziter'], QMin['template']['numerical'])
        for nac in pairs:
            string += '      { "title" : "nacme", "target" : %i , "target2" : %i, "nacmtype" : "full", "maxziter" : %i},\n' % (nac[0], nac[1], QMin['template']['maxziter'])
        string = string[:-2] + '\n      ],\n'
        string += '    "method" : [ {\n'
        if QMin['template']['method'] == 'caspt2':
            string += '    "title" : "caspt2",\n    "smith" : {\n'
            string += '      "method" : "caspt2",\n'
            for key in caspt2keys:
                string += '      "%s" : %s,\n' % (key, QMin['template'][key])
            if QMin['template']['frozen'] != -1:
                string += '      "ncore" : %i,\n' % QMin['template']['frozen']
            string = string[:-2] + '\n      },\n'
        else:
            string += '    "title" : "casscf",\n'
        for key in casscfkeys:
            string += '    "%s" : "%s",\n' % (key, QMin['template'][key])
        string += '    "nspin" : %i,\n' % (gsmult - 1)
        string += '    "print_thresh" : 1e-10,\n'
        string += '    "charge" : %i,\n' % (charge)
        string += '    "nstate" : %i,\n' % (QMin['template']['nstate'][gsmult - 1])
        string += '    "dipole" : "%s",\n' % (QMin['template']['dipole'])
        string = string[:-2] + '\n    } ]\n  },\n\n'
    #else:
    #    string += '  {\n    "title" : "casscf",\n'
    #    for key in casscfkeys:
    #        string += '    "%s" : %s,\n' % (key, QMin['template'][key])
    #    string += '    "nspin" : %i,\n' % (gsmult - 1)
    #    string += '    "print_thresh" : 1e-10,\n'
    #    string += '    "charge" : %i,\n' % (charge)
    #    string += '    "nstate" : %i,\n' % (QMin['template']['nstate'][gsmult - 1])
    #    string = string[:-2] + '\n  },\n\n'
    #    if QMin['template']['method'] == 'caspt2':
    #        string += '  { \n  '
    #        string += '  "title" : "smith",\n'
    #        string += '    "method" : "caspt2",\n'
    #        for key in caspt2keys:
    #            string += '    "%s" : "%s",\n' % (key, QMin['template'][key])
    #        string = string[:-2] + '\n   },\n\n'



    string = string[:-3] + '\n\n'
    string += ']}'

    # oldgeom=[]
    # try:
    # f=open('previous_geom')
    # oldgeom=f.readlines()
    # f.close()
    # except IOError:
    # pass


    return string



# ======================================================================= #

def shorten_DIR(string):
    maxlen = 40
    front = 12
    if len(string) > maxlen:
        return string[0:front] + '...' + string[-(maxlen - 3 - front):]
    else:
        return string + ' ' * (maxlen - len(string))

# ======================================================================= #


def runBAGEL(WORKDIR, BAGEL, ncpu, mpi, savedir, strip=False):

    prevdir = os.getcwd()
    os.chdir(WORKDIR)
    # print(ncpu)
    if mpi:
        string = 'mpirun -n ' + str(ncpu) + ' ' + os.path.join(BAGEL, 'bin', 'BAGEL') + ' '
        os.environ['OMP_NUM_THREADS'] = '1'
        os.environ['MKL_NUM_THREADS'] = '1'
        os.environ['BAGEL_NUM_THREADS'] = '1'
        print("\nRunning with MPI (mpirun -n %i)" % (ncpu))
    else:
        string = os.path.join(BAGEL, 'bin', 'BAGEL') + ' '
        os.environ['OMP_NUM_THREADS'] = str(ncpu)
        os.environ['MKL_NUM_THREADS'] = str(ncpu)
        os.environ['BAGEL_NUM_THREADS'] = str(ncpu)
        print("\nRunning with OpenMP (%i cores)" % (ncpu))
    if PRINT or DEBUG:
        starttime = datetime.datetime.now()
        sys.stdout.write('START:\t%s\t%s\t"%s"\n' % (shorten_DIR(WORKDIR), starttime, shorten_DIR(string)))
        sys.stdout.flush()
    stdoutfile = open(os.path.join(WORKDIR, 'BAGEL.out'), 'w')
    stderrfile = open(os.path.join(WORKDIR, 'BAGEL.err'), 'w')

    string += 'BAGEL.run > BAGEL.out'
    try:
        runerror = sp.call(string, shell=True, stdout=stdoutfile, stderr=stderrfile)
    except OSError:
        print('Call have had some serious problems:', OSError)
        sys.exit(57)

    stdoutfile.close()
    stderrfile.close()
    stderr = readfile(os.path.join(WORKDIR, 'BAGEL.err'))
    for line in stderr:
        if 'error' in line.lower():
            sys.stdout.write('ERROR: \t%s\t"%s"\n' % (shorten_DIR(WORKDIR), line.strip()))
            runerror += 1
            break
    if PRINT or DEBUG:
        endtime = datetime.datetime.now()
        sys.stdout.write('FINISH:\t%s\t%s\tRuntime: %s\tError Code: %i\n' % (shorten_DIR(WORKDIR), endtime, endtime - starttime, runerror))
        sys.stdout.flush()
    if DEBUG and runerror != 0:
        copydir = os.path.join(savedir, 'debug_BAGEL_stdout')
        if not os.path.isdir(copydir):
            mkdir(copydir)
        outfile = os.path.join(WORKDIR, 'BAGEL.out')
        tofile = os.path.join(copydir, "BAGEL_problems_%s.out" % (os.path.basename(WORKDIR)))
        shutil.copy(outfile, tofile)
        print('Error in %s! Copied BAGEL output to %s' % (WORKDIR, tofile))
    os.chdir(prevdir)
    # strip=False
    if strip and not DEBUG and runerror == 0:
        stripWORKDIR(WORKDIR)
    return runerror

# ======================================================================= #


def stripWORKDIR(WORKDIR):
    ls = os.listdir(WORKDIR)
    keep = ['BAGEL.run$', 'BAGEL.err$', '.out', 'archive.archive', 'orbitals.molden']
    for ifile in ls:
        delete = True
        for k in keep:
            if containsstring(k, ifile):
                delete = False
        if delete:
            rmfile = os.path.join(WORKDIR, ifile)
            if not DEBUG:
                os.remove(rmfile)

# ======================================================================= #


def moveOldFiles(QMin):
    # moves all relevant files in the savedir to old files (per job)
    if PRINT:
        print('>>>>>>>>>>>>> Moving old files')
    basenames = ['archive', 'orbitals.molden']
    if 'nooverlap' not in QMin:
        basenames.append('mos')
    for job in QMin['joblist']:
        for base in basenames:
            fromfile = os.path.join(QMin['savedir'], '%s.%i' % (base, job))
            if not os.path.isfile(fromfile):
                print('File %s not found, cannot move to OLD!' % (fromfile))
                sys.exit(58)
            if base == 'archive':
                tofile = os.path.join(QMin['savedir'], '%s.%i.old.archive' % (base, job))
            else:
                tofile = os.path.join(QMin['savedir'], '%s.%i.old' % (base, job))
            if PRINT:
                print(shorten_DIR(fromfile) + '   =>   ' + shorten_DIR(tofile))
            shutil.copy(fromfile, tofile)
    # moves all relevant files in the savedir to old files (per mult)
    basenames = []
    if 'nooverlap' not in QMin:
        basenames = ['dets']
    for job in itmult(QMin['states']):
        for base in basenames:
            fromfile = os.path.join(QMin['savedir'], '%s.%i' % (base, job))
            if not os.path.isfile(fromfile):
                print('File %s not found, cannot move to OLD!' % (fromfile))
                sys.exit(59)
            tofile = os.path.join(QMin['savedir'], '%s.%i.old' % (base, job))
            if PRINT:
                print(shorten_DIR(fromfile) + '   =>   ' + shorten_DIR(tofile))
            shutil.copy(fromfile, tofile)
    # also remove aoovl files if present
    delete = ['AO_overl', 'AO_overl.mixed']
    for f in delete:
        rmfile = os.path.join(QMin['savedir'], f)
        if os.path.isfile(rmfile):
            os.remove(rmfile)
            if PRINT:
                print('rm ' + rmfile)
    print

# ======================================================================= #


def saveFiles(WORKDIR, QMin):

    # copy the archive and molden files from master directories
    job = QMin['IJOB']
    fromfile = os.path.join(WORKDIR, 'archive.archive')
    tofile = os.path.join(QMin['savedir'], 'archive.%i' % (job))
    # print QMin['savedir']
    shutil.copy(fromfile, tofile)
    fromfile = os.path.join(WORKDIR, 'orbitals.molden.%i' % (job))
    tofile = os.path.join(QMin['savedir'], 'orbitals.molden.%i' % (job))
    shutil.copy(fromfile, tofile)
    if PRINT:
        print(shorten_DIR(tofile))

    # if necessary, extract the MOs and write them to savedir
    if 'ion' in QMin or 'nooverlap' not in QMin:
        f = os.path.join(WORKDIR, 'orbitals.molden.%i' % job)
        string = make_mos_from_Molden(f, QMin)
        mofile = os.path.join(QMin['savedir'], 'mos.%i' % job)
        writefile(mofile, string)
        if PRINT:
            print(shorten_DIR(mofile))

    # if necessary, extract the CI coefficients and write them to savedir
    if 'ion' in QMin or 'nooverlap' not in QMin:
        f = os.path.join(WORKDIR, 'BAGEL.out')
        strings = get_dets_from_cores(f, QMin)
        for f in strings:
            writefile(f, strings[f])
            if PRINT:
                print(shorten_DIR(f))

# ======================================================================= #


def convert_decimal_format(x):
    # input: float
    # returns: formatted string
    # example: input 0.0234 -> returns 0.23400000000000D-01
    count = 0
    if x == 0.0:
        s = '%0.14fD-%02d' % (0., 0)
        return s
    else:
        if abs(x) > 1.00:
            while abs(x) >= 1.00:
                count += 1
                x /= 10
            if count == 0:
                if x < 0:
                    s = '%0.14fD+%02d' % (x, count)
                    s = '-' + s[2:]
                else:
                    s = '%0.14fD+%02d' % (x, count)
            else:
                if x < 0:
                    s = '%0.14fD+%02d' % (x, count)
                    s = '-' + s[2:]
                else:
                    s = '%0.14fD+%02d' % (x, count)
            return s
        else:
            while abs(x) <= 1.00:
                count += 1
                x *= 10
            if count == 0:
                if x < 0:
                    s = '%0.14fD-%02d' % (x / 10, count)
                    s = '-' + s[2:]
                else:
                    s = '%0.14fD-%02d' % (x / 10, count)
            else:
                if x < 0:
                    s = '%0.14fD-%02d' % (x / 10, count - 1)
                    s = '-' + s[2:]
                else:
                    s = '%0.14fD-%02d' % (x / 10, count - 1)
            return s


# ======================================================================= #
def make_mos_from_Molden(moldenfile, QMin):
   # here we take the spherical MOs and transform them to cartesian ones to get the overlaps with pyQuante

    # read file
    data = readfile(moldenfile)

    # basis info
    shells = {'s': (1, 1), 'p': (3, 3), 'd': (6, 5), 'f': (10, 7), 'g': (15, 9)}
    mode = {'s': 0, 'p': 0, 'd': 0, 'f': 0, 'g': 0}    # 0: cartesian, 1: spherical
    for line in data:
        if '[5d' in line.lower():
            mode['d'] = 1
        if '7f]' in line.lower():
            mode['f'] = 1
        if '[9g]' in line.lower():
            mode['g'] = 1
    NAOs = [0, 0]          # 0: with full cart basis, 1: as in Molden
    aos = []

    # get basis info
    for iline, line in enumerate(data):
        if '[gto]' in line.lower():
            break
    else:
        print('Could not find basis set in %s!' % (moldenfile))
        sys.exit(60)
    while True:
        iline += 1
        if iline >= len(data):
            print('Could not find basis set in %s!' % (moldenfile))
            sys.exit(61)
        line = data[iline].lower()
        if '[' in line:
            break
        s = line.split()
        if len(s) >= 1 and s[0] in shells:
            NAOs[0] += shells[s[0]][0]
            NAOs[1] += shells[s[0]][mode[s[0]]]
            aos.append(s[0])

    #print(NAOs)
    #print(aos)

    for iline, line in enumerate(data):
        if '[mo]' in line.lower():
            break
    else:
        print('Could not find MO coefficients in %s!' % (moldenfile))
        sys.exit(62)

    NAO = NAOs[1]
    # get coefficients for alpha
    NMO_A = NAO
    MO_A = [[0. for i in range(NAO)] for j in range(NMO_A)]
    for imo in range(NMO_A):
        for iao in range(NAO):
            jline = iline + 4 + iao + (NAO + 3) * imo
            line = data[jline]
            MO_A[imo][iao] = float(line.split()[1])

    if any([1 == mode[i] for i in aos]):
        print('no support for spherical basis sets in overlaps')
        sys.exit(1)

    # reorder AOs from Molden order to pyscf order
    # https://www.theochem.ru.nl/molden/molden_format.html
    # https://pyscf.org/user/gto.html#ordering-of-basis-function
    reorders = {
    's': {0:0},
    'p': {0:0, 1:1, 2:2},
    'd': {0:0, 1:3, 2:5, 3:1, 4:2, 5:4},
    'f': {0:0, 1:6, 2:9, 3:3, 4:1, 5:2, 6:4, 7:8, 8:7, 9:4}
    }
    factors = {
    #'d': {1: math.sqrt(3.), 2: math.sqrt(3.), 4: math.sqrt(3.)},
    }
    MO_A_reorder = []
    for imo in range(NMO_A):
        mo = [0. for i in range(NAO)]
        ishift=0
        for ishell,s in enumerate(aos):
            for i in range(shells[s][mode[s[0]]]):
                iao_old = ishift + i
                iao_new = ishift + reorders[s][i]
                #print(ishell,s,i,iao_old,'->',iao_new)i
                if s in factors and i in factors[s]:
                    fact=factors[s][i]
                else:
                    fact=1.
                mo[iao_new] = MO_A[imo][iao_old]*math.sqrt(fact)
            ishift += shells[s][mode[s[0]]]
        MO_A_reorder.append(mo)
    MO_A = MO_A_reorder


    #d34 = math.sqrt(3. / 4.)
    #f65 = math.sqrt(6. / 5.)
    #f38 = math.sqrt(3. / 8.)
    #f58 = math.sqrt(5. / 8.)
    #f98 = math.sqrt(9. / 8.)
    #f920 = math.sqrt(9. / 20.)
    #f340 = math.sqrt(3. / 40.)
    #transform = {
    #    'd': [[(0, -0.5), (3, d34)],
    #          [(0, -0.5), (3, -d34)],
    #          [(0, 1.0)],
    #          [(4, 1.0)],
    #          [(1, 1.0)],
    #          [(2, 1.0)]],
    #    'f': [[(1, -f38), (5, f58)],
    #          [(2, -f38), (6, -f58)],
    #          [(0, 1.0)],
    #          [(1, -f340), (5, -f98)],
    #          [(2, -f340), (6, f98)],
    #          [(0, -f920), (3, d34)],
    #          [(1, f65)],
    #          [(2, f65)],
    #          [(0, -f920), (3, -d34)],
    #          [(4, 1.0)]],
    #    'g': None                         # TODO: add and check for bagel
    #}

    #tMO_A = []
    #for imo in range(NMO_A):
    #    newmo = []
    #    iao = 0
    #    for shell in aos:
    #        n_cart = shells[shell][0]
    #        if mode[shell] == 0:
    #            for i in range(n_cart):
    #                newmo.append(MO_A[imo][iao])
    #                iao += 1
    #        else:
    #            n_sph = shells[shell][1]
    #            new = [0. for i in range(n_cart)]
    #            for i in range(n_cart):
    #                for j, c in transform[shell][i]:
    #                    new[i] += c * MO_A[imo][iao + j]
    #            newmo.extend(new)
    #            iao += n_sph
    #    tMO_A.append(newmo)


    #MO_A = tMO_A
    #NAO = NAOs[0]


    # handle frozen core
    # print QMin['frozcore']
    NMO = NMO_A - QMin['frozcore']

    string = '''2mocoef
header
 1
MO-coefficients from Gaussian
 1
 %i   %i
 a
mocoef
(*)
''' % (NAO, NMO)
    x = 0
    for imo, mo in enumerate(MO_A):
        if imo < QMin['frozcore']:
            continue
        for c in mo:
            if x >= 3:
                string += '\n'
                x = 0
            string += '% 6.12e ' % c
            x += 1
        if x > 0:
            string += '\n'
            x = 0
    string += 'orbocc\n(*)\n'
    x = 0
    for i in range(NMO):
        if x >= 3:
            string += '\n'
            x = 0
        string += '% 6.12e ' % (0.0)
        x += 1

    return string

# ======================================================================= #


def get_dets_from_cores(filename, QMin):
    # obtain ci coefficients from output
    # in the case of XMS-caspt2, transform them via the rotation matrix
    cores = readfile(filename)
    job = QMin['IJOB']
    casscf = []
    totrot = []
    xms = []
    nclosed = QMin['template']['nclosed']
    nact = QMin['template']['nact']
    mult = QMin['jobs'][job]['mults'][0]
    read = False
    readrot = False
    readxms = False
    # loop over outfile to get all casscf-ci-coeffs and rot matrices
    for line in cores:
        if 'NACME' in line and read:
            casscf.append(ci_vec)
            read = False
        if 'Permanent' in line and read:
            casscf.append(ci_vec)
            read = False
        if len(line.strip()) == 0 and read:
            casscf.append(ci_vec)
            read = False
        if len(line.strip()) == 0 and readrot:
            totrot.append(rot)
            readrot = False
        if read and 'Reference calculation' not in line:
            ci_vec[line.split()[0].replace('2', '3').replace('a', '2').replace('b', '1').replace('.', '0')] = line.split()[1]

        if readrot:
            rot.append(line.split())
        if readxms:
            a += 1
            if a > 2:
                xms.append(line.split()[1:4])
                if a > 4:
                    readxms = False
        if 'ci vector' in line:
            ci_vec = {}
            read = True
        if 'MS-CASPT2 rotation matrix' in line:
            rot = []
            readrot = True
        if '* nvirt' in line:
            nvirt = int(line.split()[3])
        if 'Extended multi-state CASPT2 (XMS-CASPT2) rotation matrix' in line:
            readxms = True
            a = 0

    # rotate states to get caspt2-ci-coeffs
    if QMin['template']['method'] == 'caspt2' and len(totrot) > 0:
        # print("rot start")
        for j in range(len(totrot)):
            mats = [totrot[j]]
            newstates = []
            for k in range(len(mats)):
                rotmat = mats[k]
                for i in range(len(rotmat)):
                    ci_vec = {}
                    for j in range(len(rotmat)):
                        if abs(float(rotmat[j][i])) > 0.0:
                            for config in casscf[j]:
                                if config in ci_vec:
                                    ci_vec[config] += float(casscf[j][config]) * float(rotmat[j][i])
                                else:
                                    ci_vec[config] = float(casscf[j][config]) * float(rotmat[j][i])
                    newstates.append(ci_vec)
        casscf = deepcopy(newstates)
#        for m in casscf:
#          cisum=0
#          for entry in m:
#            cisum+=m[entry]**2
#          print(cisum)
    else:
        newstates = casscf

    # rearrange dets
    used_dets = {}
    for i in range(len(newstates)):
        for config in newstates[i]:
            if config in used_dets:
                used_dets[config][i] = float(newstates[i][config])
            else:
                used_dets[config] = [0 for k in range(len(newstates))]
                used_dets[config][i] = float(newstates[i][config])


    nclosed -= QMin['frozcore']

    # get number of alpha electrons.
    for config in used_dets:
        na = 0
        ninv = 0
        for char in config:
            if char == '1' or char == '3':
                na += 1
        break  # break after first configuration
    # construct string for det file
    nalpha = na + nclosed
    # print("hupdidupi", nclosed, nact, nvirt)
    s = '%i %i %i\n' % (len(newstates), nclosed + nact + nvirt, len(used_dets))

    # adapt sign of the configurations to be consistent with the wfoverlap convention
    for entry in used_dets:
        detstring = nclosed * 'd' + '%s' % (entry.replace('0', 'e').replace('1', 'a').replace('2', 'b').replace('3', 'd')) + nvirt * 'e' + ' '
        s += detstring
        na = 0
        ninv = 0
        for char in detstring:
            if char == 'a' or char == 'd':
                na += 1
            if char == 'b' or char == 'd':
                ninv += nalpha - na
        for contribution in used_dets[entry]:
            s += '%11.7f ' % (contribution * (-1)**ninv)
        s += '\n'
#    print s

    strings = {}
    filename = os.path.join(QMin['savedir'], 'dets.%i' % mult)
    strings[filename] = s

    return strings


# ======================================================================= #
def mkdir(DIR):
    # mkdir the DIR, or clean it if it exists
    if os.path.exists(DIR):
        if os.path.isfile(DIR):
            print('%s exists and is a file!' % (DIR))
            sys.exit(63)
        elif os.path.isdir(DIR):
            if DEBUG:
                print('Remake\t%s' % DIR)
            shutil.rmtree(DIR)
            os.makedirs(DIR)
    else:
        try:
            if DEBUG:
                print('Make\t%s' % DIR)
            os.makedirs(DIR)
        except OSError:
            print('Can not create %s\n' % (DIR))
            sys.exit(64)

# ======================================================================= #


def link(PATH, NAME, crucial=True, force=True):
    # do not create broken links
    if not os.path.exists(PATH) and crucial:
        print('Source %s does not exist, cannot create link!' % (PATH))
        sys.exit(65)
    if os.path.islink(NAME):
        if not os.path.exists(NAME):
            # NAME is a broken link, remove it so that a new link can be made
            os.remove(NAME)
        else:
            # NAME is a symlink pointing to a valid file
            if force:
                # remove the link if forced to
                os.remove(NAME)
            else:
                print('%s exists, cannot create a link of the same name!' % (NAME))
                if crucial:
                    sys.exit(66)
                else:
                    return
    elif os.path.exists(NAME):
        # NAME is not a link. The interface will not overwrite files/directories with links, even with force=True
        print('%s exists, cannot create a link of the same name!' % (NAME))
        if crucial:
            sys.exit(67)
        else:
            return
    os.symlink(PATH, NAME)

# =============================================================================================== #
# =============================================================================================== #
# =======================================  TheoDORE ============================================= #
# =============================================================================================== #
# =============================================================================================== #


def run_theodore(QMin, errorcodes):

    if 'theodore' in QMin:
        print('>>>>>>>>>>>>> Starting the TheoDORE job execution')

        for ijob in QMin['jobs']:
            if not QMin['jobs'][ijob]['restr']:
                if DEBUG:
                    print('Skipping Job %s because it is unrestricted.' % (ijob))
                continue
            WORKDIR = os.path.join(QMin['scratchdir'], 'master_%i' % ijob)
            setupWORKDIR_TH(WORKDIR, QMin)
            os.environ
            errorcodes['theodore_%i' % ijob] = runTHEODORE(WORKDIR, QMin['theodir'])

        # Error code handling
        j = 0
        string = 'Error Codes:\n'
        for i in errorcodes:
            if 'theodore' in i:
                string += '\t%s\t%i' % (i + ' ' * (10 - len(i)), errorcodes[i])
                j += 1
                if j == 4:
                    j = 0
                    string += '\n'
        print(string)
        if any((i != 0 for i in errorcodes.values())):
            print('Some subprocesses did not finish successfully!')
            sys.exit(68)

        print('')

    return errorcodes

# ======================================================================= #


def setupWORKDIR_TH(WORKDIR, QMin):
    # mkdir the WORKDIR, or clean it if it exists, then copy all necessary files from pwd and savedir

    # write dens_ana.in
    inputstring = '''rtype='BAGEL'
rfile='TAPE21'
jmol_orbitals=False
molden_orbitals=False
Om_formula=2
eh_pop=1
comp_ntos=True
print_OmFrag=True
output_file='tden_summ.txt'
prop_list=%s
at_lists=%s
''' % (str(QMin['template']['theodore_prop']), str(QMin['template']['theodore_fragment']))

    filename = os.path.join(WORKDIR, 'dens_ana.in')
    writefile(filename, inputstring)
    if DEBUG:
        print('================== DEBUG input file for WORKDIR %s =================' % (shorten_DIR(WORKDIR)))
        print(inputstring)
        print('TheoDORE input written to: %s' % (filename))
        print('====================================================================')

    return

# ======================================================================= #


def runTHEODORE(WORKDIR, THEODIR):
    prevdir = os.getcwd()
    os.chdir(WORKDIR)
    string = 'python ' + os.path.join(THEODIR, 'bin', 'analyze_tden.py')
    stdoutfile = open(os.path.join(WORKDIR, 'theodore.out'), 'w')
    stderrfile = open(os.path.join(WORKDIR, 'theodore.err'), 'w')
    if PRINT or DEBUG:
        starttime = datetime.datetime.now()
        sys.stdout.write('START:\t%s\t%s\t"%s"\n' % (shorten_DIR(WORKDIR), starttime, shorten_DIR(string)))
        sys.stdout.flush()
    try:
        runerror = sp.call(string, shell=True, stdout=stdoutfile, stderr=stderrfile)
    except OSError:
        print('Call have had some serious problems:', OSError)
        sys.exit(69)
    stdoutfile.close()
    stderrfile.close()
    if PRINT or DEBUG:
        endtime = datetime.datetime.now()
        sys.stdout.write('FINISH:\t%s\t%s\tRuntime: %s\tError Code: %i\n' % (shorten_DIR(WORKDIR), endtime, endtime - starttime, runerror))
        sys.stdout.flush()
    os.chdir(prevdir)
    return runerror

# =============================================================================================== #
# =============================================================================================== #
# =======================================  Dyson and overlap calcs ============================== #
# =============================================================================================== #
# =============================================================================================== #


def run_wfoverlap(QMin, errorcodes):

    print('>>>>>>>>>>>>> Starting the WFOVERLAP job execution')

    # do Dyson calculations
    if 'ion' in QMin:
        for ionpair in QMin['ionmap']:
            WORKDIR = os.path.join(QMin['scratchdir'], 'Dyson_%i_%i_%i_%i' % ionpair)
            files = {'aoovl': 'AO_overl',
                     'det.a': 'dets.%i' % ionpair[0],
                     'det.b': 'dets.%i' % ionpair[2],
                     'mo.a': 'mos.%i' % ionpair[1],
                     'mo.b': 'mos.%i' % ionpair[3]}
            setupWORKDIR_WF(WORKDIR, QMin, files)
            errorcodes['Dyson_%i_%i_%i_%i' % ionpair] = runWFOVERLAP(WORKDIR, QMin['wfoverlap'], memory=QMin['memory'], ncpu=QMin['ncpu'])

    # do overlap calculations
    if 'overlap' in QMin:
        get_Double_AOovl_molden(QMin)
        for m in itmult(QMin['states']):
            job = QMin['multmap'][m]
            WORKDIR = os.path.join(QMin['scratchdir'], 'WFOVL_%i_%i' % (m, job))
            files = {'aoovl': 'AO_overl.mixed',
                     'det.a': 'dets.%i.old' % m,
                     'det.b': 'dets.%i' % m,
                     'mo.a': 'mos.%i.old' % job,
                     'mo.b': 'mos.%i' % job}
            setupWORKDIR_WF(WORKDIR, QMin, files)
            errorcodes['WFOVL_%i_%i' % (m, job)] = runWFOVERLAP(WORKDIR, QMin['wfoverlap'], memory=QMin['memory'], ncpu=QMin['ncpu'])

    # Error code handling
    j = 0
    string = 'Error Codes:\n'
    for i in errorcodes:
        if 'Dyson' in i or 'WFOVL' in i:
            string += '\t%s\t%i' % (i + ' ' * (10 - len(i)), errorcodes[i])
            j += 1
            if j == 4:
                j = 0
                string += '\n'
    print(string)
    if any((i != 0 for i in errorcodes.values())):
        print('Some subprocesses did not finish successfully!')
        sys.exit(70)

    print('')

    return errorcodes

# ======================================================================= #


def setupWORKDIR_WF(WORKDIR, QMin, files):
    # mkdir the WORKDIR, or clean it if it exists, then copy all necessary files from pwd and savedir

    # setup the directory
    mkdir(WORKDIR)

    # write wfovl.inp
    inputstring = '''mix_aoovl=aoovl
a_mo=mo.a
b_mo=mo.b
a_det=det.a
b_det=det.b
a_mo_read=0
b_mo_read=0
ao_read=0
'''
    if 'ion' in QMin:
        if QMin['ndocc'] > 0:
            inputstring += 'ndocc=%i\n' % (QMin['ndocc'])
    if QMin['ncpu'] >= 8:
        inputstring += 'force_direct_dets\n'
    filename = os.path.join(WORKDIR, 'wfovl.inp')
    writefile(filename, inputstring)
    if DEBUG:
        print('================== DEBUG input file for WORKDIR %s =================' % (shorten_DIR(WORKDIR)))
        print(inputstring)
        print('wfoverlap input written to: %s' % (filename))
        print('====================================================================')

    # link input files from save
    linkfiles = ['aoovl', 'det.a', 'det.b', 'mo.a', 'mo.b']
    for f in linkfiles:
        fromfile = os.path.join(QMin['savedir'], files[f])
        tofile = os.path.join(WORKDIR, f)
        link(fromfile, tofile)

    return

# ======================================================================= #


def runWFOVERLAP(WORKDIR, WFOVERLAP, memory=100, ncpu=1):
    prevdir = os.getcwd()
    os.chdir(WORKDIR)
    string = WFOVERLAP + ' -m %i' % (memory) + ' -f wfovl.inp'
    stdoutfile = open(os.path.join(WORKDIR, 'wfovl.out'), 'w')
    stderrfile = open(os.path.join(WORKDIR, 'wfovl.err'), 'w')
    os.environ['OMP_NUM_THREADS'] = str(ncpu)
    if PRINT or DEBUG:
        starttime = datetime.datetime.now()
        sys.stdout.write('START:\t%s\t%s\t"%s"\n' % (shorten_DIR(WORKDIR), starttime, shorten_DIR(string)))
        sys.stdout.flush()
    try:
        runerror = sp.call(string, shell=True, stdout=stdoutfile, stderr=stderrfile)
    except OSError:
        print('Call have had some serious problems:', OSError)
        sys.exit(71)
    stdoutfile.close()
    stderrfile.close()
    if PRINT or DEBUG:
        endtime = datetime.datetime.now()
        sys.stdout.write('FINISH:\t%s\t%s\tRuntime: %s\tError Code: %i\n' % (shorten_DIR(WORKDIR), endtime, endtime - starttime, runerror))
        sys.stdout.flush()
    os.chdir(prevdir)
    return runerror

# ======================================================================= #


def get_Double_AOovl_molden(QMin):

    for m in itmult(QMin['states']):
        if QMin['states'] != 0:
            # selection of which multiplet does not matter for ao overlap
            filename1 = os.path.join(QMin['savedir'], 'orbitals.molden.%i.old' % m)
            filename2 = os.path.join(QMin['savedir'], 'orbitals.molden.%i' % m)

    #
    NAO, Smat = get_smat_from_Molden(filename1, filename2)

    # Smat is now full matrix NAO*NAO
    # we want the lower left quarter, but transposed
    string = '%i %i\n' % (NAO // 2, NAO // 2)
    for irow in range(NAO // 2, NAO):
        for icol in range(0, NAO // 2):
            string += '% .15e ' % (Smat[icol][irow])  # note the exchanged indices => transposition
        string += '\n'
    filename = os.path.join(QMin['savedir'], 'AO_overl.mixed')
    writefile(filename, string)
    return

# ======================================================================= #


def get_smat_from_Molden(file1, file2=''):

    try:
        from pyscf.gto import mole
        from pyscf.tools.molden import load
    except ImportError:
        print('Could not import pyscf!')
        sys.exit(72)
    # read file1
    mol, mo_energy, mo_coeff, mo_occ, irrep_labels, spins = load(file1)
    # read file2:
    mol2 = None
    if file2:
        mol2, mo_energy, mo_coeff, mo_occ, irrep_labels, spins = load(file2)
        mol = mole.conc_mol(mol, mol2)

    S = mol.intor('int1e_ovlp') #.tolist()
    
    # normalize with sqrt(diag)
    for i in range(len(S)):
       n=math.sqrt(S[i,i])
       S[i,:] /= n
       S[:,i] /= n

    return len(S), S.tolist()


# ======================================================================= #

def saveAOmatrix(WORKDIR, QMin, job):
    for i in range(5):
        test = 'master_%i' % i
        if test in job:
            break
    filename = os.path.join(WORKDIR, 'orbitals.molden.%i' % i)
    NAO, Smat = get_smat_from_Molden(filename)

    string = '%i %i\n' % (NAO, NAO)
    for irow in range(NAO):
        for icol in range(NAO):
            string += '% .7e ' % (Smat[icol][irow])
        string += '\n'
    filename = os.path.join(QMin['savedir'], 'AO_overl')
    writefile(filename, string)
    if PRINT:
        print(shorten_DIR(filename))



# ======================================================================= #



# =============================================================================================== #
# =============================================================================================== #
# ========================================= BAGEL output parsing ================================== #
# =============================================================================================== #
# =============================================================================================== #

def getQMout(QMin):

    if PRINT:
        print('>>>>>>>>>>>>> Reading output files')
    starttime = datetime.datetime.now()

    QMout = {}
    states = QMin['states']
    nstates = QMin['nstates']
    nmstates = QMin['nmstates']
    natom = QMin['natom']
    joblist = QMin['joblist']

    # Hamiltonian
    if 'h' in QMin or 'soc' in QMin:
        # make Hamiltonian
        if 'h' not in QMout:
            QMout['h'] = makecmatrix(nmstates, nmstates)
        # go through all jobs
        for job in joblist:
            # first get energies from TAPE21
            corefile = os.path.join(QMin['scratchdir'], 'master_%i/BAGEL.out' % (job))
            energies = getenergy(corefile, job, QMin)
            mults = QMin['multmap'][-job]

            for i in range(nmstates):
                for j in range(nmstates):
                    m1, s1, ms1 = tuple(QMin['statemap'][i + 1])
                    m2, s2, ms2 = tuple(QMin['statemap'][j + 1])
                    if m1 not in mults or m2 not in mults:
                        continue
                    if i == j:
                        QMout['h'][i][j] = energies[(m1, s1)]
                    elif 'soc' in QMin and QMin['jobs'][job]['restr']:
                        if m1 == m2 == 1:
                            continue
                        # TODO: invstatemap and submatrix are not defined
                        x = invstatemap[(m1, s1, ms1)]
                        y = invstatemap[(m2, s2, ms2)]
                        QMout['h'][i][j] = submatrix[x - 1][y - 1]

    # Dipole Moments
    if 'dm' in QMin:
        # make matrix
        if 'dm' not in QMout:
            QMout['dm'] = [makecmatrix(nmstates, nmstates) for i in range(3)]
        complete = {}
        # get permanent and transition dipole moments from all master files
        for job in joblist:
            dipoles = getdm(QMin, job)
            complete[job] = dipoles
        # now we use the fact, that we know the order of the calculated dm's
        # get the order of calculated properties for each multiplicity
        gradmap = QMin['gradmap']
        nacmap = QMin['nacmap']
        for job in joblist:
            jobdm = []
            for grad in gradmap:
                if grad[0] == job:
                    jobdm.append((grad[1], grad[1]))
            for nac in nacmap:
                if nac[0] == job:
                    jobdm.append((nac[1], nac[3]))
            # loop over all states, find corresponding entry in calculation
            # to map the entry from above to the extracted dm's
            for i in range(nmstates):
                m1, s1, ms1 = tuple(QMin['statemap'][i + 1])
                if m1 != job:
                    continue
                for j in range(nmstates):
                    m2, s2, ms2 = tuple(QMin['statemap'][j + 1])
                    if m2 != job:
                        continue
                    if (s1, s2) in jobdm:
                        if m1 == m2 and ms1 == ms2:
                            for ixyz in range(3):
                                QMout['dm'][ixyz][i][j] = complete[job][jobdm.index((s1, s2))][ixyz]
                                QMout['dm'][ixyz][j][i] = complete[job][jobdm.index((s1, s2))][ixyz]


    # Gradients
    if 'grad' in QMin:
        QMout['grad'] = [[[0. for i in range(3)] for j in range(natom)] for k in range(nmstates)]
        for job in joblist:
            outfile = os.path.join(QMin['scratchdir'], 'master_%i/BAGEL.out' % (job))
            grads = getgrad(outfile, QMin, job)
            for istate in QMin['statemap']:
                state = QMin['statemap'][istate]
                if state[0] == job:
                    QMout['grad'][istate - 1] = grads[state[1] - 1]

        if QMin['neglected_gradient'] != 'zero':
            for i in range(nmstates):
                m1, s1, ms1 = tuple(QMin['statemap'][i + 1])
                if not (m1, s1) in QMin['gradmap']:
                    if QMin['neglected_gradient'] == 'gs':
                        j = QMin['gsmap'][i + 1] - 1
                    elif QMin['neglected_gradient'] == 'closest':
                        e1 = QMout['h'][i][i]
                        de = 999.
                        for grad in QMin['gradmap']:
                            for k in range(nmstates):
                                m2, s2, ms2 = tuple(QMin['statemap'][k + 1])
                                if grad == (m2, s2):
                                    break
                            e2 = QMout['h'][k][k]
                            if de > abs(e1 - e2):
                                de = abs(e1 - e2)
                                j = k
                    QMout['grad'][i] = QMout['grad'][j]

    # NACdr
    if 'nacdr' in QMin:
        if 'nacdr' not in QMout:
            nacdr = [[[[0. for i in range(3)] for j in range(natom)] for k in range(nmstates)] for l in range(nmstates)]
        nacs = {}
        for job in joblist:
            jobnac = getnacdr(job, QMin)
            nacs.update(jobnac)
        for i, i1 in enumerate(QMin['statemap']):
            mult, state, ms = tuple(QMin['statemap'][i1])
            for j, j1 in enumerate(QMin['statemap']):
                if not i < j:
                    continue
                mult2, state2, ms2 = tuple(QMin['statemap'][j1])
                test = (mult, state - 1, mult2, state2 - 1)
                if not ms1 == ms2:
                    continue
                if test in nacs:
                    # TODO: check if "full" option is used as nacmtype
                    # TODO: seems to be ok without unweighting
                    # unweight the energy gap
                    # dE=QMout['h'][i][i]-QMout['h'][j][j]
                    dE = 1.
                    nacdr[i][j] = deepcopy(nacs[test])
                    nacdr[j][i] = deepcopy(nacdr[i][j])
                    for x in range(natom):
                        for y in range(3):
                            nacdr[i][j][x][y] *= +dE
                            nacdr[j][i][x][y] *= -dE
        QMout['nacdr'] = nacdr

    # Regular Overlaps
    if 'overlap' in QMin:
        if 'overlap' not in QMout:
            QMout['overlap'] = makecmatrix(nmstates, nmstates)
        for mult in itmult(QMin['states']):
            job = QMin['multmap'][mult]
            outfile = os.path.join(QMin['scratchdir'], 'WFOVL_%i_%i/wfovl.out' % (mult, job))
            out = readfile(outfile)
            if PRINT:
                print('Overlaps: ' + shorten_DIR(outfile))
            for i in range(nmstates):
                for j in range(nmstates):
                    m1, s1, ms1 = tuple(QMin['statemap'][i + 1])
                    m2, s2, ms2 = tuple(QMin['statemap'][j + 1])
                    if not m1 == m2 == mult:
                        continue
                    if not ms1 == ms2:
                        continue
                    QMout['overlap'][i][j] = getsmate(out, s1, s2)

    # Phases from overlaps
    if 'phases' in QMin:
        if 'phases' not in QMout:
            QMout['phases'] = [complex(1., 0.) for i in range(nmstates)]
        if 'overlap' in QMout:
            for i in range(nmstates):
                if QMout['overlap'][i][i].real < 0.:
                    QMout['phases'][i] = complex(-1., 0.)

    # Dyson norms
    if 'ion' in QMin:
        if 'prop' not in QMout:
            QMout['prop'] = makecmatrix(nmstates, nmstates)
        for ion in QMin['ionmap']:
            outfile = os.path.join(QMin['scratchdir'], 'Dyson_%i_%i_%i_%i/wfovl.out' % ion)
            out = readfile(outfile)
            if PRINT:
                print('Dyson:    ' + shorten_DIR(outfile))
            for i in range(nmstates):
                for j in range(nmstates):
                    m1, s1, ms1 = tuple(QMin['statemap'][i + 1])
                    m2, s2, ms2 = tuple(QMin['statemap'][j + 1])
                    if not (ion[0], ion[2]) == (m1, m2) and not (ion[0], ion[2]) == (m2, m1):
                        continue
                    if not abs(ms1 - ms2) == 0.5:
                        continue
                    # switch multiplicities such that m1 is smaller mult
                    if m1 > m2:
                        s1, s2 = s2, s1
                        m1, m2 = m2, m1
                        ms1, ms2 = ms2, ms1
                    # compute M_S overlap factor
                    if ms1 < ms2:
                        factor = (ms1 + 1. + (m1 - 1.) / 2.) / m1
                    else:
                        factor = (-ms1 + 1. + (m1 - 1.) / 2.) / m1
                    QMout['prop'][i][j] = getDyson(out, s1, s2) * factor

    # TheoDORE
    if 'theodore' in QMin:
        if 'theodore' not in QMout:
            QMout['theodore'] = makecmatrix(QMin['template']['theodore_n'], nmstates)
        for job in joblist:
            if not QMin['jobs'][job]['restr']:
                continue
            sumfile = os.path.join(QMin['scratchdir'], 'master_%i/tden_summ.txt' % job)
            omffile = os.path.join(QMin['scratchdir'], 'master_%i/OmFrag.txt' % job)
            props = get_theodore(sumfile, omffile, QMin)
            for i in range(nmstates):
                m1, s1, ms1 = tuple(QMin['statemap'][i + 1])
                if (m1, s1) in props:
                    for j in range(QMin['template']['theodore_n']):
                        QMout['theodore'][i][j] = props[(m1, s1)][j]

    # QM/MM energy terms
#    if QMin['template']['qmmm']:
#        job=QMin['joblist'][0]
#        outfile=os.path.join(QMin['scratchdir'],'master_%i/BAGEL.out' % (job))
#        QMout['qmmm_energies']=get_qmmm_energies(outfile,QMin['template']['qmmm_coupling'])

    endtime = datetime.datetime.now()
    if PRINT:
        print("Readout Runtime: %s" % (endtime - starttime))

    if DEBUG:
        copydir = os.path.join(QMin['savedir'], 'debug_BAGEL_stdout')
        if not os.path.isdir(copydir):
            mkdir(copydir)
        for job in joblist:
            outfile = os.path.join(QMin['scratchdir'], 'master_%i/BAGEL.out' % (job))
            shutil.copy(outfile, os.path.join(copydir, "BAGEL_%i.out" % job))
            if QMin['jobs'][job]['restr'] and 'theodore' in QMin:
                outfile = os.path.join(QMin['scratchdir'], 'master_%i/tden_summ.txt' % job)
                shutil.copy(outfile, os.path.join(copydir, 'THEO_%i.out' % (job)))
                outfile = os.path.join(QMin['scratchdir'], 'master_%i/OmFrag.txt' % job)
                shutil.copy(outfile, os.path.join(copydir, 'THEO_OMF_%i.out' % (job)))
        if 'grad' in QMin:
            for grad in QMin['gradmap']:
                # print QMin['jobgrad'][grad]
                # path,isgs=QMin['jobgrad'][grad]
                # outfile=os.path.join(QMin['scratchdir'],path,'BAGEL.out')
                outfile = os.path.join(QMin['scratchdir'], 'master_%i/BAGEL.out' % (grad[0]))
                shutil.copy(outfile, os.path.join(copydir, "BAGEL.out"))
        if 'overlap' in QMin:
            for mult in itmult(QMin['states']):
                job = QMin['multmap'][mult]
                outfile = os.path.join(QMin['scratchdir'], 'WFOVL_%i_%i/wfovl.out' % (mult, job))
                shutil.copy(outfile, os.path.join(copydir, 'WFOVL_%i_%i.out' % (mult, job)))
        if 'ion' in QMin:
            for ion in QMin['ionmap']:
                outfile = os.path.join(QMin['scratchdir'], 'Dyson_%i_%i_%i_%i/wfovl.out' % ion)
                shutil.copy(outfile, os.path.join(copydir, 'Dyson_%i_%i_%i_%i.out' % ion))

    return QMout

# ======================================================================= #


def getenergy(corefile, ijob, QMin):


    f = readfile(corefile)
    if PRINT:
        print('Energy:   ' + shorten_DIR(corefile))

    # figure out the excited state settings
    mults = QMin['jobs'][ijob]['mults']
    restr = QMin['jobs'][ijob]['restr']
    gsmult = mults[0]
    estates_to_extract = deepcopy(QMin['states'])
    for imult in range(len(estates_to_extract)):
        if not imult + 1 in mults:
            estates_to_extract[imult] = 0

    if QMin['template']['method'] == 'caspt2':
        if QMin['template']['ms'] == "false" and QMin['template']['xms'] == "false" or estates_to_extract[QMin['jobs'][ijob]['mults'][0] - 1] == 1:
            caspt2 = True
            mscaspt2 = False
        else:
            caspt2 = False
            mscaspt2 = True
    else:
        caspt2 = False
        mscaspt2 = False
    for imult in mults:
        nstates = estates_to_extract[imult - 1]
        if nstates > 0:
            for iline, line in enumerate(f):
                if mscaspt2:
                    if '* MS-CASPT2 energy :' in line:
                        if int(line.split()[5]) == 0:
                            gsenergy = float(line.split()[6])
                            energies = {(gsmult, 1): gsenergy}
                        else:
                            i = int(line.split()[5])
                            energies[(imult, i + (gsmult == imult))] = float(line.split()[6])

                if '* CASPT2 energy :' in line and caspt2:
                    if int(line.split()[5]) == 0:
                        gsenergy = float(line.split()[6])
                        energies = {(gsmult, 1): gsenergy}
                    else:
                        i = int(line.split()[5])
                        energies[(imult, i + (gsmult == imult))] = float(line.split()[6])

                elif not mscaspt2 and not caspt2:
                    if '* ci vector' in line and ' 0,' in line:
                        for i in range(nstates):
                            if i == 0:
                                # print("line that failes: ", line)
                                if '*' in f[iline - nstates - 1 + i]:
                                    energies = {(gsmult, 1): float(f[iline - nstates - 1 + i].split()[3])}
                                else:
                                    energies = {(gsmult, 1): float(f[iline - nstates - 1 + i].split()[2])}
                            elif '*' in f[iline - nstates - 1 + i]:
                                energies[(imult, i + (gsmult == imult))] = float(f[iline - nstates - 1 + i].split()[3])
                            else:
                                energies[(imult, i + (gsmult == imult))] = float(f[iline - nstates - 1 + i].split()[2])


    return energies


# ======================================================================= #
def getdm(QMin, job):

    cores = os.path.join(QMin['scratchdir'], 'master_%i/BAGEL.out' % (job))
    if PRINT:
        print('Dipoles:  ' + shorten_DIR(cores))
    out = readfile(cores)
    if QMin['template']['method'] == 'casscf':
        getstring = 'Permanent dipole moment: Relaxed'
    elif QMin['template']['method'] == 'caspt2' and QMin['dipolelevel'] == 0:
        getstring = 'Permanent dipole moment: CASPT2 relaxed'
    elif QMin['template']['method'] == 'caspt2' and QMin['dipolelevel'] == 1:
        getstring = 'Permanent dipole moment: Transition dipole moment between'
    dipoles = []
    for iline, line in enumerate(out):
        if getstring in line:
            dipole = out[iline + 1]
            dipole = dipole[dipole.find('(') + 1:dipole.find(')')].replace(',', '').split()
            for xyz in range(len(dipole)):
                dipole[xyz] = float(dipole[xyz])
            dipoles.append(dipole)
    # if only one state is requested in casscf, the dipole moment is not printed in the bagel output!
    if QMin['template']['method'] == 'casscf' and QMin['states'][job - 1] == 1:
        dipoles = [[0.0, 0.0, 0.0]]
        print("\nBAGEL does not compute Dipole moments in the case of a single state casscf calculation! The dipole moment was set to zero!\n")

    return dipoles

# ======================================================================= #


def getgrad(outfile, QMin, job):

    # read file
    out = readfile(outfile)
    collect = True
    natom = QMin['natom']
    grads = []
    g = []
    for iline, line in enumerate(out):
        if collect and 'Atom' in line:
            coord = [0, 0, 0]
            for xyz in range(3):
                coord[xyz] = float(out[iline + xyz + 1].split()[1])
            g.append(coord)
            if len(g) == natom:
                grads.append(g)
                g = []
        if 'NACME' in line:
            break
    finalgrads = [[[0. for i in range(3)] for j in range(natom)] for k in range(QMin['nmstates'])]
    i = 0
    for grad in QMin['gradmap']:
        if grad[0] == job:
            finalgrads[grad[1] - 1] = grads[i]
            i += 1
    grads = finalgrads

    return grads

# ======================================================================= #


def getnacdr(mult, QMin):
    outfile = os.path.join(QMin['scratchdir'], 'master_%i/BAGEL.out' % (mult))
    out = readfile(outfile)
    natom = QMin['natom']
    collect = False
    nacs = {}
    n = []
    for iline, line in enumerate(out):
        if 'NACME Target' in line:
            statepair = [int(out[iline].split()[4]), int(out[iline].split()[6])]
            collect = True
        if collect and 'Atom' in line:
            coord = [0, 0, 0]
            for xyz in range(3):
                coord[xyz] = float(out[iline + xyz + 1].split()[1])
            n.append(coord)
            if len(n) == natom:
                nacs[mult, statepair[0], mult, statepair[1]] = n
                n = []
        if 'Gradient computed' in line:
            collect = False

    return nacs


# ======================================================================= #
# def get_qmmm_energies(outfile,coupling):

    # out=readfile(outfile)
    # if PRINT:
    # print('QMMM:     '+shorten_DIR(outfile))

    # if coupling==1:
    #startstring='Q M / M M      E N E R G Y'
    # shift=2
    # elif coupling==2:
    #startstring='These results include the electrostatic interaction between QM and MM systems'
    # shift=0

    # toextract={'bond_mm':     (5,2),
    # 'angle_mm':    (6,2),
    # 'torsion_mm':  (7,2),
    # 'VdW_mm':      (10,4),
    # 'elstat_mm':   (11,2),
    # 'VdW_qmmm':    (10,5),
    # 'elstat_qmmm': (11,3)
    # }

    # iline=-1
    # while True:
    # iline+=1
    # if startstring in out[iline]:
    # break
    # iline+=shift

    # energies={}
    # for label in toextract:
    # t=toextract[label]
    # e=float(out[iline+t[0]].split()[t[1]])
    # energies[label]=e

    # return energies

# ======================================================================= #
def getsmate(out, s1, s2):
    ilines = -1
    while True:
        ilines += 1
        if ilines == len(out):
            print('Overlap of states %i - %i not found!' % (s1, s2))
            sys.exit(75)
        if containsstring('Overlap matrix <PsiA_i|PsiB_j>', out[ilines]):
            break
    ilines += 1 + s1
    f = out[ilines].split()
    return float(f[s2 + 1])

# ======================================================================= #


def getDyson(out, s1, s2):
    ilines = -1
    while True:
        ilines += 1
        if ilines == len(out):
            print('Dyson norm of states %i - %i not found!' % (s1, s2))
            sys.exit(76)
        if containsstring('Dyson norm matrix <PsiA_i|PsiB_j>', out[ilines]):
            break
    ilines += 1 + s1
    f = out[ilines].split()
    return float(f[s2 + 1])

# ======================================================================= #


def get_theodore(sumfile, omffile, QMin):
    out = readfile(sumfile)
    if PRINT:
        print('TheoDORE: ' + shorten_DIR(sumfile))
    props = {}
    for line in out[2:]:
        s = line.replace('(', ' ').replace(')', ' ').split()
        if len(s) == 0:
            continue
        n = int(s[0])
        m = int(s[1])
        props[(m, n + (m == 1))] = [theo_float(i) for i in s[5:]]

    out = readfile(omffile)
    if PRINT:
        print('TheoDORE: ' + shorten_DIR(omffile))
    for line in out[1:]:
        s = line.replace('(', ' ').replace(')', ' ').split()
        if len(s) == 0:
            continue
        n = int(s[0])
        m = int(s[1])
        props[(m, n + (m == 1))].extend([theo_float(i) for i in s[4:]])

    return props

# ======================================================================= #


def theo_float(string):
    try:
        s = float(string)
    except ValueError:
        s = 0.
    return s

# =============================================================================================== #
# =============================================================================================== #
# ========================================= Miscellaneous ======================================= #
# =============================================================================================== #
# =============================================================================================== #

# ======================================================================= #


def cleandir(directory):
    for data in os.listdir(directory):
        path = directory + '/' + data
        if os.path.isfile(path) or os.path.islink(path):
            if DEBUG:
                print('rm %s' % (path))
            try:
                os.remove(path)
            except OSError:
                print('Could not remove file from directory: %s' % (path))
        else:
            if DEBUG:
                print('')
            cleandir(path)
            os.rmdir(path)
            if DEBUG:
                print('rm %s' % (path))
    if PRINT:
        print('===> Cleaning up directory %s' % (directory))

# ======================================================================= #


def backupdata(backupdir, QMin):
    # save all files in savedir, except which have 'old' in their name
    ls = os.listdir(QMin['savedir'])
    for f in ls:
        ff = QMin['savedir'] + '/' + f
        if os.path.isfile(ff) and 'old' not in ff:
            step = int(QMin['step'][0])
            fdest = backupdir + '/' + f + '.stp' + str(step)
            shutil.copy(ff, fdest)

# =============================================================================================== #
# =============================================================================================== #
# ========================================= Main ================================================ #
# =============================================================================================== #
# =============================================================================================== #


def main():


    # Retrieve PRINT and DEBUG #delete?
    try:
        envPRINT = os.getenv('SH2BAGEL_PRINT')
        if envPRINT and envPRINT.lower() == 'false':
            global PRINT
            PRINT = False
        envDEBUG = os.getenv('SH2BAGEL_DEBUG')
        if envDEBUG and envDEBUG.lower() == 'true':
            global DEBUG
            DEBUG = True
    except ValueError:
        print('PRINT or DEBUG environment variables do not evaluate to numerical values!')
        sys.exit(77)

    # Process Command line arguments
    if len(sys.argv) != 2:
        print('Usage:\n./SHARC_BAGEL.py <QMin>\n')
        print('version:', version)
        print('date:', versiondate)
        print('changelog:\n', changelogstring)
        sys.exit(78)
    QMinfilename = sys.argv[1]

    # Print header
    printheader()


    # Read QMinfile
    QMin = readQMin(QMinfilename)
    # get the job schedule
    QMin, schedule = generate_joblist(QMin)

    printQMin(QMin)
    if DEBUG:
        pprint.pprint(schedule, depth=1)

    # run all the BAGEL jobs
    errorcodes = runjobs(schedule, QMin)

    # do all necessary overlap and Dyson calculations
    if 'ion' in QMin or 'overlap' in QMin:
        errorcodes = run_wfoverlap(QMin, errorcodes)

    # do all necessary Theodore calculations
    errorcodes = run_theodore(QMin, errorcodes)

    # read all the output files
    QMout = getQMout(QMin)


    if PRINT or DEBUG:
        printQMout(QMin, QMout)

    # backup data if requested
    if 'backup' in QMin:
        backupdata(QMin['backup'], QMin)

    # Measure time
    runtime = measuretime()
    QMout['runtime'] = runtime

    # Write QMout
    writeQMout(QMin, QMout, QMinfilename)

    # Remove Scratchfiles from SCRATCHDIR
    if not DEBUG:
        #cleandir(QMin['scratchdir'])
        if 'cleanup' in QMin:
            cleandir(QMin['savedir'])

    print
    print(datetime.datetime.now())
    print('#================ END ================#')


if __name__ == '__main__':
    main()






# kate: indent-width 4
