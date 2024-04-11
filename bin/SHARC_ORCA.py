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


# Modules:
# Operating system, isfile and related routines, move files, create directories
import os
import shutil
# External Calls to MOLCAS
import subprocess as sp
# Command line arguments
import sys
# Regular expressions
import re
# debug print for dicts and arrays
import pprint
# sqrt and other math
import math
import cmath
# runtime measurement
import datetime
# copy of arrays of arrays
from copy import deepcopy
# parallel calculations
from multiprocessing import Pool
import time
# hostname
from socket import gethostname
# write debug traces when in pool threads
import traceback
# parse Python literals from input
import ast
import struct

# ======================================================================= #

version = '3.0'
versiondate = datetime.date(2023, 4, 1)


changelogstring = '''
16.05.2018: INITIAL VERSION
- functionality as SHARC_GAUSSIAN.py, minus restricted triplets
- QM/MM capabilities in combination with TINKER
- AO overlaps computed by PyQuante (only up to f functions)

11.09.2018:
- added "basis_per_element", "basis_per_atom", and "hfexchange" keywords

03.10.2018:
Update for Orca 4.1:
- SOC for restricted singlets and triplets
- gradients for restricted triplets
- multigrad features
- orca_fragovl instead of PyQuante

16.10.2018:
Update for Orca 4.1, after revisions:
- does not work with Orca 4.0 or lower (orca_fragovl unavailable, engrad/pcgrad files)

11.10.2020:
- COBRAMM can be used for QM/MM calculations
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
    9: '9-et',
    10: '10-et',
    11: '11-et',
    12: '12-et',
    13: '13-et',
    14: '14-et',
    15: '15-et',
    16: '16-et',
    17: '17-et',
    18: '18-et',
    19: '19-et',
    20: '20-et',
    21: '21-et',
    22: '22-et',
    23: '23-et',
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
           'Li': 0, 'Be': 0, 'B': 1, 'C': 1, 'N': 1, 'O': 1, 'F': 1, 'Ne': 1,
           'Na': 1, 'Mg': 1, 'Al': 5, 'Si': 5, 'P': 5, 'S': 5, 'Cl': 5, 'Ar': 5,
           'K': 5, 'Ca': 5,
           'Sc': 5, 'Ti': 5, 'V': 5, 'Cr': 5, 'Mn': 5, 'Fe': 5, 'Co': 5, 'Ni': 5, 'Cu': 5, 'Zn': 5,
           'Ga': 9, 'Ge': 9, 'As': 9, 'Se': 9, 'Br': 9, 'Kr': 9,
           'Rb': 9, 'Sr': 9,
           'Y': 14, 'Zr': 14, 'Nb': 14, 'Mo': 14, 'Tc': 14, 'Ru': 14, 'Rh': 14, 'Pd': 14, 'Ag': 14, 'Cd': 14,
           'In': 18, 'Sn': 18, 'Sb': 18, 'Te': 18, 'I': 18, 'Xe': 18,
           'Cs': 18, 'Ba': 18,
           'La': 18,
           'Ce': 18, 'Pr': 18, 'Nd': 18, 'Pm': 18, 'Sm': 18, 'Eu': 18, 'Gd': 18, 'Tb': 18, 'Dy': 18, 'Ho': 18, 'Er': 18, 'Tm': 18, 'Yb': 18, 'Lu': 23,
           'Hf': 23, 'Ta': 23, 'W': 23, 'Re': 23, 'Os': 23, 'Ir': 23, 'Pt': 23, 'Au': 23, 'Hg': 23,
           'Tl': 34, 'Pb': 34, 'Bi': 34, 'Po': 34, 'At': 34, 'Rn': 34,
           'Fr': 34, 'Ra': 34,
           'Ac': 34,
           'Th': 34, 'Pa': 34, 'U': 34, 'Np': 34, 'Pu': 34, 'Am': 34, 'Cm': 34, 'Bk': 34, 'Cf': 34, 'Es': 34, 'Fm': 34, 'Md': 34, 'No': 34, 'Lr': 34,
           'Rf': 50, 'Db': 50, 'Sg': 50, 'Bh': 50, 'Hs': 50, 'Mt': 50, 'Ds': 50, 'Rg': 50, 'Cn': 50,
           'Nh': 50, 'Fl': 50, 'Mc': 50, 'Lv': 50, 'Ts': 50, 'Og': 50
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
au2eV = 27.2113987622
kcal_to_Eh = 0.0015936011

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


def isbinary(path):
    return (re.search(r':.* text', sp.Popen(["file", '-L', path], stdout=sp.PIPE).stdout.read()) is None)

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
# =========================================== print routines ==================================== #
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
    string += '||' + ' ' * 28 + 'SHARC - ORCA - Interface' + ' ' * 28 + '||\n'
    string += '||' + ' ' * 80 + '||\n'
    string += '||' + ' ' * 14 + 'Authors: Sebastian Mai, Lea Ibele, and Moritz Heindl' + ' ' * 14 + '||\n'
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

    string = 'Restricted:    '
    for i in itmult(QMin['states']):
        string += '%5s       ' % (QMin['jobs'][QMin['multmap'][i]]['restr'])
    print(string)

    string = 'Method: \t'
    if QMin['template']['no_tda']:
        string += 'TD-'
    else:
        string += 'TDA-'
    string += QMin['template']['functional'].split()[0].upper()
    string += '/%s' % (QMin['template']['basis'])
    parts = []
    if QMin['template']['dispersion']:
        parts.append(QMin['template']['dispersion'].split()[0].upper())
    if QMin['template']['qmmm']:
        parts.append('QM/MM')
    if len(parts) > 0:
        string += '\t('
        string += ','.join(parts)
        string += ')'
    print(string)

    string = 'Found Geo'
    if 'veloc' in QMin:
        string += ' and Veloc! '
    else:
        string += '! '
    string += 'NAtom is %i.\n' % (QMin['natom_orig'])
    print(string)

    string = 'Geometry in Bohrs (%i atoms):\n' % QMin['natom_orig']
    if DEBUG:
        for i in range(QMin['natom_orig']):
            string += '%2s ' % (QMin['geo_orig'][i][0])
            for j in range(3):
                string += '% 7.4f ' % (QMin['geo_orig'][i][j + 1])
            string += '\n'
    else:
        for i in range(min(QMin['natom_orig'], 5)):
            string += '%2s ' % (QMin['geo_orig'][i][0])
            for j in range(3):
                string += '% 7.4f ' % (QMin['geo_orig'][i][j + 1])
            string += '\n'
        if QMin['natom_orig'] > 5:
            string += '..     ...     ...     ...\n'
            string += '%2s ' % (QMin['geo_orig'][-1][0])
            for j in range(3):
                string += '% 7.4f ' % (QMin['geo_orig'][-1][j + 1])
            string += '\n'
    print(string)

    if 'veloc' in QMin and DEBUG:
        string = ''
        for i in range(QMin['natom_orig']):
            string += '%s ' % (QMin['geo_orig'][i][0])
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

    print('State map:')
    pprint.pprint(QMin['statemap'])
    print

    for i in sorted(QMin):
        if not any([i == j for j in ['h', 'dm', 'soc', 'dmdr', 'socdr', 'theodore', 'geo', 'veloc', 'states', 'comment', 'grad', 'nacdr', 'ion', 'overlap', 'template', 'statemap', 'pointcharges', 'geo_orig', 'qmmm']]):
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
    if 'socdr' in QMin:
        string += writeQMoutsocdr(QMin, QMout)
    if 'dmdr' in QMin:
        string += writeQMoutdmdr(QMin, QMout)
    if 'ion' in QMin:
        string += writeQMoutprop(QMin, QMout)
    if 'theodore' in QMin or QMin['template']['qmmm']:
        string += writeQMoutTHEODORE(QMin, QMout)
    if 'phases' in QMin:
        string += writeQmoutPhases(QMin, QMout)
    if 'grad' in QMin:
        if QMin['template']['cobramm']:
            writeQMoutgradcobramm(QMin, QMout)
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

# =================================== #


def writeQMoutgradcobramm(QMin, QMout):
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
    natom = len(QMout['pcgrad'][0])
    string = ''
    # string+='! %i Gradient Vectors (%ix%ix3, real)\n' % (3,nmstates,natom)
    i = 0
    for imult, istate, ims in itnmstates(states):
        string += '%i %i ! %i %i %i\n' % (natom, 3, imult, istate, ims)
        for atom in range(natom):
            for xyz in range(3):
                print((eformat(QMout['pcgrad'][i][atom][xyz], 9, 3)), i, atom)
                string += '%s ' % (eformat(QMout['pcgrad'][i][atom][xyz], 9, 3))
            string += '\n'
        # string+='\n'
        i += 1
    string += '\n'
    writefile("grad_charges", string)
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

    # print property matrix (flag 11) for backwards compatibility
    string = ''
    string += '! %i Property Matrix (%ix%i, complex)\n' % (11, nmstates, nmstates)
    string += '%i %i\n' % (nmstates, nmstates)
    for i in range(nmstates):
        for j in range(nmstates):
            string += '%s %s ' % (eformat(QMout['prop'][i][j].real, 12, 3), eformat(QMout['prop'][i][j].imag, 12, 3))
        string += '\n'
    string += '\n'

    # print property matrices (flag 20) in new format
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
    if QMin['template']['qmmm']:
        nprop += len(QMin['qmmm']['MMEnergy_terms'])
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
    if QMin['template']['qmmm']:
        for label in sorted(QMin['qmmm']['MMEnergy_terms']):
            descriptors.append(label)
            string += label + '\n'

    string += '! Property Vectors (%ix%i, real)\n' % (nprop, nmstates)
    if 'theodore' in QMin:
        for i in range(QMin['template']['theodore_n']):
            string += '! TheoDORE descriptor %i (%s)\n' % (i + 1, descriptors[i])
            for j in range(nmstates):
                string += '%s\n' % (eformat(QMout['theodore'][j][i].real, 12, 3))
    if QMin['template']['qmmm']:
        for label in sorted(QMin['qmmm']['MMEnergy_terms']):
            string += '! QM/MM energy contribution (%s)\n' % (label)
            for j in range(nmstates):
                string += '%s\n' % (eformat(QMin['qmmm']['MMEnergy_terms'][label], 12, 3))
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
# =========================================== QM/MM ============================================= #
# =============================================================================================== #
# =============================================================================================== #

def prepare_QMMM(QMin, table_file):
    ''' creates dictionary with:
    MM coordinates (including connectivity and atom types)
    QM coordinates (including Link atom stuff)
    point charge data (including redistribution for Link atom neighbors)
    reorder arrays (for internal processing, all QM, then all LI, then all MM)

    is only allowed to read the following keys from QMin:
    geo
    natom
    QM/MM related infos from template
    '''

    table = readfile(table_file)


    # read table file
    print('===== Running QM/MM preparation ====')
    print('Reading table file ...         ', datetime.datetime.now())
    QMMM = {}
    QMMM['qmmmtype'] = []
    QMMM['atomtype'] = []
    QMMM['connect'] = []
    allowed = ['qm', 'mm']
    # read table file
    for iline, line in enumerate(table):
        s = line.split()
        if len(s) == 0:
            continue
        if not s[0].lower() in allowed:
            print('Not allowed QMMM-type "%s" on line %i!' % (s[0], iline + 1))
            sys.exit(15)
        QMMM['qmmmtype'].append(s[0].lower())
        QMMM['atomtype'].append(s[1])
        QMMM['connect'].append(set())
        for i in s[2:]:
            QMMM['connect'][-1].add(int(i) - 1)           # internally, atom numbering starts at 0
    QMMM['natom_table'] = len(QMMM['qmmmtype'])


    # list of QM and MM atoms
    QMMM['QM_atoms'] = []
    QMMM['MM_atoms'] = []
    for iatom in range(QMMM['natom_table']):
        if QMMM['qmmmtype'][iatom] == 'qm':
            QMMM['QM_atoms'].append(iatom)
        elif QMMM['qmmmtype'][iatom] == 'mm':
            QMMM['MM_atoms'].append(iatom)

    # make connections redundant and fill bond array
    print('Checking connection table ...  ', datetime.datetime.now())
    QMMM['bonds'] = set()
    for iatom in range(QMMM['natom_table']):
        for jatom in QMMM['connect'][iatom]:
            QMMM['bonds'].add(tuple(sorted([iatom, jatom])))
            QMMM['connect'][jatom].add(iatom)
    QMMM['bonds'] = sorted(list(QMMM['bonds']))


    # find link bonds
    print('Finding link bonds ...         ', datetime.datetime.now())
    QMMM['linkbonds'] = []
    QMMM['LI_atoms'] = []
    for i, j in QMMM['bonds']:
        if QMMM['qmmmtype'][i] != QMMM['qmmmtype'][j]:
            link = {}
            if QMMM['qmmmtype'][i] == 'qm':
                link['qm'] = i
                link['mm'] = j
            elif QMMM['qmmmtype'][i] == 'mm':
                link['qm'] = j
                link['mm'] = i
            link['scaling'] = {'qm': 0.3, 'mm': 0.7}
            link['element'] = 'H'
            link['atom'] = [link['element'], 0., 0., 0.]
            for xyz in range(3):
                link['atom'][xyz + 1] += link['scaling']['mm'] * QMin['geo'][link['mm']][xyz + 1]
                link['atom'][xyz + 1] += link['scaling']['qm'] * QMin['geo'][link['qm']][xyz + 1]
            QMMM['linkbonds'].append(link)
            QMMM['LI_atoms'].append(QMMM['natom_table'] - 1 + len(QMMM['linkbonds']))
            QMMM['atomtype'].append('999')
            QMMM['connect'].append(set([link['qm'], link['mm']]))


    # check link bonds
    mm_in_links = []
    qm_in_links = []
    mm_in_link_neighbors = []
    for link in QMMM['linkbonds']:
        mm_in_links.append(link['mm'])
        qm_in_links.append(link['qm'])
        for j in QMMM['connect'][link['mm']]:
            if QMMM['qmmmtype'][j] == 'mm':
                mm_in_link_neighbors.append(j)
    mm_in_link_neighbors.extend(mm_in_links)
    # no QM atom is allowed to be bonded to two MM atoms
    if not len(qm_in_links) == len(set(qm_in_links)):
        print('Some QM atom is involved in more than one link bond!')
        sys.exit(16)
    # no MM atom is allowed to be bonded to two QM atoms
    if not len(mm_in_links) == len(set(mm_in_links)):
        print('Some MM atom is involved in more than one link bond!')
        sys.exit(17)
    # no neighboring MM atoms are allowed to be involved in link bonds
    if not len(mm_in_link_neighbors) == len(set(mm_in_link_neighbors)):
        print('An MM-link atom is bonded to another MM-link atom!')
        sys.exit(18)


    # check geometry and connection table
    if not QMMM['natom_table'] == QMin['natom']:
        print('Number of atoms in table file does not match number of atoms in QMin!')
        sys.exit(19)


    # process MM geometry (and convert to angstrom!)
    QMMM['MM_coords'] = []
    for atom in QMin['geo']:
        QMMM['MM_coords'].append([atom[0]] + [i * au2a for i in atom[1:4]])
    for ilink, link in enumerate(QMMM['linkbonds']):
        QMMM['MM_coords'].append(['HLA'] + link['atom'][1:4])


    # create reordering dicts
    print('Creating reorder mappings ...  ', datetime.datetime.now())
    QMMM['reorder_input_MM'] = {}
    QMMM['reorder_MM_input'] = {}
    j = -1
    for i, t in enumerate(QMMM['qmmmtype']):
        if t == 'qm':
            j += 1
            QMMM['reorder_MM_input'][j] = i
    for ilink, link in enumerate(QMMM['linkbonds']):
        j += 1
        QMMM['reorder_MM_input'][j] = QMMM['natom_table'] + ilink
    for i, t in enumerate(QMMM['qmmmtype']):
        if t == 'mm':
            j += 1
            QMMM['reorder_MM_input'][j] = i
    for i in QMMM['reorder_MM_input']:
        QMMM['reorder_input_MM'][QMMM['reorder_MM_input'][i]] = i


    # process QM geometry (including link atoms), QM coords in bohr!
    QMMM['QM_coords'] = []
    QMMM['reorder_input_QM'] = {}
    QMMM['reorder_QM_input'] = {}
    j = -1
    for iatom in range(QMMM['natom_table']):
        if QMMM['qmmmtype'][iatom] == 'qm':
            QMMM['QM_coords'].append(deepcopy(QMin['geo'][iatom]))
            j += 1
            QMMM['reorder_input_QM'][iatom] = j
            QMMM['reorder_QM_input'][j] = iatom
    for ilink, link in enumerate(QMMM['linkbonds']):
        QMMM['QM_coords'].append(link['atom'])
        j += 1
        QMMM['reorder_input_QM'][-(ilink + 1)] = j
        QMMM['reorder_QM_input'][j] = -(ilink + 1)


    # process charge redistribution around link bonds
    # point charges are in input geometry ordering
    print('Charge redistribution ...      ', datetime.datetime.now())
    QMMM['charge_distr'] = []
    for iatom in range(QMMM['natom_table']):
        if QMMM['qmmmtype'][iatom] == 'qm':
            QMMM['charge_distr'].append([(0., 0)])
        elif QMMM['qmmmtype'][iatom] == 'mm':
            if iatom in mm_in_links:
                QMMM['charge_distr'].append([(0., 0)])
            else:
                QMMM['charge_distr'].append([(1., iatom)])
    for link in QMMM['linkbonds']:
        mm_neighbors = []
        for j in QMMM['connect'][link['mm']]:
            if QMMM['qmmmtype'][j] == 'mm':
                mm_neighbors.append(j)
        if len(mm_neighbors) > 0:
            factor = 1. / len(mm_neighbors)
            for j in QMMM['connect'][link['mm']]:
                if QMMM['qmmmtype'][j] == 'mm':
                    QMMM['charge_distr'][j].append((factor, link['mm']))

    # pprint.pprint(QMMM)
    return QMMM

# ======================================================================= #


def execute_tinker(QMin, ff_file_path):
    '''
    run tinker to get:
    * MM energy
    * MM gradient
    * point charges

    is only allowed to read the following keys from QMin:
    qmmm
    scratchdir
    savedir
    tinker
    '''

    QMMM = QMin['qmmm']

    # prepare Workdir
    WORKDIR = os.path.join(QMin['scratchdir'], 'TINKER')
    mkdir(WORKDIR)


    # key file
    # string='parameters %s\nQMMM %i\nQM %s\n' % (
    # ff_file_path,
    # QMMM['natom_table']+len(QMMM['linkbonds']),
    # ' '.join( [ str(QMMM['reorder_input_MM'][i]+1) for i in QMMM['QM_atoms'] ] )
    # )
    # if len(QMMM['linkbonds'])>0:
    # string+='LA %s\n' % (
    # ' '.join( [ str(QMMM['reorder_input_MM'][i]+1) for i in QMMM['LI_atoms'] ] ) )
    # string+='MM %s\n' % (
    # ' '.join( [ str(QMMM['reorder_input_MM'][i]+1) for i in QMMM['MM_atoms'] ] )  )
    # string+='\nDEBUG\n'
    # if len(QMMM['linkbonds'])>0:
    # string+='atom    999    99    HLA     "Hydrogen Link Atom"        1      1.008     0\n'
    # string+='\n'
    # filename=os.path.join(WORKDIR,'TINKER.key')
    # writefile(filename,string)


    print('Writing TINKER inputs ...      ', datetime.datetime.now())
    # key file
    string = 'parameters %s\nQMMM %i\n' % (ff_file_path, QMMM['natom_table'] + len(QMMM['linkbonds']))
    string += 'QM %i %i\n' % (-1, len(QMMM['QM_atoms']))
    if len(QMMM['linkbonds']) > 0:
        string += 'LA %s\n' % (
            ' '.join([str(QMMM['reorder_input_MM'][i] + 1) for i in QMMM['LI_atoms']]))
    string += 'MM %i %i\n' % (-(1 + len(QMMM['QM_atoms']) + len(QMMM['linkbonds'])),
                              QMMM['natom_table'] + len(QMMM['linkbonds']))
    # if DEBUG:
    # string+='\nDEBUG\n'
    if QMin['ncpu'] > 1:
        string += '\nOPENMP-THREADS %i\n' % QMin['ncpu']
    if len(QMMM['linkbonds']) > 0:
        string += 'atom    999    99    HLA     "Hydrogen Link Atom"        1      1.008     0\n'
    # string+='CUTOFF 1.0\n'
    string += '\n'
    filename = os.path.join(WORKDIR, 'TINKER.key')
    writefile(filename, string)


    # xyz/type/connection file
    string = '%i\n' % (len(QMMM['MM_coords']))
    for iatom_MM in range(len(QMMM['MM_coords'])):
        iatom_input = QMMM['reorder_MM_input'][iatom_MM]
        string += '% 5i  %3s  % 16.12f % 16.12f % 16.12f  %4s  %s\n' % (
            iatom_MM + 1,
            QMMM['MM_coords'][iatom_input][0],
            QMMM['MM_coords'][iatom_input][1],
            QMMM['MM_coords'][iatom_input][2],
            QMMM['MM_coords'][iatom_input][3],
            QMMM['atomtype'][iatom_input],
            ' '.join([str(QMMM['reorder_input_MM'][i] + 1) for i in sorted(QMMM['connect'][iatom_input])])
        )
    filename = os.path.join(WORKDIR, 'TINKER.xyz')
    writefile(filename, string)


    # communication file
    string = 'SHARC 0 -1\n'
    for iatom_MM in range(len(QMMM['MM_coords'])):
        iatom_input = QMMM['reorder_MM_input'][iatom_MM]
        string += '% 16.12f % 16.12f % 16.12f\n' % tuple(QMMM['MM_coords'][iatom_input][1:4])
    filename = os.path.join(WORKDIR, 'TINKER.qmmm')
    writefile(filename, string)


    # standard input file
    string = 'TINKER.xyz'
    filename = os.path.join(WORKDIR, 'TINKER.in')
    writefile(filename, string)


    # run TINKER
    runTINKER(WORKDIR, QMin['tinker'], QMin['savedir'], strip=False, ncpu=QMin['ncpu'])


    # read out TINKER
    filename = os.path.join(WORKDIR, 'TINKER.qmmm')
    output = readfile(filename)

    # check success
    if 'MMisOK' not in output[0]:
        print('TINKER run not successful!')
        sys.exit(20)

    # get MM energy (convert from kcal to Hartree)
    print('Searching MMEnergy ...         ', datetime.datetime.now())
    QMMM['MMEnergy'] = float(output[1].split()[-1]) * kcal_to_Eh

    # get MM gradient (convert from kcal/mole/A to Eh/bohr)
    print('Searching MMGradient ...       ', datetime.datetime.now())
    QMMM['MMGradient'] = {}
    for line in output:
        if 'MMGradient' in line:
            s = line.split()
            iatom_MM = int(s[1]) - 1
            iatom_input = QMMM['reorder_MM_input'][iatom_MM]
            grad = [float(i) * kcal_to_Eh * au2a for i in s[2:5]]
            QMMM['MMGradient'][iatom_input] = grad
        if 'MMq' in line:
            break

    # get MM point charges
    print('Searching MMpc_raw ...         ', datetime.datetime.now())
    QMMM['MMpc_raw'] = {}
    for i in range(QMMM['natom_table']):
        QMMM['MMpc_raw'][i] = 0.
    iline = 0
    while True:
        iline += 1
        line = output[iline]
        if 'MMq' in line:
            break
    iatom_MM = len(QMMM['QM_atoms']) + len(QMMM['LI_atoms']) - 1
    while True:
        iline += 1
        iatom_MM += 1
        line = output[iline]
        if 'NMM' in line:
            break
        s = line.split()
        q = float(s[-1])
        QMMM['MMpc_raw'][QMMM['reorder_MM_input'][iatom_MM]] = q

    # compute actual charges (including redistribution)
    print('Redistributing charges ...     ', datetime.datetime.now())
    QMMM['MMpc'] = {}
    for i in range(QMMM['natom_table']):
        s = 0.
        for factor, iatom in QMMM['charge_distr'][i]:
            s += factor * QMMM['MMpc_raw'][iatom]
        QMMM['MMpc'][i] = s

    # make list of pointcharges without QM atoms
    print('Finalizing charges ...         ', datetime.datetime.now())
    QMMM['pointcharges'] = []
    QMMM['reorder_pc_input'] = {}
    ipc = 0
    for iatom_input in QMMM['MM_atoms']:
        atom = QMMM['MM_coords'][iatom_input]
        q = QMMM['MMpc'][iatom_input]
        QMMM['pointcharges'].append(atom[1:4] + [q])
        QMMM['reorder_pc_input'][ipc] = iatom
        ipc += 1






    # Get energy components from standard out (debug print)
    filename = os.path.join(WORKDIR, 'TINKER.out')
    output = readfile(filename)
    QMMM['MMEnergy_terms'] = {}
    for line in output:
        if 'kcal/mol' in line:
            s = line.split()
            # QMMM['MMEnergy_terms'][s[0]]=float(s[-2])*kcal_to_Eh
            QMMM['MMEnergy_terms'][s[0]] = float(s[2])

    print('====================================')
    print('\n')

    # DONE! Final results:
    # QMMM['MMEnergy']
    # QMMM['MMGradient']
    # QMMM['MMpc']
    # QMMM['QM_coords']
    # QMMM['reorder_input_QM']
    # QMMM['reorder_QM_input']

    # print('='*60)
    # print('E:',QMMM['MMEnergy'])
    # print('Grad:')
    # pprint.pprint( QMMM['MMGradient'] )
    # print('MM coord:')
    # print str(QMMM['natom_table']) +'\n'
    # for atom in QMMM['MM_coords']:
    # print atom[0], atom[1], atom[2], atom[3]
    # print('QM coord:')
    # print str(len(QMMM['QM_coords'])) +'\n'
    # for atom in QMMM['QM_coords']:
    # print atom[0], atom[1]*au2a, atom[2]*au2a, atom[3]*au2a
    # print('MM pc:')
    # for iatom,atom in enumerate(QMMM['MM_coords']):
    # q=QMMM['MMpc'][iatom]
    # if q!=0.:
    # print atom[1],atom[2],atom[3],q
    # pprint.pprint(QMMM['MMEnergy_terms'])
    # print('='*60)

    # pprint.pprint(QMMM)
    # sys.exit(21)
    return QMMM

# ======================================================================= #


def coords_same(coord1, coord2):
    thres = 1e-5
    s = 0.
    for i in range(3):
        s += (coord1[i] - coord2[i])**2
    s = math.sqrt(s)
    return s <= thres

# ======================================================================= #


def runTINKER(WORKDIR, tinker, savedir, strip=False, ncpu=1):
    prevdir = os.getcwd()
    os.chdir(WORKDIR)
    string = os.path.join(tinker, 'bin', 'tkr2qm_s') + ' '
    string += ' < TINKER.in'
    os.environ['OMP_NUM_THREADS'] = str(ncpu)
    if PRINT or DEBUG:
        starttime = datetime.datetime.now()
        sys.stdout.write('START:\t%s\t%s\t"%s"\n' % (shorten_DIR(WORKDIR), starttime, shorten_DIR(string)))
        sys.stdout.flush()
    stdoutfile = open(os.path.join(WORKDIR, 'TINKER.out'), 'w')
    stderrfile = open(os.path.join(WORKDIR, 'TINKER.err'), 'w')
    try:
        runerror = sp.call(string, shell=True, stdout=stdoutfile, stderr=stderrfile)
    except OSError:
        print('Call have had some serious problems:', OSError)
        sys.exit(22)
    stdoutfile.close()
    stderrfile.close()
    if PRINT or DEBUG:
        endtime = datetime.datetime.now()
        sys.stdout.write('FINISH:\t%s\t%s\tRuntime: %s\tError Code: %i\n' % (shorten_DIR(WORKDIR), endtime, endtime - starttime, runerror))
        sys.stdout.flush()
    if DEBUG and runerror != 0:
        copydir = os.path.join(savedir, 'debug_TINKER_stdout')
        if not os.path.isdir(copydir):
            mkdir(copydir)
        outfile = os.path.join(WORKDIR, 'TINKER.out')
        tofile = os.path.join(copydir, "TINKER_problems.out")
        shutil.copy(outfile, tofile)
        print('Error in %s! Copied TINKER output to %s' % (WORKDIR, tofile))
    os.chdir(prevdir)
    if strip and not DEBUG and runerror == 0:
        stripWORKDIR(WORKDIR)
    return runerror

# ======================================================================= #


def transform_QM_QMMM(QMin, QMout):

    # Meta data
    QMin['natom'] = QMin['natom_orig']
    QMin['geo'] = QMin['geo_orig']

    # Hamiltonian
    if 'h' in QMout:
        for i in range(QMin['nmstates']):
            QMout['h'][i][i] += QMin['qmmm']['MMEnergy']

    # Gradients
    if 'grad' in QMout:
        nmstates = QMin['nmstates']
        natom = QMin['natom_orig']
        grad = [[[0. for i in range(3)] for j in range(natom)] for k in range(nmstates)]
        # QM gradient
        for iqm in QMin['qmmm']['reorder_QM_input']:
            iqmmm = QMin['qmmm']['reorder_QM_input'][iqm]
            if iqmmm < 0:
                ilink = -iqmmm - 1
                link = QMin['qmmm']['linkbonds'][ilink]
                for istate in range(nmstates):
                    for ixyz in range(3):
                        grad[istate][link['qm']][ixyz] += QMout['grad'][istate][iqm][ixyz] * link['scaling']['qm']
                        grad[istate][link['mm']][ixyz] += QMout['grad'][istate][iqm][ixyz] * link['scaling']['mm']
            else:
                for istate in range(nmstates):
                    for ixyz in range(3):
                        grad[istate][iqmmm][ixyz] += QMout['grad'][istate][iqm][ixyz]
        # PC gradient
        for iqm, iqmmm in enumerate(QMin['qmmm']['MM_atoms']):
            for istate in range(nmstates):
                for ixyz in range(3):
                    grad[istate][iqmmm][ixyz] += QMout['pcgrad'][istate][iqm][ixyz]
        # MM gradient
        for iqmmm in range(QMin['qmmm']['natom_table']):
            for istate in range(nmstates):
                for ixyz in range(3):
                    grad[istate][iqmmm][ixyz] += QMin['qmmm']['MMGradient'][iqmmm][ixyz]
        QMout['grad'] = grad

    # pprint.pprint(QMout)
    return QMin, QMout


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
            sys.exit(23)
    else:
        try:
            os.makedirs(SCRATCHDIR)
        except OSError:
            print('Can not create SCRATCHDIR=%s\n' % (SCRATCHDIR))
            sys.exit(24)

# ======================================================================= #


def removequotes(string):
    if string.startswith("'") and string.endswith("'"):
        return string[1:-1]
    elif string.startswith('"') and string.endswith('"'):
        return string[1:-1]
    else:
        return string

# ======================================================================= #


def getOrcaVersion(path):
    # run orca with nonexisting file
    string = os.path.join(path, 'orca') + ' nonexisting'
    try:
        proc = sp.Popen(string, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    except OSError:
        print('Call have had some serious problems:', OSError)
        sys.exit(25)
    comm = proc.communicate()[0].decode()
    # find version string
    for line in comm.split('\n'):
        if 'Program Version' in line:
            s = line.split('-')[0].split()[2].split('.')
            s = tuple([int(i) for i in s])
            return s
    print('Could not find Orca version!')
    sys.exit(26)

# ======================================================================= #


def getsh2Orcakey(sh2Orca, key):
    i = -1
    while True:
        i += 1
        try:
            line = re.sub('#.*$', '', sh2Orca[i])
        except IndexError:
            break
        line = line.split(None, 1)
        if line == []:
            continue
        if key.lower() in line[0].lower():
            return line
    return ['', '']

# ======================================================================= #


def get_sh2Orca_environ(sh2Orca, key, environ=True, crucial=True):
    line = getsh2Orcakey(sh2Orca, key)
    if line[0]:
        LINE = line[1]
        LINE = removequotes(LINE).strip()
    else:
        if environ:
            LINE = os.getenv(key.upper())
            if not LINE:
                if crucial:
                    print('Either set $%s or give path to %s in ORCA.resources!' % (key.upper(), key.upper()))
                    sys.exit(27)
                else:
                    return None
        else:
            if crucial:
                print('Give path to %s in ORCA.resources!' % (key.upper()))
                sys.exit(28)
            else:
                return None
    LINE = os.path.expandvars(LINE)
    LINE = os.path.expanduser(LINE)
    if containsstring(';', LINE):
        print("$%s contains a semicolon. Do you probably want to execute another command after %s? I can't do that for you..." % (key.upper(), key.upper()))
        sys.exit(29)
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
            sys.exit(30)
        if 'end' in line:
            break
        fields = line.split()
        try:
            nacpairs.append([int(fields[0]), int(fields[1])])
        except ValueError:
            print('"nacdr select" is followed by pairs of state indices, each pair on a new line!')
            sys.exit(31)
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
        sys.exit(32)
    QMin['natom'] = natom
    if len(QMinlines) < natom + 4:
        print('Input file must contain at least:\nnatom\ncomment\ngeometry\nkeyword "states"\nat least one task')
        sys.exit(33)

    # Save Comment line
    QMin['comment'] = QMinlines[1]

    # Get geometry and possibly velocity (for backup-analytical non-adiabatic couplings)
    QMin['geo'] = []
    QMin['veloc'] = []
    hasveloc = True
    QMin['frozcore'] = 0
    QMin['Atomcharge'] = 0
    for i in range(2, natom + 2):
        # only check line formatting for first 1000 atoms
        if i < 1000 and not containsstring('[a-zA-Z][a-zA-Z]?[0-9]*.*[-]?[0-9]+[.][0-9]*.*[-]?[0-9]+[.][0-9]*.*[-]?[0-9]+[.][0-9]*', QMinlines[i]):
            print('Input file does not comply to xyz file format! Maybe natom is just wrong.')
            sys.exit(34)
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
            sys.exit(35)
    else:
        factor = 1. / au2a

    for iatom in range(len(QMin['geo'])):
        for ixyz in range(3):
            QMin['geo'][iatom][ixyz + 1] *= factor


    if 'states' not in QMin:
        print('Keyword "states" not given!')
        sys.exit(36)
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
    if 'states' not in QMin:
        print('Number of states not given in QM input file %s!' % (QMinfilename))
        sys.exit(37)

    # only singlets TODO
    # if len(QMin['states'])>1 and sum(QMin['states'][1:])>0:
        # print('Currently, only singlet states are allowed!')
        # sys.exit(38)

    possibletasks = ['h', 'soc', 'dm', 'grad', 'overlap', 'dmdr', 'socdr', 'ion', 'theodore', 'phases']
    if not any([i in QMin for i in possibletasks]):
        print('No tasks found! Tasks are %s.' % possibletasks)
        sys.exit(39)

    if 'h' not in QMin and 'soc' not in QMin:
        QMin['h'] = []

    # if  'soc' in QMin:
        # print('Spin-orbit couplings cannot be computed currently!')
        # sys.exit(40)

    if 'soc' in QMin and (len(QMin['states']) < 3 or QMin['states'][2] <= 0):
        QMin = removekey(QMin, 'soc')
        QMin['h'] = []
        print('HINT: No triplet states requested, turning off SOC request.')

    if 'samestep' in QMin and 'init' in QMin:
        print('"Init" and "Samestep" cannot be both present in QM.in!')
        sys.exit(41)

    if 'restart' in QMin and 'init' in QMin:
        print('"Init" and "Samestep" cannot be both present in QM.in!')
        sys.exit(42)

    if 'phases' in QMin:
        QMin['overlap'] = []

    if 'overlap' in QMin and 'init' in QMin:
        print('"overlap" and "phases" cannot be calculated in the first timestep! Delete either "overlap" or "init"')
        sys.exit(43)

    if 'init' not in QMin and 'samestep' not in QMin and 'restart' not in QMin:
        QMin['newstep'] = []

    if not any([i in QMin for i in ['h', 'soc', 'dm', 'grad']]) and 'overlap' in QMin:
        QMin['h'] = []

    if 'nacdt' in QMin or 'nacdr' in QMin:
        print('Within the SHARC-ORCA interface couplings can only be calculated via the overlap method. "nacdr" and "nacdt" are not supported.')
        sys.exit(44)

    if 'dmdr' in QMin:
        print('Dipole derivatives ("dmdr") not currently supported')
        sys.exit(45)

    if 'socdr' in QMin:
        print('Spin-orbit coupling derivatives ("socdr") are not implemented')
        sys.exit(46)


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
                    sys.exit(47)
                if QMin['grad'][i] > nmstates:
                    print('State for requested gradient does not correspond to any state in QM input file state list!')
                    sys.exit(48)


    # shorten states array TODO
    # QMin['states']=[QMin['states'][0]]



# --------------------------------------------- ORCA.resources ----------------------------------

    QMin['pwd'] = os.getcwd()

    # open ORCA.resources
    filename = 'ORCA.resources'
    if os.path.isfile(filename):
        sh2Orca = readfile(filename)
    else:
        print('HINT: reading resources from SH2Orc.inp')
        sh2Orca = readfile('SH2Orc.inp')


    # Set up scratchdir
    line = get_sh2Orca_environ(sh2Orca, 'scratchdir', False, False)
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
        line = get_sh2Orca_environ(sh2Orca, 'savedir', False, False)
        if line is None:
            line = QMin['pwd'] + '/SAVEDIR/'
    line = os.path.expandvars(line)
    line = os.path.expanduser(line)
    line = os.path.abspath(line)
    if 'init' in QMin:
        checkscratch(line)
    QMin['savedir'] = line
    link(QMin['savedir'], os.path.join(QMin['pwd'], 'SAVE'), False, False)


    # setup environment for Orca
    QMin['orcadir'] = get_sh2Orca_environ(sh2Orca, 'orcadir')
    if 'LD_LIBRARY_PATH' in os.environ:
        os.environ['LD_LIBRARY_PATH'] = '%s:' % (QMin['orcadir']) + os.environ['LD_LIBRARY_PATH']
    else:
        os.environ['LD_LIBRARY_PATH'] = '%s' % (QMin['orcadir'])
    QMin['OrcaVersion'] = getOrcaVersion(QMin['orcadir'], )
    print('Detected ORCA version %s' % (str(QMin['OrcaVersion'])))
    os.environ['PATH'] = '%s:' % (QMin['orcadir']) + os.environ['PATH']
    if QMin['OrcaVersion'] < (4, 1):
        print('This version of the SHARC-ORCA interface is only compatible to Orca 4.1 or higher!')
        sys.exit(49)


    # debug option
    line = getsh2Orcakey(sh2Orca, 'debug')
    if line[0]:
        if len(line) <= 1 or 'true' in line[1].lower():
            global DEBUG
            DEBUG = True


    # save_stuff option
    QMin['save_stuff'] = False
    line = getsh2Orcakey(sh2Orca, 'save_stuff')
    if line[0]:
        if len(line) <= 1 or 'true' in line[1].lower():
            QMin['save_stuff'] = True



    # debug option
    line = getsh2Orcakey(sh2Orca, 'no_print')
    if line[0]:
        if len(line) <= 1 or 'true' in line[1].lower():
            global PRINT
            PRINT = False



    # resources
    QMin['ncpu'] = 1
    line = getsh2Orcakey(sh2Orca, 'ncpu')
    if line[0]:
        try:
            QMin['ncpu'] = int(line[1])
        except ValueError:
            print('Number of CPUs does not evaluate to numerical value!')
            sys.exit(50)
    if os.environ.get('NSLOTS') is not None:
        QMin['ncpu'] = int(os.environ.get('NSLOTS'))
        print('Detected $NSLOTS variable. Will use ncpu=%i' % (QMin['ncpu']))
    elif os.environ.get('SLURM_NTASKS_PER_NODE') is not None:
        QMin['ncpu'] = int(os.environ.get('SLURM_NTASKS_PER_NODE'))
        print('Detected $SLURM_NTASKS_PER_NODE variable. Will use ncpu=%i' % (QMin['ncpu']))
    QMin['ncpu'] = max(1, QMin['ncpu'])

    QMin['delay'] = 0.0
    line = getsh2Orcakey(sh2Orca, 'delay')
    if line[0]:
        try:
            QMin['delay'] = float(line[1])
        except ValueError:
            print('Submit delay does not evaluate to numerical value!')
            sys.exit(51)

    QMin['schedule_scaling'] = 0.9
    line = getsh2Orcakey(sh2Orca, 'schedule_scaling')
    if line[0]:
        try:
            x = float(line[1])
            if 0 < x <= 1.:
                QMin['schedule_scaling'] = x
        except ValueError:
            print('"schedule_scaling" does not evaluate to numerical value!')
            sys.exit(52)


    # initial MO guess settings
    # if neither keyword is present, the interface will reuse MOs from savedir, or let ADF generate a guess
    line = getsh2Orcakey(sh2Orca, 'always_orb_init')
    if line[0]:
        QMin['always_orb_init'] = []
    line = getsh2Orcakey(sh2Orca, 'always_guess')
    if line[0]:
        QMin['always_guess'] = []
    if 'always_orb_init' in QMin and 'always_guess' in QMin:
        print('Keywords "always_orb_init" and "always_guess" cannot be used together!')
        sys.exit(53)


    # wfoverlap settings
    if 'overlap' in QMin or 'ion' in QMin:
        # WFoverlap
        QMin['wfoverlap'] = get_sh2Orca_environ(sh2Orca, 'wfoverlap', False, False)
        if QMin['wfoverlap'] is None:
            ciopath = os.path.join(os.path.expandvars(os.path.expanduser('$SHARC')), 'wfoverlap.x')
            if os.path.isfile(ciopath):
                QMin['wfoverlap'] = ciopath
            else:
                print('Give path to wfoverlap.x in ORCA.resources!')
                sys.exit(54)
        # PyQuante
        # QMin['pyquante']=get_sh2Orca_environ(sh2Orca,'pyquante',False,False)
        # if QMin['pyquante']==None or not os.path.isdir(QMin['pyquante']):
            # print('Give path to the PyQuante installation directory in ORCA.resources!')
            # sys.exit(55)
        # if 'PYTHONPATH' in os.environ:
            # os.environ['PYTHONPATH']+=os.pathsep + QMin['pyquante']
        # else:
            # os.environ['PYTHONPATH']=QMin['pyquante']
        # sys.path.append(QMin['pyquante'])

    # memory
    QMin['memory'] = 100
    line = getsh2Orcakey(sh2Orca, 'memory')
    if line[0]:
        QMin['memory'] = float(line[1])

    # truncation threshold
    QMin['wfthres'] = 0.99
    line = getsh2Orcakey(sh2Orca, 'wfthres')
    if line[0]:
        QMin['wfthres'] = float(line[1])

    # get the nooverlap keyword: no dets will be extracted if present
    line = getsh2Orcakey(sh2Orca, 'nooverlap')
    if line[0]:
        QMin['nooverlap'] = []


    # TheoDORE settings
    if 'theodore' in QMin:
        QMin['theodir'] = get_sh2Orca_environ(sh2Orca, 'theodir', False, False)
        if QMin['theodir'] is None or not os.path.isdir(QMin['theodir']):
            print('Give path to the TheoDORE installation directory in ORCA.resources!')
            sys.exit(56)
        os.environ['THEODIR'] = QMin['theodir']
        os.environ['THEODIR'] = QMin['theodir']
        if 'PYTHONPATH' in os.environ:
            os.environ['PYTHONPATH'] = os.path.join(QMin['theodir'], 'lib') + os.pathsep + QMin['theodir'] + os.pathsep + os.environ['PYTHONPATH']
            # print os.environ['PYTHONPATH']
        else:
            os.environ['PYTHONPATH'] = os.path.join(QMin['theodir'], 'lib') + os.pathsep + QMin['theodir']


    # neglected gradients
    QMin['neglected_gradient'] = 'zero'
    if 'grad' in QMin:
        line = getsh2Orcakey(sh2Orca, 'neglected_gradient')
        if line[0]:
            if line[1].lower().strip() == 'zero':
                QMin['neglected_gradient'] = 'zero'
            elif line[1].lower().strip() == 'gs':
                QMin['neglected_gradient'] = 'gs'
            elif line[1].lower().strip() == 'closest':
                QMin['neglected_gradient'] = 'closest'
            else:
                print('Unknown argument to "neglected_gradient"!')
                sys.exit(57)



# --------------------------------------------- ORCA.template ----------------------------------

    # define classes and defaults
    bools = {'no_tda': False,
             'unrestricted_triplets': False,
             'qmmm': False,
             'cobramm': False,
             'picture_change': False
             }
    strings = {'basis': '6-31G',
               'auxbasis': '',
               'functional': 'PBE',
               'dispersion': '',
               'grid': '',
               'gridx': '',
               'gridxc': '',
               'ri': '',
               'scf': '',
               'qmmm_table': 'ORCA.qmmm.table',
               'qmmm_ff_file': 'ORCA.ff',
               'keys': '',
               'paste_input_file': ''
               }
    integers = {
        'frozen': -1,
        'maxiter': 700
    }
    floats = {
        'hfexchange': -1.,
        'intacc': -1.
    }
    special = {'paddingstates': [0 for i in QMin['states']],
               'charge': [i % 2 for i in range(len(QMin['states']))],
               'theodore_prop': ['Om', 'PRNTO', 'S_HE', 'Z_HE', 'RMSeh'],
               'theodore_fragment': [],
               'basis_per_element': {},
               'ecp_per_element': {},
               'basis_per_atom': {},
               'range_sep_settings': {'do': False, 'mu': 0.14, 'scal': 1.0, 'ACM1': 0.0, 'ACM2': 0.0, 'ACM3': 1.0}
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
    template = readfile('ORCA.template')

    # go through template
    for line in template:
        orig = re.sub('#.*$', '', line).strip()
        line = orig.lower().split()
        if len(line) == 0:
            continue
        elif line[0] in bools:
            QMin['template'][line[0]] = True
        elif line[0] in strings:
            QMin['template'][line[0]] = orig.split(None, 1)[1]
        elif line[0] in integers:
            QMin['template'][line[0]] = int(float(line[1]))
        elif line[0] in floats:
            QMin['template'][line[0]] = float(line[1])
        elif line[0] in special:

            # paddingstates needs to be autoexpanded and checked
            if line[0] == 'paddingstates':
                if len(line) == 2:
                    QMin['template']['paddingstates'] = [int(line[1]) for i in range(len(QMin['states']))]
                elif len(line) - 1 >= len(QMin['states']):
                    QMin['template']['paddingstates'] = [int(line[1 + i]) for i in range(len(QMin['states']))]
                else:
                    print('Length of "paddingstates" does not match length of "states"!')
                    sys.exit(58)
                for i in range(len(QMin['template']['paddingstates'])):
                    if QMin['template']['paddingstates'][i] < 0:
                        QMin['template']['paddingstates'][i] = 0

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
                        # sys.exit(59)
                        # print('Charges from template not compatible with multiplicities!')
                        # sys.exit(60)
                else:
                    print('Length of "charge" does not match length of "states"!')
                    sys.exit(61)

            # basis_per_element can occur several times
            elif line[0] == 'basis_per_element':
                line2 = orig.split(None, 2)
                QMin['template']['basis_per_element'][line2[1]] = line2[2]

            # ecp_per_element can occur several times
            elif line[0] == 'ecp_per_element':
                line2 = orig.split(None, 2)
                QMin['template']['ecp_per_element'][line2[1]] = line2[2]

            # basis_per_atom can occur several times
            elif line[0] == 'basis_per_atom':
                line2 = orig.split(None, 2)
                QMin['template']['basis_per_atom'][int(line2[1]) - 1] = line2[2]

            # list of floats must be parsed separately
            elif line[0] == 'range_sep_settings':
                QMin['template']['range_sep_settings']['do'] = True
                QMin['template']['range_sep_settings']['mu'] = float(line[1])
                QMin['template']['range_sep_settings']['scal'] = float(line[2])
                QMin['template']['range_sep_settings']['ACM1'] = float(line[3])
                QMin['template']['range_sep_settings']['ACM2'] = float(line[4])
                QMin['template']['range_sep_settings']['ACM3'] = float(line[5])


    # go through sh2Orca for the theodore settings and QM/MM file names
    for line in sh2Orca:
        orig = re.sub('#.*$', '', line).strip()
        line = orig.lower().split()
        if len(line) == 0:
            continue
        elif line[0] in special:

            # TheoDORE properties need to be parsed in a special way
            if line[0] == 'theodore_prop':
                if '[' in orig:
                    string = orig.split(None, 1)[1]
                    QMin['template']['theodore_prop'] = ast.literal_eval(string)
                else:
                    QMin['template']['theodore_prop'] = []
                    s = orig.split(None)[1:]
                    for i in s:
                        QMin['template']['theodore_prop'].append(i)
                theodore_spelling = ['Om',
                                     'PRNTO',
                                     'Z_HE', 'S_HE', 'RMSeh',
                                     'POSi', 'POSf', 'POS',
                                     'PRi', 'PRf', 'PR', 'PRh',
                                     'CT', 'CT2', 'CTnt',
                                     'MC', 'LC', 'MLCT', 'LMCT', 'LLCT',
                                     'DEL', 'COH', 'COHh']
                for i in range(len(QMin['template']['theodore_prop'])):
                    for j in theodore_spelling:
                        if QMin['template']['theodore_prop'][i].lower() == j.lower():
                            QMin['template']['theodore_prop'][i] = j

            # TheoDORE fragments need to be parsed in a special way
            elif line[0] == 'theodore_fragment':
                if '[' in orig:
                    string = orig.split(None, 1)[1]
                    QMin['template']['theodore_fragment'] = ast.literal_eval(string)
                else:
                    s = orig.split(None)[1:]
                    l = []
                    for i in s:
                        l.append(int(i))
                    QMin['template']['theodore_fragment'].append(l)



    if QMin['template']['paste_input_file']:
        path = os.path.expandvars(os.path.expanduser(QMin['template']['paste_input_file']))
        if os.path.isfile(path):
            QMin['template']['paste_input_file'] = readfile(path)
        else:
            print('Additional input file %s not found!' % path)
            sys.exit(62)


    # do logic checks
    if not QMin['template']['unrestricted_triplets']:
        if QMin['OrcaVersion'] < (4, 1):
            if len(QMin['states']) >= 3 and QMin['states'][2] > 0:
                print('With Orca v<4.1, triplets can only be computed with the unrestricted_triplets option!')
                sys.exit(62)
        if len(QMin['template']['charge']) >= 3 and QMin['template']['charge'][0] != QMin['template']['charge'][2]:
            print('Charges of singlets and triplets differ. Please enable the "unrestricted_triplets" option!')
            sys.exit(63)
    if QMin['template']['unrestricted_triplets'] and 'soc' in QMin:
        if len(QMin['states']) >= 3 and QMin['states'][2] > 0:
            print('Request "SOC" is not compatible with "unrestricted_triplets"!')
            sys.exit(64)
    # if QMin['template']['cosmo'] and 'grad' in QMin: TODO
        # print('COSMO is not compatible with gradient calculations!')
        # sys.exit(65)
    if QMin['template']['ecp_per_element'] and 'soc' in QMin:
        if len(QMin['states']) >= 3 and QMin['states'][2] > 0:
            print('Request "SOC" is not compatible with using ECPs!')
            sys.exit(64)




# --------------------------------------------- QM/MM ----------------------------------

    # qmmm keyword
    QMin['qmmm'] = QMin['template']['qmmm']

    # prepare everything
    QMin['geo_orig'] = QMin['geo']
    QMin['natom_orig'] = QMin['natom']
    QMin['frozcore_orig'] = QMin['frozcore']
    QMin['Atomcharge_orig'] = QMin['Atomcharge']
    if QMin['qmmm']:
        # charge HINT
        print('HINT: For QM/MM calculations, you have to specify in the template file the charge for the *QM region only*!')
        print('The automatic assignment of total charge might not work if the MM part is not neutral!\n')

        # get settings from ORCA.resources
        # Tinker
        line = getsh2Orcakey(sh2Orca, 'tinker')
        if not line[0]:
            print('TINKER path not given!')
            sys.exit(66)
        line = os.path.expandvars(line[1].strip())
        line = os.path.expanduser(line)
        line = os.path.abspath(line)
        QMin['tinker'] = line
        if not os.path.isfile(os.path.join(QMin['tinker'], 'bin', 'tkr2qm_s')):
            print('TINKER executable at "%s" not found!' % os.path.join(QMin['tinker'], 'bin', 'tkr2qm_s'))
            sys.exit(67)

        # table and ff files
        for line in sh2Orca:
            orig = re.sub('#.*$', '', line).strip()
            line = orig.lower().split()
            if len(line) == 0:
                continue
            elif line[0] == 'qmmm_table':
                line2 = orig.split(None, 1)
                if len(line2) < 2:
                    print('Please specify a connection table file after "qmmm_table"!')
                    sys.exit(68)
                filename = os.path.abspath(os.path.expandvars(os.path.expanduser(line2[1])))
                QMin['template']['qmmm_table'] = filename
            elif line[0] == 'qmmm_ff_file':
                line2 = orig.split(None, 1)
                if len(line2) < 2:
                    print('Please specify a force field file after "qmmm_ff_file"!')
                    sys.exit(69)
                filename = os.path.abspath(os.path.expandvars(os.path.expanduser(line2[1])))
                QMin['template']['qmmm_ff_file'] = filename

        # prepare data structures and run Tinker
        QMin['qmmm'] = prepare_QMMM(QMin, QMin['template']['qmmm_table'])
        execute_tinker(QMin, QMin['template']['qmmm_ff_file'])   # modifies QMin['qmmm'] in place !

        # modify QMin dict
        # QMin['geo_orig']=QMin['geo']
        QMin['geo'] = QMin['qmmm']['QM_coords']
        # QMin['natom_orig']=QMin['natom']
        QMin['natom'] = len(QMin['geo'])
        # QMin['frozcore_orig']=QMin['frozcore']
        QMin['frozcore'] = 0
        # QMin['Atomcharge_orig']=QMin['Atomcharge']
        QMin['Atomcharge'] = 0
        for atom in QMin['geo']:
            symb = atom[0]
            QMin['frozcore'] += FROZENS[symb]
            QMin['Atomcharge'] += ATOMCHARGE[symb]
        QMin['pointcharges'] = QMin['qmmm']['pointcharges']






    # number of frozen core orbitals for wfoverlap
    # this code is down here because we need to check for QM/MM in template first
    line = getsh2Orcakey(sh2Orca, 'numfrozcore')
    if line[0]:
        numfroz = int(line[1])
        if numfroz == 0:
            QMin['frozcore'] = 0
        elif numfroz > 0:
            QMin['frozcore'] = numfroz
        elif numfroz < 0:
            pass        # here we take frozcore from above



    # number of doubly occupied orbitals for Dyson
    line = getsh2Orcakey(sh2Orca, 'numocc')
    if line[0]:
        numfroz = int(line[1])
        if numfroz <= 0:
            QMin['ndocc'] = 0
        elif numfroz > 0:
            QMin['ndocc'] = max(0, numfroz - QMin['frozcore'])
    else:
        QMin['ndocc'] = 0

# --------------------------------------------- Logic ----------------------------------

    # obtain the statemap
    statemap = {}
    i = 1
    for imult, istate, ims in itnmstates(QMin['states']):
        statemap[i] = [imult, istate, ims]
        i += 1
    QMin['statemap'] = statemap

    # obtain the states to actually compute
    states_to_do = deepcopy(QMin['states'])
    for i in range(len(QMin['states'])):
        if states_to_do[i] > 0:
            states_to_do[i] += QMin['template']['paddingstates'][i]
    if not QMin['template']['unrestricted_triplets']:
        if len(QMin['states']) >= 3 and QMin['states'][2] > 0:
            states_to_do[0] = max(QMin['states'][0], 1)
            req = max(QMin['states'][0] - 1, QMin['states'][2])
            states_to_do[0] = req + 1
            states_to_do[2] = req
    QMin['states_to_do'] = states_to_do

    # make the jobs
    jobs = {}
    if QMin['states_to_do'][0] > 0:
        jobs[1] = {'mults': [1], 'restr': True}
    if len(QMin['states_to_do']) >= 2 and QMin['states_to_do'][1] > 0:
        jobs[2] = {'mults': [2], 'restr': False}
    if len(QMin['states_to_do']) >= 3 and QMin['states_to_do'][2] > 0:
        if not QMin['template']['unrestricted_triplets'] and QMin['states_to_do'][0] > 0:
            if QMin['OrcaVersion'] >= (4, 1):
                jobs[1]['mults'].append(3)
            else:
                jobs[3] = {'mults': [1, 3], 'restr': True}
        else:
            jobs[3] = {'mults': [3], 'restr': False}
    if len(QMin['states_to_do']) >= 4:
        for imult, nstate in enumerate(QMin['states_to_do'][3:]):
            if nstate > 0:
                # jobs[len(jobs)+1]={'mults':[imult+4],'restr':False}
                jobs[imult + 4] = {'mults': [imult + 4], 'restr': False}
    QMin['jobs'] = jobs

    # make the multmap (mapping between multiplicity and job)
    # multmap[imult]=ijob
    # multmap[-ijob]=[imults]
    multmap = {}
    for ijob in jobs:
        job = jobs[ijob]
        for imult in job['mults']:
            multmap[imult] = ijob
        multmap[-(ijob)] = job['mults']
    multmap[1] = 1
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

    # get the set of states for which gradients actually need to be calculated
    gradmap = set()
    if 'grad' in QMin:
        for i in QMin['grad']:
            gradmap.add(tuple(statemap[i][0:2]))
            # gradmap.add(      (statemap[i][0],1) )
    gradmap = list(gradmap)
    gradmap.sort()
    QMin['gradmap'] = gradmap

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
                # print m1,job1,el1,m2,job2,el2
                if abs(m1 - m2) == 1 and abs(el1 - el2) == 1:
                    ionmap.append((m1, job1, m2, job2))
        QMin['ionmap'] = ionmap

    # number of properties/entries calculated by TheoDORE
    if 'theodore' in QMin:
        QMin['template']['theodore_n'] = len(QMin['template']['theodore_prop']) + len(QMin['template']['theodore_fragment'])**2
    else:
        QMin['template']['theodore_n'] = 0

# --------------------------------------------- File setup ----------------------------------

    # check for initial orbitals
    initorbs = {}
    if 'always_guess' in QMin:
        QMin['initorbs'] = {}
    elif 'init' in QMin or 'always_orb_init' in QMin:
        for job in QMin['joblist']:
            filename = os.path.join(QMin['pwd'], 'ORCA.gbw.init')
            if os.path.isfile(filename):
                initorbs[job] = filename
        for job in QMin['joblist']:
            filename = os.path.join(QMin['pwd'], 'ORCA.gbw.%i.init' % (job))
            if os.path.isfile(filename):
                initorbs[job] = filename
        if 'always_orb_init' in QMin and len(initorbs) < njobs:
            print('Initial orbitals missing for some jobs!')
            sys.exit(70)
        QMin['initorbs'] = initorbs
    elif 'newstep' in QMin:
        for job in QMin['joblist']:
            filename = os.path.join(QMin['savedir'], 'ORCA.gbw.%i' % (job))
            if os.path.isfile(filename):
                initorbs[job] = filename + '.old'     # file will be moved to .old
            else:
                print('File %s missing in savedir!' % (filename))
                sys.exit(71)
        QMin['initorbs'] = initorbs
    elif 'samestep' in QMin:
        for job in QMin['joblist']:
            filename = os.path.join(QMin['savedir'], 'ORCA.gbw.%i' % (job))
            if os.path.isfile(filename):
                initorbs[job] = filename
            else:
                print('File %s missing in savedir!' % (filename))
                sys.exit(72)
        QMin['initorbs'] = initorbs
    elif 'restart' in QMin:
        for job in QMin['joblist']:
            filename = os.path.join(QMin['savedir'], 'ORCA.gbw.%i.old' % (job))
            if os.path.isfile(filename):
                initorbs[job] = filename
            else:
                print('File %s missing in savedir!' % (filename))
                sys.exit(73)
        QMin['initorbs'] = initorbs


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
        print('======= DEBUG print for QMin =======')
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
        optimal[i] = nrounds // parallel_speedup(ncores, scaling)
    # print optimal
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
    # print nrounds,nslots,cpu_per_run
    return nrounds, nslots, cpu_per_run

# =============================================================================================== #


def generate_joblist(QMin):

    # pprint.pprint(QMin)

    # sort the gradients into the different jobs
    gradjob = {}
    for ijob in QMin['joblist']:
        gradjob['master_%i' % ijob] = {}
    for grad in QMin['gradmap']:
        ijob = QMin['multmap'][grad[0]]
        isgs = False
        if not QMin['jobs'][ijob]['restr']:
            if grad[1] == 1:
                isgs = True
        else:
            if grad == (1, 1):
                isgs = True
        istates = QMin['states_to_do'][grad[0] - 1]
        if QMin['OrcaVersion'] < (4, 1):
            if isgs and istates > 1:
                gradjob['grad_%i_%i' % grad] = {}
                gradjob['grad_%i_%i' % grad][grad] = {'gs': True}
            else:
                n = 0
                for gradx in gradjob['master_%i' % ijob]:
                    n += 1
                if n > 0:
                    gradjob['grad_%i_%i' % grad] = {}
                    gradjob['grad_%i_%i' % grad][grad] = {'gs': False}
                else:
                    gradjob['master_%i' % ijob][grad] = {'gs': False}
        else:
            gradjob['master_%i' % ijob][grad] = {'gs': isgs}
    # pprint.pprint(gradjob)

    # make map for states onto gradjobs
    jobgrad = {}
    for job in gradjob:
        for state in gradjob[job]:
            jobgrad[state] = (job, gradjob[job][state]['gs'])
    QMin['jobgrad'] = jobgrad
    # print gradjob
    # print
    # print jobgrad
    # sys.exit(74)

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
            remove = ['gradmap', 'ncpu']
            for r in remove:
                QMin1 = removekey(QMin1, r)
            QMin1['gradmap'] = list(gradjob[i])
            QMin1['ncpu'] = cpu_per_run[icount]
            if QMin['OrcaVersion'] < (4, 1):
                if 3 in QMin['multmap'][-QMin1['IJOB']] and QMin['jobs'][QMin1['IJOB']]['restr']:
                    QMin1['states'][0] = 1
                    QMin1['states_to_do'][0] = 1
            if QMin1['qmmm']:
                QMin1['qmmm'] = True
            icount += 1
            schedule[-1][i] = QMin1

    # add the gradient calculations
    ntasks = 0
    for i in gradjob:
        if 'grad' in i:
            ntasks += 1
    if ntasks > 0:
        nrounds, nslots, cpu_per_run = divide_slots(QMin['ncpu'], ntasks, QMin['schedule_scaling'])
        QMin['nslots_pool'].append(nslots)
        schedule.append({})
        icount = 0
        for i in gradjob:
            if 'grad' in i:
                QMin1 = deepcopy(QMin)
                mult = list(gradjob[i])[0][0]
                QMin1['IJOB'] = QMin['multmap'][mult]
                remove = ['gradmap', 'ncpu', 'h', 'soc', 'dm', 'overlap', 'ion', 'always_guess', 'always_orb_init', 'init']
                for r in remove:
                    QMin1 = removekey(QMin1, r)
                QMin1['gradmap'] = list(gradjob[i])
                QMin1['ncpu'] = cpu_per_run[icount]
                QMin1['gradonly'] = []
                if QMin1['qmmm']:
                    QMin1['qmmm'] = True
                icount += 1
                schedule[-1][i] = QMin1

    return QMin, schedule


# =============================================================================================== #
# =============================================================================================== #
# ======================================= Orca Job Execution ==================================== #
# =============================================================================================== #
# =============================================================================================== #

def runjobs(schedule, QMin):

    if 'newstep' in QMin:
        moveOldFiles(QMin)

    print('>>>>>>>>>>>>> Starting the ORCA job execution')

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
    print(string)
    if any((i != 0 for i in errorcodes.values())):
        print('Some subprocesses did not finish successfully!')
        print('See %s:%s for error messages in ORCA output.' % (gethostname(), QMin['scratchdir']))
        sys.exit(75)
    print

    if PRINT:
        print('>>>>>>>>>>>>> Saving files')
        starttime = datetime.datetime.now()
    for ijobset, jobset in enumerate(schedule):
        if not jobset:
            continue
        for ijob, job in enumerate(jobset):
            if 'master' in job:
                WORKDIR = os.path.join(QMin['scratchdir'], job)
                # if not 'samestep' in QMin or 'molden' in QMin:
                # if 'molden' in QMin or not 'nooverlap' in QMin:
                # saveMolden(WORKDIR,jobset[job])
                if 'samestep' not in QMin:
                    saveFiles(WORKDIR, jobset[job])
                if 'ion' in QMin and ijobset == 0 and ijob == 0:
                    saveAOmatrix(WORKDIR, QMin)
    # saveGeometry(QMin)
    if PRINT:
        endtime = datetime.datetime.now()
        print('Saving Runtime: %s' % (endtime - starttime))
    print

    return errorcodes

# ======================================================================= #


def run_calc(WORKDIR, QMin):
    try:
        setupWORKDIR(WORKDIR, QMin)
        strip = True
        err = runORCA(WORKDIR, QMin['orcadir'], strip)
        # err=0
    except Exception as problem:
        print('*' * 50 + '\nException in run_calc(%s)!' % (WORKDIR))
        traceback.print_exc()
        print('*' * 50 + '\n')
        raise problem

    return err

# ======================================================================= #


def setupWORKDIR(WORKDIR, QMin):
    # mkdir the WORKDIR, or clean it if it exists, then copy all necessary files from pwd and savedir
    # then put the ORCA.inp file

    # setup the directory
    mkdir(WORKDIR)

    # write ORCA.inp
    inputstring = writeORCAinput(QMin)
    filename = os.path.join(WORKDIR, 'ORCA.inp')
    writefile(filename, inputstring)
    if DEBUG:
        print('================== DEBUG input file for WORKDIR %s =================' % (shorten_DIR(WORKDIR)))
        print(inputstring)
        print('ORCA input written to: %s' % (filename))
        print('====================================================================')
    # write point charges
    if QMin['qmmm']:
        inputstring = write_pccoord_file(QMin)
        filename = os.path.join(WORKDIR, 'ORCA.pc')
        writefile(filename, inputstring)
        if DEBUG:
            print('================== DEBUG input file for WORKDIR %s =================' % (shorten_DIR(WORKDIR)))
            print(inputstring)
            print('Point charges written to: %s' % (filename))
            print('====================================================================')
    # if QMin['template']['cobramm']:
    #    inputstring=write_pc_cobramm(QMin)
    #    filename=os.path.join(WORKDIR,'charge.dat')
    #    writefile(filename,inputstring)
    if QMin['template']['cobramm']:
        currentDirectory = os.getcwd()
        fromfile = os.path.join(currentDirectory, 'charge.dat')
        tofile = tofile = os.path.join(WORKDIR, 'charge.dat')
        shutil.copy(fromfile, tofile)

    # wf file copying
    if 'master' in QMin:
        job = QMin['IJOB']
        if job in QMin['initorbs']:
            fromfile = QMin['initorbs'][job]
            tofile = os.path.join(WORKDIR, 'ORCA.gbw')
            shutil.copy(fromfile, tofile)
    elif 'grad' in QMin:
        job = QMin['IJOB']
        fromfile = os.path.join(QMin['scratchdir'], 'master_%i' % job, 'ORCA.gbw')
        tofile = os.path.join(WORKDIR, 'ORCA.gbw')
        shutil.copy(fromfile, tofile)

    # force field file copying
    # if QMin['template']['qmmm']:
        # fromfile=QMin['template']['qmmm_ff_file']
        # tofile=os.path.join(WORKDIR,'ADF.ff')
        # shutil.copy(fromfile,tofile)

    return

# ======================================================================= #


def writeORCAinput(QMin):
    # split gradmap into smaller chunks
    Nmax_gradlist = 255
    gradmaps = [sorted(QMin['gradmap'])[i:i + Nmax_gradlist] for i in range(0, len(QMin['gradmap']), Nmax_gradlist)]

    # make multi-job input
    string = ''
    for ichunk, chunk in enumerate(gradmaps):
        if ichunk >= 1:
            string += '\n\n$new_job\n\n%base "ORCA"\n\n'
        QMin_copy = deepcopy(QMin)
        QMin_copy['gradmap'] = chunk
        string += ORCAinput_string(QMin_copy)
    if not gradmaps:
        string += ORCAinput_string(QMin)
    return string


# ======================================================================= #
def ORCAinput_string(QMin):
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
    states_to_do[gsmult - 1] -= 1

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
        # sing=states_to_do[0]>0
        trip = (len(states_to_do) >= 3 and states_to_do[2] > 0)
    else:
        ncalc = max(states_to_do)
        # mults_td=''

    # whether to do SOC
    # sopert=False
    # gscorr=False
    # if 'soc' in QMin:
        # if restr:
        # nsing=QMin['states'][0]
        # if len(QMin['states'])>=3:
        # ntrip=QMin['states'][2]
        # else:
        # ntrip=0
        # if nsing+ntrip>=2 and ntrip>=1:
        # sopert=True
        # if nsing>=1:
        # gscorr=True

    # gradients
    multigrad = False
    if 'grad' in QMin and QMin['gradmap']:
        dograd = True
        egrad = ()
        for grad in QMin['gradmap']:
            if not (gsmult, 1) == grad:
                egrad = grad
        # if len(QMin['gradmap'])>1:
        if QMin['OrcaVersion'] >= (4, 1):
            multigrad = True
            singgrad = []
            tripgrad = []
            for grad in QMin['gradmap']:
                if grad[0] == gsmult:
                    singgrad.append(grad[1] - 1)
                if grad[0] == 3 and restr:
                    tripgrad.append(grad[1])
    else:
        dograd = False

    # construct the input string
    string = ''

    # main line
    string += '! '

    keys = ['basis',
            'auxbasis',
            'functional',
            'dispersion',
            'ri',
            'keys']
    for i in keys:
        string += '%s ' % (QMin['template'][i])
    keys = ['nousesym']
    for i in keys:
        string += '%s ' % (i)

    if QMin['template']['grid']:
      string += 'grid%s ' % QMin['template']['grid']
    if QMin['template']['gridx']:
        string += 'gridx%s ' % QMin['template']['gridx']
# In this way, one can change grid on individual atoms:
# %method
# SpecialGridAtoms 26,15,-1,-4         # for element 26 and, for atom index 1 and 4 (cannot change on atom 0!)
# SpecialGridIntAcc 7,6,5,5            # size of grid
# end

    if dograd:
        string += 'engrad'

    string += '\n'

    # cpu cores
    if QMin['ncpu'] > 1 and 'AOoverlap' not in QMin:
        string += '%%pal\n  nprocs %i\nend\n\n' % (QMin['ncpu'])
    string += '%%maxcore %i\n\n' % (QMin['memory'])

    # basis sets
    if QMin['template']['basis_per_element']:
        string += '%basis\n'
        for i in QMin['template']['basis_per_element']:
            string += 'newgto %s "%s" end\n' % (i, QMin['template']['basis_per_element'][i])
        if not QMin['template']['ecp_per_element']:
            string += 'end\n\n'

    # ECP basis sets
    if QMin['template']['ecp_per_element']:
        if QMin['template']['basis_per_element']:
            for i in QMin['template']['ecp_per_element']:
                string += 'newECP %s "%s" end\n' % (i, QMin['template']['ecp_per_element'][i])
            string += 'end\n\n'
        else:
            print("ECP defined without additional basis. Not implemented.")

    # frozen core
    if QMin['frozcore'] > 0:
        string += '%%method\nfrozencore -%i\nend\n\n' % (2 * QMin['frozcore'])
    else:
        string += '%method\nfrozencore FC_NONE\nend\n\n'

    # hf exchange
    if QMin['template']['hfexchange'] >= 0.:
        # string+='%%method\nScalHFX = %f\nScalDFX = %f\nend\n\n' % (QMin['template']['hfexchange'],1.-QMin['template']['hfexchange'])
        string += '%%method\nScalHFX = %f\nend\n\n' % (QMin['template']['hfexchange'])

    # Range separation
    if QMin['template']['range_sep_settings']['do']:
        string += '''%%method
 RangeSepEXX True
 RangeSepMu %f
 RangeSepScal %f
 ACM %f, %f, %f\nend\n\n
''' % (QMin['template']['range_sep_settings']['mu'],
            QMin['template']['range_sep_settings']['scal'],
            QMin['template']['range_sep_settings']['ACM1'],
            QMin['template']['range_sep_settings']['ACM2'],
            QMin['template']['range_sep_settings']['ACM3']
       )

    # Intacc
    if QMin['template']['intacc'] > 0.:
        string += '''%%method
  intacc %3.1f\nend\n\n''' % (QMin['template']['intacc'])


    # Gaussian point charge scheme
    if 'cpcm' in QMin['template']['keys'].lower():
        string += '''%cpcm
  surfacetype vdw_gaussian\nend\n\n'''



    # excited states
    if ncalc > 0 and 'AOoverlap' not in QMin:
        string += '%tddft\n'
        if not QMin['template']['no_tda']:
            string += 'tda true\n'
        else:
            string += 'tda false\n'
        if QMin['template']['gridxc']:
            string += 'gridxc %s\n' % (QMin['template']['gridxc'])
        if 'theodore' in QMin:
            string += 'tprint 0.0001\n'
        if restr and trip:
            string += 'triplets true\n'
        string += 'nroots %i\n' % (ncalc)
        if restr and 'soc' in QMin:
            string += 'dosoc true\n'
            string += 'printlevel 3\n'
        # string+="dotrans all\n" #TODO
        if dograd:
            if multigrad:
                if singgrad:
                    string += 'sgradlist '
                    string += ','.join([str(i) for i in sorted(singgrad)])
                    string += '\n'
                if tripgrad:
                    string += 'tgradlist '
                    string += ','.join([str(i) for i in sorted(tripgrad)])
                    string += '\n'
            elif egrad:
                string += 'iroot %i\n' % (egrad[1] - (gsmult == egrad[0]))
        string += 'end\n\n'

    # output
    string += '%output\n'
    if 'AOoverlap' in QMin or 'ion' in QMin or 'theodore' in QMin:
        string += 'Print[ P_Overlap ] 1\n'
    if 'master' in QMin or 'theodore' in QMin:
        string += 'Print[ P_MOs ] 1\n'
    string += 'end\n\n'

    # scf
    string += '%scf\n'
    if 'AOoverlap' in QMin:
        string += 'maxiter 0\n'
    else:
        string += 'maxiter %i\n' % (QMin['template']['maxiter'])
    string += 'end\n\n'

    # rel
    if QMin['template']['picture_change']:
        string += '%rel\nPictureChange true\nend\n\n'


    # TODO: workaround
    # if 'soc' in QMin and 'grad' in QMin:
        # string+='%rel\nonecenter true\nend\n\n'

    # charge mult geom
    string += '%coords\nCtyp xyz\nunits bohrs\n'
    if 'AOoverlap' in QMin:
        string += 'charge %i\n' % (2. * charge)
    else:
        string += 'charge %i\n' % (charge)
    string += 'mult %i\n' % (gsmult)
    string += 'coords\n'
    for iatom, atom in enumerate(QMin['geo']):
        label = atom[0]
        string += '%4s %16.9f %16.9f %16.9f' % (label, atom[1], atom[2], atom[3])
        if iatom in QMin['template']['basis_per_atom']:
            string += ' newgto "%s" end' % (QMin['template']['basis_per_atom'][iatom])
        string += '\n'
    string += 'end\nend\n\n'

    # point charges
    if QMin['qmmm']:
        string += '%pointcharges "ORCA.pc"\n\n'
    elif QMin['template']['cobramm']:
        string += '%pointcharges "charge.dat"\n\n'
    if QMin['template']['paste_input_file']:
        string += '\n'
        for line in QMin['template']['paste_input_file']:
            string += line
        string += '\n'
    return string

# ======================================================================= #


def write_pccoord_file(QMin):
    string = '%i\n' % len(QMin['pointcharges'])
    for atom in QMin['pointcharges']:
        string += '%f %f %f %f\n' % (atom[3], atom[0], atom[1], atom[2])
    return string

# def write_pc_cobramm(QMin):
#    cobcharges=open('charges.dat', 'r')
#    charges=cobcharges.readlines()
#    string='%i\n' % len(charges)
#    for atom in charges:
#        string+='%f %f %f %f\n' % (atom[3],atom[0],atom[1],atom[2])
#    return string

# ======================================================================= #


def shorten_DIR(string):
    maxlen = 50
    front = 12
    if len(string) > maxlen:
        return string[0:front] + '...' + string[-(maxlen - 3 - front):]
    else:
        return string + ' ' * (maxlen - len(string))

# ======================================================================= #


def runORCA(WORKDIR, orcadir, strip=False):
    prevdir = os.getcwd()
    os.chdir(WORKDIR)
    string = os.path.join(orcadir, 'orca') + ' '
    string += 'ORCA.inp'
    if PRINT or DEBUG:
        starttime = datetime.datetime.now()
        sys.stdout.write('START:\t%s\t%s\t"%s"\n' % (shorten_DIR(WORKDIR), starttime, shorten_DIR(string)))
        sys.stdout.flush()
    stdoutfile = open(os.path.join(WORKDIR, 'ORCA.log'), 'w')
    stderrfile = open(os.path.join(WORKDIR, 'ORCA.err'), 'w')
    try:
        runerror = sp.call(string, shell=True, stdout=stdoutfile, stderr=stderrfile)
    except OSError:
        print('Call have had some serious problems:', OSError)
        sys.exit(76)
    stdoutfile.close()
    stderrfile.close()
    if PRINT or DEBUG:
        endtime = datetime.datetime.now()
        sys.stdout.write('FINISH:\t%s\t%s\tRuntime: %s\tError Code: %i\n' % (shorten_DIR(WORKDIR), endtime, endtime - starttime, runerror))
        sys.stdout.flush()
    os.chdir(prevdir)
    if strip and not DEBUG and runerror == 0:
        stripWORKDIR(WORKDIR)
    return runerror

# ======================================================================= #


def stripWORKDIR(WORKDIR):
    ls = os.listdir(WORKDIR)
    keep = ['ORCA.inp$', 'ORCA.err$', 'ORCA.log$', 'ORCA.gbw', 'ORCA.cis', 'ORCA.engrad', 'ORCA.pcgrad', 'ORCA.molden.input', 'ORCA.pc']
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
    basenames = ['ORCA.gbw', 'ORCA.molden']
    if 'nooverlap' not in QMin:
        basenames.append('mos')
    for job in QMin['joblist']:
        for base in basenames:
            fromfile = os.path.join(QMin['savedir'], '%s.%i' % (base, job))
            if not os.path.isfile(fromfile):
                print('File %s not found, cannot move to OLD!' % (fromfile))
                sys.exit(77)
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
                sys.exit(78)
            tofile = os.path.join(QMin['savedir'], '%s.%i.old' % (base, job))
            if PRINT:
                print(shorten_DIR(fromfile) + '   =>   ' + shorten_DIR(tofile))
            shutil.copy(fromfile, tofile)
    # geometry file
    # fromfile=os.path.join(QMin['savedir'],'geom.dat')
    # tofile  =os.path.join(QMin['savedir'],'geom.dat.old')
    # if PRINT:
        # print(shorten_DIR(fromfile)+'   =>   '+shorten_DIR(tofile))
    # shutil.copy(fromfile,tofile)

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
# def saveGeometry(QMin):
    # string=''
    # for iatom,atom in enumerate(QMin['geo']):
    # label=atom[0]
    # string+='%4s %16.9f %16.9f %16.9f\n' % (label,atom[1],atom[2],atom[3])
    # filename=os.path.join(QMin['savedir'],'geom.dat')
    # writefile(filename,string)
    # if PRINT:
    # print shorten_DIR(filename)
    # return

# ======================================================================= #


def saveFiles(WORKDIR, QMin):

    # copy the gbw files from master directories
    job = QMin['IJOB']
    fromfile = os.path.join(WORKDIR, 'ORCA.gbw')
    tofile = os.path.join(QMin['savedir'], 'ORCA.gbw.%i' % (job))
    shutil.copy(fromfile, tofile)
    if PRINT:
        print(shorten_DIR(tofile))

    # make Molden files and copy to savedir
    saveMolden(WORKDIR, QMin)

    # if necessary, extract the MOs and write them to savedir
    if 'ion' in QMin or 'nooverlap' not in QMin:
        mofile = os.path.join(QMin['savedir'], 'mos.%i' % job)
        f = os.path.join(WORKDIR, 'ORCA.gbw')
        string = get_MO_from_gbw(f, QMin)
        # f=os.path.join(WORKDIR,'ORCA.log')
        # string=get_MO_from_stdout(f,QMin)
        # moldenfile=os.path.join(QMin['savedir'],'ORCA.molden.%i' % job)
        # string=make_mos_from_Molden(moldenfile,QMin)
        writefile(mofile, string)
        if PRINT:
            print(shorten_DIR(mofile))

    # if necessary, extract the TDDFT coefficients and write them to savedir
    if 'ion' in QMin or 'nooverlap' not in QMin:
        f = os.path.join(WORKDIR, 'ORCA.cis')
        strings = get_dets_from_cis(f, QMin)
        for f in strings:
            writefile(f, strings[f])
            if PRINT:
                print(shorten_DIR(f))

# ======================================================================= #


def saveMolden(WORKDIR, QMin):

    # run orca_2mkl
    prevdir = os.getcwd()
    os.chdir(WORKDIR)
    string = 'orca_2mkl ORCA -molden'
    stdoutfile = open(os.path.join(WORKDIR, 'orca_2mkl.out'), 'w')
    stderrfile = open(os.path.join(WORKDIR, 'orca_2mkl.err'), 'w')
    if PRINT or DEBUG:
        starttime = datetime.datetime.now()
        # sys.stdout.write('START:\t%s\t%s\t"%s"\n' % (shorten_DIR(WORKDIR),starttime,shorten_DIR(string)))
        sys.stdout.flush()
    try:
        runerror = sp.call(string, shell=True, stdout=stdoutfile, stderr=stderrfile)
    except OSError:
        print('Call have had some serious problems:', OSError)
        sys.exit(79)
    stdoutfile.close()
    stderrfile.close()
    if PRINT or DEBUG:
        endtime = datetime.datetime.now()
        # sys.stdout.write('FINISH:\t%s\t%s\tRuntime: %s\tError Code: %i\n' % (shorten_DIR(WORKDIR),endtime,endtime-starttime,runerror))
        sys.stdout.flush()
    os.chdir(prevdir)

    job = QMin['IJOB']
    fromfile = os.path.join(WORKDIR, 'ORCA.molden.input')
    tofile = os.path.join(QMin['savedir'], 'ORCA.molden.%i' % (job))
    shutil.copy(fromfile, tofile)
    if PRINT:
        print(shorten_DIR(tofile))


# ======================================================================= #
def get_MO_from_gbw(filename, QMin):

    # run orca_fragovl
    string = 'orca_fragovl %s %s' % (filename, filename)
    try:
        proc = sp.Popen(string, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    except OSError:
        print('Call have had some serious problems:', OSError)
        sys.exit(80)
    comm = proc.communicate()[0].decode()
    data = comm.split('\n')
    # get size of matrix
    for line in reversed(data):
        # print line
        s = line.split()
        if len(s) >= 1:
            NAO = int(line.split()[0]) + 1
            break

    job = QMin['IJOB']
    restr = QMin['jobs'][job]['restr']

    # find MO block
    iline = -1
    while True:
        iline += 1
        if len(data) <= iline:
            print('MOs not found!')
            sys.exit(81)
        line = data[iline]
        if 'FRAGMENT A MOs MATRIX' in line:
            break
    iline += 3

    # formatting
    nblock = 6
    npre = 11
    ndigits = 16
    # default_pos=[14,30,46,62,78,94]
    default_pos = [npre + 3 + ndigits * i for i in range(nblock)]  # does not include shift

    # get coefficients for alpha
    NMO_A = NAO
    MO_A = [[0. for i in range(NAO)] for j in range(NMO_A)]
    for imo in range(NMO_A):
        jblock = imo // nblock
        jcol = imo % nblock
        for iao in range(NAO):
            shift = max(0, len(str(iao)) - 3)
            jline = iline + jblock * (NAO + 1) + iao
            line = data[jline]
            # fix too long floats in strings
            dots = [idx for idx, item in enumerate(line.lower()) if '.' in item]
            diff = [dots[i] - default_pos[i] - shift for i in range(len(dots))]
            if jcol == 0:
                pre = 0
            else:
                pre = diff[jcol - 1]
            post = diff[jcol]
            # fixed
            val = float(line[npre + shift + jcol * ndigits + pre: npre + shift + ndigits + jcol * ndigits + post])
            MO_A[imo][iao] = val
    iline += ((NAO - 1) // nblock + 1) * (NAO + 1)

    # coefficients for beta
    if not restr:
        NMO_B = NAO
        MO_B = [[0. for i in range(NAO)] for j in range(NMO_B)]
        for imo in range(NMO_B):
            jblock = imo // nblock
            jcol = imo % nblock
            for iao in range(NAO):
                shift = max(0, len(str(iao)) - 3)
                jline = iline + jblock * (NAO + 1) + iao
                line = data[jline]
                # fix too long floats in strings
                dots = [idx for idx, item in enumerate(line.lower()) if '.' in item]
                diff = [dots[i] - default_pos[i] - shift for i in range(len(dots))]
                if jcol == 0:
                    pre = 0
                else:
                    pre = diff[jcol - 1]
                post = diff[jcol]
                # fixed
                val = float(line[npre + shift + jcol * ndigits + pre: npre + shift + ndigits + jcol * ndigits + post])
                MO_B[imo][iao] = val


    NMO = NMO_A - QMin['frozcore']
    if restr:
        NMO = NMO_A - QMin['frozcore']
    else:
        NMO = NMO_A + NMO_B - 2 * QMin['frozcore']

    # make string
    string = '''2mocoef
header
 1
MO-coefficients from Orca
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
    if not restr:
        x = 0
        for imo, mo in enumerate(MO_B):
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


def get_dets_from_cis(filename, QMin):

    # get general infos
    job = QMin['IJOB']
    restr = QMin['jobs'][job]['restr']
    mults = QMin['jobs'][job]['mults']
    gsmult = QMin['multmap'][-job][0]
    nstates_to_extract = deepcopy(QMin['states'])
    nstates_to_skip = [QMin['states_to_do'][i] - QMin['states'][i] for i in range(len(QMin['states']))]
    for i in range(len(nstates_to_extract)):
        if not i + 1 in mults:
            nstates_to_extract[i] = 0
            nstates_to_skip[i] = 0
        elif i + 1 == gsmult:
            nstates_to_extract[i] -= 1
    # print(job,restr,mults,gsmult,nstates_to_extract)

    # get infos from logfile
    logfile = os.path.join(os.path.dirname(filename), 'ORCA.log')
    data = readfile(logfile)
    infos = {}
    for iline, line in enumerate(data):
        if '# of contracted basis functions' in line:
            infos['nbsuse'] = int(line.split()[-1])
        if 'Orbital ranges used for CIS calculation:' in line:
            s = data[iline + 1].replace('.', ' ').split()
            infos['NFC'] = int(s[3])
            infos['NOA'] = int(s[4]) - int(s[3]) + 1
            infos['NVA'] = int(s[7]) - int(s[6]) + 1
            if restr:
                infos['NOB'] = infos['NOA']
                infos['NVB'] = infos['NVA']
            else:
                s = data[iline + 2].replace('.', ' ').split()
                infos['NOB'] = int(s[4]) - int(s[3]) + 1
                infos['NVB'] = int(s[7]) - int(s[6]) + 1

    if 'NOA' not in infos:
        nstates_onfile = 0
        charge = QMin['chargemap'][gsmult]
        nelec = float(QMin['Atomcharge'] - charge)
        infos['NOA'] = int(nelec / 2. + float(gsmult - 1) / 2.)
        infos['NOB'] = int(nelec / 2. - float(gsmult - 1) / 2.)
        infos['NVA'] = infos['nbsuse'] - infos['NOA']
        infos['NVB'] = infos['nbsuse'] - infos['NOB']
        infos['NFC'] = 0
    else:
        # get all info from cis file
        CCfile = open(filename, 'rb')
        nvec = struct.unpack('i', CCfile.read(4))[0]
        header = [struct.unpack('i', CCfile.read(4))[0] for i in range(8)]
        # print infos
        # print header
        if infos['NOA'] != header[1] - header[0] + 1:
            print('Number of orbitals in %s not consistent' % filename)
            sys.exit(82)
        if infos['NVA'] != header[3] - header[2] + 1:
            print('Number of orbitals in %s not consistent' % filename)
            sys.exit(83)
        if not restr:
            if infos['NOB'] != header[5] - header[4] + 1:
                print('Number of orbitals in %s not consistent' % filename)
                sys.exit(84)
            if infos['NVB'] != header[7] - header[6] + 1:
                print('Number of orbitals in %s not consistent' % filename)
                sys.exit(85)
        if QMin['template']['no_tda']:
            nstates_onfile = nvec // 2
        else:
            nstates_onfile = nvec


    # get ground state configuration
    # make step vectors (0:empty, 1:alpha, 2:beta, 3:docc)
    if restr:
        occ_A = [3 for i in range(infos['NFC'] + infos['NOA'])] + [0 for i in range(infos['NVA'])]
    if not restr:
        occ_A = [1 for i in range(infos['NFC'] + infos['NOA'])] + [0 for i in range(infos['NVA'])]
        occ_B = [2 for i in range(infos['NFC'] + infos['NOB'])] + [0 for i in range(infos['NVB'])]
    occ_A = tuple(occ_A)
    if not restr:
        occ_B = tuple(occ_B)

    # get infos
    nocc_A = infos['NOA']
    nvir_A = infos['NVA']
    nocc_B = infos['NOB']
    nvir_B = infos['NVB']

    # get eigenvectors
    eigenvectors = {}
    for imult, mult in enumerate(mults):
        eigenvectors[mult] = []
        if mult == gsmult:
            # add ground state
            if restr:
                key = tuple(occ_A[QMin['frozcore']:])
            else:
                key = tuple(occ_A[QMin['frozcore']:] + occ_B[QMin['frozcore']:])
            eigenvectors[mult].append({key: 1.0})
        for istate in range(nstates_to_extract[mult - 1]):
            CCfile.read(40)
            dets = {}
            for iocc in range(header[0], header[1] + 1):
                for ivirt in range(header[2], header[3] + 1):
                    dets[(iocc, ivirt, 1)] = struct.unpack('d', CCfile.read(8))[0]
            if not restr:
                for iocc in range(header[4], header[5] + 1):
                    for ivirt in range(header[6], header[7] + 1):
                        dets[(iocc, ivirt, 2)] = struct.unpack('d', CCfile.read(8))[0]
            if QMin['template']['no_tda']:
                CCfile.read(40)
                for iocc in range(header[0], header[1] + 1):
                    for ivirt in range(header[2], header[3] + 1):
                        dets[(iocc, ivirt, 1)] += struct.unpack('d', CCfile.read(8))[0]
                        dets[(iocc, ivirt, 1)] /= 2.
                if not restr:
                    for iocc in range(header[4], header[5] + 1):
                        for ivirt in range(header[6], header[7] + 1):
                            dets[(iocc, ivirt, 2)] += struct.unpack('d', CCfile.read(8))[0]
                            dets[(iocc, ivirt, 2)] /= 2.

            # pprint.pprint(dets)
            # truncate vectors
            norm = 0.
            for k in sorted(dets, key=lambda x: dets[x]**2, reverse=True):
                factor = 1.
                if norm > factor * QMin['wfthres']:
                    del dets[k]
                    continue
                norm += dets[k]**2
            # pprint.pprint(dets)
            # create strings and expand singlets
            dets2 = {}
            if restr:
                for iocc, ivirt, dummy in dets:
                    # singlet
                    if mult == 1:
                        # alpha excitation
                        key = list(occ_A)
                        key[iocc] = 2
                        key[ivirt] = 1
                        dets2[tuple(key)] = dets[(iocc, ivirt, dummy)] * math.sqrt(0.5)
                        # beta excitation
                        key[iocc] = 1
                        key[ivirt] = 2
                        dets2[tuple(key)] = dets[(iocc, ivirt, dummy)] * math.sqrt(0.5)
                    # triplet
                    elif mult == 3:
                        key = list(occ_A)
                        key[iocc] = 1
                        key[ivirt] = 1
                        dets2[tuple(key)] = dets[(iocc, ivirt, dummy)]
            else:
                for iocc, ivirt, dummy in dets:
                    if dummy == 1:
                        key = list(occ_A + occ_B)
                        key[iocc] = 0
                        key[ivirt] = 1
                        dets2[tuple(key)] = dets[(iocc, ivirt, dummy)]
                    elif dummy == 2:
                        key = list(occ_A + occ_B)
                        key[infos['NFC'] + nocc_A + nvir_A + iocc] = 0
                        key[infos['NFC'] + nocc_A + nvir_A + ivirt] = 2
                        dets2[tuple(key)] = dets[(iocc, ivirt, dummy)]
            # pprint.pprint(dets2)
            # remove frozen core
            dets3 = {}
            for key in dets2:
                problem = False
                if restr:
                    if any([key[i] != 3 for i in range(QMin['frozcore'])]):
                        problem = True
                else:
                    if any([key[i] != 1 for i in range(QMin['frozcore'])]):
                        problem = True
                    if any([key[i] != 2 for i in range(nocc_A + nvir_A + QMin['frozcore'], nocc_A + nvir_A + 2 * QMin['frozcore'])]):
                        problem = True
                if problem:
                    print('WARNING: Non-occupied orbital inside frozen core! Skipping ...')
                    continue
                    # sys.exit(86)
                if restr:
                    key2 = key[QMin['frozcore']:]
                else:
                    key2 = key[QMin['frozcore']:QMin['frozcore'] + nocc_A + nvir_A] + key[nocc_A + nvir_A + 2 * QMin['frozcore']:]
                dets3[key2] = dets2[key]
            # pprint.pprint(dets3)
            # append
            eigenvectors[mult].append(dets3)
        # skip extra roots
        for istate in range(nstates_to_skip[mult - 1]):
            CCfile.read(40)
            for iocc in range(header[0], header[1] + 1):
                for ivirt in range(header[2], header[3] + 1):
                    CCfile.read(8)
            if not restr:
                for iocc in range(header[4], header[5] + 1):
                    for ivirt in range(header[6], header[7] + 1):
                        CCfile.read(8)
            if QMin['template']['no_tda']:
                CCfile.read(40)
                for iocc in range(header[0], header[1] + 1):
                    for ivirt in range(header[2], header[3] + 1):
                        CCfile.read(8)
                if not restr:
                    for iocc in range(header[4], header[5] + 1):
                        for ivirt in range(header[6], header[7] + 1):
                            CCfile.read(8)


    strings = {}
    for imult, mult in enumerate(mults):
        filename = os.path.join(QMin['savedir'], 'dets.%i' % mult)
        strings[filename] = format_ci_vectors(eigenvectors[mult])

    return strings

# ======================================================================= #


def format_ci_vectors(ci_vectors):

    # get nstates, norb and ndets
    alldets = set()
    for dets in ci_vectors:
        for key in dets:
            alldets.add(key)
    ndets = len(alldets)
    nstates = len(ci_vectors)
    norb = len(next(iter(alldets)))

    string = '%i %i %i\n' % (nstates, norb, ndets)
    for det in sorted(alldets, reverse=True):
        for o in det:
            if o == 0:
                string += 'e'
            elif o == 1:
                string += 'a'
            elif o == 2:
                string += 'b'
            elif o == 3:
                string += 'd'
        for istate in range(len(ci_vectors)):
            if det in ci_vectors[istate]:
                string += ' %15.11f ' % ci_vectors[istate][det]
            else:
                string += ' %15.11f ' % 0.
        string += '\n'
    return string

# ======================================================================= #


def saveAOmatrix(WORKDIR, QMin):
    # filename=os.path.join(WORKDIR,'ORCA.log')
    # NAO,Smat=get_smat(filename)
    # filename=os.path.join(WORKDIR,'ORCA.molden.input')
    # NAO,Smat=get_smat_from_Molden(filename)
    filename = os.path.join(WORKDIR, 'ORCA.gbw')
    NAO, Smat = get_smat_from_gbw(filename)

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
# def get_smat(filename):

    # data=readfile(filename)

    # find MO block
    # iline=-1
    # NAO=0
    # while True:
        # iline+=1
        # if len(data)<=iline:
        # print('MOs not found!')
        # sys.exit(87)
        # line=data[iline]
        # if '# of contracted basis functions' in line:
        # NAO=int(line.split()[-1])
        # if 'OVERLAP MATRIX' in line:
        # break
    # if NAO==0:
        # print('Number of basis functions not found!')
        # sys.exit(88)
    # iline+=2

    # read matrix
    # nblock=6
    # ao_ovl=[ [ 0. for i in range(NAO) ] for j in range(NAO) ]
    # for i in range(NAO):
        # for j in range(NAO):
        # jline=iline + (i/nblock)*(NAO+1)+1+j
        # jcol =1+i%nblock
        # ao_ovl[i][j]=float(data[jline].split()[jcol])

    # return NAO,ao_ovl

# ======================================================================= #


def get_smat_from_gbw(file1, file2=''):

    if not file2:
        file2 = file1

    # run orca_fragovl
    string = 'orca_fragovl %s %s' % (file1, file2)
    try:
        proc = sp.Popen(string, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    except OSError:
        print('Call have had some serious problems:', OSError)
        sys.exit(89)
    comm = proc.communicate()[0].decode()
    out = comm.split('\n')

    # get size of matrix
    for line in reversed(out):
        # print line
        s = line.split()
        if len(s) >= 1:
            NAO = int(line.split()[0]) + 1
            break

    # find start of matrix
    iline = -1
    while True:
        iline += 1
        line = out[iline]
        if 'FRAGMENT-FRAGMENT OVERLAP MATRIX' in line:
            break

    # read matrix
    nblock = 6
    ao_ovl = [[0. for i in range(NAO)] for j in range(NAO)]
    for x in range(NAO):
        for y in range(NAO):
            block = x // nblock
            xoffset = x % nblock + 1
            yoffset = block * (NAO + 1) + y + 3 + iline
            ao_ovl[x][y] = float(out[yoffset].split()[xoffset])

    return NAO, ao_ovl


# ======================================================================= #
# def get_smat_from_Molden(file1, file2=''):

    # read file1
    # molecule=read_molden(file1)

    # read file2:
    # if file2:
    # molecule.extend(read_molden(file2))
    # pprint.pprint(molecule)

    # make PyQuante object
    # try:
    # from PyQuante.Ints import getS
    # from PyQuante.Basis.basis import BasisSet
    # from PyQuante.CGBF import CGBF
    # from PyQuante.PGBF import PGBF
    # from PyQuante import Molecule
    # from PyQuante.shell import Shell
    # except ImportError:
    # print('Could not import PyQuante!')
    # sys.exit(90)

    # class moldenBasisSet(BasisSet):
    # def __init__(self, molecule):
    # sym2powerlist = {
    # 'S' : [(0,0,0)],
    # 'P' : [(1,0,0),(0,1,0),(0,0,1)],
    # 'D' : [(2,0,0),(0,2,0),(0,0,2),(1,1,0),(1,0,1),(0,1,1)],
    # 'F' : [(3,0,0),(0,3,0),(0,0,3),(1,2,0),(2,1,0),(2,0,1),
    # (1,0,2),(0,1,2),(0,2,1), (1,1,1)],
    # 'G' : [(4,0,0),(0,4,0),(0,0,4),(3,1,0),(3,0,1),(1,3,0),
    # (0,3,1),(1,0,3),(0,1,3),(2,2,0),(2,0,2),(0,2,2),
    # (2,1,1),(1,2,1),(1,1,2)]
    # }


    # make molecule
    # atomlist=[]
    # for atom in molecule:
    # atomlist.append( (atom['el'], tuple(atom['coord'])) )
    # target=Molecule('default',atomlist,units='bohr')

    # self.bfs=[]
    # self.shells=[]
    # for iatom,atom in enumerate(molecule):
    # for bas in atom['basis']:
    # shell=Shell(bas[0])
    # for power in sym2powerlist[bas[0]]:
    # cgbf = CGBF(target[iatom].pos(), power, target[iatom].atid)
    # for alpha, coef in bas[1:]:
    # angular=sum(power)
    # coef*=alpha**(-(0.75+angular*0.5))  *  2**(-angular)  *  (2.0/math.pi)**(-0.75)
    # cgbf.add_primitive(alpha,coef)
    # cgbf.normalize()
    # self.bfs.append(cgbf)
    # shell.append(cgbf, len(self.bfs)-1)
    # self.shells.append(shell)

    # a=moldenBasisSet(molecule)
    # pprint.pprint(a.__dict__)

    # S=getS(a)
    # S=S.tolist()
    # return len(S),S


# ======================================================================= #
def read_molden(filename):
    data = readfile(filename)

    molecule = []
    # get geometry
    for iline, line in enumerate(data):
        if '[atoms]' in line.lower():
            break
    else:
        print('No geometry found in %s!' % (filename))
        sys.exit(91)

    if 'au' in line.lower():
        factor = 1.
    elif 'angstrom' in line.lower():
        factor = au2a

    while True:
        iline += 1
        line = data[iline]
        if '[' in line:
            break
        s = line.lower().split()
        atom = {'el': s[0], 'coord': [float(i) * factor for i in s[3:6]], 'basis': []}
        molecule.append(atom)

    # get basis set
    for iline, line in enumerate(data):
        if '[gto]' in line.lower():
            break
    else:
        print('No geometry found in %s!' % (filename))
        sys.exit(92)

    shells = {'s': 1, 'p': 3, 'd': 6, 'f': 10, 'g': 15}
    while True:
        iline += 1
        line = data[iline]
        if '[' in line:
            break
        s = line.lower().split()
        if len(s) == 0:
            continue
        if not s[0] in shells:
            iatom = int(s[0]) - 1
        else:
            newbf = [s[0].upper()]
            nprim = int(s[1])
            for iprim in range(nprim):
                iline += 1
                s = data[iline].split()
                newbf.append((float(s[0]), float(s[1])))
            molecule[iatom]['basis'].append(newbf)
    return molecule


# ======================================================================= #
def mkdir(DIR):
    # mkdir the DIR, or clean it if it exists
    if os.path.exists(DIR):
        if os.path.isfile(DIR):
            print('%s exists and is a file!' % (DIR))
            sys.exit(93)
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
            sys.exit(94)

# ======================================================================= #


def link(PATH, NAME, crucial=True, force=True):
    # do not create broken links
    if not os.path.exists(PATH) and crucial:
        print('Source %s does not exist, cannot create link!' % (PATH))
        sys.exit(95)
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
                    sys.exit(96)
                else:
                    return
    elif os.path.exists(NAME):
        # NAME is not a link. The interface will not overwrite files/directories with links, even with force=True
        print('%s exists, cannot create a link of the same name!' % (NAME))
        if crucial:
            sys.exit(97)
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
            else:
                mults = QMin['jobs'][ijob]['mults']
                gsmult = mults[0]
                ns = 0
                for i in mults:
                    ns += QMin['states'][i - 1] - (i == gsmult)
                if ns == 0:
                    if DEBUG:
                        print('Skipping Job %s because it contains no excited states.' % (ijob))
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
            sys.exit(98)

        print('')

    return errorcodes

# ======================================================================= #


def setupWORKDIR_TH(WORKDIR, QMin):
    # mkdir the WORKDIR, or clean it if it exists, then copy all necessary files from pwd and savedir

    # write dens_ana.in
    inputstring = '''rtype='cclib'
rfile='ORCA.log'
read_binary=True
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
    fromfile = os.path.join(WORKDIR, 'ORCA.cis')
    tofile = os.path.join(WORKDIR, 'orca.cis')
    link(fromfile, tofile)
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
    string = os.path.join(THEODIR, 'bin', 'analyze_tden.py')
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
        sys.exit(99)
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
        get_Double_AOovl_gbw(QMin)
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
        sys.exit(100)

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
        sys.exit(101)
    stdoutfile.close()
    stderrfile.close()
    if PRINT or DEBUG:
        endtime = datetime.datetime.now()
        sys.stdout.write('FINISH:\t%s\t%s\tRuntime: %s\tError Code: %i\n' % (shorten_DIR(WORKDIR), endtime, endtime - starttime, runerror))
        sys.stdout.flush()
    os.chdir(prevdir)
    return runerror


# ======================================================================= #
def get_Double_AOovl_gbw(QMin):

    # get geometries
    # filename1=os.path.join(QMin['savedir'],'ORCA.molden.1.old')
    # filename2=os.path.join(QMin['savedir'],'ORCA.molden.1')
    job = sorted(QMin['jobs'].keys())[0]
    filename1 = os.path.join(QMin['savedir'], 'ORCA.gbw.%i.old' % job)
    filename2 = os.path.join(QMin['savedir'], 'ORCA.gbw.%i' % job)

    #
    # NAO,Smat=get_smat_from_Molden(filename1,filename2)
    NAO, Smat = get_smat_from_gbw(filename1, filename2)

    # Smat is now full matrix NAO*NAO
    # we want the lower left quarter, but transposed
    # string='%i %i\n' % (NAO/2,NAO/2)
    # for irow in range(NAO/2,NAO):
    # for icol in range(0,NAO/2):
    # string+='% .15e ' % (Smat[icol][irow])          # note the exchanged indices => transposition
    # string+='\n'

    # Smat is already off-diagonal block matrix NAO*NAO
    # we want the lower left quarter, but transposed
    string = '%i %i\n' % (NAO, NAO)
    for irow in range(0, NAO):
        for icol in range(0, NAO):
            string += '% .15e ' % (Smat[irow][icol])          # note the exchanged indices => transposition
        string += '\n'
    filename = os.path.join(QMin['savedir'], 'AO_overl.mixed')
    writefile(filename, string)
    return













# =============================================================================================== #
# =============================================================================================== #
# ====================================== ORCA output parsing ================================ #
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
            logfile = os.path.join(QMin['scratchdir'], 'master_%i/ORCA.log' % (job))
            energies = getenergy(logfile, job, QMin)
            # print energies
            # also get SO matrix and mapping
            if 'soc' in QMin and QMin['jobs'][job]['restr']:
                submatrix, invstatemap = getsocm(logfile, job, QMin)
            mults = QMin['multmap'][-job]
            if 3 in mults and QMin['OrcaVersion'] < (4, 1):
                mults = [3]
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
                        x = invstatemap[(m1, s1, ms1)]
                        y = invstatemap[(m2, s2, ms2)]
                        QMout['h'][i][j] = submatrix[x - 1][y - 1]

    # Dipole Moments
    if 'dm' in QMin:
        # make matrix
        if 'dm' not in QMout:
            QMout['dm'] = [makecmatrix(nmstates, nmstates) for i in range(3)]
        # go through all jobs
        for job in joblist:
            logfile = os.path.join(QMin['scratchdir'], 'master_%i/ORCA.log' % (job))
            dipoles = gettdm(logfile, job, QMin)
            mults = QMin['multmap'][-job]
            if 3 in mults and QMin['OrcaVersion'] < (4, 1):
                mults = [3]
            for i in range(nmstates):
                m1, s1, ms1 = tuple(QMin['statemap'][i + 1])
                if m1 not in QMin['jobs'][job]['mults']:
                    continue
                for j in range(nmstates):
                    m2, s2, ms2 = tuple(QMin['statemap'][j + 1])
                    if m2 not in QMin['jobs'][job]['mults']:
                        continue
                    if i == j:
                        # TODO: does not work with restricted triplet
                        #isgs= (s1==1)
                        isgs = (QMin['gsmap'][i + 1] == i + 1)
                        if isgs:
                            logfile = os.path.join(QMin['scratchdir'], 'master_%i/ORCA.log' % (job))
                        elif (m1, s1) in QMin['gradmap']:
                            path, isgs = QMin['jobgrad'][(m1, s1)]
                            logfile = os.path.join(QMin['scratchdir'], path, 'ORCA.log')
                        else:
                            continue
                        dm = getdm(logfile, isgs)
                        for ixyz in range(3):
                            QMout['dm'][ixyz][i][j] = dm[ixyz]
                        continue
                    if not m1 == m2 == mults[0] or not ms1 == ms2:
                        continue
                    if s1 == 1:
                        for ixyz in range(3):
                            QMout['dm'][ixyz][i][j] = dipoles[(m2, s2)][ixyz]
                    elif s2 == 1:
                        for ixyz in range(3):
                            QMout['dm'][ixyz][i][j] = dipoles[(m1, s1)][ixyz]


    # Gradients
    if 'grad' in QMin:
        if 'grad' not in QMout:
            QMout['grad'] = [[[0. for i in range(3)] for j in range(natom)] for k in range(nmstates)]
        if QMin['qmmm'] and 'pcgrad' not in QMout:
            QMout['pcgrad'] = [[[0. for i in range(3)] for j in QMin['pointcharges']] for k in range(nmstates)]
        if QMin['template']['cobramm']:
            ncharges = len(readfile("charge.dat")) - 1
            QMout['pcgrad'] = [[[0. for i in range(3)] for j in range(ncharges)] for k in range(nmstates)]
        for grad in QMin['gradmap']:
            path, isgs = QMin['jobgrad'][grad]
            gsmult = QMin['jobs'][int(path.split('_')[1])]['mults'][0]
            restr = QMin['jobs'][int(path.split('_')[1])]['restr']
            if isgs:
                fname = '.ground'
                if QMin['states'][gsmult - 1] == 1:
                    fname = ''
            else:
                if restr:
                    fname = '.' + IToMult[grad[0]].lower() + '.root%i' % (grad[1] - (grad[0] == gsmult))
                else:
                    fname = '.singlet.root%i' % (grad[1] - (grad[0] == gsmult))
            logfile = os.path.join(QMin['scratchdir'], path, 'ORCA.engrad' + fname)
            g = getgrad(logfile, QMin)
            # print g
            if QMin['qmmm']:
                # TODO: seems ORCA.pcgrad.ground is not written in ORCA 5 
                if 'ground' in fname:
                    fname = ''
                logfile = os.path.join(QMin['scratchdir'], path, 'ORCA.pcgrad' + fname)
                gpc = getpcgrad(logfile, QMin)
            for istate in QMin['statemap']:
                state = QMin['statemap'][istate]
                # print grad,istate,state
                if (state[0], state[1]) == grad:
                    QMout['grad'][istate - 1] = g
                    if QMin['qmmm']:
                        QMout['pcgrad'][istate - 1] = gpc
            if QMin['template']['cobramm']:
                logfile = os.path.join(QMin['scratchdir'], path, 'ORCA.pcgrad' + fname)
                gpc = getpcgrad(logfile, QMin)
            for istate in QMin['statemap']:
                state = QMin['statemap'][istate]
                # print grad,istate,state
                if (state[0], state[1]) == grad:
                    QMout['grad'][istate - 1] = g
                    if QMin['template']['cobramm']:
                        QMout['pcgrad'][istate - 1] = gpc
            # if QMin['template']['cobramm']:
            # gradfromorca=os.path.join(QMin['scratchdir'],path,'ORCA.pcgrad')
            # gradforcobramm=os.path.join(QMin['scratchdir'],path,'grad_charges')
            # shutil.move(gradfromorca,gradforcobramm)
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
                    if QMin['qmmm']:
                        QMout['pcgrad'][i] = QMout['pcgrad'][j]

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
            else:
                mults = QMin['jobs'][job]['mults']
                gsmult = mults[0]
                ns = 0
                for i in mults:
                    ns += QMin['states'][i - 1] - (i == gsmult)
                if ns == 0:
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
    # if QMin['template']['qmmm']:
        # job=QMin['joblist'][0]
        #outfile=os.path.join(QMin['scratchdir'],'master_%i/ADF.out' % (job))
        # QMout['qmmm_energies']=get_qmmm_energies(outfile,QMin['template']['qmmm_coupling'])

    # transform back from QM to QM/MM
    if QMin['template']['qmmm']:
        QMin, QMout = transform_QM_QMMM(QMin, QMout)







    endtime = datetime.datetime.now()
    if PRINT:
        print("Readout Runtime: %s" % (endtime - starttime))

    if DEBUG:
        # pprint.pprint(QMout)
        copydir = os.path.join(QMin['savedir'], 'debug_ORCA_stdout')
        if not os.path.isdir(copydir):
            mkdir(copydir)
        for job in joblist:
            outfile = os.path.join(QMin['scratchdir'], 'master_%i/ORCA.log' % (job))
            shutil.copy(outfile, os.path.join(copydir, "ORCA_%i.log" % job))
            if QMin['jobs'][job]['restr'] and 'theodore' in QMin:
                outfile = os.path.join(QMin['scratchdir'], 'master_%i/tden_summ.txt' % job)
                try:
                    shutil.copy(outfile, os.path.join(copydir, 'THEO_%i.out' % (job)))
                except IOError:
                    pass
                outfile = os.path.join(QMin['scratchdir'], 'master_%i/OmFrag.txt' % job)
                try:
                    shutil.copy(outfile, os.path.join(copydir, 'THEO_OMF_%i.out' % (job)))
                except IOError:
                    pass
        if 'grad' in QMin:
            for grad in QMin['gradmap']:
                path, isgs = QMin['jobgrad'][grad]
                outfile = os.path.join(QMin['scratchdir'], path, 'ORCA.log')
                shutil.copy(outfile, os.path.join(copydir, "ORCA_GRAD_%i_%i.log" % grad))
        if 'overlap' in QMin:
            for mult in itmult(QMin['states']):
                job = QMin['multmap'][mult]
                outfile = os.path.join(QMin['scratchdir'], 'WFOVL_%i_%i/wfovl.out' % (mult, job))
                shutil.copy(outfile, os.path.join(copydir, 'WFOVL_%i_%i.out' % (mult, job)))
        if 'ion' in QMin:
            for ion in QMin['ionmap']:
                outfile = os.path.join(QMin['scratchdir'], 'Dyson_%i_%i_%i_%i/wfovl.out' % ion)
                shutil.copy(outfile, os.path.join(copydir, 'Dyson_%i_%i_%i_%i.out' % ion))

    if QMin['save_stuff']:
        copydir = os.path.join(QMin['savedir'], 'save_stuff')
        if not os.path.isdir(copydir):
            mkdir(copydir)
        for job in joblist:
            outfile = os.path.join(QMin['scratchdir'], 'master_%i/ORCA.log' % (job))
            shutil.copy(outfile, os.path.join(copydir, "ORCA_%i.log" % job))
            outfile = os.path.join(QMin['scratchdir'], 'master_%i/ORCA.gbw' % (job))
            shutil.copy(outfile, os.path.join(copydir, "ORCA_%i.gbw" % job))
            outfile = os.path.join(QMin['scratchdir'], 'master_%i/ORCA.cis' % (job))
            if os.path.isfile(outfile):
                shutil.copy(outfile, os.path.join(copydir, "ORCA_%i.cis" % job))



    return QMin, QMout

# ======================================================================= #


def getenergy(logfile, ijob, QMin):

    # open file
    f = readfile(logfile)
    if PRINT:
        print('Energy:   ' + shorten_DIR(logfile))

    # read ground state
    for iline, line in enumerate(f):
        if 'TOTAL SCF ENERGY' in line:
            gsenergy = float(f[iline + 3].split()[3])
        if 'Dispersion correction' in line:
            gsenergy += float(line.split()[-1])
            break

    # figure out the excited state settings
    mults = QMin['jobs'][ijob]['mults']
    restr = QMin['jobs'][ijob]['restr']
    gsmult = mults[0]
    estates_to_extract = deepcopy(QMin['states'])
    estates_to_extract[gsmult - 1] -= 1
    for imult in range(len(estates_to_extract)):
        if not imult + 1 in mults:
            estates_to_extract[imult] = 0
    for imult in range(len(estates_to_extract)):
        if imult + 1 in mults:
            estates_to_extract[imult] = max(estates_to_extract)
    # TODO: restricted triplet energies

    # extract excitation energies
    # loop also works if no energies should be extracted
    energies = {(gsmult, 1): gsenergy}
    for imult in mults:
        nstates = estates_to_extract[imult - 1]
        # print nstates
        if nstates > 0:
            strings = [['TD-DFT/TDA',
                        'TD-DFT',
                        'RPA',
                        'CIS'],
                       ['EXCITED STATES']
                       ]
            if QMin['OrcaVersion'] >= (4, 1):
                if restr:
                    if imult == 1:
                        strings.append(['SINGLETS'])
                    if imult == 3:
                        strings.append(['TRIPLETS'])
            for iline, line in enumerate(f):
                if all([any([i in line for i in st]) for st in strings]):
                    # print line
                    # if 'TD-DFT/TDA EXCITED STATES' in line or 'TD-DFT EXCITED STATES' in line or 'RPA EXCITED STATES' in line or 'CIS-EXCITED STATES' in line:
                    # if QMin['OrcaVersion']>=(4,1):
                    break
            finalstring = ['Entering ', '-EXCITATION SPECTRA']
            while True:
                iline += 1
                if iline >= len(f):
                    print('Error in parsing excitation energies')
                    sys.exit(102)
                line = f[iline]
                if any([i in line for i in finalstring]):
                    break
                if 'STATE' in line:
                    # print line
                    s = line.replace(':', ' ').split()
                    ikey = 0
                    while True:
                        ikey -= 1
                        if 'cm**-1' in s[ikey]:
                            break
                    e = gsenergy+float(s[ikey-1])*rcm_to_Eh
                    i = int(s[1])
                    if i > nstates:
                        break
                    energies[(imult, i + (gsmult == imult))] = e
    return energies

## ======================================================================= #


def getsocm(outfile, ijob, QMin):

    # read the standard out into memory
    out = readfile(outfile)
    if PRINT:
        print('SOC:      ' + shorten_DIR(outfile))



    # get number of states (nsing=ntrip in Orca)
    for line in out:
        if 'Number of roots to be determined' in line:
            nst = int(line.split()[-1])
            break
    nrS = nst
    nrT = nst

    # make statemap for the state ordering of the SO matrix
    inv_statemap = {}
    inv_statemap[(1, 1, 0.0)] = 1
    i = 1
    for x in range(nrS):
        i += 1
        inv_statemap[(1, x + 2, 0.0)] = i
    spin = [0.0, -1.0, +1.0]
    for y in range(3):
        for x in range(nrT):
            i += 1
            inv_statemap[(3, x + 1, spin[y])] = i
    #pprint.pprint( inv_statemap)

    # get matrix
    iline = -1
    while True:
        iline += 1
        line = out[iline]
        if 'The full SOC matrix' in line:
            break
    iline += 5
    ncol = 6
    real = [[0 + 0j for i in range(4 * nst + 1)] for j in range(4 * nst + 1)]
    for x in range(len(real)):
        for y in range(len(real[0])):
            block = x // ncol
            xoffset = 1 + x % ncol
            yoffset = block * (4 * nst + 2) + y
            # print iline,x,y,block,xoffset,yoffset
            val = float(out[iline + yoffset].split()[xoffset])
            if abs(val) > 1e-16:
                real[y][x] = val

    iline += ((4 * nst) // ncol + 1) * (4 * nst + 2) + 2
    for x in range(len(real)):
        for y in range(len(real[0])):
            block = x // ncol
            xoffset = 1 + x % ncol
            yoffset = block * (4 * nst + 2) + y
            val = float(out[iline + yoffset].split()[xoffset])
            if abs(val) > 1e-16:
                real[y][x] += (0 + 1j) * val

    #pprint.pprint( real)


    return real, inv_statemap

# ======================================================================= #


def gettdm(logfile, ijob, QMin):

    # open file
    f = readfile(logfile)
    if PRINT:
        print('Dipoles:  ' + shorten_DIR(logfile))

    # figure out the excited state settings
    mults = QMin['jobs'][ijob]['mults']
    if 3 in mults and QMin['OrcaVersion'] < (4, 1):
        mults = [3]
    restr = QMin['jobs'][ijob]['restr']
    gsmult = mults[0]
    estates_to_extract = deepcopy(QMin['states'])
    estates_to_extract[gsmult - 1] -= 1
    for imult in range(len(estates_to_extract)):
        if not imult + 1 in mults:
            estates_to_extract[imult] = 0

    # print "getting cool dipoles"
    # extract transition dipole moments
    dipoles = {}
    for imult in mults:
        if not imult == gsmult:
            continue
        nstates = estates_to_extract[imult - 1]
        if nstates > 0:
            for iline, line in enumerate(f):
                if '  ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS' in line:
                    # print line
                    for istate in range(nstates):
                        shift = 5 + istate
                        s = f[iline + shift].split()
                        dm = [float(i) for i in s[5:8]]
                        dipoles[(imult, istate + 1 + (gsmult == imult))] = dm
    # print dipoles
    return dipoles

# ======================================================================= #


def getdm(logfile, isgs):

    # open file
    f = readfile(logfile)
    if PRINT:
        print('Dipoles:  ' + shorten_DIR(logfile))

    if isgs:
        findstring = 'ORCA ELECTRIC PROPERTIES CALCULATION'
    else:
        findstring = '*** CIS RELAXED DENSITY ***'
    for iline, line in enumerate(f):
        if findstring in line:
            break
    while True:
        iline += 1
        line = f[iline]
        if 'Total Dipole Moment' in line:
            s = line.split()
            dmx = float(s[4])
            dmy = float(s[5])
            dmz = float(s[6])
            dm = [dmx, dmy, dmz]
            return dm

# ======================================================================= #


def getgrad(logfile, QMin):

    # initialize
    natom = QMin['natom']
    g = [[0. for i in range(3)] for j in range(natom)]

    # read file
    if os.path.isfile(logfile):
        out = readfile(logfile)
        if PRINT:
            print('Gradient: ' + shorten_DIR(logfile))

        # get gradient
        string = 'The current gradient in Eh/bohr'
        shift = 2
        for iline, line in enumerate(out):
            if string in line:
                for iatom in range(natom):
                    for ixyz in range(3):
                        s = out[iline + shift + 3 * iatom + ixyz]
                        g[iatom][ixyz] = float(s)

    # read binary file otherwise
    else:
        logfile += '.grad.tmp'
        Gfile = open(logfile, 'rb')
        if PRINT:
            print('Gradient: ' + shorten_DIR(logfile))

        # get gradient
        Gfile.read(8 + 28 * natom)    # skip header
        for iatom in range(natom):
            for ixyz in range(3):
                f = struct.unpack('d', Gfile.read(8))[0]
                g[iatom][ixyz] = f

    return g

# ======================================================================= #


def getgrad_from_log(logfile, QMin):

    # read file
    out = readfile(logfile)
    if PRINT:
        print('Gradient: ' + shorten_DIR(logfile))

    # initialize
    natom = QMin['natom']
    g = [[0. for i in range(3)] for j in range(natom)]

    # find gradients
    iline = -1
    while True:
        iline += 1
        line = out[iline]
        if 'ORCA SCF GRADIENT CALCULATION' in line:
            break


    return g

# ======================================================================= #


def getpcgrad(logfile, QMin):

    # read file
    out = readfile(logfile)
    if PRINT:
        print('Gradient: ' + shorten_DIR(logfile))

    # initialize
    # natom=len(QMin['pointcharges'])
    # g=[ [ 0. for i in range(3) ] for j in range(natom) ]

    # get gradient
    # for iatom in range(natom):
    #     for ixyz in range(3):
    #       s=out[iatom+1].split()
    #       g[iatom][ixyz]=float(s[ixyz])
    # return g
    g = []
    for iatom in range(len(out) - 1):
        atom_grad = [0. for i in range(3)]
        s = out[iatom + 1].split()
        for ixyz in range(3):
            atom_grad[ixyz] = float(s[ixyz])
        g.append(atom_grad)
    return g

## ======================================================================= #
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
            sys.exit(103)
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
            sys.exit(104)
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

    # Retrieve PRINT and DEBUG
    try:
        envPRINT = os.getenv('SH2Orc_PRINT')
        if envPRINT and envPRINT.lower() == 'false':
            global PRINT
            PRINT = False
        envDEBUG = os.getenv('SH2Orc_DEBUG')
        if envDEBUG and envDEBUG.lower() == 'true':
            global DEBUG
            DEBUG = True
    except ValueError:
        print('PRINT or DEBUG environment variables do not evaluate to numerical values!')
        sys.exit(105)

    # Process Command line arguments
    if len(sys.argv) != 2:
        print('Usage:\n./SHARC_ORCA.py <QMin>\n')
        print('version:', version)
        print('date:', versiondate)
        print('changelog:\n', changelogstring)
        sys.exit(106)
    QMinfilename = sys.argv[1]

    # Print header
    printheader()

    # Read QMinfile
    QMin = readQMin(QMinfilename)

    # get the job schedule
    QMin, schedule = generate_joblist(QMin)
    printQMin(QMin)
    if DEBUG:
        pprint.pprint(schedule, depth=3)

    # run all the ADF jobs
    errorcodes = runjobs(schedule, QMin)

    # do all necessary overlap and Dyson calculations
    errorcodes = run_wfoverlap(QMin, errorcodes)

    # do all necessary Theodore calculations
    errorcodes = run_theodore(QMin, errorcodes)

    # read all the output files
    QMin, QMout = getQMout(QMin)
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
        cleandir(QMin['scratchdir'])
        if 'cleanup' in QMin:
            cleandir(QMin['savedir'])

    print
    print(datetime.datetime.now())
    print('#================ END ================#')


if __name__ == '__main__':
    main()






# kate: indent-width 2
