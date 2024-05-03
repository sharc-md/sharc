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

#    ====================================================================
# ||                                                                       ||
# ||                General Remarks                                        ||
# ||                                                                       ||
#    ====================================================================
#
# This script uses several different specification for the electronic states under consideration.
# Generally, the input specs are like "3 Singlets, 0 Doublets and 3 Triplets"
# Based on this information, the states may be denoted by different systems.
#
# The most comprehensive denotation is: (mult, state, MS).
# In this case, states of higher multiplicity show up several times and the total number of states may be quite large.
# This total number of states is called nmstates in the script (equals to 12 for the above example)
#
# Since the MS components of a multiplet often share several expectation values, these need to be calculated only once.
# In this case, states can be safely denoted by (mult,state).
# The total number of these states is called nstates.
#
# In both systems, the states can be given indices. In this script and in SHARC, one first counts up the state quantum number,
# then the MS quantum number, and finally the multiplicity, like in this code snippet:
#
# i=0
# for mult in range(len(states)):
#     for MS in range(mult):
#         for state in range(states[mult]):
#             i+=1
#             print i, mult+1, state+1, MS-i/2
#
# more about this below in the docstrings of the iterator functions

# ======================================================================= #

# IMPLEMENTATION OF ADDITIONAL TASKS KEYWORDS, JOBS, ETC:
#
# A new task keyword in QMin has to be added to:
#             - readQMin (for consistency check)
#             - gettasks (planning the MOLCAS calculation)
#             - print QMin (optional)
#
# A new task in the Tasks array needs changes in:
#             - gettasks
#             - writeMOLCASinput
#             - redotasks
#             - printtasks

# Modules:
# Operating system, isfile and related routines, move files, create directories
import os
# External Calls to MOLPRO
import subprocess as sp
# Command line arguments
import sys
# shell utilities like copy
import shutil
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
# gethostname routine
from socket import gethostname
# reading binary files
import struct
import copy
import ast

# ======================================================================= #

version = '3.0'
versiondate = datetime.date(2023, 4, 1)


changelogstring = '''
12.02.2016:     Initial version 0.1
- CC2 and ADC(2) from Turbomole
- SOC from Orca call
- only singlets and triplets
  => Doublets and Dyson could be added later

15.03.2016:     0.1.1
- ridft can be used for the SCF calculation, but dscf has to be run afterwards anyways (for SOC and overlaps).
- Laplace-transformed SOS-CC2/ADC(2) can be used ("spin-scaling lt-sos"), but does not work for soc or trans-dm
- if ricc2 does not converge, it is rerun with different combinations of npre and nstart (until convergence or too many tries)

07.03.2017:
- wfthres is now interpreted as in other interfaces (default is 0.99 now)

25.04.2017:
- can use external basis set libraries

23.08.2017
- Resource file is now called "RICC2.resources" instead of "SH2CC2.inp" (old filename still works)

24.08.2017:
- numfrozcore in resources file can now be used to override number of frozen cores for overlaps
- added Theodore capabilities (compute descriptors, OmFrag, and NTOs (also activate MOLDEN key for that))

11.11.2020:
- COBRAMM can be used for QM/MM calculations
'''

# ======================================================================= #
# holds the system time when the script was started
starttime = datetime.datetime.now()

# global variables for printing (PRINT gives formatted output, DEBUG gives raw output)
DEBUG = False
PRINT = True

# hash table for conversion of multiplicity to the keywords used in COLUMBUS
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

# hash table for conversion of polarisations to the keywords used in COLUMBUS
IToPol = {
    0: 'X',
    1: 'Y',
    2: 'Z',
    'X': 0,
    'Y': 1,
    'Z': 2
}


NUMBERS = {'H': 1, 'He': 2,
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

BASISSETS = [
    'SV',
    'SVP',
    'SV(P)',
    'def-SVP',
    'def2-SVP',
    'dhf-SVP',
    'dhf-SVP-2c',
    'def-SV(P)',
    'def2-SV(P)',
    'dhf-SV(P)',
    'dhf-SV(P)-2c',
    'DZ',
    'DZP',
    'TZ',
    'TZP',
    'TZV',
    'TZVP',
    'def-TZVP',
    'TZVE',
    'TZVEP',
    'TZVPP',
    'def-TZVPP',
    'def2-TZVP',
    'dhf-TZVP',
    'dhf-TZVP-2c',
    'def2-TZVPP',
    'dhf-TZVPP',
    'dhf-TZVPP-2c',
    'TZVPPP',
    'QZV',
    'def-QZV',
    'def2-QZV',
    'QZVP',
    'def-QZVP',
    'def2-QZVP',
    'dhf-QZVP',
    'dhf-QZVP-2c',
    'QZVPP',
    'def-QZVPP',
    'def2-QZVPP',
    'dhf-QZVPP',
    'dhf-QZVPP-2c',
    'minix',
    'sto-3ghondo',
    '4-31ghondo',
    '6-31ghondo',
    '3-21ghondo',
    'dzphondo',
    'tzphondo',
    '6-31G',
    '6-31G*',
    '6-31G**',
    '6-311G',
    '6-311G*',
    '6-311G**',
    '6-311++G**',
    'cc-pVDZ',
    'cc-pV(D+d)Z',
    'aug-cc-pVDZ',
    'aug-cc-pV(D+d)Z',
    'YP-aug-cc-pVDZ',
    'cc-pwCVDZ',
    'aug-cc-pwCVDZ',
    'cc-pVDZ-sp',
    'cc-pVTZ',
    'cc-pV(T+d)Z',
    'aug-cc-pVTZ',
    'aug-cc-pV(T+d)Z',
    'YP-aug-cc-pVTZ',
    'cc-pwCVTZ',
    'aug-cc-pwCVTZ',
    'cc-pVTZ-sp',
    'cc-pVQZ',
    'cc-pV(Q+d)Z',
    'aug-cc-pVQZ',
    'aug-cc-pV(Q+d)Z',
    'YP-aug-cc-pVQZ',
    'cc-pwCVQZ',
    'aug-cc-pwCVQZ',
    'cc-pVQZ-sp',
    'cc-pV5Z',
    'cc-pV(5+d)Z',
    'aug-cc-pV5Z',
    'aug-cc-pV(5+d)Z',
    'YP-aug-cc-pV5Z',
    'cc-pwCV5Z',
    'aug-cc-pwCV5Z',
    'cc-pV5Z-large-s',
    'cc-pV6Z',
    'cc-pV(6+d)Z',
    'aug-cc-pV6Z',
    'aug-cc-pV(6+d)Z',
    'cc-pV6Z-sp',
    'cc-pVDZ-F12',
    'cc-pVTZ-F12',
    'cc-pVQZ-F12',
    'def2-SVPD',
    'def2-TZVPPD',
    'def2-TZVPD',
    'def2-QZVPPD',
    'def2-QZVPD'
]

# conversion factors
au2a = 0.529177211
rcm_to_Eh = 4.556335e-6
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
        f.close()
    except IOError:
        print('Could not write to file %s!' % (filename))
        sys.exit(13)

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
    if PRINT or DEBUG:
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
    string += '||' + ' ' * 27 + 'SHARC - RICC2 - Interface' + ' ' * 28 + '||\n'
    string += '||' + ' ' * 80 + '||\n'
    string += '||' + ' ' * 29 + 'Author: Sebastian Mai' + ' ' * 30 + '||\n'
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

    if DEBUG:
        pprint.pprint(QMin)
    if not PRINT:
        return
    print('==> QMin Job description for:\n%s' % (QMin['comment']))

    string = 'Tasks:  '
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
        string += '\tOverlaps'
    if 'angular' in QMin:
        string += '\tAngular'
    if 'ion' in QMin:
        string += '\tDyson norms'
    if 'dmdr' in QMin:
        string += '\tDM-Grad'
    if 'socdr' in QMin:
        string += '\tSOC-Grad'
    if 'theodore' in QMin:
        string += '\tTheoDORE'
    if 'phases' in QMin:
        string += '\tPhases'
    print(string)

    string = 'States: '
    for i in itmult(QMin['states']):
        string += '\t%i %s' % (QMin['states'][i - 1], IToMult[i])
    print(string)



    string = 'Method: \t'
    if QMin['template']['spin-scaling'] != 'none':
        string += QMin['template']['spin-scaling'].upper() + '-'
    string += QMin['template']['method'].upper()
    string += '/%s' % (QMin['template']['basis'])
    parts = []
    if QMin['template']['douglas-kroll']:
        parts.append('Douglas-Kroll')
    if QMin['template']['frozen'] != -1:
        parts.append('RICC2 frozen orbitals=%i' % (QMin['template']['frozen']))
    if QMin['template']['scf'] == 'ridft':
        parts.append('RI-SCF')
    if len(parts) > 0:
        string += '\t('
        string += ','.join(parts)
        string += ')'
    print(string)


    # if 'dm' in QMin and QMin['template']['method']=='adc(2)':
    # print('WARNING: excited-to-excited transition dipole moments in ADC(2) are zero!')


    if 'dm' in QMin and QMin['template']['method'] == 'cc2' and (QMin['states'][0] > 1 and len(QMin['states']) > 2 and QMin['states'][2] > 0):
        print('WARNING: will not calculate transition dipole moments! For CC2, please use only singlet states or ground state + triplets.')


    string = 'Found Geo'
    if 'veloc' in QMin:
        string += ' and Veloc! '
    else:
        string += '! '
    string += 'NAtom is %i.\n' % (QMin['natom'])
    print(string)

    string = '\nGeometry in Bohrs:\n'
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

    # if 'overlap' in QMin:
        # string='Overlaps:\n'
        # for i in range(1,QMin['nmstates']+1):
        # for j in range(1,QMin['nmstates']+1):
        # if [i,j] in QMin['overlap'] or [j,i] in QMin['overlap']:
        # string+='X '
        # else:
        # string+='. '
        # string+='\n'
        # print(string)

    for i in QMin:
        if not any([i == j for j in ['h', 'dm', 'soc', 'dmdr', 'socdr', 'theodore', 'geo', 'veloc', 'states', 'comment', 'LD_LIBRARY_PATH', 'grad', 'nacdr', 'ion', 'overlap', 'template', 'qmmm', 'cobramm', 'geo_orig', 'pointcharges']]):
            if not any([i == j for j in ['ionlist', 'ionmap']]) or DEBUG:
                string = i + ': '
                string += str(QMin[i])
                print(string)
        else:
            string = i + ': ...'
            print(string)
    print('\n')
    sys.stdout.flush()

# ======================================================================= #


def printcomplexmatrix(matrix, states):
    '''Prints a formatted matrix. Zero elements are not printed, blocks of different mult and MS are delimited by dashes. Also prints a matrix with the imaginary parts, if any one element has non-zero imaginary part.

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
            string += '% .5f\t' % (grad[atom][xyz])
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
    print('===> Results:\n')
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
# =========================================== output extraction ================================= #
# =============================================================================================== #
# =============================================================================================== #

# ======================================================================= #
def get_RICC2out(QMin, QMout, job):
    # job contains: 'tmexc_soc','tmexc_dm', 'spectrum','exprop_dm','static_dm','gsgrad','exgrad', 'E'
    # reads ricc2.out and adds matrix elements to QMout
    ricc2 = readfile(QMin['scratchdir'] + '/JOB/ricc2.out')

    states = QMin['states']
    nstates = QMin['nstates']
    nmstates = QMin['nmstates']
    natom = QMin['natom']

    # Hamiltonian
    if 'h' in QMin or 'soc' in QMin:
        # make Hamiltonian
        if 'h' not in QMout:
            QMout['h'] = makecmatrix(nmstates, nmstates)
        # read the matrix elements from ricc2.out
        for i in range(nmstates):
            for j in range(nmstates):
                m1, s1, ms1 = tuple(QMin['statemap'][i + 1])
                m2, s2, ms2 = tuple(QMin['statemap'][j + 1])
                if 'soc' in QMin and 'tmexc_soc' in job:
                    # read SOC elements
                    QMout['h'][i][j] = getsocme(ricc2, QMin, i + 1, j + 1)
                elif 'E' in job:
                    # read only diagonal elements
                    if i == j:
                        QMout['h'][i][j] = getenergy(ricc2, QMin, i + 1)

    # Dipole moments
    if 'dm' in QMin:
        # make dipole matrices
        if 'dm' not in QMout:
            QMout['dm'] = []
            for i in range(3):
                QMout['dm'].append(makecmatrix(nmstates, nmstates))
        # read the elements from ricc2.out
        for i in range(nmstates):
            for j in range(nmstates):
                m1, s1, ms1 = tuple(QMin['statemap'][i + 1])
                m2, s2, ms2 = tuple(QMin['statemap'][j + 1])
                if m1 != m2:
                    continue
                if ms1 != ms2:
                    continue
                if (m1, s1) == (m2, s2) == (1, 1):
                    if 'static_dm' in job or ('gsgrad', 1, 1) in job:
                        for xyz in range(3):
                            QMout['dm'][xyz][i][j] = getdiagdm(ricc2, QMin, i + 1, xyz)
                elif i == j:
                    if 'exprop_dm' in job or ('exgrad', m1, s1) in job:
                        for xyz in range(3):
                            QMout['dm'][xyz][i][j] = getdiagdm(ricc2, QMin, i + 1, xyz)
                elif (m1, s1) == (1, 1) or (m2, s2) == (1, 1):
                    if 'spectrum' in job:
                        for xyz in range(3):
                            QMout['dm'][xyz][i][j] = gettransdm(ricc2, QMin, i + 1, j + 1, xyz)
                else:
                    if 'tmexc_dm' in job:
                        for xyz in range(3):
                            QMout['dm'][xyz][i][j] = gettransdm(ricc2, QMin, i + 1, j + 1, xyz)

    # Gradients
    if 'grad' in QMin:
        # make vectors
        if 'grad' not in QMout:
            QMout['grad'] = [[[0., 0., 0.] for i in range(natom)] for j in range(nmstates)]
        if QMin['qmmm'] and 'pcgrad' not in QMout:
            QMout['pcgrad'] = [[[0. for i in range(3)] for j in QMin['pointcharges']] for k in range(nmstates)]
        if QMin['cobramm'] and 'pcgrad' not in QMout:
            ncharges = len(readfile('charge.dat'))  # -4
            # print(ncharges,"tot")
            QMout['pcgrad'] = [[[0. for i in range(3)] for j in range(ncharges)] for k in range(nmstates)]
            print(QMout['pcgrad'])
#      pcgrad=os.path.join(QMin['scratchdir'],'JOB','pc_grad')
#      specify_state=os.path.join(QMin['scratchdir'],'JOB','pc_grad.old.%s') % (nmstates)#% (mult,nexc)
        # shutil.copy(pcbrad,specify_state)
#      shutil.move(pcgrad,specify_state)
#    if QMin['cobramm'] and not 'pcgrad' in QMout: ##21.09.20
#      QMout['pcgrad']=[ [ [ 0. for i in range(3) ] for j in QMin['pointcharges'] ] for k in range(nmstates) ]
        # read the elements from ricc2.out
        for i in range(nmstates):
            m1, s1, ms1 = tuple(QMin['statemap'][i + 1])
            if (m1, s1) == (1, 1):
                tup = ('gsgrad', 1, 1)
            else:
                tup = ('exgrad', m1, s1)
            if tup in job:
                QMout['grad'][i] = getgrad(ricc2, QMin, i + 1)
                if QMin['qmmm']:
                    QMout['pcgrad'][i - 1] = getpcgrad(QMin)
                if QMin['cobramm']:
                    logfile = os.path.join(QMin['scratchdir'], 'JOB', 'pc_grad')
                    getcobrammpcgrad(logfile, QMin)
                    # gpc=getcobrammpcgrad(logfile,QMin)
                    # QMout['pcgrad'][i]=gpc
                    # print(QMout['pcgrad'][i])
                    pcgradold = os.path.join(QMin['scratchdir'], 'JOB', 'pc_grad')
                    specify_state = os.path.join(QMin['scratchdir'], 'JOB', 'pc_grad.%s.%s') % (m1, s1)  # % (mult,nexc)
        # shutil.copy(pcbrad,specify_state)
                    shutil.copy(pcgradold, specify_state)
        if QMin['neglected_gradient'] != 'zero' and 'samestep' not in QMin:
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



    return QMout

# ======================================================================= #


def getenergy(ricc2, QMin, istate):
    mult, state, ms = tuple(QMin['statemap'][istate])

    # ground state energy
    if QMin['template']['method'] == 'cc2':
        string = 'Final CC2 energy'
    elif QMin['template']['method'] == 'adc(2)':
        string = 'Final MP2 energy'
    for line in ricc2:
        if string in line:
            e = float(line.split()[5])
            break
    else:
        print('"%s" not found in ricc2.out' % (string))
        sys.exit(14)

    # return gs energy if requested
    if mult == 1 and state == 1:
        return e

    # excited state
    if QMin['template']['method'] == 'cc2':
        string = '| sym | multi | state |          CC2 excitation energies       |  %t1   |  %t2   |'
    elif QMin['template']['method'] == 'adc(2)':
        string = '| sym | multi | state |          ADC(2) excitation energies    |  %t1   |  %t2   |'

    # find correct table
    iline = -1
    while True:
        iline += 1
        if iline == len(ricc2):
            print('"%s" not found in ricc2.out' % (string))
            sys.exit(15)
        line = ricc2[iline]
        if string in line:
            break
    iline += 3

    # find correct line
    if mult == 1:
        iline += state - 1
    elif mult == 3:
        if QMin['states'][0] >= 2:
            iline += QMin['states'][0]
        iline += state

    # get energy
    line = ricc2[iline]
    e = e + float(line.split()[7])
    return complex(e, 0.)

# ======================================================================= #


def getsocme(ricc2, QMin, istate, jstate):

    if istate == jstate:
        return getenergy(ricc2, QMin, istate)

    if istate not in QMin['statemap'] or jstate not in QMin['statemap']:
        print('States %i or %i are not in statemap!' % (istate, jstate))
        sys.exit(16)

    m1, s1, ms1 = tuple(QMin['statemap'][istate])
    m2, s2, ms2 = tuple(QMin['statemap'][jstate])

    # RICC2 does not calculate triplet-triplet SOC, singlet-singlet SOC is zero
    if m1 == m2:
        return complex(0., 0.)

    # RICC2 does not calculate SOC with S0
    if (m1, s1) == (1, 1) or (m2, s2) == (1, 1):
        return complex(0., 0.)

    # find the correct table
    string = '|        States           Operator    Excitation                 Transition             |'
    iline = 0
    while True:
        iline += 1
        if iline == len(ricc2):
            print('"%s" not found in ricc2.out' % (string))
            sys.exit(17)
        line = ricc2[iline]
        if string in line:
            break
    iline += 7

    # find the correct line
    if istate > jstate:
        m1, s1, ms1, m2, s2, ms2 = m2, s2, ms2, m1, s1, ms1
    nsing = QMin['states'][0] - 1
    ntrip = QMin['states'][2]
    ntot = nsing + ntrip
    if m1 == 1:
        x1 = s1 - 1
    else:
        x1 = nsing + s1
    if m2 == 1:
        x2 = s2 - 1
    else:
        x2 = nsing + s2
    for i in range(1, 1 + ntot):
        for j in range(i + 1, 1 + ntot):
            iline += 1
            if i == x1 and j == x2:
                try:
                    s = ricc2[iline].split()
                    idx = int(10 + 2 * ms2)
                    soc = float(s[idx]) * rcm_to_Eh
                except IndexError:
                    print('Could not find SOC matrix element with istate=%i, jstate=%i, line=%i' % (istate, jstate, iline))
                return complex(soc, 0.)

# ======================================================================= #


def getdiagdm(ricc2, QMin, istate, pol):
    # finds and returns state dipole moments in ricc2.out
    m1, s1, ms1 = tuple(QMin['statemap'][istate])

    if (m1, s1) == (1, 1):
        start1string = ''
        start2string = '<<<<<<<<<<  GROUND STATE FIRST-ORDER PROPERTIES  >>>>>>>>>>>'
    else:
        start1string = '<<<<<<<<<<<<<<<  EXCITED STATE PROPERTIES  >>>>>>>>>>>>>>>>'
        start2string = 'number, symmetry, multiplicity:  % 3i a    %i' % (s1 - (m1 == 1), m1)
    findstring = '     dipole moment:'
    stopstring = 'Analysis of unrelaxed properties'

    # find correct section
    iline = -1
    while True:
        iline += 1
        if iline == len(ricc2):
            print('Could not find dipole moment of istate=%i, Fail=0' % (istate))
            sys.exit(18)
        line = ricc2[iline]
        if start1string in line:
            break

    # find correct state
    while True:
        iline += 1
        if iline == len(ricc2):
            print('Could not find dipole moment of istate=%i, Fail=1' % (istate))
            sys.exit(19)
        line = ricc2[iline]
        if start2string in line:
            break

    # find correct line
    while True:
        iline += 1
        if iline == len(ricc2):
            print('Could not find dipole moment of istate=%i, Fail=2' % (istate))
            sys.exit(20)
        line = ricc2[iline]
        if stopstring in line:
            print('Could not find dipole moment of istate=%i, Fail=3' % (istate))
            sys.exit(21)
        if findstring in line:
            break

    iline += 3 + pol
    s = ricc2[iline].split()
    dm = float(s[1])
    return complex(dm, 0.)

# ======================================================================= #


def gettransdm(ricc2, QMin, istate, jstate, pol):
    # finds and returns transition dipole moments in ricc2.out
    m1, s1, ms1 = tuple(QMin['statemap'][istate])
    m2, s2, ms2 = tuple(QMin['statemap'][jstate])

    if istate > jstate:
        m1, s1, ms1, m2, s2, ms2 = m2, s2, ms2, m1, s1, ms1

    # ground state to excited state
    if (m1, s1) == (1, 1):
        start1string = '<<<<<<<<<<<<  ONE-PHOTON ABSORPTION STRENGTHS  >>>>>>>>>>>>>'
        start2string = 'number, symmetry, multiplicity:  % 3i a    %i' % (s2 - (m2 == 1), m2)
        stopstring = '<<<<<<<<<<<<<<<  EXCITED STATE PROPERTIES  >>>>>>>>>>>>>>>>'

        # find correct section
        iline = -1
        while True:
            iline += 1
            if iline == len(ricc2):
                print('Could not find transition dipole moment of istate=%i,jstate=%i, Fail=0' % (istate, jstate))
                sys.exit(22)
            line = ricc2[iline]
            if start1string in line:
                break

        # find correct state
        while True:
            iline += 1
            if iline == len(ricc2):
                print('Could not find transition dipole moment of istate=%i,jstate=%i, Fail=1' % (istate, jstate))
                sys.exit(23)
            line = ricc2[iline]
            if stopstring in line:
                print('Could not find transition dipole moment of istate=%i,jstate=%i, Fail=2' % (istate, jstate))
                sys.exit(24)
            if start2string in line:
                break

        # find element
        iline += 7 + pol
        s = ricc2[iline].split()
        dm = 0.5 * (float(s[3]) + float(s[5]))
        return dm

    # excited to excited state
    else:
        start1string = '<<<<<<<<<<<  EXCITED STATE TRANSITION MOMENTS  >>>>>>>>>>>>'
        start2string = 'Transition moments for pair  % 3i a      % 3i a' % (s1 - (m1 == 1), s2 - (m2 == 1))
        stopstring = 'Model:'
        nostring = 'Transition and Operator of different multiplicity.'

        # invert search for triplets
        if m1 == 3:
            start1string, stopstring = stopstring, start1string

        # find correct section
        iline = -1
        while True:
            iline += 1
            if iline == len(ricc2):
                print('Could not find transition dipole moment of istate=%i,jstate=%i, Fail=4' % (istate, jstate))
                sys.exit(25)
            line = ricc2[iline]
            if start1string in line:
                break

        # find correct state
        while True:
            if m1 == 1:
                # search forward for singlet-singlet transitions
                iline += 1
            elif m1 == 3:
                # search backward from the end for triplet-triplet transitions
                iline += -1
            if iline + 2 == len(ricc2) or iline == -1:
                print('Could not find transition dipole moment of istate=%i,jstate=%i, Fail=5' % (istate, jstate))
                sys.exit(26)
            line = ricc2[iline]
            line2 = ricc2[iline + 2]
            if stopstring in line:
                print('Could not find transition dipole moment of istate=%i,jstate=%i, Fail=6' % (istate, jstate))
                sys.exit(27)
            if start2string in line and nostring not in line2:
                break

        # find element
        iline += 2 + pol
        s = ricc2[iline].split()
        dm = 0.5 * (float(s[1]) + float(s[2]))
        return dm

    print('Could not find transition dipole moment of istate=%i,jstate=%i, Fail=7' % (istate, jstate))
    sys.exit(28)


# ======================================================================= #
def getgrad(ricc2, QMin, istate):
    m1, s1, ms1 = tuple(QMin['statemap'][istate])
    natom = QMin['natom']

    if (m1, s1) == (1, 1):
        start1string = '<<<<<<<<<<  GROUND STATE FIRST-ORDER PROPERTIES  >>>>>>>>>>>'
        stop1string = '<<<<<<<<<<<<<<<  EXCITED STATE PROPERTIES  >>>>>>>>>>>>>>>>'
    else:
        start1string = '<<<<<<<<<<<<<<<  EXCITED STATE PROPERTIES  >>>>>>>>>>>>>>>>'
        stop1string = 'total wall-time'
    findstring = 'cartesian gradient of the energy (hartree/bohr)'

    # find correct section
    iline = -1
    while True:
        iline += 1
        if iline == len(ricc2):
            print('Could not find gradient of istate=%i, Fail=0' % (istate))
            sys.exit(29)
        line = ricc2[iline]
        if start1string in line:
            break
        if stop1string in line:
            print('Could not find gradient of istate=%i, Fail=1' % (istate))
            sys.exit(30)

    # find gradient
    while True:
        iline += 1
        if iline == len(ricc2):
            print('Could not find gradient of istate=%i, Fail=2' % (istate))
            sys.exit(31)
        line = ricc2[iline]
        if findstring in line:
            break
    iline += 3

    # get grad
    grad = []
    col = 0
    row = 0
    for iatom in range(natom):
        atom = []
        col += 1
        if col > 5:
            col = 1
            row += 1
        for xyz in range(3):
            line = ricc2[iline + 5 * row + xyz + 1]
            el = float(line.split()[col].replace('D', 'E'))
            atom.append(el)
        grad.append(atom)
    return grad

# ======================================================================= #


def getpcgrad(QMin):
    pcgrad = readfile(os.path.join(QMin['scratchdir'], 'JOB', 'pc_grad'))

    grad = []
    iline = 0
    for ipc, pc in enumerate(QMin['pointcharges']):
        if pc[-1] != 0:
            g = []
            iline += 1
            s = pcgrad[iline].replace('D', 'E').split()
            for i in range(3):
                g.append(float(s[i]))
            grad.append(g)
    return grad
# ======================================================================= #


def getcobrammpcgrad(logfile, QMin):

    #  # read file
    #  out=readfile(logfile)
    #  if PRINT:
    #      print('Gradient: '+shorten_DIR(logfile))

    #   gradpc=[]
    #   iline=0
    #   ncharges=(readfile(('charge.dat') ))
    # for ipc,pc in enumerate('charge.dat'):
    #   for iatom in range(len('charge.dat')):
    #     if pc[-1]!=0:
    #       g=[]
    #       iline+=1
    #       s=out[iline].replace('D','E').split()
    #       for i in range(3):
    #         g.append(float(s[i]) )
    #       gradpc.append(g)
    #   print(gradpc)
    #   return gradpc
    # #  gradpc=[]
    #  # for iatom in range(len(out)-1):
    # #      atom_grad=[0. for i in range(3)]
    # s=out[iatom+1].replace('D','E').split()
    # #      for ixyz in range(3):
    # if ixyz != '$end':
    # continue
    # atom_grad[ixyz]=float(s[ixyz])
    # gradpc.append(atom_grad)
    # print(gradpc)
    # return gradpc

    out = readfile(logfile)

    ncharges = len(out)
    gradpc = []
    out.pop(0)
    out.pop(-1)
    ncharges = len(out)
    string = ''
    # icharge=0
    for pc in range(ncharges):
        q = []
        xyz = out[pc].replace('D', 'E').split()
        # icharge+=1
        for i in range(3):
            q.append(float(xyz[i]))
        gradpc.append(q)
    filecharges = open("grad_charges", "a")
    string += '%i %i !\n' % (ncharges, 3)
    # string+='%i %i ! %i %i %i\n' % (natom,3,imult,istate,ims)
    for atom in range(ncharges):
        for xyz in range(3):
            string += '%s ' % (eformat((gradpc[atom][xyz]), 9, 3))
        string += "\n"
    filecharges.write(string)

# ======================================================================= #


def writeQMoutgradcobramm(QMin, QMout):
    '''Generates a string with the Gradient vectors in SHARC format.

    The string starts with a ! followed by a flag specifying the type of data. On the next line, natom and 3 are written, followed by      the gradient, with one line per atom and a blank line at the end. Each MS component shows up (nmstates gradients are written).

    Arguments:
    1 dictionary: QMin
    2 dictionary: QMout

    Returns:
    1 string: multiline string with the Gradient vectors'''
    ncharges = len(readfile(os.path.join(QMin['scratchdir'], 'JOB', 'pc_grad'))) - 2
    states = QMin['states']
    nstates = QMin['nstates']
    nmstates = QMin['nmstates']
    natom = len(QMout['pcgrad'][0])
    print(QMout['pcgrad'][1])
    string = ''
    print(natom)
    # string+='! %i Gradient Vectors (%ix%ix3, real)\n' % (3,nmstates,natom)
    i = 0
    for imult, istate, ims in itnmstates(states):
        string += '%i %i ! %i %i %i\n' % (natom, 3, imult, istate, ims)
        for atom in range(natom):
            for xyz in range(3):
                print((QMout['pcgrad'][i][atom], 9, 3), i, atom)
                string += '%s ' % (eformat(QMout['pcgrad'][i][atom][xyz], 9, 3))
            string += '\n'
        # string+='\n'
        i += 1
    string += '\n'
    writefile("grad_charges", string)


# ======================================================================= #
def get_smatel(out, s1, s2):
    ilines = -1
    while True:
        ilines += 1
        if ilines == len(out):
            print('Overlap of states %i - %i not found!' % (s1, s2))
            sys.exit(32)
        if containsstring('Overlap matrix <PsiA_i|PsiB_j>', out[ilines]):
            break
    ilines += 1 + s1
    f = out[ilines].split()
    return float(f[s2 + 1])

# ======================================================================= #


def get_wfovlout(QMin, QMout, path, mult):
    outfile = os.path.join(path, 'wfovl.out')
    out = readfile(outfile)

    if 'overlap' in QMin:
        nmstates = QMin['nmstates']
        if 'overlap' not in QMout:
            QMout['overlap'] = makecmatrix(nmstates, nmstates)
        # read the overlap matrix
        for i in range(nmstates):
            for j in range(nmstates):
                m1, s1, ms1 = tuple(QMin['statemap'][i + 1])
                m2, s2, ms2 = tuple(QMin['statemap'][j + 1])
                if not m1 == m2 == mult:
                    continue
                if not ms1 == ms2:
                    continue
                QMout['overlap'][i][j] = get_smatel(out, s1, s2)

    return QMout
# ======================================================================= #
# def get_dysonel(out,s1,s2):
    # ilines=-1
    # while True:
    # ilines+=1
    # if ilines==len(out):
    # print('Overlap of states %i - %i not found!' % (s1,s2))
    # sys.exit(33)
    # if containsstring('Dyson norm matrix |<PsiA_i|PsiB_j>|^2', out[ilines]):
    # break
    # ilines+=1+s1
    # f=out[ilines].split()
    # return float(f[s2+1])

# ======================================================================= #
# def get_DYSONoutput(QMin,QMout,imult,jmult):
    # with open(QMin['scratchdir']+'/OVERLAP/cioverlap.out') as f:
    # out=f.readlines()
    # out=readfile(QMin['scratchdir']+'/DYSON/dyson.out')

    # if 'ion' in QMin:
    # nmstates=QMin['nmstates']
    # if not 'prop' in QMout:
    # QMout['prop']=makecmatrix(nmstates,nmstates)
    # read the Dyson norm matrix
    # for i in range(nmstates):
    # for j in range(nmstates):
    # m1,s1,ms1=tuple(QMin['statemap'][i+1])
    # m2,s2,ms2=tuple(QMin['statemap'][j+1])
    # if not (imult,jmult)==(m1,m2) and not (imult,jmult)==(m2,m1):
    # continue
    # if not abs(ms1-ms2)==0.5:
    # continue
    # if m1>m2:
    # s1,s2=s2,s1
    # QMout['prop'][i][j]=get_dysonel(out,s1,s2)

    # return QMout


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
    # TODO: writeQMoutsocdr is not defined
    if 'socdr' in QMin:
        string += writeQMoutsocdr(QMin, QMout)
    # TODO: writeQMoutdmdr is not defined
    if 'dmdr' in QMin:
        string += writeQMoutdmdr(QMin, QMout)
    if 'ion' in QMin:
        string += writeQMoutprop(QMin, QMout)
    if 'theodore' in QMin or QMin['template']['qmmm']:
        string += writeQMoutTHEODORE(QMin, QMout)
    if 'phases' in QMin:
        string += writeQmoutPhases(QMin, QMout)
    # if  QMin['template']['cobramm']:
    #  writeQMoutgradcobramm(QMin,QMout)
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
            string += '%s %s ' % (eformat(QMout['h'][i][j].real, 12, 3), eformat(QMout['h'][i][j].imag, 12, 3))
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
                string += '%s %s ' % (eformat(QMout['dm'][xyz][i][j].real, 12, 3), eformat(QMout['dm'][xyz][i][j].imag, 12, 3))
            string += '\n'
        string += ''
    return string

# ======================================================================= #


def writeQMoutang(QMin, QMout):
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
    string += '! %i Angular Momentum Matrices (3x%ix%i, complex)\n' % (9, nmstates, nmstates)
    for xyz in range(3):
        string += '%i %i\n' % (nmstates, nmstates)
        for i in range(nmstates):
            for j in range(nmstates):
                string += '%s %s ' % (eformat(QMout['angular'][xyz][i][j].real, 12, 3), eformat(QMout['angular'][xyz][i][j].imag, 12, 3))
            string += '\n'
        string += ''
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
        string += '%i %i ! m1 %i s1 %i ms1 %i\n' % (natom, 3, imult, istate, ims)
        for atom in range(natom):
            for xyz in range(3):
                string += '%s ' % (eformat(QMout['grad'][i][atom][xyz], 12, 3))
            string += '\n'
        string += ''
        i += 1
    return string

# ======================================================================= #


def writeQMoutnacnum(QMin, QMout):
    '''Generates a string with the NAC matrix in SHARC format.

    The string starts with a ! followed by a flag specifying the type of data. In the next line, the dimensions of the matrix are given, followed by nmstates blocks of nmstates elements. Blocks are separated by a blank line.

    Arguments:
    1 dictionary: QMin
    2 dictionary: QMout

    Returns:
    1 string: multiline string with the NAC matrix'''

    states = QMin['states']
    nstates = QMin['nstates']
    nmstates = QMin['nmstates']
    natom = QMin['natom']
    string = ''
    string += '! %i Non-adiabatic couplings (ddt) (%ix%i, complex)\n' % (4, nmstates, nmstates)
    string += '%i %i\n' % (nmstates, nmstates)
    for i in range(nmstates):
        for j in range(nmstates):
            string += '%s %s ' % (eformat(QMout['nacdt'][i][j].real, 12, 3), eformat(QMout['nacdt'][i][j].imag, 12, 3))
        string += '\n'
    string += ''
    # also write wavefunction phases
    string += '! %i Wavefunction phases (%i, complex)\n' % (7, nmstates)
    for i in range(nmstates):
        string += '%s %s\n' % (eformat(QMout['phases'][i], 12, 3), eformat(0., 12, 3))
    string += '\n\n'
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
            # string+='%i %i ! %i %i %i %i %i %i\n' % (natom,3,imult,istate,ims,jmult,jstate,jms)
            string += '%i %i ! m1 %i s1 %i ms1 %i   m2 %i s2 %i ms2 %i\n' % (natom, 3, imult, istate, ims, jmult, jstate, jms)
            for atom in range(natom):
                for xyz in range(3):
                    string += '%s ' % (eformat(QMout['nacdr'][i][j][atom][xyz], 12, 3))
                string += '\n'
            string += ''
            j += 1
        i += 1
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
            string += '%s %s ' % (eformat(QMout['overlap'][j][i].real, 12, 3), eformat(QMout['overlap'][j][i].imag, 12, 3))
        string += '\n'
    string += '\n'
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

# ======================================================================= #


def writeQMoutprop(QMin, QMout):
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
            sys.exit(34)
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
        sys.exit(35)
    # no MM atom is allowed to be bonded to two QM atoms
    if not len(mm_in_links) == len(set(mm_in_links)):
        print('Some MM atom is involved in more than one link bond!')
        sys.exit(36)
    # no neighboring MM atoms are allowed to be involved in link bonds
    if not len(mm_in_link_neighbors) == len(set(mm_in_link_neighbors)):
        print('An MM-link atom is bonded to another MM-link atom!')
        sys.exit(37)


    # check geometry and connection table
    if not QMMM['natom_table'] == QMin['natom']:
        print('Number of atoms in table file does not match number of atoms in QMin!')
        sys.exit(38)


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
        sys.exit(39)

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

    # make list of pointcharges without QM atoms and zero-charge MM atoms
    print('Finalizing charges ...         ', datetime.datetime.now())
    QMMM['pointcharges'] = []
    QMMM['reorder_pc_input'] = {}
    ipc = 0
    for iatom_input in QMMM['MMpc']:
        q = QMMM['MMpc'][iatom_input]
        if q != 0:
            atom = QMMM['MM_coords'][iatom_input]
            QMMM['pointcharges'].append(atom[1:4] + [q])
            QMMM['reorder_pc_input'][ipc] = iatom_input
            ipc += 1
    # for iatom_input in QMMM['MM_atoms']:
            # atom=QMMM['MM_coords'][iatom_input]
            # q=QMMM['MMpc'][iatom_input]
            # QMMM['pointcharges'].append( atom[1:4]+[q] )
            # QMMM['reorder_pc_input'][ipc]=iatom
            # ipc+=1






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
    # print( QMMM['MMGradient'] )
    # print('MM coord:')
    # print(str(QMMM['natom_table']) +'\n')
    # for atom in QMMM['MM_coords']:
    # print(atom[0], atom[1], atom[2], atom[3])
    # print('QM coord:')
    # print(str(len(QMMM['QM_coords'])) +'\n')
    # for atom in QMMM['QM_coords']:
    # print(atom[0], atom[1]*au2a, atom[2]*au2a, atom[3]*au2a)
    # print('MM pc:')
    # for iatom,atom in enumerate(QMMM['MM_coords']):
    # q=QMMM['MMpc'][iatom]
    # if q!=0.:
    # print(atom[1],atom[2],atom[3],q)
    # pprint.pprint(QMMM['MMEnergy_terms'])
    # print('='*60)

    # pprint.pprint(QMMM)
    # sys.exit(40)
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
        sys.exit(41)
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
    # if strip and not DEBUG and runerror==0:
    # stripWORKDIR(WORKDIR)
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
        # for iqm,iqmmm in enumerate(QMin['qmmm']['MM_atoms']):
        for iqm in QMin['qmmm']['reorder_pc_input']:
            iqmmm = QMin['qmmm']['reorder_pc_input'][iqm]
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

def get_turbomole_version(turbodir):
    ls = os.listdir(turbodir)
    v = 0.0
    for f in ls:
        if "TURBOMOLE_" in f:
            s=f.split('_')
            v = int(s[1])/10
            break
    print('Turbomole version:',v)
    return v


def get_arch(turbodir):
    os.environ['TURBODIR'] = turbodir
    string = os.path.join(turbodir, 'scripts', 'sysname')
    proc = sp.Popen([string], stdout=sp.PIPE)
    output = proc.communicate()[0].strip().decode("utf-8")
    print('Architecture: %s' % output)
    return output

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
            sys.exit(42)
    else:
        try:
            os.makedirs(SCRATCHDIR)
        except OSError:
            print('Can not create SCRATCHDIR=%s\n' % (SCRATCHDIR))
            sys.exit(43)

# ======================================================================= #


def removequotes(string):
    if string.startswith("'") and string.endswith("'"):
        return string[1:-1]
    elif string.startswith('"') and string.endswith('"'):
        return string[1:-1]
    else:
        return string

# ======================================================================= #


def getsh2cc2key(sh2cc2, key):
    i = -1
    while True:
        i += 1
        try:
            line = re.sub(r'#.*$', '', sh2cc2[i])
        except IndexError:
            break
        line = line.strip().split(None, 1)
        if line == []:
            continue
        if key.lower() in line[0].lower():
            return line
    return ['', '']

# ======================================================================= #


def get_sh2cc2_environ(sh2cc2, key, environ=True, crucial=True):
    line = getsh2cc2key(sh2cc2, key)
    if line[0]:
        LINE = line[1]
    else:
        if environ:
            LINE = os.getenv(key.upper())
            if not LINE:
                print('Either set $%s or give path to %s in SH2COL.inp!' % (key.upper(), key.upper()))
                if crucial:
                    sys.exit(44)
                else:
                    return None
        else:
            print('Give path to %s in SH2COL.inp!' % (key.upper()))
            if crucial:
                sys.exit(45)
            else:
                return None
    LINE = os.path.expandvars(LINE)
    LINE = os.path.expanduser(LINE)
    LINE = os.path.abspath(LINE)
    LINE = removequotes(LINE).strip()
    if containsstring(';', LINE):
        print("$%s contains a semicolon. Do you probably want to execute another command after %s? I can't do that for you..." % (key.upper(), key.upper()))
        sys.exit(46)
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
            sys.exit(47)
        if 'end' in line:
            break
        fields = line.split()
        try:
            nacpairs.append([int(fields[0]), int(fields[1])])
        except ValueError:
            print('"nacdr select" is followed by pairs of state indices, each pair on a new line!')
            sys.exit(48)
    return nacpairs, i

# =============================================================================================== #
# =============================================================================================== #
# =========================================== readQMin and gettasks ============================= #
# =============================================================================================== #
# =============================================================================================== #


# ======================================================================= #     OK
def readQMin(QMinfilename):
    '''Reads the time-step dependent information from QMinfilename. This file contains all information from the current SHARC job: geometry, velocity, number of states, requested quantities along with additional information. The routine also checks this input and obtains a number of environment variables necessary to run COLUMBUS.

    Reads also the information from SH2COL

    Steps are:
    - open and read QMinfilename
    - Obtain natom, comment, geometry (, velocity)
    - parse remaining keywords from QMinfile
    - check keywords for consistency, calculate nstates, nmstates
    - obtain environment variables for path to COLUMBUS and scratch directory, and for error handling

    Arguments:
    1 string: name of the QMin file

    Returns:
    1 dictionary: QMin'''

    # read QMinfile
    QMinlines = readfile(QMinfilename)
    QMin = {}

    # Get natom
    try:
        natom = int(QMinlines[0])
    except ValueError:
        print('first line must contain the number of atoms!')
        sys.exit(49)
    QMin['natom'] = natom
    if len(QMinlines) < natom + 4:
        print('Input file must contain at least:\nnatom\ncomment\ngeometry\nkeyword "states"\nat least one task')
        sys.exit(50)

    # Save Comment line
    QMin['comment'] = QMinlines[1]

    # Get geometry and possibly velocity (for backup-analytical non-adiabatic couplings)
    QMin['geo'] = []
    QMin['veloc'] = []
    hasveloc = True
    for i in range(2, natom + 2):
        if not containsstring('[a-zA-Z][a-zA-Z]?[0-9]*.*[-]?[0-9]+[.][0-9]*.*[-]?[0-9]+[.][0-9]*.*[-]?[0-9]+[.][0-9]*', QMinlines[i]):
            print('Input file does not comply to xyz file format! Maybe natom is just wrong.')
            sys.exit(51)
        fields = QMinlines[i].split()
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
        if len(args) >= 1 and 'select' in args[0]:
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
            sys.exit(52)
    else:
        factor = 1. / au2a

    for iatom in range(len(QMin['geo'])):
        for ixyz in range(3):
            QMin['geo'][iatom][ixyz + 1] *= factor


    if 'states' not in QMin:
        print('Keyword "states" not given!')
        sys.exit(53)

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
        sys.exit(54)

    possibletasks = ['h', 'soc', 'dm', 'grad', 'overlap', 'dmdr', 'socdr', 'ion', 'theodore', 'phases']
    if not any([i in QMin for i in possibletasks]):
        print('No tasks found! Tasks are "h", "soc", "dm", "grad","dmdr", "socdr", "overlap" and "ion".')
        sys.exit(55)

    if 'samestep' in QMin and 'init' in QMin:
        print('"Init" and "Samestep" cannot be both present in QM.in!')
        sys.exit(56)

    if 'phases' in QMin:
        QMin['overlap'] = []

    if 'overlap' in QMin and 'init' in QMin:
        print('"overlap" and "phases" cannot be calculated in the first timestep! Delete either "overlap" or "init"')
        sys.exit(57)

    if 'init' not in QMin and 'samestep' not in QMin and 'restart' not in QMin:
        QMin['newstep'] = []

    if not any([i in QMin for i in ['h', 'soc', 'dm', 'grad']]) and 'overlap' in QMin:
        QMin['h'] = []

    if len(QMin['states']) > 3:
        print('Higher multiplicities than triplets are not supported!')
        sys.exit(58)

    if len(QMin['states']) > 1 and QMin['states'][1] > 0:
        print('No doublet states supported currently!')
        sys.exit(59)

    if len(QMin['states']) == 1 and 'soc' in QMin:
        QMin = removekey(QMin, 'soc')
        QMin['h'] = []

    if 'h' in QMin and 'soc' in QMin:
        QMin = removekey(QMin, 'h')

    if 'nacdt' in QMin or 'nacdr' in QMin:
        print('Within the SHARC-RICC2 interface, couplings can only be calculated via the overlap method. "nacdr" and "nacdt" are not supported.')
        sys.exit(60)

    if 'socdr' in QMin or 'dmdr' in QMin:
        print('Within the SHARC-RICC2 interface, "dmdr" and "socdr" are not supported.')
        sys.exit(61)

    if 'ion' in QMin:
        print('Ionization probabilities not implemented!')
        sys.exit(62)

    if 'step' not in QMin:
        QMin['step'] = [0]

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
                    sys.exit(63)
                if QMin['grad'][i] > nmstates:
                    print('State for requested gradient does not correspond to any state in QM input file state list!')
                    sys.exit(64)

    # Process the overlap requests
    # identically to the nac requests
    if 'overlap' in QMin:
        if len(QMin['overlap']) >= 1:
            nacpairs = QMin['overlap']
            for i in range(len(nacpairs)):
                if nacpairs[i][0] > nmstates or nacpairs[i][1] > nmstates:
                    print('State for requested non-adiabatic couplings does not correspond to any state in QM input file state list!')
                    sys.exit(65)
        else:
            QMin['overlap'] = [[j + 1, i + 1] for i in range(nmstates) for j in range(i + 1)]

    # obtain the statemap
    statemap = {}
    i = 1
    for imult, istate, ims in itnmstates(QMin['states']):
        statemap[i] = [imult, istate, ims]
        i += 1
    QMin['statemap'] = statemap

    # get the set of states for which gradients actually need to be calculated
    gradmap = set()
    if 'grad' in QMin:
        for i in QMin['grad']:
            gradmap.add(tuple(statemap[i][0:2]))
    gradmap = sorted(gradmap)
    QMin['gradmap'] = gradmap








    # --------------------------------------------- Resources ----------------------------------


    # environment setup

    QMin['pwd'] = os.getcwd()

    # open RICC2.resources
    filename = 'RICC2.resources'
    if os.path.isfile(filename):
        sh2cc2 = readfile(filename)
    else:
        print('HINT: reading resources from SH2CC2.inp')
        sh2cc2 = readfile('SH2CC2.inp')

    # ncpus for SMP-parallel turbomole and wfoverlap
    # this comes before the turbomole path determination
    QMin['ncpu'] = 1
    line = getsh2cc2key(sh2cc2, 'ncpu')
    if line[0]:
        try:
            QMin['ncpu'] = int(line[1])
            QMin['ncpu'] = max(1, QMin['ncpu'])
        except ValueError:
            print('Number of CPUs does not evaluate to numerical value!')
            sys.exit(66)
    os.environ['OMP_NUM_THREADS'] = str(QMin['ncpu'])
    if QMin['ncpu'] > 1:
        os.environ['PARA_ARCH'] = 'SMP'
        os.environ['PARNODES'] = str(QMin['ncpu'])


    # set TURBOMOLE paths
    QMin['turbodir'] = get_sh2cc2_environ(sh2cc2, 'turbodir')
    os.environ['TURBODIR'] = QMin['turbodir']
    arch = get_arch(QMin['turbodir'])
    version = get_turbomole_version(QMin['turbodir'])
    QMin['arch'] = arch
    QMin['tmversion'] = version
    #print('Turbomole version detected: %f    Turbomole architecture detected: %s' % (version, arch) )
    os.environ['PATH'] = '%s/scripts:%s/bin/%s:' % (QMin['turbodir'], QMin['turbodir'], arch) + os.environ['PATH']
    # print('Added to PATH:', '%s/scripts:%s/bin/%s:' % (QMin['turbodir'], QMin['turbodir'], arch))

    # set ORCA paths
    if 'soc' in QMin:
        QMin['orcadir'] = get_sh2cc2_environ(sh2cc2, 'orcadir')
        os.environ['PATH'] = '%s:' % (QMin['orcadir']) + os.environ['PATH']
        os.environ['LD_LIBRARY_PATH'] = '%s:' % (QMin['orcadir']) + os.environ['LD_LIBRARY_PATH']


    # Set up scratchdir
    line = get_sh2cc2_environ(sh2cc2, 'scratchdir', False, False)
    if line is None:
        line = QMin['pwd'] + '/SCRATCHDIR/'
    line = os.path.expandvars(line)
    line = os.path.expanduser(line)
    line = os.path.abspath(line)
    # checkscratch(line)
    QMin['scratchdir'] = line


    # Set up savedir
    if 'savedir' in QMin:
        # savedir may be read from QM.in file
        line = QMin['savedir'][0]
    else:
        line = get_sh2cc2_environ(sh2cc2, 'savedir', False, False)
        if line is None:
            line = QMin['pwd'] + '/SAVEDIR/'
    line = os.path.expandvars(line)
    line = os.path.expanduser(line)
    line = os.path.abspath(line)
    if 'init' in QMin:
        checkscratch(line)
    QMin['savedir'] = line


    # debug keyword in SH2CC2
    line = getsh2cc2key(sh2cc2, 'debug')
    if line[0]:
        if len(line) <= 1 or 'true' in line[1].lower():
            global DEBUG
            DEBUG = True

    line = getsh2cc2key(sh2cc2, 'no_print')
    if line[0]:
        if len(line) <= 1 or 'true' in line[1].lower():
            global PRINT
            PRINT = False


    # memory for Turbomole, Orca, and wfoverlap
    QMin['memory'] = 100
    line = getsh2cc2key(sh2cc2, 'memory')
    if line[0]:
        try:
            QMin['memory'] = int(line[1])
            QMin['memory'] = max(100, QMin['memory'])
        except ValueError:
            print('Run memory does not evaluate to numerical value!')
            sys.exit(67)
    else:
        print('WARNING: Please set memory in RICC2.resources (in MB)! Using 100 MB default value!')


    # initial MO guess settings
    # if neither keyword is present, the interface will reuse MOs from savedir, or use the EHT guess
    line = getsh2cc2key(sh2cc2, 'always_orb_init')
    if line[0]:
        QMin['always_orb_init'] = []
    line = getsh2cc2key(sh2cc2, 'always_guess')
    if line[0]:
        QMin['always_guess'] = []
    if 'always_orb_init' in QMin and 'always_guess' in QMin:
        print('Keywords "always_orb_init" and "always_guess" cannot be used together!')
        sys.exit(68)


    # get the nooverlap keyword: no dets will be extracted if present
    line = getsh2cc2key(sh2cc2, 'nooverlap')
    if line[0]:
        if 'overlap' in QMin or 'phases' in QMin or 'ion' in QMin:
            print('"nooverlap" is incompatible with "overlap" or "phases"!')
            sys.exit(69)
        QMin['nooverlap'] = []


    # dipole moment calculation level
    QMin['dipolelevel'] = 2
    line = getsh2cc2key(sh2cc2, 'dipolelevel')
    if line[0]:
        try:
            QMin['dipolelevel'] = int(line[1])
        except ValueError:
            print('Run memory does not evaluate to numerical value!')
            sys.exit(70)


    # wfoverlaps setting
    QMin['wfthres'] = 0.99
    line = getsh2cc2key(sh2cc2, 'wfthres')
    if line[0]:
        QMin['wfthres'] = float(line[1])
    if 'overlap' in QMin:
        # QMin['wfoverlap']=get_sh2cc2_environ(sh2cc2,'wfoverlap')
        QMin['wfoverlap'] = get_sh2cc2_environ(sh2cc2, 'wfoverlap', False, False)
        if QMin['wfoverlap'] is None:
            ciopath = os.path.join(os.path.expandvars(os.path.expanduser('$SHARC')), 'wfoverlap.x')
            if os.path.isfile(ciopath):
                QMin['wfoverlap'] = ciopath
            else:
                print('Give path to wfoverlap.x in RICC2.resources!')
        line = getsh2cc2key(sh2cc2, 'numfrozcore')
        if line[0]:
            numfroz = int(line[1])
            if numfroz == 0:
                QMin['ncore'] = 0
            elif numfroz > 0:
                QMin['ncore'] = numfroz
            elif numfroz < 0:
                pass        # here we rely on the frozen key from the template below


    # TheoDORE settings
    if 'theodore' in QMin:
        QMin['theodir'] = get_sh2cc2_environ(sh2cc2, 'theodir', False, False)
        if QMin['theodir'] is None or not os.path.isdir(QMin['theodir']):
            print('Give path to the TheoDORE installation directory in TURBOMOLE.resources!')
            sys.exit(71)
        os.environ['THEODIR'] = QMin['theodir']
        if 'PYTHONPATH' in os.environ:
            os.environ['PYTHONPATH'] += os.pathsep + os.path.join(QMin['theodir'], 'lib') + os.pathsep + QMin['theodir']
        else:
            os.environ['PYTHONPATH'] = os.path.join(QMin['theodir'], 'lib') + os.pathsep + QMin['theodir']


    # norestart setting
    line = getsh2cc2key(sh2cc2, 'no_ricc2_restart')
    if line[0]:
        QMin['no_ricc2_restart'] = []


    # neglected gradients
    QMin['neglected_gradient'] = 'zero'
    if 'grad' in QMin:
        line = getsh2cc2key(sh2cc2, 'neglected_gradient')
        if line[0]:
            if line[1].lower().strip() == 'zero':
                QMin['neglected_gradient'] = 'zero'
            elif line[1].lower().strip() == 'gs':
                QMin['neglected_gradient'] = 'gs'
            elif line[1].lower().strip() == 'closest':
                QMin['neglected_gradient'] = 'closest'
            else:
                print('Unknown argument to "neglected_gradient"!')
                sys.exit(72)





    # --------------------------------------------- Template ----------------------------------

    # open template
    template = readfile('RICC2.template')

    QMin['template'] = {}
    integers = ['frozen', 'charge']
    strings = ['basis', 'auxbasis', 'method', 'scf', 'spin-scaling', 'basislib']
    floats = []
    booleans = ['douglas-kroll']
    for i in booleans:
        QMin['template'][i] = False
    QMin['template']['method'] = 'adc(2)'
    QMin['template']['scf'] = 'dscf'
    QMin['template']['spin-scaling'] = 'none'
    QMin['template']['basislib'] = ''
    QMin['template']['charge'] = 0
    QMin['template']['frozen'] = -1

    QMin['template']['theodore_prop'] = ['Om', 'PRNTO', 'S_HE', 'Z_HE', 'RMSeh']
    QMin['template']['theodore_fragment'] = []

    for line in template:
        line = re.sub('#.*$', '', line).lower().split()
        if len(line) == 0:
            continue
        elif line[0] in integers:
            QMin['template'][line[0]] = int(line[1])
        elif line[0] in booleans:
            QMin['template'][line[0]] = True
        elif line[0] in strings:
            QMin['template'][line[0]] = line[1]
        elif line[0] in floats:
            QMin['template'][line[0]] = float(line[1])

    necessary = ['basis']
    for i in necessary:
        if i not in QMin['template']:
            print('Key %s missing in template file!' % (i))
            sys.exit(73)

    # make basis set name in correct case, so that Turbomole recognizes them
    for basis in BASISSETS:
        if QMin['template']['basis'].lower() == basis.lower():
            QMin['template']['basis'] = basis
            break
    if 'auxbasis' in QMin['template']:
        for basis in BASISSETS:
            if QMin['template']['auxbasis'].lower() == basis.lower():
                QMin['template']['auxbasis'] = basis
                break
        if QMin['template']['basislib']:
            print('Keywords "basislib" and "auxbasis" cannot be used together in template!\nInstead, create a file for the auxbasis in /basislib/cbasen/')
            sys.exit(74)

    # go through sh2cc2 for the theodore settings
    for line in sh2cc2:
        orig = re.sub('#.*$', '', line).strip()
        line = orig.lower().split()
        if len(line) == 0:
            continue

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




# --------------------------------------------- QM/MM ----------------------------------

    # qmmm keyword
    QMin['qmmm'] = False
    QMin['template']['qmmm'] = False
    i = 0
    for line in template:
        line = re.sub('#.*$', '', line).lower().split()
        if len(line) < 1:
            continue
        if line[0] == 'qmmm':
            QMin['qmmm'] = True

    QMin['cobramm'] = False
    QMin['template']['cobramm'] = False
    i = 0
    for line in template:
        line = re.sub('#.*$', '', line).lower().split()
        if len(line) < 1:
            continue
        if line[0] == 'cobramm':
            QMin['cobramm'] = True
    if QMin['cobramm']:
        QMin['template']['cobramm'] = True

    # prepare everything
    if QMin['qmmm']:
        QMin['template']['qmmm'] = True

        # get settings from RICC2.resources
        # Tinker
        line = getsh2cc2key(sh2cc2, 'tinker')
        if not line[0]:
            print('TINKER path not given!')
            sys.exit(75)
        line = os.path.expandvars(line[1].strip())
        line = os.path.expanduser(line)
        line = os.path.abspath(line)
        QMin['tinker'] = line
        if not os.path.isfile(os.path.join(QMin['tinker'], 'bin', 'tkr2qm_s')):
            print('TINKER executable at "%s" not found!' % os.path.join(QMin['tinker'], 'bin', 'tkr2qm_s'))
            sys.exit(76)

        # table and ff files
        for line in sh2cc2:
            orig = re.sub('#.*$', '', line).strip()
            line = orig.lower().split()
            if len(line) == 0:
                continue
            elif line[0] == 'qmmm_table':
                line2 = orig.split(None, 1)
                if len(line2) < 2:
                    print('Please specify a connection table file after "qmmm_table"!')
                    sys.exit(77)
                filename = os.path.abspath(os.path.expandvars(os.path.expanduser(line2[1])))
                QMin['template']['qmmm_table'] = filename
            elif line[0] == 'qmmm_ff_file':
                line2 = orig.split(None, 1)
                if len(line2) < 2:
                    print('Please specify a force field file after "qmmm_ff_file"!')
                    sys.exit(78)
                filename = os.path.abspath(os.path.expandvars(os.path.expanduser(line2[1])))
                QMin['template']['qmmm_ff_file'] = filename

        # prepare data structures and run Tinker
        QMin['qmmm'] = prepare_QMMM(QMin, QMin['template']['qmmm_table'])
        QMMMout = execute_tinker(QMin, QMin['template']['qmmm_ff_file'])

        # modify QMin dict
        QMin['geo_orig'] = QMin['geo']
        QMin['geo'] = QMin['qmmm']['QM_coords']
        QMin['natom_orig'] = QMin['natom']
        QMin['natom'] = len(QMin['geo'])
        QMin['pointcharges'] = deepcopy(QMin['qmmm']['pointcharges'])
        # QMin['pointcharges']=[]
        # for iatom in range(QMin['qmmm']['natom_table']):
        # atom=QMin['qmmm']['MM_coords'][iatom]
        # for iatom,atom in enumerate(QMin['qmmm']['MM_coords']):
        # QMin['pointcharges'].append( [atom[1],atom[2],atom[3],QMin['qmmm']['MMpc'][iatom]] )



# --------------------------------------------- logic checks ----------------------------------





    # logic checks:

    # find method
    allowed_methods = ['cc2', 'adc(2)']
    for m in allowed_methods:
        if QMin['template']['method'] == m:
            QMin['method'] = m
            break
    else:
        print('Unknown method "%s" given in RICC2.template' % (QMin['template']['method']))
        sys.exit(79)

    # find spin-scaling
    allowed_methods = ['scs', 'sos', 'lt-sos', 'none']
    if not any([QMin['template']['spin-scaling'] == i for i in allowed_methods]):
        print('Unknown spin-scaling "%s" given in RICC2.template' % (QMin['template']['spin-scaling']))
        sys.exit(80)

    # find SCF program
    allowed_methods = ['dscf', 'ridft']
    if not any([QMin['template']['scf'] == i for i in allowed_methods]):
        print('Unknown SCF program "%s" given in RICC2.template' % (QMin['template']['scf']))
        sys.exit(81)

    # get number of electrons
    nelec = 0
    for atom in QMin['geo']:
        nelec += NUMBERS[atom[0].title()]
    nelec -= QMin['template']['charge']
    QMin['nelec'] = nelec
    if nelec % 2 != 0:
        print('Currently, only even-electronic systems are possible in SHARC_RICC2.py!')
        sys.exit(82)

    # no soc for elements beyond Kr due to ECP
    if 'soc' in QMin and any([NUMBERS[atom[0].title()] > 36 for atom in QMin['geo']]):
        print('Spin-orbit couplings for elements beyond Kr do not work due to default ECP usage!')
        sys.exit(83)

    # soc and cc2 do not work together
    if 'soc' in QMin and 'cc2' in QMin['template']['method']:
        print('Currently, spin-orbit coupling is not possible at CC2 level. Please use ADC(2)!')
        sys.exit(84)

    # lt-sos-CC2/ADC(2) does not work in certain cases
    if 'lt-sos' in QMin['template']['spin-scaling']:
        if QMin['ncpu'] > 1:
            print('NOTE: Laplace-transformed SOS-%s is not fully SMP parallelized.' % (QMin['template']['method'].upper()))
            # sys.exit(85)
        if 'soc' in QMin:
            print('Laplace-transformed SOS-%s is not compatible with SOC calculation!' % (QMin['template']['method'].upper()))
            sys.exit(86)
        if QMin['dipolelevel'] == 2:
            print('Laplace-transformed SOS-%s is not compatible with dipolelevel=2!' % (QMin['template']['method'].upper()))
            sys.exit(87)

    # number of properties/entries calculated by TheoDORE
    if 'theodore' in QMin:
        QMin['template']['theodore_n'] = len(QMin['template']['theodore_prop']) + len(QMin['template']['theodore_fragment'])**2
    else:
        QMin['template']['theodore_n'] = 0


    # Check the save directory
    try:
        ls = os.listdir(QMin['savedir'])
        err = 0
    except OSError:
        err = 1
    if 'init' in QMin:
        err = 0
    elif 'overlap' in QMin:
        if 'newstep' in QMin:
            if 'mos' not in ls:
                print('File "mos" missing in SAVEDIR!')
                err += 1
            if 'coord' not in ls:
                print('File "coord" missing in SAVEDIR!')
                err += 1
            for imult, nstates in enumerate(QMin['states']):
                if nstates < 1:
                    continue
                if not 'dets.%i' % (imult + 1) in ls:
                    print('File "dets.%i.old" missing in SAVEDIR!' % (imult + 1))
                    err += 1
        elif 'samestep' in QMin or 'restart' in QMin:
            if 'mos.old' not in ls:
                print('File "mos" missing in SAVEDIR!')
                err += 1
            if 'coord.old' not in ls:
                print('File "coord.old" missing in SAVEDIR!')
                err += 1
            for imult, nstates in enumerate(QMin['states']):
                if nstates < 1:
                    continue
                if not 'dets.%i.old' % (imult + 1) in ls:
                    print('File "dets.%i.old" missing in SAVEDIR!' % (imult + 1))
                    err += 1
    if err > 0:
        print('%i files missing in SAVEDIR=%s' % (err, QMin['savedir']))
        sys.exit(88)

    if PRINT:
        printQMin(QMin)

    return QMin

# ======================================================================= #


# =============================================================================================== #
# =============================================================================================== #
# =========================================== gettasks and setup routines ======================= #
# =============================================================================================== #
# =============================================================================================== #

def get_jobs(QMin):
    # returns a list with the setup for the different ricc2 jobs
    # first, find the properties we need to calculate
    prop = set()
    prop.add('E')
    if 'soc' in QMin:
        prop.add('tmexc_soc')
    if 'grad' in QMin:
        for i in QMin['gradmap']:
            if i == (1, 1):
                prop.add(tuple(['gsgrad'] + list(i)))
            else:
                prop.add(tuple(['exgrad'] + list(i)))
    if 'dm' in QMin:
        if QMin['dipolelevel'] >= 0:
            # make only dipole moments which are for free
            if 'soc' in QMin:
                prop.add('tmexc_dm')
        if QMin['dipolelevel'] >= 1:
            # <S0|dm|T> only works for ADC(2)
            if QMin['states'][0] >= 1 and (QMin['template']['method'] == 'adc(2)' or len(QMin['states']) == 1):
                prop.add('spectrum')
        if QMin['dipolelevel'] >= 2:
            prop.add('exprop_dm')
            prop.add('static_dm')
            # tmexc does not work for CC2 if excited singlets and triples are present
            if not (QMin['template']['method'] == 'cc2' and len(QMin['states']) > 1 and QMin['states'][0] > 1):
                prop.add('tmexc_dm')

    # define the rules for distribution of jobs
    forbidden = {'E': [],
                 'tmexc_soc': ['gsgrad', 'exgrad', 'exprop_dm', 'static_dm'],
                 'tmexc_dm': [],
                 'spectrum': [],
                 'gsgrad': ['tmexc_soc'],
                 'exgrad': ['tmexc_soc', 'exgrad'],
                 'static_dm': ['tmexc_soc'],
                 'exprop_dm': ['tmexc_soc']
                 }
    if QMin['qmmm']:
        forbidden['gsgrad'].append('exgrad')
        forbidden['exgrad'].append('gsgrad')
    if QMin['cobramm']:  # 21.09.20#
        forbidden['gsgrad'].append('exgrad')
        forbidden['exgrad'].append('gsgrad')

    priority = ['E',
                'tmexc_soc',
                'tmexc_dm',
                'spectrum',
                'gsgrad',
                'exgrad',
                'static_dm',
                'exprop_dm']

    # print prop

    # second, distribute the properties into jobs
    jobs = []
    # iterate until prop is empty
    while True:
        job = set()
        # process according to priority order
        for prior in priority:
            # print('Prior:',prior)
            # print('Job before iter: ',job)
            # print('Prop pefore iter:',prop)
            # cannot delete from prop during iteration, therefore make copy
            prop2 = deepcopy(prop)
            # check if prior is in prop
            for p in prop:
                if prior in p:
                    # check if any forbidden task already in job
                    fine = True
                    for forb in forbidden[prior]:
                        for j in job:
                            if forb in j:
                                fine = False
                    # if allowed, put task into job, delete from prop
                    if fine:
                        job.add(p)
                        prop2.remove(p)
            # copy back prop after iteration
            prop = deepcopy(prop2)
            # print('Job after iter: ',job)
            # print('Prop after iter:',prop)
            # print
        jobs.append(job)
        if len(prop) == 0:
            break

    # pprint.pprint(jobs)

    return jobs


# ======================================================================= #

def gettasks(QMin):
    ''''''

    states = QMin['states']
    nstates = QMin['nstates']
    nmstates = QMin['nmstates']

    # Currently implemented keywords: h, soc, dm, grad, overlap, samestep, restart, init
    tasks = []
    # During initialization, create all temporary directories
    # and link them appropriately
    tasks.append(['mkdir', QMin['scratchdir']])
    tasks.append(['link', QMin['scratchdir'], QMin['pwd'] + '/SCRATCH', False])
    tasks.append(['mkdir', QMin['scratchdir'] + '/JOB'])
    if 'overlap' in QMin:
        tasks.append(['mkdir', QMin['scratchdir'] + '/OVERLAP'])
        tasks.append(['mkdir', QMin['scratchdir'] + '/AO_OVL'])

    if 'init' in QMin:
        tasks.append(['mkdir', QMin['savedir']])
        tasks.append(['link', QMin['savedir'], QMin['pwd'] + '/SAVE', False])

    if 'molden' in QMin and not os.path.isdir(QMin['savedir'] + '/MOLDEN'):
        tasks.append(['mkdir', QMin['savedir'] + '/MOLDEN'])

    if 'samestep' not in QMin and 'init' not in QMin and 'restart' not in QMin:
        tasks.append(['movetoold'])

    if 'backup' in QMin:
        tasks.append(['mkdir', QMin['savedir'] + '/backup/'])

    # do all TURBOMOLE calculations
    if 'overlaponly' not in QMin:

        tasks.append(['cleanup', QMin['scratchdir'] + '/JOB'])
        tasks.append(['writegeom', QMin['scratchdir'] + '/JOB'])
        tasks.append(['define', QMin['scratchdir'] + '/JOB'])
        tasks.append(['modify_control'])

        # initial orbitals
        if 'always_guess' in QMin:
            # no further action, define creates MOs
            pass
        elif 'always_orb_init' in QMin:
            # always take "mos.init" from the main directory
            tryfile = os.path.join(QMin['pwd'], 'mos.init')
            if os.path.isfile(tryfile):
                tasks.append(['getmo', tryfile])
        else:
            if 'init' in QMin:
                tryfile = os.path.join(QMin['pwd'], 'mos.init')
                if os.path.isfile(tryfile):
                    tasks.append(['getmo', tryfile])
            elif 'samestep' in QMin:
                tryfile = os.path.join(QMin['savedir'], 'mos')
                if os.path.isfile(tryfile):
                    tasks.append(['getmo', tryfile])
            elif 'restart' in QMin:
                tryfile = os.path.join(QMin['savedir'], 'mos.old')
                if os.path.isfile(tryfile):
                    tasks.append
            elif 'newstep' in QMin:
                tryfile = os.path.join(QMin['savedir'], 'mos.old')
                if os.path.isfile(tryfile):
                    tasks.append(['getmo', tryfile])

        # SCF calculation
        if QMin['template']['scf'] == 'ridft':
            tasks.append(['ridft'])
        tasks.append(['dscf'])
        if 'molden' in QMin or 'theodore' in QMin:
            tasks.append(['copymolden'])

        # Orca call for SOCs
        if 'soc' in QMin:
            tasks.append(['orca_soc'])

        # ricc2 calls:
        # prop={'tmexc_soc','tmexc_dm', 'spectrum','exprop_dm','static_dm','gs_grad','exgrad'}
        jobs = get_jobs(QMin)
        for ijob, job in enumerate(jobs):
            tasks.append(['prep_control', job])
            tasks.append(['ricc2'])
            tasks.append(['get_RICC2out', job])
            if ijob == 0:
                tasks.append(['save_data'])
                if 'nooverlap' not in QMin:
                    mults = [i + 1 for i in range(len(QMin['states'])) if QMin['states'][i] > 0]
                    tasks.append(['get_dets', QMin['scratchdir'] + '/JOB', mults])
                # if 'molden' in QMin:
                    # tasks.append(['molden'])

                if 'theodore' in QMin:
                    tasks.append(['run_theodore'])
                    tasks.append(['get_theodore'])
                    if 'molden' in QMin:
                        tasks.append(['copy_ntos'])

        if 'overlap' in QMin:
            # get mixed AO overlap
            tasks.append(['cleanup', QMin['scratchdir'] + '/AO_OVL'])
            tasks.append(['get_AO_OVL', QMin['scratchdir'] + '/AO_OVL'])

            for imult in range(len(QMin['states'])):
                if QMin['states'][imult] == 0:
                    continue
                tasks.append(['cleanup', QMin['scratchdir'] + '/OVERLAP'])
                tasks.append(['wfoverlap', QMin['scratchdir'] + '/OVERLAP', imult + 1])
                tasks.append(['get_wfovlout', QMin['scratchdir'] + '/OVERLAP', imult + 1])

    if 'backup' in QMin:
        tasks.append(['backupdata', QMin['backup']])

    if 'cleanup' in QMin:
        tasks.append(['cleanup', QMin['savedir']])
    # if not DEBUG:
        # tasks.append(['cleanup',QMin['scratchdir']])

    return tasks


# =============================================================================================== #
# =============================================================================================== #
# =========================================== SUBROUTINES TO RUNEVERYTING ======================= #
# =============================================================================================== #
# =============================================================================================== #

def mkdir(DIR):
    # mkdir the DIR, or clean it if it exists
    if os.path.exists(DIR):
        if os.path.isfile(DIR):
            print('%s exists and is a file!' % (DIR))
            sys.exit(89)
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
            sys.exit(90)

# ======================================================================= #


def link(PATH, NAME, crucial=True, force=True):
    # do not create broken links
    if not os.path.exists(PATH):
        print('Source %s does not exist, cannot create link!' % (PATH))
        sys.exit(91)
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
                    sys.exit(92)
                else:
                    return
    elif os.path.exists(NAME):
        # NAME is not a link. The interface will not overwrite files/directories with links, even with force=True
        print('%s exists, cannot create a link of the same name!' % (NAME))
        if crucial:
            sys.exit(93)
        else:
            return
    os.symlink(PATH, NAME)

# ======================================================================= #


def shorten_DIR(string):
    maxlen = 40
    front = 12
    if len(string) > maxlen:
        return string[0:front] + '...' + string[-(maxlen - 3 - front):]
    else:
        return string + ' ' * (maxlen - len(string))

# ======================================================================= #


def cleandir(directory):
    if DEBUG:
        print('===> Cleaning up directory %s\n' % (directory))
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
    if DEBUG:
        print('\n')

# ======================================================================= #


def movetoold(QMin):
    # rename all files in savedir
    saveable = ['dets', 'mos', 'coord']
    savedir = QMin['savedir']
    ls = os.listdir(savedir)
    if ls == []:
        return
    for f in ls:
        f2 = savedir + '/' + f
        if os.path.isfile(f2):
            if any([i in f for i in saveable]):
                if 'old' not in f:
                    fdest = f2 + '.old'
                    shutil.move(f2, fdest)

# ======================================================================= #


def save_data(QMin):
    # copy files to savedir
    saveable = ['mos', 'coord']
    for i in saveable:
        fromfile = os.path.join(QMin['scratchdir'], 'JOB', i)
        tofile = os.path.join(QMin['savedir'], i)
        shutil.copy(fromfile, tofile)

# ======================================================================= #


def getmo(mofile, QMin):
    if os.path.exists(mofile):
        tofile = os.path.join(QMin['scratchdir'], 'JOB', 'mos')
        shutil.copy(mofile, tofile)
    else:
        print('Could not find mocoef-file %s!' % (mofile))
        sys.exit(94)

# ======================================================================= #


def writegeom(QMin):
    factor = au2a
    fname = QMin['scratchdir'] + '/JOB/geom.xyz'
    string = '%i\n\n' % (QMin['natom'])
    for atom in QMin['geo']:
        string += atom[0]
        for xyz in range(1, 4):
            string += '  %f' % (atom[xyz] * factor)
        string += '\n'
    writefile(fname, string)

    os.chdir(QMin['scratchdir'] + '/JOB')
    error = sp.call('x2t geom.xyz > coord', shell=True)
    if error != 0:
        print('xyz2col call failed!')
        sys.exit(95)
    os.chdir(QMin['pwd'])

    # QM/MM
    if QMin['qmmm']:
        string = '$point_charges nocheck\n'
        for atom in QMin['pointcharges']:
            string += '%16.12f %16.12f %16.12f %12.9f\n' % (atom[0] / au2a, atom[1] / au2a, atom[2] / au2a, atom[3])
        string += '$end\n'
        filename = QMin['scratchdir'] + '/JOB/pc'
        writefile(filename, string)

    # COBRAMM
    if QMin['cobramm']:
        # chargefiles='charge.dat'
        # tocharge=os.path.join(QMin['scratchdir']+'/JOB/point_charges')
        # shutil.copy(chargefiles,tocharge)
        cobcharges = open('charge.dat', 'r')
        charges = cobcharges.read()
        only_atom = charges.split()
        only_atom.pop(0)
        filename = QMin['scratchdir'] + '/JOB/point_charges'
        string = '$point_charges nocheck\n'
        string += charges
        # counter=0
        # for atom in only_atom:
        #   	string+=atom
        #    string+=' '
        #    counter+=1
        #    if counter == 4:
        #      string+='\n'
        #      counter=0
        #    #string+='\n'
        string += '$end'
        writefile(filename, string)

# ======================================================================= #


def runProgram(string, workdir, outfile, errfile=''):
    prevdir = os.getcwd()
    if DEBUG:
        print(workdir)
    os.chdir(workdir)
    if PRINT or DEBUG:
        starttime = datetime.datetime.now()
        sys.stdout.write('%s\n\t%s' % (string, starttime))
        sys.stdout.flush()
    stdoutfile = open(os.path.join(workdir, outfile), 'w')
    if errfile:
        stderrfile = open(os.path.join(workdir, errfile), 'w')
    else:
        stderrfile = sp.STDOUT
    try:
        runerror = sp.call(string, shell=True, stdout=stdoutfile, stderr=stderrfile)
    except OSError:
        print('Call have had some serious problems:', OSError)
        sys.exit(96)
    stdoutfile.close()
    if errfile:
        stderrfile.close()
    if PRINT or DEBUG:
        endtime = datetime.datetime.now()
        sys.stdout.write('\t%s\t\tRuntime: %s\t\tError Code: %i\n\n' % (endtime, endtime - starttime, runerror))
    os.chdir(prevdir)
    return runerror

# ======================================================================= #


def define(path, QMin, ricc2=True):

    # first three sections
    if QMin['template']['basislib']:
        # write definrc
        string = '''basis=%s/basen
basis=%s/cbasen
''' % (QMin['template']['basislib'], QMin['template']['basislib'])
        infile = os.path.join(path, '.definerc')
        writefile(infile, string)

    # write define input
    string = '''
title: SHARC-RICC2 run
a coord
*
no
'''
    if QMin['template']['basislib']:
        string += '''lib
3
'''
    string += '''b
all %s
''' % (QMin['template']['basis'])
    if not ricc2:
        string += 'c\nall 0.\n'
    string += '''*
eht
y
%i
y
''' % (
       QMin['template']['charge']
       )

    if ricc2:
        # cc section
        string += 'cc\n'

        # frozen orbitals
        if QMin['template']['frozen'] == 0:
            pass
        elif QMin['template']['frozen'] < 0:
            string += 'freeze\n*\n'
        elif QMin['template']['frozen'] > 0:
            string += 'freeze\ncore %i\n*\n' % (QMin['template']['frozen'])

        # auxilliary basis set: cbas
        # this is mandatory for ricc2, so if the user does not give an auxbasis, we take the default one
        if 'auxbasis' not in QMin['template']:
            if QMin['template']['basislib']:
                string += 'cbas\n'
                # skip error messages in define
                elements = set()
                for atom in QMin['geo']:
                    elements.add(atom[0])
                string += '\n\n' * len(elements)
                string += '''lib
4
b
all %s
*
''' % (QMin['template']['basis'])
            else:
                string += 'cbas\n*\n'
        else:
            string += 'cbas\nb\nall %s\n*\n' % (QMin['template']['auxbasis'])

        # memory
        string += 'memory %i\n' % (QMin['memory'])

        # ricc2 section (here we set only the model)
        string += 'ricc2\n%s\n' % (QMin['template']['method'])
        if QMin['template']['spin-scaling'] == 'none':
            pass
        elif QMin['template']['spin-scaling'] == 'scs':
            string += 'scs\n'
        elif QMin['template']['spin-scaling'] == 'sos':
            string += 'sos\n'
        elif QMin['template']['spin-scaling'] == 'lt-sos':
            string += 'sos\n'
        # number of DIIS vectors
        ndiis = max(10, 5 * max(QMin['states']))
        string += 'mxdiis = %i\n' % (ndiis)
        string += '*\n'

        # leave cc input
        string += '*\n'

    # leave define
    string += '*\n'

    # string contains the input for define, now call it
    infile = os.path.join(path, 'define.input')
    writefile(infile, string)
    string = 'define < define.input'
    runerror = runProgram(string, path, 'define.output')

    if runerror != 0:
        print('define call failed! Error Code=%i Path=%s' % (runerror, path))
        sys.exit(97)


    return

# ======================================================================= #


def add_section_to_control(path, section):
    # adds a section keyword to the control file, before $end
    # if section does not start with $, $ will be prepended
    if not section[0] == '$':
        section = '$' + section
    infile = readfile(path)

    # get iline of $end
    iline = -1
    while True:
        iline += 1
        line = infile[iline]
        if '$end' in line:
            break

    outfile = infile[:iline]
    outfile.append(section + '\n')
    outfile.extend(infile[iline:])
    writefile(path, outfile)
    return

# ======================================================================= #


def add_option_to_control_section(path, section, newline):
    if not section[0] == '$':
        section = '$' + section
    infile = readfile(path)
    newline = '  ' + newline

    # get iline where section starts
    iline = -1
    while True:
        iline += 1
        if iline == len(infile):
            return
        line = infile[iline]
        if section in line:
            break

    # get jline where section ends
    jline = iline
    while True:
        jline += 1
        line = infile[jline]
        if '$' in line:
            break
        # do nothing if line is already there
        if newline + '\n' == line:
            return

    # splice together new file
    outfile = infile[:jline]
    outfile.append(newline + '\n')
    outfile.extend(infile[jline:])
    writefile(path, outfile)
    return

# ======================================================================= #


def remove_section_in_control(path, section):
    # removes a keyword and its options from control file
    if not section[0] == '$':
        section = '$' + section
    infile = readfile(path)

    # get iline where section starts
    iline = -1
    while True:
        iline += 1
        if iline == len(infile):
            return
        line = infile[iline]
        if section in line:
            break

    # get jline where section ends
    jline = iline
    while True:
        jline += 1
        line = infile[jline]
        if '$' in line:
            break

    # splice together new file
    outfile = infile[:iline] + infile[jline:]
    writefile(path, outfile)
    return

# ======================================================================= #


def modify_control(QMin):
    # this adjusts the control file for the main JOB calculations
    control = os.path.join(QMin['scratchdir'], 'JOB/control')

    if 'soc' in QMin:
        add_section_to_control(control, '$mkl')
    if QMin['template']['douglas-kroll']:
        add_section_to_control(control, '$rdkh')

    # add laplace keyword for LT-SOS
    if QMin['template']['spin-scaling'] == 'lt-sos':
        add_section_to_control(control, '$laplace')
        add_option_to_control_section(control, '$laplace', 'conv=5')

    # remove_section_in_control(control,'$optimize')
    # add_option_to_control_section(control,'$ricc2','scs')
    remove_section_in_control(control, '$scfiterlimit')
    add_section_to_control(control, '$scfiterlimit 100')

    # QM/MM point charges
    if QMin['qmmm']:
        add_option_to_control_section(control, '$drvopt', 'point charges')
        add_section_to_control(control, '$point_charges file=pc')
        add_section_to_control(control, '$point_charge_gradients file=pc_grad')
    return

    # COBRAMM
    if QMin['cobramm']:
        add_option_to_control_section(control, '$drvopt', 'point charges')
        add_section_to_control(control, '$point_charges file=point_charges')  # inserire nome file quando deciso
        add_section_to_control(control, '$point_charge_gradients file=pc_grad')
    return


# ======================================================================= #
def prep_control(QMin, job):
    # prepares the control file to calculate grad, soc, dm
    # job contains: 'tmexc_soc','tmexc_dm', 'spectrum','exprop_dm','static_dm','gsgrad','exgrad', 'E'

    control = os.path.join(QMin['scratchdir'], 'JOB/control')

    # remove sections to cleanly rewrite them
    remove_section_in_control(control, '$response')
    remove_section_in_control(control, '$excitations')

    # add number of states
    add_section_to_control(control, '$excitations')
    add_option_to_control_section(control, '$ricc2', 'maxiter 45')
    nst = QMin['states'][0] - 1       # exclude ground state here
    if nst >= 1:
        string = 'irrep=a multiplicity=1 nexc=%i npre=%i nstart=%i' % (nst, nst + 1, nst + 1)
        add_option_to_control_section(control, '$excitations', string)
    if len(QMin['states']) >= 3:
        nst = QMin['states'][2]
        if nst >= 1:
            string = 'irrep=a multiplicity=3 nexc=%i npre=%i nstart=%i' % (nst, nst + 1, nst + 1)
            add_option_to_control_section(control, '$excitations', string)

    # add response section
    if 'static_dm' or 'gsgrad' in job:
        add_section_to_control(control, '$response')

    # add property lines
    if 'tmexc_soc' in job or 'tmexc_dm' in job:
        string = 'tmexc istates=all fstates=all operators='
        prop = []
        if 'tmexc_soc' in job:
            prop.append('soc')
        if 'tmexc_dm' in job:
            prop.append('diplen')
        string += ','.join(prop)
        add_option_to_control_section(control, '$excitations', string)
    if 'spectrum' in job:
        string = 'spectrum states=all operators=diplen'
        add_option_to_control_section(control, '$excitations', string)
    if 'exprop_dm' in job:
        string = 'exprop states=all relaxed operators=diplen'
        add_option_to_control_section(control, '$excitations', string)
    if 'static_dm' in job:
        string = 'static relaxed operators=diplen'
        add_option_to_control_section(control, '$response', string)

    if QMin['cobramm']:
        add_option_to_control_section(control, '$drvopt', ' point charges')
        add_section_to_control(control, '$point_charges file=point_charges')  # inserire nome file quando deciso
        add_section_to_control(control, '$point_charge_gradients file=pc_grad')

    # add gradients
    for j in job:
        if 'gsgrad' in j:
            string = 'gradient'
            add_option_to_control_section(control, '$response', string)
        if 'exgrad' in j:
            string = 'xgrad states=(a{%i} %i)' % (j[1], j[2] - (j[1] == 1))
            add_option_to_control_section(control, '$excitations', string)


    # ricc2 restart
    if 'E' not in job and 'no_ricc2_restart' not in QMin:
        if QMin['tmversion'] >= 7.7:
            add_option_to_control_section(control, '$ricc2', 'restart on')
        else:
            add_option_to_control_section(control, '$ricc2', 'restart')
        restartfile = os.path.join(QMin['scratchdir'], 'JOB/restart.cc')
        try:
            os.remove(restartfile)
        except OSError:
            pass
    else:
        if QMin['tmversion'] >= 7.7:
            add_option_to_control_section(control, '$ricc2', 'restart off')
        else:
            add_option_to_control_section(control, '$ricc2', 'norestart')

    # D1 and D2 diagnostic
    if DEBUG and 'E' in job:
        add_option_to_control_section(control, '$ricc2', 'd1diag')
        add_section_to_control(control, '$D2-diagnostic')

    return

# ======================================================================= #


def get_dets(path, mults, QMin):
    # read all determinant expansions from working directory and put them into the savedir

    for imult in mults:
        ca = civfl_ana(path, imult, maxsqnorm=QMin['wfthres'], filestr='CCRE0')
        for istate in range(1, 1 + QMin['states'][imult - 1]):
            ca.get_state_dets(istate)
        writename = os.path.join(QMin['savedir'], 'dets.%i' % (imult))
        ca.write_det_file(QMin['states'][imult - 1], wname=writename)

        # for CC2, also save the left eigenvectors
        if QMin['template']['method'] == 'cc2':
            ca = civfl_ana(path, imult, maxsqnorm=QMin['wfthres'], filestr='CCLE0')
            for istate in range(1, 1 + QMin['states'][imult - 1]):
                ca.get_state_dets(istate)
            writename = os.path.join(QMin['savedir'], 'dets_left.%i' % (imult))
            ca.write_det_file(QMin['states'][imult - 1], wname=writename)

        if 'frozenmap' not in QMin:
            QMin['frozenmap'] = {}
        QMin['frozenmap'][imult] = ca.nfrz
    return QMin

# ======================================================================= #


def get_AO_OVL(path, QMin):
    # get double geometry
    oldgeom = readfile(os.path.join(QMin['savedir'], 'coord.old'))
    newgeom = readfile(os.path.join(QMin['savedir'], 'coord'))
    string = '$coord\n'
    wrt = False
    for line in oldgeom:
        if '$coord' in line:
            wrt = True
            continue
        elif '$' in line:
            wrt = False
            continue
        if wrt:
            string += line
    for line in newgeom:
        if '$coord' in line:
            wrt = True
            continue
        elif '$' in line:
            wrt = False
            continue
        if wrt:
            string += line
    string += '$end\n'
    tofile = os.path.join(path, 'coord')
    writefile(tofile, string)

    # call define and then add commands to control file
    define(path, QMin, ricc2=False)
    controlfile = os.path.join(path, 'control')
    remove_section_in_control(controlfile, '$scfiterlimit')
    add_section_to_control(controlfile, '$scfiterlimit 0')
    if QMin['tmversion'] >= 7.7:
        add_section_to_control(controlfile, '$intsdebug 1 sao')
    else:
        add_section_to_control(controlfile, '$intsdebug sao')
    add_section_to_control(controlfile, '$closed shells')
    add_option_to_control_section(controlfile, '$closed shells', 'a 1-2')
    add_section_to_control(controlfile, '$scfmo none')

    # write geometry again because define tries to be too clever with the double geometry
    tofile = os.path.join(path, 'coord')
    writefile(tofile, string)

    # call dscf
    string = 'dscf'
    runerror = runProgram(string, path, 'dscf.out')

    # get AO overlap matrix from dscf.out
    dscf = readfile(os.path.join(path, 'dscf.out'))
    for line in dscf:
        if ' number of basis functions   :' in line:
            nbas = int(line.split()[-1])
            break
    else:
        print('Could not find number of basis functions in dscf.out!')
        sys.exit(98)
    iline = -1
    while True:
        iline += 1
        line = dscf[iline]
        if 'OVERLAP(SAO)' in line:
            break
    iline += 1
    ao_ovl = makermatrix(nbas, nbas)
    x = 0
    y = 0
    while True:
        iline += 1
        line = dscf[iline].split()
        for el in line:
            ao_ovl[x][y] = float(el)
            ao_ovl[y][x] = float(el)
            x += 1
            if x > y:
                x = 0
                y += 1
        if y >= nbas:
            break
    # the SAO overlap in dscf output is a LOWER triangular matrix
    # hence, the off-diagonal block must be transposed

    # write AO overlap matrix to savedir
    string = '%i %i\n' % (nbas // 2, nbas // 2)
    for irow in range(nbas // 2, nbas):
        for icol in range(0, nbas // 2):
            string += '% .15e ' % (ao_ovl[icol][irow])          # note the exchanged indices => transposition
        string += '\n'
    filename = os.path.join(QMin['savedir'], 'ao_ovl')
    writefile(filename, string)

# ======================================================================= #


def wfoverlap(QMin, scradir, mult):
    # link all input files for wfoverlap
    savedir = QMin['savedir']
    link(os.path.join(savedir, 'ao_ovl'), os.path.join(scradir, 'ao_ovl'), crucial=True, force=True)
    link(os.path.join(savedir, 'mos.old'), os.path.join(scradir, 'mos.a'), crucial=True, force=True)
    link(os.path.join(savedir, 'mos'), os.path.join(scradir, 'mos.b'), crucial=True, force=True)
    if QMin['template']['method'] == 'cc2':
        link(os.path.join(savedir, 'dets_left.%i.old' % (mult)), os.path.join(scradir, 'dets.a'), crucial=True, force=True)
    else:
        link(os.path.join(savedir, 'dets.%i.old' % (mult)), os.path.join(scradir, 'dets.a'), crucial=True, force=True)
    link(os.path.join(savedir, 'dets.%i' % (mult)), os.path.join(scradir, 'dets.b'), crucial=True, force=True)

    # write input file for wfoverlap
    string = '''mix_aoovl=ao_ovl
a_mo=mos.a
b_mo=mos.b
a_det=dets.a
b_det=dets.b
a_mo_read=2
b_mo_read=2
'''
    if 'ncore' in QMin:
        icore = QMin['ncore']
    elif 'frozenmap' in QMin:
        icore = QMin['frozenmap'][mult]
    else:
        icore = 0
    string += 'ncore=%i' % (icore)
    writefile(os.path.join(scradir, 'wfovl.inp'), string)

    # run wfoverlap
    string = '%s -f wfovl.inp -m %i' % (QMin['wfoverlap'], QMin['memory'])
    runProgram(string, scradir, 'wfovl.out')

# ======================================================================= #


def run_dscf(QMin):
    workdir = os.path.join(QMin['scratchdir'], 'JOB')
    if QMin['ncpu'] > 1:
        string = 'dscf_omp'
        # TODO: check when to use OMP and when SMP, or make option
    else:
        string = 'dscf'
    runerror = runProgram(string, workdir, 'dscf.out')
    if runerror != 0:
        print('DSCF calculation crashed! Error code=%i' % (runerror))
        sys.exit(99)

    return

# ======================================================================= #


def run_ridft(QMin):
    workdir = os.path.join(QMin['scratchdir'], 'JOB')

    # add RI settings to control file
    controlfile = os.path.join(workdir, 'control')
    remove_section_in_control(controlfile, '$maxcor')
    add_section_to_control(controlfile, '$maxcor %i' % (int(QMin['memory'] * 0.6)))
    add_section_to_control(controlfile, '$ricore %i' % (int(QMin['memory'] * 0.4)))
    add_section_to_control(controlfile, '$jkbas file=auxbasis')
    add_section_to_control(controlfile, '$rij')
    add_section_to_control(controlfile, '$rik')

    if QMin['ncpu'] > 1:
        string = 'ridft_smp'
    else:
        string = 'ridft'
    runerror = runProgram(string, workdir, 'ridft.out')
    if runerror != 0:
        print('RIDFT calculation crashed! Error code=%i' % (runerror))
        sys.exit(100)

    # remove RI settings from control file
    controlfile = os.path.join(workdir, 'control')
    remove_section_in_control(controlfile, '$maxcor')
    add_section_to_control(controlfile, '$maxcor %i' % (QMin['memory']))
    remove_section_in_control(controlfile, '$rij')
    remove_section_in_control(controlfile, '$rik')

    return

# ======================================================================= #


def run_orca(QMin):
    workdir = os.path.join(QMin['scratchdir'], 'JOB')
    string = 'orca_2mkl soc -gbw'
    runerror = runProgram(string, workdir, 'orca_2mkl.out', 'orca_2mkl.err')
    if runerror != 0:
        print('orca_2mkl calculation crashed! Error code=%i' % (runerror))
        sys.exit(101)

    string = '''soc.gbw
soc.psoc
soc.soc
3
1 2 3 0 4 0 0 4
0
'''
    writefile(os.path.join(workdir, 'soc.socinp'), string)

    string = 'orca_soc soc.socinp -gbw'
    runerror = runProgram(string, workdir, 'orca_soc.out', 'orca_soc.err')
    if runerror != 0:
        print('orca_soc calculation crashed! Error code=%i' % (runerror))
        sys.exit(102)

    return


# ======================================================================= #
shift_mask = {1: (+1, +1),
              2: (-2, -1),
              3: (+1, +1),
              4: (+1, +1),
              5: (+1, +1),
              6: (+1, +1)}


def change_pre_states(workdir, itrials):
    filename = os.path.join(workdir, 'control')
    data = readfile(filename)
    data2 = copy.deepcopy(data)
    for i, line in enumerate(data):
        if 'irrep' in line:
            s = line.replace('=', ' ').split()
            mult = s[3]
            nexc = s[5]
            npre = str(int(s[7]) + shift_mask[itrials][0])
            nstart = str(int(s[9]) + shift_mask[itrials][1])
            data2[i] = 'irrep=a multiplicity=%s nexc=%s npre=%s nstart=%s\n' % (mult, nexc, npre, nstart)
        if 'maxiter 45' in line:
            data2[i] = 'maxiter 100\n'
    writefile(filename, data2)
    # TODO: check how these can be added nicely and whether they improve convergence robustly
    # string = 'thrdiis = 5'
    # add_option_to_control_section(filename, '$excitations', string)
    # string = 'thrpreopt = 6'
    # add_option_to_control_section(filename, '$excitations', string)

# ======================================================================= #


def run_ricc2(QMin):
    workdir = os.path.join(QMin['scratchdir'], 'JOB')

    # enter loop until convergence of CC2/ADC(2)
    itrials = 0
    while True:
        if QMin['ncpu'] > 1:
            string = 'ricc2_omp'
        else:
            string = 'ricc2'
        runerror = runProgram(string, workdir, 'ricc2.out')
        if runerror != 0:
            print('RICC2 calculation crashed! Error code=%i' % (runerror))
            ok = False
        # check for convergence in output file
        filename = os.path.join(workdir, 'ricc2.out')
        data = readfile(filename)
        ok = True
        for line in data:
            if 'NO CONVERGENCE' in line:
                ok = False
                break
        if ok:
            break
        # go only here if no convergence
        itrials += 1
        if itrials > max(shift_mask):
            print('Not able to obtain convergence in RICC2. Aborting...')
            sys.exit(103)
        print('No convergence of excited-state calculations! Restarting with modified number of preoptimization states...')
        change_pre_states(workdir, itrials)

    return

# ======================================================================= #


def copymolden(QMin):
    # run tm2molden in scratchdir
    string = 'molden.input\nY\n'
    filename = os.path.join(QMin['scratchdir'], 'JOB', 'tm2molden.input')
    writefile(filename, string)
    string = 'tm2molden < tm2molden.input'
    path = os.path.join(QMin['scratchdir'], 'JOB')
    runProgram(string, path, 'tm2molden.output')

    if 'molden' in QMin:
        # create directory
        moldendir = QMin['savedir'] + '/MOLDEN/'
        if not os.path.isdir(moldendir):
            mkdir(moldendir)

        # save the molden.input file
        f = QMin['scratchdir'] + '/JOB/molden.input'
        fdest = moldendir + '/step_%s.molden' % (QMin['step'][0])
        shutil.copy(f, fdest)

# ======================================================================= #


def backupdata(backupdir, QMin):
    # save all files in savedir, except which have 'old' in their name
    ls = os.listdir(QMin['savedir'])
    for f in ls:
        ff = QMin['savedir'] + '/' + f
        if os.path.isfile(ff) and 'old' not in ff:
            fdest = backupdir + '/' + f
            shutil.copy(ff, fdest)
    # save molden files
    if 'molden' in QMin:
        ff = os.path.join(QMin['savedir'], 'MOLDEN', 'step_%s.molden' % (QMin['step'][0]))
        fdest = os.path.join(backupdir, 'step_%s.molden' % (QMin['step'][0]))
        shutil.copy(ff, fdest)

# ======================================================================= #
# def runcioverlap(QMin,mult):
    # workdir=QMin['scratchdir']+'/OVERLAP'
    # os.environ['OMP_NUM_THREADS']=str(QMin['ncpu'])
    # string='%s -f cioverlap.in -m %i &> cioverlap.out' % (QMin['wfoverlap'],QMin['colmem'])
    # runerror=runProgram(string,workdir)
    # if runerror!=0:
        # print('cioverlap call not successful!')
        # sys.exit(104)

    # backup dyson.out
    # if 'backup' in QMin:
        # f=QMin['scratchdir']+'/OVERLAP/cioverlap.out'
        # fdest=QMin['backup']+'/cioverlap.out.m%i' % (mult)
        # shutil.copy(f,fdest)


# ======================================================================= #
def runeverything(tasks, QMin):

    if PRINT or DEBUG:
        print('=============> Entering RUN section <=============\n\n')

    QMout = {}
    for task in tasks:
        if DEBUG:
            print(task)
        if task[0] == 'movetoold':
            movetoold(QMin)
        if task[0] == 'mkdir':
            mkdir(task[1])
        if task[0] == 'link':
            if len(task) == 4:
                link(task[1], task[2], task[3])
            else:
                link(task[1], task[2])
        if task[0] == 'getmo':
            getmo(task[1], QMin)
        if task[0] == 'define':
            define(task[1], QMin)
        if task[0] == 'modify_control':
            modify_control(QMin)
        if task[0] == 'writegeom':
            writegeom(QMin)
        if task[0] == 'save_data':
            save_data(QMin)
        # if task[0]=='backupdata':
            # backupdata(task[1],QMin)
        if task[0] == 'copymolden':
            copymolden(QMin)
        if task[0] == 'get_dets':
            QMin = get_dets(task[1], task[2], QMin)
        if task[0] == 'cleanup':
            cleandir(task[1])
        if task[0] == 'dscf':
            run_dscf(QMin)
        if task[0] == 'ridft':
            run_ridft(QMin)
        if task[0] == 'orca_soc':
            run_orca(QMin)
        if task[0] == 'prep_control':
            prep_control(QMin, task[1])
        if task[0] == 'ricc2':
            run_ricc2(QMin)
        if task[0] == 'get_RICC2out':
            QMout = get_RICC2out(QMin, QMout, task[1])
        if task[0] == 'get_AO_OVL':
            get_AO_OVL(task[1], QMin)
        if task[0] == 'wfoverlap':
            wfoverlap(QMin, task[1], task[2])
        if task[0] == 'get_wfovlout':
            QMout = get_wfovlout(QMin, QMout, task[1], task[2])
        if task[0] == 'run_theodore':
            setupWORKDIR_TH(QMin)
            run_theodore(QMin)
        if task[0] == 'get_theodore':
            QMout = get_theodore(QMin, QMout)
        if task[0] == 'copy_ntos':
            copy_ntos(QMin)

    # if no dyson pairs were calculated because of selection rules, put an empty matrix
    if 'prop' not in QMout and 'ion' in QMin:
        QMout['prop'] = makecmatrix(QMin['nmstates'], QMin['nmstates'])

    # Phases from overlaps
    if 'phases' in QMin:
        if 'phases' not in QMout:
            QMout['phases'] = [complex(1., 0.) for i in range(QMin['nmstates'])]
        if 'overlap' in QMout:
            for i in range(QMin['nmstates']):
                if QMout['overlap'][i][i].real < 0.:
                    QMout['phases'][i] = complex(-1., 0.)

    # transform back from QM to QM/MM
    if QMin['template']['qmmm']:
        QMin, QMout = transform_QM_QMMM(QMin, QMout)


    return QMin, QMout


# =============================================================================================== #
# =============================================================================================== #
# =======================================  TheoDORE ============================================= #
# =============================================================================================== #
# =============================================================================================== #

# ======================================================================= #
def run_theodore(QMin):
    workdir = os.path.join(QMin['scratchdir'], 'JOB')
    string = 'python %s/bin/analyze_tden.py' % (QMin['theodir'])
    runerror = runProgram(string, workdir, 'theodore.out')
    if runerror != 0:
        print('Theodore calculation crashed! Error code=%i' % (runerror))
        sys.exit(105)
    return

# ======================================================================= #


def setupWORKDIR_TH(QMin):
    # mkdir the WORKDIR, or clean it if it exists, then copy all necessary files from pwd and savedir

    WORKDIR = os.path.join(QMin['scratchdir'], 'JOB')
    # write dens_ana.in
    inputstring = '''rtype='ricc2'
rfile='ricc2.out'
mo_file='molden.input'
jmol_orbitals=False
molden_orbitals=%s
read_binary=True
comp_ntos=True
alphabeta=False
Om_formula=2
eh_pop=1
print_OmFrag=True
output_file='tden_summ.txt'
prop_list=%s
at_lists=%s
''' % (('molden' in QMin),
        str(QMin['template']['theodore_prop']),
        str(QMin['template']['theodore_fragment']))

    filename = os.path.join(WORKDIR, 'dens_ana.in')
    writefile(filename, inputstring)
    return

# ======================================================================= #


def get_theodore(QMin, QMout):
    if 'theodore' not in QMout:
        QMout['theodore'] = makecmatrix(QMin['template']['theodore_n'], QMin['nmstates'])
        sumfile = os.path.join(QMin['scratchdir'], 'JOB/tden_summ.txt')
        omffile = os.path.join(QMin['scratchdir'], 'JOB/OmFrag.txt')
        props = get_props(sumfile, omffile, QMin)
        for i in range(QMin['nmstates']):
            m1, s1, ms1 = tuple(QMin['statemap'][i + 1])
            if (m1, s1) in props:
                for j in range(QMin['template']['theodore_n']):
                    QMout['theodore'][i][j] = props[(m1, s1)][j]
    return QMout

# ======================================================================= #


def get_props(sumfile, omffile, QMin):
    out = readfile(sumfile)
    props = {}
    for line in out[2:]:
        s = line.replace('(', ' ').replace(')', ' ').split()
        if len(s) == 0:
            continue
        n = int(s[0])
        m = int(s[1])
        props[(m, n + (m == 1))] = [theo_float(i) for i in s[5:]]

    out = readfile(omffile)
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

# ======================================================================= #


def copy_ntos(QMin):

    # create directory
    moldendir = QMin['savedir'] + '/MOLDEN/'
    if not os.path.isdir(moldendir):
        mkdir(moldendir)

    # save the nto_x-x.a.mld files
    for i in QMin['statemap']:
        m, s, ms = QMin['statemap'][i]
        if m == 1 and s == 1:
            continue
        if m > 1 and ms != float(m - 1) / 2:
            continue
        f = os.path.join(QMin['scratchdir'], 'JOB', 'nto_%i-%i-a.mld' % (s - (m == 1), m))
        fdest = moldendir + '/step_%s__nto_%i_%i.molden' % (QMin['step'][0], m, s)
        shutil.copy(f, fdest)

# ======================================================================= #
# ======================================================================= #
# ======================================================================= #


class civfl_ana:
    def __init__(self, path, imult, maxsqnorm=1.0, debug=False, filestr='CCRE0'):
        self.det_dict = {}  # dictionary with determinant strings and cicoefficient information
        self.sqcinorms = {}  # CI-norms
        self.path = path
        if imult not in [1, 3]:
            print('CCR* file readout implemented only for singlets and triplets!')
            sys.exit(106)
        self.mult = imult
        self.maxsqnorm = maxsqnorm
        self.debug = debug
        self.nmos = -1  # number of MOs
        self.nfrz = 0  # number of frozen orbs
        self.nocc = -1  # number of occupied orbs (including frozen)
        self.nvir = -1  # number of virtuals
        self.filestr = filestr
        self.read_control()
# ================================================== #

    def read_control(self):
        '''
        Reads nmos, nfrz, nvir from control file
        '''
        controlfile = os.path.join(self.path, 'control')
        control = readfile(controlfile)
        for iline, line in enumerate(control):
            if 'nbf(AO)' in line:
                s = line.split('=')
                self.nmos = int(s[-1])
            if '$closed shells' in line:
                s = control[iline + 1].split()[1].split('-')
                self.nocc = int(s[-1])
            if 'implicit core' in line:
                s = line.split()
                self.nfrz = int(s[2])
        if self.nmos == -1:
            mosfile = os.path.join(self.path, 'mos')
            mos = readfile(mosfile)
            for line in mos:
                if "eigenvalue" in line:
                    self.nmos = int(line.split()[0])
        if any([self.nmos == -1, self.nfrz == -1, self.nocc == -1]):
            print('Number of orbitals not found: nmos=%i, nfrz=%i, nocc=%i' % (self.nmos, self.nfrz, self.nocc))
            sys.exit(107)
        self.nvir = self.nmos - self.nocc
# ================================================== #

    def get_state_dets(self, state):
        """
        Get the transition matrix from CCR* file and add to det_dict.
        """
        if (self.mult, state) == (1, 1):
            det = self.det_string(0, self.nocc, 'de')
            self.det_dict[det] = {1: 1.}
            return
        try:
            filename = ('%s% 2i% 3i% 4i' % (self.filestr, 1, self.mult, state - (self.mult == 1))).replace(' ', '-')
            filename = os.path.join(self.path, filename)
            CCfile = open(filename, 'rb')
        except IOError:
            # if the files are not there, use the right eigenvectors
            filename = ('%s% 2i% 3i% 4i' % ('CCRE0', 1, self.mult, state - (self.mult == 1))).replace(' ', '-')
            filename = os.path.join(self.path, filename)
            CCfile = open(filename, 'rb')
        # skip 8 byte
        CCfile.read(8)
        # read method from 8 byte
        method = str(struct.unpack('8s', CCfile.read(8))[0])
        # skip 8 byte
        CCfile.read(8)
        # read number of CSFs from 4 byte
        nentry = struct.unpack('i', CCfile.read(4))[0]
        # skip 4 byte
        CCfile.read(4)
        # read 8 byte as long int
        versioncheck = struct.unpack('l', CCfile.read(8))[0]
        if versioncheck == 0:
            # skip 16 byte in Turbomole >=7.1
            CCfile.read(16)
        else:
            # skip 8 byte in Turbomole <=7.0
            CCfile.read(8)
        # checks
        if 'CCS' in method:
            print('ERROR: preoptimization vector found in file: %s' % (filename))
            sys.exit(108)
        if not nentry == self.nvir * (self.nocc - self.nfrz):
            print('ERROR: wrong number of entries found in file: %s' % (filename))
        # get data
        state_dict = {}
        nact = self.nocc - self.nfrz
        for iocc in range(nact):
            for ivirt in range(self.nvir):
                coef = struct.unpack('d', CCfile.read(8))[0]
                if self.mult == 1:
                    det = self.det_string(iocc + self.nfrz, self.nocc + ivirt, 'ab')
                    state_dict[det] = coef
                elif self.mult == 3:
                    det = self.det_string(iocc + self.nfrz, self.nocc + ivirt, 'aa')
                    state_dict[det] = coef
        # renormalize
        vnorm = 0.
        for i in state_dict:
            vnorm += state_dict[i]**2
        vnorm = math.sqrt(vnorm)
        # truncate data
        state_dict2 = {}
        norm = 0.
        for i in sorted(state_dict, key=lambda x: state_dict[x]**2, reverse=True):
            state_dict2[i] = state_dict[i] / vnorm
            norm += state_dict2[i]**2
            if norm > self.maxsqnorm:
                break
        # put into general det_dict, also adding the b->a excitation for singlets
        if self.mult == 1:
            for i in state_dict2:
                coef = state_dict2[i] / math.sqrt(2.)
                j = i.replace('a', 't').replace('b', 'a').replace('t', 'b')
                if i in self.det_dict:
                    self.det_dict[i][state] = coef
                else:
                    self.det_dict[i] = {state: coef}
                if j in self.det_dict:
                    self.det_dict[j][state] = -coef
                else:
                    self.det_dict[j] = {state: -coef}
        elif self.mult == 3:
            for i in state_dict2:
                coef = state_dict2[i]
                if i in self.det_dict:
                    self.det_dict[i][state] = coef
                else:
                    self.det_dict[i] = {state: coef}
# ================================================== #

    def det_string(self, fromorb, toorb, spin):
        if fromorb >= self.nocc or toorb < self.nocc or fromorb >= self.nmos or toorb >= self.nmos:
            print('Error generating determinant string!')
            sys.exit(109)
        string = 'd' * self.nocc + 'e' * (self.nmos - self.nocc)
        string = string[:fromorb] + spin[0] + string[fromorb + 1:toorb] + spin[1] + string[toorb + 1:]
        return string
# ================================================== #

    def sort_key(self, key):
        """
        For specifying the sorting order of the determinants.
        """
        return key.replace('d', '0').replace('a', '1').replace('b', '1')
# ================================================== #

    def sort_key2(self, key):
        """
        For specifying the sorting order of the determinants.
        """
        return key.replace('d', '0').replace('a', '0').replace('b', '1').replace('e', '1')
# ================================================== #

    def write_det_file(self, nstate, wname='dets', wform=' % 14.10f'):
        string = '%i %i %i\n' % (nstate, self.nmos, len(self.det_dict))
        for det in sorted(sorted(self.det_dict, key=self.sort_key2), key=self.sort_key):
            string += det
            for istate in range(1, nstate + 1):
                try:
                    string += wform % (self.det_dict[det][istate])
                except KeyError:
                    string += wform % (0.)
            string += '\n'
        writefile(wname, string)


# ========================== Main Code =============================== #
def main():

    # Retrieve PRINT and DEBUG
    try:
        envPRINT = os.getenv('SH2CC2_PRINT')
        if envPRINT and envPRINT.lower() == 'false':
            global PRINT
            PRINT = False
        envDEBUG = os.getenv('SH2CC2_DEBUG')
        if envDEBUG and envDEBUG.lower() == 'true':
            global DEBUG
            DEBUG = True
    except ValueError:
        print('PRINT or DEBUG environment variables do not evaluate to logical values!')
        sys.exit(110)

    # Process Command line arguments
    if len(sys.argv) != 2:
        print('Usage:\n./SHARC_RICC2.py <QMin>\n')
        print('version:', version)
        print('date:', versiondate)
        print('changelog:\n', changelogstring)
        sys.exit(111)
    QMinfilename = sys.argv[1]

    # Print header
    printheader()

    # Read QMinfile
    QMin = readQMin(QMinfilename)

    # Process Tasks
    Tasks = gettasks(QMin)
    if DEBUG:
        pprint.pprint(Tasks)

    # do all runs
    QMin, QMout = runeverything(Tasks, QMin)

    printQMout(QMin, QMout)

    # Measure time
    runtime = measuretime()
    QMout['runtime'] = runtime

    # Write QMout
    writeQMout(QMin, QMout, QMinfilename)

    if PRINT or DEBUG:
        print(datetime.datetime.now())
        print('#================ END ================#')


if __name__ == '__main__':
    main()
