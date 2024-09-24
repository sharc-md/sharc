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

# ======================================================================= #
# Modules:
# Operating system, isfile and related routines, move files, create directories
import os
import shutil
# External Calls to MOLPRO
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
# gethostname routine
from socket import gethostname
# combinatorics
import itertools
# parallel calculations
from multiprocessing import Pool
import time
import traceback


# =========================================================0

version = '3.0'
versiondate = datetime.date(2023, 4, 1)


changelogstring = '''
06.02.:
- changed the records for cpmcscf to 5xxx.9

07.02.:
- added angular keyword (angular momentum evaulation and extraction)
- extraction of best obtained accuracy in cpmcscf
- nogprint,orbitals,civectors added in MOLPRO input

08.02.:
- added removal of SCRATCHDIR after job finished successfully
- added expansion of environment variables and ~ in paths

13.03.:
- added facility for selective analytical NACs
- added input for nac ana select
- added environment variables for PRINT and DEBUG

08.05.:
- added CI:pspace task

22.05.:
- Smat is written transposed now (to be compatible with Fortran)

06.06.:
- MCSCF convergence thresholds increased by default, to help cpmcscf convergence

11.10.:
- changed keyword "nac" to "nacdr", "nacdt" and "overlap"
=>NOT COMPATIBLE WITH OLDER VERSIONS OF SHARC!

19.02.2014:
- modified qmout write routines to be compatible with the new SHARC

11.03.2014:
- changed keyword "restart" to "samestep" to avoid ambiguity
=>NOT COMPATIBLE WITH OLDER VERSIONS OF SHARC!

11.06.2014:
- "grad" can now have no arguments. "grad" is equivalent to "grad all"

16.07.2014:
- savedir from QM.in or SH2PRO.inp

22.09.2014:
- improved calculation of overlap matrices

08.10.2014:     1.0
- official release version, no changes to 0.2

18.12.2014:
- fixed a bug where CPMCSCF solutions converging in zero iterations (full CI case) are not treated properly
- fixed a bug where gradients are not read out if "grad" is given without specifying "all" or the states

04.07.2016:
- COMPLETE REWORK
- multi-job capabilities (several independent CASSCF calcs for different multiplicities)
- overlaps through wfoverlap
- Dyson norms for ionization
- NACdt (plus checknacs stuff), angular keyword decrepated
- phase adjustment between CASSCF and MRCI calculations
- parallelization (multi-job, multi-grad)
- new Input file structure (like MOLCAS interface)
=> NOT BACKWARDS-COMPATIBLE WITH OLDER MOLPRO INTERFACE INPUTS

27.09.2016:
- added "basis_external" keyword for MOLPRO.template, which allows specifying a file whose content is taken as the basis set definition

23.08.2017:
- Resource file is now called "MOLPRO.resources" instead of "SH2PRO.inp" (old filename still works)

24.08.2017:
- added the numfrozcore and numocc keywords for Dyson norm and overlap calculations
- gradaccumax and gradaccudefault are now keywords in the template, not in the resources file (not backwards compatible)
'''

# ======================================================================= #
# holds the system time when the script was started
starttime = datetime.datetime.now()

# global variables for printing (PRINT gives formatted output, DEBUG gives raw output)
DEBUG = False
PRINT = True

# hash table for conversion of multiplicity to the keywords used in MOLPRO
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

# hash table for conversion of polarisations to the keywords used in MOLPRO
IToPol = {
    0: 'X',
    1: 'Y',
    2: 'Z',
    'X': 0,
    'Y': 1,
    'Z': 2
}

# conversion factors
au2a = 0.529177211
rcm_to_Eh = 4.556335e-6

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
           'La': 23, 'Hf': 23, 'Ta': 23, 'W': 23, 'Re': 23, 'Os': 23, 'Ir': 23, 'Pt': 23, 'Au': 23, 'Hg': 23,
           'Tl': 23, 'Pb': 23, 'Bi': 23, 'Po': 23, 'At': 23, 'Rn': 23
           }

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


def link(PATH, NAME, crucial=True, force=True):
    # do not create broken links
    if not os.path.exists(PATH):
        print('Source %s does not exist, cannot create link!' % (PATH))
        sys.exit(14)
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
                    sys.exit(15)
                else:
                    return
    elif os.path.exists(NAME):
        # NAME is not a link. The interface will not overwrite files/directories with links, even with force=True
        print('%s exists, cannot create a link of the same name!' % (NAME))
        if crucial:
            sys.exit(16)
        else:
            return
    os.symlink(PATH, NAME)

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
    string += '||' + ' ' * 25 + 'SHARC - MOLPRO2012 - Interface' + ' ' * 25 + '||\n'
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
    '''If PRINT, prints a formatted Summary of the control variables read from the input file.

    Arguments:
    1 dictionary: QMin'''

    if DEBUG:
        pprint.pprint(QMin)
    if not PRINT:
        return
    print('==> QMin Job description for:\n%s' % (QMin['comment']))

    string = 'Tasks:'
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
    if 'phases' in QMin:
        string += '\tPhases'
    print(string)

    string = 'States: '
    for i in itmult(QMin['states']):
        string += '\t%i %s' % (QMin['states'][i - 1], IToMult[i])
    print(string)

    string = '\nMethod:\n'
    string += 'Electrons per Multiplicity:'
    for m in range(QMin['maxmult']):
        if m + 1 in QMin['multmap']:
            job = QMin['multmap'][m + 1]
            string += '\t%i' % (QMin['template']['nelec'][job][m])
        else:
            string += '\t/'
    string += '\nJob ID per Multiplicity:'
    for m in range(QMin['maxmult']):
        string += '\t%i' % (QMin['jobs'][m])
    string += '\nLevel of Theory per Job ID:\n'
    for j in QMin['joblist']:
        string += '%i:\tSA(%i' % (j, QMin['template']['roots'][j][0])
        for i in QMin['template']['roots'][j][1:]:
            string += '+%i' % (i)
        nelec = QMin['template']['nelec'][j][QMin['multmap'][-j][0] - 1] - 2 * QMin['template']['closed'][j - 1]
        norb = QMin['template']['occ'][j - 1] - QMin['template']['closed'][j - 1]
        string += ')-CASSCF(%i,%i)' % (nelec, norb)
        if 'basis_external' in QMin['template']:
            string += '/%s\n' % ('CUSTOM BASIS')
        else:
            string += '/%s\n' % (QMin['template']['basis'])
    if QMin['template']['dkho'] > 0:
        string += 'Using Douglas-Kroll-Hess Hamiltonian\n'
    print(string)

    string = 'Found Geo'
    if 'veloc' in QMin:
        string += ' and Veloc! '
    else:
        string += '! '
    string += 'NAtom is %i.\n' % (QMin['natom'])
    print(string)

    string = 'Geometry in Bohrs:\n'
    for i in range(QMin['natom']):
        string += '%s ' % (QMin['geo'][i][0])
        for j in range(3):
            string += '% 7.4f ' % (QMin['geo'][i][j + 1])
        string += '\n'
    print(string)

    if 'veloc' in QMin:
        string = ''
        for i in range(QMin['natom']):
            string += '%s ' % (QMin['geo'][i][0])
            for j in range(3):
                string += '% 7.4f ' % (QMin['veloc'][i][j])
            string += '\n'
        print(string)

    if 'grad' in QMin:
        string = 'Gradients:   '
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

    # if 'overlap' in QMin:
        # string='Overlaps:\n'
        # for i in range(1,QMin['nmstates']+1):
        # for j in range(1,QMin['nmstates']+1):
        # if [i,j] in QMin['overlap'] or [j,i] in QMin['overlap']:
        #string+='X '
        # else:
        #string+='. '
        # string+='\n'
        # print(string)

    do_not_print = ['h', 'soc', 'dm', 'geo', 'veloc', 'states', 'comment', 'grad', 'nacdr', 'ion', 'overlap', 'nacdt']
    print_if_debug = ['ionlist']
    for i in QMin:
        if i not in do_not_print:
            if i not in print_if_debug or DEBUG:
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


def printQMout(QMin, QMout):
    '''If PRINT, prints a summary of all requested QM output values. Matrices are formatted using printcomplexmatrix, vectors using printgrad.

    Arguments:
    1 dictionary: QMin
    2 dictionary: QMout'''

    if DEBUG:
        pprint.pprint(QMout)
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
    # overlaps
    if 'overlap' in QMin:
        print('=> Overlap matrix:\n')
        matrix = QMout['overlap']
        printcomplexmatrix(matrix, states)
        if 'phases' in QMout:
            print('=> Wavefunction Phases:\n%i\n' % (nmstates))
            for i in range(nmstates):
                print('% 3.1f % 3.1f' % (QMout['phases'][i].real, QMout['phases'][i].imag))
            print('\n')
    # Property matrix (dyson norms)
    if 'ion' in QMin and 'prop' in QMout:
        print('=> Property matrix:\n')
        matrix = QMout['prop']
        printcomplexmatrix(matrix, states)
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


# ======================================================================= #     OK
def nextblock(data, program='*', occ=1):
    '''Scans the list of strings data for the next occurence of MOLPRO program block for program. Returns the line number where the block ends and the block itself.

    Arguments:
    1 list of strings: data
    2 string: MOLPRO program name (like "MULTI" in "1PROGRAM * MULTI" )
    3 integer: return the i-th block of the specified program

    Returns:
    1 integer: line number in data where the specified block ends
    2 list of strings: The specified block'''

    progdata = []
    i = -1
    while True:
        try:
            i += 1
            line = data[i].split()
            if line == []:
                continue
            if containsstring('1PROGRAM', line[0]) and containsstring(program, line[2]):
                occ -= 1
                if occ == 0:
                    break
        except IndexError:
            print('Block %s not found in routine nextblock! Probably MOLPRO encountered an error not anticipated in this script. Check the MOLPRO output!' % (program))
            sys.exit(17)
    progdata.append(data[i])
    i += 1
    while i < len(data) and not containsstring('1PROGRAM', data[i]):
        progdata.append(data[i])
        i += 1
    return i, progdata


# ======================================================================= #
def getcienergy(out, mult, state):

    ilines = 0
    # look for CI program block
    while ilines < len(out):
        if containsstring(r'PROGRAM \* CI', out[ilines]):
            # look for multiplicity
            while ilines < len(out):
                if containsstring('Reference symmetry', out[ilines]):
                    if containsstring(IToMult[mult], out[ilines]):
                        # look for energy
                        while ilines < len(out):
                            # if '********************************************************' in out[ilines]:
                            # break
                            if containsstring(r'!(MRCI|CI\(SD\)) STATE ?[0-9]+\.1 Energy', out[ilines]):
                                kstate = int(out[ilines].replace('.', ' ').replace('E', ' ').split()[2])
                                if kstate == state:
                                    return float(out[ilines].split()[-1])
                            ilines += 1
                    else:
                        break
                ilines += 1
        ilines += 1
    print('CI energy of state %i in mult %i not found!' % (state, mult))
    sys.exit(18)

# ======================================================================= #


def getcidm(out, mult, state1, state2, pol):
    '''Searches a complete MOLPRO output file for a cartesian component of a dipole moment between the two specified states.

    Only takes one multiplicity, since in this script, only non-relativistic dipole moments are calculated.
    If state1==state2, then this returns a state dipole moment, otherwise a transition dipole moment.

    Arguments:
    1 list of strings: MOLPRO output
    2 integer: mult
    3 integer: state1
    4 integer: state2
    5 integer (0,1,2) or character (X,Y,Z): Polarisation

    Returns:
    1 float: cartesian dipole moment component in atomic units'''

    if pol == 'X' or pol == 'Y' or pol == 'Z':
        pol = IToPol[pol]
    ilines = 0
    while ilines < len(out):
        if containsstring(r'PROGRAM \* CI', out[ilines]):
            while ilines < len(out):
                if containsstring('Reference symmetry', out[ilines]):
                    if containsstring(IToMult[mult], out[ilines]):
                        # expectation values are in the results section, transition moments seperately
                        if state1 == state2:
                            while not containsstring(r'\*\*\*', out[ilines]):
                                if containsstring(r'!.* STATE ?[0-9]+\.1 Dipole moment', out[ilines]):
                                    kstate = int(out[ilines].replace('.', ' ').replace('E', ' ').split()[2])
                                    if kstate == state1:
                                        return float(out[ilines].split()[-3 + pol])
                                ilines += 1
                        else:
                            while not containsstring(r'\*\*\*', out[ilines]):
                                if containsstring(r'MRCI trans.*<.*\|DM.\|.*>', out[ilines]):
                                    braket = out[ilines].replace('<', ' ').replace('>', ' ').replace('|', ' ').replace('.', ' ').split()
                                    s1 = int(braket[2])
                                    s2 = int(braket[5])
                                    p = IToPol[braket[4][2]]
                                    if p == pol and ((s1 == state1 and s2 == state2) or (s1 == state2 and s2 == state1)):
                                        return float(out[ilines].split()[3])
                                ilines += 1
                            # return 0.
                    else:
                        break
                ilines += 1
        ilines += 1
    if state1 == state2:
        print('Dipole moment of state %i, mult %i, not found!' % (state1, mult))
        sys.exit(19)
    else:
        # print('CI dipole moment of states %i and %i in mult %i not found!' % (state1,state2,mult))
        return 0.

# ======================================================================= #


def getciang(out, mult, state1, state2, pol):
    '''Searches a complete MOLPRO output file for a cartesian component of a angular momentum between the two specified states.

    Only takes one multiplicity, since in this script, only non-relativistic angular momenta are calculated.
    If state1==state2, then this returns zero, otherwise a transition angular momentum.

    Arguments:
    1 list of strings: MOLPRO output
    2 integer: mult
    3 integer: state1
    4 integer: state2
    5 integer (0,1,2) or character (X,Y,Z): Polarisation

    Returns:
    1 complex: cartesian angular momentum component in atomic units'''

    if state1 == state2:
        return complex(0., 0.)
    if pol == 'X' or pol == 'Y' or pol == 'Z':
        pol = IToPol[pol]
    ilines = 0
    while ilines < len(out):
        if containsstring(r'PROGRAM \* CI', out[ilines]):
            while ilines < len(out):
                if containsstring('Reference symmetry', out[ilines]):
                    if containsstring(IToMult[mult], out[ilines]):

                        while not containsstring(r'\*\*\*', out[ilines]):
                            if containsstring(r'MRCI trans.*<.*\|L.\|.*>', out[ilines]):
                                braket = out[ilines].replace('<', ' ').replace('>', ' ').replace('|', ' ').replace('.', ' ').split()
                                s1 = int(braket[2])
                                s2 = int(braket[5])
                                p = IToPol[braket[4][1]]
                                if p == pol and s1 == state1 and s2 == state2:
                                    return complex(out[ilines].split()[3].replace('i', 'j'))
                                if p == pol and s1 == state2 and s2 == state1:
                                    return -complex(out[ilines].split()[3].replace('i', 'j')).conjugate()
                            ilines += 1
                        return complex(0., 0.)
                    else:
                        break
                ilines += 1
        ilines += 1
    print('CI angular momentum of states %i and %i in mult %i not found!' % (state1, state2, mult))
    sys.exit(20)

# ======================================================================= #


def istate_in_job(m1, s1, ms1, states):
    k = -1
    for im in range(m1):
        for ims in range(im + 1):
            for ist in range(states[im]):
                k += 1
                if (im + 1, ims - im / 2., ist + 1) == (m1, float(ms1), s1):
                    return k

# ======================================================================= #
# def getsocme(out,mstate1,mstate2,states):


def getsocme(out, istate, jstate, QMin):

    istate += 1
    jstate += 1
    if istate not in QMin['statemap'] or jstate not in QMin['statemap']:
        print('States %i or %i are not in statemap!' % (istate, jstate))
        sys.exit(21)

    m1, s1, ms1 = tuple(QMin['statemap'][istate])
    m2, s2, ms2 = tuple(QMin['statemap'][jstate])
    job1 = QMin['multmap'][m1]
    job2 = QMin['multmap'][m2]

    if job1 != job2:
        return complex(0., 0.)
    job = job1

    # get numbers of states within job
    mults = QMin['multmap'][-job]
    states_in_job = [0 for imult in range(max(mults))]
    for imult in mults:
        states_in_job[imult - 1] = QMin['states'][imult - 1]
    nmstates = sum([(imult + 1) * ix for imult, ix in enumerate(states_in_job)])

    i = istate_in_job(m1, s1, ms1, states_in_job)
    j = istate_in_job(m2, s2, ms2, states_in_job)
    # print istate,jstate,i,j

    iline = -1
    while iline < len(out):
        iline += 1
        line = out[iline]
        if '   Spin-orbit calculation in the basis of zeroth order wave functions' in line:
            break
    else:
        print('SOC matrix not found in master_%i/MOLPRO.out!' % (job))
        sys.exit(22)
    iline += 3
    eref = float(out[iline].replace('=', ' ').split()[-1])

    iline = -1
    while iline < len(out):
        iline += 1
        line = out[iline]
        if ' Spin-Orbit Matrix (CM-1)' in line:
            # if '  State Sym Spin    / Nr.' in line:
            break
    else:
        print('SOC matrix not found in master_%i/MOLPRO.out!' % (job))
        sys.exit(23)
    iline += 5
    # iline+=2

    rcm_to_Eh = 4.556335e-6
    # get a single matrix element
    block = (j) // 10
    yoffset = (i) * 3 + block * (3 * nmstates + 3)
    xoffset = (j) % 10
    # block=(j)/8
    #yoffset=(i)*3 + block*(3*nmstates+2)
    # xoffset=(j)%8
    # print iline
    # print block,xoffset,yoffset
    # print out[iline+yoffset]
    # print out[iline+yoffset+1]
    real = float(out[iline + yoffset].split()[4 + xoffset]) * rcm_to_Eh
    imag = float(out[iline + yoffset + 1].split()[xoffset]) * rcm_to_Eh
    if istate == jstate:
        real += eref
    return complex(real, imag)




    # sys.exit(24)

    # iline=-1
    # while True:
    # while True:
    # iline+=1
    # if iline==len(out):
    # print('No Spin-Orbit CI output found for multiplicities %i and %i!' % (mult1,mult2))
    # sys.exit(25)
    # if 'SEWLS' in out[iline]:
    # break
    # iline+=8
    # found=[False,False]
    # while True:
    # iline+=1
    # if not 'Wavefunction restored' in out[iline]:
    # break
    # mult=int(float(out[iline].split()[7])*2+1)
    # if mult==mult1:
    # found[0]=True
    # if mult==mult2:
    # found[1]=True
    # if all(found):
    # break
    # while True:
    # iline+=1
    # if 'Lowest unperturbed energy E0=' in out[iline]:
    # eref=complex(float(out[iline].split()[4]),0)
    # break
    # while True:
    # iline+=1
    # if r'Spin-Orbit Matrix (CM-1)' in out[iline]:
    # break
    # iline+=5

    # rcm_to_Eh=4.556335e-6

    # mstate1=0
    # mstate2=0
    # i=0
    # for imult,istate,ims in itnmstates(states):
    # i+=1
    # if (imult,istate,ims)==(mult1,state1,ms1):
    # mstate1=i
    # if (imult,istate,ims)==(mult2,state2,ms2):
    # mstate2=i
    # nmstates=i

    # if mstate1==0:
    # print('Mult %i, State %i, MS %i not in SOC matrix after line %i!' % (mult1,state1,ms1,iline))
    # sys.exit(26)
    # if mstate2==0:
    # print('Mult %i, State %i, MS %i not in SOC matrix after line %i!' % (mult2,state2,ms2,iline))
    # sys.exit(27)

    # get a single matrix element
    # block=(mstate2-1)/10
    #yoffset=(mstate1-1)*3 + block*(3*nmstates+3)
    # xoffset=(mstate2-1)%10
    # real=float(out[iline+yoffset].split()[4+xoffset])*rcm_to_Eh
    # imag=float(out[iline+yoffset+1].split()[xoffset])*rcm_to_Eh
    # if mstate1==mstate2:
    # real+=eref
    # return complex(real,imag)

# ======================================================================= #
def getgrad(out, mult, state, natom):
    '''Searches a MOLPRO output for a SA-MCSCF gradient of a specified state.

    Arguments:
    1 list of strings: MOLPRO output
    2 integer: mult
    3 integer: state
    4 integer: natom

    Returns:
    1 list of list of floats: gradient vector (natom x 3) in atomic units'''

    ilines = 0
    grad = []
    multfound = False
    statefound = False
    # look for FORCE program block
    while ilines < len(out):
        if containsstring(r'PROGRAM \* ALASKA', out[ilines]):
            # look for multiplicity and state
            jlines = ilines
            while not containsstring(r'\*\*\*', out[jlines]):
                if containsstring(IToMult[mult], out[jlines]):
                    multfound = True
                    break
                jlines += 1
            jlines = ilines
            while not containsstring(r'\*\*\*', out[jlines]):
                if containsstring('SA-MC GRADIENT FOR STATE', out[jlines]):
                    line = out[jlines].replace('E', ' ').replace('.', ' ').split()
                    if state == int(line[5]):
                        statefound = True
                    break
                jlines += 1
            if multfound and statefound:
                jlines += 4
                for i in range(natom):
                    line = out[jlines + i].split()
                    for j in range(3):
                        try:
                            line[j + 1] = float(line[j + 1])
                        except ValueError:
                            print('Bad float in gradient in line %i! Maybe natom is wrong.' % (ilines + i))
                    grad.append(line[1:])
                return grad
            else:
                multfound = False
                statefound = False
                ilines += 1
                continue
        ilines += 1
    print('Gradient of state %i in mult %i not found!' % (state, mult))
    sys.exit(28)

# ======================================================================= #


def getnacana(out, mult, state1, state2, natom):
    '''Searches a MOLPRO output file for an analytical non-adiabatic coupling vector from SA-MCSCF.

    Arguments:
    1 list of strings: MOLPRO output
    2 integer: mult
    3 integer: state1
    4 integer: state2
    5 integer: natom

    Returns:
    1 list of list of floats: non-adiabatic coupling vector (natom x 3) in atomic units'''

    ilines = 0
    grad = []
    multfound = False
    statefound = False
    # diagonal couplings are zero
    if state1 == state2:
        for i in range(natom):
            grad.append([0., 0., 0.])
        return grad
    # look for FORCE program block
    while ilines < len(out):
        if containsstring(r'PROGRAM \* ALASKA', out[ilines]):
            # look for multiplicity and state
            jlines = ilines
            while not containsstring(r'\*\*\*', out[jlines]):
                if containsstring(IToMult[mult], out[jlines]):
                    multfound = True
                    break
                jlines += 1
            jlines = ilines
            while not containsstring(r'\*\*\*', out[jlines]):
                if containsstring('SA-MC NACME FOR STATES', out[jlines]):
                    line = out[jlines].replace('.', ' ').replace('-', ' ').split()
                    # make sure the NACs are antisymmetric
                    if state1 == int(line[5]) and state2 == int(line[7]):
                        statefound = True
                        factor = 1.
                    if state1 == int(line[7]) and state2 == int(line[5]):
                        statefound = True
                        factor = -1.
                    break
                jlines += 1
            if multfound and statefound:
                jlines += 4
                for i in range(natom):
                    line = out[jlines + i].split()
                    for j in range(3):
                        try:
                            line[j + 1] = factor * float(line[j + 1])
                        except ValueError:
                            print('Bad float in gradient in line %i! Maybe natom is wrong.' % (ilines + i))
                    grad.append(line[1:])
                return grad
            else:
                multfound = False
                statefound = False
                ilines += 1
                continue
        ilines += 1
    print('Non-adiatic coupling of states %i - %i in mult %i not found!' % (state1, state2, mult))
    sys.exit(29)

# ======================================================================= #


def getmrcioverlap(out, mult, state1, state2):
    '''Searches a MOLPRO output for a single MRCI overlap (from a CI trans calculation).

    Arguments:
    1 list of strings: MOLPRO output
    2 integer: mult
    3 integer: state1
    4 integer: state2

    Returns:
    1 float: MRCI overlap (THIS MATRIX IS NOT SYMMETRIC!)'''

    ilines = 0
    while ilines < len(out):
        if containsstring('Ket wavefunction restored from record .*\\.3', out[ilines]):
            line = out[ilines].replace('.', ' ').split()
            if mult == int(line[5]) - 6000:
                break
        ilines += 1
    while not containsstring(r'\*\*\*', out[ilines]):
        if containsstring('!MRCI overlap', out[ilines]):
            braket = out[ilines].replace('<', ' ').replace('>', ' ').replace('|', ' ').replace('.', ' ').split()
            s1 = int(braket[2])
            s2 = int(braket[4])
            # overlap matrix is NOT symmetric!
            if s1 == state1 and s2 == state2:
                return float(out[ilines].split()[3])
        ilines += 1
    return 0.

# ======================================================================= #


def getnacnum(out, mult, state1, state2):
    '''Searches a MOLPRO output for a single non-adiabatic coupling matrix element from a DDR calculation.

    Arguments:
    1 list of strings: MOLPRO output
    2 integer: mult
    3 integer: state1
    4 integer: state2

    Returns:
    1 float: NAC matrix element'''

    # diagonal couplings are zero
    if state1 == state2:
        return 0.
    ilines = 0
    multfound = False
    statefound = False
    while ilines < len(out):
        if containsstring('Construct non-adiabatic coupling elements by finite difference method', out[ilines]):
            jlines = ilines
            while not containsstring('\\*\\*\\*', out[jlines]):
                if containsstring('Transition density \\(R\\|R\\+DR\\)', out[jlines]):
                    line = out[jlines].replace('.', ' ').replace('-', ' ').split()
                    if mult == int(line[4]) - 8000:
                        multfound = True
                    if state1 == int(line[8]) and state2 == int(line[10]):
                        statefound = True
                        factor = 1.
                    if state1 == int(line[10]) and state2 == int(line[8]):
                        statefound = True
                        factor = -1.
                    if multfound and statefound:
                        jlines += 5
                        return factor * float(out[jlines].split()[2])
                    else:
                        multfound = False
                        statefound = False
                        ilines += 1
                jlines += 1
        ilines += 1
    print('Non-adiatic coupling of states %i - %i in mult %i not found!' % (state1, state2, mult))
    sys.exit(30)

# ======================================================================= #


def getsmate(out, s1, s2):
    ilines = -1
    while True:
        ilines += 1
        if ilines == len(out):
            print('Overlap of states %i - %i not found!' % (s1, s2))
            sys.exit(31)
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
            sys.exit(32)
        if containsstring('Dyson norm matrix <PsiA_i|PsiB_j>', out[ilines]):
            break
    ilines += 1 + s1
    f = out[ilines].split()
    return float(f[s2 + 1])





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

    # add header info
    string += '! 0 Basic information\nstates '
    for i in QMin['states']:
        string += '%i ' % i
    string += '\nnmstates %i\n' % QMin['nmstates']
    string += 'natom %i\n' % QMin['natom']
    string += 'npc 0\n'
    string += 'charges '
    for i in QMin['states']:
        string += '%i ' % 0
    string += '\n\n'

    if 'h' in QMin or 'soc' in QMin:
        string += writeQMoutsoc(QMin, QMout)
    if 'dm' in QMin:
        string += writeQMoutdm(QMin, QMout)
    if 'grad' in QMin:
        string += writeQMoutgrad(QMin, QMout)
    if 'nacdr' in QMin:
        string += writeQMoutnacana(QMin, QMout)
    if 'overlap' in QMin:
        string += writeQMoutnacsmat(QMin, QMout)
    if 'ion' in QMin:
        string += writeQMoutprop(QMin, QMout)
    if 'phases' in QMin:
        string += writeQmoutPhases(QMin, QMout)
    string += writeQMouttime(QMin, QMout)
    writefile(outfilename, string)
    # try:
    # outfile=open(outfilename,'w')
    # outfile.write(string)
    # outfile.close()
    # except IOError:
    # print('Could not write QM output!')
    # sys.exit(33)
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
        string += '%i %i ! %i %i %i\n' % (natom, 3, imult, istate, ims)
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
    string += '! %i Wavefunction phases (%i, complex)\n%i\n' % (7, nmstates, nmstates)
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

    string = '! 8 Runtime\n%s\n' % (eformat(QMout['runtime'], 12, 3))
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
    return string

# ======================================================================= #


def writeQmoutPhases(QMin, QMout):

    string = '! 7 Phases\n%i ! for all nmstates\n' % (QMin['nmstates'])
    for i in range(QMin['nmstates']):
        string += '%s %s\n' % (eformat(QMout['phases'][i].real, 9, 3), eformat(QMout['phases'][i].imag, 9, 3))
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
            sys.exit(34)
    else:
        try:
            os.makedirs(SCRATCHDIR)
        except OSError:
            print('Can not create SCRATCHDIR=%s\n' % (SCRATCHDIR))
            sys.exit(35)

# ======================================================================= #


def removequotes(string):
    if string.startswith("'") and string.endswith("'"):
        return string[1:-1]
    elif string.startswith('"') and string.endswith('"'):
        return string[1:-1]
    else:
        return string

# ======================================================================= #


def getsh2prokey(sh2pro, key):
    i = -1
    while True:
        i += 1
        try:
            line = re.sub('#.*$', '', sh2pro[i])
        except IndexError:
            break
        line = line.split(None, 1)
        if line == []:
            continue
        if key.lower() in line[0].lower():
            return line
    return ['', '']

# ======================================================================= #


def get_sh2pro_environ(sh2pro, key, environ=True, crucial=True):
    line = getsh2prokey(sh2pro, key)
    if line[0]:
        LINE = line[1]
        LINE = removequotes(LINE).strip()
    else:
        if environ:
            LINE = os.getenv(key.upper())
            if not LINE:
                if crucial:
                    print('Either set $%s or give path to %s in SH2CAS.inp!' % (key.upper(), key.upper()))
                    sys.exit(36)
                else:
                    return ''
        else:
            if crucial:
                print('Give path to %s in SH2CAS.inp!' % (key.upper()))
                sys.exit(37)
            else:
                return ''
    LINE = os.path.expandvars(LINE)
    LINE = os.path.expanduser(LINE)
    if containsstring(';', LINE):
        print("$%s contains a semicolon. Do you probably want to execute another command after %s? I can't do that for you..." % (key.upper(), key.upper()))
        sys.exit(38)
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
            sys.exit(39)
        if 'end' in line:
            break
        fields = line.split()
        try:
            nacpairs.append([int(fields[0]), int(fields[1])])
        except ValueError:
            print('"nacdr select" is followed by pairs of state indices, each pair on a new line!')
            sys.exit(40)
    return nacpairs, i

# ======================================================================= #     OK


def readQMin(QMinfilename):
    '''Reads the time-step dependent information from QMinfilename. This file contains all information from the current SHARC job: geometry, velocity, number of states, requested quantities along with additional information. The routine also checks this input and obtains a number of environment variables necessary to run MOLPRO.

    Steps are:
    - open and read QMinfilename
    - Obtain natom, comment, geometry (, velocity)
    - parse remaining keywords from QMinfile
    - check keywords for consistency, calculate nstates, nmstates
    - obtain environment variables for path to MOLPRO and scratch directory, and for error handling

    Arguments:
    1 string: name of the QMin file

    Returns:
    1 dictionary: QMin'''


    # read QMinfile
    QMinlines = readfile(QMinfilename)
    QMin = {}



    # Get natom
    try:
        natom = int(QMinlines[0].split()[0])
    except ValueError:
        print('first line must contain the number of atoms!')
        sys.exit(41)
    QMin['natom'] = natom
    if len(QMinlines) < natom + 4:
        print('Input file must contain at least:\nnatom\ncomment\ngeometry\nkeyword "states"\nat least one task')
        sys.exit(42)



    # Save Comment line
    QMin['comment'] = QMinlines[1]



    # Get geometry and possibly velocity
    QMin['geo'] = []
    QMin['veloc'] = []
    hasveloc = True
    for i in range(2, natom + 2):
        if not containsstring('[a-zA-Z][a-zA-Z]?[0-9]*.*[-]?[0-9]+[.][0-9]*.*[-]?[0-9]+[.][0-9]*.*[-]?[0-9]+[.][0-9]*', QMinlines[i]):
            print('Input file does not comply to xyz file format! Maybe natom is just wrong.')
            sys.exit(43)
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


    # Units conversion
    if 'unit' in QMin:
        if QMin['unit'][0] == 'angstrom':
            factor = 1. / au2a
        elif QMin['unit'][0] == 'bohr':
            factor = 1.
        else:
            print('Dont know input unit %s!' % (QMin['unit'][0]))
            sys.exit(44)
    else:
        factor = 1. / au2a
    for iatom in range(len(QMin['geo'])):
        for ixyz in range(3):
            QMin['geo'][iatom][ixyz + 1] *= factor


    # State number
    if 'states' not in QMin:
        print('Keyword "states" not given!')
        sys.exit(45)
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
    possibletasks = ['h', 'soc', 'dm', 'grad', 'nacdr', 'overlap', 'ion', 'molden', 'phases']
    if not any([i in QMin for i in possibletasks]):
        print('No tasks found! Tasks are "h", "soc", "dm", "grad", "nacdr", "overlap", "ion", and "molden".')
        sys.exit(46)

    if 'samestep' in QMin and 'init' in QMin:
        print('"Init" and "Samestep" cannot be both present in QM.in!')
        sys.exit(47)

    if 'phases' in QMin:
        QMin['overlap'] = []

    if 'overlap' in QMin and 'init' in QMin:
        print('"overlap" and "phases" cannot be calculated in the first timestep! Delete either "overlap" or "init"')
        sys.exit(48)

    if 'init' not in QMin and 'samestep' not in QMin and 'restart' not in QMin:
        QMin['newstep'] = []

    if len(QMin['states']) > 8:
        print('Higher multiplicities than octets are not supported!')
        sys.exit(49)

    if 'h' in QMin and 'soc' in QMin:
        QMin = removekey(QMin, 'h')

    if 'nacdt' in QMin:
        print('Within the SHARC-MOLPRO interface, couplings can only be calculated via "nacdr" or "overlap".')
        sys.exit(50)

    if 'dmdr' in QMin:
        print('Dipole moment gradients not available!')
        sys.exit(51)

    if 'socdr' in QMin:
        print('Spin-orbit gradients not available!')
        sys.exit(52)

    if 'nacdr' in QMin:
        QMin['docicas'] = True

    # Process the gradient requests
    if 'grad' in QMin:
        if len(QMin['grad']) == 0 or QMin['grad'][0] == 'all':
            QMin['grad'] = [i + 1 for i in range(nmstates)]
        else:
            for i in range(len(QMin['grad'])):
                try:
                    QMin['grad'][i] = int(QMin['grad'][i])
                except ValueError:
                    print('Arguments to keyword "grad" must be "all" or a list of integers!')
                    sys.exit(53)
                if QMin['grad'][i] > nmstates:
                    print('State for requested gradient does not correspond to any state in QM input file state list!')
                    sys.exit(54)

    # Process the non-adiabatic coupling requests
    # type conversion has already been done
    if 'nacdr' in QMin:
        if len(QMin['nacdr']) >= 1:
            nacpairs = QMin['nacdr']
            for i in range(len(nacpairs)):
                if nacpairs[i][0] > nmstates or nacpairs[i][1] > nmstates:
                    print('State for requested non-adiabatic couplings does not correspond to any state in QM input file state list!')
                    sys.exit(55)
        else:
            QMin['nacdr'] = [[j + 1, i + 1] for i in range(nmstates) for j in range(i)]

    # obtain the statemap
    statemap = {}
    i = 1
    for imult, istate, ims in itnmstates(QMin['states']):
        statemap[i] = [imult, istate, ims]
        i += 1
    QMin['statemap'] = statemap
    QMin['maxmult'] = max([i[0] for i in QMin['statemap'].values()])

    # get the set of states for which gradients actually need to be calculated
    gradmap = set()
    if 'grad' in QMin:
        for i in QMin['grad']:
            gradmap.add(tuple(statemap[i][0:2]))
    gradmap = sorted(gradmap)
    QMin['gradmap'] = gradmap

    # get the list of statepairs for NACdr calculation
    nacmap = set()
    if 'nacdr' in QMin:
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








    # environment setup

    QMin['pwd'] = os.getcwd()


    # open MOLPRO.resources
    filename = 'MOLPRO.resources'
    if os.path.isfile(filename):
        sh2pro = readfile(filename)
    else:
        print('HINT: reading resources from SH2PRO.inp')
        sh2pro = readfile('SH2PRO.inp')


    # ncpus
    QMin['ncpu'] = 1
    line = getsh2prokey(sh2pro, 'ncpu')
    if line[0]:
        try:
            QMin['ncpu'] = int(line[1])
            QMin['ncpu'] = max(1, QMin['ncpu'])
        except ValueError:
            print('Number of CPUs does not evaluate to numerical value!')
            sys.exit(56)
    os.environ['OMP_NUM_THREADS'] = str(QMin['ncpu'])


    # MOLPRO path
    QMin['molpro'] = get_sh2pro_environ(sh2pro, 'molpro', False)

    # MOLPRO extra arguments
    QMin['molpro_arguments'] = ''
    line = getsh2prokey(sh2pro, 'molpro_arguments')
    if line[0]:
        QMin['molpro_arguments'] = line[1].strip()


    # Set up scratchdir
    line = get_sh2pro_environ(sh2pro, 'scratchdir', False, False)
    if not line:
        line = os.path.join(QMin['pwd'], 'SCRATCHDIR/')
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
        line = get_sh2pro_environ(sh2pro, 'savedir', False, False)
        if not line:
            line = os.path.join(QMin['pwd'], 'SAVEDIR/')
    line = os.path.expandvars(line)
    line = os.path.expanduser(line)
    line = os.path.abspath(line)
    if 'init' in QMin:
        checkscratch(line)
    QMin['savedir'] = line


    # debug keyword in SH2PRO
    line = getsh2prokey(sh2pro, 'debug')
    if line[0]:
        if len(line) <= 1 or 'true' in line[1].lower():
            global DEBUG
            DEBUG = True

    line = getsh2prokey(sh2pro, 'no_print')
    if line[0]:
        if len(line) <= 1 or 'true' in line[1].lower():
            global PRINT
            PRINT = False


    # memory for MOLPRO and wfoverlap
    QMin['memory'] = 100
    line = getsh2prokey(sh2pro, 'memory')
    if line[0]:
        try:
            QMin['memory'] = int(line[1])
            QMin['memory'] = max(100, QMin['memory'])
        except ValueError:
            print('Run memory does not evaluate to numerical value!')
            sys.exit(57)
    else:
        print('WARNING: Please set memory in SH2PRO.inp (in MB)! Using 100 MB default value!')


    # initial MO guess settings
    # if neither keyword is present, the interface will reuse MOs from savedir, or use the EHT guess
    line = getsh2prokey(sh2pro, 'always_orb_init')
    if line[0]:
        QMin['always_orb_init'] = []
    line = getsh2prokey(sh2pro, 'always_guess')
    if line[0]:
        QMin['always_guess'] = []
    if 'always_orb_init' in QMin and 'always_guess' in QMin:
        print('Keywords "always_orb_init" and "always_guess" cannot be used together!')
        sys.exit(58)


    # get the nooverlap keyword: no dets will be extracted if present
    line = getsh2prokey(sh2pro, 'nooverlap')
    if line[0]:
        QMin['nooverlap'] = []


    # wfoverlaps setting
    if 'overlap' in QMin or 'ion' in QMin or 'docicas' in QMin:
        # QMin['wfoverlap']=get_sh2pro_environ(sh2pro,'wfoverlap')
        QMin['wfoverlap'] = get_sh2pro_environ(sh2pro, 'wfoverlap', False, False)
        if not QMin['wfoverlap']:
            ciopath = os.path.join(os.path.expandvars(os.path.expanduser('$SHARC')), 'wfoverlap.x')
            if os.path.isfile(ciopath):
                QMin['wfoverlap'] = ciopath
            else:
                print('Give path to wfoverlap.x in MOLPRO.resources!')
                sys.exit(59)
        # get ncore and ndocc
        line = getsh2prokey(sh2pro, 'numfrozcore')
        if line[0]:
            QMin['ncore'] = int(line[1])
        else:
            QMin['ncore'] = 0
        line = getsh2prokey(sh2pro, 'numocc')
        if line[0]:
            QMin['ndocc'] = int(line[1])
        else:
            QMin['ndocc'] = 0
        if QMin['ncore'] < 0:
            QMin['ncore'] = sum([FROZENS[atom[0]] for atom in QMin['geo']])
        if QMin['ndocc'] < 0:
            pass              # ndocc is dependent on job
        if 'nooverlap' in QMin:
            print('"nooverlap" keyword present, but overlap calculation required.')
            sys.exit(60)
        QMin['wfthres'] = 1.0
        line = getsh2prokey(sh2pro, 'wfthres')
        if line[0]:
            QMin['wfthres'] = float(line[1])



    # job delay
    QMin['delay'] = 0.0
    line = getsh2prokey(sh2pro, 'delay')
    if line[0]:
        try:
            QMin['delay'] = float(line[1])
        except ValueError:
            print('Submit delay does not evaluate to numerical value!')
            sys.exit(61)



    # open template
    template = readfile('MOLPRO.template')
    temp = []
    for line in template:
        line = re.sub('#.*$', '', line).split()
        if len(line) == 0:
            continue
        temp.append(line)
    QMin['template'] = {}

    # first collect the "simple" inputs
    integers = ['dkho']
    strings = ['basis', 'basis_external']
    floats = ['gradaccudefault', 'gradaccumax']
    booleans = []
    for i in booleans:
        QMin['template'][i] = False
    QMin['template']['dkho'] = 0
    QMin['template']['gradaccudefault'] = 1e-7
    QMin['template']['gradaccumax'] = 1e-2
    for line in temp:
        if line[0].lower() in integers:
            QMin['template'][line[0]] = int(line[1])
        elif line[0].lower() in booleans:
            QMin['template'][line[0]] = True
        elif line[0].lower() in strings:
            QMin['template'][line[0]] = line[1]
        elif line[0].lower() in floats:
            QMin['template'][line[0]] = float(line[1])

    # get external basis set block
    if 'basis_external' in QMin['template']:
        if os.path.isfile(QMin['template']['basis_external']):
            QMin['template']['basis_block'] = readfile(QMin['template']['basis_external'])
        else:
            print('"basis_external" key not readable: "%s"!' % QMin['template']['basis_external'])
            sys.exit(62)
    else:
        if 'basis' not in QMin['template']:
            print('Key "basis" missing in template file!' % (i))
            sys.exit(63)

    # check for completeness
    # necessary=['basis']
    # for i in necessary:
        # if not i in QMin['template']:
            # print('Key %s missing in template file!' % (i))
            # sys.exit(64)

    # now collect the casscf settings
    # jobs keyword
    jobs = [1 for i in range(QMin['maxmult'])]
    for line in temp:
        if line[0] == 'jobs':
            jobs = [int(i) for i in line[1:]]
    njobs = max(jobs)
    if any([i <= 0 for i in jobs]):
        print('Job ID numbers must be positive! Jobs: %s' % (jobs))
        sys.exit(65)
    if len(jobs) < QMin['maxmult']:
        print('No jobs for multiplicities larger than %i!' % (len(jobs)))
        sys.exit(66)
    QMin['jobs'] = jobs
    QMin['njobs'] = njobs
    # print('jobs:',QMin['jobs'])

    # get multmap
    multmap = {}
    for mult in itmult(QMin['states']):
        job = jobs[mult - 1]
        multmap[mult] = job
        if -job in multmap:
            multmap[-job].append(mult)
        else:
            multmap[-job] = [mult]
    for i in range(njobs):
        if not -(i + 1) in multmap:
            multmap[-(i + 1)] = []
    QMin['multmap'] = multmap
    # print('multmap:',multmap)

    # get the joblist
    joblist = set()
    for i in multmap:
        if i > 0:
            joblist.add(multmap[i])
    joblist = list(joblist)
    joblist.sort()
    QMin['joblist'] = joblist

    # orbital settings
    for line in temp:
        if line[0] == 'occ':
            QMin['template']['occ'] = [int(i) for i in line[1:]]
        if line[0] == 'closed':
            QMin['template']['closed'] = [int(i) for i in line[1:]]
    if 'occ' not in QMin['template']:
        print('No "occ" statement given!')
        sys.exit(67)
    if 'closed' not in QMin['template']:
        print('No "closed" statement given!')
        sys.exit(68)
    if len(QMin['template']['occ']) == 1:
        QMin['template']['occ'] = QMin['template']['occ'] * njobs
    if len(QMin['template']['closed']) == 1:
        QMin['template']['closed'] = QMin['template']['closed'] * njobs
    if len(QMin['template']['occ']) != njobs:
        print('Invalid "occ" specification! Give either 1 or %i values after "occ"!' % (njobs))
        sys.exit(69)
    if len(QMin['template']['closed']) != njobs:
        print('Invalid "closed" specification! Give either 1 or %i values after "closed"!' % (njobs))
        sys.exit(70)
    if any([i < 0 for i in QMin['template']['closed']]):
        print('Number of closed-shell orbitals must be positive! closed=%s' % (QMin['template']['closed']))
        sys.exit(71)
    if any([QMin['template']['closed'][i] >= QMin['template']['occ'][i] for i in range(njobs)]):
        print('Number of occupied orbitals must be larger than number of closed orbitals!')
        sys.exit(72)
    # print('closed:',QMin['template']['closed'])
    # print('occ:',QMin['template']['occ'])

    # wavefunction settings
    roots = {}
    rootpad = {}
    nelec = {}

    # roots
    i = 0
    for line in temp:
        if line[0] == 'roots':
            i += 1
            if i > QMin['njobs']:
                print('Too many "roots" statements (at least %i statements, but only %i jobs).' % (i, QMin['njobs']))
                sys.exit(73)
            f = [int(j) for j in line[1:]]
            if len(f) < QMin['maxmult']:
                f = f + [0] * (QMin['maxmult'] - len(f))
            # if len(f)>QMin['maxmult']:
                # f=f[:QMin['maxmult']]
            if any([j < 0 for j in f]):
                print('Number of roots must be positive! %s' % (f))
                sys.exit(74)
            for m in QMin['multmap'][-i]:
                if QMin['states'][m - 1] > f[m - 1]:
                    print('Not enough roots in job %i in multiplicity %i!' % (i, m))
                    sys.exit(75)
            roots[i] = f
    if len(roots) != njobs:
        print('Invalid number of "roots" statements! Please, give exactly %i "roots" statements.' % (njobs))
        sys.exit(76)
    # print('roots:',roots)

    # rootpad
    i = 0
    for line in temp:
        if line[0] == 'rootpad':
            i += 1
            f = [int(j) for j in line[1:]]
            if len(f) < len(roots[i]):
                f = f + [0] * (len(roots[i]) - len(f))
            # if len(f)>QMin['maxmult']:
                # f=f[:QMin['maxmult']]
            if any([j < 0 for j in f]):
                print('Values for "rootpad" must be positive! %s' % (f))
                sys.exit(77)
            rootpad[i] = f
    if len(rootpad) < njobs:
        for i in range(len(rootpad) + 1, njobs + 1):
            rootpad[i] = [0] * len(roots[i])
    # print('rootpad:',rootpad)

    # nelec
    i = 0
    for line in temp:
        if line[0] == 'nelec':
            i += 1
            f = [int(j) for j in line[1:]]
            if len(f) == 1:
                x = f[0]
                f = [x - ((x + m + 1) % 2 == 0) for m in range(len(roots[i]))]
            if len(f) != len(roots[i]):
                print('Invalid specification of "nelec"! Give either 1 or %i values after "nelec"!' % (len(roots[i])))
                sys.exit(78)
            if any([j <= 0 for j in f]):
                print('Number of electrons must be positive! %s' % (f))
                sys.exit(79)
            if any([(f[m] + m + 1) % 2 == 0 for m in range(len(roots[i]))]):
                print('Number of electrons incompatible with multiplicity! %s' % (f))
                sys.exit(80)
            nelec[i] = f
    if len(nelec) == 0:
        print('No "nelec" statements given!')
        sys.exit(81)
    if len(nelec) == 1:
        for i in range(njobs - 1):
            nelec[i + 2] = nelec[1]
    if len(nelec) != njobs:
        print('Invalid number of "nelec" statements! Give either 1 or %i "nelec" statements' % (njobs))
        sys.exit(82)
    # print('nelec:',nelec)

    QMin['template']['nelec'] = nelec
    QMin['template']['roots'] = roots
    QMin['template']['rootpad'] = rootpad

    # make the ionmap
    if 'ion' in QMin:
        ionmap = []
        for m1 in itmult(QMin['states']):
            job1 = QMin['multmap'][m1]
            el1 = QMin['template']['nelec'][job1][m1 - 1]
            for m2 in itmult(QMin['states']):
                if m1 >= m2:
                    continue
                job2 = QMin['multmap'][m2]
                el2 = QMin['template']['nelec'][job2][m2 - 1]
                # print m1,job1,el1,m2,job2,el2
                if abs(m1 - m2) == 1 and abs(el1 - el2) == 1:
                    ionmap.append((m1, job1, m2, job2))
        QMin['ionmap'] = ionmap

    # check the template
    # for i,job in enumerate(QMin['template']['casscf']):
        # headline=job['head'][0].rstrip()
        # if not 'maxit' in headline:
        # headline=headline+',maxit=40'
        # if not 'energy' in headline:
        # headline=headline+',energy=0.1e-7'
        # if not 'gradient' in headline:
        # headline=headline+',gradient=0.1e-6'
        # if not 'step' in headline:
        # headline=headline+',step=0.1e-2'
        # job['head'][0]=headline+'\n'

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
        QMin['backup'] = backupdir1



    # check for initial orbitals
    initorbs = {}
    if 'always_guess' in QMin:
        QMin['initorbs'] = {}
    elif 'init' in QMin or 'always_orb_init' in QMin:
        for job in QMin['joblist']:
            filename = os.path.join(QMin['pwd'], 'wf.%i.init' % (job))
            if os.path.isfile(filename):
                initorbs[job] = filename
            else:
                filename = os.path.join(QMin['pwd'], 'wf.init')
                if os.path.isfile(filename):
                    initorbs[job] = filename
        if 'always_orb_init' in QMin and len(initorbs) < njobs:
            print('Initial orbitals missing for some jobs!')
            sys.exit(83)
        QMin['initorbs'] = initorbs
    elif 'newstep' in QMin:
        for job in QMin['joblist']:
            filename = os.path.join(QMin['savedir'], 'wf.%i' % (job))
            if os.path.isfile(filename):
                initorbs[job] = filename + '.old'   # file will be moved to .old
            else:
                print('File %s missing in savedir!' % (filename))
                sys.exit(84)
        QMin['initorbs'] = initorbs
    elif 'samestep' in QMin:
        for job in QMin['joblist']:
            filename = os.path.join(QMin['savedir'], 'wf.%i' % (job))
            if os.path.isfile(filename):
                initorbs[job] = filename
            else:
                print('File %s missing in savedir!' % (filename))
                sys.exit(85)
        QMin['initorbs'] = initorbs
    elif 'restart' in QMin:
        for job in QMin['joblist']:
            filename = os.path.join(QMin['savedir'], 'wf.%i.old' % (job))
            if os.path.isfile(filename):
                initorbs[job] = filename
            else:
                print('File %s missing in savedir!' % (filename))
                sys.exit(86)
        QMin['initorbs'] = initorbs

    return QMin

# =============================================================================================== #
# =============================================================================================== #
# =========================================== Job Scheduling ==================================== #
# =============================================================================================== #
# =============================================================================================== #


def generate_joblist(QMin):

    joblist = []

    # add the initial integral computation
    QMin1 = deepcopy(QMin)
    QMin1['integrals'] = []
    remove = ['nacdr', 'grad', 'h', 'soc', 'dm', 'overlap', 'ion']
    for r in remove:
        QMin1 = removekey(QMin1, r)
    joblist.append({'integrals': QMin1})

    # add the master calculations for each job
    joblist.append({})
    for job in QMin['joblist']:
        QMin1 = deepcopy(QMin)
        QMin1['master'] = True
        QMin1['JOB'] = job
        remove = ['grad', 'nacdr', 'comment', 'ncpu']
        for r in remove:
            QMin1 = removekey(QMin1, r)
        joblist[-1]['master_%i' % (job)] = QMin1

    # add the gradient calculations
    joblist.append({})
    for grad in QMin['gradmap']:
        QMin1 = deepcopy(QMin)
        QMin1['gradmap'] = [grad]
        QMin1['JOB'] = QMin['multmap'][grad[0]]
        remove = ['h', 'soc', 'dm', 'overlap', 'ion', 'nacdr', 'always_guess', 'always_orb_init', 'comment', 'ncpu', 'init', 'docicas']
        for r in remove:
            QMin1 = removekey(QMin1, r)
        joblist[-1]['grad_%i_%i' % tuple(grad)] = QMin1

    # add the nac calculations
    for nac in QMin['nacmap']:
        QMin1 = deepcopy(QMin)
        QMin1['nacmap'] = [nac]
        QMin1['JOB'] = QMin['multmap'][nac[0]]
        remove = ['h', 'soc', 'dm', 'overlap', 'ion', 'grad', 'always_guess', 'always_orb_init', 'comment', 'ncpu', 'init', 'docicas']
        for r in remove:
            QMin1 = removekey(QMin1, r)
        joblist[-1]['nac_%i_%i_%i_%i' % tuple(nac)] = QMin1

    return joblist

# =============================================================================================== #
# =============================================================================================== #
# =========================================== Job Execution ===================================== #
# =============================================================================================== #
# =============================================================================================== #


def runjobs(joblist, QMin):

    if 'newstep' in QMin:
        moveOldFiles(QMin)

    # job='integrals'
    # QMin1=joblist[0][job]
    # job='master_1'
    # QMin1=joblist[1][job]
    # job='grad_1_1'
    # QMin1=joblist[2][job]
    # job='nac_1_1_1_2'
    # QMin1=joblist[2][job]
    # WORKDIR=os.path.join(QMin['scratchdir'],job)
    # run_calc(WORKDIR,QMin1)

    # sys.exit(87)

    print('>>>>>>>>>>>>> Starting the MOLPRO job execution')

    for jobset in joblist:
        if not jobset:
            continue
        errorcodes = {}
        pool = Pool(processes=QMin['ncpu'])
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
            sys.exit(88)





        for job in jobset:
            if 'master' in job:
                WORKDIR = os.path.join(QMin['scratchdir'], job)
                saveFiles(WORKDIR, jobset[job])
            if 'integrals' in job:
                if 'ion' in QMin or 'docicas' in QMin:
                    WORKDIR = os.path.join(QMin['scratchdir'], job)
                    saveAOovl(WORKDIR, jobset[job])

        print('')


    if PRINT:
        string = '  ' + '=' * 40 + '\n'
        string += '||' + ' ' * 40 + '||\n'
        string += '||' + ' ' * 10 + 'All Tasks completed!' + ' ' * 10 + '||\n'
        string += '||' + ' ' * 40 + '||\n'
        string += '  ' + '=' * 40 + '\n'

    return errorcodes

# ======================================================================= #


def run_calc(WORKDIR, QMin):
    try:
        for irun in range(5):
            Tasks = gettasks(QMin)
            setupWORKDIR(WORKDIR, Tasks, QMin)
            strip = 'integrals' not in QMin
            err = runMOLPRO(WORKDIR, QMin['molpro'], strip, QMin['molpro_arguments'])

            if 'grad' in QMin or 'nacdr' in QMin:
                out = readfile(os.path.join(WORKDIR, 'MOLPRO.out'))
                conv = getCPstatus(out)
                if conv == 0.:
                    break
                elif conv == -1:
                    return 97
                elif conv <= QMin['template']['gradaccumax']:
                    QMin['template']['gradaccudefault'] = 1.1 * conv
                elif conv > QMin['template']['gradaccumax']:
                    print('CRASHED:\t%s\tCP-MCSCF did not converge.' % (WORKDIR))
                    return 96

            elif 'master' in QMin:
                out = readfile(os.path.join(WORKDIR, 'MOLPRO.out'))
                errors = findCICASerrors(out)
                if 'excessgrad' in errors:
                    if 'pspace_mc' in QMin:
                        QMin['pspace_mc'] *= 4.
                    else:
                        QMin['pspace_mc'] = 1.
                if 'convrefci' in errors:
                    if 'pspace_ci' in QMin:
                        QMin['pspace_ci'] *= 4.
                    else:
                        QMin['pspace_ci'] = 1.
                if len(errors) == 0:
                    break

            else:
                break
        else:
            print('CRASHED:\t%s\tDid not converge.' % (WORKDIR))
            return 96
        return err
    except Exception as problem:
        print('*' * 50 + '\nException in run_calc(%s)!' % (WORKDIR))
        traceback.print_exc()
        print('*' * 50 + '\n')
        raise problem

# ======================================================================= #


def gettasks(QMin):

    states = QMin['states']
    nstates = QMin['nstates']
    nmstates = QMin['nmstates']

    tasks = []
    # pprint.pprint(QMin)

    # get job ID
    if 'integrals' in QMin:
        job = None
    else:
        job = QMin['JOB']

    # set header
    files = {1: 'integrals'}
    if 'integrals' not in QMin:
        if 'master' in QMin:
            files[2] = 'wf,new'
            if job in QMin['initorbs']:
                files[3] = 'wf.guess,old'
        elif 'grad' in QMin or 'nacdr' in QMin:
            files[2] = 'wf'
    tasks.append(['header', files])

    # AOoverlap
    if 'AOoverlap' in QMin:
        tasks.append(['coarseINT'])

    # integration
    if 'integrals' in QMin:
        tasks.append(['integrate'])
        tasks.append(['matrop', 'S'])
        return tasks

    # AO overlap and det output
    if 'ion' in QMin or 'docicas' in QMin or 'nooverlap' not in QMin:
        tasks.append(['gprint', 0.00000005])

    # master
    if 'master' in QMin:
        # MCSCF
        if 'pspace_mc' in QMin:
            tasks.append(['mcscf', job, QMin['pspace_mc']])
        else:
            tasks.append(['mcscf', job, -1.])
        tasks.append(['molden'])
        # CI
        mults = QMin['multmap'][-job]
        for m in mults:
            if 'pspace_ci' in QMin:
                tasks.append(['ci', job, m, states[m - 1], QMin['pspace_ci']])
            else:
                tasks.append(['ci', job, m, states[m - 1], -1.])
        # SOC
        if 'soc' in QMin:
            tupel = ['cihlsmat']
            for m in mults:
                tupel.append(m)
            tasks.append(tupel)

    # gradient calculations
    if 'grad' in QMin:
        for grad in QMin['gradmap']:
            tasks.append(['cpgrad', job, grad, QMin['template']['gradaccudefault']])
            tasks.append(['forcegrad', grad])

    # NAC calculations
    if 'nacdr' in QMin:
        for nac in QMin['nacmap']:
            tasks.append(['cpnac', job, nac, QMin['template']['gradaccudefault']])
            tasks.append(['forcenac', nac])

    if DEBUG:
        print(tasks)
    return tasks

# ======================================================================= #


def setupWORKDIR(WORKDIR, tasks, QMin):
    # mkdir the WORKDIR, or clean it if it exists, then copy all necessary files from pwd and savedir
    # then put the geom.xyz and MOLPRO.input files

    # setup the directory
    if os.path.exists(WORKDIR):
        if os.path.isfile(WORKDIR):
            print('%s exists and is a file!' % (WORKDIR))
            sys.exit(89)
        elif os.path.isdir(WORKDIR):
            if DEBUG:
                print('Remake\t%s' % WORKDIR)
            shutil.rmtree(WORKDIR)
            os.makedirs(WORKDIR)
    else:
        try:
            if DEBUG:
                print('Make\t%s' % WORKDIR)
            os.makedirs(WORKDIR)
        except OSError:
            print('Can not create %s\n' % (WORKDIR))
            sys.exit(90)

    # write MOLPRO.input
    inputstring = writeMOLPROinput(tasks, QMin)
    filename = os.path.join(WORKDIR, 'MOLPRO.input')
    writefile(filename, inputstring)
    if DEBUG:
        print(inputstring)
        print('MOLPRO input written to: %s' % (filename))

    # integral copying
    if 'integrals' not in QMin:
        fromfile = os.path.join(QMin['scratchdir'], 'integrals', 'integrals')
        tofile = os.path.join(WORKDIR, 'integrals')
        shutil.copy(fromfile, tofile)

    # wf file copying
    if 'integrals' in QMin:
        pass
    elif 'master' in QMin:
        job = QMin['JOB']
        if job in QMin['initorbs']:
            fromfile = QMin['initorbs'][job]
            tofile = os.path.join(WORKDIR, 'wf.guess')
            shutil.copy(fromfile, tofile)
    elif 'grad' in QMin or 'nacdr' in QMin:
        job = QMin['JOB']
        fromfile = os.path.join(QMin['scratchdir'], 'master_%i' % job, 'wf')
        tofile = os.path.join(WORKDIR, 'wf')
        shutil.copy(fromfile, tofile)

    return

# ======================================================================= #


def writeMOLPROinput(tasks, QMin):

    if 'integrals' in QMin:
        job = None
    else:
        job = QMin['JOB']

    # pprint.pprint(tasks)
    string = ''

    for task in tasks:

        # header ======================================================================================== #
        if task[0] == 'header':
            string += '***,MOLPRO input from SHARC-MOLPRO interface %s\n' % version
            string += 'memory,%i,k\n\n' % (QMin['memory'] * 125)
            for i in sorted(task[1]):
                string += 'file,%i,./%s\n' % (i, task[1][i])
            string += '\n'
            if QMin['template']['dkho'] > 0:
                string += 'dkho=%i\n' % (QMin['template']['dkho'])
            if 'basis_block' in QMin['template']:
                string += 'basis={\n%s\n}\n\n' % (''.join(QMin['template']['basis_block']))
            else:
                string += 'basis=%s\n\n' % (QMin['template']['basis'])
            string += 'nosym\nbohr\ngeometry={\n'
            for iatom, atom in enumerate(QMin['geo']):
                string += '%s%i %16.9f %16.9f %16.9f\n' % (atom[0], iatom + 1, atom[1], atom[2], atom[3])
            string += '}\n\n'

        # make Twoelectron integrals cheap when calculating AO overlaps
        elif task[0] == 'coarseINT':
            string += 'GTHRESH,THROVL=-1e6,TWOINT=1.d9,PREFAC=1.d9\nGDIRECT\n\n'

        # integrate ======================================================================================== #
        elif task[0] == 'integrate':
            string += 'int\n\n'

        # matrop ======================================================================================== #
        elif task[0] == 'matrop':
            string += '{matrop\nload,%s\nprint,%s\n}\n\n' % (task[1], task[1])

        # gprint ======================================================================================== #
        elif task[0] == 'gprint':
            string += 'gprint,orbitals,civectors;\n'
            string += 'gthresh,thrprint=0.,printci=%.10f;\n\n' % (task[1])

        # gprint ======================================================================================== #
        elif task[0] == 'molden':
            string += 'PUT,MOLDEN,orbs.molden;\n'

        # mcscf ======================================================================================== #
        elif task[0] == 'mcscf':
            # this mcscf runs are only for energies, not for gradients etc.
            maxit = 40
            energy = 1e-8
            gradient = 1e-7
            step = 1e-3
            string += '{casscf,maxit=%i,energy=%.10f,gradient=%.10f,step=%.10f\n' % (maxit, energy, gradient, step)
            string += 'frozen,0\nclosed,%i\nocc,%i\n' % (QMin['template']['closed'][task[1] - 1], QMin['template']['occ'][task[1] - 1])
            if task[1] in QMin['initorbs']:
                string += 'start,2140.3\n'
            string += 'orbital,2140.2\n'
            if task[2] > 0.:
                string += 'pspace,%.2f\n' % (task[2])
            for m in range(len(QMin['template']['roots'][task[1]])):
                nelec = QMin['template']['nelec'][task[1]][m]
                roots = QMin['template']['roots'][task[1]][m]
                rootpad = QMin['template']['rootpad'][task[1]][m]
                if roots + rootpad == 0:
                    continue
                string += 'wf,%i,1,%i\nstate,%i\nweight' % (nelec, m, roots + rootpad)
                for i in range(roots):
                    string += ',1'
                for i in range(rootpad):
                    string += ',0'
                string += '\n'
            string += '};\n\n'

        # ci ======================================================================================== #
        elif task[0] == 'ci':
            nelec = QMin['template']['nelec'][task[1]][task[2] - 1]
            string += '{ci\nmaxiter,250,1000\norbital,2140.2\nsave,%i.2\nnoexc\n' % (6000 + task[2])
            if task[4] > 0.:
                string += 'pspace,%.2f\n' % (task[4])
            string += 'core,%i\nwf,%i,1,%i\nstate,%i\n};\n\n' % (
                QMin['template']['closed'][task[1] - 1],
                nelec,
                task[2] - 1,
                task[3])

        # spin-orbit ======================================================================================== #
        elif task[0] == 'cihlsmat':
            string += '{ci\nhlsmat,amfi'
            for i in range(len(task) - 1):
                string += ',%i.2' % (6000 + task[i + 1])
            string += '\nprint,hls=1\n};\n\n'

        # mcscf with cp ======================================================================================== #
        elif task[0] == 'cpgrad' or task[0] == 'cpnac':
            # this mcscf runs are for gradients
            maxit = 40
            energy = 1e-8
            gradient = 1e-7
            step = 1e-3
            string += '{casscf,maxit=%i,energy=%.10f,gradient=%.10f,step=%.10f\n' % (maxit, energy, gradient, step)
            string += 'frozen,0\nclosed,%i\nocc,%i\n' % (QMin['template']['closed'][task[1] - 1], QMin['template']['occ'][task[1] - 1])
            if task[1] in QMin['initorbs']:
                string += 'start,2140.2\n'
            string += 'orbital,2200.2\n'
            # string+='dont,orbital\n'         # TODO: why can't MOLPRO restart properly?
            for m in range(QMin['maxmult']):
                nelec = QMin['template']['nelec'][task[1]][m]
                roots = QMin['template']['roots'][task[1]][m]
                rootpad = QMin['template']['rootpad'][task[1]][m]
                if roots + rootpad == 0:
                    continue
                string += 'wf,%i,1,%i\nstate,%i\nweight' % (nelec, m, roots + rootpad)
                for i in range(roots):
                    string += ',1'
                for i in range(rootpad):
                    string += ',0'
                string += '\n'
            string += 'print,micro\n'
            if task[0] == 'cpgrad':
                string += 'cpmcscf,grad,state=%i.1,ms2=%i,record=%i.2,accu=%18.15f\n' % (
                    task[2][1],
                    task[2][0] - 1,
                    5000 + task[2][0] * 100 + task[2][1],
                    task[3])
            elif task[0] == 'cpnac':
                string += 'cpmcscf,nacm,state1=%i.1,state2=%i.1,ms2=%i,record=%i.2,accu=%18.15f\n' % (
                    task[2][1],
                    task[2][3],
                    task[2][0] - 1,
                    5020 + task[2][0] * 100 + 10 * task[2][1] + task[2][3],
                    task[3])
            string += '};\n\n'

        # force for grad =============================================================================== #
        elif task[0] == 'forcegrad':
            string += '{force;\nsamc,%i.2\n};\n\n' % (5000 + task[1][0] * 100 + task[1][1])

        # force for nac ================================================================================== #
        elif task[0] == 'forcenac':
            string += '{force;\nsamc,%i.2\n};\n\n' % (5020 + task[1][0] * 100 + 10 * task[1][1] + task[1][3])

        else:  # ========================================================================================== #
            print('Unknown task %s found in writeMOLPROinput!' % task)
            sys.exit(91)

    return string





    ## move old wavefunction files ====================================================================================== #
    # if not tasks[0][0]=='initstep' and not tasks[0][0]=='samestep' and not tasks[0][0]=='restartstep':
    # fromfile=os.path.join(QMin['savedir'],'wf.last')
    #tofile  =os.path.join(QMin['savedir'],'wf.prelast')
    # try:
    # os.rename(fromfile,tofile)
    # except OSError:
    # pass
    # fromfile=os.path.join(QMin['savedir'],'wf.current')
    #tofile  =os.path.join(QMin['savedir'],'wf.last')
    # try:
    # os.rename(fromfile,tofile)
    # except OSError:
    # pass
    # rmfiles=['integrals','wf.current','wf.last']
    # for f in rmfiles:
    # exfile=os.path.join(QMin['scratchdir'],f)
    # if os.path.exists(exfile):
    # os.remove(exfile)

    # if tasks[0][0]=='samestep' or tasks[0][0]=='restartstep':
    # inp.write('file,2,./wf.current\n')
    # exfile=os.path.join(QMin['scratchdir'],'wf.current')
    # if not os.path.exists(exfile):
    # exfile=os.path.join(QMin['savedir'],'wf.current')
    # if not os.path.exists(exfile):
    # print('No "wf.current" found in scratchdir or savedir!')
    #sys.exit (92)
    # else:
    # tofile=os.path.join(QMin['scratchdir'],'wf.current')
    # shutil.copy(exfile,tofile)
    # else:
    # inp.write('file,2,./wf.current,new\n')


    # if tasks[0][0]=='initstep':
    # inp.write('file,3,./wf.init\n\n')
    # exfile=os.path.join(QMin['pwd'],'wf.init')
    # if os.path.exists(exfile):
    # tofile=os.path.join(QMin['scratchdir'],'wf.init')
    # shutil.copy(exfile,tofile)
    # else:
    # inp.write('file,3,./wf.last\n\n')
    # exfile=os.path.join(QMin['scratchdir'],'wf.last')
    # if os.path.exists(exfile):
    # exfile=os.path.join(QMin['savedir'],'wf.last')
    # if not os.path.exists(exfile):
    # print('No "wf.last" found in scratchdir or savedir!')
    #sys.exit (93)
    # else:
    # tofile=os.path.join(QMin['scratchdir'],'wf.last')
    # shutil.copy(exfile,tofile)

    ## mcscf:pspace: create the same casscf block as above, but include a pspace statement (convergence helper...) ====== #
    # if task[0]=='mcscf:pspace':
    # string+='{'+''.join(QMin['template']['casscf'][task[1]-1]['head'])
    # if 'init' in QMin:
    # string+='start,2140.3\n'
    # elif 'newstep' in QMin:
    #string+='start,%i.3\n' % (2140+task[1]-1)
    #string+='orbital,%i.2\n' % (2140+task[1]-1)
    #string+='pspace,%.2f\n' % (task[2])
    # string+=mult_to_WFblock(QMin['template']['casscf'][task[1]-1]['mult'])
    # string+='};\n\n'

    ## same as above, but with nstati statement (convergence helper) ==================================================== #
    # elif task[0]=='ci:nstati':
    # print('ci:nstati is not yet implemented!')
    # sys.exit(94)

    ## same as above, but with pspace statement (convergence helper) ==================================================== #
    # elif task[0]=='ci:pspace':
    # print header
    # string+='{ci\nmaxiter,250,1000\norbital,%i.2\nsave,%i.2\nnoexc\npspace,%i\ncore,%i\n' % (
    # 2140+task[1]-1,
    # 6000+task[2],
    # task[4],
    # QMin['template']['casscf'][task[1]-1]['closed'])
    # nelec=QMin['template']['casscf'][task[1]-1]['mult'][task[2]]['nelec']
    #string+='wf,%i,%i,%i\nstate,%i\n};\n\n' % (nelec,1,task[2]-1,task[3])


# ======================================================================= #
def runMOLPRO(WORKDIR, MOLPRO, strip=False, arguments=''):
    prevdir = os.getcwd()
    os.chdir(WORKDIR)
    os.environ['WorkDir'] = WORKDIR
    string = os.path.join(MOLPRO, 'molpro') + ' '
    if arguments:
        string += ' ' + arguments + ' '
    string += '-W./ -I./ -d./ MOLPRO.input'
    stdoutfile = open(os.path.join(WORKDIR, 'MOLPRO.out'), 'w')
    stderrfile = open(os.path.join(WORKDIR, 'MOLPRO.err'), 'w')
    if PRINT or DEBUG:
        starttime = datetime.datetime.now()
        sys.stdout.write('START:\t%s\t%s\t"%s"\n' % (WORKDIR, starttime, string))
        sys.stdout.flush()
    try:
        runerror = sp.call(string, shell=True, stdout=stdoutfile, stderr=stderrfile)
    except OSError:
        print('Call have had some serious problems:', OSError)
        sys.exit(95)
    stdoutfile.close()
    stderrfile.close()
    if PRINT or DEBUG:
        endtime = datetime.datetime.now()
        sys.stdout.write('FINISH:\t%s\t%s\tRuntime: %s\tError Code: %i\n' % (WORKDIR, endtime, endtime - starttime, runerror))
        sys.stdout.flush()
    os.chdir(prevdir)
    if strip and not DEBUG:
        stripWORKDIR(WORKDIR)
    return runerror

# ======================================================================= #


def stripWORKDIR(WORKDIR):
    ls = os.listdir(WORKDIR)
    keep = ['MOLPRO.out$', 'wf', 'orbs.molden']
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


def cleandir(directory):
    # ''''''
    # if PRINT:
    # print('===> Removing SCRATCHDIR=%s\n' % (SCRATCHDIR))
    # for data in os.listdir(SCRATCHDIR):
    # path=SCRATCHDIR+'/'+data
    # try:
    # if DEBUG or PRINT:
    # print('rm %s' % (path))
    # os.remove(path)
    # except OSError:
    # print('Could not remove file from SCRATCHDIR: %s' % (path))
    # try:
    # if DEBUG or PRINT:
    # print('rm %s\n\n' % (SCRATCHDIR))
    # os.rmdir(SCRATCHDIR)
    # except OSError:
    # print('Could not remove SCRATCHDIR=%s' % (SCRATCHDIR))
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


def moveOldFiles(QMin):
    # move geom file if it exists
    if 'nooverlap' not in QMin:
        fromfile = os.path.join(QMin['savedir'], 'geom')
        if not os.path.isfile(fromfile):
            print('File %s not found, cannot move to OLD!' % (fromfile))
            sys.exit(96)
        tofile = os.path.join(QMin['savedir'], 'geom.old')
        if DEBUG:
            print('Copy:\t%s\n\t==>\n\t%s' % (fromfile, tofile))
        shutil.copy(fromfile, tofile)
    # moves all relevant files in the savedir to old files
    basenames = ['wf']
    if 'nooverlap' not in QMin:
        basenames.append('mo')
    for job in QMin['joblist']:
        for base in basenames:
            fromfile = os.path.join(QMin['savedir'], '%s.%i' % (base, job))
            if not os.path.isfile(fromfile):
                print('File %s not found, cannot move to OLD!' % (fromfile))
                sys.exit(97)
            tofile = os.path.join(QMin['savedir'], '%s.%i.old' % (base, job))
            if DEBUG:
                print('Copy:\t%s\n\t==>\n\t%s' % (fromfile, tofile))
            shutil.copy(fromfile, tofile)
    # moves all relevant files in the savedir to old files
    if 'nooverlap' not in QMin:
        basenames = ['det_ci']
        for job in itmult(QMin['states']):
            for base in basenames:
                fromfile = os.path.join(QMin['savedir'], '%s.%i' % (base, job))
                if not os.path.isfile(fromfile):
                    print('File %s not found, cannot move to OLD!' % (fromfile))
                    sys.exit(98)
                tofile = os.path.join(QMin['savedir'], '%s.%i.old' % (base, job))
                if DEBUG:
                    print('Copy:\t%s\n\t==>\n\t%s' % (fromfile, tofile))
                shutil.copy(fromfile, tofile)
    # also remove aoovl files if present
    delete = ['aoovl', 'aoovl_double']
    for f in delete:
        rmfile = os.path.join(QMin['savedir'], f)
        if os.path.isfile(rmfile):
            os.remove(rmfile)

# ======================================================================= #


def saveFiles(WORKDIR, QMin):
    # see https://www.ibm.com/support/pages/apar/IJ29942
    shutil._USE_CP_SENDFILE = False
    # copy the wf files from master directories
    job = QMin['JOB']
    fromfile = os.path.join(WORKDIR, 'wf')
    tofile = os.path.join(QMin['savedir'], 'wf.%i' % (job))
    shutil.copy(fromfile, tofile)

    # if necessary, extract the MOs and write them to savedir
    if 'ion' in QMin or 'docicas' in QMin or 'nooverlap' not in QMin:
        out = readfile(os.path.join(WORKDIR, 'MOLPRO.out'))
        string = get_MO_from_out(out)
        mofile = os.path.join(QMin['savedir'], 'mo.%i' % job)
        writefile(mofile, string)

    # if necessary, extract the CASSCF coefficients and write them to savedir
    # TODO: actually, the CASSCF det files are not needed in savedir, could be placed in a keepdir
    if 'docicas' in QMin:
        out = readfile(os.path.join(WORKDIR, 'MOLPRO.out'))
        for im, m in enumerate(QMin['multmap'][-job]):
            roots = QMin['template']['roots'][job]
            n = 0
            for i in roots:
                if i > 0:
                    n += 1
            if n == 1:
                string = get_CASdet_from_out(out, 0, QMin['states'][m - 1])
            else:
                n = 0
                for i, j in enumerate(roots):
                    if j > 0:
                        n += 1
                    if i + 1 == m:
                        break
                string = get_CASdet_from_out(out, m, QMin['states'][m - 1])
            detfile = os.path.join(QMin['savedir'], 'det_cas.%i' % m)
            writefile(detfile, string)

    # if necessary, extract the MRCI coefficients, convert to determinants, and write them to savedir
    if 'docicas' in QMin or 'ion' in QMin or 'nooverlap' not in QMin:
        out = readfile(os.path.join(WORKDIR, 'MOLPRO.out'))
        for im, m in enumerate(QMin['multmap'][-job]):
            string = get_CIdet_from_out(out, m)
            detfile = os.path.join(QMin['savedir'], 'det_ci.%i' % m)
            writefile(detfile, string)

    # if necessary, save the current geometry in savedir
    if 'nooverlap' not in QMin:
        string = ''
        for iatom, atom in enumerate(QMin['geo']):
            string += '%s %16.9f %16.9f %16.9f\n' % (atom[0], atom[1], atom[2], atom[3])
        geofile = os.path.join(QMin['savedir'], 'geom')
        writefile(geofile, string)

    # save molden files
    if 'molden' in QMin:
        fromfile = os.path.join(WORKDIR, 'orbs.molden')
        tofile = os.path.join(QMin['savedir'], 'orbs.%i.molden' % (job))
        shutil.copy(fromfile, tofile)

# ======================================================================= #


def saveAOovl(WORKDIR, QMin):
    # saveAOovl files from integrals directory

    outfile = os.path.join(WORKDIR, 'MOLPRO.out')
    out = readfile(outfile)

    # get number of AOs
    for line in out:
        if 'NUMBER OF CONTRACTIONS:' in line:
            nao = int(line.split()[3])
            break

    # find AO overlap matrix
    for iline, line in enumerate(out):
        if ' SYMMETRY BLOCK 1.1' in line:
            break
    else:
        print('Did not find AO overlap matrix!')
        sys.exit(99)

    # detect format
    line = out[iline + 1]
    if 'E' in line:
        formatting = 2
        q = 24
    elif len(line.split()) == 0:
        formatting = 1
    else:
        formatting = 0
        q = len(line.split())

    # get matrix
    AOovl = []
    for irow in range(nao):
        AOovl.append([])
        if formatting == 2:
            iline += 1
            line = out[iline]
            i = 0
            s = []
            while True:
                x = line[1 + i * q:1 + (i + 1) * q]
                s.append(x)
                i += 1
                if 1 + (i + 1) * q > len(line):
                    break
            for y in s:
                AOovl[-1].append(float(y))
        elif formatting == 1:
            iline += 1
            for x in range((nao - 1) // 10 + 1):
                iline += 1
                line = out[iline]
                s = line.split()
                for y in s:
                    AOovl[-1].append(float(y))
        elif formatting == 0:
            for x in range((nao - 1) // q + 1):
                iline += 1
                line = out[iline]
                s = line.split()
                for y in s:
                    AOovl[-1].append(float(y))

    # format string
    string = '%i %i\n' % (nao, nao)
    for irow in range(nao):
        for icol in range(nao):
            string += '%11.8f ' % (AOovl[irow][icol])
        string += '\n'

    # write to file
    aofile = os.path.join(QMin['savedir'], 'aoovl')
    writefile(aofile, string)

    return

# ======================================================================= #


def getCPstatus(out):
    for iline, line in enumerate(out):
        if 'for CP-MCSCF' in line:
            break
    else:
        return -1.
    iline += 1
    conv = 1e6
    while True:
        iline += 1
        line = out[iline]
        if 'Convergence reached' in line:
            return 0.
    return conv

# ======================================================================= #


def findCICASerrors(out):
    errors = set()
    for line in out:
        if 'EXCESSIVE GRADIENT IN CI' in line:
            errors.add('excessgrad')
        if 'NO CONVERGENCE IN REFERENCE CI' in line or 'Sometimes it helps to redefine P space' in line:
            errors.add('convrefci')
    return errors

# =============================================================================================== #
# =============================================================================================== #
# =======================================  MO and Det extraction  =============================== #
# =============================================================================================== #
# =============================================================================================== #


def fortraneformat(number, prec=12, exp_digits=2):
    neg = number < 0.
    if neg:
        s = "%.*e" % (prec, -number)
    else:
        s = "%.*e" % (prec, number)
    mantissa, exp = s.split('e')
    mantissa = mantissa.replace('.', '').replace('-', '')[:prec]
    if number != 0:
        exp = int(exp) + 1
    if neg:
        prefix = '-.'
    else:
        prefix = '0.'
    return "%s%sE%+0*d" % (prefix, mantissa, exp_digits + 1, int(exp))

# ======================================================================= #


def get_MO_from_out(out):
    # extracts the first occurence of orbitals and formats them into lumorb format for savedir
    # find number of orbitals
    found = 0
    for line in out:
        if 'Number of closed-shell orbitals:' in line:
            nclosed = int(line.replace('(',' ').replace(')',' ').split()[-2])
            found += 1
        elif 'Number of active  orbitals:' in line:
            nact = int(line.replace('(',' ').replace(')',' ').split()[-2])
            found += 1
        elif 'Number of external orbitals:' in line:
            nex = int(line.replace('(',' ').replace(')',' ').split()[-2])
            found += 1
            if found == 2:
                # this means that no closed-shell orbitals are present in output
                found += 1
                nclosed = 0
        elif '1PROGRAM * CI ' in line:
            found -= 1000
        if found == 3:
            break
    else:
        print('Did not find number of orbital specifications in output!')
        sys.exit(100)
    norb = nclosed + nact + nex
    # print norb

    # find the orbitals
    for iline, line in enumerate(out):
        if 'NATURAL ORBITALS' in line:
            break
    iline += 7 + (norb) // 10

    # extract the closed and active orbitals (the others are never occupied anyways)
    orbitals = []
    for imo in range(nclosed + nact):
        icol = 3
        s = out[iline].split()
        orb = []
        for iao in range(norb):
            orb.append(float(s[icol]))
            icol += 1
            if (iao + 1) % 10 == 0 and not (iao + 1) == norb:
                iline += 1
                s = out[iline].split()
                icol = 0
        orbitals.append(orb)
        iline += 2

    # format the lumorb string
    string = '#INPORB 1.1\n#INFO\n* MOs from SHARC-MOLPRO interface in lumorb format\n0 1 0\n%i\n%i\n#ORB\n' % (norb, nclosed + nact)
    for imo in range(nclosed + nact):
        string += '* ORBITAL 1 %i\n' % (imo + 1)
        for iao in range(norb):
            string += fortraneformat(orbitals[imo][iao])
            if (iao + 1) % 4 == 0:
                string += '\n'
        if not (iao + 1) % 4 == 0:
            string += '\n'

    return string

# ======================================================================= #


def get_CASdet_from_out(out, imult, nstates):
    # extracts the first occurence of determinants from CASSCF and formats them for savedir
    # find number of orbitals
    found = 0
    for line in out:
        if 'Number of closed-shell orbitals:' in line:
            nclosed = int(line.replace('(',' ').replace(')',' ').split()[-2])
            found += 1
        elif 'Number of active  orbitals:' in line:
            nact = int(line.replace('(',' ').replace(')',' ').split()[-2])
            found += 1
            if found == 1:
                found += 1
                nclosed = 0
        elif '1PROGRAM * CI ' in line:
            found -= 1000
        if found == 2:
            break
    else:
        print('Did not find number of orbital specifications in output!')
        sys.exit(101)
    norb = nclosed + nact

    # find the determinants
    for iline, line in enumerate(out):
        if imult == 0:
            string = 'CI Coefficients'
            multi = ' '
        else:
            string = 'CI Coefficients of symmetry'
            multi = IToMult[imult]
        if string in line and multi in line: 
            break
    iline += 2

    # get the determinants
    ci_vectors = {'ndocc': nclosed, 'nvirt': 0}
    while True:
        iline += 1
        line = out[iline]
        s = line.split()
        if len(s) == 0:
            break
        if line[1:2] != ' ':
            # new det entry
            key = s[0].replace('2', '3').replace('a', '1').replace('b', '2')
            key = tuple([int(i) for i in key])
            coeff = [float(i) for i in s[1:]]
            ci_vectors[key] = coeff
        else:
            # more coefficients for previous entry
            coeff = [float(i) for i in s]
            ci_vectors[key].extend(coeff)

    # truncate to nstates
    for key in ci_vectors:
        if key == 'ndocc' or key == 'nvirt':
            continue
        ci_vectors[key] = ci_vectors[key][0:nstates]

    return format_ci_vectors(ci_vectors)

# ======================================================================= #


def get_CIdet_from_out(out, imult):
    # get CI coefficients for states of imult multiplicity
    # find CI block
    for iline, line in enumerate(out):
        if 'Reference symmetry:' in line and IToMult[imult] in line:
            break

    # find number of orbitals
    found = 0
    while True:
        iline += 1
        if iline >= len(out):
            print('Did not find number of orbital specifications in output!')
            sys.exit(102)
        line = out[iline]
        if 'Number of core orbitals:' in line:
            nclosed = int(line.split()[-4])
            found += 1
        elif 'Number of active  orbitals:' in line:
            nact = int(line.split()[-4])
            found += 1
            if found == 1:
                found += 1
                nclosed = 0
        elif '1PROGRAM * CI ' in line:
            found += -1000
        if found == 2:
            break
    norb = nclosed + nact

    # find coefficients
    while True:
        iline += 1
        if iline >= len(out):
            print('Did not find number of orbital specifications in output!')
            sys.exit(103)
        line = out[iline]
        if ' Reference coefficients greater than' in line:
            break
        elif '1PROGRAM * CI ' in line:
            print('Did not find number of orbital specifications in output!')
            sys.exit(104)
    iline += 1

    # get csfs
    csf_vectors = {}
    while True:
        iline += 1
        line = out[iline]
        s = line.split()
        if len(s) == 0:
            break
        if line[1:2] != ' ':
            # new det entry
            key = s[0].replace('2', '3').replace('/', '1').replace('\\', '2')
            key = tuple([int(i) for i in key])
            coeff = [float(i) for i in s[1:]]
            csf_vectors[key] = coeff
        else:
            # more coefficients for previous entry
            coeff = [float(i) for i in s]
            csf_vectors[key].extend(coeff)

    # convert to determinants
    ci_vectors = {'ndocc': nclosed, 'nvirt': 0}
    for csf in csf_vectors:
        dets = decompose_csf((imult - 1) / 2., list(csf))
        coeff = csf_vectors[csf]
        for det in dets:
            c = [dets[det] * i for i in coeff]
            if det in ci_vectors:
                for istate in range(len(coeff)):
                    ci_vectors[det][istate] += c[istate]
            else:
                ci_vectors[det] = c

    return format_ci_vectors(ci_vectors)

# ======================================================================= #


def decompose_csf(ms2, step):
    # ms2 is M_S value
    # step is step vector for CSF (e.g. 3333012021000)

    def powmin1(x):
        a = [1, -1]
        return a[x % 2]

    # calculate key numbers
    nopen = sum([i == 1 or i == 2 for i in step])
    nalpha = int(nopen / 2. + ms2)
    norb = len(step)

    # make reference determinant
    refdet = deepcopy(step)
    for i in range(len(refdet)):
        if refdet[i] == 1:
            refdet[i] = 2

    # get the b vector and the set of open shell orbitals
    bval = []
    openorbs = []
    b = 0
    for i in range(norb):
        if step[i] == 1:
            b += 1
        elif step[i] == 2:
            b -= 1
        bval.append(b)
        if refdet[i] == 2:
            openorbs.append(i)

    # loop over the possible determinants
    dets = {}
    # get all possible combinations of nalpha orbitals from the openorbs set
    for localpha in itertools.combinations(openorbs, nalpha):
        # make determinant string
        det = deepcopy(refdet)
        for i in localpha:
            det[i] = 1

        # get coefficient
        coeff = 1.
        sign = +1
        m2 = 0
        for k in range(norb):
            if step[k] == 0:
                pass
                # sign*=powmin1(bval[k])
            elif step[k] == 1:
                m2 += powmin1(det[k] + 1)
                num = bval[k] + powmin1(det[k] + 1) * m2
                if num == 0.:
                    break
                denom = 2. * bval[k]
                sign *= powmin1(bval[k])
                coeff *= 1. * num / denom
            elif step[k] == 2:
                m2 += powmin1(det[k] - 1)
                num = bval[k] + 2 + powmin1(det[k]) * m2
                if num == 0.:
                    break
                denom = 2. * (bval[k] + 2)
                sign *= powmin1(bval[k] - det[k] + 1)
                coeff *= 1. * num / denom
            elif step[k] == 3:
                # sign*=powmin1(bval[k])
                num = 1.

        # add determinant to dict if coefficient non-zero
        if num != 0.:
            dets[tuple(det)] = -1. * sign * math.sqrt(coeff)

    return dets

# ======================================================================= #


def format_ci_vectors(ci_vectors):
    # get nstates, norb and ndets
    for key in ci_vectors:
        if key != 'ndocc' and key != 'nvirt':
            nstates = len(ci_vectors[key])
            norb = len(key)
            break
    ndets = len(ci_vectors) - 2
    ndocc = ci_vectors['ndocc']
    nvirt = ci_vectors['nvirt']

    # sort determinant strings
    dets = []
    for key in ci_vectors:
        if key != 'ndocc' and key != 'nvirt':
            dets.append(key)
    dets.sort(reverse=True)

    string = '%i %i %i\n' % (nstates, norb + ndocc + nvirt, ndets)
    for det in dets:
        for i in range(ndocc):
            string += 'd'
        for o in det:
            if o == 0:
                string += 'e'
            elif o == 1:
                string += 'a'
            elif o == 2:
                string += 'b'
            elif o == 3:
                string += 'd'
        for i in range(nvirt):
            string += 'e'
        for c in ci_vectors[det]:
            string += ' %16.12f ' % c
        string += '\n'
    return string

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
            files = {'aoovl': 'aoovl',
                     'det.a': 'det_ci.%i' % ionpair[0],
                     'det.b': 'det_ci.%i' % ionpair[2],
                     'mo.a': 'mo.%i' % ionpair[1],
                     'mo.b': 'mo.%i' % ionpair[3]}
            setupWORKDIR_WF(WORKDIR, QMin, files)
            errorcodes['Dyson_%i_%i_%i_%i' % ionpair] = runWFOVERLAP(WORKDIR, QMin['wfoverlap'], memory=QMin['memory'], ncpu=QMin['ncpu'])

    # do CI-CAS overlaps
    if 'docicas' in QMin:
        for m in itmult(QMin['states']):
            job = QMin['multmap'][m]
            WORKDIR = os.path.join(QMin['scratchdir'], 'CI_CAS_%i_%i' % (m, job))
            files = {'aoovl': 'aoovl',
                     'det.a': 'det_ci.%i' % m,
                     'det.b': 'det_cas.%i' % m,
                     'mo.a': 'mo.%i' % job,
                     'mo.b': 'mo.%i' % job}
            setupWORKDIR_WF(WORKDIR, QMin, files)
            errorcodes['CI_CAS_%i_%i' % (m, job)] = runWFOVERLAP(WORKDIR, QMin['wfoverlap'], memory=QMin['memory'], ncpu=QMin['ncpu'])

    # do overlap calculations
    if 'overlap' in QMin:
        get_Double_AOovl(QMin)
        for m in itmult(QMin['states']):
            job = QMin['multmap'][m]
            WORKDIR = os.path.join(QMin['scratchdir'], 'WFOVL_%i_%i' % (m, job))
            files = {'aoovl': 'aoovl_double',
                     'det.b': 'det_ci.%i' % m,
                     'det.a': 'det_ci.%i.old' % m,
                     'mo.b': 'mo.%i' % job,
                     'mo.a': 'mo.%i.old' % job}
            setupWORKDIR_WF(WORKDIR, QMin, files)
            errorcodes['WFOVL_%i_%i' % (m, job)] = runWFOVERLAP(WORKDIR, QMin['wfoverlap'], memory=QMin['memory'], ncpu=QMin['ncpu'])

    # Error code handling
    j = 0
    string = 'Error Codes:\n'
    for i in errorcodes:
        if 'Dyson' in i or 'CI_CAS' in i or 'WFOVL' in i:
            string += '\t%s\t%i' % (i + ' ' * (10 - len(i)), errorcodes[i])
            j += 1
            if j == 4:
                j = 0
                string += '\n'
    print(string)
    if any((i != 0 for i in errorcodes.values())):
        print('Some subprocesses did not finish successfully!')
        sys.exit(105)

    print('')

    if PRINT:
        string = '  ' + '=' * 40 + '\n'
        string += '||' + ' ' * 40 + '||\n'
        string += '||' + ' ' * 10 + 'All Tasks completed!' + ' ' * 10 + '||\n'
        string += '||' + ' ' * 40 + '||\n'
        string += '  ' + '=' * 40 + '\n'

    return errorcodes

# ======================================================================= #


def runWFOVERLAP(WORKDIR, WFOVERLAP, memory=100, ncpu=1):
    prevdir = os.getcwd()
    os.chdir(WORKDIR)
    os.environ['WorkDir'] = WORKDIR
    string = WFOVERLAP + ' -m %i' % (memory) + ' -f wfovl.inp'
    stdoutfile = open(os.path.join(WORKDIR, 'wfovl.out'), 'w')
    stderrfile = open(os.path.join(WORKDIR, 'wfovl.err'), 'w')
    os.environ['OMP_NUM_THREADS'] = str(ncpu)
    if PRINT or DEBUG:
        starttime = datetime.datetime.now()
        sys.stdout.write('START:\t%s\t%s\t"%s"\n' % (WORKDIR, starttime, string))
        sys.stdout.flush()
    try:
        runerror = sp.call(string, shell=True, stdout=stdoutfile, stderr=stderrfile)
    except OSError:
        print('Call have had some serious problems:', OSError)
        sys.exit(106)
    stdoutfile.close()
    stderrfile.close()
    if PRINT or DEBUG:
        endtime = datetime.datetime.now()
        sys.stdout.write('FINISH:\t%s\t%s\tRuntime: %s\tError Code: %i\n' % (WORKDIR, endtime, endtime - starttime, runerror))
        sys.stdout.flush()
    os.chdir(prevdir)
    # if strip and not DEBUG:
    # stripWORKDIR(WORKDIR)
    return runerror

# ======================================================================= #


def setupWORKDIR_WF(WORKDIR, QMin, files):
    # mkdir the WORKDIR, or clean it if it exists, then copy all necessary files from pwd and savedir
    # then put the geom.xyz and MOLPRO.input files

    # setup the directory
    if os.path.exists(WORKDIR):
        if os.path.isfile(WORKDIR):
            print('%s exists and is a file!' % (WORKDIR))
            sys.exit(107)
        elif os.path.isdir(WORKDIR):
            if DEBUG:
                print('Remake\t%s' % WORKDIR)
            shutil.rmtree(WORKDIR)
            os.makedirs(WORKDIR)
    else:
        try:
            if DEBUG:
                print('Make\t%s' % WORKDIR)
            os.makedirs(WORKDIR)
        except OSError:
            print('Can not create %s\n' % (WORKDIR))
            sys.exit(108)

    # write wfovl.inp
    inputstring = '''mix_aoovl=aoovl
a_mo=mo.a
b_mo=mo.b
a_det=det.a
b_det=det.b
a_mo_read=1
b_mo_read=1
ao_read=0
'''
    inputstring += 'ncore=%i\n' % (QMin['ncore'])
    if 'ion' in QMin:
        if QMin['ndocc'] < 0:
            ndocc = min(QMin['template']['closed'])
        else:
            ndocc = QMin['ndocc']
        inputstring += 'ndocc=%i\n' % (ndocc)
    filename = os.path.join(WORKDIR, 'wfovl.inp')
    writefile(filename, inputstring)
    if DEBUG:
        print(inputstring)
        print('MOLPRO input written to: %s' % (filename))

    # link input files from save
    linkfiles = ['aoovl', 'det.a', 'det.b', 'mo.a', 'mo.b']
    for f in linkfiles:
        fromfile = os.path.join(QMin['savedir'], files[f])
        tofile = os.path.join(WORKDIR, f)
        link(fromfile, tofile)

    return

# ======================================================================= #


def get_Double_AOovl(QMin):
    # get old geometry
    old = readfile(os.path.join(QMin['savedir'], 'geom.old'))
    oldgeom = []
    shift = 1.1e-5
    for line in old:
        s = line.split()
        oldgeom.append([s[0], float(s[1]) + shift, float(s[2]) + shift, float(s[3]) + shift])

    # build QMin
    QMin1 = deepcopy(QMin)
    QMin1['geo'] = oldgeom + QMin1['geo']
    QMin1['integrals'] = []
    QMin1['AOoverlap'] = []
    remove = ['nacdr', 'grad', 'h', 'soc', 'dm', 'overlap', 'ion']
    for r in remove:
        QMin1 = removekey(QMin1, r)

    # run the calculation
    tasks = gettasks(QMin1)
    WORKDIR = os.path.join(QMin['scratchdir'], 'AOoverlap')
    setupWORKDIR(WORKDIR, tasks, QMin1)
    err = runMOLPRO(WORKDIR, QMin['molpro'], True)

    # get the AO matrix
    outfile = os.path.join(WORKDIR, 'MOLPRO.out')
    out = readfile(outfile)

    # get number of AOs
    for line in out:
        if 'NUMBER OF CONTRACTIONS:' in line:
            nao = int(line.split()[3])
            break

    # find AO overlap matrix
    for iline, line in enumerate(out):
        if ' SYMMETRY BLOCK 1.1' in line:
            break
    else:
        print('Did not find AO overlap matrix!')
        sys.exit(109)

    # detect format
    line = out[iline + 1]
    if 'E' in line:
        formatting = 2
        q = 24
    elif len(line.split()) == 0:
        formatting = 1
    else:
        formatting = 0
        q = len(line.split())

    # get matrix
    AOovl = []
    for irow in range(nao):
        AOovl.append([])
        if formatting == 2:
            iline += 1
            line = out[iline]
            i = 0
            s = []
            while True:
                x = line[1 + i * q:1 + (i + 1) * q]
                s.append(x)
                i += 1
                if 1 + (i + 1) * q > len(line):
                    break
            for y in s:
                AOovl[-1].append(float(y))
        elif formatting == 1:
            iline += 1
            for x in range((nao - 1) // 10 + 1):
                iline += 1
                line = out[iline]
                s = line.split()
                for y in s:
                    AOovl[-1].append(float(y))
        elif formatting == 0:
            for x in range((nao - 1) // q + 1):
                iline += 1
                line = out[iline]
                s = line.split()
                for y in s:
                    AOovl[-1].append(float(y))


    # get off-diagonal block of AO matrix
    AOovl2 = []
    for irow in range(nao // 2):
        # AOovl2.append( AOovl[nao/2+irow][:nao/2])  # lower left block
        AOovl2.append(AOovl[irow][nao // 2:])         # upper right block

    # format string
    # IMPORTANT: upper right block should not be transposed
    string = '%i %i\n' % (nao // 2, nao // 2)
    for irow in range(nao // 2):
        for icol in range(nao // 2):
            string += '%11.8f ' % (AOovl2[irow][icol])
        string += '\n'

    # write to file
    aofile = os.path.join(QMin['savedir'], 'aoovl_double')
    writefile(aofile, string)

    return

# =============================================================================================== #
# =============================================================================================== #
# =======================================  QMout ================================================ #
# =============================================================================================== #
# =============================================================================================== #

# ======================================================================= #


def check_phases(phases):
    thres = 0.95
    newphases = []
    for p in phases:
        if abs(p) < thres:
            print('CASSCF-CASCI overlaps are not unity!')
            sys.exit(110)
        if p < 0.:
            newphases.append(-1.)
        else:
            newphases.append(+1.)
    return newphases


# ======================================================================= #
def getQMout(QMin):

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
            outfile = os.path.join(QMin['scratchdir'], 'master_%i/MOLPRO.out' % (job))
            out = readfile(outfile)
            mults = QMin['multmap'][-job]
            for i in range(nmstates):
                for j in range(nmstates):
                    m1, s1, ms1 = tuple(QMin['statemap'][i + 1])
                    m2, s2, ms2 = tuple(QMin['statemap'][j + 1])
                    if m1 not in mults or m2 not in mults:
                        continue
                    if i == j:
                        QMout['h'][i][j] = getcienergy(out, m1, s1)
                    elif 'soc' in QMin:
                        QMout['h'][i][j] = getsocme(out, i, j, QMin)

    # Dipole Moments
    if 'dm' in QMin:
        # make matrix
        if 'dm' not in QMout:
            QMout['dm'] = [makecmatrix(nmstates, nmstates) for i in range(3)]
        # go through all jobs
        for job in joblist:
            outfile = os.path.join(QMin['scratchdir'], 'master_%i/MOLPRO.out' % (job))
            out = readfile(outfile)
            mults = QMin['multmap'][-job]
            for i in range(nmstates):
                for j in range(nmstates):
                    m1, s1, ms1 = tuple(QMin['statemap'][i + 1])
                    m2, s2, ms2 = tuple(QMin['statemap'][j + 1])
                    if not m1 == m2 or not ms1 == ms2 or m1 not in mults or m2 not in mults:
                        continue
                    for ixyz in range(3):
                        QMout['dm'][ixyz][i][j] = getcidm(out, m1, s1, s2, ixyz)

    # Regular Overlaps
    if 'overlap' in QMin:
        if 'overlap' not in QMout:
            QMout['overlap'] = makecmatrix(nmstates, nmstates)
        for mult in itmult(QMin['states']):
            job = QMin['multmap'][mult]
            outfile = os.path.join(QMin['scratchdir'], 'WFOVL_%i_%i/wfovl.out' % (mult, job))
            out = readfile(outfile)
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

    # Gradients
    if 'grad' in QMin:
        if 'grad' not in QMout:
            QMout['grad'] = [[[0. for i in range(3)] for j in range(natom)] for k in range(nmstates)]
        for grad in QMin['gradmap']:
            outfile = os.path.join(QMin['scratchdir'], 'grad_%i_%i/MOLPRO.out' % grad)
            out = readfile(outfile)
            g = getgrad(out, grad[0], grad[1], natom)
            for istate in QMin['statemap']:
                state = QMin['statemap'][istate]
                if (state[0], state[1]) == grad:
                    QMout['grad'][istate - 1] = g

    # Nonadiabatic couplings
    if 'nacdr' in QMin:
        if 'nacdr' not in QMout:
            QMout['nacdr'] = [[[[0. for i in range(3)] for j in range(natom)] for k in range(nmstates)] for l in range(nmstates)]
        # get phases
        allphases = {}
        for mult in itmult(QMin['states']):
            job = QMin['multmap'][mult]
            outfile = os.path.join(QMin['scratchdir'], 'CI_CAS_%i_%i/wfovl.out' % (mult, job))
            out = readfile(outfile)
            phases = []
            for i in range(QMin['states'][mult - 1]):
                phases.append(getsmate(out, i + 1, i + 1))
            allphases[mult] = check_phases(phases)
        # get NACs
        for nac in QMin['nacmap']:
            outfile = os.path.join(QMin['scratchdir'], 'nac_%i_%i_%i_%i/MOLPRO.out' % nac)
            out = readfile(outfile)
            g = getnacana(out, nac[0], nac[1], nac[3], natom)
            phase1 = allphases[nac[0]][nac[1] - 1]
            phase2 = allphases[nac[0]][nac[3] - 1]
            # print('correcting:',nac,phase1,phase2,phase1*phase2)
            for iatom in range(natom):
                for ixyz in range(3):
                    g[iatom][ixyz] *= phase1 * phase2
            gneg = deepcopy(g)
            for iatom in range(natom):
                for ixyz in range(3):
                    gneg[iatom][ixyz] *= -1
            for istate in QMin['statemap']:
                for jstate in QMin['statemap']:
                    state1 = QMin['statemap'][istate]
                    state2 = QMin['statemap'][jstate]
                    if not state1[2] == state2[2]:
                        continue
                    if (state1[0], state1[1], state2[0], state2[1]) == nac:
                        QMout['nacdr'][istate - 1][jstate - 1] = g
                    elif (state2[0], state2[1], state1[0], state1[1]) == nac:
                        QMout['nacdr'][istate - 1][jstate - 1] = gneg

    return QMout

# ========================== Main Code =============================== #


def main():

    # Retrieve PRINT and DEBUG
    try:
        envPRINT = os.getenv('SH2PRO_PRINT')
        if envPRINT and envPRINT.lower() == 'false':
            global PRINT
            PRINT = False
        envDEBUG = os.getenv('SH2PRO_DEBUG')
        if envDEBUG and envDEBUG.lower() == 'true':
            global DEBUG
            DEBUG = True
    except ValueError:
        print('PRINT or DEBUG environment variables do not evaluate to boolean values!')
        sys.exit(111)

    # Process Command line arguments
    if len(sys.argv) != 2:
        print('Usage:\n./SHARC_MOLPRO.py <QMin>\n')
        print('version:', version)
        print('date:', versiondate)
        print('changelog:\n', changelogstring)
        sys.exit(112)
    QMinfilename = sys.argv[1]

    # Print header
    printheader()

    # Read QMinfile
    QMin = readQMin(QMinfilename)
    printQMin(QMin)
    # pprint.pprint(QMin)

    # get the job schedule
    joblist = generate_joblist(QMin)
    # pprint.pprint(joblist,depth=2)

    # run all the MOLPRO jobs
    errorcodes = runjobs(joblist, QMin)

    # do all necessary overlap and Dyson calculations
    errorcodes = run_wfoverlap(QMin, errorcodes)

    # read all the output files
    QMout = getQMout(QMin)
    if PRINT or DEBUG:
        printQMout(QMin, QMout)

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

    if PRINT or DEBUG:
        print(datetime.datetime.now())
        print('#================ END ================#')


if __name__ == '__main__':
    main()
