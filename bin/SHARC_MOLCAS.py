#!/usr/bin/env python3

# ******************************************
#
#    SHARC Program Suite
#
#    Copyright (c) 2023 University of Vienna
#    Copyright (c) 2023 University of Minnesota
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

# ======================================================================= #
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
from socket import gethostname
import itertools
# write debug traces when in pool threads
import traceback
# for diabatization
import numpy as np



# ======================================================================= #

version = '3.0'
versiondate = datetime.date(2023, 2, 3)



changelogstring = '''
07.07.2014:
- Gradients can be setup with MOLCAS in parallel fashion, using 1 core per gradient.
- QM/MM support added (using MOLCAS and TINKER).

03.10.2014:
- Changed title lines
- readQMin was harmonized with the other interfaces:
  * reads MOLCAS.template as before
  * reads SH2CAS.inp (keywords "molcas", "scratchdir", "memory")
  * Project is the comment on the second line of QM.in, stripped of whitespace and prepended with "Comment_"
  * SHQM directory is the pwd (the directory where the interface is started
  * gradaccudefault and gradaccumax as in MOLPRO interface
- for SS-CASSCF (in any given multiplicity) no MCLR is executed
- changed getgrad to also find SS-CASSCF gradients
- in writeMOLCASinput, section task[0]=="ddrdiab" fixed a bug (mixup of variables mult and i)

08.10.2014:     1.0
- official release version, no changes to 0.3

05.02.2015:
- fixed a bug in getsmate(), where the overlap matrix is read from the wrong RASSI run

06.02.2015:
- major rewrite started...

12.02.2015:
- major rewrite finished
- New features:
  * Numerical gradients for Cholesky-based methods, CASPT2 and MS-CASPT2
  * Dipole moment derivatives and spin-orbit coupling derivatives
  * Numerical differentiation parallelized
  * Project is now always "MOLCAS"
  * Restart if MCLR did not converge
- most of the code was redesigned from scratch to make it easier to maintain
- backwards compatibility to input for version 1.0

18.02.2015:
- fixed a bug with SOC matrix element readout
- added a delay time (default 0 sec) for starting parallel jobs, which can be set in SH2CAS.inp
- Default is Douglas-Kroll, to be backwards-compatible with MOLCAS interface 1.0. DKH can be turned off with option "no-douglas-kroll"
- fixed a bug with DKH-CASSCF Energy readout

20.02.2015:
- made retrieving the hostname more portable
- made molcas call more portable (pipes)
- made getversion() more portable by reading from $MOLCAS/.molcasversion
- added support for states of different electron number (number of active electrons might be one lower than given)
- added support for gradmode=1 with QM/MM
- changed: SS-CASSCF gradients always lead to gradmode=1

26.03.2015:
- MOLCAS version is used to write the AMFI keyword into the appropriate section (&GATEWAY or &SEWARD)
- MOLCAS.template keyword "cholesky_analytical" leads to analytical gradient calculations for CD-CASSCF. Without the keyword, CD-CASSCF gradients are evaluated numerically.

28.05.2015:
- added "frozen" keyword to modify the number of frozen orbitals in CASPT2

10.11.2015:
- added support for "ion" (Dyson calculations) using the wfoverlap code
- fixed gradients for MS-CASPT2 calcs with only one state
- fixed gradients for (N-1)-electron states
- added "molden" keyword

04.08.2015:
- added "angmom" keyword for GATEWAY, so that SOC can be calculated for Natom<3

23.08.2017:
- added the "rootpad" keyword to the template, which can be used to request extra, zero-weight states in SA
- added the "baslib" keyword to the template, which can be used to load custom basis sets
- Resource file is now called "MOLCAS.resources" instead of "SH2CAS.inp" (old filename still works)
- The connection table file for QM/MM is now called MOLCAS.qmmm.table (old filename still works)

24.08.2017:
- added the numfrozcore and numocc keywords for Dyson norm calculations
- gradaccumax, gradaccudefault, displ are now keywords in the template, not in the resources file (not backwards compatible)
- cleanup keyword now available

29.01.2018:     2.0
- rework for openMOLCAS

29.06.2018:
- can do nonadiabatic couplings (nacdr keyword)

09.11.2018:
- can work with MPI parallel MOLCAS
- new keyword: "mpi_parallel" in resources file activates MPI parallel mode
- no dynamics scheduling due to integral copying issues (every gradient is run with the same ncpu as the master)

15.11.2018:
- can now do PCM computations (only numerical gradients)

11.10.2020:
- COBRAMM can be used for QM/MM calculation

25.05.2021:
- numerical gradients with diabatization (diab_num_grad keyword)

06.08.2022:
- can now do CMS-PDFT with OpenMolcas

03.02.2023:
- convert to python 3
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

# conversion factors
au2a = 0.529177211
rcm_to_Eh = 4.556335e-6

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
    string += '||' + ' ' * 27 + 'SHARC - MOLCAS - Interface' + ' ' * 27 + '||\n'
    string += '||' + ' ' * 80 + '||\n'
    string += '||' + ' ' * 21 + 'Authors: Sebastian Mai, Martin Richter,' + ' ' * 20 + '||\n'
    string += '||' + ' ' * 15 + 'Severin Polonius, Philipp Marquetand and Yinan Shu' + ' ' * 15 + '||\n'
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
    if 'phases' in QMin:
        string += '\tPhases'
    print(string)

    string = 'States: '
    for i in itmult(QMin['states']):
        string += '\t%i %s' % (QMin['states'][i - 1], IToMult[i])
    print(string)

    string = 'Method: \t'
    string += 'SA(%i' % (QMin['template']['roots'][0])
    for i in QMin['template']['roots'][1:]:
        string += '|%i' % (i)
    string += ')-'
    string += QMin['template']['method'].upper()
    string += '(%s,%s)/%s' % (tuple(QMin['template']['nactel']), (QMin['template']['ras1'],QMin['template']['ras2'],QMin['template']['ras3']), QMin['template']['basis'])
    parts = []
    if QMin['template']['cholesky']:
        parts.append('RICD')
    if not QMin['template']['no-douglas-kroll']:
        parts.append('Douglas-Kroll')
    if QMin['method'] > 0 and QMin['template']['ipea'] != 0.25:
        parts.append('IPEA=%4.2f' % (QMin['template']['ipea']))
    if QMin['method'] > 0 and QMin['template']['imaginary'] != 0.00:
        parts.append('Imaginary Shift=%4.2f' % (QMin['template']['imaginary']))
    if QMin['template']['frozen'] != -1:
        parts.append('CASPT2 frozen orbitals=%i' % (QMin['template']['frozen']))
    if QMin['template']['diab_num_grad']:
        parts.append('diabatization for num. grad')
    if len(parts) > 0:
        string += '\t('
        string += ','.join(parts)
        string += ')'
    if QMin['template']['qmmm']:
        string += '\t+Tinker'
    if 'cobramm' in QMin:
        string += '\t+COBRAMM'
    print(string)
    # say, if CAS(n-1,m) is used for any multiplicity
    oddmults = False
    for i in QMin['statemap'].values():
        if (QMin['template']['nactel'][0] + i[0]) % 2 == 0:
            oddmults = True
    if oddmults:
        string = '\t\t' + ['Even ', 'Odd '][QMin['template']['nactel'][0] % 2 == 0]
        string += 'numbers of electrons are treated with CAS(%i,%i).' % (QMin['template']['nactel'][0] - 1, QMin['template']['ras2'])
        print(string)
    # CAS(2,2) does not allow for SOC calculations (bug in MOLCAS)
    # if QMin['template']['nactel']==2 and QMin['template']['ras2']==2:
        # if 'soc' in QMin or 'socdr' in QMin:
        # print('WARNING: CAS(2,2) yields zero cm^-1 for all SOC matrix elements in MOLCAS!')

    string = 'Found Geo'
    if 'veloc' in QMin:
        string += ' and Veloc! '
    else:
        string += '! '
    string += 'NAtom is %i.\n' % (QMin['natom'])
    print(string)

    string = '\nGeometry in Bohrs:\n'
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

    for i in QMin:
        if not any([i == j for j in ['h', 'dm', 'soc', 'dmdr', 'socdr', 'geo', 'veloc', 'states', 'comment', 'LD_LIBRARY_PATH', 'grad', 'nacdr', 'ion', 'overlap', 'template']]):
            if not any([i == j for j in ['ionlist', 'ionmap']]) or DEBUG:
                string = i + ': '
                string += str(QMin[i])
                print(string)
    print('\n')
    sys.stdout.flush()

# ======================================================================= #


def printtasks(tasks):
    '''If PRINT, prints a formatted table of the tasks in the tasks list.

    Arguments:
    1 list of lists: tasks list (see gettasks for specs)'''

    # if DEBUG:
    # pprint.pprint(tasks)
    if not PRINT:
        return
    print('==> Task Queue:\n')
    for i in range(len(tasks)):
        task = tasks[i]
        if task[0] == 'gateway':
            print('GATEWAY')
        elif task[0] == 'seward':
            print('SEWARD')
        elif task[0] == 'link':
            print('Link\t%s\t--> \t%s' % (task[2], task[1]))
        elif task[0] == 'copy':
            print('Copy\t%s\t==> \t%s' % (task[1], task[2]))
        elif task[0] == 'rasscf':
            print('RASSCF\tMultiplicity: %i\tStates: %i\tJOBIPH=%s\tLUMORB=%s' % (task[1], task[2], task[3], task[4]))
        elif task[0] == 'rasscf-cms':
            print('RASSCF-CMS\tMultiplicity: %i\tStates: %i\tJOBIPH=%s\tLUMORB=%s\tRLXROOT=%i' % (task[1], task[2], task[3], task[4], task[5]))
        # elif task[0]=='rasscf-rlx':
            # print('RASSCF\tMultiplicity: %i\tStates: %i\tRLXROOT=%i' % (task[1], task[2], task[3]))
        elif task[0] == 'alaska':
            print('ALASKA')
        elif task[0] == 'mclr':
            print('MCLR')
        elif task[0] == 'mclr-cms':
            print('MCLR-CMS')
        elif task[0] == 'caspt2':
            print('CASPT2\tMultiplicity: %i\tStates: %i\tMULTISTATE=%s' % (task[1], task[2], task[3]))
        elif task[0] == 'cms-pdft':
            print('CMS-PDFT\tFunctional: %s' % (task[1]))
        elif task[0] == 'rassi':
            print('RASSI\t%s\tStates: %s' % ({'soc': 'Spin-Orbit Coupling', 'dm': 'Dipole Moments', 'overlap': 'Overlaps'}[task[1]], task[2]))
        elif task[0] == 'espf':
            print('ESPF')
        elif task[0] == 'cobramm':
            print('COBRAMM')
        else:
            print(task)
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
def getversion(out, MOLCAS):
    allowedrange = [(18.0,25.999), (8.29999, 9.30001)]
    # first try to find $MOLCAS/.molcasversion
    molcasversion = os.path.join(MOLCAS, '.molcasversion')
    if os.path.isfile(molcasversion):
        vf = open(molcasversion)
        string = vf.readline()
        vf.close()
        print('Content: "%s"\n' % string)
    # otherwise try to read this from the output file
    else:
        string = ''
        for i in range(50):
            line = out[i]
            s = line.split()
            for j, el in enumerate(s):
                if 'version' in el:
                    string = s[j + 1]
                    break
            if string != '':
                break
    a = re.search('[0-9]+\\.[0-9]+', string)
    if a is None:
        print('No MOLCAS version found.\nCheck whether MOLCAS path is set correctly in MOLCAS.resources\nand whether $MOLCAS/.molcasversion exists.')
        sys.exit(17)
    v = float(a.group())
    if not any([i[0] <= v <= i[1] for i in allowedrange]):
        # allowedrange[0]<=v<=allowedrange[1]:
        print('MOLCAS version %3.1f not supported! ' % (v))
        sys.exit(18)
    if DEBUG:
        print('Found MOLCAS version %3.1f\n' % (v))
    return v

# ======================================================================= #


def getcienergy(out, mult, state, version, method, dkh):
    '''Searches a complete MOLCAS output file for the MRCI energy of (mult,state).

    Arguments:
    1 list of strings: MOLCAS output
    2 integer: mult
    3 integer: state

    Returns:
    1 float: total CI energy of specified state in hartree'''

    if method == 0:
        modulestring = '&RASSCF'
        spinstring = 'Spin quantum number'
        if dkh:
            energystring = 'Final state energy(ies)'
        else:
            energystring = '::    RASSCF root number'
            stateindex = 4
            enindex = 7
    elif method == 1:
        modulestring = '&CASPT2'
        spinstring = 'Spin quantum number'
        energystring = '::    CASPT2 Root'
        stateindex = 3
        enindex = 6
    elif method == 2:
        modulestring = '&CASPT2'
        spinstring = 'Spin quantum number'
        energystring = '::    MS-CASPT2 Root'
        stateindex = 3
        enindex = 6
    elif method == 3:
        modulestring = '&MCPDFT'
        spinstring = 'Spin quantum number'
        energystring = 'Total MC-PDFT energy for state'
        stateindex = 5
        enindex = 6
    elif method == 4:
        modulestring = '&MCPDFT'
        spinstring = 'Spin quantum number'
        energystring = '::    XMS-PDFT Root'
        stateindex = 3
        enindex = 6
    elif method == 5:
        modulestring = '&MCPDFT'
        spinstring = 'Spin quantum number'
        energystring = '::    CMS-PDFT Root'
        stateindex = 3
        enindex = 6
    elif method == 6:
        modulestring = '&CASPT2'
        spinstring = 'Spin quantum number'
        energystring = '::    XMS-CASPT2 Root'
        stateindex = 3
        enindex = 6


    module = False
    correct_mult = False
    for i, line in enumerate(out):
        if modulestring in line:
            module = True
        elif spinstring in line and module:
            spin = float(line.split()[3])
            if int(2 * spin) + 1 == mult:
                correct_mult = True
        elif energystring in line and module and correct_mult:
            if method == 0 and dkh:
                l = out[i + 4 + state].split()
                return float(l[1])
            else:
                l = line.split()
                if int(l[stateindex]) == state:
                    return float(l[enindex])
    print('CI energy of state %i in mult %i not found!' % (state, mult))
    sys.exit(19)

# ======================================================================= #


def getcidm(out, mult, state1, state2, pol, version):
    # two cases:
    # - Dipole moments are in RASSI calculation with only one JOBIPH file
    # - Dipole moments are in RASSI calculation with two JOBIPH files of same multiplicity

    if pol == 'X' or pol == 'Y' or pol == 'Z':
        pol = IToPol[pol]

    modulestring = '&RASSI'
    spinstring = 'SPIN MULTIPLICITY:'
    stopstring = 'The following data are common to all the states'
    stop2string = 'Special properties section'
    statesstring = 'Nr of states:'
    #matrixstring=' PROPERTY: MLTPL  1   COMPONENT:   %i' % (pol+1)
    matrixstring = 'PROPERTY: MLTPL  1   COMPONENT:*   %i' % (pol + 1)

    # first, find the correct RASSI output section for the given multiplicity
    module = False
    jobiphmult = []

    for iline, line in enumerate(out):
        if modulestring in line:
            module = True
            jobiphmult = []
        elif module:
            if spinstring in line:
                jobiphmult.append(int(line.split()[-1]))
            if stopstring in line:
                if all(i == mult for i in jobiphmult):
                    break
                else:
                    module = False
    else:
        print('DM element not found!', mult, state1, state2, pol)
        print('No RASSI run for multiplicity %i found!' % (mult))
        sys.exit(20)

    # Now start searching at iline, looking for the requested matrix element
    for jline, line in enumerate(out[iline + 1:]):
        if stop2string in line:
            print('DM element not found!', mult, state1, state2, pol)
            print('Found correct RASSI run, but too few matrix elements!')
            sys.exit(21)
        if statesstring in line:
            nstates = int(line.split()[-1])
            if len(jobiphmult) == 2:
                stateshift = nstates // 2
            else:
                stateshift = 0
        if matrixstring in line:
            block = (stateshift + state2 - 1) // 4
            rowshift = 3 + stateshift + state1 + (6 + nstates) * block
            colshift = 1 + (stateshift + state2 - 1) % 4

            return float(out[iline + jline + rowshift + 1].split()[colshift])


# ======================================================================= #
def getMOLCASstatenumber(mult, state, ms, states):
    statenumber = 0
    for m, mstates in enumerate(states):  # iterate over multiplicities
        if m + 1 < mult:
            statenumber += mstates * (m + 1)
        else:  # correct multiplicity found
            for nstate in range(1, mstates + 1):  # iterate over states
                if nstate < state:
                    statenumber += m + 1
                else:  # correct state found
                    statenumber += int(ms + 1 + 0.5 * (mult))
                    return statenumber
    print('getMOLCASstatenumber Error: mult=%i, state=%i, ms=%i not in' % (mult, state, ms), states)
    quit(1)

# ======================================================================= #
# SOCME_START_ILINE=-1
# SOCME_FILE_ID=-1


def getsocme(out, mult1, state1, ms1, mult2, state2, ms2, states, version, method, dkh):
    '''Searches a MOLCAS output for an element of the Spin-Orbit hamiltonian matrix. Also converts from cm^-1 to hartree and adds the diagonal shift.

    Arguments:
    1 list of strings: MOLCAS output
    2-4 integer: multiplicity, state and ms for state1
    5-7 integer: multiplicity, state and ms for state2
    8 list of integer: states specs

    Returns:
    1 complex: SO hamiltonian matrix element in hartree'''
    rcm_to_Eh = 4.556335e-6
    socstring = 'I1  S1  MS1    I2  S2  MS2    Real part    Imag part      Absolute'
    stopstring = '----------------------------------------------------------------------'

    # return diagonal elements
    if mult1 == mult2 and state1 == state2 and ms1 == ms2:
        return complex(getcienergy(out, mult1, state1, version, method, dkh), 0.0)

    if len(states) == 1:
        return complex(0.0, 0.0)

    if mult1 == mult2 == 1:
        return complex(0.0, 0.0)

    # otherwise, find state indices s1 and s2
    s1 = getMOLCASstatenumber(mult1, state1, ms1, states)
    s2 = getMOLCASstatenumber(mult2, state2, ms2, states)

    ## accelerated version of finding the SOC section if it was already found for this particular output file
    #if SOCME_START_ILINE==-1 or SOCME_FILE_ID!=hash(tuple(out)):
    #    ## look for spin-orbit section
    #    for iline,line in enumerate(out):
    #        if socstring in line:
    #            break
    #    else:
    #        print('No Spin-Orbit section found in output!')
    #        sys.exit(22)
    #    global SOCME_START_ILINE
    #    SOCME_START_ILINE=iline
    #    global SOCME_FILE_ID
    #    SOCME_FILE_ID=hash(tuple(out))
    #else:
    #    iline=SOCME_START_ILINE

    # look for spin-orbit section
    for iline, line in enumerate(out):
        if socstring in line:
            break
    else:
        print('No Spin-Orbit section found in output!')
        sys.exit(23)

    # look for matrix element
    socme = complex(0.0, 0.0)
    while True:
        iline += 1
        if stopstring in out[iline]:
            break
        l = out[iline].split()
        I1 = int(l[0])
        I2 = int(l[3])
        if (I1, I2) == (s1, s2) or (I1, I2) == (s2, s1):
            if s2 <= s1:
                socme = complex(float(l[6]), +float(l[7]))
            else:
                socme = complex(float(l[6]), -float(l[7]))
            break
    return socme * rcm_to_Eh

# ======================================================================= #


def getgrad(out, mult, state, QMin, version):

    espf = QMin['template']['qmmm']
    natom = QMin['natom']
    roots = QMin['template']['roots']
    ssgrad = roots[mult - 1] == 1          # for SS-CASSCF no MCLR is needed

    if espf:
        gradstring = 'Molecular gradients, after ESPF'
        versionshift = 8

    mclr = False
    alaska = False
    rasscf = False
    multfound = False
    statefound = False
    grad = []

    for i, line in enumerate(out):
        if '&RASSCF' in line:
            rasscf = True
            mclr = False
            alaska = False
        elif '&MCLR' in line:
            rasscf = False
            mclr = True
            alaska = False
        elif '&ALASKA' in line and multfound and statefound:
            rasscf = False
            mclr = False
            alaska = True
        elif 'Spin quantum number' in line and rasscf:
            if int(round(float(line.split()[3]) * 2)) + 1 == mult:
                multfound = True
            else:
                multfound = False
            if ssgrad:
                statefound = True
            else:
                statefound = False
        elif ' Lagrangian multipliers are calculated for state no.' in line and mclr and multfound:
            if int(line.split()[7].strip()) == state:
                statefound = True
        elif alaska and multfound and statefound:
            if espf:
                if gradstring in line:
                    for atom in range(natom):
                        if atom + 1 in QMin['active_qmmm_atoms']:
                            atomindex = QMin['active_qmmm_atoms'].index(atom + 1)
                            grad.append([float(out[i + versionshift + atomindex].split()[xyz + 1]) for xyz in range(3)])
                        else:
                            grad.append([0.0 for xyz in range(3)])
                    return grad
            else:
                if 'Molecular gradients' in line:
                    for atom in range(natom):
                        grad.append([float(out[i + 8 + atom].split()[1 + xyz]) for xyz in range(3)])
                    return grad
    else:
        print('Gradient of state %i in mult %i not found!' % (state, mult))
        sys.exit(24)


# ======================================================================= #
def getefield(out, mult, state, QMin):

    roots = QMin['template']['roots']
    ssgrad = roots[mult - 1] == 1          # for SS-CASSCF no MCLR is needed

    rasscf, mclr, alaska = False, False, False
    multfound, statefound = False, False
    efield = []

    for i, line in enumerate(out):
        if '&RASSCF' in line:
            rasscf, mclr, alaska = True, False, False
        elif '&MCLR' in line:
            rasscf, mclr, alaska = False, True, False
        elif '&ALASKA' in line and multfound and statefound:
            rasscf, mclr, alaska = False, False, True
        elif 'Spin quantum number' in line and rasscf:
            if int(round(float(line.split()[3]) * 2)) + 1 == mult:
                multfound = True
            else:
                multfound = False
            if ssgrad:
                statefound = True
            else:
                statefound = False
        elif ' Lagrangian multipliers are calculated for state no.' in line and mclr and multfound:
            if int(line.split()[7].strip()) == state:
                statefound = True
        elif mclr and multfound and statefound:
            if 'Electric field:' in line:
                natomMM, efield = 0, []
                while True:
                    try:
                        efield.append([float(out[i + 2 + natomMM].split()[1 + xyz]) for xyz in range(3)])
                    except IndexError:
                        break
                    natomMM += 1
                return efield
    else:
        print('Electric field of state %i in mult %i not found!' % (state, mult))
        sys.exit(24)

# ======================================================================= #


def getnacdr(out, mult, state, state2, QMin, version):

    # espf=QMin['template']['qmmm']
    natom = QMin['natom']
    roots = QMin['template']['roots']
    # ssgrad=roots[mult-1]==1          # for SS-CASSCF no MCLR is needed

    #if espf:
    #    if 7<=version<8:
    #        gradstring='After ESPF, gradients are'
    #        versionshift=2
    #    elif 8<=version<9:
    #        gradstring='Molecular gradients, after ESPF'
    #        versionshift=8

    mclr = False
    alaska = False
    rasscf = False
    multfound = False
    statefound = False
    grad = []

    for i, line in enumerate(out):
        if '&RASSCF' in line:
            rasscf = True
            mclr = False
            alaska = False
        elif '&MCLR' in line:
            rasscf = False
            mclr = True
            alaska = False
        elif '&ALASKA' in line and multfound and statefound:
            rasscf = False
            mclr = False
            alaska = True
        elif 'Spin quantum number' in line and rasscf:
            if int(round(float(line.split()[3]) * 2)) + 1 == mult:
                multfound = True
            else:
                multfound = False
            # if ssgrad:
                # statefound=True
            # else:
                # statefound=False
        elif ' Lagrangian multipliers are calculated for state' in line and mclr and multfound:
            if 'states' in line:
                s = line.split()
                s1 = int(s[-2].replace('/', ''))
                s2 = int(s[-1])
                if (s1, s2) == (state, state2):
                    # if int(line.split()[7].strip()) == state:
                    statefound = True
        elif alaska and multfound and statefound:
            #if espf:
            #    if gradstring in line:
            #        for atom in range(natom):
            #            if atom+1 in QMin['active_qmmm_atoms']:
            #                atomindex=QMin['active_qmmm_atoms'].index(atom+1)
            #                grad.append([ float(out[i+versionshift+atomindex].split()[xyz+1]) for xyz in range(3) ])
            #            else:
            #                grad.append([ 0.0 for xyz in range(3)])
            #        return grad
            #else:
            if 'Total derivative coupling' in line:
                for atom in range(natom):
                    grad.append([float(out[i + 8 + atom].split()[1 + xyz]) for xyz in range(3)])
                return grad
    else:
        print('NAC vector of states %i / %i in mult %i not found!' % (state, state2, mult))
        sys.exit(25)




# ======================================================================= #
def getsmate(out, mult, state1, state2, states):
    # one case:
    # - Dipole moments are in RASSI calculation with two JOBIPH files of same multiplicity

    modulestring = '&RASSI'
    spinstring = 'SPIN MULTIPLICITY:'
    stopstring = 'The following data are common to all the states'
    stop2string = 'MATRIX ELEMENTS OF 1-ELECTRON OPERATORS'
    statesstring = 'Nr of states:'
    matrixstring = 'OVERLAP MATRIX FOR THE ORIGINAL STATES:'

    # first, find the correct RASSI output section for the given multiplicity
    module = False
    jobiphmult = []

    for iline, line in enumerate(out):
        if modulestring in line:
            module = True
            jobiphmult = []
        elif module:
            if spinstring in line:
                jobiphmult.append(int(line.split()[-1]))
            if stopstring in line:
                if jobiphmult == [mult, mult]:
                    break
                else:
                    module = False
    else:
        print('Overlap element not found!', mult, state1, state2)
        print('No correct RASSI run for multiplicity %i found!' % (mult))
        sys.exit(26)

    # Now start searching at iline, looking for the requested matrix element
    for jline, line in enumerate(out[iline + 1:]):
        if stop2string in line:
            print('Overlap element not found!', mult, state1, state2)
            print('Found correct RASSI run, but too few matrix elements!')
            sys.exit(27)
        if statesstring in line:
            nstates = int(line.split()[-1])
        if matrixstring in line:
            rowshift = 1
            for i in range(nstates // 2 + state2 - 1):
                rowshift += i // 5 + 1
            rowshift += 1 + (state1 - 1) // 5
            colshift = (state1 - 1) % 5

            return float(out[iline + jline + rowshift + 1].split()[colshift])

# ======================================================================= #


def getQMout(out, QMin):
    '''Constructs the requested matrices and vectors using the get<quantity> routines.

    The dictionary QMout contains all the requested properties. Its content is dependent on the keywords in QMin:
    - 'h' in QMin:
                    QMout['h']: list(nmstates) of list(nmstates) of complex, the non-relaticistic hamiltonian
    - 'soc' in QMin:
                    QMout['h']: list(nmstates) of list(nmstates) of complex, the spin-orbit hamiltonian
    - 'dm' in QMin:
                    QMout['dm']: list(3) of list(nmstates) of list(nmstates) of complex, the three dipole moment matrices
    - 'grad' in QMin:
                    QMout['grad']: list(nmstates) of list(natom) of list(3) of float, the gradient vectors of every state (even if "grad all" was not requested, all nmstates gradients are contained here)
    - 'nac' in QMin and QMin['nac']==['num']:
                    QMout['nac']: list(nmstates) of list(nmstates) of complex, the non-adiabatic coupling matrix
                    QMout['mrcioverlap']: list(nmstates) of list(nmstates) of complex, the MRCI overlap matrix
                    QMout['h']: like with QMin['h']
    - 'nac' in QMin and QMin['nac']==['ana']:
                    QMout['nac']: list(nmstates) of list(nmstates) of list(natom) of list(3) of float, the matrix of coupling vectors
    - 'nac' in QMin and QMin['nac']==['smat']:
                    QMout['nac']: list(nmstates) of list(nmstates) of complex, the adiabatic-diabatic transformation matrix
                    QMout['mrcioverlap']: list(nmstates) of list(nmstates) of complex, the MRCI overlap matrix
                    QMout['h']: like with QMin['h']

    Arguments:
    1 list of strings: Concatenated MOLCAS output
    2 dictionary: QMin

    Returns:
    1 dictionary: QMout'''


    # get version of MOLCAS
    version = QMin['version']
    method = QMin['method']

    # Currently implemented keywords: h, soc, dm, grad, nac (num,ana,smat)
    states = QMin['states']
    nstates = QMin['nstates']
    nmstates = QMin['nmstates']
    natom = QMin['natom']
    QMout = {}
    # h: get CI energies of all ci calculations and construct hamiltonian, returns a matrix(nmstates,nmstates)
    if 'h' in QMin:
        # no spin-orbit couplings, hamilton operator diagonal, only one loop
        h = makecmatrix(nmstates, nmstates)
        for istate, i in enumerate(QMin['statemap'].values()):
            mult, state, ms = tuple(i)
            if states[mult - 1] == 1 and method in [2,6]:
                method1 = 1
            else:
                method1 = method
            h[istate][istate] = complex(getcienergy(out, mult, state, version, method1, not QMin['template']['no-douglas-kroll']))
        QMout['h'] = h
    # SOC: get SOC matrix and construct hamiltonian, returns a matrix(nmstates,nmstates)
    if 'soc' in QMin:
        # soc: matrix is not diagonal, two nested loop
        soc = makecmatrix(nmstates, nmstates)
        for istate, i in enumerate(QMin['statemap'].values()):
            for jstate, j in enumerate(QMin['statemap'].values()):
                mult1, state1, ms1 = tuple(i)
                mult2, state2, ms2 = tuple(j)
                if states[mult1 - 1] == 1 and method in [2,6]:
                    method1 = 1
                else:
                    method1 = method
                soc[istate][jstate] = getsocme(out, mult1, state1, ms1, mult2, state2, ms2, QMin['states'], version, method1, not QMin['template']['no-douglas-kroll'])
        QMout['h'] = soc
    # DM: get vector of three dipole matrices, three nested loops, returns a list of three matrices(nmstates,nmstates)
    if 'dm' in QMin:
        dm = []
        for xyz in range(3):
            dm.append(makecmatrix(nmstates, nmstates))
            for istate, i in enumerate(QMin['statemap']):
                for jstate, j in enumerate(QMin['statemap']):
                    mult1, state1, ms1 = tuple(QMin['statemap'][i])
                    mult2, state2, ms2 = tuple(QMin['statemap'][j])
                    if mult1 == mult2 and ms1 == ms2:
                        dm[xyz][istate][jstate] = complex(getcidm(out, mult1, state1, state2, xyz, version))
                    else:
                        dm[xyz][istate][jstate] = complex(0.0)
        QMout['dm'] = dm
    # Grad:  returns a list of nmstates vectors
    if 'grad' in QMin:
        grad = []
        for i in QMin['statemap']:
            mult, state, ms = tuple(QMin['statemap'][i])
            if (mult, state) in QMin['gradmap']:
                grad.append(getgrad(out, mult, state, QMin, version))
            else:
                gradatom = []
                for iatom in range(natom):
                    gradatom.append([0., 0., 0.])
                grad.append(gradatom)
        QMout['grad'] = grad
    # EField:  returns the electric field computed at given points, for COBRAMM QM/MM
    if 'grad' in QMin and QMin['template']['cobramm']:
        efield = []
        for i in QMin['statemap']:
            mult, state, ms = tuple(QMin['statemap'][i])
            if (mult, state) in QMin['gradmap']:
                efield.append(getefield(out, mult, state, QMin))
            else:
                efield.append([])
        QMout['efield'] = efield
    # NACdr:  returns a list of nmstates x nmstates vectors
    if 'nacdr' in QMin:
        nacdr = [[[[0. for i in range(3)] for j in range(natom)] for k in range(nmstates)] for l in range(nmstates)]
        for i, i1 in enumerate(QMin['statemap']):
            mult, state, ms = tuple(QMin['statemap'][i1])
            for j, j1 in enumerate(QMin['statemap']):
                if not i < j:
                    continue
                mult2, state2, ms2 = tuple(QMin['statemap'][j1])
                if (mult, state, mult2, state2) in QMin['nacmap']:
                    nacdr[i][j] = getnacdr(out, mult, state, state2, QMin, version)
                    nacdr[j][i] = deepcopy(nacdr[i][j])
                    for x in range(natom):
                        for y in range(3):
                            nacdr[j][i][x][y] *= -1.
        QMout['nacdr'] = nacdr
    # NAC: case of keyword "smat": returns a matrix(nmstates,nmstates)
    if 'overlap' in QMin:
        nac = makecmatrix(nmstates, nmstates)
        for istate, i in enumerate(QMin['statemap']):
            for jstate, j in enumerate(QMin['statemap']):
                mult1, state1, ms1 = tuple(QMin['statemap'][i])
                mult2, state2, ms2 = tuple(QMin['statemap'][j])
                if 'overlap_nacs' in QMin:
                    if mult1 not in [i[0] for i in QMin['nacmap']]:
                        continue
                    if mult2 not in [i[0] for i in QMin['nacmap']]:
                        continue
                if mult1 == mult2 and ms1 == ms2:
                    nac[istate][jstate] = complex(getsmate(out, mult1, state1, state2, states))
                else:
                    nac[istate][jstate] = complex(0.0)
        QMout['overlap'] = nac
    # Phases from overlaps
    if 'phases' in QMin:
        if 'phases' not in QMout:
            QMout['phases'] = [complex(1., 0.) for i in range(nmstates)]
        if 'overlap' in QMout:
            for i in range(nmstates):
                if QMout['overlap'][i][i].real < 0.:
                    QMout['phases'][i] = complex(-1., 0.)
    return QMout

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
    if 'nacdr' in QMin:
        string += writeQMoutnacana(QMin, QMout)
    if 'overlap' in QMin:
        string += writeQMoutnacsmat(QMin, QMout)
    if 'socdr' in QMin:
        string += writeQMoutsocdr(QMin, QMout)
    if 'dmdr' in QMin:
        string += writeQMoutdmdr(QMin, QMout)
    if 'ion' in QMin:
        string += writeQMoutprop(QMin, QMout)
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
            sys.exit(28)
    else:
        try:
            os.makedirs(SCRATCHDIR)
        except OSError:
            print('Cannot create SCRATCHDIR=%s\n' % (SCRATCHDIR))
            sys.exit(29)

# ======================================================================= #


def removequotes(string):
    if string.startswith("'") and string.endswith("'"):
        return string[1:-1]
    elif string.startswith('"') and string.endswith('"'):
        return string[1:-1]
    else:
        return string

# ======================================================================= #


def getsh2caskey(sh2cas, key):
    i = -1
    while True:
        i += 1
        try:
            line = re.sub('#.*$', '', sh2cas[i])
        except IndexError:
            break
        line = line.split(None, 1)
        if line == []:
            continue
        if key.lower() in line[0].lower():
            return line
    return ['', '']

# ======================================================================= #


def get_sh2cas_environ(sh2cas, key, environ=True, crucial=True):
    line = getsh2caskey(sh2cas, key)
    if line[0]:
        LINE = line[1]
        LINE = removequotes(LINE).strip()
    else:
        if environ:
            LINE = os.getenv(key.upper())
            if not LINE:
                if crucial:
                    print('Either set $%s or give path to %s in MOLCAS.resources!' % (key.upper(), key.upper()))
                    sys.exit(30)
                else:
                    return ''
        else:
            if crucial:
                print('Give path to %s in MOLCAS.resources!' % (key.upper()))
                sys.exit(31)
            else:
                return ''
    LINE = os.path.expandvars(LINE)
    LINE = os.path.expanduser(LINE)
    if containsstring(';', LINE):
        print("$%s contains a semicolon. Do you probably want to execute another command after %s? I can't do that for you..." % (key.upper(), key.upper()))
        sys.exit(32)
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
            sys.exit(33)
        if 'end' in line:
            break
        fields = line.split()
        try:
            nacpairs.append([int(fields[0]), int(fields[1])])
        except ValueError:
            print('"nacdr select" is followed by pairs of state indices, each pair on a new line!')
            sys.exit(34)
    return nacpairs, i

# ======================================================================= #         OK


def readQMin(QMinfilename):
    '''Reads the time-step dependent information from QMinfilename. This file contains all information from the current SHARC job: geometry, velocity, number of states, requested quantities along with additional information. The routine also checks this input and obtains a number of environment variables necessary to run MOLCAS.

    Steps are:
    - open and read QMinfilename
    - Obtain natom, comment, geometry (, velocity)
    - parse remaining keywords from QMinfile
    - check keywords for consistency, calculate nstates, nmstates
    - obtain environment variables for path to MOLCAS and scratch directory, and for error handling

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
        sys.exit(35)
    QMin['natom'] = natom
    if len(QMinlines) < natom + 4:
        print('Input file must contain at least:\nnatom\ncomment\ngeometry\nkeyword "states"\nat least one task')
        sys.exit(36)

    # Save Comment line
    QMin['comment'] = QMinlines[1]

    # Get geometry and possibly velocity (for backup-analytical non-adiabatic couplings)
    QMin['geo'] = []
    QMin['veloc'] = []
    hasveloc = True
    for i in range(2, natom + 2):
        if not containsstring('[a-zA-Z][a-zA-Z]?[0-9]*.*[-]?[0-9]+[.][0-9]*.*[-]?[0-9]+[.][0-9]*.*[-]?[0-9]+[.][0-9]*', QMinlines[i]):
            print('Input file does not comply to xyz file format! Maybe natom is just wrong.')
            sys.exit(37)
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
            sys.exit(38)
    else:
        factor = 1. / au2a

    for iatom in range(len(QMin['geo'])):
        for ixyz in range(3):
            QMin['geo'][iatom][ixyz + 1] *= factor


    if 'states' not in QMin:
        print('Keyword "states" not given!')
        sys.exit(39)
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
        sys.exit(40)

    possibletasks = ['h', 'soc', 'dm', 'grad', 'overlap', 'dmdr', 'socdr', 'ion', 'phases']
    if not any([i in QMin for i in possibletasks]):
        print('No tasks found! Tasks are "h", "soc", "dm", "grad","dmdr", "socdr", "overlap" and "ion".')
        sys.exit(41)

    if 'samestep' in QMin and 'init' in QMin:
        print('"Init" and "Samestep" cannot be both present in QM.in!')
        sys.exit(42)

    if 'phases' in QMin:
        QMin['overlap'] = []

    if 'overlap' in QMin and 'init' in QMin:
        print('"overlap" and "phases" cannot be calculated in the first timestep! Delete either "overlap" or "init"')
        sys.exit(43)

    if 'init' not in QMin and 'samestep' not in QMin:
        QMin['newstep'] = []

    if not any([i in QMin for i in ['h', 'soc', 'dm', 'grad']]) and 'overlap' in QMin:
        QMin['h'] = []

    if len(QMin['states']) > 8:
        print('Higher multiplicities than octets are not supported!')
        sys.exit(44)

    if 'h' in QMin and 'soc' in QMin:
        QMin = removekey(QMin, 'h')

    if 'nacdt' in QMin:
        print('Within the SHARC-MOLCAS interface, "nacdt" is not supported.')
        sys.exit(45)

    if 'molden' in QMin:
        os.environ['MOLCAS_MOLDEN'] = 'ON'
        if 'samestep' in QMin:
            print('HINT: Not producing Molden files in "samestep" mode!')
            del QMin['molden']

    # if 'ion' in QMin:
        # print('Ionization probabilities not implemented!')
        # sys.exit(46)

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

    # Process the overlap requests
    # identically to the nac requests
    if 'overlap' in QMin:
        if len(QMin['overlap']) >= 1:
            nacpairs = QMin['overlap']
            for i in range(len(nacpairs)):
                if nacpairs[i][0] > nmstates or nacpairs[i][1] > nmstates:
                    print('State for requested non-adiabatic couplings does not correspond to any state in QM input file state list!')
                    sys.exit(49)
        else:
            QMin['overlap'] = [[j + 1, i + 1] for i in range(nmstates) for j in range(i + 1)]

    # Process the non-adiabatic coupling requests
    # type conversion has already been done
    if 'nacdr' in QMin:
        if len(QMin['nacdr']) >= 1:
            nacpairs = QMin['nacdr']
            for i in range(len(nacpairs)):
                if nacpairs[i][0] > nmstates or nacpairs[i][1] > nmstates:
                    print('State for requested non-adiabatic couplings does not correspond to any state in QM input file state list!')
                    sys.exit(50)
        else:
            QMin['nacdr'] = [[j + 1, i + 1] for i in range(nmstates) for j in range(i)]


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







    # open MOLCAS.resources
    filename = 'MOLCAS.resources'
    if os.path.isfile(filename):
        sh2cas = readfile(filename)
    else:
        print('Warning: No MOLCAS.resources found!')
        print('Reading resources from SH2CAS.inp')
        sh2cas = readfile('SH2CAS.inp')

    QMin['pwd'] = os.getcwd()

    QMin['molcas'] = get_sh2cas_environ(sh2cas, 'molcas')
    os.environ['MOLCAS'] = QMin['molcas']

    driver = get_sh2cas_environ(sh2cas, 'driver', crucial=False)
    if driver == '':
        driver = os.path.join(QMin['molcas'], 'bin', 'pymolcas')
        if not os.path.isfile(driver):
            driver = os.path.join(QMin['molcas'], 'bin', 'molcas.exe')
            if not os.path.isfile(driver):
                print('No driver (pymolcas or molcas.exe) found in $MOLCAS/bin. Please add the path to the driver via the "driver" keyword.')
                sys.exit(52)
    QMin['driver'] = driver

    QMin['tinker'] = get_sh2cas_environ(sh2cas, 'tinker', crucial=False)
    if QMin['tinker'] == '':
        QMin = removekey(QMin, 'tinker')
    else:
        os.environ['TINKER'] = QMin['tinker']

    if 'ion' in QMin:
        QMin['wfoverlap'] = get_sh2cas_environ(sh2cas, 'wfoverlap', crucial=False)
        if not QMin['wfoverlap']:
            ciopath = os.path.join(os.path.expandvars(os.path.expanduser('$SHARC')), 'wfoverlap.x')
            if os.path.isfile(ciopath):
                QMin['wfoverlap'] = ciopath
            else:
                print('Give path to wfoverlap.x in MOLCAS.resources!')
                sys.exit(51)
        # read ncore and ndocc from resources
        line = getsh2caskey(sh2cas, 'numfrozcore')
        if line[0]:
            try:
                QMin['ncore'] = max(0, int(line[1]))
            except ValueError:
                print('numfrozcore does not evaluate to numerical value!')
                sys.exit(52)
        line = getsh2caskey(sh2cas, 'numocc')
        if line[0]:
            try:
                QMin['ndocc'] = int(line[1])
            except ValueError:
                print('numocc does not evaluate to numerical value!')
                sys.exit(53)


    # Set up scratchdir
    line = get_sh2cas_environ(sh2cas, 'scratchdir', False, False)
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
        line = get_sh2cas_environ(sh2cas, 'savedir', False, False)
        if line is None or line == '':
            line = os.path.join(QMin['pwd'], 'SAVEDIR')
    line = os.path.expandvars(line)
    line = os.path.expanduser(line)
    line = os.path.abspath(line)
    if 'init' in QMin:
        checkscratch(line)
    QMin['savedir'] = line


    line = getsh2caskey(sh2cas, 'debug')
    if line[0]:
        if len(line) <= 1 or 'true' in line[1].lower():
            global DEBUG
            DEBUG = True

    line = getsh2caskey(sh2cas, 'no_print')
    if line[0]:
        if len(line) <= 1 or 'true' in line[1].lower():
            global PRINT
            PRINT = False

    QMin['memory'] = 500
    line = getsh2caskey(sh2cas, 'memory')
    if line[0]:
        try:
            QMin['memory'] = int(line[1])
        except ValueError:
            print('MOLCAS memory does not evaluate to numerical value!')
            sys.exit(54)
    else:
        print('WARNING: Please set memory for MOLCAS in MOLCAS.resources (in MB)! Using 500 MB default value!')
    os.environ['MOLCASMEM'] = str(QMin['memory'])
    os.environ['MOLCAS_MEM'] = str(QMin['memory'])

    QMin['ncpu'] = 1
    line = getsh2caskey(sh2cas, 'ncpu')
    if line[0]:
        try:
            QMin['ncpu'] = int(line[1])
        except ValueError:
            print('Number of CPUs does not evaluate to numerical value!')
            sys.exit(55)

    QMin['mpi_parallel'] = False
    line = getsh2caskey(sh2cas, 'mpi_parallel')
    if line[0]:
        QMin['mpi_parallel'] = True


    QMin['schedule_scaling'] = 0.6
    line = getsh2caskey(sh2cas, 'schedule_scaling')
    if line[0]:
        try:
            x = float(line[1])
            if 0 < x <= 2.:
                QMin['schedule_scaling'] = x
        except ValueError:
            print('Schedule scaling does not evaluate to numerical value!')
            sys.exit(56)

    QMin['Project'] = 'MOLCAS'
    os.environ['Project'] = QMin['Project']

    QMin['delay'] = 0.0
    line = getsh2caskey(sh2cas, 'delay')
    if line[0]:
        try:
            QMin['delay'] = float(line[1])
        except ValueError:
            print('Submit delay does not evaluate to numerical value!')
            sys.exit(57)

    QMin['Project'] = 'MOLCAS'
    os.environ['Project'] = QMin['Project']
    os.environ['MOLCAS_OUTPUT'] = 'PWD'

    line = getsh2caskey(sh2cas, 'always_orb_init')
    if line[0]:
        QMin['always_orb_init'] = []
    line = getsh2caskey(sh2cas, 'always_guess')
    if line[0]:
        QMin['always_guess'] = []
    if 'always_orb_init' in QMin and 'always_guess' in QMin:
        print('Keywords "always_orb_init" and "always_guess" cannot be used together!')
        sys.exit(58)

    # open template
    template = readfile('MOLCAS.template')

    QMin['template'] = {}
    integers = ['nactel', 'inactive', 'ras2', 'frozen', 'ras1', 'ras3']
    strings = ['basis', 'method', 'baslib', 'pdft-functional']
    floats = ['ipea', 'imaginary', 'gradaccumax', 'gradaccudefault', 'displ', 'rasscf_thrs_e', 'rasscf_thrs_rot', 'rasscf_thrs_egrd', 'cholesky_accu']
    booleans = ['cholesky', 'no-douglas-kroll', 'douglas-kroll', 'qmmm', 'cholesky_analytical', 'diab_num_grad', 'cobramm', 'caspt2_ana_grad']
    for i in booleans:
        QMin['template'][i] = False
    QMin['template']['roots'] = [0 for i in range(8)]
    QMin['template']['rootpad'] = [0 for i in range(8)]
    QMin['template']['ras1'] = 0
    QMin['template']['ras3'] = 0

    QMin['template']['method'] = 'casscf'
    QMin['template']['pdft-functional'] = 't:pbe'
    QMin['template']['baslib'] = ''
    QMin['template']['ipea'] = 0.25
    QMin['template']['imaginary'] = 0.00
    QMin['template']['frozen'] = -1
    QMin['template']['iterations'] = [200, 100]
    QMin['template']['gradaccumax'] = 1.e-2
    QMin['template']['gradaccudefault'] = 1.e-4
    QMin['template']['displ'] = 0.005
    QMin['template']['cholesky_accu'] = 1e-4
    QMin['template']['rasscf_thrs_e'] = 1e-8
    QMin['template']['rasscf_thrs_rot'] = 1e-4            # TODO: apparent default in MOLCAS is 0.1
    QMin['template']['rasscf_thrs_egrd'] = 1e-4
    QMin['template']['pcmset'] = {'solvent': 'water', 'aare': 0.4, 'r-min': 1.0, 'on': False}
    QMin['template']['pcmstate'] = (QMin['statemap'][1][0], QMin['statemap'][1][1])
    QMin['template']['no-douglas-kroll'] = True
    QMin['template']['caspt2_ana_grad'] = False

    for line in template:
        orig = re.sub('#.*$', '', line).split(None, 1)
        line = re.sub('#.*$', '', line).lower().split()
        if len(line) == 0:
            continue
        if 'spin' in line[0]:
            QMin['template']['roots'][int(line[1]) - 1] = int(line[3])
        elif 'roots' in line[0]:
            for i, n in enumerate(line[1:]):
                QMin['template']['roots'][i] = int(n)
        elif 'nactel' in line[0]:
            QMin['template']['nactel'] = [int(i) for i in line[1:]]
        elif 'rootpad' in line[0]:
            for i, n in enumerate(line[1:]):
                QMin['template']['rootpad'][i] = int(n)
        elif 'baslib' in line[0]:
            QMin['template']['baslib'] = os.path.abspath(orig[1])
        elif 'iterations' in line[0]:
            if len(line) >= 3:
                QMin['template']['iterations'] = [int(i) for i in line[-2:]]
            elif len(line) == 2:
                QMin['template']['iterations'][0] = int(line[-1])
        elif line[0] in integers:
            QMin['template'][line[0]] = int(line[1])
        elif line[0] in booleans:
            QMin['template'][line[0]] = True
        elif line[0] in strings:
            QMin['template'][line[0]] = line[1]
        elif line[0] in floats:
            QMin['template'][line[0]] = float(line[1])
        elif 'pcmset' in line[0]:
            # order: solvent, aare, r-min
            QMin['template']['pcmset']['on'] = True
            QMin['template']['pcmset']['solvent'] = line[1]
            if len(line) >= 3:
                QMin['template']['pcmset']['aare'] = float(line[2])
            if len(line) >= 4:
                QMin['template']['pcmset']['r-min'] = float(line[3])
        elif 'pcmstate' in line[0]:
            QMin['template']['pcmstate'] = (int(line[1]), int(line[2]))

    # roots must be larger or equal to states
    for i, n in enumerate(QMin['template']['roots']):
        if i == len(QMin['states']):
            break
        if not n >= QMin['states'][i]:
            print('Too few states in state-averaging in multiplicity %i! %i requested, but only %i given' % (i + 1, QMin['states'][i], n))
            sys.exit(59)

    # check nactel
    if len(QMin['template']['nactel']) == 1:
        QMin['template']['nactel'] = [ QMin['template']['nactel'][0], 0, 0 ]
    elif len(QMin['template']['nactel']) == 3:
        pass
    else:
        print('You must specify either 1 or 3 numbers for "nactel"!')
        sys.exit(59)


    # check rootpad
    for i, n in enumerate(QMin['template']['rootpad']):
        if i == len(QMin['states']):
            break
        if not n >= 0:
            print('Rootpad must not be negative!')
            sys.exit(60)

    # condense roots list
    for i in range(len(QMin['template']['roots']) - 1, 0, -1):
        if QMin['template']['roots'][i] == 0:
            QMin['template']['roots'].pop(i)
        else:
            break
    QMin['template']['rootpad'] = QMin['template']['rootpad'][:len(QMin['template']['roots'])]

    # check roots versus number of electrons
    #nelec=QMin['template']['inactive']*2+QMin['template']['nactel']
    #for i,n in enumerate(QMin['states']):
    #    if n>0:
    #        if not (QMin['template']['nactel']+i)%2==0:
    #            print('Number of electrons is %i, but states of multiplicity %i are requested.' % (nelec,i+1))
    #            sys.exit(61)

    necessary = ['basis', 'nactel', 'ras2', 'inactive']
    for i in necessary:
        if i not in QMin['template']:
            print('Key %s missing in template file!' % (i))
            sys.exit(62)


    # logic checks:
    if QMin['template']['pcmset']['on']:
        if QMin['template']['qmmm']:
            print('PCM and QM/MM cannot be used together!')

    if QMin['template']['pcmset']['on']:
        if QMin['template']['cobramm']:
            print('PCM and QM/MM cannot be used together!')

    # modern OpenMolcas can do analytical Cholesky gradients/NACs
    QMin['template']['cholesky_analytical']=QMin['template']['cholesky']

    # Douglas-Kroll new treatment
    QMin['template']['no-douglas-kroll'] = not QMin['template']['douglas-kroll']

    # QM/MM mode needs Tinker
    if QMin['template']['qmmm'] and 'tinker' not in QMin:
        print('Please set $TINKER or give path to tinker in MOLCAS.resources!')
        sys.exit(63)

    # find method
    allowed_methods = ['casscf', 'caspt2', 'ms-caspt2', 'mc-pdft', 'xms-pdft', 'cms-pdft', 'xms-caspt2']
    # 0: casscf
    # 1: caspt2 (single state)
    # 2: ms-caspt2
    # 3: mc-pdft (single state)
    # 4: xms-pdft
    # 5: cms-pdft
    # 6: xms-caspt2
    for i, m in enumerate(allowed_methods):
        if QMin['template']['method'] == m:
            QMin['method'] = i
            break
    else:
        print('Unknown method "%s" given in MOLCAS.template' % (QMin['template']['method']))
        sys.exit(64)

    # find functional if it is cms-pdft
    if QMin['method'] == 5:
        allowed_functionals=['tpbe','t:pbe','ft:pbe','t:blyp','ft:blyp','t:revPBE','ft:revPBE','t:LSDA','ft:LSDA']
        for i,m in enumerate(allowed_functionals):
            if QMin['template']['pdft-functional'] == m:
                QMin['pdft-functional'] = i
                break
        else:
            print('Warning! No analytical gradients for cms-pdft and "%s" given in MOLCAS.template' % (QMin['template']['pdft-functional']))
            print('Using numerical gradients. Analytical gradients only for t:pbe, ft:pbe, t:blyp, ft:blyp, t:revPBE, ft:revPBE, t:LSDA, or ft:LSDA.')
            QMin['pdft-functional'] = -1

    # decide which type of gradients to do:
    # 0 = analytical CASSCF gradients in one MOLCAS input file (less overhead, but recommended only under certain circumstances)
    # 1 = analytical CASSCF gradients in separate MOLCAS inputs, possibly distributed over several CPUs (DEFAULT)
    # 2 = numerical gradients (CASPT2, MS-CASPT2, Cholesky-CASSCF; or for dmdr and socdr), possibly distributed over several CPUs
    # 3 = analytical CMS-PDFT gradients in one MOLCAS input file (less overhead, but recommended only under certain circumstances)
    if 'dmdr' in QMin or 'socdr' in QMin or 'grad' in QMin or 'nacdr' in QMin:
        if 'dmdr' in QMin or 'socdr' in QMin:
            print('Numerical gradients due to dmdr or socdr...')
            QMin['gradmode'] = 2
        elif QMin['template']['cholesky'] and not QMin['template']['cholesky_analytical']:
            print('Numerical gradients due to old Cholesky code...')
            QMin['gradmode'] = 2
        elif QMin['template']['pcmset']['on']:
            print('Numerical gradients due to PCM...')
            QMin['gradmode'] = 2
        elif QMin['method'] in [1, 2, 4, 6] and not QMin['template']['caspt2_ana_grad']:
            print('Numerical gradients due to SS-CASPT2, MS-CASPT2, or XMS-PDFT...')
            QMin['gradmode'] = 2
        elif QMin['method'] == 5 and QMin['pdft-functional'] == -1:
                print('Numerical gradients due to XMS-PDFT with non-compatible functional...')
                QMin['gradmode'] = 2
        else:
            if QMin['ncpu'] > 0:
                QMin['gradmode'] = 1
            else:
                # check if any gradient to be calculated is a SS-CASSCF gradient
                ss_grads = False
                for i in QMin['gradmap']:
                    if QMin['template']['roots'][i[0] - 1] == 1:
                        ss_grads = True
                if ss_grads:
                    QMin['gradmode'] = 1
                else:
                    QMin['gradmode'] = 0
        if QMin['gradmode'] == 2:
            QMin['displ'] = QMin['template']['displ'] / au2a
    else:
        QMin['gradmode'] = 0
    QMin['ncpu'] = max(1, QMin['ncpu'])

    # currently, QM/MM is only allowed with CASSCF and analytical gradients on one CPU
    # will be available in the future
    if QMin['template']['qmmm'] and QMin['gradmode'] == 2:
        print('QM/MM is only possible currently with analytical gradients.')
        sys.exit(65)
    if 'nacdr' in QMin and QMin['gradmode'] == 2:
        print('Nonadiabatic coupling vectors are only possible currently with analytical gradients.')
        sys.exit(66)
    if QMin['template']['qmmm'] and 'nacdr' in QMin:
        print('Nonadiabatic coupling vectors are currently not available with QM/MM.')
        sys.exit(67)
    if 'nacdr' in QMin and QMin['method'] not in [0,2,6]:
        print('Nonadiabatic coupling vectors are only possible currently with RASSCF and MS/XMS-CASPT2.')
        sys.exit(66)



    # QM/MM Tinker setup
    if QMin['template']['qmmm']:
        QMin['active_qmmm_atoms'] = []
        keycontent = readfile('MOLCAS.qmmm.key')
        for line in keycontent:
            if 'qmmm' in line.lower() and 'electrostatics' not in line.lower():
                QMin['total_qmmm_natom'] = int(line.split()[1])
            elif containsstring('^qm ', line.lower()) or containsstring('^mm ', line.lower()):
                line = line.split()
                if len(line) == 3 and int(line[1]) < 0:     # range definition (e.g. QM -1 4)
                    start = -int(line[1])
                    end = int(line[2])
                    QMin['active_qmmm_atoms'] += [i for i in range(start, end + 1)]
                else:                                   # list definition (e.g. MM 5 6 7 8)
                    QMin['active_qmmm_atoms'] += [int(i) for i in line[1:]]
        # print('total amount of qmmm atoms given in key file:', QMin['total_qmmm_natom'])
        # print('number of indices given in key file:', len(QMin['active_qmmm_atoms']))


    # Check the save directory
    try:
        ls = os.listdir(QMin['savedir'])
        err = 0
    except OSError:
        print('Problems reading SCRADIR=%s' % (QMin['savedir']))
        sys.exit(68)
    if 'init' in QMin:
        err = 0
    elif 'samestep' in QMin:
        for imult, nstates in enumerate(QMin['states']):
            if nstates < 1:
                continue
            if not 'MOLCAS.%i.JobIph' % (imult + 1) in ls:
                print('File "MOLCAS.%i.JobIph" missing in SAVEDIR!' % (imult + 1))
                err += 1
        if 'overlap' in QMin:
            for imult, nstates in enumerate(QMin['states']):
                if nstates < 1:
                    continue
                if not 'MOLCAS.%i.JobIph.old' % (imult + 1) in ls:
                    print('File "MOLCAS.%i.JobIph.old" missing in SAVEDIR!' % (imult + 1))
                    err += 1
    elif 'overlap' in QMin:
        for imult, nstates in enumerate(QMin['states']):
            if nstates < 1:
                continue
            if not 'MOLCAS.%i.JobIph' % (imult + 1) in ls:
                print('File "MOLCAS.%i.JobIph" missing in SAVEDIR!' % (imult + 1))
                err += 1
    if err > 0:
        print('%i files missing in SAVEDIR=%s' % (err, QMin['savedir']))
        sys.exit(69)

    QMin['version'] = getversion([''] * 50, QMin['molcas'])

    if PRINT:
        printQMin(QMin)

    return QMin


# =============================================================================================== #
# =============================================================================================== #
# =========================================== gettasks and setup routines ======================= #
# =============================================================================================== #
# =============================================================================================== #

def gettasks(QMin):

    # Currently implemented keywords: soc, dm, grad,
    tasks = []
    if 'pargrad' not in QMin:
        tasks.append(['gateway'])
        tasks.append(['seward'])

    if QMin['template']['qmmm']:
        if 'pargrad' in QMin:
            tasks.append(['gateway'])
            tasks.append(['seward'])
        tasks.append(['espf'])

    if QMin['template']['pcmset']['on']:
        pcm_mult = QMin['template']['pcmstate'][0]
        list_to_do = [(pcm_mult - 1, QMin['states'][pcm_mult - 1])]
        for i, j in enumerate(QMin['states']):
            if not i + 1 == pcm_mult:
                list_to_do.append((i, j))
    else:
        list_to_do = [(i, j) for i, j in enumerate(QMin['states'])]

    # for imult,nstates in enumerate(QMin['states']):
    for imult, nstates in list_to_do:
        if nstates == 0:
            continue

        # find the correct initial MO file
        mofile = ''
        if 'always_guess' not in QMin:
            if 'init' in QMin or 'always_orb_init' in QMin:
                ls = os.listdir(QMin['pwd'])
                for i in ls:
                    if 'MOLCAS.%i.JobIph.init' % (imult + 1) in i:
                        mofile = os.path.join(QMin['pwd'], 'MOLCAS.%i.JobIph.init' % (imult + 1))
                        break
                    elif 'MOLCAS.%i.RasOrb.init' % (imult + 1) in i:
                        mofile = os.path.join(QMin['pwd'], 'MOLCAS.%i.RasOrb.init' % (imult + 1))
                        break
            elif 'samestep' in QMin:
                mofile = os.path.join(QMin['savedir'], 'MOLCAS.%i.JobIph' % (imult + 1))
            elif 'displacement' in QMin:
                mofile = os.path.join(QMin['savedir'], 'MOLCAS.%i.JobIph.master' % (imult + 1))
            else:
                mofile = os.path.join(QMin['savedir'], 'MOLCAS.%i.JobIph.old' % (imult + 1))
        if not mofile == '':
            if 'JobIph' in mofile:
                tasks.append(['link', mofile, 'JOBOLD'])
            elif 'RasOrb' in mofile:
                tasks.append(['link', mofile, 'INPORB'])

        # RASSCF
        # ['rasscf',imult,sa-nstates,jobiph]
        if 'samestep' not in QMin or 'always_orb_init' in QMin:
            jobiph = 'JobIph' in mofile
            rasorb = 'RasOrb' in mofile
            if QMin['method']==5 and QMin['pdft-functional'] > -1:
                if not 'init' in QMin:
                    tasks.append(['copy',os.path.join(QMin['savedir'],'Do_Rotate.%i.txt' % (imult+1)),'Do_Rotate.txt'])
            tasks.append(['rasscf', imult + 1, QMin['template']['roots'][imult], jobiph, rasorb])
            if QMin['method'] == 4:
                tasks[-1].append(['XMSI'])
            if QMin['method'] == 5:
                tasks[-1].append(['CMSI'])
            # next lines make no sense ~~~
            if jobiph and not (QMin['method'] == 5 and QMin['pdft-functional'] > -1):
                tasks.append(['rm', 'JOBOLD'])
            if QMin['method'] == 0 or QMin['method'] == 3:
                tasks.append(['copy', 'MOLCAS.JobIph', 'MOLCAS.%i.JobIph' % (imult + 1)])
            elif QMin['method'] == 5 and QMin['pdft-functional'] > -1:
                tasks.append(['copy','MOLCAS.JobIph','MOLCAS.%i.JobIph' % (imult+1)])
                if not 'init' in QMin:
                    tasks.append(['rm', 'JOBOLD'])
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if 'ion' in QMin:
                tasks.append(['copy', 'MOLCAS.RasOrb', 'MOLCAS.%i.RasOrb' % (imult + 1)])
            if 'molden' in QMin:
                tasks.append(['copy', 'MOLCAS.rasscf.molden', 'MOLCAS.%i.molden' % (imult + 1)])

            if QMin['method'] in [3, 4, 5]:
                keys = ['KSDFT=%s' % QMin['template']['pdft-functional']]
                if QMin['gradmode'] == 0 and 'grad' in QMin:
                    keys.append('GRAD')
                else:
                    keys.append('noGrad')
                if QMin['method'] in [4, 5]:
                    keys.append('MSPDFT')
                    keys.append('WJOB')
                if QMin['method'] in [5]:
                    keys.append('CMMI=0')
                    keys.append('CMSS=Do_Rotate.txt')
                    keys.append('CMTH=1.0d-10')
                tasks.append(['mcpdft', keys])
                if QMin['method'] in [4, 5]:
                    tasks.append(['copy', 'MOLCAS.JobIph', 'MOLCAS.%i.JobIph' % (imult + 1)])


        if 'samestep' not in QMin:
            if QMin['method'] in [1,2,6]:
                # caspt2
                tasks.append(['caspt2', imult + 1, nstates, QMin['method'] ])
                # copy JobIphs
                tasks.append(['copy', 'MOLCAS.JobMix', 'MOLCAS.%i.JobIph' % (imult + 1)])

#            if QMin['method'] in [3, 4]:
#                # mc-pdft
#                keys = ['KSDFT=%s' % QMin['template']['pdft-functional']]
#                if QMin['method'] == 3 and 'grad' in QMin:
#                    keys.append('GRAD')
#                else:
#                    keys.append('noGrad')
#                if QMin['method'] == 4:
#                    keys.append('MSPDFT')
#                    keys.append('WJOB')
#                tasks.append(['mcpdft', keys])
#                # copy JobIphs
#                if QMin['method'] == 4:
#                    tasks.append(['copy', 'MOLCAS.JobIph', 'MOLCAS.%i.JobIph' % (imult + 1)])
#
#            if QMin['method'] == 5:
#                # mc-pdft
#                keys = ['KSDFT=%s' % QMin['template']['pdft-functional']]
#                if 'grad' in QMin:
#                    keys.append('GRAD')
#                else:
#                    keys.append('noGrad')
#                keys.append('MSPDFT')
#                keys.append('WJOB')
#                tasks.append(['mcpdft', keys])
#                # copy JobIphs
#                tasks.append(['copy', 'MOLCAS.JobIph', 'MOLCAS.%i.JobIph' % (imult + 1)])

        # Gradients
        if QMin['gradmode'] == 0:
            for i in QMin['gradmap']:
                if i[0] == imult + 1:
                    if QMin['template']['roots'][imult] == 1:
                        # SS-CASSCF
                        if 'samestep' in QMin:
                            tasks.append(['rasscf', imult + 1, QMin['template']['roots'][imult], True, False])
                            if QMin['method'] == 3:
                                tasks.append(['mcpdft', ['KSDFT=%s' % QMin['template']['pdft-functional'], 'GRAD']])
                            if QMin['method'] in [2,6]:
                                print('Error: single state gradient with MS/XMS-CASPT2')
                                sys.exit(99)
                        tasks.append(['alaska'])
                    else:
                        # SA-CASSCF
                        tasks.append(['link', 'MOLCAS.%i.JobIph' % (imult + 1), 'JOBOLD'])
                        # tasks.append(['rasscf-rlx',imult+1,QMin['template']['roots'][imult],i[1]])
                        # tasks.append(['mclr',QMin['template']['gradaccudefault']])
                        # tasks.append(['alaska'])

                        if QMin['method'] == 3:
                            tasks.append(['rasscf', imult + 1, QMin['template']['roots'][imult], True, False, ['RLXROOT=%i' % i[1]]])
                            tasks.append(['mcpdft', ['KSDFT=%s' % QMin['template']['pdft-functional'], 'GRAD']])
                            tasks.append(['alaska'])
                        elif QMin['method'] == 5:
                            if not 'init' in QMin:
                                tasks.append(['copy',os.path.join(QMin['savedir'],'Do_Rotate.%i.txt' % (imult+1)),'Do_Rotate.txt'])
                            #tasks.append(['rasscf-cms',imult+1,QMin['template']['roots'][imult],True,False,i[1]])
                            tasks.append(['rasscf', imult + 1, QMin['template']['roots'][imult], True, False, ['RLXROOT=%i' % i[1], 'CMSI' ] ])
                            tasks.append(['mcpdft', ['KSDFT=%s' % QMin['template']['pdft-functional'], 'GRAD','MSPDFT','WJOB'] ])
                            tasks.append(['alaska', i[1] ])
                        elif QMin['method'] == 0:
                            tasks.append(['rasscf', imult + 1, QMin['template']['roots'][imult], True, False])
                            tasks.append(['mclr', QMin['template']['gradaccudefault'], 'sala=%i' % i[1]])
                            tasks.append(['alaska'])
                        elif QMin['method'] in [2,6]:
                            tasks.append(['rasscf', imult + 1, QMin['template']['roots'][imult], True, False])
                            tasks.append(['caspt2', imult + 1, nstates, QMin['method'], 'GRDT\nrlxroot = %i' % i[1]   ])
                            tasks.append(['mclr', QMin['template']['gradaccudefault']])
                            tasks.append(['alaska'])


        # Nac vectors
        if QMin['gradmode'] == 0:
            for i in QMin['nacmap']:
                if i[0] == imult + 1:
                    tasks.append(['link', 'MOLCAS.%i.JobIph' % (imult + 1), 'JOBOLD'])
                    if QMin['method'] == 0:
                        tasks.append(['rasscf', imult + 1, QMin['template']['roots'][imult], True, False])
                        tasks.append(['mclr', QMin['template']['gradaccudefault'], 'nac=%i %i' % (i[1], i[3])])
                        tasks.append(['alaska'])
                    elif QMin['method'] == 5:
                            if not 'init' in QMin:
                                tasks.append(['copy',os.path.join(QMin['savedir'],'Do_Rotate.%i.txt' % (imult+1)),'Do_Rotate.txt'])
                            #tasks.append(['rasscf-cms',imult+1,QMin['template']['roots'][imult],True,False,i[1]])
                            tasks.append(['rasscf', imult + 1, QMin['template']['roots'][imult], True, False, ['CMSI' ] ])
                            tasks.append(['mcpdft', ['KSDFT=%s' % QMin['template']['pdft-functional'], 'GRAD','MSPDFT','WJOB'] ])
                            tasks.append(['mclr', QMin['template']['gradaccudefault'], 'nac=%i %i' % (i[1], i[3])])
                            tasks.append(['alaska' ])
                    elif QMin['method'] in [2,6]:
                        tasks.append(['rasscf', imult + 1, QMin['template']['roots'][imult], True, False])
                        tasks.append(['caspt2', imult + 1, nstates, QMin['method'], 'GRDT\nnac = %i %i' % (i[1],i[3])   ])
                        #tasks.append(['mclr', QMin['template']['gradaccudefault']])
                        tasks.append(['alaska', i[1], i[3] ])



#        #CMS-PDFT Gradients
#        if QMin['method'] == 5 and QMin['pdft-functional'] > -1:
#            for i in QMin['gradmap']:
#                if i[0]==imult+1:
#                    if not 'init' in QMin:
#                        tasks.append(['copy',os.path.join(QMin['savedir'],'Do_Rotate.%i.txt' % (imult+1)),'Do_Rotate.txt'])
#                    #tasks.append(['link','MOLCAS.%i.JobIph' % (imult+1),'JOBOLD'])
#                    tasks.append(['rasscf-cms',imult+1,QMin['template']['roots'][imult],True,False,i[1]])
#                    tasks.append(['cms-pdft',QMin['template']['pdft-functional']])
#                    tasks.append(['copy','Do_Rotate.txt',os.path.join(QMin['savedir'],'Do_Rotate.%i.txt' % (imult+1))])
#                    tasks.append(['mclr-cms',QMin['template']['gradaccudefault']])
#                    tasks.append(['alaska'])

        # RASSI for overlaps
        if 'overlap' in QMin:
            if 'overlap_nacs' in QMin:
                if imult+1 not in [ i[0] for i in QMin['nacmap']]:
                    continue
                tasks.append(['link', 'MOLCAS.%i.JobIph' % (imult + 1), 'JOB001'])
                tasks.append(['link', 'MOLCAS.JobIph', 'JOB002'])
            else:
                if 'displacement' in QMin:
                    tasks.append(['link', os.path.join(QMin['savedir'], 'MOLCAS.%i.JobIph.master' % (imult + 1)), 'JOB001'])
                else:
                    tasks.append(['link', os.path.join(QMin['savedir'], 'MOLCAS.%i.JobIph.old' % (imult + 1)), 'JOB001'])
                tasks.append(['link', 'MOLCAS.%i.JobIph' % (imult + 1), 'JOB002'])
            tasks.append(['rassi', 'overlap', [nstates, nstates]])

        # RASSI for Dipole moments only if overlap-RASSI is not needed
        elif 'dm' in QMin or 'ion' in QMin:
            tasks.append(['link', 'MOLCAS.%i.JobIph' % (imult + 1), 'JOB001'])
            tasks.append(['rassi', 'dm', [nstates]])


    # SOC
    if 'soc' in QMin:
        i = 0
        roots = []
        for imult, nstates in enumerate(QMin['states']):
            if nstates == 0:
                continue
            i += 1
            roots.append(nstates)
            tasks.append(['link', 'MOLCAS.%i.JobIph' % (imult + 1), 'JOB%03i' % (i)])
        tasks.append(['rassi', 'soc', roots])

    if DEBUG:
        printtasks(tasks)

    return tasks

# ======================================================================= #


def writeMOLCASinput(tasks, QMin):

    string = ''

    for task in tasks:

        if task[0] == 'gateway':
            if QMin['template']['qmmm']:
                string += '&GATEWAY\nTINKER\nGROUP=NOSYM\nBASIS=%s\n' % (QMin['template']['basis'])
            else:
                string += '&GATEWAY\nCOORD=MOLCAS.xyz\nGROUP=NOSYM\nBASIS=%s\n' % (QMin['template']['basis'])
            if 'soc' in QMin:
                string += 'AMFI\n'
                string += 'angmom\n0 0 0\n'
            if QMin['template']['baslib']:
                string += 'BASLIB\n%s\n' % QMin['template']['baslib']
            if QMin['template']['cholesky']:
                string += 'RICD\nCDTHreshold=%f\n' % (QMin['template']['cholesky_accu'])
            if QMin['template']['pcmset']['on']:
                string += 'RF-INPUT\nPCM-MODEL\nSOLVENT = %s\nAARE = %f\nR-MIN = %f\nEND OF RF-INPUT\n' % (QMin['template']['pcmset']['solvent'], QMin['template']['pcmset']['aare'], QMin['template']['pcmset']['r-min'])
            # if 'cobramm' in QMin['template']:
            if QMin['template']['cobramm']:
                with open("charge.dat", "r") as fin:
                    string += fin.read() + "\n"
            string += '\n'

        elif task[0] == 'seward':
            string += '&SEWARD\n'
            if not QMin['template']['no-douglas-kroll']:
                string+='Relativistic = R02O02\nRELINT\n'
            # if 'soc' in QMin and QMin['version']<=8.0:
                # string+='AMFI\n'
            # if 'soc' in QMin:
                # string+='AMFI\n'
            if QMin['template']['cholesky']:
                string += 'DOANA\n'
            #if it is cms-pdft, add CMSI
            if QMin['method'] in [3,4,5]:
#            if QMin['method'] == 5 and QMin['pdft-functional'] > -1:
                string += 'GRID INPUT\n'
#                string += 'GRID=ULTRAFINE\n'
                string += 'NORO\n'
                string += 'NOSC\n'
                string += 'END OF GRID INPUT\n'
            string += '\n'

        elif task[0] == 'espf':
            string += '&ESPF\nEXTERNAL=TINKER\n\n'

        elif task[0] == 'link':
            name = os.path.basename(task[1])
            string += '>> COPY %s %s\n\n' % (name, task[2])

        elif task[0] == 'copy':
            string += '>> COPY %s %s\n\n' % (task[1], task[2])

        elif task[0] == 'rm':
            string += '>> RM %s\n\n' % (task[1])

        elif task[0] == 'rasscf':
            nactel = QMin['template']['nactel'][0]
            nhole = QMin['template']['nactel'][1]
            nelec = QMin['template']['nactel'][2]
            npad = QMin['template']['rootpad'][task[1] - 1]
            if (nactel - task[1]) % 2 == 0:
                nactel -= 1
            string += '&RASSCF\nSPIN=%i\nNACTEL=%i %i %i\nINACTIVE=%i\nRAS2=%i\nITERATIONS=%i,%i\n' % (
                task[1],
                nactel,
                nhole,
                nelec,
                QMin['template']['inactive'],
                QMin['template']['ras2'],
                QMin['template']['iterations'][0],
                QMin['template']['iterations'][1])
            if QMin['template']['ras1'] > 0:
                string += 'RAS1=%i\n' % (QMin['template']['ras1'])
            if QMin['template']['ras3'] > 0:
                string += 'RAS3=%i\n' % (QMin['template']['ras3'])
            if npad == 0:
                string += 'CIROOT=%i %i 1\n' % (task[2], task[2])
            else:
                string += 'CIROOT=%i %i; ' % (task[2], task[2] + npad)
                for i in range(task[2]):
                    string += '%i ' % (i + 1)
                string += ';'
                for i in range(task[2]):
                    string += '%i ' % (1)
                string += '\n'
            if not QMin['method'] in [2,6]:
                string += 'ORBLISTING=NOTHING\nPRWF=0.1\n'
            if 'grad' in QMin and QMin['gradmode'] < 2:
                string += 'THRS=1.0e-10 1.0e-06 1.0e-06\n'
            else:
                string += 'THRS=%14.12f %14.12f %14.12f\n' % (QMin['template']['rasscf_thrs_e'], QMin['template']['rasscf_thrs_rot'], QMin['template']['rasscf_thrs_egrd'])
            if task[3]:
                string += 'JOBIPH\n'
            elif task[4]:
                string += 'LUMORB\n'
            if len(task) >= 6:
                for a in task[5]:
                    string += a + '\n'
            if QMin['template']['pcmset']['on']:
                if task[1] == QMin['template']['pcmstate'][0]:
                    string += 'RFROOT = %i\n' % QMin['template']['pcmstate'][1]
                else:
                    string += 'NONEQUILIBRIUM\n'
            string += '\n'

#   NOTE: CMS gradients are now evaluated in gradmode = 0
#        elif task[0]=='rasscf-cms':
#            nactel=QMin['template']['nactel']
#            npad=QMin['template']['rootpad'][task[1]-1]
#            if (nactel-task[1])%2==0:
#                nactel-=1
#            string+='&RASSCF\nSPIN=%i\nNACTEL=%i 0 0\nINACTIVE=%i\nRAS2=%i\n' % (
#                    task[1],
#                    nactel,
#                    QMin['template']['inactive'],
#                    QMin['template']['ras2'])
#            if npad==0:
#                string+='CIROOT=%i %i 1\n' % (task[2],task[2])
#            else:
#                string+='CIROOT=%i %i; ' % (task[2],task[2]+npad)
#                for i in range(task[2]):
#                    string+='%i ' % (i+1)
#                string+=';'
#                for i in range(task[2]):
#                    string+='%i ' % (1)
#                string+='\n'
#            string+='RLXROOT=%i\n' % (task[5])
#            string+='ORBLISTING=NOTHING\nPRWF=0.1\n'
#            if 'grad' in QMin and QMin['gradmode']<2:
#                string+='THRS=1.0e-10 1.0e-06 1.0e-06\n'
#            else:
#                string+='THRS=%14.12f %14.12f %14.12f\n' % (QMin['template']['rasscf_thrs_e'],QMin['template']['rasscf_thrs_rot'],QMin['template']['rasscf_thrs_egrd'])
#            if task[3]:
#                string+='JOBIPH\n'
#            elif task[4]:
#                string+='LUMORB\n'
#            if QMin['template']['pcmset']['on']:
#                if task[1]==QMin['template']['pcmstate'][0]:
#                    string+='RFROOT = %i\n' % QMin['template']['pcmstate'][1]
#                else:
#                    string+='NONEQUILIBRIUM\n'
#            string+='CMSI\n'
#            string+='CMMI=0\n'
#            if not 'init' in QMin:
#                string+='CMSS=Do_Rotate.txt\n'
#            string+='CMTH=1.0d-10\n'
#            string+='\n'
#
#            string+='\n'

        #elif task[0]=='rasscf-rlx':
            #nactel=QMin['template']['nactel']
            #if (nactel-task[1])%2==0:
                #nactel-=1
            #string+='&RASSCF\nSPIN=%i\nNACTEL=%i 0 0\nINACTIVE=%i\nRAS2=%i\nCIROOT=%i %i 1\nRLXROOT=%i\n' % (
                    #task[1],
                    #nactel,
                    #QMin['template']['inactive'],
                    #QMin['template']['ras2'],
                    #task[2],task[2],
                    #task[3])
            #string+='ORBLISTING=NOTHING\nPRWF=0.1\nJOBIPH\n'
            #string+='\n'

        elif task[0] == 'caspt2':
            string += '&CASPT2\nSHIFT=0.0\nIMAGINARY=%5.3f\nIPEASHIFT=%4.2f\nMAXITER=%i\n' % (
                QMin['template']['imaginary'],
                QMin['template']['ipea'],
                120)
            if QMin['template']['frozen'] != -1:
                string += 'FROZEN=%i\n' % (QMin['template']['frozen'])
            if QMin['method'] == 1:
                string += 'NOMULT\n'
            elif QMin['method'] == 2:
                string += 'MULTISTATE= %i ' % (task[2])
            elif QMin['method'] == 6:
                string += 'XMULTISTATE= %i ' % (task[2])
            if QMin['method'] in [2,6]:
                for i in range(task[2]):
                    string += '%i ' % (i + 1)
            string += '\nOUTPUT=BRIEF\nPRWF=0.1\n'
            if QMin['template']['pcmset']['on']:
                string += 'RFPERT\n'
            if len(task) == 5:
                string += task[4] 
            string += '\n'

#        elif task[0] == 'cms-pdft':
#            string += '&MCPDFT\nKSDFT=%s\n' % (task[1])
#            string += 'MSPD\n'
#            if 'grad' in QMin and QMin['gradmode'] == 3:
#                string += 'GRAD\n'
#            string += '\n'

        elif task[0] == 'mcpdft':
            string += '&MCPDFT\n'
            for i in task[1]:
                string += i + '\n'
            string += '\n\n'


        elif task[0] == 'rassi':
            string += '&RASSI\nNROFJOBIPHS\n%i' % (len(task[2]))
            for i in task[2]:
                string += ' %i' % (i)
            string += '\n'
            for i in task[2]:
                string += ' '.join([str(j) for j in range(1, i + 1)]) + '\n'
            string += 'MEIN\n'
            if QMin['method'] > 0:
                string += 'EJOB\n'
            if task[1] != 'soc' and 'master' in QMin and 'ion' in QMin:
                # smallest value printed by MOLCAS is 0.00001
                string += 'CIPR\nTHRS=0.000005d0\n'
            if task[1] == 'dm':
                pass
            elif task[1] == 'soc':
                string += 'SPINORBIT\nSOCOUPLING=0.0d0\nEJOB\n'
            elif task[1] == 'overlap':
                string += 'OVERLAPS\n'
            string += '\n'

        elif task[0] == 'mclr':
            string += '&MCLR\nTHRESHOLD=%f\n' % (task[1])
            if len(task) > 2:
                string += '%s\n' % (task[2])
            string += '\n'

#        elif task[0]=='mclr-cms':
#            string+='&MCLR\nTHRESHOLD=%f\n\n' % (task[1])

        elif task[0] == 'alaska':
            string += '&ALASKA\n'
            if len(task)==2:
                string+='root=%i\n' % (task[1])
            elif len(task)==3:
                string+='nac=%i %i\n' % (task[1], task[2])
            string+='\n'

        else:
            print('Unknown task keyword %s found in writeMOLCASinput!' % task[0])
            print(task)
            sys.exit(70)

    return string

# ======================================================================= #


def writegeomfile(QMin):
    string = ''
    if QMin['template']['qmmm']:
        try:
            geomtmpfile = open('MOLCAS.qmmm.table', 'r')
        except IOError:
            try:
                geomtmpfile = open('MOLCAS.qmmm.template', 'r')
            except IOError:
                print('Could not find file "MOLCAS.qmmm.table"!')
                sys.exit(71)
        geomtemplate = geomtmpfile.readlines()
        geomtmpfile.close()
        string += '%i\n\n' % (QMin['natom'])
        for iatom, atom in enumerate(QMin['geo']):
            tmpline = geomtemplate[iatom + 1].split()
            for xyz in range(3):
                tmpline[xyz + 2] = ' %f ' % (QMin['geo'][iatom][xyz + 1] * au2a)
            string += ' '.join(tmpline) + '\n'
    else:
        string += '%i\n\n' % (QMin['natom'])
        for iatom, atom in enumerate(QMin['geo']):
            string += '%s%i ' % (atom[0], iatom + 1)
            for xyz in range(1, 4):
                string += ' %f' % (atom[xyz] * au2a)
            string += '\n'

    return string

# ======================================================================= #


def setupWORKDIR(WORKDIR, tasks, QMin):
    # mkdir the WORKDIR, or clean it if it exists, then copy all necessary JobIph files from pwd and savedir
    # then put the geom.xyz and MOLCAS.input files

    # setup the directory
    if os.path.exists(WORKDIR):
        if os.path.isfile(WORKDIR):
            print('%s exists and is a file!' % (WORKDIR))
            sys.exit(72)
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
            sys.exit(73)

    # write geom file
    geomstring = writegeomfile(QMin)
    filename = os.path.join(WORKDIR, 'MOLCAS.xyz')
    writefile(filename, geomstring)
    if DEBUG:
        print(geomstring)
        print('Geom written to: %s' % (filename))

    # write MOLCAS.input
    inputstring = writeMOLCASinput(tasks, QMin)
    filename = os.path.join(WORKDIR, 'MOLCAS.input')
    writefile(filename, inputstring)
    if DEBUG:
        print(inputstring)
        print('MOLCAS input written to: %s' % (filename))

    # make subdirs
    if QMin['mpi_parallel']:
        for i in range(QMin['ncpu'] - 1):
            subdir = os.path.join(WORKDIR, 'tmp_%i' % (i + 1))
            os.makedirs(subdir)

    # JobIph copying
    copyfiles = set()
    for task in tasks:
        if task[0] == 'link' and task[1][0] == '/':
            copyfiles.add(task[1])
    for files in copyfiles:
        if DEBUG:
            print('Copy:\t%s\n\t==>\n\t%s' % (files, WORKDIR))
        shutil.copy(files, WORKDIR)
        if QMin['mpi_parallel']:
            for i in range(QMin['ncpu'] - 1):
                subdir = os.path.join(WORKDIR, 'tmp_%i' % (i + 1))
                if DEBUG:
                    print('Copy:\t%s\n\t==>\n\t%s' % (files, subdir))
                shutil.copy(files, subdir)

    # copy QM/MM related files
    if QMin['template']['qmmm']:
        fromfile = os.path.join(QMin['pwd'], 'MOLCAS.qmmm.key')
        tofile = os.path.join(WORKDIR, 'MOLCAS.key')
        # read in and expand environment variables
        content = readfile(fromfile)
        # print(content)
        for i in range(len(content)):
            if 'parameters' in content[i]:
                s = content[i].split()
                s[1] = os.path.abspath(os.path.expanduser(os.path.expandvars(s[1])))
                content[i] = ' '.join(s)
        if DEBUG:
            print('Copy:\t%s\n\t==>\n\t%s' % (fromfile, tofile))
            print(content)
        # shutil.copy(fromfile,tofile)
        writefile(tofile, content)
        if QMin['mpi_parallel']:
            for i in range(QMin['ncpu'] - 1):
                subdir = os.path.join(WORKDIR, 'tmp_%i' % (i + 1))
                tofile = os.path.join(subdir, 'MOLCAS.key')
                if DEBUG:
                    print('Copy:\t%s\n\t==>\n\t%s' % (fromfile, tofile))
                writefile(tofile, content)

    # link integral files
    if 'pargrad' in QMin and not QMin['template']['qmmm']:
        copyfiles = [('MOLCAS.RunFile', 'MOLCAS.RunFile')]

        # linkfiles=[('MOLCAS.OneInt','ONEINT')]
        linkfiles = []
        copyfiles.append(('MOLCAS.OneInt', 'ONEINT'))
        if QMin['template']['cholesky']:
            toappend = ['ChVec', 'QVec', 'ChRed', 'ChDiag', 'ChRst', 'ChMap']
            ls = os.listdir(os.path.join(QMin['scratchdir'], 'master'))
            for i in toappend:
                for f in ls:
                    if i in f:
                        #linkfiles.append( (f,f) )
                        copyfiles.append((f, f))
        else:
            #linkfiles.append( ('MOLCAS.OrdInt','ORDINT') )
            copyfiles.append(('MOLCAS.OrdInt', 'ORDINT'))

        if QMin['method'] in [3,4,5]:
            copyfiles.append(('GRIDFILE','GRIDFILE'))
            copyfiles.append(('MOLCAS.NqGrid','MOLCAS.NqGrid'))

        for ifile in copyfiles:
            fromfile = os.path.join(QMin['scratchdir'], 'master', ifile[0])
            tofile = os.path.join(WORKDIR, ifile[1])
            shutil.copy(fromfile, tofile)

        for ifile in linkfiles:
            fromfile = os.path.join(QMin['scratchdir'], 'master', ifile[0])
            tofile = os.path.join(WORKDIR, ifile[1])
            os.symlink(fromfile, tofile)

        if QMin['mpi_parallel']:
            for i in range(QMin['ncpu'] - 1):
                subdir = 'tmp_%i' % (i + 1)

                for ifile in copyfiles:
                    fromfile = os.path.join(QMin['scratchdir'], 'master', subdir, ifile[0])
                    tofile = os.path.join(WORKDIR, subdir, ifile[1])
                    shutil.copy(fromfile, tofile)

                for ifile in linkfiles:
                    fromfile = os.path.join(QMin['scratchdir'], 'master', subdir, ifile[0])
                    tofile = os.path.join(WORKDIR, subdir, ifile[1])
                    os.symlink(fromfile, tofile)
    return


# ======================================================================= #
def runMOLCAS(WORKDIR, MOLCAS, driver, ncpu, strip=False):
    prevdir = os.getcwd()
    os.chdir(WORKDIR)
    os.environ['WorkDir'] = WORKDIR
    os.environ['MOLCAS_NPROCS'] = str(ncpu)
    path = driver
    if not os.path.isfile(path):
        print('ERROR: could not find Molcas driver ("pymolcas" or "molcas.exe") in $MOLCAS/bin!')
        sys.exit(74)
    string = path + ' MOLCAS.input'
    stdoutfile = open(os.path.join(WORKDIR, 'MOLCAS.out'), 'w')
    stderrfile = open(os.path.join(WORKDIR, 'MOLCAS.err'), 'w')
    if PRINT or DEBUG:
        starttime = datetime.datetime.now()
        sys.stdout.write('START:\t%s\t%s\t"%s"\n' % (WORKDIR, starttime, string))
        sys.stdout.flush()
    try:
        runerror = sp.call(string, shell=True, stdout=stdoutfile, stderr=stderrfile)
        # pass
    except OSError:
        print('Call have had some serious problems:', OSError)
        sys.exit(75)
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


def doDisplacement(QMin, idir, displ):
    iatom, ixyz, isign = tuple(idir)
    QMin1 = deepcopy(QMin)
    QMin1['geo'][iatom][ixyz + 1] += isign * displ
    return QMin1


# ======================================================================= #

def parallel_speedup(N, scaling):
    # computes the parallel speedup from Amdahls law
    # with scaling being the fraction of parallelizable work and (1-scaling) being the serial part
    return 1. / ((1 - scaling) + scaling / N)

# ======================================================================= #


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
    # print(nrounds,nslots,cpu_per_run)
    return nrounds, nslots, cpu_per_run


# ======================================================================= #
def generate_joblist(QMin):
    '''split the full job into subtasks, each with a QMin dict, a WORKDIR
    structure: joblist = [ {WORKDIR: QMin, ..}, {..}, .. ]
    each element of the joblist is a set of jobs,
    and all jobs from the first set need to be completed before the second set can be processed.'''

    joblist = []
    if QMin['gradmode'] == 0 or QMin['gradmode'] == 3:
        # case of serial gradients on one cpu
        QMin1 = deepcopy(QMin)
        QMin1['master'] = []
        if 'ion' in QMin:
            QMin1['keepintegrals'] = []
        if QMin['mpi_parallel']:
            QMin1['ncpu'] = QMin['ncpu']
        else:
            QMin1['ncpu'] = 1
        QMin['nslots_pool'] = [1]
        joblist.append({'master': QMin1})


    elif QMin['gradmode'] == 1:
        # case of analytical gradients for several states on several cpus

        # we will do wavefunction and dm, soc, overlap always first
        # afterwards we will do all gradients and nacdr asynchonously
        QMin1 = deepcopy(QMin)
        QMin1['master'] = []
        QMin1['keepintegrals'] = []
        QMin1['gradmap'] = []
        QMin1['nacmap'] = []
        if 'ion' in QMin:
            QMin1['keepintegrals'] = []
        if QMin['mpi_parallel']:
            QMin1['ncpu'] = QMin['ncpu']
        else:
            QMin1['ncpu'] = 1
        QMin['nslots_pool'] = [1]
        joblist.append({'master': QMin1})

        QMin2 = deepcopy(QMin)
        remove = ['h', 'soc', 'dm', 'always_guess', 'always_orb_init', 'comment', 'ncpu', 'init', 'veloc', 'overlap', 'ion']
        for r in remove:
            QMin2 = removekey(QMin2, r)
        QMin2['gradmode'] = 0
        QMin2['pargrad'] = []
        QMin2['samestep'] = []
        ntasks = len(QMin['gradmap']) + len(QMin['nacmap'])
        if QMin['mpi_parallel']:
            # nrounds,nslots,cpu_per_run=divide_slots(QMin['ncpu'],ntasks,QMin['schedule_scaling'])
            nrounds = ntasks
            nslots = 1
            cpu_per_run = [QMin['ncpu']] * ntasks
        else:
            nrounds = 1
            nslots = QMin['ncpu']
            cpu_per_run = [1] * ntasks
        joblist.append({})
        icount = 0
        for grad in QMin['gradmap']:
            QMin3 = deepcopy(QMin2)
            QMin3['gradmap'] = [grad]
            QMin3['nacmap'] = []
            QMin3['ncpu'] = cpu_per_run[icount]
            #for i in range(len(QMin3['states'])):
            #    if not i+1 == grad[0]:
            #        QMin3['states'][i] = 0
            icount += 1
            joblist[-1]['grad_%i_%i' % grad] = QMin3
        for nac in QMin['nacmap']:
            QMin3 = deepcopy(QMin2)
            QMin3['nacmap'] = [nac]
            QMin3['gradmap'] = []
            #for i in range(len(QMin3['states'])):
            #    if not i+1 == nac[0]:
            #        QMin3['states'][i] = 0
            QMin3['overlap'] = [[j + 1, i + 1] for i in range(QMin['nmstates']) for j in range(i + 1)]
            QMin3['overlap_nacs'] = []
            QMin3['ncpu'] = cpu_per_run[icount]
            icount += 1
            joblist[-1]['nacdr_%i_%i_%i_%i' % nac] = QMin3
        QMin['nslots_pool'].append(nslots)

    elif QMin['gradmode'] == 2:
        # case of numerical gradients for ALL states, plus optionally gradients of DM and SOC
        # if only energy gradients:
        # -> do central point first, and n-1 displacements in parallel
        # -> do all other displacements afterwards
        QMin1 = deepcopy(QMin)
        QMin1['master'] = []
        if 'ion' in QMin:
            QMin1['keepintegrals'] = []
        QMin1['gradmap'] = []
        QMin1['master_displacement'] = []
        if 'h' not in QMin and 'soc' not in QMin and 'socdr' not in QMin:
            QMin1['h'] = []
        remove = ['grad', 'socdr', 'dmdr']
        for r in remove:
            QMin1 = removekey(QMin1, r)
        if 'socdr' in QMin:
            QMin1['soc'] = []
        if 'dmdr' in QMin:
            QMin1['dm'] = []
        if 'ion' in QMin:
            QMin1['keepintegrals'] = []
        if QMin['mpi_parallel']:
            QMin1['ncpu'] = QMin['ncpu']
        else:
            QMin1['ncpu'] = 1
        QMin['nslots_pool'] = [1]
        joblist.append({'master': QMin1})

        QMin2 = deepcopy(QMin)
        remove = ['comment', 'ncpu', 'veloc', 'grad', 'h', 'soc', 'dm', 'overlap', 'socdr', 'dmdr', 'ion']
        for r in remove:
            QMin2 = removekey(QMin2, r)
        QMin2['newstep'] = []
        QMin2['gradmap'] = []
        ntasks = 6 * QMin['natom']
        if QMin['mpi_parallel']:
            # nrounds,nslots,cpu_per_run=divide_slots(QMin['ncpu'],ntasks,QMin['schedule_scaling'])
            nrounds = ntasks
            nslots = 1
            cpu_per_run = [QMin['ncpu']] * ntasks
        else:
            nrounds = 1
            nslots = QMin['ncpu']
            cpu_per_run = [1] * ntasks
        QMin['nslots_pool'].append(nslots)

        #if 'socdr' in QMin or 'dmdr' in QMin:
        #    idispl=QMin['ncpu']-1
        #else:
        #    idispl=0
        icount = 0
        joblist.append({})
        for iatom in range(QMin['natom']):
            for ixyz in range(3):
                for isign in [-1., 1.]:
                    # idispl+=1
                    # if idispl==QMin['ncpu']:
                    #     joblist.append({})
                    QMin3 = deepcopy(QMin2)
                    QMin3 = doDisplacement(QMin3, [iatom, ixyz, isign], QMin['displ'])

                    # if 'socdr' in QMin or 'dmdr' in QMin or idispl>QMin['ncpu']:
                    QMin3['displacement'] = []
                    remove = ['always_guess', 'always_orb_init', 'init']
                    for r in remove:
                        QMin3 = removekey(QMin3, r)
                    if 'socdr' in QMin:
                        QMin3['soc'] = []
                    elif 'grad' in QMin:
                        QMin3['h'] = []
                    if 'dmdr' in QMin:
                        QMin3['dm'] = []
                    # if 'socdr' in QMin or 'dmdr' in QMin:
                    QMin3['overlap'] = [[j + 1, i + 1] for i in range(QMin['nmstates']) for j in range(i + 1)]
                    QMin3['ncpu'] = cpu_per_run[icount]
                    icount += 1

                    jobname = 'displ_%i_%i_%s' % (iatom, ixyz, {-1.: 'p', 1.: 'n'}[isign])
                    joblist[-1][jobname] = QMin3

    if DEBUG:
        pprint.pprint(joblist, depth=3)
    return QMin, joblist

# ======================================================================= #


def run_calc(WORKDIR, QMin):
    err = 96
    irun = -1
    while err == 96:
        irun += 1
        if 'grad' in QMin:
            if irun > 0:
                QMin['template']['gradaccudefault'] *= 10.
            if QMin['template']['gradaccudefault'] > QMin['template']['gradaccumax']:
                print('CRASHED:\t%s\tMCLR did not converge.' % (WORKDIR))
                return 96
        else:
            if irun > 0:
                print('CRASHED:\t%s\tDid not converge.' % (WORKDIR))
                return 96
        if irun >= 10:
            print('CRASHED:\t%s\tDid not converge. Giving up after 10 tries.' % (WORKDIR))
        try:
            Tasks = gettasks(QMin)
            setupWORKDIR(WORKDIR, Tasks, QMin)
            strip = 'keepintegrals' not in QMin
            err = runMOLCAS(WORKDIR, QMin['molcas'], QMin['driver'], QMin['ncpu'], strip)
        except Exception as problem:
            print('*' * 50 + '\nException in run_calc(%s)!' % (WORKDIR))
            traceback.print_exc()
            print('*' * 50 + '\n')
            raise problem
    return err

# ======================================================================= #


def runjobs(joblist, QMin):

    if 'newstep' in QMin:
        moveJobIphs(QMin)

    print('>>>>>>>>>>>>> Starting the job execution')

    errorcodes = {}
    for ijobset, jobset in enumerate(joblist):
        if not jobset:
            continue
        pool = Pool(processes=QMin['nslots_pool'][ijobset])
        for job in jobset:
            QMin1 = jobset[job]
            WORKDIR = os.path.join(QMin['scratchdir'], job)

            errorcodes[job] = pool.apply_async(run_calc, [WORKDIR, QMin1])
            # errorcodes[job]=run_calc(WORKDIR,QMin1)
            time.sleep(QMin['delay'])
        pool.close()
        pool.join()

        if 'master' in jobset:
            WORKDIR = os.path.join(QMin['scratchdir'], 'master')
            saveJobIphs(WORKDIR, jobset['master'])

        print('')

    for i in errorcodes:
        errorcodes[i] = errorcodes[i].get()

    if PRINT:
        string = '  ' + '=' * 40 + '\n'
        string += '||' + ' ' * 40 + '||\n'
        string += '||' + ' ' * 10 + 'All Tasks completed!' + ' ' * 10 + '||\n'
        string += '||' + ' ' * 40 + '||\n'
        string += '  ' + '=' * 40 + '\n'
        print(string)
        j = 0
        string = 'Error Codes:\n\n'
        for i in errorcodes:
            string += '\t%s\t%i' % (i + ' ' * (10 - len(i)), errorcodes[i])
            j += 1
            if j == 4:
                j = 0
                string += '\n'
        print(string)

    if any((i != 0 for i in errorcodes.values())):
        print('Some subprocesses did not finish successfully!')
        # sys.exit(76)

    return errorcodes

# ======================================================================= #


def collectOutputs(joblist, QMin, errorcodes):

    QMout = {}

    for jobset in joblist:
        for job in jobset:
            if errorcodes[job] == 0:
                outfile = os.path.join(QMin['scratchdir'], job, 'MOLCAS.out')
                print('Reading %s' % (outfile))
                out = readfile(outfile)
                QMout[job] = getQMout(out, jobset[job])
                if 'displacement' in jobset[job]:
                    QMout[job] = verifyQMout(QMout[job], jobset[job], out)
            else:
                if 'master' in job or 'grad' in job:
                    print('Job %s did not finish sucessfully!' % (job))
                    sys.exit(77)
                elif 'displ' in job:
                    QMout[job] = get_zeroQMout(jobset[job])

    # if DEBUG:
        # pprint.pprint(QMout,width=130)
    if DEBUG:
        for i in sorted(QMout):
            QMout1 = QMout[i]
            for j in joblist:
                if i in j:
                    QMin1 = j[i]
                    break
            print('==============================> %s <==============================' % (i))
            printQMout(QMin1, QMout1)

    return QMout

# ======================================================================= #


def phase_correction(matrix):
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


def loewdin_orthonormalization(A):
    '''
    returns loewdin orthonormalized matrix
    '''

    # S = A^T * A
    S = np.dot(A.T, A)

    # S^d = U^T * S * U
    S_diag_only, U = np.linalg.eigh(S)

    # calculate the inverse sqrt of the diagonal matrix
    S_diag_only_inverse_sqrt = [1. / (float(d) ** 0.5) for d in S_diag_only]
    S_diag_inverse_sqrt = np.diag(S_diag_only_inverse_sqrt)

    # calculate inverse sqrt of S
    S_inverse_sqrt = np.dot(np.dot(U, S_diag_inverse_sqrt), U.T)

    # calculate loewdin orthonormalized matrix
    A_lo = np.dot(A, S_inverse_sqrt)

    # normalize A_lo
    A_lo = A_lo.T
    length = len(A_lo)
    A_lon = np.zeros((length, length), dtype=complex)

    for i in range(length):
        norm_of_col = np.linalg.norm(A_lo[i])
        A_lon[i] = [e / (norm_of_col ** 0.5) for e in A_lo[i]][0]

    return A_lon.T

# ======================================================================= #


def calculate_W_dQi(H, S, e_ref):
    '''
    Return diabatized H and S
    '''

    # get diagonal of Hamiltonian
    H = np.diag([e - e_ref for e in np.diag(H)])

    # do phase correction if necessary
    if any([x for x in np.diag(S) if x < 0]):
        S = phase_correction(S)

    # do loewdin orthonorm. on overlap matrix
    U = loewdin_orthonormalization(np.matrix(S))

    return np.dot(np.dot(U.T, H), U), U



# ======================================================================= #
def overlapsign(x):
    overlapthreshold = 0.8
    if abs(x) < overlapthreshold:
        return 0.0
    else:
        return math.copysign(1, x)

# ======================================================================= #


def numdiff(enp, enn, enc, displ, o1p, o2p, o1n, o2n, iatom, idir):
    o1p = overlapsign(o1p)
    o2p = overlapsign(o2p)
    o1n = overlapsign(o1n)
    o2n = overlapsign(o2n)

    enp *= o1p * o2p
    enn *= o1n * o2n

    if (o1p == 0.0 or o2p == 0.0) and (o1n == 0.0 or o2n == 0.0):
        print('Numerical differentiation failed, both displacements have bad overlap! iatom=%i, idir=%i' % (iatom, idir))
        sys.exit(78)
    if o1p == 0.0 or o2p == 0.0:
        print('Using one-sided NumDiff for iatom=%i, idir=%i. Retaining only negative displacement.' % (iatom, idir))
        g = (enc - enn) / displ
    elif o1n == 0.0 or o2n == 0.0:
        print('Using one-sided NumDiff for iatom=%i, idir=%i. Retaining only positive displacement.' % (iatom, idir))
        g = (enp - enc) / displ
    else:
        g = (enp - enn) / 2. / displ
    return -g

# ======================================================================= #


def arrangeQMout(QMin, QMoutall, QMoutDyson):

    # sys.exit(0)

    QMout = {}
    if 'h' in QMin or 'soc' in QMin:
        QMout['h'] = QMoutall['master']['h']
    if 'dm' in QMin:
        QMout['dm'] = QMoutall['master']['dm']
    if 'overlap' in QMin:
        QMout['overlap'] = QMoutall['master']['overlap']
    # Phases from overlaps
    if 'phases' in QMin:
        if 'phases' not in QMout:
            QMout['phases'] = [complex(1., 0.) for i in range(QMin['nmstates'])]
        if 'overlap' in QMout:
            for i in range(QMin['nmstates']):
                if QMout['overlap'][i][i].real < 0.:
                    QMout['phases'][i] = complex(-1., 0.)

    if 'grad' in QMin:
        if QMin['gradmode'] == 0 or QMin['gradmode']==3:
            QMout['grad'] = QMoutall['master']['grad']

        elif QMin['gradmode'] == 1:
            zerograd = [[0.0 for xyz in range(3)] for iatom in range(QMin['natom'])]
            grad = []
            for i in sorted(QMin['statemap']):
                mult, state, ms = tuple(QMin['statemap'][i])
                if (mult, state) in QMin['gradmap']:
                    name = 'grad_%i_%i' % (mult, state)
                    grad.append(QMoutall[name]['grad'][i - 1])
                else:
                    grad.append(zerograd)
            QMout['grad'] = grad

        elif QMin['gradmode'] == 2:
            grad = [[[0.0 for xyz in range(3)] for iatom in range(QMin['natom'])] for istate in range(QMin['nmstates'])]
            for iatom in range(QMin['natom']):
                for xyz in range(3):
                    namep = 'displ_%i_%i_p' % (iatom, xyz)
                    namen = 'displ_%i_%i_n' % (iatom, xyz)
                    displ = QMin['displ']

                    # diabatization
                    if QMin['template']['diab_num_grad']:
                        Hmaster = QMoutall['master']['h']
                        Hpos = deepcopy(QMoutall[namep]['h'])
                        Spos = deepcopy(QMoutall[namep]['overlap'])
                        Hneg = deepcopy(QMoutall[namen]['h'])
                        Sneg = deepcopy(QMoutall[namen]['overlap'])
                        Hpos, Spos = calculate_W_dQi(Hpos, Spos, QMout['h'][0][0])
                        Hneg, Sneg = calculate_W_dQi(Hneg, Sneg, QMout['h'][0][0])
                        QMoutall[namep]['h_diab'] = Hpos
                        QMoutall[namep]['s_diab'] = Spos
                        QMoutall[namen]['h_diab'] = Hneg
                        QMoutall[namen]['s_diab'] = Sneg
                    else:
                        Hmaster = QMoutall['master']['h']
                        Hpos = QMoutall[namep]['h']
                        Hneg = QMoutall[namen]['h']
                        Spos = QMoutall[namep]['overlap']
                        Sneg = QMoutall[namen]['overlap']


                    for istate in range(QMin['nmstates']):

                        enc = Hmaster[istate][istate].real

                        enp = Hpos[istate][istate].real
                        ovp = Spos[istate][istate].real

                        enn = Hneg[istate][istate].real
                        ovn = Sneg[istate][istate].real

                        g = numdiff(enp, enn, enc, displ, ovp, ovp, ovn, ovn, iatom, xyz)
                        grad[istate][iatom][xyz] = g
            QMout['grad'] = grad

    if 'grad' in QMin and QMin['template']['cobramm']:
        if QMin['gradmode'] == 1:
            with open("efield.dat", "w") as f:
                for i in sorted(QMin['statemap']):
                    mult, state, ms = tuple(QMin['statemap'][i])
                    if (mult, state) in QMin['gradmap']:
                        name = 'grad_%i_%i' % (mult, state)
                        efield = QMoutall[name]['efield'][i - 1]
                        f.write("{0} 3 ! {1} {2}\n".format(len(efield), mult, state))
                        for e in efield:
                            f.write("{0:20.10f} {1:20.10f} {2:20.10f}\n".format(*e))
                    else:
                        f.write("{0} 3 ! {1} {2}\n".format(0, mult, state))

    if 'nacdr' in QMin:
        QMout['nacdr'] = [[[[0. for i in range(3)] for j in range(QMin['natom'])] for k in range(QMin['nmstates'])] for l in range(QMin['nmstates'])]
        for i in sorted(QMin['statemap']):
            for j in sorted(QMin['statemap']):
                m1, s1, ms1 = tuple(QMin['statemap'][i])
                m2, s2, ms2 = tuple(QMin['statemap'][j])
                if not m1 == m2:
                    continue
                if not ms1 == ms2:
                    continue
                if s1 == s2:
                    continue
                if (m1, s1, m2, s2) in QMin['nacmap']:
                    name = 'nacdr_%i_%i_%i_%i' % (m1, min(s1, s2), m2, max(s1, s2))
                    QMout['nacdr'][i - 1][j - 1] = deepcopy(QMoutall[name]['nacdr'][i - 1][j - 1])
                    QMout['nacdr'][j - 1][i - 1] = deepcopy(QMoutall[name]['nacdr'][j - 1][i - 1])
                    ovl = QMoutall[name]['overlap'][i - 1][i - 1].real
                    ovl *= QMoutall[name]['overlap'][j - 1][j - 1].real
                    for iatom in range(QMin['natom']):
                        for xyz in range(3):
                            QMout['nacdr'][i - 1][j - 1][iatom][xyz] *= ovl
                            QMout['nacdr'][j - 1][i - 1][iatom][xyz] *= ovl

    if 'socdr' in QMin:
        socdr = [[[[0.0 for xyz in range(3)] for iatom in range(QMin['natom'])] for istate in range(QMin['nmstates'])] for jstate in range(QMin['nmstates'])]
        displ = QMin['displ']
        for iatom in range(QMin['natom']):
            for xyz in range(3):
                namep = 'displ_%i_%i_p' % (iatom, xyz)
                namen = 'displ_%i_%i_n' % (iatom, xyz)
                for istate in range(QMin['nmstates']):
                    for jstate in range(QMin['nmstates']):
                        if istate == jstate:
                            continue

                        # diabatization
                        if QMin['template']['diab_num_grad']:
                            Hmaster = QMoutall['master']['h']

                            Hpos = deepcopy(QMoutall[namep]['h'])
                            Hpos = np.array(Hpos) - np.diag([e - QMout['h'][0][0] for e in np.diag(Hpos)])
                            Spos = QMoutall[namep]['s_diab']
                            Hpos = np.dot(np.dot(Spos.T, Hpos), Spos)

                            Hneg = deepcopy(QMoutall[namen]['h'])
                            Hneg = np.array(Hneg) - np.diag([e - QMout['h'][0][0] for e in np.diag(Hneg)])
                            Sneg = QMoutall[namen]['s_diab']
                            Hneg = np.dot(np.dot(Sneg.T, Hneg), Sneg)
                        else:
                            Hmaster = QMoutall['master']['h']
                            Hpos = QMoutall[namep]['h']
                            Hneg = QMoutall[namen]['h']
                            Spos = QMoutall[namep]['overlap']
                            Sneg = QMoutall[namen]['overlap']

                        enc = Hmaster[istate][jstate]

                        enp = Hpos[istate][jstate]
                        o1p = Spos[istate][istate].real
                        o2p = Spos[jstate][jstate].real

                        enn = Hneg[istate][jstate]
                        o1n = Sneg[istate][istate].real
                        o2n = Sneg[jstate][jstate].real

                        g = numdiff(enp, enn, enc, displ, o1p, o2p, o1n, o2n, iatom, xyz)
                        socdr[istate][jstate][iatom][xyz] = g
        QMout['socdr'] = socdr

    if 'dmdr' in QMin:
        dmdr = [[[[[0.0 for xyz in range(3)] for iatom in range(QMin['natom'])] for istate in range(QMin['nmstates'])] for jstate in range(QMin['nmstates'])] for ipol in range(3)]
        displ = QMin['displ']
        for iatom in range(QMin['natom']):
            for xyz in range(3):
                namep = 'displ_%i_%i_p' % (iatom, xyz)
                namen = 'displ_%i_%i_n' % (iatom, xyz)
                for ipol in range(3):
                    for istate in range(QMin['nmstates']):
                        for jstate in range(QMin['nmstates']):

                            # diabatization
                            if QMin['template']['diab_num_grad']:
                                Hmaster = QMoutall['master']['dm'][ipol]

                                Hpos = deepcopy(QMoutall[namep]['dm'][ipol])
                                Hpos = np.array(Hpos) - np.diag([e - QMout['h'][0][0] for e in np.diag(Hpos)])
                                Spos = QMoutall[namep]['s_diab']
                                Hpos = np.dot(np.dot(Spos.T, Hpos), Spos)

                                Hneg = deepcopy(QMoutall[namen]['dm'][ipol])
                                Hneg = np.array(Hneg) - np.diag([e - QMout['h'][0][0] for e in np.diag(Hneg)])
                                Sneg = QMoutall[namen]['s_diab']
                                Hneg = np.dot(np.dot(Sneg.T, Hneg), Sneg)
                            else:
                                Hmaster = QMoutall['master']['dm'][ipol]
                                Hpos = QMoutall[namep]['dm'][ipol]
                                Hneg = QMoutall[namen]['dm'][ipol]
                                Spos = QMoutall[namep]['overlap']
                                Sneg = QMoutall[namen]['overlap']

                            enc = Hmaster[istate][jstate].real

                            enp = Hpos[istate][jstate].real
                            o1p = Spos[istate][istate].real
                            o2p = Spos[jstate][jstate].real

                            enn = Hneg[istate][jstate].real
                            o1n = Sneg[istate][istate].real
                            o2n = Sneg[jstate][jstate].real

                            g = numdiff(enp, enn, enc, displ, o1p, o2p, o1n, o2n, iatom, xyz)
                            dmdr[ipol][istate][jstate][iatom][xyz] = g
        QMout['dmdr'] = dmdr

    if 'ion' in QMin:
        QMout['prop'] = QMoutDyson

    if PRINT:
        print('\n===================================================================')
        print('========================= Final Results ===========================')
        print('===================================================================')
        printQMout(QMin, QMout)

    return QMout

# ======================================================================= #


def getcaspt2weight(out, mult, state):
    modulestring = '&CASPT2'
    spinstring = 'Spin quantum number'
    statestring = 'Compute H0 matrices for state'
    refstring = 'Reference weight'
    stateindex = 5
    refindex = 2

    module = False
    correct_mult = False
    correct_state = False
    for i, line in enumerate(out):
        if modulestring in line:
            module = True
        elif 'Stop Module' in line:
            module = False
            correct_mult = False
            correct_state = False
        elif spinstring in line and module:
            spin = float(line.split()[3])
            if int(2 * spin) + 1 == mult:
                correct_mult = True
        elif statestring in line and module and correct_mult:
            if state == int(line.split()[stateindex]):
                correct_state = True
        elif refstring in line and module and correct_mult and correct_state:
            return float(line.split()[refindex])
    print('CASPT2 reference weight of state %i in mult %i not found!' % (state, mult))
    sys.exit(80)

# ======================================================================= #


def getcaspt2transform(out, mult):
    modulestring = '&CASPT2'
    spinstring = 'Spin quantum number'
    statestring = 'Number of CI roots used'
    matrixstring = 'Eigenvectors:'
    singlestatestring = 'This is a CASSCF or RASSCF reference function'
    stateindex = 5
    refindex = 2

    module = False
    correct_mult = False
    nstates = 0

    for i, line in enumerate(out):
        if modulestring in line:
            module = True
        elif 'Stop Module' in line:
            module = False
            correct_mult = False
        elif spinstring in line and module:
            spin = float(line.split()[3])
            if int(2 * spin) + 1 == mult:
                correct_mult = True
        elif statestring in line and module:
            nstates = int(line.split()[stateindex])
        elif singlestatestring in line and module and correct_mult:
            return [[1.]]
        elif matrixstring in line and module and correct_mult:
            t = [[0. for x in range(nstates)] for y in range(nstates)]
            for x in range(nstates):
                for y in range(nstates):
                    lineshift = i + y + 1 + x // 5 * (nstates + 1)
                    indexshift = x % 5
                    t[x][y] = float(out[lineshift].split()[indexshift])
            return t
    print('MS-CASPT2 transformation matrix in mult %i not found!' % (mult))
    sys.exit(81)

# ======================================================================= #


def verifyQMout(QMout, QMin, out):
    # checks whether a displacement calculation gave a sensible result
    # currently only checks for CASPT2 problems (reference weight)

    refweight_ratio = 0.80

    if QMin['method'] in [0, 3]:
        # CASSCF case
        pass
    elif QMin['method'] == 5 and QMin['pdft-functional'] != -1:
        pass
    elif QMin['method'] in [1, 2, 4, 5, 6]:
        # SS-CASPT2 and MS-CASPT2 cases
        refs = []
        for istate in range(QMin['nmstates']):
            mult, state, ms = tuple(QMin['statemap'][istate + 1])
            if QMin['method'] in [1, 2]:
                refs.append(getcaspt2weight(out, mult, state))
            elif QMin['method'] in [3,4,5]:
                refs.append(1.)
            # print mult,state,refs[-1]

        # MS-CASPT2: get eigenvectors and transform
        if QMin['method'] in [2,6]:
            offset = 0
            for imult, nstate in enumerate(QMin['states']):
                if nstate == 0:
                    continue
                mult = imult + 1
                for ims in range(mult):
                    t = getcaspt2transform(out, mult)
                    # pprint.pprint(t)
                    refslice = refs[offset:offset + nstate]
                    # print refslice
                    newref = [0. for i in range(nstate)]
                    for i in range(nstate):
                        for j in range(nstate):
                            newref[i] += refslice[j] * t[i][j]**2
                        refs[offset + i] = newref[i]
                    # print newref
                    offset += nstate
        for istate in range(QMin['nmstates']):
            mult, state, ms = tuple(QMin['statemap'][istate + 1])
            # print mult,state,refs[istate]

        # check the reference weights and set overlap to zero if not acceptable
            for istate in range(QMin['nmstates']):
                if refs[istate] < max(refs) * refweight_ratio:
                    QMout['overlap'][istate][istate] = complex(0., 0.)
                    # print('Set to zero:',istate)

    return QMout

# ======================================================================= #


def get_zeroQMout(QMin):
    nmstates = QMin['nmstates']
    natom = QMin['natom']
    QMout = {}
    if 'h' in QMin or 'soc' in QMin:
        QMout['h'] = [[complex(0.0) for i in range(nmstates)] for j in range(nmstates)]
    if 'dm' in QMin:
        QMout['dm'] = [[[complex(0.0) for i in range(nmstates)] for j in range(nmstates)] for xyz in range(3)]
    if 'overlap' in QMin:
        QMout['overlap'] = [[complex(0.0) for i in range(nmstates)] for j in range(nmstates)]
    if 'grad' in QMin:
        QMout['grad'] = [[[0., 0., 0.] for i in range(natom)] for j in range(nmstates)]
    return QMout







# ======================================================================= #
def cleanupSCRATCH(SCRATCHDIR):
    ''''''
    if PRINT:
        print('===> Removing directory %s\n' % (SCRATCHDIR))
    try:
        if True:
            shutil.rmtree(SCRATCHDIR)
        else:
            print('not removing anything. SCRATCHDIR is %s' % SCRATCHDIR)
    except OSError:
        print('Could not remove directory %s' % (SCRATCHDIR))

# ======================================================================= #


def saveJobIphs(WORKDIR, QMin):
    # Moves JobIph files from WORKDIR to savedir
    for imult, nstates in enumerate(QMin['states']):
        if nstates < 1:
            continue
        fromfile = os.path.join(WORKDIR, 'MOLCAS.%i.JobIph' % (imult + 1))
        if not os.path.isfile(fromfile):
            print('File %s not found, cannot move to savedir!' % (fromfile))
            sys.exit(82)
        tofile = os.path.join(QMin['savedir'], 'MOLCAS.%i.JobIph' % (imult + 1))
        if DEBUG:
            print('Copy:\t%s\n\t==>\n\t%s' % (fromfile, tofile))
        shutil.copy(fromfile, tofile)

        if 'molden' in QMin:
            # copy MOLDEN files
            path = os.path.join(QMin['savedir'], 'MOLDEN')
            if not os.path.isdir(path):
                try:
                    os.makedirs(path)
                except OSError:
                    pass
            try:
                fromfile = os.path.join(WORKDIR, 'MOLCAS.%i.molden' % (imult + 1))
                tofile = os.path.join(QMin['savedir'], 'MOLDEN', 'MOLCAS.%i.molden' % (imult + 1))
                shutil.copy(fromfile, tofile)
            except OSError:
                pass

        if 'master_displacement' in QMin:
            fromfile = os.path.join(QMin['savedir'], 'MOLCAS.%i.JobIph' % (imult + 1))
            tofile = os.path.join(QMin['savedir'], 'MOLCAS.%i.JobIph.master' % (imult + 1))
            if DEBUG:
                print('Copy:\t%s\n\t==>\n\t%s' % (fromfile, tofile))
            shutil.copy(fromfile, tofile)

# ======================================================================= #


def moveJobIphs(QMin):
    # moves all relevant JobIph files in the savedir to old-JobIph files
    # deletes also all old .master files
    for imult, nstates in enumerate(QMin['states']):
        if nstates < 1:
            continue
        fromfile = os.path.join(QMin['savedir'], 'MOLCAS.%i.JobIph' % (imult + 1))
        if not os.path.isfile(fromfile):
            print('File %s not found, cannot move to OLD!' % (fromfile))
            sys.exit(83)
        tofile = os.path.join(QMin['savedir'], 'MOLCAS.%i.JobIph.old' % (imult + 1))
        if DEBUG:
            print('Copy:\t%s\n\t==>\n\t%s' % (fromfile, tofile))
        shutil.copy(fromfile, tofile)

    ls = os.listdir(QMin['savedir'])
    for i in ls:
        if '.master' in i:
            rmfile = os.path.join(QMin['savedir'], i)
            os.remove(rmfile)

# ======================================================================= #


def stripWORKDIR(WORKDIR):
    ls = os.listdir(WORKDIR)
    keep = ['MOLCAS.out', 'MOLCAS\\.[1-9]\\.JobIph', 'MOLCAS\\.[1-9]\\.RasOrb', 'MOLCAS\\.[1-9]\\.molden']
    for ifile in ls:
        delete = True
        for k in keep:
            if containsstring(k, ifile):
                delete = False
        if delete:
            rmfile = os.path.join(WORKDIR, ifile)
            if not DEBUG:
                if os.path.isdir(rmfile):
                    cleanupSCRATCH(rmfile)
                else:
                    os.remove(rmfile)



# =============================================================================================== #
# =============================================================================================== #
# ===========================================  Dyson norms  ===================================== #
# =============================================================================================== #
# =============================================================================================== #

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
            if step[k] == 1:
                m2 += powmin1(det[k] + 1)
                num = bval[k] + powmin1(det[k] + 1) * m2
                denom = 2. * bval[k]
                if num == 0.:
                    break
                coeff *= 1. * num / denom
            elif step[k] == 2:
                m2 += powmin1(det[k] - 1)
                num = bval[k] + 2 + powmin1(det[k]) * m2
                denom = 2. * (bval[k] + 2)
                sign *= powmin1(bval[k] + 2 - det[k])
                if num == 0.:
                    break
                coeff *= 1. * num / denom
            elif step[k] == 3:
                sign *= powmin1(bval[k])
                num = 1.

        # add determinant to dict if coefficient non-zero
        if num != 0.:
            dets[tuple(det)] = 1. * sign * math.sqrt(coeff)

    #pprint.pprint( dets)
    return dets

# ======================================================================= #


def get_determinants(out, mult):

    # first, find the correct RASSI output section for the given multiplicity
    modulestring = '&RASSI'
    spinstring = 'SPIN MULTIPLICITY:'
    stopstring = 'The following data are common to all the states'
    module = False
    jobiphmult = []
    for iline, line in enumerate(out):
        if modulestring in line:
            module = True
            jobiphmult = []
        elif module:
            if spinstring in line:
                jobiphmult.append(int(line.split()[-1]))
            if stopstring in line:
                if all(i == mult for i in jobiphmult):
                    break
                else:
                    module = False
    else:
        print('Determinants not found!', mult)
        print('No RASSI run for multiplicity %i found!' % (mult))
        sys.exit(84)

    # ndocc and nvirt
    ndocc = -1
    nvirt = -1
    while True:
        iline += 1
        line = out[iline]
        if ' INACTIVE ' in line:
            ndocc = int(line.split()[-1])
        if ' SECONDARY ' in line:
            nvirt = int(line.split()[-1])
        if ndocc != -1 and nvirt != -1:
            break

    # Get number of states
    while True:
        iline += 1
        line = out[iline]
        if 'Nr of states:' in line:
            nstates = int(line.split()[-1])
            break

    # Now start searching at iline, collecting all CI vectors
    # the dict has det-string tuples as keys and lists of float as values
    ci_vectors = {}
    statesstring = 'READCI called for state'
    stopstring = '****************************'
    finalstring = 'HAMILTONIAN MATRIX'
    done = set()

    finished = False
    while True:
        while True:
            iline += 1
            line = out[iline]
            if statesstring in line:
                state = int(line.split()[-1])
                break
            if finalstring in line:
                finished = True
                break
        if finished:
            break
        if state in done:
            continue
        done.add(state)
        iline+=13
        # TODO: verify now many lines need to be skipped! in newer versions might only be 10 or 11!
        while True:
            iline += 1
            line = out[iline]
            if stopstring in line:
                break
            s = line.split()
            if s == []:
                continue
            coef = float(s[-2])
            if coef == 0.:
                continue
            csf = s[-3]
            step = []
            for i in csf:
                if i == '2':
                    step.append(3)
                elif i == 'd':
                    step.append(2)
                elif i == 'u':
                    step.append(1)
                elif i == '0':
                    step.append(0)
            dets = decompose_csf((mult - 1) / 2., step)
            # add dets to ci_vectors
            for det in dets:
                if det in ci_vectors:
                    d = state - len(ci_vectors[det])
                    if d > 0:
                        ci_vectors[det].extend([0.] * d)
                    ci_vectors[det][state - 1] += coef * dets[det]
                else:
                    ci_vectors[det] = [0.] * state
                    ci_vectors[det][state - 1] += coef * dets[det]
    for det in ci_vectors.keys():
        d = nstates - len(ci_vectors[det])
        if d > 0:
            ci_vectors[det].extend([0.] * d)
        # if all( i==0. for i in ci_vectors[det] ):
            #del ci_vectors[det]
    ci_vectors['ndocc'] = ndocc
    ci_vectors['nvirt'] = nvirt
    #pprint.pprint( ci_vectors)
    return ci_vectors

# ======================================================================= #


def format_ci_vectors(ci_vectors):
    # get nstates, norb and ndets
    for key in ci_vectors:
        if key != 'ndocc' and key != 'nvirt':
            nstates = len(ci_vectors[key])
            norb = len(key)
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

# ======================================================================= #


def runWFOVERLAPS(WORKDIR, wfoverlaps, memory=100, ncpu=1):
    prevdir = os.getcwd()
    os.chdir(WORKDIR)
    string = wfoverlaps + ' -m %i' % (memory) + ' -f dyson.in'
    stdoutfile = open(os.path.join(WORKDIR, 'dyson.out'), 'w')
    stderrfile = open(os.path.join(WORKDIR, 'dyson.err'), 'w')
    os.environ['OMP_NUM_THREADS'] = str(ncpu)
    if PRINT or DEBUG:
        starttime = datetime.datetime.now()
        sys.stdout.write('START:\t%s\t%s\t"%s"\n' % (WORKDIR, starttime, string))
        sys.stdout.flush()
    try:
        runerror = sp.call(string, shell=True, stdout=stdoutfile, stderr=stderrfile)
    except OSError:
        print('Call have had some serious problems:', OSError)
        sys.exit(85)
    stdoutfile.close()
    stderrfile.close()
    if PRINT or DEBUG:
        endtime = datetime.datetime.now()
        sys.stdout.write('FINISH:\t%s\t%s\tRuntime: %s\tError Code: %i\n' % (WORKDIR, endtime, endtime - starttime, runerror))
        sys.stdout.flush()
    os.chdir(prevdir)
    return runerror

# ======================================================================= #


def get_dysonel(out, s1, s2):
    ilines = -1
    while True:
        ilines += 1
        if ilines == len(out):
            print('Overlap of states %i - %i not found!' % (s1, s2))
            sys.exit(86)
        if containsstring('Dyson norm matrix |<PsiA_i|PsiB_j>|^2', out[ilines]):
            break
    ilines += 1 + s1
    f = out[ilines].split()
    return float(f[s2 + 1])

# ======================================================================= #


def do_Dyson(QMin):
    print('\n>>>>>>>>>>>>> Starting the Dyson calculations')

    # create list of mult-pairings
    states = QMin['states']
    mult_pairs = []
    for i in range(len(states) - 1):
        if states[i] != 0 and states[i + 1] != 0:
            mult_pairs.append((i + 1, i + 2))

    # create directories for each pair
    for pair in mult_pairs:
        path = os.path.join(QMin['scratchdir'], 'Dyson_%i_%i' % (pair[0], pair[1]))
        if os.path.exists(path):
            if os.path.isfile(path):
                print('%s exists and is a file!\n' % (path))
                sys.exit(87)
        else:
            try:
                os.makedirs(path)
            except OSError:
                print('Can not create %s\n' % (QMin['scratchdir']))
                sys.exit(88)

    # create list of multiplicities to treat and mapping to directories
    mult_map = {}
    for pair in mult_pairs:
        if pair[0] in mult_map:
            mult_map[pair[0]].append(pair)
        else:
            mult_map[pair[0]] = [pair]
        if pair[1] in mult_map:
            mult_map[pair[1]].append(pair)
        else:
            mult_map[pair[1]] = [pair]
    # print mult_map

    # obtain the determinant files from master/MOLCAS.out
    path = os.path.join(QMin['scratchdir'], 'master', 'MOLCAS.out')
    out = readfile(path)
    nfrozen = {}
    for mult in mult_map:
        ci_vectors = get_determinants(out, mult)
        nfrozen[mult] = ci_vectors['ndocc']
        for pair in mult_map[mult]:
            path = os.path.join(QMin['scratchdir'], 'Dyson_%i_%i' % (pair[0], pair[1]), 'dets.%i' % mult)
            # print path
            string = format_ci_vectors(ci_vectors)
            writefile(path, string)

    # link the RasOrb files from master to Dyson_*_*
    for mult in mult_map:
        fromfile = os.path.join(QMin['scratchdir'], 'master', 'MOLCAS.%i.RasOrb' % mult)
        for pair in mult_map[mult]:
            tofile = os.path.join(QMin['scratchdir'], 'Dyson_%i_%i' % (pair[0], pair[1]), 'MOLCAS.%i.RasOrb' % mult)
            link(fromfile, tofile)

    # link the integral files
    fromRunFile = os.path.join(QMin['scratchdir'], 'master', 'MOLCAS.RunFile')
    fromOneInt = os.path.join(QMin['scratchdir'], 'master', 'MOLCAS.OneInt')
    for pair in mult_pairs:
        path = os.path.join(QMin['scratchdir'], 'Dyson_%i_%i' % (pair[0], pair[1]))
        toRunFile = os.path.join(path, 'RUNFILE')
        link(fromRunFile, toRunFile)
        toOneInt = os.path.join(path, 'ONEINT')
        link(fromOneInt, toOneInt)

    # write input files
    for pair in mult_pairs:
        path = os.path.join(QMin['scratchdir'], 'Dyson_%i_%i' % (pair[0], pair[1]))
        inputfile = os.path.join(path, 'dyson.in')
        string = '''ao_read=1
a_mo=MOLCAS.%i.RasOrb
b_mo=MOLCAS.%i.RasOrb
a_det=dets.%i
b_det=dets.%i
a_mo_read=1
b_mo_read=1
same_aos
''' % (pair[0], pair[1], pair[0], pair[1])
        if 'ncore' in QMin:
            string += 'ncore=%i\n' % QMin['ncore']
        if 'ndocc' in QMin:
            string += 'ndocc=%i\n' % QMin['ndocc']
        #frozen=min( nfrozen[pair[0]],nfrozen[pair[1]] )
        #string+='\nndocc=%i\n' % (frozen)
        writefile(inputfile, string)

    # run the jobs subsequently with full number of CPUs
    errorcodes = {}
    for pair in mult_pairs:
        path = os.path.join(QMin['scratchdir'], 'Dyson_%i_%i' % (pair[0], pair[1]))
        runerror = runWFOVERLAPS(path, QMin['wfoverlap'], QMin['memory'], QMin['ncpu'])
        errorcodes['Dyson_%i_%i' % (pair[0], pair[1])] = runerror
    if PRINT:
        string = '\n  ' + '=' * 40 + '\n'
        string += '||' + ' ' * 40 + '||\n'
        string += '||' + ' ' * 10 + 'All Tasks completed!' + ' ' * 10 + '||\n'
        string += '||' + ' ' * 40 + '||\n'
        string += '  ' + '=' * 40 + '\n'
        print(string)
        j = 0
        string = 'Error Codes:\n\n'
        for i in errorcodes:
            string += '\t%s\t%i\n' % (i + ' ' * (10 - len(i)), errorcodes[i])
        print(string)
    if any((i != 0 for i in errorcodes.values())):
        print('Some subprocesses did not finish successfully!')
        sys.exit(89)

    # get the Dyson norms
    nmstates = QMin['nmstates']
    QMoutDyson = makecmatrix(nmstates, nmstates)
    for pair in mult_pairs:
        path = os.path.join(QMin['scratchdir'], 'Dyson_%i_%i' % (pair[0], pair[1]), 'dyson.out')
        out = readfile(path)
        imult = pair[0]
        jmult = pair[1]
        for i in range(nmstates):
            for j in range(nmstates):
                m1, s1, ms1 = tuple(QMin['statemap'][i + 1])
                m2, s2, ms2 = tuple(QMin['statemap'][j + 1])
                if not (imult, jmult) == (m1, m2) and not (imult, jmult) == (m2, m1):
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
                QMoutDyson[i][j] = get_dysonel(out, s1, s2) * factor

    return QMoutDyson




















# =============================================================================================== #
# =============================================================================================== #
# =========================================== Main routine  ===================================== #
# =============================================================================================== #
# =============================================================================================== #

# ========================== Main Code =============================== #
def main():

    # Retrieve PRINT and DEBUG
    try:
        envPRINT = os.getenv('SH2CAS_PRINT')
        if envPRINT and envPRINT.lower() == 'false':
            global PRINT
            PRINT = False
        envDEBUG = os.getenv('SH2CAS_DEBUG')
        if envDEBUG and envDEBUG.lower() == 'true':
            global DEBUG
            DEBUG = True
    except ValueError:
        print('PRINT or DEBUG environment variables do not evaluate to numerical values!')
        sys.exit(90)

    # Process Command line arguments
    if len(sys.argv) != 2:
        print('Usage:\n./SHARC_MOLCAS.py <QMin>\n')
        print('version:', version)
        print('date:', versiondate)
        print('changelog:\n', changelogstring)
        sys.exit(91)
    QMinfilename = sys.argv[1]

    # Print header
    printheader()

    # Read QMinfile
    QMin = readQMin(QMinfilename)

    # make list of jobs
    QMin, joblist = generate_joblist(QMin)

    # run all MOLCAS jobs
    errorcodes = runjobs(joblist, QMin)

    # get output
    QMoutall = collectOutputs(joblist, QMin, errorcodes)

    # extract data, perform Dyson calculations
    if 'ion' in QMin:
        QMoutDyson = do_Dyson(QMin)
    else:
        QMoutDyson = None

    # format final output
    QMout = arrangeQMout(QMin, QMoutall, QMoutDyson)

    # Measure time
    runtime = measuretime()
    QMout['runtime'] = runtime

    # Write QMout
    writeQMout(QMin, QMout, QMinfilename)

    # Remove Scratchfiles from SCRATCHDIR
    if not DEBUG:
        cleanupSCRATCH(QMin['scratchdir'])
        if 'cleanup' in QMin:
            cleanupSCRATCH(QMin['savedir'])
    if PRINT or DEBUG:
        print('#================ END ================#')


if __name__ == '__main__':
    main()
