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

# TODO: transfer plotting function from make_fitscript.py (prettier plots)

# Modules:
# Operating system, isfile and related routines, move files, create directories
import sys
import os
import shutil
# External Calls
import subprocess as sp
# Regular expressions
import re
# debug print(for dicts and arrays)
import pprint
# sqrt and other math
import math
# copy of arrays of arrays
from copy import deepcopy
# others
import time
import datetime
from optparse import OptionParser
import readline
import colorsys
import random
from socket import gethostname

import numpy as np
import scipy.integrate as spint
import scipy.optimize as spopt

try:
    import numpy
    NUMPY = True
except ImportError:
    NUMPY = False

# =========================================================0

version = '3.0'
versiondate = datetime.date(2023, 4, 1)

changelogstring = '''

'''

# ======================================================================= #
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



# ======================================================================================================================


def displaywelcome():
    string = '\n'
    string += '  ' + '=' * 80 + '\n'
    string += '||' + '{:^80}'.format('') + '||\n'
    string += '||' + '{:^80}'.format('Direct fitting for SHARC populations') + '||\n'
    string += '||' + '{:^80}'.format('') + '||\n'
    string += '||' + '{:^80}'.format('Author: Sebastian Mai') + '||\n'
    string += '||' + '{:^80}'.format('') + '||\n'
    string += '||' + '{:^80}'.format('Version:' + version) + '||\n'
    string += '||' + '{:^80}'.format(versiondate.strftime("%d.%m.%y")) + '||\n'
    string += '||' + '{:^80}'.format('') + '||\n'
    string += '  ' + '=' * 80 + '\n\n'
    string += '''
This script fits SHARC populations (as generated with populations.py)
to general kinetic models based on first-order population transfer.

'''
    print(string)

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================


class rgbcolor:
    '''an object which you initialize with a list of integers
  and whose hexcolor() routine returns a hex-coded color for a given pair (index,state)
  initialize: [6,0,3]
  - each non-empty group is allocated the same space on the colorwheel
  - each group space is divided equally between the elements of this group

  => first group gets 180deg of the colorwheel, each element gets 30deg
  => second group is empty, does not get space
  => third group gets 180deg of the colorwheel, each element gets 60deg

  - the script also allows to eliminate certain colors
  - if the Index,Element pair is invalid (e.g. (2,1) for the above input), it returns white (#FFFFFF)

  Example of usage:
  a=[6,0,3]
  R=rgbcolor(a)
  for index,num in enumerate(a):
    for el in range(num):
      print(index,el,R.hexcolor(index+1,el+1))
  print 2,1,R.hexcolor(2,1)

  Output:
  1     1         #FF0000
  1     2         #FF7F00
  1     3         #FFFF00
  1     4         #7FFF00
  1     5         #00FF00
  1     6         #00FF7F
  3     1         #00FFFF
  3     2         #0000FF
  3     2         #FF00FF
  2     1         #FFFFFF #invalid, hence white
  '''

    def __init__(self, initlist):
        # Hues:
        # 0.0     0.15       ...
        # Red     Yellow     ...
        excluded = [
            [0.12, 0.22]       # exclude yellow hues from the colorwheel
        ]
        # excluded-list must be sorted, for each element x[0]<=x[1]
        # each excluded[i][0]<=excluded[i+1][0]
        # everything between 0 and 1
        # sort the pairs
        for i, el in enumerate(excluded):
            excluded[i] = [min(el), max(el)]
        # sort the list
        excluded.sort(key=lambda x: x[0])
        # make all elements between 0 and 1
        temp1 = []
        for i, el in enumerate(excluded):
            temp1.append([min(1., max(0., el[0])), max(0., min(1., el[1]))])
        # resolve overlapping ranges
        temp2 = [[0., 0.]]
        for i, el in enumerate(temp1):
            if el[0] >= temp2[-1][1]:
                temp2.append(el)
            else:
                temp2[-1][1] = el[1]
        self.excluded = temp2
        # number of non-empty groups
        self.initlist = initlist
        n = 0
        for index, el in enumerate(initlist):
            # negative numbers in initlist are not allowed, make them to zero
            self.initlist[index] = max(0, el)
            if el > 0:
                n += 1
        self.n = n                    # number of non-empty groups
        self.m = len(initlist)        # number of groups
        # available colorspace self.a and shifts to skip excluded regions self.ex
        self.a = 1.
        self.ex = [0. for i in self.excluded]
        for i, el in enumerate(self.excluded):
            self.a -= el[1] - el[0]
            self.ex[i] = el[1] - el[0]
        # list of starting values
        self.startlist = [0. for i in range(self.m)]
        for index, el in enumerate(initlist):
            if index == 0:
                continue
            if el > 0:
                self.startlist[index] = self.a / self.n + self.startlist[index - 1]
            else:
                self.startlist[index] = self.startlist[index - 1]
        # list of increments
        self.incrlist = [0. for i in range(self.m)]
        for index, el in enumerate(initlist):
            if el > 0:
                self.incrlist[index] = self.a / self.n / el

    def rgb_to_hex(self, rgb):
        # rgb is a list with three elements, which are floats 0<=x<=1
        color = 0
        for i in range(3):
            a = max(min(rgb[i], 1.0), 0.0)
            color += int(255 * a) * 256**(2 - i)
        # color=int(255*triple[0])*256**2+int(255*triple[1])*256+int(255*triple[2])
        string = hex(color)[2:].upper()
        string = '#' + '0' * (6 - len(string)) + string
        return string

    def hexcolor(self, index, el):
        if not 1 <= index <= self.m:
            return '#FFFFFF'
        if not 1 <= el <= self.initlist[index - 1]:
            return '#FFFFFF'
        deg = self.startlist[index - 1] + self.incrlist[index - 1] * (el - 1)
        for i, el in enumerate(self.excluded):
            if deg > el[0]:
                deg += self.ex[i]
        rgbtriple = colorsys.hsv_to_rgb(deg, 1, 1)
        return self.rgb_to_hex(rgbtriple)

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================


def open_keystrokes():
    global KEYSTROKES
    KEYSTROKES = open('KEYSTROKES.tmp', 'w')


def close_keystrokes():
    KEYSTROKES.close()
    shutil.move('KEYSTROKES.tmp', 'KEYSTROKES.make_fit')

# ===================================


def question(question, typefunc, default=None, autocomplete=True, ranges=False):
    if typefunc == int or typefunc == float:
        if default is not None and not isinstance(default, list):
            print('Default to int or float question must be list!')
            quit(1)
    if typefunc == str and autocomplete:
        readline.set_completer_delims(' \t\n;')
        readline.parse_and_bind("tab: complete")    # activate autocomplete
    else:
        readline.parse_and_bind("tab: ")            # deactivate autocomplete

    while True:
        s = question
        if default is not None:
            if typefunc == bool or typefunc == str:
                s += ' [%s]' % (str(default))
            elif typefunc == int or typefunc == float:
                s += ' ['
                for i in default:
                    s += str(i) + ' '
                s = s[:-1] + ']'
        if typefunc == str and autocomplete:
            s += ' (autocomplete enabled)'
        if typefunc == int and ranges:
            s += ' (range comprehension enabled)'
        s += ' '

        line = input(s)
        line = re.sub(r'#.*$', '', line).strip()
        if not typefunc == str:
            line = line.lower()

        if line == '' or line == '\n':
            if default is not None:
                KEYSTROKES.write(line + ' ' * (40 - len(line)) + ' #' + s + '\n')
                return default
            else:
                continue

        if typefunc == bool:
            posresponse = ['y', 'yes', 'true', 't', 'ja', 'si', 'yea', 'yeah', 'aye', 'sure', 'definitely']
            negresponse = ['n', 'no', 'false', 'f', 'nein', 'nope']
            if line in posresponse:
                KEYSTROKES.write(line + ' ' * (40 - len(line)) + ' #' + s + '\n')
                return True
            elif line in negresponse:
                KEYSTROKES.write(line + ' ' * (40 - len(line)) + ' #' + s + '\n')
                return False
            else:
                print('I didn''t understand you.')
                continue

        if typefunc == str:
            KEYSTROKES.write(line + ' ' * (40 - len(line)) + ' #' + s + '\n')
            return line

        if typefunc == float:
            # float will be returned as a list
            f = line.split()
            try:
                for i in range(len(f)):
                    f[i] = typefunc(f[i])
                KEYSTROKES.write(line + ' ' * (40 - len(line)) + ' #' + s + '\n')
                return f
            except ValueError:
                print('Please enter floats!')
                continue

        if typefunc == int:
            # int will be returned as a list
            f = line.split()
            out = []
            try:
                for i in f:
                    if ranges and '~' in i:
                        q = i.split('~')
                        for j in range(int(q[0]), int(q[1]) + 1):
                            out.append(j)
                    else:
                        out.append(int(i))
                KEYSTROKES.write(line + ' ' * (40 - len(line)) + ' #' + s + '\n')
                return out
            except ValueError:
                if ranges:
                    print('Please enter integers or ranges of integers (e.g. "-3~-1  2  5~7")!')
                else:
                    print('Please enter integers!')
                continue

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================


def label_valid(label):
    # allowed format: letter followed by letters and numbers and underscore
    if '__' in label:
        return False
    if label == 'F' or label == 'x':
        return False
    if re.match(r'^[a-zA-Z][a-zA-Z0-9_]*$', label) is None:
        return False
    else:
        return True

# ===================================================


def print_reactions(rate_matrix, specmap):
    # get longest label
    n = 3
    for i in range(len(rate_matrix)):
        if len(specmap[i]) > n:
            n = len(specmap[i])
        for j in range(len(rate_matrix)):
            if len(rate_matrix[i][j]) > n:
                n = len(rate_matrix[i][j])
    formatstring = ' %%%is ' % (n)
    # construct string
    string = ' ' * (n + 2) + '|'
    for i in range(len(rate_matrix)):
        string += formatstring % (specmap[i])
    string += '\n'
    string += '-' * (n + 2) + '+' + '-' * ((n + 2) * len(rate_matrix)) + '\n'
    for i in range(len(rate_matrix)):
        string += formatstring % (specmap[i]) + '|'
        for j in range(len(rate_matrix)):
            label = rate_matrix[i][j]
            if label == '':
                label = '.'
            string += formatstring % (label)
        string += '\n'
    string += '\n(Initial species: rows; Final species: columns)\n'
    return string

# ===================================================


def nullspace(a, rtol=1e-5):
    u, s, v = numpy.linalg.svd(a)
    rank = (s > rtol * s[0]).sum()
    return rank, v[rank:].T.copy()

# ===================================================


def get_cycles(rate_matrix):
    # get nspec and nreact
    nspec = len(rate_matrix)
    nreact = 0
    for i in range(nspec):
        for j in range(nspec):
            if rate_matrix[i][j] != '':
                nreact += 1
    if nreact == 0:
        return -1
    # construct stochiometry matrix
    A = [[0 for i in range(nreact)] for j in range(nspec)]
    ireact = 0
    for i in range(nspec):
        for j in range(nspec):
            if rate_matrix[i][j] != '':
                A[i][ireact] = -1
                A[j][ireact] = +1
                ireact += 1
    # pprint.pprint(A)
    # get nullspace of A
    if NUMPY:
        rank, null = nullspace(A)
        nullrank = len(A[0]) - rank
        if nullrank > 0:
            print('  The reaction network contains %i %s!' % (nullrank, ['cycles', 'cycle'][nullrank == 1]))
        return nullrank
    else:
        print('  Hint: Cannot check for cycles without NUMPY!')
        return -1

# ===================================================


def check_pop_file(content):
    # checks whether the content of a populations file is valid.
    # checks number of columns, time range and ordering, whether all rows have same number of columns
    # also extracts the numerical data
    data = []
    ncol = -1
    maxtime = -1.
    for line in content:
        line = re.sub(r'#.*$', '', line)
        if line == '\n':
            continue
        s = line.split()
        # check time
        time = float(s[0])
        if maxtime == -1. and time != 0.:
            print('  Time does not start at zero!')
            return False, 0, 0, []
        if time < 0.:
            print('  Negative times detected!')
            return False, 0, 0, []
        if time < maxtime:
            print('  Times not ordered!')
            return False, 0, 0, []               # TODO: maybe this check is not necessary
        maxtime = time
        # check data
        d = [float(i) for i in s]
        # if any( [i<0. for i in d] ):
        # print('  Negative populations detected!')
        # return False,0,0,[]
        col = len(d)
        if ncol == -1:
            ncol = col
        elif ncol != col:
            print('  Inconsistent number of columns detected!')
            if ncol > col:
                ncol = col
        data.append(d)
    return True, maxtime, ncol, data

# ===================================================


def get_infos():
    '''This routine asks for definitions of the kinetic model and the populations file.'''

    INFOS = {}

    print('{:#^60}'.format(''))
    print('{:#^60}'.format(' Kinetics Model '))
    print('{:#^60}'.format('') + '\n\n')
    # =========================== Define the kinetic model species ==================================
    print('{:-^60}'.format('Model Species') + '\n')
    print('''First, please specify the set of species used in your model kinetics.

Possible input:
+ <label> <label> ...   Adds one or several species to the set
- <label> <label> ...   Removes one or several species from the set
show                    Show the currently defined set of species
end                     Finish species input

Each label must be unique. Enter the labels without quotes.
''')
    species = []
    while True:
        line = question('Input:', str, 'end', False)
        s = line.split()
        if 'end' in s[0].lower():
            if len(species) == 0:
                print('  No species added yet!')
            else:
                break
        elif 'show' in s[0].lower():
            print('  Current set:  %s\n' % (species))
        elif '+' in s[0]:
            for i in s[1:]:
                if i in species:
                    print('  Species \'%s\' already in set!' % (i))
                else:
                    if label_valid(i):
                        species.append(i)
                        print('  Species \'%s\' added!' % (i))
                    else:
                        print('  Invalid label \'%s\'! Labels must be a letter followed by letters, \n  numbers and single underscores! "F" and "x" are reserved!' % (i))
        elif '-' in s[0]:
            for i in s[1:]:
                if i in species:
                    species.remove(i)
                    print('  Species \'%s\' removed!' % (i))
                else:
                    print('  Species \'%s\' not in set!' % (i))
        else:
            print('  I did not understand you.')
    print('\nFinal species set:  %s\n' % (species))
    nspec = len(species)
    specmap = {}
    for i in range(len(species)):
        specmap[species[i]] = i
        specmap[i] = species[i]
    INFOS['nspec'] = nspec
    INFOS['specmap'] = specmap


    # =========================== Define the kinetic model reactions ==================================
    print('{:-^60}'.format('Model Elementary Reactions') + '\n')
    print('''Second, please specify the set of elementary reactions in your model kinetics.

Possible input:
+ <species1> <species2> <rate_label>       Add a reaction from species1 to species2 with labelled rate constant
- <rate_label>                             Remove the reaction(s) with the given rate constant
show                                       Show the currently defined set of reactions (as directed adjacency matrix)
end                                        Finish reaction input

Each rate label must be unique.
''')
    rate_matrix = [['' for i in range(nspec)] for j in range(nspec)]
    rateset = set()
    while True:
        line = question('Input:', str, 'end', False)
        s = line.split()
        if 'end' in s[0].lower():
            if len(rateset) == 0:
                print('  No reactions added yet!')
            else:
                break
        elif 'show' in s[0].lower():
            print(print_reactions(rate_matrix, specmap))
        elif '+' in s[0]:
            if len(s) < 4:
                print('Please write "+ species1 species2 ratelabel"!')
                continue
            if s[1] == s[2]:
                print('  Species labels identical! No reaction added.')
                continue
            if s[1] in specmap and s[2] in specmap and not s[3] in specmap:
                if rate_matrix[specmap[s[1]]][specmap[s[2]]] != '':
                    print('Please remove rate constant %s first!' % (rate_matrix[specmap[s[1]]][specmap[s[2]]]))
                    continue
                if not s[3] in rateset:
                    if label_valid(s[3]):
                        rateset.add(s[3])
                        rate_matrix[specmap[s[1]]][specmap[s[2]]] = s[3]
                        rank = get_cycles(rate_matrix)
                        print('  Reaction from \'%s\' to \'%s\' with rate label \'%s\' added!' % (s[1], s[2], s[3]))
                    else:
                        print('  Invalid label \'%s\'! Labels must be a letter followed by letters, numbers and single underscores!' % (s[3]))
                else:
                    print('  Rate label \'%s\' already defined!' % (s[3]))
                    #anyways=question('Do you want to add it anyways (i.e., use two reactions with same rate constant)?',bool,False)
                    # if anyways:
                    rate_matrix[specmap[s[1]]][specmap[s[2]]] = s[3]
                    rank = get_cycles(rate_matrix)
                    print('  Reaction from \'%s\' to \'%s\' with rate label \'%s\' added!' % (s[1], s[2], s[3]))
            else:
                if not s[1] in specmap:
                    print('  Species \'%s\' not defined!' % (s[1]))
                if not s[2] in specmap:
                    print('  Species \'%s\' not defined!' % (s[2]))
                if s[3] in specmap:
                    print('  Label \'%s\' already used for a species!' % (s[3]))
        elif '-' in s[0]:
            if s[1] in rateset:
                rateset.remove(s[1])
                for i in range(nspec):
                    for j in range(nspec):
                        if rate_matrix[i][j] == s[1]:
                            rate_matrix[i][j] = ''
                            rank = get_cycles(rate_matrix)
            else:
                print('  Rate label \'%s\' not defined!' % (s[1]))
        else:
            print('  I did not understand you.')
    print('\nFinal reaction network:')
    print(print_reactions(rate_matrix, specmap))
    INFOS['rateset'] = rateset
    INFOS['rate_matrix'] = rate_matrix
    INFOS['rank'] = rank



    ratemap = {}
    for i, x in enumerate(sorted(rateset)):
        ratemap[i] = x
        ratemap[x] = i
    INFOS['nrates'] = len(rateset)
    INFOS['ratemap'] = ratemap


    # make list of reactions
    rates = []
    for irate in range(INFOS['nrates']):
        el = INFOS['ratemap'][irate]
        r = []
        for ir, row in enumerate(INFOS['rate_matrix']):
            for ic, element in enumerate(row):
                if element == el:
                    r.append((ir, ic))
        rates.append(r)
    INFOS['rates'] = rates



    # =========================== Define the kinetic model initial conditions ==================================
    print('{:-^60}'.format('Model Initial Conditions') + '\n')
    print('''Third, please specify species with non-zero initial populations.

Possible input:
+ <species>       Declare species to have non-zero initial population
- <species>       Remove species from the set of non-zero initial populations
show              Show the currently defined non-zero initial populations
end               Finish initial condition input
''')
    initset = set()
    while True:
        line = question('Input:', str, 'end', False)
        s = line.split()
        if 'end' in s[0].lower():
            if len(initset) == 0:
                print('  No species with non-zero initial population yet!')
            else:
                break
        elif 'show' in s[0].lower():
            print('  Current set:  %s\n' % (list(initset)))
        elif '+' in s[0]:
            for i in s[1:]:
                if i not in specmap:
                    print('  Species \'%s\' not defined!' % (i))
                elif i in initset:
                    print('  Species \'%s\' already in set!' % (i))
                else:
                    initset.add(i)
                    print('  Species \'%s\' added!' % (i))
        elif '-' in s[0]:
            for i in s[1:]:
                if i in initset:
                    initset.remove(i)
                    print('  Species \'%s\' removed!' % (i))
                else:
                    print('  Species \'%s\' not in set!' % (i))
        else:
            print('  I did not understand you.')
    print('\nFinal initial species set:  %s\n' % (list(initset)))
    INFOS['initset'] = initset

    #
    INFOS['initial'] = []
    for i in sorted(INFOS['initset']):
        INFOS['initial'].append(INFOS['specmap'][i])
    INFOS['ninitial'] = len(INFOS['initial'])





    print('{:#^60}'.format(''))
    print('{:#^60}'.format(' Fitting Data '))
    print('{:#^60}'.format('') + '\n\n')

    # =========================== Bootstrapping or not ==================================
    print('{:-^60}'.format('Operation mode') + '\n')
    print('''This script can work with the following output:
* pop.out (file from populations.py)
* bootstrap_data/ (directory from populations.py)
Using only the pop.out allows fitting and obtaining time constants.
Using the bootstrap data instead additionally allows for realistic error estimates.
''')
    INFOS['do_bootstrap'] = question('Do you want to use bootstrap data?', bool, False)

    if INFOS['do_bootstrap']:
        INFOS['bootstrap_cycles'] = question('How many bootstrap samples?', int, [100])[0]



    # =========================== Define the data file ==================================
    print('\n' + '{:-^60}'.format('Population data file') + '\n')
    if INFOS['do_bootstrap']:
        print('''Please specify the path to the bootstrap data directory (as generated by populations.py).\n''')
        while True:
            bsdir = question('Bootstrap data directory:', str, 'bootstrap_data/', True)
            if not os.path.isdir(bsdir):
                print('  Directory not found!')
                continue
            ls = os.listdir(bsdir)
            VALID = True
            MAXTIME = []
            NCOL = []
            DATA = []
            for i in ls:
                if 'pop_' not in i:
                    continue
                popfile = os.path.join(bsdir, i)
                content = readfile(popfile)
                valid, maxtime, ncol, data = check_pop_file(content)
                if not valid:
                    print('  File format not valid (%s)!' % popfile)
                    VALID = False
                else:
                    MAXTIME.append(maxtime)
                    NCOL.append(ncol)
                    DATA.append(data)
            if not VALID:
                continue
            s = set(MAXTIME)
            if len(s) > 1:
                print('Bootstrap files have different maximum time!')
                continue
            s = set(NCOL)
            if len(s) > 1:
                print('Bootstrap files have different number of columns')
                continue
            s = set([len(i) for i in DATA])
            if len(s) > 1:
                print('Bootstrap files have different number of time steps!')
                continue
            maxtime = MAXTIME[0]
            ncol = NCOL[0]
            data = DATA
            popfile = bsdir
            INFOS['ntraj'] = len(DATA)
            print('  Detected maximal time of %7.1f fs and %i columns (time plus %i data columns).' % (maxtime, ncol, ncol - 1))
            break

        print
        INFOS['write_bootstrap_fits'] = question('Do you want to write fitting curves for all bootstrap cycles?', bool, False)
    else:
        print('''Please specify the path to the population data file (as generated by populations.py).\n''')
        while True:
            popfile = question('Populations file:', str, 'pop.out', True)
            if not os.path.isfile(popfile):
                print('  File not found!')
                continue
            content = readfile(popfile)
            valid, maxtime, ncol, data = check_pop_file(content)
            if not valid:
                print('  File format not valid!')
                continue
            else:
                print('  Detected maximal time of %7.1f fs and %i columns (time plus %i data columns).' % (maxtime, ncol, ncol - 1))
                break
    INFOS['maxtime'] = maxtime
    INFOS['ncol'] = ncol
    INFOS['data'] = data
    INFOS['popfile'] = os.path.abspath(popfile)

    # =========================== Define the data -- species mapping ==================================
    print('\n' + '{:-^60}'.format('Population-to-Species Mapping for Fit') + '\n')
    print('''Please specify which model species should be fitted to which data file columns.
For example, you can fit the label 'S0' to column 2:
  S0 = 2
You can also fit the sum of two species to a column:
  T1 T2 = 5
You can also fit a species to the sum of several columns:
  T_all = 5 6 7
You can even fit a sum of species to a sum of columns:
  T1 T2 = 5 6 7

On the right side, "~" can be used to indicate ranges:
  T1 T2 = 5~9

Possible input:
<species1> <species2> ... = <columns1> <column2> ...            Set one mapping
show                                                            Show mapping
end                                                             Finish mapping input
reset                                                           Redo the mapping input

Each species label must be used at most once.
Each column number (except for \'1\', which denotes the time) must be used at most once.
''')
    print('Set of species:        %s' % (species))
    columns = [i for i in range(2, ncol + 1)]
    print('Set of column numbers: %s' % (columns))

    species_groups = []
    columns_groups = []
    ngroups = 0
    while True:
        line = question('Input:', str, 'end', False)
        s = line.split()
        if 'end' in s[0].lower():
            if ngroups == 0:
                print('  No valid input yet!')
            else:
                break
        elif 'show' in s[0].lower():
            print('  Current mapping groups:')
            for i in range(ngroups):
                string = '    '
                for j in species_groups[i]:
                    string += ' %s ' % (j)
                string += ' = '
                for j in columns_groups[i]:
                    string += ' %i ' % (j)
                print(string)
        elif ' = ' in line:
            if s[0] == '=' or s[-1] == '=':
                print('  Invalid input! Put species labels to the left of \'=\' and column number to the right!')
                continue
            do_species = True
            valid = True
            new_species_group = []
            new_columns_group = []
            for i in s:
                if i == '=':
                    if not do_species:
                        print('More than 1 "=" used!')
                        valid = False
                        break
                    else:
                        do_species = False
                        continue
                if do_species:
                    if i not in species:
                        print('  Species label \'%s\' not defined!' % (i))
                        valid = False
                        break
                    if any([i in j for j in species_groups]):
                        print('  Species label \'%s\' already assigned!' % (i))
                        valid = False
                        break
                    if i in new_species_group:
                        print('  Species label \'%s\' used twice!' % (i))
                        valid = False
                        break
                    new_species_group.append(i)
                else:
                    try:
                        if '~' in i:
                            ii = []
                            q = i.split('~')
                            for j in range(int(q[0]), int(q[1]) + 1):
                                ii.append(j)
                        else:
                            ii = [int(i)]
                    except ValueError:
                        print('  Could not understand!')
                        valid = False
                        break
                    for i in ii:
                        if i not in columns:
                            print('  Column number %i not in data file!' % (i))
                            valid = False
                            break
                        if any([i in j for j in columns_groups]):
                            print('  Column number %i already assigned!' % (i))
                            valid = False
                            break
                        if i in new_columns_group:
                            print('  Columns number %i used twice!' % (i))
                            valid = False
                            break
                        new_columns_group.append(i)
            if valid:
                species_groups.append(new_species_group)
                columns_groups.append(new_columns_group)
                ngroups += 1
        elif 'reset' in s[0].lower():
            species_groups = []
            columns_groups = []
            ngroups = 0
            print('  Mappings reset! Please repeat input!')
    print('Final mappings:')
    for i in range(ngroups):
        string = '    '
        for j in species_groups[i]:
            string += ' %s ' % (j)
        string += ' = '
        for j in columns_groups[i]:
            string += ' %i ' % (j)
        print(string)
    INFOS['species_groups'] = species_groups
    INFOS['columns_groups'] = columns_groups
    INFOS['ngroups'] = ngroups

    # get the summation defs
    summation = []
    for i in INFOS['species_groups']:
        s = []
        for j in i:
            ind = INFOS['specmap'][j]
            s.append(ind)
        summation.append(s)
    INFOS['summation'] = summation


    print('\n')
    print('{:#^60}'.format(''))
    print('{:#^60}'.format(' Fitting procedure '))
    print('{:#^60}'.format(''))

    # =========================== Initial guesses ==================================

    print('\n\n' + '{:-^40}'.format('Initial guesses') + '\n')

    def print_guesses(INFOS, y0, p0):
        for i in range(INFOS['nrates']):
            name = INFOS['ratemap'][i]
            t = p0[i]
            print('  time constant ( %-12s ): %12.4f fs' % (name, 1. / t))
        for i in range(INFOS['ninitial']):
            name = INFOS['specmap'][INFOS['initial'][i]]
            t = y0[i]
            print('  initial pop   ( %-12s ): %12.4f' % (name, t))
        print

    # rate guesses
    p0 = [1. / (100 + 20 * i) for i in range(INFOS['nrates'])]

    # initial guesses
    y0 = [1. for i in range(INFOS['ninitial'])]

    print('''Please check the initial guesses for the parameters

Possible input:
label = value     Set an initial guess (detects type automatically and computes k=1/t for rates)
show              Show the currently defined non-zero initial populations
end               Finish initial condition input
''')

    print_guesses(INFOS, y0, p0)
    while True:
        line = question('Input:', str, 'end', False)
        s = line.split()
        if 'show' in s[0].lower():
            print_guesses(INFOS, y0, p0)
        elif 'end' in s[0].lower():
            break
        elif '=' in line:
            if len(s) != 3:
                print('  Format must be "label = value" (including spaces)')
                continue
            if s[0] in INFOS['ratemap']:
                val = 1. / float(s[2])
                ind = INFOS['ratemap'][s[0]]
                p0[ind] = val
            elif s[0] in INFOS['specmap']:
                if not s[0] in INFOS['initset']:
                    print('  Initial population of "%s" cannot be non-zero.' % (s[0]))
                    continue
                val = float(s[2])
                ind = INFOS['specmap'][s[0]]
                ind = INFOS['initial'].index(ind)
                y0[ind] = val
            else:
                print('  Unknown label "%s"' % (s[0]))
                continue
        else:
            print('  Could not understand!')
    INFOS['p0'] = p0
    INFOS['y0'] = y0
    print('Final guess parameters:')
    print_guesses(INFOS, y0, p0)

    # =========================== Optimize initial pops ==================================
    print('\n\n' + '{:-^40}'.format('Optimize initial populations') + '\n')
    INFOS['opt_init'] = question('Do you want to optimize the initial populations (otherwise only the rates)?', bool, True)

    # =========================== Positive rates ==================================

    print('\n\n' + '{:-^40}'.format('Constrained optimization') + '\n')
    INFOS['bounds'] = question('Do you want to restrict all rates/initial populations to be non-negative?', bool, True)



    return INFOS

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================


class globalfunction():

    # ------------------------------------------------
    def __init__(self, nspecies, ratedefs, initialdefs, sumdefs, Tarray, p0, y0):

        # Number of species
        self.nspecies = nspecies

        # Setup the summation
        for j in sumdefs:
            for i in j:
                if not 0 <= i < self.nspecies:
                    print('Illegal summation definition: %s' % (ratedefs))
                    sys.exit(1)
        self.sumdefs = sumdefs

        # check the rate definitions
        for j in ratedefs:
            for i in j:
                if not 0 <= i[0] < self.nspecies:
                    print('Illegal rate definition: %s' % (ratedefs))
                    sys.exit(1)
                if not 0 <= i[1] < self.nspecies:
                    print('Illegal rate definition: %s' % (ratedefs))
                    sys.exit(1)
        self.ratedefs = ratedefs
        self.nrates = len(ratedefs)

        # make initial rate matrix
        if not len(p0) == self.nrates:
            print('Initial parameters must have same length as rate definitions!')
            sys.exit(1)
        self.p = p0
        self.set_ratematrix(p0)

        # check initial definitions
        for i in initialdefs:
            if not 0 <= i < self.nspecies:
                print('Illegal initial definition: %s' % (ratedefs))
                sys.exit(1)
        self.initdefs = initialdefs
        self.ninit = len(self.initdefs)

        # check initial data
        if not self.ninit == len(y0):
            print('Initial data must have same length as initial definitions!')
            sys.exit(1)
        self.y0 = y0
        self.set_initvector(self.y0)

        # make time array
        self.T = Tarray
        self.tmax = max(Tarray)

        # fill dictionary with function values
        self.vals = {0.: self.y}
        self.fill_vals(self.T)

    # ------------------------------------------------
    def set_ratematrix(self, p):

        self.p = p
        ratematrix = [[0. for i in range(self.nspecies)] for j in range(self.nspecies)]
        for i, ind1 in enumerate(self.ratedefs):
            for ind in ind1:
                ratematrix[ind[1]][ind[0]] += p[i]
                ratematrix[ind[0]][ind[0]] -= p[i]
        self.rates = np.array(ratematrix)

    # ------------------------------------------------
    def set_initvector(self, y):

        self.y = [0. for i in range(self.nspecies)]
        for i, ind in enumerate(self.initdefs):
            self.y[ind] = self.y0[i]

    # ------------------------------------------------
    def fun(self, t, y):
        return np.dot(self.rates, y)

    # ------------------------------------------------
    def fill_vals(self, T):

        # initialize Runge-Kutta
        RK = spint.RK45(self.fun, 0., self.y, self.tmax + 0.1)

        # perform initial step
        RK.step()
        interpol = RK.dense_output()

        # perform the steps in order
        for it, t in enumerate(sorted(T)):
            if t in self.vals:
                continue
            while True:
                if not interpol.t_min <= t < interpol.t_max:
                    try:
                        RK.step()
                        interpol = RK.dense_output()
                    except RuntimeError:
                        if not t == self.tmax:
                            print('Error 1')
                            sys.exit(1)
                        else:
                            break
                else:
                    break
            y = interpol(t)
            self.vals[t] = y

    # ------------------------------------------------
    def __call__(self, t, *params):
        # check for modified parameters and reinitialize
        p1 = list(params[:self.nrates])
        if len(params) >= self.nrates + self.ninit:
            y1 = list(params[-self.ninit:])
        else:
            y1 = self.y0
        if not y1 == self.y0 or not p1 == self.p:
            if not y1 == self.y0:
                self.set_initvector(y1)
            if not p1 == self.p:
                self.set_ratematrix(p1)
            self.vals = {0.: self.y}
            self.fill_vals(self.T)
        # evaluate
        state = int(t // self.tmax)
        tnew = t % self.tmax
        if tnew not in self.vals:
            self.fill_vals([tnew])
        s = 0.
        for i in self.sumdefs[state]:
            s += self.vals[tnew][i]
        return s

    # ------------------------------------------------
    def call_array(self, T, *params):
        # print(parameters)
        #string='Parameters: '
        # for i in list(params):
        #string+='%12.9f  ' % i
        # string+='\n'
        # sys.stderr.write(string)
        # check for modified parameters and reinitialize
        p1 = list(params[:self.nrates])
        if len(params) >= self.nrates + self.ninit:
            y1 = list(params[-self.ninit:])
        else:
            y1 = self.y0
        if not y1 == self.y0 or not p1 == self.p:
            if not y1 == self.y0:
                self.set_initvector(y1)
            if not p1 == self.p:
                self.set_ratematrix(p1)
            self.vals = {0.: self.y}
            self.fill_vals(self.T)
        # evaluate
        missingT = []
        for t in T:
            state = int(t // self.tmax)
            tnew = t % self.tmax
            if tnew not in self.vals:
                missingT.append(tnew)
        if missingT:
            self.fill_vals(missingT)
        output = []
        for t in T:
            state = int(t // self.tmax)
            tnew = t % self.tmax
            s = 0.
            for i in self.sumdefs[state]:
                s += self.vals[tnew][i]
            output.append(s)
        return output

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================


def create_bootstrap_data(Tdata1, RNGarray, INFOS):

    # get the fitting data
    nsteps = len(Tdata1) - 1
    Tdata = [0. for i in range(nsteps * INFOS['ngroups'])]
    Y1 = [0. for i in range(nsteps * INFOS['ngroups'])]
    Y2 = [0. for i in range(nsteps * INFOS['ngroups'])]
    for igroup in range(INFOS['ngroups']):
        cols = INFOS['columns_groups'][igroup]
        for istep in range(nsteps):
            for itraj in RNGarray:
                s = 0.
                for col in cols:
                    s += INFOS['data'][itraj][istep][col - 1]
                Y1[istep + igroup * nsteps] += s
                Y2[istep + igroup * nsteps] += s**2
            Tdata[istep + igroup * nsteps] = Tdata1[istep] + igroup * INFOS['maxtime']

    Ydata = [Y1[i] / INFOS['ntraj'] for i in range(nsteps * INFOS['ngroups'])]
    Yerr = [math.sqrt(Y2[i] / INFOS['ntraj'] - Ydata[i]**2) for i in range(nsteps * INFOS['ngroups'])]
    for i in range(nsteps * INFOS['ngroups']):
        if Yerr[i] == 0:
            Yerr[i] = 0.001

    return Tdata, Ydata, Yerr

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================


def make_fit(INFOS):

    print('\n' + '{:#^60}'.format(' Fitting ') + '\n')

    # rate guesses
    p0 = deepcopy(INFOS['p0'])

    # initial guesses
    y0 = deepcopy(INFOS['y0'])

    # bounds
    if INFOS['bounds']:
        bounds = (1e-6, np.inf)
    else:
        bounds = (-np.inf, np.inf)



    # ------------------------------------------------

    # PREPARE the data
    if INFOS['do_bootstrap']:
        # get time data
        Tdata1 = []
        for line in INFOS['data'][0]:
            Tdata1.append(line[0])
        # check the time data for consistency
        for data in INFOS['data'][1:]:
            Tdata2 = []
            for line in data:
                Tdata2.append(line[0])
            if not Tdata1 == Tdata2:
                print('Time data inconsistent!')
                sys.exit(1)

        RNGarray = [i for i in range(INFOS['ntraj'])]
        Tdata, Ydata, Yerr = create_bootstrap_data(Tdata1, RNGarray, INFOS)
        Yerr_absol = True

        # TODO: Could use the actual Yerr, but in some cases this gives very bad fits, so we deactivate it here:
        Yerr = None
        Yerr_absol = False


    # ------------------------------------------------
    else:
        # get time data
        Tdata1 = []
        for line in INFOS['data']:
            Tdata1.append(line[0])


        # get the fitting data
        nsteps = len(Tdata1) - 1
        Tdata = [0. for i in range(nsteps * INFOS['ngroups'])]
        Ydata = [0. for i in range(nsteps * INFOS['ngroups'])]
        for igroup in range(INFOS['ngroups']):
            cols = INFOS['columns_groups'][igroup]
            for istep in range(nsteps):
                # print(istep+igroup*nsteps)
                for col in cols:
                    Ydata[istep + igroup * nsteps] += INFOS['data'][istep][col - 1]
                Tdata[istep + igroup * nsteps] = Tdata1[istep] + igroup * INFOS['maxtime']
        Yerr = None
        Yerr_absol = False

    # ------------------------------------------------

    # initialize fitting function
    F = globalfunction(INFOS['nspec'], INFOS['rates'], INFOS['initial'], INFOS['summation'], Tdata1, p0, y0)



    print('\n' + '{:-^40}'.format(' Iterations ') + '\n')
    # get optimal parameters
    if INFOS['opt_init']:
        OPT = spopt.curve_fit(F.call_array, Tdata, Ydata, p0=p0 + y0, bounds=bounds, sigma=Yerr, absolute_sigma=Yerr_absol, verbose=2)
    else:
        OPT = spopt.curve_fit(F.call_array, Tdata, Ydata, p0=p0, bounds=bounds, sigma=Yerr, absolute_sigma=Yerr_absol, verbose=2)
    popt = OPT[0].tolist()
    OPT_orig = deepcopy(OPT)
    # print(popt)



    print('\n\n' + '{:-^40}'.format(' Final parameters ') + '\n')

    const_names = []
    const_values = []
    for i in range(INFOS['nrates']):
        name = INFOS['ratemap'][i]
        const_names.append(name)
        t = 1. / OPT[0][i]
        const_values.append(t)
        dt = np.sqrt(OPT[1][i][i]) / OPT[0][i]**2
        perc = dt / t * 100.
        print('time constant ( %-12s ): %12.4f fs +/- %12.4f fs (%7.2f %%)' % (name, t, dt, perc))
    for i in range(INFOS['ninitial']):
        name = INFOS['specmap'][INFOS['initial'][i]]
        const_names.append(name)
        if INFOS['opt_init']:
            ind = i + INFOS['nrates']
            t = OPT[0][ind]
            dt = np.sqrt(OPT[1][ind][ind])
            perc = dt / t * 100.
            print('initial pop   ( %-12s ): %12.4f    +/- %12.4f    (%7.2f %%)' % (name, t, dt, perc))
        else:
            t = y0[i]
            print('initial pop   ( %-12s ): %12.4f' % (name, t))
        const_values.append(t)
    print
    if INFOS['do_bootstrap']:
        print('''These time constants include errors that assume that the population data
is free of uncertainty. Bootstrapping analysis is following now.''')
        # print('''These time constants include errors estimated from the standard deviation
# of the original population data. Bootstrapping analysis is following now.''')
    else:
        print('''These time constants include errors that assume that the population data
is free of uncertainty. If you want to compute a more accurate error estimate, provide bootstrapping
data to this script (bootstrapping data can be prepared with populations.py).''')


    # print(function values and data together)
    string = ''
    for it, T in enumerate(Tdata):
        string += '%12.9f %12.9f %12.9f' % (T, Ydata[it], F(T, *popt))
        if Yerr:
            string += '   %12.9f' % (Yerr[it])
        string += '\n'
    print('\nRaw data and fitted functions written to "fit_results.txt".')
    writefile('fit_results.txt', string)


    # make gnuplot script
    string = gnuplot_string(INFOS, const_values)
    print('\nGNUPLOT script written to "fit_results.gp".')
    writefile('fit_results.gp', string)

    sys.stdout.flush()


    # ------------------------------------------------

    if INFOS['do_bootstrap']:
        verbose = False
        print('\n' + '{:#^60}'.format(' Bootstrapping ') + '\n')
        if INFOS['write_bootstrap_fits']:
            print('Writing individual results to %s/fit_results_%%i.txt ...\n' % (INFOS['popfile']))

        p0 = popt[:INFOS['nrates']]
        if INFOS['opt_init']:
            y0 = popt[INFOS['nrates']:]
        else:
            y0 = INFOS['y0']
        F = globalfunction(INFOS['nspec'], INFOS['rates'], INFOS['initial'], INFOS['summation'], Tdata1, p0, y0)

        results = []
        if not verbose:
            string = ' Cycle '
            for i in const_names:
                string += '%12s ' % i
            string += ' Time'
            print(string)
        begintime = datetime.datetime.now()
        for iboot in range(INFOS['bootstrap_cycles']):
            try:
                if verbose:
                    print('\n{:.^30}'.format(' Cycle %i ' % (iboot + 1)) + '\n')
                    verb = 2
                else:
                    verb = 0
                RNGarray = [random.randint(0, INFOS['ntraj'] - 1) for i in range(INFOS['ntraj'])]
                Tdata, Ydata, Yerr = create_bootstrap_data(Tdata1, RNGarray, INFOS)
                Yerr = None
                Yerr_absol = False
                y1 = deepcopy(y0)
                p1 = deepcopy(p0)
                if INFOS['opt_init']:
                    OPT = spopt.curve_fit(F.call_array, Tdata, Ydata, p0=p1 + y1, bounds=bounds, sigma=Yerr, absolute_sigma=Yerr_absol, verbose=verb)
                else:
                    OPT = spopt.curve_fit(F.call_array, Tdata, Ydata, p0=p1, bounds=bounds, sigma=Yerr, absolute_sigma=Yerr_absol, verbose=verb)
                pboot = OPT[0].tolist()

                if verbose:
                    print
                R = {}
                string = '%6i ' % (iboot + 1)
                for i in range(INFOS['nrates']):
                    name = INFOS['ratemap'][i]
                    t = 1. / OPT[0][i]
                    dt = np.sqrt(OPT[1][i][i]) / OPT[0][i]**2
                    perc = dt / t * 100.
                    if verbose:
                        print('time constant ( %-12s ): %12.4f fs +/- %12.4f fs (%7.2f %%)' % (name, t, dt, perc))
                    else:
                        string += '%12.4f ' % (t)
                    R[name] = t
                for i in range(INFOS['ninitial']):
                    name = INFOS['specmap'][INFOS['initial'][i]]
                    if INFOS['opt_init']:
                        ind = i + INFOS['nrates']
                        t = OPT[0][ind]
                        dt = np.sqrt(OPT[1][ind][ind])
                        perc = dt / t * 100.
                        if verbose:
                            print('initial pop   ( %-12s ): %12.4f    +/- %12.4f    (%7.2f %%)' % (name, t, dt, perc))
                        else:
                            string += '%12.4f ' % (t)
                    else:
                        t = y0[i]
                        if verbose:
                            print('initial pop   ( %-12s ): %12.4f' % (name, t))
                        else:
                            string += '%12.4f ' % (t)
                    R[name] = t
                results.append(R)
                deltatime = datetime.datetime.now() - begintime
                begintime = datetime.datetime.now()
                if not verbose:
                    string += ' %s' % deltatime
                    print(string)
                else:
                    print('Time: ', deltatime)

                if INFOS['write_bootstrap_fits']:
                    string = ''
                    for it, T in enumerate(Tdata):
                        string += '%12.9f %12.9f %12.9f' % (T, Ydata[it], F(T, *pboot))
                        if Yerr:
                            string += '   %12.9f' % (Yerr[it])
                        string += '\n'
                    filename = os.path.join(INFOS['popfile'], 'fit_results_%i.txt' % iboot)
                    writefile(filename, string)

                sys.stdout.flush()
            except KeyboardInterrupt:
                print('Aborted, going to final analysis...')
                time.sleep(0.5)
                break


        # final analysis
        print('\n>>>>>>>>>>>>> Finished the bootstrapping cycles ...')
        string_all = ''

        final_errors = {}
        for ikey, key in enumerate(const_names):
            if key in INFOS['ratemap']:
                string = '\n{:-^110}'.format(' Analysis for time constant "%s" ' % key) + '\n'
            elif key in INFOS['specmap']:
                string = '\n{:-^105}'.format(' Analysis for initial population "%s" ' % key) + '\n'

            data = []
            for i in results:
                data.append(i[key])
            mini = min(data)
            maxi = max(data)
            mean_a = mean_arith(data)
            stdev_a = stdev_arith(data, mean_a)
            mean_g = mean_geom(data)
            stdev_g = stdev_geom(data, mean_g)
            final_errors[key] = stdev_a

            string += '''
    Arithmetic analysis:          %12.4f +/- %12.4f
                                            ( +/-   %8.2f %%)

    Geometric analysis:           %12.4f  +  %12.4f  -  %12.4f
                                            (  +    %8.2f %%  -    %8.2f %%)

    Minimum and maximum:          %12.4f      and       %12.4f

    Histogram:
    ==========''' % (mean_a,
                     stdev_a,
                     stdev_a / mean_a * 100.,
                     mean_g,
                     mean_g * (stdev_g - 1.),
                     -mean_g * (1. / stdev_g - 1.),
                     (stdev_g - 1.) * 100.,
                     -(1. / stdev_g - 1.) * 100.,
                     mini,
                     maxi
                     )
            print(string)
            string_all += string + '\n'

            # make histogram
            nhisto = 11
            whisto = 3
            bins = [mean_g * stdev_g**((float(i) - nhisto / 2) / (nhisto / 2) * whisto) for i in range(nhisto - 1)]
            h = histogram(bins)
            hout = [0 for i in range(nhisto)]
            for x in data:
                i = h.put(x)
                hout[i] += 1

            # draw histogram
            height = 12
            dx = float(height) / max(hout)
            string = ''
            for ih in range(height):
                for ib in hout:
                    x = ib * dx
                    if x >= (height - ih):
                        s = '#'
                    else:
                        s = ' '
                    string += ' %5s' % s
                string += '| %i\n' % (int((height - ih) / dx))
            string += '     '
            for i in bins:
                string += '   |  '
            string += '\n   '
            for i in bins:
                if i < 10.:
                    string += ' %5.3f' % (i)
                elif i < 100.:
                    string += ' %5.2f' % (i)
                elif i < 1000.:
                    string += ' %5.1f' % (i)
                else:
                    string += ' %5i' % (i)
            print(string)
            string_all += string + '\n'



        print('\n\n' + '{:-^40}'.format(' Final parameters ') + '\n')

        const_names = []
        for i in range(INFOS['nrates']):
            name = INFOS['ratemap'][i]
            const_names.append(name)
            t = 1. / OPT_orig[0][i]
            dt = final_errors[name]
            perc = dt / t * 100.
            print('time constant ( %-12s ): %12.4f fs +/- %12.4f fs (%7.2f %%)' % (name, t, dt, perc))
        for i in range(INFOS['ninitial']):
            name = INFOS['specmap'][INFOS['initial'][i]]
            const_names.append(name)
            if INFOS['opt_init']:
                ind = i + INFOS['nrates']
                t = OPT_orig[0][ind]
                dt = final_errors[name]
                perc = dt / t * 100.
                print('initial pop   ( %-12s ): %12.4f    +/- %12.4f    (%7.2f %%)' % (name, t, dt, perc))
            else:
                t = y0[i]
                print('initial pop   ( %-12s ): %12.4f' % (name, t))
        print



        # write results to file
        string = '\n\nFull data:\n\n'
        string += '%10s ' % 'Sample'
        for key in const_names:
            string += '%12s ' % key
        string += '\n'
        for i, c in enumerate(results):
            string += '%10i ' % (i + 1)
            for key in const_names:
                string += '%12.4f ' % c[key]
            string += '\n'
        string_all += string

        print('\nOutput (analysis and full fitted data) written to "fit_bootstrap.txt".')
        writefile('fit_bootstrap.txt', string_all)






# ======================================== #
def mean_arith(data):
    s = 0.
    for i in data:
        s += i
    return s / len(data)

# ======================================== #


def stdev_arith(data, mean=-9999):
    if mean == -9999:
        m = mean_arith(data)
    else:
        m = mean
    s = 0.
    for i in data:
        s += (i - m)**2
    s = s / (len(data) - 1)
    return math.sqrt(s)

# ======================================== #


def mean_geom(data):
    s = 0.
    for i in data:
        s += math.log(i)
    s = s / len(data)
    return math.exp(s)

# ======================================== #


def stdev_geom(data, mean=-9999):
    if mean == -9999:
        m = math.log(mean_geom(data))
    else:
        m = math.log(mean)
    s = 0.
    for i in data:
        s += (math.log(i) - m)**2
    s = s / (len(data) - 1)
    return math.exp(math.sqrt(s))


# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

class histogram:
    def __init__(self, binlist):
        '''binlist must be a list of floats
    Later, all floats x with binlist[i-1]<x<=binlist[i] will return i'''
        self.binlist = sorted(binlist)
        self.len = len(binlist) + 1

    def put(self, x):
        i = 0
        for el in self.binlist:
            if x <= el:
                return i
            else:
                i += 1
        return i

    def __repr__(self):
        s = 'Histogram object: '
        for i in self.binlist:
            s += '%f ' % (i)
        return s

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================


def gnuplot_string(INFOS, const_values):
    colorpalette = rgbcolor([INFOS['ngroups']])

    # header
    string = '#\n'
    string += '# +' + '-' * 60 + '+\n'
    string += '# |' + '{: ^60}'.format('Fit plotting script') + '|\n'
    string += '# +' + '-' * 60 + '+\n#\n#\n'

    # add as comment the definitions of the kinetic model
    string += '# *** Definition of the kinetic model: ***\n'
    s = print_reactions(INFOS['rate_matrix'], INFOS['specmap']).splitlines()
    for line in s:
        string += '#' + line + '\n'
    string += '#\n#\n'
    string += '# *** Species and initial value: ***\n'
    for i in range(INFOS['nspec']):
        if INFOS['specmap'][i] in INFOS['initset']:
            initvalue = INFOS['specmap'][i] + '__0'
        else:
            initvalue = '0'
        string += '# %s' % (INFOS['specmap'][i]) + ' ' * (10 - len(INFOS['specmap'][i])) + initvalue + '\n'
    string += '#\n#\n'
    string += '# *** Reaction rates: ***\n'
    for i in INFOS['rateset']:
        string += '# %s\n' % (i)
    string += '#\n#\n'
    string += '\n\n# ========================================================\n'

    # add gnuplot global options
    string += '#<<\n'
    string += '# *** Gnuplot general options: ***\n'
    string += '''set xlabel "Global Fit Time Axis (fs)"
set ylabel "Population"
set xrange [0:%.2f]
set yrange [0:1]
unset key
set tmargin 1
set bmargin 4
set lmargin 8
set rmargin 22
''' % (INFOS['maxtime'] * INFOS['ngroups'])
    string += '\n'

    # add borders between subplots
    for i in range(INFOS['ngroups'] - 1):
        string += 'set arrow from %f,0 to %f,1 nohead front\n' % (INFOS['maxtime'] * (i + 1), INFOS['maxtime'] * (i + 1))

    # add label with rates
    string += '# *** Label with time constants: ***\n'
    string += 'set label 1 "'
    j = 0
    for i in range(INFOS['nrates']):
        name = INFOS['ratemap'][i]
        t = const_values[i]
        j += 1
        string += 't(%s) = %7.1f fs\\n' % (name, t)
    for i in range(INFOS['ninitial']):
        name = INFOS['specmap'][INFOS['initial'][i]]
        t = const_values[j]
        j += 1
        string += 'p0(%s) = %7.1f\\n' % (name, t)
    string += '" at %.2f,0.8 left\n' % (INFOS['maxtime'] * (0.1 + INFOS['ngroups']))
    string += '\n'

    # add labels for groups
    string += '# *** Label with groups: ***\n'
    for i in range(INFOS['ngroups']):
        string += 'set label %i "' % (i + 2)
        string += '+'.join(INFOS['species_groups'][i])
        string += '\\n'
        string += '+'.join(['$%i' % q for q in INFOS['columns_groups'][i]])
        string += '" at %.2f,0.95 center\n' % ((0.5 + i) * INFOS['maxtime'])

    # plot
    string += 'p "fit_results.txt" u 1:2 w p pt 7 ps 0.4 lc rgb "black", "" u 1:3 w l lw 2 lc rgb "red"\n\npause -1'

    # Put time stamp infos at the end
    string += '\n\n# *** Infos: ***\n'
    string += '# %s@%s\n' % (os.environ['USER'], gethostname())
    string += '# Date: %s\n' % (datetime.datetime.now())
    string += '# Current directory: %s\n\n' % (os.getcwd())

    return string


# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def main():
    '''Main routine'''

    usage = '''
python make_fit.py

This interactive script lets the user specify a kinetic model and a populations.py output file and produces a
GNUPLOT script which allows to fit the model parameters to the populations
'''
    description = ''
    parser = OptionParser(usage=usage, description=description)

    displaywelcome()
    open_keystrokes()

    # get input
    INFOS = get_infos()

    # echo input
    print('\n\n' + '{:#^60}'.format('Full input') + '\n')
    for item in sorted(INFOS):
        if not item == 'data':
            print(item, ' ' * (25 - len(item)), INFOS[item])
        elif item == 'data':
            print(item, ' ' * (25 - len(item)), '[ ... ]')
    print('')
    go_on = question('Do you want to continue?', bool, True)
    if not go_on:
        quit(0)
    print('')

    # do work
    # functionstring=get_functions_from_maxima(INFOS)
    # write_gnuscript(INFOS,functionstring)
    # write_fitting_data(INFOS)
    make_fit(INFOS)


    # finalize
    close_keystrokes()

# ======================================================================================================================


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print('\nCtrl+C makes me a sad SHARC ;-(\n')
        quit(0)
