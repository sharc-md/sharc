#!/usr/bin/env python2

#******************************************
#
#    SHARC Program Suite
#
#    Copyright (c) 2019 University of Vienna
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
#******************************************

#!/usr/bin/env python2

# TODO: G16 writes X+Y/X-Y or X/Y ??


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

# =========================================================0
# compatibility stuff

if sys.version_info[0]!=2:
    print 'This is a script for Python 2!'
    sys.exit(0)

if sys.version_info[1]<5:
    def any(iterable):
        for element in iterable:
            if element:
                return True
        return False

    def all(iterable):
        for element in iterable:
            if not element:
                return False
        return True

# ======================================================================= #

version='2.1'
versiondate=datetime.date(2019,9,1)



changelogstring='''
30.11.2017: INITIAL VERSION
- only singlets
- h, dm, grad, overlap are working

01.12.2017:
- all features of ADF interface, minus SOCs

08.01.2019:
- added "basis_external" keyword
'''

# ======================================================================= #
# holds the system time when the script was started
starttime=datetime.datetime.now()

# global variables for printing (PRINT gives formatted output, DEBUG gives raw output)
DEBUG=False
PRINT=True

# hash table for conversion of multiplicity to the keywords used in MOLCAS
IToMult={
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
IToPol={
                0: 'X', 
                1: 'Y', 
                2: 'Z', 
                'X': 0, 
                'Y': 1, 
                'Z': 2
                }

#Number of frozen core orbitals
FROZENS = {'H':  0, 'He': 0,
'Li': 1, 'Be': 1, 'B':  1, 'C':  1,  'N': 1,  'O': 1, 'F':  1, 'Ne':1,
'Na':1, 'Mg':1, 'Al':5, 'Si':5,  'P':5,  'S':5, 'Cl':5, 'Ar':5,
'K': 5, 'Ca':5,
'Sc':5, 'Ti':5, 'V': 5, 'Cr':5, 'Mn':5, 'Fe':5, 'Co':5, 'Ni':5, 'Cu':5, 'Zn':5,
'Ga':9, 'Ge':9, 'As':9, 'Se':9, 'Br':9, 'Kr':9,
'Rb':9, 'Sr':9,
'Y':14,  'Zr':14, 'Nb':14, 'Mo':14, 'Tc':14, 'Ru':14, 'Rh':14, 'Pd':14, 'Ag':14, 'Cd':14,
'In':18, 'Sn':18, 'Sb':18, 'Te':18,  'I':18, 'Xe':18,
'Cs':18, 'Ba':18,
'La':23, 
'Ce':23, 'Pr':23, 'Nd':23, 'Pm':23, 'Sm':23, 'Eu':23, 'Gd':23, 'Tb':23, 'Dy':23, 'Ho':23, 'Er':23, 'Tm':23, 'Yb':23, 'Lu':23,
         'Hf':23, 'Ta':23,  'W':23, 'Re':23, 'Os':23, 'Ir':23, 'Pt':23, 'Au':23, 'Hg':23,
'Tl':23, 'Pb':23, 'Bi':23, 'Po':23, 'At':23, 'Rn':23, 
'Fr':30, 'Ra':30,
'Ac':30, 
'Th':30, 'Pa':30,  'U':30, 'Np':30, 'Pu':30, 'Am':30, 'Cm':30, 'Bk':30, 'Cf':30, 'Es':30, 'Fm':30, 'Md':30, 'No':30, 'Lr':30,
         'Rf':30, 'Db':30, 'Sg':30, 'Bh':30, 'Hs':30, 'Mt':30, 'Ds':30, 'Rg':30, 'Cn':30,
'Nh':39, 'Fl':39, 'Mc':39, 'Lv':39, 'Ts':39, 'Og':39
}


ATOMCHARGE = {'H':1, 'He':2,
'Li':3, 'Be':4, 'B':5, 'C':6,  'N':7,  'O':8, 'F':9, 'Ne':10,
'Na':11, 'Mg':12, 'Al':13, 'Si':14,  'P':15,  'S':16, 'Cl':17, 'Ar':18,
'K':19, 'Ca':20,
'Sc':21, 'Ti':22, 'V':23, 'Cr':24, 'Mn':25, 'Fe':26, 'Co':27, 'Ni':28, 'Cu':29, 'Zn':30,
'Ga':31, 'Ge':32, 'As':33, 'Se':34, 'Br':35, 'Kr':36,
'Rb':37, 'Sr':38,
'Y':39,  'Zr':40, 'Nb':41, 'Mo':42, 'Tc':43, 'Ru':44, 'Rh':45, 'Pd':46, 'Ag':47, 'Cd':48,
'In':49, 'Sn':50, 'Sb':51, 'Te':52,  'I':53, 'Xe':54,
'Cs':55, 'Ba':56,
'La':57, 
'Ce':58, 'Pr':59, 'Nd':60, 'Pm':61, 'Sm':62, 'Eu':63, 'Gd':64, 'Tb':65, 'Dy':66, 'Ho':67, 'Er':68, 'Tm':69, 'Yb':70, 'Lu':71,
         'Hf':72, 'Ta':73,  'W':74, 'Re':75, 'Os':76, 'Ir':77, 'Pt':78, 'Au':79, 'Hg':80,
'Tl':81, 'Pb':82, 'Bi':83, 'Po':84, 'At':85, 'Rn':86, 
'Fr':87, 'Ra':88,
'Ac':89, 
'Th':90, 'Pa':91,  'U':92, 'Np':93, 'Pu':94, 'Am':95, 'Cm':96, 'Bk':97, 'Cf':98, 'Es':99,'Fm':100,'Md':101,'No':102,'Lr':103,
        'Rf':104,'Db':105,'Sg':106,'Bh':107,'Hs':108,'Mt':109,'Ds':110,'Rg':111,'Cn':112,
'Nh':113,'Fl':114,'Mc':115,'Lv':116,'Ts':117,'Og':118
}

# conversion factors
au2a=0.529177211
rcm_to_Eh=4.556335e-6
D2au=0.393430307
au2eV=27.2113987622

# =============================================================================================== #
# =============================================================================================== #
# =========================================== general routines ================================== #
# =============================================================================================== #
# =============================================================================================== #

# ======================================================================= #
def readfile(filename):
  try:
    f=open(filename)
    out=f.readlines()
    f.close()
  except IOError:
    print 'File %s does not exist!' % (filename)
    sys.exit(12)
  return out

# ======================================================================= #
def writefile(filename,content):
  # content can be either a string or a list of strings
  try:
    f=open(filename,'w')
    if isinstance(content,list):
      for line in content:
        f.write(line)
    elif isinstance(content,str):
      f.write(content)
    else:
      print 'Content %s cannot be written to file!' % (content)
      sys.exit(13)
    f.close()
  except IOError:
    print 'Could not write to file %s!' % (filename)
    sys.exit(14)

# ======================================================================= #
def isbinary(path):
  return (re.search(r':.* text',sp.Popen(["file", '-L', path], stdout=sp.PIPE).stdout.read())is None)

# ======================================================================= #
def eformat(f, prec, exp_digits):
    '''Formats a float f into scientific notation with prec number of decimals and exp_digits number of exponent digits.

    String looks like:
    [ -][0-9]\.[0-9]*E[+-][0-9]*

    Arguments:
    1 float: Number to format
    2 integer: Number of decimals
    3 integer: Number of exponent digits

    Returns:
    1 string: formatted number'''

    s = "% .*e"%(prec, f)
    mantissa, exp = s.split('e')
    return "%sE%+0*d"%(mantissa, exp_digits+1, int(exp))

# ======================================================================= #
def measuretime():
    '''Calculates the time difference between global variable starttime and the time of the call of measuretime.

    Prints the Runtime, if PRINT or DEBUG are enabled.

    Arguments:
    none

    Returns:
    1 float: runtime in seconds'''

    endtime=datetime.datetime.now()
    runtime=endtime-starttime
    hours=runtime.seconds/3600
    minutes=runtime.seconds/60-hours*60
    seconds=runtime.seconds%60
    print '==> Runtime:\n%i Days\t%i Hours\t%i Minutes\t%i Seconds\n\n' % (runtime.days,hours,minutes,seconds)
    total_seconds=runtime.days*24*3600+runtime.seconds+runtime.microseconds/1.e6
    return total_seconds

# ======================================================================= #
def removekey(d,key):
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
def containsstring(string,line):
    '''Takes a string (regular expression) and another string. Returns True if the first string is contained in the second string.

    Arguments:
    1 string: Look for this string
    2 string: within this string

    Returns:
    1 boolean'''

    a=re.search(string,line)
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
        if states[i]<1:
            continue
        yield i+1
    return

# ======================================================================= #
def itnmstates(states):

    for i in range(len(states)):
        if states[i]<1:
            continue
        for k in range(i+1):
            for j in range(states[i]):
                yield i+1,j+1,k-i/2.
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

    print starttime,gethostname(),os.getcwd()
    if not PRINT:
        return
    string='\n'
    string+='  '+'='*80+'\n'
    string+='||'+' '*80+'||\n'
    string+='||'+' '*26+'SHARC - GAUSSIAN - Interface'+' '*26+'||\n'
    string+='||'+' '*80+'||\n'
    string+='||'+' '*15+'Authors: Sebastian Mai and Maximilian F.S.J. Menger'+' '*14+'||\n'
    string+='||'+' '*80+'||\n'
    string+='||'+' '*(36-(len(version)+1)/2)+'Version: %s' % (version)+' '*(35-(len(version))/2)+'||\n'
    lens=len(versiondate.strftime("%d.%m.%y"))
    string+='||'+' '*(37-lens/2)+'Date: %s' % (versiondate.strftime("%d.%m.%y"))+' '*(37-(lens+1)/2)+'||\n'
    string+='||'+' '*80+'||\n'
    string+='  '+'='*80+'\n\n'
    print string
    if DEBUG:
        print changelogstring

# ======================================================================= #
def printQMin(QMin):

  if not PRINT:
    return
  print '==> QMin Job description for:\n%s' % (QMin['comment'])

  string='Mode:   '
  if 'init' in QMin:
      string+='\tINIT'
  if 'restart' in QMin:
      string+='\tRESTART'
  if 'samestep' in QMin:
      string+='\tSAMESTEP'
  if 'newstep' in QMin:
      string+='\tNEWSTEP'

  string+='\nTasks:  '
  if 'h' in QMin:
    string+='\tH'
  if 'soc' in QMin:
    string+='\tSOC'
  if 'dm' in QMin:
    string+='\tDM'
  if 'grad' in QMin:
    string+='\tGrad'
  if 'nacdr' in QMin:
    string+='\tNac(ddr)'
  if 'nacdt' in QMin:
    string+='\tNac(ddt)'
  if 'overlap' in QMin:
    string+='\tOverlap'
  if 'angular' in QMin:
    string+='\tAngular'
  if 'ion' in QMin:
    string+='\tDyson'
  if 'dmdr' in QMin:
    string+='\tDM-Grad'
  if 'socdr' in QMin:
    string+='\tSOC-Grad'
  if 'theodore' in QMin:
    string+='\tTheoDORE'
  if 'phases' in QMin:
    string+='\tPhases'
  print string

  string='States:        '
  for i in itmult(QMin['states']):
    string+='% 2i %7s  ' % (QMin['states'][i-1],IToMult[i])
  print string

  string='Charges:       '
  for i in itmult(QMin['states']):
    string+='%+2i %7s  ' % (QMin['chargemap'][i],'')
  print string

  string='Restricted:    '
  for i in itmult(QMin['states']):
    string+='%5s       ' % (QMin['jobs'][QMin['multmap'][i]]['restr'])
  print string

  string='Method: \t'
  if QMin['template']['no_tda']:
      string+='TD-'
  else:
      string+='TDA-'
  string+=QMin['template']['functional'].split()[0].upper()
  string+='/%s' % (QMin['template']['basis'])
  parts=[]
  if QMin['template']['dispersion']:
      parts.append(QMin['template']['dispersion'].split()[0].upper())
  #if QMin['template']['qmmm']:
      #parts.append('QM/MM')
  if QMin['template']['scrf']:
      parts.append('SCRF')
  if len(parts)>0:
      string+='\t('
      string+=','.join(parts)
      string+=')'
  print string

  string='Found Geo'
  if 'veloc' in QMin:
    string+=' and Veloc! '
  else:
    string+='! '
  string+='NAtom is %i.\n' % (QMin['natom'])
  print string

  string='Geometry in Bohrs (%i atoms):\n' % QMin['natom']
  if DEBUG:
    for i in range(QMin['natom']):
      string+='%2s ' % (QMin['geo'][i][0])
      for j in range(3):
        string+='% 7.4f ' % (QMin['geo'][i][j+1])
      string+='\n'
  else:
    for i in range(min(QMin['natom'],5)):
      string+='%2s ' % (QMin['geo'][i][0])
      for j in range(3):
        string+='% 7.4f ' % (QMin['geo'][i][j+1])
      string+='\n'
    if QMin['natom']>5:
      string+='..     ...     ...     ...\n'
      string+='%2s ' % (QMin['geo'][-1][0])
      for j in range(3):
        string+='% 7.4f ' % (QMin['geo'][-1][j+1])
      string+='\n'
  print string

  if 'veloc' in QMin and DEBUG:
    string=''
    for i in range(QMin['natom']):
      string+='%s ' % (QMin['geo'][i][0])
      for j in range(3):
        string+='% 7.4f ' % (QMin['veloc'][i][j])
      string+='\n'
    print string

  if 'grad' in QMin:
    string='Gradients requested:   '
    for i in range(1,QMin['nmstates']+1):
      if i in QMin['grad']:
        string+='X '
      else:
        string+='. '
    string+='\n'
    print string

  print 'State map:'
  pprint.pprint(QMin['statemap'])
  print

  for i in sorted(QMin):
    if not any( [i==j for j in ['h','dm','soc','dmdr','socdr','theodore','geo','veloc','states','comment', 'grad','nacdr','ion','overlap','template','statemap'] ] ):
      if not any( [i==j for j in ['ionlist'] ] ) or DEBUG:
        string=i+': '
        string+=str(QMin[i])
        print string
  print '\n'
  sys.stdout.flush()


# ======================================================================= #
def printcomplexmatrix(matrix,states):
    '''Prints a formatted matrix. Zero elements are not printed, blocks of different mult and MS are delimited by dashes. Also prints a matrix with the imaginary parts, of any one element has non-zero imaginary part.

    Arguments:
    1 list of list of complex: the matrix
    2 list of integers: states specs'''

    nmstates=0
    for i in range(len(states)):
        nmstates+=states[i]*(i+1)
    string='Real Part:\n'
    string+='-'*(11*nmstates+nmstates/3)
    string+='\n'
    istate=0
    for imult,i,ms in itnmstates(states):
        jstate=0
        string+='|'
        for jmult,j,ms2 in itnmstates(states):
            if matrix[istate][jstate].real==0.:
                string+=' '*11
            else:
                string+='% .3e ' % (matrix[istate][jstate].real)
            if j==states[jmult-1]:
                string+='|'
            jstate+=1
        string+='\n'
        if i==states[imult-1]:
            string+='-'*(11*nmstates+nmstates/3)
            string+='\n'
        istate+=1
    print string
    imag=False
    string='Imaginary Part:\n'
    string+='-'*(11*nmstates+nmstates/3)
    string+='\n'
    istate=0
    for imult,i,ms in itnmstates(states):
        jstate=0
        string+='|'
        for jmult,j,ms2 in itnmstates(states):
            if matrix[istate][jstate].imag==0.:
                string+=' '*11
            else:
                imag=True
                string+='% .3e ' % (matrix[istate][jstate].imag)
            if j==states[jmult-1]:
                string+='|'
            jstate+=1
        string+='\n'
        if i==states[imult-1]:
            string+='-'*(11*nmstates+nmstates/3)
            string+='\n'
        istate+=1
    string+='\n'
    if imag:
        print string

# ======================================================================= #
def printgrad(grad,natom,geo):
    '''Prints a gradient or nac vector. Also prints the atom elements. If the gradient is identical zero, just prints one line.

    Arguments:
    1 list of list of float: gradient
    2 integer: natom
    3 list of list: geometry specs'''

    string=''
    iszero=True
    for atom in range(natom):
        if not DEBUG:
            if atom==5:
                string+='...\t...\t     ...\t     ...\t     ...\n'
            if 5<=atom<natom-1:
                continue
        string+='%i\t%s\t' % (atom+1,geo[atom][0])
        for xyz in range(3):
            if grad[atom][xyz]!=0:
                iszero=False
            g=grad[atom][xyz]
            if isinstance(g,float):
                string+='% .5f\t' % (g)
            elif isinstance(g,complex):
                string+='% .5f\t% .5f\t\t' % (g.real,g.imag)
        string+='\n'
    if iszero:
        print '\t\t...is identical zero...\n'
    else:
        print string

# ======================================================================= #
def printtheodore(matrix,QMin):
    string='%6s ' % 'State'
    for i in QMin['template']['theodore_prop']:
        string+='%6s ' % i
    for i in range(len(QMin['template']['theodore_fragment'])):
        for j in range(len(QMin['template']['theodore_fragment'])):
            string+='  Om%1i%1i ' % (i+1,j+1)
    string+='\n'+'-------'*(1+QMin['template']['theodore_n'])+'\n'
    istate=0
    for imult,i,ms in itnmstates(QMin['states']):
        istate+=1
        string+='%6i ' % istate
        for i in matrix[istate-1]:
            string+='%6.4f ' % i.real
        string+='\n'
    print string

# ======================================================================= #
def printQMout(QMin,QMout):
    '''If PRINT, prints a summary of all requested QM output values. Matrices are formatted using printcomplexmatrix, vectors using printgrad. 

    Arguments:
    1 dictionary: QMin
    2 dictionary: QMout'''

    #if DEBUG:
        #pprint.pprint(QMout)
    if not PRINT:
        return
    states=QMin['states']
    nstates=QMin['nstates']
    nmstates=QMin['nmstates']
    natom=QMin['natom']
    print '\n\n>>>>>>>>>>>>> Results\n'
    # Hamiltonian matrix, real or complex
    if 'h' in QMin or 'soc' in QMin:
        eshift=math.ceil(QMout['h'][0][0].real)
        print '=> Hamiltonian Matrix:\nDiagonal Shift: %9.2f' % (eshift)
        matrix=deepcopy(QMout['h'])
        for i in range(nmstates):
            matrix[i][i]-=eshift
        printcomplexmatrix(matrix,states)
    # Dipole moment matrices
    if 'dm' in QMin:
        print '=> Dipole Moment Matrices:\n'
        for xyz in range(3):
            print 'Polarisation %s:' % (IToPol[xyz])
            matrix=QMout['dm'][xyz]
            printcomplexmatrix(matrix,states)
    # Gradients
    if 'grad' in QMin:
        print '=> Gradient Vectors:\n'
        istate=0
        for imult,i,ms in itnmstates(states):
            print '%s\t%i\tMs= % .1f:' % (IToMult[imult],i,ms)
            printgrad(QMout['grad'][istate],natom,QMin['geo'])
            istate+=1
    # Overlaps
    if 'overlap' in QMin:
        print '=> Overlap matrix:\n'
        matrix=QMout['overlap']
        printcomplexmatrix(matrix,states)
        if 'phases' in QMout:
            print '=> Wavefunction Phases:\n'
            for i in range(nmstates):
                print '% 3.1f % 3.1f' % (QMout['phases'][i].real,QMout['phases'][i].imag)
            print '\n'
    # Spin-orbit coupling derivatives
    if 'socdr' in QMin:
        print '=> Spin-Orbit Gradient Vectors:\n'
        istate=0
        for imult,i,ims in itnmstates(states):
            jstate=0
            for jmult,j,jms in itnmstates(states):
                print '%s\t%i\tMs= % .1f -- %s\t%i\tMs= % .1f:' % (IToMult[imult],i,ims,IToMult[jmult],j,jms)
                printgrad(QMout['socdr'][istate][jstate],natom,QMin['geo'])
                jstate+=1
            istate+=1
    # Dipole moment derivatives
    if 'dmdr' in QMin:
        print '=> Dipole moment derivative vectors:\n'
        istate=0
        for imult,i,msi in itnmstates(states):
            jstate=0
            for jmult,j,msj in itnmstates(states):
                if imult==jmult and msi==msj:
                    for ipol in range(3):
                        print '%s\tStates %i - %i\tMs= % .1f\tPolarization %s:' % (IToMult[imult],i,j,msi,IToPol[ipol])
                        printgrad(QMout['dmdr'][ipol][istate][jstate],natom,QMin['geo'])
                jstate+=1
            istate+=1
  # Property matrix (dyson norms)
    if 'ion' in QMin and 'prop' in QMout:
        print '=> Property matrix:\n'
        matrix=QMout['prop']
        printcomplexmatrix(matrix,states)
  # TheoDORE
    if 'theodore' in QMin:
        print '=> TheoDORE results:\n'
        matrix=QMout['theodore']
        printtheodore(matrix,QMin)
    sys.stdout.flush()


# =============================================================================================== #
# =============================================================================================== #
# ======================================= Matrix initialization ================================= #
# =============================================================================================== #
# =============================================================================================== #

# ======================================================================= #         OK
def makecmatrix(a,b):
    '''Initialises a complex axb matrix.

    Arguments:
    1 integer: first dimension
    2 integer: second dimension

    Returns;
    1 list of list of complex'''

    mat=[ [ complex(0.,0.) for i in range(a) ] for j in range(b) ]
    return mat

# ======================================================================= #         OK
def makermatrix(a,b):
    '''Initialises a real axb matrix.

    Arguments:
    1 integer: first dimension
    2 integer: second dimension

    Returns;
    1 list of list of real'''

    mat=[ [ 0. for i in range(a) ] for j in range(b) ]
    return mat


# =============================================================================================== #
# =============================================================================================== #
# =========================================== QMout writing ===================================== #
# =============================================================================================== #
# =============================================================================================== #


# ======================================================================= #
def writeQMout(QMin,QMout,QMinfilename):
    '''Writes the requested quantities to the file which SHARC reads in. The filename is QMinfilename with everything after the first dot replaced by "out". 

    Arguments:
    1 dictionary: QMin
    2 dictionary: QMout
    3 string: QMinfilename'''

    k=QMinfilename.find('.')
    if k==-1:
        outfilename=QMinfilename+'.out'
    else:
        outfilename=QMinfilename[:k]+'.out'
    if PRINT:
        print '===> Writing output to file %s in SHARC Format\n' % (outfilename)
    string=''
    if 'h' in QMin or 'soc' in QMin:
        string+=writeQMoutsoc(QMin,QMout)
    if 'dm' in QMin:
        string+=writeQMoutdm(QMin,QMout)
    if 'grad' in QMin:
        string+=writeQMoutgrad(QMin,QMout)
    if 'overlap' in QMin:
        string+=writeQMoutnacsmat(QMin,QMout)
    if 'socdr' in QMin:
        string+=writeQMoutsocdr(QMin,QMout)
    if 'dmdr' in QMin:
        string+=writeQMoutdmdr(QMin,QMout)
    if 'ion' in QMin:
        string+=writeQMoutprop(QMin,QMout)
    if 'theodore' in QMin: # or QMin['template']['qmmm']:
        string+=writeQMoutTHEODORE(QMin,QMout)
    if 'phases' in QMin:
        string+=writeQmoutPhases(QMin,QMout)
    string+=writeQMouttime(QMin,QMout)
    outfile=os.path.join(QMin['pwd'],outfilename)
    writefile(outfile,string)
    return

# ======================================================================= #
def writeQMoutsoc(QMin,QMout):
    '''Generates a string with the Spin-Orbit Hamiltonian in SHARC format.

    The string starts with a ! followed by a flag specifying the type of data. In the next line, the dimensions of the matrix are given, followed by nmstates blocks of nmstates elements. Blocks are separated by a blank line.

    Arguments:
    1 dictionary: QMin
    2 dictionary: QMout

    Returns:
    1 string: multiline string with the SOC matrix'''

    states=QMin['states']
    nstates=QMin['nstates']
    nmstates=QMin['nmstates']
    natom=QMin['natom']
    string=''
    string+='! %i Hamiltonian Matrix (%ix%i, complex)\n' % (1,nmstates,nmstates)
    string+='%i %i\n' % (nmstates,nmstates)
    for i in range(nmstates):
        for j in range(nmstates):
            string+='%s %s ' % (eformat(QMout['h'][i][j].real,9,3),eformat(QMout['h'][i][j].imag,9,3))
        string+='\n'
    string+='\n'
    return string

# ======================================================================= #
def writeQMoutdm(QMin,QMout):
    '''Generates a string with the Dipole moment matrices in SHARC format.

    The string starts with a ! followed by a flag specifying the type of data. In the next line, the dimensions of the matrix are given, followed by nmstates blocks of nmstates elements. Blocks are separated by a blank line. The string contains three such matrices.

    Arguments:
    1 dictionary: QMin
    2 dictionary: QMout

    Returns:
    1 string: multiline string with the DM matrices'''

    states=QMin['states']
    nstates=QMin['nstates']
    nmstates=QMin['nmstates']
    natom=QMin['natom']
    string=''
    string+='! %i Dipole Moment Matrices (3x%ix%i, complex)\n' % (2,nmstates,nmstates)
    for xyz in range(3):
        string+='%i %i\n' % (nmstates,nmstates)
        for i in range(nmstates):
            for j in range(nmstates):
                string+='%s %s ' % (eformat(QMout['dm'][xyz][i][j].real,9,3),eformat(QMout['dm'][xyz][i][j].imag,9,3))
            string+='\n'
        #string+='\n'
    string+='\n'
    return string

# ======================================================================= #
def writeQMoutgrad(QMin,QMout):
    '''Generates a string with the Gradient vectors in SHARC format.

    The string starts with a ! followed by a flag specifying the type of data. On the next line, natom and 3 are written, followed by the gradient, with one line per atom and a blank line at the end. Each MS component shows up (nmstates gradients are written).

    Arguments:
    1 dictionary: QMin
    2 dictionary: QMout

    Returns:
    1 string: multiline string with the Gradient vectors'''

    states=QMin['states']
    nstates=QMin['nstates']
    nmstates=QMin['nmstates']
    natom=QMin['natom']
    string=''
    string+='! %i Gradient Vectors (%ix%ix3, real)\n' % (3,nmstates,natom)
    i=0
    for imult,istate,ims in itnmstates(states):
        string+='%i %i ! %i %i %i\n' % (natom,3,imult,istate,ims)
        for atom in range(natom):
            for xyz in range(3):
                string+='%s ' % (eformat(QMout['grad'][i][atom][xyz],9,3))
            string+='\n'
        #string+='\n'
        i+=1
    string+='\n'
    return string

# ======================================================================= #
def writeQMoutnacsmat(QMin,QMout):
    '''Generates a string with the adiabatic-diabatic transformation matrix in SHARC format.

    The string starts with a ! followed by a flag specifying the type of data. In the next line, the dimensions of the matrix are given, followed by nmstates blocks of nmstates elements. Blocks are separated by a blank line.

    Arguments:
    1 dictionary: QMin
    2 dictionary: QMout

    Returns:
    1 string: multiline string with the transformation matrix'''

    states=QMin['states']
    nstates=QMin['nstates']
    nmstates=QMin['nmstates']
    natom=QMin['natom']
    string=''
    string+='! %i Overlap matrix (%ix%i, complex)\n' % (6,nmstates,nmstates)
    string+='%i %i\n' % (nmstates,nmstates)
    for j in range(nmstates):
        for i in range(nmstates):
            string+='%s %s ' % (eformat(QMout['overlap'][j][i].real,9,3),eformat(QMout['overlap'][j][i].imag,9,3))
        string+='\n'
    string+='\n'
    return string

# ======================================================================= #
def writeQMoutdmdr(QMin,QMout):

  states=QMin['states']
  nmstates=QMin['nmstates']
  natom=QMin['natom']
  string=''
  string+='! %i Dipole moment derivatives (%ix%ix3x%ix3, real)\n' % (12,nmstates,nmstates,natom)
  i=0
  for imult,istate,ims in itnmstates(states):
    j=0
    for jmult,jstate,jms in itnmstates(states):
      for ipol in range(3):
        string+='%i %i ! m1 %i s1 %i ms1 %i   m2 %i s2 %i ms2 %i   pol %i\n' % (natom,3,imult,istate,ims,jmult,jstate,jms,ipol)
        for atom in range(natom):
          for xyz in range(3):
            string+='%s ' % (eformat(QMout['dmdr'][ipol][i][j][atom][xyz],12,3))
          string+='\n'
        string+=''
      j+=1
    i+=1
  string+='\n'
  return string

# ======================================================================= #
def writeQMoutsocdr(QMin,QMout):

  states=QMin['states']
  nmstates=QMin['nmstates']
  natom=QMin['natom']
  string=''
  string+='! %i Spin-Orbit coupling derivatives (%ix%ix3x%ix3, complex)\n' % (13,nmstates,nmstates,natom)
  i=0
  for imult,istate,ims in itnmstates(states):
    j=0
    for jmult,jstate,jms in itnmstates(states):
        string+='%i %i ! m1 %i s1 %i ms1 %i   m2 %i s2 %i ms2 %i\n' % (natom,3,imult,istate,ims,jmult,jstate,jms)
        for atom in range(natom):
            for xyz in range(3):
                string+='%s %s ' % (eformat(QMout['socdr'][i][j][atom][xyz].real,12,3),eformat(QMout['socdr'][i][j][atom][xyz].imag,12,3))
        string+='\n'
        string+=''
        j+=1
    i+=1
  string+='\n'
  return string

# ======================================================================= #
def writeQMoutprop(QMin,QMout):

    nmstates=QMin['nmstates']

    # print property matrix (flag 11) for backwards compatibility
    string=''
    string+='! %i Property Matrix (%ix%i, complex)\n' % (11,nmstates,nmstates)
    string+='%i %i\n' % (nmstates,nmstates)
    for i in range(nmstates):
        for j in range(nmstates):
            string+='%s %s ' % (eformat(QMout['prop'][i][j].real,12,3),eformat(QMout['prop'][i][j].imag,12,3))
        string+='\n'
    string+='\n'

    # print property matrices (flag 20) in new format
    string+='! %i Property Matrices\n' % (20)
    string+='%i    ! number of property matrices\n' % (1)

    string+='! Property Matrix Labels (%i strings)\n' % (1)
    string+='Dyson norms\n'

    string+='! Property Matrices (%ix%ix%i, complex)\n' % (1,nmstates,nmstates)
    string+='%i %i   ! Dyson norms\n' % (nmstates,nmstates)
    for i in range(nmstates):
        for j in range(nmstates):
            string+='%s %s ' % (eformat(QMout['prop'][i][j].real,12,3),eformat(QMout['prop'][i][j].imag,12,3))
        string+='\n'
    string+='\n'

    return string

# ======================================================================= #
def writeQMoutTHEODORE(QMin,QMout):

    nmstates=QMin['nmstates']
    nprop=QMin['template']['theodore_n']
    #if QMin['template']['qmmm']:
        #nprop+=7
    if nprop<=0:
        return '\n'

    string=''

    string+='! %i Property Vectors\n' % (21)
    string+='%i    ! number of property vectors\n' % (nprop)

    string+='! Property Vector Labels (%i strings)\n' % (nprop)
    descriptors=[]
    if 'theodore' in QMin:
        for i in QMin['template']['theodore_prop']:
            descriptors.append('%s' % i)
            string+=descriptors[-1]+'\n'
        for i in range(len(QMin['template']['theodore_fragment'])):
            for j in range(len(QMin['template']['theodore_fragment'])):
                descriptors.append('Om_{%i,%i}' % (i+1,j+1))
                string+=descriptors[-1]+'\n'
    #if QMin['template']['qmmm']:
        #for label in QMout['qmmm_energies']:
            #descriptors.append(label)
            #string+=label+'\n'

    string+='! Property Vectors (%ix%i, real)\n' % (nprop,nmstates)
    if 'theodore' in QMin:
        for i in range(QMin['template']['theodore_n']):
            string+='! TheoDORE descriptor %i (%s)\n' % (i+1,descriptors[i])
            for j in range(nmstates):
                string+='%s\n' % (eformat(QMout['theodore'][j][i].real,12,3))
    #if QMin['template']['qmmm']:
        #for label in QMout['qmmm_energies']:
            #string+='! QM/MM energy contribution (%s)\n' % (label)
            #for j in range(nmstates):
                #string+='%s\n' % (eformat(QMout['qmmm_energies'][label],12,3))
    string+='\n'

    return string

# ======================================================================= #
def writeQmoutPhases(QMin,QMout):

    string='! 7 Phases\n%i ! for all nmstates\n' % (QMin['nmstates'])
    for i in range(QMin['nmstates']):
        string+='%s %s\n' % (eformat(QMout['phases'][i].real,9,3),eformat(QMout['phases'][i].imag,9,3))
    return string

# ======================================================================= #
def writeQMouttime(QMin,QMout):
    '''Generates a string with the quantum mechanics total runtime in SHARC format.

    The string starts with a ! followed by a flag specifying the type of data. In the next line, the runtime is given

    Arguments:
    1 dictionary: QMin
    2 dictionary: QMout

    Returns:
    1 string: multiline string with the runtime'''

    string='! 8 Runtime\n%s\n' % (eformat(QMout['runtime'],9,3))
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

    exist=os.path.exists(SCRATCHDIR)
    if exist:
        isfile=os.path.isfile(SCRATCHDIR)
        if isfile:
            print '$SCRATCHDIR=%s exists and is a file!' % (SCRATCHDIR)
            sys.exit(15)
    else:
        try:
            os.makedirs(SCRATCHDIR)
        except OSError:
            print 'Can not create SCRATCHDIR=%s\n' % (SCRATCHDIR)
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
def getGaussianVersion(groot):
  tries={'g09': '09',
         'g16': '16'}
  ls=os.listdir(groot)
  for i in tries:
    if i in ls:
      return tries[i]
  else:
    print 'Found no executable (possible names: %s) in $groot!' % (list(tries))
    sys.exit(17)

# ======================================================================= #
def getsh2Gaukey(sh2Gau,key):
  i=-1
  while True:
    i+=1
    try:
      line=re.sub('#.*$','',sh2Gau[i])
    except IndexError:
      break
    line=line.split(None,1)
    if line==[]:
      continue
    if key.lower() in line[0].lower():
      return line
  return ['','']

# ======================================================================= #
def get_sh2Gau_environ(sh2Gau,key,environ=True,crucial=True):
  line=getsh2Gaukey(sh2Gau,key)
  if line[0]:
    LINE=line[1]
    LINE=removequotes(LINE).strip()
  else:
    if environ:
      LINE=os.getenv(key.upper())
      if not LINE:
        if crucial:
          print 'Either set $%s or give path to %s in GAUSSIAN.resources!' % (key.upper(),key.upper())
          sys.exit(18)
        else:
          return None
    else:
      if crucial:
        print 'Give path to %s in GAUSSIAN.resources!' % (key.upper())
        sys.exit(19)
      else:
        return None
  LINE=os.path.expandvars(LINE)
  LINE=os.path.expanduser(LINE)
  if containsstring(';',LINE):
    print "$%s contains a semicolon. Do you probably want to execute another command after %s? I can't do that for you..." % (key.upper(),key.upper())
    sys.exit(20)
  return LINE

# ======================================================================= #
def get_pairs(QMinlines,i):
  nacpairs=[]
  while True:
    i+=1
    try:
      line=QMinlines[i].lower()
    except IndexError:
      print '"keyword select" has to be completed with an "end" on another line!'
      sys.exit(21)
    if 'end' in line:
      break
    fields=line.split()
    try:
      nacpairs.append([int(fields[0]),int(fields[1])])
    except ValueError:
      print '"nacdr select" is followed by pairs of state indices, each pair on a new line!'
      sys.exit(22)
  return nacpairs,i

# ======================================================================= #         OK
def readQMin(QMinfilename):
    '''Reads the time-step dependent information from QMinfilename. 

    Arguments:
    1 string: name of the QMin file

    Returns:
    1 dictionary: QMin'''


# --------------------------------------------- QM.in ----------------------------------

    QMinlines=readfile(QMinfilename)
    QMin={}

    # Get natom
    try:
        natom=int(QMinlines[0])
    except ValueError:
        print 'first line must contain the number of atoms!'
        sys.exit(23)
    QMin['natom']=natom
    if len(QMinlines)<natom+4:
        print 'Input file must contain at least:\nnatom\ncomment\ngeometry\nkeyword "states"\nat least one task'
        sys.exit(24)

    # Save Comment line
    QMin['comment']=QMinlines[1]

    # Get geometry and possibly velocity (for backup-analytical non-adiabatic couplings)
    QMin['geo']=[]
    QMin['veloc']=[]
    hasveloc=True
    QMin['frozcore']=0
    QMin['Atomcharge']=0
    for i in range(2,natom+2):
        if not containsstring('[a-zA-Z][a-zA-Z]?[0-9]*.*[-]?[0-9]+[.][0-9]*.*[-]?[0-9]+[.][0-9]*.*[-]?[0-9]+[.][0-9]*', QMinlines[i]):
            print 'Input file does not comply to xyz file format! Maybe natom is just wrong.'
            sys.exit(25)
        fields=QMinlines[i].split()
        fields[0]=fields[0].title()
        symb = fields[0]
        QMin['frozcore']+=FROZENS[symb]
        QMin['Atomcharge']+=ATOMCHARGE[symb]
        for j in range(1,4):
            fields[j]=float(fields[j])
        QMin['geo'].append(fields[0:4])
        if len(fields)>=7:
            for j in range(4,7):
                fields[j]=float(fields[j])
            QMin['veloc'].append(fields[4:7])
        else:
            hasveloc=False
    if not hasveloc:
        QMin=removekey(QMin,'veloc')


    # Parse remaining file
    i=natom+1
    while i+1<len(QMinlines):
        i+=1
        line=QMinlines[i]
        line=re.sub('#.*$','',line)
        if len(line.split())==0:
            continue
        key=line.lower().split()[0]
        if 'savedir' in key:
            args=line.split()[1:]
        else:
            args=line.lower().split()[1:]
        if key in QMin:
            print 'Repeated keyword %s in line %i in input file! Check your input!' % (key,i+1)
            continue  # only first instance of key in QM.in takes effect
        if len(args)>=1 and args[0]=='select':
            pairs,i=get_pairs(QMinlines,i)
            QMin[key]=pairs
        else:
            QMin[key]=args

    if 'unit' in QMin:
        if QMin['unit'][0]=='angstrom':
            factor=1./au2a
        elif QMin['unit'][0]=='bohr':
            factor=1.
        else:
            print 'Dont know input unit %s!' % (QMin['unit'][0])
            sys.exit(26)
    else:
        factor=1./au2a

    for iatom in range(len(QMin['geo'])):
        for ixyz in range(3):
            QMin['geo'][iatom][ixyz+1]*=factor


    if not 'states' in QMin:
        print 'Keyword "states" not given!'
        sys.exit(27)
    # Calculate states, nstates, nmstates
    for i in range(len(QMin['states'])):
        QMin['states'][i]=int(QMin['states'][i])
    reduc=0
    for i in reversed(QMin['states']):
        if i==0:
            reduc+=1
        else:
            break
    for i in range(reduc):
        del QMin['states'][-1]
    nstates=0
    nmstates=0
    for i in range(len(QMin['states'])):
        nstates+=QMin['states'][i]
        nmstates+=QMin['states'][i]*(i+1)
    QMin['nstates']=nstates
    QMin['nmstates']=nmstates


    # Various logical checks
    if not 'states' in QMin:
        print 'Number of states not given in QM input file %s!' % (QMinfilename)
        sys.exit(28)

    ## only singlets TODO
    #if len(QMin['states'])>1 and sum(QMin['states'][1:])>0:
      #print 'Currently, only singlet states are allowed!'
      #sys.exit(29)

    possibletasks=['h','soc','dm','grad','overlap','dmdr','socdr','ion','theodore','phases']
    if not any([i in QMin for i in possibletasks]):
        print 'No tasks found! Tasks are %s.' % possibletasks
        sys.exit(30)

    if not 'h' in QMin and not 'soc' in QMin:
        QMin['h']=[]

    if  'soc' in QMin:
        print 'Spin-orbit couplings cannot be computed currently!'
        sys.exit(31)

    #if 'soc' in QMin and (len(QMin['states'])<3 or QMin['states'][2]<=0):
        #QMin=removekey(QMin,'soc')
        #QMin['h']=[]
        #print 'HINT: No triplet states requested, turning off SOC request.'

    if 'samestep' in QMin and 'init' in QMin:
        print '"Init" and "Samestep" cannot be both present in QM.in!'
        sys.exit(32)

    if 'restart' in QMin and 'init' in QMin:
        print '"Init" and "Samestep" cannot be both present in QM.in!'
        sys.exit(33)

    if 'phases' in QMin:
        QMin['overlap']=[]

    if 'overlap' in QMin and 'init' in QMin:
        print '"overlap" and "phases" cannot be calculated in the first timestep! Delete either "overlap" or "init"'
        sys.exit(34)

    if not 'init' in QMin and not 'samestep' in QMin and not 'restart'in QMin:
        QMin['newstep']=[]

    if not any([i in QMin for i in ['h','soc','dm','grad']]) and 'overlap' in QMin:
        QMin['h']=[]

    if 'nacdt' in QMin or 'nacdr' in QMin:
        print 'Within the SHARC-GAUSSIAN interface couplings can only be calculated via the overlap method. "nacdr" and "nacdt" are not supported.'
        sys.exit(35)

    if 'dmdr' in QMin:
        print 'Dipole derivatives ("dmdr") not currently supported'
        sys.exit(36)

    if 'socdr' in QMin:
        print 'Spin-orbit coupling derivatives ("socdr") are not implemented'
        sys.exit(37)


    # Check for correct gradient list
    if 'grad' in QMin:
        if len(QMin['grad'])==0 or QMin['grad'][0]=='all':
            QMin['grad']=[ i+1 for i in range(nmstates)]
            #pass
        else:
            for i in range(len(QMin['grad'])):
                try:
                    QMin['grad'][i]=int(QMin['grad'][i])
                except ValueError:
                    print 'Arguments to keyword "grad" must be "all" or a list of integers!'
                    sys.exit(38)
                if QMin['grad'][i]>nmstates:
                    print 'State for requested gradient does not correspond to any state in QM input file state list!'
                    sys.exit(39)


    ## shorten states array TODO
    #QMin['states']=[QMin['states'][0]]



# --------------------------------------------- GAUSSIAN.resources ----------------------------------

    QMin['pwd']=os.getcwd()

    # open GAUSSIAN.resources
    filename='GAUSSIAN.resources'
    if os.path.isfile(filename):
        sh2Gau=readfile(filename)
    else:
        print 'HINT: reading resources from SH2Gau.inp'
        sh2Gau=readfile('SH2Gau.inp')


    # Set up scratchdir
    line=get_sh2Gau_environ(sh2Gau,'scratchdir',False,False)
    if line==None:
        line=QMin['pwd']+'/SCRATCHDIR/'
    line=os.path.expandvars(line)
    line=os.path.expanduser(line)
    line=os.path.abspath(line)
    #checkscratch(line)
    QMin['scratchdir']=line
    link(QMin['scratchdir'], os.path.join(QMin['pwd'],'SCRATCH'), False,False)


    # Set up savedir
    if 'savedir' in QMin:
        # savedir may be read from QM.in file
        line=QMin['savedir'][0]
    else:
        line=get_sh2Gau_environ(sh2Gau,'savedir',False,False)
        if line==None:
            line=QMin['pwd']+'/SAVEDIR/'
    line=os.path.expandvars(line)
    line=os.path.expanduser(line)
    line=os.path.abspath(line)
    if 'init' in QMin:
        checkscratch(line)
    QMin['savedir']=line
    link(QMin['savedir'],os.path.join(QMin['pwd'],'SAVE'),False,False)


    # setup environment for Gaussian
    QMin['groot']=get_sh2Gau_environ(sh2Gau,'groot')
    QMin['Gversion']=getGaussianVersion(QMin['groot'])
    print 'Detected GAUSSIAN version %s' % QMin['Gversion']
    os.environ['g%sroot' % QMin['Gversion']]=QMin['groot']
    os.environ['GAUSS_EXEDIR']=QMin['groot']
    os.environ['GAUSS_SCRDIR']='.'
    os.environ['PATH']='$GAUSS_EXEDIR:'+os.environ['PATH']
    QMin['GAUSS_EXEDIR']=QMin['groot']
    QMin['GAUSS_EXE']=os.path.join(QMin['groot'],'g%s'%QMin['Gversion'])


    # debug option
    line=getsh2Gaukey(sh2Gau,'debug')
    if line[0]:
        if len(line)<=1 or 'true' in line[1].lower():
            global DEBUG
            DEBUG=True


    # debug option
    line=getsh2Gaukey(sh2Gau,'no_print')
    if line[0]:
        if len(line)<=1 or 'true' in line[1].lower():
            global PRINT
            PRINT=False



    # resources
    QMin['ncpu']=1
    line=getsh2Gaukey(sh2Gau,'ncpu')
    if line[0]:
        try:
            QMin['ncpu']=int(line[1])
        except ValueError:
            print 'Number of CPUs does not evaluate to numerical value!'
            sys.exit(40)
    if os.environ.get('NSLOTS') != None:
        QMin['ncpu']=int(os.environ.get('NSLOTS'))
        print 'Detected $NSLOTS variable. Will use ncpu=%i' % (QMin['ncpu'])
    elif os.environ.get('SLURM_NTASKS_PER_NODE') != None:
        QMin['ncpu']=int(os.environ.get('SLURM_NTASKS_PER_NODE'))
        print 'Detected $SLURM_NTASKS_PER_NODE variable. Will use ncpu=%i' % (QMin['ncpu'])
    QMin['ncpu']=max(1,QMin['ncpu'])

    QMin['delay']=0.0
    line=getsh2Gaukey(sh2Gau,'delay')
    if line[0]:
        try:
            QMin['delay']=float(line[1])
        except ValueError:
            print 'Submit delay does not evaluate to numerical value!'
            sys.exit(41)

    QMin['schedule_scaling']=0.9
    line=getsh2Gaukey(sh2Gau,'schedule_scaling')
    if line[0]:
        try:
            x=float(line[1])
            if 0<x<=1.:
                QMin['schedule_scaling']=x
        except ValueError:
            print '"schedule_scaling" does not evaluate to numerical value!'
            sys.exit(42)


    # initial MO guess settings
    # if neither keyword is present, the interface will reuse MOs from savedir, or let ADF generate a guess
    line=getsh2Gaukey(sh2Gau,'always_orb_init')
    if line[0]:
        QMin['always_orb_init']=[]
    line=getsh2Gaukey(sh2Gau,'always_guess')
    if line[0]:
        QMin['always_guess']=[]
    if 'always_orb_init' in QMin and 'always_guess' in QMin:
        print 'Keywords "always_orb_init" and "always_guess" cannot be used together!'
        sys.exit(43)


    # wfoverlap settings
    if 'overlap' in QMin or 'ion' in QMin:
        QMin['wfoverlap']=get_sh2Gau_environ(sh2Gau,'wfoverlap',False,False)
        if QMin['wfoverlap']==None:
            ciopath=os.path.join(os.path.expandvars(os.path.expanduser('$SHARC')),'wfoverlap.x')
            if os.path.isfile(ciopath):
                QMin['wfoverlap']=ciopath
            else:
                print 'Give path to wfoverlap.x in GAUSSIAN.resources!'
                sys.exit(44)

    # memory
    QMin['memory']=100
    line=getsh2Gaukey(sh2Gau,'memory')
    if line[0]:
        QMin['memory']= float(line[1])

    # truncation threshold
    QMin['wfthres']=0.99 
    line=getsh2Gaukey(sh2Gau,'wfthres')
    if line[0]:
        QMin['wfthres']= float(line[1])

    # get the nooverlap keyword: no dets will be extracted if present
    line=getsh2Gaukey(sh2Gau,'nooverlap')
    if line[0]:
        QMin['nooverlap']=[]


    # TheoDORE settings
    if 'theodore' in QMin:
        QMin['theodir']=get_sh2Gau_environ(sh2Gau,'theodir',False,False)
        if QMin['theodir']==None or not os.path.isdir(QMin['theodir']):
            print 'Give path to the TheoDORE installation directory in GAUSSIAN.resources!'
            sys.exit(45)
        os.environ['THEODIR']=QMin['theodir']
        if 'PYTHONPATH' in os.environ:
          os.environ['PYTHONPATH']+=os.pathsep + os.path.join(QMin['theodir'],'lib') + os.pathsep + QMin['theodir']
        else:
          os.environ['PYTHONPATH']=os.path.join(QMin['theodir'],'lib') + os.pathsep + QMin['theodir']


    # neglected gradients
    QMin['neglected_gradient']='zero'
    if 'grad' in QMin:
        line=getsh2Gaukey(sh2Gau,'neglected_gradient')
        if line[0]:
            if line[1].lower().strip()=='zero':
                QMin['neglected_gradient']='zero'
            elif line[1].lower().strip()=='gs':
                QMin['neglected_gradient']='gs'
            elif line[1].lower().strip()=='closest':
                QMin['neglected_gradient']='closest'
            else:
                print 'Unknown argument to "neglected_gradient"!'
                sys.exit(46)



# --------------------------------------------- GAUSSIAN.template ----------------------------------

    # define classes and defaults
    bools   ={'denfit'                  :False,
              'no_tda'                  :False,
              'unrestricted_triplets'   :False
              }
    strings ={'basis'                   :'6-31G',
              'functional'              :'PBEPBE',
              'dispersion'              :'',
              'grid'                    :'finegrid',
              'scrf'                    :'',
              'scf'                     :'',
              'qmmm_table'              :'GAUSSIAN.qmmm.table',
              'qmmm_ff_file'            :'GAUSSIAN.ff',
              'iop'                     :'',
              'keys'                    :'',
              'basis_external'          :''
              }
    integers={
              }
    floats  ={
              }
    special ={'paddingstates'           :[0 for i in QMin['states']],
              'charge'                  :[i%2 for i in range(len(QMin['states']))],
              'theodore_prop'           :['Om','PRNTO','S_HE','Z_HE','RMSeh'],
              'theodore_fragment'       :[]
              }

    # create QMin subdictionary
    QMin['template']={}
    for i in bools:
        QMin['template'][i]=bools[i]
    for i in strings:
        QMin['template'][i]=strings[i]
    for i in integers:
        QMin['template'][i]=integers[i]
    for i in floats:
        QMin['template'][i]=floats[i]
    for i in special:
        QMin['template'][i]=special[i]

    # open template
    template=readfile('GAUSSIAN.template')

    # go through template
    for line in template:
        orig=re.sub('#.*$','',line).strip()
        line=orig.lower().split()
        if len(line)==0:
            continue
        elif line[0] in bools:
            QMin['template'][line[0]]=True
        elif line[0] in strings:
            QMin['template'][line[0]]=orig.split(None,1)[1]
        elif line[0] in integers:
            QMin['template'][line[0]]=int(float(line[1]))
        elif line[0] in floats:
            QMin['template'][line[0]]=float(line[1])
        elif line[0] in special:

            # paddingstates needs to be autoexpanded and checked
            if line[0]=='paddingstates':
                if len(line)==2:
                    QMin['template']['paddingstates']=[int(line[1])   for i in range(len(QMin['states']))]
                elif len(line)-1>=len(QMin['states']):
                    QMin['template']['paddingstates']=[int(line[1+i]) for i in range(len(QMin['states']))]
                else:
                    print 'Length of "paddingstates" does not match length of "states"!'
                    sys.exit(47)
                for i in range(len(QMin['template']['paddingstates'])):
                    if QMin['template']['paddingstates'][i]<0:
                        QMin['template']['paddingstates'][i]=0

            # charge needs to be autoexpanded, checked and assigned to the multiplicities
            elif line[0]=='charge':
                if len(line)==2:
                    charge=int(float(line[1]))
                    if (QMin['Atomcharge']+charge)%2==1 and len(QMin['states'])>1:
                        print 'HINT: Charge shifted by -1 to be compatible with multiplicities.'
                        charge-=1
                    QMin['template']['charge']=[i%2+charge for i in range(len(QMin['states']))]
                    print 'HINT: total charge per multiplicity automatically assigned, please check (%s).' % QMin['template']['charge']
                    print 'You can set the charge in the template manually for each multiplicity ("charge 0 +1 0 ...")\n'
                elif len(line)-1>=len(QMin['states']):
                    QMin['template']['charge']=[int(float(line[1+i])) for i in range(len(QMin['states']))]
                    compatible=True
                    for imult,cha in enumerate(QMin['template']['charge']):
                        if not (QMin['Atomcharge']+cha+imult)%2==0:
                            compatible=False
                    if not compatible:
                        print 'Charges from template not compatible with multiplicities!'
                        sys.exit(48)
                else:
                    print 'Length of "charge" does not match length of "states"!'
                    sys.exit(49)


    # go through sh2Gau for the theodore settings and QM/MM file names
    for line in sh2Gau:
        orig=re.sub('#.*$','',line).strip()
        line=orig.lower().split()
        if len(line)==0:
            continue
        elif line[0] in special:

            # TheoDORE properties need to be parsed in a special way
            if line[0]=='theodore_prop':
                if '[' in orig:
                    string=orig.split(None,1)[1]
                    QMin['template']['theodore_prop']=ast.literal_eval(string)
                else:
                    QMin['template']['theodore_prop']=[]
                    s=orig.split(None)[1:]
                    for i in s:
                        QMin['template']['theodore_prop'].append(i)
                theodore_spelling=['Om', 
                                   'PRNTO', 
                                   'Z_HE', 'S_HE', 'RMSeh',
                                   'POSi', 'POSf', 'POS', 
                                   'PRi', 'PRf', 'PR', 'PRh',
                                   'CT', 'CT2', 'CTnt',
                                   'MC', 'LC', 'MLCT', 'LMCT', 'LLCT', 
                                   'DEL', 'COH', 'COHh']
                for i in range(len(QMin['template']['theodore_prop'])):
                    for j in theodore_spelling:
                        if QMin['template']['theodore_prop'][i].lower()==j.lower():
                            QMin['template']['theodore_prop'][i]=j

            # TheoDORE fragments need to be parsed in a special way
            elif line[0]=='theodore_fragment':
                if '[' in orig:
                    string=orig.split(None,1)[1]
                    QMin['template']['theodore_fragment']=ast.literal_eval(string)
                else:
                    s=orig.split(None)[1:]
                    l=[]
                    for i in s:
                        l.append(int(i))
                    QMin['template']['theodore_fragment'].append(l)

            ## qmmm_table is a filename which needs to be checked
            #elif line[0]=='qmmm_table':
                #line2=orig.split(None,1)
                #if len(line2)<2:
                    #print 'Please specify a connection table file after "qmmm_table"!'
                    #sys.exit(50)
                #filename=os.path.abspath(os.path.expandvars(os.path.expanduser(line2[1])))
                #QMin['template']['qmmm_table']=filename

            ## qmmm_ff_file is a filename which needs to be checked
            #elif line[0]=='qmmm_ff_file':
                #line2=orig.split(None,1)
                #if len(line2)<2:
                    #print 'Please specify a force field file after "qmmm_ff_file"!'
                    #sys.exit(51)
                #filename=os.path.abspath(os.path.expandvars(os.path.expanduser(line2[1])))
                #QMin['template']['qmmm_ff_file']=filename


    ## read external basis set
    if QMin['template']['basis_external']:
      QMin['template']['basis']='gen'
      QMin['template']['basis_external']=readfile(QMin['template']['basis_external'])
      


    ##do logic checks
    #if QMin['template']['qmmm']:
        #filename=QMin['template']['qmmm_table']
        #if not os.path.isfile(filename):
            #print 'Connection table file "%s" does not exist!' % filename
            #sys.exit(52)
        #filename=QMin['template']['qmmm_ff_file']
        #if not os.path.isfile(filename):
            #print 'Force field file "%s" does not exist!' % filename
            #sys.exit(53)
        #print 'HINT: For QM/MM calculations, you have to specify in the template file the charge for the *QM region only*!'
        #print 'The automatic assignment of total charge might not work if the MM part is not neutral!\n'
    if not QMin['template']['unrestricted_triplets']:
        if len(QMin['template']['charge'])>=3 and QMin['template']['charge'][0]!=QMin['template']['charge'][2]:
            print 'Charges of singlets and triplets differ. Please enable the "unrestricted_triplets" option!'
            sys.exit(54)


    ## check the connection table file and frozen core
    #if QMin['template']['qmmm']:
        #out=readfile(QMin['template']['qmmm_table'])
        #link_atoms={}
        #qm_atoms={}
        #mm_atoms={}
        #found_mm=False
        #for iline,line in enumerate(out):
            #if 'subend' in line:
                #break
            #s=line.split()
            #if 'qm' in s[2].lower():
                #if found_mm:
                    #print 'In %s, all QM/LI atoms must occur consecutively at the beginning!' % QMin['template']['qmmm_table']
                    #sys.exit(55)
                #qm_atoms[iline]=QMin['geo'][iline][0]
            #elif 'li' in s[2].lower():
                #if found_mm:
                    #print 'In %s, all QM/LI atoms must occur consecutively at the beginning!' % QMin['template']['qmmm_table']
                    #sys.exit(56)
                #link_atoms[iline]=''
            #elif 'mm' in s[2].lower():
                #found_mm=True
                #mm_atoms[iline]=s[1]
        #nlink=len(link_atoms)
        #nqmatom=len(qm_atoms)+nlink
        #natom_table=nqmatom+len(mm_atoms)
        #if natom_table!=QMin['natom']:
            #print 'Number of atoms in connection table (%i) is inconsistent with %s (%i)!' % (natom_table,QMinfilename,QMin['natom'])
            #sys.exit(57)
        #if nlink>0:
            #links_found=False
            #for iline,line in enumerate(out):
                #if 'link_bonds' in line.lower() and not '!' in line:
                    #links_found=True
                    #break
            #if not links_found:
                #print 'Please add a "link_bonds" block to %s!' % (QMin['template']['qmmm_table'])
                #print '''Example:
#...
#6       H1      QM      4
#7       CT      LI      4       8       9       10
#8       HC      MM      7
#...
  #subend

  #link_bonds
#7 - 4 1.4 H H1'''
                #sys.exit(58)
            #while True:
                #iline+=1
                #if iline>=len(out):
                    #break
                #line=out[iline]
                #s=line.split()
                #n=int(s[0])
                #link_atoms[n-1]=s[4]
        #QMin['frozcore']=0 
        #atom_frags=set()
        #for i in qm_atoms:
            #QMin['frozcore']+=FROZENS[qm_atoms[i]]
            #atom_frags.add(qm_atoms[i])
        #for i in link_atoms:
            #QMin['frozcore']+=FROZENS[link_atoms[i]]
            #atom_frags.add(link_atoms[i])
        #QMin['atom_frags']=atom_frags




    # number of frozen core orbitals for wfoverlap
    # this call is down here because we might need to check for QM/MM in template first
    line=getsh2Gaukey(sh2Gau,'numfrozcore')
    if line[0]:
        numfroz=int(line[1])
        if numfroz==0:
            QMin['frozcore']=0 
        elif numfroz>0:
            QMin['frozcore']=numfroz
        elif numfroz<0:
            pass        # here we take frozcore from above



    # number of doubly occupied orbitals for Dyson
    line=getsh2Gaukey(sh2Gau,'numocc')
    if line[0]:
        numfroz=int(line[1])
        if numfroz<=0:
            QMin['ndocc']=0 
        elif numfroz>0:
            QMin['ndocc']=max(0,numfroz-QMin['frozcore'])
    else:
        QMin['ndocc']=0 

# --------------------------------------------- Logic ----------------------------------

    # obtain the statemap 
    statemap={}
    i=1
    for imult,istate,ims in itnmstates(QMin['states']):
        statemap[i]=[imult,istate,ims]
        i+=1
    QMin['statemap']=statemap

    # obtain the states to actually compute
    states_to_do=deepcopy(QMin['states'])
    for i in range(len(QMin['states'])):
        if states_to_do[i]>0:
            states_to_do[i]+=QMin['template']['paddingstates'][i]
    if not QMin['template']['unrestricted_triplets']:
        if len(QMin['states'])>=3 and QMin['states'][2]>0 and QMin['states'][0]<=1 :
            if 'soc' in QMin:
                states_to_do[0]=2
            else:
                states_to_do[0]=1
    QMin['states_to_do']=states_to_do

    # make the jobs
    jobs={}
    if QMin['states_to_do'][0]>0:
        jobs[1]={'mults':[1],'restr':True}
    if len(QMin['states_to_do'])>=2 and QMin['states_to_do'][1]>0:
        jobs[2]={'mults':[2],'restr':False}
    if len(QMin['states_to_do'])>=3 and QMin['states_to_do'][2]>0:
        if not QMin['template']['unrestricted_triplets'] and QMin['states_to_do'][0]>0:
            #jobs[1]['mults'].append(3)
            jobs[3]={'mults':[1,3],'restr':True}
        else:
            jobs[3]={'mults':[3],'restr':False}
    if len(QMin['states_to_do'])>=4:
        for imult,nstate in enumerate(QMin['states_to_do'][3:]):
            if nstate>0:
                jobs[len(jobs)+1]={'mults':[imult+4],'restr':False}
    QMin['jobs']=jobs

    # make the multmap (mapping between multiplicity and job)
    # multmap[imult]=ijob
    # multmap[-ijob]=[imults]
    multmap={}
    for ijob in jobs:
        job=jobs[ijob]
        for imult in job['mults']:
            multmap[imult]=ijob
        multmap[-(ijob)]=job['mults']
    multmap[1]=1
    QMin['multmap']=multmap

    # get the joblist
    joblist=set()
    for i in jobs:
        joblist.add(i)
    joblist=list(joblist)
    joblist.sort()
    QMin['joblist']=joblist
    njobs=len(joblist)
    QMin['njobs']=njobs

    # make the gsmap
    gsmap={}
    for i in range(QMin['nmstates']):
        m1,s1,ms1=tuple(QMin['statemap'][i+1])
        gs=(m1,1,ms1)
        job=QMin['multmap'][m1]
        if m1==3 and QMin['jobs'][job]['restr']:
            gs=(1,1,0.0)
        for j in range(QMin['nmstates']):
            m2,s2,ms2=tuple(QMin['statemap'][j+1])
            if (m2,s2,ms2)==gs:
                break
        gsmap[i+1]=j+1
    QMin['gsmap']=gsmap

    # get the set of states for which gradients actually need to be calculated
    gradmap=set()
    if 'grad' in QMin:
        for i in QMin['grad']:
            gradmap.add( tuple(statemap[i][0:2]) )
            #gradmap.add(      (statemap[i][0],1) )
    gradmap=list(gradmap)
    gradmap.sort()
    QMin['gradmap']=gradmap

    # make the chargemap
    QMin['chargemap']={}
    for i,c in enumerate(QMin['template']['charge']):
        QMin['chargemap'][i+1]=c

    # make the ionmap
    if 'ion' in QMin:
        ionmap=[]
        for m1 in itmult(QMin['states']):
            job1=QMin['multmap'][m1]
            el1=QMin['chargemap'][m1]
            for m2 in itmult(QMin['states']):
                if m1>=m2:
                    continue
                job2=QMin['multmap'][m2]
                el2=QMin['chargemap'][m2]
                #print m1,job1,el1,m2,job2,el2
                if abs(m1-m2)==1 and abs(el1-el2)==1:
                    ionmap.append( (m1,job1,m2,job2) )
        QMin['ionmap']=ionmap

    # number of properties/entries calculated by TheoDORE
    if 'theodore' in QMin:
        QMin['template']['theodore_n']=len(QMin['template']['theodore_prop']) + len(QMin['template']['theodore_fragment'])**2
    else:
        QMin['template']['theodore_n']=0

# --------------------------------------------- File setup ----------------------------------

    # check for initial orbitals
    initorbs={}
    if 'always_guess' in QMin:
        QMin['initorbs']={}
    elif 'init' in QMin or 'always_orb_init' in QMin:
        for job in QMin['joblist']:
            filename=os.path.join(QMin['pwd'],'GAUSSIAN.chk.init')
            if os.path.isfile(filename):
                initorbs[job]=filename
        for job in QMin['joblist']:
            filename=os.path.join(QMin['pwd'],'GAUSSIAN.chk.%i.init' % (job))
            if os.path.isfile(filename):
                initorbs[job]=filename
        if 'always_orb_init' in QMin and len(initorbs)<njobs:
            print 'Initial orbitals missing for some jobs!'
            sys.exit(59)
        QMin['initorbs']=initorbs
    elif 'newstep' in QMin:
        for job in QMin['joblist']:
            filename=os.path.join(QMin['savedir'],'GAUSSIAN.chk.%i' % (job))
            if os.path.isfile(filename):
                initorbs[job]=filename+'.old'     # file will be moved to .old
            else:
                print 'File %s missing in savedir!' % (filename)
                sys.exit(60)
        QMin['initorbs']=initorbs
    elif 'samestep' in QMin:
        for job in QMin['joblist']:
            filename=os.path.join(QMin['savedir'],'GAUSSIAN.chk.%i' % (job))
            if os.path.isfile(filename):
                initorbs[job]=filename
            else:
                print 'File %s missing in savedir!' % (filename)
                sys.exit(61)
        QMin['initorbs']=initorbs
    elif 'restart' in QMin:
        for job in QMin['joblist']:
            filename=os.path.join(QMin['savedir'],'GAUSSIAN.chk.%i.old' % (job))
            if os.path.isfile(filename):
                initorbs[job]=filename
            else:
                print 'File %s missing in savedir!' % (filename)
                sys.exit(62)
        QMin['initorbs']=initorbs


    # make name for backup directory
    if 'backup' in QMin:
      backupdir=QMin['savedir']+'/backup'
      backupdir1=backupdir
      i=0
      while os.path.isdir(backupdir1):
        i+=1
        if 'step' in QMin:
          backupdir1=backupdir+'/step%s_%i' % (QMin['step'][0],i)
        else:
          backupdir1=backupdir+'/calc_%i' % (i)
      QMin['backup']=backupdir

    if DEBUG:
        print '======= DEBUG print for QMin ======='
        pprint.pprint(QMin)
        print '===================================='
    return QMin

# =============================================================================================== #
# =============================================================================================== #
# =========================================== Job Scheduling ==================================== #
# =============================================================================================== #
# =============================================================================================== #

def parallel_speedup(N,scaling):
    # computes the parallel speedup from Amdahls law
    # with scaling being the fraction of parallelizable work and (1-scaling) being the serial part
    return 1./((1-scaling)+scaling/N)

def divide_slots(ncpu,ntasks,scaling):
    # this routine figures out the optimal distribution of the tasks over the CPU cores
    #   returns the number of rounds (how many jobs each CPU core will contribute to),
    #   the number of slots which should be set in the Pool,
    #   and the number of cores for each job.
    minpar=1
    ntasks_per_round=ncpu/minpar
    if ncpu==1:
        ntasks_per_round=1
    ntasks_per_round=min(ntasks_per_round,ntasks)
    optimal={}
    for i in range(1,1+ntasks_per_round):
        nrounds=int(math.ceil(float(ntasks)/i))
        ncores=ncpu/i
        optimal[i]=nrounds/parallel_speedup(ncores,scaling)
    #print optimal
    best=min(optimal,key=optimal.get)
    nrounds=int(math.ceil(float(ntasks)/best))
    ncores=ncpu/best

    cpu_per_run=[0 for i in range(ntasks)]
    if nrounds==1:
        itask=0
        for icpu in range(ncpu):
            cpu_per_run[itask]+=1
            itask+=1
            if itask>=ntasks:
                itask=0
        nslots=ntasks
    else:
        for itask in range(ntasks):
            cpu_per_run[itask]=ncores
        nslots=ncpu/ncores
    #print nrounds,nslots,cpu_per_run
    return nrounds,nslots,cpu_per_run

# =============================================================================================== #

def generate_joblist(QMin):

    #pprint.pprint(QMin)

    # sort the gradients into the different jobs
    gradjob={}
    for ijob in QMin['joblist']:
        gradjob['master_%i'%ijob]={}
    for grad in QMin['gradmap']:
        ijob=QMin['multmap'][grad[0]]
        isgs=False
        istates=QMin['states_to_do'][grad[0]-1]
        if not QMin['jobs'][ijob]['restr']:
            if grad[1]==1:
                isgs=True
        else:
            if grad==(1,1):
                isgs=True
        if isgs and istates>1:
            gradjob['grad_%i_%i'%grad]={}
            gradjob['grad_%i_%i'%grad][grad]={'gs':True}
        else:
            n=0
            for gradx in gradjob['master_%i'%ijob]:
                n+=1
            if n>0:
                gradjob['grad_%i_%i'%grad]={}
                gradjob['grad_%i_%i'%grad][grad]={'gs':False}
            else:
                gradjob['master_%i'%ijob][grad]={'gs':False}
    #pprint.pprint(gradjob)

    # make map for states onto gradjobs
    jobgrad={}
    for job in gradjob:
        for state in gradjob[job]:
            jobgrad[state]=(job,gradjob[job][state]['gs'])
    QMin['jobgrad']=jobgrad
    #print gradjob
    #print 
    #print jobgrad
    #sys.exit(63)

    schedule=[]
    QMin['nslots_pool']=[]

    # add the master calculations
    ntasks=0
    for i in gradjob:
        if 'master' in i:
            ntasks+=1
    nrounds,nslots,cpu_per_run=divide_slots(QMin['ncpu'],ntasks,QMin['schedule_scaling'])
    QMin['nslots_pool'].append(nslots)
    schedule.append({})
    icount=0
    for i in sorted(gradjob):
        if 'master' in i:
            QMin1=deepcopy(QMin)
            QMin1['master']=True
            QMin1['IJOB']=int(i.split('_')[1])
            remove=['gradmap','ncpu']
            for r in remove:
                QMin1=removekey(QMin1,r)
            QMin1['gradmap']=list(gradjob[i])
            QMin1['ncpu']=cpu_per_run[icount]
            if 3 in QMin['multmap'][-QMin1['IJOB']] and QMin['jobs'][QMin1['IJOB']]['restr']:
              QMin1['states'][0]=1
              QMin1['states_to_do'][0]=1
            icount+=1
            schedule[-1][i]=QMin1

    # add the gradient calculations
    ntasks=0
    for i in gradjob:
        if 'grad' in i:
            ntasks+=1
    if ntasks>0:
        nrounds,nslots,cpu_per_run=divide_slots(QMin['ncpu'],ntasks,QMin['schedule_scaling'])
        QMin['nslots_pool'].append(nslots)
        schedule.append({})
        icount=0
        for i in gradjob:
            if 'grad' in i:
                QMin1=deepcopy(QMin)
                mult=list(gradjob[i])[0][0]
                QMin1['IJOB']=QMin['multmap'][mult]
                remove=['gradmap','ncpu','h','soc','dm','overlap','ion','always_guess','always_orb_init','init']
                for r in remove:
                    QMin1=removekey(QMin1,r)
                QMin1['gradmap']=list(gradjob[i])
                QMin1['ncpu']=cpu_per_run[icount]
                QMin1['gradonly']=[]
                icount+=1
                schedule[-1][i]=QMin1

    return QMin,schedule


# =============================================================================================== #
# =============================================================================================== #
# ======================================= ADF Job Execution ===================================== #
# =============================================================================================== #
# =============================================================================================== #

def runjobs(schedule,QMin):

    if 'newstep' in QMin:
        moveOldFiles(QMin)

    print '>>>>>>>>>>>>> Starting the GAUSSIAN job execution'

    errorcodes={}
    for ijobset,jobset in enumerate(schedule):
        if not jobset:
            continue
        pool = Pool(processes=QMin['nslots_pool'][ijobset])
        for job in jobset:
            QMin1=jobset[job]
            WORKDIR=os.path.join(QMin['scratchdir'],job)

            errorcodes[job]=pool.apply_async(run_calc , [WORKDIR,QMin1])
            time.sleep(QMin['delay'])
        pool.close()
        pool.join()

    for i in errorcodes:
        errorcodes[i]=errorcodes[i].get()
    j=0
    string='Error Codes:\n'
    for i in errorcodes:
        string+='\t%s\t%i' % (i+' '*(10-len(i)),errorcodes[i])
        j+=1
        if j==4:
            j=0
            string+='\n'
    print string
    if any((i!=0 for i in errorcodes.values())):
        print 'Some subprocesses did not finish successfully!'
        print 'See %s:%s for error messages in GAUSSIAN output.' % (gethostname(),QMin['scratchdir'])
        sys.exit(64)
    print

    if PRINT:
        print '>>>>>>>>>>>>> Saving files'
        starttime=datetime.datetime.now()
    for ijobset,jobset in enumerate(schedule):
        if not jobset:
            continue
        for job in jobset:
            if 'master' in job:
                WORKDIR=os.path.join(QMin['scratchdir'],job)
                if not 'samestep' in QMin:
                    saveFiles(WORKDIR,jobset[job])
                if 'ion' in QMin and ijobset==0:
                    saveAOmatrix(WORKDIR,QMin)
    saveGeometry(QMin)
    if PRINT:
        endtime=datetime.datetime.now()
        print 'Saving Runtime: %s' % (endtime-starttime)
    print

    return errorcodes

# ======================================================================= #
def run_calc(WORKDIR,QMin):
    try:
        setupWORKDIR(WORKDIR,QMin)
        strip=True
        err=runGaussian(WORKDIR,QMin['GAUSS_EXE'],strip)
        err=0
    except Exception, problem:
        print '*'*50+'\nException in run_calc(%s)!' % (WORKDIR)
        traceback.print_exc()
        print '*'*50+'\n'
        raise problem

    return err

# ======================================================================= #
def setupWORKDIR(WORKDIR,QMin):
    # mkdir the WORKDIR, or clean it if it exists, then copy all necessary files from pwd and savedir
    # then put the GAUSSIAN.com file

    # setup the directory
    mkdir(WORKDIR)

    # write GAUSSIAN.com
    inputstring=writeGAUSSIANinput(QMin)
    filename=os.path.join(WORKDIR,'GAUSSIAN.com')
    writefile(filename,inputstring)
    if DEBUG:
        print '================== DEBUG input file for WORKDIR %s =================' % (shorten_DIR(WORKDIR))
        print inputstring
        print 'GAUSSIAN input written to: %s' % (filename)
        print '===================================================================='

    # wf file copying
    if 'master' in QMin:
        job=QMin['IJOB']
        if job in QMin['initorbs']:
            fromfile=QMin['initorbs'][job]
            tofile=os.path.join(WORKDIR,'GAUSSIAN.chk')
            shutil.copy(fromfile,tofile)
    elif 'grad' in QMin:
        job=QMin['IJOB']
        fromfile=os.path.join(QMin['scratchdir'],'master_%i' % job, 'GAUSSIAN.chk')
        tofile=os.path.join(WORKDIR,'GAUSSIAN.chk')
        shutil.copy(fromfile,tofile)

    ## force field file copying
    #if QMin['template']['qmmm']:
        #fromfile=QMin['template']['qmmm_ff_file']
        #tofile=os.path.join(WORKDIR,'ADF.ff')
        #shutil.copy(fromfile,tofile)

    return

# ======================================================================= #
def writeGAUSSIANinput(QMin):

    #pprint.pprint(QMin)

    # general setup
    job=QMin['IJOB']
    gsmult=QMin['multmap'][-job][0]
    restr=QMin['jobs'][job]['restr']
    charge=QMin['chargemap'][gsmult]

    # excited states to calculate
    states_to_do=QMin['states_to_do']
    for imult in range(len(states_to_do)):
        if not imult+1 in QMin['multmap'][-job]:
            states_to_do[imult]=0
    states_to_do[gsmult-1]-=1

    # do minimum number of states for gradient jobs
    if 'gradonly' in QMin:
        gradmult=QMin['gradmap'][0][0]
        gradstat=QMin['gradmap'][0][1]
        for imult in range(len(states_to_do)):
            if imult+1==gradmult:
                states_to_do[imult]=gradstat-(gradmult==gsmult)
            else:
                states_to_do[imult]=0

    # number of states to calculate
    if restr:
        ncalc=max(states_to_do)
        sing=states_to_do[0]>0
        trip=(len(states_to_do)>=3 and states_to_do[2]>0)
        if sing and trip:
          mults_td=',50-50'
        elif sing and not trip:
          mults_td=',singlets'
        elif trip and not sing:
          mults_td=',triplets'
    else:
        ncalc=max(states_to_do)
        mults_td=''

    ## whether to do SOC
    #sopert=False
    #gscorr=False
    #if 'soc' in QMin:
        #if restr:
            #nsing=QMin['states'][0]
            #if len(QMin['states'])>=3:
                #ntrip=QMin['states'][2]
            #else:
                #ntrip=0
            #if nsing+ntrip>=2 and ntrip>=1:
                #sopert=True
                #if nsing>=1:
                    #gscorr=True

    # gradients
    if QMin['gradmap']:
        dograd=True
        egrad=()
        for grad in QMin['gradmap']:
            if not (gsmult,1)==grad:
                egrad=grad
        #if multigrad_possible:
            #singgrad=[]
            #tripgrad=[]
            #for grad in QMin['gradmap']:
                #if not (gsmult,1)==grad:
                    #if grad[0]==gsmult:
                        #singgrad.append(grad[1]-1)
                    #if grad[0]==3 and restr:
                        #tripgrad.append(grad[1])
    else:
        dograd=False

    # construct the input string TODO
    string=''

    # link 0
    string+='%%MEM=%iMB\n' % (QMin['memory'])
    string+='%%NProcShared=%i\n' % (QMin['ncpu'])
    string+='%%Chk=%s\n' % ('GAUSSIAN.chk')
    if 'AOoverlap' in QMin or 'ion' in QMin:
        string+='%%Rwf=%s\n' % ('GAUSSIAN.rwf')
        if 'AOoverlap' in QMin:
            string+='%KJob l302\n'
    string+='\n'

    # Route section
    data=['p',
          'nosym',
          'unit=AU',
          QMin['template']['functional']
          ]
    if not QMin['template']['functional'].lower()=='dftba':
        data.append(QMin['template']['basis'])
    if dograd:
        data.append('force')
    if 'AOoverlap' in QMin:
        data.append('IOP(2/12=3)')
    if QMin['template']['dispersion']:
        data.append('EmpiricalDispersion=%s' % QMin['template']['dispersion'])
    if QMin['template']['grid']:
        data.append('int(grid=%s)' % QMin['template']['grid'])
    if QMin['template']['denfit']:
        data.append('denfit')
    if ncalc>0:
        if not QMin['template']['no_tda']:
            s='tda'
        else:
            s='td'
        s+='(nstates=%i%s' % (ncalc,mults_td)
        if dograd and egrad:
            s+=',root=%i' % (egrad[1]-(gsmult==egrad[0]))
        if not 'master' in QMin and 'grad' in QMin:
            s+=',read'
        s+=')'
        data.append(s)
    if QMin['template']['scrf']:
        s=','.join(QMin['template']['scrf'].split())
        data.append('scrf(%s)' % s)
    if QMin['template']['scf']:
        s=','.join(QMin['template']['scf'].split())
        data.append('scf(%s)' % s)
    if QMin['template']['iop']:
        s=','.join(QMin['template']['iop'].split())
        data.append('iop(%s)' % s)
    if QMin['template']['keys']:
        data.extend([QMin['template']['keys']])
    if 'theodore' in QMin:
        data.append('pop=full')
        data.append('IOP(9/40=3)')
        data.append('GFPRINT')
    string+='#'
    for i in data:
        string+=i+'\n'


    # title
    string+='\nSHARC-GAUSSIAN job\n\n'

    # charge/mult and geometry
    if 'AOoverlap' in QMin:
      string+='%i %i\n' % (2.*charge,gsmult)
    else:
      string+='%i %i\n' % (charge,gsmult)
    for iatom,atom in enumerate(QMin['geo']):
        label=atom[0]
        string+='%4s %16.9f %16.9f %16.9f\n' % (label,atom[1],atom[2],atom[3])
    string+='\n'
    if QMin['template']['functional'].lower()=='dftba':
      string+='@GAUSS_EXEDIR:dftba.prm\n'
    if QMin['template']['basis_external']:
      for line in QMin['template']['basis_external']:
        string+=line
    string+='\n\n'


    return string

# ======================================================================= #
def shorten_DIR(string):
    maxlen=40
    front=12
    if len(string)>maxlen:
        return string[0:front]+'...'+string[-(maxlen-3-front):]
    else:
        return string+' '*(maxlen-len(string))

# ======================================================================= #
def runGaussian(WORKDIR,GAUSS_EXE,strip=False):
    prevdir=os.getcwd()
    os.chdir(WORKDIR)
    string=GAUSS_EXE+' '
    string+='< GAUSSIAN.com'
    if PRINT or DEBUG:
        starttime=datetime.datetime.now()
        sys.stdout.write('START:\t%s\t%s\t"%s"\n' % (shorten_DIR(WORKDIR),starttime,shorten_DIR(string)))
        sys.stdout.flush()
    stdoutfile=open(os.path.join(WORKDIR,'GAUSSIAN.log'),'w')
    stderrfile=open(os.path.join(WORKDIR,'GAUSSIAN.err'),'w')
    try:
        runerror=sp.call(string,shell=True,stdout=stdoutfile,stderr=stderrfile)
    except OSError:
        print 'Call have had some serious problems:',OSError
        sys.exit(65)
    stdoutfile.close()
    stderrfile.close()
    #stderr=readfile(os.path.join(WORKDIR,'ADF.err'))
    #for line in stderr:
        #if 'error' in line.lower():
            #sys.stdout.write('ERROR: \t%s\t"%s"\n' % (shorten_DIR(WORKDIR),line.strip()))
            #runerror+=1
            #break
    if PRINT or DEBUG:
        endtime=datetime.datetime.now()
        sys.stdout.write('FINISH:\t%s\t%s\tRuntime: %s\tError Code: %i\n' % (shorten_DIR(WORKDIR),endtime,endtime-starttime,runerror))
        sys.stdout.flush()
    os.chdir(prevdir)
    if strip and not DEBUG and runerror==0:
        stripWORKDIR(WORKDIR)
    return runerror

# ======================================================================= #
def stripWORKDIR(WORKDIR):
    ls=os.listdir(WORKDIR)
    keep=['GAUSSIAN.com$','GAUSSIAN.err$','GAUSSIAN.log$','GAUSSIAN.chk','GAUSSIAN.fchk','GAUSSIAN.rwf']
    for ifile in ls:
        delete=True
        for k in keep:
            if containsstring(k,ifile):
                delete=False
        if delete:
            rmfile=os.path.join(WORKDIR,ifile)
            if not DEBUG:
                os.remove(rmfile)

# ======================================================================= #
def moveOldFiles(QMin):
    # moves all relevant files in the savedir to old files (per job)
    if PRINT:
        print '>>>>>>>>>>>>> Moving old files'
    basenames=['GAUSSIAN.chk']
    if not 'nooverlap' in QMin:
        basenames.append('mos')
    for job in QMin['joblist']:
        for base in basenames:
            fromfile=os.path.join(QMin['savedir'],'%s.%i' % (base,job))
            if not os.path.isfile(fromfile):
                print 'File %s not found, cannot move to OLD!' % (fromfile)
                sys.exit(66)
            tofile    =os.path.join(QMin['savedir'],'%s.%i.old' % (base,job))
            if PRINT:
                print shorten_DIR(fromfile)+'   =>   '+shorten_DIR(tofile)
            shutil.copy(fromfile,tofile)
    # moves all relevant files in the savedir to old files (per mult)
    basenames=[]
    if not 'nooverlap' in QMin:
        basenames=['dets']
    for job in itmult(QMin['states']):
        for base in basenames:
            fromfile=os.path.join(QMin['savedir'],'%s.%i' % (base,job))
            if not os.path.isfile(fromfile):
                print 'File %s not found, cannot move to OLD!' % (fromfile)
                sys.exit(67)
            tofile    =os.path.join(QMin['savedir'],'%s.%i.old' % (base,job))
            if PRINT:
                print shorten_DIR(fromfile)+'   =>   '+shorten_DIR(tofile)
            shutil.copy(fromfile,tofile)
    # geometry file
    fromfile=os.path.join(QMin['savedir'],'geom.dat')
    tofile  =os.path.join(QMin['savedir'],'geom.dat.old')
    if PRINT:
      print shorten_DIR(fromfile)+'   =>   '+shorten_DIR(tofile)
    shutil.copy(fromfile,tofile)
    # also remove aoovl files if present
    delete=['AO_overl','AO_overl.mixed']
    for f in delete:
        rmfile=os.path.join(QMin['savedir'],f)
        if os.path.isfile(rmfile):
            os.remove(rmfile)
            if PRINT:
                print 'rm '+rmfile
    print

# ======================================================================= #
def saveGeometry(QMin):
  string=''
  for iatom,atom in enumerate(QMin['geo']):
    label=atom[0]
    string+='%4s %16.9f %16.9f %16.9f\n' % (label,atom[1],atom[2],atom[3])
  filename=os.path.join(QMin['savedir'],'geom.dat')
  writefile(filename,string)
  if PRINT:
    print shorten_DIR(filename)
  return

# ======================================================================= #
def saveFiles(WORKDIR,QMin):

    # copy the TAPE21 from master directories
    job=QMin['IJOB']
    fromfile=os.path.join(WORKDIR,'GAUSSIAN.chk')
    tofile=os.path.join(QMin['savedir'],'GAUSSIAN.chk.%i' % (job))
    shutil.copy(fromfile,tofile)
    if PRINT:
        print shorten_DIR(tofile)

    # if necessary, extract the MOs and write them to savedir
    if 'ion' in QMin or not 'nooverlap' in QMin:
        f=os.path.join(WORKDIR,'GAUSSIAN.chk')
        string=get_MO_from_chk(f,QMin)
        mofile=os.path.join(QMin['savedir'],'mos.%i' % job)
        writefile(mofile,string)
        if PRINT:
            print shorten_DIR(mofile)

    # if necessary, extract the TDDFT coefficients and write them to savedir
    if 'ion' in QMin or not 'nooverlap' in QMin:
        f=os.path.join(WORKDIR,'GAUSSIAN.chk')
        strings=get_dets_from_chk(f,QMin)
        for f in strings:
            writefile(f,strings[f])
            if PRINT:
                print shorten_DIR(f)

# ======================================================================= #
def get_rwfdump(groot,filename,number):
  WORKDIR=os.path.dirname(filename)
  prevdir=os.getcwd()
  os.chdir(WORKDIR)
  dumpname='rwfdump.txt'
  string='%s/rwfdump %s %s %s' % (groot,os.path.basename(filename),dumpname,number)
  #print string
  try_shells=['sh','bash','csh','tcsh']
  for shell in try_shells:
    try:
      runerror=sp.call(string,shell=True,executable=shell)
    except OSError:
      print 'Gaussian rwfdump has serious problems:',OSError
      sys.exit(68)
  string=readfile(dumpname)
  os.chdir(prevdir)
  return string


# ======================================================================= #
def get_MO_from_chk(filename,QMin):

    job=QMin['IJOB']
    restr=QMin['jobs'][job]['restr']

    # extract alpha orbitals
    data=get_rwfdump(QMin['groot'],filename,'524R')
    for iline,line in enumerate(data):
      if "Dump of file" in line:
        break
    mocoef_A=[]
    while True:
      iline+=1
      if iline>=len(data):
        break
      s=data[iline].split()
      for i in s:
        mocoef_A.append(float(i.replace('D','E')))
    NAO=int(math.sqrt(len(mocoef_A)))
    NMO_A=NAO
    MO_A=[ mocoef_A[NAO*i:NAO*(i+1)] for i in range(NAO) ]

    # extract beta orbitals
    if not restr:
      data=get_rwfdump(QMin['groot'],filename,'526R')
      for iline,line in enumerate(data):
        if "Dump of file" in line:
          break
      mocoef_B=[]
      while True:
        iline+=1
        if iline>=len(data):
          break
        s=data[iline].split()
        for i in s:
          mocoef_B.append(float(i.replace('D','E')))
      if not NAO==int(math.sqrt(len(mocoef_B))):
        print 'Problem in orbital reading!'
        sys.exit(69)
      NMO_B=NAO
      MO_B=[ mocoef_B[NAO*i:NAO*(i+1)] for i in range(NAO) ]


    NMO=NMO_A      -  QMin['frozcore']
    if restr:
        NMO=NMO_A      -  QMin['frozcore']
    else:
        NMO=NMO_A+NMO_B-2*QMin['frozcore']

    # make string
    string='''2mocoef
header
 1
MO-coefficients from Gaussian
 1
 %i   %i
 a
mocoef
(*)
''' % (NAO,NMO)
    x=0
    for imo,mo in enumerate(MO_A):
        if imo<QMin['frozcore']:
            continue
        for c in mo:
            if x>=3:
                string+='\n'
                x=0
            string+='% 6.12e ' % c
            x+=1
        if x>0:
            string+='\n'
            x=0
    if not restr:
        x=0
        for imo,mo in enumerate(MO_B):
            if imo<QMin['frozcore']:
                continue
            for c in mo:
                if x>=3:
                    string+='\n'
                    x=0
                string+='% 6.12e ' % c
                x+=1
            if x>0:
                string+='\n'
                x=0
    string+='orbocc\n(*)\n'
    x=0
    for i in range(NMO):
        if x>=3:
            string+='\n'
            x=0
        string+='% 6.12e ' % (0.0)
        x+=1

    return string

# ======================================================================= #
def get_dets_from_chk(filename,QMin):

    # get general infos
    job=QMin['IJOB']
    restr=QMin['jobs'][job]['restr']
    mults=QMin['jobs'][job]['mults']
    if 3 in mults:
        mults=[3]
    gsmult=QMin['multmap'][-job][0]
    nstates_to_extract=deepcopy(QMin['states'])
    for i in range(len(nstates_to_extract)):
        if not i+1 in mults:
            nstates_to_extract[i]=0
        elif i+1==gsmult:
            nstates_to_extract[i]-=1

    # get infos from logfile
    logfile=os.path.join(os.path.dirname(filename),'GAUSSIAN.log')
    data=readfile(logfile)
    infos={}
    for iline,line in enumerate(data):
      if 'NBsUse=' in line:
        s=line.split()
        infos['nbsuse']=int(s[1])
      if 'Range of M.O.s used for correlation:' in line:
        for i in [1,2]:
          s=data[iline+i].replace('=',' ').split()
          for j in range(5):
            infos[s[2*j]]=int(s[2*j+1])

    if not 'NOA' in  infos:
      nstates_onfile=0
      charge=QMin['chargemap'][gsmult]
      nelec=float(QMin['Atomcharge']-charge)
      infos['NOA']=int(nelec/2. + float(gsmult-1)/2. )
      infos['NOB']=int(nelec/2. - float(gsmult-1)/2. )
      infos['NVA']=infos['nbsuse']-infos['NOA']
      infos['NVB']=infos['nbsuse']-infos['NOB']
      infos['NFC']=0
    else:
      # get all info from checkpoint
      data=get_rwfdump(QMin['groot'],filename,'635R')
      for iline,line in enumerate(data):
        if "Dump of file" in line:
          break
      eigenvectors_array=[]
      while True:
        iline+=1
        if iline>=len(data):
          break
        s=data[iline].split()
        for i in s:
          try:
            eigenvectors_array.append(float(i.replace('D','E')))
          except ValueError:
            eigenvectors_array.append(float('NaN'))
      nstates_onfile=(len(eigenvectors_array)-12)/(4+8*(infos['NOA']*infos['NVA']+infos['NOB']*infos['NVB']))
    #print nstates_onfile
    #print len(eigenvectors_array)
    #print infos

    # get ground state configuration
    # make step vectors (0:empty, 1:alpha, 2:beta, 3:docc)
    if restr:
        occ_A=[ 3 for i in range(infos['NFC']+infos['NOA']) ]+[ 0 for i in range(infos['NVA']) ]
    if not restr:
        occ_A=[ 1 for i in range(infos['NFC']+infos['NOA']) ]+[ 0 for i in range(infos['NVA']) ]
        occ_B=[ 2 for i in range(infos['NFC']+infos['NOB']) ]+[ 0 for i in range(infos['NVB']) ]
    occ_A=tuple(occ_A)
    if not restr:
        occ_B=tuple(occ_B)

    # get infos
    nocc_A=infos['NOA']
    nvir_A=infos['NVA']
    nocc_B=infos['NOB']
    nvir_B=infos['NVB']

    # get eigenvectors
    eigenvectors={}
    for imult,mult in enumerate(mults):
        eigenvectors[mult]=[]
        if mult==gsmult:
            # add ground state
            if restr:
                key=tuple(occ_A[QMin['frozcore']:])
            else:
                key=tuple(occ_A[QMin['frozcore']:]+occ_B[QMin['frozcore']:])
            eigenvectors[mult].append( {key:1.0} )
        for istate in range(nstates_to_extract[mult-1]):
            # get X+Y vector
            startindex=12+istate*(nvir_A*nocc_A+nvir_B*nocc_B)
            endindex  =startindex+nvir_A*nocc_A+nvir_B*nocc_B
            eig=[ i for i in eigenvectors_array[startindex:endindex] ]
            # get X-Y vector
            startindex=12+istate*(nvir_A*nocc_A+nvir_B*nocc_B)+4*nstates_onfile*(nvir_A*nocc_A+nvir_B*nocc_B)
            endindex  =startindex+nvir_A*nocc_A+nvir_B*nocc_B
            eigl=[ i for i in eigenvectors_array[startindex:endindex] ]
            # get X vector
            for i in range(len(eig)):
              eig[i]=(eig[i]+eigl[i])/2.
            # make dictionary
            dets={}
            if restr:
                for iocc in range(nocc_A):
                    for ivirt in range(nvir_A):
                        index=iocc*nvir_A+ivirt
                        dets[ (iocc,ivirt,1) ]=eig[index]
            else:
                for iocc in range(nocc_A):
                    for ivirt in range(nvir_A):
                        index=iocc*nvir_A+ivirt
                        dets[ (iocc,ivirt,1) ]=eig[index]
                for iocc in range(nocc_B):
                    for ivirt in range(nvir_B):
                        index=iocc*nvir_B+ivirt + nvir_A*nocc_A
                        dets[ (iocc,ivirt,2) ]=eig[index]
            #pprint.pprint(dets)
            # truncate vectors
            norm=0.
            for k in sorted(dets,key=lambda x: dets[x]**2,reverse=True):
                if restr:
                    factor=0.5
                else:
                    factor=1.
                if norm>factor*QMin['wfthres']:
                    del dets[k]
                    continue
                norm+=dets[k]**2
            #pprint.pprint(dets)
            # create strings and expand singlets
            dets2={}
            if restr:
                for iocc,ivirt,dummy in dets:
                    # singlet
                    if mult==1:
                        # alpha excitation
                        key=list(occ_A)
                        key[infos['NFC']+iocc]=2
                        key[infos['NFC']+nocc_A+ivirt]=1
                        dets2[tuple(key)]=dets[ (iocc,ivirt,dummy) ]
                        # beta excitation
                        key[infos['NFC']+iocc]=1
                        key[infos['NFC']+nocc_A+ivirt]=2
                        dets2[tuple(key)]=dets[ (iocc,ivirt,dummy) ]
                    # triplet
                    elif mult==3:
                        key=list(occ_A)
                        key[infos['NFC']+iocc]=1
                        key[infos['NFC']+nocc_A+ivirt]=1
                        dets2[tuple(key)]=dets[ (iocc,ivirt,dummy) ]*math.sqrt(2.)
            else:
                for iocc,ivirt,dummy in dets:
                    if dummy==1:
                        key=list(occ_A+occ_B)
                        key[infos['NFC']+iocc]=0
                        key[infos['NFC']+nocc_A+ivirt]=1
                        dets2[tuple(key)]=dets[ (iocc,ivirt,dummy) ]
                    elif dummy==2:
                        key=list(occ_A+occ_B)
                        key[2*infos['NFC']+nocc_A+nvir_A + iocc]=0
                        key[2*infos['NFC']+nocc_A+nvir_A + nocc_B+ivirt]=2
                        dets2[tuple(key)]=dets[ (iocc,ivirt,dummy) ]
            #pprint.pprint(dets2)
            # remove frozen core
            dets3={}
            for key in dets2:
                problem=False
                if restr:
                    if any( [key[i]!=3 for i in range(QMin['frozcore']) ] ):
                        problem=True
                else:
                    if any( [key[i]!=1 for i in range(QMin['frozcore']) ] ):
                        problem=True
                    if any( [key[i]!=2 for i in range(nocc_A+nvir_A+QMin['frozcore'], nocc_A+nvir_A + 2*QMin['frozcore']) ] ):
                        problem=True
                if problem:
                    print 'WARNING: Non-occupied orbital inside frozen core! Skipping ...'
                    continue
                    #sys.exit(70)
                if restr:
                    key2=key[QMin['frozcore']:]
                else:
                    key2=key[QMin['frozcore']:QMin['frozcore']+nocc_A+nvir_A] + key[nocc_A+nvir_A+2*QMin['frozcore']:]
                dets3[key2]=dets2[key]
            #pprint.pprint(dets3)
            # append
            eigenvectors[mult].append(dets3)

    strings={}
    for imult,mult in enumerate(mults):
        filename=os.path.join(QMin['savedir'],'dets.%i' % mult)
        strings[filename]=format_ci_vectors(eigenvectors[mult])

    return strings

# ======================================================================= #
def format_ci_vectors(ci_vectors):

    # get nstates, norb and ndets
    alldets=set()
    for dets in ci_vectors:
        for key in dets:
            alldets.add(key)
    ndets=len(alldets)
    nstates=len(ci_vectors)
    norb=len(next(iter(alldets)))

    string='%i %i %i\n' % (nstates,norb,ndets)
    for det in sorted(alldets,reverse=True):
        for o in det:
            if o==0:
                string+='e'
            elif o==1:
                string+='a'
            elif o==2:
                string+='b'
            elif o==3:
                string+='d'
        for istate in range(len(ci_vectors)):
            if det in ci_vectors[istate]:
                string+=' %11.7f ' % ci_vectors[istate][det]
            else:
                string+=' %11.7f ' % 0.
        string+='\n'
    return string

# ======================================================================= #
def saveAOmatrix(WORKDIR,QMin):
    filename=os.path.join(WORKDIR,'GAUSSIAN.rwf')
    NAO,Smat=get_smat(filename,QMin['groot'])

    string='%i %i\n' % (NAO,NAO)
    for irow in range(NAO):
        for icol in range(NAO):
            string+='% .15e ' % (Smat[icol][irow])
        string+='\n'
    filename=os.path.join(QMin['savedir'],'AO_overl')
    writefile(filename,string)
    if PRINT:
        print shorten_DIR(filename)

# ======================================================================= #
def get_smat(filename,groot):

    # get all info from checkpoint
    data=get_rwfdump(groot,filename,'514R')

    # extract matrix
    for iline,line in enumerate(data):
      if "Dump of file" in line:
        break
    Smat=[]
    while True:
      iline+=1
      if iline>=len(data):
        break
      s=data[iline].split()
      for i in s:
        Smat.append(float(i.replace('D','E')))
    NAO=int(math.sqrt(2.*len(Smat)+0.25)-0.5)

    # Smat is lower triangular matrix, len is NAO*(NAO+1)/2
    ao_ovl=makermatrix(NAO,NAO)
    x=0
    y=0
    for el in Smat:
        ao_ovl[x][y]=el
        ao_ovl[y][x]=el
        x+=1
        if x>y:
            x=0
            y+=1
    return NAO,ao_ovl

# ======================================================================= #
def mkdir(DIR):
    # mkdir the DIR, or clean it if it exists
    if os.path.exists(DIR):
        if os.path.isfile(DIR):
            print '%s exists and is a file!' % (DIR)
            sys.exit(71)
        elif os.path.isdir(DIR):
            if DEBUG:
                print 'Remake\t%s' % DIR
            shutil.rmtree(DIR)
            os.makedirs(DIR)
    else:
        try:
            if DEBUG:
                print 'Make\t%s' % DIR
            os.makedirs(DIR)
        except OSError:
            print 'Can not create %s\n' % (DIR)
            sys.exit(72)

# ======================================================================= #

def link(PATH,NAME,crucial=True,force=True):
    # do not create broken links
    if not os.path.exists(PATH) and crucial:
        print 'Source %s does not exist, cannot create link!' % (PATH)
        sys.exit(73)
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
                print '%s exists, cannot create a link of the same name!' % (NAME)
                if crucial:
                    sys.exit(74)
                else:
                    return
    elif os.path.exists(NAME):
        # NAME is not a link. The interface will not overwrite files/directories with links, even with force=True
        print '%s exists, cannot create a link of the same name!' % (NAME)
        if crucial:
            sys.exit(75)
        else:
            return
    os.symlink(PATH, NAME)

# =============================================================================================== #
# =============================================================================================== #
# =======================================  TheoDORE ============================================= #
# =============================================================================================== #
# =============================================================================================== #

def run_theodore(QMin,errorcodes):

    if 'theodore' in QMin:
        print '>>>>>>>>>>>>> Starting the TheoDORE job execution'

        for ijob in QMin['jobs']:
            if not QMin['jobs'][ijob]['restr']:
                if DEBUG:
                    print 'Skipping Job %s because it is unrestricted.' % (ijob)
                continue
            else:
                mults=QMin['jobs'][ijob]['mults']
                gsmult=mults[0]
                ns=0
                for i in mults:
                    ns+=QMin['states'][i-1]-(i==gsmult)
                if ns==0:
                    if DEBUG:
                        print 'Skipping Job %s because it contains no excited states.' % (ijob)
                    continue
            WORKDIR=os.path.join(QMin['scratchdir'],'master_%i' % ijob)
            setupWORKDIR_TH(WORKDIR,QMin)
            os.environ
            errorcodes['theodore_%i' % ijob]=runTHEODORE(WORKDIR,QMin['theodir'])

        # Error code handling
        j=0
        string='Error Codes:\n'
        for i in errorcodes:
            if 'theodore' in i:
                string+='\t%s\t%i' % (i+' '*(10-len(i)),errorcodes[i])
                j+=1
                if j==4:
                    j=0
                    string+='\n'
        print string
        if any((i!=0 for i in errorcodes.values())):
            print 'Some subprocesses did not finish successfully!'
            sys.exit(76)

        print ''

    return errorcodes

# ======================================================================= #
def setupWORKDIR_TH(WORKDIR,QMin):
    # mkdir the WORKDIR, or clean it if it exists, then copy all necessary files from pwd and savedir

    # write dens_ana.in
    inputstring='''rtype='cclib'
rfile='GAUSSIAN.log'
jmol_orbitals=False
molden_orbitals=False
Om_formula=2
eh_pop=1
comp_ntos=True
print_OmFrag=True
output_file='tden_summ.txt'
prop_list=%s
at_lists=%s
''' % (str(QMin['template']['theodore_prop']),str(QMin['template']['theodore_fragment']))

    filename=os.path.join(WORKDIR,'dens_ana.in')
    writefile(filename,inputstring)
    if DEBUG:
        print '================== DEBUG input file for WORKDIR %s =================' % (shorten_DIR(WORKDIR))
        print inputstring
        print 'TheoDORE input written to: %s' % (filename)
        print '===================================================================='

    return

# ======================================================================= #
def runTHEODORE(WORKDIR,THEODIR):
    prevdir=os.getcwd()
    os.chdir(WORKDIR)
    string='python2 '+os.path.join(THEODIR,'bin','analyze_tden.py')
    stdoutfile=open(os.path.join(WORKDIR,'theodore.out'),'w')
    stderrfile=open(os.path.join(WORKDIR,'theodore.err'),'w')
    if PRINT or DEBUG:
        starttime=datetime.datetime.now()
        sys.stdout.write('START:\t%s\t%s\t"%s"\n' % (shorten_DIR(WORKDIR),starttime,shorten_DIR(string)))
        sys.stdout.flush()
    try:
        runerror=sp.call(string,shell=True,stdout=stdoutfile,stderr=stderrfile)
    except OSError:
        print 'Call have had some serious problems:',OSError
        sys.exit(77)
    stdoutfile.close()
    stderrfile.close()
    if PRINT or DEBUG:
        endtime=datetime.datetime.now()
        sys.stdout.write('FINISH:\t%s\t%s\tRuntime: %s\tError Code: %i\n' % (shorten_DIR(WORKDIR),endtime,endtime-starttime,runerror))
        sys.stdout.flush()
    os.chdir(prevdir)
    return runerror

# =============================================================================================== #
# =============================================================================================== #
# =======================================  Dyson and overlap calcs ============================== #
# =============================================================================================== #
# =============================================================================================== #

def run_wfoverlap(QMin,errorcodes):

    print '>>>>>>>>>>>>> Starting the WFOVERLAP job execution'

    # do Dyson calculations
    if 'ion' in QMin:
        for ionpair in QMin['ionmap']:
            WORKDIR=os.path.join(QMin['scratchdir'],'Dyson_%i_%i_%i_%i' % ionpair)
            files={'aoovl':'AO_overl', 
                   'det.a': 'dets.%i' % ionpair[0],
                   'det.b': 'dets.%i' % ionpair[2],
                   'mo.a':    'mos.%i' % ionpair[1],
                   'mo.b':    'mos.%i' % ionpair[3] }
            setupWORKDIR_WF(WORKDIR,QMin,files)
            errorcodes['Dyson_%i_%i_%i_%i' % ionpair]=runWFOVERLAP(WORKDIR,QMin['wfoverlap'],memory=QMin['memory'],ncpu=QMin['ncpu'])

    # do overlap calculations
    if 'overlap' in QMin:
        get_Double_AOovl(QMin)
        for m in itmult(QMin['states']):
            job=QMin['multmap'][m]
            WORKDIR=os.path.join(QMin['scratchdir'],'WFOVL_%i_%i' % (m,job))
            files={'aoovl':'AO_overl.mixed', 
                         'det.a': 'dets.%i.old' % m,
                         'det.b': 'dets.%i' % m,
                         'mo.a':    'mos.%i.old' % job,
                         'mo.b':    'mos.%i' % job }
            setupWORKDIR_WF(WORKDIR,QMin,files)
            errorcodes['WFOVL_%i_%i' % (m,job)]=runWFOVERLAP(WORKDIR,QMin['wfoverlap'],memory=QMin['memory'],ncpu=QMin['ncpu'])

    # Error code handling
    j=0
    string='Error Codes:\n'
    for i in errorcodes:
        if 'Dyson' in i or 'WFOVL' in i:
            string+='\t%s\t%i' % (i+' '*(10-len(i)),errorcodes[i])
            j+=1
            if j==4:
                j=0
                string+='\n'
    print string
    if any((i!=0 for i in errorcodes.values())):
        print 'Some subprocesses did not finish successfully!'
        sys.exit(78)

    print ''

    return errorcodes

# ======================================================================= #
def setupWORKDIR_WF(WORKDIR,QMin,files):
    # mkdir the WORKDIR, or clean it if it exists, then copy all necessary files from pwd and savedir

    # setup the directory
    mkdir(WORKDIR)

    # write wfovl.inp
    inputstring='''mix_aoovl=aoovl
a_mo=mo.a
b_mo=mo.b
a_det=det.a
b_det=det.b
a_mo_read=0
b_mo_read=0
ao_read=0
'''
    if 'ion' in QMin:
        if QMin['ndocc']>0:
            inputstring+='ndocc=%i\n' % (QMin['ndocc'])
    if QMin['ncpu']>=8:
        inputstring+='force_direct_dets\n'
    filename=os.path.join(WORKDIR,'wfovl.inp')
    writefile(filename,inputstring)
    if DEBUG:
        print '================== DEBUG input file for WORKDIR %s =================' % (shorten_DIR(WORKDIR))
        print inputstring
        print 'wfoverlap input written to: %s' % (filename)
        print '===================================================================='

    # link input files from save
    linkfiles=[ 'aoovl', 'det.a', 'det.b', 'mo.a', 'mo.b' ]
    for f in linkfiles:
        fromfile=os.path.join(QMin['savedir'],files[f])
        tofile    =os.path.join(WORKDIR,f)
        link(fromfile,tofile)

    return

# ======================================================================= #
def runWFOVERLAP(WORKDIR,WFOVERLAP,memory=100,ncpu=1):
    prevdir=os.getcwd()
    os.chdir(WORKDIR)
    string=WFOVERLAP+' -m %i' % (memory)+' -f wfovl.inp'
    stdoutfile=open(os.path.join(WORKDIR,'wfovl.out'),'w')
    stderrfile=open(os.path.join(WORKDIR,'wfovl.err'),'w')
    os.environ['OMP_NUM_THREADS']=str(ncpu)
    if PRINT or DEBUG:
        starttime=datetime.datetime.now()
        sys.stdout.write('START:\t%s\t%s\t"%s"\n' % (shorten_DIR(WORKDIR),starttime,shorten_DIR(string)))
        sys.stdout.flush()
    try:
        runerror=sp.call(string,shell=True,stdout=stdoutfile,stderr=stderrfile)
    except OSError:
        print 'Call have had some serious problems:',OSError
        sys.exit(79)
    stdoutfile.close()
    stderrfile.close()
    if PRINT or DEBUG:
        endtime=datetime.datetime.now()
        sys.stdout.write('FINISH:\t%s\t%s\tRuntime: %s\tError Code: %i\n' % (shorten_DIR(WORKDIR),endtime,endtime-starttime,runerror))
        sys.stdout.flush()
    os.chdir(prevdir)
    return runerror


# ======================================================================= #
def get_Double_AOovl(QMin):


    # get geometries
    filename1=os.path.join(QMin['savedir'],'geom.dat')
    oldgeo=get_geometry(filename1)
    filename2=os.path.join(QMin['savedir'],'geom.dat.old')
    newgeo=get_geometry(filename2)

    ## apply shift
    #shift=1e-5
    #for iatom in range(len(oldgeo)):
        #for ixyz in range(3):
            #oldgeo[iatom][1+ixyz]+=shift

    # build QMin   # TODO: always singlet for AOoverlaps
    QMin1=deepcopy(QMin)
    QMin1['geo']=oldgeo+newgeo
    QMin1['AOoverlap']=[filename1,filename2]
    QMin1['IJOB']=QMin['joblist'][0]
    QMin1['natom']=len(newgeo)
    remove=['nacdr','grad','h','soc','dm','overlap','ion']
    for r in remove:
        QMin1=removekey(QMin1,r)

    # run the calculation
    WORKDIR=os.path.join(QMin['scratchdir'],'AOoverlap')
    err=run_calc(WORKDIR,QMin1)

    # get output
    filename=os.path.join(WORKDIR,'GAUSSIAN.rwf')
    NAO,Smat=get_smat(filename,QMin['groot'])

    # adjust the diagonal blocks for DFTB-A
    if QMin['template']['functional']=='dftba':
      Smat=adjust_DFTB_Smat(Smat,NAO,QMin)

    ## Smat is now full matrix NAO*NAO
    ## we want the lower left quarter, but transposed
    string='%i %i\n' % (NAO/2,NAO/2)
    for irow in range(NAO/2,NAO):
        for icol in range(0,NAO/2):
            string+='% .15e ' % (Smat[icol][irow])          # note the exchanged indices => transposition
        string+='\n'
    filename=os.path.join(QMin['savedir'],'AO_overl.mixed')
    writefile(filename,string)

    return

# ======================================================================= #
def get_geometry(filename):
    data=readfile(filename)
    geometry=[]
    for line in data:
      s=line.split()
      geometry.append( [ s[0], float(s[1]), float(s[2]), float(s[3]) ] )
    return geometry

# ======================================================================= #
def adjust_DFTB_Smat(Smat,NAO,QMin):
  # list with the number of basis functions for basis set VSTO-6G* (used for DFTBA in Gaussian)
  nbasis={1: ['h','he'],
          2: ['li','be'],
          4: ['b','c','n','o','f','ne']}
  nbs={}
  for i in nbasis:
    for el in nbasis[i]:
      nbs[el]=i
  Nb=0
  itot=0
  mapping={}
  for ii,i in enumerate(QMin['geo']):
    try:
      Nb+=nbs[i[0].lower()]
    except KeyError:
      print 'Error: Overlaps with DFTB need further testing!'
      sys.exit(80)
    for j in range(nbs[i[0].lower()]):
      mapping[itot]=ii
      itot+=1
  #print mapping
  #sys.exit(81)
  # make interatomic overlap blocks unit matrices
  for i in range(Nb):
    ii=mapping[i]
    for j in range(Nb):
      jj=mapping[j]
      if ii!=jj:
        continue
      if i==j:
        Smat[i][j+Nb]=1.
        Smat[i+Nb][j]=1.
      else:
        Smat[i][j+Nb]=0.
        Smat[i+Nb][j]=0.
  return Smat

# =============================================================================================== #
# =============================================================================================== #
# ====================================== GAUSSIAN output parsing ================================ #
# =============================================================================================== #
# =============================================================================================== #

def getQMout(QMin):

    if PRINT:
        print '>>>>>>>>>>>>> Reading output files'
    starttime=datetime.datetime.now()

    QMout={}
    states=QMin['states']
    nstates=QMin['nstates']
    nmstates=QMin['nmstates']
    natom=QMin['natom']
    joblist=QMin['joblist']

    # TODO:
    # excited state energies and transition moments could be read from rwfdump "770R"
    # KS orbital energies: "522R"
    # geometry SEEMS TO BE in "507R"
    # 1TDM might be in "633R"


    # Hamiltonian
    if 'h' in QMin:     # or 'soc' in QMin:
        # make Hamiltonian
        if not 'h' in QMout:
            QMout['h']=makecmatrix(nmstates,nmstates)
        # go through all jobs
        for job in joblist:
            # first get energies from TAPE21
            logfile=os.path.join(QMin['scratchdir'],'master_%i/GAUSSIAN.log' % (job))
            energies=getenergy(logfile,job,QMin)
            ## also get SO matrix and mapping
            #if 'soc' in QMin and QMin['jobs'][job]['restr']:
                #outfile=os.path.join(QMin['scratchdir'],'master_%i/ADF.out' % (job))
                #submatrix,invstatemap=getsocm(outfile,t21file,job,QMin)
            mults=QMin['multmap'][-job]
            if 3 in mults:
              mults=[3]
            for i in range(nmstates):
                for j in range(nmstates):
                    m1,s1,ms1=tuple(QMin['statemap'][i+1])
                    m2,s2,ms2=tuple(QMin['statemap'][j+1])
                    if not m1 in mults or not m2 in mults:
                        continue
                    if i==j:
                        QMout['h'][i][j]=energies[(m1,s1)]
                    #elif 'soc' in QMin and QMin['jobs'][job]['restr']:
                        #if m1==m2==1:
                            #continue
                        #x=invstatemap[(m1,s1,ms1)]
                        #y=invstatemap[(m2,s2,ms2)]
                        #QMout['h'][i][j]=submatrix[x-1][y-1]

    # Dipole Moments
    if 'dm' in QMin:
        # make matrix
        if not 'dm' in QMout:
            QMout['dm']=[makecmatrix(nmstates,nmstates) for i in range(3)]
        # go through all jobs
        for job in joblist:
            logfile=os.path.join(QMin['scratchdir'],'master_%i/GAUSSIAN.log' % (job))
            dipoles=gettdm(logfile,job,QMin)
            mults=QMin['multmap'][-job]
            mults=QMin['multmap'][-job]
            if 3 in mults:
              mults=[3]
            for i in range(nmstates):
                m1,s1,ms1=tuple(QMin['statemap'][i+1])
                if not m1 in QMin['jobs'][job]['mults']:
                    continue
                for j in range(nmstates):
                    m2,s2,ms2=tuple(QMin['statemap'][j+1])
                    if not m2 in QMin['jobs'][job]['mults']:
                        continue
                    if i==j and (m1,s1) in QMin['gradmap']:
                        path,isgs=QMin['jobgrad'][(m1,s1)]
                        logfile=os.path.join(QMin['scratchdir'],path,'GAUSSIAN.log')
                        dm=getdm(logfile)
                        for ixyz in range(3):
                            QMout['dm'][ixyz][i][j]=dm[ixyz]
                    if i==j:
                        continue
                    if not m1==m2==mults[0] or not ms1==ms2:
                        continue
                    if s1==1:
                        for ixyz in range(3):
                            QMout['dm'][ixyz][i][j]=dipoles[(m2,s2)][ixyz]
                    elif s2==1:
                        for ixyz in range(3):
                            QMout['dm'][ixyz][i][j]=dipoles[(m1,s1)][ixyz]


    # Gradients
    if 'grad' in QMin:
        if not 'grad' in QMout:
            QMout['grad']=[ [ [ 0. for i in range(3) ] for j in range(natom) ] for k in range(nmstates) ]
        for grad in QMin['gradmap']:
            path,isgs=QMin['jobgrad'][grad]
            logfile=os.path.join(QMin['scratchdir'],path,'GAUSSIAN.log')
            g=getgrad(logfile,QMin)
            for istate in QMin['statemap']:
                state=QMin['statemap'][istate]
                if (state[0],state[1])==grad:
                    QMout['grad'][istate-1]=g
        if QMin['neglected_gradient']!='zero':
            for i in range(nmstates):
                m1,s1,ms1=tuple(QMin['statemap'][i+1])
                if not (m1,s1) in QMin['gradmap']:
                    if QMin['neglected_gradient']=='gs':
                        j=QMin['gsmap'][i+1]-1
                    elif QMin['neglected_gradient']=='closest':
                        e1=QMout['h'][i][i]
                        de=999.
                        for grad in QMin['gradmap']:
                            for k in range(nmstates):
                                m2,s2,ms2=tuple(QMin['statemap'][k+1])
                                if grad==(m2,s2):
                                    break
                            e2=QMout['h'][k][k]
                            if de>abs(e1-e2):
                                de=abs(e1-e2)
                                j=k
                    QMout['grad'][i]=QMout['grad'][j]

    # Regular Overlaps
    if 'overlap' in QMin:
        if not 'overlap' in QMout:
            QMout['overlap']=makecmatrix(nmstates,nmstates)
        for mult in itmult(QMin['states']):
            job=QMin['multmap'][mult]
            outfile=os.path.join(QMin['scratchdir'],'WFOVL_%i_%i/wfovl.out' % (mult,job))
            out=readfile(outfile)
            if PRINT:
                print 'Overlaps: '+shorten_DIR(outfile)
            for i in range(nmstates):
                for j in range(nmstates):
                    m1,s1,ms1=tuple(QMin['statemap'][i+1])
                    m2,s2,ms2=tuple(QMin['statemap'][j+1])
                    if not m1==m2==mult:
                        continue
                    if not ms1==ms2:
                        continue
                    QMout['overlap'][i][j]=getsmate(out,s1,s2)

    # Phases from overlaps
    if 'phases' in QMin:
        if not 'phases' in QMout:
            QMout['phases']=[ complex(1.,0.) for i in range(nmstates) ]
        if 'overlap' in QMout:
            for i in range(nmstates):
                if QMout['overlap'][i][i].real<0.:
                    QMout['phases'][i]=complex(-1.,0.)

    # Dyson norms
    if 'ion' in QMin:
        if not 'prop' in QMout:
            QMout['prop']=makecmatrix(nmstates,nmstates)
        for ion in QMin['ionmap']:
            outfile=os.path.join(QMin['scratchdir'],'Dyson_%i_%i_%i_%i/wfovl.out' % ion)
            out=readfile(outfile)
            if PRINT:
                print 'Dyson:    '+shorten_DIR(outfile)
            for i in range(nmstates):
                for j in range(nmstates):
                    m1,s1,ms1=tuple(QMin['statemap'][i+1])
                    m2,s2,ms2=tuple(QMin['statemap'][j+1])
                    if not (ion[0],ion[2])==(m1,m2) and not (ion[0],ion[2])==(m2,m1):
                        continue
                    if not abs(ms1-ms2)==0.5:
                        continue
                    # switch multiplicities such that m1 is smaller mult
                    if m1>m2:
                        s1,s2=s2,s1
                        m1,m2=m2,m1
                        ms1,ms2=ms2,ms1
                    # compute M_S overlap factor
                    if ms1<ms2:
                        factor=( ms1+1.+(m1-1.)/2. )/m1
                    else:
                        factor=( -ms1+1.+(m1-1.)/2. )/m1
                    QMout['prop'][i][j]=getDyson(out,s1,s2)*factor

    # TheoDORE
    if 'theodore' in QMin:
        if not 'theodore' in QMout:
            QMout['theodore']=makecmatrix(QMin['template']['theodore_n'],nmstates)
        for job in joblist:
            if not QMin['jobs'][job]['restr']:
                continue
            else:
                mults=QMin['jobs'][job]['mults']
                gsmult=mults[0]
                ns=0
                for i in mults:
                    ns+=QMin['states'][i-1]-(i==gsmult)
                if ns==0:
                    continue
            sumfile=os.path.join(QMin['scratchdir'],'master_%i/tden_summ.txt' % job)
            omffile=os.path.join(QMin['scratchdir'],'master_%i/OmFrag.txt' % job)
            props=get_theodore(sumfile,omffile,QMin)
            for i in range(nmstates):
                m1,s1,ms1=tuple(QMin['statemap'][i+1])
                if (m1,s1) in props:
                    for j in range(QMin['template']['theodore_n']):
                        QMout['theodore'][i][j]=props[(m1,s1)][j]

    ## QM/MM energy terms
    #if QMin['template']['qmmm']:
        #job=QMin['joblist'][0]
        #outfile=os.path.join(QMin['scratchdir'],'master_%i/ADF.out' % (job))
        #QMout['qmmm_energies']=get_qmmm_energies(outfile,QMin['template']['qmmm_coupling'])

    endtime=datetime.datetime.now()
    if PRINT:
        print "Readout Runtime: %s" % (endtime-starttime)

    if DEBUG:
        copydir=os.path.join(QMin['savedir'],'debug_GAUSSIAN_stdout')
        if not os.path.isdir(copydir):
            mkdir(copydir)
        for job in joblist:
            outfile=os.path.join(QMin['scratchdir'],'master_%i/GAUSSIAN.log' % (job))
            shutil.copy(outfile,os.path.join(copydir,"GAUSSIAN_%i.log" % job))
            if QMin['jobs'][job]['restr'] and 'theodore' in QMin:
                outfile=os.path.join(QMin['scratchdir'],'master_%i/tden_summ.txt' % job)
                try:
                    shutil.copy(outfile,os.path.join(copydir,'THEO_%i.out' % (job)))
                except IOError:
                    pass
                outfile=os.path.join(QMin['scratchdir'],'master_%i/OmFrag.txt' % job)
                try:
                    shutil.copy(outfile,os.path.join(copydir,'THEO_OMF_%i.out' % (job)))
                except IOError:
                    pass
        if 'grad' in QMin:
            for grad in QMin['gradmap']:
                path,isgs=QMin['jobgrad'][grad]
                outfile=os.path.join(QMin['scratchdir'],path,'GAUSSIAN.log')
                shutil.copy(outfile,os.path.join(copydir,"GAUSSIAN_GRAD_%i_%i.log" % grad))
        if 'overlap' in QMin:
            for mult in itmult(QMin['states']):
                job=QMin['multmap'][mult]
                outfile=os.path.join(QMin['scratchdir'],'WFOVL_%i_%i/wfovl.out' % (mult,job))
                shutil.copy(outfile,os.path.join(copydir,'WFOVL_%i_%i.out' % (mult,job)))
        if 'ion' in QMin:
            for ion in QMin['ionmap']:
                outfile=os.path.join(QMin['scratchdir'],'Dyson_%i_%i_%i_%i/wfovl.out' % ion)
                shutil.copy(outfile,os.path.join(copydir,'Dyson_%i_%i_%i_%i.out' % ion))

    return QMout

# ======================================================================= #
def getenergy(logfile,ijob,QMin):

    # open file
    f=readfile(logfile)
    if PRINT:
      print 'Energy:   '+shorten_DIR(logfile)

    # read ground state
    for line in f:
      if ' SCF Done:' in line:
        gsenergy=float(line.split()[4])

    # figure out the excited state settings
    mults=QMin['jobs'][ijob]['mults']
    restr=QMin['jobs'][ijob]['restr']
    gsmult=mults[0]
    estates_to_extract=deepcopy(QMin['states'])
    estates_to_extract[gsmult-1]-=1
    for imult in range(len(estates_to_extract)):
      if not imult+1 in mults:
        estates_to_extract[imult]=0
    for imult in range(len(estates_to_extract)):
      if imult+1 in mults:
        estates_to_extract[imult]=max(estates_to_extract)

    # extract excitation energies
    # loop also works if no energies should be extracted
    energies={(gsmult,1): gsenergy}
    for imult in mults:
      nstates=estates_to_extract[imult-1]
      if nstates>0:
        istate=0
        for line in f:
          if 'Excited State' in line:
            if restr:
              if not IToMult[imult] in line:
                continue
            energies[ (imult,istate+1+(gsmult==imult)) ] =float(line.split()[4])/au2eV+gsenergy
            istate+=1
            if istate>=nstates:
              break
    return energies

## ======================================================================= #
#def getsocm(outfile,t21file,ijob,QMin):

    ## read the ADF standard out into memory
    #out=readfile(outfile)
    #if PRINT:
        #print 'SOC:      '+shorten_DIR(outfile)

    ## open TAPE21
    #import kf
    #f=kf.kffile(t21file)
    #if PRINT:
        #print 'SOC:      '+shorten_DIR(t21file)


    ## get number of excitations from TAPE21
    #if QMin['ADFversion']>=(2017,208):
        #nrS=int(f.read('Excitations SS A','nr of excenergies')[0])
        #nrT=int(f.read('Excitations ST A','nr of excenergies')[0])
    #else:
        #nrS=int(f.read('All excitations','nr excitations')[0])
        #nrT=nrS
    #nrexci=nrS+3*nrT

    ## get GSCORR variable
    #GSCORR=False
    #for line in out:
        #if 'gscorr' in line.lower():
            #GSCORR=True
            #break
    #if GSCORR:
        #nrexci+=1

    #if QMin['ADFversion'] >= (2017,208):
        ## read SOC matrix from TAPE21
        #real_tri=f.read('Excitations SO A','SOmat-R')
        #imag_tri=f.read('Excitations SO A','SOmat-I')
        #real=[ [ 0+0j for i in range(nrexci) ] for j in range(nrexci) ]
        #x=0
        #y=0
        #for i in range(len(real_tri)):
            #if abs(real_tri[i])<1e-15:
                #real_tri[i]=0.
            #if abs(imag_tri[i])<1e-15:
                #imag_tri[i]=0.
            #real[x][y]=real_tri[i] + (0+1j)*imag_tri[i]
            #real[y][x]=real_tri[i] + (0-1j)*imag_tri[i]
            #x+=1
            #if x>y:
                #y+=1
                #x=0
    #else:
        ## read real matrix from stdout and make Hermitian
        #real=readSOC(out,'======  SO matrix real part',nrexci)
        #imag=readSOC(out,'======  SO matrix imaginary part',nrexci)
        #for x in range(len(real)):
            #for y in range(len(real[0])):
                #if x<y:
                    #real[x][y]+=(0+1j)*imag[x][y]
                #else:
                    #real[x][y]+=(0-1j)*imag[x][y]

    ## make statemap for the state ordering of the SO matrix
    #inv_statemap={}
    #if GSCORR:
        #inv_statemap[(1,1,0.0)]=nrexci
    #i=0
    #for x in range(nrS):
        #i+=1
        #inv_statemap[(1,x+2,0.0)] =i
    #for x in range(nrT):
        #i+=3
        #inv_statemap[(3,x+1, 0.0)]=i-2
        #inv_statemap[(3,x+1,+1.0)]=i-1
        #inv_statemap[(3,x+1,-1.0)]=i-0

    #return real,inv_statemap

## ======================================================================= #
#def readSOC(out,string,nrexci):

    ## find starting string
    #iline=-1
    #while True:
        #iline+=1
        #line=out[iline]
        #if string in line:
            #break

    ## read
    #matrix=[ [ 0. for i in range(nrexci) ] for j in range(nrexci) ]
    #for x in range(nrexci):
        #for y in range(nrexci):
            #if x>y:
                #x1,y1=y,x
            #else:
                #x1,y1=x,y
            #block=  x1/4
            #xoffset=x1%4+1
            #yoffset=block*(3+nrexci)+4+y1-block*(block+1)*2
            #matrix[x][y]=float(out[iline+yoffset].split()[xoffset])
            #if abs(matrix[x][y])<1e-15:
                #matrix[x][y]=0.

    #return matrix

# ======================================================================= #
def gettdm(logfile,ijob,QMin):

    # open file
    f=readfile(logfile)
    if PRINT:
      print 'Dipoles:  '+shorten_DIR(logfile)

    # figure out the excited state settings
    mults=QMin['jobs'][ijob]['mults']
    if 3 in mults:
      mults=[3]
    restr=QMin['jobs'][ijob]['restr']
    gsmult=mults[0]
    estates_to_extract=deepcopy(QMin['states'])
    estates_to_extract[gsmult-1]-=1
    for imult in range(len(estates_to_extract)):
      if not imult+1 in mults:
        estates_to_extract[imult]=0

    # get ordering of states in Gaussian output
    istate=[ int(i+1==gsmult) for i in range(len(QMin['states'])) ]
    index=0
    gaustatemap={}
    for iline,line in enumerate(f):
      if 'Excited State' in line:
        if restr:
          s=line.replace('-',' ').split()
          imult=IToMult[s[3]]
          istate[imult-1]+=1
          gaustatemap[(imult,istate[imult-1])]=index
          index+=1
        else:
          imult=gsmult
          istate[imult-1]+=1
          gaustatemap[(imult,istate[imult-1])]=index
          index+=1
    #print gaustatemap


    # extract transition dipole moments
    dipoles={}
    for imult in mults:
      if not imult==gsmult:
        continue
      nstates=estates_to_extract[imult-1]
      if nstates>0:
        for iline,line in enumerate(f):
          if 'Ground to excited state transition electric dipole moments ' in line:
            for istate in range(nstates):
              shift=gaustatemap[(imult,istate+1+(gsmult==imult))]
              s=f[iline+2+shift].split()
              dipoles[ (imult,istate+1+(gsmult==imult)) ]=[ float(i) for i in s[1:4]]
    #print dipoles
    return dipoles

# ======================================================================= #
def getdm(logfile):

    # open file
    f=readfile(logfile)
    if PRINT:
      print 'Dipoles:  '+shorten_DIR(logfile)

    for iline,line in enumerate(f):
        if 'Forces (Hartrees/Bohr)' in line:
            s=f[iline-2].split('=')[1].replace('D','E')
            dmx=float(s[0:15])
            dmy=float(s[15:30])
            dmz=float(s[30:45])
            dm=[dmx,dmy,dmz]
            return dm

# ======================================================================= #
def getgrad(logfile,QMin):

    # read file and check if ego is active
    out=readfile(logfile)
    if PRINT:
        print 'Gradient: '+shorten_DIR(logfile)

    # initialize
    natom=QMin['natom']
    g=[ [ 0. for i in range(3) ] for j in range(natom) ]

    # get gradient
    string='Forces (Hartrees/Bohr)'
    shift=3
    for iline,line in enumerate(out):
        if string in line:
            for iatom in range(natom):
                s=out[iline+shift+iatom].split()
                for i in range(3):
                    g[iatom][i]=-float(s[2+i])

    return g


## ======================================================================= #
#def get_qmmm_energies(outfile,coupling):

    #out=readfile(outfile)
    #if PRINT:
        #print 'QMMM:     '+shorten_DIR(outfile)

    #if coupling==1:
        #startstring='Q M / M M      E N E R G Y'
        #shift=2
    #elif coupling==2:
        #startstring='These results include the electrostatic interaction between QM and MM systems'
        #shift=0

    #toextract={'bond_mm':     (5,2),
               #'angle_mm':    (6,2),
               #'torsion_mm':  (7,2),
               #'VdW_mm':      (10,4),
               #'elstat_mm':   (11,2),
               #'VdW_qmmm':    (10,5),
               #'elstat_qmmm': (11,3)
              #}

    #iline=-1
    #while True:
        #iline+=1
        #if startstring in out[iline]:
            #break
    #iline+=shift

    #energies={}
    #for label in toextract:
        #t=toextract[label]
        #e=float(out[iline+t[0]].split()[t[1]])
        #energies[label]=e

    #return energies

# ======================================================================= #
def getsmate(out,s1,s2):
    ilines=-1
    while True:
        ilines+=1
        if ilines==len(out):
            print 'Overlap of states %i - %i not found!' % (s1,s2)
            sys.exit(82)
        if containsstring('Overlap matrix <PsiA_i|PsiB_j>', out[ilines]):
            break
    ilines+=1+s1
    f=out[ilines].split()
    return float(f[s2+1])

# ======================================================================= #
def getDyson(out,s1,s2):
    ilines=-1
    while True:
        ilines+=1
        if ilines==len(out):
            print 'Dyson norm of states %i - %i not found!' % (s1,s2)
            sys.exit(83)
        if containsstring('Dyson norm matrix <PsiA_i|PsiB_j>', out[ilines]):
            break
    ilines+=1+s1
    f=out[ilines].split()
    return float(f[s2+1])

# ======================================================================= #
def get_theodore(sumfile,omffile,QMin):  # TODO!
    out=readfile(sumfile)
    if PRINT:
        print 'TheoDORE: '+shorten_DIR(sumfile)
    props={}
    for line in out[2:]:
        s=line.replace('?',' ').split()
        if len(s)==0:
            continue
        n=int(re.search('([0-9]+)',s[0]).groups()[0])
        m=re.search('([a-zA-Z]+)',s[0]).groups()[0]
        for i in IToMult:
            if isinstance(i,str) and m in i:
                m=IToMult[i]
                break
        props[(m,n+(m==1))]=[ theo_float(i) for i in s[4:]]

    out=readfile(omffile)
    if PRINT:
        print 'TheoDORE: '+shorten_DIR(omffile)
    for line in out[1:]:
        s=line.replace('(',' ').replace(')',' ').split()
        if len(s)==0:
            continue
        n=int(re.search('([0-9]+)',s[0]).groups()[0])
        m=re.search('([a-zA-Z]+)',s[0]).groups()[0]
        for i in IToMult:
            if isinstance(i,str) and m in i:
                m=IToMult[i]
                break
        props[(m,n+(m==1))].extend([ theo_float(i) for i in s[2:]])

    return props

# ======================================================================= #
def theo_float(string):
    try:
        s=float(string)
    except ValueError:
        s=0.
    return s

# =============================================================================================== #
# =============================================================================================== #
# ========================================= Miscellaneous ======================================= #
# =============================================================================================== #
# =============================================================================================== #

# ======================================================================= #
def cleandir(directory):
    for data in os.listdir(directory):
        path=directory+'/'+data
        if os.path.isfile(path) or os.path.islink(path):
            if DEBUG:
                print 'rm %s' % (path)
            try:
                os.remove(path)
            except OSError:
                print 'Could not remove file from directory: %s' % (path)
        else:
            if DEBUG:
                print ''
            cleandir(path)
            os.rmdir(path)
            if DEBUG:
                print 'rm %s' % (path)
    if PRINT:
        print '===> Cleaning up directory %s' % (directory)

# ======================================================================= #
def backupdata(backupdir,QMin):
    # save all files in savedir, except which have 'old' in their name
    ls=os.listdir(QMin['savedir'])
    for f in ls:
        ff=QMin['savedir']+'/'+f
        if os.path.isfile(ff) and not 'old' in ff:
            step = int(QMin['step'][0])
            fdest=backupdir+'/'+f+'.stp'+str(step)
            shutil.copy(ff,fdest)

# =============================================================================================== #
# =============================================================================================== #
# ========================================= Main ================================================ #
# =============================================================================================== #
# =============================================================================================== #

def main():

    # Retrieve PRINT and DEBUG
    try:
        envPRINT=os.getenv('SH2Gau_PRINT')
        if envPRINT and envPRINT.lower()=='false':
            global PRINT
            PRINT=False
        envDEBUG=os.getenv('SH2Gau_DEBUG')
        if envDEBUG and envDEBUG.lower()=='true':
            global DEBUG
            DEBUG=True
    except ValueError:
        print 'PRINT or DEBUG environment variables do not evaluate to numerical values!'
        sys.exit(84)

    # Process Command line arguments
    if len(sys.argv)!=2:
        print 'Usage:\n./SHARC_GAUSSIAN.py <QMin>\n'
        print 'version:',version
        print 'date:',versiondate
        print 'changelog:\n',changelogstring
        sys.exit(85)
    QMinfilename=sys.argv[1]

    # Print header
    printheader()

    # Read QMinfile
    QMin=readQMin(QMinfilename)

    # get the job schedule
    QMin,schedule=generate_joblist(QMin)
    printQMin(QMin)
    if DEBUG:
        pprint.pprint(schedule,depth=1)

    # run all the ADF jobs
    errorcodes=runjobs(schedule,QMin)

    ## do all necessary overlap and Dyson calculations
    errorcodes=run_wfoverlap(QMin,errorcodes)

    ## do all necessary Theodore calculations
    errorcodes=run_theodore(QMin,errorcodes)

      # read all the output files
    QMout=getQMout(QMin)
    if PRINT or DEBUG:
        printQMout(QMin,QMout)

    # backup data if requested
    if 'backup' in QMin:
        backupdata(QMin['backup'],QMin)

    # Measure time
    runtime=measuretime()
    QMout['runtime']=runtime

    # Write QMout
    writeQMout(QMin,QMout,QMinfilename)

    # Remove Scratchfiles from SCRATCHDIR
    if not DEBUG:
        cleandir(QMin['scratchdir'])
        if 'cleanup' in QMin:
            cleandir(QMin['savedir'])

    print
    print datetime.datetime.now()
    print '#================ END ================#'

if __name__ == '__main__':
    main()






# kate: indent-width 2
