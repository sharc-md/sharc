#!/usr/bin/env python2

#******************************************
#
#    SHARC Program Suite
#
#    Copyright (c) 2018 University of Vienna
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

# This script calculates QC results for a model system
# 
# Reads QM.in
# Calculates SOC matrix, dipole moments, gradients, nacs and overlaps
# Writes these back to QM.out

from copy import deepcopy 
import math
import sys
import re
import os
import stat
import shutil
import datetime
import pprint
try:
  import numpy
  NONUMPY=False
except ImportError:
  import subprocess as sp
  NONUMPY=True


# =========================================================
# compatibility stuff

if sys.version_info[0]!=2:
  print 'This is a script for Python 2!'
  sys.exit(0)

if sys.version_info[1]<4:
  print 'INFO: Script is not tested for Python <2.4! Proceed at own risk!'

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

# =========================================================
# some constants
PRINT = True
DEBUG = False
CM_TO_HARTREE = 1./219474.6
HARTREE_TO_EV = 27.211396132
U_TO_AMU = 1./5.4857990943e-4
BOHR_TO_ANG=0.529177211
PI = math.pi

version='2.0'
versiondate=datetime.date(2018,2,1)


# hash table for conversion of multiplicity to the keywords used in MOLPRO
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

# hash table for conversion of polarisations to the keywords used in MOLPRO
IToPol={
        0: 'X', 
        1: 'Y', 
        2: 'Z', 
        'X': 0, 
        'Y': 1, 
        'Z': 2
        }

# =========================================================
# =========================================================
# =========================================================

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
  3 integer: MS value'''

  for i in range(len(states)):
    if states[i]<1:
      continue
    for k in range(i+1):
      for j in range(states[i]):
        yield i+1,j+1,k-i/2.
  return

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
    string+='%i\t%s\t' % (atom+1,geo[atom][0])
    for xyz in range(3):
      if grad[atom][xyz]!=0:
        iszero=False
      string+='% .5f\t' % (grad[atom][xyz])
    string+='\n'
  if iszero:
    print '\t\t...is identical zero...\n'
  else:
    print string

# ======================================================================= #

def printQMout(QMin,QMout):
  '''If PRINT, prints a summary of all requested QM output values. Matrices are formatted using printcomplexmatrix, vectors using printgrad. 

  Arguments:
  1 dictionary: QMin
  2 dictionary: QMout'''

  if DEBUG:
    pprint.pprint(QMout)
  if not PRINT:
    return
  states=QMin['states']
  nstates=QMin['nstates']
  nmstates=QMin['nmstates']
  natom=QMin['natom']
  print '===> Results:\n'
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
      printgrad(QMout['grad'][istate],natom,QMin['geom'])
      istate+=1
  # Non-adiabatic couplings
  if 'nacdt' in QMin:
    print '=> Numerical Non-adiabatic couplings:\n'
    matrix=QMout['nacdt']
    printcomplexmatrix(matrix,states)
    matrix=deepcopy(QMout['mrcioverlap'])
    for i in range(nmstates):
      for j in range(nmstates):
        matrix[i][j]=complex(matrix[i][j])
    print '=> MRCI overlaps:\n'
    printcomplexmatrix(matrix,states)
    if 'phases' in QMout:
      print '=> Wavefunction Phases:\n%i\n' % (nmstates)
      for i in range(nmstates):
        print '% 3.1f % 3.1f' % (QMout['phases'][i].real,QMout['phases'][i].imag)
      print '\n'
  if 'nacdr' in QMin:
      print '=> Analytical Non-adiabatic coupling vectors:\n'
      istate=0
      for imult,i,msi in itnmstates(states):
        jstate=0
        for jmult,j,msj in itnmstates(states):
          if imult==jmult and msi==msj:
            print '%s\tStates %i - %i\tMs= % .1f:' % (IToMult[imult],i,j,msi)
            printgrad(QMout['nacdr'][istate][jstate],natom,QMin['geom'])
          jstate+=1
        istate+=1
  if 'dmdr' in QMin:
      print '=> Dipole moment derivative vectors:\n'
      istate=0
      for imult,i,msi in itnmstates(states):
        jstate=0
        for jmult,j,msj in itnmstates(states):
          if imult==jmult and msi==msj:
            for ipol in range(3):
              print '%s\tStates %i - %i\tMs= % .1f\tPolarization %s:' % (IToMult[imult],i,j,msi,IToPol[ipol])
              printgrad(QMout['dmdr'][istate][jstate][ipol],natom,QMin['geom'])
          jstate+=1
        istate+=1
  if 'overlap' in QMin:
    print '=> Overlap matrix:\n'
    matrix=QMout['overlap']
    printcomplexmatrix(matrix,states)
    if 'phases' in QMout:
      print '=> Wavefunction Phases:\n%i\n' % (nmstates)
      for i in range(nmstates):
        print '% 3.1f % 3.1f' % (QMout['phases'][i].real,QMout['phases'][i].imag)
      print '\n'
  # Angular momentum matrices
  if 'angular' in QMin:
    print '=> Angular Momentum Matrices:\n'
    for xyz in range(3):
      print 'Polarisation %s:' % (IToPol[xyz])
      matrix=QMout['angular'][xyz]
      printcomplexmatrix(matrix,states)
  sys.stdout.flush()

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
            sys.exit(12)
    else:
        try:
            os.makedirs(SCRATCHDIR)
        except OSError:
            print 'Can not create SCRATCHDIR=%s\n' % (SCRATCHDIR)
            sys.exit(13)

# =========================================================
# =========================================================
# =========================================================
# diagonalization and transformation stuff
class diagonalizer:
  def __init__(self):
    exe=os.getenv('SHARC')
    exe=os.path.expanduser(os.path.expandvars(exe))+'/diagonalizer.x'
    if not os.path.isfile(exe):
      print 'SHARC auxilliary diagonalizer not found at %s!' % (exe)
      sys.exit(14)
    self.exe=exe
  def eigh(self,H):
    STDIN='C %i %i\nTitle\n' % (len(H),len(H))
    for x in H:
      for y in x:
        STDIN+='%20.13f %20.13f ' % (y.real,y.imag)
      STDIN+='\n'
    proc=sp.Popen(self.exe,stdin=sp.PIPE,stdout=sp.PIPE)
    STDOUT=proc.communicate(input=STDIN)[0].split('\n')
    shift=1
    for ix in range(len(H)):
      line=STDOUT[shift+ix].split()
      for iy in range(len(H)):
        H[ix][iy]=complex(float(line[0+2*iy]),float(line[1+2*iy]))
    U=[ [ 0. for i in range(len(H)) ] for j in range(len(H)) ]
    shift=2+len(H)
    for ix in range(len(H)):
      line=STDOUT[shift+ix].split()
      for iy in range(len(H)):
        U[ix][iy]=complex(float(line[0+2*iy]),float(line[1+2*iy]))
    return H,U

# =========================================================
def diagonalize(A):
  if NONUMPY:
    diagon=diagonalizer()
  else:
    diagon=numpy.linalg

  # diagonalize Hamiltonian
  Hd,U=diagon.eigh(A)
  if not NONUMPY:
    Hd=numpy.diag(Hd)

  return Hd,U



# =========================================================
def transform(A,U):
  '''returns U^T.A.U'''

  if NONUMPY:
    temp=[ [ 0. for i in range(len(U)) ] for j in range(len(U)) ]
    B=[ [ 0. for i in range(len(U)) ] for j in range(len(U)) ]
    for a in range(len(U)):
      for b in range(len(U)):
        for i in range(len(U)):
          temp[a][b]+=U[i][a].conjugate()*A[i][b]
    for a in range(len(U)):
      for b in range(len(U)):
        for i in range(len(U)):
          B[a][b]+=temp[a][i]*U[i][b]
  else:
    Ucon=[ [ U[i][j].conjugate() for i in range(len(U)) ] for j in range(len(U)) ]
    B=numpy.dot(Ucon,numpy.dot(A,U))

  return B

# =========================================================
def matmult(A,B,transA=False,transB=False):
  if transA:
    tempA=[ [ A[i][j].conjugate() for i in range(len(A)) ] for j in range(len(A)) ]
  else:
    tempA=A
  if transB:
    tempB=[ [ B[i][j].conjugate() for i in range(len(B)) ] for j in range(len(B)) ]
  else:
    tempB=B
  if NONUMPY:
    Res=[ [ 0. for i in range(len(A)) ] for j in range(len(A)) ]
    for a in range(len(A)):
      for b in range(len(A)):
        for i in range(len(A)):
          Res[a][b]+=tempA[a][i]*tempB[i][b]
  else:
    Res=numpy.dot(tempA,tempB)
  return Res

# =============================================================================================== #
# =============================================================================================== #
# =========================================== QMout writing ===================================== #
# =============================================================================================== #
# =============================================================================================== #


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
def writeQMout(QMin,QMout,QMinfilename):
  '''Writes the requested quantities to the file which SHARC reads in. The filename is QMinfilename with everything after the first dot replaced by "out". 

  Arguments:
  1 dictionary: QMin
  2 dictionary: QMout
  3 string: QMinfilename'''

  outfilename=os.path.join(QMin['savedir'],'geom.out')
  string=''
  for iatom in range(QMin['natom']):
    string+=QMin['geom'][iatom][0]
    for i in range(3):
      string+=' %12.9f ' % (QMin['geom'][iatom][i+1])
    string+='\n'
  string+='\n'
  try:
    outfile=open(outfilename,'w')
    outfile.write(string)
    outfile.close()
  except IOError:
    print 'Could not write to savedir!'
    sys.exit(15)

  os.chdir(QMin['pwd'])
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
  if 'angular' in QMin:
    string+=writeQMoutang(QMin,QMout)
  if 'grad' in QMin:
    string+=writeQMoutgrad(QMin,QMout)
  if 'nacdr' in QMin:
    string+=writeQMoutnacana(QMin,QMout)
  if 'overlap' in QMin:
    string+=writeQMoutnacsmat(QMin,QMout)
  if 'ion' in QMin:
    string+=writeQMoutprop(QMin,QMout)
  if 'dmdr' in QMin:
    string+=writeQMoutDMgrad(QMin,QMout)
  if 'phases' in QMin:
    string+=writeQmoutPhases(QMin,QMout)
  string+=writeQMouttime(QMin,QMout)
  try:
    outfile=open(outfilename,'w')
    outfile.write(string)
    outfile.close()
  except IOError:
    print 'Could not write QM output!'
    sys.exit(16)
  if 'backup' in QMin:
    try:
      outfile=open(QMin['backup']+'/'+outfilename,'w')
      outfile.write(string)
      outfile.close()
    except IOError:
      print 'WARNING: Could not write QM output backup!'
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
      string+='%s %s ' % (eformat(QMout['h'][i][j].real,12,3),eformat(QMout['h'][i][j].imag,12,3))
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
        string+='%s %s ' % (eformat(QMout['dm'][xyz][i][j].real,12,3),eformat(QMout['dm'][xyz][i][j].imag,12,3))
      string+='\n'
    string+=''
  return string

# ======================================================================= #
def writeQMoutang(QMin,QMout):
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
  string+='! %i Angular Momentum Matrices (3x%ix%i, complex)\n' % (9,nmstates,nmstates)
  for xyz in range(3):
    string+='%i %i\n' % (nmstates,nmstates)
    for i in range(nmstates):
      for j in range(nmstates):
        string+='%s %s ' % (eformat(QMout['angular'][xyz][i][j].real,12,3),eformat(QMout['angular'][xyz][i][j].imag,12,3))
      string+='\n'
    string+=''
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
    string+='%i %i ! m1 %i s1 %i ms1 %i\n' % (natom,3,imult,istate,ims)
    for atom in range(natom):
      for xyz in range(3):
        string+='%s ' % (eformat(QMout['grad'][i][atom][xyz],12,3))
      string+='\n'
    string+=''
    i+=1
  return string

# ======================================================================= #
def writeQMoutnacnum(QMin,QMout):
  '''Generates a string with the NAC matrix in SHARC format.

  The string starts with a ! followed by a flag specifying the type of data. In the next line, the dimensions of the matrix are given, followed by nmstates blocks of nmstates elements. Blocks are separated by a blank line.

  Arguments:
  1 dictionary: QMin
  2 dictionary: QMout

  Returns:
  1 string: multiline string with the NAC matrix'''

  states=QMin['states']
  nstates=QMin['nstates']
  nmstates=QMin['nmstates']
  natom=QMin['natom']
  string=''
  string+='! %i Non-adiabatic couplings (ddt) (%ix%i, complex)\n' % (4,nmstates,nmstates)
  string+='%i %i\n' % (nmstates,nmstates)
  for i in range(nmstates):
    for j in range(nmstates):
      string+='%s %s ' % (eformat(QMout['nacdt'][i][j].real,12,3),eformat(QMout['nacdt'][i][j].imag,12,3))
    string+='\n'
  string+=''
  # also write wavefunction phases
  string+='! %i Wavefunction phases (%i, complex)\n' % (7,nmstates)
  for i in range(nmstates):
    string+='%s %s\n' % (eformat(QMout['phases'][i],12,3),eformat(0.,12,3))
  string+='\n\n'
  return string

# ======================================================================= #
def writeQMoutnacana(QMin,QMout):
  '''Generates a string with the NAC vectors in SHARC format.

  The string starts with a ! followed by a flag specifying the type of data. On the next line, natom and 3 are written, followed by the gradient, with one line per atom and a blank line at the end. Each MS component shows up (nmstates x nmstates vectors are written).

  Arguments:
  1 dictionary: QMin
  2 dictionary: QMout

  Returns:
  1 string: multiline string with the NAC vectors'''

  states=QMin['states']
  nstates=QMin['nstates']
  nmstates=QMin['nmstates']
  natom=QMin['natom']
  string=''
  string+='! %i Non-adiabatic couplings (ddr) (%ix%ix%ix3, real)\n' % (5,nmstates,nmstates,natom)
  i=0
  for imult,istate,ims in itnmstates(states):
    j=0
    for jmult,jstate,jms in itnmstates(states):
      string+='%i %i ! m1 %i s1 %i ms1 %i   m2 %i s2 %i ms2 %i\n' % (natom,3,imult,istate,ims,jmult,jstate,jms)
      for atom in range(natom):
        for xyz in range(3):
          string+='%s ' % (eformat(QMout['nacdr'][i][j][atom][xyz],12,3))
        string+='\n'
      string+=''
      j+=1
    i+=1
  return string

# ======================================================================= #
def writeQMoutDMgrad(QMin,QMout):

  states=QMin['states']
  nstates=QMin['nstates']
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
            string+='%s ' % (eformat(QMout['dmdr'][i][j][ipol][atom][xyz],12,3))
          string+='\n'
        string+=''
      j+=1
    i+=1
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
      string+='%s %s ' % (eformat(QMout['overlap'][j][i].real,12,3),eformat(QMout['overlap'][j][i].imag,12,3))
    string+='\n'
  string+='\n'
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

# ======================================================================= #
def writeQMoutprop(QMin,QMout):
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
  string+='! %i Property Matrix (%ix%i, complex)\n' % (11,nmstates,nmstates)
  string+='%i %i\n' % (nmstates,nmstates)
  for i in range(nmstates):
    for j in range(nmstates):
      string+='%s %s ' % (eformat(QMout['prop'][i][j].real,12,3),eformat(QMout['prop'][i][j].imag,12,3))
    string+='\n'
  string+='\n'
  return string

# ======================================================================= #
def writeQmoutPhases(QMin,QMout):

  string='! 7 Phases\n%i ! for all nmstates\n' % (QMin['nmstates'])
  for i in range(QMin['nmstates']):
    string+='%s %s\n' % (eformat(QMout['phases'][i].real,9,3),eformat(QMout['phases'][i].imag,9,3))
  return string

# =========================================================
# =========================================================
# =========================================================

# =========================================================
# class for real function of three coordinates
# includes derivatives
#
# using eval is probably risky, because it can execute any code...
class func_mat:
  def split_strings(self,n,strings):
    a=[]
    for i in range(n):
      s=strings[i].strip().split(',')
      if any([j.strip()=='' for j in s[0:i+1]]):
        print 'Matrix elements missing in definition!'
        sys.exit(17)
      a.append(s)
    return a

  def __init__(self, rstrings, istrings=None):
    # rstring: list of strings defining the matrix elements for the real part
    # istrings: for the imaginary part
    # both matrices should be defined with lower triangle matrices
    self.n=int(len(rstrings))
    self.cmpx=False
    self.r=self.split_strings(self.n,rstrings)

    if istrings!=None:
      self.cmpx=True
      self.i=self.split_strings(self.n,istrings)

  def mat(self,_geom,_var):
    # geom is list of list of atom coordinates
    # e.g. [['I', 0.6, 0.0, 0.0], ['Br', 2.4, 0.0, 0.0]]
    # var is dictionary of variable mappings
    # e.g. {'y': (1, 0), 'x': (0, 0)}

    # set the variables
    for _v in _var:
      if isinstance(_var[_v],list):
        _i,_j=tuple(_var[_v])
        exec('%s=%f' % (_v,_geom[_i][_j+1]) )
      else:
        exec('%s=%f' % (_v,_var[_v] ) )

    # build the real matrix
    _R=[ [ 0. for _i in range(self.n) ] for _j in range(self.n) ]
    for _i in range(self.n):
      for _j in range(_i+1):
        _R[_i][_j]=eval(self.r[_i][_j])

    # build the imaginary matrix
    _I=[ [ 0. for _i in range(self.n) ] for _j in range(self.n) ]
    if self.cmpx:
      for _i in range(self.n):
        for _j in range(_i+1):
          _I[_i][_j]=eval(self.i[_i][_j])

    # build the Hermitian matrix
    _M=[ [ 0. for _i in range(self.n) ] for _j in range(self.n) ]
    for _i in range(self.n):
      for _j in range(self.n):
        if _i<_j:
          _M[_i][_j]=complex(_R[_j][_i],-_I[_j][_i])
        elif _i==_j:
          _M[_i][_j]=complex(_R[_i][_j],0)
        elif _i>_j:
          _M[_i][_j]=complex(_R[_i][_j],_I[_i][_j])

    # return
    return _M

# =========================================================
def read_QMin():
  # reads the geometry, unit keyword, nstates keyword
  # does not read the request keywords, since it calculates by default all quantities
  QMin={}
  f=open('QM.in')
  qmin=f.readlines()
  f.close()

  QMin['natom']=int(qmin[0])
  QMin['comment']=qmin[1]

  # get geometry
  line=qmin[2].split()[1:4]
  geom=[ [ float(line[i]) for i in range(3) ] ]

  geom=[]
  for i in range(2,QMin['natom']+2):
    line=qmin[i].split()
    for j in range(3):
      line[j+1]=float(line[j+1])
    geom.append(line)

  # find states keyword
  for line in qmin:
    s=line.split()
    if len(s)==0:
      continue
    if 'states' in s[0].lower():
      states=[]
      for iatom in range(len(s)-1):
        states.append(int(s[iatom+1]))
      break
  else:
    print 'No state keyword given!'
    sys.exit(18)
  nstates=0
  nmstates=0
  for mult,i in enumerate(states):
    nstates+=i
    nmstates+=(mult+1)*i
  QMin['states']=states
  QMin['nstates']=nstates
  QMin['nmstates']=nmstates
  statemap={}
  i=1
  for imult,nstates in enumerate(states):
    if nstates==0:
      continue
    for ims in range(imult+1):
      ms=ims-imult/2.
      for istate in range(nstates):
        statemap[i]=[imult+1,istate+1,ms]
        i+=1
  QMin['statemap']=statemap

  # find unit keyword
  factor=1.
  for line in qmin:
    s=line.split()
    if len(s)==0:
      continue
    if 'unit' in s[0].lower():
      if not 'bohr' in s[1].lower():
        factor=BOHR_TO_ANG
  for i in range(QMin['natom']):
    for j in range(3):
      geom[i][j+1]/=factor

  # find init, samestep, restart
  for line in qmin:
    line=line.split('#')[0]
    s=line.split()
    if len(s)==0:
      continue
    if 'init' in s[0].lower():
      QMin['init']=[]
    if 'samestep' in s[0].lower():
      QMin['samestep']=[]
    if 'restart' in s[0].lower():
      QMin['restart']=[]

  # find savedir
  QMin['savedir']='./SAVEDIR/'
  for line in qmin:
    s=line.split()
    if len(s)==0:
      continue
    if 'savedir' in s[0].lower():
      QMin['savedir']=s[1]
  QMin['savedir']=os.path.abspath(os.path.expanduser(os.path.expandvars(QMin['savedir'])))

  if 'init' in QMin:
    checkscratch(QMin['savedir'])
  if not 'init' in QMin and not 'samestep' in QMin and not 'restart' in QMin:
    fromfile=os.path.join(QMin['savedir'],'geom.out')
    tofile=os.path.join(QMin['savedir'],'geom_old.out')
    shutil.copy(fromfile,tofile)

  # read old geometry from savedir/geom_old.out
  try:
    filename=os.path.join(QMin['savedir'],'geom_old.out')
    f=open(filename)
    #f=open('QM.out')
    qmin=f.readlines()
    f.close()

    geomold=[]
    for i in range(QMin['natom']):
      line=qmin[i].split()
      for j in range(3):
        line[j+1]=float(line[j+1])
      geomold.append(line)
  except IOError:
    geomold=geom

  QMin['geom']=geom
  QMin['geomold']=geomold

  # find forbidden keywords and optional keywords
  for line in qmin:
    s=line.lower().split()
    if len(s)==0:
      continue
    if 'nacdt' in s[0] or 'nacdr' in s[0]:
      print 'NACDR and NACDT are not supported!'
      sys.exit(19)
    if 'dmdr' in s[0]:
      QMin['dmdr']=[]


  # add request keywords
  QMin['soc']=[]
  QMin['dm']=[]
  QMin['grad']=[]
  QMin['overlap']=[]
  QMin['phases']=[]
  QMin['pwd']=os.getcwd()
  return QMin

# =========================================================
def find_lines(nlines,match,strings):
  smatch=match.lower().split()
  iline=-1
  while True:
    iline+=1
    if iline==len(strings):
      return []
    line=strings[iline].lower().split()
    if tuple(line)==tuple(smatch):
      return strings[iline+1:iline+nlines+1]


# =========================================================
def read_SH2Ana(QMin):
  # reads SH2Ana.inp, deletes comments and blank lines
  SH2ANA={}
  if os.path.isfile('Analytical.template'):
    f=open('Analytical.template')
  else:
    f=open('SH2Ana.inp')
  sh2ana=f.readlines()
  f.close()

  # check natoms
  natom=int(sh2ana[0])
  if not natom==QMin['natom']:
    print 'Natom from QM.in and from SH2Ana.inp are inconsistent!'
    sys.exit(20)

  # check nstates
  nstates=int(sh2ana[1])
  if not nstates==QMin['nmstates']:
    print 'NMstates from QM.in and nstates from SH2Ana.inp are inconsistent!'
    sys.exit(21)

  # read the coordinate <-> variable mapping
  gvar=[]
  for i in range(natom):
    s=sh2ana[i+2].lower().split()
    if s[0]!=QMin['geom'][i][0].lower():
      print 'Inconsistent atom labels in QM.in and SH2Ana.inp!'
      sys.exit(22)
    gvar.append(s[1:])
  var={}
  for i in range(natom):
    for j in range(3):
      v=gvar[i][j]
      if v=='0':
        continue
      if v[0:1]=='_':
        print 'Variable names must not start with an underscore!'
        sys.exit(23)
      if v in var:
        print 'Repeated variable in geom<->var mapping in SH2Ana.inp!'
        sys.exit(24)
      var[v]=[i,j]

  # read additional variables
  iline=-1
  while True:
    iline+=1
    if iline==len(sh2ana):
      break
    line=re.sub('#.*$','',sh2ana[iline])
    s=line.lower().split()
    if s==[]:
      continue
    if 'variables' in s[0]:
      while True:
        iline+=1
        line=re.sub('#.*$','',sh2ana[iline])
        s=line.split()
        if s==[]:
          continue
        if 'end' in s[0].lower():
          break
        if s[0][0:1]=='_':
          print 'Variable names must not start with an underscore!'
          sys.exit(25)
        if s[0] in var:
          print 'Repeated variable in additional variables in SH2Ana.inp!'
          sys.exit(26)
        var[s[0]]=float(s[1])

  SH2ANA['var']=var
  SH2ANA['gvar']=gvar


  # look out for the nodiag keyword
  iline=-1
  while True:
    iline+=1
    if iline==len(sh2ana):
      break
    line=re.sub('#.*$','',sh2ana[iline])
    s=line.lower().split()
    if 'nodiag' in s:
      QMin['nodiag']=[]


  # obtain the Hamiltonian
  Hstring=find_lines(nstates,'Hamiltonian',sh2ana)
  if Hstring==[]:
    print 'No Hamiltonian defined in SH2Ana.inp!'
    sys.exit(27)
  fmat=func_mat(Hstring)
  SH2ANA['H']=fmat.mat(QMin['geom'],SH2ANA['var'])

  # obtain the old Hamiltonian
  SH2ANA['Hold']=fmat.mat(QMin['geomold'],SH2ANA['var'])

  # obtain the derivatives
  SH2ANA['deriv']={}
  for v in var:
    if isinstance(var[v],list):
      Dstring=find_lines(nstates,'Derivatives %s' % (v),sh2ana)
      if Dstring==[]:
        print 'No derivative matrix for variable %s defined in SH2Ana.inp!' % (v)
        sys.exit(28)
      fmat=func_mat(Dstring)
      SH2ANA['deriv'][v]=fmat.mat(QMin['geom'],SH2ANA['var'])

  # obtain the dipole matrices
  SH2ANA['dipole']={}
  for idir in range(1,4):
    Dstring=find_lines(nstates,'Dipole %s' % (idir),sh2ana)
    if Dstring==[]:
      SH2ANA['dipole'][idir]=[ [ complex(0.,0.) for i in range(nstates) ] for j in range(nstates) ]
    else:
      fmat=func_mat(Dstring)
      SH2ANA['dipole'][idir]=fmat.mat(QMin['geom'],SH2ANA['var'])

  # obtain the dipole derivative matrices
  SH2ANA['dipolederiv']={}
  for idir in range(1,4):
    SH2ANA['dipolederiv'][idir]={}
    for v in var:
      if isinstance(var[v],list):
        Dstring=find_lines(nstates,'Dipolederivatives %s %s' % (idir,v),sh2ana)
        if Dstring==[]:
          SH2ANA['dipolederiv'][idir][v]=[ [ complex(0.,0.) for i in range(nstates) ] for j in range(nstates) ]
        else:
          fmat=func_mat(Dstring)
          SH2ANA['dipolederiv'][idir][v]=fmat.mat(QMin['geom'],SH2ANA['var'])

  # obtain the SO matrix
  Rstring=find_lines(nstates,'SpinOrbit R',sh2ana)
  if Rstring==[]:
    SH2ANA['soc']=[ [ complex(0.,0.) for i in range(nstates) ] for j in range(nstates) ]
  else:
    Istring=find_lines(nstates,'SpinOrbit I',sh2ana)
    if Istring==[]:
      fmat=func_mat(Rstring)
    else:
      fmat=func_mat(Rstring,Istring)
    SH2ANA['soc']=fmat.mat(QMin['geom'],SH2ANA['var'])

  return SH2ANA, QMin


# ============================================================================
# ============================================================================
# ============================================================================
def getQMout(QMin,SH2ANA):
  '''Calculates the MCH Hamiltonian, SOC matrix ,overlap matrix, gradients, DM'''

  QMout={}

  #pprint.pprint(SH2ANA,width=192)

  if NONUMPY:
    diagon=diagonalizer()
  else:
    diagon=numpy.linalg

  # diagonalize Hamiltonian
  if not 'nodiag' in QMin:
    Hd,U=diagonalize(SH2ANA['H'])
  else:
    Hd=SH2ANA['H']
    U=[ [ i==j for i in range(QMin['nmstates']) ] for j in range(QMin['nmstates']) ]

  #pprint.pprint(SH2ANA['H'])
  #pprint.pprint(U)

  # diagonalize old Hamiltonian
  if not 'nodiag' in QMin:
    Uold=diagonalize(SH2ANA['Hold'])[1]
  else:
    Uold=[ [ i==j for i in range(QMin['nmstates']) ] for j in range(QMin['nmstates']) ]

  # get gradients
  ggrad=[]
  for iatom in range(QMin['natom']):
    ggrad.append([])
    for idir in range(3):
      if SH2ANA['gvar'][iatom][idir]=='0':
        ggrad[-1].append( [ 0. for i in range(QMin['nmstates']) ] )
      else:
        v=SH2ANA['gvar'][iatom][idir]
        Dmatrix=transform(SH2ANA['deriv'][v],U)
        ggrad[-1].append( [ Dmatrix[i][i].real for i in range(QMin['nmstates']) ] )
  # rearrange gradients 
  grad=[ [ [ ggrad[iatom][idir][istate] for idir in range(3) ] for iatom in range(QMin['natom']) ] for istate in range(QMin['nmstates']) ]

  # transform dipole matrices
  dipole=[]
  for idir in range(3):
    Dmatrix=transform(SH2ANA['dipole'][idir+1],U)
    dipole.append(Dmatrix)

  # get dipole derivatives
  gdmdr=[]
  for ipol in range(1,4):
    gdmdr.append([])
    for iatom in range(QMin['natom']):
      gdmdr[-1].append([])
      for idir in range(3):
        if SH2ANA['gvar'][iatom][idir]=='0':
          gdmdr[-1][-1].append( [ [ 0. for i in range(QMin['nmstates']) ] for j in range(QMin['nmstates']) ] )
        else:
          v=SH2ANA['gvar'][iatom][idir]
          Dmatrix=transform(SH2ANA['dipolederiv'][ipol][v],U)
          gdmdr[-1][-1].append( [ [ Dmatrix[i][j].real for i in range(QMin['nmstates']) ] for j in range(QMin['nmstates']) ] )
  # rearrange dipole gradients 
  #pprint.pprint(gdmdr,width=192)
  dmdr= [ 
          [ 
            [ 
              [ 
                [ 
                  gdmdr[ipol][iatom][idir][istate][jstate] 
                for idir in range(3) ] 
              for iatom in range(QMin['natom']) ] 
            for ipol in range(3) ] 
          for istate in range(QMin['nmstates']) ] 
        for jstate in range(QMin['nmstates']) ]

  # get overlap matrix
  overlap=matmult(Uold,U,transA=True)

  # transform SOC matrix
  SO=transform(SH2ANA['soc'],U)
  for i in range(QMin['nmstates']):
    SO[i][i]=complex(0.,0.)
  Hfull=[ [ Hd[i][j]+SO[i][j] for i in range(QMin['nmstates']) ] for j in range(QMin['nmstates']) ]


  # assign QMout elements
  QMout['h']=Hfull
  QMout['dm']=dipole
  QMout['grad']=grad
  QMout['dmdr']=dmdr
  QMout['overlap']=overlap
  QMout['runtime']=0.
  # Phases from overlaps
  if 'phases' in QMin:
    if not 'phases' in QMout:
      QMout['phases']=[ complex(1.,0.) for i in range(QMin['nmstates']) ]
    if 'overlap' in QMout:
      for i in range(QMin['nmstates']):
        if QMout['overlap'][i][i].real<0.:
          QMout['phases'][i]=complex(-1.,0.)

  #pprint.pprint(QMout,width=192)






  return QMout

# ============================================================================
def main():

  QMin=read_QMin()
  SH2ANA,QMin=read_SH2Ana(QMin)
  #pprint.pprint( QMin)
  #pprint.pprint( SH2ANA)

  QMout=getQMout(QMin,SH2ANA)

  printQMout(QMin,QMout)

  # Write QMout
  writeQMout(QMin,QMout,'QM.in')

  print '#================ END ================#'

# ============================================================================
if __name__ == '__main__':
  try:
    main()
  except KeyboardInterrupt:
    print '\nCtrl+C makes me a sad SHARC ;-(\n'
    sys.exit(0)
