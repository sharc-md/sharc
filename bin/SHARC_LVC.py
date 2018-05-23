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

# This script calculates QC results for a system described by the LVC model
#
# Reads QM.in
# Calculates SOC matrix, dipole moments, gradients, nacs and overlaps
# Writes these back to QM.out

import time
(tc, tt) = (time.clock(), time.time())

import math
import sys
import os
import datetime
import shutil
from copy import deepcopy
try:
  # Importing numpy takes about 100 ms, which is ~50% of the execution time!
  import numpy
  NONUMPY=False
except ImportError:
  import subprocess as sp
  NONUMPY=True

print "Import: CPU time: % .3f s, wall time: %.3f s"%(time.clock() - tc, time.time() - tt)

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
def printgrad(grad,natom,geo,prnorm=False):
  '''Prints a gradient or nac vector. Also prints the atom elements. If the gradient is identical zero, just prints one line.

  Arguments:
  1 list of list of float: gradient
  2 integer: natom
  3 list of list: geometry specs'''

  norm = 0.
  string=''
  iszero=True
  for atom in range(natom):
    string+='%i\t%s\t' % (atom+1,geo[atom][0])
    for xyz in range(3):
      norm += grad[atom][xyz] * grad[atom][xyz]
      if grad[atom][xyz]!=0:
        iszero=False
      string+='% .5f\t' % (grad[atom][xyz])
    string+='\n'
  if prnorm:
    print 'Norm: %.6f'%norm
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
            printgrad(QMout['nacdr'][istate][jstate],natom,QMin['geom'], True)
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
        assert abs(float(line[1+2*iy])) <= 1.e-10
        H[ix][iy]=float(line[0+2*iy])
    U=[ [ 0. for i in range(len(H)) ] for j in range(len(H)) ]
    shift=2+len(H)
    for ix in range(len(H)):
      line=STDOUT[shift+ix].split()
      for iy in range(len(H)):
        assert abs(float(line[1+2*iy])) <= 1.e-10
        U[ix][iy]=float(line[0+2*iy])
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
    #Ucon=[ [ U[i][j].conjugate() for i in range(len(U)) ] for j in range(len(U)) ]
    B=numpy.dot(numpy.array(U).T,numpy.dot(A,U))

  return B

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
  try:
    mantissa, exp = s.split('e')
  except:
    print f, s
    raise
  return "%sE%+0*d"%(mantissa, exp_digits+1, int(exp))


# ======================================================================= #
def writeQMout(QMin,QMout,QMinfilename):
  '''Writes the requested quantities to the file which SHARC reads in. The filename is QMinfilename with everything after the first dot replaced by "out".

  Arguments:
  1 dictionary: QMin
  2 dictionary: QMout
  3 string: QMinfilename'''

  os.chdir(QMin['pwd'])
  k=QMinfilename.find('.')
  if k==-1:
    outfilename=QMinfilename+'.out'
  else:
    outfilename=QMinfilename[:k]+'.out'
  if PRINT:
    print '===> Writing output to file %s in SHARC Format\n' % (outfilename)
  string=''
  for iatom in range(QMin['natom']):
    string+=QMin['geom'][iatom][0]
    for i in range(3):
      string+=' %12.9f ' % (QMin['geom'][iatom][i+1])
    string+='\n'
  string+='\n'
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
  string+=writeQMouttime(QMin,QMout)
  try:
    outfile=open(outfilename,'w')
    outfile.write(string)
    outfile.close()
  except IOError:
    print 'Could not write QM output!'
    sys.exit(15)
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
    sys.exit(16)
  nstates=0
  nmstates=0
  for mult,i in enumerate(states):
    nstates+=i
    nmstates+=(mult+1)*i
  QMin['states']=states
  QMin['nstates']=nstates
  QMin['nmstates']=nmstates
  QMin['nmult'] = 0
  statemap={}
  i=1
  for imult,nstates in enumerate(states):
    if nstates==0:
      continue
    QMin['nmult'] += 1
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
  QMin['geom']=geom

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
    fromfile=os.path.join(QMin['savedir'],'U.out')
    if not os.path.isfile(fromfile):
      print 'ERROR: savedir does not contain U.out! Maybe you need to add "init" to QM.in.'
      sys.exit(17)
    tofile=os.path.join(QMin['savedir'],'Uold.out')
    shutil.copy(fromfile,tofile)


  # find forbidden keywords and optional keywords
  for line in qmin:
    s=line.lower().split()
    if len(s)==0:
      continue
    for t in ['h','soc', 'nacdr', 'dm', 'grad', 'overlap']:
        if s[0] in t:
            QMin[s[0]] = []
    if 'nacdt' in s[0]:
      print 'NACDT is not supported!'
      sys.exit(18)
    if 'dmdr' in s[0]:
      print 'DMDR is not supported!'
      sys.exit(19)

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
def read_LVC_mat(nmstates, header, rfile):
  mat=[ [ complex(0.,0.) for i in range(nmstates) ] for j in range(nmstates) ]

  # real part
  tmp=find_lines(nmstates,header+' R',rfile)
  if not tmp==[]:
    for i, line in enumerate(tmp):
      for j, val in enumerate(line.split()):
        mat[i][j] += float(val)

  # imaginary part
  tmp=find_lines(nmstates,header+' I',rfile)
  if not tmp==[]:
    for i, line in enumerate(tmp):
      for j, val in enumerate(line.split()):
        mat[i][j] += float(val) * 1j

  return mat
# =========================================================
def read_V0(QMin, SH2LVC, fname='V0.txt'):
  """"
  Reads information about the ground-state potential from V0.txt.
  Returns the displacement vector.
  """
  try:
    f=open(fname)
  except IOError:
    print 'Input file %s not found.'%fname
    sys.exit(20)
  v0=f.readlines()
  f.close()

  # read the coordinates and compute Cartesian displacements
  disp=[] # displacement as one 3N-vector
  SH2LVC['Ms']=[] # Squareroot masses in a.u.
  geom = QMin['geom']
  tmp = find_lines(QMin['natom'], 'Geometry',v0)
  for i in range(QMin['natom']):
    s=tmp[i].lower().split()
    if s[0]!=geom[i][0].lower():
      print s[0], geom[i][0]
      print 'Inconsistent atom labels in QM.in and %s!'%fname
      sys.exit(21)
    disp += [geom[i][1] - float(s[2]), geom[i][2] - float(s[3]), geom[i][3] - float(s[4])]
    SH2LVC['Ms'] += 3*[(float(s[5])*U_TO_AMU)**.5]

  # Frequencies (a.u.)
  tmp = find_lines(1, 'Frequencies',v0)
  if tmp==[]:
    print 'No Frequencies defined in %s!'%fname
    sys.exit(22)
  SH2LVC['Om'] = [float(o) for o in tmp[0].split()]

  # Normal modes in mass-weighted coordinates
  tmp = find_lines(len(SH2LVC['Om']), 'Mass-weighted normal modes', v0)
  if tmp==[]:
    print 'No normal modes given in %s!'%fname
    sys.exit(23)
  SH2LVC['V']  = [map(float,line.split()) for line in tmp] # transformation matrix

  return disp

# =========================================================

def read_SH2LVC(QMin, fname='LVC.template'):
  # reads LVC.template, deletes comments and blank lines
  SH2LVC={}
  try:
    f=open(fname)
  except IOError:
    try:
      f=open('SH2LVC.inp')
    except IOError:
      print 'Input file "LVC.template" not found.'
      sys.exit(24)
  sh2lvc=f.readlines()
  f.close()

  disp = read_V0(QMin, SH2LVC, sh2lvc[0].strip())

  # check nstates
  states=[int(s) for s in sh2lvc[1].split()]
  if not states==QMin['states']:
    print 'states from QM.in and nstates from LVC.template are inconsistent!', QMin['states'], states
    sys.exit(25)
  nstates = QMin['nstates']
  nmstates = QMin['nmstates']
  nmult = len(states)
  r3N = range(3*QMin['natom'])
  Om = SH2LVC['Om']

  # Transform the coordinates to dimensionless mass-weighted normal modes
  MR = [SH2LVC['Ms'][i] * disp[i] for i in r3N]
  MRV = [0. for i in r3N]
  for i in r3N:
    MRV[i] = sum(MR[j] * SH2LVC['V'][j][i] for j in r3N)
  Q =  [MRV[i] * Om[i]**0.5 for i in r3N]

  # Compute the ground state potential and gradient
  V0 = sum(0.5 * Om[i] * Q[i]*Q[i] for i in r3N)
  HMCH =  [[[0. for istate in range(states[imult])] for jstate in range(states[imult])] for imult in range(nmult)]
  for imult in range(nmult):
    for istate in range(states[imult]):
      HMCH[imult][istate][istate] = V0

  dHMCH = [[[[0. for istate in range(states[imult])] for jstate in range(states[imult])] for imult in range(nmult)] for i in r3N]
  for i in r3N:
    for imult in range(nmult):
      for istate in range(states[imult]):
        dHMCH[i][imult][istate][istate] = Om[i] * Q[i]

  # Add the vertical energies (epsilon)
  # Enter in separate lines as:
  # <n_epsilon>
  # <mult> <state> <epsilon>
  # <mult> <state> <epsilon>

  tmp = find_lines(1, 'epsilon',sh2lvc)
  if not tmp==[]:
    eps = []
    neps = int(tmp[0])
    tmp = find_lines(neps+1, 'epsilon', sh2lvc)
    for line in tmp[1:]:
      words = line.split()
      eps.append((int(words[0])-1, int(words[1])-1, float(words[-1])))

    for e in eps:
      (imult, istate, val) = e
      HMCH[imult][istate][istate] += val

  #for imult in range(nmult): print numpy.array(HMCH[imult])

  # Add the intrastate LVC constants (kappa)
  # Enter in separate lines as:
  # <n_kappa>
  # <mult> <state> <mode> <kappa>
  # <mult> <state> <mode> <kappa>

  tmp = find_lines(1, 'kappa', sh2lvc)
  if not tmp==[]:
    kappa = []
    nkappa = int(tmp[0])
    tmp = find_lines(nkappa+1, 'kappa', sh2lvc)
    for line in tmp[1:]:
      words = line.split()
      kappa.append((int(words[0])-1, int(words[1])-1, int(words[2])-1, float(words[-1])))

    for k in kappa:
      (imult, istate, i, val) = k
      HMCH[imult][istate][istate]  += val * Q[i]
      dHMCH[i][imult][istate][istate] += val

  # Add the interstate LVC constants (lambda)
  # Enter in separate lines as:
  # <n_lambda>
  # <mult> <state1> <state2> <mode> <lambda>
  # <mult> <state1> <state2> <mode> <lambda>

  tmp = find_lines(1, 'lambda', sh2lvc)
  if not tmp==[]:
    lam = []
    nlam = int(tmp[0])
    tmp = find_lines(nlam+1, 'lambda', sh2lvc)
    for line in tmp[1:]:
      words = line.split()
      lam.append((int(words[0])-1, int(words[1])-1, int(words[2])-1, int(words[3])-1, float(words[-1])))

    for l in lam:
      (imult, istate, jstate, i, val) = l
      HMCH[imult][istate][jstate]  += val * Q[i]
      HMCH[imult][jstate][istate]  += val * Q[i]
      dHMCH[i][imult][istate][jstate] += val
      dHMCH[i][imult][jstate][istate] += val

  SH2LVC['H']  = HMCH
  SH2LVC['dH'] = dHMCH

  SH2LVC['dipole'] = {}
  SH2LVC['dipole'][1]= read_LVC_mat(nmstates, 'DMX', sh2lvc)
  SH2LVC['dipole'][2]= read_LVC_mat(nmstates, 'DMY', sh2lvc)
  SH2LVC['dipole'][3]= read_LVC_mat(nmstates, 'DMZ', sh2lvc)

  # obtain the SOC matrix
  SH2LVC['soc'] = read_LVC_mat(nmstates, 'SOC', sh2lvc)

  return SH2LVC, QMin

# ============================================================================
# ============================================================================
# ============================================================================
def getQMout(QMin,SH2LVC):
  '''Calculates the MCH Hamiltonian, SOC matrix ,overlap matrix, gradients, DM'''

  QMout={}

  nmult = len(QMin['states'])
  r3N = range(3*QMin['natom'])

  # Diagonalize Hamiltonian and expand to the full ms-basis
  U  = [[ 0. for i in range(QMin['nmstates']) ] for j in range(QMin['nmstates']) ]
  Hd = [[ 0. for i in range(QMin['nmstates']) ] for j in range(QMin['nmstates']) ]
  dHfull = [[[ 0. for i in range(QMin['nmstates']) ] for j in range(QMin['nmstates']) ] for iQ in r3N]
  offs = 0
  for imult in range(nmult):
    dim = QMin['states'][imult]
    if not dim == 0:
      Hdtmp,Utmp=diagonalize(SH2LVC['H'][imult])
      for ms in range(imult+1):
        for i in range(dim):
          Hd[i+offs][i+offs] = Hdtmp[i][i]
          U[i+offs][offs:offs+dim] = Utmp[i]
          for iQ in r3N:
            dHfull[iQ][i+offs][offs:offs+dim] = SH2LVC['dH'][iQ][imult][i]
        offs += dim
#  print "QMout1: CPU time: % .3f s, wall time: %.3f s"%(time.clock() - tc, time.time() - tt)

  # Transform the gradients to the MCH basis
  dE = [[[0. for iQ in range(3*QMin['natom'])] for istate in range(QMin['nmstates'])] for jstate in range(QMin['nmstates'])]
  for iQ in r3N:
    dEmat = transform(dHfull[iQ], U)
    for istate in range(QMin['nmstates']):
      for jstate in range(QMin['nmstates']):
        dE[istate][jstate][iQ] = dEmat[istate][jstate]
#  print "QMout2: CPU time: % .3f s, wall time: %.3f s"%(time.clock() - tc, time.time() - tt)

  # Convert the gradient to Cartesian coordinates
  #   -> It would be more efficent to do this only for unique Ms values
  VOdE = [0. for i in r3N]
  grad = []
  for istate in range(QMin['nmstates']):
    OdE = [0. for iQ in r3N]
    for iQ in r3N:
      if abs(SH2LVC['Om'][iQ]) > 1.e-8:
        OdE[iQ] = dE[istate][istate][iQ] * SH2LVC['Om'][iQ]**0.5
    if NONUMPY:
      for iQ in r3N:
        VOdE[iQ] = sum(SH2LVC['V'][iQ][jQ] * OdE[jQ] for jQ in r3N)
    else:
      VOdE = numpy.dot(SH2LVC['V'], OdE)

    grad.append([])
    for iat in range(QMin['natom']):
      grad[-1].append([VOdE[3*iat] * SH2LVC['Ms'][3*iat], VOdE[3*iat+1] * SH2LVC['Ms'][3*iat+1], VOdE[3*iat+2] * SH2LVC['Ms'][3*iat+2]])
#  print "QMout3: CPU time: % .3f s, wall time: %.3f s"%(time.clock() - tc, time.time() - tt)

  if 'nacdr' in QMin:
    nonac = [ [0., 0., 0.] for iat in range(QMin['natom'])]
    QMout['nacdr'] = [[ nonac for istate in range(QMin['nmstates'])] for jstate in range(QMin['nmstates'])]
    istate =- 1
    for imult,ist,ims in itnmstates(QMin['states']):
      istate += 1

      jstate =- 1
      for jmult,jst,jms in itnmstates(QMin['states']):
        jstate += 1

        if imult == jmult and ims == jms and istate < jstate:
          OdE = [0. for iQ in r3N]
          for iQ in r3N:
            if abs(SH2LVC['Om'][iQ]) > 1.e-8:
              OdE[iQ] = dE[istate][jstate][iQ] * SH2LVC['Om'][iQ]**0.5
          if NONUMPY:
            for iQ in r3N:
              VOdE[iQ] = sum(SH2LVC['V'][iQ][jQ] * OdE[jQ] for jQ in r3N)
          else:
            VOdE = numpy.dot(SH2LVC['V'], OdE)

          deriv = []
          for iat in range(QMin['natom']):
            deriv.append([VOdE[3*iat] * SH2LVC['Ms'][3*iat], VOdE[3*iat+1] * SH2LVC['Ms'][3*iat+1], VOdE[3*iat+2] * SH2LVC['Ms'][3*iat+2]])

          Einv = (Hd[jstate][jstate] - Hd[istate][istate]) ** (-1.)

          QMout['nacdr'][istate][jstate] = [ [ c * Einv for c in d] for d in deriv]
          QMout['nacdr'][jstate][istate] = [ [-c * Einv for c in d] for d in deriv]

  # transform dipole matrices
  dipole=[]
  for idir in range(3):
    Dmatrix=transform(SH2LVC['dipole'][idir+1],U)
    dipole.append(Dmatrix)

  # get overlap matrix
  Uoldfile=os.path.join(QMin['savedir'],'Uold.out')
  if 'init' in QMin:
    overlap = [ [ float(i==j) for i in range(QMin['nmstates']) ] for j in range(QMin['nmstates']) ]
  else:
    Uold = [[float(v) for v in line.split()] for line in open(Uoldfile, 'r').readlines()]
    if NONUMPY:
      overlap = [ [ 0. for i in range(QMin['nmstates']) ] for j in range(QMin['nmstates']) ]
      rS = range(QMin['nmstates'])
      for a in rS:
        for b in rS:
          for i in rS:
            overlap[a][b]+=Uold[i][a]*U[i][b]
    else:
      overlap = numpy.dot(numpy.array(Uold).T,U)

  Ufile=os.path.join(QMin['savedir'],'U.out')
  f = open(Ufile, 'w')
  for line in U:
    for c in line:
      f.write(str(c) + ' ')
    f.write('\n')
  f.close()

  # transform SOC matrix
  SO=transform(SH2LVC['soc'],U)
  for i in range(QMin['nmstates']):
    SO[i][i]=complex(0.,0.)
  Hfull=[ [ Hd[i][j]+SO[i][j] for i in range(QMin['nmstates']) ] for j in range(QMin['nmstates']) ]

  # assign QMout elements
  QMout['h']=Hfull
  QMout['dm']=dipole
  QMout['grad']=grad
  #QMout['dmdr']=dmdr
  QMout['overlap']=overlap
  QMout['runtime']=0.

  #pprint.pprint(QMout,width=192)

  return QMout

# ============================================================================
def main():

  QMin=read_QMin()
  SH2LVC,QMin=read_SH2LVC(QMin)
  print "SH2LVC: CPU time: % .3f s, wall time: %.3f s"%(time.clock() - tc, time.time() - tt)

  QMout=getQMout(QMin,SH2LVC)
  print "QMout:  CPU time: % .3f s, wall time: %.3f s"%(time.clock() - tc, time.time() - tt)

  # This print routine takes about 10 ms, i.e. ~5% of the total execution time
  printQMout(QMin,QMout)
  #print "Print:  CPU time: % .3f s, wall time: %.3f s"%(time.clock() - tc, time.time() - tt)

  # Write QMout
  writeQMout(QMin,QMout,'QM.in')
  print "Write:  CPU time: % .3f s, wall time: %.3f s"%(time.clock() - tc, time.time() - tt)

  print "Final:  CPU time: % .3f s, wall time: %.3f s"%(time.clock() - tc, time.time() - tt)
  print '#================ END ================#'

# ============================================================================
if __name__ == '__main__':
  try:
    main()
  except KeyboardInterrupt:
    print '\nCtrl+C makes me a sad SHARC ;-(\n'
    sys.exit(0)
