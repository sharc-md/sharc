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

# Script for printing excitation energies, oscillator strengths and other quantities from QM.out file
# 
# usage python QMout_print.py [options] <QM.out>

import math
import cmath
import sys
import os
import pprint
from optparse import OptionParser

try:
  import numpy
  NONUMPY=False
except ImportError:
  import subprocess as sp
  NONUMPY=True
  
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

# some constants
DEBUG = False
CM_TO_HARTREE = 1./219474.6     #4.556335252e-6 # conversion factor from cm-1 to Hartree
HARTREE_TO_EV = 27.211396132    # conversion factor from Hartree to eV
U_TO_AMU = 1./5.4857990943e-4            # conversion from g/mol to amu
ANG_TO_BOHR = 1./0.529177211    #1.889725989      # conversion from Angstrom to bohr
PI = math.pi

# hash table for conversion of multiplicity to the keywords used in COLUMBUS
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

# ======================================================================= #
def readfile(filename):
  try:
    f=open(filename)
    out=f.readlines()
    f.close()
  except IOError:
    print 'File %s could not be read!' % (filename)
    sys.exit(12)
  return out

# ======================================================================= #
def itnmstates(states):
  for i in range(len(states)):
    if states[i]<1:
      continue
    for k in range(i+1):
      for j in range(states[i]):
        yield i+1,j+1,k-i/2.
  return

# ======================================================================= #
def read_QMout(path,nstates,natom,request):
  targets={'h':         {'flag': 1,
                         'type': complex,
                         'dim':  (1,nstates,nstates)},
           'dm':        {'flag': 2,
                         'type': complex,
                         'dim':  (3,nstates,nstates)},
           'grad':      {'flag': 3,
                         'type': float,
                         'dim':  (nstates,natom,3)}
          }

  # read QM.out
  lines=readfile(path)

  # obtain all targets
  QMout={}
  for t in targets:
    if t in request:
      iline=-1
      while True:
        iline+=1
        if iline>=len(lines):
          print 'Could not find target %s with flag %i in file %s!' % (t,targets[t]['flag'],f)
          sys.exit(11)
        line=lines[iline]
        if '! %i' % (targets[t]['flag']) in line:
          break
      values=[]
      for iblocks in range(targets[t]['dim'][0]):
        iline+=1
        block=[]
        for irow in range(targets[t]['dim'][1]):
          iline+=1
          line=lines[iline].split()
          if targets[t]['type']==complex:
            row=[ complex(float(line[2*i]),float(line[2*i+1])) for i in range(targets[t]['dim'][2]) ]
          elif targets[t]['type']==float:
            row=[ float(line[i]) for i in range(targets[t]['dim'][2]) ]
          block.append(row)
        values.append(block)
      QMout[t]=values

  #pprint.pprint(QMout)
  return QMout

# ======================================================================= #
def read_QMin(path,request):
  QMin={}
  qminlines=readfile(path)
  natom=int(qminlines[0])

  for r in request:
    for iline,line in enumerate(qminlines):
      if iline<=natom+2:
        continue
      s=line.split(None,1)
      if r in s[0]:
        QMin[s[0]]=s[1]
  if 'natom' in request:
    QMin['natom']=natom
  #pprint.pprint(QMin)
  return QMin

# ======================================================================================================================

class diagonalizer:
  def __init__(self):
    exe=os.getenv('SHARC')
    exe=os.path.expanduser(os.path.expandvars(exe))+'/diagonalizer.x'
    if not os.path.isfile(exe):
      print 'SHARC auxilliary diagonalizer not found at %s!' % (exe)
      sys.exit(1)
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

# ======================================================================================================================

def transform(H,DM,P):
  '''transforms the H and DM matrices in the representation where H is diagonal.'''

  if NONUMPY:
    diagon=diagonalizer()
    H,U=diagon.eigh(H)
    UDMU=[ [ [ 0. for i in range(len(H)) ] for j in range(len(H)) ] for k in range(3) ]
    for xyz in range(3):
      temp=[ [ 0. for i in range(len(H)) ] for j in range(len(H)) ]
      for a in range(len(H)):
        for b in range(len(H)):
          for i in range(len(H)):
            temp[a][b]+=U[i][a].conjugate()*DM[xyz][i][b]
      for a in range(len(H)):
        for b in range(len(H)):
          for i in range(len(H)):
            UDMU[xyz][a][b]+=temp[a][i]*U[i][b]
    DM=UDMU

    if P!=None:
      UPU=[ [ 0. for i in range(len(H)) ] for j in range(len(H)) ]
      for a in range(len(H)):
        for b in range(len(H)):
          for i in range(len(H)):
            UPU[a][b]+=U[i][a].conjugate()*P[i][b]
      P=[ [ 0. for i in range(len(H)) ] for j in range(len(H)) ]
      for a in range(len(H)):
        for b in range(len(H)):
          for i in range(len(H)):
            P[a][b]+=temp[a][i]*U[i][b]

  else:
    eig,U=numpy.linalg.eigh(H)
    Ucon=[ [ 0. for i in range(len(H)) ] for j in range(len(H)) ]
    for ix in range(len(U)):
      for iy in range(len(U)):
        Ucon[ix][iy]=U[iy][ix].conjugate()
        if ix==iy:
          H[ix][iy]=complex(eig[ix])
        else:
          H[ix][iy]=complex(0)
    UDMU=[0,0,0]
    for xyz in range(3):
      UDMU[xyz]=numpy.dot(Ucon,numpy.dot(DM[xyz],U))
    DM=UDMU

    if P!=None:
      UPU=numpy.dot(Ucon,numpy.dot(P,U))
      P=UPU

  return H,DM,U

# ========================== Main Code =============================== #
def main():

  usage='''
QMout_print.py [options] QM.out

This script reads a QM.out file from a SHARC interface and prints
excitation energies and oscillator strengths.
'''

  description=''

  parser = OptionParser(usage=usage, description=description)
  parser.add_option('-i', dest='i', type=str, nargs=1, default='', help="QM.in file (to read number of states)")
  parser.add_option('-e', dest='e', type=float, nargs=1, default=0.0, help="Absolute energy shift (float, default=0.0)")
  parser.add_option('-s', dest='s', type=str, nargs=1, default='', help="Number of states (in quotes separated by whitespace)")
  parser.add_option('-n', dest='n', type=int, nargs=1, default=1, help="Number of atoms")
  parser.add_option('-D', dest='D', action='store_true',help="Diagonalize")

  #parser.add_option('-n', dest='n', type=int, nargs=1, default=3, help="Number of geometries to be generated (integer, default=3)")
  #parser.add_option('-r', dest='r', type=int, nargs=1, default=16661, help="Seed for the random number generator (integer, default=16661)")
  #parser.add_option('--MOLPRO', dest='M', action='store_true',help="Assume a MOLPRO frequency file (default=assume MOLDEN file)")
  #parser.add_option('-m', dest='m', action='store_true',help="Enter non-default atom masses")
  #parser.add_option('--keep_trans_rot', dest='KTR', action='store_true',help="Keep translational and rotational components")

  (options, args) = parser.parse_args()
  ezero=options.e
  qminfile=options.i
  qmoutfile=args[0]

  QMin={'states': 1, 'natom': options.n}
  if qminfile!='':
    QMin=read_QMin(qminfile,['states','natom'])
    sstates=QMin['states']
  elif options.s!='':
    sstates=options.s
  else:
    sstates='1'
  states=[]
  for i in sstates.split():
    states.append(int(i))
  QMin['states']=states
  nmstates=0
  for i in range(len(states)):
    nmstates+=states[i]*(i+1)
  QMin['nmstates']=nmstates

  # obtain the statemap 
  statemap={}
  i=1
  for imult,istate,ims in itnmstates(QMin['states']):
    statemap[i]=[imult,istate,ims]
    i+=1
  QMin['statemap']=statemap

  QMout=read_QMout(qmoutfile,QMin['nmstates'],QMin['natom'],['h','dm'])

  sys.stderr.write( 'Number of states: %s\n' % (states) )
  sys.stderr.write( '%5s  %11s %16s %12s %12s   %6s\n' % ('State','Label','E (E_h)','dE (eV)','f_osc','Spin') )

  if options.D:
    h,dm,U=transform(QMout['h'][0],QMout['dm'],None)
    QMout['h']=[h]
    QMout['dm']=dm

  #pprint.pprint(QMin)
  #pprint.pprint(QMout)


  energies=[]
  fosc=[]
  if options.D:
    for istate in range(QMin['nmstates']):
      e=QMout['h'][0][istate][istate].real
      energies.append(e)
      # spin
      spin=0.
      ist=0
      imax=0.
      for jstate in range(QMin['nmstates']):
        m,s,ms=QMin['statemap'][jstate+1]
        c=(U[jstate][istate]*U[jstate][istate].conjugate()).real
        spin+=m*c
        if c>imax:
          ist=(m,s)
          imax=c
      # fosc
      dmx=QMout['dm'][0][istate][0].real
      dmy=QMout['dm'][1][istate][0].real
      dmz=QMout['dm'][2][istate][0].real
      f=2./3.*(e-energies[0])*(dmx**2+dmy**2+dmz**2)
      fosc.append(f)
      #else:
        #dmx=dmy=dmz=0.
        #fosc.append(0.)
      if ezero!=0.0:
        de=(e-ezero)*27.21
      else:
        de=(e-energies[0])*27.21
      string='%5i %10s%02i %16.10f %12.8f %12.8f   %6.4f' % (istate+1,IToMult[ist[0]][0],ist[1]-(ist[0]<=2),e,de,fosc[-1],spin)
      print string
  else:
    for istate in range(QMin['nmstates']):
      m,s,ms=QMin['statemap'][istate+1]
      if 2*ms+1!=m:
        continue
      e=QMout['h'][0][istate][istate].real
      energies.append(e)
      if m==1 and s>0:
        dmx=QMout['dm'][0][istate][0].real
        dmy=QMout['dm'][1][istate][0].real
        dmz=QMout['dm'][2][istate][0].real
        f=2./3.*(e-energies[0])*(dmx**2+dmy**2+dmz**2)
        fosc.append(f)
      else:
        dmx=dmy=dmz=0.
        fosc.append(0.)
      if ezero!=0.0:
        de=(e-ezero)*27.21
      else:
        de=(e-energies[0])*27.21
      string='%5i %10s%02i %16.10f %12.8f %12.8f   %6.4f' % (istate+1,IToMult[m][0],s-(m<=2),e,de,fosc[-1],m)
      print string










if __name__ == '__main__':
    main()