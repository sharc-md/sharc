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

import sys
import os
#import re
#import string
#import math
#import imp
try:
  import numpy
except ImportError:
  print 'The kf module required to read ADF binary files needs numpy. Please install numpy and then try again'
  sys.exit()

adf=os.path.expandvars('$ADFHOME')
sys.path.append(adf+'/scripting')
import kf

# ================================================

def read_t21(filename):
  f = kf.kffile(filename)

  try:
    f.sections()
  except:
    print 'File does not seem to be a TAPE21 file...'
    return None

  Freq = f.read("Freq","Frequencies").tolist()

  natom = int(f.read("Freq","nr of atoms"))

  fragtype=f.read("Geometry","fragmenttype").tolist()
  atomtype_index=f.read("Geometry","fragment and atomtype index").tolist()[natom:]
  Atomsymbs=[ fragtype[i-1] for i in atomtype_index ]

  xyz = f.read("Freq","xyz").tolist()
  FreqCoord=[]
  x=0
  for i in range(natom):
    atom=[Atomsymbs[i]]+xyz[x:x+3]
    FreqCoord.append(atom)
    x+=3

  Normalmodes= f.read("Freq","Normalmodes").tolist()
  Modes={}
  for i in range(3*natom):
    m=[]
    for j in range(natom):
      n=Normalmodes[3*natom*i+3*j:3*natom*i+3*j+3]
      m.append(n)
    Modes[i]=m

  Int={}

  F={'FREQ':            Freq,
     'FR-COORD':        FreqCoord,
     'FR-NORM-COORD':   Modes,
     'INT':             Int,
     'natom':           natom}
  return F

# ================================================

def read_ADFout(filename):
  data=open(filename).readlines()

  iline=-1
  while True:
    iline+=1
    if iline>=len(data):
      print 'Could not find Frequencies output!'
      return None
    line=data[iline]
    if 'F R E Q U E N C' in line:
      break
  iline+=10
  FreqCoord=[]
  natom=0
  while True:
    iline+=1
    line=data[iline]
    if '----' in line:
      break
    s=line.split()
    try:
      atom=[ s[1], float(s[2]), float(s[3]), float(s[4]) ]
    except IndexError:
      print 'Could not find optimized coordinates!'
      return None
    FreqCoord.append(atom)
    natom+=1

  while True:
    iline+=1
    line=data[iline]
    if 'Vibrations and Normal Modes  ***  (cartesian coordinates, NOT mass-weighted)  ***' in line:
      break
  iline+=6
  Modes={}
  x=0
  y=0
  while True:
    iline+=1
    line=data[iline]
    if 'List of All Frequencies:' in line:
      break
    s=line.split()
    if '----' in line:
      y=len(s)
      x+=y
      for i in range(y):
        Modes[x-y+i]=[]
    if len(s)<=3:
      continue
    for i in range(y):
      m=[ float(s[1+3*i]),float(s[2+3*i]),float(s[3+3*i]) ]
      Modes[x-y+i].append(m)

  iline+=8
  Int={}
  Freq=[]
  x=0
  while True:
    iline+=1
    line=data[iline]
    s=line.split()
    if len(s)==0:
      break
    Freq.append(float(s[0]))
    Int[x]=float(s[2])
    x+=1

  F={'FREQ':            Freq,
     'FR-COORD':        FreqCoord,
     'FR-NORM-COORD':   Modes,
     'INT':             Int,
     'natom':           natom}
  return F

# ================================================

def format_molden(F):
  string='[MOLDEN FORMAT]\n[FREQ]\n'
  for i in F['FREQ']:
    string+='%8.2f\n' % i
  string+='[FR-COORD]\n'
  for atom in F['FR-COORD']:
    string+='%4s  %12.9f  %12.9f  %12.9f\n' % tuple(atom)
  string+='[FR-NORM-COORD]\n'
  for i in range(len(F['FR-NORM-COORD'])):
    mode=F['FR-NORM-COORD'][i]
    string+='Vibration %4i\n' % (i+1)
    for m in mode:
      string+='  %12.9f  %12.9f  %12.9f\n' % tuple(m)
  if len(F['INT'])>0:
    string+='[INT]\n'
    for i in range(len(F['FR-NORM-COORD'])):
      if i in F['INT']:
        string+='%12.6f\n' % F['INT'][i]
      else:
        string+='%12.6f\n' % 0.
  return string






# ================================================

filename = sys.argv[1]

F=read_t21(filename)
if not F:
  F=read_ADFout(filename)
if not F:
  print 'Could not parse file %s!' % filename
  sys.exit(1)

string=format_molden(F)
outfile = open(filename+'.molden','w')
outfile.write(string)
outfile.close()



# save the shell command
command='python '+' '.join(sys.argv)
f=open('KEYSTROKES.ADF_freq','w')
f.write(command)
f.close()
