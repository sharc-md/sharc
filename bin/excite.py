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

# Script for the calculation of Wigner distributions from molden frequency files
# 
# usage 
import sys
if sys.version_info[0]!=2:
  sys.stdout.write('The SHARC suite is not compatible with Python 3! Use Python 2 (>2.6)!')
  sys.exit(0)


import copy
import math
import re
import os
import stat
import shutil
import random
import datetime
from optparse import OptionParser
import readline
import time

try:
  import numpy
  NONUMPY=False
except ImportError:
  import subprocess as sp
  NONUMPY=True

# =========================================================0
# compatibility stuff

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
BOHR_TO_ANG=0.529177211
PI = math.pi

version='2.0'
versionneeded=[0.2, 1.0, 2.0, float(version)]
versiondate=datetime.date(2018,2,1)

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def try_read(l,index,typefunc,default):
  try:
    if typefunc==bool:
      return 'True'==l[index]
    else:
      return typefunc(l[index])
  except IndexError:
    return typefunc(default)
  except ValueError:
    print 'Could not initialize object!'
    quit(1)

# ======================================================================================================================

class ATOM:
  def __init__(self,symb='??',num=0.,coord=[0.,0.,0.],m=0.,veloc=[0.,0.,0.]):
    self.symb  = symb
    self.num   = num
    self.coord = coord
    self.mass  = m
    self.veloc = veloc
    self.Ekin=0.5*self.mass * sum( [ self.veloc[i]**2 for i in range(3) ] )

  def init_from_str(self,initstring=''):
    f=initstring.split()
    self.symb  =   try_read(f,0,str,  '??')
    self.num   =   try_read(f,1,float,0.)
    self.coord = [ try_read(f,i,float,0.) for i in range(2,5) ]
    self.mass  =   try_read(f,5,float,0.)*U_TO_AMU
    self.veloc = [ try_read(f,i,float,0.) for i in range(6,9) ]
    self.Ekin=0.5*self.mass * sum( [ self.veloc[i]**2 for i in range(3) ] )

  def __str__(self):
    s ='%2s % 5.1f '               % (self.symb, self.num)
    s+='% 12.8f % 12.8f % 12.8f '  % tuple(self.coord)
    s+='% 12.8f '                  % (self.mass/U_TO_AMU)
    s+='% 12.8f % 12.8f % 12.8f'   % tuple(self.veloc)
    return s

  def EKIN(self):
    self.Ekin=0.5*self.mass * sum( [ self.veloc[i]**2 for i in range(3) ] )
    return self.Ekin

  def geomstring(self):
    s='  %2s % 5.1f % 12.8f % 12.8f % 12.8f % 12.8f' % (self.symb,self.num,self.coord[0],self.coord[1],self.coord[2],self.mass/U_TO_AMU)
    return s

  def velocstring(self):
    s=' '*11+'% 12.8f % 12.8f % 12.8f' % tuple(self.veloc)
    return s

# ======================================================================================================================

class STATE:
  def __init__(self,i=0,e=0.,eref=0.,dip=[0.,0.,0.]):
    self.i       = i
    self.e       = e.real
    self.eref    = eref.real
    self.dip     = dip
    self.Excited = False
    self.Eexc    = self.e-self.eref
    self.Fosc    = (2./3.*self.Eexc*sum( [i*i.conjugate() for i in self.dip] ) ).real
    if self.Eexc==0.:
      self.Prob  = 0.
    else:
      self.Prob  = self.Fosc/self.Eexc**2

  def init_from_str(self,initstring):
    f=initstring.split()
    self.i       =   try_read(f,0,int,  0 )
    self.e       =   try_read(f,1,float,0.)
    self.eref    =   try_read(f,2,float,0.)
    self.dip     = [ complex( try_read(f,i,float,0.),try_read(f,i+1,float,0.) ) for i in [3,5,7] ]
    self.Excited =   try_read(f,11,bool, False)
    self.Eexc    = self.e-self.eref
    self.Fosc    = (2./3.*self.Eexc*sum( [i*i.conjugate() for i in self.dip] ) ).real
    if self.Eexc==0.:
      self.Prob  = 0.
    else:
      self.Prob  = self.Fosc/self.Eexc**2

  def __str__(self):
    s ='%03i % 18.10f % 18.10f ' % (self.i,self.e,self.eref)
    for i in range(3):
      s+='% 12.8f % 12.8f ' % (self.dip[i].real,self.dip[i].imag)
    s+='% 12.8f % 12.8f %s' % (self.Eexc*HARTREE_TO_EV,self.Fosc,self.Excited)
    return s

  def Excite(self,max_Prob,erange):
    try:
      Prob=self.Prob/max_Prob
    except ZeroDivisionError:
      Prob=-1.
    if not (erange[0] <= self.Eexc <= erange[1]):
      Prob=-1.
    self.Excited=(random.random() < Prob)

# ======================================================================================================================

class INITCOND:
  def __init__(self,atomlist=[],eref=0.,epot_harm=0.):
    self.atomlist=atomlist
    self.eref=eref
    self.Epot_harm=epot_harm
    self.natom=len(atomlist)
    self.Ekin=sum( [atom.Ekin for atom in self.atomlist] )
    self.statelist=[]
    self.nstate=0
    self.Epot=epot_harm

  def addstates(self,statelist):
    self.statelist=statelist
    self.nstate=len(statelist)
    self.Epot=self.statelist[0].e-self.eref

  def init_from_file(self,f,eref,index):
    while True: 
      line=f.readline()
      #if 'Index     %i' % (index) in line:
      if re.search('Index\s+%i' % (index),line):
        break
      if line=='\n':
        continue
      if line=='':
        print 'Initial condition %i not found in file %s' % (index,f.name)
        quit(1)
    f.readline()        # skip one line, where "Atoms" stands
    atomlist=[]
    while True:
      line=f.readline()
      if 'States' in line:
        break
      atom=ATOM()
      atom.init_from_str(line)
      atomlist.append(atom)
    statelist=[]
    while True:
      line=f.readline()
      if 'Ekin' in line:
        break
      state=STATE()
      state.init_from_str(line)
      statelist.append(state)
    epot_harm=0.
    while not line=='\n' and not line=='':
      line=f.readline()
      if 'epot_harm' in line.lower():
        epot_harm=float(line.split()[1])
        break
    self.atomlist=atomlist
    self.eref=eref
    self.Epot_harm=epot_harm
    self.natom=len(atomlist)
    self.Ekin=sum( [atom.Ekin for atom in self.atomlist] )
    self.statelist=statelist
    self.nstate=len(statelist)
    if self.nstate>0:
      self.Epot=self.statelist[0].e-self.eref
    else:
      self.Epot=epot_harm

  def __str__(self):
    s='Atoms\n'
    for atom in self.atomlist:
      s+=str(atom)+'\n'
    s+='States\n'
    for state in self.statelist:
      s+=str(state)+'\n'
    s+='Ekin      % 16.12f a.u.\n' % (self.Ekin)
    s+='Epot_harm % 16.12f a.u.\n' % (self.Epot_harm)
    s+='Epot      % 16.12f a.u.\n' % (self.Epot)
    s+='Etot_harm % 16.12f a.u.\n' % (self.Epot_harm+self.Ekin)
    s+='Etot      % 16.12f a.u.\n' % (self.Epot+self.Ekin)
    s+='\n\n'
    return s

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def itnmstates(states):
  for i in range(len(states)):
    if states[i]<1:
      continue
    for k in range(i+1):
      for j in range(states[i]):
        yield i+1,j+1,k-i/2.
  return

def get_statemap(states):
  statemap={}
  i=1
  for imult,istate,ims in itnmstates(states):
    statemap[i]=[imult,istate,ims]
    i+=1
  return statemap

def print_statemap(statemap,diag=False):
  n=len(statemap)
  if diag:
    s='# State map for diagonal states:\n#State\tQuant\n'
    for i in range(1,n+1):
      s+='%i\t%i\n' % (i,i)
  else:
    s='# State map for MCH states:\n#State\tMult\tM_s\tQuant\n'
    for i in range(1,n+1):
      (mult,state,ms)=statemap[i]
      s+='%i\t%i\t%+3.1f\t%i\n' % (i,mult,ms,state)
  return s

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def check_initcond_version(string,must_be_excited=False):
  if not 'sharc initial conditions file' in string.lower():
    return False
  f=string.split()
  for i,field in enumerate(f):
    if 'version' in field.lower():
      try:
        v=float(f[i+1])
        if not v in versionneeded:
          return False
      except IndexError:
        return False
  if must_be_excited:
    if not 'excited' in string.lower():
      return False
  return True


# ======================================================================================================================

def centerstring(string,n,pad=' '):
  l=len(string)
  if l>=n:
    return string
  else:
    return  pad*((n-l+1)/2)+string+pad*((n-l)/2)

# ======================================================================================================================

def displaywelcome():
  string='\n'
  string+='  '+'='*80+'\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Excite initial conditions for SHARC',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Author: Sebastian Mai',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Version:'+version,80)+'||\n'
  string+='||'+centerstring(versiondate.strftime("%d.%m.%y"),80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  string+='''
This script automatizes to read-out the results of initial excited-state calculations for SHARC.
It calculates oscillator strength (in MCH and diagonal basis) and stochastically 
determines initial states for trajectories.
  '''
  print string

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def open_keystrokes():
  global KEYSTROKES
  KEYSTROKES=open('KEYSTROKES.tmp','w')

def close_keystrokes():
  KEYSTROKES.close()
  shutil.move('KEYSTROKES.tmp','KEYSTROKES.excite')

# ===================================

def question(question,typefunc,default=None,autocomplete=True,ranges=False):
  if typefunc==int or typefunc==float:
    if not default==None and not isinstance(default,list):
      print 'Default to int or float question must be list!'
      quit(1)
  if typefunc==str and autocomplete:
    readline.set_completer_delims(' \t\n;')
    readline.parse_and_bind("tab: complete")    # activate autocomplete
  else:
    readline.parse_and_bind("tab: ")            # deactivate autocomplete

  while True:
    s=question
    if default!=None:
      if typefunc==bool or typefunc==str:
        s+= ' [%s]' % (str(default))
      elif typefunc==int or typefunc==float:
        s+= ' ['
        for i in default:
          s+=str(i)+' '
        s=s[:-1]+']'
    if typefunc==str and autocomplete:
      s+=' (autocomplete enabled)'
    if typefunc==int and ranges:
      s+=' (range comprehension enabled)'
    s+=' '

    line=raw_input(s)
    line=re.sub('#.*$','',line).strip()
    if not typefunc==str:
      line=line.lower()

    if line=='' or line=='\n':
      if default!=None:
        KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
        return default
      else:
        continue

    if typefunc==bool:
      posresponse=['y','yes','true', 't', 'ja',  'si','yea','yeah','aye','sure','definitely']
      negresponse=['n','no', 'false', 'f', 'nein', 'nope']
      if line in posresponse:
        KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
        return True
      elif line in negresponse:
        KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
        return False
      else:
        print 'I didn''t understand you.'
        continue

    if typefunc==str:
      KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
      return line

    if typefunc==float:
      # float will be returned as a list
      f=line.split()
      try:
        for i in range(len(f)):
          f[i]=typefunc(f[i])
        KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
        return f
      except ValueError:
        print 'Please enter floats!'
        continue

    if typefunc==int:
      # int will be returned as a list
      f=line.split()
      out=[]
      try:
        for i in f:
          if ranges and '~' in i:
            q=i.split('~')
            for j in range(int(q[0]),int(q[1])+1):
              out.append(j)
          else:
            out.append(int(i))
        KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
        return out
      except ValueError:
        if ranges:
          print 'Please enter integers or ranges of integers (e.g. "-3~-1  2  5~7")!'
        else:
          print 'Please enter integers!'
        continue


# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def read_matrix(qmout,i,filename):
  line=qmout[i].split()
  try:
    n1=int(line[0])
    n2=int(line[1])
  except ValueError:
    print 'Could not read number of states of matrix in %s' % (filename)
    return None
  if not n1==n2:
    print 'Non-square matrix in %s!' % (filename)
    return None
  i+=1
  A=[]
  try:
    for irow in range(n1):
      line=qmout[i].split()
      a=[]
      for icol in range(n1):
        re=float(line[0+2*icol])
        im=float(line[1+2*icol])
        a.append(complex(re,im))
      A.append(a)
      i+=1
  except (ValueError,IndexError):
    print 'Matrix malformatted in %s' % (filename)
    return None
  return A

# =======================================

def find_flag(qmout,flag,filename):
  i=0
  try:
    while not '! %i' % (flag) in qmout[i]:
      i+=1
  except IndexError:
    print 'No matrix with flag %i in %s!' % (flag,filename)
    return None
  return i

# =======================================

def extractQMout(filename,readP=False,readS=False):
  '''Takes the path to a QM.out file and returns the Hamiltonian, the Dipole matrices and the property matrix'''
  try:
    qmoutf=open(filename,'r')
    qmout=qmoutf.readlines()
    qmoutf.close()
  except IOError:
    print 'Could not find %s!' % (filename)
    return None,None,None,None

  i=find_flag(qmout,1,filename)
  if i==None:
    return None,None,None,None
  H=read_matrix(qmout,i+1,filename)

  i=find_flag(qmout,2,filename)
  if i==None:
    DM=None
  else:
    #return H,None,None
    DM=[]
    for idir in range(3):
      DM.append(read_matrix(qmout,i+1,filename))
      if DM[-1]==None:
        DM=None
        break
      i+=len(DM[-1])+1

  if readP:
    i=find_flag(qmout,11,filename)
    if i==None:
      return H,DM,None,None
    P=read_matrix(qmout,i+1,filename)
  else:
    P=None

  if readS:
    i=find_flag(qmout,6,filename)
    if i==None:
      return H,DM,P,None
    S=read_matrix(qmout,i+1,filename)
  else:
    S=None

  return H,DM,P,S

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

  return H,DM,P

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def get_infos(INFOS):
  '''This routine asks for the paths of the initconds file and ICONDS directory, for energy window and the representation.'''

  print centerstring('Initial conditions file',60,'-')+'\n'
  # open the initconds file
  try:
    initfile='initconds'
    initf=open(initfile)
    line=initf.readline()
    if check_initcond_version(line):
      print 'Initial conditions file "initconds" detected. Do you want to use this?'
      if not question('Use file "initconds"?',bool,True):
        initf.close()
        raise IOError
    else:
      initf.close()
      raise IOError
  except IOError:
    print '\nIf you do not have an initial conditions file, prepare one with wigner.py!\n'
    print 'Please enter the filename of the initial conditions file.'
    while True:
      initfile=question('Initial conditions filename:',str,'initconds')
      initfile=os.path.expanduser(os.path.expandvars(initfile))
      if os.path.isdir(initfile):
        print 'Is a directory: %s' % (initfile)
        continue
      if not os.path.isfile(initfile):
        print 'File does not exist: %s' % (initfile)
        continue
      try:
        initf=open(initfile,'r')
      except IOError:
        print 'Could not open: %s' % (initfile)
        continue
      line=initf.readline()
      if check_initcond_version(line):
        break
      else:
        print 'File does not contain initial conditions!'
        continue
  # read the header
  INFOS['ninit']=int(initf.readline().split()[1])
  INFOS['natom']=int(initf.readline().split()[1])
  INFOS['repr']=initf.readline().split()[1]
  INFOS['eref']=float(initf.readline().split()[1])
  INFOS['eharm']=float(initf.readline().split()[1])

  # get guess for number of states
  line=initf.readline()
  if 'states' in line.lower():
    states=[]
    l=line.split()
    for i in range(1,len(l)):
      states.append(int(l[i]))
    INFOS['states']=states
  else:
    INFOS['states']=None

  while True:
    line=initf.readline()
    if 'Equilibrium' in line:
      break
    if line=='':
      print 'File malformatted! No equilibrium geometry!'
      quit(1)
  equi=[]
  for i in range(INFOS['natom']):
    line=initf.readline()
    atom=ATOM()
    atom.init_from_str(line)
    equi.append(atom)
  INFOS['equi']=equi
  initf.seek(0)                 # rewind the initf file
  INFOS['initf']=initf
  print '\nFile "%s" contains %i initial conditions.' % (initfile,INFOS['ninit'])
  print 'Number of atoms is %i\n' % (INFOS['natom'])


  print centerstring('Generate excited state lists',60,'-')+'\n'
  print '''Using the following options, excited state lists can be added to the initial conditions:

1       Generate a list of dummy states
2       Read excited-state information from ab initio calculations (from setup_init.py)'''
  allowed=[1,2]
  guess_gen=[2]
  if any([i in INFOS['repr'].lower()   for i in ['mch','diag']]):
    allowed.append(3)
    print '3       Keep existing excited-state information'
    guess_gen=[3]
  print ''
  while True:
    INFOS['gen_list']=question('How should the excited-state lists be generated?',int,guess_gen)[0]
    if not INFOS['gen_list'] in allowed:
      print 'Please give one of the following integer: %s' % (allowed)
      continue
    break

  if INFOS['gen_list']==1:
    INFOS['read_QMout']=False
    INFOS['make_list']=True
  elif INFOS['gen_list']==2:
    INFOS['read_QMout']=True
    INFOS['make_list']=False
  elif INFOS['gen_list']==3:
    INFOS['read_QMout']=False
    INFOS['make_list']=False


  if INFOS['read_QMout']:
    print 'Please enter the path to the directory containing the ICOND subdirectories.'
    while True:
      path=question('Path to ICOND directories:',str)
      path=os.path.expanduser(os.path.expandvars(path))
      path=os.path.abspath(path)
      if not os.path.isdir(path):
        print 'Is not a directory or does not exist: %s' % (path)
        continue
      else:
        ls=os.listdir(path)
        n=0
        for i in ls:
          if 'ICOND' in i:
            n+=1
        if n==0:
          print 'Does not contain any ICOND directories: %s' % (path)
          continue
        else:
          break
    print '\n%s\nDirectory contains %i subdirectories.' % (path,n)
    if n<INFOS['ninit']+1:
      print 'There are more initial conditions in %s.' % (initfile)
    INFOS['iconddir']=path
    INFOS['ncond']=n
    print ''


  if INFOS['make_list']:
    print '\nPlease enter the number of states as a list of integers\ne.g. 3 0 3 for three singlets, zero doublets and three triplets.'
    while True:
      states=question('Number of states:',int)
      if len(states)==0:
        continue
      if any(i<0 for i in states):
        print 'Number of states must be positive!'
        continue
      break
    print ''
    nstates=0
    for mult,i in enumerate(states):
      nstates+=(mult+1)*i
    print 'Number of states: '+str(states)
    print 'Total number of states: %i\n' % (nstates)
    print ''
    INFOS['states']=states
    INFOS['nstates']=nstates


  if INFOS['make_list'] or INFOS['read_QMout']:
    print centerstring('Excited-state representation',60,'-')
    if INFOS['read_QMout']:
      print '''\nThis script can calculate the excited-state energies and oscillator strengths in two representations.
These representations are:
- MCH representation: Only the diagonal elements of the Hamiltonian are taken into account. The states are the spin-free states as calculated in the quantum chemistry code. This option is usually sufficient for systems with small SOC (below 300 cm^-1).
- diagonal representation: The Hamiltonian including spin-orbit coupling is diagonalized. The states are spin-corrected, fully adiabatic. Note that for this the excited-state calculations have to include spin-orbit couplings. This is usually not necessary for systems with small SOC.
'''
    else:
      print '''\nThis script needs to set the electronic state representation.
There are two representations:
- MCH representation: Only the diagonal elements of the Hamiltonian are taken into account. The states are the spin-free states as calculated in the quantum chemistry code. This option is usually sufficient for systems with small SOC (below 300 cm^-1).
- diagonal representation: The Hamiltonian including spin-orbit coupling is diagonalized. The states are spin-corrected, fully adiabatic. Note that for this the excited-state calculations have to include spin-orbit couplings. This is usually not necessary for systems with small SOC.
'''
    INFOS['diag']=question('Do you want to use the diagonal representation (yes=diag, no=MCH)?',bool)
    if INFOS['diag'] and INFOS['read_QMout']:
      qmfilename=INFOS['iconddir']+'/ICOND_00000/QM.in'
      if os.path.isfile(qmfilename):
        soc_there=False
        qmfile=open(qmfilename,'r')
        for line in qmfile:
          if 'soc' in line.lower():
            soc_there=True
        qmfile.close()
        if not soc_there:
          print '\nDiagonal representation specified, but \n%s\n says there are no SOCs in the QM.out files.\nUsing MCH representation.' % (qmfilename)
          INFOS['diag']=False
          time.sleep(2)
      else:
        print 'Could not determine whether calculations include SOC.'
    print ''
    if INFOS['diag']:
      INFOS['repr']='diag'
    else:
      INFOS['repr']='MCH'


    if INFOS['read_QMout']:
      qmfilename=INFOS['iconddir']+'/ICOND_00000/QM.in'
      if os.path.isfile(qmfilename):
        qmfile=open(qmfilename,'r')
        for line in qmfile:
          if re.search('^\s?ion\s?',line.lower()):
            INFOS['ion']=question('Use ionization probabilities instead of dipole moments?',bool,False)
          if 'states' in line.lower():
            states=[]
            l=line.split()
            for i in range(1,len(l)):
              states.append(int(l[i]))
            INFOS['states']=states
        qmfile.close()
      if not 'ion' in INFOS:
        INFOS['ion']=False


    print '\n'+centerstring('Reference energy',60,'-')+'\n'
    if INFOS['read_QMout']:
      qmfilename=INFOS['iconddir']+'/ICOND_00000/QM.out'
    if INFOS['make_list']:
      eref_from_file=question('Do you have conducted an ab initio calculation at the equilibrium geometry?',bool)
      if eref_from_file:
        while True:
          qmfilename=question('Path to the QM.out file of the calculation:',str)
          if not os.path.isfile(qmfilename):
            print 'File %s does not exist!' % (qmfilename)
            continue
          break
      else:
        qmfilename=''
    if os.path.isfile(qmfilename):
      H,DM,P,Smat=extractQMout(qmfilename)
      if H!=None:
        if INFOS['diag']:
          H,DM,P=transform(H,DM,P)
        INFOS['eref']=H[0][0].real
        print 'Reference energy read from file \n%s' % (qmfilename)
        print 'E_ref= %16.12f' % (INFOS['eref'])
    else:
      print '\nPlease enter the ground state equilibrium energy in hartree.'
      INFOS['eref']=question('Reference energy (hartree): ',float)[0]
    print ''




  print '\n'+centerstring('Excited-state selection',60,'-')+'\n'
  print '''Using the following options, the excited states can be flagged as valid initial states for dynamics:

1       Unselect all initial states
2       Provide a list of desired initial states'''
  allowed=[1,2]
  guess_gen=[2]
  if not INFOS['make_list']:
    print '3       Simulate delta-pulse excitation based on excitation energies and oscillator strengths'
    allowed.append(3)
    guess_gen=[3]
  if not INFOS['make_list'] and not INFOS['read_QMout']:
    print '4       Keep selection (i.e., only print statistics on the excited states and exit)'
    allowed.append(4)
  print ''
  while True:
    INFOS['excite']=question('How should the excited states be flagged?',int,guess_gen)[0]
    if not INFOS['excite'] in allowed:
      print 'Please give one of the following integer: %s' % (allowed)
      continue
    break
  print ''



  if INFOS['excite']==1:
    INFOS['allowed']=set()
    INFOS['erange']=[-2.,-1.]


  if INFOS['excite']==3 or (INFOS['excite']==2 and not INFOS['make_list']) or INFOS['excite']==4:
    print '\n'+centerstring('Excitation window',60,'-')
    if INFOS['excite']==4:
      print '\nEnter the energy window for counting.'
    else:
      print '\nEnter the energy window for exciting the trajectories.'
    while True:
      erange=question('Range (eV):',float,[0.,10.])
      if erange[0]>=erange[1]:
        print 'Range empty!'
        continue
      break
    print '\nScript will allow excitations only between %f eV and %f eV.\n' % (erange[0],erange[1])
    erange[0]/=HARTREE_TO_EV
    erange[1]/=HARTREE_TO_EV
    INFOS['erange']=erange




  INFOS['diabatize']=False
  if INFOS['excite']==2:
    print '\n'+centerstring('Considered states',60,'-')

    if INFOS['read_QMout'] and INFOS['repr']=='MCH':
      qmfilename=INFOS['iconddir']+'/ICOND_00001/QM.in'
      if os.path.isfile(qmfilename):
        qmfile=open(qmfilename,'r')
        for line in qmfile:
          if re.search('^\s?overlap\s?',line.lower()):
            print '\nThe vertical excitation calculations were done with overlaps with a reference.\nReference overlaps can be used to obtain diabatic states.\n'
            INFOS['diabatize']=question('Do you want to specify the initial states in a diabatic picture?',bool,False)
        qmfile.close()

    print '''\nPlease give a list of all states which should be 
flagged as valid initial states for the dynamics.
Note that this is applied to all initial conditions.'''
    if INFOS['diabatize']:
      print '\nNOTE: These numbers are interpreted as diabatic states.\nThe diabatic basis is the set of states computed in ICOND_00000/.\nPlease carefully analyze these states to decide which diabatic states to request.'
      #print 'NOTE: You can only enter one initial state.'
    if 'states' in INFOS:
      diago=(INFOS['repr']=='diag')
      print print_statemap(get_statemap(INFOS['states']),diag=diago)

    while True:
      allowed_states=question('List of initial states:',int,ranges=True)
      if any([i<=0 for i in allowed_states]):
        print 'State indices must be positive!'
        continue
      #if INFOS['diabatize']:
        #if len(allowed_states)>1:
          #print 'Only one initial state allowed!'
          #continue
      break
    INFOS['allowed']=set(allowed_states)
    if not 'erange' in INFOS:
      INFOS['erange']=[float('-inf'),float('inf')]

  if INFOS['read_QMout']:
    print centerstring('Considered states',60,'-')+'\n'
    print 'From which state should the excitation originate (for computation of excitation energies and oscillator strength)?'
    INFOS['initstate']=question('Lower state for excitation?',int,[1])[0]-1
  else:
    INFOS['initstate']=0

  if INFOS['excite']==3:
    if 'states' in INFOS:
      diago=(INFOS['repr']=='diag')
      print print_statemap(get_statemap(INFOS['states']),diag=diago)
    allstates=question('Do you want to include all states in the selection?',bool,True)
    if allstates:
      INFOS['allowed']=set()
    else:
      print '\nPlease enter the states which you want to EXCLUDE from the selection procedure.'
      a=question('Excluded states:',int,ranges=True)
      INFOS['allowed']=set([-i for i in a])
    print ''

  if INFOS['excite']==3:
    print centerstring('Random number seed',60,'-')+'\n'
    print 'Please enter a random number generator seed (type "!" to initialize the RNG from the system time).'
    while True:
      line=question('RNG Seed: ',str,'!',False)
      if line=='!':
        random.seed()
        break
      try:
        rngseed=int(line)
        random.seed(rngseed)
      except ValueError:
        print 'Please enter an integer or "!".'
        continue
      break
    print ''


  if INFOS['excite']==4:
    INFOS['allowed']=set()


  return INFOS

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def get_initconds(INFOS):
  ''''''

  print 'Reading initial condition file ...'
  if not INFOS['read_QMout'] and not INFOS['make_list']:
    INFOS['initf'].seek(0)
    while True:
      line=INFOS['initf'].readline()
      if 'Repr' in line:
        INFOS['diag']=line.split()[1].lower()=='diag'
        INFOS['repr']=line.split()[1]
      if 'Eref' in line:
        INFOS['eref']=float(line.split()[1])
        break

  initlist=[]
  width_bar=50
  for icond in range(1,INFOS['ninit']+1):
    initcond=INITCOND()
    initcond.init_from_file(INFOS['initf'],INFOS['eref'],icond)
    initlist.append(initcond)
    done=width_bar*(icond)/INFOS['ninit']
    sys.stdout.write('\r  Progress: ['+'='*done+' '*(width_bar-done)+'] %3i%%' % (done*100/width_bar))
  print '\nNumber of initial conditions in file:       %5i' % (INFOS['ninit'])
  return initlist

# ======================================================================================================================

def make_list(INFOS,initlist):
  print '\nMaking dummy states ...'
  width_bar=50
  for icond in range(1,INFOS['ninit']+1):
    estates=[]
    for istate in range(INFOS['nstates']):
      estates.append(STATE(i=istate+1))
    initlist[icond-1].addstates(estates)
    done=width_bar*(icond)/INFOS['ninit']
    sys.stdout.write('\r  Progress: ['+'='*done+' '*(width_bar-done)+'] %3i%%' % (done*100/width_bar))
  print '\nNumber of initial conditions where states were added:   %5i' % (INFOS['ninit'])
  return initlist

# ======================================================================================================================

def get_QMout(INFOS,initlist):
  ''''''

  print '\nReading QM.out data ...'
  if NONUMPY and  INFOS['diag']:
    print 'NUMPY not found, will use external SHARC diagonalizer...'
    global diagon
    diagon=diagonalizer()
  ncond=0
  initstate=INFOS['initstate']
  width_bar=50
  for icond in range(1,INFOS['ninit']+1):
    # look for a QM.out file
    qmfilename=INFOS['iconddir']+'/ICOND_%05i/QM.out' % (icond)
    done=width_bar*(icond)/INFOS['ninit']
    sys.stdout.write('\r  Progress: ['+'='*done+' '*(width_bar-done)+'] %3i%%' % (done*100/width_bar))
    if not os.path.isfile(qmfilename):
      #print 'No QM.out for ICOND_%05i!' % (icond)
      continue
    ncond+=1
    H,DM,P,Smat=extractQMout(qmfilename,INFOS['ion'],INFOS['diabatize'])
    if INFOS['diag']:
      H,DM,P=transform(H,DM,P)
    if INFOS['diabatize']:
      thres=0.5
      #string=''
      N=Smat[0][0].real**2
      for i in range(len(Smat)):
        for j in range(len(Smat[0])):
          Smat[i][j]=Smat[i][j].real**2/N
          #string+='%5.3f  ' % Smat[i][j]
        #string+='\n'
      #print string
      Diabmap={}
      for i in range(len(Smat)):
        j=Smat[i].index(max(Smat[i]))
        if Smat[i][j]>=thres:
          Diabmap[i]=j
      #print icond,Diabmap
    # generate list of excited states
    estates=[]
    for istate in range(len(H)):
      if INFOS['ion']:
        dip=[math.sqrt(abs(P[initstate][istate])),0,0]
      else:
        dip=[DM[i][initstate][istate] for i in range(3)]
      estate=STATE(len(estates)+1,  H[istate][istate],  H[initstate][initstate],   dip)
      estates.append(estate)
    initlist[icond-1].addstates(estates)
    if INFOS['diabatize']:
      initlist[icond-1].Diabmap=Diabmap
  print '\nNumber of initial conditions with QM.out:   %5i' % (ncond)
  return initlist

# ======================================================================================================================

def excite(INFOS,initlist):
  emin=INFOS['erange'][0]
  emax=INFOS['erange'][1]
  if not INFOS['excite']==4:
    if INFOS['excite']==3:
      # get the maximum oscillator strength
      maxprob=0
      for i,icond in enumerate(initlist):
        if icond.statelist==[]:
          continue
        for j,jstate in enumerate(icond.statelist):
          if emin <= jstate.Eexc <= emax:
            if -(j+1) not in INFOS['allowed']:
              if jstate.Prob>maxprob:
                maxprob=jstate.Prob
    # set the excitation flags
    print '\nSelecting initial states ...'
    width_bar=50
    nselected=0
    for i,icond in enumerate(initlist):
      done=width_bar*(i+1)/len(initlist)
      sys.stdout.write('\r  Progress: ['+'='*done+' '*(width_bar-done)+'] %3i%%' % (done*100/width_bar))
      if icond.statelist==[]:
        continue
      else:
        if INFOS['excite']==1:
          for jstate in icond.statelist:
            jstate.Excited=False
        elif INFOS['excite']==2:
          if INFOS['diabatize']:
            Diabmap=icond.Diabmap
            #print i,Diabmap
            allowed=[]
            for q in INFOS['allowed']:
              if q-1 in Diabmap:
                allowed.append(Diabmap[q-1]+1)
          else:
            allowed=INFOS['allowed']
          for j,jstate in enumerate(icond.statelist):
            if emin <= jstate.Eexc <= emax and j+1 in allowed:
              jstate.Excited=True
              nselected+=1
            else:
              jstate.Excited=False
        elif INFOS['excite']==3:
          # and excite
          for j,jstate in enumerate(icond.statelist):
            if emin <= jstate.Eexc <= emax:
              if -(j+1) not in INFOS['allowed']:
                jstate.Excite(maxprob,INFOS['erange'])
                if jstate.Excited:
                  nselected+=1
              else:
                jstate.Excited=False
            else:
              jstate.Excited=False
    print '\nNumber of initial states:                   %5i' % (nselected)

  # statistics
  maxprob=0.
  nexc=[0]
  ninrange=[0]
  ntotal=[0]
  for i,icond in enumerate(initlist):
    if icond.statelist==[]:
      continue
    else:
      for j,jstate in enumerate(icond.statelist):
        if j+1>len(ntotal):
          ntotal.append(0)
        if j+1>len(ninrange):
          ninrange.append(0)
        if j+1>len(nexc):
          nexc.append(0)
        ntotal[j]+=1
        if emin <= jstate.Eexc <= emax:
          ninrange[j]+=1
        if jstate.Excited:
          nexc[j]+=1
  print '\nNumber of initial conditions excited:'
  print 'State   Selected   InRange   Total'
  for i in range(len(ntotal)):
    print '  % 3i       % 4i      % 4i    % 4i' % (i+1,nexc[i],ninrange[i],ntotal[i])
  return initlist

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def writeoutput(initlist,INFOS):
  outfilename=INFOS['initf'].name+'.excited'
  if os.path.isfile(outfilename):
    overw=question('Overwrite %s? ' % (outfilename),bool,False)
    print ''
    if overw:
      try:
        outf=open(outfilename,'w')
      except IOError:
        print 'Could not open: %s' % (outfilename)
        outf=None
    else:
      outf=None
    if not outf:
      while True:
        outfilename=question('Please enter the output filename: ',str)
        try:
          outf=open(outfilename,'w')
        except IOError:
          print 'Could not open: %s' % (outfilename)
          continue
        break
  else:
    outf=open(outfilename,'w')

  print 'Writing output to %s ...' % (outfilename)

  string='''SHARC Initial conditions file, version %s   <Excited>
Ninit     %i
Natom     %i
Repr      %s
Eref      %18.10f
Eharm     %18.10f
''' % (version,INFOS['ninit'],INFOS['natom'],INFOS['repr'],INFOS['eref'],INFOS['eharm'])
  if INFOS['states']:
    string+='States    '
    for n in INFOS['states']:
      string+='%i ' % (n)
  string+='\n\n\nEquilibrium\n'

  for atom in INFOS['equi']:
    string+=str(atom)+'\n'
  string+='\n\n'

  for i,icond in enumerate(initlist):
    string+= 'Index     %i\n%s' % (i+1, str(icond))
  outf.write(string)
  outf.close()

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def main():
  '''Main routine'''

  usage='''
python excite.py

This interactive script reads out initconds files and QM.out files from excitation calculations and combines these
information to determine which initial conditions are bright enough for a dynamics simulation.
'''
  description=''
  parser = OptionParser(usage=usage, description=description)
  #parser.add_option('--no-excitation', dest='E', action='store_true',default=False,help="Sets all excitations to false.")
  #parser.add_option('--ground-state-only', dest='G', action='store_true',default=False,help="Selects the ground state of all initial conditions, and no excited states (e.g., for dynamics with laser excitation).")
  #(options, args) = parser.parse_args()

  displaywelcome()
  open_keystrokes()


  #INFOS={'do_excitations': not options.E, 'ground_state_only': options.G}
  INFOS={}
  INFOS=get_infos(INFOS)

  print '\n\n'+centerstring('Full input',60,'#')+'\n'
  for item in INFOS:
    if not item=='equi':
      print item, ' '*(25-len(item)), INFOS[item]
  print ''
  go_on=question('Do you want to continue?',bool,True)
  if not go_on:
    quit(0)
  print ''

  initlist=get_initconds(INFOS)

  if INFOS['read_QMout']:
    initlist=get_QMout(INFOS,initlist)
  if INFOS['make_list']:
    initlist=make_list(INFOS,initlist)
  initlist=excite(INFOS,initlist)

  if not INFOS['excite']==4:
    writeoutput(initlist,INFOS)
  else:
    print 'Nothing done, will not write output.'

  close_keystrokes()

# ======================================================================================================================

if __name__ == '__main__':
  try:
    main()
  except KeyboardInterrupt:
    print '\nCtrl+C makes me a sad SHARC ;-(\n'
    quit(0)
