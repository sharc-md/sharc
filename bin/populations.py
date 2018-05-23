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

# Interactive script for the setup of dynamics calculations for SHARC
# 
# usage: python setup_traj.py

import copy
import math
import sys
import re
import os
import stat
import shutil
import subprocess as sp
import datetime
import random
from optparse import OptionParser
import readline
import time
import colorsys
import pprint

try:
  import numpy
  NONUMPY=False
except ImportError:
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
BOHR_TO_ANG=0.529177211
AU_TO_FS=0.024188843
PI = math.pi

version='2.0'
versiondate=datetime.date(2018,2,1)


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
def itnmstates(states):

  x=0
  for i in range(len(states)):
    if states[i]<1:
      continue
    for k in range(i+1):
      for j in range(states[i]):
        x+=1
        yield i+1,j+1,k-i/2.,x
      x-=states[i]
    x+=states[i]
  return


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
    f.close()
  except IOError:
    print 'Could not write to file %s!' % (filename)
    sys.exit(13)

# ======================================================================= #
def mkdir(DIR):
    # mkdir the DIR, or clean it if it exists
    if os.path.exists(DIR):
        if os.path.isfile(DIR):
            print '%s exists and is a file!' % (DIR)
            sys.exit(69)
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
            sys.exit(70)

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

class output_dat:
  def __init__(self,filename):
    self.data=readfile(filename)
    self.filename=filename
    # get number of states
    for line in self.data:
      if 'nstates_m' in line:
        try:
          s=line.split()[0:-2]
          self.states=[ int(i) for i in s ]
        except ValueError:
          s=line.split()[1:-1]
          self.states=[ int(i) for i in s ]
        break
    nm=0
    for i,n in enumerate(self.states):
      nm+=n*(i+1)
    self.nmstates=nm
    # get line numbers where new timesteps start
    self.startlines=[]
    iline=-1
    while True:
      iline+=1
      if iline==len(self.data):
        break
      if 'Step' in self.data[iline]:
        self.startlines.append(iline)
    self.current=0
    #print self.states
    #print self.nmstates
    #print self.startlines
    #print self.current

  def __iter__(self):
    return self

  def next(self):
    # returns time step, U matrix and diagonal state
    # step
    current=self.current
    self.current+=1
    if current+1>len(self.startlines):
      raise StopIteration
    # U matrix starts at startlines[current]+5+nmstates
    U=[ [ 0 for i in range(self.nmstates) ] for j in range(self.nmstates) ]
    for iline in range(self.nmstates):
      index=self.startlines[current]+4+self.nmstates+iline
      line=self.data[index]
      s=line.split()
      for j in range(self.nmstates):
        r=float(s[2*j])
        i=float(s[2*j+1])
        U[iline][j]=complex(r,i)
    # diagonal state, has to search linearly
    while True:
      index+=1
      if index>len(self.data) or index==self.startlines[iline+1]:
        print 'Error reading timestep %i in file %s' % (current,self.filename)
        sys.exit(11)
      line=self.data[index]
      if 'states (diag, MCH)' in line:
        state_diag=int(self.data[index+1].split()[0])
        break
    return current,U,state_diag




















# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def centerstring(string,n,pad=' '):
  l=len(string)
  if l>=n:
    return string
  else:
    return  pad*((n-l+1)/2)+string+pad*((n-l)/2)

def displaywelcome():
  print 'Script for population computation started...\n'
  string='\n'
  string+='  '+'='*80+'\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Reading populations from SHARC dynamics',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Author: Sebastian Mai',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Version:'+version,80)+'||\n'
  string+='||'+centerstring(versiondate.strftime("%d.%m.%y"),80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  string+='''
This script calculates ensemble populations 
(e.g. based on the classically occupied state or based on the quantum amplitudes).
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
  shutil.move('KEYSTROKES.tmp','KEYSTROKES.populations')

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

class histogram:
  def __init__(self,binlist):
    '''binlist must be a list of floats
Later, all floats x with binlist[i-1]<x<=binlist[i] will return i'''
    self.binlist=sorted(binlist)
    self.len=len(binlist)+1
  def put(self,x):
    i=0
    for el in self.binlist:
      if x<=el:
        return i
      else:
        i+=1
    return i
  def __repr__(self):
    s='Histogram object: '
    for i in self.binlist:
      s+='%f ' % (i)
    return s

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def get_general():
  ''''''

  INFOS={}

  print centerstring('Paths to trajectories',60,'-')
  print '\nPlease enter the paths to all directories containing the "TRAJ_0XXXX" directories.\nE.g. Sing_2/ and Sing_3/. \nPlease enter one path at a time, and type "end" to finish the list.'
  count=0
  paths=[]
  while True:
    path=question('Path: ',str,'end')
    if path=='end':
      if len(paths)==0:
        print 'No path yet!'
        continue
      print ''
      break
    path=os.path.expanduser(os.path.expandvars(path))
    if not os.path.isdir(path):
      print 'Does not exist or is not a directory: %s' % (path)
      continue
    if path in paths:
      print 'Already included.'
      continue
    ls=os.listdir(path)
    print ls
    for i in ls:
      if 'TRAJ' in i:
        count+=1
    print 'Found %i subdirectories in total.\n' % count
    paths.append(path)
  INFOS['paths']=paths
  print 'Total number of subdirectories: %i\n' % (count)


  # get guessstates from SHARC input of first subdirectory
  ls=os.listdir(INFOS['paths'][0])
  for i in ls:
    if 'TRAJ' in i:
      break
  inputfilename=INFOS['paths'][0]+'/'+i+'/input'
  guessstates=None
  LD_dynamics=False
  if os.path.isfile(inputfilename):
    inputfile=open(inputfilename)
    for line in inputfile:
      if 'nstates' in line.lower():
        guessstates=[]
        l=re.sub('#.*$','',line).strip().split()
        for i in range(1,len(l)):
          guessstates.append(int(l[i]))
      if 'coupling' in line.lower():
        if 'overlap' in line.lower():
          LD_dynamics=True


  allowed=[i for i in range(1,12)]
  print centerstring('Analyze Mode',60,'-')
  print '''\nThis script can analyze the classical populations in different ways:
1       Number of trajectories in each diagonal state                                   from output.lis
2       Number of trajectories in each MCH state                                        from output.lis
3       Number of trajectories in each MCH state (multiplets summed up)                 from output.lis
4       Number of trajectories whose total spin value falls into certain intervals      from output.lis
5       Number of trajectories whose dipole moment falls into certain intervals         from output.lis
6       Number of trajectories whose oscillator strength falls into certain intervals   from output_data/fosc.out

It can also sum the quantum amplitudes:
7       Quantum amplitudes in diagonal picture                                          from output_data/coeff_diag.out
8       Quantum amplitudes in MCH picture                                               from output_data/coeff_MCH.out
9       Quantum amplitudes in MCH picture (multiplets summed up)                        from output_data/coeff_MCH.out

It can also transform the classical diagonal populations to MCH basis (might take long):
10      Transform diagonal populations to MCH states                                    from output.dat
11      Transform diagonal populations to MCH states (multiplets summed up)             from output.dat'''
  if LD_dynamics:
    print '20      Quantum amplitudes in diabatic picture                                          from output_data/coeff_diab.out'
    allowed.append(20)
  while True:
    num=question('Analyze mode:',int)[0]
    if not num in allowed:
      print 'Please enter one of the following integers: %s!' % (allowed)
      continue
    if guessstates!=None and len(guessstates)==1 and num==4:
      print 'Only singlet states, analysis unnecessary.'
      continue
    break
  INFOS['mode']=num
  print ''



  if INFOS['mode'] in [6,7,8,9,20]:
    print 'Run data_extractor.x for each trajectory prior to performing the analysis?\nFor many or long trajectories, this might take some time.'
    run_extractor=question('Run data_extractor.x?',bool,True)
    if run_extractor:
      run_full=not question('Run data_extractor.x only if output.dat newer than output_data/',bool,True)
    else:
      run_full=False
  else:
    run_extractor=False
    run_full=False
  INFOS['run_extractor']=run_extractor
  INFOS['run_extractor_full']=run_full



  if INFOS['mode'] in [1,2,3,7,8,9,20]:

    print centerstring('Number of states',60,'-')
    print '\nPlease enter the number of states as a list of integers\ne.g. 3 0 3 for three singlets, zero doublets and three triplets.'
    while True:
      states=question('Number of states:',int,guessstates)
      if len(states)==0:
        continue
      if any(i<0 for i in states):
        print 'Number of states must be positive!'
        continue
      break
    print ''
    nstates=0
    nmstates=0
    for mult,i in enumerate(states):
      nstates+=i
      nmstates+=(mult+1)*i
    INFOS['states']=states
    INFOS['nstates']=nstates
    INFOS['nmstates']=nmstates
    # obtain the statemap 
    statemap={}
    i=1
    for imult,istate,ims,instate in itnmstates(INFOS['states']):
      statemap[i]=[imult,istate,ims,instate]
      i+=1
    INFOS['statemap']=statemap


  if INFOS['mode'] in [4,5,6]:
    print centerstring('Intervals',60,'-')
    print '\nPlease enter the interval limits, all on one line.'
    if INFOS['mode'] in [4]:
      guess=[]
      for i in range(len(guessstates)-1):
        guess.append(i+0.5)
    else:
      guess=None
    limits=question('Interval limits: ',float,guess)
    INFOS['histo']=histogram(limits)
  print ''

  print centerstring('Normalization',60,'-')+'\n'
  INFOS['normalize']=question('Normalize the populations?',bool,True)
  print ''

  print centerstring('Simulation time',60,'-')
  print '\nUp to which simulation time should the analysis be performed? (Trajectories which are shorter are continued with their last values.)'
  while True:
    time=question('Simulation time (in fs): ',float,[1000.])[0]
    if time <0.:
      print 'Time must be positive!'
      continue
    break
  INFOS['maxtime']=time
  print ''


  print centerstring('Setup for bootstrapping?',60,'-')+'\n'
  print '\nThe population data can be analyzed by fitting with a kinetic model (via make_fitscript.py). In order to estimate errors for these time constants (via bootstrapping), additional data needs to be saved here.'
  INFOS['bootstrap']=question('Save data for bootstrapping?',bool,False)
  if INFOS['bootstrap']:
    INFOS['bootstrap_dir']=question('Directory for data?',str,'bootstrap_data/')


  print centerstring('Gnuplot script',60,'-')+'\n'
  INFOS['gnuplot']=question('Gnuplot script?',bool,False)
  if INFOS['gnuplot']:
    INFOS['gnuplot_out']=question('Gnuplot script filename?',str,'populations.gp')
  print ''

  return INFOS

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def do_calc(INFOS):

  forbidden=['crashed','running','dead','dont_analyze']

  #run the data extractor, if necessary
  if INFOS['run_extractor']:
    # first check whether $SHARC contains the exctractor
    print 'Running data_extractor...'
    sharcpath=os.getenv('SHARC')
    if sharcpath==None:
      print 'Please set $SHARC to the directory containing the SHARC executables!'
      sys.exit(1)
    else:
      if not os.path.isfile(sharcpath+'/data_extractor.x'):
        print '$SHARC does not contain data_extractor.x!'
        sys.exit(1)
      else:
        cwd=os.getcwd()
        for idir in INFOS['paths']:
          ls=os.listdir(idir)
          for itraj in ls:
            if not 'TRAJ_' in itraj:
              continue
            path=idir+'/'+itraj
            print path
            # check whether output_data/expec.out is newer than output.dat
            update=False
            if not os.path.isfile(path+'/output_data/expec.out'):
              update=True
            if not update:
              time_dat=os.path.getmtime(path+'/output.dat')
              time_expec=os.path.getmtime(path+'/output_data/expec.out')
              if time_dat > time_expec or INFOS['run_extractor_full']:
                update=True
            if update:
              os.chdir(path)
              io=sp.call(sharcpath+'/data_extractor.x -s -f output.dat > /dev/null 2> /dev/null',shell=True)
              if io!=0:
                print 'WARNING: extractor call failed for %s! Exit code %i' % (path,io)
              os.chdir(cwd)
            else:
              pass
    print 'Extraction finished!\n'

  width=30
  # prepare the list of output.lis files
  files=[]
  ntraj=0
  print 'Checking the directories...'
  for idir in INFOS['paths']:
    ls=os.listdir(idir)
    for itraj in ls:
      if not 'TRAJ_' in itraj:
        continue
      path=idir+'/'+itraj
      s=path+' '*(width-len(path))
      if INFOS['mode'] in [1,2,3,4,5]:
        pathfile=path+'/output.lis'
      elif INFOS['mode'] in [6]:
        pathfile=path+'/output_data/fosc.out'
      elif INFOS['mode'] in [7]:
        pathfile=path+'/output_data/coeff_diag.out'
      elif INFOS['mode'] in [8,9]:
        pathfile=path+'/output_data/coeff_MCH.out'
      elif INFOS['mode'] in [20]:
        pathfile=path+'/output_data/coeff_diab.out'
      elif INFOS['mode'] in [10,11]:
        pathfile=path+'/output.dat'
      if not os.path.isfile(pathfile):
        s+='%s NOT FOUND' % (pathfile)
        print s
        continue
      lstraj=os.listdir(path)
      valid=True
      for i in lstraj:
        if i.lower() in forbidden:
          s+='DETECTED FILE %s' % (i.lower())
          print s
          valid=False
          break
      if not valid:
        continue
      s+='OK'
      print s
      ntraj+=1
      files.append(pathfile)
  print 'Number of trajectories: %i' % (ntraj)
  if ntraj==0:
    print 'No valid trajectories found, exiting...'
    sys.exit(0)

  # get timestep
  if INFOS['mode'] in [1,2,3,4,5,6,7,8,9,20]:
    for ifile in files:
      lisf=open(ifile)
      file_valid=True
      while True:
        line=lisf.readline()
        if line=='':
          file_valid=False
          break
        if line[0]=='#':
          continue
        break
      if not file_valid:
        lisf.close()
        continue
      f=line.split()
      if INFOS['mode'] in [1,2,3,4,5]:
        t0=float(f[1])
      elif INFOS['mode'] in [6,7,8,9,20]:
        t0=float(f[0])
      N=0
      while True:
        line=lisf.readline()
        if len(line)==0:
          break
        if line[0]=='#':
          continue
        f=line.split()
        l2=line
        N+=1
      if N==0:
        lisf.close()
        continue
      f=l2.split()
      if INFOS['mode'] in [1,2,3,4,5]:
        dt=(float(f[1])-t0)/N
      elif INFOS['mode'] in [6,7,8,9,20]:
        dt=(float(f[0])-t0)/N
      if dt==0.:
        print 'ERROR: Timestep is zero.'
        quit(1)
      lisf.close()
      break
  elif INFOS['mode'] in [10,11]:
    for ifile in files:
      lisf=open(ifile)
      for line in lisf:
        if 'dtstep' in line:
          try:
            dt=float(line.split()[0])*AU_TO_FS
          except ValueError:
            dt=float(line.split()[-1])*AU_TO_FS
          break
      else:
        lisf.close()
        continue
      lisf.close()
      break

  # get number of steps
  nsteps=int(INFOS['maxtime']/dt)+1

  # get nstates
  if INFOS['mode'] in [1,2,7,8,20]:
    nstates=INFOS['nmstates']
  elif INFOS['mode'] in [3,9]:
    nstates=INFOS['nstates']
  elif INFOS['mode'] in [4,5,6]:
    nstates=len(INFOS['histo'].binlist)+1
  elif INFOS['mode'] in [10,11]:
    output_first=output_dat(files[0])
    INFOS['nmstates']=output_first.nmstates
    INFOS['states']=output_first.states
    nstates=0
    for i in INFOS['states']:
      nstates+=i
    # obtain the statemap 
    statemap={}
    i=1
    for imult,istate,ims,instate in itnmstates(INFOS['states']):
      statemap[i]=[imult,istate,ims,instate]
      i+=1
    INFOS['statemap']=statemap
  print 'Found dt=%f, nsteps=%i, nstates=%i\n' % (dt,nsteps,nstates)
  INFOS['nstates']=nstates

  # get populations
  width=60
  pop_full=[ [ [0. for j in range(nstates) ] for i in range(nsteps) ] for ifile in files ]
  shortest=9999999.
  longest=0.
  for fileindex,ifile in enumerate(files):
    if INFOS['mode'] in [10,11]:
      output_current=output_dat(ifile)
      istep=-1
      for istep,U,state_diag in output_current:
        #print istep,state_diag
        vec2=[ U[i][state_diag-1] for i in range(len(U)) ]
        vec=[ 0. for i in range(nstates)]
        if INFOS['mode'] in [10]:
          for i in range(nstates):
            vec[i]=vec2[i].real**2+vec2[i].imag**2
        elif INFOS['mode'] in [11]:
          for i in range(INFOS['nmstates']):
            state=INFOS['statemap'][i+1][3]-1
            vec[state]+=vec2[i].real**2+vec2[i].imag**2
        for istate in range(nstates):
          pop_full[fileindex][istep][istate]+=vec[istate]
      if dt*istep<shortest:
        shortest=dt*istep
      if dt*istep>longest:
        longest=dt*istep
      if istep==-1:
        print '%s' % (ifile)+' '*(width-len(ifile))+' %i\tZero Timesteps found!' % (t)
        ntraj-=1
        continue
      else:
        print '%s' % (ifile)+' '*(width-len(ifile))+' %i' % (istep)
      while istep+1<nsteps:
        istep+=1
        if INFOS['mode'] in [10,11]:
          for i in range(nstates):
            pop_full[fileindex][istep][i]+=vec[i]
    else:
      lisf=open(ifile)
      t=-1
      for line in lisf:
        if line[0]=='#':
          continue
        f=line.split()
        t+=1
        if t>=nsteps:
          break

        if INFOS['mode'] in [1,2,3,4,5,6]:
          if INFOS['mode']==1:
            state=int(f[2])-1
          elif INFOS['mode']==2:
            state=int(f[3])-1
          elif INFOS['mode']==3:
            state=int(f[3])
            # state in nm scheme to state in n scheme
            state=INFOS['statemap'][state][3]-1
          elif INFOS['mode']==4:
            state=INFOS['histo'].put(float(f[9]))
          elif INFOS['mode']==5:
            state=INFOS['histo'].put(float(f[8]))
          elif INFOS['mode']==6:
            state=INFOS['histo'].put(float(f[1]))
          pop_full[fileindex][t][state]+=1
        elif INFOS['mode'] in [7,8,9,20]:
          vec=[ 0. for i in range(nstates)]
          if INFOS['mode'] in [7,8,20]:
            for i in range(nstates):
              vec[i]=float(f[2+2*i])**2+float(f[3+2*i])**2
          if INFOS['mode']==9:
            for i in range(INFOS['nmstates']):
              state=INFOS['statemap'][i+1][3]-1
              #imult,istate,ims=IstateToMultState(i+1,INFOS['states'])
              #state=MultStateToIstate(imult,istate,INFOS['states'])-1
              vec[state]+=float(f[2+2*i])**2+float(f[3+2*i])**2
          for i in range(nstates):
            pop_full[fileindex][t][i]+=vec[i]
      lisf.close()
      if dt*t<shortest:
        shortest=dt*t
      if dt*t>longest:
        longest=dt*t
      if t==-1:
        print '%s' % (ifile)+' '*(width-len(ifile))+'%i\tZero Timesteps found!' % (t)
        ntraj-=1
        continue
      else:
        print '%s' % (ifile)+' '*(width-len(ifile))+'%i' % (t)
      while t+1<nsteps:
        t+=1
        if INFOS['mode'] in [1,2,3,4,5,6]:
          pop_full[fileindex][t][state]+=1
        elif INFOS['mode'] in [7,8,9,20]:
          for i in range(nstates):
            pop_full[fileindex][t][i]+=vec[i]
  print 'Shortest trajectory: %f' % (shortest)
  print 'Longest trajectory: %f' % (longest)
  print 'Number of trajectories: %i' % (ntraj)
  INFOS['shortest']=shortest
  INFOS['longest']=longest

  # make pop array
  pop=[ [0. for j in range(nstates) ] for i in range(nsteps) ]        # first index is time, second is state
  for fileindex,ifile in enumerate(files):
    for i in range(nsteps):
      for j in range(nstates):
        pop[i][j]+=pop_full[fileindex][i][j]

  # write populations
  s='#Mode: %i\n' % INFOS['mode']
  s+='#%15i ' % (1)
  for i in range(nstates):
    s+='%16i ' % (i+2)
  s+='\n'
  s+='#%15s ' % ('Time (fs)')
  for i in range(nstates):

    if INFOS['mode'] in [1,7]:
      s+='%16s ' % ('X%i' % (i+1))
    elif INFOS['mode'] in [2,8,20,10]:
      mult,state,ms=tuple(INFOS['statemap'][i+1][0:3])
      #IstateToMultState(i+1,INFOS['states'])
      string='%s %i %i' % (IToMult[mult][0:3],state,ms)
      s+='%16s ' % (string)
    elif INFOS['mode'] in [3,9,11]:
      mult,state=tuple(INFOS['statemap'][i+1][0:2])
      #INstateToMultState(i+1,INFOS['states'])
      string='%s %i' % (IToMult[mult][0:3],state)
      s+='%16s ' % (string)
    elif INFOS['mode'] in [4,5,6]:
      if i<len(INFOS['histo'].binlist):
        string='< %.2e' % (INFOS['histo'].binlist[i])
      else:
        string='> %.2e' % (INFOS['histo'].binlist[-1])
      s+='%16s ' % (string)

  s+='\n'
  for i,line in enumerate(pop):
    s+='%16.9f ' % (i*dt)
    for el in line:
      if INFOS['normalize']:
        x=float(el)/ntraj
      else:
        x=float(el)
      s+='%16.9f ' % (x)
    s+='\n'
  #print s

  print ''
  outfilename='pop.out'
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

  print 'Writing to %s ...' % (outfilename)
  outf.write(s)
  outf.close()


  # save bootstrapping data
  if INFOS['bootstrap']:
    print 'Writing to %s ...' % (INFOS['bootstrap_dir'])
    mkdir(INFOS['bootstrap_dir'])
    for fileindex,ifile in enumerate(files):
      filename=os.path.join(INFOS['bootstrap_dir'],'pop_%i.dat' % (fileindex+1))
      s='#Mode: %i\n' % INFOS['mode']
      s+='#%15i ' % (1)
      for i in range(nstates):
        s+='%16i ' % (i+2)
      s+='\n'
      s+='#%15s ' % ('Time (fs)')
      for i in range(nstates):

        if INFOS['mode'] in [1,7]:
          s+='%16s ' % ('X%i' % (i+1))
        elif INFOS['mode'] in [2,8,20,10]:
          mult,state,ms=tuple(INFOS['statemap'][i+1][0:3])
          #IstateToMultState(i+1,INFOS['states'])
          string='%s %i %i' % (IToMult[mult][0:3],state,ms)
          s+='%16s ' % (string)
        elif INFOS['mode'] in [3,9,11]:
          mult,state=tuple(INFOS['statemap'][i+1][0:2])
          #INstateToMultState(i+1,INFOS['states'])
          string='%s %i' % (IToMult[mult][0:3],state)
          s+='%16s ' % (string)
        elif INFOS['mode'] in [4,5,6]:
          if i<len(INFOS['histo'].binlist):
            string='< %.2e' % (INFOS['histo'].binlist[i])
          else:
            string='> %.2e' % (INFOS['histo'].binlist[-1])
          s+='%16s ' % (string)

      s+='\n'
      for i,line in enumerate(pop_full[fileindex]):
        s+='%16.9f ' % (i*dt)
        for el in line:
          s+='%16.9f ' % (float(el))
        s+='\n'
      writefile(filename,s)
    #for i in range(nsteps):






  INFOS['outputfile']=outfilename
  INFOS['ntraj']=ntraj
  return INFOS

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
    print index,el,R.hexcolor(index+1,el+1)
print 2,1,R.hexcolor(2,1)

                            .--------- running index over the groups
                            | .------- number of groups
                            | | .----- running index over elements
                            | | | .--- number of groups
                            | | | | .- number of elements in group
                            | | | | |
Output:                     v v v v v
1     1         #FF0000  =  0/2+0/2/6
1     2         #FF7F00  =  0/2+1/2/6
1     3         #FFFF00  =  0/2+2/2/6
1     4         #7FFF00  =  0/2+3/2/6
1     5         #00FF00  =  0/2+4/2/6
1     6         #00FF7F  =  0/2+5/2/6
3     1         #00FFFF  =  1/2+1/2/3
3     2         #0000FF  =  1/2+2/2/3
3     2         #FF00FF  =  1/2+3/2/3
2     1         #FFFFFF #invalid, hence white
'''
  def __init__(self,initlist):
    # Hues:
    # 0.0     0.15       ...
    # Red     Yellow     ...
    excluded=[
      [0.12,0.22]       # exclude yellow hues from the colorwheel
      ]
    # excluded-list must be sorted, for each element x[0]<=x[1]
    # each excluded[i][0]<=excluded[i+1][0]
    # everything between 0 and 1
    # sort the pairs
    for i,el in enumerate(excluded):
      excluded[i]=[min(el),max(el)]
    # sort the list
    excluded.sort(key=lambda x:x[0])
    #make all elements between 0 and 1
    temp1=[]
    for i,el in enumerate(excluded):
      temp1.append([ min(1.,max(0.,el[0])), max(0.,min(1.,el[1])) ])
    # resolve overlapping ranges
    temp2=[[0.,0.]]
    for i,el in enumerate(temp1):
      if el[0]>=temp2[-1][1]:
        temp2.append(el)
      else:
        temp2[-1][1]=el[1]
    self.excluded=temp2
    # number of non-empty groups
    self.initlist=initlist
    n=0
    for index,el in enumerate(initlist):
      # negative numbers in initlist are not allowed, make them to zero
      self.initlist[index]=max(0,el)
      if el>0:
        n+=1
    self.n=n                    # number of non-empty groups
    self.m=len(initlist)        # number of groups
    # available colorspace self.a and shifts to skip excluded regions self.ex
    self.a=1.
    self.ex=[ 0. for i in self.excluded ]
    for i,el in enumerate(self.excluded):
      self.a-=el[1]-el[0]
      self.ex[i]=el[1]-el[0]
    # list of starting values
    self.startlist=[ 0. for i in range(self.m) ]
    for index,el in enumerate(initlist):
      if index==0:
        continue
      if el>0:
        self.startlist[index]=self.a/self.n+self.startlist[index-1]
      else:
        self.startlist[index]=self.startlist[index-1]
    # list of increments
    self.incrlist=[ 0. for i in range(self.m) ]
    for index,el in enumerate(initlist):
      if el>0:
        self.incrlist[index]=self.a/self.n/el
  def rgb_to_hex(self,rgb):
    #rgb is a list with three elements, which are floats 0<=x<=1
    color=0
    for i in range(3):
      a=max(min(rgb[i],1.0),0.0)
      color+=int(255*a)*256**(2-i)
    #color=int(255*triple[0])*256**2+int(255*triple[1])*256+int(255*triple[2])
    string=hex(color)[2:].upper()
    string='#'+'0'*(6-len(string))+string
    return string
  def hexcolor(self,index,el):
    if not 1<=index<=self.m:
      return '#FFFFFF'
    if not 1<=el<=self.initlist[index-1]:
      return '#FFFFFF'
    deg=self.startlist[index-1]+self.incrlist[index-1]*(el-1)
    for i,el in enumerate(self.excluded):
      if deg>el[0]:
        deg+=self.ex[i]
    rgbtriple=colorsys.hsv_to_rgb(deg,1,1)
    return self.rgb_to_hex(rgbtriple)




# =============================================

def make_gnuplot(INFOS):
  
  #if not INFOS['diag'] and 'states' in INFOS:
    #sortedMCH=True
    #statemap={}
    #i=1
    #for imult,istate,ims,instate in itnmstates(INFOS['states']):
      #statemap[i]=[imult,istate,ims,instate]
      #i+=1
    #nstates=sum(INFOS['states'])
    #nmstates=INFOS['nstate']
    #maxmult=len(INFOS['states'])
    #angdiff=1./maxmult/max(INFOS['states'])
  #else:
    #sortedMCH=False
    #nstates=INFOS['nstate']
    #nmstates=INFOS['nstate']
    #angdiff=math.ceil(nstates/2+0.5)/nstates

  # calculate maxmult and angdiff
  if 'states' in INFOS:
    if INFOS['mode'] in [1,7]:
      init=[INFOS['nmstates']]
      #maxmult=1
      #angdiff=1./maxmult/max(INFOS['states'])
    else:
      init=INFOS['states']
      #maxmult=len(INFOS['states'])
      #angdiff=1./maxmult/max(INFOS['states'])
  elif 'histo' in INFOS:
    init=[INFOS['histo'].len]
    #maxmult=INFOS['histo'].len
    #angdiff=1./maxmult
  nstates=INFOS['nstates']
  R=rgbcolor(init)


  # gnustring
  title={1: 'Classical populations (diagonal)',
         2: 'Classical populations (MCH)',
         3: 'Classical populations (MCH, multiplets)',
         4: 'Populations (spin expectation value)',
         5: 'Populations (state dipole moment)',
         6: 'Populations (oscillator strength)',
         7: 'Quantum amplitudes (diagonal)',
         8: 'Quantum amplitudes (MCH)',
         9: 'Quantum amplitudes (MCH, multiplets)',
        10: 'Transformed classical populations (MCH)',
        11: 'Transformed classical populations (MCH)',
        20: 'Quantum amplitudes (diabatic)'
        }

  gnustring='''set title "%s\\n%i Trajectories (Shortest %.1f fs, Longest %.1f fs)"

set xrange [%f:%f]
set yrange [%f:%f]
set xlabel 'Time (fs)'
set ylabel 'Population'

set term pngcairo size 640,480
set out '%s.png'
''' % (title[INFOS['mode']],
       INFOS['ntraj'],
       INFOS['shortest'],
       INFOS['longest'],
       0.,
       INFOS['maxtime'],
       0.,
       1.,
       INFOS['gnuplot_out']
       )

  for istate in range(1,1+nstates):
    if istate==1:
      gnustring+='plot "%s" ' % (INFOS['outputfile'])
    else:
      gnustring+='     ""        '

    if INFOS['mode'] in [1,7]:
      gnustring+='u 1:%i w l tit "State %i" lw 2.5 lc rgbcolor "%s"' % (istate+1,istate,R.hexcolor(1,istate) )
    elif INFOS['mode'] in [2,8,10,20]:
      mult,state,ms=tuple(INFOS['statemap'][istate][0:3])
      name=IToMult[mult]+' %i' % (state-(mult==1 or mult==2))
      gnustring+='u 1:%i w l tit "%s" lw 2.5 lc rgbcolor "%s"' % (istate+1,name,R.hexcolor(mult,state))
    elif INFOS['mode'] in [3,9,11]:
      gnustring+='u 1:%i' % (istate+1)
      for i in INFOS['statemap']:
        if istate==INFOS['statemap'][i][3]:
          mult=INFOS['statemap'][i][0]
          state=INFOS['statemap'][i][1]
      name=IToMult[mult]+' %i' % (state-(mult==1 or mult==2))
      gnustring+=' w l tit "%s" lw 2.5 lc rgbcolor "%s"' % (
      name,
      R.hexcolor(mult,state)
      )
    elif INFOS['mode'] in [4,5,6]:
      if istate-1<len(INFOS['histo'].binlist):
        name='< %.2e' % (INFOS['histo'].binlist[istate-1])
      else:
        name='> %.2e' % (INFOS['histo'].binlist[-1])
      gnustring+='u 1:%i w l tit "%s" lw 2.5 lc rgbcolor "%s"' % (istate+1,name,R.hexcolor(1,istate))

    if istate!=nstates:
      gnustring+=', \\\n'
  gnustring+='\n'

  out=open(INFOS['gnuplot_out'],'w')
  out.write(gnustring)
  out.close()
  print 'Gnuplot script written to "%s"' % (INFOS['gnuplot_out'])


# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def main():
  '''Main routine'''

  usage='''
python populations.py

This interactive program reads information from output.lis files or output_data/ and calculates ensemble populations.
'''

  description=''
  displaywelcome()
  open_keystrokes()

  INFOS=get_general()

  print centerstring('Full input',60,'#')+'\n'
  for item in INFOS:
    print item, ' '*(25-len(item)), INFOS[item]
  print ''
  calc=question('Do you want to do the specified analysis?',bool,True)
  print ''

  if calc:
    INFOS=do_calc(INFOS)

    if INFOS['gnuplot']:
      make_gnuplot(INFOS)

  close_keystrokes()


# ======================================================================================================================

if __name__ == '__main__':
  try:
    main()
  except KeyboardInterrupt:
    print '\nCtrl+C makes me a sad SHARC ;-(\n'
    quit(0)
