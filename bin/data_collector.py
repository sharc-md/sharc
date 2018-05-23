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

PIPELINE='''
**************************************************************************************************
Pipeline:
=========
            Collecting
                 |
                 1-------\ 
                 |       |
                 |   Smoothing
                 |       |
                 |<------/
                 |
  /--------------2--------\ 
  |                       |
  |                 Synchronizing
  |                       |
  |          /------------3--------------\ 
  |          |                           |
  |          |                     Convoluting(X)
  |          |                           |
  |          4-------\                   6---------\ 
  |          |       |                   |         |
  |          |   Averaging               |     Summing(Y)
  |          |       |                   |         |
  |          |<------/                   |<--------/
  |          |                           |
  |          5-------\                   7---------\ 
  |          |       |                   |         |
  |          |   Statistics              |   Integrating(X)
  |          |       |                   |         |
  |          |<------/                   |<--------/
  |          |                           |
  |          |                           8---------\ 
  |          |                           |         |
  |          |               /-----------9         |
  |          |               |           |         |
  |          |         Integrating(T)    |   Convoluting(T)
  |          |               |           |         |
  |          |               \---------->|         |
  |          |                           |         |
  |          |                           |<--------/
  |          |                           |
  |          |<-------------------------10
  |          |                           |
Type1      Type2                       Type3

Procedure Explanations:
=======================
- Collecting:
  extracts the requested columns from the files, creating a Type1 dataset

- Smoothing:
  applies some smoothing procedure to each trajectory independently, creating a new Type1 dataset

- Synchronizing:
  merges all trajectories into one dataset with a common time axis, creating a Type2 dataset

- Averaging:
  for each X and Y column, compute average and standard deviation across all trajectories, creating a new Type2 dataset

- Statistics:
  for each X and Y column, compute average and standard deviation across all trajectories and all time steps, creating a new Type2 dataset

- Convoluting(X):
  merges all trajectories into a dataset with common time and X axes, using a convolution along the X axis; creates a Type3 dataset

- Sum(Y):
  if multiple Y values are present, sum them up for each T and X, creating a new Type3 dataset

- Integrating(X):
  integrates along the X axis within some specified bounds, creating a new Type3 dataset

- Convoluting(T):
  applies a convolution along the time axis, creating a new Type3 dataset

- Integrating(T):
  performs a cumulative summation along the time axis (such that the final integral is stored in the last time step), creating a new Type3 dataset

Dataset Explanations:
====================
- Type1 dataset:
  Independent trajectories with possibly different time axes
  (not intended for plotting)
  ***
##i  path                 time   x1   x2  ...   y1   y1  ...
  0  TRAJ_00001/filename   0.0  1.6  3.1  ...  0.1  6.1  ...
  0  TRAJ_00001/filename   0.5  1.7  2.7  ...  0.1  6.0  ...
  0  TRAJ_00001/filename   1.0  1.9  2.2  ...  0.2  6.1  ...
  ...

  1  TRAJ_00002/filename   0.0  1.7  2.9  ...  0.2  6.2  ...
  1  TRAJ_00002/filename   1.0  1.8  2.3  ...  0.1  6.3  ...
  1  TRAJ_00002/filename   2.0  1.7  1.6  ...  0.1  6.1  ...
  ...
  ***


- Type2 dataset:
  Data with a common time axis, with possibly missing entries for some trajectories
  (can be plotted as hair figures, etc)
  ***
##       <-- TRAJ_00001 -->  ...  <-- TRAJ_00002 -->  ...
##time   x1   y1   x2   y2  ...   x1   y1   x2   y2  ...
  0.0    1.6  0.1  3.1  6.1  ...  1.7  0.2  2.9  6.2  ...
  0.5    1.7  0.1  2.7  6.0  ...  nan  nan  nan  nan  ...
  1.0    1.9  0.2  2.2  6.1  ...  1.8  0.1  2.3  6.3  ...
  ...
  ***


- Type3 dataset:
  Common time and X axes, Y values obtained by convolution
  (can be plotted as 3D plots)
  ***
##time   x   y1  y2   ...
  0.0  1.2  0.2  0.0  ...
  0.0  1.3  0.3  0.0  ...
  0.0  1.4  0.5  0.0  ...
  0.0  1.5  0.7  0.0  ...
  ...

  0.5  1.2  0.2  0.0  ...
  0.5  1.3  0.3  0.0  ...
  0.5  1.4  0.5  0.0  ...
  0.5  1.5  0.7  0.0  ...
  ...
  ***
**************************************************************************************************
'''





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
class gauss:
  def __init__(self,fwhm):
    self.fwhm=fwhm
    self.c=-4.*math.log(2.)/fwhm**2             # this factor needs to be evaluated only once
  def ev(self,A,x0,x):
    return A*math.exp( self.c*(x-x0)**2)        # this routine does only the necessary calculations

class lorentz:
  def __init__(self,fwhm):
    self.fwhm=fwhm
    self.c=0.25*fwhm**2
  def ev(self,A,x0,x):
    return A/( (x-x0)**2/self.c+1)

class boxfunction:
  def __init__(self,fwhm):
    self.fwhm=fwhm
    self.w=0.5*fwhm
  def ev(self,A,x0,x):
    if abs(x-x0)<self.w:
      return A
    else:
      return 0.

class lognormal:
  def __init__(self,fwhm):
    self.fwhm=fwhm
  def ev(self,A,x0,x):
    if x<=0 or x0<=0:
      return 0.
    # for lognormal distribution, the factor for the exponent depends on x0
    c=(math.log( (self.fwhm+math.sqrt(self.fwhm**2+4.*x0**2))/(2.*x0)))**2
    # note that the function does not take a value of A at x0
    # instead, the function is normalized such that its maximum will have a value of A (at x<=x0)
    return A*x0/x*math.exp( -c/(4.*math.log(2.)) -math.log(2.)*(math.log(x)-math.log(x0))**2/c)

kernels={1: {'f': gauss,       'description': 'Gaussian function',          'factor':1.5},
         2: {'f': lorentz,     'description': 'Lorentzian function',        'factor':2.5},
         3: {'f': boxfunction, 'description': 'Rectangular window function','factor':0.6},
         4: {'f': lognormal,   'description': 'Log-normal function',        'factor':1.5}}

class spectrum:
  def __init__(self,npts,emin,emax,fwhm,lineshape):
    self.npts=npts
    if lineshape==1:
      self.f=gauss(fwhm)
    elif lineshape==2:
      self.f=lorentz(fwhm)
    elif lineshape==3:
      self.f=boxfunction(fwhm)
    elif lineshape==4:
      self.f=lognormal(fwhm)
    self.en=[ emin + float(i)/self.npts*(emax-emin) for i in range(self.npts+1) ]       # the energy grid needs to be calculated only once
    self.spec=[ 0. for i in range(self.npts+1) ]
  def add(self,A,x0):
    if A==0.:
      return
    for i in range(self.npts+1):
      self.spec[i]+=self.f.ev(A,x0,self.en[i])

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
  print 'Script for data collecting started...\n'
  string='\n'
  string+='  '+'='*80+'\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Reading table data from SHARC dynamics',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Author: Sebastian Mai',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Version:'+version,80)+'||\n'
  string+='||'+centerstring(versiondate.strftime("%d.%m.%y"),80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  string+='''
This script collects table data from SHARC trajectories, smooths them, synchronizes them,
convolutes them, and computes averages and similar statistics.
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
  shutil.move('KEYSTROKES.tmp','KEYSTROKES.data_collector')

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

# ======================================================================= #         OK
def containsstring(string,line):

    a=re.search(string,line)
    if a:
        return True
    else:
        return False

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def get_general():
  ''''''

  INFOS={}

  # -------------------------------- Running data_extractor or geo.py ----------------------------
  # TODO

  # ---------------------------------------- File selection --------------------------------------

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
      if 'TRAJ' in i or 'ICOND' in i:
        count+=1
    print 'Found %i subdirectories in total.\n' % count
    paths.append(path)
  INFOS['paths']=paths
  print 'Total number of subdirectories: %i\n' % (count)

  # make list of TRAJ paths
  width=50
  forbidden=['crashed','running','dead','dont_analyze']
  dirs=[]
  ntraj=0
  print 'Checking the directories...'
  for idir in INFOS['paths']:
    ls=os.listdir(idir)
    for itraj in sorted(ls):
      if not 'TRAJ_' in itraj and not 'ICOND_' in itraj:
        continue
      path=idir+'/'+itraj
      if not os.path.isdir(path):
        continue
      s=path+' '*(width-len(path))
      lstraj=os.listdir(path)
      valid=True
      for i in lstraj:
        if i.lower() in forbidden:
          s+='DETECTED FILE %s' % (i.lower())
          #print s
          valid=False
          break
      if not valid:
        continue
      s+='OK'
      #print s
      ntraj+=1
      dirs.append(path)
  print 'Number of trajectories: %i' % (ntraj)
  if ntraj==0:
    print 'No valid trajectories found, exiting...'
    sys.exit(0)

  # check the dirs
  print 'Checking for common files...'
  allfiles={}
  for d in dirs:
    for dirpath, dirnames, filenames in os.walk(d):
      for f in filenames:
        line=os.path.join( os.path.relpath(dirpath,d), f)
        if line in allfiles:
          allfiles[line]+=1
        else:
          allfiles[line]=1
  exclude=['QM/ADF\.template',
           'QM/ADF\.resources',
           'QM/Analytical\.template',
           'QM/COLUMBUS\.resources',
           'QM/GAUSSIAN\.template',
           'QM/GAUSSIAN\.resources',
           'QM/LCV\.template',
           'QM/MOLPRO\.template',
           'QM/MOLPRO\.resources',
           'QM/MOLCAS\.template',
           'QM/MOLCAS\.resources',
           'QM/RICC2\.template',
           'QM/RICC2\.resources',
           'QM/.*init',
           'QM/SCRATCH',
           'QM/SAVE',
           'QM/runQM\.sh',
           'QM/QM\.in',
           'QM/QM\.out',
           'QM/QM\.log',
           'QM/QM\.err',
           '\./output\.dat',
           '\./output\.log',
           '\./output\.xyz',
           '\./output\.dat\.ext',
           '\./input',
           '\./geom',
           '\./veloc',
           '\./coeff',
           '\./restart\.traj',
           '\./restart\.ctrl',
           '\./run\.sh',
           'restart/.*',
           '\./.*init',
           '\./STOP',
           '\./CRASHED',
           '\./RUNNING',
           '\./DONT_ANALYZE'
           ]
  allfiles2={}
  for line in allfiles:
    if allfiles[line]>=2 and not any( [containsstring(i,line) for i in exclude] ):
      allfiles2[line]=allfiles[line]
  print '\nList of files common to the trajectory directories:\n'
  print '%6s %20s   %s' % ('Index','Number of appearance','Relative file path')
  print '-'*58
  allfiles_index={}
  for iline,line in enumerate(sorted(allfiles2)):
    allfiles_index[iline]=line
    print '%6i %20i   %s' % (iline,allfiles2[line],line)

  # choose one of these files
  print '\nPlease give the relative file path of the file you want to collect:'
  while True:
    string=question('File path or index:',str,'0',False)
    try:
      string=allfiles_index[int(string)]
    except ValueError:
      pass
    if string in allfiles2:
      INFOS['filepath']=string
      break
    else:
      print 'I did not understand %s' % string

  # make list of files
  allfiles=[]
  for d in dirs:
    f=os.path.join(d,INFOS['filepath'])
    if os.path.isfile(f):
      allfiles.append( f )
  INFOS['allfiles']=allfiles

  # ---------------------------------------- Columns --------------------------------------

  print '\n'+centerstring('Data columns',60,'-')+'\n'
  # get number of columns
  filename=allfiles[0]
  testfile=readfile(filename)
  for line in testfile:
    if not '#' in line:
      ncol=len(line.split())
      break
  print 'Number of columns in the file:   %i' % (ncol)
  INFOS['ncol']=ncol

  # select columns
  print '\nPlease select the data columns for the analysis:'
  print 'For T column: \n  only enter one (positive) column index. \n  If 0, the line number will be used instead.'
  print 'For X column: \n  enter one or more column indices. \n  If 0, all entries of that column will be set to 1. \n  If negative, the read numbers will be multiplied by -1.'
  print 'For Y column: \n  enter as many column indices as for X. \n  If 0, all entries of that column will be set to 1. \n  If negative, the read numbers will be multiplied by -1.'
  print ''
  while True:
    INFOS['colT']=question('T column (time):',int,[1])[0]
    if 0<=INFOS['colT']<=ncol:
      # 0:   use line number (neglecting commented or too short lines)
      # 1-n: use that line for time data
      break
    else:
      print 'Please enter a number between 0 and %i!' % ncol
  while True:
    INFOS['colX']=question('X columns:',int,[2],ranges=True)
    if all( [-ncol<=x<=ncol for x in INFOS['colX'] ] ):
      INFOS['nX']=len(INFOS['colX'])
      break
    else:
      print 'Please enter a set of numbers between %i and %i!' % (-ncol,ncol)
  while True:
    default=[0 for i in INFOS['colX']]
    INFOS['colY']=question('Y columns:',int,default,ranges=True)
    if all( [-ncol<=x<=ncol for x in INFOS['colY'] ] ) and len(INFOS['colY'])==len(INFOS['colX']):
      INFOS['nY']=len(INFOS['colY'])
      break
    else:
      print 'Please enter a set of %i numbers between %i and %i!' % (len(INFOS['colX']),-ncol,ncol)

  print 'Selected columns:'
  print 'T: %s     X: %s    Y: %s\n' % (str(INFOS['colT']),str(INFOS['colX']),str(INFOS['colY']))

  # ---------------------------------------- Analysis procedure --------------------------------------

  print centerstring('Analysis procedure',60,'-')+'\n'
  show=question('Show possible workflow options?',bool,True)
  if show:
    print '\nThe following diagram shows which workflows are possible with this script:'
    print PIPELINE

  # Question 1
  print '\n'+centerstring('1 Smoothing',40,'-')+'\n'
  if question('Do you want to apply smoothing to the individual trajectories?',bool,False):
    print '\nChoose one of the following smoothing functions:'
    for i in sorted(kernels):
      print '%i  %s' % (i,kernels[i]['description'])
    while True:
      i=question('Choose one of the functions:',int,[1])[0]
      if i in kernels:
        break
      else:
        print 'Choose one of the following: %s' % (list(kernels))
    w=question('Choose width of the smoothing function (in units of column %i):' % (INFOS['colT']),float,[10.0])[0]
    INFOS['smoothing']={'function': kernels[i]['f'](w)}
  else:
    INFOS['smoothing']={}

  # Question 2
  print '\n'+centerstring('2 Synchronizing',40,'-')+'\n'
  if question('Do you want to synchronize the data?',bool,True):
    INFOS['synchronizing']=True
  else:
    INFOS['synchronizing']=False

  # first branching
  INFOS['averaging']={}
  INFOS['statistics']={}
  INFOS['convolute_X']={}
  INFOS['sum_Y']=False
  INFOS['integrate_X']={}
  INFOS['convolute_T']={}
  INFOS['integrate_T']={}
  INFOS['type3_to_type2']=False

  # Question 3
  if INFOS['synchronizing']:
    print '\n'+centerstring('3 Convoluting along X',40,'-')+'\n'
    if question('Do you want to apply convolution in X direction?',bool,False):
      print '\nChoose one of the following convolution kernels:'
      for i in sorted(kernels):
        print '%i  %s' % (i,kernels[i]['description'])
      while True:
        kern=question('Choose one of the functions:',int,[1])[0]
        if kern in kernels:
          break
        else:
          print 'Choose one of the following: %s' % (list(kernels))
      w=question('Choose width of the smoothing function (in units of the X columns):',float,[1.0])[0]
      INFOS['convolute_X']={'function': kernels[kern]['f'](w)}
      #print 'Choose the size of the grid along X:'
      INFOS['convolute_X']['npoints']=question('Size of the grid along X:',int,[25])[0]
      print '\nChoose minimum and maximum of the grid along X:'
      print 'Enter either a single number a (X grid from  xmin-a*width  to  xmax+a*width)'
      print '        or two numbers a and b (X grid from  a  to  b)'
      INFOS['convolute_X']['xrange']=question('Xrange:',float,[kernels[kern]['factor']])
      if len(INFOS['convolute_X']['xrange'])>2:
        INFOS['convolute_X']['xrange']=INFOS['convolute_X']['xrange'][:2]

  # Question 4
  if INFOS['synchronizing'] and not INFOS['convolute_X']:
    print '\n'+centerstring('4 Averaging',40,'-')+'\n'
    if question('Do you want to average the data columns across all trajectories?',bool,False):
      print 'Choose one of the following options:'
      print '%i  %s' % (1,'Arithmetic average and standard deviation')
      print '%i  %s' % (2,'Geometric average and standard deviation')
      while True:
        av=question('Choose one of the options:',int,[1])[0]
        if av in [1,2]:
          break
        else:
          print 'Choose one of the following: %s' % ([1,2])
      if av==1:
        INFOS['averaging']={'mean': mean_arith, 'stdev': stdev_arith}
      elif av==2:
        INFOS['averaging']={'mean': mean_geom, 'stdev': stdev_geom}

  # Question 4
  if INFOS['synchronizing'] and not INFOS['convolute_X']:
    print '\n'+centerstring('5 Total statistics',40,'-')+'\n'
    if question('Do you want to compute the total mean and standard deviation over all time steps?',bool,False):
      print 'Choose one of the following options:'
      print '%i  %s' % (1,'Arithmetic average and standard deviation')
      print '%i  %s' % (2,'Geometric average and standard deviation')
      while True:
        av=question('Choose one of the options:',int,[1])[0]
        if av in [1,2]:
          break
        else:
          print 'Choose one of the following: %s' % ([1,2])
      if av==1:
        INFOS['statistics']={'mean': mean_arith, 'stdev': stdev_arith}
      elif av==2:
        INFOS['statistics']={'mean': mean_geom, 'stdev': stdev_geom}

  # Question 6
  if INFOS['synchronizing'] and INFOS['convolute_X']:
    print '\n'+centerstring('6 Sum over all Y',40,'-')+'\n'
    INFOS['sum_Y']=question('Do you want to sum up all Y values?',bool,False)

  # Question 7
  if INFOS['synchronizing'] and INFOS['convolute_X']:
    print '\n'+centerstring('7 Integrate along X',40,'-')+'\n'
    if question('Do you want to integrate in X direction?',bool,False):
      print 'Please specify the lower and upper bounds for the integration:'
      while True:
        INFOS['integrate_X']['xrange']=question('Xmin and Xmax:',float,[0.0,10.0])
        if len(INFOS['integrate_X']['xrange'])>=2:
          INFOS['integrate_X']['xrange']=INFOS['integrate_X']['xrange'][:2]
          break

  # Question 8
  if INFOS['synchronizing'] and INFOS['convolute_X']:
    print '\n'+centerstring('8 Convoluting along T',40,'-')+'\n'
    if question('Do you want to apply convolution in T direction?',bool,False):
      print 'Choose one of the following convolution kernels:'
      for i in sorted(kernels):
        print '%i  %s' % (i,kernels[i]['description'])
      while True:
        kern=question('Choose one of the functions:',int,[1])[0]
        if kern in kernels:
          break
        else:
          print 'Choose one of the following: %s' % (list(kernels))
      w=question('Choose width of the smoothing function (in units of the X columns):',float,[25.0])[0]
      INFOS['convolute_T']={'function': kernels[kern]['f'](w)}
      #print 'Choose the size of the grid along X:'
      INFOS['convolute_T']['npoints']=question('Size of the grid along T:',int,[200])[0]
      print '\nChoose minimum and maximum of the grid along X:'
      print 'Enter either a single number a (X grid from  xmin-a*width  to  xmax+a*width)'
      print '        or two numbers a and b (X grid from  a  to  b)'
      INFOS['convolute_T']['xrange']=question('Xrange:',float,[kernels[kern]['factor']])
      if len(INFOS['convolute_T']['xrange'])>2:
        INFOS['convolute_T']['xrange']=INFOS['convolute_T']['xrange'][:2]

  # Question 9
  if INFOS['synchronizing'] and INFOS['convolute_X'] and not INFOS['convolute_T']:
    print '\n'+centerstring('9 Integrating along T',40,'-')+'\n'
    INFOS['integrate_T']=question('Do you want to integrate in T direction?',bool,False)

  # Question 10
  if INFOS['synchronizing'] and INFOS['convolute_X']:
    print '\n'+centerstring('10 Convert to Type2 dataset',40,'-')+'\n'
    print 'If you performed integration along X, the data might be better formatted as Type2 dataset.'
    recommend=bool(INFOS['integrate_X'])
    INFOS['type3_to_type2']=question('Do you want to output as Type2 dataset?',bool,recommend)





  #pprint.pprint(INFOS)
  #sys.exit(1)




  return INFOS

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================



def do_calc(INFOS):

  outindex=0
  outstring=''

  # TODO: 

  print '\n\n>>>>>>>>>>>>>>>>>>>>>> Started data analysis\n'

  # ---------------------- collect data -------------------------------
  if True:
    print 'Collecting the data ...'
    data1=collect_data(INFOS)
    outindex=1
    filename=make_filename(outindex,INFOS,outstring)
    print '>>>> Writing output to file "%s"...\n' % filename
    writefile(filename,stringType1(data1,INFOS))

  # ---------------------- apply temporal smoothing -------------------------------
  if INFOS['smoothing']:
    print 'Applying temporal smoothing ...'
    data1=smoothing_xy(INFOS,data1)
    outindex=1
    outstring+='_sm'
    filename=make_filename(outindex,INFOS,outstring)
    print '>>>> Writing output to file "%s"...\n' % filename
    writefile(filename,stringType1(data1,INFOS))

  # ---------------------- apply synchronization -------------------------------
  if INFOS['synchronizing']:
    print 'Synchronizing temporal data ...'
    data2 = synchronize(INFOS,data1)
    outindex=2
    outstring+='_sy'
    filename=make_filename(outindex,INFOS,outstring)
    print '>>>> Writing output to file "%s"...\n' % filename
    writefile(filename,stringType2(data2))

  # ---------------------- compute averages --------------------
  if INFOS['averaging']:
    print 'Computing averages ...'
    data2=calc_average(INFOS,data2)
    outindex=2
    outstring+='_av'
    filename=make_filename(outindex,INFOS,outstring)
    print '>>>> Writing output to file "%s"...\n' % filename
    writefile(filename,stringType2(data2))

  # ---------------------- compute averages --------------------
  if INFOS['statistics']:
    print 'Computing total statistics ...'
    data2=calc_statistics(INFOS,data2)
    outindex=2
    outstring+='_st'
    filename=make_filename(outindex,INFOS,outstring)
    print '>>>> Writing output to file "%s"...\n' % filename
    writefile(filename,stringType2(data2))

  # ---------------------- convoluting X --------------------
  if INFOS['convolute_X']:
    print 'Convoluting data (along X column) ...'
    data3=do_x_convolution(INFOS,data2)
    outindex=3
    outstring+='_cX'
    filename=make_filename(outindex,INFOS,outstring)
    print '>>>> Writing output to file "%s"...\n' % filename
    writefile(filename,stringType3(data3))

  # ---------------------- convoluting X --------------------
  if INFOS['sum_Y']:
    print 'Summing all Y values ...'
    data3=do_y_summation(INFOS,data3)
    outindex=3
    outstring+='_sY'
    filename=make_filename(outindex,INFOS,outstring)
    print '>>>> Writing output to file "%s"...\n' % filename
    writefile(filename,stringType3(data3))

  # ---------------------- integrating X --------------------
  if INFOS['convolute_X'] and INFOS['integrate_X']:
    print 'Integrating data (along X column) ...'
    data3=integrate_X(INFOS,data3)
    outindex=3
    outstring+='_iX'
    filename=make_filename(outindex,INFOS,outstring)
    print '>>>> Writing output to file "%s"...\n' % filename
    writefile(filename,stringType3(data3))

  # ---------------------- convoluting T --------------------
  if INFOS['convolute_X'] and INFOS['convolute_T']:
    print 'Convoluting data (along T column) ...'
    data3=do_t_convolution(INFOS,data3)
    outindex=3
    outstring+='_cT'
    filename=make_filename(outindex,INFOS,outstring)
    print '>>>> Writing output to file "%s"...\n' % filename
    writefile(filename,stringType3(data3))

  # ---------------------- integrating T --------------------
  if INFOS['convolute_X'] and INFOS['integrate_T']:
    print 'Integrating data (along T column) ...'
    data3=integrate_T(INFOS,data3)
    outindex=3
    outstring+='_iT'
    filename=make_filename(outindex,INFOS,outstring)
    print '>>>> Writing output to file "%s"...\n' % filename
    writefile(filename,stringType3(data3))

  # ---------------------- integrating T --------------------
  if INFOS['convolute_X'] and INFOS['type3_to_type2']:
    print 'Converting to Type2 dataset ...'
    data2=type3_to_type2(INFOS,data3)
    outindex=2
    outstring+='_cv'
    filename=make_filename(outindex,INFOS,outstring)
    print '>>>> Writing output to file "%s"...\n' % filename
    writefile(filename,stringType2(data2))


  # ---------------------- write -------------------------------
  #filename=make_filename(outindex,INFOS,outstring)

  print
  print 'Finished!'
  #print 'Writing output to file "%s"...\n' % filename
  #if outindex==1:
    #writefile(filename,stringType1(data1,INFOS))
  #elif outindex==2:
    #writefile(filename,stringType2(data2))
  #elif outindex==3:
    #writefile(filename,stringType3(data3))

  return INFOS


# ===============================================

def make_filename(outindex,INFOS,outstring):
  filename='collected_data_%i_' % (INFOS['colT'])
  for i in INFOS['colX']:
    filename+='%i' % i
  filename+='_'
  for i in INFOS['colY']:
    filename+='%i' % i
  filename+=outstring+'.type%i.txt' % (outindex)
  if len(filename)>=255:
    filename=filename[:15]+'...'+filename[-35:]
  return filename

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def collect_data(INFOS):
  data1={}
  maxcol=max( [INFOS['colT']]+[abs(i) for i in INFOS['colX']]+[abs(i) for i in INFOS['colY']] )
  width_bar=50
  for it1,f in enumerate(INFOS['allfiles']):
    done=width_bar*(it1+1)/len(INFOS['allfiles'])
    sys.stdout.write('\r  Progress: ['+'='*done+' '*(width_bar-done)+'] %3i%%' % (done*100/width_bar))
    #print '  ... %s' % f
    data1[f]=[]
    text=readfile(f)
    iline=-1
    for line in text:
      if '#' in line:
        continue
      s=line.split()
      if len(s) < maxcol:
        continue
      iline+=1
      if INFOS['colT']==0:
        t=iline
      else:
        t=float( s[INFOS['colT']-1] )
      x=[]
      for i in INFOS['colX']:
        if i==0:
          x.append(1.)
        elif i<0:
          x.append(-float( s[-i-1] ))
        elif i>0:
          x.append(+float( s[i-1] ))
      y=[]
      for i in INFOS['colY']:
        if i==0:
          y.append(1.)
        elif i<0:
          y.append(-float( s[i-1] ))
        elif i>0:
          y.append(+float( s[i-1] ))
      data1[f].append( tuple([t]+x+y) )
    data1[f].sort(key=lambda x: x[0])
  print
  return data1

# ===========================================
def smoothing_xy( INFOS, data1 ):
  data2={}
  f=INFOS['smoothing']['function']
  width_bar=50
  for it1,key in enumerate(data1):
    done=width_bar*(it1+1)/len(data1)
    sys.stdout.write('\r  Progress: ['+'='*done+' '*(width_bar-done)+'] %3i%%' % (done*100/width_bar))
    sys.stdout.flush()
    #print '  ... %s' % key
    data2[key]=[]
    for T in data1[key]:
      t=T[0]
      out=[t]
      for i in range(1,len(T)):
        n=0.
        s=0.
        for T1 in data1[key]:
          t1=T1[0]
          w=f.ev(1.,t,t1)
          if w>0.:
            n+=w
            s+=w*T1[i]
        out.append(s/n)
      data2[key].append(out)
  print
  return data2

# ===========================================
def synchronize( INFOS, data1 ):
  # get all times
  times=set()
  for traj in data1:
    for T in data1[traj]:
      times.add(T[0])
  times=list(times)
  times.sort()
  # order data and add NaNs
  data2=[ [] for i in times ]
  width_bar=50
  for ik,key in enumerate(sorted(data1)):
    done=width_bar*(ik+1)/len(data1)
    sys.stdout.write('\r  Progress: ['+'='*done+' '*(width_bar-done)+'] %3i%%' % (done*100/width_bar))
    #print '  ... %s' % traj
    iterator=iter(data1[key])
    t=min(times)-1.
    for it1,t1 in enumerate(times):
      if t<t1:
        try:
          T=next(iterator)
          t=T[0]
        except StopIteration:
          t=None
      if t1==t:
        d=tuple(T[1:])
      else:
        d=tuple( [ float('NaN') for i in T[1:] ] )
      data2[it1].append(d)
  print
  # convert to dict
  data3={'times':times,
         'data':data2}
  # find extrema of data
  data3['tmin']=min(times)
  data3['tmax']=max(times)
  nx=len(data2[0][0])/2
  xmin=data2[0][0][0]
  xmax=xmin
  ymin=data2[0][0][nx]
  ymax=ymin
  for T in data2:
    for X in T:
      xmin=min( [xmin]+list(X[:nx]) )
      xmax=max( [xmax]+list(X[:nx]) )
      ymin=min( [ymin]+list(X[nx:]) )
      ymax=max( [ymax]+list(X[nx:]) )
  # remember which columns to print and make labels
  mask=[]
  labels=[]
  for i in INFOS['colX']:
    if i!=0:
      mask.append(True)
      labels.append('X Column %3i' % i)
    else:
      mask.append(False)
  for i in INFOS['colY']:
    if i!=0:
      mask.append(True)
      labels.append('Y Column %3i' % i)
    else:
      mask.append(False)
  toprint=[]
  labels=['Time']+labels*len(data1)
  for traj in data1:
    toprint.append(mask)
  data3['xmin']=xmin
  data3['xmax']=xmax
  data3['ymin']=ymin
  data3['ymax']=ymax
  data3['toprint']=toprint
  data3['labels']=labels
  return data3

# ===========================================
def calc_average(INFOS,data2):
  data3=[]
  times=[]
  ndata=[]
  width_bar=50
  for it1,t1 in enumerate(data2['times']):
    done=width_bar*(it1+1)/len(data2['times'])
    sys.stdout.write('\r  Progress: ['+'='*done+' '*(width_bar-done)+'] %3i%%' % (done*100/width_bar))
    T=data2['data'][it1]
    means=[]
    stdevs=[]
    for ix in range(len(T[0])):
      array=[]
      for traj in T:
        if traj[ix]==traj[ix]:
          array.append(traj[ix])
      xmean=INFOS['averaging']['mean'](array)
      xstdv=INFOS['averaging']['stdev'](array,xmean)
      means.append(xmean)
      stdevs.append(xstdv)
    ndata.append(len(array))
    times.append(t1)
    data3.append([tuple(means+stdevs)])
  print
  # remember which columns to print
  toprint=copy.copy(data2['toprint'][0])
  # make labels
  labels=['Time']
  for i in INFOS['colX']:
    if i!=0:
      labels.append('X Mean %3i' % i)
  for i in INFOS['colY']:
    if i!=0:
      labels.append('Y Mean %3i' % i)
  for i in INFOS['colX']:
    if i!=0:
      labels.append('X Stdev %3i' % i)
  for i in INFOS['colY']:
    if i!=0:
      labels.append('Y Stdev %3i' % i)
  labels.append('Count')
  # convert to type2
  data4={}
  data4['ndata']=ndata
  data4['times']=times
  data4['data']=data3
  data4['tmin']=min(times)
  data4['tmax']=max(times)
  data4['xmin']=data2['xmin']
  data4['xmax']=data2['xmax']
  data4['ymin']=data2['ymin']
  data4['ymax']=data2['ymax']
  data4['toprint']=[toprint+toprint]
  data4['labels']=labels
  return data4

# ===========================================
def calc_statistics(INFOS,data2):
  data3=[]
  times=[]
  ndata=[]
  arrays=[ [] for i in data2['data'][0][0] ]
  width_bar=50
  for it1,t1 in enumerate(data2['times']):
    done=width_bar*(it1+1)/len(data2['times'])
    sys.stdout.write('\r  Progress: ['+'='*done+' '*(width_bar-done)+'] %3i%%' % (done*100/width_bar))
    T=data2['data'][it1]
    means=[]
    stdevs=[]
    for ix in range(len(T[0])):
      #array=[]
      for traj in T:
        if traj[ix]==traj[ix]:
          arrays[ix].append(traj[ix])
      xmean=INFOS['statistics']['mean'](arrays[ix])
      xstdv=INFOS['statistics']['stdev'](arrays[ix],xmean)
      means.append(xmean)
      stdevs.append(xstdv)
    ndata.append(len(arrays[ix]))
    times.append(t1)
    data3.append([tuple(means+stdevs)])
  print
  # remember which columns to print
  toprint=copy.copy(data2['toprint'][0])
  # make labels
  labels=['Time']
  for i in INFOS['colX']:
    if i!=0:
      labels.append('X Mean %3i' % i)
  for i in INFOS['colY']:
    if i!=0:
      labels.append('Y Mean %3i' % i)
  for i in INFOS['colX']:
    if i!=0:
      labels.append('X Stdev %3i' % i)
  for i in INFOS['colY']:
    if i!=0:
      labels.append('Y Stdev %3i' % i)
  labels.append('Count')
  # convert to type2
  data4={}
  data4['ndata']=ndata
  data4['times']=times
  data4['data']=data3
  data4['tmin']=min(times)
  data4['tmax']=max(times)
  data4['xmin']=data2['xmin']
  data4['xmax']=data2['xmax']
  data4['ymin']=data2['ymin']
  data4['ymax']=data2['ymax']
  data4['toprint']=[toprint+toprint]
  data4['labels']=labels
  return data4

# ===========================================
def do_x_convolution(INFOS,data2):
  # set up xrange
  width=INFOS['convolute_X']['function'].fwhm
  xmin=data2['xmin']
  xmax=data2['xmax']
  if not INFOS['convolute_X']['xrange']:
    xmin=xmin-2.*width
    xmax=xmax+2.*width
  elif len(INFOS['convolute_X']['xrange'])==2:
    xmin=INFOS['convolute_X']['xrange'][0]
    xmax=INFOS['convolute_X']['xrange'][1]
  elif len(INFOS['convolute_X']['xrange'])==1:
    xmin=xmin-INFOS['convolute_X']['xrange'][0]*width
    xmax=xmax+INFOS['convolute_X']['xrange'][0]*width
  # do convolution
  data3=[]
  width_bar=50
  for it1,t1 in enumerate(data2['times']):
    done=width_bar*(it1+1)/len(data2['times'])
    sys.stdout.write('\r  Progress: ['+'='*done+' '*(width_bar-done)+'] %3i%%' % (done*100/width_bar))
    ny=len(data2['data'][it1][0])/2
    spec=[ spectrum(INFOS['convolute_X']['npoints']-1,xmin,xmax,1.0,1) for i in range(ny) ]
    for i in range(ny):
      spec[i].f=INFOS['convolute_X']['function']
    for T in data2['data'][it1]:
      for i in range(ny):
        if T[ny+i]==T[ny+i] and T[i]==T[i]:
          spec[i].add(T[ny+i],T[i])
    d=[]
    for ix in range(INFOS['convolute_X']['npoints']):
      d.append( [spec[i].spec[ix] for i in range(ny)] )
    data3.append(d)
  print
  # make labels
  labels=['Time','X_axis']
  for i,x in enumerate(INFOS['colX']):
    labels.append('Conv(%i,%i)' % (x,INFOS['colY'][i]))
  # make type3 dictionary:
  data4={}
  data4['data']=data3
  data4['times']=copy.copy(data2['times'])
  data4['xvalues']=copy.copy(spec[0].en)
  data4['labels']=labels
  data4['tmin']=min(data4['times'])
  data4['tmax']=max(data4['times'])
  data4['xmin']=min(data4['xvalues'])
  data4['xmax']=max(data4['xvalues'])
  return data4

# ===========================================
def do_t_convolution(INFOS,data3):
  # set up trange
  width=INFOS['convolute_T']['function'].fwhm
  tmin=data3['tmin']
  tmax=data3['tmax']
  if not INFOS['convolute_T']['xrange']:
    tmin=tmin-2.*width
    tmax=tmax+2.*width
  elif len(INFOS['convolute_T']['xrange'])==2:
    tmin=INFOS['convolute_T']['xrange'][0]
    tmax=INFOS['convolute_T']['xrange'][1]
  elif len(INFOS['convolute_T']['xrange'])==1:
    tmin=tmin-INFOS['convolute_T']['xrange'][0]*width
    tmax=tmax+INFOS['convolute_T']['xrange'][0]*width
  # do convolution
  allspec=[]
  for ix1,x1 in enumerate(data3['xvalues']):
    ny=len(data3['data'][0][ix1])
    spec=[ spectrum(INFOS['convolute_T']['npoints']-1,tmin,tmax,1.0,1) for i in range(ny) ]
    for i in range(ny):
      spec[i].f=INFOS['convolute_T']['function']
    allspec.append(spec)
  width_bar=50
  for it1,t1 in enumerate(data3['times']):
    done=width_bar*(it1+1)/len(data3['times'])
    sys.stdout.write('\r  Progress: ['+'='*done+' '*(width_bar-done)+'] %3i%%' % (done*100/width_bar))
    for ix1,x1 in enumerate(data3['xvalues']):
      for i in range(ny):
        allspec[ix1][i].add(data3['data'][it1][ix1][i],t1)
  print
  times=copy.copy(allspec[0][0].en)
  data4=[]
  for it1,t1 in enumerate(times):
    d=[]
    for ix1,x1 in enumerate(data3['xvalues']):
      d.append( [allspec[ix1][i].spec[it1] for i in range(ny)] )
    data4.append(d)
  # make type3 dictionary:
  data5={}
  data5['data']=data4
  data5['times']=times
  data5['xvalues']=copy.copy(data3['xvalues'])
  data5['labels']=data3['labels']
  data5['tmin']=min(data5['times'])
  data5['tmax']=max(data5['times'])
  data5['xmin']=min(data5['xvalues'])
  data5['xmax']=max(data5['xvalues'])
  return data5

# ===========================================
def integrate_T(INFOS,data3):
  # do cumulative sum for all x values
  data4=copy.deepcopy(data3)
  width_bar=50
  for it1,t1 in enumerate(data3['times']):
    done=width_bar*(it1+1)/len(data3['times'])
    sys.stdout.write('\r  Progress: ['+'='*done+' '*(width_bar-done)+'] %3i%%' % (done*100/width_bar))
    if it1==0:
      continue
    for ix1,x1 in enumerate(data3['xvalues']):
      ny=len(data4['data'][it1][ix1])
      for i in range(ny):
        data4['data'][it1][ix1][i]+=data4['data'][it1-1][ix1][i]
  print
  return data4

# ===========================================
def integrate_X(INFOS,data3):
  # sum up for all x values below a, between a and b, and above b
  xmin=INFOS['integrate_X']['xrange'][0]
  xmax=INFOS['integrate_X']['xrange'][1]
  xvalues=[-1,0,1]
  width_bar=50
  data=[]
  for it1,t1 in enumerate(data3['times']):
    done=width_bar*(it1+1)/len(data3['times'])
    sys.stdout.write('\r  Progress: ['+'='*done+' '*(width_bar-done)+'] %3i%%' % (done*100/width_bar))
    ny=len(data3['data'][it1][0])
    d=[ [ 0. for i in range(ny) ] for j in range(3) ]
    for ix1,x1 in enumerate(data3['xvalues']):
      if x1<xmin:
        out=0
      elif xmin<=x1<=xmax:
        out=1
      elif xmax<x1:
        out=2
      for i in range(ny):
        d[out][i]+=data3['data'][it1][ix1][i]
    data.append(d)
  print
  # make type3 dictionary:
  data5={}
  data5['data']=data
  data5['times']=copy.copy(data3['times'])
  data5['xvalues']=xvalues
  data5['labels']=data3['labels']
  data5['tmin']=min(data5['times'])
  data5['tmax']=max(data5['times'])
  data5['xmin']=min(data5['xvalues'])
  data5['xmax']=max(data5['xvalues'])
  return data5

# ===========================================
def do_y_summation(INFOS,data3):
  width_bar=50
  data=[]
  for it1,t1 in enumerate(data3['times']):
    done=width_bar*(it1+1)/len(data3['times'])
    sys.stdout.write('\r  Progress: ['+'='*done+' '*(width_bar-done)+'] %3i%%' % (done*100/width_bar))
    ny=len(data3['data'][it1][0])
    d=[ [ 0. ] for j in data3['xvalues'] ]
    for ix1,x1 in enumerate(data3['xvalues']):
      for i in range(ny):
        d[ix1][0]+=data3['data'][it1][ix1][i]
    data.append(d)
  print
  # make type3 dictionary:
  data5={}
  data5['data']=data
  data5['times']=copy.copy(data3['times'])
  data5['xvalues']=copy.copy(data3['xvalues'])
  data5['labels']=['Time','X_axis','Y_sum']
  data5['tmin']=min(data5['times'])
  data5['tmax']=max(data5['times'])
  data5['xmin']=min(data5['xvalues'])
  data5['xmax']=max(data5['xvalues'])
  return data5

# ===========================================
def type3_to_type2(INFOS,data3):
  # copy data
  data4=copy.deepcopy(data3)
  # adjust to type2
  del data4['xvalues']
  data4['tmin']=min(data4['times'])
  data4['tmax']=max(data4['times'])
  data4['xmin']=0
  data4['xmax']=0
  labels=['Time']
  for ix1,x1 in enumerate(data3['xvalues']):
    labels.append('X=%8f' % x1)
  data4['labels']=labels
  return data4


# ======================================================================================================================

def mean_arith(data):
  if len(data)<1:
    return float('NaN')
  s=0.
  for i in data:
    s+=i
  return s/len(data)

# ======================================== #
def stdev_arith(data,mean=None):
  if len(data)<2:
    return float('NaN')
  if mean==None:
    m=mean_arith(data)
  else:
    m=mean
  s=0.
  for i in data:
    s+=(i-m)**2
  s=s/(len(data)-1)
  return math.sqrt(s)

# ======================================== #
def mean_geom(data):
  if len(data)<1:
    return float('NaN')
  s=0.
  for i in data:
    s+=math.log(i)
  s=s/len(data)
  return math.exp(s)

# ======================================== #
def stdev_geom(data,mean=None):
  if len(data)<2:
    return float('NaN')
  if mean==None:
    m=math.log(mean_geom(data))
  else:
    m=math.log(mean)
  s=0.
  for i in data:
    s+=(math.log(i)-m)**2
  s=s/(len(data)-1)
  return math.exp(math.sqrt(s))

# ======================================== #



# ======================================================================================================================

def stringType1(type1,INFOS):
  # make header
  longest=max( [len(key) for key in type1] )
  string=   '#    1 '+' '*(longest-1)+'2'+' '*15+'3'
  for i in range(2*len(INFOS['colX'])):
    string+='             %3i' % (i+4)
  string+='\n#Index '+' '*(longest-8)+'Filename'+' '*11+' Time'
  for i in INFOS['colX']:
    string+='    X Column %3i' % (i)
  for i in INFOS['colY']:
    string+='    Y Column %3i' % (i)
  string+='\n'
  # make data string
  for ikey,key in enumerate(sorted(type1)):
    for d in type1[key]:
      string+='%6i %s ' % (ikey,key)
      for i in d:
        string+='% .8E ' % i
      string+='\n'
    string+='\n'
  return string

# ======================================================================================================================

def stringType2(type2):
  string=''
  # make header
  if 'labels' in type2:
    string+='#'
    for i,label in enumerate(type2['labels']):
      if i==0:
        string+='%14i ' % (i+1)
      else:
        string+='%15i ' % (i+1)
    string+='\n#'
    for i,label in enumerate(type2['labels']):
      if i==0:
        string+='%14s ' % label
      else:
        string+='%15s ' % label
    string+='\n'
  # make data string
  for it,t in enumerate(type2['times']):
    string+='% .8E' % t
    for iT,T in enumerate(type2['data'][it]):
      for ix,x in enumerate(T):
        if not 'toprint' in type2 or type2['toprint'][iT][ix]:
          string+=' % .8E' % (x)
    if 'ndata' in type2:
      string+=' % .8E' % (type2['ndata'][it])
    string+='\n'
  return string

# ======================================================================================================================

def stringType3(type3):
  # make header
  string=''
  if 'labels' in type3:
    string+='#'
    for i,label in enumerate(type3['labels']):
      if i==0:
        string+='%14i ' % (i+1)
      else:
        string+='%15i ' % (i+1)
    string+='\n#'
    for i,label in enumerate(type3['labels']):
      if i==0:
        string+='%14s ' % label
      else:
        string+='%15s ' % label
    string+='\n'
  # make data string
  for it,t in enumerate(type3['times']):
    for ix,x in enumerate(type3['xvalues']):
      string+='% .8E % .8E ' % (t,x)
      for y in type3['data'][it][ix]:
        string+='% .8E ' % y
      string+='\n'
    string+='\n'
  return string

# ======================================================================================================================

def readType1(strings):
  print 'Type1 cannot be read currently!'
  sys.exit(1)
  #data1={}
  #for line in strings:
    #s=line.split()
    #if len(s)<1:
      #continue
    #key=s[1]
    #values=tuple( [ float(i) for i in s[2:] ] )
    #if not key in data1:
      #data1[key]=[]
    #data1[key].append(values)
  #for key in data:
    #data1[key].sort(key=lambda x: x[0])
  #return data1

# ======================================================================================================================

def readType3(strings):
  print 'Type2 cannot be read currently!'
  sys.exit(1)

# ======================================================================================================================

def readType3(strings):
  print 'Type3 cannot be read currently!'
  sys.exit(1)

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def main():
  '''Main routine'''

  usage='''
python data_collector.py

This interactive program reads table information from SHARC trajectories.
'''

  description=''
  displaywelcome()
  open_keystrokes()

  INFOS=get_general()

  print '\n\n'+centerstring('Full input',60,'#')+'\n'
  for item in INFOS:
    print item, ' '*(25-len(item)), INFOS[item]
  print ''
  calc=question('Do you want to do the specified analysis?',bool,True)
  print ''

  if calc:
    INFOS=do_calc(INFOS)

  close_keystrokes()


# ======================================================================================================================

if __name__ == '__main__':
  try:
    main()
  except KeyboardInterrupt:
    print '\nCtrl+C makes me a sad SHARC ;-(\n'
    quit(0)
