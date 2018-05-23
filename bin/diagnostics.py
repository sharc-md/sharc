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

from copy import deepcopy 
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
        s=line.split()[0:-2]
        break
    self.states=[ int(i) for i in s ]
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
  #print 'Script for setup of initial conditions started...\n'
  string='\n'
  string+='  '+'='*80+'\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Diagnostic tool for trajectories from SHARC dynamics',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Author: Sebastian Mai and Moritz Heindl',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Version:'+version,80)+'||\n'
  string+='||'+centerstring(versiondate.strftime("%d.%m.%y"),80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  string+='''
This script reads output.dat files from SHARC trajectories and checks:
* missing files
* normal termination
* total energy conservation
* total population conservation
* discontinuities in potential and kinetic energy
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
  shutil.move('KEYSTROKES.tmp','KEYSTROKES.diagnostics')

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

def print_settings(settings,header='Current settings:'):
  order=[
    'missing_output',
    'missing_restart',
    'normal_termination',
    'always_update',
    'etot_window',
    'etot_step',
    'epot_step',
    'ekin_step',
    'pop_window',
    'hop_energy',
    'intruders'
  ]
  print header
  for i in order:
    print '%22s : %s' % (i,settings[i])
  return




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
  print INFOS
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
        print l
        for i in range(1,len(l)):
          guessstates.append(int(l[i]))
      if 'coupling' in line.lower():
        if 'overlap' in line.lower():
          LD_dynamics=True


  # default diagnostics settings
  defaults={
    'normal_termination':True,
    'missing_output':True,
    'missing_restart':True,
    'etot_window':0.2,
    'etot_step':0.1,
    'epot_step':0.7,
    'ekin_step':0.7,
    'pop_window':1e-7,
    'hop_energy':1.0,
    'intruders':False,
    'always_update':False
  }
  helptext={
    'normal_termination':'Checks for exit status of trajectory (RUNNING, CRASHED, FINISHED).',
    'missing_output':'Checks if "output.lis", "output.log", "output.xyz", "output.dat" are existing.',
    'missing_restart':'Checks if "restart.ctrl", "restart.traj", "restart/" are existing.',
    'etot_window':'Maximum permissible drift in total energy (in eV).',
    'etot_step':'Maximum permissible total energy difference between two successive timesteps (in eV).',
    'epot_step':'Maximum permissible active state potential energy difference between two successive timesteps (in eV). Not checked for timesteps where a hop occurred.',
    'ekin_step':'Maximum permissible kinetic energy difference between two successive timesteps (in eV).',
    'pop_window':'Maximum permissible drift in total population.',
    'hop_energy':'Maximum permissible change in active state energy difference during a surface hop (in eV).',
    'intruders':'Checks if intruder state messages in "output.log" refer to active state.',
    'always_update':'Run data_extractor.x for all trajectories, even if all files have up-to-date time stamps.'
  }
  if LD_dynamics:
    defaults['intruders']=True

  # get settings
  print centerstring('Diagnostic settings',60,'-')
  print '\nPlease, adjust the diagnostic settings according to your preferences.'
  print 'You can use the following commands:\nshow\t\tPrints the current settings\nhelp\t\tPrints explanations for the keys\nend\t\tSave and continue\n<key> <value>\tAdjust setting.\n'
  INFOS['settings']=deepcopy(defaults)
  print_settings(INFOS['settings'])
  print ''
  while True:
    line=question('? ',str,'end',False).lower()
    if line=='end':
      break
    if 'show' in line:
      print_settings(INFOS['settings'])
      continue
    if 'help' in line:
      s=line.split()
      if len(s)==1:
        print_settings(helptext,header='Explanation for keys:')
      elif len(s)==2:
        if s[1] in helptext:
          print '%22s : %s' % (s[1],helptext[s[1]])
      continue
    s=line.split()
    if len(s)!=2:
      print 'Please enter "<key> <setting>".'
      continue
    try:
      f=int(s[1])
      INFOS['settings'][s[0]]=f
      continue
    except ValueError:
      pass
    try:
      f=float(s[1])
      INFOS['settings'][s[0]]=f
      continue
    except ValueError:
      pass
    if 'true' in s[1]:
      INFOS['settings'][s[0]]=True
      continue
    elif 'false' in s[1]:
      INFOS['settings'][s[0]]=False
      continue
    else:
      INFOS['settings'][s[0]]=s[1]



  if not LD_dynamics:
    print 'HINT: Intruder state check only possible if trajectories were propagated with "coupling overlap".'
    defaults['settings']=False


  return INFOS

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
      if x<el:
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

def look_for_files(filelist,path,s,trajectories):
      files=filelist
      s=s
      for ifile in files:
        f=os.path.join(path,ifile)
        s+=ifile[-4:]
        if os.path.isfile(f):
          trajectories[path]['files'][ifile]=True
          s+=' .. '
        else:
          trajectories[path]['files'][ifile]=False
          s+=' !! '
      return s, trajectories

# ======================================================================================================================

def check_files(path,trajectories,INFOS):
  # check if files are there
  trajectories[path]['files']={}
  s='    Output files:     '
  files=['output.lis','output.log','output.dat','output.xyz']
  s, trajectories = look_for_files(files,path,s,trajectories)
  if all(trajectories[path]['files'].values()):
    s+='OK'
    missing = False
    if INFOS['settings']['missing_output']:
      print s
  else:
    s+='Files missing!'
    trajectories[path]['maxsteps']=0
    trajectories[path]['tana']=0.
    missing = True
    print s
    #continue
  #print INFOS['settings']['missing_output']


  # check for restart files
  if INFOS['settings']['missing_restart']:
    files=['restart.ctrl','restart.traj']
    s='    Restart files:    '
    s, trajectories = look_for_files(files,path,s,trajectories)
    ls2=os.path.join(path,'restart')
    if not os.path.isdir(ls2) or len(os.listdir(ls2))==0:
      s+='restart/ !! '
      trajectories[path]['files']['restart']=False
    else:
      s+='restart/ .. '
      trajectories[path]['files']['restart']=True
    if all(trajectories[path]['files'].values()):
      s+='    OK'
    else:
      s+='    Restart might not be possible.'
    if INFOS['settings']['missing_restart']:
      print s
  return trajectories, missing

# ======================================================================================================================

def check_runtime(path, trajectories,INFOS):
  # get maximum run time
  f=os.path.join(path,'output.log')
  f=readfile(f)
  if not check_printlevel(f):
    pass
  trajectories[path]['tana'] = 0
  for line in reversed(f):
    trajectories[path]['laststep']=0
    trajectories[path]['maxsteps']=1
    if 'entering timestep' in line.lower():
      trajectories[path]['laststep']=int(line.split()[3])
      break
  for line in reversed(f):
    if 'found nsteps=' in line.lower():
      trajectories[path]['maxsteps']=int(line.split()[2])
      trajectories[path]['dtstep']=float(line.split()[5])
  s='    Progress:         ['
  progress=float(trajectories[path]['laststep'])/trajectories[path]['maxsteps']
  s+='='*int(25*progress) + ' '*(25-int(25*progress))+']     %.1f of %.1f fs' % (trajectories[path]['laststep']*trajectories[path]['dtstep'], trajectories[path]['maxsteps']*trajectories[path]['dtstep'])
  #s+='\n'
  if INFOS['settings']['normal_termination']:
    print s
  return trajectories, f

# ======================================================================================================================

def check_termination(path, trajectories,INFOS,f):
  # check for normal termination
  trajectories[path]['terminated']=False
  trajectories[path]['crashed']=False
  trajectories[path]['stopped']=False
  trajectories[path]['stuck']=False
  for line in reversed(f[-30:]):
    if 'total wallclock time' in line.lower():
      trajectories[path]['terminated']=True
    elif 'file stop detected' in line.lower():
      trajectories[path]['stopped']=True
    elif 'qm call was not successful' in line.lower():
      trajectories[path]['crashed']=True
  s='    Status:                                           '
  if trajectories[path]['terminated']:
    if trajectories[path]['crashed']:
      s+='CRASHED'
    elif trajectories[path]['stopped']:
      s+='FINISHED (stopped by user)'
    else:
      s+='FINISHED'
  else:
    #check how much time passed since last QM call and compare to average
    #calculation time. Label trajectory STUCK if too much time passed (5x).
    #f = reversed(f)
    timesteps = []
    count = 0
    countmax = min(10,trajectories[path]['laststep'])
    for line in range(len(f)):
      if 'ntering timestep' in f[line].lower():
        check_old=f[line+2].split()
        if check_old[0] == "Start":
          timesteps.append(f[line+2].strip()[12:])
        else:
          timesteps.append(f[line+2].strip())
        count += 1
      if count == countmax: 
        break
    total = datetime.timedelta()
    for entry in range(len(timesteps)-2):
      tstart = datetime.datetime.strptime(timesteps[entry+1], '%a %b %d %H:%M:%S %Y') 
      tend = datetime.datetime.strptime(timesteps[entry+2], '%a %b %d %H:%M:%S %Y')
      total += tend-tstart
    for line in range(len(f)):
      if 'ntering timestep' in f[-line].lower():
        check_old=f[-line+2].split()
        if check_old[0] == "Start":
          tstart = datetime.datetime.strptime(f[-line+2].strip()[12:], '%a %b %d %H:%M:%S %Y')
          break
        else:
          tstart = datetime.datetime.strptime(f[-line+2].strip(), '%a %b %d %H:%M:%S %Y')
          break
    tend=datetime.datetime.now()
    tdiff = tend-tstart
    stuck=False
    if 5*total/(count-2) < tdiff:
      trajectories[path]['stuck']=True
      stuck = True
    if stuck: 
      s+='STUCK (%s since last QM call)' % (str(tdiff)[:-7])
    else:
      s+='RUNNING'
  if INFOS['settings']['normal_termination']:
    print s
  
  return trajectories

# ======================================================================================================================

def check_length(path,trajectories,filelength,filename):
  #checks if the number of entries in a file corresponds to the number 
  #of time steps in the output.log file
  if filelength != trajectories[path]['laststep'] and filelength != trajectories[path]['laststep']+1:
    return 'Wrong step nr in %s ' % (filename)
  else:
    return ''

# ======================================================================================================================

def check_printlevel(logdata):
  for line in logdata:
    if "Print level:" in line:
      level=int(line.split()[-1])
      if level<2:
        return False
      else:
        return True
  return False

# ======================================================================================================================

def check_consistency(path,trajectories,data,filename):
  #checks if no timesteps are omitted in a given file
  problem = ''
  deltatime = float(trajectories[path]['dtstep'])
  if filename == 'output.lis':
    deltatime = 1
  for line in data:
    if not '#' in line:
      x = line.split()
      if float(x[0]) == 0.:
        prevtime = 0
        tana = 0
      elif float(x[0]) - deltatime - prevtime == 0.:
        prevtime = float(x[0])
        tana = prevtime
        pass
      else:
        problem = 'Missing steps in %s' % (filename) #(int(prevtime/deltatime))
        if filename == 'output.lis':
          tana = float(prevtime*float(trajectories[path]['dtstep']))
        else:
          tana = float(prevtime)
        break
    
  return problem, tana

# ======================================================================================================================

def check_energies(path,trajectories,INFOS,hops):
  #look for large changes in the total, kinetic, and potential energy
  # inbetween time steps. Check also if a large change in energy was observed 
  #during a hop.
  f=os.path.join(path,'output_data','energy.out')
  if os.path.isfile(f):
    f=readfile(f)
    try:
      problem, tana_length = check_consistency(path,trajectories,f,'energy.out')
    except:
      print '\n    An error occured while trying to check energy.out for consistency.' 
      trajectories[path]['error'] = True
    if problem == '':
      problem = check_length(path,trajectories,len(f)-3,'energy.out')
    for line in f:
      if '#' in line:
        continue
      x=line.split()
      t=float(x[0])
      e=[ float(i) for i in x[1:] ]
      if t==0.:
        eold=e
        etotmin=e[2]
        etotmax=e[2]
      elif t > tana_length:
        tana = tana_length
        break
      hop = False
      currstep = int(t/trajectories[path]['dtstep']) 
      if currstep in hops:
        hop = True
      ok=True
      tana=t
      if etotmin>e[2]:
        etotmin=e[2]
      if etotmax<e[2]:
        etotmax=e[2]
      if abs(etotmax-etotmin)>INFOS['settings']['etot_window']:
        ok=False
        problem='Large fluctuation in Etot'
      if not hop:
        if abs(e[0]-eold[0]) > INFOS['settings']['ekin_step']:
          ok=False
          problem='Large step in Ekin'
        if abs(e[1]-eold[1]) > INFOS['settings']['epot_step']:
          ok=False
          problem='Large step in Epot'
      else:
        if abs(e[1]-eold[1]) > INFOS['settings']['hop_energy']:
          ok=False
          problem='Large dE during hop'
      if abs(e[2]-eold[2]) > INFOS['settings']['etot_step']:
        ok=False
        problem='Large step in Etot'
      if not ok:
        break
      eold=e
    if len(f) <= 3:
      tana = 0. 
      problem='Empty energy.out file'
    trajectories[path]['tana']=tana
    trajectories[path]['problem']=problem
    s='    Energy:           ' + problem + ' '*(32-len(problem))
    if problem:
      s+='at %.2f fs' % tana
    else:
      s+='OK'
    print s
  else:
    problem='"energy.out" missing'
    s='    Energy:           ' + problem + ' '*(32-len(problem))+'!!'
    trajectories[path]['tana']=0.
    trajectories[path]['problem']=problem
    print s
  return trajectories

# ======================================================================================================================

def check_populations(path,trajectories,INFOS):
  #look for large changes in the total population inbetween time steps
  f=os.path.join(path,'output_data','coeff_diag.out')
  if os.path.isfile(f):
    f=readfile(f)
    tana = 0
    problem=''
    try:
      problem, tana_length = check_consistency(path,trajectories,f,'coeff_diag.out')
    except:
      print '\n    An error occured while trying to check coeff_diag.out for consistency.'
      trajectories[path]['error'] = True
    if problem == '':
      problem = check_length(path,trajectories,len(f)-3,'coeff_diag.out')
    for line in f:
      if '#' in line:
        continue
      x=line.split()
      t=float(x[0])
      pop=float(x[1])
      if t==0.:
        popmin=pop
        popmax=pop
      elif t > tana_length:
        tana = tana_length
        break
      ok=True
      tana=t
      if popmin>pop:
        popmin=pop
      if popmax<pop:
        popmax=pop
      if abs(popmax-popmin)>INFOS['settings']['pop_window']:
        ok=False
        problem='Fluctuation in Population'
      if not ok:
        break
    if len(f) <= 3:
      tana = 0. 
      problem='Empty coeff_diag.out file'
    trajectories[path]['tana']=min(tana,trajectories[path]['tana'])
    trajectories[path]['problem']=problem
    s='    Population:       ' + problem + ' '*(32-len(problem))
    if problem:
      s+='at %.2f fs' % tana
    else:
      s+='OK'
    print s
  else:
    problem='"coeff_diag.out" missing'
    s='    Population:       ' + problem + ' '*(32-len(problem))+'!!'
    trajectories[path]['tana']=0.
    trajectories[path]['problem']=problem
    print s
  return trajectories

# ======================================================================================================================

def check_intruders(path,trajectories,INFOS,lis,tana,problem_length):
  # control for intruder states by comparing detected intruder states in the
  #output.log  to active states in output.lis
  if INFOS['settings']['intruders']:
    f=os.path.join(path,'output.log')
    f=readfile(f)
    ok=True
    problem=''
    notpossible=False
    if not check_printlevel(f):
      notpossible=True
    prevstep=0
    for line in f:
      if 'ntering timestep' in line:
        tstep=int(line.split()[3])
        if tstep == 0:
          prevstep = 0
        elif tstep - 1 != prevstep:
          tana = prevstep * trajectories[path]['dtstep']
          problem = 'Missing steps in output.log'
          break
        prevstep = tstep
        if tstep > tana:
          ok = False
          problem = problem_length
          break
      if 'RESTART requested.' in line:
        prevstep-=1
      if 'State: ' in line:
        intruder=int(line.split()[1])
        if not notpossible:
          state = lis[tstep][2]
          if state==intruder:
            problem='Intruder state found'
            ok=False
            tana=tstep*trajectories[path]['dtstep']
        else:
          problem='Intruder state found (cannot determine time step)'
          ok=False
          tana=0.
        if not ok:
          break
    else:
      tana = trajectories[path]['laststep']*trajectories[path]['dtstep']
    trajectories[path]['tana']=min(tana,trajectories[path]['tana'])
    s='    Intruder states:  ' + problem + ' '*(32-len(problem))
    trajectories[path]['problem']=problem
    if problem:
      s+='at %.2f fs' % tana
    else:
      s+='OK'
    print s
  return trajectories

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def do_calc(INFOS):

  sharcpath=os.getenv('SHARC')
  if sharcpath==None:
    print 'Please set $SHARC to the directory containing the SHARC executables!'
    sys.exit(1)
  cwd=os.getcwd()

  # go through directories
  trajectories={}
  ntraj=0
  print 'Checking the directories...'
  for idir in INFOS['paths']:
    ls=os.listdir(idir)
    ls.sort()
    for itraj in ls:
      if not 'TRAJ_' in itraj:
        continue
      path=os.path.join(idir,itraj)
      trajectories[path]={}
      print centerstring(' '+path+' ',80,'~')+'\n'
      trajectories[path]['error'] = False
      trajectories[path]['filelength'] = ''

      try:
        trajectories, missing = check_files(path,trajectories,INFOS)
      except:
        print '\n    An error occured while trying to look for the files.\n'
        trajectories[path]['error'] = True
        continue
      if missing:
        continue

      try:
        trajectories, f = check_runtime(path, trajectories,INFOS)
      except:
        print '\n    An error occured while trying to extract the runtime.\n \
   Files may be corrupted.\n'
        trajectories[path]['error'] = True
        continue

      try:
        trajectories = check_termination(path, trajectories,INFOS,f)
      except:
        print '\n    An error occured while trying to extract the status.\n \
   Files may be corrupted.\n'
        trajectories[path]['error'] = True   
        continue
 
      # run data extractor
      update=False
      if not os.path.isfile(os.path.join(path,'output_data','expec.out')):
        update=True
      if not update:
        time_dat=os.path.getmtime(os.path.join(path,'output.dat'))
        try:
          time_expec=os.path.getmtime(os.path.join(path,'output_data','expec.out'))
          if time_dat > time_expec:
            update=True
          time_expec=os.path.getmtime(os.path.join(path,'output_data','energy.out'))
          if time_dat > time_expec:
            update=True
          time_expec=os.path.getmtime(os.path.join(path,'output_data','coeff_diag.out'))
          if time_dat > time_expec:
            update=True
        except OSError:
          update=True
      if INFOS['settings']['always_update']:
        update=True

      # run extractor
      if update:
        sys.stdout.write('    Data extractor...                                 ')
        sys.stdout.flush()
        os.chdir(path)
        io=sp.call(sharcpath+'/data_extractor.x -a output.dat > /dev/null 2> /dev/null',shell=True)
        if io!=0:
          print 'WARNING: extractor call failed for %s! Exit code %i' % (path,io)
        os.chdir(cwd)
        sys.stdout.write('OK\n')
      else:
        pass

      #save output.lis as dict with timestep as keys
      lis = {}
      f=os.path.join(path,'output.lis')
      f=readfile(f)
      hops = []
      step = 0
      for line in f:
        if '#' in line:
          if 'Surface Hop' in line:
            hops.append(step)
          continue
        x=line.split()
        lis[float(x[0])]=x[1:]
        step += 1
      try:
        problem, tana = check_consistency(path,trajectories,f,'output.lis')
      except:
        print '\n    An error occured while trying to check output.lis for consistency.'
        trajectories[path]['error'] = True 
      if problem == '':
        problem = check_length(path,trajectories,len(lis),'output.lis')
      try:
        trajectories = check_energies(path,trajectories,INFOS,hops)
      except:
        print '\n    An error occured while trying to extract the energies.\n \
   Files may be corrupted.\n'
        trajectories[path]['error'] = True   
      try:
        trajectories = check_populations(path,trajectories,INFOS)
      except:
        print '\n    An error occured while trying to extract the populations.\n \
   Files may be corrupted.\n'
        trajectories[path]['error'] = True   
      try:
        trajectories = check_intruders(path,trajectories,INFOS,lis,tana,problem)
      except:
        print '\n    An error occured while trying to extract possible intruder states.\n \
   Files may be corrupted.\n'
        trajectories[path]['error'] = True   


      #sys.stdout.write(s)

      if  trajectories[path]['filelength'] != '':
        print trajectories[path]['filelength']
      print '\n\n'



  # statistics
  #pprint.pprint(trajectories)
  #print a summarizing table 
  trajsorted=sorted(trajectories,key=lambda x: trajectories[x]['tana'])
  print '\n\n'+centerstring(' Summary ',80,'=')+'\n'

  print '%30s %6s %6s %6s %6s' % ('Trajectory','Files?','Status','Length','T_use')
  print '%30s %6s %6s %6s %6s\n' % ('','','','(fs)','(fs)')

  maxtime=0.
  for i in trajectories:
    if 'laststep' in trajectories[i] and 'dtstep' in trajectories[i]:
      maxtime=max(maxtime,trajectories[i]['laststep']*trajectories[i]['dtstep'])
  maxtime=100.*int(maxtime/100.+0.999)
  nhisto=5
  hist=histogram( [ maxtime/nhisto*i for i in range(1,1+nhisto) ] )
  hist_data=[0]*(nhisto+1)
  #print hist

  for itraj in trajsorted:
    if all( [trajectories[itraj]['files'][i] for i in ['output.log','output.dat','output.lis','output.xyz'] ] ) and not trajectories[itraj]['error']:
      complete='OK'
    else:
      complete='!!'
    if trajectories[itraj]['filelength'] != '':
      complete='**'
    if complete=='OK':
      if trajectories[itraj]['crashed']:
        status='CRASH'
      elif trajectories[itraj]['stopped']:
        status='STOP'
      elif trajectories[itraj]['stuck']:
        status='STUCK'
      elif trajectories[itraj]['terminated']:
        status='FINISH'
      else:
        status='RUN'
    else:
      status=''
    if complete=='OK':
      full=trajectories[itraj]['maxsteps']*trajectories[itraj]['dtstep']
      length=trajectories[itraj]['laststep']*trajectories[itraj]['dtstep']

      t_use=trajectories[itraj]['tana']
    else:
      full=1000.
      length=0.
      t_use=0.
    diagram='['
    diagram+='='*(int(25.*t_use/full))
    diagram+='-'*(int(25.*length/full)-int(25.*t_use/full))
    diagram+=' '*(25-int(25.*length/full))
    diagram+=']'
    print '%30s %6s %6s %6s %6s   %s' % (itraj,complete,status,length,t_use,diagram)
    hist_data[hist.put(t_use)]+=1

  s='\nThis many trajectories can be used for an analysis up to the given time:\n'
  total=sum(hist_data)
  for i in range(nhisto):
    s+='up to % 5.1f fs:    % 3i  trajectories\n' % (hist.binlist[i],total-sum(hist_data[:i+1]))
  print s

  # get a guess for the T_use threshold
  t_uses=[ trajectories[itraj]['tana'] for itraj in trajsorted ]
  topt=0.
  nopt=0
  for i,t in enumerate(t_uses):
    if (total-i)*t > topt*nopt:
      topt=t
      nopt=total-i

  print centerstring(' Trajectory Flagging ',60,'-')
  print '\nYou can now flag the trajectories according to their maximum usable time.\nIn this way, you can restrict the analysis tools to the set of trajectories with sufficient simulation time.\n'
  flag=question('Do you want to flag the trajectories?',bool,True)
  if flag:
    t_thres=question('Threshold for T_use (fs):',float,[topt])[0]
    nanalyze=0

    for itraj in trajectories:
      tana=trajectories[itraj]['tana']
      f=os.path.join(itraj,'DONT_ANALYZE')
      if tana>=t_thres:
        nanalyze+=1
        if os.path.isfile(f):
          os.remove(f)
      else:
        if not os.path.isfile(f):
          open(f,'w').close()

    print
    print 'Flagged  %i trajectories for analysis.' % (nanalyze)
    print 'Excluded %i trajectories from analysis.' % (total-nanalyze)


  return INFOS

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def main():
  '''Main routine'''

  usage='''
python diagnostics.py

This interactive program reads trajectory files and checks their validity.
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

  close_keystrokes()


# ======================================================================================================================

if __name__ == '__main__':
  try:
    main()
  except KeyboardInterrupt:
    print '\nCtrl+C makes me a sad SHARC ;-(\n'
    quit(0)
