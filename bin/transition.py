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

# Script for counting hopping events
# 
# usage: python transition.py

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
PI = math.pi

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

version='2.0'
versiondate=datetime.date(2018,2,1)

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
  print 'Script for hop counting started...\n'
  string='\n'
  string+='  '+'='*80+'\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Counting hopping events from SHARC dynamics',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Author: Sebastian Mai',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Version:'+version,80)+'||\n'
  string+='||'+centerstring(versiondate.strftime("%d.%m.%y"),80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  string+='''
This script reads output.lis files files and counts all hopping events
to produce a matrix with the transition counts.
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
  shutil.move('KEYSTROKES.tmp','KEYSTROKES.transition')

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

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def get_general():
  ''''''

  INFOS={}

  # Path to trajectories
  print centerstring('Paths to trajectories',60,'-')
  print '\nPlease enter the paths to all directories containing the "TRAJ_0XXXX" directories.\nE.g. S_2 and S_3. \nPlease enter one path at a time, and type "end" to finish the list.'
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
  if os.path.isfile(inputfilename):
    inputfile=open(inputfilename)
    for line in inputfile:
      if 'nstates' in line.lower():
        guessstates=[]
        l=line.split()
        for i in range(1,len(l)):
          guessstates.append(int(l[i]))


  # Analyze mode
  print centerstring('Analyze Mode',60,'-')
  print '''\nThis script finds the transition matrix:
1        In MCH basis                                                    from output.lis
2        In MCH basis (ignoring hops within one multiplet)               from output.lis

This script can also print the transition matrix for each timestep:
3        In MCH basis                                                    from output.lis
4        In MCH basis (ignoring hops within one multiplet)               from output.lis
'''
  while True:
    num=question('Analyze mode:',int)[0]
    if not 1<=num<=4:
      print 'Please enter an integer between 1 and 2!'
      continue
    break
  INFOS['mode']=num
  print ''


  # Run data extractor
  if INFOS['mode'] in []:
    print 'Run data_extractor.x for each trajectory prior to performing the analysis?\nFor many or long trajectories, this might take some time.'
    run_extractor=yesnoquestion('Run data_extractor.x?')
  else:
    run_extractor=False
  INFOS['run_extractor']=run_extractor


  # Number of states
  if INFOS['mode'] in [1,2,3,4]:
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

  # Simulation time
  if INFOS['mode'] in [1,2,3,4]:
    print centerstring('Simulation time',60,'-')
    print '\nUp to which simulation time should the analysis be performed?'
    while True:
      time=question('Simulation time (in fs): ',float,[1000.])[0]
      if time <0.:
        print 'Time must be positive!'
        continue
      break
    INFOS['maxtime']=time
    print ''

  # Intervals
  if INFOS['mode'] in []:
    print centerstring('Intervals',60,'-')
    print '\nPlease enter the interval limits, all on one line.'
    while True:
      nst=raw_input('Interval limits: ')
      nst=re.sub('#.*$','',nst)
      nst=nst.split()
      if len(nst)==0:
        continue
      limits=[]
      try:
        for i in nst:
          limits.append(float(i))
      except ValueError:
        print 'Please enter a list of floats!'
        continue
      break
    INFOS['histo']=histogram(limits)
  print ''


  # States involved in hopping and direction
  #if INFOS['mode'] in [1]:
    #INFOS['fromstates']=[]
    #print centerstring('States involved in surface hop',60,'-')+'\n'
    #print 'In this analysis mode, all geometries are fetched where a trajectory switches from a given MCH state to another given MCH state.\n\nPlease enter the old MCH state involved (mult state):'
    #while True:
      #rmult,rstate=tuple(question('State 1:',int)[0:2])
      #if rmult>len(INFOS['states']):
        #print '%i is larger than maxmult (%i)!' % (rmult,len(INFOS['states']))
        #continue
      #if rstate>INFOS['states'][rmult-1]:
        #print 'Only %i states of mult %i' % (INFOS['states'][rmult-1],rmult)
        #continue
      #break
    #INFOS['fromstates'].append([rmult,rstate])

    #INFOS['tostates']=[]
    #print '\nPlease enter the new MCH state involved (mult state):'
    #while True:
      #rmult,rstate=tuple(question('State 2:',int)[0:2])
      #if rmult>len(INFOS['states']):
        #print '%i is larger than maxmult (%i)!' % (rmult,len(INFOS['states']))
        #continue
      #if rstate>INFOS['states'][rmult-1]:
        #print 'Only %i states of mult %i' % (INFOS['states'][rmult-1],rmult)
        #continue
      #break
    #INFOS['tostates'].append([rmult,rstate])

    #print '''\nDirection:
#1       Forwards
#2       Backwards
#3       Two-way
#'''
    #while True:
      #num=question('Direction mode:',int,[3])[0]
      #if not 1<=num<=3:
        #print 'Please enter an integer between 1 and 3!'
        #continue
      #break
    #INFOS['dirmode']=num
    #print ''
    #if num==1:
      #pass
    #elif num==2:
      #INFOS['fromstates'],INFOS['tostates']=INFOS['tostates'],INFOS['fromstates']
    #elif num==3:
      #INFOS['fromstates'].extend(INFOS['tostates'])
      #INFOS['tostates']=INFOS['fromstates']

  return INFOS

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def print_transition_matrix(transition,labels):
  string=' '*8+'|'
  for i in range(len(transition)):
    string+='%8s' % (labels[i].rstrip())
  string+='\n'
  string+='-'*8+'+'+'-'*(8*len(transition))+'\n'
  for i in range(len(transition)):
    string+='%8s' % (labels[i])+'|'
    for j in range(len(transition)):
      string+=' %7i' % (transition[i][j])
    string+='\n'
  return string

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
            os.chdir(path)
            io=sp.call(sharcpath+'/data_extractor.x output.dat > /dev/null 2> /dev/null',shell=True)
            if io!=0:
              print 'WARNING: extractor call failed for %s!' % (path)
            os.chdir(cwd)
    print 'Extraction finished!\n'

  width=30
  # prepare the list of output.lis files
  files=[]
  ntraj=0
  print 'Checking the directories...'
  for idir in INFOS['paths']:
    ls=os.listdir(idir)
    ls.sort()
    for itraj in ls:
      if not 'TRAJ_' in itraj:
        continue
      path=idir+'/'+itraj
      s=path+' '*(width-len(path))
      if INFOS['mode'] in [1,2,3,4]:
        pathfile=path+'/output.lis'
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
    print 'No trajectories found, exiting...'
    sys.exit(0)

  # get timestep
  if INFOS['mode'] in [1,2,3,4]:
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
      t0=float(f[1])
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
      dt=(float(f[1])-t0)/N
      if dt==0.:
        print 'ERROR: Timestep is zero.'
        quit(1)
      lisf.close()
      break

    # get number of steps
    nsteps=int(INFOS['maxtime']/dt)+1
    print 'Number of steps: %i' % (nsteps)


  # make empty transition matrices
  if INFOS['mode']==1:
    transition=[ [ 0 for i in range(INFOS['nmstates']) ] for j in range(INFOS['nmstates']) ]
  elif INFOS['mode']==2:
    transition=[ [ 0 for i in range(INFOS['nstates']) ] for j in range(INFOS['nstates']) ]
  elif INFOS['mode']==3:
    transition=[ [ [ 0 for i in range(INFOS['nmstates']) ] for j in range(INFOS['nmstates']) ] for t in range(nsteps) ]
  elif INFOS['mode']==4:
    transition=[ [ [ 0 for i in range(INFOS['nstates']) ] for j in range(INFOS['nstates']) ] for t in range(nsteps) ]

  # make state mapping and labels
  mapping={}
  labels={}
  if INFOS['mode'] in [1,3]:
    for i in range(INFOS['nmstates']):
      mapping[i]=i
      mult,state,ms=tuple(INFOS['statemap'][i+1][0:3])
      label='%1s%i%+3.1f' % (IToMult[mult][0:1],state-(mult<=2),ms)
      labels[i]=label
  elif INFOS['mode'] in [2,4]:
    for i in range(INFOS['nmstates']):
      mult,state,ms,j=tuple(INFOS['statemap'][i+1])
      mapping[i]=j-1
      label='%1s%i    ' % (IToMult[mult][0:1],state-(mult<=2))
      labels[i]=label

  # loop over the permissible trajectories
  string=''
  for ipath in files:
    f=open(ipath)
    lis=f.readlines()
    f.close()
    # go through the lis file line by line, skipping commented lines
    oldstate=-1
    istep=0
    for line in lis:
      if '#' in line:
        continue
      s=line.split()
      state=int(s[3])
      if oldstate==-1:
        oldstate=state
        continue
      if INFOS['mode'] in [1,2]:
        transition[mapping[state-1]][mapping[oldstate-1]]+=1
        istep+=1
        if istep>nsteps:
          break
      elif INFOS['mode'] in [3,4]:
        transition[istep][mapping[state-1]][mapping[oldstate-1]]+=1
        istep+=1
        if istep==len(transition) or istep>nsteps:
          break
      oldstate=state
  #print string

  if INFOS['mode'] in [1,2]:
    print '\n'
    print centerstring('Results',60,'*')

    string=print_transition_matrix(transition,labels)
    print 'Full transition matrix:'
    print string

    sumtrans=copy.deepcopy(transition)
    for i in range(len(transition)):
      for j in range(len(transition)):
        if i<j:
          sumtrans[i][j]=transition[i][j]+transition[j][i]
        elif i==j:
          sumtrans[i][j]=transition[i][j]
        elif i>j:
          sumtrans[i][j]=0
    string=print_transition_matrix(sumtrans,labels)
    print 'Sum transition matrix:'
    print string

    length=len(transition)+1
    difftrans=[ [ 0 for i in range(length) ] for j in range(length)]
    for i in range(len(transition)):
      for j in range(len(transition)):
        difftrans[i][j]=transition[i][j]-transition[j][i]
    for i in range(len(transition)):
      difftrans[i][length-1]=sum( [difftrans[i][j] for j in range(len(transition)) ])
      difftrans[length-1][i]=sum( [difftrans[j][i] for j in range(len(transition)) ])
    difftrans[length-1][length-1]=0
    labels[length-1]='Sum   '
    string=print_transition_matrix(difftrans,labels)
    print 'Difference transition matrix:'
    print string

  elif INFOS['mode'] in [3,4]:
    print '\n'
    print centerstring('Results',60,'*')
    print '\n'

    # make header
    s='#%7i ' % (1)
    k=1
    for i in range(len(transition[0])):
      for j in range(len(transition[0])):
        k+=1
        s+='%5i ' % (k)
    s+='\n'
    s+='#%7s ' % ('Time')
    for i in range(len(transition[0])):
      for j in range(len(transition[0])):
        s+='%2s-%2s ' % (labels[i][:2],labels[j][:2])
    s+='\n'

    for istep,tran in enumerate(transition):
      s+='%8.1f ' % (dt*istep)
      for row in tran:
        for col in row:
          s+='%+5i ' % (col)
      s+='\n'

    # write to file
    outfilename='transition_full.out'
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





# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================


# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def main():
  '''Main routine'''

  usage='''
python transition.py

This interactive program counts hopping events in an ensemble and presents the results as
a transition matrix.
'''

  description=''
  displaywelcome()
  open_keystrokes()

  INFOS=get_general()

  print centerstring('Full input',60,'#')+'\n'
  for item in INFOS:
    if not item=='statemap':
      print item, ' '*(25-len(item)), INFOS[item]
  print ''
  calc=question('Do you want to do the specified analysis?',bool,True)
  print ''

  if calc:
    do_calc(INFOS)

  close_keystrokes()


# ======================================================================================================================

if __name__ == '__main__':
  try:
    main()
  except KeyboardInterrupt:
    print '\nCtrl+C makes me a sad SHARC ;-(\n'
    quit(0)
