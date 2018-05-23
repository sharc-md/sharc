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

# Interactive script for the calculation of time constant errors
# 
# usage: python bootstrap.py

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

# parallel calculations
from multiprocessing import Pool
#from multiprocessing.pool import ThreadPool as Pool
#from multiprocessing.dummy import Pool

# write debug traces when in pool threads
import traceback
import signal

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

def centerstring(string,n,pad=' '):
  l=len(string)
  if l>=n:
    return string
  else:
    return  pad*((n-l+1)/2)+string+pad*((n-l)/2)

def displaywelcome():
  print 'Script for bootstrap analysis of population fits started...\n'
  string='\n'
  string+='  '+'='*80+'\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Bootstrap time constants from SHARC dynamics',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Author: Sebastian Mai',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Version:'+version,80)+'||\n'
  string+='||'+centerstring(versiondate.strftime("%d.%m.%y"),80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  string+='''
This script reads ensemble populations (from populations.py) and a kinetic model (from make_fitscript.py)
and computes statistics on the kinetic model fitting.
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
  shutil.move('KEYSTROKES.tmp','KEYSTROKES.bootstrap')

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
  print centerstring('Paths to bootstrap data',60,'-')
  print '\nPlease enter the path to the directory containing the raw bootstrap data.\nThis data can be generated with populations.py\n'
  while True:
    path=question('Path: ',str,'bootstrap_data/')
    path=os.path.abspath(os.path.expanduser(os.path.expandvars(path)))
    if not os.path.isdir(path):
      print 'Does not exist or is not a directory: %s' % (path)
      continue
    ls=os.listdir(path)
    print ls
    count=0
    for i in ls:
      if 'pop_' in i:
        count+=1
    if count==0:
      print 'Directory does not contain bootstrap data!'
      continue
    print 'Found %i bootstrap data files.' % count
    break
  INFOS['bootstrap_dir']=path
  print

  # detect number, step, length, ncols
  data=readfile(os.path.join(path,ls[0]))
  dt=-1.
  for iline,line in enumerate(data):
    if '#' in line:
      continue
    break
  s=line.split()
  t0=float(s[0])
  ncol=len(s)-1
  s=data[iline+1].split()
  t1=float(s[0])
  dt=t1-t0
  for iline,line in reversed(list(enumerate(data))):
    if '#' in line:
      continue
    if line.strip()=='':
      continue
    break
  tn=float(line.split()[0])
  steps=int(tn/dt)
  INFOS['ntraj']=count
  INFOS['dt']=dt
  INFOS['steps']=steps
  INFOS['ncol']=ncol

  # nboot
  print centerstring('Number of bootstrap cycles',60,'-')
  print '\nPlease enter the number of bootstrapping cycles to be performed.'
  INFOS['nboot']=question('Number of bootstrap cycles: ',int,[10])[0]

  # Random number seed
  print '\nPlease enter a random number generator seed (type "!" to initialize the RNG from the system time).'
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

  # fitting script
  print centerstring('Fitting script',60,'-')
  print '\nPlease provide the path to the desired fitting script.\nThis script can be generated with make_fitscript.py, and should subsequently me adjusted (suitable guesses for fitted constants).'
  while True:
    path=question('Path: ',str,'model_fit.gp')
    path=os.path.abspath(os.path.expanduser(os.path.expandvars(path)))
    if not os.path.isfile(path):
      print 'Does not exist: %s' % (path)
      continue
    break
  INFOS['fitfile']=path

  print '\nPlease provide the command to execute gnuplot.'
  INFOS['gnuplot']=question('Command: ',str,'gnuplot')

  #parallel runs       TODO: currently is very inefficient due to subprocess creation overhead
  print centerstring('Parallel computation',60,'-')
  print '\nbootstrap.py can run the fitting (which is done through Gnuplot) on multiple CPU).\n(however, the overhead is currently large, so speedup is limited)'
  INFOS['ncpu']=question('Number of CPUs to use: ',int,[1])[0]
  #INFOS['ncpu']=1
  # how many runs to do per cycle
  npercycle=int(math.sqrt(INFOS['nboot']))
  npercycle-=npercycle%INFOS['ncpu']
  #if INFOS['ncpu']==1:
    #npercycle=1
  INFOS['npercycle']=max(npercycle,INFOS['ncpu'])

  return INFOS

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def do_calc(INFOS):

  # read and adjust the fit gnuplot script
  fitscript=readfile(INFOS['fitfile'])
  #print fitscript
  fit2=''
  output_mode=False
  for line in fitscript:
    if '>>' in line:
      output_mode=True
    elif '<<' in line:
      output_mode=False
    elif '<><> 1' in line:
      fit2+='set term png\nset out "model_fit.png"\n'
    elif output_mode:
      fit2+=line
  fit2=fit2.replace('model_fit.dat','data')
  fit2=fit2.replace('model_fit','fit')
  #print fit2

  # read in the ensemble data
  print '>>>>>>>>>>>>> Reading the ensemble data ...'
  ls=os.listdir(INFOS['bootstrap_dir'])
  pop_full=[ [ [0. for j in range(INFOS['ncol']) ] for i in range(INFOS['steps']) ] for k in range(INFOS['ntraj']) ]
  for ifile,filename in enumerate(sorted(ls)):
    if not 'pop_' in filename:
      continue
    filename2=os.path.join(INFOS['bootstrap_dir'],filename)
    data=readfile(filename2)
    print ifile,filename2
    istep=-1
    for line in data:
      if '#' in line:
        continue
      istep+=1
      if istep>=INFOS['steps']:
        break
      s=line.split()
      if len(s)==0:
        continue
      for istate in range(len(s)-1):
        #print ifile,istep,istate
        pop_full[ifile][istep][istate]=float(s[istate+1])

  # make temporary directories
  tmproot=os.path.join(INFOS['bootstrap_dir'],'tmp')
  mkdir(tmproot)



  # enter the main loop: generating files for gnuplot, running gnuplot, reading out...
  print '\n>>>>>>>>>>>>> Starting the bootstrapping cycles ...'
  if INFOS['ncpu']>1:
    print '             (do not use Ctrl-C in parallel mode)\n'
  else:
    print '             (use Ctrl-C to skip the remaining cycles and go to the final analysis)\n'
  outfile=open('bootstrap_cycles.out','w')
  string='Step | %14s | %6s: %8s     %8s | ...\n' % ('Run Time','Key','g.Mean','g.Stdv')
  string+='='*59
  print string
  outfile.write(string+'\n')
  starttime=datetime.datetime.now()
  begintime=datetime.datetime.now()
  idone=0
  constants_all=[]
  prevdir=os.getcwd()
  while True:
    try:
      idone_step=0
      indices=[ [ 0 for i in range(INFOS['ntraj']) ] for j in range(INFOS['npercycle']) ]
      for icpu in range(INFOS['npercycle']):
        idone_step+=1
        idone+=1

        # sample data
        for itraj in range(INFOS['ntraj']):
          r=random.randint(0,INFOS['ntraj']-1)
          indices[icpu][itraj]=r
        #print idone,indices[icpu]

        if idone>=INFOS['nboot']:
          break

      if INFOS['ncpu']>1:
        for icpu in range(INFOS['npercycle']):
          tmpdir=os.path.join(tmproot,'cpu_%i' % icpu)
          mkdir(tmpdir)
      else:
        tmpdir=os.path.join(tmproot,'cpu_%i' % 0)
        mkdir(tmpdir)

      constants_step=[]
      if INFOS['ncpu']>1:
        pool = Pool(processes=INFOS['ncpu'])
        try:
          for icpu in range(idone_step):
            directory=os.path.join(tmproot,'cpu_%i' % icpu)
            constants=pool.apply_async(make_job , [pop_full,indices[icpu],INFOS['dt'],directory,fit2,INFOS['gnuplot']])
            constants_step.append(constants)
          pool.close()
          pool.join()
        except Exception, e:
          pool.close()
          pool.join()
          os.chdir(prevdir)
          raise KeyboardInterrupt
        for i in range(len(constants_step)):
          constants_step[i]=constants_step[i].get()
      else:
        try:
          for icpu in range(idone_step):
            directory=os.path.join(tmproot,'cpu_%i' % 0)
            constants=make_job(pop_full,indices[icpu],INFOS['dt'],directory,fit2,INFOS['gnuplot'])
            constants_step.append(constants)
        except Exception, e:
          os.chdir(prevdir)
          raise KeyboardInterrupt

      for i in constants_step:
        constants_all.append(i)

      if len(constants_all)>=3:
        string=print_intermediate_statistics(constants_all)
        s='%4i | %s %s' % (idone,datetime.datetime.now()-starttime,string)
        print s
        outfile.write(s+'\n')
        outfile.flush()
        starttime=datetime.datetime.now()

      if idone>=INFOS['nboot']:
        break
    except KeyboardInterrupt:
      print 'Aborted, going to final analysis...'
      time.sleep(0.5)
      break
    except ValueError:
      print 'Value Error (e.g., negative time constant), ignoring results of this cycle...'
      continue
  print 'Run time =',datetime.datetime.now()-begintime
  outfile.close()


  allconsts=set()
  for i in constants_all:
    for key in i:
      allconsts.add(key)


  # final analysis
  print '\n>>>>>>>>>>>>> Finished the bootstrapping cycles ...'
  string_all=''

  for key in allconsts:
    string='\n'+centerstring(' Analysis for time constant "%s" ' % key,110,'-')+'\n'

    # calculate
    data=[]
    for i in constants_all:
      if key in i:
        data.append(i[key])
    mini=min(data)
    maxi=max(data)
    mean_a=mean_arith(data)
    stdev_a=stdev_arith(data,mean_a)
    mean_g=mean_geom(data)
    stdev_g=stdev_geom(data,mean_g)

    string+='''
  Arithmetic analysis:          %12.6f +/- %12.6f
                                           ( +/-   %8.2f %%)

  Geometric analysis:           %12.6f  +  %12.6f  -  %12.6f
                                           (  +    %8.2f %%  -    %8.2f %%)

  Minimum and maximum:          %12.6f      and       %12.6f

  Histogram:
  ==========''' % (mean_a,
       stdev_a,
       stdev_a/mean_a*100.,
       mean_g,
       mean_g*(stdev_g-1.),
       -mean_g*(1./stdev_g-1.),
       (stdev_g-1.)*100.,
       -(1./stdev_g-1.)*100.,
       mini,
       maxi
       )
    print string
    string_all+=string+'\n'

    # make histogram
    nhisto=11
    whisto=3
    bins=[ mean_g*stdev_g**( (float(i)-nhisto/2)/(nhisto/2)*whisto ) for i in range(nhisto-1) ]
    h=histogram(bins)
    hout=[ 0 for i in range(nhisto) ]
    for x in data:
      i=h.put(x)
      hout[i]+=1

    # draw histogram
    height=12
    dx=float(height)/max(hout)
    string=''
    for ih in range(height):
      for ib in hout:
        x=ib*dx
        if x>=(height-ih):
          s='#'
        else:
          s=' '
        string+=' %5s' % s
      string+='| %i\n' % ( int((height-ih)/dx))
    string+='     '
    for i in bins:
      string+='   |  '
    string+='\n   '
    for i in bins:
      if i<10.:
        string+=' %5.3f' % (i)
      elif i<100.:
        string+=' %5.2f' % (i)
      elif i<1000.:
        string+=' %5.1f' % (i)
      else:
        string+=' %5i' % (i)
    print string
    string_all+=string+'\n'

  #write results to file
  string='\n\nFull data:\n\n'
  keylist=[]
  for key in constants_all[0]:
    keylist.append(key)
  string+='%10s ' % 'Sample'
  for key in keylist:
    string+='%12s ' % key
  string+='\n'
  for i,c in enumerate(constants_all):
    string+='%10i ' % (i+1)
    for key in keylist:
      if key in c:
        string+='%12.6f ' % c[key]
      else:
        string+='%12s ' % 'NaN'
    string+='\n'
  string_all+=string

  print '\nOutput (analysis and full fitted data) written to "bootstrap.out".'
  writefile('bootstrap.out',string_all)




  return INFOS

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

class KeyboardInterruptError(Exception): pass

def make_job(pop_full,indices,dt,directory,fit2,gnuplot):
  sys.tracebacklimit=0
  #signal.signal(signal.SIGINT, signal.SIG_IGN)
  try:
    # write data
    string=make_data_string(pop_full,indices,dt)
    filename=os.path.join(directory,'data')
    writefile(filename,string)
    # write fit file
    filename=os.path.join(directory,'fit.gp')
    writefile(filename,fit2)
    # run gnuplot
    constants=run_gnuplot(directory,'fit.gp',gnuplot)
  except KeyboardInterrupt:
    raise KeyboardInterruptError()
  except Exception, problem:
    print '*'*50+'\nException in make_job(%s)!' % (directory)
    traceback.print_exc()
    print '*'*50+'\n'
    raise problem

  return constants

# ======================================================================= #
def make_data_string(pop_full,indices,dt):
  ntraj=len(pop_full)
  steps=len(pop_full[0])
  ncol=len(pop_full[0][0])
  pop=[ [ 0. for i in range(ncol) ] for j in range(steps) ]
  for istep in range(steps):
    for icol in range(ncol):
      for x in indices:
        pop[istep][icol]+=pop_full[x][istep][icol]
  string=''
  for istep in range(ncol*steps):
    string+='%12.6f ' % (istep*dt)
    for icol in range(ncol):
      string+='%12.6f ' % (pop[istep%steps][icol]/ntraj)
    string+='\n'
  return string

# ======================================================================= #
def run_gnuplot(directory,filename,gnuplot):
  prevdir=os.getcwd()
  os.chdir(directory)
  string='%s %s' % (gnuplot,filename)
  ps = sp.Popen(string, shell=True, stdout=sp.PIPE,stderr=sp.STDOUT)
  output = ps.communicate()[0]
  if ps.returncode!=0:
    print 'ERROR: gnuplot call not successful!'
  output=output.split('\n')
  constants={}
  for line in output:
    if '&&&' in line:
      #print line
      s=line.split()
      constants[s[1]]=float(s[-1])
      if constants[s[1]]<=0.:
        del constants[s[1]]
  os.chdir(prevdir)
  return constants

# ======================================================================= #
def print_intermediate_statistics(constants_all):
  string=''
  allconsts=set()
  for i in constants_all:
    for key in i:
      allconsts.add(key)
  for key in allconsts:
    data=[]
    for i in constants_all:
      if key in i:
        data.append(i[key])
    #mean_a=mean_arith(data)
    #stdev_a=stdev_arith(data,mean_a)
    mean_g=mean_geom(data)
    stdev_g=stdev_geom(data,mean_g)
    string+='| %6s: %8.2f +/- %8.2f ' % (key,mean_g,mean_g*(stdev_g-1.))
    #string='%s%6s : %8.2f +/- %8.2f    |    %8.2f + %8.2f - %8.2f' % (' '*15,
                                                                    #key,
                                                                    #mean_a,
                                                                    #stdev_a,
                                                                    #mean_g,
                                                                    #mean_g*(stdev_g-1.),
                                                                    #-mean_g*(1./stdev_g-1.)
                                                                    #)
    #print string
  return string




# ======================================== #
def mean_arith(data):
  s=0.
  for i in data:
    s+=i
  return s/len(data)

# ======================================== #
def stdev_arith(data,mean=-9999):
  if mean==-9999:
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
  s=0.
  for i in data:
    s+=math.log(i)
  s=s/len(data)
  return math.exp(s)

# ======================================== #
def stdev_geom(data,mean=-9999):
  if mean==-9999:
    m=math.log(mean_geom(data))
  else:
    m=math.log(mean)
  s=0.
  for i in data:
    s+=(math.log(i)-m)**2
  s=s/(len(data)-1)
  return math.exp(math.sqrt(s))




# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def main():
  '''Main routine'''

  usage='''
python bootstrap.py

This interactive program combines a bootstrap_data directory and a model-fit gnuplot script
and computes error statistics for the model fit.
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

  close_keystrokes()
  if calc:
    INFOS=do_calc(INFOS)



# ======================================================================================================================

if __name__ == '__main__':
  try:
    main()
  except KeyboardInterrupt:
    print '\nCtrl+C makes me a sad SHARC ;-(\n'
    quit(0)
