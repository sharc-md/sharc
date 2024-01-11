#!/usr/bin/env python3

#******************************************
#
#    SHARC Program Suite
#
#    Copyright (c) 2023 University of Vienna
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


# Modules:
# Operating system, isfile and related routines, move files, create directories
import sys

import os
import shutil
# External Calls
import subprocess as sp
# Regular expressions
import re
# debug print for dicts and arrays
import pprint
# sqrt and other math
import math
# copy of arrays of arrays
from copy import deepcopy
# others
import time
import datetime
from optparse import OptionParser
import readline
import colorsys
from socket import gethostname

try:
  import numpy
  NUMPY=True
except ImportError:
  NUMPY=False

# =========================================================0
# compatibility stuff

if sys.version_info[0]!=2:
    print('This is a script for Python 2!')
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


# ======================================================================= #

version = '3.0'
versiondate = datetime.date(2023, 4, 1)

changelogstring='''

'''

# ======================================================================= #
# hash table for conversion of multiplicity to the keywords used in MOLCAS
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

# conversion factors
au2a=0.529177211
rcm_to_Eh=4.556335e-6

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
    print('File %s does not exist!' % (filename))
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
      print('Content %s cannot be written to file!' % (content))
    f.close()
  except IOError:
    print('Could not write to file %s!' % (filename))
    sys.exit(13)

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
  string+='||'+centerstring('Generate GNUPLOT fitting scripts for SHARC populations',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Author: Sebastian Mai',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Version:'+version,80)+'||\n'
  string+='||'+centerstring(versiondate.strftime("%d.%m.%y"),80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  string+='''
This script automatizes the generation of GNUPLOT scripts for fitting of SHARC populations 
(as generated with populations.py) to general kinetic models based on first-order population transfer.

'''
  print(string)

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

Output:
1     1         #FF0000
1     2         #FF7F00
1     3         #FFFF00
1     4         #7FFF00
1     5         #00FF00
1     6         #00FF7F
3     1         #00FFFF
3     2         #0000FF
3     2         #FF00FF
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

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def open_keystrokes():
  global KEYSTROKES
  KEYSTROKES=open('KEYSTROKES.tmp','w')

def close_keystrokes():
  KEYSTROKES.close()
  shutil.move('KEYSTROKES.tmp','KEYSTROKES.make_fitscript')

# ===================================

def question(question,typefunc,default=None,autocomplete=True,ranges=False):
  if typefunc==int or typefunc==float:
    if not default==None and not isinstance(default,list):
      print('Default to int or float question must be list!')
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

    line=input(s)
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
        print('I didn''t understand you.')
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
        print('Please enter floats!')
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
          print('Please enter integers or ranges of integers (e.g. "-3~-1  2  5~7")!')
        else:
          print('Please enter integers!')
        continue

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def label_valid(label):
  # allowed format: letter followed by letters and numbers and underscore
  if '__' in label:
    return False
  if label=='F' or label=='x':
    return False
  if re.match('^[a-zA-Z][a-zA-Z0-9_]*$',label)==None:
    return False
  else:
    return True

# ===================================================

def print_reactions(rate_matrix,specmap):
  # get longest label
  n=3
  for i in range(len(rate_matrix)):
    if len(specmap[i])>n:
      n=len(specmap[i])
    for j in range(len(rate_matrix)):
      if len(rate_matrix[i][j])>n:
        n=len(rate_matrix[i][j])
  formatstring=' %%%is ' % (n)
  # construct string
  string=' '*(n+2)+'|'
  for i in range(len(rate_matrix)):
    string+=formatstring % (specmap[i])
  string+='\n'
  string+='-'*(n+2)+'+'+'-'*((n+2)*len(rate_matrix))+'\n'
  for i in range(len(rate_matrix)):
    string+=formatstring % (specmap[i])+'|'
    for j in range(len(rate_matrix)):
      label=rate_matrix[i][j]
      if label=='':
        label='.'
      string+=formatstring % (label)
    string+='\n'
  string+='\n(Initial species: rows; Final species: columns)\n'
  return string

# ===================================================

def nullspace(a, rtol=1e-5):
    u, s, v = numpy.linalg.svd(a)
    rank = (s > rtol*s[0]).sum()
    return rank,v[rank:].T.copy()

# ===================================================

def get_cycles(rate_matrix):
  # get nspec and nreact
  nspec=len(rate_matrix)
  nreact=0
  for i in range(nspec):
    for j in range(nspec):
      if rate_matrix[i][j]!='':
        nreact+=1
  if nreact==0:
    return -1
  # construct stochiometry matrix
  A=[ [ 0 for i in range(nreact) ] for j in range(nspec) ]
  ireact=0
  for i in range(nspec):
    for j in range(nspec):
      if rate_matrix[i][j]!='':
        A[i][ireact]=-1
        A[j][ireact]=+1
        ireact+=1
  #pprint.pprint(A)
  # get nullspace of A
  if NUMPY:
    rank,null=nullspace(A)
    nullrank=len(A[0])-rank
    if nullrank>0:
      print('  The reaction network contains %i %s!' % (nullrank,['cycles','cycle'][nullrank==1]))
    return nullrank
  else:
    print('  Hint: Cannot check for cycles without NUMPY!')
    return -1

# ===================================================

def check_pop_file(content):
  # checks whether the content of a populations file is valid.
  # checks number of columns, time range and ordering, whether all rows have same number of columns
  # also extracts the numerical data
  data=[]
  ncol=-1
  maxtime=-1.
  for line in content:
    line=re.sub('#.*$','',line)
    if line=='\n':
      continue
    s=line.split()
    # check time
    time=float(s[0])
    if maxtime==-1. and time!=0.:
      print('  Time does not start at zero!')
      return False,0,0,[]
    if time<0.:
      print('  Negative times detected!')
      return False,0,0,[]
    if time <maxtime:
      print('  Times not ordered!')
      return False,0,0,[]               # TODO: maybe this check is not necessary
    maxtime=time
    # check data
    d=[ float(i) for i in s ]
    if any( [i<0. for i in d] ):
      print('  Negative populations detected!')
      return False,0,0,[]
    col=len(d)
    if ncol==-1:
      ncol=col
    elif ncol!=col:
      print('  Inconsistent number of columns detected!')
      if ncol>col:
        ncol=col
    data.append(d)
  return True,maxtime,ncol,data

# ===================================================

def get_infos():
  '''This routine asks for definitions of the kinetic model and the populations file.'''

  INFOS={}

  print(centerstring('',60,'#'))
  print(centerstring('Kinetics Model',60,'#'))
  print(centerstring('',60,'#')+'\n\n')
  # =========================== Define the kinetic model species ==================================
  print(centerstring('Model Species',60,'-')+'\n')
  print('''First, please specify the set of species used in your model kinetics.

Possible input:
+ <label> <label> ...   Adds one or several species to the set
- <label> <label> ...   Removes one or several species from the set
show                    Show the currently defined set of species
end                     Finish species input

Each label must be unique. Enter the labels without quotes.
''')
  species=[]
  while True:
    line=question('Input:',str,'end',False)
    s=line.split()
    if 'end' in s[0].lower():
      if len(species)==0:
        print('  No species added yet!')
      else:
        break
    elif 'show' in s[0].lower():
      print('  Current set:  %s\n' % (species))
    elif '+' in s[0]:
      for i in s[1:]:
        if i in species:
          print('  Species \'%s\' already in set!' % (i))
        else:
          if label_valid(i):
            species.append(i)
            print('  Species \'%s\' added!' % (i))
          else:
            print('  Invalid label \'%s\'! Labels must be a letter followed by letters, \n  numbers and single underscores! "F" and "x" are reserved!' % (i))
    elif '-' in s[0]:
      for i in s[1:]:
        if i in species:
          species.remove(i)
          print('  Species \'%s\' removed!' % (i))
        else:
          print('  Species \'%s\' not in set!' % (i))
    else:
      print('  I did not understand you.')
  print('\nFinal species set:  %s\n' % (species))
  nspec=len(species)
  specmap={}
  for i in range(len(species)):
    specmap[species[i]]=i
    specmap[i]=species[i]
  INFOS['nspec']=nspec
  INFOS['specmap']=specmap

  # =========================== Define the kinetic model reactions ==================================
  print(centerstring('Model Elementary Reactions',60,'-')+'\n')
  print('''Second, please specify the set of elementary reactions in your model kinetics.

Possible input:
+ <species1> <species2> <rate_label>       Add a reaction from species1 to species2 with labelled rate constant
- <rate_label>                             Remove the reaction(s) with the given rate constant
show                                       Show the currently defined set of reactions (as directed adjacency matrix)
end                                        Finish reaction input

Each rate label must be unique.
''')
  rate_matrix=[ [ '' for i in range(nspec) ] for j in range(nspec) ]
  rateset=set()
  while True:
    line=question('Input:',str,'end',False)
    s=line.split()
    if 'end' in s[0].lower():
      if len(rateset)==0:
        print('  No reactions added yet!')
      else:
        break
    elif 'show' in s[0].lower():
      print(print_reactions(rate_matrix,specmap))
    elif '+' in s[0]:
      if len(s)<4:
        print('Please write "+ species1 species2 ratelabel"!')
        continue
      if s[1]==s[2]:
        print('  Species labels identical! No reaction added.')
        continue
      if s[1] in specmap and s[2] in specmap and not s[3] in specmap:
        if rate_matrix[specmap[s[1]]][specmap[s[2]]]!='':
          print('Please remove rate constant %s first!' % (rate_matrix[specmap[s[1]]][specmap[s[2]]]))
          continue
        if not s[3] in rateset:
          if label_valid(s[3]):
            rateset.add(s[3])
            rate_matrix[specmap[s[1]]][specmap[s[2]]]=s[3]
            rank=get_cycles(rate_matrix)
            print('  Reaction from \'%s\' to \'%s\' with rate label \'%s\' added!' % (s[1],s[2],s[3]))
          else:
            print('  Invalid label \'%s\'! Labels must be a letter followed by letters, numbers and single underscores!' % (s[3]))
        else:
          print('  Rate label \'%s\' already defined!' % (s[3]))
          anyways=question('Do you want to add it anyways (i.e., use two reactions with same rate constant)?',bool,False)
          if anyways:
            rate_matrix[specmap[s[1]]][specmap[s[2]]]=s[3]
            rank=get_cycles(rate_matrix)
            print('  Reaction from \'%s\' to \'%s\' with rate label \'%s\' added!' % (s[1],s[2],s[3]))
      else:
        if not s[1] in specmap:
          print('  Species \'%s\' not defined!' % (s[1]))
        if not s[2] in specmap:
          print('  Species \'%s\' not defined!' % (s[2]))
        if s[3] in specmap:
          print('  Label \'%s\' already used for a species!' % (s[3]))
    elif '-' in s[0]:
      if s[1] in rateset:
        rateset.remove(s[1])
        for i in range(nspec):
          for j in range(nspec):
            if rate_matrix[i][j]==s[1]:
              rate_matrix[i][j]=''
              rank=get_cycles(rate_matrix)
      else:
        print('  Rate label \'%s\' not defined!' % (s[1]))
    else:
      print('  I did not understand you.')
  print('\nFinal reaction network:')
  print(print_reactions(rate_matrix,specmap))
  INFOS['rateset']=rateset
  INFOS['rate_matrix']=rate_matrix
  INFOS['rank']=rank

  # =========================== Define the kinetic model initial conditions ==================================
  print(centerstring('Model Initial Conditions',60,'-')+'\n')
  print('''Third, please specify species with non-zero initial populations.

Possible input:
+ <species>       Declare species to have non-zero initial population
- <species>       Remove species from the set of non-zero initial populations
show              Show the currently defined non-zero initial populations
end               Finish initial condition input
''')
  initset=set()
  while True:
    line=question('Input:',str,'end',False)
    s=line.split()
    if 'end' in s[0].lower():
      if len(initset)==0:
        print('  No species with non-zero initial population yet!')
      else:
        break
    elif 'show' in s[0].lower():
      print('  Current set:  %s\n' % (list(initset)))
    elif '+' in s[0]:
      for i in s[1:]:
        if i not in specmap:
          print('  Species \'%s\' not defined!' % (i))
        elif i in initset:
          print('  Species \'%s\' already in set!' % (i))
        else:
          initset.add(i)
          print('  Species \'%s\' added!' % (i))
    elif '-' in s[0]:
      for i in s[1:]:
        if i in initset:
          initset.remove(i)
          print('  Species \'%s\' removed!' % (i))
        else:
          print('  Species \'%s\' not in set!' % (i))
    else:
      print('  I did not understand you.')
  print('\nFinal initial species set:  %s\n' % (list(initset)))
  INFOS['initset']=initset

  print(centerstring('',60,'#'))
  print(centerstring('Fitting Data',60,'#'))
  print(centerstring('',60,'#')+'\n\n')

  # =========================== Define the data file ==================================
  print(centerstring('Population data file',60,'-')+'\n')
  print('''Please specify the path to the population data file (as generated by populations.py).
''')
  while True:
    popfile=question('Populations file:',str,'pop.out',True)
    if not os.path.isfile(popfile):
      print('  File not found!')
      continue
    content=readfile(popfile)
    valid,maxtime,ncol,data=check_pop_file(content)
    if not valid:
      print('  File format not valid!')
      continue
    else:
      print('  Detected maximal time of %7.1f fs and %i columns (time plus %i data columns).' % (maxtime,ncol,ncol-1))
      break
  INFOS['maxtime']=maxtime
  INFOS['ncol']=ncol
  INFOS['data']=data
  INFOS['popfile']=os.path.abspath(popfile)

  # =========================== Define the data -- species mapping ==================================
  print('\n'+centerstring('Population-to-Species Mapping for Fit',60,'-')+'\n')
  print('''Please specify which model species should be fitted to which data file columns.
For example, you can fit the label 'S0' to column 2:
  S0 = 2
You can also fit the sum of two species to a column:
  T1 T2 = 5
You can also fit a species to the sum of several columns:
  T_all = 5 6 7
You can even fit a sum of species to a sum of columns:
  T1 T2 = 5 6 7

On the right side, "~" can be used to indicate ranges:
  T1 T2 = 5~9

Possible input:
<species1> <species2> ... = <columns1> <column2> ...            Set one mapping
show                                                            Show mapping
end                                                             Finish mapping input
reset                                                           Redo the mapping input

Each species label must be used at most once.
Each column number (except for \'1\', which denotes the time) must be used at most once.
''')
  print('Set of species:        %s' % (species))
  columns=[i for i in range(2,ncol+1)]
  print('Set of column numbers: %s' % (columns))

  species_groups=[]
  columns_groups=[]
  ngroups=0
  while True:
    line=question('Input:',str,'end',False)
    s=line.split()
    if 'end' in s[0].lower():
      if ngroups==0:
        print('  No valid input yet!')
      else:
        break
    elif 'show' in s[0].lower():
      print('  Current mapping groups:')
      for i in range(ngroups):
        string='    '
        for j in species_groups[i]:
          string+=' %s ' % (j)
        string+=' = '
        for j in columns_groups[i]:
          string+=' %i ' % (j)
        print(string)
    elif ' = ' in line:
      if s[0]=='=' or s[-1]=='=':
        print('  Invalid input! Put species labels to the left of \'=\' and column number to the right!')
        continue
      do_species=True
      valid=True
      new_species_group=[]
      new_columns_group=[]
      for i in s:
        if i=='=':
          if not do_species:
            print('More than 1 "=" used!')
            valid=False
            break
          else:
            do_species=False
            continue
        if do_species:
          if not i in species:
            print('  Species label \'%s\' not defined!' % (i))
            valid=False
            break
          if any( [ i in j for j in species_groups ] ):
            print('  Species label \'%s\' already assigned!' % (i))
            valid=False
            break
          if i in new_species_group:
            print('  Species label \'%s\' used twice!' % (i))
            valid=False
            break
          new_species_group.append(i)
        else:
          try:
            if '~' in i:
              ii=[]
              q=i.split('~')
              for j in range(int(q[0]),int(q[1])+1):
                ii.append(j)
            else:
              ii=[int(i)]
          except ValueError:
            print('  Could not understand!')
            valid=False
            break
          for i in ii:
            if not i in columns:
              print('  Column number %i not in data file!' % (i))
              valid=False
              break
            if any( [ i in j for j in columns_groups ] ):
              print('  Column number %i already assigned!' % (i))
              valid=False
              break
            if i in new_columns_group:
              print('  Columns number %i used twice!' % (i))
              valid=False
              break
            new_columns_group.append(i)
      if valid:
        species_groups.append(new_species_group)
        columns_groups.append(new_columns_group)
        ngroups+=1
    elif 'reset' in s[0].lower():
      species_groups=[]
      columns_groups=[]
      ngroups=0
      print('  Mappings reset! Please repeat input!')
  print('Final mappings:')
  for i in range(ngroups):
    string='    '
    for j in species_groups[i]:
      string+=' %s ' % (j)
    string+=' = '
    for j in columns_groups[i]:
      string+=' %i ' % (j)
    print(string)
  INFOS['species_groups']=species_groups
  INFOS['columns_groups']=columns_groups
  INFOS['ngroups']=ngroups

  return INFOS

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def wrap_line(line,lenline=150):
  allowed_after='()/*-+='
  if line=='' or line=='\n':
    return []
  elif len(line)<lenline:
    return [line.strip()]
  else:
    for i in range(lenline,0,-1):
      if line[i:i+1] in allowed_after:
        return [line[:i+1].strip()]+wrap_line(line[i+1:],lenline)

# ===================================================

def get_functions_from_maxima(INFOS):
  # write MAXIMA input file
  # system of differential equations
  string='display2d:false;\neq:[\n'
  for i in range(INFOS['nspec']):
    string+='\'diff(%s(t),t)= ' % (INFOS['specmap'][i])
    nterms=0
    for j in range(INFOS['nspec']):
      if INFOS['rate_matrix'][i][j]!='':
        string+=' -%s*%s(t) ' % (INFOS['rate_matrix'][i][j],INFOS['specmap'][i])
        nterms+=1
      if INFOS['rate_matrix'][j][i]!='':
        string+=' +%s*%s(t) ' % (INFOS['rate_matrix'][j][i],INFOS['specmap'][j])
        nterms+=1
    if nterms==0:
      string+='0'
    if i+1==INFOS['nspec']:
      string+='\n];\n'
    else:
      string+=',\n'
  string+='\n'
  # initial values
  for i in range(INFOS['nspec']):
    if INFOS['specmap'][i] in INFOS['initset']:
      initvalue=INFOS['specmap'][i]+'__0'
    else:
      initvalue='0'
    string+='atvalue(%s(t),t=0,%s);\n' % (INFOS['specmap'][i],initvalue)
  string+='\n'
  # assume all rates positive
  for i in INFOS['rateset']:
    string+='assume(%s>0);\n' % (i)
  string+='\n'
  # solve equation system
  string+='A:desolve(eq,['
  for i in range(INFOS['nspec']):
    string+='%s(t)' % (INFOS['specmap'][i])
    if i+1==INFOS['nspec']:
      string+=']);\n'
    else:
      string+=','
  string+='\n'
  # format output
  string+='load(f90);\n\n'
  for i in range(INFOS['nspec']):
    string+='print("###")$\nf90(A[%i])$\n' % (i+1)
  # write MAXIMA input
  cwd=os.getcwd()
  infile=os.path.join(cwd,'MAXIMA.input')
  writefile(infile,string)

  # call MAXIMA in interactive batch mode
  string='maxima -b %s' % (infile)
  while True:
    print('Calling MAXIMA computer algebra system with following command:\n  %s' % (string))
    go_on=question('Run this command?',bool,True)
    if go_on:
      break
    else:
      print('Please enter the path to the MAXIMA computer algebra system:')
      path=question('Path to MAXIMA:',str,None,True)
      string=path+' -b %s' % (infile)
  print('\n')
  outfile=os.path.join(cwd,'MAXIMA.output')
  errfile=os.path.join(cwd,'MAXIMA.err')
  stdoutfile=open(outfile,'wb')
  stderrfile=open(errfile,'w')
  print('\n\nRunning MAXIMA interactively...')
  print()
  print('(If MAXIMA takes long, its output will be shown and STDIN switched to MAXIMA.\nPlease try to answer any questions that MAXIMA is prompting,\ne.g., "Is <some constant> positive, negative, or zero?"\nIn this case, answer with "pos;", "neg;", or "zero;")')
  print()
  print('*'*100)
  p=sp.Popen(string,shell=True,stdout=sp.PIPE,stderr=stderrfile,bufsize=1,stdin=sp.PIPE)
  isleep=0
  while True:
    isleep+=1
    try:
      time.sleep(0.5)
      if isleep>=3:
        while True:
          byte=p.stdout.read(1)
          if byte:
            sys.stdout.write(byte)
            sys.stdout.flush()
            stdoutfile.write(byte)
          else:
            break
      else:
        print('waiting ...')
    except KeyboardInterrupt:
      p.kill()
      print('*'*100)
      print('*** MAXIMA execution halted ***\n\n')
      raise
    if p.poll() is None:
      pass
    else:
      break
  while True:
    byte=p.stdout.read(1)
    if byte:
      stdoutfile.write(byte)
    else:
      break
  p.communicate()
  stdoutfile.close()
  stderrfile.close()
  print('*'*100)
  print('\n*** Done! ***')


  # extract function definitions from MAXIMA output
  out=readfile(outfile)
  functions=[]
  iline=-1
  active=False
  while True:
    iline+=1
    if iline>=len(out):
      if len(functions)==INFOS['nspec']:
        break
      else:
        print('*** Output does not contain all function definitions!')
        sys.exit(1)
    line=out[iline]
    if '###' in line[0:4]:
      for j in range(1,4):
        if '=' in out[iline+j]:
          iline+=j
          break
      string=''
      while True:
        line=out[iline]
        if line[0]=='&':
          line=line[1:]
        if not '&' in line:
          string+=line
          break
        else:
          line=line[:-2]
          string+=line
          iline+=1
      functions.append(string)
  #print functions

  # format as one string with proper newlines
  lenline=80
  string=''
  for function in functions:
    s=wrap_line(function,lenline)
    for i in range(len(s)):
      string+=s[i]
      if i+1==len(s):
        string+='\n\n'
      else:
        string+='\\\n'
  return string

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def get_global_function(species_groups,maxtime,timeshift):
  if len(species_groups)==1:
    string=''
    for i,spec in enumerate(species_groups[0]):
      string+='%s(x-%.2f)' % (spec,timeshift)
      if i+1<len(species_groups[0]):
        string+='+'
    return string
  elif len(species_groups)==2:
    string='x<%.2f ? ' % (maxtime+timeshift)
    for i,spec in enumerate(species_groups[0]):
      string+='%s(x-%.2f)' % (spec,timeshift)
      if i+1<len(species_groups[0]):
        string+='+'
    string+=' : '
    for i,spec in enumerate(species_groups[1]):
      string+='%s(x-%.2f)' % (spec,timeshift+maxtime)
      if i+1<len(species_groups[1]):
        string+='+'
    return string
  else:
    string='x<%.2f ? ' % (maxtime+timeshift)
    for i,spec in enumerate(species_groups[0]):
      string+='%s(x-%.2f)' % (spec,timeshift)
      if i+1<len(species_groups[0]):
        string+='+'
    string+=' : ( '
    string+=get_global_function(species_groups[1:],maxtime,timeshift+maxtime)
    string+=' ) '
    return string

# ===================================================

def get_global_using(columns_groups,maxtime,timeshift):
  if len(columns_groups)==1:
    string=''
    for i,col in enumerate(columns_groups[0]):
      string+='$%i' % (col)
      if i+1<len(columns_groups[0]):
        string+='+'
    return string
  elif len(columns_groups)==2:
    string='$1<%.2f ? ' % (maxtime+timeshift)
    for i,col in enumerate(columns_groups[0]):
      string+='$%i' % (col)
      if i+1<len(columns_groups[0]):
        string+='+'
    string+=' : '
    for i,col in enumerate(columns_groups[1]):
      string+='$%i' % (col)
      if i+1<len(columns_groups[1]):
        string+='+'
    return string
  else:
    string='$1<%.2f ? ' % (maxtime+timeshift)
    for i,col in enumerate(columns_groups[0]):
      string+='$%i' % (col)
      if i+1<len(columns_groups[0]):
        string+='+'
    string+=' : ( '
    string+=get_global_using(columns_groups[1:],maxtime,timeshift+maxtime)
    string+=' ) '
    return string

# ===================================================

def write_gnuscript(INFOS,functionstring):
  colorpalette=rgbcolor([INFOS['ngroups']])

  # header
  string ='#>>\n'
  string+='# +'+'-'*60+'+\n'
  string+='# |'+centerstring('Fitting script',60,' ')+'|\n'
  string+='# +'+'-'*60+'+\n#\n#\n'

  # add as comment the definitions of the kinetic model
  string+='# *** Definition of the kinetic model: ***\n'
  s=print_reactions(INFOS['rate_matrix'],INFOS['specmap']).splitlines()
  for line in s:
    string+='#' + line + '\n'
  string+='#\n#\n'
  string+='# *** Species and initial value: ***\n'
  for i in range(INFOS['nspec']):
    if INFOS['specmap'][i] in INFOS['initset']:
      initvalue=INFOS['specmap'][i]+'__0'
    else:
      initvalue='0'
    string+='# %s' % (INFOS['specmap'][i])+' '*(10-len(INFOS['specmap'][i]))+initvalue+'\n'
  string+='#\n#\n'
  string+='# *** Reaction rates: ***\n'
  for i in INFOS['rateset']:
    string+='# %s\n' % (i)
  string+='#\n#\n'

  # add the function definitions
  string+='# *** Species population function definitions: ***\n'
  string+=functionstring
  string+='\n\n# ========================================================\n'

  # add the rate constant guesses
  string+='# *** Reaction rates initial guesses: ***\n# These are given in units of inverse fs (e.g., time constant of 100 fs is written as 1./100.).\n# TODO: Please change to some suitable values!\n'
  for i,rate in enumerate(list(INFOS['rateset'])):
    string+='%s = 1./%i.\n' % (rate,100+i)
  string+='\n'

  # add the initial populations
  string+='# *** Species population initial guesses: ***\n# TODO: Please change to some suitable values!\n'
  for i in INFOS['initset']:
    string+='%s__0 = 1./%i.\n' % (i,len(INFOS['initset']))
  string+='\n'

  # add gnuplot global options
  string+='#<<\n'
  string+='# *** Gnuplot general options: ***\n'
  string+='''set xlabel "Time (fs)"
set ylabel "Population"
set xrange [0:%.2f]
set yrange [0:1]
set key at %.2f,1.00 top right
''' % (INFOS['maxtime'],INFOS['maxtime'])
  string+='\n'

  # add label with rates
  string+='# *** Label with time constants: ***\n'
  string+='set label 1 sprintf("'
  for i in INFOS['rateset']:
    string+='t_%s = %%7.1f fs\\n' % (i)
  for i in INFOS['initset']:
    string+='%s__0 = %%7.1f\\n' % (i)
  string+='"'
  for i in INFOS['rateset']:
    string+=',1./%s' % (i)
  for i in INFOS['initset']:
    string+=',%s__0' % (i)
  string+=') at %.2f,0.95\n' % (0.3*INFOS['maxtime'])
  string+='\n'

  # Initial plot
  string+='# *** Initial plot before fit: ***\n'
  string+='set title "Initial plot before fit\\nPress ENTER to setup global fit."\n'
  string+='p \\\n'
  for igroup,group in enumerate(INFOS['species_groups']):
    string+='  '
    for j in range(len(group)):
      string+='%s(x)' % (group[j])
      if j+1<len(group):
        string+='+'
      else:
        string+=' '
    string+='t "'
    for j in range(len(group)):
      string+='%s' % (group[j])
      if j+1<len(group):
        string+='+'
    string+='" w l lw 2 lc rgbcolor "%s", \\\n' % (colorpalette.hexcolor(1,igroup+1))
  for igroup,group in enumerate(INFOS['columns_groups']):
    string+='  "model_fit.dat" u 1:('
    for j in range(len(group)):
      string+='$%i' % (group[j])
      if j+1<len(group):
        string+='+'
    string+=') t "'
    for j in range(len(group)):
      string+='$%i' % (group[j])
      if j+1<len(group):
        string+='+'
    string+= '" w l lw 0.5 lc rgbcolor "%s"' % (colorpalette.hexcolor(1,igroup+1))
    if igroup+1<INFOS['ngroups']:
      string+=', \\\n'
    else:
      string+='\n'
  string+='\npause -1\n\n'

  # defining the global fit
  string+='#>>\n'
  string+='# *** Global fit setup: ***\n'
  string+='set xrange [0:%.2f]\nunset label\n' % (INFOS['ngroups']*INFOS['maxtime'])
  string+='set yrange [0:1]\n'
  string+='set key at %.2f,1.00 top right\n' % (INFOS['ngroups']*INFOS['maxtime'])
  string+='set fit logfile "model_fit.log"\n'

  # global fitting function
  string+='F(x)= '
  string+=get_global_function(INFOS['species_groups'],INFOS['maxtime'],0)
  string+='\n\n'
  usingstring=get_global_using(INFOS['columns_groups'],INFOS['maxtime'],0)

  # global plot
  string+='#<<\n'
  string+='# *** Global plot before fitting: ***\n'
  string+='set title "Global plot before fitting\\nPress ENTER to perform fit."\n'
  string+='p "model_fit.dat" u 1:( %s ) t "Data" w p pt 6 ps 0.5 lc rgbcolor "black", F(x) t "Fitting function" w l lw 2 lc rgbcolor "red"\npause -1\n' % (usingstring)
  string+='\n'

  # global fit
  string+='#>>\n'
  string+='# *** Execute global fit: ***\n# TODO: remove the "#" on the next line to also fit the initial populations.\n'
  viastring=''
  for i,rate in enumerate(list(INFOS['rateset'])):
    viastring+='%s' % (rate)
    if i+1<len(INFOS['rateset']):
      viastring+=','
  viastring+='#'
  for i,init in enumerate(list(INFOS['initset'])):
    viastring+=',%s__0' % (init)
  string+='fit F(x) "model_fit.dat" u 1:( %s ) via %s' % (usingstring,viastring)
  string+='\n\n'

  # global plot after fitting
  string+='#<<\n'
  string+='# *** Global plot after fitting: ***\n'
  string+='set title "Global plot after fitting\\nPress ENTER to display final results."\n'
  string+='p "model_fit.dat" u 1:( %s ) t "Data" w p pt 6 ps 0.5 lc rgbcolor "black", F(x) t "Fitting function" w l lw 2 lc rgbcolor "red"\npause -1\n' % (usingstring)
  string+='\n'

  # add gnuplot global options
  string+='#>>\n'
  string+='# *** Gnuplot general options: ***\n'
  string+='''set xlabel "Time (fs)"
set ylabel "Population"
set xrange [0:%.2f]
set yrange [0:1]
set key at %.2f,1.00 top right
''' % (INFOS['maxtime'],INFOS['maxtime'])
  string+='\n'

  # add label with rates
  string+='# *** Label with time constants: ***\n'
  string+='set label 1 sprintf("'
  for i in INFOS['rateset']:
    string+='t_%s = %%7.1f fs\\n' % (i)
  for i in INFOS['initset']:
    string+='%s__0 = %%7.1f\\n' % (i)
  string+='"'
  for i in INFOS['rateset']:
    string+=',1./%s' % (i)
  for i in INFOS['initset']:
    string+=',%s__0' % (i)
  string+=') at %.2f,0.95\n' % (0.3*INFOS['maxtime'])
  string+='\n'

  # Final plot
  string+='#<><> 1\n'
  string+='# *** Final plot after fit: ***\n'
  string+='set title "Final plot after fitting\\nPress ENTER to save as PNG and TXT."\n'
  string+='p \\\n'
  for igroup,group in enumerate(INFOS['species_groups']):
    string+='  '
    for j in range(len(group)):
      string+='%s(x)' % (group[j])
      if j+1<len(group):
        string+='+'
      else:
        string+=' '
    string+='t "'
    for j in range(len(group)):
      string+='%s' % (group[j])
      if j+1<len(group):
        string+='+'
    string+='" w l lw 2 lc rgbcolor "%s", \\\n' % (colorpalette.hexcolor(1,igroup+1))
  for igroup,group in enumerate(INFOS['columns_groups']):
    string+='  "model_fit.dat" u 1:('
    for j in range(len(group)):
      string+='$%i' % (group[j])
      if j+1<len(group):
        string+='+'
    string+=') t "'
    for j in range(len(group)):
      string+='$%i' % (group[j])
      if j+1<len(group):
        string+='+'
    string+= '" w l lw 0.5 lc rgbcolor "%s"' % (colorpalette.hexcolor(1,igroup+1))
    if igroup+1<INFOS['ngroups']:
      string+=', \\\n'
    else:
      string+='\n'
  string+='#<<\n'
  string+='\npause -1\n\n'
  string+='#>>\n'

  # For saving the final plot to a file and a table
  string+='# *** File output: ***\n'
  string+='#<<\n'
  string+='set title "%s@%s, %s\\nFile:%s"\n' % (os.environ['USER'],gethostname(),datetime.datetime.now(),INFOS['popfile'])
  string+='set term pngcairo size 800,480\nset out "model_fit.png"\nreplot\n'
  string+='#>>\n'
  string+='set table "model_fit.txt"\nreplot\n\n'

  # Writing the time constants and initial populations
  for i in INFOS['rateset']:
    string+='print "&&& %s : ", 1./%s\n' % (i,i)
  for i in INFOS['initset']:
    string+='print "&&& %s__0 : ", %s__0\n' % (i,i)

  # Put time stamp infos at the end
  string+='# *** Infos: ***\n'
  string+='# %s@%s\n' % (os.environ['USER'],gethostname())
  string+='# Date: %s\n' % (datetime.datetime.now())
  string+='# Current directory: %s\n\n' % (os.getcwd())

  # Write gnuplot script
  outfilename='model_fit.gp'
  print('Writing GNUPLOT script to %s ...' % (outfilename))
  writefile(outfilename,string)

  return

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def write_fitting_data(INFOS):
  string=''
  for igroup in range(INFOS['ngroups']):
    for step in INFOS['data']:
      line='%16.12f ' % (step[0]+igroup*INFOS['maxtime'])
      for d in step[1:]:
        line+='%16.12f ' % (d)
      line+='\n'
      string+=line

  # Write file
  outfilename='model_fit.dat'
  print('Writing GNUPLOT populations data to %s ...' % (outfilename))
  writefile(outfilename,string)

  return

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def print_messages(INFOS):
  string=''
  # decoration
  string+='\n'+centerstring('',62,'*')+'\n'
  string+='*'+centerstring(' Final Instructions and Hints ',60,' ')+'*\n'
  string+=centerstring('',62,'*')+'\n'

  # how to run the fitting script
  string+='\n'
  string+='''* How to run the fit *
  The fitting script is a script for GNUPLOT. Execute it with:
    user@host> gnuplot model_fit.gp
  Note that the file "model_fit.dat" must be present for the GNUPLOT script to work.
  By pressing ENTER, you can proceed through these stages:
    1) Initial plot of the population data and the fitting functions using the initial fitting parameters
    2) Plot of the data and functions (transformed for global fit) before fitting
    3) Perform the global fit and plot of the data and functions
    4) Final plot of the population data and the fitting functions, also presenting the final fit parameters
    5) Write the final plot to "model_fit.png" and write the data to "model_fit.txt"
'''

  # adjust parameters
  string+='\n'
  string+='* Adjust initial parameters *\n  Please adjust the starting guess values of the following fitting parameters:\n'
  for i in INFOS['rateset']:
    string+='    %s\n' % (i)
  for i in INFOS['initset']:
    string+='    %s__0\n' % (i)
  string+='  in file "model_fit.gp" before running the fitting procedure.\n'

  # whether initial populations should be fitted
  string+='\n'
  string+='* Fitting of initial populations *\n  The fitting script was setup to only optimize the rate constants in the global fit.\n  If you intend to also optimize the initial populations, please modify the fit command like this:\n'
  string+='    fit F(x) [...] via '
  for i in INFOS['rateset']:
    string+='%s,' % (i)
  for j,i in enumerate(list(INFOS['initset'])):
    string+='%s__0' % (i)
    if j+1<len(INFOS['initset']):
      string+=','
  string+='\n'

  # warning if not all species or columns are used
  unused_spec=[]
  for i in range(INFOS['nspec']):
    spec=INFOS['specmap'][i]
    used=False
    for j in INFOS['species_groups']:
      if spec in j:
        used=True
        break
    if not used:
      unused_spec.append(spec)
  unused_col=[]
  for i in range(2,1+INFOS['ncol']):
    used=False
    for j in INFOS['columns_groups']:
      if i in j:
        used=True
        break
    if not used:
      unused_col.append(i)
  if len(unused_col)>0 or len(unused_spec)>0:
    string+='\n'
    string+='**** Unused species or data columns ****\n'
    if len(unused_spec)>0:
      string+='  You defined species in your kinetic model which are not used in the global fit:\n'
      for i in unused_spec:
        string+='    %s\n' % (i)
    if len(unused_col)>0:
      string+='  The population data file contains columns which are not used in the global fit:\n'
      for i in unused_col:
        string+='    %i\n' % (i)
    string+='  Please ensure that this is in accord with your intentions.\n'

  # if cycles are in the model
  if INFOS['rank']==-1:
    string+='\n'
    string+='**** Cycles in reaction network ****\n'
    string+='  Without NUMPY, this script cannot check the reaction network for cycles.\n  Cycles might make the fit hard to converge or lead to large uncertainties for some fitting parameters.\n  Please carefully check your results.\n'
  if INFOS['rank']>0:
    string+='\n'
    string+='**** Cycles in reaction network ****\n'
    string+='  The reaction network defined contains %i cycles.\n  Cycles might make the fit hard to converge or lead to large uncertainties for some fitting parameters.\n  Please carefully check your results.\n' % (INFOS['rank'])

  string+=centerstring('',62,'*')+'\n'
  print(string)

  return

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def main():
  '''Main routine'''

  usage='''
python make_fitscript.py

This interactive script lets the user specify a kinetic model and a populations.py output file and produces a 
GNUPLOT script which allows to fit the model parameters to the populations
'''
  description=''
  parser = OptionParser(usage=usage, description=description)

  displaywelcome()
  open_keystrokes()

  # get input
  INFOS=get_infos()

  # echo input
  print('\n\n'+centerstring('Full input',60,'#')+'\n')
  for item in INFOS:
    if not item=='data':
      print(item, ' '*(25-len(item)), INFOS[item])
    elif item=='data':
      print(item, ' '*(25-len(item)), '[ ... ]')
  print('')
  go_on=question('Do you want to continue?',bool,True)
  if not go_on:
    quit(0)
  print('')

  # do work
  functionstring=get_functions_from_maxima(INFOS)
  write_gnuscript(INFOS,functionstring)
  write_fitting_data(INFOS)

  # display final messages
  print_messages(INFOS)

  # finalize
  close_keystrokes()

# ======================================================================================================================

if __name__ == '__main__':
  try:
    main()
  except KeyboardInterrupt:
    print('\nCtrl+C makes me a sad SHARC ;-(\n')
    quit(0)














