#!/usr/bin/env python2

#******************************************
#
#    SHARC Program Suite
#
#    Copyright (c) 2019 University of Vienna
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

# Interactive script for the setup of displacement calculations for SHARC
#
# usage: python setup_displacement.py

import sys
import re
import os
import stat
import shutil
from optparse import OptionParser

# for ordered dictionary output
from collections import OrderedDict
# to easily write/read data structure to/from file
#import pickle
import json


# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================


## checks for python 2 compatibility, defines any()/all() & sets version, versiondate, versionneeded
## contains dict & routines for the QM interfaces
## IO stuff (only make_directory)
## stdin/out stuff
## constants
#from modules import compatibility, interfaces, IO, output, constants


import sys 
import datetime

if sys.version_info[0] != 2:
  print 'This is a script for Python 2!'
  sys.exit(0)

if sys.version_info[1] < 5:
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

version = '2.1'
versionneeded = [0.2, 1.0, 2.0, 2.1, float(version)]
versiondate = datetime.date(2019,9,1)

# ======================================================================================================================

import os
#from output import question


def make_directory(displacement_dir):
  '''Creates a directory'''

  if os.path.isfile(displacement_dir):
    print '\nWARNING: %s is a file!' % (displacement_dir)
    return -1
  
  if os.path.isdir(displacement_dir):
    if len(os.listdir(displacement_dir)) == 0: return 0
    else:
      print '\nWARNING: %s/ is not empty!' % (displacement_dir)
      
      if not 'overwrite' in globals():
        global overwrite
        overwrite = question('Do you want to overwrite files in this folder? ', bool, False)
      
      if overwrite: return 0
      else: return -1
  else:
    try:
      os.mkdir(displacement_dir)
    except OSError:
      print '\nWARNING: %s cannot be created!' % (displacement_dir)
      return -1
    return 0

# ======================================================================================================================


import readline
import re

#from compatibility import *


def centerstring(string, n, pad=' '):
  l = len(string)
  if l >= n:
    return string
  else:
    return  pad * ((n - l + 1) / 2) + string + pad * ((n - l) / 2)


def displaywelcome():
  print 'Script for setup of LVC parametrization started...\n'
  string='\n'
  string+='  '+'='*80+'\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('LVC parametrization for SHARC dynamics',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Author: Simon Kropf and Sebastian Mai',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Version:'+version,80)+'||\n'
  string+='||'+centerstring(versiondate.strftime("%d.%m.%y"),80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  string+='This script automatizes the setup of excited-state calculations\nin order to parametrize LVC models for SHARC dynamics.'
  print string

def open_keystrokes():
  global KEYSTROKES
  KEYSTROKES=open('KEYSTROKES.tmp','w')

def close_keystrokes():
  KEYSTROKES.close()
  shutil.move('KEYSTROKES.tmp','KEYSTROKES.setup_LVCparam')


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


def print_INFOS(INFOS):
  print '\n' + centerstring('Full input', 60, '#') + '\n'
  
  for item in INFOS:
    i = 0 # counter for very long lists we do not want full output)

    if isinstance(INFOS[item], list) or isinstance(INFOS[item], tuple):
      first = True
      for elem in INFOS[item]:
        if i >= 10: break
        if first:
          print item, ' ' * (25 - len(item) - 1), elem
          first = False
        else:
          print ' ' * 25, elem

        i += 1


    elif isinstance(INFOS[item], dict):
      first = True
      for k, v in INFOS[item].items():
        if i >= 10: break
        if first:
          print item, ' ' * (25 - len(item)) + '%s: %s' % (k, v)
          first = False
        else:
          print ' ' * 25 + ' %s: %s' % (k, v)

        i += 1
    else: print item, ' ' * (25 - len(item) - 1), INFOS[item]

    if i >= 10:
      print ' ' * 25, '.'
      print ' ' * 25, '(' + str(len(INFOS[item]) - i) + ' more)'
      print ' ' * 25, '.'
  return



def reduce_big_list_to_short_str(big_list):
  '''
  Takes possibly big lists of numbers (e. g. list of all normal modes)
  and reduces them to short string for clean output to user

  e. g.: [7 8 9 12 13 14 17 20 21 22 23 25 26 28 29 30] => '[7~9 12~14 17 20~23 25 26 28~30]'

  returns shortened list string
  '''

  # if empty return [None]
  if not big_list: return [None]

  # start list_string, sort bit_list, set vars
  short_list_str = '('
  big_list = sorted(big_list)
  i_start, i_current, i_lastadded = 0, 0, 0
  
  # while index is within big_list
  while i_current < len(big_list) - 1:
    # check if next element is within range & continue
    if big_list[i_current] + 1 == big_list[i_current + 1]:
      i_current += 1
      continue
    # range ended - create shortened string
    else:
      # no range just one number alone
      if i_current == i_start:
        short_list_str += '%i ' % big_list[i_current]
      # range detected - shorten it
      else:
        # special case for 2 neighbour numbers
        if big_list[i_start] + 1 == big_list[i_current]:
          short_list_str += '%i %i ' % (big_list[i_start], big_list[i_current])
        else:
          short_list_str += '%i~%i ' % (big_list[i_start], big_list[i_current])

      # set vars accordingly for next run
      i_current += 1
      i_start = i_current
      i_lastadded = i_current

  # code above will always leave out last range/number - add it here (that's why we need i_lastadded) 
  if i_lastadded != i_current:
    # special case again for 2 neighbouring numbers
    if big_list[i_start] + 1 == big_list[i_current]:
      short_list_str += '%i %i' % (big_list[i_start], big_list[i_current])
    else:
      short_list_str += '%i~%i' % (big_list[i_start], big_list[i_current])
  else: short_list_str += str(big_list[i_current])
  
  # close bracket & return
  short_list_str += ')'
  return short_list_str


def reduce_displacement_dictionary_to_output_str(big_dictionary):
  '''
  used for shortening displacement dict therefore can never be empty

  e. g.
  for:
      OrderedDict({ 7: 0.05,  8: 0.05,  9: 0.05, 12: 0.05, 13: 0.05,
                   14: 0.04, 15: 0.05, 16: 0.05, 19: 0.03, 21: 0.03,
                   22: 0.03, 23: 0.03, 27: 0.01 })

  returns:
      [7~9 12 13 15 16]: 0.05
      [14]: 0.04
      [19 21~23]: 0.03
      [27]: 0.01

  returns string with ranges for same displacements
  '''
  output_str = ''

  # getting all different displacements
  displacement_list = []
  for normal_mode, displacement in big_dictionary.items():
    if not displacement in displacement_list: displacement_list.append(displacement)

  # running through all different displacements
  normal_mode_list_big = []
  for displacement in displacement_list:
    # adding all normal modes with same displacement to list
    for normal_mode, disp in big_dictionary.items():
      if displacement == disp:
        normal_mode_list_big.append(normal_mode)

    # reduce all normal modes with same displacement and add it to nice output
    output_str += '%s: %g\n' % (reduce_big_list_to_short_str(normal_mode_list_big), displacement)
    normal_mode_list_big = []

  return output_str

# ======================================================================================================================


from math import pi

CM_TO_HARTREE = 1./219474.6     # 4.556335252e-6 # conversion factor from cm-1 to Hartree
HARTREE_TO_EV = 27.211396132    # conversion factor from Hartree to eV
U_TO_AMU = 1./5.4857990943e-4   # conversion from g/mol to amu
ANG_TO_BOHR = 1./0.529177211    # 1.889725989    # conversion from Angstrom to bohr
PI = pi

pthresh = 1.e-5**2

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

NUMBERS = {'H':1, 'He':2,
'Li':3, 'Be':4, 'B':5, 'C':6,  'N':7,  'O':8, 'F':9, 'Ne':10,
'Na':11, 'Mg':12, 'Al':13, 'Si':14,  'P':15,  'S':16, 'Cl':17, 'Ar':18,
'K':19, 'Ca':20,
'Sc':21, 'Ti':22, 'V':23, 'Cr':24, 'Mn':25, 'Fe':26, 'Co':27, 'Ni':28, 'Cu':29, 'Zn':30,
'Ga':31, 'Ge':32, 'As':33, 'Se':34, 'Br':35, 'Kr':36,
'Rb':37, 'Sr':38,
'Y':39,  'Zr':40, 'Nb':41, 'Mo':42, 'Tc':43, 'Ru':44, 'Rh':45, 'Pd':46, 'Ag':47, 'Cd':48,
'In':49, 'Sn':50, 'Sb':51, 'Te':52,  'I':53, 'Xe':54,
'Cs':55, 'Ba':56,
'La':57, 
'Ce':58, 'Pr':59, 'Nd':60, 'Pm':61, 'Sm':62, 'Eu':63, 'Gd':64, 'Tb':65, 'Dy':66, 'Ho':67, 'Er':68, 'Tm':69, 'Yb':70, 'Lu':71,
         'Hf':72, 'Ta':73,  'W':74, 'Re':75, 'Os':76, 'Ir':77, 'Pt':78, 'Au':79, 'Hg':80,
'Tl':81, 'Pb':82, 'Bi':83, 'Po':84, 'At':85, 'Rn':86, 
'Fr':87, 'Ra':88,
'Ac':89, 
'Th':90, 'Pa':91,  'U':92, 'Np':93, 'Pu':94, 'Am':95, 'Cm':96, 'Bk':97, 'Cf':98, 'Es':99,'Fm':100,'Md':101,'No':102,'Lr':103,
        'Rf':104,'Db':105,'Sg':106,'Bh':107,'Hs':108,'Mt':109,'Ds':110,'Rg':111,'Cn':112,
'Nh':113,'Fl':114,'Mc':115,'Lv':116,'Ts':117,'Og':118
}


# Atomic Weights of the most common isotopes
# From https://chemistry.sciences.ncsu.edu/msf/pdf/IsotopicMass_NaturalAbundance.pdf
MASSES = {'H' :   1.007825 * U_TO_AMU,
          'He':   4.002603 * U_TO_AMU,
          'Li':   7.016004 * U_TO_AMU,
          'Be':   9.012182 * U_TO_AMU,
          'B' :  11.009305 * U_TO_AMU,
          'C' :  12.000000 * U_TO_AMU,
          'N' :  14.003074 * U_TO_AMU,
          'O' :  15.994915 * U_TO_AMU,
          'F' :  18.998403 * U_TO_AMU,
          'Ne':  19.992440 * U_TO_AMU,
          'Na':  22.989770 * U_TO_AMU,
          'Mg':  23.985042 * U_TO_AMU,
          'Al':  26.981538 * U_TO_AMU,
          'Si':  27.976927 * U_TO_AMU,
          'P' :  30.973762 * U_TO_AMU,
          'S' :  31.972071 * U_TO_AMU,
          'Cl':  34.968853 * U_TO_AMU,
          'Ar':  39.962383 * U_TO_AMU,
          'K' :  38.963707 * U_TO_AMU,
          'Ca':  39.962591 * U_TO_AMU,
          'Sc':  44.955910 * U_TO_AMU,
          'Ti':  47.947947 * U_TO_AMU,
          'V' :  50.943964 * U_TO_AMU,
          'Cr':  51.940512 * U_TO_AMU,
          'Mn':  54.938050 * U_TO_AMU,
          'Fe':  55.934942 * U_TO_AMU,
          'Co':  58.933200 * U_TO_AMU,
          'Ni':  57.935348 * U_TO_AMU,
          'Cu':  62.929601 * U_TO_AMU,
          'Zn':  63.929147 * U_TO_AMU,
          'Ga':  68.925581 * U_TO_AMU,
          'Ge':  73.921178 * U_TO_AMU,
          'As':  74.921596 * U_TO_AMU,
          'Se':  79.916522 * U_TO_AMU,
          'Br':  78.918338 * U_TO_AMU,
          'Kr':  83.911507 * U_TO_AMU,
          'Rb':  84.911789 * U_TO_AMU,
          'Sr':  87.905614 * U_TO_AMU,
          'Y' :  88.905848 * U_TO_AMU,
          'Zr':  89.904704 * U_TO_AMU,
          'Nb':  92.906378 * U_TO_AMU,
          'Mo':  97.905408 * U_TO_AMU,
          'Tc':  98.907216 * U_TO_AMU,
          'Ru': 101.904350 * U_TO_AMU,
          'Rh': 102.905504 * U_TO_AMU,
          'Pd': 105.903483 * U_TO_AMU,
          'Ag': 106.905093 * U_TO_AMU,
          'Cd': 113.903358 * U_TO_AMU,
          'In': 114.903878 * U_TO_AMU,
          'Sn': 119.902197 * U_TO_AMU,
          'Sb': 120.903818 * U_TO_AMU,
          'Te': 129.906223 * U_TO_AMU,
          'I' : 126.904468 * U_TO_AMU,
          'Xe': 131.904154 * U_TO_AMU,
          'Cs': 132.905447 * U_TO_AMU,
          'Ba': 137.905241 * U_TO_AMU,
          'La': 138.906348 * U_TO_AMU,
          'Ce': 139.905435 * U_TO_AMU,
          'Pr': 140.907648 * U_TO_AMU,
          'Nd': 141.907719 * U_TO_AMU,
          'Pm': 144.912744 * U_TO_AMU,
          'Sm': 151.919729 * U_TO_AMU,
          'Eu': 152.921227 * U_TO_AMU,
          'Gd': 157.924101 * U_TO_AMU,
          'Tb': 158.925343 * U_TO_AMU,
          'Dy': 163.929171 * U_TO_AMU,
          'Ho': 164.930319 * U_TO_AMU,
          'Er': 165.930290 * U_TO_AMU,
          'Tm': 168.934211 * U_TO_AMU,
          'Yb': 173.938858 * U_TO_AMU,
          'Lu': 174.940768 * U_TO_AMU,
          'Hf': 179.946549 * U_TO_AMU,
          'Ta': 180.947996 * U_TO_AMU,
          'W' : 183.950933 * U_TO_AMU,
          'Re': 186.955751 * U_TO_AMU,
          'Os': 191.961479 * U_TO_AMU,
          'Ir': 192.962924 * U_TO_AMU,
          'Pt': 194.964774 * U_TO_AMU,
          'Au': 196.966552 * U_TO_AMU,
          'Hg': 201.970626 * U_TO_AMU,
          'Tl': 204.974412 * U_TO_AMU,
          'Pb': 207.976636 * U_TO_AMU,
          'Bi': 208.980383 * U_TO_AMU,
          'Po': 208.982416 * U_TO_AMU,
          'At': 209.987131 * U_TO_AMU,
          'Rn': 222.017570 * U_TO_AMU,
          'Fr': 223.019731 * U_TO_AMU, 
          'Ra': 226.025403 * U_TO_AMU,
          'Ac': 227.027747 * U_TO_AMU, 
          'Th': 232.038050 * U_TO_AMU, 
          'Pa': 231.035879 * U_TO_AMU, 
          'U' : 238.050783 * U_TO_AMU, 
          'Np': 237.048167 * U_TO_AMU,
          'Pu': 244.064198 * U_TO_AMU, 
          'Am': 243.061373 * U_TO_AMU, 
          'Cm': 247.070347 * U_TO_AMU, 
          'Bk': 247.070299 * U_TO_AMU, 
          'Cf': 251.079580 * U_TO_AMU, 
          'Es': 252.082972 * U_TO_AMU,
          'Fm': 257.095099 * U_TO_AMU,
          'Md': 258.098425 * U_TO_AMU,
          'No': 259.101024 * U_TO_AMU,
          'Lr': 262.109692 * U_TO_AMU,
          'Rf': 267. * U_TO_AMU,
          'Db': 268. * U_TO_AMU,
          'Sg': 269. * U_TO_AMU,
          'Bh': 270. * U_TO_AMU,
          'Hs': 270. * U_TO_AMU,
          'Mt': 278. * U_TO_AMU,
          'Ds': 281. * U_TO_AMU,
          'Rg': 282. * U_TO_AMU,
          'Cn': 285. * U_TO_AMU,
          'Nh': 286. * U_TO_AMU,
          'Fl': 289. * U_TO_AMU,
          'Mc': 290. * U_TO_AMU,
          'Lv': 293. * U_TO_AMU,
          'Ts': 294. * U_TO_AMU,
          'Og': 294. * U_TO_AMU
}


# Isotopes used for the masses
ISOTOPES={'H' : 'H-1' ,
          'He': 'He-4',
          'Li': 'Li-7',
          'Be': 'Be-9*',
          'B' : 'B_11' ,
          'C' : 'C-12' ,
          'N' : 'N-14' ,
          'O' : 'O-16' ,
          'F' : 'F-19*' ,
          'Ne': 'Ne-20',
          'Na': 'Na-23*',
          'Mg': 'Mg-24',
          'Al': 'Al-27*',
          'Si': 'Si-28',
          'P' : 'P-31*' ,
          'S' : 'S-32' ,
          'Cl': 'Cl-35',
          'Ar': 'Ar-40',
          'K' : 'K-39' ,
          'Ca': 'Ca-40',
          'Sc': 'Sc-45*',
          'Ti': 'Ti-48',
          'V' : 'V-51' ,
          'Cr': 'Cr-52',
          'Mn': 'Mn-55*',
          'Fe': 'Fe-56',
          'Co': 'Co-59*',
          'Ni': 'Ni-58',
          'Cu': 'Cu-63',
          'Zn': 'Zn-64',
          'Ga': 'Ga-69',
          'Ge': 'Ge-74',
          'As': 'As-75*',
          'Se': 'Se-80',
          'Br': 'Br-79',
          'Kr': 'Kr-84',
          'Rb': 'Rb-85',
          'Sr': 'Sr-88',
          'Y' : 'Y-89*' ,
          'Zr': 'Zr-90',
          'Nb': 'Nb-93*',
          'Mo': 'Mo-98',
          'Tc': 'Tc-98',
          'Ru': 'Ru-102',
          'Rh': 'Rh-103*',
          'Pd': 'Pd-106',
          'Ag': 'Ag-107',
          'Cd': 'Cd-114',
          'In': 'In-115',
          'Sn': 'Sn-120',
          'Sb': 'Sb-121',
          'Te': 'Te-130',
          'I' : 'I-127*' ,
          'Xe': 'Xe-132',
          'Cs': 'Cs-133*',
          'Ba': 'Ba-138',
          'La': 'La-139',
          'Ce': 'Ce-140',
          'Pr': 'Pr-141*',
          'Nd': 'Nd-142',
          'Pm': 'Pm-145',
          'Sm': 'Sm-152',
          'Eu': 'Eu-153',
          'Gd': 'Gd-158',
          'Tb': 'Tb-159*',
          'Dy': 'Dy-164',
          'Ho': 'Ho-165*',
          'Er': 'Er-166',
          'Tm': 'Tm-169*',
          'Yb': 'Yb-174',
          'Lu': 'Lu-175',
          'Hf': 'Hf-180',
          'Ta': 'Ta-181',
          'W' : 'W-184' ,
          'Re': 'Re-187',
          'Os': 'Os-192',
          'Ir': 'Ir-193',
          'Pt': 'Pt-195',
          'Au': 'Au-197*',
          'Hg': 'Hg-202',
          'Tl': 'Tl-205',
          'Pb': 'Pb-208',
          'Bi': 'Bi-209*',
          'Po': 'Po-209',
          'At': 'At-210',
          'Rn': 'Rn-222',
          'Fr': 'Fr-223', 
          'Ra': 'Ra-226',
          'Ac': 'Ac-227', 
          'Th': 'Th-232*', 
          'Pa': 'Pa-231*', 
          'U' : 'U-238' , 
          'Np': 'Np-237',
          'Pu': 'Pu-244', 
          'Am': 'Am-243', 
          'Cm': 'Cm-247', 
          'Bk': 'Bk-247', 
          'Cf': 'Cf-251', 
          'Es': 'Es-252',
          'Fm': 'Fm-257',
          'Md': 'Md-258',
          'No': 'No-259',
          'Lr': 'Lr-262',
          'Rf': 'Rf-267',
          'Db': 'Db-268',
          'Sg': 'Sg-269',
          'Bh': 'Bh-270',
          'Hs': 'Hs-270',
          'Mt': 'Mt-278',
          'Ds': 'Ds-281',
          'Rg': 'Rg-282',
          'Cn': 'Cn-285',
          'Nh': 'Nh-286',
          'Fl': 'Fl-289',
          'Mc': 'Mc-290',
          'Lv': 'Lv-293',
          'Ts': 'Ts-294',
          'Og': 'Og-294'
}

# ======================================================================================================================


import re
import os
import ast
import time
import shutil

#from output import *


Interfaces={
  1: {'script':          'SHARC_MOLPRO.py',
      'name':            'molpro',
      'description':     'MOLPRO (only CASSCF)',
      'get_routine':     'get_MOLPRO',
      'prepare_routine': 'prepare_MOLPRO',
      'features':        {'overlap': ['wfoverlap'],
                          'dyson':   ['wfoverlap'],
                          'nacdr':   ['wfoverlap'],
                          'phases':  ['wfoverlap'],
                          'soc':     []             },
      'pysharc':          False
     },
  2: {'script':          'SHARC_COLUMBUS.py',
      'name':            'columbus',
      'description':     'COLUMBUS (CASSCF, RASSCF and MRCISD), using SEWARD integrals',
      'get_routine':     'get_COLUMBUS',
      'prepare_routine': 'prepare_COLUMBUS',
      'features':        {'overlap': ['wfoverlap'],
                          'dyson':   ['wfoverlap'],
                          'phases':  ['wfoverlap'],
                          'nacdr':   [],
                          'soc':     []               },
      'pysharc':          False
     },
  3: {'script':          'SHARC_Analytical.py',
      'name':            'analytical',
      'description':     'Analytical PESs',
      'get_routine':     'get_Analytical',
      'prepare_routine': 'prepare_Analytical',
      'features':        {'overlap': [],
                          'dipolegrad':[],
                          'phases':  [],
                          'soc':     []             },
      'pysharc':          False
     },
  4: {'script':          'SHARC_MOLCAS.py',
      'name':            'molcas',
      'description':     'MOLCAS (CASSCF, CASPT2, MS-CASPT2)',
      'get_routine':     'get_MOLCAS',
      'prepare_routine': 'prepare_MOLCAS',
      'features':        {'overlap': [],
                          'dyson':   ['wfoverlap'],
                          'dipolegrad':[],
                          'phases':  [],
                          'nacdr':   [],
                          'soc':     []             },
      'pysharc':          False
     },
  5: {'script':          'SHARC_ADF.py',
      'name':            'adf',
      'description':     'ADF (DFT, TD-DFT)',
      'get_routine':     'get_ADF',
      'prepare_routine': 'prepare_ADF',
      'features':        {'overlap': ['wfoverlap'],
                          'dyson':   ['wfoverlap'],
                          'theodore':['theodore'],
                          'phases':  ['wfoverlap'],
                          'soc':     []                 },
      'pysharc':          False
     },
  6: {'script':          'SHARC_RICC2.py',
      'name':            'ricc2',
      'description':     'TURBOMOLE (ricc2 with CC2 and ADC(2))',
      'get_routine':     'get_RICC2',
      'prepare_routine': 'prepare_RICC2',
      'features':        {'overlap': ['wfoverlap'],
                          'theodore':['theodore'],
                          'phases':  ['wfoverlap'],
                          'soc':     []                 },
      'pysharc':          False
     },
  7: {'script':          'SHARC_LVC.py',
      'name':            'lvc',
      'description':     'LVC Hamiltonian',
      'get_routine':     'get_LVC',
      'prepare_routine': 'prepare_LVC',
      'features':        {'overlap': [],
                          'nacdr':   [],
                          'phases':  [],
                          'soc':     []                 },
      'pysharc':          True,
      'pysharc_driver':   'pysharc_lvc.py'
     },
  8: {'script':          'SHARC_GAUSSIAN.py',
      'name':            'gaussian',
      'description':     'GAUSSIAN (DFT, TD-DFT)',
      'get_routine':     'get_GAUSSIAN',
      'prepare_routine': 'prepare_GAUSSIAN',
      'features':        {'overlap': ['wfoverlap'],
                          'dyson':   ['wfoverlap'],
                          'theodore':['theodore'],
                          'phases':  ['wfoverlap']        },
      'pysharc':          False
     },
  9: {'script':          'SHARC_ORCA.py',
      'name':            'orca',
      'description':     'ORCA (DFT, TD-DFT, HF, CIS)',
      'get_routine':     'get_ORCA',
      'prepare_routine': 'prepare_ORCA',
      'features':        {'overlap': ['wfoverlap'],
                          'dyson':   ['wfoverlap'],
                          'theodore':['theodore'],
                          'phases':  ['wfoverlap'],
                          'soc':     []},
      'pysharc':          False
     },
  10:{'script':          'SHARC_BAGEL.py',
      'name':            'bagel',
      'description':     'BAGEL (CASSCF, CASPT2, (X)MS-CASPT2)',
      'get_routine':     'get_BAGEL',
      'prepare_routine': 'prepare_BAGEL',
      'features':        {'overlap': ['wfoverlap'],
                          'dyson':   ['wfoverlap'],
                          'nacdr' : [],
                          'dipolegrad':[],
                          'phases':  [],                     },
      'pysharc':          False
     }
  }



def checktemplate_MOLPRO(filename):
  necessary=['basis','closed','occ','nelec','roots']
  try:
    f=open(filename)
    data=f.readlines()
    f.close()
  except IOError:
    print 'Could not open template file %s' % (filename)
    return False
  i=0
  for l in data:
    if necessary[i] in l:
      i+=1
      if i+1==len(necessary):
        return True
  print 'The template %s seems to be incomplete! It should contain: ' % (filename) +str(necessary)
  return False

# =================================================

def get_MOLPRO(INFOS):
  '''This routine asks for all questions specific to MOLPRO:
  - path to molpro
  - scratch directory
  - MOLPRO.template
  - wf.init
  '''

  string='\n  '+'='*80+'\n'
  string+='||'+centerstring('MOLPRO Interface setup',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  print string

  print centerstring('Path to MOLPRO',60,'-')+'\n'
  path=os.getenv('MOLPRO')
  path=os.path.expanduser(os.path.expandvars(path))
  if not path=='':
    path='$MOLPRO/'
  else:
    path=None
  #if path!='':
    #print 'Environment variable $MOLPRO detected:\n$MOLPRO=%s\n' % (path)
    #if question('Do you want to use this MOLPRO installation?',bool,True):
      #INFOS['molpro']=path
  #if not 'molpro' in INFOS:
  print '\nPlease specify path to MOLPRO directory (SHELL variables and ~ can be used, will be expanded when interface is started).\n'
  INFOS['molpro']=question('Path to MOLPRO executable:',str,path)
  print ''


  print centerstring('Scratch directory',60,'-')+'\n'
  print 'Please specify an appropriate scratch directory. This will be used to temporally store the integrals. The scratch directory will be deleted after the calculation. Remember that this script cannot check whether the path is valid, since you may run the calculations on a different machine. The path will not be expanded by this script.'
  INFOS['scratchdir']=question('Path to scratch directory:',str)
  print ''


  print centerstring('MOLPRO input template file',60,'-')+'\n'
  print '''Please specify the path to the MOLPRO.template file. This file must be a valid MOLPRO input file for a CASSCF calculation. It should contain the following settings:
- memory settings
- Basis set (possibly also Douglas-Kroll settings etc.)
- CASSCF calculation with:
  * Number of frozen, closed and occupied orbitals
  * wf and state cards for the specification of the wavefunction
MOLPRO.template files can easily be created using molpro_input.py (Open a second shell if you need to create one now).

The MOLPRO interface will generate the remaining MOLPRO input automatically.
'''
  if os.path.isfile('MOLPRO.template'):
    if checktemplate_MOLPRO('MOLPRO.template'):
      print 'Valid file "MOLPRO.template" detected. '
      usethisone=question('Use this template file?',bool,True)
      if usethisone:
        INFOS['molpro.template']='MOLPRO.template'
  if not 'molpro.template' in INFOS:
    while True:
      filename=question('Template filename:',str)
      if not os.path.isfile(filename):
        print 'File %s does not exist!' % (filename)
        continue
      if checktemplate_MOLPRO(filename):
        break
    INFOS['molpro.template']=filename
  print ''


  print centerstring('Initial wavefunction: MO Guess',60,'-')+'\n'
  print '''Please specify the path to a MOLPRO wavefunction file containing suitable starting MOs for the CASSCF calculation. Please note that this script cannot check whether the wavefunction file and the Input template are consistent!

If you optimized your geometry with MOLPRO/CASSCF you can reuse the "wf" file from the optimization.
'''
  if question('Do you have an initial wavefunction file?',bool,True):
    while True:
      filename=question('Initial wavefunction file:',str,'wf.init')
      if os.path.isfile(filename):
        break
      else:
        print 'File not found!'
    INFOS['molpro.guess']=filename
  else:
    print 'WARNING: Remember that CASSCF calculations may run very long and/or yield wrong results without proper starting MOs.'
    time.sleep(2)
    INFOS['molpro.guess']=False


  print centerstring('MOLPRO Ressource usage',60,'-')+'\n'
  print '''Please specify the amount of memory available to MOLPRO (in MB). For calculations including moderately-sized CASSCF calculations and less than 150 basis functions, around 2000 MB should be sufficient.
'''
  INFOS['molpro.mem']=abs(question('MOLPRO memory:',int,[500])[0])
  print '''Please specify the number of CPUs to be used by EACH calculation.
'''
  INFOS['molpro.ncpu']=abs(question('Number of CPUs:',int,[1])[0])



  if 'wfoverlap' in INFOS['needed']:
    print '\n'+centerstring('WFoverlap setup',60,'-')+'\n'
    INFOS['molpro.wfpath']=question('Path to wavefunction overlap executable:',str,'$SHARC/wfoverlap.x')
    # TODO: not asked for: numfrozcore, numocc

  return INFOS

# =================================================

def prepare_MOLPRO(INFOS,iconddir):
  # write MOLPRO.resources
  try:
    sh2pro=open('%s/MOLPRO.resources' % (iconddir), 'w')
  except IOError:
    print 'IOError during prepareMOLPRO, iconddir=%s' % (iconddir)
    quit(1)
  string='''molpro %s
scratchdir %s/%s/

memory %i
ncpu %i
''' % (INFOS['molpro'],INFOS['scratchdir'],iconddir,INFOS['molpro.mem'],INFOS['molpro.ncpu'])
  if 'wfoverlap' in INFOS['needed']:
    string+='wfoverlap %s\n' % (INFOS['molpro.wfpath'])
  else:
    string+='\nnooverlap\n'
  sh2pro.write(string)
  sh2pro.close()

  # copy MOs and template
  cpfrom=INFOS['molpro.template']
  cpto='%s/MOLPRO.template' % (iconddir)
  shutil.copy(cpfrom,cpto)
  if INFOS['molpro.guess']:
    cpfrom=INFOS['molpro.guess']
    cpto='%s/wf.init' % (iconddir)
    shutil.copy(cpfrom,cpto)

  return

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def checktemplate_COLUMBUS(TEMPLATE, mult):
  '''Checks whether TEMPLATE is a file or directory. If a file or does not exist, it quits with exit code 1, if it is a directory, it checks whether all important input files are there. Does not check for all input files, since runc does this, too.

  Arguments:
  1 string: path to TEMPLATE

  returns whether input is for isc keyword or socinr keyword
  and returns the DRT of the given multiplicity'''

  exist=os.path.exists(TEMPLATE)
  if exist:
    isfile=os.path.isfile(TEMPLATE)
    if isfile:
      #print 'TEMPLATE=%s exists and is a file!' % (TEMPLATE)
      return None,None,None
    necessary=['control.run','mcscfin','tranin','propin']
    lof=os.listdir(TEMPLATE)
    for i in necessary:
      if not i in lof:
        #print 'Did not find input file %s! Did you prepare the input according to the instructions?' % (i)
        return None,None,None
    cidrtinthere=False
    ciudginthere=False
    for i in lof:
      if 'cidrtin' in i:
        cidrtinthere=True
      if 'ciudgin' in i:
        ciudginthere=True
    if not cidrtinthere or not ciudginthere:
      #print 'Did not find input file %s.*! Did you prepare the input according to the instructions?' % (i)
      return None,None,None
  else:
    #print 'Directory %s does not exist!' % (TEMPLATE)
    return None,None,None

  # get integral program
  try:
    intprog=open(TEMPLATE+'/intprogram')
    line=intprog.readline()
    if 'hermit' in line:
      INTPROG='dalton'
    elif 'seward' in line:
      INTPROG='seward'
    else:
      return None,None,None
  except IOError:
    return None,None,None

  # check cidrtin and cidrtin* for the multiplicity
  try:
    cidrtin=open(TEMPLATE+'/cidrtin')
    line=cidrtin.readline().split()
    if line[0].lower()=='y':
      maxmult=int(cidrtin.readline().split()[0])
      cidrtin.readline()
      nelec=int(cidrtin.readline().split()[0])
      if mult<=maxmult and (mult+nelec)%2!=0:
        return 1, (mult+1)/2,INTPROG    # socinr=1, single=-1, isc=0
      else:
        return None,None,None
    else:
      mult2=int(cidrtin.readline().split()[0])
      if mult!=mult2:
        #print 'Multiplicity %i cannot be treated in directory %s (single DRT)!'  % (mult,TEMPLATE)
        return None,None,None
      return -1,1,INTPROG
  except IOError:
    # find out in which DRT the requested multiplicity is
    for i in range(1,9):        # COLUMBUS can treat at most 8 DRTs
      try:
        cidrtin=open(TEMPLATE+'/cidrtin.%i' % i)
      except IOError:
        return None,None,None
      cidrtin.readline()
      mult2=int(cidrtin.readline().split()[0])
      if mult==mult2:
        return 0,i,INTPROG
      cidrtin.close()

# =================================================

def get_COLUMBUS(INFOS):
  '''This routine asks for all questions specific to COLUMBUS:
  - path to COLUMBUS
  - scratchdir
  - path to template directory
  - mocoef
  - memory
  '''

  string='\n  '+'='*80+'\n'
  string+='||'+centerstring('COLUMBUS Interface setup',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  print string


  print centerstring('Path to COLUMBUS',60,'-')+'\n'
  path=os.getenv('COLUMBUS')
  if path=='':
    path=None
  else:
    path='$COLUMBUS/'
  #path=os.path.expanduser(os.path.expandvars(path))
  #if path!='':
    #print 'Environment variable $COLUMBUS detected:\n$COLUMBUS=%s\n' % (path)
    #if question('Do you want to use this COLUMBUS installation?',bool,True):
      #INFOS['columbus']=path
  #if not 'columbus' in INFOS:
  print '\nPlease specify path to COLUMBUS directory (SHELL variables and ~ can be used, will be expanded when interface is started).\n'
  INFOS['columbus']=question('Path to COLUMBUS:',str,path)
  print ''


  print centerstring('Scratch directory',60,'-')+'\n'
  print 'Please specify an appropriate scratch directory. This will be used to temporally store all COLUMBUS files. The scratch directory will be deleted after the calculation. Remember that this script cannot check whether the path is valid, since you may run the calculations on a different machine. The path will not be expanded by this script.'
  INFOS['scratchdir']=question('Path to scratch directory:',str)
  print ''


  print centerstring('COLUMBUS input template directory',60,'-')+'\n'
  print '''Please specify the path to the COLUMBUS template directory.
The directory must contain subdirectories with complete COLUMBUS input file sets for the following steps:
- Integrals with SEWARD/MOLCAS
- SCF
- MCSCF
- SO-MRCI (even if no Spin-Orbit couplings will be calculated)
The COLUMBUS interface will generate the remaining COLUMBUS input automatically, depending on the number of states.

In order to setup the COLUMBUS input, use COLUMBUS' input facility colinp. For further information, see the Spin-orbit tutorial for COLUMBUS [1].

[1] http://www.univie.ac.at/columbus/docs_COL70/tutorial-SO.pdf
'''
  while True:
    path=question('Path to templates:',str)
    path=os.path.expanduser(os.path.expandvars(path))
    path=os.path.abspath(path)
    if not os.path.isdir(path):
      print 'Directory %s does not exist!' % (path)
      continue

    content=os.listdir(path)
    multmap={}
    allOK=True
    for mult in range(1,1+len(INFOS['states'])):
      if INFOS['states'][mult-1]==0:
        continue
      found=False
      for d in content:
        template=path+'/'+d
        socitype,drt,intprog=checktemplate_COLUMBUS(template,mult)
        if socitype==None:
          continue
        if not d[-1]=='/':
          d+='/'
        multmap[mult]=d
        found=True
        break
      if not found:
        print 'No input directory for multiplicity %i!' % (mult)
        allOK=False
        continue
    if allOK:
      break
  print '\nAccepted path: %s\n' % (path)

  print '''Check whether the jobs are assigned correctly to the multiplicities. Use the following commands:
  mult job        make <mult> use the input in <job>
  show            show the mapping of multiplicities to jobs
  end             confirm this mapping
'''
  for i in multmap:
    print '%i ==> %s' % (i,multmap[i])
  while True:
    line=question('Adjust job mapping:',str,'end',False)
    if 'show' in line.lower():
      for i in multmap:
        print '%i ==> %s' % (i,multmap[i])
      continue
    elif 'end' in line.lower():
      break
    else:
      f=line.split()
      try:
        m=int(f[0])
        j=f[1]
      except (ValueError,IndexError):
        continue
      if not m in multmap:
        print 'Multiplicity %i not necessary!' % (m)
        continue
      if not os.path.isdir(path+'/'+j):
        print 'No template subdirectory %s!' % (j)
        continue
      if not j[-1]=='/':
        j+='/'
      multmap[m]=j
  print ''

  mocoefmap={}
  for job in set([ multmap[i] for i in multmap]):
    mocoefmap[job]=multmap[min(multmap)]
  print '''Check whether the mocoeffiles are assigned correctly to the jobs. Use the following commands:
  job mocoefjob   make <job> use the mocoeffiles from <mocoefjob>
  show            show the mapping of multiplicities to jobs
  end             confirm this mapping
'''
  width=max([ len(i) for i in mocoefmap] )
  for i in mocoefmap:
    print '%s' % (i) +' '*(width-len(i))+ ' <== %s' % (mocoefmap[i])
  while True:
    line=question('Adjust mocoef mapping:',str,'end',False)
    if 'show' in line.lower():
      for i in mocoefmap:
        print '%s <== %s' % (i,mocoefmap[i])
      continue
    elif 'end' in line.lower():
      break
    else:
      f=line.split()
      try:
        j=f[0]
        m=f[1]
      except (ValueError,IndexError):
        continue
      if not m[-1]=='/':
        m+='/'
      if not j[-1]=='/':
        j+='/'
      mocoefmap[j]=m
  print ''

  INFOS['columbus.template']=path
  INFOS['columbus.multmap']=multmap
  INFOS['columbus.mocoefmap']=mocoefmap
  INFOS['columbus.intprog']=intprog

  INFOS['columbus.copy_template']=question('Do you want to copy the template directory to each trajectory (Otherwise it will be linked)?',bool,False)
  if INFOS['columbus.copy_template']:
    INFOS['columbus.copy_template_from']=INFOS['columbus.template']
    INFOS['columbus.template']='./COLUMBUS.template/'


  print centerstring('Initial wavefunction: MO Guess',60,'-')+'\n'
  print '''Please specify the path to a COLUMBUS mocoef file containing suitable starting MOs for the CASSCF calculation.
'''
  init=question('Do you have an initial mocoef file?',bool,True)
  if init:
    while True:
      line=question('Mocoef filename:',str,'mocoef_mc.init')
      line=os.path.expanduser(os.path.expandvars(line))
      if os.path.isfile(line):
          break
      else:
        print 'File not found!'
        continue
    INFOS['columbus.guess']=line
  else:
    print 'WARNING: Remember that CASSCF calculations may run very long and/or yield wrong results without proper starting MOs.'
    time.sleep(2)
    INFOS['columbus.guess']=False
  print ''


  print centerstring('COLUMBUS Memory usage',60,'-')+'\n'
  print '''Please specify the amount of memory available to COLUMBUS (in MB). For calculations including moderately-sized CASSCF calculations and less than 150 basis functions, around 2000 MB should be sufficient.
'''
  INFOS['columbus.mem']=abs(question('COLUMBUS memory:',int)[0])


  # Ionization
  #print '\n'+centerstring('Ionization probability by Dyson norms',60,'-')+'\n'
  #INFOS['ion']=question('Dyson norms?',bool,False)
  #if INFOS['ion']:
  if 'wfoverlap' in INFOS['needed']:
    print '\n'+centerstring('WFoverlap setup',60,'-')+'\n'
    INFOS['columbus.dysonpath']=question('Path to wavefunction overlap executable:',str,'$SHARC/wfoverlap.x')
    INFOS['columbus.ciothres']=question('Determinant screening threshold:',float,[0.97])[0]
    INFOS['columbus.numfrozcore']=question('Number of frozen core orbitals for overlaps (-1=as in template):',int,[-1])[0]
    if 'ion' in INFOS and INFOS['ion']:
      INFOS['columbus.numocc']=question('Number of doubly occupied orbitals for Dyson:',int,[0])[0]

  return INFOS

# =================================================

def prepare_COLUMBUS(INFOS,iconddir):
  # write COLUMBUS.resources
  try:
    sh2col=open('%s/COLUMBUS.resources' % (iconddir), 'w')
  except IOError:
    print 'IOError during prepareCOLUMBUS, directory=%i' % (iconddir)
    quit(1)
  string= 'columbus %s\nscratchdir %s/%s/WORK\n' % (INFOS['columbus'],INFOS['scratchdir'],iconddir)
  string+='savedir %s/%s/savedir\ntemplate %s\nmemory %i\n\n' % (INFOS['scratchdir'],iconddir, INFOS['columbus.template'],INFOS['columbus.mem'])
  string+='integrals %s\n' % (INFOS['columbus.intprog'])
  for mult in INFOS['columbus.multmap']:
    string+='DIR %i %s\n' % (mult,INFOS['columbus.multmap'][mult])
  string+='\n'
  for job in INFOS['columbus.mocoefmap']:
    string+='MOCOEF %s %s\n' % (job,INFOS['columbus.mocoefmap'][job])
  if 'wfoverlap' in INFOS['needed']:
    string+='wfoverlap %s\n' % (INFOS['columbus.dysonpath'])
    string+='wfthres %s\n' % (INFOS['columbus.ciothres'])
    if INFOS['columbus.numfrozcore']>=0:
      string+='numfrozcore %i\n' % (INFOS['columbus.numfrozcore'])
    if 'columbus.numocc' in INFOS:
      string+='numocc %i\n' % (INFOS['columbus.numocc'])
  else:
    string+='\nnooverlap\n'
  sh2col.write(string)
  sh2col.close()

  # copy MOs and template
  if INFOS['columbus.guess']:
    cpfrom=INFOS['columbus.guess']
    cpto='%s/mocoef_mc.init' % (iconddir)
    shutil.copy(cpfrom,cpto)

  if INFOS['columbus.copy_template']:
    copy_from=INFOS['columbus.copy_template_from']
    copy_to=iconddir+'/COLUMBUS.template/'
    if os.path.exists(copy_to):
      shutil.rmtree(copy_to)
    shutil.copytree(copy_from,copy_to)

  return

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def check_Analytical_block(data,identifier,nstates,eMsg):
  iline=-1
  while True:
    iline+=1
    if iline==len(data):
      if eMsg:
        print 'No matrix %s defined!' % (identifier)
      return False
    line=re.sub('#.*$','',data[iline]).split()
    if line==[]:
      continue
    ident=identifier.split()
    fits=True
    for i,el in enumerate(ident):
      if not el.lower() in line[i].lower():
        fits=False
        break
    if fits:
      break
  strings=data[iline+1:iline+1+nstates]
  for i,el in enumerate(strings):
    a=el.strip().split(',')
    if len(a)<i+1:
      if eMsg:
        print '%s matrix is not a lower triangular matrix with n=%i!' % (identifier,nstates)
      return False
  return True

# =================================================

def checktemplate_Analytical(filename,req_nstates,eMsg=True):
  f=open(filename)
  data=f.readlines()
  f.close()

  # check whether first two lines are positive integers
  try:
    natom=int(data[0])
    nstates=int(data[1])
  except ValueError:
    if eMsg:
      print 'First two lines must contain natom and nstates!'
    return False
  if natom<1 or nstates<1:
    if eMsg:
      print 'natom and nstates must be positive!'
    return False
  if nstates!=req_nstates:
    if eMsg:
      print 'Template file is for %i states!' % (nstates)
    return False

  # check the next natom lines
  variables=set()
  for i in range(2,2+natom):
    line=data[i]
    match=re.match('\s*[a-zA-Z]*\s+[a-zA-Z0][a-zA-Z0-9_]*\s+[a-zA-Z0][a-zA-Z0-9_]*\s+[a-zA-Z0][a-zA-Z0-9_]*',line)
    if not match:
      if eMsg:
        print 'Line %i malformatted!' % (i+1)
      return False
    else:
      a=line.split()
      for j in range(3):
        match=re.match('\s*[a-zA-Z][a-zA-Z0-9_]*',a[j+1])
        if match:
          variables.add(a[j+1])

  # check variable blocks
  iline=-1
  while True:
    iline+=1
    if iline==len(data):
      break
    line=re.sub('#.*$','',data[iline]).split()
    if line==[]:
      continue
    if 'variables' in line[0].lower():
      while True:
        iline+=1
        if iline==len(data):
          if eMsg:
            print 'Non-terminated variables block!'
          return False
        line=re.sub('#.*$','',data[iline]).split()
        if line==[]:
          continue
        if 'end' in line[0].lower():
          break
        match=re.match('[a-zA-Z][a-zA-Z0-9_]*',line[0])
        if not match:
          if eMsg:
            print 'Invalid variable name: %s' % (line[0])
          return False
        try:
          a=float(line[1])
        except ValueError:
          if eMsg:
            print 'Non-numeric value for variable %s' % (line[0])
          return False
        except IndexError:
          if eMsg:
            print 'No value for variable %s' % (line[0])
          return False

  # check hamiltonian block
  line='hamiltonian'
  a=check_Analytical_block(data,line,nstates,eMsg)
  if not a:
    return False

  # check derivatives of each variable
  for v in variables:
    line='derivatives %s' % (v)
    a=check_Analytical_block(data,line,nstates,eMsg)
    if not a:
      return False

  return True

# =================================================

def get_Analytical(INFOS):

  string='\n  '+'='*80+'\n'
  string+='||'+centerstring('Analytical PES Interface setup',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  print string

  if os.path.isfile('Analytical.template'):
    if checktemplate_Analytical('Analytical.template',INFOS['nstates'],eMsg=True):
      print 'Valid file "Analytical.template" detected. '
      usethisone=question('Use this template file?',bool,True)
      if usethisone:
        INFOS['analytical.template']='Analytical.template'
  if not 'analytical.template' in INFOS:
    while True:
      filename=question('Template filename:',str)
      if not os.path.isfile(filename):
        print 'File %s does not exist!' % (filename)
        continue
      if checktemplate_Analytical(filename,INFOS['nstates']):
        break
    INFOS['analytical.template']=filename
  print ''

  return INFOS

# =================================================

def prepare_Analytical(INFOS,iconddir):
  # copy Analytical.template

  # copy MOs and template
  cpfrom=INFOS['analytical.template']
  cpto='%s/Analytical.template' % (iconddir)
  shutil.copy(cpfrom,cpto)

  return

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def get_LVC(INFOS):

  string='\n  '+'='*80+'\n'
  string+='||'+centerstring('LVC Interface setup',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  print string

  if os.path.isfile('LVC.template'):
    print 'File "LVC.template" detected. '
    usethisone=question('Use this template file?',bool,True)
    if usethisone:
      INFOS['LVC.template']='LVC.template'
  if not 'LVC.template' in INFOS:
    while True:
      filename=question('Template filename:',str)
      if not os.path.isfile(filename):
        print 'File %s does not exist!' % (filename)
        continue

      break
    INFOS['LVC.template']=filename
  print ''

  return INFOS

# =================================================

def prepare_LVC(INFOS,iconddir):
  # copy LVC.template

  # copy MOs and template
  cpfrom=INFOS['LVC.template']
  cpto='%s/LVC.template' % (iconddir)
  shutil.copy(cpfrom,cpto)

  return

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def check_MOLCAS_qmmm(filename):
  f=open(filename)
  data=f.readlines()
  f.close()
  for line in data:
    if 'qmmm' in line.lower():
      return True
  return False

# =================================================

def checktemplate_MOLCAS(filename,INFOS):
  necessary=['basis','ras2','nactel','inactive']
  try:
    f=open(filename)
    data=f.readlines()
    f.close()
  except IOError:
    print 'Could not open template file %s' % (filename)
    return False
  valid=[]
  for i in necessary:
    for l in data:
      if i in re.sub('#.*$','',l):
        valid.append(True)
        break
    else:
      valid.append(False)
  if not all(valid):
    print 'The template %s seems to be incomplete! It should contain: ' % (filename) +str(necessary)
    return False
  roots_there=False
  for l in data:
    l=re.sub('#.*$','',l).lower().split()
    if len(l)==0:
      continue
    if 'roots' in l[0]:
      roots_there=True
  if not roots_there:
    for mult,state in enumerate(INFOS['states']):
      if state<=0:
        continue
      valid=[]
      for l in data:
        if 'spin' in re.sub('#.*$','',l).lower():
          f=l.split()
          if int(f[1])==mult+1:
            valid.append(True)
            break
      else:
        valid.append(False)
  if not all(valid):
    string='The template %s seems to be incomplete! It should contain the keyword "spin" for ' % (filename)
    for mult,state in enumerate(INFOS['states']):
      if state<=0:
        continue
      string+='%s, ' % (IToMult[mult+1])
    string=string[:-2]+'!'
    print string
    return False
  return True

# =================================================

def get_MOLCAS(INFOS):
  '''This routine asks for all questions specific to MOLPRO:
  - path to molpro
  - scratch directory
  - MOLPRO.template
  - wf.init
  '''

  string='\n  '+'='*80+'\n'
  string+='||'+centerstring('MOLCAS Interface setup',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  print string

  print centerstring('Path to MOLCAS',60,'-')+'\n'
  path=os.getenv('MOLCAS')
  #path=os.path.expanduser(os.path.expandvars(path))
  if path=='':
    path=None
  else:
    path='$MOLCAS/'
      #print 'Environment variable $MOLCAS detected:\n$MOLCAS=%s\n' % (path)
      #if question('Do you want to use this MOLCAS installation?',bool,True):
        #INFOS['molcas']=path
    #if not 'molcas' in INFOS:
  print '\nPlease specify path to MOLCAS directory (SHELL variables and ~ can be used, will be expanded when interface is started).\n'
  INFOS['molcas']=question('Path to MOLCAS:',str,path)
  print ''


  print centerstring('Scratch directory',60,'-')+'\n'
  print 'Please specify an appropriate scratch directory. This will be used to temporally store the integrals. The scratch directory will be deleted after the calculation. Remember that this script cannot check whether the path is valid, since you may run the calculations on a different machine. The path will not be expanded by this script.'
  INFOS['scratchdir']=question('Path to scratch directory:',str)
  print ''


  print centerstring('MOLCAS input template file',60,'-')+'\n'
  print '''Please specify the path to the MOLCAS.template file. This file must contain the following settings:

basis <Basis set>
ras2 <Number of active orbitals>
nactel <Number of active electrons>
inactive <Number of doubly occupied orbitals>
roots <Number of roots for state-averaging>

The MOLCAS interface will generate the appropriate MOLCAS input automatically.
'''
  if os.path.isfile('MOLCAS.template'):
    if checktemplate_MOLCAS('MOLCAS.template',INFOS):
      print 'Valid file "MOLCAS.template" detected. '
      usethisone=question('Use this template file?',bool,True)
      if usethisone:
        INFOS['molcas.template']='MOLCAS.template'
  if not 'molcas.template' in INFOS:
    while True:
      filename=question('Template filename:',str)
      if not os.path.isfile(filename):
        print 'File %s does not exist!' % (filename)
        continue
      if checktemplate_MOLCAS(filename,INFOS):
        break
    INFOS['molcas.template']=filename
  print ''


  # QMMM 
  if check_MOLCAS_qmmm(INFOS['molcas.template']):
    print centerstring('MOLCAS+TINKER QM/MM setup',60,'-')+'\n'
    print 'Your template specifies a QM/MM calculation. Please specify the path to TINKER.' 
    path=os.getenv('TINKER')
    if path=='':
      path=None
    else:
      path='$TINKER/'
    print '\nPlease specify path to TINKER bin/ directory (SHELL variables and ~ can be used, will be expanded when interface is started).\n'
    INFOS['tinker']=question('Path to TINKER/bin:',str,path)
    print 'Please give the key and connection table files.'
    while True:
      filename=question('Key file:',str)
      if not os.path.isfile(filename):
        print 'File %s does not exist!' % (filename)
      else:
        break
    INFOS['MOLCAS.fffile']=filename
    while True:
      filename=question('Connection table file:',str)
      if not os.path.isfile(filename):
        print 'File %s does not exist!' % (filename)
      else:
        break
    INFOS['MOLCAS.ctfile']=filename


  print centerstring('Initial wavefunction: MO Guess',60,'-')+'\n'
  print '''Please specify the path to a MOLCAS JobIph file containing suitable starting MOs for the CASSCF calculation. Please note that this script cannot check whether the wavefunction file and the Input template are consistent!
'''
  string='Do you have initial wavefunction files for '
  for mult,state in enumerate(INFOS['states']):
    if state<=0:
      continue
    string+='%s, ' % (IToMult[mult+1])
  string=string[:-2]+'?'
  if question(string,bool,True):
    while True:
      jobiph_or_rasorb=question('JobIph files (1) or RasOrb files (2)?',int)[0]
      if jobiph_or_rasorb in [1,2]:
        break
    INFOS['molcas.jobiph_or_rasorb']=jobiph_or_rasorb
    INFOS['molcas.guess']={}
    for mult,state in enumerate(INFOS['states']):
      if state<=0:
        continue
      while True:
        if jobiph_or_rasorb==1:
          guess_file='MOLCAS.%i.JobIph.init' % (mult+1)
        else:
          guess_file='MOLCAS.%i.RasOrb.init' % (mult+1)
        filename=question('Initial wavefunction file for %ss:' % (IToMult[mult+1]),str,guess_file)
        if os.path.isfile(filename):
          INFOS['molcas.guess'][mult+1]=filename
          break
        else:
          print 'File not found!'
  else:
    print 'WARNING: Remember that CASSCF calculations may run very long and/or yield wrong results without proper starting MOs.'
    time.sleep(2)
    INFOS['molcas.guess']={}


  print centerstring('MOLCAS Ressource usage',60,'-')+'\n'
  print '''Please specify the amount of memory available to MOLCAS (in MB). For calculations including moderately-sized CASSCF calculations and less than 150 basis functions, around 2000 MB should be sufficient.
'''
  INFOS['molcas.mem']=abs(question('MOLCAS memory:',int,[1000])[0])
  print '''Please specify the number of CPUs to be used by EACH calculation.
'''
  INFOS['molcas.ncpu']=abs(question('Number of CPUs:',int,[1])[0])



  # Ionization
  #print '\n'+centerstring('Ionization probability by Dyson norms',60,'-')+'\n'
  #INFOS['ion']=question('Dyson norms?',bool,False)
  #if INFOS['ion']:
  if 'wfoverlap' in INFOS['needed']:
    print '\n'+centerstring('WFoverlap setup',60,'-')+'\n'
    INFOS['molcas.wfoverlap']=question('Path to wavefunction overlap executable:',str,'$SHARC/wfoverlap.x')
    # TODO not asked for: numfrozcore, numocc

  return INFOS

# =================================================

def prepare_MOLCAS(INFOS,iconddir):
  # write MOLCAS.resources
  try:
    sh2cas=open('%s/MOLCAS.resources' % (iconddir), 'w')
  except IOError:
    print 'IOError during prepareMOLCAS, iconddir=%s' % (iconddir)
    quit(1)
  project='MOLCAS'
  string='molcas %s\nscratchdir %s/%s/\nmemory %i\nncpu %i\nproject %s' % (INFOS['molcas'],INFOS['scratchdir'],iconddir,INFOS['molcas.mem'],INFOS['molcas.ncpu'],project)
  if 'wfoverlap' in INFOS['needed']:
    string+='\nwfoverlap %s\n' % INFOS['molcas.wfoverlap']
  else:
    string+='\nnooverlap\n'
  if 'tinker' in INFOS:
    string+='tinker %s' % (INFOS['tinker'])
  sh2cas.write(string)
  sh2cas.close()

  # copy MOs and template
  cpfrom=INFOS['molcas.template']
  cpto='%s/MOLCAS.template' % (iconddir)
  shutil.copy(cpfrom,cpto)
  if not INFOS['molcas.guess']=={}:
    for i in INFOS['molcas.guess']:
      if INFOS['molcas.jobiph_or_rasorb']==1:
        cpfrom=INFOS['molcas.guess'][i]
        cpto='%s/%s.%i.JobIph.init' % (iconddir,project,i)
      else:
        cpfrom=INFOS['molcas.guess'][i]
        cpto='%s/%s.%i.RasOrb.init' % (iconddir,project,i)
      shutil.copy(cpfrom,cpto)

  if 'MOLCAS.fffile' in INFOS:
    cpfrom1=INFOS['MOLCAS.fffile']
    cpto1='%s/MOLCAS.qmmm.key' % (iconddir)
    shutil.copy(cpfrom1,cpto1)

  if 'MOLCAS.ctfile' in INFOS:
    cpfrom1=INFOS['MOLCAS.ctfile']
    cpto1='%s/MOLCAS.qmmm.table' % (iconddir)
    shutil.copy(cpfrom1,cpto1)
  return

#======================================================================================================================
#======================================================================================================================
#======================================================================================================================

def checktemplate_ADF(filename,INFOS):
  necessary=['basis','functional','charge']
  try:
    f=open(filename)
    data=f.readlines()
    f.close()
  except IOError:
    print 'Could not open template file %s' % (filename)
    return False
  valid=[]
  for i in necessary:
    for l in data:
      line=l.lower().split()
      if len(line)==0:
        continue
      line=line[0]
      if i==re.sub('#.*$','',line):
        valid.append(True)
        break
    else:
      valid.append(False)
  if not all(valid):
    print 'The template %s seems to be incomplete! It should contain: ' % (filename) +str(necessary)
    return False
  return True

# =================================================

def qmmm_job(filename,INFOS):
  necessary=['qmmm']
  try:
    f=open(filename)
    data=f.readlines()
    f.close()
  except IOError:
    print 'Could not open template file %s' % (filename)
    return False
  valid=[]
  for i in necessary:
    for l in data:
      line=l.lower().split()
      if len(line)==0:
        continue
      line=line[0]
      if i==re.sub('#.*$','',line):
        valid.append(True)
        break
    else:
      valid.append(False)
  if not all(valid):
    return False
  return True

# =================================================

def get_ADF(INFOS):
  '''This routine asks for all questions specific to ADF:
  - path to ADF
  - scratch directory
  - ADF.template
  - TAPE21
  '''

  string='\n  '+'='*80+'\n'
  string+='||'+centerstring('ADF Interface setup',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  print string

  print centerstring('Path to ADF',60,'-')+'\n'
  path=os.getenv('ADFHOME')
  if path:
    path='$ADFHOME/'
  adfrc=question('Setup from adfrc.sh file?',bool,True)
  if adfrc:
    if path:
      path='$ADFHOME/adfrc.sh'
    print '\nPlease specify path to the adfrc.sh file (SHELL variables and ~ can be used, will be expanded when interface is started).\n'
    path=question('Path to adfrc.sh file:',str,path)
    INFOS['adfrc']=os.path.abspath(os.path.expanduser(os.path.expandvars(path)))
    print 'Will use adfrc= %s' % INFOS['adfrc']
    INFOS['adf']='$ADFHOME'
    INFOS['scmlicense']='$SCMLICENSE'
    print ''
  else:
    print '\nPlease specify path to ADF directory (SHELL variables and ~ can be used, will be expanded when interface is started).\n'
    INFOS['adf']=question('Path to ADF:',str,path)
    print ''
    print centerstring('Path to ADF license file',60,'-')+'\n'
    path=os.getenv('SCMLICENSE')
    #path=os.path.expanduser(os.path.expandvars(path))
    if path=='':
      path=None
    else:
      path='$SCMLICENSE'
    print'\nPlease specify path to ADF license.txt\n'
    INFOS['scmlicense']=question('Path to license:',str,path)
    print ''





  # scratch
  print centerstring('Scratch directory',60,'-')+'\n'
  print 'Please specify an appropriate scratch directory. This will be used to run the ADF calculations. The scratch directory will be deleted after the calculation. Remember that this script cannot check whether the path is valid, since you may run the calculations on a different machine. The path will not be expanded by this script.'
  INFOS['scratchdir']=question('Path to scratch directory:',str)
  print ''


  # template file
  print centerstring('ADF input template file',60,'-')+'\n'
  print '''Please specify the path to the ADF.template file. This file must contain the following keywords:

basis <basis>
functional <type> <name>
charge <x> [ <x2> [ <x3> ...] ]

The ADF interface will generate the appropriate ADF input automatically.
'''
  if os.path.isfile('ADF.template'):
    if checktemplate_ADF('ADF.template',INFOS):
      print 'Valid file "ADF.template" detected. '
      usethisone=question('Use this template file?',bool,True)
      if usethisone:
        INFOS['ADF.template']='ADF.template'
  if not 'ADF.template' in INFOS:
    while True:
      filename=question('Template filename:',str)
      if not os.path.isfile(filename):
        print 'File %s does not exist!' % (filename)
        continue
      if checktemplate_ADF(filename,INFOS):
        break
    INFOS['ADF.template']=filename
  print ''


  # QMMM
  if qmmm_job(INFOS['ADF.template'],INFOS):
    print centerstring('ADF QM/MM setup',60,'-')+'\n'
    print 'Your template specifies a QM/MM calculation. Please give the force field and connection table files.'
    while True:
      filename=question('Force field file:',str)
      if not os.path.isfile(filename):
        print 'File %s does not exist!' % (filename)
      else:
        break
    INFOS['ADF.fffile']=filename
    while True:
      filename=question('Connection table file:',str)
      if not os.path.isfile(filename):
        print 'File %s does not exist!' % (filename)
      else:
        break
    INFOS['ADF.ctfile']=filename


  # initial MOs
  print centerstring('Initial restart: MO Guess',60,'-')+'\n'
  print '''Please specify the path to an ADF TAPE21 file containing suitable starting MOs for the ADF calculation. Please note that this script cannot check whether the wavefunction file and the Input template are consistent!
'''
  if question('Do you have a restart file?',bool,True):
    if True:
      while True:
        filename=question('Restart file:',str,'ADF.t21.init')
        if os.path.isfile(filename):
          INFOS['adf.guess']=filename
          break
        else:
          print 'Could not find file "%s"!' % (filename)
  else:
    INFOS['adf.guess']={}


  # Resources
  print centerstring('ADF Ressource usage',60,'-')+'\n'
  print '''Please specify the number of CPUs to be used by EACH calculation.
'''
  INFOS['adf.ncpu']=abs(question('Number of CPUs:',int)[0])

  if INFOS['adf.ncpu']>1:
    print '''Please specify how well your job will parallelize.
A value of 0 means that running in parallel will not make the calculation faster, a value of 1 means that the speedup scales perfectly with the number of cores.
Typical values for ADF are 0.90-0.98 for LDA/GGA functionals and 0.50-0.80 for hybrids (better if RIHartreeFock is used).'''
    INFOS['adf.scaling']=min(1.0,max(0.0,question('Parallel scaling:',float,[0.8])[0] ))
  else:
    INFOS['adf.scaling']=0.9


  # Ionization
  #print '\n'+centerstring('Ionization probability by Dyson norms',60,'-')+'\n'
  #INFOS['ion']=question('Dyson norms?',bool,False)
  #if INFOS['ion']:
  if 'wfoverlap' in INFOS['needed']:
    print '\n'+centerstring('WFoverlap setup',60,'-')+'\n'
    INFOS['adf.wfoverlap']=question('Path to wavefunction overlap executable:',str,'$SHARC/wfoverlap.x')
    print ''
    print 'State threshold for choosing determinants to include in the overlaps'
    print 'For hybrids (and without TDA) one should consider that the eigenvector X may have a norm larger than 1'
    INFOS['adf.ciothres']=question('Threshold:',float,[0.99])[0]
    print ''
    INFOS['adf.mem']=question('Memory for wfoverlap (MB):',int,[1000])[0]
    # TODO not asked: numfrozcore and numocc

    #print 'Please state the number of core orbitals you wish to freeze for the overlaps (recommended to use for at least the 1s orbital and a negative number uses default values)?'
    #print 'A value of -1 will use the defaults used by ADF for a small frozen core and 0 will turn off the use of frozen cores'
    #INFOS['frozcore_number']=question('How many orbital to freeze?',int,[-1])[0]


  # TheoDORE
  theodore_spelling=['Om',
                    'PRNTO',
                    'Z_HE', 'S_HE', 'RMSeh',
                    'POSi', 'POSf', 'POS',
                    'PRi', 'PRf', 'PR', 'PRh',
                    'CT', 'CT2', 'CTnt',
                    'MC', 'LC', 'MLCT', 'LMCT', 'LLCT',
                    'DEL', 'COH', 'COHh']
  #INFOS['theodore']=question('TheoDORE analysis?',bool,False)
  if 'theodore' in INFOS['needed']:
    print '\n'+centerstring('Wave function analysis by TheoDORE',60,'-')+'\n'

    INFOS['adf.theodore']=question('Path to TheoDORE directory:',str,'$THEODIR')
    print ''

    print 'Please give a list of the properties to calculate by TheoDORE.\nPossible properties:'
    string=''
    for i,p in enumerate(theodore_spelling):
      string+='%s ' % (p)
      if (i+1)%8==0:
        string+='\n'
    print string
    l=question('TheoDORE properties:',str,'Om  PRNTO  S_HE  Z_HE  RMSeh')
    if '[' in l:
      INFOS['theodore.prop']=ast.literal_eval(l)
    else:
      INFOS['theodore.prop']=l.split()
    print ''

    print 'Please give a list of the fragments used for TheoDORE analysis.'
    print 'You can use the list-of-lists from dens_ana.in'
    print 'Alternatively, enter all atom numbers for one fragment in one line. After defining all fragments, type "end".'
    if qmmm_job(INFOS['ADF.template'],INFOS):
      print 'You should only include the atom numbers of QM and link atoms.'
    INFOS['theodore.frag']=[]
    while True:
      l=question('TheoDORE fragment:',str,'end')
      if 'end' in l.lower():
        break
      if '[' in l:
        try:
          INFOS['theodore.frag']=ast.literal_eval(l)
          break
        except ValueError:
          continue
      f=[ int(i) for i in l.split() ]
      INFOS['theodore.frag'].append(f)
    INFOS['theodore.count']=len(INFOS['theodore.prop'])+len(INFOS['theodore.frag'])**2


  return INFOS

# =================================================

def prepare_ADF(INFOS,iconddir):
  # write ADF.resources
  try:
    sh2cas=open('%s/ADF.resources' % (iconddir), 'w')
  except IOError:
    print 'IOError during prepareADF, iconddir=%s' % (iconddir)
    quit(1)
#  project='ADF'
  string='adfhome %s\nscmlicense %s\nscratchdir %s/%s/\nncpu %i\nschedule_scaling %f\n' % (INFOS['adf'],INFOS['scmlicense'],INFOS['scratchdir'],iconddir,INFOS['adf.ncpu'],INFOS['adf.scaling'])
  if 'wfoverlap' in INFOS['needed']:
    string+='wfoverlap %s\nwfthres %f\n' % (INFOS['adf.wfoverlap'],INFOS['adf.ciothres'])
    string+='memory %i\n' % (INFOS['adf.mem'])
    #string+='numfrozcore %i\n' %(INFOS['frozcore_number'])
  else:
    string+='nooverlap\n'
  if 'theodore' in INFOS['needed']:
    string+='theodir %s\n' % (INFOS['adf.theodore'])
    string+='theodore_prop %s\n' % (INFOS['theodore.prop'])
    string+='theodore_fragment %s\n' % (INFOS['theodore.frag'])
  if 'ADF.fffile' in INFOS:
    string+='qmmm_ff_file ADF.qmmm.ff\n'
  if 'ADF.ctfile' in INFOS:
    string+='qmmm_table ADF.qmmm.table\n'
  sh2cas.write(string)
  sh2cas.close()

  # copy MOs and template
  cpfrom=INFOS['ADF.template']
  cpto='%s/ADF.template' % (iconddir)
  shutil.copy(cpfrom,cpto)

  if INFOS['adf.guess']:
    cpfrom1=INFOS['adf.guess']
    cpto1='%s/ADF.t21_init' % (iconddir)
    shutil.copy(cpfrom1,cpto1)

  if 'ADF.fffile' in INFOS:
    cpfrom1=INFOS['ADF.fffile']
    cpto1='%s/ADF.qmmm.ff' % (iconddir)
    shutil.copy(cpfrom1,cpto1)

  if 'ADF.ctfile' in INFOS:
    cpfrom1=INFOS['ADF.ctfile']
    cpto1='%s/ADF.qmmm.table' % (iconddir)
    shutil.copy(cpfrom1,cpto1)

  return

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def checktemplate_RICC2(filename,INFOS):
  necessary=['basis']
  try:
    f=open(filename)
    data=f.readlines()
    f.close()
  except IOError:
    print 'Could not open template file %s' % (filename)
    return False
  valid=[]
  for i in necessary:
    for l in data:
      line=l.lower()
      if i in re.sub('#.*$','',line):
        valid.append(True)
        break
    else:
      valid.append(False)
  if not all(valid):
    print 'The template %s seems to be incomplete! It should contain: ' % (filename) +str(necessary)
    return False
  return True

# =================================================

def get_RICC2(INFOS):
  string='\n  '+'='*80+'\n'
  string+='||'+centerstring('Turbomole RICC2 Interface setup',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  print string

  print centerstring('Path to TURBOMOLE',60,'-')+'\n'
  path=os.getenv('TURBODIR')
  if path=='':
    path=None
  else:
    path='$TURBODIR/'
  print '\nPlease specify path to TURBOMOLE directory (SHELL variables and ~ can be used, will be expanded when interface is started).\n'
  INFOS['turbomole']=question('Path to TURBOMOLE:',str,path)
  print ''

  if INFOS['soc']:
    print centerstring('Path to ORCA',60,'-')+'\n'
    path=os.getenv('ORCADIR')
    if path=='':
      path=None
    else:
      path='$ORCADIR/'
    print '\nPlease specify path to ORCA directory (SHELL variables and ~ can be used, will be expanded when interface is started).\n\nORCA is necessary for the calculation of spin-orbit couplings with ricc2.\n'
    INFOS['orca']=question('Path to ORCA:',str,path)
    print ''


  print centerstring('Scratch directory',60,'-')+'\n'
  print 'Please specify an appropriate scratch directory. This will be used to temporally store the integrals. The scratch directory will be deleted after the calculation. Remember that this script cannot check whether the path is valid, since you may run the calculations on a different machine. The path will not be expanded by this script.'
  INFOS['scratchdir']=question('Path to scratch directory:',str)
  print ''


  print centerstring('RICC2 input template file',60,'-')+'\n'
  print '''Please specify the path to the RICC2.template file. This file must contain the following settings:

basis <Basis set>

In addition, it can contain the following:

auxbasis <Basis set>
charge <integer>
method <"ADC(2)" or "CC2">                      # only ADC(2) can calculate spin-orbit couplings
frozen <number of frozen core orbitals>
spin-scaling <"none", "SCS", or "SOS">
douglas-kroll                                   # DKH is only used if this keyword is given

'''
  if os.path.isfile('RICC2.template'):
    if checktemplate_RICC2('RICC2.template',INFOS):
      print 'Valid file "RICC2.template" detected. '
      usethisone=question('Use this template file?',bool,True)
      if usethisone:
        INFOS['ricc2.template']='RICC2.template'
  if not 'ricc2.template' in INFOS:
    while True:
      filename=question('Template filename:',str)
      if not os.path.isfile(filename):
        print 'File %s does not exist!' % (filename)
        continue
      if checktemplate_RICC2(filename,INFOS):
        break
    INFOS['ricc2.template']=filename
  print ''


  print centerstring('Initial wavefunction: MO Guess',60,'-')+'\n'
  print '''Please specify the path to a Turbomole "mos" file containing suitable starting MOs for the calculation. Please note that this script cannot check whether the file and the input template are consistent!
'''
  string='Do you have an initial orbitals file?'
  if question(string,bool,True):
    while True:
      guess_file='mos'
      filename=question('Initial wavefunction file:',str,guess_file)
      if os.path.isfile(filename):
        INFOS['ricc2.guess']=filename
        break
      else:
        print 'File not found!'
  else:
    INFOS['ricc2.guess']=[]


  print centerstring('RICC2 Ressource usage',60,'-')+'\n'
  print '''Please specify the amount of memory available to Turbomole (in MB).
'''
  INFOS['ricc2.mem']=abs(question('RICC2 memory:',int,[1000])[0])
  print '''Please specify the number of CPUs to be used by EACH trajectory.
'''
  INFOS['ricc2.ncpu']=abs(question('Number of CPUs:',int,[1])[0])



  if 'wfoverlap' in INFOS['needed']:
    print '\n'+centerstring('WFoverlap setup',60,'-')+'\n'
    INFOS['ricc2.wfoverlap']=question('Path to wavefunction overlap executable:',str,'$SHARC/wfoverlap.x')
    print ''
    print 'State threshold for choosing determinants to include in the overlaps'
    #print 'For hybrids (and without TDA) one should consider that the eigenvector X may have a norm larger than 1'
    INFOS['ricc2.ciothres']=question('Threshold:',float,[0.99])[0]



  # TheoDORE
  theodore_spelling=['Om',
                    'PRNTO',
                    'Z_HE', 'S_HE', 'RMSeh',
                    'POSi', 'POSf', 'POS',
                    'PRi', 'PRf', 'PR', 'PRh',
                    'CT', 'CT2', 'CTnt',
                    'MC', 'LC', 'MLCT', 'LMCT', 'LLCT',
                    'DEL', 'COH', 'COHh']
  #INFOS['theodore']=question('TheoDORE analysis?',bool,False)
  if 'theodore' in INFOS['needed']:
    print '\n'+centerstring('Wave function analysis by TheoDORE',60,'-')+'\n'

    INFOS['ricc2.theodore']=question('Path to TheoDORE directory:',str,'$THEODIR')
    print ''

    print 'Please give a list of the properties to calculate by TheoDORE.\nPossible properties:'
    string=''
    for i,p in enumerate(theodore_spelling):
      string+='%s ' % (p)
      if (i+1)%8==0:
        string+='\n'
    print string
    l=question('TheoDORE properties:',str,'Om  PRNTO  S_HE  Z_HE  RMSeh')
    if '[' in l:
      INFOS['theodore.prop']=ast.literal_eval(l)
    else:
      INFOS['theodore.prop']=l.split()
    print ''

    print 'Please give a list of the fragments used for TheoDORE analysis.'
    print 'You can use the list-of-lists from dens_ana.in'
    print 'Alternatively, enter all atom numbers for one fragment in one line. After defining all fragments, type "end".'
    INFOS['theodore.frag']=[]
    while True:
      l=question('TheoDORE fragment:',str,'end')
      if 'end' in l.lower():
        break
      if '[' in l:
        try:
          INFOS['theodore.frag']=ast.literal_eval(l)
          break
        except ValueError:
          continue
      f=[ int(i) for i in l.split() ]
      INFOS['theodore.frag'].append(f)
    INFOS['theodore.count']=len(INFOS['theodore.prop'])+len(INFOS['theodore.frag'])**2




  return INFOS

# =================================================

def prepare_RICC2(INFOS,iconddir):
  # write RICC2.resources
  try:
    sh2cas=open('%s/RICC2.resources' % (iconddir), 'w')
  except IOError:
    print 'IOError during prepare_RICC2, iconddir=%s' % (iconddir)
    quit(1)
  string='''turbodir %s
scratchdir %s/%s
memory %i
ncpu %i
dipolelevel 1
''' % (INFOS['turbomole'],
       INFOS['scratchdir'],
       iconddir,
       INFOS['ricc2.mem'],
       INFOS['ricc2.ncpu'])
  if INFOS['soc']:
    string+='orcadir %s\n' % (INFOS['orca'])
  if 'wfoverlap' in INFOS['needed']:
    string+='wfoverlap %s\nwfthres %f\n' % (INFOS['ricc2.wfoverlap'],INFOS['ricc2.ciothres'])
    #string+='numfrozcore %i\n' %(INFOS['frozcore_number'])
  else:
    string+='nooverlap\n'
  if 'theodore' in INFOS['needed']:
    string+='theodir %s\n' % (INFOS['ricc2.theodore'])
    string+='theodore_prop %s\n' % (INFOS['theodore.prop'])
    string+='theodore_fragment %s\n' % (INFOS['theodore.frag'])
  sh2cas.write(string)
  sh2cas.close()

  # copy MOs and template
  cpfrom=INFOS['ricc2.template']
  cpto='%s/RICC2.template' % (iconddir)
  shutil.copy(cpfrom,cpto)
  if INFOS['ricc2.guess']:
    cpfrom1=INFOS['ricc2.guess']
    cpto1='%s/mos.init' % (iconddir)
    shutil.copy(cpfrom1,cpto1)
  return

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def checktemplate_GAUSSIAN(filename,INFOS):
  necessary=['basis','functional','charge']
  try:
    f=open(filename)
    data=f.readlines()
    f.close()
  except IOError:
    print 'Could not open template file %s' % (filename)
    return False
  valid=[]
  for i in necessary:
    for l in data:
      line=l.lower().split()
      if len(line)==0:
        continue
      line=line[0]
      if i==re.sub('#.*$','',line):
        valid.append(True)
        break
    else:
      valid.append(False)
  if not all(valid):
    print 'The template %s seems to be incomplete! It should contain: ' % (filename) +str(necessary)
    return False
  return True

# =================================================

def get_GAUSSIAN(INFOS):
  '''This routine asks for all questions specific to GAUSSIAN:
  - path to GAUSSIAN
  - scratch directory
  - GAUSSIAN.template
  - TAPE21
  '''

  string='\n  '+'='*80+'\n'
  string+='||'+centerstring('GAUSSIAN Interface setup',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  print string

  print centerstring('Path to GAUSSIAN',60,'-')+'\n'
  tries=['g16root','g09root','g03root']
  for i in tries:
    path=os.getenv(i)
    if path:
      path='$%s/' % i
      break
  #gaussianprofile=question('Setup from gaussian.profile file?',bool,True)
  #if gaussianprofile:
    #if path:
      #path='%s/gaussian.profile' % path
    #print '\nPlease specify path to the gaussian.profile file (SHELL variables and ~ can be used, will be expanded when interface is started).\n'
    #path=question('Path to GAUSSIAN:',str,path)
    #INFOS['gaussianprofile']=os.path.abspath(os.path.expanduser(os.path.expandvars(path)))
    #print 'Will use gaussianprofile= %s' % INFOS['gaussianprofile']
    #INFOS['gaussian']='$GAUSSIANHOME'
    #print ''
  #else:
  print '\nPlease specify path to GAUSSIAN directory (SHELL variables and ~ can be used, will be expanded when interface is started).\n'
  INFOS['groot']=question('Path to GAUSSIAN:',str,path)
  print ''


  # scratch
  print centerstring('Scratch directory',60,'-')+'\n'
  print 'Please specify an appropriate scratch directory. This will be used to run the GAUSSIAN calculations. The scratch directory will be deleted after the calculation. Remember that this script cannot check whether the path is valid, since you may run the calculations on a different machine. The path will not be expanded by this script.'
  INFOS['scratchdir']=question('Path to scratch directory:',str)
  print ''


  # template file
  print centerstring('GAUSSIAN input template file',60,'-')+'\n'
  print '''Please specify the path to the GAUSSIAN.template file. This file must contain the following keywords:
  
basis <basis>
functional <type> <name>
charge <x> [ <x2> [ <x3> ...] ] 

The GAUSSIAN interface will generate the appropriate GAUSSIAN input automatically.
'''
  if os.path.isfile('GAUSSIAN.template'):
    if checktemplate_GAUSSIAN('GAUSSIAN.template',INFOS):
      print 'Valid file "GAUSSIAN.template" detected. '
      usethisone=question('Use this template file?',bool,True)
      if usethisone:
        INFOS['GAUSSIAN.template']='GAUSSIAN.template'
  if not 'GAUSSIAN.template' in INFOS:
    while True:
      filename=question('Template filename:',str)
      if not os.path.isfile(filename):
        print 'File %s does not exist!' % (filename)
        continue
      if checktemplate_GAUSSIAN(filename,INFOS):
        break
    INFOS['GAUSSIAN.template']=filename
  print ''



  # initial MOs
  print centerstring('Initial restart: MO Guess',60,'-')+'\n'
  print '''Please specify the path to an GAUSSIAN chk file containing suitable starting MOs for the GAUSSIAN calculation. Please note that this script cannot check whether the wavefunction file and the Input template are consistent!
'''
  if question('Do you have a restart file?',bool,True):
    if True:
      while True:
        filename=question('Restart file:',str,'GAUSSIAN.chk.init')
        if os.path.isfile(filename):
          INFOS['gaussian.guess']=filename
          break
        else:
          print 'Could not find file "%s"!' % (filename)
  else:
    INFOS['gaussian.guess']={}


  # Resources
  print centerstring('GAUSSIAN Ressource usage',60,'-')+'\n'
  print '''Please specify the number of CPUs to be used by EACH calculation.
'''
  INFOS['gaussian.ncpu']=abs(question('Number of CPUs:',int)[0])

  if INFOS['gaussian.ncpu']>1:
    print '''Please specify how well your job will parallelize.
A value of 0 means that running in parallel will not make the calculation faster, a value of 1 means that the speedup scales perfectly with the number of cores.
Typical values for GAUSSIAN are 0.90-0.98.'''
    INFOS['gaussian.scaling']=min(1.0,max(0.0,question('Parallel scaling:',float,[0.9])[0] ))
  else:
    INFOS['gaussian.scaling']=0.9

  INFOS['gaussian.mem']=question('Memory (MB):',int,[1000])[0]

  # Ionization
  #print '\n'+centerstring('Ionization probability by Dyson norms',60,'-')+'\n'
  #INFOS['ion']=question('Dyson norms?',bool,False)
  #if INFOS['ion']:
  if 'wfoverlap' in INFOS['needed']:
    print '\n'+centerstring('WFoverlap setup',60,'-')+'\n'
    INFOS['gaussian.wfoverlap']=question('Path to wavefunction overlap executable:',str,'$SHARC/wfoverlap.x')
    print ''
    print 'State threshold for choosing determinants to include in the overlaps'
    print 'For hybrids without TDA one should consider that the eigenvector X may have a norm larger than 1'
    INFOS['gaussian.ciothres']=question('Threshold:',float,[0.99])[0]
    print ''
    # TODO not asked: numfrozcore and numocc

    #print 'Please state the number of core orbitals you wish to freeze for the overlaps (recommended to use for at least the 1s orbital and a negative number uses default values)?'
    #print 'A value of -1 will use the defaults used by GAUSSIAN for a small frozen core and 0 will turn off the use of frozen cores'
    #INFOS['frozcore_number']=question('How many orbital to freeze?',int,[-1])[0]


  # TheoDORE
  theodore_spelling=['Om', 
                    'PRNTO', 
                    'Z_HE', 'S_HE', 'RMSeh',
                    'POSi', 'POSf', 'POS', 
                    'PRi', 'PRf', 'PR', 'PRh',
                    'CT', 'CT2', 'CTnt',
                    'MC', 'LC', 'MLCT', 'LMCT', 'LLCT', 
                    'DEL', 'COH', 'COHh']
  #INFOS['theodore']=question('TheoDORE analysis?',bool,False)
  if 'theodore' in INFOS['needed']:
    print '\n'+centerstring('Wave function analysis by TheoDORE',60,'-')+'\n'

    INFOS['gaussian.theodore']=question('Path to TheoDORE directory:',str,'$THEODIR')
    print ''

    print 'Please give a list of the properties to calculate by TheoDORE.\nPossible properties:'
    string=''
    for i,p in enumerate(theodore_spelling):
      string+='%s ' % (p)
      if (i+1)%8==0:
        string+='\n'
    print string
    l=question('TheoDORE properties:',str,'Om  PRNTO  S_HE  Z_HE  RMSeh')
    if '[' in l:
      INFOS['theodore.prop']=ast.literal_eval(l)
    else:
      INFOS['theodore.prop']=l.split()
    print ''

    print 'Please give a list of the fragments used for TheoDORE analysis.'
    print 'You can use the list-of-lists from dens_ana.in'
    print 'Alternatively, enter all atom numbers for one fragment in one line. After defining all fragments, type "end".'
    if qmmm_job(INFOS['GAUSSIAN.template'],INFOS):
      print 'You should only include the atom numbers of QM and link atoms.'
    INFOS['theodore.frag']=[]
    while True:
      l=question('TheoDORE fragment:',str,'end')
      if 'end' in l.lower():
        break
      if '[' in l:
        try:
          INFOS['theodore.frag']=ast.literal_eval(l)
          break
        except ValueError:
          continue
      f=[ int(i) for i in l.split() ]
      INFOS['theodore.frag'].append(f)
    INFOS['theodore.count']=len(INFOS['theodore.prop'])+len(INFOS['theodore.frag'])**2


  return INFOS

# =================================================

def prepare_GAUSSIAN(INFOS,iconddir):
  # write GAUSSIAN.resources
  try:
    sh2cas=open('%s/GAUSSIAN.resources' % (iconddir), 'w')
  except IOError:
    print 'IOError during prepareGAUSSIAN, iconddir=%s' % (iconddir)
    quit(1)
#  project='GAUSSIAN'
  string='groot %s\nscratchdir %s/%s/\nncpu %i\nschedule_scaling %f\n' % (INFOS['groot'],INFOS['scratchdir'],iconddir,INFOS['gaussian.ncpu'],INFOS['gaussian.scaling'])
  string+='memory %i\n' % (INFOS['gaussian.mem'])
  if 'wfoverlap' in INFOS['needed']:
    string+='wfoverlap %s\nwfthres %f\n' % (INFOS['gaussian.wfoverlap'],INFOS['gaussian.ciothres'])
    #string+='numfrozcore %i\n' %(INFOS['frozcore_number'])
  else:
    string+='nooverlap\n'
  if 'theodore' in INFOS['needed']:
    string+='theodir %s\n' % (INFOS['gaussian.theodore'])
    string+='theodore_prop %s\n' % (INFOS['theodore.prop'])
    string+='theodore_fragment %s\n' % (INFOS['theodore.frag'])
  sh2cas.write(string)
  sh2cas.close()

  # copy MOs and template
  cpfrom=INFOS['GAUSSIAN.template']
  cpto='%s/GAUSSIAN.template' % (iconddir)
  shutil.copy(cpfrom,cpto)

  if INFOS['gaussian.guess']:
    cpfrom1=INFOS['gaussian.guess']
    cpto1='%s/GAUSSIAN.chk.init' % (iconddir)
    shutil.copy(cpfrom1,cpto1)

  return


#======================================================================================================================
#======================================================================================================================
#======================================================================================================================

def checktemplate_ORCA(filename,INFOS):
  necessary=['basis','functional','charge']
  try:
    f=open(filename)
    data=f.readlines()
    f.close()
  except IOError:
    print 'Could not open template file %s' % (filename)
    return False
  valid=[]
  for i in necessary:
    for l in data:
      line=l.lower().split()
      if len(line)==0:
        continue
      line=line[0]
      if i==re.sub('#.*$','',line):
        valid.append(True)
        break
    else:
      valid.append(False)
  if not all(valid):
    print 'The template %s seems to be incomplete! It should contain: ' % (filename) +str(necessary)
    return False
  return True

# =================================================

def qmmm_job(filename,INFOS):
  necessary=['qmmm']
  try:
    f=open(filename)
    data=f.readlines()
    f.close()
  except IOError:
    print 'Could not open template file %s' % (filename)
    return False
  valid=[]
  for i in necessary:
    for l in data:
      line=l.lower().split()
      if len(line)==0:
        continue
      line=line[0]
      if i==re.sub('#.*$','',line):
        valid.append(True)
        break
    else:
      valid.append(False)
  if not all(valid):
    return False
  return True

# =================================================

def get_ORCA(INFOS):
  '''This routine asks for all questions specific to ORCA:
  - path to ORCA
  - scratch directory
  - ORCA.template
  - initial gbw file
  '''

  string='\n  '+'='*80+'\n'
  string+='||'+centerstring('ORCA Interface setup',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  print string

  print centerstring('Path to ORCA',60,'-')+'\n'
  print '\nPlease specify path to ORCA directory (SHELL variables and ~ can be used, will be expanded when interface is started).\n'
  INFOS['orcadir']=question('Path to ORCA:',str,'$ORCADIR')
  print ''




  # scratch
  print centerstring('Scratch directory',60,'-')+'\n'
  print 'Please specify an appropriate scratch directory. This will be used to run the ORCA calculations. The scratch directory will be deleted after the calculation. Remember that this script cannot check whether the path is valid, since you may run the calculations on a different machine. The path will not be expanded by this script.'
  INFOS['scratchdir']=question('Path to scratch directory:',str)
  print ''


  # template file
  print centerstring('ORCA input template file',60,'-')+'\n'
  print '''Please specify the path to the ORCA.template file. This file must contain the following keywords:

basis <basis>
functional <type> <name>
charge <x> [ <x2> [ <x3> ...] ]

The ORCA interface will generate the appropriate ORCA input automatically.
'''
  if os.path.isfile('ORCA.template'):
    if checktemplate_ORCA('ORCA.template',INFOS):
      print 'Valid file "ORCA.template" detected. '
      usethisone=question('Use this template file?',bool,True)
      if usethisone:
        INFOS['ORCA.template']='ORCA.template'
  if not 'ORCA.template' in INFOS:
    while True:
      filename=question('Template filename:',str)
      if not os.path.isfile(filename):
        print 'File %s does not exist!' % (filename)
        continue
      if checktemplate_ORCA(filename,INFOS):
        break
    INFOS['ORCA.template']=filename
  print ''


  # QMMM
  if qmmm_job(INFOS['ORCA.template'],INFOS):
    print centerstring('ORCA+TINKER QM/MM setup',60,'-')+'\n'
    print 'Your template specifies a QM/MM calculation. Please specify the path to TINKER.' 
    path=os.getenv('TINKER')
    if path=='':
      path=None
    else:
      path='$TINKER/'
    print '\nPlease specify path to TINKER bin/ directory (SHELL variables and ~ can be used, will be expanded when interface is started).\n'
    INFOS['tinker']=question('Path to TINKER/bin:',str,path)
    while True:
      filename=question('Force field file:',str)
      if not os.path.isfile(filename):
        print 'File %s does not exist!' % (filename)
      else:
        break
    INFOS['ORCA.fffile']=filename
    while True:
      filename=question('Connection table file:',str)
      if not os.path.isfile(filename):
        print 'File %s does not exist!' % (filename)
      else:
        break
    INFOS['ORCA.ctfile']=filename


  # initial MOs
  print centerstring('Initial restart: MO Guess',60,'-')+'\n'
  print '''Please specify the path to an ORCA gbw file containing suitable starting MOs for the ORCA calculation. Please note that this script cannot check whether the wavefunction file and the Input template are consistent!
'''
  if question('Do you have a restart file?',bool,True):
    if True:
      while True:
        filename=question('Restart file:',str,'ORCA.gbw')
        if os.path.isfile(filename):
          INFOS['orca.guess']=filename
          break
        else:
          print 'Could not find file "%s"!' % (filename)
  else:
    INFOS['orca.guess']={}


  # Resources
  print centerstring('ORCA Ressource usage',60,'-')+'\n'
  print '''Please specify the number of CPUs to be used by EACH calculation.
'''
  INFOS['orca.ncpu']=abs(question('Number of CPUs:',int)[0])

  if INFOS['orca.ncpu']>1:
    print '''Please specify how well your job will parallelize.
A value of 0 means that running in parallel will not make the calculation faster, a value of 1 means that the speedup scales perfectly with the number of cores.'''
    INFOS['orca.scaling']=min(1.0,max(0.0,question('Parallel scaling:',float,[0.8])[0] ))
  else:
    INFOS['orca.scaling']=0.9
  INFOS['orca.mem']=question('Memory (MB):',int,[1000])[0]


  # Ionization
  #print '\n'+centerstring('Ionization probability by Dyson norms',60,'-')+'\n'
  #INFOS['ion']=question('Dyson norms?',bool,False)
  #if INFOS['ion']:
  if 'wfoverlap' in INFOS['needed']:
    print '\n'+centerstring('WFoverlap setup',60,'-')+'\n'
    INFOS['orca.wfoverlap']=question('Path to wavefunction overlap executable:',str,'$SHARC/wfoverlap.x')
    print ''
    print 'State threshold for choosing determinants to include in the overlaps'
    print 'For hybrids (and without TDA) one should consider that the eigenvector X may have a norm larger than 1'
    INFOS['orca.ciothres']=question('Threshold:',float,[0.99])[0]
    print ''

    ## PyQuante
    #print '\n'+centerstring('PyQuante setup',60,'-')+'\n'
    #INFOS['orca.pyquante']=question('Path to PyQuante lib directory:',str,'$PYQUANTE')
    #print ''


  # TheoDORE
  theodore_spelling=['Om',
                    'PRNTO',
                    'Z_HE', 'S_HE', 'RMSeh',
                    'POSi', 'POSf', 'POS',
                    'PRi', 'PRf', 'PR', 'PRh',
                    'CT', 'CT2', 'CTnt',
                    'MC', 'LC', 'MLCT', 'LMCT', 'LLCT',
                    'DEL', 'COH', 'COHh']
  #INFOS['theodore']=question('TheoDORE analysis?',bool,False)
  if 'theodore' in INFOS['needed']:
    print '\n'+centerstring('Wave function analysis by TheoDORE',60,'-')+'\n'

    INFOS['orca.theodore']=question('Path to TheoDORE directory:',str,'$THEODIR')
    print ''

    print 'Please give a list of the properties to calculate by TheoDORE.\nPossible properties:'
    string=''
    for i,p in enumerate(theodore_spelling):
      string+='%s ' % (p)
      if (i+1)%8==0:
        string+='\n'
    print string
    l=question('TheoDORE properties:',str,'Om  PRNTO  S_HE  Z_HE  RMSeh')
    if '[' in l:
      INFOS['theodore.prop']=ast.literal_eval(l)
    else:
      INFOS['theodore.prop']=l.split()
    print ''

    print 'Please give a list of the fragments used for TheoDORE analysis.'
    print 'You can use the list-of-lists from dens_ana.in'
    print 'Alternatively, enter all atom numbers for one fragment in one line. After defining all fragments, type "end".'
    if qmmm_job(INFOS['ORCA.template'],INFOS):
      print 'You should only include the atom numbers of QM and link atoms.'
    INFOS['theodore.frag']=[]
    while True:
      l=question('TheoDORE fragment:',str,'end')
      if 'end' in l.lower():
        break
      if '[' in l:
        try:
          INFOS['theodore.frag']=ast.literal_eval(l)
          break
        except ValueError:
          continue
      f=[ int(i) for i in l.split() ]
      INFOS['theodore.frag'].append(f)
    INFOS['theodore.count']=len(INFOS['theodore.prop'])+len(INFOS['theodore.frag'])**2


  return INFOS

# =================================================

def prepare_ORCA(INFOS,iconddir):
  # write ORCA.resources
  try:
    sh2cas=open('%s/ORCA.resources' % (iconddir), 'w')
  except IOError:
    print 'IOError during prepareORCA, iconddir=%s' % (iconddir)
    quit(1)
#  project='ORCA'
  string='orcadir %s\nscratchdir %s/%s/\nncpu %i\nschedule_scaling %f\n' % (INFOS['orcadir'],INFOS['scratchdir'],iconddir,INFOS['orca.ncpu'],INFOS['orca.scaling'])
  string+='memory %i\n' % (INFOS['orca.mem'])
  if 'wfoverlap' in INFOS['needed']:
    string+='wfoverlap %s\nwfthres %f\n' % (INFOS['orca.wfoverlap'],INFOS['orca.ciothres'])
    #string+='numfrozcore %i\n' %(INFOS['frozcore_number'])
  else:
    string+='nooverlap\n'
  if 'theodore' in INFOS['needed']:
    string+='theodir %s\n' % (INFOS['orca.theodore'])
    string+='theodore_prop %s\n' % (INFOS['theodore.prop'])
    string+='theodore_fragment %s\n' % (INFOS['theodore.frag'])
  if 'tinker' in INFOS:
    string+='tinker %s\n' % (INFOS['tinker'])
  if 'ORCA.fffile' in INFOS:
    string+='qmmm_ff_file ORCA.qmmm.ff\n'
  if 'ORCA.ctfile' in INFOS:
    string+='qmmm_table ORCA.qmmm.table\n'
  sh2cas.write(string)
  sh2cas.close()

  # copy MOs and template
  cpfrom=INFOS['ORCA.template']
  cpto='%s/ORCA.template' % (iconddir)
  shutil.copy(cpfrom,cpto)

  if INFOS['orca.guess']:
    cpfrom1=INFOS['orca.guess']
    cpto1='%s/ORCA.gbw.init' % (iconddir)
    shutil.copy(cpfrom1,cpto1)

  if 'ORCA.fffile' in INFOS:
    cpfrom1=INFOS['ORCA.fffile']
    cpto1='%s/ORCA.qmmm.ff' % (iconddir)
    shutil.copy(cpfrom1,cpto1)

  if 'ORCA.ctfile' in INFOS:
    cpfrom1=INFOS['ORCA.ctfile']
    cpto1='%s/ORCA.qmmm.table' % (iconddir)
    shutil.copy(cpfrom1,cpto1)

  return




# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def checktemplate_BAGEL(filename,INFOS):
  necessary=['basis','df_basis','nact','nclosed']
  try:
    f=open(filename)
    data=f.readlines()
    f.close()
  except IOError:
    print 'Could not open template file %s' % (filename)
    return False
  valid=[]
  for i in necessary:
    for l in data:
      if i in re.sub('#.*$','',l):
        valid.append(True)
        break
    else:
      valid.append(False)
  if not all(valid):
    print 'The template %s seems to be incomplete! It should contain: ' % (filename) +str(necessary)
    return False
  roots_there=False
  for l in data:
    l=re.sub('#.*$','',l).lower().split()
    if len(l)==0:
      continue
    if 'nstate' in l[0]:
      roots_there=True
  if not roots_there:
    for mult,state in enumerate(INFOS['states']):
      if state<=0:
        continue
      valid=[]
      for l in data:
        if 'spin' in re.sub('#.*$','',l).lower():
          f=l.split()
          if int(f[1])==mult+1:
            valid.append(True)
            break
      else:
        valid.append(False)
  if not all(valid):
    string='The template %s seems to be incomplete! It should contain the keyword "spin" for ' % (filename)
    for mult,state in enumerate(INFOS['states']):
      if state<=0:
        continue
      string+='%s, ' % (IToMult[mult+1])
    string=string[:-2]+'!'
    print string
    return False
  return True

# =================================================

def get_BAGEL(INFOS):
  '''This routine asks for all questions specific to BAGEL:
  - path to bagel
  - scratch directory
  - BAGEL.template
  - wf.init
  '''

  string='\n  '+'='*80+'\n'
  string+='||'+centerstring('BAGEL Interface setup',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  print string

  print centerstring('Path to BAGEL',60,'-')+'\n'
  path=os.getenv('BAGEL')
  #path=os.path.expanduser(os.path.expandvars(path))
  if path=='':
    path=None
  else:
    path='$BAGEL/'
      #print 'Environment variable $MOLCAS detected:\n$MOLCAS=%s\n' % (path)
      #if question('Do you want to use this MOLCAS installation?',bool,True):
        #INFOS['molcas']=path
    #if not 'molcas' in INFOS:
  print '\nPlease specify path to BAGEL directory (SHELL variables and ~ can be used, will be expanded when interface is started).\n'
  INFOS['bagel']=question('Path to BAGEL:',str,path)
  print ''


  print centerstring('Scratch directory',60,'-')+'\n'
  print 'Please specify an appropriate scratch directory. This will be used to temporally store the integrals. The scratch directory will be deleted after the calculation. Remember that this script cannot check whether the path is valid, since you may run the calculations on a different machine. The path will not be expanded by this script.'
  INFOS['scratchdir']=question('Path to scratch directory:',str)
  print ''


  print centerstring('BAGEL input template file',60,'-')+'\n'
  print '''Please specify the path to the BAGEL.template file. This file must contain the following settings:

basis <Basis set>
df_basis <Density fitting basis set>
nact <Number of active orbitals>
nclosed <Number of doubly occupied orbitals>
nstate <Number of states for state-averaging>

The BAGEL interface will generate the appropriate BAGEL input automatically.
'''
  if os.path.isfile('BAGEL.template'):
    if checktemplate_BAGEL('BAGEL.template',INFOS):
      print 'Valid file "BAGEL.template" detected. '
      usethisone=question('Use this template file?',bool,True)
      if usethisone:
        INFOS['bagel.template']='BAGEL.template'
  if not 'bagel.template' in INFOS:
    while True:
      filename=question('Template filename:',str)
      if not os.path.isfile(filename):
        print 'File %s does not exist!' % (filename)
        continue
      if checktemplate_BAGEL(filename,INFOS):
        break
    INFOS['bagel.template']=filename
  print ''

  print centerstring('Dipole level',60,'-')+'\n'
  print 'Please specify the desired amount of calculated dipole moments:\n0 -only dipole moments that are for free are calculated\n1 -calculate all transition dipole moments between the (singlet) ground state and all singlet states for absorption spectra\n2 -calculate all dipole moments'
  INFOS['dipolelevel']=question('Requested dipole level:',int,[0])[0]
  print ''


  # QMMM 
#  if check_MOLCAS_qmmm(INFOS['molcas.template']):
#    print centerstring('MOLCAS+TINKER QM/MM setup',60,'-')+'\n'
#    print 'Your template specifies a QM/MM calculation. Please specify the path to TINKER.' 
#    path=os.getenv('TINKER')
#    if path=='':
#      path=None
#    else:
#      path='$TINKER/'
#    print '\nPlease specify path to TINKER bin/ directory (SHELL variables and ~ can be used, will be expanded when interface is started).\n'
#    INFOS['tinker']=question('Path to TINKER/bin:',str,path)
#    print 'Please give the key and connection table files.'
 #   while True:
 ##     filename=question('Key file:',str)
#      if not os.path.isfile(filename):
#        print 'File %s does not exist!' % (filename)
#      else:
#        break
#    INFOS['MOLCAS.fffile']=filename
#    while True:
#      filename=question('Connection table file:',str)
#      if not os.path.isfile(filename):
#        print 'File %s does not exist!' % (filename)
#      else:
 #       break
#    INFOS['MOLCAS.ctfile']=filename



  print centerstring('Initial wavefunction: MO Guess',60,'-')+'\n'
  print '''Please specify the path to a MOLCAS JobIph file containing suitable starting MOs for the CASSCF calculation. Please note that this script cannot check whether the wavefunction file and the Input template are consistent!
'''
  INFOS['bagel.guess']={}
  string='Do you have initial wavefunction files for '
  for mult,state in enumerate(INFOS['states']):
    if state<=0:
      continue
    string+='%s, ' % (IToMult[mult+1])
  string=string[:-2]+'?'
  if question(string,bool,True):
    for mult,state in enumerate(INFOS['states']):
      if state<=0:
        continue
      while True:
        guess_file='archive.%i.init' % (mult+1)
        filename=question('Initial wavefunction file for %ss:' % (IToMult[mult+1]),str,guess_file)
        if os.path.isfile(filename):
          INFOS['bagel.guess'][mult+1]=filename
          break
        else:
          print 'File not found!'
  else:
    print 'WARNING: Remember that CASSCF calculations may run very long and/or yield wrong results without proper starting MOs.'
    time.sleep(1)

  print centerstring('BAGEL Ressource usage',60,'-')+'\n'#TODO

  print '''Please specify the number of CPUs to be used by EACH calculation.
'''
  INFOS['bagel.ncpu']=abs(question('Number of CPUs:',int,[1])[0])
  
  if INFOS['bagel.ncpu']>1:
    INFOS['bagel.mpi']=question('Use MPI mode (no=OpenMP)?',bool,False)
  else:
    INFOS['bagel.mpi']=False





  ## Ionization
  #need_wfoverlap=False
  #print centerstring('Ionization probability by Dyson norms',60,'-')+'\n'
  #INFOS['ion']=question('Dyson norms?',bool,False)
  #if 'ion' in INFOS and INFOS['ion']:
    #need_wfoverlap=True

  # wfoverlap
  if 'wfoverlap' in INFOS['needed']:
    print '\n'+centerstring('WFoverlap setup',60,'-')+'\n'
    INFOS['bagel.wfoverlap']=question('Path to wavefunction overlap executable:',str,'$SHARC/wfoverlap.x')
    # TODO not asked for: numfrozcore, numocc
    print '''Please specify the path to the PyQuante directory.
'''
    INFOS['bagel.pyq']=question('PyQuante path:',str)
    print '''Please specify the amount of memory available to wfoverlap.x (in MB). \n (Note that BAGEL's memory cannot be controlled)
'''
    INFOS['bagel.mem']=abs(question('wfoverlap.x memory:',int,[1000])[0])
  else:
    INFOS['bagel.mem']=1000
    INFOS['bagel.pyq']=''

  return INFOS

# =================================================

def prepare_BAGEL(INFOS,iconddir):
  # write BAGEL.resources
  try:
    sh2cas=open('%s/BAGEL.resources' % (iconddir), 'w')
  except IOError:
    print 'IOError during prepareBAGEL, iconddir=%s' % (iconddir)
    quit(1)
  project='BAGEL'
  string='bagel %s\npyquante %s\nscratchdir %s/%s/\nmemory %i\nncpu %i\ndipolelevel %i\nproject %s' % (INFOS['bagel'],INFOS['bagel.pyq'],INFOS['scratchdir'],iconddir,INFOS['bagel.mem'],INFOS['bagel.ncpu'],INFOS['dipolelevel'],project)
  
  if INFOS['bagel.mpi']:
    string+='mpi_parallel\n'
  if 'wfoverlap' in INFOS['needed']:
    string+='\nwfoverlap %s\n' % INFOS['bagel.wfoverlap']
  else:
    string+='\nnooverlap\n'
#  if 'tinker' in INFOS:
#    string+='tinker %s' % (INFOS['tinker'])
  sh2cas.write(string)
  sh2cas.close()

  # copy MOs and template
  cpfrom=INFOS['bagel.template']
  cpto='%s/BAGEL.template' % (iconddir)
  shutil.copy(cpfrom,cpto)
  if not INFOS['bagel.guess']=={}:
    for i in INFOS['bagel.guess']:
      cpfrom=INFOS['bagel.guess'][i]
      cpto='%s/%s.%i.init' % (iconddir,'archive',i)
      shutil.copy(cpfrom,cpto)



  return






# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def get_general():
  '''
  This routine questions from the user some general information:
  - V0.txt file
  - # of displacements
  - one/two-sided derivation
  - number of states
  - QM package
  - SOCs
  
  returns INFOS dictionary
  '''
  INFOS = {}


  ## -------------------- getting & reading V0.txt -------------------- ##
  print centerstring('V0.txt file', 60, '-') + '\n'

  # check if V0 exists
  v0file = 'V0.txt'
  try:
    if os.path.isfile(v0file):
      print 'Ground-state file "V0.txt" detected. Do you want to use this?'
      if not question('Use file "V0.txt"?', bool, True): raise IOError
    else: raise IOError
  except IOError:
    print '\nIf you do not have an ground-state file, prepare one with \'wigner.py -l <molden-file>\'!\n'
    print 'Please enter the filename of the ground-state file.'
    while True:
      v0file = question('Ground-state filename:', str, 'V0.txt')
      v0file = os.path.expanduser(os.path.expandvars(v0file))
      if os.path.isdir(v0file):
        print 'Is a directory: %s' % (v0file)
        continue
      if not os.path.isfile(v0file):
        print 'File does not exist: %s' % (v0file)
        continue
      if os.path.isfile(v0file): break

  INFOS['v0f'] = os.path.abspath(v0file)

  # read V0.txt
  with open(v0file, 'r') as v0f:
    content = v0f.readlines()
    v0f.close()

  INFOS = read_V0(INFOS, content)

  # output to user
  print '\nFile "%s" contains %i atoms and we will use %i frequencies/normal modes.\n(others are zero)\n' % (v0file, len(INFOS['atoms']), len(INFOS['freqencies']))


  ## -------------------- number of states -------------------- ##
  print centerstring('Number of states', 60, '-')
  print '\nPlease enter the number of states as a list of integers\ne.g. 3 0 3 for three singlets, zero doublets and three triplets.'
  
  while True:
    states = question('Number of states:', int)

    if not states: continue
    if any(i < 0 for i in states):
      print 'Number of states must be positive!'
      continue
    break
    
  nstates = 0
  for mult, i in enumerate(states):
    nstates += (mult + 1) * i

  # output to user
  print '\nNumber of states: %s' % (states)
  print 'Total number of states: %i\n' % (nstates)

  # saving input
  INFOS['states'] = states
  INFOS['nstates'] = nstates


  ## -------------------- QM package -------------------- ##
  print centerstring('Choose the quantum chemistry interface', 60, '-')
  print '\nPlease specify the quantum chemistry interface (enter any of the following numbers):'
  
  # list available interfaces
  for i in Interfaces:
    print '%i\t%s' % (i, Interfaces[i]['description'])
  print ''
  
  # getting interface
  while True:
    num = question('Interface number:', int)[0]
    if num in Interfaces: break
    else:
      print 'Please input one of the following: %s!' % ([i for i in Interfaces])
  
  # output to user
  print 'The used interface will be: ' + Interfaces[num]['description'] + '\n'

  # save input
  INFOS['interface'] = num
  INFOS['needed'] = []


  ## -------------------- Setup SOCs -------------------- ##
  print centerstring('Spin-orbit couplings (SOCs)',60,'-')+'\n'
  if len(states) > 1:
    if 'soc' in Interfaces[num]['features']:
      print 'Do you want to compute spin-orbit couplings?\n'
      soc = question('Spin-Orbit calculation?', bool, True)
      if soc:
        print 'Will calculate spin-orbit matrix.'
    else:
      print 'Interface cannot provide SOCs: not calculating spin-orbit matrix.'
      soc = False
  else:
    print 'Only singlets specified: not calculating spin-orbit matrix.'
    soc = False
  print ''
  
  # save input
  INFOS['soc'] = soc
  if INFOS['soc']:
    INFOS['needed'].extend(Interfaces[num]['features']['soc'])


  ## -------------------- whether to do gradients or numerical -------------------- ##
  print centerstring('Analytical gradients', 60, '-')+'\n'

  INFOS['ana_grad'] = question('Do you want to use analytical gradients for kappa terms?', bool, True)

  print '\nAnalytical gradients for kappas: %r\n' % INFOS['ana_grad']


  ## -------------------- whether to do gradients or numerical -------------------- ##
  print centerstring('Analytical nonadiabatic coupling vectors', 60, '-')+'\n'
  if 'nacdr' in Interfaces[num]['features']:

    INFOS['ana_nac'] = question('Do you want to use analytical nonadiabatic coupling vectors for lambda terms?', bool, False)

  else:
    INFOS['ana_nac'] = False

  print 'Do you want to use analytical nonadiabatic coupling vectors for lambdas: %r\n' % INFOS['ana_nac']
  if INFOS['ana_nac']:
    INFOS['needed'].extend(Interfaces[num]['features']['nacdr'])

  ## -------------------- Whether to do overlaps -------------------- ##
  if (not INFOS['ana_grad']) or (not INFOS['ana_nac']):
    do_overlaps=True
  else:
    do_overlaps=False

  if do_overlaps:
    if not 'overlap' in Interfaces[num]['features']:
      print 'Interface cannot provide overlaps, numerical mode not possible!'
      sys.exit(1)

    INFOS['do_overlaps'] = True
    INFOS['needed'].extend(Interfaces[num]['features']['overlap'])
  else:
    INFOS['do_overlaps'] = False




  ## -------------------- select normal modes -------------------- ##
  print centerstring('Normal modes', 60, '-') + '\n'

  if not question('Do you want to make LVC parameters for all normal modes?', bool, True):
    not_using = question('Which of the normal modes do you NOT want to use? %s' % (reduce_big_list_to_short_str(INFOS['normal_modes'].keys())), int, ranges = True)

    # delete selected normal modes & jump over non-existing ones
    for k in not_using:
      if not (k in INFOS['normal_modes'].keys()): print 'Normal mode %i doesn\'t exist. Going on with the rest.' % (k)
      else: del INFOS['normal_modes'][k]

  # output to user
  print '\nWe will use the following normal modes: %s\n' % (reduce_big_list_to_short_str(INFOS['normal_modes'].keys()))






  ## ----------------------------------------------------------------- ##
  ## -------------------- Infos for displacements -------------------- ##
  if INFOS['do_overlaps']:


    ## -------------------- reading displacements -------------------- ##
    print centerstring('Displacements', 60, '-') + '\n'
    displacement_magnitudes = OrderedDict(sorted({ i: 0.05 for i in INFOS['normal_modes'].keys() }.items(), key = lambda t: t[0]))

    # getting user input
    if question('Do you want to use other displacements than the default of [0.05]?', bool, False):
      while True:
        print ''
        displacement = abs(question('Displacement magnitude:', float, [0.05])[0])
        
        if displacement == 0:
          print 'Only non-zero displacement magnitudes are allowed.'
          continue
        
        indexes = question('For which normal modes do you want to use the displacement of %g? %s:' % (displacement, reduce_big_list_to_short_str(INFOS['normal_modes'].keys())), int, ranges = True)
        
        # check if normal mode of input exists and add it to displacement_magnitudes dict
        for k in indexes:
          if k in INFOS['normal_modes'].keys():
            displacement_magnitudes[k] = displacement
          else: print 'Normal mode %i doesn\'t exist. Therefore can\'t set displacement. Going on with the rest.' % (k)
          
        print ''
        if not question('Do you want to define more displacements?', bool, False): break
    
    # saving input
    INFOS['displacement_magnitudes'] = displacement_magnitudes

    # output to user
    print '\nScript will use displacement magnitudes of:\n\n%s \n' % (reduce_displacement_dictionary_to_output_str(displacement_magnitudes))



    ## -------------------- ignore problematic states -------------------- ##
    if not (INFOS['ana_grad'] and INFOS['ana_nac']):
      print centerstring('Intruder states', 60, '-')

      print '''
  Intruder states can be detected by small overlap matrix elements.
  Affected numerical kappa/lambda terms will be ignored and not written to the parameter file.
'''
      INFOS['ignore_problematic_states'] = not question('Do you want to check for intruder states?', bool, True)

      print '\nIgnore problematic states: %r\n' % INFOS['ignore_problematic_states']

    else:
      INFOS['ignore_problematic_states'] = False



    ## -------------------- one/two-sided derivation -------------------- ##
    print centerstring('One-/Two-sided derivation', 60, '-')
    print '\nOne-/Two-sided derivation of normal modes.'

    # getting user input
    one_sided_derivations = {}
    n_normal_modes = len(INFOS['normal_modes'])
      
    derivations = question('Choose for which normal modes you want to use one-sided derivation %s:' % (reduce_big_list_to_short_str(INFOS['normal_modes'])), int, [None], ranges = True)

    # check if normal mode of input exists and add it to one_sided_derivation dict
    if derivations != [None]:
      for k in derivations:
        if k in INFOS['normal_modes'].keys():
          one_sided_derivations[k] = True
        else: print 'Normal mode %i doesn\'t exist. Therefore can\'t set one-sided derivation. Going on with the rest.' % (k)

    #saving input
    INFOS['one-sided_derivations'] = OrderedDict(sorted(one_sided_derivations.items(), key = lambda t: t[0]))

    # output to user
    print '\nOne-sided derivation will be used on: %s\n' % (reduce_big_list_to_short_str(one_sided_derivations.keys()))






  ### -------------------- number of frozen cores -------------------- ##
  #print centerstring('Number of frozen cores', 60, '-')
  #if question('\nDo you want to set frozen cores?', bool, False):
    #print '\nPlease enter the amount of frozen cores you want to have.'
    
    #while True:
      #frozen_cores = question('Number of frozen cores:', int)

      #if not frozen_cores: continue
      #if frozen_cores < 0:
        #print 'Number of frozen cores must be equal/bigger than zero!'
        #continue
      #break

    ## output to user
    #print '\nNumber of frozen cores: %s\n' % (frozen_cores[0])

    ## saving input
    #INFOS['frozen_cores'] = frozen_cores[0]
  #else:
    #INFOS['frozen_cores'] = -1




  ## -------------------- Calculate displacements -------------------- ##
  INFOS = calculate_displacements(INFOS)


  ## -------------------- Create path information -------------------- ##
  INFOS['result_path'] = 'DSPL_RESULTS'
  INFOS['paths'] = {}
  INFOS['paths']['0eq'] = 'DSPL_%03i_eq' % (0)
  if INFOS['do_overlaps']:
    INFOS['paths'].update({ k: 'DSPL_%03i_%s' % (int(k[:-1]), k[-1]) for k, v in INFOS['displacements'].items() })
  # sort dict
  INFOS['paths'] = OrderedDict(sorted(INFOS['paths'].items(), key = lambda t: t[0]))




  INFOS['cwd'] = os.getcwd()
  return INFOS

# ======================================================================================================================

def read_V0(INFOS, content):
  '''
  reads V0.txt and puts the data into the INFOS dictionary

  returns INFOS dictionary
  '''
  # current header
  header = ''
  # set headers of V0.txt file
  headers = { 'geo': 'Geometry\n', 'freq': 'Frequencies\n', 'mwnmodes': 'Mass-weighted normal modes\n' }
  # init list/dicts
  INFOS['atoms'], freqencies, normal_modes = [], {}, {}

  for line in content:
    # check if within atom lines
    if header == headers['geo'] and line != headers['freq']:
      elements = line.strip().split()
      
      # add atom
      INFOS['atoms'].append({'atom': elements[0],
        'e-': int(float(elements[1])),
        'coords [bohr]': (float(elements[2]), float(elements[3]), float(elements[4])),
        'mass (amu)': float(elements[5]) * U_TO_AMU
       })
    
    # check if within frequencies and add them 
    if header == headers['freq'] and line != headers['mwnmodes']:
      freqencies = { i + 1: freq for i, freq in zip(range(len(line.split())), map(float, line.strip().split())) }
      
      # init number of normal_modes with lists, as to be able to append them later
      for i in range(len(freqencies)):
        normal_modes[i + 1] = []
    
    # within normal modes
    if header == headers['mwnmodes']:
      elements = line.strip().split()
      
      # every column is a normal mode - assign accordingly
      for i in range(len(elements)):
        normal_modes[i + 1].append(float(elements[i]))
    
    # change current header
    if line in headers.values():
      header = line

  # save as ordered dict for nice output later on
  INFOS['normal_modes'] = OrderedDict(sorted({ i: normal_mode for i, normal_mode in normal_modes.items() if freqencies[i] != 0 }.items(), key = lambda t: t[0]))
  INFOS['freqencies'] = OrderedDict(sorted({ i: freq for i, freq in freqencies.items() if freq != 0 }.items(), key = lambda t: t[0]))

  return INFOS

# ======================================================================================================================

def get_runscript_info(INFOS):
  '''
  Gets all the necessary information from the user for the runscripts

  returns INFOS dictionary
  '''
  
  string =  '\n  ' + '=' * 80
  string += '\n||' + centerstring('Run mode setup', 80) + '||'
  string += '\n  ' + '=' * 80 + '\n\n'
  print string

  print centerstring('Run script', 60, '-') + '\n'
  print '''  This script can generate the run scripts for each initial condition in two modes:

    - In mode 1, the calculation is run in subdirectories of the current directory.

    - In mode 2, the input files are transferred to another directory (e.g. a local
      scratch directory), the calculation is run there, results are copied back and
      the temporary directory is deleted. Note that this temporary directory is not
      the same as the "scratchdir" employed by the interfaces.

  Note that in any case this script will create the input subdirectories in the
  current working directory.'''

  print '\n  In case of mode 1, the calculations will be run in: \'%s\'\n' % (INFOS['cwd'])
  INFOS['here'] = question('Use mode 1 (i.e., calculate here)?', bool, True)
  
  if not INFOS['here']:
    print '\nWhere do you want to perform the calculations? Note that this script cannot check\nwhether the path is valid.'
    INFOS['copydir'] = question('Run directory?', str)
  print ''

  print centerstring('Submission script', 60, '-') + '\n'
  print 'During the setup, a script for running all initial conditions sequentially in batch\nmode is generated. Additionally, a queue submission script can be generated for all\ninitial conditions.'
  INFOS['qsub'] = question('Generate submission script?', bool, False)
  
  if INFOS['qsub']:
    INFOS['qsub'] = True
    print '\nPlease enter a queue submission command, including possibly options to the queueing\nsystem, e.g. for SGE: "qsub -q queue.q -S /bin/bash -cwd" (Do not type quotes!).'
    INFOS['qsubcommand'] = question('Submission command?', str, None, False)
    INFOS['proj'] = question('Project Name:', str, None, False)

  print ''
  return INFOS

# ======================================================================================================================

def calculate_displacements(INFOS):
  '''
  Calculates all displacements to set up the calculations
  and the full transformation matrix of dimensionless mass-weighted coordinates
  
  returns INFOS dictionary
  '''

  # dividing normal modes by sqrt(frequency)
  fw_normal_modes = {}
  for i, normal_mode in INFOS['normal_modes'].items():
    fw_normal_modes[i] = ([nm / (INFOS['freqencies'][i] ** 0.5) for nm in normal_mode])

  # dividing the normal modes by sqrt(atom_mass)
  fmw_normal_modes = {}
  for i, fw_normal_mode in fw_normal_modes.items():
    j = 0
    fmw_normal_mode = []
    for atom in INFOS['atoms']:
      fmw_normal_mode.append(fw_normal_mode[j] / (atom['mass (amu)'] ** 0.5))
      fmw_normal_mode.append(fw_normal_mode[j + 1] / (atom['mass (amu)'] ** 0.5))
      fmw_normal_mode.append(fw_normal_mode[j + 2] / (atom['mass (amu)'] ** 0.5))
      j += 3
      
    fmw_normal_modes[i] = fmw_normal_mode
  
  # writing frequency and mass weighted normal modes to dict
  INFOS['fmw_normal_modes'] = fmw_normal_modes

  # writing displacements by multiplying normal modes with displacement magnitudes
  displacements = {}
  if INFOS['do_overlaps']:
    for k, normal_mode in fmw_normal_modes.items():
      displacements[str(k) + 'p'] = ([nm * INFOS['displacement_magnitudes'][k] for nm in normal_mode])

      # for two sided derivation
      if not k in INFOS['one-sided_derivations']:
        displacements[str(k) + 'n'] = ([nm * INFOS['displacement_magnitudes'][k] * (-1) for nm in normal_mode])

  # sort displacements and return
  INFOS['displacements'] = OrderedDict(sorted(displacements.items(), key = lambda t: t[0]))
  return INFOS

# ======================================================================================================================

def write_QM_in(INFOS, displacement_key, displacement_value, displacement_dir):
  '''
  Writes QM.in file for displacement calculations
  '''

  # open writable QM.in file
  try:
    qmin = open('%s/QM.in' % (displacement_dir), 'w')
  except IOError:
    print 'IOError during write_QM_in, displacement_dir=%s' % (displacement_dir)
    quit(1)

  # number of atoms
  string = '%i\n' % (len(INFOS['atoms']))

  # number of current initial condition
  string += 'Displacement %s\n' % (displacement_key)

  # add eq
  if displacement_key == None:
    for atom in INFOS['atoms']:
      string += '%s %f %f %f\n' % (atom['atom'], atom['coords [bohr]'][0], atom['coords [bohr]'][1], atom['coords [bohr]'][2])
  
  # for non eq add displacements
  else:
    i = 0
    for atom in INFOS['atoms']:
      string += '%s %f %f %f\n' % (atom['atom'],
                                   atom['coords [bohr]'][0] + INFOS['displacements'][displacement_key][i],
                                   atom['coords [bohr]'][1] + INFOS['displacements'][displacement_key][i + 1],
                                   atom['coords [bohr]'][2] + INFOS['displacements'][displacement_key][i + 2])
      i += 3

  # unit def
  string += 'unit bohr\n'

  # states def
  string += 'states '
  for i in INFOS['states']: string += '%i ' % (i)
  string += '\n'

  # eq: init ; displacement: overlap + cleanup
  if displacement_key == None: string += 'init\n'
  else:
    string += 'overlap\n'
    string += 'cleanup\n'

  # write molden file
  string += 'molden\n'
  # set savedir
  string += 'savedir ./SAVE/\n'

  

  # spin orbit coupling
  if INFOS['soc']: string += '\nSOC\n'
  else: string += '\nH\n'

  # dipole moment
  string += 'DM\n'

  # gradient
  if displacement_key is None and INFOS['ana_grad']:
    string += 'GRAD\n'
  if displacement_key is None and INFOS['ana_nac']:
    string += 'NACDR\n'

  # if theodore set it
  if displacement_key is None and 'theodore' in INFOS and INFOS['theodore']: 
    string += 'theodore\n'

  qmin.write(string)
  qmin.close()

# ======================================================================================================================

def write_runscript(INFOS, displacement_dir):
  '''
  writes the runscript in each subdirectory
  '''

  filename = '%s/run.sh' % (displacement_dir)
  try:
    runscript = open(filename, 'w')
  except IOError:
    print 'IOError during write_runscript, displacement_dir = %s' % (displacement_dir)
    quit(1)
    
  if 'proj' in INFOS:
    projname='%4s_%5s' % (INFOS['proj'][0:4], displacement_dir[-6:-1])
  else:
    projname='init_%5s' % (displacement_dir[-6:-1])


  intstring=''
  if 'adfrc' in INFOS:
    intstring = '. %s\nexport PYTHONPATH=$ADFHOME/scripting:$PYTHONPATH' % (INFOS['adfrc'])


  refstring=''
  if displacement_dir != INFOS['result_path'] + '/' + INFOS['paths']['0eq']:
    refstring = 'if [ -d ../../%s/SAVE ];\nthen\n  if [ -d ./SAVE ];\n  then\n    rm -r ./SAVE\n  fi\n  cp -r ../../%s/SAVE ./\nelse\n  echo "Should do a reference overlap calculation, but the reference data in ../../%s/ seems not OK."\n  exit 1\nfi' % (INFOS['result_path'] + '/' + INFOS['paths']['0eq'], INFOS['result_path'] + '/' + INFOS['paths']['0eq'], INFOS['result_path'] + '/' + INFOS['paths']['0eq'])
    
  ## generate run scripts here
  # for here mode
  if INFOS['here']:
    string = '#!/bin/bash\n\n#$-N %s\n\n%s\n\nPRIMARY_DIR=%s/%s/\n\ncd $PRIMARY_DIR\n%s\n\n$SHARC/%s QM.in >> QM.log 2>> QM.err\n' % (projname, intstring, INFOS['cwd'], displacement_dir, refstring, Interfaces[INFOS['interface']]['script'])
  # for remote mode
  else:
    string = '#!/bin/bash\n\n#$-N %s\n\n%s\n\nPRIMARY_DIR=%s/%s/\nCOPY_DIR=%s/%s/\n\ncd $PRIMARY_DIR\n%s\n\nmkdir -p $COPY_DIR\ncp -r $PRIMARY_DIR/* $COPY_DIR\ncd $COPY_DIR\n\n$SHARC/%s QM.in >> QM.log 2>> QM.err\n\ncp -r $COPY_DIR/QM.* $COPY_DIR/SAVE/ $PRIMARY_DIR\nrm -r $COPY_DIR\n' % (projname,intstring,INFOS['cwd'], displacement_dir, INFOS['copydir'], displacement_dir, refstring, Interfaces[INFOS['interface']]['script'])

  
  # run, close & make executable
  runscript.write(string)
  runscript.close()
  os.chmod(filename, os.stat(filename).st_mode | stat.S_IXUSR)

  return

# ======================================================================================================================

def write_displacement_info(INFOS):
  '''
  Writes displacement info to log-file and INFOS to info-file
  '''

  displacement_log_filename = '%s/displacements.log' % (INFOS['result_path'])
  displacement_info_filename = '%s/displacements.json' % (INFOS['result_path'])

  # open writable displacements.info file
  try:
    displacement_log = open(displacement_log_filename, 'w')
    displacement_info = open(displacement_info_filename, 'w')
  except IOError:
    print 'IOError during opening writeable %s or %s - file. Quitting.' % (displacement_log_filename, displacement_info_filename)
    quit(1)

  # write INFOS to info file
  #pickle.dump(INFOS, displacement_info)
  json.dump(INFOS, displacement_info, sort_keys=True, indent=4)

  # writing header
  displacement_log.write('natoms: %i\nstates: %s\nnstates: %i\n\nnormal_mode displacement_magnitude p/n path\n' % (len(INFOS['atoms']), str(INFOS['states'])[1:-1].replace(',', ''), INFOS['nstates']))

  # writing eq
  displacement_log.write('-%s-%s-%s%s/\n' % (' ' * 11, ' ' * 22, ' ' * 3, INFOS['result_path'] + '/' + INFOS['paths']['0eq']))

  # writing all displacements to log file
  if INFOS['do_overlaps']:
    for displacement_key, v in INFOS['displacements'].items():
      normal_mode = int(displacement_key[:-1])

      displacement_log.write('%i%s%g%s%c%s%s\n' % (normal_mode,
                                                  ' ' * (12 - len(str(normal_mode))),
                                                  INFOS['displacement_magnitudes'][normal_mode],
                                                  ' ' * (23 - len(str(INFOS['displacement_magnitudes'][normal_mode]))),
                                                  displacement_key[-1],
                                                  ' ' * 3,
                                                  INFOS['result_path'] + '/' + INFOS['paths'][displacement_key]))

  displacement_log.close()
  displacement_info.close()

# ======================================================================================================================

#def write_frozen_cores(INFOS):
  #for path in INFOS['paths'].values():
    #relative_path = INFOS['result_path'] + '/' + path

    #for file_name in os.listdir(relative_path):
      #if file_name.endswith('.resources'):
        #try:
          #file = open(relative_path + '/' + file_name, 'a')
          #file.write('numfrozcore %i\n' % INFOS['frozen_cores'])
          #file.close()
        #except IOError:
          #print 'Interface resource file was not found in %s\nCouldn\'t write frozen cores!' % (path)

# ======================================================================================================================

def setup_equilibrium(INFOS):
  '''
  Sets up the eq condition
  '''
  
  displacement_dir = INFOS['result_path'] + '/' + INFOS['paths']['0eq']
  io = make_directory(displacement_dir)
  
  # eq condition already exists
  if io != 0:
    print 'Skipping equlibrium %s!' % (displacement_dir)
    return True

  # create eq condition
  write_QM_in(INFOS, None, None, displacement_dir)
  #getattr(interfaces, [Interfaces[INFOS['interface']]['prepare_routine']][0])(INFOS, displacement_dir)
  globals()[Interfaces[ INFOS['interface']]['prepare_routine'] ](INFOS, displacement_dir)
  write_runscript(INFOS, displacement_dir)
  
  return False

# ======================================================================================================================

def setup_all(INFOS):
  '''
  This routine sets up the directories for the initial calculations.
  '''

  if make_directory(INFOS['result_path']) == -1:
    print 'Results folder will not be created or overwritten. Quitting.'
    quit(1)

  string =  '\n  ' + '=' * 80
  string += '\n||' + centerstring('Setting up directories...', 80) + '||'
  string += '\n  ' + '=' * 80 + '\n\n'
  print string

  # define local variables
  all_run_filename = '%s/all_run_dspl.sh' % (INFOS['result_path'])
  all_qsub_filename = '%s/all_qsub_dspl.sh' % (INFOS['result_path'])


  # write current working directory to all_run_dspl.sh & all_qsub_dspl.sh
  all_run = open(all_run_filename, 'w')
  string = '#/bin/bash\n\nCWD=%s\n\n' % (INFOS['cwd'])
  all_run.write(string)
  
  # add queueing script if wanted
  if INFOS['qsub']:
    all_qsub = open(all_qsub_filename, 'w')
    string = '#/bin/bash\n\nCWD=%s\n\n' % (INFOS['cwd'])
    all_qsub.write(string)

  # if eq doesn't exist yet, set it up & add it to the run_all* files
  if not setup_equilibrium(INFOS):
    eq_dir = INFOS['result_path'] + '/' + INFOS['paths']['0eq']

    string = 'cd $CWD/%s/\nbash run.sh\ncd $CWD\necho %s >> %s/DONE\n' % (eq_dir, eq_dir, INFOS['result_path'])
    all_run.write(string)
    
    if INFOS['qsub']:
      string = 'cd $CWD/%s/\n%s run.sh\ncd $CWD\n' % (eq_dir, INFOS['qsubcommand'])
      all_qsub.write(string)


  ## set up all displacement calculations
  if INFOS['do_overlaps']:
    dispacements_done = 0
    width_progressbar = 50
    number_of_displacements = len(INFOS['displacements'].keys())

    # iterating through displacements
    for displacement_key, displacement_value in INFOS['displacements'].items():
      displacement_dir = INFOS['result_path'] + '/' + INFOS['paths'][displacement_key]
      
      dispacements_done += 1
      done = dispacements_done * width_progressbar / number_of_displacements
      
      sys.stdout.write('\rProgress: [' + '=' * done + ' ' * (width_progressbar - done) + '] %3i%%' % (done * 100 / width_progressbar))
      sys.stdout.flush()
      
      if make_directory(displacement_dir) != 0:
        print 'Skipping displacement %s!' % (displacement_dir)
        continue
      
      # write QM.in, interfaces & runscript for displacement
      write_QM_in(INFOS, displacement_key, displacement_value, displacement_dir)
      #getattr(interfaces, [Interfaces[INFOS['interface']]['prepare_routine']][0])(INFOS, displacement_dir)
      globals()[Interfaces[ INFOS['interface']]['prepare_routine'] ](INFOS, displacement_dir)
      write_runscript(INFOS, displacement_dir)
      
      string = 'cd $CWD/%s/\nbash run.sh\ncd $CWD\necho %s >> %s/DONE\n' % (displacement_dir, displacement_dir, INFOS['result_path'])
      all_run.write(string)
      
      if INFOS['qsub']:
        string = 'cd $CWD/%s/\n%s run.sh\ncd $CWD\n' % (displacement_dir, INFOS['qsubcommand'])
        all_qsub.write(string)

  
  # close filehandlers & make executable
  all_run.close()
  os.chmod(all_run_filename, os.stat(all_run_filename).st_mode | stat.S_IXUSR)
  
  if INFOS['qsub']:
    all_qsub.close()
    os.chmod(all_qsub_filename, os.stat(all_qsub_filename).st_mode | stat.S_IXUSR)


  # write displacement info
  write_displacement_info(INFOS)

  ## write frozen cores
  #if INFOS['frozen_cores'] != -1:
    #write_frozen_cores(INFOS)

  print '\n'

# ======================================================================================================================

def main():
  '''
  Main routine
  '''

  usage='''python setup_displacement.py'''

  parser = OptionParser(usage = usage, description = '')

  displaywelcome()
  open_keystrokes()

  # get general information for calcultion
  INFOS = get_general()
  # get interface info - use reflection to get chosen routine
  #INFOS = getattr(interfaces, [Interfaces[INFOS['interface']]['get_routine']][0])(INFOS)
  INFOS = globals()[Interfaces[INFOS['interface']]['get_routine']](INFOS)
  INFOS = get_runscript_info(INFOS)

  print_INFOS(INFOS)

  if question('Do you want to setup the specified calculations?', bool, True):
    setup_all(INFOS)
    #shutil.move('KEYSTROKES.tmp', 'KEYSTROKES.setup_displacement')

  close_keystrokes()


if __name__ == '__main__':
  try:
    main()
  except KeyboardInterrupt:
    print '\nCtrl+C makes me a sad SHARC ;-(\n'
    quit(0)
