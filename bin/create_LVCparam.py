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

import os
import sys
#import pickle
import json
import itertools
import numpy as np
from optparse import OptionParser



def json_load_byteified(file_handle):
    return _byteify(
        json.load(file_handle, object_hook=_byteify),
        ignore_dicts=True
    )

def json_loads_byteified(json_text):
    return _byteify(
        json.loads(json_text, object_hook=_byteify),
        ignore_dicts=True
    )

def _byteify(data, ignore_dicts = False):
    # if this is a unicode string, return its string representation
    if isinstance(data, unicode):
        return data.encode('utf-8')
    # if this is a list of values, return list of byteified values
    if isinstance(data, list):
        return [ _byteify(item, ignore_dicts=True) for item in data ]
    # if this is a dictionary, return dictionary of byteified keys and values
    # but only if we haven't already byteified it
    if isinstance(data, dict) and not ignore_dicts:
        return {
            _byteify(key, ignore_dicts=True): _byteify(value, ignore_dicts=True)
            for key, value in data.iteritems()
        }
    # if it's anything else, return it in its original form
    return data


#from modules import constants
#from modules import output
# replace by import from official source below (QMout2LVC needs to be extended with flag 6 for overlaps to work)
# 'overlap':      {'flag': 6,
#                  'type': complex,
#                  'dim':  (nstates,nstates)}
#from modules.QMout2LVC import read_QMout, LVC_complex_mat
#from modules.SHARC_LVC import itnmstates
#import imp
#from imp.load_source('', os.environ['SHARC'] + 'QMout2LVC.py') import read_QMout
#from imp.load_source('', os.environ['SHARC'] + 'SHARC_LVC.py') import itnmstates



#import sys 
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


# ======================================================================= #

#from math import pi

#CM_TO_HARTREE = 1./219474.6     # 4.556335252e-6 # conversion factor from cm-1 to Hartree
#HARTREE_TO_EV = 27.211396132    # conversion factor from Hartree to eV
U_TO_AMU = 1./5.4857990943e-4   # conversion from g/mol to amu
#ANG_TO_BOHR = 1./0.529177211    # 1.889725989    # conversion from Angstrom to bohr
#PI = pi

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

#NUMBERS = {'H':1, 'He':2,
#'Li':3, 'Be':4, 'B':5, 'C':6,  'N':7,  'O':8, 'F':9, 'Ne':10,
#'Na':11, 'Mg':12, 'Al':13, 'Si':14,  'P':15,  'S':16, 'Cl':17, 'Ar':18,
#'K':19, 'Ca':20,
#'Sc':21, 'Ti':22, 'V':23, 'Cr':24, 'Mn':25, 'Fe':26, 'Co':27, 'Ni':28, 'Cu':29, 'Zn':30,
#'Ga':31, 'Ge':32, 'As':33, 'Se':34, 'Br':35, 'Kr':36,
#'Rb':37, 'Sr':38,
#'Y':39,  'Zr':40, 'Nb':41, 'Mo':42, 'Tc':43, 'Ru':44, 'Rh':45, 'Pd':46, 'Ag':47, 'Cd':48,
#'In':49, 'Sn':50, 'Sb':51, 'Te':52,  'I':53, 'Xe':54,
#'Cs':55, 'Ba':56,
#'La':57, 
#'Ce':58, 'Pr':59, 'Nd':60, 'Pm':61, 'Sm':62, 'Eu':63, 'Gd':64, 'Tb':65, 'Dy':66, 'Ho':67, 'Er':68, 'Tm':69, 'Yb':70, 'Lu':71,
         #'Hf':72, 'Ta':73,  'W':74, 'Re':75, 'Os':76, 'Ir':77, 'Pt':78, 'Au':79, 'Hg':80,
#'Tl':81, 'Pb':82, 'Bi':83, 'Po':84, 'At':85, 'Rn':86, 
#'Fr':87, 'Ra':88,
#'Ac':89, 
#'Th':90, 'Pa':91,  'U':92, 'Np':93, 'Pu':94, 'Am':95, 'Cm':96, 'Bk':97, 'Cf':98, 'Es':99,'Fm':100,'Md':101,'No':102,'Lr':103,
        #'Rf':104,'Db':105,'Sg':106,'Bh':107,'Hs':108,'Mt':109,'Ds':110,'Rg':111,'Cn':112,
#'Nh':113,'Fl':114,'Mc':115,'Lv':116,'Ts':117,'Og':118
#}


# Atomic Weights of the most common isotopes
# From https://chemistry.sciences.ncsu.edu/msf/pdf/IsotopicMass_NaturalAbundance.pdf
#MASSES = {'H' :   1.007825 * U_TO_AMU,
          #'He':   4.002603 * U_TO_AMU,
          #'Li':   7.016004 * U_TO_AMU,
          #'Be':   9.012182 * U_TO_AMU,
          #'B' :  11.009305 * U_TO_AMU,
          #'C' :  12.000000 * U_TO_AMU,
          #'N' :  14.003074 * U_TO_AMU,
          #'O' :  15.994915 * U_TO_AMU,
          #'F' :  18.998403 * U_TO_AMU,
          #'Ne':  19.992440 * U_TO_AMU,
          #'Na':  22.989770 * U_TO_AMU,
          #'Mg':  23.985042 * U_TO_AMU,
          #'Al':  26.981538 * U_TO_AMU,
          #'Si':  27.976927 * U_TO_AMU,
          #'P' :  30.973762 * U_TO_AMU,
          #'S' :  31.972071 * U_TO_AMU,
          #'Cl':  34.968853 * U_TO_AMU,
          #'Ar':  39.962383 * U_TO_AMU,
          #'K' :  38.963707 * U_TO_AMU,
          #'Ca':  39.962591 * U_TO_AMU,
          #'Sc':  44.955910 * U_TO_AMU,
          #'Ti':  47.947947 * U_TO_AMU,
          #'V' :  50.943964 * U_TO_AMU,
          #'Cr':  51.940512 * U_TO_AMU,
          #'Mn':  54.938050 * U_TO_AMU,
          #'Fe':  55.934942 * U_TO_AMU,
          #'Co':  58.933200 * U_TO_AMU,
          #'Ni':  57.935348 * U_TO_AMU,
          #'Cu':  62.929601 * U_TO_AMU,
          #'Zn':  63.929147 * U_TO_AMU,
          #'Ga':  68.925581 * U_TO_AMU,
          #'Ge':  73.921178 * U_TO_AMU,
          #'As':  74.921596 * U_TO_AMU,
          #'Se':  79.916522 * U_TO_AMU,
          #'Br':  78.918338 * U_TO_AMU,
          #'Kr':  83.911507 * U_TO_AMU,
          #'Rb':  84.911789 * U_TO_AMU,
          #'Sr':  87.905614 * U_TO_AMU,
          #'Y' :  88.905848 * U_TO_AMU,
          #'Zr':  89.904704 * U_TO_AMU,
          #'Nb':  92.906378 * U_TO_AMU,
          #'Mo':  97.905408 * U_TO_AMU,
          #'Tc':  98.907216 * U_TO_AMU,
          #'Ru': 101.904350 * U_TO_AMU,
          #'Rh': 102.905504 * U_TO_AMU,
          #'Pd': 105.903483 * U_TO_AMU,
          #'Ag': 106.905093 * U_TO_AMU,
          #'Cd': 113.903358 * U_TO_AMU,
          #'In': 114.903878 * U_TO_AMU,
          #'Sn': 119.902197 * U_TO_AMU,
          #'Sb': 120.903818 * U_TO_AMU,
          #'Te': 129.906223 * U_TO_AMU,
          #'I' : 126.904468 * U_TO_AMU,
          #'Xe': 131.904154 * U_TO_AMU,
          #'Cs': 132.905447 * U_TO_AMU,
          #'Ba': 137.905241 * U_TO_AMU,
          #'La': 138.906348 * U_TO_AMU,
          #'Ce': 139.905435 * U_TO_AMU,
          #'Pr': 140.907648 * U_TO_AMU,
          #'Nd': 141.907719 * U_TO_AMU,
          #'Pm': 144.912744 * U_TO_AMU,
          #'Sm': 151.919729 * U_TO_AMU,
          #'Eu': 152.921227 * U_TO_AMU,
          #'Gd': 157.924101 * U_TO_AMU,
          #'Tb': 158.925343 * U_TO_AMU,
          #'Dy': 163.929171 * U_TO_AMU,
          #'Ho': 164.930319 * U_TO_AMU,
          #'Er': 165.930290 * U_TO_AMU,
          #'Tm': 168.934211 * U_TO_AMU,
          #'Yb': 173.938858 * U_TO_AMU,
          #'Lu': 174.940768 * U_TO_AMU,
          #'Hf': 179.946549 * U_TO_AMU,
          #'Ta': 180.947996 * U_TO_AMU,
          #'W' : 183.950933 * U_TO_AMU,
          #'Re': 186.955751 * U_TO_AMU,
          #'Os': 191.961479 * U_TO_AMU,
          #'Ir': 192.962924 * U_TO_AMU,
          #'Pt': 194.964774 * U_TO_AMU,
          #'Au': 196.966552 * U_TO_AMU,
          #'Hg': 201.970626 * U_TO_AMU,
          #'Tl': 204.974412 * U_TO_AMU,
          #'Pb': 207.976636 * U_TO_AMU,
          #'Bi': 208.980383 * U_TO_AMU,
          #'Po': 208.982416 * U_TO_AMU,
          #'At': 209.987131 * U_TO_AMU,
          #'Rn': 222.017570 * U_TO_AMU,
          #'Fr': 223.019731 * U_TO_AMU, 
          #'Ra': 226.025403 * U_TO_AMU,
          #'Ac': 227.027747 * U_TO_AMU, 
          #'Th': 232.038050 * U_TO_AMU, 
          #'Pa': 231.035879 * U_TO_AMU, 
          #'U' : 238.050783 * U_TO_AMU, 
          #'Np': 237.048167 * U_TO_AMU,
          #'Pu': 244.064198 * U_TO_AMU, 
          #'Am': 243.061373 * U_TO_AMU, 
          #'Cm': 247.070347 * U_TO_AMU, 
          #'Bk': 247.070299 * U_TO_AMU, 
          #'Cf': 251.079580 * U_TO_AMU, 
          #'Es': 252.082972 * U_TO_AMU,
          #'Fm': 257.095099 * U_TO_AMU,
          #'Md': 258.098425 * U_TO_AMU,
          #'No': 259.101024 * U_TO_AMU,
          #'Lr': 262.109692 * U_TO_AMU,
          #'Rf': 267. * U_TO_AMU,
          #'Db': 268. * U_TO_AMU,
          #'Sg': 269. * U_TO_AMU,
          #'Bh': 270. * U_TO_AMU,
          #'Hs': 270. * U_TO_AMU,
          #'Mt': 278. * U_TO_AMU,
          #'Ds': 281. * U_TO_AMU,
          #'Rg': 282. * U_TO_AMU,
          #'Cn': 285. * U_TO_AMU,
          #'Nh': 286. * U_TO_AMU,
          #'Fl': 289. * U_TO_AMU,
          #'Mc': 290. * U_TO_AMU,
          #'Lv': 293. * U_TO_AMU,
          #'Ts': 294. * U_TO_AMU,
          #'Og': 294. * U_TO_AMU
#}


# Isotopes used for the masses
#ISOTOPES={'H' : 'H-1' ,
          #'He': 'He-4',
          #'Li': 'Li-7',
          #'Be': 'Be-9*',
          #'B' : 'B_11' ,
          #'C' : 'C-12' ,
          #'N' : 'N-14' ,
          #'O' : 'O-16' ,
          #'F' : 'F-19*' ,
          #'Ne': 'Ne-20',
          #'Na': 'Na-23*',
          #'Mg': 'Mg-24',
          #'Al': 'Al-27*',
          #'Si': 'Si-28',
          #'P' : 'P-31*' ,
          #'S' : 'S-32' ,
          #'Cl': 'Cl-35',
          #'Ar': 'Ar-40',
          #'K' : 'K-39' ,
          #'Ca': 'Ca-40',
          #'Sc': 'Sc-45*',
          #'Ti': 'Ti-48',
          #'V' : 'V-51' ,
          #'Cr': 'Cr-52',
          #'Mn': 'Mn-55*',
          #'Fe': 'Fe-56',
          #'Co': 'Co-59*',
          #'Ni': 'Ni-58',
          #'Cu': 'Cu-63',
          #'Zn': 'Zn-64',
          #'Ga': 'Ga-69',
          #'Ge': 'Ge-74',
          #'As': 'As-75*',
          #'Se': 'Se-80',
          #'Br': 'Br-79',
          #'Kr': 'Kr-84',
          #'Rb': 'Rb-85',
          #'Sr': 'Sr-88',
          #'Y' : 'Y-89*' ,
          #'Zr': 'Zr-90',
          #'Nb': 'Nb-93*',
          #'Mo': 'Mo-98',
          #'Tc': 'Tc-98',
          #'Ru': 'Ru-102',
          #'Rh': 'Rh-103*',
          #'Pd': 'Pd-106',
          #'Ag': 'Ag-107',
          #'Cd': 'Cd-114',
          #'In': 'In-115',
          #'Sn': 'Sn-120',
          #'Sb': 'Sb-121',
          #'Te': 'Te-130',
          #'I' : 'I-127*' ,
          #'Xe': 'Xe-132',
          #'Cs': 'Cs-133*',
          #'Ba': 'Ba-138',
          #'La': 'La-139',
          #'Ce': 'Ce-140',
          #'Pr': 'Pr-141*',
          #'Nd': 'Nd-142',
          #'Pm': 'Pm-145',
          #'Sm': 'Sm-152',
          #'Eu': 'Eu-153',
          #'Gd': 'Gd-158',
          #'Tb': 'Tb-159*',
          #'Dy': 'Dy-164',
          #'Ho': 'Ho-165*',
          #'Er': 'Er-166',
          #'Tm': 'Tm-169*',
          #'Yb': 'Yb-174',
          #'Lu': 'Lu-175',
          #'Hf': 'Hf-180',
          #'Ta': 'Ta-181',
          #'W' : 'W-184' ,
          #'Re': 'Re-187',
          #'Os': 'Os-192',
          #'Ir': 'Ir-193',
          #'Pt': 'Pt-195',
          #'Au': 'Au-197*',
          #'Hg': 'Hg-202',
          #'Tl': 'Tl-205',
          #'Pb': 'Pb-208',
          #'Bi': 'Bi-209*',
          #'Po': 'Po-209',
          #'At': 'At-210',
          #'Rn': 'Rn-222',
          #'Fr': 'Fr-223', 
          #'Ra': 'Ra-226',
          #'Ac': 'Ac-227', 
          #'Th': 'Th-232*', 
          #'Pa': 'Pa-231*', 
          #'U' : 'U-238' , 
          #'Np': 'Np-237',
          #'Pu': 'Pu-244', 
          #'Am': 'Am-243', 
          #'Cm': 'Cm-247', 
          #'Bk': 'Bk-247', 
          #'Cf': 'Cf-251', 
          #'Es': 'Es-252',
          #'Fm': 'Fm-257',
          #'Md': 'Md-258',
          #'No': 'No-259',
          #'Lr': 'Lr-262',
          #'Rf': 'Rf-267',
          #'Db': 'Db-268',
          #'Sg': 'Sg-269',
          #'Bh': 'Bh-270',
          #'Hs': 'Hs-270',
          #'Mt': 'Mt-278',
          #'Ds': 'Ds-281',
          #'Rg': 'Rg-282',
          #'Cn': 'Cn-285',
          #'Nh': 'Nh-286',
          #'Fl': 'Fl-289',
          #'Mc': 'Mc-290',
          #'Lv': 'Lv-293',
          #'Ts': 'Ts-294',
          #'Og': 'Og-294'
#}


# ======================================================================= #

#import readline
#import re

##from compatibility import *


def centerstring(string, n, pad=' '):
  l = len(string)
  if l >= n:
    return string
  else:
    return  pad * ((n - l + 1) / 2) + string + pad * ((n - l) / 2)


def displaywelcome():
  print 'Script for setup of displacements started...\n'
  string='\n'
  string+='  '+'='*80+'\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Compute LVC parameters',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Author: Simon Kropf and Sebastian Mai',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Version:'+version,80)+'||\n'
  string+='||'+centerstring(versiondate.strftime("%d.%m.%y"),80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  string+='This script automatizes the setup of excited-state calculations for displacements\nfor SHARC dynamics.'
  print string


#def question(question,typefunc,default=None,autocomplete=True,ranges=False):
  #KEYSTROKES=open('KEYSTROKES.tmp','a')
  #if typefunc==int or typefunc==float:
    #if not default==None and not isinstance(default,list):
      #print 'Default to int or float question must be list!'
      #quit(1)
  #if typefunc==str and autocomplete:
    #readline.set_completer_delims(' \t\n;')
    #readline.parse_and_bind("tab: complete")    # activate autocomplete
  #else:
    #readline.parse_and_bind("tab: ")            # deactivate autocomplete

  #while True:
    #s=question
    #if default!=None:
      #if typefunc==bool or typefunc==str:
        #s+= ' [%s]' % (str(default))
      #elif typefunc==int or typefunc==float:
        #s+= ' ['
        #for i in default:
          #s+=str(i)+' '
        #s=s[:-1]+']'
    #if typefunc==str and autocomplete:
      #s+=' (autocomplete enabled)'
    #if typefunc==int and ranges:
      #s+=' (range comprehension enabled)'
    #s+=' '

    #line=raw_input(s)
    #line=re.sub('#.*$','',line).strip()
    #if not typefunc==str:
      #line=line.lower()

    #if line=='' or line=='\n':
      #if default!=None:
        #KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
        #return default
      #else:
        #continue

    #if typefunc==bool:
      #posresponse=['y','yes','true', 't', 'ja',  'si','yea','yeah','aye','sure','definitely']
      #negresponse=['n','no', 'false', 'f', 'nein', 'nope']
      #if line in posresponse:
        #KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
        #return True
      #elif line in negresponse:
        #KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
        #return False
      #else:
        #print 'I didn''t understand you.'
        #continue

    #if typefunc==str:
      #KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
      #return line

    #if typefunc==float:
      ## float will be returned as a list
      #f=line.split()
      #try:
        #for i in range(len(f)):
          #f[i]=typefunc(f[i])
        #KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
        #return f
      #except ValueError:
        #print 'Please enter floats!'
        #continue

    #if typefunc==int:
      ## int will be returned as a list
      #f=line.split()
      #out=[]
      #try:
        #for i in f:
          #if ranges and '~' in i:
            #q=i.split('~')
            #for j in range(int(q[0]),int(q[1])+1):
              #out.append(j)
          #else:
            #out.append(int(i))
        #KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
        #return out
      #except ValueError:
        #if ranges:
          #print 'Please enter integers or ranges of integers (e.g. "-3~-1  2  5~7")!'
        #else:
          #print 'Please enter integers!'
        #continue
  
  #KEYSTROKES.close()


#def print_INFOS(INFOS):
  #print '\n' + centerstring('Full input', 60, '#') + '\n'
  
  #for item in INFOS:
    #i = 0 # counter for very long lists we do not want full output)

    #if isinstance(INFOS[item], list) or isinstance(INFOS[item], tuple):
      #first = True
      #for elem in INFOS[item]:
        #if i >= 10: break
        #if first:
          #print item, ' ' * (25 - len(item) - 1), elem
          #first = False
        #else:
          #print ' ' * 25, elem

        #i += 1


    #elif isinstance(INFOS[item], dict):
      #first = True
      #for k, v in INFOS[item].items():
        #if i >= 10: break
        #if first:
          #print item, ' ' * (25 - len(item)) + '%s: %s' % (k, v)
          #first = False
        #else:
          #print ' ' * 25 + ' %s: %s' % (k, v)

        #i += 1
    #else: print item, ' ' * (25 - len(item) - 1), INFOS[item]

    #if i >= 10:
      #print ' ' * 25, '.'
      #print ' ' * 25, '(' + str(len(INFOS[item]) - i) + ' more)'
      #print ' ' * 25, '.'
  #return



#def reduce_big_list_to_short_str(big_list):
  #'''
  #Takes possibly big lists of numbers (e. g. list of all normal modes)
  #and reduces them to short string for clean output to user

  #e. g.: [7 8 9 12 13 14 17 20 21 22 23 25 26 28 29 30] => '[7~9 12~14 17 20~23 25 26 28~30]'

  #returns shortened list string
  #'''

  ## if empty return [None]
  #if not big_list: return [None]

  ## start list_string, sort bit_list, set vars
  #short_list_str = '('
  #big_list = sorted(big_list)
  #i_start, i_current, i_lastadded = 0, 0, 0
  
  ## while index is within big_list
  #while i_current < len(big_list) - 1:
    ## check if next element is within range & continue
    #if big_list[i_current] + 1 == big_list[i_current + 1]:
      #i_current += 1
      #continue
    ## range ended - create shortened string
    #else:
      ## no range just one number alone
      #if i_current == i_start:
        #short_list_str += '%i ' % big_list[i_current]
      ## range detected - shorten it
      #else:
        ## special case for 2 neighbour numbers
        #if big_list[i_start] + 1 == big_list[i_current]:
          #short_list_str += '%i %i ' % (big_list[i_start], big_list[i_current])
        #else:
          #short_list_str += '%i~%i ' % (big_list[i_start], big_list[i_current])

      ## set vars accordingly for next run
      #i_current += 1
      #i_start = i_current
      #i_lastadded = i_current

  ## code above will always leave out last range/number - add it here (that's why we need i_lastadded) 
  #if i_lastadded != i_current:
    ## special case again for 2 neighbouring numbers
    #if big_list[i_start] + 1 == big_list[i_current]:
      #short_list_str += '%i %i' % (big_list[i_start], big_list[i_current])
    #else:
      #short_list_str += '%i~%i' % (big_list[i_start], big_list[i_current])
  #else: short_list_str += str(big_list[i_current])
  
  ## close bracket & return
  #short_list_str += ')'
  #return short_list_str


#def reduce_displacement_dictionary_to_output_str(big_dictionary):
  #'''
  #used for shortening displacement dict therefore can never be empty

  #e. g.
  #for:
      #OrderedDict({ 7: 0.05,  8: 0.05,  9: 0.05, 12: 0.05, 13: 0.05,
                   #14: 0.04, 15: 0.05, 16: 0.05, 19: 0.03, 21: 0.03,
                   #22: 0.03, 23: 0.03, 27: 0.01 })

  #returns:
      #[7~9 12 13 15 16]: 0.05
      #[14]: 0.04
      #[19 21~23]: 0.03
      #[27]: 0.01

  #returns string with ranges for same displacements
  #'''
  #output_str = ''

  ## getting all different displacements
  #displacement_list = []
  #for normal_mode, displacement in big_dictionary.items():
    #if not displacement in displacement_list: displacement_list.append(displacement)

  ## running through all different displacements
  #normal_mode_list_big = []
  #for displacement in displacement_list:
    ## adding all normal modes with same displacement to list
    #for normal_mode, disp in big_dictionary.items():
      #if displacement == disp:
        #normal_mode_list_big.append(normal_mode)

    ## reduce all normal modes with same displacement and add it to nice output
    #output_str += '%s: %g\n' % (reduce_big_list_to_short_str(normal_mode_list_big), displacement)
    #normal_mode_list_big = []

  #return output_str

# ======================================================================= #
def itnmstates(states):
  '''Takes an array of the number of states in each multiplicity and generates an iterator over all states specified. Iterates also over all MS values of all states.

  Example:
  [3,0,3] yields 12 iterations with
  1,1,0
  1,2,0
  1,3,0
  3,1,-1
  3,2,-1
  3,3,-1
  3,1,0
  3,2,0
  3,3,0
  3,1,1
  3,2,1
  3,3,1

  Arguments:
  1 list of integers: States specification

  Returns:
  1 integer: multiplicity
  2 integer: state
  3 float: MS value'''

  for i in range(len(states)):
    if states[i]<1:
      continue
    for k in range(i+1):
      for j in range(states[i]):
        yield i+1,j+1,k-i/2.
  return


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
def read_QMout(path,nstates,natom,request):
  targets={'h':         {'flag': 1,
                         'type': complex,
                         'dim':  (nstates,nstates)},
           'dm':        {'flag': 2,
                         'type': complex,
                         'dim':  (3,nstates,nstates)},
           'grad':      {'flag': 3,
                         'type': float,
                         'dim':  (nstates,natom,3)},
           'nacdr':      {'flag': 5,
                         'type': float,
                         'dim':  (nstates,nstates,natom,3)},
            'overlap':   {'flag': 6,
                          'type': complex,
                          'dim':  (nstates,nstates)}
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
          print 'Could not find "%s" (flag "%i") in file %s!' % (t,targets[t]['flag'],path)
          sys.exit(11)
        line=lines[iline]
        if '! %i' % (targets[t]['flag']) in line:
          break
      values=[]
      # =========== single matrix
      if len(targets[t]['dim'])==2:
        iline+=1
        for irow in range(targets[t]['dim'][0]):
          iline+=1
          line=lines[iline].split()
          if targets[t]['type']==complex:
            row=[ complex(float(line[2*i]),float(line[2*i+1])) for i in range(targets[t]['dim'][1]) ]
          elif targets[t]['type']==float:
            row=[ float(line[i]) for i in range(targets[t]['dim'][1]) ]
          values.append(row)
      # =========== list of matrices
      elif len(targets[t]['dim'])==3:
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
      # =========== matrix of matrices
      elif len(targets[t]['dim'])==4:
        for iblocks in range(targets[t]['dim'][0]):
          sblock=[]
          for jblocks in range(targets[t]['dim'][1]):
            iline+=1
            block=[]
            for irow in range(targets[t]['dim'][2]):
              iline+=1
              line=lines[iline].split()
              if targets[t]['type']==complex:
                row=[ complex(float(line[2*i]),float(line[2*i+1])) for i in range(targets[t]['dim'][3]) ]
              elif targets[t]['type']==float:
                row=[ float(line[i]) for i in range(targets[t]['dim'][3]) ]
              block.append(row)
            sblock.append(block)
          values.append(sblock)
      QMout[t]=values

  #pprint.pprint(QMout)
  return QMout



# ======================================================================= #

def LVC_complex_mat(header, mat, deldiag=False, oformat=' % .7e'):
    rnonzero = False
    inonzero = False

    rstr = header + ' R\n'
    istr = header + ' I\n'
    for i in range(len(mat)):
        for j in range(len(mat)):
            val = mat[i][j].real
            if deldiag and i==j:
                val = 0.
            rstr += oformat%val
            if val*val > pthresh: rnonzero = True

            val = mat[i][j].imag
            if deldiag and i==j:
                val = 0.
            istr += oformat%val
            if val*val > pthresh: inonzero = True

        rstr += '\n'
        istr += '\n'

    retstr = ''
    if rnonzero: retstr += rstr
    if inonzero: retstr += istr

    return retstr





































# ======================================================================= #


def loewdin_orthonormalization(A):
  '''
  returns loewdin orthonormalized matrix
  '''

  # S = A^T * A
  S = np.dot(A.T, A)
  
  # S^d = U^T * S * U
  S_diag_only, U = np.linalg.eigh(S)
  
  # calculate the inverse sqrt of the diagonal matrix
  S_diag_only_inverse_sqrt = [1. / (float(d) ** 0.5) for d in S_diag_only]
  S_diag_inverse_sqrt = np.diag(S_diag_only_inverse_sqrt)
  
  # calculate inverse sqrt of S
  S_inverse_sqrt = np.dot(np.dot(U, S_diag_inverse_sqrt), U.T)
  
  # calculate loewdin orthonormalized matrix
  A_lo = np.dot(A, S_inverse_sqrt)

  # normalize A_lo
  A_lo = A_lo.T
  length = len(A_lo)
  A_lon = np.zeros((length, length), dtype = np.complex)

  for i in range(length):
    norm_of_col = np.linalg.norm(A_lo[i])
    A_lon[i] = [e / (norm_of_col ** 0.5) for e in A_lo[i]][0]

  return A_lon.T

# ======================================================================= #

def partition_matrix(matrix, multiplicity, states):
  '''
  return the first partitioned matrix of the given multiplicity

  e. g.: (3 0 2) states

    [111, 121, 131,   0,   0,   0,   0,   0,   0]       returns for multiplicity of 1:
    [112, 122, 132,   0,   0,   0,   0,   0,   0]             [111, 121, 131]
    [113, 123, 133,   0,   0,   0,   0,   0,   0]             [112, 122, 132]
    [  0,   0,   0, 311, 321,   0,   0,   0,   0]             [113, 123, 133]
    [  0,   0,   0, 312, 322,   0,   0,   0,   0] ====>
    [  0,   0,   0,   0,   0, 311, 321,   0,   0]       returns for multiplicity of 3:
    [  0,   0,   0,   0,   0, 312, 322,   0,   0]               [311, 321]
    [  0,   0,   0,   0,   0,   0,   0, 311, 321]               [312, 322]
    [  0,   0,   0,   0,   0,   0,   0, 312, 322]

    123 ^= 1...multiplicity
           2...istate
           3...jstate
  '''
  # get start index based on given multiplicity
  start_index = 0
  for i, state in enumerate(states):
    if (i + 1) == multiplicity: break
    else: start_index += state

  # size of the partition ^= state for given multiplicity
  size = states[multiplicity - 1]

  # create empty partition
  partition = np.zeros((size, size), dtype = complex)

  # get the partition out of the matrix
  for i in range(start_index, start_index + size):
    for j in range(start_index, start_index + size):
      partition[i - start_index][j - start_index] = matrix[i][j]

  return partition

# ======================================================================= #

def phase_correction(matrix):
  length = len(matrix)
  phase_corrected_matrix = [[.0 for x in range(length)] for x in range(length)]

  for i in range(length):
    diag = matrix[i][i].real

    # look if diag is significant and negative & switch phase
    if diag ** 2 > 0.5 and diag < 0:
      for j in range(length):
        phase_corrected_matrix[j][i] = matrix[j][i] * -1
    # otherwise leave values as is
    else:
      for j in range(length):
        phase_corrected_matrix[j][i] = matrix[j][i]

  return phase_corrected_matrix

# ======================================================================= #

def check_overlap_diagonal(matrix, states, normal_mode, displacement, ignore_problematic_states):
  '''
  Checks for problematic states (diagonals**2 of overlap matrix smaller than 0.5)
  '''
  problematic_states = {}

  for imult in range(len(states)):
    part_matrix = partition_matrix(matrix, imult + 1, states)

    for state in range(len(part_matrix)):
      sum_column = sum([part_matrix[j][state] ** 2 for j in range(len(part_matrix))])
      if sum_column < 0.5:
        print '* Problematic state %i in %i%s: %s' % (state + 1, int(normal_mode), displacement, IToMult[imult + 1])
        problematic_states[str(normal_mode) + displacement] = imult + 1

  return problematic_states

# ======================================================================= #

def calculate_W_dQi(H, S, e_ref, normal_mode, displ):
  '''
  Calculates the displacement matrix
  '''

  # get diagonalised hamiltonian
  H = np.diag([e - e_ref for e in np.diag(H)])

  # do phase correction if necessary
  if any([x for x in np.diag(S) if x < 0]):
    S = phase_correction(S)

  # do loewdin orthonorm. on overlap matrix
  U = loewdin_orthonormalization(np.matrix(S))

  return np.dot(np.dot(U.T, H), U)

# ======================================================================= #


def write_LVC_template(INFOS):
  lvc_template_content = '%s\n' % (INFOS['v0f'])
  lvc_template_content += str(INFOS['states'])[1:-1].replace(',', '') + '\n'

  #print INFOS

  # print some infos
  print '\nData extraction started ...'
  print 'Number of states:',INFOS['nstates']
  print 'Number of atoms:',len(INFOS['atoms'])
  print 'Kappas:', ['numerical','analytical'][INFOS['ana_grad']]
  print 'Lambdas:', ['numerical','analytical'][INFOS['ana_nac']]
  print
  print 'Reading files ...'
  print 

  # extract data from central point
  requests=['h','dm']
  if INFOS['ana_grad']:
    requests.append('grad')
  if INFOS['ana_nac']:
    requests.append('nacdr')
  path=os.path.join(INFOS['paths']['0eq'] , 'QM.out')
  print path, requests
  QMout_eq = read_QMout(path , INFOS['nstates'], len(INFOS['atoms']), requests)

  # ------------------ epsilon ----------------------
  epsilon_str_list = []

  i = 0
  e_ref = QMout_eq['h'][0][0]

  # run through all multiplicities
  for imult in range(len(INFOS['states'])):
    # partition matrix for every multiplicity
    partition = partition_matrix(QMout_eq['h'], imult + 1, INFOS['states'])

    # run over diagonal and get epsilon values
    for istate in range(len(partition)):
      epsilon_str_list.append('%3i %3i % .10f\n' % (imult + 1, istate + 1, (partition[istate][istate] - e_ref).real))

  # add results to template string
  lvc_template_content += 'epsilon\n'
  lvc_template_content += '%i\n' % (len(epsilon_str_list))
  lvc_template_content += ''.join(sorted(epsilon_str_list))

  # ------------------- kappa -----------------------
  nkappa = 0
  kappa_str_list = []
  r3N = [i for i in range(3 * len(INFOS['atoms'])) ]

  # run through all possible states
  if INFOS['ana_grad']:
    for i, sti in enumerate(itnmstates(INFOS['states'])):
      imult, istate, ims = sti

      if ims == (imult - 1) / 2.:
        # puts the gradient matrix into a list, has form: [ x, y, z, x, y, z, x, y, z]
        gradient = list(itertools.chain(*QMout_eq['grad'][i]))

        # runs through normal modes
        for normal_mode in INFOS['fmw_normal_modes'].keys():

          # calculates kappa from normal modes and grad
          kappa = sum([INFOS['fmw_normal_modes'][normal_mode][ixyz] * gradient[ixyz] for ixyz in r3N])

          # writes kappa to result string
          if kappa ** 2 > pthresh:
            kappa_str_list.append('%3i %3i %5i % .5e\n' % (imult, istate, int(normal_mode), kappa))
            nkappa += 1

  # ------------------------ lambda --------------------------
  lam = 0
  nlambda = 0
  lambda_str_list = []

  if INFOS['ana_nac']:

    for i, sti in enumerate(itnmstates(INFOS['states'])):
      imult, istate, ims = sti

      if ims != (imult - 1) / 2.:
        continue

      for j, stj in enumerate(itnmstates(INFOS['states'])):
        jmult, jstate, jms = stj

        if jms != (jmult - 1) / 2.:
          continue

        if i>=j:
          continue

        if imult!=jmult:
          continue

        if ims!=jms:
          continue

        nacvector = list(itertools.chain(*QMout_eq['nacdr'][i][j]))

        # runs through normal modes
        for normal_mode in INFOS['fmw_normal_modes'].keys():

          # calculates lambd from normal modes and grad
          dE = (QMout_eq['h'][j][j]-QMout_eq['h'][i][i]).real
          lambd = sum([INFOS['fmw_normal_modes'][normal_mode][ixyz] * nacvector[ixyz] for ixyz in r3N]) * dE

          # writes lambd to result string
          if lambd ** 2 > pthresh:
            #lambda_str_list.append('%3i %3i %5i % .5e\n' % (imult, istate, int(normal_mode), lambd))
            lambda_str_list.append('%3i %3i %3i %3i % .5e\n' % (imult, istate, jstate, int(normal_mode), lambd))
            nlambda += 1


  # ------------------------ numerical kappas and lambdas --------------------------

  if not ( INFOS['ana_nac'] and INFOS['ana_grad'] ):
    if not 'displacements' in INFOS:
      print 'No displacement info found in "displacements.json"!'
      sys.exit(1)

    if not INFOS['ana_nac'] and not INFOS['ana_grad']:
      whatstring='kappas and lambdas'
    elif not INFOS['ana_grad']:
      whatstring='kappas'
    elif not INFOS['ana_nac']:
      whatstring='lambdas'

    # running through all normal modes
    for normal_mode, v in INFOS['normal_modes'].items():

      twosided=False

      # get pos displacement
      pos_displ_mag = INFOS['displacement_magnitudes'][normal_mode]

      # get hamiltonian & overlap matrix from QM.out
      path=os.path.join(INFOS['paths'][str(normal_mode) + 'p'] , 'QM.out')
      requests=['h', 'overlap']
      print path, requests
      pos_H, pos_S = read_QMout(path , INFOS['nstates'], len(INFOS['atoms']), requests).values()

      # check diagonal of S & print warning
      INFOS['problematic_mults'] = check_overlap_diagonal(pos_S, INFOS['states'], normal_mode, 'p', INFOS['ignore_problematic_states'])

      # calculate displacement matrix
      pos_W_dQi = calculate_W_dQi(pos_H, pos_S, e_ref, normal_mode, 'p')


      # Check for two-sided differentiation
      if str(normal_mode) + 'n' in INFOS['displacements']:
        twosided=True
        # get neg displacement
        neg_displ_mag = INFOS['displacement_magnitudes'][normal_mode]

        # get hamiltonian & overlap matrix from QM.out
        path=os.path.join(INFOS['paths'][str(normal_mode) + 'n'] , 'QM.out')
        requests=['h', 'overlap']
        print path, requests
        neg_H, neg_S = read_QMout(path , INFOS['nstates'], len(INFOS['atoms']), requests).values()

        # check diagonal of S & print warning if wanted
        INFOS['problematic_mults'].update(check_overlap_diagonal(neg_S, INFOS['states'], normal_mode, 'n', INFOS['ignore_problematic_states']))

        # calculate displacement matrix
        neg_W_dQi = calculate_W_dQi(neg_H, neg_S, e_ref, normal_mode, 'n')


      # Loop over multiplicities to get kappas and lambdas
      for imult in range(len(INFOS['states'])):

        # checking problematic states
        if INFOS['ignore_problematic_states']:
          if str(normal_mode) + 'p' in INFOS['problematic_mults']:
            if INFOS['problematic_mults'][str(normal_mode) + 'p'] == imult + 1:
              print 'Not producing %s for normal mode: %s' % (whatstring,normal_mode)
              continue
          if str(normal_mode) + 'n' in INFOS['problematic_mults']:
            if twosided and INFOS['problematic_mults'][str(normal_mode) + 'n'] == imult + 1:
              print '! Not producing %s for multiplicity %i for normal mode: %s' % (whatstring,imult+1,normal_mode)
              continue

        # partition matrices
        pos_partition = partition_matrix(pos_W_dQi, imult + 1, INFOS['states'])
        if twosided:
          neg_partition = partition_matrix(neg_W_dQi, imult + 1, INFOS['states'])
        partition_length = len(pos_partition)

        # get lambdas and kappas
        for i in range(partition_length):
          if not INFOS['ana_grad']:
            if not twosided:
              kappa = pos_partition[i][i].real / pos_displ_mag
            else:
              kappa =  (pos_partition[i][i] - neg_partition[i][i]).real / (pos_displ_mag + neg_displ_mag)
            if kappa ** 2 > pthresh:
              kappa_str_list.append('%3i %3i %5i % .5e\n' % (imult+1, i+1, int(normal_mode), kappa))
              nkappa += 1

          if not INFOS['ana_nac']:
            for j in range(partition_length):
              if i>=j:
                continue
              if not twosided:
                lam = pos_partition[i][j].real / pos_displ_mag
              else:
                lam = (pos_partition[i][j] - neg_partition[i][j]).real / (pos_displ_mag + neg_displ_mag)
              if lam ** 2 > pthresh:
                lambda_str_list.append('%3i %3i %3i %3i % .5e\n' % (imult + 1, i + 1, j + 1, int(normal_mode), lam))
                nlambda += 1






  # add results to template string
  lvc_template_content += 'kappa\n'
  lvc_template_content += '%i\n' % (nkappa)
  lvc_template_content += ''.join(sorted(kappa_str_list))

  lvc_template_content += 'lambda\n'
  lvc_template_content += '%i\n' % (nlambda)
  lvc_template_content += ''.join(sorted(lambda_str_list))


  # ----------------------- matrices ------------------------------
  lvc_template_content += LVC_complex_mat('SOC', QMout_eq['h'], deldiag=True)
  lvc_template_content += LVC_complex_mat('DMX', QMout_eq['dm'][0])
  lvc_template_content += LVC_complex_mat('DMY', QMout_eq['dm'][1])
  lvc_template_content += LVC_complex_mat('DMZ', QMout_eq['dm'][2])


  # -------------------- write to file ----------------------------
  print '\nFinished!\nLVC parameters written to file: LVC.template\n'
  lvc_template = open('LVC.template', 'w')
  lvc_template.write(lvc_template_content)
  lvc_template.close()



# ======================================================================= #
# ======================================================================= #
# ======================================================================= #

def main():
  '''Main routine'''
  script_name = sys.argv[0].split('/')[-1]
  
  usage='''python %s''' % (script_name)

  parser = OptionParser(usage = usage, description = '')

  displaywelcome()

  # load INFOS object from file
  displacement_info_filename = 'displacements.json'
  try:
    with open(displacement_info_filename, 'r') as displacement_info:
      INFOS = json_load_byteified(displacement_info)
      displacement_info.close()
  except IOError:
    print 'IOError during opening readable %s - file. Quitting.' % (displacement_info_filename)
    quit(1)

  # set manually for old calcs
  # INFOS['ignore_problematic_states'] = True

  # write LVC.template
  write_LVC_template(INFOS)


# ======================================================================= #


if __name__ == '__main__':
  try:
    main()
  except KeyboardInterrupt:
    print '\nCtrl+C occured. Exiting.\n'
    sys.exit()
