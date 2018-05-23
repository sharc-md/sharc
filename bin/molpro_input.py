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

import copy
import math
import sys
import re
import os
import stat
import shutil
import datetime
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

version='2.0'
versiondate=datetime.date(2018,2,1)


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
  string='\n'
  string+='  '+'='*80+'\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('MOLPRO Input file generator',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Author: Sebastian Mai',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Version:'+version,80)+'||\n'
  string+='||'+centerstring(versiondate.strftime("%d.%m.%y"),80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  string+='''
This script allows to quickly create MOLPRO input files for single-points calculations,
ground state optimizations, frequency calculations and SA-CASSCF calculations. 
It also generates MOLPRO.template files to be used with the SHARC-MOLPRO Interface.
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
  shutil.move('KEYSTROKES.tmp','KEYSTROKES.molpro_input')

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

def show_massses(masslist):
  s='Number\tType\tMass\n'
  for i,atom in enumerate(masslist):
    s+='%i\t%2s\t%12.9f %s\n' % (i+1,atom[0],atom[1], ['','*'][atom[1]!=MASSES[atom[0]]])
  print s

def ask_for_masses(masslist):
  print '''
Please enter non-default masses:
+ number mass           use non-default mass <mass> for atom <number>
- number                remove non-default mass for atom <number> (default mass will reinstated)
show                    show atom masses
end                     finish input for non-default masses
'''
  show_massses(masslist)
  while True:
    line=question('Change an atoms mass:',str,'end',False)
    if 'end' in line:
      break
    if 'show' in line:
      show_massses(masslist)
      continue
    if '+' in line:
      f=line.split()
      if len(f)<3:
        continue
      try:
        num=int(f[1])
        mass=float(f[2])
      except ValueError:
        continue
      if not 0<=num<=len(masslist):
        print 'Atom %i does not exist!' % (num)
        continue
      masslist[num-1][1]=mass
      continue
    if '-' in line:
      f=line.split()
      if len(f)<2:
        continue
      try:
        num=int(f[1])
      except ValueError:
        continue
      if not 0<=num<=len(masslist):
        print 'Atom %i does not exist!' % (num)
        continue
      masslist[num-1][1]=MASSES[masslist[num-1][0]]
      continue
  return masslist

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def get_infos():
  '''Asks for the settings of the calculation:
- type (single point, optimization+freq or MOLPRO.template
- level of theory
- basis set
- douglas kroll
- memory
- geometry

specific:
- opt: freq?
- CASSCF: docc, act'''

  INFOS={}

  # Type of calculation
  print centerstring('Type of calculation',60,'-')
  print '''\nThis script generates input for the following types of calculations:
  1       Single point calculations (HF, DFT, MP2, SS/SA-CASSCF, EOM-CCSD)
  2       Optimizations & Frequency calculations (HF, DFT, MP2, SS/SA-CASSCF)
  3       MOLPRO.template file for dynamics (SA-CASSCF)
  4       Crossing point optimization: CI and MXP (SA-CASSCF)
Please enter the number corresponding to the type of calculation.
'''
  while True:
    ctype=question('Type of calculation:',int)[0]
    if not ctype in [1,2,3,4]:
      print 'Enter an integer (1-4)!'
      continue
    break
  INFOS['ctype']=ctype
  freq=False
  if ctype==2:
    freq=question('Frequency calculation?',bool,True)
  INFOS['freq']=freq
  print ''


  guessnact=None
  guessnorb=None
  guessnelec=None
  guessbase=None
  guessstates=None
  guessmem=None


  # Geometry
  print centerstring('Geometry',60,'-')
  if ctype==3:
    print '\nNo geometry necessary for MOLPRO.template generation\n'
    INFOS['geom']=None
    # see whether a MOLPRO.input file is there, where we can take the number of electrons from
    nelec=0
    try:
      molproinput=open('MOLPRO.input','r')
      for line in molproinput:
        if 'wf,' in line and not './wf,' in line:
          guessnelec=[int(line.split(',')[1])]
          mult=int(line.split(',')[3])
        if 'state,' in line:
          if guessstates==None:
            guessstates=[]
          nstate=int(line.split(',')[1])
          for i in range(mult-len(guessstates)):
            guessstates.append(0)
          guessstates.append(nstate)
        if 'closed,' in line:
          nclosed=int(line.split(',')[1])
        if 'occ,' in line:
          nocc=int(line.split(',')[1])
        if 'basis=' in line:
          guessbase=line.split('=')[1].strip()
        if 'memory' in line:
          guessmem=[int(line.split(',')[1])/125]
      try:
        guessnorb=[nocc-nclosed]
        guessnact=[guessnelec[0]-2*nclosed]
      except:
        pass
    except (IOError,ValueError):
      pass
    # continue with asking for number of electrons
    while True:
      nelec=question('Number of electrons: ',int,guessnelec,False)[0]
      if nelec<=0:
        print 'Enter a positive number!'
        continue
      break
    INFOS['nelec']=nelec
  else:
    print '\nPlease specify the geometry file (xyz format, Angstroms):'
    while True:
      path=question('Geometry filename:',str,'geom.xyz')
      try:
        gf=open(path,'r')
      except IOError:
        print 'Could not open: %s' % (path)
        continue
      g=gf.readlines()
      gf.close()
      try:
        natom=int(g[0])
      except ValueError:
        print 'Malformatted: %s' % (path)
        continue
      geom=[]
      ncharge=0
      fine=True
      for i in range(natom):
        try:
          line=g[i+2].split()
        except IndexError:
          print 'Malformatted: %s' % (path)
          fine=False
        try:
          atom=[line[0],float(line[1]),float(line[2]),float(line[3])]
        except (IndexError,ValueError):
          print 'Malformatted: %s' % (path)
          fine=False
          continue
        geom.append(atom)
        try:
          ncharge+=NUMBERS[atom[0]]
        except KeyError:
          print 'Atom type %s not supported!' % (atom[0])
          fine=False
      if not fine:
        continue
      else:
        break
    print 'Number of atoms: %i\nNuclear charge: %i\n' % (natom,ncharge)
    INFOS['geom']=geom
    INFOS['ncharge']=ncharge
    INFOS['natom']=natom
    print 'Enter the total (net) molecular charge:'
    while True:
      charge=question('Charge:',int,[0])[0]
      break
    INFOS['nelec']=ncharge-charge
    print 'Number of electrons: %i\n' % (ncharge-charge)

  # Masses
  if INFOS['freq']:
    # make default mass list
    masslist=[]
    for atom in geom:
      masslist.append( [atom[0],MASSES[atom[0]]] )
    # ask
    INFOS['nondefmass']=not question('Use standard masses (most common isotope)?',bool,True)
    if INFOS['nondefmass']:
      INFOS['masslist']=ask_for_masses(masslist)
    else:
      INFOS['masslist']=masslist

  # Level of theory
  print '\n'+centerstring('Level of theory',60,'-')
  allowed=[1,2,3,4,5]
  s='''\nSupported by this script are:
  1       HF
  2       DFT %s
  3       MP2 %s
  4       SS-CASSCF
  5       SA-CASSCF %s
''' % tuple(3*[['','(Only numerical frequencies)'][INFOS['freq']]])
  if ctype==1 and INFOS['nelec']%2==0:
    s+='  6       EOM-CCSD\n'
    allowed.append(6)
  else:
    s+='EOM-CCSD is only possible single-point calculations of singlet (even-electron) states.\n'
  print s

  if ctype==3 or ctype==4:
    ltype=5
    print 'Choosing SA-CASSCF for MOLPRO.template generation.'
  else:
    while True:
      ltype=question('Level of theory:',int)[0]
      if not ltype in allowed:
        print 'Enter an integer in %s!' % (allowed)
        continue
      break
  INFOS['ltype']=ltype

  # DFT
  if ltype==2:
    print '''Commons functionals and their names in MOLPRO:
  B3LYP      B3LYP, B3LYP3, B3LYP5
  BP86       B-P
  PBE        PBE
  PBE0       PBE0
'''
    func=question('Functional:',str,None,False)
    INFOS['dft.func']=func
    disp=question('Dispersion correction? ',bool)
    INFOS['dft.disp']=disp

  # basis set
  print '\nPlease enter the basis set.'
  cadpac=(ctype==2 and ltype+freq>=5) or ctype==3 or ctype==4
  if (ctype==2 and ltype+freq>=5) or ctype==4:
    print 'For SA-CASSCF Optimizations/Frequencies and SS-CASSCF Frequencies,\nonly segmented basis sets are allowed.'
  if ctype==3:
    print 'For MOLPRO.template generation, only segmented basis sets are allowed.'
  print '''Common available basis sets:
  Pople:     6-31G**, 6-311G, 6-31+G, 6-31G(d,p), ...
  Dunning:   cc-pVXZ, aug-cc-pVXZ, cc-pVXZ-DK, ...    %s
  Turbomole: def2-SV(P), def2-SVP, def2-TZVP, ...
  ANO:       ROOS                                     %s''' % (['','not available'][cadpac],['','not available'][cadpac])
  basis=question('Basis set:',str,guessbase,False)
  INFOS['basis']=basis

  # douglas kroll
  dk=question('Douglas-Kroll scalar-relativistic integrals?',bool,True)
  INFOS['DK']=dk

  # CASSCF
  if ltype==4 or ltype==5:
    print '\n'+centerstring('CASSCF Settings',60,'-')+'\n'
    while True:
      nact=question('Number of active electrons:',int,guessnact)[0]
      if nact<=0:
        print 'Enter a positive number!'
        continue
      if (INFOS['nelec']-nact)%2!=0:
        print 'nelec-nact must be even!'
        continue
      if INFOS['nelec']<nact:
        print 'Number of active electrons cannot be larger than total number of electrons!'
        continue
      break
    INFOS['cas.nact']=nact
    while True:
      norb=question('Number of active orbitals:',int,guessnorb)[0]
      if norb<=0:
        print 'Enter a positive number!'
        continue
      if 2*norb<nact:
        print 'norb must be larger than nact/2!'
        continue
      break
    INFOS['cas.norb']=norb

  if ltype<5:
    print '\nPlease enter the multiplicity (1=singlet, 2=doublet, 3=triplet, ...)'
    while True:
      mult=question('Multiplicity:',int,[1])[0]
      if mult<=0:
        print 'Enter a positive number!'
        continue
      if (INFOS['nelec']-mult-1)%2!=0:
        print 'Nelec is %i, so mult cannot be %i' % (INFOS['nelec'],mult)
        continue
      break
    INFOS['mult']=mult
    if ltype==4:
      INFOS['cas.nstates']=[0 for i in range(mult)]
      INFOS['cas.nstates'][mult-1]=1
      INFOS['maxmult']=mult
  elif ltype==5:
    print 'Please enter the number of states as a list of integers\ne.g. 3 0 3 for three singlets, zero doublets and three triplets.'
    while True:
      states=question('Number of states:',int,guessstates)
      maxmult=len(states)
      for i in range(maxmult):
        n=states[i]
        if not ctype==3:
          if (not i%2==INFOS['nelec']%2) and int(n)>0:
            print 'Nelec is %i. Ignoring states with mult=%i!' % (INFOS['nelec'], i+1)
            states[i]=0
        if n<0:
          states[i]=0
      if sum(states)==0:
        print 'No states!'
        continue
      break
    s='Accepted number of states:'
    for i in states:
      s+=' %i' % (i)
    print s
    INFOS['maxmult']=len(states)
    INFOS['cas.nstates']=states
    if ctype==2:
      print '\nPlease specify the state to optimize\ne.g. 3 2 for the second triplet state.'
      while True:
        rmult,rstate=tuple(question('Root:',int,[1,1]))
        if not 1<=rmult<=INFOS['maxmult']:
          print '%i must be between 1 and %i!' % (rmult,INFOS['maxmult'])
          continue
        if not 1<=rstate<=states[rmult-1]:
          print 'Only %i states of mult %i' % (states[rmult-1],rmult)
          continue
        break
      INFOS['cas.root']=[rmult,rstate]
    if ctype==4:
      print '\nPlease specify the first state involved in the optimization\ne.g. 3 2 for the second triplet state.'
      while True:
        rmult,rstate=tuple(question('Root:',int,[1,1]))
        if not 1<=rmult<=INFOS['maxmult']:
          print '%i must be between 1 and %i!' % (rmult,INFOS['maxmult'])
          continue
        if not 1<=rstate<=states[rmult-1]:
          print 'Only %i states of mult %i' % (states[rmult-1],rmult)
          continue
        break
      INFOS['cas.root1']=[rmult,rstate]
      print '\nPlease specify the second state involved in the optimization\ne.g. 3 2 for the second triplet state.'
      while True:
        rmult,rstate=tuple(question('Root:',int,[1,2]))
        if not 1<=rmult<=INFOS['maxmult']:
          print '%i must be between 1 and %i!' % (rmult,INFOS['maxmult'])
          continue
        if not 1<=rstate<=states[rmult-1]:
          print 'Only %i states of mult %i' % (states[rmult-1],rmult)
          continue
        break
      INFOS['cas.root2']=[rmult,rstate]
      if INFOS['cas.root1']==INFOS['cas.root2']:
        print 'Both states are identical, please use calculation type 2 for optimizations of state minima.'
        quit(1)
      if INFOS['cas.root1'][0]==INFOS['cas.root2'][0]:
        print 'Multiplicities of both states identical, optimizing a conical intersection.'
        INFOS['cas.opt_ci']=True
      else:
        print 'Multiplicities of both states different, optimizing a minimum crossing point.'
        INFOS['cas.opt_ci']=False
    if ctype==1 and maxmult>1:
      INFOS['soci']=question('Do Spin-Orbit CASCI after CASSCF?',bool,False)
  elif ltype==6:
    print '\nPlease enter the number of singlet states (1=only CCSD, >1=EOM-CCSD)'
    while True:
      nstates=question('Number of states:',int,[1])[0]
      if nstates<=0:
        print 'Enter a positive number!'
        continue
      break
    if nstates==1:
      INFOS['ccsd.eom']=False
    else:
      INFOS['ccsd.eom']=True
    INFOS['ccsd.states']=nstates
    INFOS['mult']=1
    if nstates>1:
      INFOS['ccsd.trans']=question('Calculate oscillator strength?',bool,False)

  if not INFOS['ctype']==3:
    print '\n'+centerstring('Memory',60,'-')
    print '\nRecommendation: for small systems: 100-300 MB, for medium-sized systems: 1000-2000 MB\n'
    mem=abs(question('Memory in MB: ',int,guessmem)[0])
    mem=max(mem,50)
    INFOS['mem']=mem

  print ''

  return INFOS

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def setup_input(INFOS):
  ''''''

  if INFOS['ctype']==3:
    inpf='MOLPRO.template'
  else:
    inpf='MOLPRO.input'
  print 'Writing input to %s' % (inpf)
  try:
    inp=open(inpf,'w')
  except IOError:
    print 'Could not open %s for write!' % (inpf)
    quit(1)

# ======================================================================
  # go to routine for new interface template
  if INFOS['ctype']==3:
    s=get_template(INFOS)
    s+='\n\n---\n'
    s+='#Infos:\n'
    s+='#%s@%s\n' % (os.environ['USER'],os.environ['HOSTNAME'])
    s+='#Date: %s\n' % (datetime.datetime.now())
    s+='#Current directory: %s\n\n' % (os.getcwd())
    inp.write(s)
    return
# ======================================================================

  s='***,%s generated by molpro_input.py Version %s\n' % (inpf,version)
  s+='memory,%i,k\n\n' % (INFOS['mem']*125)     # convert to Mega-Words
  if INFOS['ctype']!=3:
    s+='file,1,./integrals,scratch\n'
    s+='file,2,./wf,new   ! remove ",new" if you want to restart\n'
  s+='\n\nprint,orbitals,civectors;\n\n'
  if INFOS['DK']:
    s+='dkroll=1\ndkho=2\n'
  s+='basis=%s\n\n' % (INFOS['basis'])

  if INFOS['geom']:
    s+='nosym\nangstrom\ngeometry={\n'
    for iatom,atom in enumerate(INFOS['geom']):
      s+='%s%i % 16.9f % 16.9f % 16.9f\n' % (atom[0],iatom+1,atom[1],atom[2],atom[3])
    s+='}\n\n'
  if INFOS['freq']:
    s+='mass,isotope\n'
    for iatom,atom in enumerate(INFOS['geom']):
      s+='mass,init,%s%i=%f\n' % (atom[0],iatom+1,INFOS['masslist'][iatom][1])
    s+='mass,print\n\n'

  # ============================ HF
  if INFOS['ltype']==1:
    if INFOS['nelec']%2==0:
      s+='{hf'
    else:
      s+='{uhf'
    if INFOS['mult']!=1 or INFOS['ncharge']!=INFOS['nelec']:
      s+='\nwf,%i,1,%i\n' % (INFOS['nelec'],INFOS['mult']-1)
    s+='};\n\n'
  # ============================ DFT
  elif INFOS['ltype']==2:
    if INFOS['nelec']%2==0:
      s+='{ks'
    else:
      s+='{uks'
    s+=',%s' % (INFOS['dft.func'])
    if INFOS['dft.disp']:
      s+=';disp'
    if INFOS['mult']!=1 or INFOS['ncharge']!=INFOS['nelec']:
      s+='\nwf,%i,1,%i\n' % (INFOS['nelec'],INFOS['mult']-1)
    s+='};\n\n'
  # ============================ MP2
  elif INFOS['ltype']==3:
    if INFOS['nelec']%2==0:
      s+='{hf'
    else:
      s+='{uhf'
    if INFOS['mult']!=1 or INFOS['ncharge']!=INFOS['nelec']:
      s+='\nwf,%i,1,%i\n' % (INFOS['nelec'],INFOS['mult']-1)
    if INFOS['nelec']%2==0:
      s+='};\n{mp2};\n\n'
    else:
      s+='};\n{ump2};\n\n'
  # ============================ CCSD
  elif INFOS['ltype']==6:
    if INFOS['nelec']%2==0:
      s+='{hf'
    else:
      s+='{uhf'
    if INFOS['mult']!=1 or INFOS['ncharge']!=INFOS['nelec']:
      s+='\nwf,%i,1,%i\n' % (INFOS['nelec'],INFOS['mult']-1)
    if INFOS['nelec']%2==0:
      s+='};\n{ccsd'
    else:
      s+='};\n{uccsd'
    if INFOS['ccsd.eom']:
      s+='\neom,-%i.1%s' % (INFOS['ccsd.states'],['',',trans=1'][INFOS['ccsd.trans']])
    s+='};\n\n'
  # ============================ CASSCF
  elif INFOS['ltype']==4 or INFOS['ltype']==5:
    s+='{casscf\n'
    s+='frozen,0\nclosed,%i\n' % ((INFOS['nelec']-INFOS['cas.nact'])/2)
    s+='occ,%i\n' % (INFOS['cas.norb']+(INFOS['nelec']-INFOS['cas.nact'])/2)
    if INFOS['ctype']<3:
      s+='!start,2140.2       ! uncomment if restarting\n'
      s+='orbital,2140.2\n'
    if INFOS['ctype']==1:
      s+='!rotate,-1.1,-1.1   ! uncomment if rotating orbitals\n'
    for i,n in enumerate(INFOS['cas.nstates']):
      if n==0:
        continue
      s+='wf,%i,1,%i\n' % (INFOS['nelec'],i)
      s+='state,%i\n' % (n)
      s+='weight'+',1'*n+'\n'

    if INFOS['ctype']==2:
      if INFOS['ltype']==5:
        s+='\ncpmcscf,grad,state=%i.1,ms2=%i,record=5001.1,accu=1e-7\n' % (INFOS['cas.root'][1],INFOS['cas.root'][0]-1)
      if INFOS['ltype']==4 and INFOS['freq']:
        s+='\ncpmcscf,hess,accu=1e-4\n'
    if INFOS['ctype']==4:
      if INFOS['cas.opt_ci']:
        if INFOS['cas.root1'][1]>INFOS['cas.root2'][1]:
          INFOS['cas.root1'],INFOS['cas.root2']=INFOS['cas.root2'],INFOS['cas.root1']
        s+='\ncpmcscf,nacm,state1=%i.1,state2=%i.1,ms2=%i,record=5001.1,accu=1e-7\n' % (
          INFOS['cas.root1'][1],
          INFOS['cas.root2'][1],
          INFOS['cas.root1'][0]-1)
      s+='\ncpmcscf,grad,state=%i.1,            ms2=%i,record=5002.1,accu=1e-7\n' % (INFOS['cas.root1'][1],INFOS['cas.root1'][0]-1)
      s+='\ncpmcscf,grad,state=%i.1,            ms2=%i,record=5003.1,accu=1e-7\n' % (INFOS['cas.root2'][1],INFOS['cas.root2'][0]-1)

    s+='};\n\n'

  if INFOS['ctype']==2:
    s+='{optg,maxit=50};\n'
    if INFOS['freq']:
      s+='{frequencies};\n'
    s+='\n'

  if INFOS['ctype']==4:
    if INFOS['cas.opt_ci']:
      recs=[1,2,3]
    else:
      recs=[2,3]
    for irec in recs:
      s+='{force\nsamc,%i.1\nconical,6100.1%s}\n\n' % (
        5000+irec,
        [',nodc',''][INFOS['cas.opt_ci']])
    s+='{optg,maxit=50,startcmd=casscf,gradient=1e-4};\n\n'
    s+='{casscf\n'
    s+='frozen,0\nclosed,%i\n' % ((INFOS['nelec']-INFOS['cas.nact'])/2)
    s+='occ,%i\n' % (INFOS['cas.norb']+(INFOS['nelec']-INFOS['cas.nact'])/2)
    for i,n in enumerate(INFOS['cas.nstates']):
      if n==0:
        continue
      s+='wf,%i,1,%i\n' % (INFOS['nelec'],i)
      s+='state,%i\n' % (n)
      s+='weight'+',1'*n+'\n'
    s+='};\n\n'

  if INFOS['ctype']==1:
    s+='PUT,MOLDEN,geom.molden\n'
  elif INFOS['ctype']==2 or INFOS['ctype']==4:
    if INFOS['freq']:
      s+='PUT,MOLDEN,freq.molden\n'
    else:
      s+='PUT,MOLDEN,opt.molden\n'

  if 'soci' in INFOS and INFOS['soci']:
    s+='\n\n'
    for i,n in enumerate(INFOS['cas.nstates']):
      if n==0:
        continue
      s+='{ci\nmaxiter,250,1000\norbital,2140.2\nsave,%i.2\nnoexc\ncore,%i\n' % (6001+i,(INFOS['nelec']-INFOS['cas.nact'])/2)
      s+='wf,%i,%i,%i\nstate,%i\n}\n\n' % (INFOS['nelec'],1,i,n)
    s+='{ci\nhlsmat,amfi'
    for i,n in enumerate(INFOS['cas.nstates']):
      if n==0:
        continue
      s+=',%i.2' % (6001+i)
    s+='\nprint,hls=1\n}\n\n'

  s+='\n\n---\n'
  s+='!Infos:\n'
  s+='!%s@%s\n' % (os.environ['USER'],os.environ['HOSTNAME'])
  s+='!Date: %s\n' % (datetime.datetime.now())
  s+='!Current directory: %s\n\n' % (os.getcwd())

  inp.write(s)

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def get_template(INFOS):
  s='basis %s\n' % (INFOS['basis'])
  if INFOS['DK']:
    s+='dkho 2\n'
  s+='''closed %i
occ %i
nelec %i
roots''' % ( (INFOS['nelec']-INFOS['cas.nact'])/2,
             (INFOS['nelec']-INFOS['cas.nact'])/2+INFOS['cas.norb'],
             INFOS['nelec'] )
  for i in INFOS['cas.nstates']:
    s+=' %i' % (i)
  s+='\nrootpad'
  for i in INFOS['cas.nstates']:
    s+=' 0'
  s+='\n'
  return s



# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def set_runscript(INFOS):

  if INFOS['ctype']==3:
    return

  print ''
  if not question('Runscript?',bool,True):
    return
  print ''

  # MOLPRO executable
  print centerstring('Path to MOLPRO',60,'-')+'\n'
  path=os.getenv('MOLPRO')
  path=os.path.expanduser(os.path.expandvars(path))
  if not path.endswith('/molpro'):
    path+='/molpro'
  if path!='':
    print 'Environment variable $MOLPRO detected:\n$MOLPRO=%s\n' % (path)
    if question('Do you want to use this MOLPRO installation?',bool,True):
      INFOS['molpro']=path
  if not 'molpro' in INFOS:
    print '\nPlease specify path to MOLPRO directory (SHELL variables and ~ can be used, will be expanded when interface is started).\n'
    INFOS['molpro']=question('Path to MOLPRO:',str)
  print ''


  # Scratch directory
  print centerstring('Scratch directory',60,'-')+'\n'
  print 'Please specify an appropriate scratch directory. This will be used to temporally store the integrals. The scratch directory will be deleted after the calculation. Remember that this script cannot check whether the path is valid, since you may run the calculations on a different machine. The path will not be expanded by this script.'
  INFOS['scratchdir']=question('Path to scratch directory:',str)+'/WORK'
  print ''

  runscript='run_MOLPRO.sh'
  print 'Writing run script %s' % (runscript)
  try:
    runf=open(runscript,'w')
  except IOError:
    print 'Could not write %s' (runscript)
    return

  string='''#!/bin/bash

PRIMARY_DIR=%s
SCRATCH_DIR=%s
cd $PRIMARY_DIR
mkdir -p $SCRATCH_DIR

%s MOLPRO.input -W$PRIMARY_DIR -I$SCRATCH_DIR -d$SCRATCH_DIR

rm -r $SCRATCH_DIR  ''' % (os.getcwd(), INFOS['scratchdir'], INFOS['molpro'])

  runf.write(string)
  runf.close()
  os.chmod(runscript, os.stat(runscript).st_mode | stat.S_IXUSR)

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def main():
  '''Main routine'''

  usage='''
python molpro_input.py

This interactive program prepares a MOLPRO input file for ground state optimizations and frequency calculations with HF, DFT, MP2 and CASSCF. It also generates input for SA-CASSCF excited-state calculations (MOLPRO.template files to be used with the SHARC-MOLPRO interface).
'''

  description=''
  parser = OptionParser(usage=usage, description=description)

  displaywelcome()
  open_keystrokes()

  INFOS=get_infos()

  print centerstring('Full input',60,'#')+'\n'
  for item in INFOS:
    print item, ' '*(15-len(item)), INFOS[item]
  print ''

  setup_input(INFOS)
  set_runscript(INFOS)
  print '\nFinished\n'

  close_keystrokes()

# ======================================================================================================================

if __name__ == '__main__':
  try:
    main()
  except KeyboardInterrupt:
    print '\nCtrl+C makes me a sad SHARC ;-(\n'
    quit(0)