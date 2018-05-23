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
# usage python wigner.py [-n <NUMBER>] <MOLDEN-FILE>

import copy
import math
import cmath
import random
import sys
import datetime
from optparse import OptionParser
import re

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

np = True
try:
    import numpy
except ImportError:
    print 'numpy package not installed'
    np = False



# some constants
DEBUG = False
CM_TO_HARTREE = 1./219474.6     #4.556335252e-6 # conversion factor from cm-1 to Hartree
HARTREE_TO_EV = 27.211396132    # conversion factor from Hartree to eV
U_TO_AMU = 1./5.4857990943e-4            # conversion from g/mol to amu
ANG_TO_BOHR = 1./0.529177211    #1.889725989      # conversion from Angstrom to bohr
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



# thresholds
LOW_FREQ = 10.0 # threshold in cm^-1 for ignoring rotational and translational low frequencies

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
def ask_for_masses():
  print '''
Option -m used, please enter non-default masses:
+ number mass           add non-default mass <mass> for atom <number>
- number                remove non-default mass for atom <number> (default mass will be used)
show                    show non-default atom masses
end                     finish input for non-default masses
'''
  MASS_LIST={}
  while True:
    line=raw_input()
    if 'end' in line:
      break
    if 'show' in line:
      s='-----------------------\nAtom               Mass\n'
      for i in MASS_LIST:
        s+='% 4i %18.12f\n' % (i,MASS_LIST[i])
      s+='-----------------------'
      print s
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
      MASS_LIST[num]=mass*U_TO_AMU
      continue
    if '-' in line:
      f=line.split()
      if len(f)<2:
        continue
      try:
        num=int(f[1])
      except ValueError:
        continue
      del MASS_LIST[num]
      continue
  return MASS_LIST


# ======================================================================================================================

def get_mass(symb,number):
  if 'MASS_LIST' in globals() and number in MASS_LIST:
    return MASS_LIST[number]
  else:
    try:
      return MASSES[symb]
    except KeyError:
      print 'No default mass for atom %s' % (symb)
      quit(1)


# ======================================================================================================================

def import_from_molden(filename,scaling,flag,lvc=False):
  '''This function imports atomic coordinates and normal modes from a MOLDEN
file. Returns molecule and modes as the other function does.
'''
  f=open(filename)
  data=f.readlines()
  f.close()

  # find coordinate block
  iline=0
  while not 'FR-COORD' in data[iline]:
    iline+=1
    if iline==len(data):
      print 'Could not find coordinates in %s!' % (filename)
      quit(1)
  # get atoms
  iline+=1
  natom=0
  molecule=[]
  while not '[' in data[iline]:
    f=data[iline].split()
    symb=f[0].lower().title()
    num=NUMBERS[symb]
    coord=[ float(f[i+1]) for i in range(3) ]
    natom+=1
    mass=get_mass(symb,natom)
    whichatoms.append(symb)
    molecule.append(ATOM(symb,num,coord,mass))
    iline+=1

  # find number of frequencies
  iline=-1
  nmodes=-1
  while True:
    iline+=1
    if iline==len(data):
      nmodes=3*natom
      break
    line=data[iline]
    if 'N_FREQ' in line:
      nmodes=int(data[iline+1])
      break

  # warn, if too few normal modes were found
  if nmodes<3*natom:
    print '*'*51+'\nWARNING: Less than 3*N_atom normal modes extracted!\n'+'*'*51+'\n'

  # obtain all frequencies, including low ones
  iline=0
  modes=[]
  while not '[FREQ]' in data[iline]:
    iline+=1
  iline+=1
  for imode in range(nmodes):
    try:
      mode={'freq':float(data[iline+imode])*CM_TO_HARTREE * scaling}
      modes.append(mode)
    except ValueError:
      print '*'*51+'\nWARNING: Less than 3*N_atom normal modes, but no [N_FREQ] keyword!\n'+'*'*51+'\n'
      nmodes=imode
      break

  # obtain normal coordinates
  iline=0
  while not 'FR-NORM-COORD' in data[iline]:
    iline+=1
  iline+=1
  for imode in range(nmodes):
    iline+=1
    move=[]
    for iatom in range(natom):
      f=data[iline].split()
      move.append([ float(f[i]) for i in range(3) ])
      iline+=1
    modes[imode]['move']=move
    # normalization stuff
    norm = 0.0
    for j, atom in enumerate(molecule):
      for xyz in range(3):
        norm += modes[imode]['move'][j][xyz]**2
    norm = math.sqrt(norm)
    if norm ==0.0 and modes[imode]['freq']>=LOW_FREQ*CM_TO_HARTREE:
      print 'WARNING: Displacement vector of mode %i is null vector. Ignoring this mode!' % (imode+1)
      modes[imode]['freq']=0.

  if not lvc:
    # delete low modes
    newmodes=[]
    for imode in range(nmodes):
      if modes[imode]['freq']<0.:
        print 'Detected negative frequency!'
      if modes[imode]['freq']>=LOW_FREQ*CM_TO_HARTREE:
        newmodes.append(modes[imode])
    modes=newmodes

  nmodes = len(modes)
  modes = determine_normal_modes_format(modes,molecule,nmodes,flag)

  return molecule, modes

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================


#def factorial(n,start):
    #"""This function calculates the factorial of n."""
    #p=1.
    #for i in range(n):
        #if i >= start:
            #p*=i+1
    #return p

#def ana_laguerre(n, x):
    #"""This function analytically calculates the value of the nth order
#Laguerre polynomial at point x. Computational limitations restrict very high 
#excited vibrational states (above 170). However, the harmonic approximation 
#is no good approximation at these high-lying vibrational states!"""
    #total = 0.
    #for m in range(n+1):
      #entry = (-1.)**m*factorial(n,m)/(factorial(n-m,1)*(factorial(m,1)))*(x**m)
      #total += entry
    #return total

def facfac_loop(n):
  yield 0, 1.
  r=1.
  for m in range(1,n+1):
    r*=float(n-m+1)/m**2
    yield m,r
  return

def ana_laguerre(n, x):
    """This function analytically calculates the value of the nth order
Laguerre polynomial at point x. Computational limitations restrict very high 
excited vibrational states (above 170). However, the harmonic approximation 
is no good approximation at these high-lying vibrational states!"""
    total = 0.
    for m,r in facfac_loop(n):
      entry = (-1.)**m * r * x**m
      total += entry
    return total



def determine_state(mode):
  """This function determines the vibrational state of the
mode for the system at a certain temperature."""
  # every state has a finite probability of being populated
  # at a finite temperature. we restrict only to so many
  # states that the sum of populations is "thresh"
  thresh=0.9999
  # also consider that with higher excited states harmonic
  # approximation is probably worse.
  freq = mode['freq']/CM_TO_HARTREE
  exponent = freq/(0.695035*temperature)  # factor for conversion cm-1 to K
  if exponent > 800:
    exponent = 600
    print '''The partition function is too close to zero due to very low temperature or very high frequency! It was set to %e''' % (math.exp(-exponent/2.) /\
 ( 1. - math.exp(-exponent) ))
  partition_function = math.exp(-exponent/2.) / \
                       ( 1. - math.exp(-exponent) )
  n=-1
  sum_p=0.
  prob=[]
  # calculate probabilities until sum is larger than threshold
  while True:
    n += 1
    p = math.exp(-exponent*(n+1./2.))/partition_function
    prob.append(p)
    sum_p += prob[n]
    if sum_p >= thresh:
      break
  n=-1
  probability = 0.
  # generate random number that is smaller than threshold
  while True:
    random_state = random.random()
    if random_state < sum_p:
      break
  # determine state number by comparing with random number
  while True:
    n += 1
    probability += prob[n]
    if probability >= random_state:
      return n
      break


def wigner(Q, P, mode):
    """This function calculates the Wigner distribution for
a single one-dimensional harmonic oscillator.
Q contains the dimensionless coordinate of the
oscillator and P contains the corresponding momentum.
n is the number of the vibrational state (default 0).
The function returns a probability for this set of parameters."""
    if temperature == 0:
      n = 0
    else:
      n = determine_state(mode)
      #square of the factorial becomes to large to handle. Keep in mind,
      #that the harmonic approximation is most likely not valid at these
      #excited states
      if n > 500:
        if high_temp:
          n = -1
          print 'Highest considered vibrational state reached! Discarding this probability.'
        else:
          print 'The calculated excited vibrational state for this normal mode exceeds the limit of the calculation.\nThe harmonic approximation is not valid for high vibrational states of low-frequency normal modes. The vibrational state ',n,' was set to 150. If you want to discard these states instead (due to oversampling of state nr 150), use the -T option.'
          n = 500
    if n == 0: # vibrational ground state
        return (math.exp(-Q**2) * math.exp(-P**2), 0.)
    # TODO: what about n==-1 ??
    else: # vibrational excited state
        rhosquare = 2.0 * (P**2 + Q**2)
        W = (-1.0)**n * ana_laguerre(n,rhosquare) \
          * math.exp(-rhosquare/2.0)
        n=float(n)
        return (W,n)

def plot_wigner_functions():
    """This function creates several output files for the first 11
Wigner functions for Q and P in the interval [-3, +3].
It also plots the wigner distribution for Q from -6 to +6
for the first 11 vibrational states."""
    grid = [0.06*i for i in range(-100, 101)]
    outstring = ''
    for Q in grid:
        outstring += '%6.2f' % Q
        for n in range(11):
            # integrate over P
            W_tot = 0.0
            for P in grid:
                W_tot += wigner(Q, P, n)
            outstring += ' % 12.8e' % W_tot
        outstring += '\n'
    outfile = open('wignerplot.out', 'w')
    outfile.write(outstring)
    outfile.close()

    grid = [0.06*i for i in range(-100, 101)]
    outstring = ''
    for Q in grid:
        outstring += '%6.2f' % Q
        for n in range(11):
            outstring += ' % 12.8e' % ana_laguerre(n, Q)
        outstring += '\n'
    outfile = open('laguerreplot.out', 'w')
    outfile.write(outstring)
    outfile.close()

    grid = [0.06*i for i in range(-100, 101)]
    for n in range(11):
        outstring = ''
        for Q in grid:
            for P in grid:
                W = wigner(Q, P, n)
                outstring += '%6.2f %6.2f % 12.8e\n' % (Q, P, W)
            outstring += '\n'
        filename = 'wigner%02i.out' % n
        outfile = open(filename, 'w')
        outfile.write(outstring)
        outfile.close()

def get_center_of_mass(molecule):
    """This function returns a list containing the center of mass
of a molecule."""
    mass = 0.0
    for atom in molecule:
        mass += atom.mass
    com = [0.0 for xyz in range(3)]
    for atom in molecule:
        for xyz in range(3):
            com[xyz] += atom.coord[xyz] * atom.mass / mass
    return com

def restore_center_of_mass(molecule, ic):
    """This function restores the center of mass for the distorted
geometry of an initial condition."""
    # calculate original center of mass
    com = get_center_of_mass(molecule)
    # caluclate center of mass for initial condition of molecule
    com_distorted = get_center_of_mass(ic)
    # get difference vector and restore original center of mass
    diff = [com[xyz] - com_distorted[xyz] for xyz in range(3)]
    for atom in ic:
        for xyz in range(3):
            atom.coord[xyz] += diff[xyz]

def remove_translations(ic):
    """This function calculates the movement of the center of mass
of an initial condition for a small timestep and removes this vector
from the initial condition's velocities."""
    # get center of mass at t = 0.0
    com = get_center_of_mass(ic)
    # get center of mass at t = dt = 0.01
    ic2 = copy.deepcopy(ic)
    dt = 0.01
    for atom in ic2:
        for xyz in range(3):
            atom.coord[xyz] += dt*atom.veloc[xyz]
    com2 = get_center_of_mass(ic2)
    # calculate velocity of center of mass and remove it
    v_com = [ (com2[xyz]-com[xyz])/dt for xyz in range(3) ]
    for atom in ic:
        for xyz in range(3):
            atom.veloc[xyz] -= v_com[xyz]
    if DEBUG:
        # check if v_com now is really zero
        # get center of mass at t = 0.0
        com = get_center_of_mass(ic)
        # get center of mass at t = dt = 1.0
        ic2 = copy.deepcopy(ic)
        dt = 1.0
        for atom in ic2:
            for xyz in range(3):
                atom.coord[xyz] += dt*atom.veloc[xyz]
        com2 = get_center_of_mass(ic2)
        # calculate velocity of center of mass and remove it
        v_com = [ (com2[xyz]-com[xyz])/dt for xyz in range(3) ]
        print v_com


def det(m):
    """This function calculates the determinant of a 3x3 matrix."""
    return m[0][0]*m[1][1]*m[2][2] + m[0][1]*m[1][2]*m[2][0] \
         + m[0][2]*m[1][0]*m[2][1] - m[0][0]*m[1][2]*m[2][1] \
         - m[0][1]*m[1][0]*m[2][2] - m[0][2]*m[1][1]*m[2][0]

def inverted(m):
    """This function calculates the inverse of a 3x3 matrix."""
    norm = m[0][0] * (m[1][1]*m[2][2] - m[1][2]*m[2][1]) \
         + m[0][1] * (m[1][2]*m[2][0] - m[1][0]*m[2][2]) \
         + m[0][2] * (m[1][0]*m[2][1] - m[1][1]*m[2][0])
    m_inv = [[0.0 for i in range(3)] for j in range(3)]
    m_inv[0][0] = (m[1][1]*m[2][2] - m[1][2]*m[2][1]) / norm
    m_inv[0][1] = (m[0][2]*m[2][1] - m[0][1]*m[2][2]) / norm
    m_inv[0][2] = (m[0][1]*m[1][2] - m[0][2]*m[1][1]) / norm
    m_inv[1][0] = (m[1][2]*m[2][0] - m[1][0]*m[2][2]) / norm
    m_inv[1][1] = (m[0][0]*m[2][2] - m[0][2]*m[2][0]) / norm
    m_inv[1][2] = (m[0][2]*m[1][0] - m[0][0]*m[1][2]) / norm
    m_inv[2][0] = (m[1][0]*m[2][1] - m[1][1]*m[2][0]) / norm
    m_inv[2][2] = (m[0][1]*m[2][0] - m[0][0]*m[2][1]) / norm
    m_inv[2][2] = (m[0][0]*m[1][1] - m[0][1]*m[1][0]) / norm
    return m_inv

def matmul(m1, m2):
    """This function multiplies two NxN matrices m1 and m2."""
    # get dimensions of resulting matrix
    n = len(m1)
    # calculate product
    result = [[0.0 for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                result[i][j] += m1[i][k]*m2[k][j]
    return result

def cross_prod(a, b):
    """This function calculates the cross product of two
3 dimensional vectors."""
    result = [0.0 for i in range(3)]
    result[0] = a[1]*b[2] - b[1]*a[2]
    result[1] = a[2]*b[0] - a[0]*b[2]
    result[2] = a[0]*b[1] - b[0]*a[1]
    return result

def linmapping(lm, y):
    z = [0.0 for i in range(3)]
    z[0] = lm[0][0]*y[0] + lm[0][1]*y[1] + lm[0][2]*y[2]
    z[1] = lm[1][0]*y[0] + lm[1][1]*y[1] + lm[1][2]*y[2]
    z[2] = lm[2][0]*y[0] + lm[2][1]*y[1] + lm[2][2]*y[2]
    return z

def remove_rotations(ic):
    # copy initial condition object
    ictmp = copy.deepcopy(ic)
    # move center of mass to coordinates (0, 0, 0)
    com = get_center_of_mass(ic)
    for atom in ictmp:
        for xyz in range(3):
            atom.coord[xyz] -= com[xyz]
    # calculate moment of inertia tensor
    I = [[0.0 for i in range(3)] for j in range(3)]
    for atom in ictmp:
        I[0][0] += atom.mass*(atom.coord[1]**2 + atom.coord[2]**2)
        I[1][1] += atom.mass*(atom.coord[0]**2 + atom.coord[2]**2)
        I[2][2] += atom.mass*(atom.coord[0]**2 + atom.coord[1]**2)
        I[0][1] -= atom.mass * atom.coord[0] * atom.coord[1]
        I[0][2] -= atom.mass * atom.coord[0] * atom.coord[2]
        I[1][2] -= atom.mass * atom.coord[1] * atom.coord[2]
    I[1][0] = I[0][1]
    I[2][0] = I[0][2]
    I[2][1] = I[1][2]
    if det(I) > 0.01: # checks if I is invertible
        ch = matmul(I, inverted(I))
        # calculate angular momentum
        ang_mom = [0.0 for i in range(3)]
        for atom in ictmp:
            mv = [0.0 for i in range(3)]
            for xyz in range(3):
                mv[xyz] = atom.mass * atom.veloc[xyz]
            L = cross_prod(mv, atom.coord)
            for xyz in range(3):
                ang_mom[xyz] -= L[xyz]
        # calculate angular velocity
        ang_vel = linmapping(inverted(I), ang_mom)
        for i,atom in enumerate(ictmp):
            v_rot = cross_prod(ang_vel, atom.coord) # calculate rotational velocity
            for xyz in range(3):
                ic[i].veloc[xyz] -= v_rot[xyz] # remove rotational velocity
    else:
        print 'WARNING: moment of inertia tensor is not invertible'

def constrain_displacement(molecule, ic, threshold=0.5):
    """This function ensures, that each atom of a generated initial
condition is not displaced further, than a given threshold from its
original position. Threshold is given in bohr."""
    for i, atom in enumerate(molecule):
        diff_vector = [ic[i].coord[xyz]-atom.coord[xyz]
                       for xyz in range(3)]
        displacement = 0.0
        for xyz in range(3):
            displacement += diff_vector[xyz]**2
        displacement = math.sqrt(displacement)
        if displacement > threshold:
            if DEBUG:
                print 'displacment for atom %i %s is %f' % (i, atom.symb, displacement)
            # shorten diff_vector to length of threshold
            for xyz in range(3):
                diff_vector[xyz] /= displacement/threshold
            # apply changes to initial condition
            for xyz in range(3):
                ic[i]['coords'][xyz] = atom.coord[xyz] \
                                       + diff_vector[xyz]

def determine_normal_modes_format(modes, molecule, nmodes, flag):
  '''This function determines the input format of the normal modes by trying to
transform them to mass-weighted coordinates and seeing which of the four methods
was able to do so via checking if the normal modes are now orthogonal. The mass-
weighted normal coordinates are then returned'''

  print '\nStarting normal mode format determination...'

  #generate different set of modes that each undergo a different transformation
  #modes_1, modes_2, modes_3 and modes are represented by the numbers 1, 2, 3
  #and 4, where 1 stands for gaussian-type coordinates, 2 for cartesian coordinates,
  #3 for Colombus-type coordinates and 4 for already mass-weighted coordinates.
  modes_1 = copy.deepcopy(modes)
  modes_2 = copy.deepcopy(modes)
  modes_3 = copy.deepcopy(modes)
  allmodes = [modes_1,modes_2,modes_3,modes]
  normformat = ["gaussian-type (Gaussian, Turbomole, Q-Chem, ADF, Orca)","cartesian (Molpro, Molcas)","columbus-type (Columbus)","mass-weighted"]

  #apply transformations to normal modes
  for imode in range(nmodes):
    norm = 0.0
    for j, atom in enumerate(molecule):
      for xyz in range(3):
        norm += modes_2[imode]['move'][j][xyz]**2*atom.mass/U_TO_AMU
    norm = math.sqrt(norm)
    if norm == 0.0 and modes[imode]['freq']>=LOW_FREQ*CM_TO_HARTREE:
      print 'WARNING: Displacement vector of mode %i is null vector. Ignoring this mode!' % (imode+1)
      for normmodes in allmodes:
        normmodes[imode]['freq']=0.0
    for j, atom in enumerate(molecule):
      for xyz in range(3):
        modes_1[imode]['move'][j][xyz] /= norm/math.sqrt(atom.mass/U_TO_AMU)
        modes_2[imode]['move'][j][xyz] *= math.sqrt(atom.mass/U_TO_AMU)
        modes_3[imode]['move'][j][xyz] *= math.sqrt(atom.mass/U_TO_AMU)/math.sqrt(ANG_TO_BOHR)
  if flag != 0:
    print "Using input flag",flag, "for", normformat[flag-1],"coordinates. Skipping normal mode analysis. "
    return allmodes[flag-1]

  elif int(flag) <= 4:
    #create dotproduct matrices of the normal mode multiplication
    #for all three transformations.
    if np:
      matrix = [[] for i in range(4)]
      for coord in range(len(molecule)):
        for xyz in range(3):
          displacement = [[] for i in range(4)]
          for mode in range(nmodes):
            for nr in range(len(allmodes)):
              displacement[nr].append(allmodes[nr][mode]['move'][coord][xyz])
          for nr in range(len(allmodes)):
            matrix[nr].append(displacement[nr])
      newmatrix = [[] for i in range(4)]
      results = [[] for i in range(4)]
      for nr in range(len(allmodes)):
        newmatrix[nr] = numpy.array(matrix[nr])
        results[nr] = numpy.dot(newmatrix[nr].T,newmatrix[nr])

    else:
      #do the slow matrix multiplication of every normal mode with every other
      #this approach is approximately 25 times slower than the numpy approach
      #at a 190 atom test-system but works always
      results = [[[0 for j in range(nmodes)]for i in range(nmodes)] for k in range(len(allmodes))]
      for mode1 in range(nmodes):
        for mode2 in range(nmodes):
          dotproduct = [0,0,0,0]
          for coord, atom in enumerate(molecule):
            for xyz in range(3):
              for i in range(len(allmodes)):
                dotproduct[i] += allmodes[i][mode1]['move'][coord][xyz]*allmodes[i][mode2]['move'][coord][xyz]
          for i in range(len(dotproduct)):
            results[i][mode1][mode2] = dotproduct[i]


    #check for orthogonal matrices
    diagonalcheck = [[],[]]
    thresh = 0.05
    for result in results:
      trace = 0
      for i in range(len(result)):
        trace += result[i][i]
        result[i][i] -= 1
      diagonalcheck[0].append(trace)
      #print all matrices
      #for row in result:
        #string = ''
        #for entry in row:
        #  string += "%4.1f" % (float(entry))
        #print string
      if any( [abs(i) > thresh for j in result for i in j ] ):
        diagonalcheck[1].append(0)
      else:
        diagonalcheck[1].append(1)
    possibleflags = []
    for i in range(4):
      if diagonalcheck[0][i] > nmodes-1 and diagonalcheck[0][i]/nmodes-1 < thresh and diagonalcheck[1][i] == 1:
        possibleflags.append(i+1)
        #this means that previous flag is overwritten if multiple checks were positive.
        #However ordering of the checks is made in a way that it gives the right result
        #for all QM programs tested so far.
        nm_flag = i
    #check for input flag
    try:
      print "Final format specifier: %s [%s]" % (nm_flag+1, normformat[nm_flag])
    except UnboundLocalError:
      print "The normal mode analysis was unable to diagonalize the normal modes."
      print "Input is therefore neither in cartesian, gaussian-type, Columbus-type, or mass weighted coordinates."
      exit(1)
    if len(possibleflags) != 1:
      string = '\n'
      for entry in possibleflags:
        string += '  %s \n' % (normformat[entry-1])
      print "Multiple possible flags have been identified: %s" % (string[:-2])
      print "The most likely assumption is %s coordinates."  % (normformat[nm_flag])
      print "These have been used in the creation of inital conditions."
      print "\nYou can override this behavior by setting the -f [int] flag in the command line:"
      string = ''
      for mode in range(len(normformat)):
         string += "  "+str(mode+1) + "\t" + (normformat[mode]) +"\n"
      print string
    else:
      print "The normal modes input format was determined to be %s coordinates." % (normformat[nm_flag])
    #return the set of transformed normal modes that resulted in an orthogonal matrix (mass-weighted)
    return allmodes[nm_flag]
  else:
    print "Wrong input, please specify a valid flag [0,1,2,3,4]!"
    quit(1)



def sample_initial_condition(molecule, modes):
  """This function samples a single initial condition from the
modes and atomic coordinates by the use of a Wigner distribution.
The first atomic dictionary in the molecule list contains also
additional information like kinetic energy and total harmonic
energy of the sampled initial condition.
Method is based on L. Sun, W. L. Hase J. Chem. Phys. 133, 044313
(2010) nonfixed energy, independent mode sampling."""
  # copy the molecule in equilibrium geometry
  atomlist = copy.deepcopy(molecule) # initialising initial condition object
  Epot = 0.0
  for atom in atomlist:
    atom.veloc = [0.0, 0.0, 0.0] # initialise velocity lists
  for mode in modes: # for each uncoupled harmonatomlist oscillator
    while True:
      # get random Q and P in the interval [-3,+3]
      # this interval is good for vibrational ground state
      # should be increased for higher states
      random_Q = random.random()*10.0 - 5.0
      random_P = random.random()*10.0 - 5.0
      # calculate probability for this set of P and Q with Wigner distr.
      probability = wigner(random_Q, random_P, mode)
      if probability[0]>1. or probability[0]<0.:
        if temperature == 0:
          print 'WARNING: wrong probability %f detected!' % (probability[0])
      elif probability[0] > random.random():
        break # coordinates accepted
    # now transform the dimensionless coordinate into a real one
    # paper says, that freq_factor is sqrt(2*PI*freq)
    # QM programs directly give angular frequency (2*PI is not needed)
    freq_factor = math.sqrt(mode['freq'])
    # Higher frequencies give lower displacements and higher momentum.
    # Therefore scale random_Q and random_P accordingly:
    random_Q /= freq_factor
    random_P *= freq_factor
    # add potential energy of this mode to total potential energy
    Epot += 0.5 * mode['freq']**2 * random_Q**2
    for i, atom in enumerate(atomlist): # for each atom
      for xyz in range(3): # and each direction
        # distort geometry according to normal mode movement
        # and unweigh mass-weighted normal modes
        if not UEG:
          atom.coord[xyz] += random_Q * mode['move'][i][xyz] * math.sqrt(1./atom.mass)
        # add velocity
        if not UZV:
          atom.veloc[xyz] += random_P * mode['move'][i][xyz] * math.sqrt(1./atom.mass)
      atom.EKIN()
  if not KTR:
    restore_center_of_mass(molecule, atomlist)
    remove_translations(atomlist)
    remove_rotations(atomlist)

  ic = INITCOND(atomlist,0.,Epot)
  return ic

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

#def initial_condition_to_string(ic):
    #"""This function converts an initial condition into a formatted
#string and returns it."""
    #outstring = 'Geometry of molecule (in bohr):\n'
    #for atom in ic.atomlist:
        #outstring += ' %2s %5.1f %12.8f %12.8f %12.8f %12.8f\n' \
                     #% (atom.symb, NUMBERS[atom.symb],
                        #atom.coord[0], atom.coord[1],
                        #atom.coord[2], atom.mass/U_TO_AMU)
    #outstring += 'Velocities of the single atoms (a.u.):\n'
    #for atom in ic.atomlist:
        #outstring += '%12.8f %12.8f %12.8f\n' % (atom.veloc[0],
                                #atom.veloc[1], atom.veloc[2])
    #for key in sorted(ic[0].keys(), key=str.lower):
        #if key not in ('veloc', 'symbol', 'coords', 'mass') and not key.startswith('Excited state'):
            #outstring += '%s: %12.8f a.u. ( %12.8f eV)\n' \
                          #% (key, ic[0][key], ic[0][key]*HARTREE_TO_EV)
        #elif key.startswith('Excited state'):
            #outstring += '%s:\n' % key
            #for exckey in sorted(ic[0][key].keys(), key=str.lower):
                #if exckey in ('E_final', 'E_exc'):
                    #outstring += '    %s: %12.8f a.u. ( %12.8f eV)\n' \
                     #% (exckey, ic[0][key][exckey], ic[0][key][exckey]*HARTREE_TO_EV)
                #elif exckey in ('Osc'):
                    #outstring += '    %s: %12.8f\n' % (exckey, ic[0][key][exckey])
                #elif exckey in ('Excitation'):
                    #outstring += '    %s: %s\n' % (exckey, ic[0][key][exckey])
                #else:
                    #print exckey, key, ic[0][key][exckey]
    #return outstring

# ======================================================================================================================

def create_initial_conditions_string(molecule, modes, ic_list, eref=0.0):
  """This function converts an list of initial conditions into a string."""
  ninit=len(ic_list)
  natom=ic_list[0].natom
  representation='None'
  #eref
  eharm=0.
  for mode in modes:
    eharm+=mode['freq']*0.5
  string='''SHARC Initial conditions file, version %s
Ninit     %i
Natom     %i
Repr      %s
Temp      %18.10f
Eref      %18.10f
Eharm     %18.10f

Equilibrium
''' % (version,ninit,natom,representation,temperature,eref,eharm)
  for atom in molecule:
    string+=str(atom)+'\n'
  string+='\n\n'

  for i, ic in enumerate(ic_list):
    string += 'Index     %i\n%s' % (i+1, str(ic))
  return string

# ======================================================================================================================


def create_initial_conditions_list(amount, molecule, modes):
    """This function creates 'amount' initial conditions from the
data given in 'molecule' and 'modes'. Output is returned
as a list containing all initial condition objects."""
    print 'Sampling initial conditions'
    ic_list = []
    width = 50
    idone = 0
    for i in range(1,amount+1): # for each requested initial condition
        # sample the initial condition
        ic = sample_initial_condition(molecule, modes)
        ic_list.append(ic)
        idone += 1
        done = idone*width/(amount)
        sys.stdout.write('\rProgress: ['+'='*done+' '*(width-done)+'] %3i%%' % (done*100/width))
        sys.stdout.flush()
    print '\n'
    return ic_list

# ======================================================================================================================

def make_dyn_file(ic_list,filename):
  #if not os.path.exists('init_geoms'):
    #os.mkdir('init_geoms')
  #for state in range(states):
  fl=open(filename,'w')
  string=''
  for i,ic in enumerate(ic_list):
    string+='%i\n%i\n' % (ic.natom,i)
    for atom in ic.atomlist:
      string+='%s' % (atom.symb)
      for j in range(3):
        string+=' %f' % (atom.coord[j]/ANG_TO_BOHR)
      string+='\n'
  fl.write(string)
  fl.close()

# ======================================================================================================================

def lvc_input(molecule, modes):
    """
    Print an input file for for SHARC_LVC.py
    """
    if not len(modes)==len(modes[0]['move']*3):
      print 'Error: Less than 3N normal modes in file!'
      print 'For LVC model setup, all 3N modes need to be present.'
      sys.exit(1)
    print "Creating V0.txt for SHARC_LVC.py ..."

    wf = open('V0.txt', 'w')
    wf.write('Geometry\n')
    for atom in molecule:
      wf.write(str(atom)[:61]+'\n')
    wf.write('Frequencies\n')
    for mode in modes:
      wf.write(' % .10f'%(mode['freq']))
    wf.write('\n')
    wf.write('Mass-weighted normal modes\n')
    for j, atom in enumerate(molecule):
      for xyz in range(3):
          for imode in range(len(modes)):
              wf.write(' % .10f'%modes[imode]['move'][j][xyz])
          wf.write('\n')
    wf.close()

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def main():
  '''Main routine'''

  usage='''
Wigner.py [options] filename.molden

This script reads a MOLDEN file containing frequencies and normal modes [1]
and generates a Wigner distribution of geometries and velocities.

The creation of the geometries and velocities is based on the
sampling of the Wigner distribution of a quantum harmonic oscillator,
as described in [2] (non-fixed energy, independent mode sampling).

[1] http://www.cmbi.ru.nl/molden/molden_format.html
[2] L. Sun, W. L. Hase: J. Chem. Phys. 133, 044313 (2010)
'''

  description=''

  parser = OptionParser(usage=usage, description=description)
  parser.add_option('-n', dest='n', type=int, nargs=1, default=3, help="Number of geometries to be generated (integer, default=3)")
  parser.add_option('-m', dest='m', action='store_true',help="Enter non-default atom masses")
  parser.add_option('-s', dest='s', type=float, nargs=1, default=1.0, help="Scaling factor for the energies (float, default=1.0)")
  parser.add_option('-t', dest='t', type=float, nargs=1, default=0., help="Temperature (float, default=0.0)")
  parser.add_option('-T', dest='T', action='store_true', help="Discard high vibrational states in the temperature sampling ")
  parser.add_option('-L', dest='L', type=float, nargs=1, default=10.0, help="Discard frequencies below this value in cm-1 (float, default=10.)")

  parser.add_option('-o', dest='o', type=str, nargs=1, default='initconds', help="Output filename (string, default=""initconds"")")
  parser.add_option('-x', dest='X', action='store_true',help="Generate a xyz file with the sampled geometries in addition to the initconds file")
  parser.add_option('-l', dest='lvc', action='store_true', help='Generate input for SHARC_LVC.py (V0.txt) rather than initconds')

  parser.add_option('-r', dest='r', type=int, nargs=1, default=16661, help="Seed for the random number generator (integer, default=16661)")
  parser.add_option('-f', dest='f', type=int, nargs=1, default='0', help="Define the type of read normal modes. 0 for automatic assignement, 1 for gaussian-type normal modes (Gaussian, Turbomole, Q-Chem, ADF, Orca), 2 for cartesian normal modes (Molcas, Molpro), 3 for Columbus-type (Columbus), or 4 for mass-weighted. (integer, default=0)")
  
  parser.add_option('--keep_trans_rot', dest='KTR', action='store_true',help="Keep translational and rotational components")
  parser.add_option('--use_eq_geom',    dest='UEG', action='store_true',help="For all samples, use the equilibrium geometry (only sample velocities)")
  parser.add_option('--use_zero_veloc', dest='UZV', action='store_true',help="For all samples, set velocities to zero")
  
  (options, args) = parser.parse_args()

  random.seed(options.r)
  amount=options.n
  if len(args)==0:
    print usage
    quit(1)
  filename=args[0]
  outfile=options.o
  nondefmass=options.m
  scaling=options.s
  flag=options.f
  lvc = options.lvc
  global LOW_FREQ
  LOW_FREQ=max(0.0000001,options.L)

  print '''Initial condition generation started...
INPUT  file                  = "%s"
OUTPUT file                  = "%s"
Number of geometries         = %i
Random number generator seed = %i
Temperature                  = %f''' % (filename, outfile, options.n, options.r,options.t)
  if nondefmass:
    global MASS_LIST
    MASS_LIST = ask_for_masses()
  else:
    print ''
  if scaling!=1.0:
    print 'Scaling factor               = %f\n' % (scaling)

  global KTR
  KTR=options.KTR
  global UEG
  UEG=options.UEG
  global UZV
  UZV=options.UZV

  global temperature
  temperature=options.t
  if temperature!=0:
    print 'Using temperature-dependent sampling'
  #   print 'Script only for finite temperature. Sorry...'
  #   sys.exit()

  global high_temp
  if options.T:
    high_temp = True
  else:
    high_temp = False


  global whichatoms
  whichatoms=[]

  molecule, modes = import_from_molden(filename, scaling,flag,options.lvc)

  string='\nGeometry:\n'
  for atom in molecule:
    string+=str(atom)[:61]+'\n'
  string+='Assumed Isotopes: '
  for i in set(whichatoms):
    string+=ISOTOPES[i]+' '
  string+='\nIsotopes with * are pure isotopes.\n'
  print string

  string='Frequencies (cm^-1) used in the calculation:\n'
  for i,mode in enumerate(modes):
    string+='%4i %12.4f\n' % (i+1,mode['freq']/CM_TO_HARTREE)
  print string

  if options.lvc:
      lvc_input(molecule, modes)
  else:
      #print 'Generating %i initial conditions' % amount
      ic_list = create_initial_conditions_list(amount, molecule, modes)
      #print 'Writing output to initconds'
      outfile = open(outfile, 'w')
      outstring = create_initial_conditions_string(molecule, modes, ic_list)
      outfile.write(outstring)
      outfile.close()

  if options.X and not options.lvc:
    make_dyn_file(ic_list,options.o+'.xyz')

  # save the shell command
  command='python '+' '.join(sys.argv)
  f=open('KEYSTROKES.wigner','w')
  f.write(command)
  f.close()

# ======================================================================================================================

if __name__ == '__main__':
    main()
