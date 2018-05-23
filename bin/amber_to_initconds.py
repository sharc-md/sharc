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
import re
from optparse import OptionParser


# =========================================================
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


# =========================================================

# some constants
DEBUG = False
CM_TO_HARTREE = 1./219474.6     #4.556335252e-6 # conversion factor from cm-1 to Hartree
HARTREE_TO_EV = 27.211396132    # conversion factor from Hartree to eV
U_TO_AMU = 1./5.4857990943e-4            # conversion from g/mol to amu
ANG_TO_BOHR = 1./0.529177211    #1.889725989      # conversion from Angstrom to bohr
AMBERVEL_TO_AU = 0.0009350161
PI = math.pi

version='2.0'
versiondate=datetime.date(2018,2,1)

# =========================================================

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
    sys.exit(13)
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
      sys.exit(14)
    f.close()
  except IOError:
    print 'Could not write to file %s!' % (filename)
    sys.exit(15)

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

def restore_center_of_mass(ic):
    """This function restores the center of mass for the distorted
geometry of an initial condition."""
    # calculate original center of mass
    com = [0.0 for xyz in range(3)]
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
        atom.EKIN()
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

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def ask_for_masses():
  print '''
Option -m used, please enter non-default masses:
+ number mass           add non-default mass <mass> for atom <number> (counting starts at 1)
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

def read_mass_from_prmtop(filename,masses):
  data=readfile(filename)
  iline=-1
  while True:
    iline+=1
    if iline>=len(data):
      print 'Could not find masses in %s' % (filename)
      sys.exit(1)
    line=data[iline]
    if 'FLAG MASS' in line:
      break
  iline+=1
  iatom=0
  #masses={}
  while True:
    iline+=1
    line=data[iline]
    if 'FLAG' in line:
      break
    s=line.split()
    for i in s:
      iatom+=1
      if not iatom in masses:
        masses[iatom]=float(i)*U_TO_AMU
  return masses

# ======================================================================================================================

def get_mass(symb,number,MASSLIST):
  if number in MASSLIST:
    return MASSLIST[number]
  else:
    try:
      return MASSES[symb]
    except KeyError:
      print 'No default mass for atom %s' % (symb)
      sys.exit(1)

# ======================================================================================================================

def get_atoms_from_prmtop(filename):
  data=readfile(filename)
  iline=-1
  while True:
    iline+=1
    if iline>=len(data):
      print 'Could not find atoms in %s' % (filename)
      sys.exit(1)
    line=data[iline]
    if 'FLAG ATOMIC_NUMBER' in line:
      break
  atoms=[]
  iline+=1
  iatom=0
  while True:
    iline+=1
    line=data[iline]
    if 'FLAG' in line:
      break
    s=line.split()
    for i in s:
      iatom+=1
      n=int(i)
      if n<=0:
        el=raw_input('\nElement for atom number %i:\n' % (iatom))[0:2].title()
      else:
        for q in NUMBERS:
          if n==NUMBERS[q]:
            el=q
            break
      atoms.append(el)
      #print iatom,el
  print
  return atoms

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def get_coords(INFOS):

  # open the prmtop file
  MASSLIST=read_mass_from_prmtop(INFOS['filename_prmtop'], INFOS['masslist'])
  ATOMS=get_atoms_from_prmtop(INFOS['filename_prmtop'])
  natom=len(ATOMS)



  dt=INFOS['timestep']  # timestep of Amber in fs
  dt*=41.341373         # convert to atomic units
  dt*=0.5               # only take half a timestep to get correct geometry




  # initialize arrays
  ic_list=[]
  igeom=0

  # go through the data
  for filename in INFOS['filename_rsts']:
    data=readfile(filename)
    iline=2
    colcount=0
    atomlist=[]

    for iatom in range(natom):
      if colcount==6:
        iline+=1
        colcount=0
      line=data[iline]
      s=line.split()
      xyz=[ float(i)*ANG_TO_BOHR    for i in s[colcount:colcount+3] ]
      colcount+=3
      symb=ATOMS[iatom]
      num=NUMBERS[symb]
      vel=[0.,0.,0.]
      mass=get_mass(symb,iatom+1,MASSLIST)
      #print iatom,symb,num,xyz,mass
      atomlist.append( ATOM(symb,num,xyz,mass,vel) )
    iline+=1
    colcount=0

    for iatom in range(natom):
      if colcount==6:
        iline+=1
        colcount=0
      line=data[iline]
      s=line.split()
      vel=[ float(i)*AMBERVEL_TO_AU    for i in s[colcount:colcount+3] ]
      colcount+=3
      #print iatom,vel
      for i in range(3):
        atomlist[iatom].coord[i]-=dt*vel[i]
      atomlist[iatom].veloc=vel
      atomlist[iatom].EKIN()

    igeom+=1
    if not INFOS['KTR']:
      restore_center_of_mass(atomlist)
      remove_translations(atomlist)
      remove_rotations(atomlist)
    sys.stdout.write('Structure % 5i: %s  ' % (igeom,INFOS['filename_rsts'][igeom-1]))
    if igeom==1:
      sys.stdout.write('(Reference geometry)')
      molecule=INITCOND(atomlist,0.,0.)
    else:
      sys.stdout.write('(Saved for initconds)')
      ic_list.append( INITCOND(atomlist,0.,0.) )
    print ''



  return molecule,ic_list


# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def create_initial_conditions_string(molecule, ic_list, eref=0.0):
  """This function converts an list of initial conditions into a string."""
  ninit=len(ic_list)
  natom=ic_list[0].natom
  representation='None'
  #eref
  eharm=0.
  #for mode in modes:
    #eharm+=mode['freq']*0.5
  string='''SHARC Initial conditions file, version %s
Ninit     %i
Natom     %i
Repr      %s
Eref      %18.10f
Eharm     %18.10f

Equilibrium
''' % (version,ninit,natom,representation,eref,eharm)
  for atom in molecule.atomlist:
    string+=str(atom)+'\n'
  string+='\n\n'

  for i, ic in enumerate(ic_list):
    string += 'Index     %i\n%s' % (i+1, str(ic))
  return string

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
# ======================================================================================================================
# ======================================================================================================================

def main():
  '''Main routine'''

  # command line option setup
  usage='''
amber_to_initconds.py [options] prmtop rst7 [rst7 [ rst7 ...] ]

This script reads a prmtop and any number of restart files.
The data is then transformed and written to initconds format.
'''
  description=''
  parser = OptionParser(usage=usage, description=description)
  parser.add_option('-t', dest='t', type=float, nargs=1, default=None, help="Time step employed by Amber in femtosecond")
  parser.add_option('-o', dest='o', type=str, nargs=1, default='initconds', help="Output filename (string, default=""initconds"")")
  parser.add_option('-x', dest='X', action='store_true',help="Generate a xyz file with the sampled geometries in addition to the initconds file")
  parser.add_option('-m', dest='m', action='store_true',help="Enter non-default atom masses")
  parser.add_option('--keep_trans_rot', dest='KTR', action='store_true',help="Keep translational and rotational components")
  parser.add_option('--use_zero_veloc', dest='UZV', action='store_true',help="For all samples, set velocities to zero")

  # arg processing
  (options, args) = parser.parse_args()
  if len(args)==0:
    print usage
    quit(1)

  # options
  INFOS={}
  INFOS['filename_prmtop']=args[0]
  INFOS['filename_rsts']=args[1:]
  if not options.t:
    print 'ERROR: please specify the length of the time step employed in Amber with the -t option!\nThis is necessary to correctly convert the Amber data (leapfrog-style) to SHARC data (velocity-Verlet style).'
    sys.exit(1)
  else:
    INFOS['timestep']=options.t
  INFOS['outfile']=options.o
  INFOS['masslist']={}
  if options.m:
    INFOS['masslist']=ask_for_masses()
  INFOS['KTR']=options.KTR
  INFOS['UZV']=options.UZV




  print '''Initial condition generation started...
prmtop file                  = "%s"
rst files                    = "%s"
OUTPUT file                  = "%s"
Number of geometries         = %i''' % (INFOS['filename_prmtop'], 
                                        INFOS['filename_rsts'], 
                                        INFOS['outfile'], 
                                        len(INFOS['filename_rsts'])
                                        )


  #print 'Generating %i initial conditions' % amount
  molecule, ic_list = get_coords(INFOS)
  #print 'Writing output to initconds'
  outfile = open(INFOS['outfile'], 'w')
  outstring = create_initial_conditions_string(molecule, ic_list)
  outfile.write(outstring)
  outfile.close()

  if options.X:
    make_dyn_file(ic_list,options.o+'.xyz')

  # save the shell command
  command='python '+' '.join(sys.argv)
  f=open('KEYSTROKES.amber_to_initconds','w')
  f.write(command)
  f.close()

# ======================================================================================================================

if __name__ == '__main__':
    main()
