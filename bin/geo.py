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

import os
import sys
import math
import copy
import re
import datetime
from optparse import OptionParser

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

version='2.0'
versiondate=datetime.date(2018,2,1)

allowedreq=['a','d','r','p','q','x','y','z','5','6','c','i','j','k','l']

# This array contains all Boeyens classification symbols for 5-membered rings.
# (phi): symbol
BOEYENS_5 = {
  (  0.): 'E_1',
  ( 36.): '2^E',
  ( 72.): 'E_3',
  (108.): '4^E',
  (144.): 'E_5',
  (180.): '1^E',
  (216.): 'E_2',
  (252.): '3^E',
  (288.): 'E_4',
  (324.): '5^E',
  ( 18.): '^2H_1',
  ( 54.): '^2H_3',
  ( 90.): '^4H_3',
  (126.): '^4H_5',
  (162.): '^1H_5',
  (198.): '^1H_2',
  (234.): '^3H_2',
  (270.): '^3H_4',
  (306.): '^5H_4',
  (342.): '^5H_1'
}

# This array contains all Boeyens classification symbols for 6-membered rings.
# (phi,theta): symbol
BOEYENS_6 = {
  (   0.,   0.): '^1C_4',
  (   0., 180.): '^4C_1',
  (   0., 54.7): '^1E',
  (  60., 54.7): 'E_2',
  ( 120., 54.7): '^3E',
  ( 180., 54.7): 'E_4',
  ( 240., 54.7): '^5E',
  ( 300., 54.7): 'E_6',
  (   0.,125.3): '^4E',
  (  60.,125.3): 'E_5',
  ( 120.,125.3): '^6E',
  ( 180.,125.3): 'E_1',
  ( 240.,125.3): '^2E',
  ( 300.,125.3): 'E_3',
  (  30., 50.8): '^1H_2',
  (  90., 50.8): '^3H_2',
  ( 150., 50.8): '^3H_4',
  ( 210., 50.8): '^5H_4',
  ( 270., 50.8): '^5H_6',
  ( 330., 50.8): '^1H_6',
  (  30.,129.2): '^4H_5',
  (  90.,129.2): '^6H_5',
  ( 150.,129.2): '^6H_1',
  ( 210.,129.2): '^2H_1',
  ( 270.,129.2): '^2H_3',
  ( 330.,129.2): '^4H_3',
  (  30., 67.5): '^1S_2',
  (  90., 67.5): '^3S_2',
  ( 150., 67.5): '^3S_4',
  ( 210., 67.5): '^5S_4',
  ( 270., 67.5): '^5S_6',
  ( 330., 67.5): '^1S_6',
  (  30.,112.5): '^4S_5',
  (  90.,112.5): '^6S_5',
  ( 150.,112.5): '^6S_1',
  ( 210.,112.5): '^2S_1',
  ( 270.,112.5): '^2S_3',
  ( 330.,112.5): '^4S_3',
  (   0.,  90.): '^{1,4}B',
  (  60.,  90.): 'B_{2,5}',
  ( 120.,  90.): '^{3,6}B',
  ( 180.,  90.): 'B_{4,1}',
  ( 240.,  90.): '^{2,5}B',
  ( 300.,  90.): 'B_{3,6}',
  (  30.,  90.): '^4T_2',
  (  90.,  90.): '^6T_2',
  ( 150.,  90.): '^3T_1',
  ( 210.,  90.): '^2T_4',
  ( 270.,  90.): '^2T_6',
  ( 330.,  90.): '^1T_3',
}

# ================== global variables for options ================= #

p=4     # default number of decimals
f=20    # default field width
Bohrs=False
Radians=False
ang2bohr=1.88972612
deg2rad=math.pi/180.0

# =================== general 3D vector operations ================== #

def rscalar3d(a,b):
  '''Scalar product of 3-dimensional vectors a and b.'''
  s=0
  for i in range(3):
    s+=a[i]*b[i]
  return s

def rcross3d(a,b):
  '''Vector product of 3-dimensional vectors a and b.'''
  s=[0,0,0]
  s[0]=a[1]*b[2]-a[2]*b[1]
  s[1]=a[2]*b[0]-a[0]*b[2]
  s[2]=a[0]*b[1]-a[1]*b[0]
  return s

def rnorm3d(a):
  '''Norm of a 3-dimensional vector.'''
  return math.sqrt(rscalar3d(a,a))

def rangle3d(a,b):
  '''Angle between two 3-dimensional vectors.'''
  x = rscalar3d(a,b) / ( rnorm3d(a)*rnorm3d(b) )
  angle=math.acos( max( -1.0, min( 1.0, x) ) )
  if Radians:
    return angle
  else:
    return angle/deg2rad

def rmult3d(scal,vec):
  '''Product of a scalar number with a 3d vector.'''
  for i in range(3):
    vec[i]*=scal
  return vec

# ===================== Cremer Pople parameters ==================== #

def center(atoms):
  '''Calculates the geometric center of N atoms.'''
  N=len(atoms)
  r=[0.,0.,0.]
  for j in range(N):
    for i in range(3):
      r[i]+=atoms[j][i]/N
  return r

def translate(atoms):
  '''Shifts atom coordinates so that their geometric center is the origin of the coordinate system.'''
  r=center(atoms)
  N=len(atoms)
  for j in range(N):
    for i in range(3):
      atoms[j][i]-=r[i]
  return atoms

def normvec(atoms):
  '''Calculates a vector which is orthogonal to the mean plane defined by the atoms.'''
  N=len(atoms)
  r1=[0,0,0]
  for j in range(N):
    factor=math.sin(2.*math.pi*j/N)
    for i in range(3):
      r1[i]+=factor*atoms[j][i]
  r2=[0,0,0]
  for j in range(N):
    factor=math.cos(2.*math.pi*j/N)
    for i in range(3):
      r2[i]+=factor*atoms[j][i]
  n=rcross3d(r1,r2)
  norm=rscalar3d(n,n)
  for i in range(3):
    n[i]/=math.sqrt(norm)
  return n

def project(atoms,n):
  '''Calculates the z coordinates relative to the mean plane defined by the atoms.'''
  N=len(atoms)
  z=[]
  for j in range(N):
    z.append(rscalar3d(atoms[j],n))
  return z

def getz(atoms):
  '''Uses the above routines to calculate the relative z coordinates relative to a 
mean plane of the atoms.'''
  atoms=translate(atoms)
  #N=len(atoms)
  n=normvec(atoms)
  z=project(atoms,n)
  return z

# ================= coordinate calculation routines ===================== #

def dist(a,b):
  '''Bond distance between atoms a and b.'''
  s=[0,0,0]
  for i in range(3):
    s[i]=b[i]-a[i]
  if Bohrs:
    return rnorm3d(s)*ang2bohr
  else:
    return rnorm3d(s)

def ang(a,b,c):
  '''Bond angle of a-b and b-c.'''
  r1=[0,0,0]
  r2=[0,0,0]
  for i in range(3):
    r1[i]=a[i]-b[i]
    r2[i]=c[i]-b[i]
  return rangle3d(r1,r2)

def dih(a,b,c,d):
  '''Dihedral angle between the a-b-c and b-c-d planes.'''
  r1=[0,0,0]
  r2=[0,0,0]
  r3=[0,0,0]
  for i in range(3):
    r1[i]=b[i]-a[i]
    r2[i]=c[i]-b[i]
    r3[i]=d[i]-c[i]
  q1=rcross3d(r1,r2)
  q2=rcross3d(r2,r3)
  if q1==[0.,0.,0.] or q2==[0.,0.,0.]:
    sys.stderr.write('Undefined dihedral angle!')
    return float('NaN')
  Q=rcross3d(q1,q2)
  if Q==[0.,0.,0.]:
    sign=1.
  elif rangle3d(Q,r2) < 90.*[deg2rad,1.][Radians==None]:
    sign=1.
  else:
    sign=-1.
  return rangle3d(q1,q2)*sign

def pyr(a,b,c,d):
  '''Pyramidalization angle between the a-b bond and the b-c-d plane.'''
  r1=[0,0,0]
  r2=[0,0,0]
  r3=[0,0,0]
  for i in range(3):
    r1[i]=a[i]-b[i]
    r2[i]=c[i]-b[i]
    r3[i]=d[i]-b[i]
  q1=rcross3d(r2,r3)
  if q1==[0.,0.,0.]:
    sys.stderr.write('Undefined pyramidalization angle!')
    return float('NaN')
  if Radians:
    return 90.0*deg2rad-rangle3d(q1,r1)
  else:
    return 90.0-rangle3d(q1,r1)

def pyr_new(a,b,c,d):
  '''Angle between the a-b bond and average of the b-c and b-d bonds.'''
  r1=[0,0,0]
  r2=[0,0,0]
  for i in range(3):
    r1[i]=a[i]-b[i]
    r2[i]=(d[i]+c[i])/2.-b[i]
  if Radians:
    phi=180.0*deg2rad-rangle3d(r1,r2)
  else:
    phi=180.0-rangle3d(r1,r2)
  if pyr(a,b,c,d)<0.:
    phi*=-1.
  return phi

def x(a):
  '''X coordinate of atom a.'''
  return a[0]

def y(a):
  '''Y coordinate of atom a.'''
  return a[1]

def z(a):
  '''Z coordinate of atom a.'''
  return a[2]

def CP5(a,b,c,d,e):
  '''Calculates the Cremer-Pople parameters for a 5-membered ring, q and phi.'''
  atoms2=[a,b,c,d,e]
  atoms=copy.deepcopy(atoms2)
  # getz() changes the atom positions, so we have to copy the positions first! 
  z=getz(atoms)
  h1=0.
  h2=0.
  for j in range(5):
    h1+=z[j]*math.cos(4.*math.pi*j/5.)
    h2-=z[j]*math.sin(4.*math.pi*j/5.)
  h1*=math.sqrt(2./5.)
  h2*=math.sqrt(2./5.)
  ph=math.atan2(h2,h1)
  if ph<0:
    ph=2*math.pi+ph
  q=math.sqrt(h1**2+h2**2)
  if Bohrs:
    q*=ang2bohr
  if Radians:
    return q,ph
  else:
    return q,abs(ph)/deg2rad

def CP6(a,b,c,d,e,f):
  '''Calculates the Cremer-Pople parameters for a 6-membered ring, Q, phi and theta.'''
  atoms2=[a,b,c,d,e,f]
  atoms=copy.deepcopy(atoms2)
  # getz() changes the atom positions, so we have to copy the positions first! 
  z=getz(atoms)
  h1=0.
  h2=0.
  for j in range(6):
    h1+=z[j]*math.cos(4.*math.pi*j/6.)
    h2-=z[j]*math.sin(4.*math.pi*j/6.)
  h1*=math.sqrt(2./6.)
  h2*=math.sqrt(2./6.)
  ph=math.atan2(h2,h1)+math.pi
  if ph<0:
    ph=2*math.pi+ph
  q2=math.sqrt(h1**2+h2**2)
  q3=0.
  for j in range(6):
    q3+=(-1)**(j)*z[j]
  q3/=math.sqrt(6.)
  Q=math.sqrt(q2**2+q3**2)
  if Bohrs:
    Q*=ang2bohr
  th=math.atan2(q3,q2)+math.pi/2
  if Radians:
    return Q,ph,th
  else:
    return Q,ph/deg2rad,th/deg2rad

def orthodromic_distance(phi1_degree, theta1_degree, phi2_degree, theta2_degree):
  '''Calculates the orthodromic (great-circle) distance on a sphere between two points (phi1,theta1) and (phi2,theta2)'''
  breite1 = 90. - theta1_degree
  breite2 = 90. - theta2_degree
  if phi1_degree > 180.:
    laenge1 = phi1_degree - 360.
  else:
    laenge1 = phi1_degree
  if phi2_degree > 180.:
    laenge2 = phi2_degree - 360.
  else:
    laenge2 = phi2_degree
  phi1 = math.radians(laenge1)
  theta1 = math.radians(breite1)
  phi2 = math.radians(laenge2)
  theta2 = math.radians(breite2)
  delta_phi = abs(phi2 - phi1)
  delta_theta = abs(theta2 - theta1)
  return math.acos(math.sin(theta1)*math.sin(theta2) + math.cos(theta1)*math.cos(theta2)*math.cos(delta_phi))

def Boeyens6(ph,th):
  '''Determines the Boeyens symbol for a given phi and theta.'''
  dist_list = []
  if Radians:
    ph/=deg2rad
    th/=deg2rad
  for angles, conformation in BOEYENS_6.items():
      dist_list.append([conformation, orthodromic_distance(angles[0], angles[1], ph, th)])
  dist_list.sort(key=lambda d: d[1])
  return dist_list[0][0]

def angular_distance(a,b):
  return min( (a-b)%360., (b-a)%360.)

def Boeyens5(ph):
  '''Determines the Boeyens symbol for a given phi.'''
  dist_list = []
  if Radians:
    ph/=deg2rad
  for angles,conformation in BOEYENS_5.items():
      dist_list.append([conformation, angular_distance(angles, ph)])
  dist_list.sort(key=lambda d: d[1])
  return dist_list[0][0]

def Angle_between_3rings(a1,a2,a3,b1,b2,b3):
  '''Calculates the angle between the normal vectors of two 3-rings (a1,a2,a3) and (b1,b2,b3).'''
  ring1=[a1,a2,a3]
  ring2=[b1,b2,b3]
  n1=normvec(ring1)
  n2=normvec(ring2)
  return rangle3d(n1,n2)

def Angle_between_4rings(a1,a2,a3,a4,b1,b2,b3,b4):
  '''Calculates the angle between the normal vectors of two 4-rings (a1,a2,a3,a4) and (b1,b2,b3,b4).'''
  ring1=[a1,a2,a3,a4]
  ring2=[b1,b2,b3,b4]
  n1=normvec(ring1)
  n2=normvec(ring2)
  return rangle3d(n1,n2)

def Angle_between_5rings(a1,a2,a3,a4,a5,b1,b2,b3,b4,b5):
  '''Calculates the angle between the normal vectors of two 5-rings (a1,a2,a3,a4,a5) and (b1,b2,b3,b4,b5).'''
  ring1=[a1,a2,a3,a4,a5]
  ring2=[b1,b2,b3,b4,b5]
  n1=normvec(ring1)
  n2=normvec(ring2)
  return rangle3d(n1,n2)

def Angle_between_6rings(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6):
  '''Calculates the angle between the normal vectors of two 6-rings (a1,a2,a3,a4,a5,a6) and (b1,b2,b3,b4,b5,b6).'''
  ring1=[a1,a2,a3,a4,a5,a6]
  ring2=[b1,b2,b3,b4,b5,b6]
  n1=normvec(ring1)
  n2=normvec(ring2)
  return rangle3d(n1,n2)

# ================================================================= #

def checkreq(s,natom):
  '''Checks one line of input.
Checks:
- whether line is empty
- whether first string has 1 character
- whether first string is one of the valid requests
- whether the number of atoms complies with the request
- whether all atom numbers are <= natom
- whether any atom number is repeated

Also converts the atom numbers to integer.'''
  if s==[]:
    return False
  if len(s[0])!=1:
    sys.stderr.write('Requests start with one-letter keys!\n')
    return False
  s[0]=s[0].lower()
  if not s[0] in allowedreq:
    sys.stderr.write('''Valid requests are:\n
    r\tbond distance (2 atoms)\n
    a\tbond angle (3 atoms)\n
    d\tdihedral angle (4 atoms)\n
    p\tpyramidalization angle (4 atoms)\n
    q\tpyramidalization angle, alternative definition (4 atoms)\n
    x\tx coordinate (1 atom)\n
    y\ty coordinate (1 atom)\n
    z\tz coordinate (1 atom)\n
    5\tCremer-Pople parameters q2, phi2, and Boeyens-like term (5 atoms)\n
    6\tCremer-Pople parameters Q, phi, theta and Boeyens term (6 atoms)\n
    i\tangle between normal vectors of two 3-rings (6 atoms)\n
    j\tangle between normal vectors of two 4-rings (8 atoms)\n
    k\tangle between normal vectors of two 5-rings (10 atoms)\n
    l\tangle between normal vectors of two 6-rings (12 atoms)\n
    c\tComment line (shortened to the first 20 characters)\n''')
    return False
  elif s[0]=='x' or s[0]=='y' or s[0]=='z':
    if not len(s)==2:
      sys.stderr.write('Coordinate needs one atom as argument!\n')
      return False
  elif s[0]=='r':
    if not len(s)==3:
      sys.stderr.write('Distance needs two atoms as arguments!\n')
      return False
  elif s[0]=='a':
    if not len(s)==4:
      sys.stderr.write('Angle needs three atoms as arguments!\n')
      return False
  elif s[0]=='d' or s[0]=='p' or s[0]=='q':
    if not len(s)==5:
      sys.stderr.write('Dihedral/Pyramidalization angle needs four atoms as arguments!\n')
      return False
  elif s[0]=='5':
    if not len(s)==6:
      sys.stderr.write('Five-Ring Cremer-Pople parameters needs five atoms as arguments!\n')
      return False
  elif s[0]=='6':
    if not len(s)==7:
      sys.stderr.write('Six-Ring Cremer-Pople parameters needs six atoms as arguments!\n')
      return False
  elif s[0]=='i':
    if not len(s)==7:
      sys.stderr.write('Angle between two 3-rings needs six atoms as arguments!\n')
      return False
  elif s[0]=='j':
    if not len(s)==9:
      sys.stderr.write('Angle between two 4-rings needs eight atoms as arguments!\n')
      return False
  elif s[0]=='k':
    if not len(s)==11:
      sys.stderr.write('Angle between two 5-rings needs ten atoms as arguments!\n')
      return False
  elif s[0]=='l':
    if not len(s)==13:
      sys.stderr.write('Angle between two 6-rings needs twelve atoms as arguments!\n')
      return False
  for i in range(len(s)-1):
    s[i+1]=int(s[i+1])
  atoms=True
  for i in range(len(s)-1):
    if s[i+1]>natom or s[i+1]<=0:
      atoms=False
  if not atoms:
    sys.stderr.write('There are only %i atoms in this molecule!\n' % (natom))
    return False
  atoms=True
  for i in range(len(s)-1):
    for j in range(len(s)-i-2):
      if s[i+1]==s[i+j+2]:
        atoms=False
  if not atoms:
    sys.stderr.write('Repeated atom!\n')
    return False
  return True

def tableheader(req):
  '''Creates a two-line string with the column number and a column label.'''
  comment_bonus=50
  s='#'+' '*(f-6) + '% 5i|' % (1)
  i=0
  for r in req:
    i+=1
    if r[0]=='c':
      s+=' '*comment_bonus
    if r[0]=='i':
      s+=' '*(1)
    if r[0]=='j':
      s+=' '*(max(26,f)-20+1)
    if r[0]=='k':
      s+=' '*(max(32,f)-20+1)
    if r[0]=='l':
      s+=' '*(max(38,f)-20+1)
    s+=' '*(f-5) + '% 5i|' % (i+1)
    if r[0]=='5':
      for j in range(2):
        i+=1
        s+=' '*(f-5) + '% 5i|' % (i+1)
    if r[0]=='6':
      for j in range(3):
        i+=1
        s+=' '*(f-5) + '% 5i|' % (i+1)
  s+='\n'
  s+='#'+' '*(f-5)+'time|'
  for r in req:
    if r[0]=='x' or r[0]=='y' or r[0]=='z':
      s+=' '*(f-4)+r[0]+'%3i|' % (r[1])
    elif r[0]=='r':
      s+=' '*(f-7)+r[0]+'%3i%3i|' % (r[1],r[2])
    elif r[0]=='a':
      s+=' '*(f-10)+r[0]+'%3i%3i%3i|' % (r[1],r[2],r[3])
    elif r[0]=='d' or r[0]=='p' or r[0]=='q':
      s+=' '*(f-13)+r[0]+'%3i%3i%3i%3i|' % (r[1],r[2],r[3],r[4])
    elif r[0]=='5':
      s+=' '*(f-16)+'q'+'%3i%3i%3i%3i%3i|'  % (r[1],r[2],r[3],r[4],r[5])
      s+=' '*(f-17)+'ph'+'%3i%3i%3i%3i%3i|' % (r[1],r[2],r[3],r[4],r[5])
      s+=' '*(f-17)+'By'+'%3i%3i%3i%3i%3i|' % (r[1],r[2],r[3],r[4],r[5])
    elif r[0]=='6':
      s+=' '*(f-20)+' Q'+'%3i%3i%3i%3i%3i%3i|'  % (r[1],r[2],r[3],r[4],r[5],r[6])
      s+=' '*(f-20)+'ph'+'%3i%3i%3i%3i%3i%3i|' % (r[1],r[2],r[3],r[4],r[5],r[6])
      s+=' '*(f-20)+'th'+'%3i%3i%3i%3i%3i%3i|' % (r[1],r[2],r[3],r[4],r[5],r[6])
      s+=' '*(f-20)+'By'+'%3i%3i%3i%3i%3i%3i|' % (r[1],r[2],r[3],r[4],r[5],r[6])
    elif r[0]=='i':
      s+=' '*(f-20)+'A3'+'%3i%3i%3i,%3i%3i%3i|'  % (r[1],r[2],r[3],r[4],r[5],r[6])
    elif r[0]=='j':
      s+=' '*(f-20)+'A4'+'%3i%3i%3i%3i,%3i%3i%3i%3i|'  % (r[1],r[2],r[3],r[4],r[5],r[6],r[7],r[8])
    elif r[0]=='k':
      s+=' '*(f-20)+'A5'+'%3i%3i%3i%3i%3i,%3i%3i%3i%3i%3i|'  % (r[1],r[2],r[3],r[4],r[5],r[6],r[7],r[8],r[9],r[10])
    elif r[0]=='l':
      s+=' '*(f-20)+'A6'+'%3i%3i%3i%3i%3i%3i,%3i%3i%3i%3i%3i%3i|'  % (r[1],r[2],r[3],r[4],r[5],r[6],r[7],r[8],r[9],r[10],r[11],r[12])
    elif r[0]=='c':
      s+=' '*(f-7+comment_bonus)+'Comment '
  return s

def ang_or_bohr(a):
  if Bohrs:
    return a*ang2bohr
  else:
    return a

def calculate(g,req,comm):
  '''Creates a one-line string containing all requested internal coordinates for geometry g.'''
  comment_bonus=50
  s=''
  formatstring='%%%i.%if ' % (f,p)
  stringstring='%%%is ' % (f)
  commentstring='%%%is ' % (f+comment_bonus)
  for r in req:
    if r[0]=='x':
      s+=formatstring % (ang_or_bohr(g[r[1]-1][0]))
    elif r[0]=='y':
      s+=formatstring % (ang_or_bohr(g[r[1]-1][1]))
    elif r[0]=='z':
      s+=formatstring % (ang_or_bohr(g[r[1]-1][2]))
    elif r[0]=='r':
      s+=formatstring % (dist(g[r[1]-1],g[r[2]-1]))
    elif r[0]=='a':
      s+=formatstring % (ang(g[r[1]-1],g[r[2]-1],g[r[3]-1]))
    elif r[0]=='p':
      s+=formatstring % (pyr(g[r[1]-1],g[r[2]-1],g[r[3]-1],g[r[4]-1]))
    elif r[0]=='q':
      s+=formatstring % (pyr_new(g[r[1]-1],g[r[2]-1],g[r[3]-1],g[r[4]-1]))
    elif r[0]=='d':
      s+=formatstring % (dih(g[r[1]-1],g[r[2]-1],g[r[3]-1],g[r[4]-1]))
    elif r[0]=='5':
      q,ph=CP5(g[r[1]-1],g[r[2]-1],g[r[3]-1],g[r[4]-1],g[r[5]-1])
      s+=formatstring % (q)
      s+=formatstring % (ph)
      s+=stringstring % ('$'+Boeyens5(ph)+'$')
    elif r[0]=='6':
      Q,ph,th=CP6(g[r[1]-1],g[r[2]-1],g[r[3]-1],g[r[4]-1],g[r[5]-1],g[r[6]-1])
      s+=formatstring % (Q)
      s+=formatstring % (ph)
      s+=formatstring % (th)
      s+=stringstring % ('$'+Boeyens6(ph,th)+'$')
    elif r[0]=='i':
      a=Angle_between_3rings(g[r[1]-1],g[r[2]-1],g[r[3]-1],g[r[4]-1],g[r[5]-1],g[r[6]-1])
      s+=' '*(1)
      s+=formatstring % (a)
    elif r[0]=='j':
      a=Angle_between_4rings(g[r[1]-1],g[r[2]-1],g[r[3]-1],g[r[4]-1],g[r[5]-1],g[r[6]-1],g[r[7]-1],g[r[8]-1])
      s+=' '*(max(26,f)-20+1)
      s+=formatstring % (a)
    elif r[0]=='k':
      a=Angle_between_5rings(g[r[1]-1],g[r[2]-1],g[r[3]-1],g[r[4]-1],g[r[5]-1],g[r[6]-1],g[r[7]-1],g[r[8]-1],g[r[9]-1],g[r[10]-1])
      s+=' '*(max(32,f)-20+1)
      s+=formatstring % (a)
    elif r[0]=='l':
      a=Angle_between_6rings(g[r[1]-1],g[r[2]-1],g[r[3]-1],g[r[4]-1],g[r[5]-1],g[r[6]-1],g[r[7]-1],g[r[8]-1],g[r[9]-1],g[r[10]-1],g[r[11]-1],g[r[12]-1])
      s+=' '*(max(38,f)-20+1)
      s+=formatstring % (a)
    elif r[0]=='c':
      if comm[0:f+comment_bonus].strip()=='':
        comm=' '*(f-14+comment_bonus)+'<EMPTY_STRING>'
      s+=commentstring % (comm[0:f+comment_bonus].strip())
  return s
# ================================================================= #

def main():
  '''Main routine.

Parts:
- command line option parsing
- reading the geometry file (creates array geo)
- reading the requested internal coordinates (creates array req)
- writes the table header
- loop over all geometries, calculation of internal coordinates'''

  usage='''
Geo.py [options] < Geo.inp > Geo.out

Reads requests from STDIN
Writes formatted table to STDOUT

Calculates internal coordinates from xyz files.


Geo.py Version %s Date %s

This command line tool calculates internal coordinates (see below) for
xyz files with several consecutive geometries, like from a molecular 
dynamics simulation.

The program can calculate bond lengths, bond angles, dihedral angles 
and pyramidalization angles. For example, a pyramidalization angle for 
a C-NH2 group would be defined as the angle between the C-N bond and 
the NH2 plane. Inter-ring angles, angles between the mean plane of two 
rings, can be computed for pairs of 3-, 4-, 5-, and 6-membered rings.

Additionally, Cremer-Pople parameters for 5- and 6-membered rings can be 
calculated [1]. For 6-membered rings, also the Boeyens classification 
symbols [2] can be generated. For 5-membered rings, the classification
symbols are inspired by the Boeyens scheme.
The program can also output x, y or z coordinates of single atoms. 

Internal coordinates are specified on STDIN, with one coordinate per line.
Each line consists of a one-letter key followed by a number of atom indices,
e.g.
r 1 2
is the bond length between atoms 1 and 2.
The one-letter keys are:
    r\tbond distance (2 atoms)
    a\tbond angle (3 atoms)
    d\tdihedral angle (4 atoms)
    p\tpyramidalization angle (4 atoms)
    q\tpyramidalization angle, alternative definition (4 atoms)
    i\tangle between normal vectors of two 3-rings (6 atoms)
    j\tangle between normal vectors of two 4-rings (8 atoms)
    k\tangle between normal vectors of two 5-rings (10 atoms)
    l\tangle between normal vectors of two 6-rings (12 atoms)
    x\tx coordinate (1 atom)
    y\ty coordinate (1 atom)
    z\tz coordinate (1 atom)
    5\tCremer-Pople parameters q2, phi2 and conformation terms (5 atoms)
    6\tCremer-Pople parameters Q, phi, theta and Boeyens terms (6 atoms)
    c\tComment line from xyz file (truncated to 20 characters)

[1] D. Cremer and J. A. Pople: "A General Definition of Ring Puckering Coordinates",
J. Am. Chem. Soc., 1975, 97, 1354-1358.

[2] J. C. A. Boeyens: "The Conformation of Six-Membered Rings",
J. Cryst. Mol. Struct., 1977, 8, 317-320.
''' % (version,versiondate)

  description=''

  parser = OptionParser(usage=usage, description=description)
  parser.add_option('-p', dest='p', type=int, nargs=1, default=4,help="number of decimals (default=4)")
  parser.add_option('-w', dest='f', type=int, nargs=1, default=20,help="field width (default=20)")
  parser.add_option('-b', dest='b', action='store_true',help="switch to bohrs (default is angstrom)")
  parser.add_option('-r', dest='r', action='store_true',help="switch to radians (default is degrees)")
  parser.add_option('-g', dest='g', type="string", nargs=1, default="output.xyz",help="geometry file in xyz format (default=output.xyz)")
  parser.add_option('-t', dest='t', type=float, nargs=1, default=1.0,help="timestep between successive geometries is fs (default=1.0 fs)")
  (options, args) = parser.parse_args()
  global p,f,Bohrs,Radians
  if options.f>=20:
    f=options.f
  else:
    sys.stderr.write('INFO: Field length f cannot be smaller than 20!')
  if 0<=options.p<=f-3:
    p=options.p
  else:
    p=f-3
  Bohrs=options.b
  Radians=options.r
  dt=options.t

  geofilename=options.g
  try:
    geofile=open(geofilename,'r')
  except IOError:
    sys.stderr.write('ERROR: Geometry file %s does not exist!\n' % (geofilename))
    sys.exit(1)
  geo=geofile.readlines()
  geofile.close()
  natom=int(geo[0].split()[0])

  sys.stderr.write('Enter the internal coordinate specifications:\n')
  answered=False
  req=[]
  iline=0
  while not answered:
    try:
      s=raw_input()
      iline+=1
      s=re.sub('#.*$','',s)
      if 'end' in s:
        answered=True
        continue
      s=s.split()
      s[0]=s[0][0:1]
      if checkreq(s,natom):
        req.append(s)
      elif not sys.stdin.isatty():
        sys.stderr.write('... error on line %i\n' % (iline))
    except EOFError:
      answered=True
    except IndexError:
      pass

  print tableheader(req)
  sys.stderr.write('Number of internal coordinate requests: % 3i\n' % (len(req)) )

  line=0
  t=0
  while line<len(geo):
    try:
      n=int(geo[line].split()[0])
    except IndexError:
      sys.stderr.write('ERROR: did not find number of atoms! Line= %i, step= %i' % (line,t) )
      sys.exit(1)
    if not n==natom:
      sys.stderr.write('ERROR: Number of atoms inconsistent! Line= %i, step= %i' % (line,t) )
      sys.exit(1)
    line+=1
    comm=geo[line]
    line+=1
    g=[]
    try:
      for i in range(natom):
        geoline=geo[line].split()
        a=[]
        for j in range(3):
          a.append(float(geoline[j+1]))
        g.append(a)
        line+=1
    except ValueError:
      sys.stderr.write('ERROR: Error while reading geometry! Line= %i\n' % (line) )
      sys.exit(1)

    formatstring='%%%i.%if ' % (f,p)
    s=calculate(g,req,comm)
    print formatstring % (t*dt) +s
    t+=1
    sys.stderr.write('\rNumber of geometries: % 6i' % (t))

  sys.stderr.write('\nFINISHED!\n\n')
  sys.exit(0)

if __name__ == '__main__':
  try:
    main()
  except KeyboardInterrupt:
    print '\nExited without writing...\n'
    quit(0)