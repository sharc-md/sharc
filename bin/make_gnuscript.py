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

# ======================================================================= #
# Modules:
# Operating system, isfile and related routines, move files, create directories
import os
# External Calls to MOLPRO
import subprocess as sp
# Command line arguments
import sys
# Regular expressions
import re
# debug print for dicts and arrays
import pprint
# sqrt and other math
import math
# runtime measurement
import datetime
# copy of arrays of arrays
from copy import deepcopy
# colorstuff
import colorsys

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

# ======================================================================= #

version='2.0'
versiondate=datetime.date(2018,2,1)

# hash table for conversion of multiplicity to the keywords used in MOLPRO
IToMult={
         1: 'Singlet', 
         2: 'Doublet', 
         3: 'Triplet', 
         4: 'Quartet', 
         5: 'Quintet', 
         6: 'Sextet ', 
         7: 'Septet ', 
         8: 'Octet  ', 
         'Singlet': 1, 
         'Doublet': 2, 
         'Triplet': 3, 
         'Quartet': 4, 
         'Quintet': 5, 
         'Sextet ': 6, 
         'Septet ': 7, 
         'Octet  ': 8
         }

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
  3 integer: MS value'''

  n=1
  for i in range(len(states)):
    if states[i]<1:
      continue
    for k in range(i+1):
      for j in range(states[i]):
        yield i+1,j+1,k-i/2.,n
        n+=1
  return


# ======================================================================= #
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
    self.excluded=[
      [0.12,0.22]       # exclude yellow hues from the colorwheel
      ]
    # number of non-empty groups
    self.initlist=initlist
    n=0
    for index,el in enumerate(initlist):
      self.initlist[index]=max(0,el)
      if el>0:
        n+=1
    self.n=n
    self.m=len(initlist)
    # available colorspace
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
    triple=[0,0,0]
    for i in range(3):
      triple[i]=max(min(rgb[i],1.0),0.0)
    color=int(255*triple[0])*256**2+int(255*triple[1])*256+int(255*triple[2])
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

# ========================== Main Code =============================== #
def main():

  if len(sys.argv)==1:
    print 'Usage:\n./make_gnuscript.py <S> <D> <T> <Q> <5> <6> <7> <8>\n'
    quit(1)
  if len(sys.argv)>9:
    print 'Only multiplicities up to octets are supported!'
    quit(1)
  states=sys.argv[1:]
  for i in range(len(states)):
    states[i]=int(states[i])
  nstates=0
  nmstates=0
  for mult,i in enumerate(states):
    nstates+=i
    nmstates+=(mult+1)*i
  maxmult=len(states)

  gnustring=''

  # first plot: energies in diag picture, with spin expectation value and dipole expectation value
  # needs to have the expec.out file from data_extractor
  # write header for first plot
  gnustring+='unset key\nunset colorbox\nset title "Energies in diagonal basis"\nset xlabel "Time t in fs"\nset ylabel "Energy in eV"\nset cbrange [0:17]\n'
  gnustring+='set palette defined (0.0 "gray90", 1e-5 "gray60", 1e-4 "gray30", 1e-3 "orange", 1e-2 "red", 1e-1 "magenta", 1e-0 "blue", 10 "blue", 11 "green", 12 "red", 13 "turquoise", 14 "orange", 15 "cyan", 16 "brown", 17 "skyblue")\n\n'
  gnustring+=  'plot "output_data/expec.out" u 1:($%i)               title "Total Energy" lw % 6.2f lc rgbcolor "#000000" w l, \\\n' % (4,0.5)
  for i in range(nmstates):
    gnustring+='""               u 1:($%i):(abs($%i)+10) title "State %i"     lw % 6.2f pal w l, \\\n'                   % (5+i,5+1*nmstates+i,i+1,4.5)
  #for i in range(nmstates):
    gnustring+='""               u 1:($%i):(abs($%i))    title "State %i"     lw % 6.2f pal w l, \\\n'                   % (5+i,5+2*nmstates+i,i+1,3.5)
  gnustring+=  '""               u 1:($%i)               title "Trajectory"   lw % 6.2f lc rgbcolor "#000000" pt 6 w p\n\n' % (3,1.0)
  gnustring+='pause -1\n\n'

  # second plot: coefficients in MCH picture
  #angdiff=1./maxmult/max(states)
  R=rgbcolor(states)
  gnustring+='set key\nset yrange [0:1]\nset ylabel "Wavefunction Amplitude"\nset title "MCH Quantum Amplitudes"\n'
  gnustring+='plot "output_data/coeff_MCH.out" \tu 1:2  \ttitle "Sum of Amplitudes" \tlw % 6.2f \tlc rgbcolor "#000000" \tw l, \\\n' % (3.0)
  for mult,state,ms,n in itnmstates(states):
    gnustring+='"" \t\t\tu 1:($%i**2+$%i**2)  \ttitle "%s %i, %i" \tlw % 6.2f \tlc rgbcolor "%s" \tw l' %    (1+2*n,2+2*n,IToMult[mult][0],state-(mult==1),ms,mult,R.hexcolor(mult,state) )
    if n<nmstates:
      gnustring+=', \\\n'
  gnustring+='\n\npause -1\n\n'

  # third plot: coefficients in diag picture
  #angdiff=math.ceil(nmstates/2+0.5)/nmstates
  R=rgbcolor([nmstates])
  gnustring+='set yrange [0:1]\nset ylabel "Wavefunction Amplitude"\nset title "Diag Quantum Amplitudes"\n'
  gnustring+='plot "output_data/coeff_diag.out" \tu 1:2  \ttitle "Sum of Amplitudes" \tlw % 6.2f \tlc rgbcolor "#000000" \tw l, \\\n' % (3.0)
  i=1
  for state in range(nmstates):
    gnustring+='"" \t\t\tu 1:($%i**2+$%i**2)  \ttitle "State %i" \tlw % 6.2f \tlc rgbcolor "%s" \tw l' %    (3+2*state,4+2*state,state+1,2.0,R.hexcolor(1,state+1))
    if i<nmstates:
      gnustring+=', \\\n'
    i+=1
  gnustring+='\n\npause -1\n\n'

  # fourth plot: hopping probabilities in diag picture
  R=rgbcolor([nmstates])
  gnustring+='set xrange[0:*]\nset ylabel "Hopping Probability"\nset title "Diag Hopping Probablities and Random Number"\nset style fill solid 0.25 border\n'
  for i in range(nmstates,0,-1):
    if i==nmstates:
      gnustring+='plot "output_data/prob.out"\t'
    else:
      gnustring+='""\t\t'
    gnustring+='u 1:%i  \ttitle "State %i" \tlw % 6.2f \tlc rgbcolor "%s" \tw boxes, \\\n' %    (2+i,i,1.0,R.hexcolor(1,i))
  gnustring+='""\t\tu 1:($1>0 ? $2 : 1/0)\ttitle "Random number" \t\t\tlc rgbcolor "black"\tw lp'
  gnustring+='\n\npause -1\n\n'

  print gnustring

# ======================================================================= #

if __name__ == '__main__':
    main()