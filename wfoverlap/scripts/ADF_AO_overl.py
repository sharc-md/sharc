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

import sys
import os
import re
import string
import math
import imp
import numpy as np
sys.path.append('/usr/license/adf/adf2016.101/scripting')
import kf



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




filename = 'TAPE15'
filename = sys.argv[1]
file1 = kf.kffile(filename)
NAO = int(file1.read('Basis','naos'))
Smat = file1.read("Matrices","Smat")
npart_a = file1.read("A","npart")
npart = npart_a.tolist()


Full_Smat=[ [ 0. for i in range(NAO) ] for j in range(NAO) ]
x=0
y=0
for el in Smat:
  x1 = npart.index(x+1)
  y1 = npart.index(y+1)
  Full_Smat[x1][y1]=el
  Full_Smat[y1][x1]=el
  x+=1
  if x>y:
    y+=1
    x=0

string='%i %i\n' % (NAO,NAO)
x=0
for i in range(NAO):
  for j  in range(NAO):
    string+='%6.12e ' % Full_Smat[i][j]
    x+=1
  string+='\n'
writefile('AO_overl',string)




