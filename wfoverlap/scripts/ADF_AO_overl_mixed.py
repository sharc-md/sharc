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
from numpy import *
sys.path.append('/usr/license/adf/adf2016.101/scripting')
from kf import *

filename = sys.argv[1]
file1 = kffile(filename)

NAO = file1.read('Basis','naos')
Smat = file1.read("Matrices","Smat")
Square_SMAT = []
Diag_SMAT=[]
i=0
j=[]
k=0
for a in range(0,int(NAO)):
    i=i+1
    j.append(i)
    Smat_column=[]
    l=math.fsum(j)
    for b in range(int(k),int(l)):
       Smat_column.append(Smat[b])
    k=l
    Diag_SMAT.append(Smat_column)

for a in range(0,int(NAO)):
    SMAT_column = []
    for b in range(0,int(NAO)):
        n=len(Diag_SMAT[a])
        if b == a:
           SMAT_column.append(Diag_SMAT[a][a])
        elif b < n : 
           SMAT_column.append(Diag_SMAT[a][b])
        else:
           SMAT_column.append(Diag_SMAT[b][a])
    Square_SMAT.append(SMAT_column)

NAO=int(NAO)

outfile= open('AO_overl','w')


outfile.write(' %i %i \n' %(NAO/2,NAO/2))
for b in range(NAO/2,NAO):
    for a in range(0,NAO/2):


#outfile.write(' %i %i \n' %(NAO,NAO))
#for b in range(0,NAO):
    #for a in range(0,NAO):



        outfile.write('  %6.12e  ' % (float(Square_SMAT[a][b])))
        if a == int(NAO/2)-1:
           outfile.write(' \n')
