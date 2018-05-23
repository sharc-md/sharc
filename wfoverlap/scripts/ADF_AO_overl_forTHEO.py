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
# runtime measurement
import datetime



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


starttime=datetime.datetime.now()


filename = 'TAPE15'
file1 = kf.kffile(filename)
NAO = int(file1.read('Basis','naos'))
Smat = file1.read("Matrices","Smat")
npart_a = file1.read("A","npart")
npart = npart_a.tolist()


endtime=datetime.datetime.now()
print 'TAPE15 read! \tTook',endtime-starttime,'seconds'
starttime=datetime.datetime.now()


filename = 'TAPE21'
file1 = kf.kffile(filename)
NSO = int(file1.read('SFOs','number'))
SFOmat = file1.read('A','SFO')

endtime=datetime.datetime.now()
print 'TAPE21 read! \tTook',endtime-starttime,'seconds'
starttime=datetime.datetime.now()

#print NAO,NSO
#print npart


Full_Smat=[ [ 0. for i in range(NAO) ] for j in range(NAO) ]
x=0
y=0
anti_npart={}
for i,el in enumerate(npart):
  anti_npart[el-1]=i

for el in Smat:
  #x1 = npart[x]-1
  #y1 = npart[y]-1
  #x1 = npart.index(x+1)
  #y1 = npart.index(y+1)
  x1=anti_npart[x]
  y1=anti_npart[y]
  #x1=x
  #y1=y
  Full_Smat[x1][y1]=el
  Full_Smat[y1][x1]=el
  x+=1
  if x>y:
    y+=1
    x=0

string=''
x=0
for i in range(NAO):
  for j  in range(NAO):
    string+='%12.9f ' % Full_Smat[i][j]
    x+=1
  string+='\n'
writefile('Smat',string)

endtime=datetime.datetime.now()
print 'Smat rearranged! \tTook',endtime-starttime,'seconds'
starttime=datetime.datetime.now()




Full_SFOmat=[ [ 0. for i in range(NAO) ] for j in range(NSO) ]
x=0
y=0
for el in SFOmat:
  Full_SFOmat[x][y]=el
  x+=1
  if x>=NSO:
    y+=1
    x=0

string=''
x=0
for i in range(NSO):
  for j in range(NAO):
    string+='%12.9f ' % Full_SFOmat[i][j]
    x+=1
  string+='\n'
writefile('SFOmat',string)


endtime=datetime.datetime.now()
print 'SFOmat rearranged! \tTook',endtime-starttime,'seconds'
starttime=datetime.datetime.now()



Full_SFOmat=np.array(Full_SFOmat)
Full_Smat  =np.array(Full_Smat)

endtime=datetime.datetime.now()
print 'Numpy arrays! \tTook',endtime-starttime,'seconds'
starttime=datetime.datetime.now()

prodmat=np.dot( Full_SFOmat, np.dot(Full_Smat,Full_SFOmat.T))

endtime=datetime.datetime.now()
print 'Matmul done! \tTook',endtime-starttime,'seconds'
starttime=datetime.datetime.now()

if len(sys.argv)>1:
  print 'Off-diagonal block!'
  string='%i %i\n' % (NSO/2,NSO/2)
  for i in range(NSO/2,NSO):
    for j in range(0,NSO/2):
      string+='%6.12e ' % prodmat[i][j]
    string+='\n'
else:
  print 'Full block!'
  string='%i %i\n' % (NSO,NSO)
  for i in range(NSO):
    for j in range(NSO):
      string+='%6.12e ' % prodmat[i][j]
    string+='\n'
writefile('SFO_overl',string)
print 'DONE!'



