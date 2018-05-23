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
#import re
#import string
import math
#import imp
#from numpy import *
from copy import deepcopy
import datetime
import subprocess as sp

path=os.path.expanduser(os.path.expandvars('$ADFHOME/scripting'))
sys.path.append(path)
import kf

PRINT=True
DEBUG=True

# ======================================================================= #
def shorten_DIR(string):
    maxlen=40
    front=12
    if len(string)>maxlen:
        return string[0:front]+'...'+string[-(maxlen-3-front):]
    else:
        return string+' '*(maxlen-len(string))
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
        f.write(line+'\n')
    elif isinstance(content,str):
      f.write(content)
    else:
      print 'Content %s cannot be written to file!' % (content)
      sys.exit(14)
    f.close()
  except IOError:
    print 'Could not write to file %s!' % (filename)
    sys.exit(15)


# ======================================================================= #         OK
def makermatrix(a,b):
    '''Initialises a real axb matrix.

    Arguments:
    1 integer: first dimension
    2 integer: second dimension

    Returns;
    1 list of list of real'''

    mat=[ [ 0. for i in range(a) ] for j in range(b) ]
    return mat


# ======================================================================= #
def get_smat(filename):

    import kf
    f = kf.kffile(filename)
    NAO = f.read('Basis','naos')
    Smat = f.read('Matrices','Smat')
    f.close()

    # Smat is lower triangular matrix, len is NAO*(NAO+1)/2
    ao_ovl=makermatrix(NAO,NAO)
    x=0
    y=0
    for el in Smat:
        ao_ovl[x][y]=el
        ao_ovl[y][x]=el
        x+=1
        if x>y:
            x=0
            y+=1
    return NAO,ao_ovl

# ======================================================================= #
def runADF(WORKDIR,ADF,ncpu,strip=False):
    prevdir=os.getcwd()
    os.chdir(WORKDIR)
    string=os.path.join(ADF,'bin','adf')+' '
    string+='-n %i < ADF.run' % (ncpu)
    if PRINT or DEBUG:
        starttime=datetime.datetime.now()
        sys.stdout.write('START:\t%s\t%s\t"%s"\n' % (shorten_DIR(WORKDIR),starttime,shorten_DIR(string)))
        sys.stdout.flush()
    stdoutfile=open(os.path.join(WORKDIR,'ADF.out'),'w')
    stderrfile=open(os.path.join(WORKDIR,'ADF.err'),'w')
    try:
        runerror=sp.call(string,shell=True,stdout=stdoutfile,stderr=stderrfile)
    except OSError:
        print 'Call have had some serious problems:',OSError
        sys.exit(66)
    stdoutfile.close()
    stderrfile.close()
    if os.path.isfile(os.path.join(WORKDIR,'TAPE13')):
        runerror=1
    if PRINT or DEBUG:
        endtime=datetime.datetime.now()
        sys.stdout.write('FINISH:\t%s\t%s\tRuntime: %s\tError Code: %i\n' % (shorten_DIR(WORKDIR),endtime,endtime-starttime,runerror))
        sys.stdout.flush()
    os.chdir(prevdir)
    if strip and not DEBUG and runerror==0:
        stripWORKDIR(WORKDIR)
    return runerror

# ======================================================================= #
def stripWORKDIR(WORKDIR):
    ls=os.listdir(WORKDIR)
    keep=['ADF.run$','ADF.err$','ADF.out$','TAPE21','TAPE15']
    for ifile in ls:
        delete=True
        for k in keep:
            if containsstring(k,ifile):
                delete=False
        if delete:
            rmfile=os.path.join(WORKDIR,ifile)
            if not DEBUG:
                os.remove(rmfile)

# ======================================================================= #
def get_Double_AOovl(f1,f2):


    # get old geometry
    oldgeo,old_nao=get_geometry(f1)
    newgeo,new_nao=get_geometry(f2)

    # apply shift
    #shift=1e-5
    #for iatom in range(len(oldgeo)):
        #for ixyz in range(3):
            #oldgeo[iatom][1+ixyz]+=shift

    # get input from t21
    inp=get_input(f1)

    # add atoms, save t15, sharcoverlap
    inp2=[]
    inp2.extend( ['units','length bohr','end'] )
    inp2.extend( ['atoms'] )
    i=1
    for atom in oldgeo:
        line='%i %s %f %f %f f=f1' % (i,atom[0],atom[1],atom[2],atom[3])
        inp2.append(line)
        i+=1
    for atom in newgeo:
        line='%i %s %f %f %f f=f2' % (i,atom[0],atom[1],atom[2],atom[3])
        inp2.append(line)
        i+=1
    inp2.extend( ['end','fragments', 'f1 %s' % (f1), 'f2 %s' % (f2), 'end', 'save TAPE15', 'sharcoverlap', 'ignoreoverlap','CALCOVERLAPONLY'] )
    inp2.extend(inp)

    # write input
    writefile('ADF.run',inp2)

    # run the calculation
    WORKDIR=os.getcwd()
    ADF=os.path.expanduser(os.path.expandvars('$ADFHOME/'))
    err=runADF(WORKDIR,ADF,1,strip=True)

    # get output
    filename=os.path.join(WORKDIR,'TAPE15')
    NAO,Smat=get_smat(filename)

    ## Smat is now full matrix NAO*NAO
    ## we want the lower left quarter, but transposed
    string='%i %i\n' % (old_nao,new_nao)
    for icol in range(0,old_nao):
        for irow in range(old_nao,NAO):
            string+='% .15e ' % (Smat[icol][irow])          # note the exchanged indices => transposition
        string+='\n'

    return string

# ======================================================================= #

def get_input(t21file):
    import kf
    f=kf.kffile(t21file)
    inp=f.read('General','Input').tolist()
    #print inp

    d_sections=['units','atoms','fragments','qmmm']
    d_keys=[]

    do=True
    newinp=[]
    for line in inp:
        if not line.split():
            continue
        if not do:
            if 'end'==line.lower().split()[0]:
                do=True
            continue
        if line.lower().split()[0] in d_sections:
            do=False
            continue
        if line.lower().split()[0] in d_keys:
            continue
        newinp.append(line)
    #print newinp

    return newinp


# ======================================================================= #
def get_geometry(t21file):

    import kf
    f=kf.kffile(t21file)
    geom=f.read('Geometry','xyz InputOrder').tolist()
    atomtype=f.read('Geometry','atomtype').tolist()
    atomindex=f.read('Geometry','fragment and atomtype index').tolist()
    atomorder=f.read('Geometry','atom order index').tolist()
    natom=int(f.read('Geometry','nr of atoms'))
    NAO = int(f.read('Basis','naos'))

    geometry=[]
    for iatom in range(natom):
        index=atomorder[iatom]-1
        atype=atomindex[natom+index]-1
        el=atomtype[atype].strip()
        geometry.append( [el, 
                          geom[3*iatom+0],
                          geom[3*iatom+1],
                          geom[3*iatom+2] ] )
    return geometry,NAO


# ======================================================================= #
# ======================================================================= #
# ======================================================================= #

filename1 = sys.argv[1]
filename2 = sys.argv[2]

string=get_Double_AOovl(filename1,filename2)
writefile('AOovl',string)




















# kate: indent-width 4