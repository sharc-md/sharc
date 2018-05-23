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
#from copy import copy
import copy
from optparse import OptionParser
import datetime

path=os.path.expanduser(os.path.expandvars('$ADFHOME/scripting'))
sys.path.append(path)
import kf

starttime=datetime.datetime.now()


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

# ======================================================================= #
def get_MO_from_tape21(filename,nfrozen):

    print '%-40s'%'  importing ...',datetime.datetime.now()-starttime

    # get all info from TAPE21
    import kf
    f = kf.kffile(filename)

    print '%-40s'%'  getting info ...',datetime.datetime.now()-starttime

    nspin=int(f.read('ActiveFrag','nspin'))
    if nspin==1:
      restr=True
    else:
      restr=False
    NAO=int(f.read('Basis','naos'))
    npart=f.read('A','npart').tolist()
    NMO_A=int(f.read('A','nmo_A'))
    mocoef_A=f.read('A','Eigen-Bas_A').tolist()
    if not restr:
        NMO_B=int(f.read('A','nmo_B'))
        mocoef_B=f.read('A','Eigen-Bas_B').tolist()

    # prepare npart
    for i in range(len(npart)):
        npart[i]-=1

    print '%-40s'%'  building matrices ...',datetime.datetime.now()-starttime

    # build MO matrices
    MO_A=[ [ 0. for iao in range(NAO) ] for imo in range(NMO_A) ]
    iao=0
    imo=0
    for i,el in enumerate(mocoef_A):
        iao1=npart[iao]
        MO_A[imo][iao1]=el
        iao+=1
        if iao>=NAO:
            iao=0
            imo+=1

    if not restr:
        MO_B=[ [ 0. for iao in range(NAO) ] for imo in range(NMO_B) ]
        iao=0
        imo=0
        for i,el in enumerate(mocoef_B):
            iao1=npart[iao]
            MO_B[imo][iao1]=el
            iao+=1
            if iao>=NAO:
                iao=0
                imo+=1

    if restr:
        NMO=NMO_A      -  nfrozen
    else:
        NMO=NMO_A+NMO_B-2*nfrozen

    print '%-40s'%'  formatting ...',datetime.datetime.now()-starttime

    # make string
    string='''2mocoef
header
 1
MO-coefficients from ADF
 1
 %i   %i
 a
mocoef
(*)
''' % (NAO,NMO)
    x=0
    for imo,mo in enumerate(MO_A):
        if imo<nfrozen:
            continue
        for c in mo:
            if x>=3:
                string+='\n'
                x=0
            string+='% 6.12e ' % c
            x+=1
        if x>0:
            string+='\n'
            x=0
    if not restr:
        x=0
        for imo,mo in enumerate(MO_B):
            if imo<nfrozen:
                continue
            for c in mo:
                if x>=3:
                    string+='\n'
                    x=0
                string+='% 6.12e ' % c
                x+=1
            if x>0:
                string+='\n'
                x=0
    string+='orbocc\n(*)\n'
    x=0
    for i in range(NMO):
        if x>=3:
            string+='\n'
            x=0
        string+='% 6.12e ' % (0.0)
        x+=1

    return string


# ======================================================================= #
def get_dets_from_tape21(filename,mults,wfthres,nfrozen):


    print '%-40s'%'  importing ...',datetime.datetime.now()-starttime

    # get all info from TAPE21
    import kf
    f = kf.kffile(filename)

    nspin=int(f.read('ActiveFrag','nspin'))
    if nspin==1:
      restr=True
    else:
      restr=False

    print '%-40s'%'  getting general info ...',datetime.datetime.now()-starttime
    # get general infos
    if restr:
        if 1 in mults:
            extrmults=['S','T']
        else:
            extrmults=['T']
        gsmult=1
    else:
        extrmults=['S']
        #get gsmult
        occ_A=f.read('A','froc_A').tolist()
        occ_B=f.read('A','froc_B').tolist()
        nA=sum(occ_A)
        nB=sum(occ_B)
        ndiff=int(nA-nB)
        gsmult=ndiff+1
    nstates_to_extract=[ 0 for i in range(max(mults)) ]
    if not restr or 1 in mults:
        n=int(f.read('Excitations SS A','nr of excenergies'))
        nstates_to_extract[gsmult-1]=n
    if restr and 3 in mults:
      n=int(f.read('Excitations ST A','nr of excenergies'))
      if n:
        nstates_to_extract[2]=n
    inp=f.read('General','Input')
    tda=False
    for line in inp:
      if 'tda' in line.lower():
        tda=True




    lhybrid=f.read('General','lhybrid')[0]
    occ_A=f.read('A','froc_A').tolist()
    occ_A=[ int(i) for i in occ_A ]
    if restr:
        nocc_A=sum(occ_A)/2
        nvir_A=len(occ_A)-nocc_A
    else:
        nocc_A=sum(occ_A)
        nvir_A=len(occ_A)-nocc_A
        NMO_B=int(f.read('A','nmo_B'))
        occ_B=f.read('A','froc_B').tolist()
        occ_B=[ int(i) for i in occ_B ]
        nocc_B=sum(occ_B)
        nvir_B=len(occ_B)-nocc_B

    # make step vectors (0:empty, 1:alpha, 2:beta, 3:docc)
    if restr:
        m={0:0,2:3}
        occ_A=[ m[i] for i in occ_A ]
    else:
        m={0:0,1:2}
        occ_B=[ m[i] for i in occ_B ]

    occ_A=tuple(occ_A)
    if not restr:
        occ_B=tuple(occ_B)

    #print occ_A
    #print nocc_A,nvir_A
    #if not restr:
        #print occ_B
        #print nocc_B,nvir_B

    print '%-40s'%'  processing eigenvectors ...',datetime.datetime.now()-starttime
    # get eigenvectors
    eigenvectors={}
    for imult,mult in enumerate(mults):
        print '%-40s'%('    Multiplicity: %i' % (mult)),datetime.datetime.now()-starttime
        eigenvectors[mult]=[]
        if mult==gsmult:
            # add ground state
            if restr:
                key=tuple(occ_A[nfrozen:])
            else:
                key=tuple(occ_A[nfrozen:]+occ_B[nfrozen:])
            eigenvectors[mult].append( {key:1.0} )
        for istate in range(nstates_to_extract[mult-1]):
            print '%-40s'%('      State: %i' % (istate+1)),datetime.datetime.now()-starttime
            print '%-40s'%'        Reading ...',datetime.datetime.now()-starttime
            section='Excitations S%s A' % extrmults[imult]
            key='eigenvector %i' % (istate+1)
            try:
                eig=f.read(section,key).tolist()
            except AttributeError:
                print 'No eigenvectors found in file %s!' % (filename)
                sys.exit(11)
            if lhybrid and not tda:
                key='left eigenvector %i' % (istate+1)
                eigl=f.read(section,key).tolist()
                for i in range(len(eig)):
                    eig[i]=(eig[i]+eigl[i])/2.
            # make dictionary
            print '%-40s'%'        Converting ...',datetime.datetime.now()-starttime
            dets={}
            if restr:
                for iocc in range(nocc_A):
                    for ivirt in range(nvir_A):
                        index=iocc*nvir_A+ivirt
                        dets[ (iocc,ivirt,1) ]=eig[index]
            else:
                for iocc in range(nocc_A):
                    for ivirt in range(nvir_A):
                        index=iocc*nvir_A+ivirt
                        dets[ (iocc,ivirt,1) ]=eig[index]
                # beta excitation
                for iocc in range(nocc_B):
                    for ivirt in range(nvir_B):
                        index=iocc*nvir_B+ivirt + max(nvir_A,nvir_B)*max(nocc_A,nocc_B)
                        dets[ (iocc,ivirt,2) ]=eig[index]
            # truncate vector
            print '%-40s'%'        Truncating vector ...',datetime.datetime.now()-starttime
            norm=0.
            for k in sorted(dets,key=lambda x: dets[x]**2,reverse=True):
                if norm>wfthres:
                    del dets[k]
                    continue
                norm+=dets[k]**2
            # create strings
            print '%-40s'%'        Making determinants ...',datetime.datetime.now()-starttime
            dets2={}
            if restr:
                for iocc,ivirt,dummy in dets:
                    if mult==1:
                        # alpha excitation
                        key=list(occ_A)
                        key[iocc]=2
                        key[nocc_A+ivirt]=1
                        dets2[tuple(key)]=dets[ (iocc,ivirt,dummy) ]/math.sqrt(2.)
                        # beta excitation
                        key[iocc]=1
                        key[nocc_A+ivirt]=2
                        dets2[tuple(key)]=dets[ (iocc,ivirt,dummy) ]/math.sqrt(2.)
                    elif mult==3:
                        key=list(occ_A)
                        key[iocc]=1
                        key[nocc_A+ivirt]=1
                        dets2[tuple(key)]=dets[ (iocc,ivirt,dummy) ]
            else:
                for iocc,ivirt,dummy in dets:
                    if dummy==1:
                        key=list(occ_A+occ_B)
                        key[iocc]=0
                        key[nocc_A+ivirt]=1
                        dets2[tuple(key)]=dets[ (iocc,ivirt,dummy) ]
                    elif dummy==2:
                        key=list(occ_A+occ_B)
                        key[nocc_A+nvir_A + iocc]=0
                        key[nocc_A+nvir_A + nocc_B+ivirt]=2
                        dets2[tuple(key)]=dets[ (iocc,ivirt,dummy) ]
            # remove frozen core
            print '%-40s'%'        Removing frozen core ...',datetime.datetime.now()-starttime
            dets3={}
            for key in dets2:
                problem=False
                if restr:
                    if any( [key[i]!=3 for i in range(nfrozen) ] ):
                        problem=True
                else:
                    if any( [key[i]!=1 for i in range(nfrozen) ] ):
                        problem=True
                    if any( [key[i]!=2 for i in range(nocc_A+nvir_A, nocc_A+nvir_A + nfrozen) ] ):
                        problem=True
                if problem:
                    print 'WARNING: Non-occupied orbital inside frozen core!'
                    continue
                if restr:
                    key2=key[nfrozen:]
                else:
                    key2=key[nfrozen:nocc_A+nvir_A] + key[nocc_A+nvir_A+nfrozen:]
                dets3[key2]=dets2[key]
            # append
            eigenvectors[mult].append(dets3)

    print '%-40s'%'  formatting eigenvectors ...',datetime.datetime.now()-starttime
    strings={}
    for imult,mult in enumerate(mults):
        filename=os.path.join('dets.%i' % mult)
        strings[filename]=format_ci_vectors(eigenvectors[mult])

    return strings

# ======================================================================= #
def format_ci_vectors(ci_vectors):

    # get nstates, norb and ndets
    alldets=set()
    for dets in ci_vectors:
        for key in dets:
            alldets.add(key)
    ndets=len(alldets)
    nstates=len(ci_vectors)
    norb=len(next(iter(alldets)))

    string='%i %i %i\n' % (nstates,norb,ndets)
    for det in sorted(alldets,reverse=True):
        for o in det:
            if o==0:
                string+='e'
            elif o==1:
                string+='a'
            elif o==2:
                string+='b'
            elif o==3:
                string+='d'
        for istate in range(len(ci_vectors)):
            if det in ci_vectors[istate]:
                string+=' %11.7f ' % ci_vectors[istate][det]
            else:
                string+=' %11.7f ' % 0.
        string+='\n'
    return string



# ======================================================================= #
# ======================================================================= #
# ======================================================================= #

parser = OptionParser(usage='', description='')
parser.add_option('-m', dest='m', type=str,   nargs=1, default='1', help="Multiplicities, e.g., '1 3' ")
parser.add_option('-f', dest='f', type=int,   nargs=1, default='0', help="N frozen core, e.g., '5'")
parser.add_option('-t', dest='t', type=float, nargs=1, default='0.99', help="Truncation threshold, e.g., '0.99'")
(options, args) = parser.parse_args()




filename = args[0]
mults=[ int(i) for i in options.m.split() ]
nfrozen=max(0,int(options.f))
wfthres=max(0.,float(options.t))

print '%-40s'%'Starting MO processing ...',datetime.datetime.now()-starttime
string=get_MO_from_tape21(filename,nfrozen)
writefile('mos',string)
print '%-40s'%'Finished MO processing.',datetime.datetime.now()-starttime

print '%-40s'%'Starting determinant processing ...',datetime.datetime.now()-starttime
strings=get_dets_from_tape21(filename,mults,wfthres,nfrozen)
for i in strings:
  writefile(i,strings[i])
print '%-40s'%'Finished determinant processing.',datetime.datetime.now()-starttime



# kate: indent-width 4
