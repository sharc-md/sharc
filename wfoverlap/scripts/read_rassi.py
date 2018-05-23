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
"""
This program reads RASSI output and converts it to the native ASCII format
used by cioverlap.x
Code adapted from SHARC (written by S. Mai)
"""

import itertools, math
from copy import deepcopy

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
    sys.exit(12)
  return out

# ======================================================================= #
def writefile(filename,content,lvprt=1):
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
    if lvprt>=1:
        print 'File %s written.'%filename
  except IOError:
    print 'Could not write to file %s!' % (filename)
    sys.exit(13)    

# ======================================================================= #

# =============================================================================================== #
# =============================================================================================== #
# ===========================================  Overlaps  ======================================== #
# =============================================================================================== #
# =============================================================================================== #

def decompose_csf(ms2,step):
  # ms2 is M_S value
  # step is step vector for CSF (e.g. 3333012021000)

  def powmin1(x):
    a=[1,-1]
    return a[x%2]

  # calculate key numbers
  nopen=sum( [ i==1 or i==2 for i in step] )
  nalpha=int(nopen/2.+ms2)
  norb=len(step)

  # make reference determinant
  refdet=deepcopy(step)
  for i in range(len(refdet)):
    if refdet[i]==1:
      refdet[i]=2

  # get the b vector and the set of open shell orbitals
  bval=[]
  openorbs=[]
  b=0
  for i in range(norb):
    if step[i]==1:
      b+=1
    elif step[i]==2:
      b-=1
    bval.append(b)
    if refdet[i]==2:
      openorbs.append(i)

  # loop over the possible determinants
  dets={}
  # get all possible combinations of nalpha orbitals from the openorbs set
  for localpha in itertools.combinations(openorbs, nalpha):
    # make determinant string
    det=deepcopy(refdet)
    for i in localpha:
      det[i]=1

    # get coefficient
    coeff=1.
    sign=+1
    m2=0
    for k in range(norb):
      if step[k]==0:
        pass
        #sign*=powmin1(bval[k])
      elif step[k]==1:
        m2+=powmin1(det[k]+1)
        num=bval[k]+powmin1(det[k]+1)*m2
        if num==0.:
          break
        denom=2.*bval[k]
        sign*=powmin1(bval[k])
        coeff*=1.*num/denom
      elif step[k]==2:
        m2+=powmin1(det[k]-1)
        num=bval[k]+2+powmin1(det[k])*m2
        if num==0.:
          break
        denom=2.*(bval[k]+2)
        sign*=powmin1(bval[k]-det[k]+1)
        coeff*=1.*num/denom
      elif step[k]==3:
        #sign*=powmin1(bval[k])
        num=1.

    # add determinant to dict if coefficient non-zero
    if num!=0.:
      dets[tuple(det)]=1.*sign*math.sqrt(coeff)

  #pprint.pprint( dets)
  return dets

# ======================================================================= #
def get_determinants(out,mult):

    # first, find the correct RASSI output section for the given multiplicity
    modulestring='MOLCAS executing module RASSI'
    #modulestring='     &RASSI'
    spinstring='SPIN MULTIPLICITY:'
    stopstring='The following data are common to all the states'
    module=False
    jobiphmult=[]
    for iline, line in enumerate(out):
        if modulestring in line:
            module=True
            jobiphmult=[]
        elif module:
            if spinstring in line:
                jobiphmult.append(int(line.split()[-1]))
            if stopstring in line:
                print 'jobiphmult: ', jobiphmult
                if all(i==mult for i in jobiphmult):
                    break
                else:
                    module=False
    else:
        print 'Determinants not found!', mult
        print 'No RASSI run for multiplicity %i found!' % (mult)
        sys.exit(15)

    # ndocc and nvirt
    ndocc=-1
    nvirt=-1
    while True:
        iline+=1
        line=out[iline]
        if ' INACTIVE ' in line:
            ndocc=int(line.split()[1])
        if ' SECONDARY ' in line:
            nvirt=int(line.split()[1])
        if ndocc!=-1 and nvirt!=-1:
            break

    # Get number of states
    while True:
        iline+=1
        line=out[iline]
        if 'Nr of states:' in line:
            nstates=int(line.split()[-1])
            break

    # Now start searching at iline, collecting all CI vectors
    # the dict has det-string tuples as keys and lists of float as values
    ci_vectors={}
    statesstring='READCI called for state'
    stopstring='****************************'

    for istate in range(nstates):
        while True:
            iline+=1
            line=out[iline]
            if statesstring in line:
                state=int(line.split()[-1])
                if state==istate+1: break
        iline+=10
        while True:
            iline+=1
            line=out[iline]
            if stopstring in line:
                break
            s=line.split()
            if s==[]:
                continue
            try:
              coef=float(s[-2])
            except ValueError:
              print "WARNING cannot read line for state %i:"%(istate+1)
              print line.rstrip()
              coef=0.
            if coef==0.:
                continue
            csf=s[-3]
            step=[]
            for i in csf:
                if i=='2':
                    step.append(3)
                elif i=='d':
                    step.append(2)
                elif i=='u':
                    step.append(1)
                elif i=='0':
                    step.append(0)
            dets=decompose_csf( (mult-1)/2. ,step)
            # add dets to ci_vectors
            for det in dets:
                if det in ci_vectors:
                    d=state-len(ci_vectors[det])
                    if d>0:
                        ci_vectors[det].extend( [0.]*d )
                    ci_vectors[det][state-1]+=coef*dets[det]
                else:
                    ci_vectors[det]=[0.]*state
                    ci_vectors[det][state-1]+=coef*dets[det]
    for det in ci_vectors.keys():
        d=nstates-len(ci_vectors[det])
        if d>0:
            ci_vectors[det].extend( [0.]*d )
        #if all( i==0. for i in ci_vectors[det] ):
            #del ci_vectors[det]
    ci_vectors['ndocc']=ndocc
    ci_vectors['nvirt']=nvirt
    #pprint.pprint( ci_vectors)
    return ci_vectors

# ======================================================================= #
def format_ci_vectors(ci_vectors):
    # get nstates, norb and ndets
    for key in ci_vectors:
        if key!='ndocc' and key!='nvirt':
            nstates=len(ci_vectors[key])
            norb=len(key)
    ndets=len(ci_vectors)-2
    ndocc=ci_vectors['ndocc']
    nvirt=ci_vectors['nvirt']

    # sort determinant strings
    dets=[]
    for key in ci_vectors:
        if key!='ndocc' and key!='nvirt':
            dets.append(key)
    dets.sort(reverse=True)

    string='%i %i %i\n' % (nstates,norb+ndocc+nvirt,ndets)
    for det in dets:
        for i in range(ndocc):
            string+='d'
        for o in det:
            if o==0:
                string+='e'
            elif o==1:
                string+='a'
            elif o==2:
                string+='b'
            elif o==3:
                string+='d'
        for i in range(nvirt):
            string+='e'
        for c in ci_vectors[det]:
            string+=' %16.12f ' % c
        string+='\n'
    return string

if __name__ == '__main__':
    import sys
    
    print "read_rassi.py <molcas_log> <mult> [<det_file>]"
    
    if len(sys.argv) <= 2:
        print "Enter at least 2 argmuments"
        sys.exit()
    mlog = sys.argv[1]
    mult = int(sys.argv[2])    
    try:
        wname = sys.argv[3]
    except IndexError:
        wname = 'dets'
    
    print "Reading %s ..."%mlog
    out = readfile(mlog)
    ci_vectors = get_determinants(out, mult)
    nfrozen = ci_vectors['ndocc']
    string=format_ci_vectors(ci_vectors)
    writefile(wname,string)    