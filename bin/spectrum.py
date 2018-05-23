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
import colorsys
import re
import pprint


# =========================================================0
# compatibility stuff

if sys.version_info[0]!=2:
  print 'This is a script for Python 2!'
  quit(0)

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

# some constants
DEBUG = False
CM_TO_HARTREE = 1./219474.6     #4.556335252e-6 # conversion factor from cm-1 to Hartree
HARTREE_TO_EV = 27.211396132    # conversion factor from Hartree to eV
U_TO_AMU = 1./5.4857990943e-4            # conversion from g/mol to amu
ANG_TO_BOHR = 1./0.529177211    #1.889725989      # conversion from Angstrom to bohr
PI = math.pi

version='2.0'
versionneeded=[0.2, 1.0, 2.0, float(version)]
versiondate=datetime.date(2018,2,1)


IToMult={
         1: 'Singlet', 
         2: 'Doublet', 
         3: 'Triplet', 
         4: 'Quartet', 
         5: 'Quintet', 
         6: 'Sextet', 
         7: 'Septet', 
         8: 'Octet', 
         'Singlet': 1, 
         'Doublet': 2, 
         'Triplet': 3, 
         'Quartet': 4, 
         'Quintet': 5, 
         'Sextet': 6, 
         'Septet': 7, 
         'Octet': 8
         }

# =========================================================0

class gauss:
  def __init__(self,fwhm):
    self.c=-4.*math.log(2.)/fwhm**2             # this factor needs to be evaluated only once
    self.f=fwhm
    self.norm=self.f/2.*math.sqrt(math.pi/math.log(2.))
  def ev(self,A,x0,x):
    return A*math.exp( self.c*(x-x0)**2)        # this routine does only the necessary calculations

class lorentz:
  def __init__(self,fwhm):
    self.f=fwhm
    self.c=0.25*fwhm**2
    self.norm=math.pi*self.f/2.
  def ev(self,A,x0,x):
    return A/( (x-x0)**2/self.c+1)

class lognormal:
  def __init__(self,fwhm):
    self.f=fwhm
    self.norm=1.       # TODO: currently not implemented 
  def ev(self,A,x0,x):
    if x<=0 or x0<=0:
      return 0.
    # for lognormal distribution, the factor for the exponent depends on x0
    c=(math.log( (self.f+math.sqrt(self.f**2+4.*x0**2))/(2.*x0)))**2
    # note that the function does not take a value of A at x0
    # instead, the function is normalized such that its maximum will have a value of A (at x<=x0)
    return A*x0/x*math.exp( -c/(4.*math.log(2.)) -math.log(2.)*(math.log(x)-math.log(x0))**2/c)

class spectrum:
  def __init__(self,npts,emin,emax,fwhm,lineshape):
    self.npts=npts
    if lineshape==1:
      self.f=gauss(fwhm)
    elif lineshape==2:
      self.f=lorentz(fwhm)
    elif lineshape==4:
      self.f=lognormal(fwhm)
    self.en=[ emin + float(i)/self.npts*(emax-emin) for i in range(self.npts+1) ]       # the energy grid needs to be calculated only once
    self.spec=[ 0. for i in range(self.npts+1) ]
    self.fnorm=self.f.norm
  def add(self,A,x0):
    if A==0.:
      return
    for i in range(self.npts+1):
      self.spec[i]+=self.f.ev(A,x0,self.en[i])

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
      if line=='\n':
        continue
      if line=='':
        print 'Initial condition %i not found in file %s' % (index,f.name)
        quit(1)
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
    epot_harm=0.
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

def check_initcond_version(string,must_be_excited=False):
  if not 'sharc initial conditions file' in string.lower():
    return False
  f=string.split()
  for i,field in enumerate(f):
    if 'version' in field.lower():
      try:
        v=float(f[i+1])
        if not v in versionneeded:
          return False
      except IndexError:
        return False
  if must_be_excited:
    if not 'excited' in string.lower():
      return False
  return True



# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def get_initconds(INFOS):
  ''''''

  try:
    initf=open(INFOS['filename'],'r')
  except IOError:
    print 'Could not open file %s!' % (INFOS['filename'])
    quit(1)

  line=initf.readline()
  if not check_initcond_version(line):
    print 'File malformatted!'
    quit(1)
  if not 'Excited' in line:
    print 'File is no output of excite.py!'
    quit(1)
  try:
    INFOS['ninit']=int(initf.readline().split()[1])
  except ValueError:
    print 'Could not read number of initial conditions!'
    quit(1)
  initf.readline()    # skip natom
  INFOS['repr']=initf.readline().split()[1]
  if INFOS['repr']=='MCH':
    INFOS['diag']=False
  else:
    INFOS['diag']=True
  try:
    INFOS['eref']=float(initf.readline().split()[1])
  except ValueError:
    print 'Could not read reference energy!'
    quit(1)
  if INFOS['eref']==0.:
    print 'WARNING: Reference energy is zero.'

  sys.stdout.write( 'Number of initial conditions: %i\n' % (INFOS['ninit']))
  sys.stdout.write( 'Reference energy %16.10f\n' % (INFOS['eref']))
  sys.stdout.write( 'Representation: %s\n' % (['MCH','diag'][INFOS['diag']]))

  initf.readline()    # skip eharm
  line=initf.readline()
  if 'states' in line.lower():
    f=line.split()
    try:
      INFOS['states']=[ int(i) for i in f[1:] ]
      INFOS['nstate']=0
      for mult,n in enumerate(INFOS['states']):
        INFOS['nstate']+=(mult+1)*n
    except:
      pass


  if INFOS['irange'][0]<=0:
    INFOS['irange'][0]=1
  if INFOS['irange'][1]>INFOS['ninit']:
    INFOS['irange'][1]=INFOS['ninit']
  sys.stdout.write( 'Reading initial conditions %i to %i\n\n' % (INFOS['irange'][0],INFOS['irange'][1]))


  statelist=[]
  for icond in range(INFOS['irange'][0],INFOS['irange'][1]+1):
    # get list of excited states
    initcond=INITCOND()
    initcond.init_from_file(initf,INFOS['eref'],icond)
    if len(initcond.statelist)==0:
      continue
    for i,state in enumerate(initcond.statelist):
      if i+1>len(statelist):
        statelist.append([])
      statelist[i].append(state)

  if len(statelist)!=INFOS['nstate']:
    del INFOS['states']
  INFOS['nstate']=len(statelist)
  sys.stdout.write( 'Number of states: %i\n' % (len(statelist)))
  if len(statelist)==0:
    sys.stdout.write('No excited-state information found! \nPerhaps this is an output file of wigner.py. In this case, \nplease first perform the excited-state calculations using setup_init.py and excite.py!\n\n')
    quit(0)

  sys.stdout.write('Number of initial conditions with excited-state information (per state):\n')
  for i in statelist:
    sys.stdout.write('%i ' % (len(i)))
  sys.stdout.write('\n\n')

  return statelist,INFOS

# ======================================================================================================================

def make_spectra(statelist,INFOS):
  speclist=[ spectrum(INFOS['npts'],INFOS['erange'][0],INFOS['erange'][1],INFOS['fwhm'],INFOS['lineshape']) for i in range(INFOS['nstate']) ]

  width=50
  idone=0
  imax=sum( len(states) for states in statelist)
  done=0

  for istate,states in enumerate(statelist):
    for icond,cond in enumerate(states):
      idone+=1
      if done<idone*width/imax:
        done=idone*width/imax
        sys.stdout.write('\rProgress: ['+'='*done+' '*(width-done)+'] %3i%%' % (done*100/width))
        sys.stdout.flush()

      if not INFOS['selected'] or cond.Excited:
        if INFOS['dos_switch']:
          speclist[istate].add(1.,cond.Eexc)
        else:
          speclist[istate].add(cond.Fosc,cond.Eexc)
  sys.stdout.write('\n')

  return speclist

# ======================================================================================================================

def make_spectra_bootstrap(statelist,INFOS):

  # check for rectangular statelist
  ncond=len(statelist[0])
  for states in statelist:
    if not ncond==len(states):
      print 'Error: Bootstrapping not possible for non-rectangular initconds file!'
      return

  # make list of admissible initial conditions
  admiss=[ i for i in range(ncond) ]

  width=50
  idone=0
  imax=INFOS['bootstraps']
  done=0

  #pprint.pprint(statelist)

  allspec=[]
  for ibootstrap in range(INFOS['bootstraps']):
    # generate list
    use=[]
    for icond in range(len(states)):
      r=random.randint(0,ncond-1)
      use.append(admiss[r])
    #print use
    spec=spectrum(INFOS['npts'],INFOS['erange'][0],INFOS['erange'][1],INFOS['fwhm'],INFOS['lineshape'])
    for istate,states in enumerate(statelist):
      for icond in use:
        cond=states[icond]
        if not INFOS['selected'] or cond.Excited:
          if INFOS['dos_switch']:
            spec.add(1.,cond.Eexc)
          else:
            spec.add(cond.Fosc,cond.Eexc)
    allspec.append(spec)
    idone+=1
    if done<idone*width/imax:
      done=idone*width/imax
      sys.stdout.write('\rProgress: ['+'='*done+' '*(width-done)+'] %3i%%' % (done*100/width))
      sys.stdout.flush()
  sys.stdout.write('\n')

  mean_spec =spectrum(INFOS['npts'],INFOS['erange'][0],INFOS['erange'][1],INFOS['fwhm'],INFOS['lineshape'])
  stdev_specp=spectrum(INFOS['npts'],INFOS['erange'][0],INFOS['erange'][1],INFOS['fwhm'],INFOS['lineshape'])
  stdev_specm=spectrum(INFOS['npts'],INFOS['erange'][0],INFOS['erange'][1],INFOS['fwhm'],INFOS['lineshape'])

  for ipt in range(INFOS['npts']):
    data=[]
    for j in allspec:
      data.append(j.spec[ipt])
    mean_spec.spec[ipt]=mean_geom(data)
    stdev=stdev_geom(data,mean_spec.spec[ipt])
    stdev_specp.spec[ipt]=mean_spec.spec[ipt]*(stdev**3-1.)
    stdev_specm.spec[ipt]=mean_spec.spec[ipt]*(1./stdev**3-1.)

  allspec=[mean_spec,stdev_specp,stdev_specm]+allspec


  print_spectra(allspec,INFOS['bootstrapfile'])



  #return speclist

# ======================================== #
def mean_geom(data):
  s=0.
  try:
    for i in data:
      s+=math.log(i)
  except ValueError:
    return float('nan')
  s=s/len(data)
  return math.exp(s)

# ======================================== #
def stdev_geom(data,mean=-9999):
  if not mean==mean:
    return float('nan')
  if mean==-9999:
    m=math.log(mean_geom(data))
  else:
    m=math.log(mean)
  s=0.
  for i in data:
    s+=(math.log(i)-m)**2
  s=s/(len(data)-1)
  return math.exp(math.sqrt(s))


# ======================================================================================================================

def print_spectra(speclist,outputfile):
  s='#%15i ' % (1)
  for i in range(len(speclist)+1):
    s+='%16i ' % (i+2)
  s+='\n'
  s+='#%15s ' % ('Energy (eV)')
  for i in range(len(speclist)):
    s+='%16i ' % (i+1)
  s+='%16s\n' % ('Total')
  npts=speclist[0].npts
  maxsum=0.
  for i in range(npts+1):
    s+='%16.12f ' % (HARTREE_TO_EV*speclist[0].en[i])
    summe=0.
    for j in range(len(speclist)):
      s+='%16.12f ' % (speclist[j].spec[i])
      summe+=speclist[j].spec[i]
    s+='%16.12f\n' % (summe)
    if summe>maxsum:
      maxsum=summe

  out=open(outputfile,'w')
  out.write(s)
  out.close()
  return maxsum

# ======================================================================================================================

def print_line_spectra(statelist,outputfile,INFOS):
  s='#%15i ' % (1)
  for i in range(len(statelist)+1):
    s+='%16i ' % (i+2)
  s+='\n'
  s+='#%15s ' % ('Energy (eV)')
  for i in range(len(statelist)):
    s+='%16i ' % (i+1)
  s+='%16s\n' % ('Total')
  nstate=len(statelist)
  for i,state in enumerate(statelist):
    for ex in state:
      s+='%16.12f ' % (ex.Eexc*HARTREE_TO_EV)
      for j in range(i):
        s+='%16.12f ' % (0.)

      if ex.Excited or not INFOS['selected']:
        s+='%16.12f ' % (ex.Fosc)
      else:
        s+='%16.12f ' % (0.)

      for j in range(nstate-i-1):
        s+='%16.12f ' % (0.)

      if ex.Excited or not INFOS['selected']:
        s+='%16.12f ' % (ex.Fosc)
      else:
        s+='%16.12f ' % (0.)
      s+='\n'
  out=open(outputfile,'w')
  out.write(s)
  out.close()

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def itnmstates(states):

  x=0
  for i in range(len(states)):
    if states[i]<1:
      continue
    for k in range(i+1):
      for j in range(states[i]):
        x+=1
        yield i+1,j+1,k-i/2.,x
      x-=states[i]
    x+=states[i]
  return

# =============================================

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
    # Hues:
    # 0.0     0.15       ...
    # Red     Yellow     ...
    excluded=[
      [0.12,0.22]       # exclude yellow hues from the colorwheel
      ]
    # excluded-list must be sorted, for each element x[0]<=x[1]
    # each excluded[i][0]<=excluded[i+1][0]
    # everything between 0 and 1
    # sort the pairs
    for i,el in enumerate(excluded):
      excluded[i]=[min(el),max(el)]
    # sort the list
    excluded.sort(key=lambda x:x[0])
    #make all elements between 0 and 1
    temp1=[]
    for i,el in enumerate(excluded):
      temp1.append([ min(1.,max(0.,el[0])), max(0.,min(1.,el[1])) ])
    # resolve overlapping ranges
    temp2=[[0.,0.]]
    for i,el in enumerate(temp1):
      if el[0]>=temp2[-1][1]:
        temp2.append(el)
      else:
        temp2[-1][1]=el[1]
    self.excluded=temp2
    # number of non-empty groups
    self.initlist=initlist
    n=0
    for index,el in enumerate(initlist):
      # negative numbers in initlist are not allowed, make them to zero
      self.initlist[index]=max(0,el)
      if el>0:
        n+=1
    self.n=n                    # number of non-empty groups
    self.m=len(initlist)        # number of groups
    # available colorspace self.a and shifts to skip excluded regions self.ex
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
    #rgb is a list with three elements, which are floats 0<=x<=1
    color=0
    for i in range(3):
      a=max(min(rgb[i],1.0),0.0)
      color+=int(255*a)*256**(2-i)
    #color=int(255*triple[0])*256**2+int(255*triple[1])*256+int(255*triple[2])
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

# =============================================

def make_gnuplot(outputfile,INFOS):
  if not INFOS['diag'] and 'states' in INFOS:
    sortedMCH=True
    statemap={}
    i=1
    for imult,istate,ims,instate in itnmstates(INFOS['states']):
      statemap[i]=[imult,istate,ims,instate]
      i+=1
    nstates=sum(INFOS['states'])
    nmstates=INFOS['nstate']
    init=INFOS['states']
  else:
    sortedMCH=False
    nstates=INFOS['nstate']
    nmstates=INFOS['nstate']
    #angdiff=math.ceil(nstates/2+0.5)/nstates
    init=[INFOS['nstate']]

  R=rgbcolor(init)

  if INFOS['dos_switch']:
    title='Density-of-states spectrum'
  else:
    title='Absorption spectrum'

  gnustring='''set title "%s (%s%s)\\n%i Initial conditions, %s representation"

set xrange [%f:%f]
set yrange [%f:%f]
set xlabel 'Energy (eV)'
set ylabel 'Absorption spectrum (normalized)'

set style fill transparent solid 0.25 border
set term pngcairo size 640,480
set out '%s.png'

''' % (title,
       ['Gaussian','Lorentzian','Lines','Log-normal'][INFOS['lineshape']-1],
       [', FWHM=%f eV' % (INFOS['fwhm']*HARTREE_TO_EV),''][INFOS['lineshape']==3],
       INFOS['ninit'],
       ['MCH','diagonal'][INFOS['diag']],
       INFOS['erange'][0]*HARTREE_TO_EV,
       INFOS['erange'][1]*HARTREE_TO_EV,
       0.,
       1.,
       INFOS['outputfile']
       )

  if 'maxsum' in INFOS:
    gnustring+='N=%f\n\n' % (INFOS['maxsum'])
  else:
    gnustring+='N=1.0\n\n'

  for istate in range(1,1+nstates):
    if istate==1:
      gnustring+='plot "%s" ' % (INFOS['outputfile'])
    else:
      gnustring+='     ""             '
    if sortedMCH:
      gnustring+='u 1:((0'
      for i in statemap:
        if istate==statemap[i][3]:
          gnustring+='+$%i' % (i+1)
          mult=statemap[i][0]
          state=statemap[i][1]
      gnustring+=')/N) w %s tit "%s" lw 0.5 lc rgbcolor "%s"' % (
      ['filledcu','i'][INFOS['lineshape']==3],
      IToMult[mult]+'_%i' % (state-(mult==1 or mult==2)),
      R.hexcolor(mult,state)
      )
    else:
      gnustring+='u 1:(($%i)/N) w %s tit "%s" lw 0.5 lc rgbcolor "%s"''' % (
      istate,
      ['filledcu','i'][INFOS['lineshape']==3],
      'State %i' % (istate),
      R.hexcolor(1,istate)
      )
    if not INFOS['lineshape']==3 or istate!=nstates:
      gnustring+=', \\\n'

  if not INFOS['lineshape']==3:
    gnustring+='     ""             u 1:(($%i)/N)   w l        tit "Sum"       lw 3.0 lc rgbcolor "black"\n\n' % (nmstates+2)

  out=open(outputfile,'w')
  out.write(gnustring)
  out.close()
  print 'Gnuplot script written to "%s"' % (outputfile)

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def main():
  '''Main routine'''

  usage='''
spectrum.py [options] initconds.excited

This script reads a SHARC initconds file containing excited-state information of the
initial conditions and generates a spectrum. The spectrum is written to file.

version %s
date %s
''' % (version,versiondate)

  description=''

  parser = OptionParser(usage=usage, description=description)
  parser.add_option('-n', dest='n', type=int, nargs=1, default=500, help="Grid points for the spectrum (integer, default=500)")
  parser.add_option('-e', dest='e', type=float, nargs=2, default=(0.,10.), help="Energy range (eV) for the spectrum grid (float, default=0 to 10 eV)")
  parser.add_option('-i', dest='i', type=int, nargs=2, default=(1,1000), help="Index range for the evaluation of the spectrum (default=1 to 1000)")
  parser.add_option('-f', dest='f', type=float, nargs=1, default=0.1, help="FWHM (eV) for the spectrum (float, default=0.1 eV)")
  parser.add_option('-o', dest='o', type=str, nargs=1, default='', help="Output filename")
  parser.add_option('-G', dest='G', action='store_true',default=False,help="Use Gaussian convolution (default)")
  parser.add_option('-L', dest='L', action='store_true',default=False,help="Use Lorentzian convolution")
  parser.add_option('-N', dest='N', action='store_true',default=False,help="Use Log-normal convolution")
  parser.add_option('-s', dest='s', action='store_true',default=False,help="Use only selected initial conditions")
  parser.add_option('-l', dest='l', action='store_true',default=False,help="Make a line spectrum instead of a convolution")
  parser.add_option('-D', dest='D', action='store_true',default=False,help="Calculate density of states instead of absorption spectrum")
  parser.add_option('--gnuplot', dest='gp', type=str, nargs=1, default='', help="Write a gnuplot script to this file")

  parser.add_option('-b', dest='b', type=str, nargs=1, default='spectrum_bootstrap.out', help="Output file for bootstrap analysis of total spectrum")
  parser.add_option('-B', dest='B', type=int, nargs=1, default=0, help="Number of bootstrap cycles")
  parser.add_option('-r', dest='r', type=int, nargs=1, default=16661, help="Seed for the random number generator (integer, default=16661)")

  (options, args) = parser.parse_args()


  if len(args)<=0:
    print 'Please give the filename of the initconds.excited file!\n'+usage
    quit(1)
  filename=args[0]
  INFOS={}
  INFOS['npts']=options.n
  sys.stdout.write('Number of grid points: %i\n' % (INFOS['npts']))
  INFOS['erange']=[options.e[0]/HARTREE_TO_EV,options.e[1]/HARTREE_TO_EV]
  sys.stdout.write('Energy range: %.3f to %.3f eV\n' % (INFOS['erange'][0]*HARTREE_TO_EV,INFOS['erange'][1]*HARTREE_TO_EV))
  INFOS['irange']=list(options.i)
  INFOS['fwhm']=options.f/HARTREE_TO_EV
  INFOS['lineshape']=1
  if options.L:
    INFOS['lineshape']=2
  if options.N:
    INFOS['lineshape']=4
  if options.L and options.G:
    sys.stdout.write('WARNING: Both Gaussian and Lorentzian convolution specified, will use Lorentzian!')
  if options.N and options.G:
    sys.stdout.write('WARNING: Both Gaussian and Log-normal convolution specified, will use Log-normal!')
  if options.l:
    INFOS['lineshape']=3
  if INFOS['lineshape']==1:
    sys.stdout.write('Lineshape: Gaussian (FWHM=%.3f eV)\n' % (options.f))
  elif INFOS['lineshape']==2:
    sys.stdout.write('Lineshape: Lorentzian (FWHM=%.3f eV)\n' % (options.f))
  elif INFOS['lineshape']==4:
    sys.stdout.write('Lineshape: Log-normal (FWHM=%.3f eV)\n' % (options.f))
  elif INFOS['lineshape']==3:
    sys.stdout.write('Lineshape: Line Spectrum\n')
  INFOS['selected']=options.s
  INFOS['filename']=filename
  INFOS['dos_switch']=options.D
  random.seed(options.r)
  INFOS['bootstrapfile']=options.b
  INFOS['bootstraps']=options.B

  if options.o=='':
    if INFOS['dos_switch']:
      outputfile='density_of_states.out'
    else:
      outputfile=['spectrum.out','spectrum_line.out'][options.l]
  else:
    outputfile=options.o

  statelist,INFOS=get_initconds(INFOS)
  if options.l:
    print_line_spectra(statelist,outputfile,INFOS)
  else:
    speclist=make_spectra(statelist,INFOS)
    maxsum=print_spectra(speclist,outputfile)
    INFOS['maxsum']=maxsum
    sys.stdout.write('\nMaximum of the absorption spectrum: %f\n' % (maxsum))

    if INFOS['bootstraps']>0:
      make_spectra_bootstrap(statelist,INFOS)
  print '\nOutput spectrum written to "%s".' % (outputfile)

  INFOS['outputfile']=outputfile
  if options.gp!='':
    make_gnuplot(options.gp,INFOS)

  print ''

  # save the shell command
  command='python '+' '.join(sys.argv)
  f=open('KEYSTROKES.spectrum','w')
  f.write(command)
  f.close()

# ======================================================================================================================

if __name__ == '__main__':
    main()
