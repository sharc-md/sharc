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
import re

numargerror="""Usage:
python seward_double-mol.py <dir1> <dir2> [run]
"""

atomlabel=re.compile('[a-zA-Z][a-zA-Z]?[1-9][0-9]*')
numberlabel=re.compile('[-]?[0-9]+[.][0-9]*')
basislabel=re.compile('[bB]asis')
basisinfo=re.compile('[a-zA-Z][a-zA-Z]?[.]')

def isatom(line):
  parts=line.split()
  if len(parts)!=4:
    return False
  match=atomlabel.match(parts[0])
  for i in range(1,3):
    match=match and numberlabel.match(parts[i])
  if match:
    return True
  else:
    return False

def isinfo(line):
  match=basislabel.search(line)
  if match:
    return True
  match=basisinfo.search(line)
  if match:
    return True
  else:
    return False

def isdelete(line):
  if isatom(line) or isinfo(line):
    return False
  else:
    return True

def write_molcasin(olddir, newdir, molcasin_name):
  geomold=open('%s/WORK/molcas.input.seward' % olddir,'r')
  geomnew=open('%s/WORK/molcas.input.seward' % newdir,'r')
  molcasin=open(molcasin_name,'w')
  
  # Write molcas input header
  
  molcastext="""&SEWARD &END
ONEOnly
EXPErt
NOGUessorb
NODElete
NODKroll
NOAMfi
MULTipoles
0

"""
  molcasin.write(molcastext)
  
  # Write old geometry including basis set infos and with increasing labels
  
  data=geomold.readlines()
  geomold.close()
  atom=1
  for ln,line in enumerate(data):
    if isdelete(line):
      continue
    if isinfo(line):
      molcasin.write(line)
    if isatom(line):
      parts=line.split()
      parts[0] = re.sub("\d+", "", parts[0])+str(atom)
      atom+=1
      line=' '.join(parts)
      molcasin.write(line+'\n')
  
  # Write new geometry including basis set infos and with further increasing labels
  
  datanew=geomnew.readlines()
  geomnew.close()
  for ln,line in enumerate(datanew):
    if isdelete(line):
      continue
    if isinfo(line):
      if line!=data[ln]:
        print ' WARNING: Inconsistent basis set information in line %i!' % (ln+1)
        #sys.exit(1)
      molcasin.write(line)
    if isatom(line):
      parts=line.split()
      parts[0] = re.sub("\d+", "", parts[0])
      oldatom=re.sub("\d+", "", data[ln].split()[0])
      if parts[0]!=oldatom:
        print ' WARNING: Different atoms in line %i!' % (ln+1)
        #sys.exit(2)
      parts[0]+='%s' % (atom)
      atom+=1
      line=' '.join(parts)
      molcasin.write(line+'\n')
  
  molcasin.write('End of Input\n')
  molcasin.close()
  print "File %s written."%molcasin_name
  
def run_seward(molcasin_name):
    import subprocess
    
    print "Running seward ..."
    os.environ['Project']='double'
    os.environ['WorkDir']=os.getcwd()
    os.environ['ThisDir']=os.getcwd()
    
    mout = open('molcas.output.seward', 'w')
    retval = subprocess.call(['molcas', molcasin_name], stdout=mout)
    mout.close()
    
    print "   Finished with returncode: %i\n"%retval
    if not retval==0:
      sys.exit(retval)
    
    os.symlink('%s/%s.RunFile'%(os.environ['WorkDir'], os.environ['Project']), 'RUNFILE')
    os.symlink('%s/%s.OneInt'%(os.environ['WorkDir'], os.environ['Project']),  'ONEINT')
    
def write_ciovin(olddir, newdir, ciovin='cioverlap.input'):
  print "Creating file %s"%ciovin
  if os.path.exists(ciovin):
    print "File already exists, exiting ..."
    return
    
  cin = open(ciovin, 'w')
  
  cin.write("a_mo=%s/MOCOEFS/mocoef_mc.sp\n"%olddir)
  cin.write("b_mo=%s/MOCOEFS/mocoef_mc.sp\n"%newdir)
  cin.write("a_det=%s/WORK/dets\n"%olddir)
  cin.write("b_det=%s/WORK/dets\n"%newdir)
  cin.write("ao_read=1\n")

  cin.close()
  print "File %s written."%ciovin

# ========================== Start of Code =============================== #

if __name__ == '__main__':

  if len(sys.argv)<3 or len(sys.argv)>4:
    print numargerror
    quit()
  
  olddir=sys.argv[1]
  newdir=sys.argv[2]
  molcasin_name='molcas.input.double'
  
  write_molcasin(olddir, newdir, molcasin_name)
  
  if (len(sys.argv)==4) and (sys.argv[3]=='run'):
    run_seward(molcasin_name)
    
  write_ciovin(olddir, newdir)
  
  sys.exit(0)