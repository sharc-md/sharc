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

numargerror="""Usage:
python dalton_double-mol.py <dir1> <dir2> [run]
"""

def write_daltcomm(maxpri=10):
    dc = open('daltcomm', 'w')
    
    dc.write("""**DALTONINPUT
.INTEGRALS
.PRINT
    2
**INTEGRALS
.PRINT
    2
.NOSUP
.NOTWO
*READIN
.MAXPRI
   %i
**END OF INPUTS\n"""%maxpri)
    
    dc.close()
    
    print "File daltcomm written."
    
def write_daltaoin(olddir, newdir):
    dout=open('daltaoin', 'w')
    dold=open('%s/WORK/daltaoin' % olddir,'r')
    
    for i, line in enumerate(dold):
        if not i == 3:
            dout.write(line)
        else:
            words = line.split()
            num = int(words[1]) * 2
            dout.write("s  %2i    0           0.10D-14\n"%num)
    
    dold.close()
    
    for i, line in enumerate(open('%s/WORK/daltaoin' % newdir,'r')):
        if i > 3: dout.write(line)
    dout.close()
    print "File daltaoin written."
    
def run_dalton(mem=1000):
    import subprocess
    
    command = ["%s/dalton.x"%os.environ['COLUMBUS'], "-m", "%i"%mem]
    print "command:", command
    
    hout = open('double_hermitls', 'w')
    retval = subprocess.call(command, stdout=hout)
    hout.close()
    
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
    cin.write("ao_read=2\n")
  
    cin.close()
    print "File %s written."%ciovin
    
# ========================== Start of Code =============================== #

if __name__ == '__main__':
    import sys

    if len(sys.argv)<3 or len(sys.argv)>4:
      print numargerror
      quit()
    
    olddir=sys.argv[1]
    newdir=sys.argv[2]
  
    write_daltcomm(maxpri=10)
    write_daltaoin(olddir, newdir)
    
    if (len(sys.argv)==4) and (sys.argv[3]=='run'):
        run_dalton()

    
    write_ciovin(olddir, newdir)