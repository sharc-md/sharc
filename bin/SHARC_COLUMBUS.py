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

#  ====================================================================
#||                                                                    ||
#||                             General Remarks                        ||
#||                                                                    ||
#  ====================================================================
#
# This script uses several different specification for the electronic states under consideration.
# Generally, the input specs are like "3 Singlets, 0 Doublets and 3 Triplets"
# Based on this information, the states may be denoted by different systems.
#
# The most comprehensive denotation is: (mult, state, MS).
# In this case, states of higher multiplicity show up several times and the total number of states may be quite large.
# This total number of states is called nmstates in the script (equals to 12 for the above example)
#
# Since the MS components of a multiplet often share several expectation values, these need to be calculated only once.
# In this case, states can be safely denoted by (mult,state).
# The total number of these states is called nstates.
#
# In both systems, the states can be given indices. In this script and in SHARC, one first counts up the state quantum number,
# then the MS quantum number, and finally the multiplicity, like in this code snippet:
#
# i=0
# for mult in range(len(states)):
#   for MS in range(mult):
#     for state in range(states[mult]):
#       i+=1
#       print i, mult+1, state+1, MS-i/2
#
# more about this below in the docstrings of the iterator functions

# ======================================================================= #

# IMPLEMENTATION OF ADDITIONAL TASKS KEYWORDS, JOBS, ETC:
#
# A new task keyword in QMin has to be added to:
#       - readQMin (for consistency check)
#       - gettasks (planning the MOLPRO calculation)
#       - print QMin (optional)
#
# A new task in the Tasks array needs changes in:
#       - gettasks
#       - writeMOLPROinput
#       - redotasks
#       - printtasks

# ======================================================================= #
# Modules:
# Operating system, isfile and related routines, move files, create directories
import os
# External Calls to MOLPRO
import subprocess as sp
# Command line arguments
import sys
# shell utilities like copy
import shutil
# Regular expressions
import re
# debug print for dicts and arrays
import pprint
# sqrt and other math
import math
import cmath
# runtime measurement
import datetime
# copy of arrays of arrays
from copy import deepcopy
from socket import gethostname



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


changelogstring='''
28.05.2014:     0.3 MAJOR UPDATE
- added ionization via Dyson norms
- can use several COLUMBUS input directories
- storage is divided into SAVE and KEEP directories
- restartable
- molden keyword

03.06.2014:     0.3.01
- fixed bug when doing socinr calculations with e.g. quintets but triplets missing
- dyson output is saved with backup keyword

16.07.2014:     0.3.02
- fixed a bug where the savedir was read in lowercase from the QM.in file

25.08.2014:     0.3.03
- fixed a bug where the order of states from socinr jobs is not correct in the SO-Hamiltonian

26.08.2014:     0.3.04
- savedir and scratchdir now have defaults ($PWD/SAVEDIR/ and $PWD/SCRATCHDIR/)
- added keyword "always_orb_init" to SH2COL.inp, to force the use of the provided initial MOs in all timesteps
- the maximum multiplicity actually needed for SOCI calculations is now written to cidrtin
  (the value in cidrtin could be higher than needed, which will break runc)
- links to savedir and scratchdir are set now

09.09.2014:     0.3.05
- multiple mocoef_mc.init files can be used (mocoef_mc.init.<job> is used for <job>, if available, otherwise mocoef_mc.init)

23.09.2014:     0.3.06
- data is not saved for "samestep" runs. This would make the wavefunction phases in the determinant files inconsistent with previous calculations at the same geometry.

08.10.2014:     1.0
- official release version, no changes to 0.3.06

20.02.2015:
- fixed a bug with obtaining the hostname, and updated the debug facilities

01.04.2015:     1.1
- moved to new cioverlap code of Felix Plasser, Matthias Ruckenbauer, Sebastian Mai
- Dalton integrals can be used to run jobs without SOC

16.04.2015:
- cioverlap.x default path is $SHARC/cioverlap.x

10.05.2015:
- keyword "rasscf" in runc can now be used (keyword "molcas_rasscf" in SH2COL.inp) to do the orbital optimization with MOLCAS.

11.05.2015:     1.1.1
- New Dyson code uses same determinant format as cioverlap now. Both codes employ the same wfthres value.
- cioverlap, dyson output, molden files and runls are backupped now

22.06.2015:     1.1.2
- the same code (wfoverlap) is now used for overlaps and Dyson norms
- ciudgin thresholds are now used unaltered (except that the last threshold is reused for additional states)

01.11.2015:
- fixed "restart" keyword
- added support for BASLIB keyword in Seward

07.03.2017:
- wfthres is now interpreted as in other interfaces (default is 0.97 now)

23.08.2017:
- Dyson norms are now correctly scaled depending on Ms value
- Resource file is now called "COLUMBUS.resources" instead of "SH2COL.inp" (old filename still works)

24.08.2017:
- "always_orb_init" and "always_guess" are now keywords in resources file, "always_mocoef_init" and "always_scf" are deprecated
- "numfrozcore" and "numocc" are now keywords, "ncore" is deprecated ("ndocc" was not functional previously)

01.02.2018:    2.0
- RELEASE version
'''

# ======================================================================= #
# holds the system time when the script was started
starttime=datetime.datetime.now()

# global variables for printing (PRINT gives formatted output, DEBUG gives raw output)
DEBUG=False
PRINT=True

# hash table for conversion of multiplicity to the keywords used in COLUMBUS
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

# hash table for conversion of polarisations to the keywords used in COLUMBUS
IToPol={
        0: 'X',
        1: 'Y',
        2: 'Z',
        'X': 0,
        'Y': 1,
        'Z': 2
        }

# conversion factors
au2a=0.529177211
rcm_to_Eh=4.556335e-6

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
def isbinary(path):
  return (re.search(r':.* text',sp.Popen(["file", '-L', path], stdout=sp.PIPE).stdout.read())is None)




# ======================================================================= #
def eformat(f, prec, exp_digits):
  '''Formats a float f into scientific notation with prec number of decimals and exp_digits number of exponent digits.

  String looks like:
  [ -][0-9]\.[0-9]*E[+-][0-9]*

  Arguments:
  1 float: Number to format
  2 integer: Number of decimals
  3 integer: Number of exponent digits

  Returns:
  1 string: formatted number'''

  s = "% .*e"%(prec, f)
  mantissa, exp = s.split('e')
  return "%sE%+0*d"%(mantissa, exp_digits+1, int(exp))

# ======================================================================= #
def measuretime():
  '''Calculates the time difference between global variable starttime and the time of the call of measuretime.

  Prints the Runtime, if PRINT or DEBUG are enabled.

  Arguments:
  none

  Returns:
  1 float: runtime in seconds'''

  endtime=datetime.datetime.now()
  runtime=endtime-starttime
  if PRINT or DEBUG:
    hours=runtime.seconds/3600
    minutes=runtime.seconds/60-hours*60
    seconds=runtime.seconds%60
    print '==> Runtime:\n%i Days\t%i Hours\t%i Minutes\t%i Seconds\n\n' % (runtime.days,hours,minutes,seconds)
  total_seconds=runtime.days*24*3600+runtime.seconds+runtime.microseconds/1.e6
  return total_seconds





# ======================================================================= #
def itmult(states):
  '''Takes an array of the number of states in each multiplicity and generates an iterator over all multiplicities with non-zero states.

  Example:
  [3,0,3] yields two iterations with
  1
  3

  Arguments:
  1 list of integers: States specification

  Returns:
  1 integer: multiplicity'''

  for i in range(len(states)):
    if states[i]<1:
      continue
    yield i+1
  return

# ======================================================================= #
def itnstates(states):
  '''Takes an array of the number of states in each multiplicity and generates an iterator over all states specified. Different values of MS for each state are not taken into account.

  Example:
  [3,0,3] yields six iterations with
  1,1
  1,2
  1,3
  3,1
  3,2
  3,3

  Arguments:
  1 list of integers: States specification

  Returns:
  1 integer: multiplicity
  2 integer: state'''

  for i in range(len(states)):
    if states[i]<1:
      continue
    for j in range(states[i]):
      yield i+1,j+1
  return

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

  for i in range(len(states)):
    if states[i]<1:
      continue
    for k in range(i+1):
      for j in range(states[i]):
        yield i+1,j+1,k-i/2.
  return

# ======================================================================= #
def ittwostates(states):
  '''Takes an array of the number of states in each multiplicity and generates an iterator over all pairs of states (s1/=s2 and s1<s2), which have the same multiplicity. Different values of MS for each state are not taken into account.

  Example:
  [3,0,3] yields six iterations with
  1 1 2
  1 1 3
  1 2 3
  3 1 2
  3 1 3
  3 2 3

  Arguments:
  1 list of integers: States specification

  Returns:
  1 integer: multiplicity
  2 integer: state 1
  3 integer: state 2'''

  for i in range(len(states)):
    if states[i]<2:
      continue
    for j1 in range(states[i]):
      for j2 in range(j1+1,states[i]):
        yield i+1,j1+1,j2+1
  return

# ======================================================================= #
def ittwostatesfull(states):
  '''Takes an array of the number of states in each multiplicity and generates an iterator over all pairs of states (all combinations), which have the same multiplicity. Different values of MS for each state are not taken into account.

  Example:
  [3,0,3] yields 18 iterations with
  1 1 1
  1 1 2
  1 1 3
  1 2 1
  1 2 2
  1 2 3
  1 3 1
  1 3 2
  1 3 3
  3 1 1
  3 1 2
  3 1 3
  3 2 1
  3 2 2
  3 2 3
  3 3 1
  3 3 2
  3 3 3

  Arguments:
  1 list of integers: States specification

  Returns:
  1 integer: multiplicity
  2 integer: state 1
  3 integer: state 2'''

  for i in itmult(states):
    for j in itnstates([states[i-1]]):
      for k in itnstates([states[i-1]]):
        yield i,j[1],k[1]
  return

# ======================================================================= #
def IstateToMultState(i,states):
  '''Takes a state index in nmstates counting scheme and converts it into (mult,state).

  Arguments:
  1 integer: state index
  2 list of integers: states specification

  Returns:
  1 integer: mult
  2 integer: state'''

  j=i
  for mult,state,ms in itnmstates(states):
    i-=1
    if i==0:
      return mult,state
  print 'state %i is not in states:' % (j),states
  sys.exit(14)

# ======================================================================= #
def IstateToMultStateMS(i,states):
  '''Takes a state index in nmstates counting scheme and converts it into (mult,state).

  Arguments:
  1 integer: state index
  2 list of integers: states specification

  Returns:
  1 integer: mult
  2 integer: state'''

  for mult,state,ms in itnmstates(states):
    i-=1
    if i==0:
      return mult,state,ms
  print 'state %i is not in states:',states
  sys.exit(15)

# ======================================================================= #
def MultStateMSToIstateCOL(mult,state,ms,states):
  '''Takes a tuple (mult,state) and returns all indices i in nmstates scheme, which correspond to this mult and state.

  Arguments:
  1 integer: mult
  2 integer: state
  3 list of integers: states specification

  Returns:
  1 integer: state index in nmstates scheme'''

  if mult-1>len(states) or state>states[mult-1]:
    print 'No state %i, mult %i in states:' % (state,mult),states
    sys.exit(16)
  n=1
  for i in range(len(states)):
    if states[i]<1:
      continue
    for j in range(states[i]):
      for k in range(i+1):
        if i+1==mult and j+1==state and k-i/2.==ms:
          return n
        n+=1

# ======================================================================= #
def MultStateToIstate(mult,state,states):
  '''Takes a tuple (mult,state) and returns all indices i in nmstates scheme, which correspond to this mult and state.

  Arguments:
  1 integer: mult
  2 integer: state
  3 list of integers: states specification

  Returns:
  1 integer: state index in nmstates scheme'''

  if mult-1>len(states) or state>states[mult-1]:
    print 'No state %i, mult %i in states:' % (state,mult),states
    sys.exit(17)
  i=1
  for imult,istate,ms in itnmstates(states):
    if imult==mult and istate==state:
      yield i
    i+=1
  return

# ======================================================================= #
def MultStateToIstateJstate(mult,state1,state2,states):
  '''Takes (mult,state1,state2) and returns all index tuples (i,j) in nmstates scheme, which correspond to this mult and pair of states. Only returns combinations, where both states have the same MS value.

  Arguments:
  1 integer: mult
  2 integer: state1
  3 integer: state2
  4 list of integers: states specification

  Returns:
  1 integer: state1 index in nmstates scheme
  2 integer: state2 index in nmstates scheme'''

  if mult-1>len(states) or state1>states[mult-1] or state2>states[mult-1]:
    print 'No states %i, %i mult %i in states:' % (state1,state2,mult),states
    sys.exit(18)
  i=1
  k=-1
  for imult,istate,ms in itnmstates(states):
    if imult==mult and istate==state1:
      j=1
      for jmult,jstate,ms2 in itnmstates(states):
        if jmult==mult and jstate==state2:
          k+=1
          if k%(mult+1)==0:
            yield i,j
        j+=1
    i+=1
  return

# ======================================================================= #
def removekey(d,key):
  '''Removes an entry from a dictionary and returns the dictionary.

  Arguments:
  1 dictionary
  2 anything which can be a dictionary keyword

  Returns:
  1 dictionary'''

  if key in d:
    r = dict(d)
    del r[key]
    return r
  return d

# ======================================================================= #     OK
def containsstring(string,line):
  '''Takes a string (regular expression) and another string. Returns True if the first string is contained in the second string.

  Arguments:
  1 string: Look for this string
  2 string: within this string

  Returns:
  1 boolean'''

  a=re.search(string,line)
  if a:
    return True
  else:
    return False


# =============================================================================================== #
# =============================================================================================== #
# =========================================== print routines ==================================== #
# =============================================================================================== #
# =============================================================================================== #

# ======================================================================= #
def printheader():
  '''Prints the formatted header of the log file. Prints version number and version date

  Takes nothing, returns nothing.'''

  print starttime,gethostname(),os.getcwd()
  if not PRINT:
    return
  string='\n'
  string+='  '+'='*80+'\n'
  string+='||'+' '*80+'||\n'
  string+='||'+' '*25+'SHARC - COLUMBUS 7 - Interface'+' '*25+'||\n'
  string+='||'+' '*80+'||\n'
  string+='||'+' '*29+'Author: Sebastian Mai'+' '*30+'||\n'
  string+='||'+' '*80+'||\n'
  string+='||'+' '*(36-(len(version)+1)/2)+'Version: %s' % (version)+' '*(35-(len(version))/2)+'||\n'
  lens=len(versiondate.strftime("%d.%m.%y"))
  string+='||'+' '*(37-lens/2)+'Date: %s' % (versiondate.strftime("%d.%m.%y"))+' '*(37-(lens+1)/2)+'||\n'
  string+='||'+' '*80+'||\n'
  string+='  '+'='*80+'\n\n'
  print string
  if DEBUG:
    print changelogstring

# ======================================================================= #
def printQMin(QMin):
  '''If PRINT, prints a formatted Summary of the control variables read from the input file.

  Arguments:
  1 dictionary: QMin'''

  if DEBUG:
    pprint.pprint(QMin)
  if not PRINT:
    return
  print '==> QMin Job description for:\n%s' % (QMin['comment'])
  string='Tasks:'
  if 'h' in QMin:
    string+='\tH'
  if 'soc' in QMin:
    string+='\tSOC'
  if 'dm' in QMin:
    string+='\tDM'
  if 'grad' in QMin:
    string+='\tGrad'
  if 'nacdr' in QMin:
    string+='\tNac(ddr)'
  if 'nacdt' in QMin:
    string+='\tNac(ddt)'
  if 'overlap' in QMin:
    string+='\tOverlaps'
  if 'angular' in QMin:
    string+='\tAngular'
  if 'ion' in QMin:
    string+='\tDyson norms'
  if 'phases' in QMin:
    string+='\tPhases'
  print string
  string='States: '
  for i in itmult(QMin['states']):
    string+='\t%i %s' % (QMin['states'][i-1],IToMult[i])
  print string
  string='Found Geo'
  if 'veloc' in QMin:
    string+=' and Veloc! '
  else:
    string+='! '
  string+='NAtom is %i.\n' % (QMin['natom'])
  print string
  string=''
  for i in range(QMin['natom']):
    string+='%s ' % (QMin['geo'][i][0])
    for j in range(3):
      string+='% 7.4f ' % (QMin['geo'][i][j+1])
    string+='\n'
  print string
  if 'veloc' in QMin:
    string=''
    for i in range(QMin['natom']):
      string+='%s ' % (QMin['geo'][i][0])
      for j in range(3):
        string+='% 7.4f ' % (QMin['veloc'][i][j])
      string+='\n'
    print string
  if 'grad' in QMin:
    string='Gradients:   '
    for i in range(1,QMin['nmstates']+1):
      if i in QMin['grad']:
        string+='X '
      else:
        string+='. '
    string+='\n'
    print string
  if 'nacdr' in QMin:
    string='Non-adiabatic couplings:\n'
    for i in range(1,QMin['nmstates']+1):
      for j in range(1,QMin['nmstates']+1):
        if [i,j] in QMin['nacdr'] or [j,i] in QMin['nacdr']:
          string+='X '
        else:
          string+='. '
      string+='\n'
    print string
  #if 'overlap' in QMin:
    #string='Overlaps:\n'
    #for i in range(1,QMin['nmstates']+1):
      #for j in range(1,QMin['nmstates']+1):
        #if [i,j] in QMin['overlap'] or [j,i] in QMin['overlap']:
          #string+='X '
        #else:
          #string+='. '
      #string+='\n'
    #print string
  #if 'ion' in QMin:
    #string='Dyson norms: '
    #for i in range(1,QMin['nmstates']+1):
      #if i in QMin['ion']:
        #string+='X '
      #else:
        #string+='. '
    #string+='\n'
    #print string
  for i in QMin:
    if not any( [i==j for j in ['h','dm','soc','geo','veloc','states','comment','LD_LIBRARY_PATH', 'grad','nacdr','ion','overlap'] ] ):
      if not any( [i==j for j in ['ionlist','ionmap'] ] ) or DEBUG:
        string=i+': '
        string+=str(QMin[i])
        print string
  print '\n'
  sys.stdout.flush()

# ======================================================================= #
def printtasks(tasks):
  '''If PRINT, prints a formatted table of the tasks in the tasks list.

  Arguments:
  1 list of lists: tasks list (see gettasks for specs)'''

  if DEBUG:
    pprint.pprint(tasks)
  if not PRINT:
    return
  print '==> Task Queue:\n'
  for i in range(len(tasks)):
    task=tasks[i]
    if task[0]=='movetoold':
      print 'Move SAVE to old'
    if task[0]=='backupdata':
      print 'Backup data\t%s' % (task[1])
    elif task[0]=='save_data':
      print 'Save data\t%s' % (task[1])
    elif task[0]=='keep_data':
      print 'Keep data\t%s' % (task[1])
    elif task[0]=='get_COLout':
      print 'Parse Output\t%s' % (task[1])
    elif task[0]=='cleanup':
      print 'Clean directory\t%s' % (task[1])
    elif task[0]=='mkdir':
      print 'Make directory\t%s' % (task[1])
    elif task[0]=='link':
      print 'Link\t\t%s\n\t--> \t%s' % (task[2],task[1])
    elif task[0]=='write':
      print 'Write file\t%s' % (task[1])
    elif task[0]=='getmo':
      print 'Get mocoef file\t%s' % (task[1])
    elif task[0]=='getinput':
      print 'Get input files\t%s' % (task[1])
    elif task[0]=='getmcinput':
      print 'Get MCSCF input\t%s' % (task[1])
    elif task[0]=='checkinput':
      print 'Check input files'
    elif task[0]=='writegeom':
      print 'Write geometry\nCall xyz2col.x'
    elif task[0]=='runc':
      print 'Run COLUMBUS'
    elif task[0]=='writemolcas':
      print 'Write MOLCAS AO overlap integral input\n\tfrom\t%s\n\tand\t%s\n\tto\t%s' % (task[1],task[2],task[3])
    elif task[0]=='runmolcas':
      print 'Run MOLCAS'
    elif task[0]=='runmkciovinp':
      print 'Run mkciovinp'
    elif task[0]=='runcioverlap':
      print 'Run cioverlap'
    elif task[0]=='saveexcitlistfile':
      print 'Save excitlistfile\tDRT%i' % (task[1])
    elif task[0]=='keep_excitlists':
      print 'Get excitlistfiles from \n\t\t%s' % (task[1])
    elif task[0]=='get_CIOoutput':
      print 'Parse Output \t%s, DRT %i' % (task[1],task[2])
    elif task[0]=='make_dets':
      print 'Run civecconsolidate'
    elif task[0]=='dyson':
      print 'Dyson norm\tMult %i, State %i\tand\tMult %i, State %i' % tuple(task[1])
    else:
      print task
  print '\n'
  sys.stdout.flush()

# ======================================================================= #
def printcomplexmatrix(matrix,states):
  '''Prints a formatted matrix. Zero elements are not printed, blocks of different mult and MS are delimited by dashes. Also prints a matrix with the imaginary parts, if any one element has non-zero imaginary part.

  Arguments:
  1 list of list of complex: the matrix
  2 list of integers: states specs'''

  nmstates=0
  for i in range(len(states)):
    nmstates+=states[i]*(i+1)
  string='Real Part:\n'
  string+='-'*(11*nmstates+nmstates/3)
  string+='\n'
  istate=0
  for imult,i,ms in itnmstates(states):
    jstate=0
    string+='|'
    for jmult,j,ms2 in itnmstates(states):
      if matrix[istate][jstate].real==0.:
        string+=' '*11
      else:
        string+='% .3e ' % (matrix[istate][jstate].real)
      if j==states[jmult-1]:
        string+='|'
      jstate+=1
    string+='\n'
    if i==states[imult-1]:
      string+='-'*(11*nmstates+nmstates/3)
      string+='\n'
    istate+=1
  print string
  imag=False
  string='Imaginary Part:\n'
  string+='-'*(11*nmstates+nmstates/3)
  string+='\n'
  istate=0
  for imult,i,ms in itnmstates(states):
    jstate=0
    string+='|'
    for jmult,j,ms2 in itnmstates(states):
      if matrix[istate][jstate].imag==0.:
        string+=' '*11
      else:
        imag=True
        string+='% .3e ' % (matrix[istate][jstate].imag)
      if j==states[jmult-1]:
        string+='|'
      jstate+=1
    string+='\n'
    if i==states[imult-1]:
      string+='-'*(11*nmstates+nmstates/3)
      string+='\n'
    istate+=1
  string+='\n'
  if imag:
    print string

# ======================================================================= #
def printgrad(grad,natom,geo):
  '''Prints a gradient or nac vector. Also prints the atom elements. If the gradient is identical zero, just prints one line.

  Arguments:
  1 list of list of float: gradient
  2 integer: natom
  3 list of list: geometry specs'''

  string=''
  iszero=True
  for atom in range(natom):
    string+='%i\t%s\t' % (atom+1,geo[atom][0])
    for xyz in range(3):
      if grad[atom][xyz]!=0:
        iszero=False
      string+='% .5f\t' % (grad[atom][xyz])
    string+='\n'
  if iszero:
    print '\t\t...is identical zero...\n'
  else:
    print string

# ======================================================================= #
def printQMout(QMin,QMout):
  '''If PRINT, prints a summary of all requested QM output values. Matrices are formatted using printcomplexmatrix, vectors using printgrad.

  Arguments:
  1 dictionary: QMin
  2 dictionary: QMout'''

  if DEBUG:
    pprint.pprint(QMout)
  if not PRINT:
    return
  states=QMin['states']
  nstates=QMin['nstates']
  nmstates=QMin['nmstates']
  natom=QMin['natom']
  print '===> Results:\n'
  # Hamiltonian matrix, real or complex
  if 'h' in QMin or 'soc' in QMin:
    eshift=math.ceil(QMout['h'][0][0].real)
    print '=> Hamiltonian Matrix:\nDiagonal Shift: %9.2f' % (eshift)
    matrix=deepcopy(QMout['h'])
    for i in range(nmstates):
      matrix[i][i]-=eshift
    printcomplexmatrix(matrix,states)
  # Dipole moment matrices
  if 'dm' in QMin:
    print '=> Dipole Moment Matrices:\n'
    for xyz in range(3):
      print 'Polarisation %s:' % (IToPol[xyz])
      matrix=QMout['dm'][xyz]
      printcomplexmatrix(matrix,states)
  # Gradients
  if 'grad' in QMin:
    print '=> Gradient Vectors:\n'
    istate=0
    for imult,i,ms in itnmstates(states):
      print '%s\t%i\tMs= % .1f:' % (IToMult[imult],i,ms)
      printgrad(QMout['grad'][istate],natom,QMin['geo'])
      istate+=1
  # Non-adiabatic couplings
  if 'nacdt' in QMin:
    print '=> Numerical Non-adiabatic couplings:\n'
    matrix=QMout['nacdt']
    printcomplexmatrix(matrix,states)
    matrix=deepcopy(QMout['mrcioverlap'])
    for i in range(nmstates):
      for j in range(nmstates):
        matrix[i][j]=complex(matrix[i][j])
    print '=> MRCI overlaps:\n'
    printcomplexmatrix(matrix,states)
    if 'phases' in QMout:
      print '=> Wavefunction Phases:\n'
      for i in range(nmstates):
        print '% 3.1f % 3.1f' % (QMout['phases'][i].real,QMout['phases'][i].imag)
      print '\n'
  if 'nacdr' in QMin:
    print '=> Analytical Non-adiabatic coupling vectors:\n'
    istate=0
    for imult,i,msi in itnmstates(states):
      jstate=0
      for jmult,j,msj in itnmstates(states):
        if imult==jmult and msi==msj:
          print '%s\tStates %i - %i\tMs= % .1f:' % (IToMult[imult],i,j,msi)
          printgrad(QMout['nacdr'][istate][jstate],natom,QMin['geo'])
        jstate+=1
      istate+=1
  if 'overlap' in QMin:
    print '=> Overlap matrix:\n'
    matrix=QMout['overlap']
    printcomplexmatrix(matrix,states)
    if 'phases' in QMout:
      print '=> Wavefunction Phases:\n'
      for i in range(nmstates):
        print '% 3.1f % 3.1f' % (QMout['phases'][i].real,QMout['phases'][i].imag)
      print '\n'
  # Angular momentum matrices
  if 'angular' in QMin:
    print '=> Angular Momentum Matrices:\n'
    for xyz in range(3):
      print 'Polarisation %s:' % (IToPol[xyz])
      matrix=QMout['angular'][xyz]
      printcomplexmatrix(matrix,states)
  # Property matrix (dyson norms)
  if 'ion' in QMin:
    print '=> Property matrix:\n'
    matrix=QMout['prop']
    printcomplexmatrix(matrix,states)
  sys.stdout.flush()



# =============================================================================================== #
# =============================================================================================== #
# =========================================== ????????????????? ================================= #
# =============================================================================================== #
# =============================================================================================== #




# ======================================================================= #     OK
def makecmatrix(a,b):
  '''Initialises a complex axb matrix.

  Arguments:
  1 integer: first dimension
  2 integer: second dimension

  Returns;
  1 list of list of complex'''

  mat=[ [ complex(0.,0.) for i in range(a) ] for j in range(b) ]
  return mat

# ======================================================================= #     OK
def makermatrix(a,b):
  '''Initialises a real axb matrix.

  Arguments:
  1 integer: first dimension
  2 integer: second dimension

  Returns;
  1 list of list of real'''

  mat=[ [ 0. for i in range(a) ] for j in range(b) ]
  return mat





# =============================================================================================== #
# =============================================================================================== #
# =========================================== output extraction ================================= #
# =============================================================================================== #
# =============================================================================================== #

# ======================================================================= #
def getcienergy(runls,QMin,istate):
  mult,state,ms=tuple(QMin['statemap'][istate])
  job,drt=tuple(QMin['multmap'][mult])
  socimode=QMin['socimap'][job]

  ilines=0
  if socimode in [0,1]:         # for socinr or isc keywords
    while ilines<len(runls):
      if 'Starting ciudgav with' in runls[ilines]:
        s=runls[ilines].split()
        if int(s[5])==drt:
          istate=0
          while not 'all zwalks' in runls[ilines]:
            if 'total mr-sdci energy' in runls[ilines]:
              istate+=1
              if istate==state:
                return float(runls[ilines].split()[3])
            ilines+=1
      ilines+=1
  elif socimode==-1:            # for single-drt jobs without extra keyword
    istate=0
    while not 'all zwalks' in runls[ilines]:
      if 'total mr-sdci energy' in runls[ilines]:
        istate+=1
        if istate==state:
          return float(runls[ilines].split()[3])
      ilines+=1
  print 'CI energy of state %i in drt %i not found!' % (state,drt)
  sys.exit(19)

# ======================================================================= #
def istate_in_job(m1,s1,ms1,states):
  k=-1
  for im in range(m1):
    for ist in range(states[im]):
      for ims in range(im+1):
        k+=1
        if (im+1,ims-im/2+1/2,ist+1)==(m1,ms1,s1):
          return k

# ======================================================================= #
def getsocme(runls,QMin,istate,jstate):
  if not istate in QMin['statemap'] or not jstate in QMin['statemap']:
    print 'States %i or %i are not in statemap!' % (istate,jstate)
    sys.exit(20)

  m1,s1,ms1=tuple(QMin['statemap'][istate])
  m2,s2,ms2=tuple(QMin['statemap'][jstate])
  job1=QMin['multmap'][m1][0]
  job2=QMin['multmap'][m2][0]

  if job1!=job2:
    return complex(0.,0.)
  job=job1

  if not QMin['socimap'][job]==1:
    print 'Job %s is not a socinr job!' % (job)
    return complex(0.,0.)

  ## get numbers of the states within this job
  #i=istate-1
  #for mult in range(1,m1):
    #if QMin['states'][mult-1]==0:
      #continue
    #if QMin['multmap'][mult][0]!=job:
      #i-=mult*QMin['states'][mult-1]
  #j=jstate-1
  #for mult in range(1,m2):
    #if QMin['states'][mult-1]==0:
      #continue
    #if QMin['multmap'][mult][0]!=job:
      #j-=mult*QMin['states'][mult-1]

  # get numbers of states within job
  mults=QMin['multmap'][job]
  states_in_job=[ 0 for imult in range(max(mults)) ]
  for imult in mults:
    states_in_job[imult-1]=QMin['states'][imult-1]

  ## get index of state 1
  #k=-1
  #for im in range(m1):
    #for ist in range(states_in_job[im]):
      #for ims in range(im+1):
        #k+=1
        #print im+1,ims-im/2+1/2,ist+1
        #if (im+1,ims-im/2+1/2,ist+1)==(m1,ms1,s1):
          #i=k
          #break
  ## get index of state 2
  #k=-1
  #for im in range(m1):
    #for ist in range(states_in_job[im]):
      #for ims in range(im+1):
        #k+=1
        #print im+1,ims-im/2+1/2,ist+1
        #if (im+1,ims-im/2+1/2,ist+1)==(m2,ms2,s2):
          #j=k
          #break

  i=istate_in_job(m1,s1,ms1,states_in_job)
  j=istate_in_job(m2,s2,ms2,states_in_job)

  # find reference energy (simply the first mrci energy encountered)
  ilines=0
  while not ' total mr-sdci energy'in runls[ilines]:
    ilines+=1
  eref=float(runls[ilines].split()[3])

  # find the matrix
  while not '-------------- HT matrix in cm-1 -------------------' in runls[ilines]:
    ilines+=1
    if ilines==len(runls):
      print 'SO Matrix not found!'
      sys.exit(21)
  ilines+=2

  # find a matrix element
  if i>=j:
    # use j as col and i as row
    f=runls[ilines+i].split()
    real=float(f[j])*rcm_to_Eh
  else:
    # use i as col and j as row
    f=runls[ilines+j].split()
    real=float(f[i])*rcm_to_Eh
  imag=0.
  if i==j:
    real+=eref
  return complex(real,imag)

## ======================================================================= #
def getcidm(runls,QMin,istate,jstate,idir):
  m1,s1,ms1=tuple(QMin['statemap'][istate])
  m2,s2,ms2=tuple(QMin['statemap'][jstate])
  job1=QMin['multmap'][m1][0]
  job2=QMin['multmap'][m2][0]

  if job1!=job2:
    return 0.
  job=job1

  if m1!=m2:
    return 0.

  if ms1!=ms2:
    return 0.

  if idir=='X' or idir=='Y' or idir=='Z':
    idir=IToPol[idir]

  drt=QMin['multmap'][m1][1]

  # get a matrix element
  ilines=0
  if s1==s2:
    while ilines<len(runls):
      if ' --- ci properties for state ' in runls[ilines]:
        f=runls[ilines].split()
        if drt==int(f[8]) and s1==int(f[5]):
          while not '-------------------------' in runls[ilines]:
            if 'total' in runls[ilines]:
              f=runls[ilines].split()
              return float(f[idir+1])
            ilines+=1
      ilines+=1

  else:
    while ilines<len(runls):
      if ' Calculating transition moment for' in runls[ilines]:
        f=runls[ilines].replace('(',' ').replace(')',' ').replace(',',' ').split()
        if drt==int(f[4]):
          x1=int(f[5])
          x2=int(f[8])
          if (s1==x1 and s2==x2) or (s2==x1 and s1==x2):
            while not '-------------------------' in runls[ilines]:
              if 'total (elec)' in runls[ilines]:
                f=runls[ilines].split()
                return float(f[idir+2])
              ilines+=1
      ilines+=1

  print 'CI dipole moment of states %i and %i in drt %i not found!' % (s1,s2,drt)
  sys.exit(22)


# ======================================================================= #
def getgrad(runls,QMin,istate):
  mult,state,ms=tuple(QMin['statemap'][istate])
  job,drt=tuple(QMin['multmap'][mult])
  natom=QMin['natom']

  ilines=0
  grad=[]
  while ilines<len(runls):
    if ' === Gradient of DRT' in runls[ilines]:
      f=runls[ilines].replace(',',' ').split()
      if drt==int(f[4]) and state==int(f[6]):
        ilines+=1
        for atom in range(natom):
          line=runls[ilines+atom].replace('D','E').split()
          for i in range(3):
            line[i]=float(line[i])
          grad.append(line)
        return grad
    ilines+=1
  print 'Gradient of state %i in drt %i not found!' % (state,drt)
  sys.exit(23)

# ======================================================================= #
def getnacana(runls,QMin,istate,jstate):
  m1,s1,ms1=tuple(QMin['statemap'][istate])
  m2,s2,ms2=tuple(QMin['statemap'][jstate])
  job1=QMin['multmap'][m1][0]
  job2=QMin['multmap'][m2][0]
  natom=QMin['natom']

  nac=[ [0.,0.,0.] for i in range(natom) ]

  if job1!=job2:
    return nac

  if m1!=m2:
    return nac

  if ms1!=ms2:
    return nac

  if s1==s2:
    return nac

  drt=QMin['multmap'][m1][1]

  # get a non-adiabatic coupling vector
  ilines=0
  nac=[]
  interst = False
  while ilines<len(runls):
    if 'deltae=' in runls[ilines]:
      ideltae = -1./float(runls[ilines].replace('=', ' ').split()[-1])
    if 'coupling for DRT' in runls[ilines]:
      if 'Interstate' in  runls[ilines]:
        interst = True
      else:
        interst = False
      f=runls[ilines].replace(':',' ').split()
      if drt==int(f[5]):
        if (s1==int(f[7]) and s2==int(f[12])) or (s2==int(f[7]) and s1==int(f[12])):
          ilines+=1
          for atom in range(natom):
            line=runls[ilines+atom].replace('D','E').split()
            for i in range(3):
              line[i]=float(line[i])
              if interst:
                line[i]*=ideltae
              if s1>s2:
                line[i]*=-1.
            nac.append(line)
          return nac
    ilines+=1
  print 'Interstate coupling of states %i and %i in drt %i not found!' % (s1,s2,drt)
  sys.exit(24)

# ======================================================================= #
def get_COLout(QMin,QMout,job):
  '''Updates QMout by reading the runls of job'''

  #with open(QMin['scratchdir']+'/JOB/runls') as f:
    #runls=f.readlines()
  runls=readfile(QMin['scratchdir']+'/JOB/runls')

  if 'backup' in QMin:
    f=QMin['scratchdir']+'/JOB/runls'
    fdest=QMin['backup']+'/runls.%s' % (job.replace('/','_'))
    shutil.copy(f,fdest)

  states=QMin['states']
  nstates=QMin['nstates']
  nmstates=QMin['nmstates']
  natom=QMin['natom']

  # Hamiltonian
  if 'h' in QMin or 'soc' in QMin:
    # make Hamiltonian
    if not 'h' in QMout:
      QMout['h']=makecmatrix(nmstates,nmstates)
    # read the matrix elements from runls
    for i in range(nmstates):
      for j in range(nmstates):
        m1,s1,ms1=tuple(QMin['statemap'][i+1])
        m2,s2,ms2=tuple(QMin['statemap'][j+1])
        if not QMin['multmap'][m1][0]==QMin['multmap'][m2][0]==job:
          continue
        if 'soc' in QMin and QMin['socimap'][job]==1:
          # read SOC elements
          QMout['h'][i][j]=getsocme(runls,QMin,i+1,j+1)
        else:
          # read only diagonal elements
          if i==j:
            QMout['h'][i][j]=getcienergy(runls,QMin,i+1)

  # Dipole matrices
  if 'dm' in QMin:
    # make matrix
    if not 'dm' in QMout:
      QMout['dm']=[]
      for xyz in range(3):
        QMout['dm'].append(makecmatrix(nmstates,nmstates))
    # read the matrix elements from runls
    for i in range(nmstates):
      for j in range(nmstates):
        m1,s1,ms1=tuple(QMin['statemap'][i+1])
        m2,s2,ms2=tuple(QMin['statemap'][j+1])
        if not QMin['multmap'][m1][0]==QMin['multmap'][m2][0]==job:
          continue
        for idir in range(3):
          QMout['dm'][idir][i][j]=getcidm(runls,QMin,i+1,j+1,idir)

  # Gradients
  # full list, even with selection of gradients
  if 'grad' in QMin:
    # initialize array
    if not 'grad' in QMout:
      QMout['grad']=[ [ [0.,0.,0.] for i in range(natom) ] for j in range(nmstates) ]
    # read gradients of current job
    for i in range(nmstates):
      m,s,ms=tuple(QMin['statemap'][i+1])
      if not QMin['multmap'][m][0]==job:
        continue
      requested=False
      for j in QMin['grad']:
        m2,s2,ms2=tuple(QMin['statemap'][j])
        if m==m2 and s==s2:
          requested=True
      if requested:
        QMout['grad'][i]=getgrad(runls,QMin,i+1)

  # Non-adiabatic couplings
  # full list, even with selection
  if 'nacdr' in QMin:
    # initialize
    if not 'nacdr' in QMout:
      QMout['nacdr']=[ [ [ [0.,0.,0.] for i in range(natom) ] for j in range(nmstates) ] for k in range(nmstates) ]
    # read nac vectors
    for i in range(nmstates):
      for j in range(nmstates):
        m1,s1,ms1=tuple(QMin['statemap'][i+1])
        m2,s2,ms2=tuple(QMin['statemap'][j+1])
        if not ms1==ms2:
          continue
        if not QMin['multmap'][m1][0]==QMin['multmap'][m2][0]==job:
          continue
        requested=False
        for sel in QMin['nacdr']:
          k,l=tuple(sel)
          m3,s3,ms3=tuple(QMin['statemap'][k])
          m4,s4,ms4=tuple(QMin['statemap'][l])
          if m1==m2==m3==m4:
            if (s1,s2)==(s3,s4) or (s1,s2)==(s4,s3):
              requested=True
        if requested:
          QMout['nacdr'][i][j]=getnacana(runls,QMin,i+1,j+1)

  return QMout
# ======================================================================= #
def get_smatel(out,s1,s2):
  ilines=-1
  while True:
    ilines+=1
    if ilines==len(out):
      print 'Overlap of states %i - %i not found!' % (s1,s2)
      sys.exit(25)
    if containsstring('Overlap matrix <PsiA_i|PsiB_j>', out[ilines]):
      break
  ilines+=1+s1
  f=out[ilines].split()
  return float(f[s2+1])

# ======================================================================= #
def get_CIOoutput(QMin,QMout,job,drt):
  #with open(QMin['scratchdir']+'/OVERLAP/cioverlap.out') as f:
    #out=f.readlines()
  out=readfile(QMin['scratchdir']+'/OVERLAP/cioverlap.out')

  if 'overlap' in QMin:
    nmstates=QMin['nmstates']
    if not 'overlap' in QMout:
      QMout['overlap']=makecmatrix(nmstates,nmstates)
    # read the overlap matrix
    mult=QMin['multmap'][job][drt-1]
    for i in range(nmstates):
      for j in range(nmstates):
        m1,s1,ms1=tuple(QMin['statemap'][i+1])
        m2,s2,ms2=tuple(QMin['statemap'][j+1])
        if not m1==m2==mult:
          continue
        if not ms1==ms2:
          continue
        QMout['overlap'][i][j]=get_smatel(out,s1,s2)

  return QMout
# ======================================================================= #
def get_dysonel(out,s1,s2):
  ilines=-1
  while True:
    ilines+=1
    if ilines==len(out):
      print 'Overlap of states %i - %i not found!' % (s1,s2)
      sys.exit(26)
    if containsstring('Dyson norm matrix |<PsiA_i|PsiB_j>|^2', out[ilines]):
      break
  ilines+=1+s1
  f=out[ilines].split()
  return float(f[s2+1])

# ======================================================================= #
def get_DYSONoutput(QMin,QMout,imult,jmult):
  #with open(QMin['scratchdir']+'/OVERLAP/cioverlap.out') as f:
    #out=f.readlines()
  out=readfile(QMin['scratchdir']+'/DYSON/dyson.out')

  if 'ion' in QMin:
    nmstates=QMin['nmstates']
    if not 'prop' in QMout:
      QMout['prop']=makecmatrix(nmstates,nmstates)
    # read the Dyson norm matrix
    for i in range(nmstates):
      for j in range(nmstates):
        m1,s1,ms1=tuple(QMin['statemap'][i+1])
        m2,s2,ms2=tuple(QMin['statemap'][j+1])
        if not (imult,jmult)==(m1,m2) and not (imult,jmult)==(m2,m1):
          continue
        if not abs(ms1-ms2)==0.5:
          continue
        # switch multiplicities such that m1 is smaller mult
        if m1>m2:
          s1,s2=s2,s1
          m1,m2=m2,m1
          ms1,ms2=ms2,ms1
        # compute M_S overlap factor
        if ms1<ms2:
          factor=( ms1+1.+(m1-1.)/2. )/m1
        else:
          factor=( -ms1+1.+(m1-1.)/2. )/m1
        QMout['prop'][i][j]=get_dysonel(out,s1,s2)*factor

  return QMout




# =============================================================================================== #
# =============================================================================================== #
# =========================================== QMout writing ===================================== #
# =============================================================================================== #
# =============================================================================================== #


# ======================================================================= #
def writeQMout(QMin,QMout,QMinfilename):
  '''Writes the requested quantities to the file which SHARC reads in. The filename is QMinfilename with everything after the first dot replaced by "out".

  Arguments:
  1 dictionary: QMin
  2 dictionary: QMout
  3 string: QMinfilename'''

  os.chdir(QMin['pwd'])
  k=QMinfilename.find('.')
  if k==-1:
    outfilename=QMinfilename+'.out'
  else:
    outfilename=QMinfilename[:k]+'.out'
  if PRINT:
    print '===> Writing output to file %s in SHARC Format\n' % (outfilename)
  string=''
  if 'h' in QMin or 'soc' in QMin:
    string+=writeQMoutsoc(QMin,QMout)
  if 'dm' in QMin:
    string+=writeQMoutdm(QMin,QMout)
  if 'angular' in QMin:
    string+=writeQMoutang(QMin,QMout)
  if 'grad' in QMin:
    string+=writeQMoutgrad(QMin,QMout)
  if 'nacdr' in QMin:
    string+=writeQMoutnacana(QMin,QMout)
  if 'overlap' in QMin:
    string+=writeQMoutnacsmat(QMin,QMout)
  if 'ion' in QMin:
    string+=writeQMoutprop(QMin,QMout)
  if 'phases' in QMin:
    string+=writeQmoutPhases(QMin,QMout)
  string+=writeQMouttime(QMin,QMout)
  writefile(outfilename,string)
  #try:
    #outfile=open(outfilename,'w')
    #outfile.write(string)
    #outfile.close()
  #except IOError:
    #print 'Could not write QM output!'
    #sys.exit(27)
  if 'backup' in QMin:
    writefile(QMin['backup']+'/'+outfilename,string)
    #try:
      #outfile=open(QMin['backup']+'/'+outfilename,'w')
      #outfile.write(string)
      #outfile.close()
    #except IOError:
      #print 'WARNING: Could not write QM output backup!'
  return

# ======================================================================= #
def writeQMoutsoc(QMin,QMout):
  '''Generates a string with the Spin-Orbit Hamiltonian in SHARC format.

  The string starts with a ! followed by a flag specifying the type of data. In the next line, the dimensions of the matrix are given, followed by nmstates blocks of nmstates elements. Blocks are separated by a blank line.

  Arguments:
  1 dictionary: QMin
  2 dictionary: QMout

  Returns:
  1 string: multiline string with the SOC matrix'''

  states=QMin['states']
  nstates=QMin['nstates']
  nmstates=QMin['nmstates']
  natom=QMin['natom']
  string=''
  string+='! %i Hamiltonian Matrix (%ix%i, complex)\n' % (1,nmstates,nmstates)
  string+='%i %i\n' % (nmstates,nmstates)
  for i in range(nmstates):
    for j in range(nmstates):
      string+='%s %s ' % (eformat(QMout['h'][i][j].real,12,3),eformat(QMout['h'][i][j].imag,12,3))
    string+='\n'
  string+='\n'
  return string

# ======================================================================= #
def writeQMoutdm(QMin,QMout):
  '''Generates a string with the Dipole moment matrices in SHARC format.

  The string starts with a ! followed by a flag specifying the type of data. In the next line, the dimensions of the matrix are given, followed by nmstates blocks of nmstates elements. Blocks are separated by a blank line. The string contains three such matrices.

  Arguments:
  1 dictionary: QMin
  2 dictionary: QMout

  Returns:
  1 string: multiline string with the DM matrices'''

  states=QMin['states']
  nstates=QMin['nstates']
  nmstates=QMin['nmstates']
  natom=QMin['natom']
  string=''
  string+='! %i Dipole Moment Matrices (3x%ix%i, complex)\n' % (2,nmstates,nmstates)
  for xyz in range(3):
    string+='%i %i\n' % (nmstates,nmstates)
    for i in range(nmstates):
      for j in range(nmstates):
        string+='%s %s ' % (eformat(QMout['dm'][xyz][i][j].real,12,3),eformat(QMout['dm'][xyz][i][j].imag,12,3))
      string+='\n'
    string+=''
  return string

# ======================================================================= #
def writeQMoutang(QMin,QMout):
  '''Generates a string with the Dipole moment matrices in SHARC format.

  The string starts with a ! followed by a flag specifying the type of data. In the next line, the dimensions of the matrix are given, followed by nmstates blocks of nmstates elements. Blocks are separated by a blank line. The string contains three such matrices.

  Arguments:
  1 dictionary: QMin
  2 dictionary: QMout

  Returns:
  1 string: multiline string with the DM matrices'''

  states=QMin['states']
  nstates=QMin['nstates']
  nmstates=QMin['nmstates']
  natom=QMin['natom']
  string=''
  string+='! %i Angular Momentum Matrices (3x%ix%i, complex)\n' % (9,nmstates,nmstates)
  for xyz in range(3):
    string+='%i %i\n' % (nmstates,nmstates)
    for i in range(nmstates):
      for j in range(nmstates):
        string+='%s %s ' % (eformat(QMout['angular'][xyz][i][j].real,12,3),eformat(QMout['angular'][xyz][i][j].imag,12,3))
      string+='\n'
    string+=''
  return string

# ======================================================================= #
def writeQMoutgrad(QMin,QMout):
  '''Generates a string with the Gradient vectors in SHARC format.

  The string starts with a ! followed by a flag specifying the type of data. On the next line, natom and 3 are written, followed by the gradient, with one line per atom and a blank line at the end. Each MS component shows up (nmstates gradients are written).

  Arguments:
  1 dictionary: QMin
  2 dictionary: QMout

  Returns:
  1 string: multiline string with the Gradient vectors'''

  states=QMin['states']
  nstates=QMin['nstates']
  nmstates=QMin['nmstates']
  natom=QMin['natom']
  string=''
  string+='! %i Gradient Vectors (%ix%ix3, real)\n' % (3,nmstates,natom)
  i=0
  for imult,istate,ims in itnmstates(states):
    string+='%i %i ! m1 %i s1 %i ms1 %i\n' % (natom,3,imult,istate,ims)
    for atom in range(natom):
      for xyz in range(3):
        string+='%s ' % (eformat(QMout['grad'][i][atom][xyz],12,3))
      string+='\n'
    string+=''
    i+=1
  return string

# ======================================================================= #
def writeQMoutnacnum(QMin,QMout):
  '''Generates a string with the NAC matrix in SHARC format.

  The string starts with a ! followed by a flag specifying the type of data. In the next line, the dimensions of the matrix are given, followed by nmstates blocks of nmstates elements. Blocks are separated by a blank line.

  Arguments:
  1 dictionary: QMin
  2 dictionary: QMout

  Returns:
  1 string: multiline string with the NAC matrix'''

  states=QMin['states']
  nstates=QMin['nstates']
  nmstates=QMin['nmstates']
  natom=QMin['natom']
  string=''
  string+='! %i Non-adiabatic couplings (ddt) (%ix%i, complex)\n' % (4,nmstates,nmstates)
  string+='%i %i\n' % (nmstates,nmstates)
  for i in range(nmstates):
    for j in range(nmstates):
      string+='%s %s ' % (eformat(QMout['nacdt'][i][j].real,12,3),eformat(QMout['nacdt'][i][j].imag,12,3))
    string+='\n'
  string+=''
  # also write wavefunction phases
  string+='! %i Wavefunction phases (%i, complex)\n' % (7,nmstates)
  for i in range(nmstates):
    string+='%s %s\n' % (eformat(QMout['phases'][i],12,3),eformat(0.,12,3))
  string+='\n\n'
  return string

# ======================================================================= #
def writeQMoutnacana(QMin,QMout):
  '''Generates a string with the NAC vectors in SHARC format.

  The string starts with a ! followed by a flag specifying the type of data. On the next line, natom and 3 are written, followed by the gradient, with one line per atom and a blank line at the end. Each MS component shows up (nmstates x nmstates vectors are written).

  Arguments:
  1 dictionary: QMin
  2 dictionary: QMout

  Returns:
  1 string: multiline string with the NAC vectors'''

  states=QMin['states']
  nstates=QMin['nstates']
  nmstates=QMin['nmstates']
  natom=QMin['natom']
  string=''
  string+='! %i Non-adiabatic couplings (ddr) (%ix%ix%ix3, real)\n' % (5,nmstates,nmstates,natom)
  i=0
  for imult,istate,ims in itnmstates(states):
    j=0
    for jmult,jstate,jms in itnmstates(states):
      #string+='%i %i ! %i %i %i %i %i %i\n' % (natom,3,imult,istate,ims,jmult,jstate,jms)
      string+='%i %i ! m1 %i s1 %i ms1 %i   m2 %i s2 %i ms2 %i\n' % (natom,3,imult,istate,ims,jmult,jstate,jms)
      for atom in range(natom):
        for xyz in range(3):
          string+='%s ' % (eformat(QMout['nacdr'][i][j][atom][xyz],12,3))
        string+='\n'
      string+=''
      j+=1
    i+=1
  return string

# ======================================================================= #
def writeQMoutnacsmat(QMin,QMout):
  '''Generates a string with the adiabatic-diabatic transformation matrix in SHARC format.

  The string starts with a ! followed by a flag specifying the type of data. In the next line, the dimensions of the matrix are given, followed by nmstates blocks of nmstates elements. Blocks are separated by a blank line.

  Arguments:
  1 dictionary: QMin
  2 dictionary: QMout

  Returns:
  1 string: multiline string with the transformation matrix'''

  states=QMin['states']
  nstates=QMin['nstates']
  nmstates=QMin['nmstates']
  natom=QMin['natom']
  string=''
  string+='! %i Overlap matrix (%ix%i, complex)\n' % (6,nmstates,nmstates)
  string+='%i %i\n' % (nmstates,nmstates)
  for j in range(nmstates):
    for i in range(nmstates):
      string+='%s %s ' % (eformat(QMout['overlap'][j][i].real,12,3),eformat(QMout['overlap'][j][i].imag,12,3))
    string+='\n'
  string+='\n'
  return string

# ======================================================================= #
def writeQMouttime(QMin,QMout):
  '''Generates a string with the quantum mechanics total runtime in SHARC format.

  The string starts with a ! followed by a flag specifying the type of data. In the next line, the runtime is given

  Arguments:
  1 dictionary: QMin
  2 dictionary: QMout

  Returns:
  1 string: multiline string with the runtime'''

  string='! 8 Runtime\n%s\n' % (eformat(QMout['runtime'],9,3))
  return string

# ======================================================================= #
def writeQMoutprop(QMin,QMout):
  '''Generates a string with the Spin-Orbit Hamiltonian in SHARC format.

  The string starts with a ! followed by a flag specifying the type of data. In the next line, the dimensions of the matrix are given, followed by nmstates blocks of nmstates elements. Blocks are separated by a blank line.

  Arguments:
  1 dictionary: QMin
  2 dictionary: QMout

  Returns:
  1 string: multiline string with the SOC matrix'''

  states=QMin['states']
  nstates=QMin['nstates']
  nmstates=QMin['nmstates']
  natom=QMin['natom']
  string=''
  string+='! %i Property Matrix (%ix%i, complex)\n' % (11,nmstates,nmstates)
  string+='%i %i\n' % (nmstates,nmstates)
  for i in range(nmstates):
    for j in range(nmstates):
      string+='%s %s ' % (eformat(QMout['prop'][i][j].real,12,3),eformat(QMout['prop'][i][j].imag,12,3))
    string+='\n'
  string+='\n'
  return string

# ======================================================================= #
def writeQmoutPhases(QMin,QMout):

  string='! 7 Phases\n%i ! for all nmstates\n' % (QMin['nmstates'])
  for i in range(QMin['nmstates']):
    string+='%s %s\n' % (eformat(QMout['phases'][i].real,9,3),eformat(QMout['phases'][i].imag,9,3))
  return string




# =============================================================================================== #
# =============================================================================================== #
# =========================================== SUBROUTINES TO readQMin =========================== #
# =============================================================================================== #
# =============================================================================================== #

def checkscratch(SCRATCHDIR):
  '''Checks whether SCRATCHDIR is a file or directory. If a file, it quits with exit code 1, if its a directory, it passes. If SCRATCHDIR does not exist, tries to create it.

  Arguments:
  1 string: path to SCRATCHDIR'''

  exist=os.path.exists(SCRATCHDIR)
  if exist:
    isfile=os.path.isfile(SCRATCHDIR)
    if isfile:
      print '$SCRATCHDIR=%s exists and is a file!' % (SCRATCHDIR)
      sys.exit(28)
  else:
    try:
      os.makedirs(SCRATCHDIR)
    except OSError:
      print 'Can not create SCRATCHDIR=%s\n' % (SCRATCHDIR)
      sys.exit(29)

# ======================================================================= #
def checktemplate(TEMPLATE, mult,states, integrals):
  '''Checks whether TEMPLATE is a file or directory. If a file or does not exist, it quits with exit code 1, if it is a directory, it checks whether all important input files are there. Does not check for all input files, since runc does this, too.

  Arguments:
  1 string: path to TEMPLATE

  returns whether input is for isc keyword or socinr keyword
  and returns the DRT of the given multiplicity'''

  exist=os.path.exists(TEMPLATE)
  if exist:
    isfile=os.path.isfile(TEMPLATE)
    if isfile:
      print 'TEMPLATE=%s exists and is a file!' % (TEMPLATE)
      sys.exit(30)
    necessary=['control.run','mcscfin','tranin']
    if integrals=='seward':
      necessary.append('molcas.input')
    elif integrals=='dalton':
      necessary.append('daltaoin')
    lof=os.listdir(TEMPLATE)
    for i in necessary:
      if not i in lof:
        print 'Did not find input file %s! Did you prepare the input according to the instructions?' % (i)
        sys.exit(31)
    cidrtinthere=False
    ciudginthere=False
    for i in lof:
      if 'cidrtin' in i:
        cidrtinthere=True
      if 'ciudgin' in i:
        ciudginthere=True
    if not cidrtinthere or not ciudginthere:
      print 'Did not find input file %s.*! Did you prepare the input according to the instructions?' % (i)
      sys.exit(32)
  else:
    print 'Directory %s does not exist!' % (TEMPLATE)
    sys.exit(33)

  # check cidrtin and cidrtin* for the multiplicity
  try:
    cidrtin=open(TEMPLATE+'/cidrtin')
    line=cidrtin.readline().split()
    if line[0].lower()=='y':
      maxmult=int(cidrtin.readline().split()[0])
      cidrtin.readline()
      nelec=int(cidrtin.readline().split()[0])
      if mult<=maxmult and (mult+nelec)%2!=0:
        drt=0
        for m1,n in enumerate(states):
          if m1+1>mult:
            break
          if (m1+1+nelec)%2!=0 and n>0:
            drt+=1
        return 1, drt    # socinr=1, single=-1, isc=0
      else:
        print 'Multiplicity %i cannot be treated in directory %s (socinr)!'  % (mult,TEMPLATE)
        sys.exit(34)
    else:
      mult2=int(cidrtin.readline().split()[0])
      if mult!=mult2:
        print 'Multiplicity %i cannot be treated in directory %s (single DRT)!'  % (mult,TEMPLATE)
        sys.exit(35)
      # not a soci calculation, but cidrtin file -> only one DRT
      return -1,1
  except IOError:
    # find out in which DRT the requested multiplicity is
    for i in range(1,9):        # COLUMBUS can treat at most 8 DRTs
      try:
        cidrtin=open(TEMPLATE+'/cidrtin.%i' % i)
      except IOError:
        print 'Multiplicity %i cannot be treated in directory %s (isc)!'  % (mult,TEMPLATE)
        sys.exit(36)
      cidrtin.readline()
      mult2=int(cidrtin.readline().split()[0])
      if mult==mult2:
        return 0,i
      cidrtin.close()

# ======================================================================= #
def removequotes(string):
  if string.startswith("'") and string.endswith("'"):
    return string[1:-1]
  elif string.startswith('"') and string.endswith('"'):
    return string[1:-1]
  else:
    return string

# ======================================================================= #
def getsh2colkey(sh2col,key):
  i=-1
  while True:
    i+=1
    try:
      line=re.sub('#.*$','',sh2col[i])
    except IndexError:
      break
    line=line.split(None,1)
    if line==[]:
      continue
    if key.lower() in line[0].lower():
      return line
  return ['','']

# ======================================================================= #
def get_sh2col_environ(sh2col,key,environ=True,crucial=True):
  line=getsh2colkey(sh2col,key)
  if line[0]:
    LINE=line[1]
  else:
    if environ:
      LINE=os.getenv(key.upper())
      if not LINE:
        print 'Either set $%s or give path to %s in COLUMBUS.resources!' % (key.upper(),key.upper())
        if crucial:
          sys.exit(37)
        else:
          return None
    else:
      print 'Give path to %s in COLUMBUS.resources!' % (key.upper())
      if crucial:
        sys.exit(38)
      else:
        return None
  LINE=os.path.expandvars(LINE)
  LINE=os.path.expanduser(LINE)
  LINE=os.path.abspath(LINE)
  LINE=removequotes(LINE).strip()
  if containsstring(';',LINE):
    print "$%s contains a semicolon. Do you probably want to execute another command after %s? I can't do that for you..." % (key.upper(),key.upper())
    sys.exit(39)
  return LINE

# ======================================================================= #
def get_pairs(QMinlines,i):
  nacpairs=[]
  while True:
    i+=1
    try:
      line=QMinlines[i].lower()
    except IndexError:
      print '"keyword select" has to be completed with an "end" on another line!'
      sys.exit(40)
    if 'end' in line:
      break
    fields=line.split()
    try:
      nacpairs.append([int(fields[0]),int(fields[1])])
    except ValueError:
      print '"nacdr select" is followed by pairs of state indices, each pair on a new line!'
      sys.exit(41)
  return nacpairs,i

# =============================================================================================== #
# =============================================================================================== #
# =========================================== readQMin and gettasks ============================= #
# =============================================================================================== #
# =============================================================================================== #


# ======================================================================= #     OK
def readQMin(QMinfilename):
  '''Reads the time-step dependent information from QMinfilename. This file contains all information from the current SHARC job: geometry, velocity, number of states, requested quantities along with additional information. The routine also checks this input and obtains a number of environment variables necessary to run COLUMBUS.

  Reads also the information from SH2COL

  Steps are:
  - open and read QMinfilename
  - Obtain natom, comment, geometry (, velocity)
  - parse remaining keywords from QMinfile
  - check keywords for consistency, calculate nstates, nmstates
  - obtain environment variables for path to COLUMBUS and scratch directory, and for error handling

  Arguments:
  1 string: name of the QMin file

  Returns:
  1 dictionary: QMin'''

  # read QMinfile
  try:
    QMinfile=open(QMinfilename,'r')
  except IOError:
    print 'QM input file "%s" not found!' % (QMinfilename)
    sys.exit(42)
  QMinlines=QMinfile.readlines()
  QMinfile.close()
  QMin={}



  # Get natom
  try:
    natom=int(QMinlines[0].split()[0])
  except ValueError:
    print 'first line must contain the number of atoms!'
    sys.exit(43)
  QMin['natom']=natom
  if len(QMinlines)<natom+4:
    print 'Input file must contain at least:\nnatom\ncomment\ngeometry\nkeyword "states"\nat least one task'
    sys.exit(44)



  # Save Comment line
  QMin['comment']=QMinlines[1]



  # Get geometry and possibly velocity
  QMin['geo']=[]
  QMin['veloc']=[]
  hasveloc=True
  for i in range(2,natom+2):
    if not containsstring('[a-zA-Z][a-zA-Z]?[0-9]*.*[-]?[0-9]+[.][0-9]*.*[-]?[0-9]+[.][0-9]*.*[-]?[0-9]+[.][0-9]*', QMinlines[i]):
      print 'Input file does not comply to xyz file format! Maybe natom is just wrong.'
      sys.exit(45)
    fields=QMinlines[i].split()
    for j in range(1,4):
      fields[j]=float(fields[j])
    QMin['geo'].append(fields[0:4])
    if len(fields)>=7:
      for j in range(4,7):
        fields[j]=float(fields[j])
      QMin['veloc'].append(fields[4:7])
    else:
      hasveloc=False
  if not hasveloc:
    QMin=removekey(QMin,'veloc')



  # Parse remaining file
  i=natom+1
  while i+1<len(QMinlines):
    i+=1
    line=QMinlines[i]
    line=re.sub('#.*$','',line)
    if len(line.split())==0:
      continue
    key=line.lower().split()[0]
    if 'savedir' in key:
      args=line.split()[1:]
    else:
      args=line.lower().split()[1:]
    if key in QMin:
      print 'Repeated keyword %s in line %i in input file! Check your input!' % (linekeyword(line),i+1)
      continue  # only first instance of key in QM.in takes effect
    if len(args)>=1 and 'select' in args[0]:
      pairs,i=get_pairs(QMinlines,i)
      QMin[key]=pairs
    else:
      QMin[key]=args



  # Calculate states, nstates, nmstates
  for i in range(len(QMin['states'])):
    QMin['states'][i]=int(QMin['states'][i])

  nstates=0
  nmstates=0
  for i in range(len(QMin['states'])):
    nstates+=QMin['states'][i]
    nmstates+=QMin['states'][i]*(i+1)
  QMin['nstates']=nstates
  QMin['nmstates']=nmstates



  # Various logical checks
  if not 'states' in QMin:
    print 'Number of states not given in QM input file %s!' % (QMinfilename)
    sys.exit(46)

  possibletasks=['h','soc','dm','grad','nacdr','nacdt','overlap','angular','phases']
  if not any([i in QMin for i in possibletasks]):
    print 'No tasks found! Tasks are "h", "soc", "dm", "grad", "nacdt", "nacdr" and "overlap".'
    sys.exit(47)

  if ('samestep' in QMin and 'init' in QMin) or ('restart' in QMin and 'init' in QMin):
    print '"Init" and "Samestep" cannot be both present in QM.in!'
    sys.exit(48)

  if 'phases' in QMin:
    QMin['overlap']=[]

  if 'overlap' in QMin and 'init' in QMin:
    print '"overlap" and "phases" cannot be calculated in the first timestep! Delete either "overlap" or "init"'
    sys.exit(49)

  if not 'init' in QMin and not 'samestep' in QMin and not 'restart' in QMin:
    QMin['newstep']=[]

  if not any([i in QMin for i in ['h','soc','dm','grad','nacdt','nacdr']]) and ('overlap' in QMin or 'ion' in QMin):
    QMin['h']=[]

  if len(QMin['states'])>8:
    print 'Higher multiplicities than octets are not supported!'
    sys.exit(50)

  if 'h' in QMin and 'soc' in QMin:
    QMin=removekey(QMin,'h')

  if 'nacdt' in QMin:
    print 'Numerical non-adiabatic couplings not available! Use "nacdr" instead!'
    sys.exit(51)

  if 'dmdr' in QMin:
    print 'Dipole moment gradients not available!'
    sys.exit(52)

  if 'socdr' in QMin:
    print 'Spin-orbit coupling gradients not available!'
    sys.exit(53)

  if not 'step' in QMin:
    QMin['step']=['0']


  # Process the gradient requests
  if 'grad' in QMin:
    if len(QMin['grad'])==0 or QMin['grad'][0]=='all':
      QMin['grad']=[ i+1 for i in range(nmstates)]
    else:
      for i in range(len(QMin['grad'])):
        try:
          QMin['grad'][i]=int(QMin['grad'][i])
        except ValueError:
          print 'Arguments to keyword "grad" must be "all" or a list of integers!'
          sys.exit(54)
        if QMin['grad'][i]>nmstates:
          print 'State for requested gradient does not correspond to any state in QM input file state list!'
          sys.exit(55)

  # Process the non-adiabatic coupling requests
  # type conversion has already been done
  if 'nacdr' in QMin:
    if len(QMin['nacdr'])>=1:
      nacpairs=QMin['nacdr']
      for i in range(len(nacpairs)):
        if nacpairs[i][0]>nmstates or nacpairs[i][1]>nmstates:
          print 'State for requested non-adiabatic couplings does not correspond to any state in QM input file state list!'
          sys.exit(56)
    else:
      QMin['nacdr']=[ [j+1,i+1] for i in range(nmstates) for j in range(i)]

  # Process the overlap requests
  # identically to the nac requests
  #if 'overlap' in QMin:
    #if len(QMin['overlap'])>=1:
      #nacpairs=QMin['overlap']
      #for i in range(len(nacpairs)):
        #if nacpairs[i][0]>nmstates or nacpairs[i][1]>nmstates:
          #print 'State for requested overlap does not correspond to any state in QM input file state list!'
          #sys.exit(57)
    #else:
      #QMin['overlap']=[ [j+1,i+1] for i in range(nmstates) for j in range(i)]

  # Process the ionization requests
  #if 'ion' in QMin:
    #if len(QMin['ion'])>=1:
      #nacpairs=QMin['ion']
      #for i in range(len(nacpairs)):
        #if nacpairs[i][0]>nmstates or nacpairs[i][1]>nmstates:
          #print 'State for requested Dyson norm does not correspond to any state in QM input file state list!'
          #sys.exit(58)
    #else:
      #QMin['ion']=[ [j+1,i+1] for i in range(nmstates) for j in range(i)]
    #if len(QMin['ion'])==0 or QMin['ion'][0]=='all':
      #QMin['ion']=[i+1 for i in range(nmstates)]
    #else:
      #for i in range(len(QMin['ion'])):
        #try:
          #QMin['ion'][i]=int(QMin['ion'][i])
        #except ValueError:
          #print 'Arguments to keyword "ion" must be "all" or a list of integers!'
          #sys.exit(59)
        #if QMin['ion'][i]>nmstates:
          #print 'State for requested ionization yield does not correspond to any state in QM input file state list!'
          #sys.exit(60)



  # open COLUMBUS.resources
  if os.path.isfile('COLUMBUS.resources'):
    sh2colf=open('COLUMBUS.resources','r')
  else:
    sh2colf=open('SH2COL.inp','r')
  sh2col=sh2colf.readlines()
  sh2colf.close()



  # Set up environment variables: $COLUMBUS, $MOLCAS, TEMPLATE, Dyson
  QMin['pwd']=os.getcwd()

  QMin['columbus']=get_sh2col_environ(sh2col,'columbus')
  os.environ['COLUMBUS']=QMin['columbus']

  if 'ion' in QMin or 'overlap' in QMin:
    QMin['molcas']=get_sh2col_environ(sh2col,'molcas')
    os.environ['MOLCAS']=QMin['molcas']

  QMin['template']=get_sh2col_environ(sh2col,'template',False)

  if 'ion' in QMin or 'overlap' in QMin:
    QMin['wfoverlap']=get_sh2col_environ(sh2col,'wfoverlap',False,False)
    if QMin['wfoverlap']==None:
      ciopath=os.path.join(os.path.expandvars(os.path.expanduser('$SHARC')),'wfoverlap.x')
      if os.path.isfile(ciopath):
        QMin['wfoverlap']=ciopath
      else:
        print 'Give path to wfoverlap.x in COLUMBUS.resources!'
        sys.exit(61)
  #if 'overlap' in QMin:
    #QMin['cioverlap']=get_sh2col_environ(sh2col,'cioverlap',False,False)
    #if QMin['cioverlap']==None:
      #ciopath=os.path.join(os.path.expandvars(os.path.expanduser('$SHARC')),'cioverlap.x')
      #if os.path.isfile(ciopath):
        #QMin['cioverlap']=ciopath
      #else:
        #print 'Give path to cioverlap.x in COLUMBUS.resources!'
        #sys.exit(62)

  # Set up scratchdir
  QMin['scratchdir']=get_sh2col_environ(sh2col,'scratchdir',False,False)
  if QMin['scratchdir']==None:
    QMin['scratchdir']=QMin['pwd']+'/SCRATCHDIR/'
  checkscratch(QMin['scratchdir'])

  ## Set up savedir
  #if not 'savedir' in QMin:
    ## savedir may be read from QM.in file
    #QMin['savedir']=get_sh2col_environ(sh2col,'savedir',False,False)
    #if QMin['savedir']==None:
      #QMin['savedir']=QMin['pwd']+'/SAVEDIR/'
  #else:
    #QMin['savedir']=QMin['savedir'][0]
  #checkscratch(QMin['savedir'])

  # Set up savedir
  if 'savedir' in QMin:
    # savedir may be read from QM.in file
    line=QMin['savedir'][0]
  else:
    line=get_sh2col_environ(sh2col,'savedir',False,False)
    if line==None:
      line=QMin['pwd']+'/SAVEDIR/'
  line=os.path.expandvars(line)
  line=os.path.expanduser(line)
  line=os.path.abspath(line)
  if 'init' in QMin:
    checkscratch(line)
  QMin['savedir']=line
  link(QMin['savedir'],os.path.join(QMin['pwd'],'SAVE'),False,False)

  line=getsh2colkey(sh2col,'debug')
  if line[0]:
    #if line[1].lower().strip()=='true':
    global DEBUG
    DEBUG=True

  line=getsh2colkey(sh2col,'no_print')
  if line[0]:
    #if line[1].lower().strip()=='true':
    global PRINT
    PRINT=True


  # Set default memory for runc and default screening threshold for cioverlaps
  QMin['colmem']=10
  line=getsh2colkey(sh2col,'memory')
  if line[0]:
    try:
      QMin['colmem']=int(line[1])
    except ValueError:
      print 'COLUMBUS memory does not evaluate to numerical value!'
      sys.exit(63)
  else:
    print 'WARNING: Please set memory for COLUMBUS in COLUMBUS.resources (in MB)! Using 10 MB default value!'

  line=getsh2colkey(sh2col,'integrals')
  if line[0]:
    arg=line[1].strip()
    allowed=['dalton','seward']
    if not arg in allowed:
      print 'Do not know integral program "%s".' % (arg)
      sys.exit(64)
    QMin['integrals']=arg
  else:
    QMin['integrals']='seward'

  if QMin['integrals']!='seward' and 'soc' in QMin:
    print 'Cannot calculate spin-orbit couplings with integral program "%s"!' % (QMin['integrals'])
    sys.exit(65)
  #if QMin['integrals']!='seward' and 'ion' in QMin:
    #print 'Cannot calculate Dyson norms with integral program "%s"!' % (QMin['integrals'])
    #sys.exit(66)

  line=getsh2colkey(sh2col,'molcas_rasscf')
  if line[0]:
    QMin['molcas_rasscf']=True
  else:
    QMin['molcas_rasscf']=False
  if QMin['molcas_rasscf']:
    if 'grad' in QMin or 'nacdr' in QMin:
      print 'Cannot calculate gradients/NAC with orbitals from MOLCAS rasscf!'
      sys.exit(67)
    if QMin['integrals']!='seward':
      print 'Orbitals from MOLCAS rasscf only works with seward integrals!'
      sys.exit(68)
    if 'molden' in QMin:
      print 'Cannot make molden file with orbitals from MOLCAS rasscf!'
      sys.exit(69)

  line=getsh2colkey(sh2col,'nooverlap')
  if line[0]:
    QMin['nooverlap']=[]
  if 'nooverlap' in QMin and 'overlap' in QMin:
    print 'keyword "overlap" in QM.in, but keyword "nooverlap" in COLUMBUS.resources!'
    sys.exit(70)
  if 'nooverlap' in QMin and 'ion' in QMin:
    print 'keyword "ion" in QM.in, ignoring keyword "nooverlap"!'
    del QMin['nooverlap']


  if not 'nooverlap' in QMin:
    QMin['wfthres']=0.97
    line=getsh2colkey(sh2col,'wfthres')
    if line[0]:
      try:
        QMin['wfthres']=float(line[1])
      except ValueError:
        print 'WFoverlaps threshold variable does not evaluate to numerical value!'
        sys.exit(71)
    else:
      print 'WARNING: Please set wfthres to some appropriate value (floating point number)! Using 0.97 default value!'

  #if 'ion' in QMin:
    ## dyson.x in principal does not need to screen anymore, since this is done by read_civfl
    #QMin['dysonthres']=1e-12
    #line=getsh2colkey(sh2col,'dysonthres')
    #if line[0]:
      #try:
        #QMin['dysonthres']=float(line[1])
      #except ValueError:
        #print 'Dyson threshold variable does not evaluate to numerical value!'
        #sys.exit(72)

  line=getsh2colkey(sh2col,'always_guess')
  if line[0]:
    QMin['always_guess']=[]
  else:
    line=getsh2colkey(sh2col,'always_scf')
    if line[0]:
      QMin['always_orb_init']=[]

  line=getsh2colkey(sh2col,'always_orb_init')
  if line[0]:
    QMin['always_orb_init']=[]
  else:
    line=getsh2colkey(sh2col,'always_mocoef_init')
    if line[0]:
      QMin['always_orb_init']=[]

  QMin['ncpu']=1
  line=getsh2colkey(sh2col,'ncpu')
  if line[0]:
    try:
      QMin['ncpu']=int(line[1])
    except ValueError:
      print 'Number of CPUs does not evaluate to numerical value!'
      sys.exit(73)

  # path to modified runc
  line=getsh2colkey(sh2col,'runc')
  if line[0]:
    RUNC=line[1]
    RUNC=os.path.expandvars(RUNC)
    RUNC=os.path.expanduser(RUNC)
    RUNC=removequotes(RUNC).strip()
    if containsstring(';',RUNC):
      print "runc path contains a semicolon. Do you probably want to execute another command after runc? I can't do that for you..."
      sys.exit(74)
    QMin['runc']=RUNC
  else:
    QMin['runc']=os.path.join(QMin['columbus'],'runc')


  # read ncore and ndocc from SH2COL
  QMin['ncore']=-1
  line=getsh2colkey(sh2col,'numfrozcore')
  if line[0]:
    try:
      QMin['ncore']=max(0,int(line[1]))
    except ValueError:
      print 'numfrozcore does not evaluate to numerical value!'
      sys.exit(75)
  line=getsh2colkey(sh2col,'numocc')
  if line[0]:
    try:
      QMin['ndocc']=int(line[1])
    except ValueError:
      print 'numocc does not evaluate to numerical value!'
      sys.exit(76)


  # get the necessary template mappings
  # Map: state -> mult,state,ms                         statemap
  # Map: mult -> dir, DRT and dir, DRT -> mult          multmap
  # Map: dir -> mocoefdir                               mocoefmap
  # Map: dir -> socinr/isc                              socimap
  # List: dir   (dependency-resolved order)             joblist

  # obtain the statemap
  statemap={}
  i=1
  for imult,istate,ims in itnmstates(QMin['states']):
    statemap[i]=[imult,istate,ims]
    i+=1
  QMin['statemap']=statemap

  # get the multmap
  multmap={}
  for mult in itmult(QMin['states']):
    # find the appropriate line in COLUMBUS.resources
    i=-1
    while True:
      i+=1
      try:
        line=re.sub('#.*$','',sh2col[i])
      except IndexError:
        print 'Multiplicity %i has no template directory given in COLUMBUS.resources!' % (mult)
        sys.exit(77)
      line=line.split()
      if len(line)==0:
        continue
      if line[0].lower()=='dir':
        if int(line[1])==mult:
          break
    # lines look like ['dir', '1', '1_3/']
    # put into multmap
    if line[2][-1]!='/':
      line[2]+='/'
    multmap[mult]=line[2]
    if not line[2] in multmap:
      multmap[line[2]]=[mult]
    else:
      multmap[line[2]].append(mult)


  # get the mocoefmap
  mocoefmap={}
  # first get all jobs
  for mult in itmult(QMin['states']):
    if not multmap[mult] in mocoefmap:
      mocoefmap[multmap[mult]]=None
  # now look up for all jobs the mocoefdir
  for job in mocoefmap:
    # find the line in COLUMBUS.resources
    i=-1
    while True:
      i+=1
      try:
        line=re.sub('#.*$','',sh2col[i])
      except IndexError:
        print 'WARNING: no mocoef directory specified for %s, will use its own mocoefs' % (job)
        line=['mocoef',job,job]
        break
      line=line.split()
      if len(line)==0:
        continue
      if line[0].lower()=='mocoef':
        if line[1][-1]!='/':
          line[1]+='/'
        if line[1] not in mocoefmap:
          continue
        if line[1]==job:
          break
    # line looks like ['mocoef', '1_3/', '1_3/']
    # put into mocoefmap
    if line[2][-1]!='/':
      line[2]+='/'
    mocoefmap[job]=line[2]
  # do a topological sort of the mocoefmap
  joblist=toposort(mocoefmap)



  # now put the DRTs inside the multmap and create the socimap
  socimap={}
  for mult in itmult(QMin['states']):
    job=multmap[mult]
    socimode,drt=checktemplate(QMin['template']+'/'+job,mult,QMin['states'],QMin['integrals'])
    multmap[mult]=[job,drt]
    multmap[job][drt-1]=mult
    socimap[job]=socimode

  # make the multpairlist
  # { state: set([states]) }
  if 'ion' in QMin:
    multset=set()
    for istate in QMin['statemap']:
      m1,s1,ms1=tuple(QMin['statemap'][istate])
      multset.add(m1)
    multpairset=set()
    for i in multset:
      for j in multset:
        if j-i==1:
          multpairset.add( (i,j) )
    QMin['multpairlist']=list(multpairset)

    #ionmap={}
    #for istate in QMin['ion']:
      #m1,s1,ms1=tuple(QMin['statemap'][istate])
      #ionmap[istate]=set()
      #for jstate in range(1,nmstates+1):
        #m2,s2,ms2=tuple(QMin['statemap'][jstate])
        #if abs(m1-m2)==1 and abs(ms1-ms2)==0.5:
          #ionmap[istate].add(jstate)
    ## the ionlist contains only the pairs of states which need to be calculated
    #ionlist=[]
    #for istate in ionmap:
      #for jstate in ionmap[istate]:
        #m1,s1,ms1=tuple(QMin['statemap'][istate])
        #m2,s2,ms2=tuple(QMin['statemap'][jstate])
        ## check if already in ionlist
        #isthere=False
        #for i in ionlist:
          #if (m1,s1,m2,s2)==tuple(i) or (m2,s2,m1,s1)==tuple(i):
            #isthere=True
            #break
        #if not isthere:
          #ionlist.append([m1,s1,m2,s2])
    ##QMin['ionmap']=ionmap
    #QMin['ionlist']=ionlist


  # check the savedir
  if not 'init' in QMin:
    required=[]
    for job in joblist:
      jobdir=job.replace('/','_')
      if 'overlap' in QMin:
        if 'samestep' in QMin or 'restart' in QMin:
          for drt,mult in enumerate(multmap[job]):
            required.append('dets.%i.old' % (mult))
          if QMin['integrals']=='seward':
            required.append('molcas.input.seward.%s.old' % (jobdir))
          elif QMin['integrals']=='dalton':
            required.append('daltaoin.%s.old' % (jobdir))
          if QMin['molcas_rasscf']:
            required.append('molcas.RasOrb.%s.old' % (jobdir))
          else:
            required.append('mocoef.%s.old' % (jobdir))
        else:
          for drt,mult in enumerate(multmap[job]):
            required.append('dets.%i' % (mult))
          if QMin['integrals']=='seward':
            required.append('molcas.input.seward.%s' % (jobdir))
          elif QMin['integrals']=='dalton':
            required.append('daltaoin.%s' % (jobdir))
          if QMin['molcas_rasscf']:
            required.append('molcas.RasOrb.%s' % (jobdir))
          else:
            required.append('mocoef.%s' % (jobdir))
      if job==mocoefmap[job]:
        if 'restart' in QMin:
          if QMin['molcas_rasscf']:
            required.append('molcas.RasOrb.%s.old' % (jobdir))
          else:
            required.append('mocoef.%s.old' % (jobdir))
        else:
          if QMin['molcas_rasscf']:
            required.append('molcas.RasOrb.%s' % (jobdir))
          else:
            required.append('mocoef.%s' % (jobdir))
    ls=os.listdir(QMin['savedir'])
    for f in required:
      if not f in ls:
        print 'File %s missing in savedir %s!' % (f,QMin['savedir'])
        sys.exit(78)


  QMin['multmap']=multmap
  QMin['mocoefmap']=mocoefmap
  QMin['joblist']=joblist
  QMin['socimap']=socimap



  if 'backup' in QMin:
    backupdir=QMin['savedir']+'/backup'
    backupdir1=backupdir
    i=0
    while os.path.isdir(backupdir1):
      i+=1
      if 'step' in QMin:
        backupdir1=backupdir+'/step%s_%i' % (QMin['step'][0],i)
      else:
        backupdir1=backupdir+'/calc_%i' % (i)
    QMin['backup']=backupdir1


  print '\n\n\n'
  return QMin

## ======================================================================= #

def toposort(jobmap):
  '''does a topological sort of the list jobmap according to the Kahn (1962) doi:10.1145/368996.369025 algorithm'''

  #sortedjobmap=[]
  ## first, make a useful representation of the graph
  ## gather the set of nodes without incoming edges in edgefree
  #graph={}
  #edgefree=set([])
  #for j in jobmap:
    #if j==jobmap[j]:
      #edgefree.add(j)
    #else:
      #graph[j]=jobmap[j]

  ## do the loop
  #while len(edgefree)!=0:
    ## take a random element from edgefree
    #n=edgefree.pop()
    ## put it into the sortedlist
    #sortedjobmap.append(n)
    ## remove nodes pointing to n
    #toremove=[]
    #for m in graph:
      #if graph[m]==n:
        #toremove.append(m)
    #for m in toremove:
      #del graph[m]
      #edgefree.add(m)
  ## if there is leftover in graph, then there was a cyclic depencency
  #if len(graph)>0:
    #print 'Cyclic dependency in the jobmap!'
    #print jobmap
    #sys.exit(79)

  # actually, the situation is easier: there must only be self-loops and paths of length 1
  sortedjobmap=[]
  edgefree=set([])
  for j in jobmap:
    if jobmap[j]==j:
      edgefree.add(j)
  for j in edgefree:
    sortedjobmap.append(j)
  for j in jobmap:
    if not j in edgefree:
      if jobmap[j] in edgefree:
        sortedjobmap.append(j)
      else:
        print 'Dependecy violation in mocoefdir specifications!'
        sys.exit(80)

  return sortedjobmap

## ======================================================================= #








































# ======================================================================= #
def gettasks(QMin):
  ''''''

  states=QMin['states']
  nstates=QMin['nstates']
  nmstates=QMin['nmstates']

  # Currently implemented keywords: soc, dm, grad, nac, samestep, init
  tasks=[]
  # During initialization, create all temporary directories
  # and link them appropriately
  tasks.append(['mkdir',QMin['scratchdir']])
  tasks.append(['link', QMin['scratchdir'],QMin['pwd']+'/SCRATCH',False])
  tasks.append(['mkdir',QMin['scratchdir']+'/JOB'])
  tasks.append(['mkdir',QMin['scratchdir']+'/OVERLAP'])
  tasks.append(['mkdir',QMin['scratchdir']+'/DYSON'])
  tasks.append(['mkdir',QMin['scratchdir']+'/MOLCAS'])
  if not os.path.isdir(QMin['scratchdir']+'/KEEP'):
    tasks.append(['mkdir',QMin['scratchdir']+'/KEEP'])

  if 'init' in QMin:
    tasks.append(['mkdir',QMin['savedir']])
    tasks.append(['link', QMin['savedir'],QMin['pwd']+'/SAVE',False])

  if 'molden' in QMin and not os.path.isdir(QMin['savedir']+'/MOLDEN'):
    tasks.append(['mkdir',QMin['savedir']+'/MOLDEN'])

  if not 'samestep' in QMin and not 'init' in QMin and not 'restart' in QMin:
    tasks.append(['movetoold'])

  if 'backup' in QMin:
    tasks.append(['mkdir',QMin['savedir']+'/backup/'])
    tasks.append(['mkdir',QMin['backup']])

  # do all COLUMBUS calculations
  if not 'overlaponly' in QMin:

    for job in QMin['joblist']:

      tasks.append(['cleanup',QMin['scratchdir']+'/JOB'])
      tasks.append(['getinput',QMin['template']+'/'+job])
      tasks.append(['getmcinput',QMin['template']+'/'+QMin['mocoefmap'][job]])
      tasks.append(['checkinput'])
      tasks.append(['writegeom'])

      # control.run ==========================================================
      string='niter=1\n'
      if 'dm' in QMin or 'grad' in QMin or 'nacdr' in QMin:
        string+='ciudgmom\nciprop\n'
      if QMin['molcas_rasscf']:
        string+='rasscf\n'
      else:
        string+='mcscf\n'
      if QMin['integrals']=='seward':
        string+='seward\n'
      elif QMin['integrals']=='dalton':
        string+='hermit\n'
      # check in socimap whether isc keyword or socinr keyword is necessary
      if QMin['socimap'][job]==1:
        # socinr
        string+='ciudgav\nsocinr='
        for mult in range(len(states)):
          if states[mult]<=0:
            string+='0'
          elif QMin['multmap'][mult+1][0]==job:
            string+='%i' % (states[mult])
          else:
            string+='0'
          if mult+1<len(states):
            string+=':'
        string+='\n'
      elif QMin['socimap'][job]==0:
        # isc
        string+='ciudgav\nisc\n'
      elif QMin['socimap'][job]==-1:
        # single drt
        string+='ciudg\n'
        # in this case, the ciudgin file has to be adjusted for the number of states
      if 'grad' in QMin or 'nacdr' in QMin:
        string+='nadcoupl\n'
      if 'always_guess' in QMin:
        string+='scf\n'
      elif 'init' in QMin and QMin['mocoefmap'][job]==job:
        if QMin['molcas_rasscf']:
          filestring='molcas.RasOrb.init'
        else:
          filestring='mocoef_mc.init'
        tryfile=os.path.join(QMin['pwd'],filestring)
        if not os.path.isfile(tryfile) and not os.path.isfile(tryfile+'.'+job.replace('/','')):
          string+='scf\n'
      if 'molden' in QMin:
        string+='molden\n'
      tasks.append(['write',QMin['scratchdir']+'/JOB/control.run',string])
      if DEBUG:
        print string

      # transmomin ==========================================================
      transmomstring='CI\n'
      # use the list of gradients from QM.in or produce one
      if 'grad' in QMin:
        #if QMin['grad'][0]=='all':
          #gargs=[i+1 for i in range(nmstates)]
        #else:
        gargs=QMin['grad']
      else:
        gargs=[]
      # use the list of gradients from QM.in or produce one
      if 'nacdr' in QMin:
        #if len(QMin['nacdr'])>=2:
        nargs=QMin['nacdr']
        #else:
          #nargs=[ [i+1,j+1] for i in range(nmstates) for j in range(nmstates)]
      else:
        nargs=[]
      # go over all multiplicities of the current job
      for drt,mult in enumerate(QMin['multmap'][job]):
        for j1 in range(1,1+states[mult-1]):
          for j2 in range(j1,1+states[mult-1]):
            string='%i %i %i %i ' % (drt+1,j1,drt+1,j2)
            # figure out whether the gradient/nac should be calculated
            if j1==j2:
              # gradient
              requested=False
              for iarg in gargs:
                # go through all requests
                imult,istate,ims=tuple(QMin['statemap'][iarg])
                if imult==mult and istate==j1:
                  requested=True
              if not requested:
                string+='T'
            else:
              # nac
              requested=False
              for iarg in nargs:
                # go through all requests
                imult1,istate1,ims1=tuple(QMin['statemap'][iarg[0]])
                imult2,istate2,ims2=tuple(QMin['statemap'][iarg[1]])
                if mult==imult1==imult2 and ims1==ims2 and ( (j1,j2)==(istate1,istate2) or (j1,j2)==(istate2,istate1) ):
                #if mult==imult1==imult2 and ( (j1==istate1 and j2==istate2) or (j1==istate2 and j2==istate1) ) and ims1==ims2:
                  requested=True
              if not requested:
                string+='T'
            transmomstring+=string+'\n'
      tasks.append(['write',QMin['scratchdir']+'/JOB/transmomin',transmomstring])
      if DEBUG:
        print transmomstring

      # ciudgin ==========================================================
      # if socinr is not used, the ciudgin file has to be adjusted for the correct number of states
      if not QMin['socimap'][job]==1:
        for drt,mult in enumerate(QMin['multmap'][job]):
          numberstates=states[mult-1]
          if QMin['socimap'][job]==0:
            filename='ciudgin.drt%i' % (drt+1)
          elif QMin['socimap'][job]==-1:
            filename='ciudgin'
          ciudgin=readfile(QMin['template']+'/'+job+'/'+filename)
          #f=open(QMin['template']+'/'+job+'/'+filename)
          #ciudgin=f.readlines()
          #f.close()
          for i in range(len(ciudgin)):
            if 'nroot' in ciudgin[i].lower():
              ciudgin[i]=' NROOT = %i\n' % (numberstates)
            if 'rtolci' in ciudgin[i].lower():
              tol=ciudgin[i].split('=')[1].split(',')
              tol2=[]
              for j in tol:
                try:
                  tol2.append(float(j))
                except ValueError:
                  break
              ciudgin[i]=' RTOLCI ='
              for j in range(numberstates):
                if j<len(tol2):
                  ciudgin[i]+='%9.6f,' % (tol2[j])
                else:
                  ciudgin[i]+='%9.6f,' % (tol2[-1])
              ciudgin[i]+='\n'
          string=' '.join(ciudgin)
          tasks.append(['write',QMin['scratchdir']+'/JOB/'+filename,string])

      # cidrtin ==========================================================
      # if socinr is used, adjust the maximum multiplicity in cidrtin to the one actually needed
      if QMin['socimap'][job]==1:
        max_needed=max(QMin['multmap'][job])
        f=open(QMin['template']+'/'+job+'/cidrtin')
        cidrtin=f.readlines()
        f.close()
        for i in range(len(cidrtin)):
          if 'maximal spin multiplicity' in cidrtin[i]:
            cidrtin[i]='%i / maximal spin multiplicity\n' % (max_needed)
        string=' '.join(cidrtin)
        tasks.append(['write',QMin['scratchdir']+'/JOB/cidrtin',string])

      # prepare the mocoef guess ============================================
      if not 'always_guess' in QMin:
        if 'always_orb_init' in QMin:
          if QMin['molcas_rasscf']:
            filestring='molcas.RasOrb.init'
          else:
            filestring='mocoef_mc.init'
          tryfile=os.path.join(QMin['pwd'],filestring+'.'+job.replace('/',''))
          if os.path.isfile(tryfile):
            tasks.append(['getmo',tryfile])
          else:
            tryfile=os.path.join(QMin['pwd'],filestring)
            tasks.append(['getmo',tryfile])
        else:
          if QMin['mocoefmap'][job]==job:
            if 'init' in QMin:
              if QMin['molcas_rasscf']:
                filestring='molcas.RasOrb.init'
              else:
                filestring='mocoef_mc.init'
              tryfile=os.path.join(QMin['pwd'],filestring+'.'+job.replace('/',''))
              if os.path.isfile(tryfile):
                tasks.append(['getmo',tryfile])
              else:
                tryfile=os.path.join(QMin['pwd'],filestring)
                if os.path.isfile(tryfile):
                  tasks.append(['getmo',tryfile])
                # if this file is not there, an scf calculation is performed (see above)
            else:
              if QMin['molcas_rasscf']:
                filestring='molcas.RasOrb'
              else:
                filestring='mocoef'
              mocoefjob=QMin['mocoefmap'][job].replace('/','_')
              if 'samestep' in QMin:
                tryfile=os.path.join(QMin['savedir'],filestring+'.'+mocoefjob)
              else:
                tryfile=os.path.join(QMin['savedir'],filestring+'.'+mocoefjob+'.old')
              tasks.append(['getmo',tryfile])
          else:
            mocoefjob=QMin['mocoefmap'][job].replace('/','_')
            if QMin['molcas_rasscf']:
              filestring='molcas.RasOrb'
            else:
              filestring='mocoef'
            tryfile=os.path.join(QMin['savedir'],filestring+'.'+mocoefjob)
            tasks.append(['getmo',tryfile])

      # run COLUMBUS ========================================================
      tasks.append(['runc'])

      # get output ==========================================================
      tasks.append(['get_COLout',job])

      # do the civecconsolidate runs, keep data =============================
      # make_dets before save_data, because save_data moves the eivector files
      #if 'ion' in QMin:
        #tasks.append(['make_dets',job])
      if not 'nooverlap' in QMin or 'ion' in QMin:
        tasks.append(['make_dets_new',job])

      # save mocoef, eigenvectors, molcas,input =============================
      if not 'samestep' in QMin:
        tasks.append(['save_data',job])

      # keep slaterfiles ====================================================
      tasks.append(['keep_data',job])

      # keep molden file ====================================================
      if 'molden' in QMin:
        tasks.append(['copymolden',job])

      # remove temporary files ==============================================
      #if not DEBUG:
        #tasks.append(['cleanup',QMin['scratchdir']+'/JOB/WORK'])

  # End of loop over COLUMBUS jobs
  # now comes the dyson and cioverlaps parts

  # Dyson part
  #if 'ion' in QMin:
    #for pair in QMin['ionlist']:
      ## pair is [m1,s1,m2,s2]
      #tasks.append(['dyson',pair])


  ## Overlap part
  if 'overlap' in QMin:
    tasks.append(['cleanup',QMin['scratchdir']+'/OVERLAP'])

    # prepare the overlap integrals
    # IMPORTANT: we make the restriction that all jobs have the same AO basis set
    jobdir=QMin['joblist'][0].replace('/','_')
    if QMin['integrals']=='seward':
      tasks.append(['mkdir',QMin['scratchdir']+'/OVERLAP/molcas'])
      tasks.append(
        ['writemolcas',
        QMin['savedir']+'/molcas.input.seward.%s.old' % (jobdir),
        QMin['savedir']+'/molcas.input.seward.%s' % (jobdir),
        QMin['scratchdir']+'/OVERLAP/molcas.input']
      )
      tasks.append(['runmolcas'])
      tasks.append(['link',QMin['scratchdir']+'/OVERLAP/molcas/molcas.RunFile',QMin['scratchdir']+'/OVERLAP/RUNFILE'])
      tasks.append(['link',QMin['scratchdir']+'/OVERLAP/molcas/molcas.OneInt', QMin['scratchdir']+'/OVERLAP/ONEINT'])
      tasks.append(['link',QMin['scratchdir']+'/OVERLAP/molcas/molcas.env', QMin['scratchdir']+'/OVERLAP/molcas.env'])
    elif QMin['integrals']=='dalton':
      tasks.append(
        ['writedalton',
        QMin['savedir']+'/daltaoin.%s.old' % (jobdir),
        QMin['savedir']+'/daltaoin.%s' % (jobdir),
        QMin['scratchdir']+'/OVERLAP/daltcomm',
        QMin['scratchdir']+'/OVERLAP/daltaoin',
        QMin['template']+'/'+QMin['joblist'][0]]
      )
      tasks.append(['rundalton'])

    for job in QMin['joblist']:
      jobdir=job.replace('/','_')
      if QMin['molcas_rasscf']:
        tasks.append(['mo_convert',QMin['savedir']+'/molcas.RasOrb.%s.old' % ( jobdir),QMin['savedir']+'/mocoef.%s.old' % (jobdir),QMin['savedir']])
        tasks.append(['mo_convert',QMin['savedir']+'/molcas.RasOrb.%s' % ( jobdir),QMin['savedir']+'/mocoef.%s' % (jobdir),QMin['savedir']])
      tasks.append(['link',QMin['savedir']+'/mocoef.%s.old' % (jobdir),QMin['scratchdir']+'/OVERLAP/mocoef1'])
      tasks.append(['link',QMin['savedir']+'/mocoef.%s'     % (jobdir),QMin['scratchdir']+'/OVERLAP/mocoef2'])

      for drt,mult in enumerate(QMin['multmap'][job]):
        tasks.append(['link',QMin['savedir']+'/dets.%i.old' % (mult), QMin['scratchdir']+'/OVERLAP/dets1'])
        tasks.append(['link',QMin['savedir']+'/dets.%i' % (mult),     QMin['scratchdir']+'/OVERLAP/dets2'])

        tasks.append(['write_ciovin',mult])
        tasks.append(['runcioverlap',mult])
        tasks.append(['get_CIOoutput',job,drt+1])




  ## Dyson part
  if 'ion' in QMin:
    tasks.append(['cleanup',QMin['scratchdir']+'/DYSON'])

    # IMPORTANT: we make the restriction that all jobs have the same AO basis set
    # we always use the AO overlap matrix of the first job in the joblist
    jobdir=QMin['joblist'][0].replace('/','_')
    if QMin['integrals']=='seward':
      tasks.append(['link',QMin['scratchdir']+'/KEEP/RUNFILE.%s' % (jobdir),   QMin['scratchdir']+'/DYSON/RUNFILE'])
      tasks.append(['link',QMin['scratchdir']+'/KEEP/ONEINT.%s' % (jobdir),    QMin['scratchdir']+'/DYSON/ONEINT'])
      tasks.append(['link',QMin['scratchdir']+'/KEEP/molcas.env.%s' % (jobdir),QMin['scratchdir']+'/DYSON/molcas.env'])
    elif QMin['integrals']=='dalton':
      tasks.append(['link',QMin['scratchdir']+'/KEEP/aoints.%s' % (jobdir),   QMin['scratchdir']+'/DYSON/aoints'])

    for multpair in QMin['multpairlist']:
      imult,jmult=tuple(multpair)

      # put the MO files into DYSON/
      ijobdir=QMin['multmap'][imult][0].replace('/','_')
      if QMin['molcas_rasscf']:
        tasks.append(['mo_convert',QMin['savedir']+'/molcas.RasOrb.%s' % ( ijobdir),QMin['savedir']+'/mocoef.%s' % (ijobdir),QMin['savedir']])
      tasks.append(['link',QMin['savedir']+'/mocoef.%s'     % (ijobdir),QMin['scratchdir']+'/DYSON/mocoef1'])

      jjobdir=QMin['multmap'][jmult][0].replace('/','_')
      if QMin['molcas_rasscf']:
        tasks.append(['mo_convert',QMin['savedir']+'/molcas.RasOrb.%s' % ( jjobdir),QMin['savedir']+'/mocoef.%s' % (jjobdir),QMin['savedir']])
      tasks.append(['link',QMin['savedir']+'/mocoef.%s'     % (jjobdir),QMin['scratchdir']+'/DYSON/mocoef2'])

      # link the determinants
      tasks.append(['link',QMin['savedir']+'/dets.%i' % (imult), QMin['scratchdir']+'/DYSON/dets1'])
      tasks.append(['link',QMin['savedir']+'/dets.%i' % (jmult), QMin['scratchdir']+'/DYSON/dets2'])

      # Run the Dyson calculation
      tasks.append(['write_dysonin',imult,jmult])
      tasks.append(['rundyson',imult,jmult])
      tasks.append(['get_DYSONoutput',imult,jmult])











  if 'backup' in QMin:
    tasks.append(['backupdata',QMin['backup']])

  if 'cleanup' in QMin:
    tasks.append(['cleanup',QMin['savedir']])
    tasks.append(['cleanup',QMin['scratchdir']])

  return tasks



# =============================================================================================== #
# =============================================================================================== #
# =========================================== SUBROUTINES TO RUNEVERYTING ======================= #
# =============================================================================================== #
# =============================================================================================== #

def mkdir(PATH):
  if os.path.exists(PATH):
    if os.path.isfile(PATH):
      print '%s exists and is a file!' % (PATH)
      sys.exit(81)
    else:
      ls=os.listdir(PATH)
      if not ls==[]:
        print 'INFO: %s exists and is a non-empty directory!' % (PATH)
        #sys.exit(82)
  else:
    os.makedirs(PATH)

# ======================================================================= #

def movetoold(QMin):
  # rename all eivectors, mocoef, molcas.input
  saveable=['dets','eivectors','mocoef','molcas.input','molcas.RasOrb','daltaoin']
  savedir=QMin['savedir']
  ls=os.listdir(savedir)
  if ls==[]:
    return
  for f in ls:
    f2=savedir+'/'+f
    if os.path.isfile(f2):
      if any( [ i in f for i in saveable ] ):
        if not 'old' in f:
          fdest=f2+'.old'
          shutil.move(f2,fdest)

# ======================================================================= #

#def link(PATH, NAME,crucial=True,force=False):
  ## do not create broken links
  #if not os.path.exists(PATH):
    #print 'Source %s does not exist, cannot create link!' % (PATH)
    #sys.exit(83)
  ## do ln -f only if NAME is already a link
  #if os.path.exists(NAME):
    #if os.path.islink(NAME) or force:
      #os.remove(NAME)
    #else:
      #print '%s exists, cannot create a link of the same name!' % (NAME)
      #if crucial:
        #sys.exit(84)
      #else:
        #return
  #if not os.path.exists(os.path.realpath(NAME)):
    #if os.path.islink(NAME):
    ## NAME is already a broken link
      #os.remove(NAME)
    #else:
      #return
  #os.symlink(PATH, NAME)

# ======================================================================= #

def link(PATH,NAME,crucial=True,force=True):
  # do not create broken links
  if not os.path.exists(PATH):
    print 'Source %s does not exist, cannot create link!' % (PATH)
    sys.exit(85)
  if os.path.islink(NAME):
    if not os.path.exists(NAME):
      # NAME is a broken link, remove it so that a new link can be made
      os.remove(NAME)
    else:
      # NAME is a symlink pointing to a valid file
      if force:
        # remove the link if forced to
        os.remove(NAME)
      else:
        print '%s exists, cannot create a link of the same name!' % (NAME)
        if crucial:
          sys.exit(86)
        else:
          return
  elif os.path.exists(NAME):
    # NAME is not a link. The interface will not overwrite files/directories with links, even with force=True
    print '%s exists, cannot create a link of the same name!' % (NAME)
    if crucial:
      sys.exit(87)
    else:
      return
  os.symlink(PATH, NAME)


# ======================================================================= #
def cleandir(directory):
  if DEBUG:
    print '===> Cleaning up directory %s\n' % (directory)
  for data in os.listdir(directory):
    path=directory+'/'+data
    if os.path.isfile(path) or os.path.islink(path):
      if DEBUG:
        print 'rm %s' % (path)
      try:
        os.remove(path)
      except OSError:
        print 'Could not remove file from directory: %s' % (path)
    else:
      if DEBUG:
        print ''
      cleandir(path)
      os.rmdir(path)
      if DEBUG:
        print 'rm %s' % (path)
  if DEBUG:
    print '\n'

# ======================================================================= #
def getmo(mofile,QMin):
  if os.path.exists(mofile):
    if QMin['molcas_rasscf']:
      tofile=os.path.join(QMin['scratchdir'],'JOB','molcas.RasOrb')
    else:
      tofile=os.path.join(QMin['scratchdir'],'JOB','mocoef')
    shutil.copy(mofile,tofile)
  else:
    print 'Could not find mocoef-file %s!' % (mofile)
    sys.exit(88)

# ======================================================================= #
def getinput(path,QMin):
  if os.path.exists(path) and os.path.isdir(path):
    ls=os.listdir(path)
    if not ls==[]:
      for f in ls:
        if os.path.isfile(path+'/'+f):
          shutil.copy(path+'/'+f,QMin['scratchdir']+'/JOB')
    else:
      print 'Template directory %s is empty!' % (path)
      sys.exit(89)
  else:
    print 'Template directory %s does not exist or is a file!' % (path)
    sys.exit(90)

# ======================================================================= #
def getmcinput(path,QMin):
  if os.path.exists(path) and os.path.isdir(path):
    ls=os.listdir(path)
    if not ls==[]:
      for f in ls:
        if 'mcdrtin' in f or 'mcscfin' in f or 'molcas.input' in f:
          if os.path.isfile(path+'/'+f):
            shutil.copy(path+'/'+f,QMin['scratchdir']+'/JOB')
    else:
      print 'Template directory %s is empty!' % (path)
      sys.exit(91)
  else:
    print 'Template directory %s does not exist or is a file!' % (path)
    sys.exit(92)

# ======================================================================= #
def checkinput(QMin):
  # check molcas.input for relevant keywords and necessary basis set information
  if QMin['integrals']=='seward':
    fname=QMin['scratchdir']+'/JOB/molcas.input'
    molcas=readfile(fname)
    necessary=['angmom']
    if 'soc' in QMin:
      necessary.append('amfi')
    if QMin['molcas_rasscf']:
      necessary.append('&rasscf')
    allthere=[False for i in necessary]
    for line in molcas:
      for i,nec in enumerate(necessary):
        if nec in line.lower():
          allthere[i]=True
    if not all(allthere):
      print '%s is missing one of the keywords %s!' % (fname,necessary)
      sys.exit(93)
    iline=0
    for atom in QMin['geo']:
      while True:
        if containsstring('%s\..*' % (atom[0].lower()),molcas[iline].lower()):
          break
        iline+=1
        if iline==len(molcas):
          print 'Basis set for atom %s not in molcas.input!' % (atom[0])
          sys.exit(94)
  # check cigrdin and mcscfin compatibility
  if 'grad' in QMin or 'nacdr' in QMin:
    # look for explicit HMC option in mcscfin
    fname=QMin['scratchdir']+'/JOB/mcscfin'
    mcscfin=readfile(fname)
    for line in mcscfin:
      if 'npath' in line.lower():
        f=line.split('=',1)[-1].split(',')
        explHMC=False
        for fi in f:
          if fi=='11':
            explHMC=True
    # adjust cigrdin according to explicit HMC option
    fname=QMin['scratchdir']+'/JOB/cigrdin'
    cigrdin=readfile(fname)
    for iline in range(len(cigrdin)):
      line=cigrdin[iline]
      if 'mdir' in line.lower():
        f=line.split(',')
        for fi in range(len(f)):
          if 'mdir' in f[fi].lower():
            f[fi]='mdir=%i' % ([1,0][explHMC])
        cigrdin[iline]=','.join(f)
    #with open(fname,'w') as f:
      #for line in cigrdin:
        #f.write(line)
    writefile(fname,cigrdin)

# ======================================================================= #
def writegeom(QMin):
  if 'unit' in QMin:
    if QMin['unit'][0]=='angstrom':
      factor=1.
    elif QMin['unit'][0]=='bohr':
      factor=au2a
    else:
      print 'Dont know input unit %s!' % (QMin['unit'][0])
      sys.exit(95)
  else:
    factor=1.

  fname=QMin['scratchdir']+'/JOB/geom.xyz'
  string='%i\n\n' % (QMin['natom'])
  for atom in QMin['geo']:
    string+=atom[0]
    for xyz in range(1,4):
      string+='  %f' % (atom[xyz]*factor)
      #g.write('  %f' % (atom[xyz]*factor))
    #g.write('\n')
    string+='\n'
  writefile(fname,string)
  #try:
    #with open(fname,'w') as g:
      #g.write('%i\n\n' % (QMin['natom']) )
      #for atom in QMin['geo']:
        #g.write(atom[0])
        #for xyz in range(1,4):
          #g.write('  %f' % (atom[xyz]*factor))
        #g.write('\n')
  #except IOError:
    #print 'Could not open geometry file %s!' % (fname)
    #sys.exit(96)

  os.chdir(QMin['scratchdir']+'/JOB')
  error=sp.call('%s/xyz2col.x < %s' % (QMin['columbus'],'geom.xyz'),shell=True)
  if error!=0:
    print 'xyz2col call failed!'
    sys.exit(97)
  os.chdir(QMin['pwd'])

# ======================================================================= #
def runProgram(string,workdir):
  prevdir=os.getcwd()
  if DEBUG:
    print workdir
  os.chdir(workdir)
  if PRINT or DEBUG:
    starttime=datetime.datetime.now()
    sys.stdout.write('%s\n\t%s' % (string,starttime))
    sys.stdout.flush()
  try:
    runerror=sp.call(string,shell=True)
  except OSError:
    print 'Call have had some serious problems:',OSError
    sys.exit(98)
  if PRINT or DEBUG:
    endtime=datetime.datetime.now()
    sys.stdout.write('\t%s\t\tRuntime: %s\t\tError Code: %i\n\n' % (endtime,endtime-starttime,runerror))
  os.chdir(prevdir)
  return runerror

# ======================================================================= #
def runCOLUMBUS(QMin):
  # setup
  workdir=QMin['scratchdir']+'/JOB'
  if DEBUG:
    string='%s -m %i -debug > runls' % (QMin['runc'],QMin['colmem'])
  else:
    string='%s -m %i > runls' % (QMin['runc'],QMin['colmem'])
  runerror=runProgram(string,workdir)

  # copy debug infos if crashed
  if runerror!=0:
    print 'COLUMBUS calculation crashed! Error code=%i' % (runerror)
    s=QMin['savedir']+'/COLUMBUS-debug/'
    s1=s
    i=0
    while os.path.exists(s1):
      i+=1
      s1=s+'%i/' % (i)
    if PRINT or DEBUG:
      print '=> Saving all text files from WORK directory to %s\n' % (s1)
    os.mkdir(s1)

    dirs=[
        os.path.join(QMin['scratchdir'],'JOB'),
        os.path.join(QMin['scratchdir'],'JOB/WORK'),
        os.path.join(QMin['scratchdir'],'JOB/LISTINGS')]
    maxsize=3*1024**2
    for d in dirs:
      ls=os.listdir(d)
      if DEBUG:
        print d
        print ls
      for i in ls:
        f=os.path.join(d,i)
        if os.stat(f)[6]<=maxsize and not os.path.isdir(f):
          try:
            shutil.copy(f,s1)
          except OSError:
            pass
          if PRINT or DEBUG:
            print i
    sys.exit(99)

# ======================================================================= #
def copymolden(job,QMin):
  # create directory
  jobdir=job.replace('/','_')
  moldendir=QMin['savedir']+'/MOLDEN/'+job
  if not os.path.isdir(moldendir):
    mkdir(moldendir)
  # save the molcas.input file
  f=QMin['scratchdir']+'/JOB/MOLDEN/molden_mo_mc.sp'
  fdest=moldendir+'/step_%s.molden' % (QMin['step'][0])
  shutil.move(f,fdest)

# ======================================================================= #
def save_data(job,QMin):
  jobdir=job.replace('/','_')
  # save all mocoef files, even if not necessary for other jobs
  if QMin['molcas_rasscf']:
    f=QMin['scratchdir']+'/JOB/MOCOEFS/molcas.RasOrb'
    fdest=QMin['savedir']+'/molcas.RasOrb.%s' % (jobdir)
  else:
    f=QMin['scratchdir']+'/JOB/MOCOEFS/mocoef_mc.sp'
    fdest=QMin['savedir']+'/mocoef.%s' % (jobdir)
  if os.path.isfile(f):
    shutil.copy(f,fdest)
  else:
    print 'Could not find %s!' % (f)
    sys.exit(100)
  if QMin['integrals']=='seward':
    # save the molcas.input file
    f=QMin['scratchdir']+'/JOB/WORK/molcas.input.seward'
    fdest=QMin['savedir']+'/molcas.input.seward.%s' % (jobdir)
    shutil.copy(f,fdest)
  elif QMin['integrals']=='dalton':
    f=QMin['scratchdir']+'/JOB/WORK/daltaoin'
    fdest=QMin['savedir']+'/daltaoin.%s' % (jobdir)
    shutil.copy(f,fdest)

# ======================================================================= #
def backupdata(backupdir,QMin):
  # save all files in savedir, except which have 'old' in their name
  ls=os.listdir(QMin['savedir'])
  for f in ls:
    ff=QMin['savedir']+'/'+f
    if os.path.isfile(ff) and not 'old' in ff:
      fdest=backupdir+'/'+f
      shutil.copy(ff,fdest)
  # save molden files
  if 'molden' in QMin:
    for job in QMin['joblist']:
      ff=os.path.join(QMin['savedir'],'MOLDEN',job,'step_%s.molden' % (QMin['step'][0]))
      fdest=os.path.join(backupdir,'step_%s.%s.molden' % (QMin['step'][0],job.replace('/','_')))
      shutil.copy(ff,fdest)


# ======================================================================= #
def keep_data(job,QMin):
  jobdir=job.replace('/','_')

  if 'ion' in QMin:
    if QMin['integrals']=='seward':
      # move the ONEINT, RUNFILE and molcas.env files
      f=QMin['scratchdir']+'/JOB/WORK/molcas.OneInt'
      fdest=QMin['scratchdir']+'/KEEP/ONEINT.%s' % (jobdir)
      shutil.move(f,fdest)
      f=QMin['scratchdir']+'/JOB/WORK/molcas.RunFile'
      fdest=QMin['scratchdir']+'/KEEP/RUNFILE.%s' % (jobdir)
      shutil.move(f,fdest)
      f=QMin['scratchdir']+'/JOB/WORK/molcas.env'
      fdest=QMin['scratchdir']+'/KEEP/molcas.env.%s' % (jobdir)
      shutil.move(f,fdest)
  elif QMin['integrals']=='dalton':
      f=QMin['scratchdir']+'/JOB/WORK/aoints'
      fdest=QMin['scratchdir']+'/KEEP/aoints.%s' % (jobdir)
      shutil.move(f,fdest)

# ======================================================================= #
def writemolcas(oldgeom,newgeom,molcasinfile):
  # =================
  atomlabel=re.compile('[a-zA-Z][a-zA-Z]?[1-9][0-9]*')
  numberlabel=re.compile('[-]?[0-9]+[.][0-9]*')
  #basislabel=re.compile('[bB]asis')
  basislabel=re.compile('^[bB]asis [sS]et|^[Ee]nd of [bB]asis')
  basisinfo=re.compile('^[a-zA-Z][a-zA-Z]?[.]')
  # =================
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
  # =================
  def isinfo(line):
    match=basislabel.search(line)
    if match:
      return True
    match=basisinfo.search(line)
    if match:
      return True
    else:
      return False

  geomold=readfile(oldgeom)
  geomnew=readfile(newgeom)

  string='&SEWARD &END\nONEOnly\nEXPErt\nNOGUessorb\nNODElete\nNODKroll\nNOAMfi\nMULTipoles\n0\n\n'

  # check for BASLIB keyword
  for i,line in enumerate(geomold):
    if 'baslib' in line.lower():
      string+='BASLIB\n'
      string+=geomold[i+1]+'\n'
      break

  # go through old geometry file
  atom=1
  for line in geomold:
    if isinfo(line):
      string+=line
    if isatom(line):
      parts=line.split()
      parts[0] = re.sub("\d+", "", parts[0])+str(atom)
      atom+=1
      line=' '.join(parts)
      string+=line+'\n'
  # go through the new geometry file
  for ln,line in enumerate(geomnew):
    if isinfo(line):
      if line!=geomold[ln]:
        print 'Inconsistent basis set information in line %i!' % (ln+1)
        #sys.exit(101)
      string+=line
    if isatom(line):
      parts=line.split()
      parts[0] = re.sub("\d+", "", parts[0])
      oldatom=re.sub("\d+", "", geomold[ln].split()[0])
      if parts[0]!=oldatom:
        print 'Different atoms in line %i!' % (ln+1)
        #sys.exit(102)
      parts[0]+=str(atom)
      atom+=1
      line=' '.join(parts)
      string+=line+'\n'
  string+='End of Input'

  # write to file
  #with open(molcasinfile,'w') as f:
    #f.write(string)
  writefile(molcasinfile,string)

# ======================================================================= #
def runmolcas(QMin):
  os.environ['MOLCAS_PROJECT']='molcas'
  os.environ['MOLCAS_WORKDIR']='molcas'
  workdir=QMin['scratchdir']+'/OVERLAP'
  prog=os.path.join(QMin['molcas'],'bin','molcas.exe')
  if not os.path.isfile(prog):
    prog=os.path.join(QMin['molcas'],'bin','pymolcas')
    if not os.path.isfile(prog):
      print 'No MOLCAS driver (molcas.exe or pymolcas) found in %s/bin/!' % (QMin['molcas'])
      sys.exit(103)
  string='%s %s &> %s' % (prog,'molcas.input','molcas.output')
  runerror=runProgram(string,workdir)
  if runerror!=0:
    print 'MOLCAS call not successful!'
    sys.exit(104)

# ======================================================================= #
def writedalton(oldgeom,newgeom,daltcomm,daltaoin,template):
  # get maxpri from template/daltcomm
  dcomm=readfile(os.path.join(template,'daltcomm'))
  for iline,line in enumerate(dcomm):
    if 'MAXPRI' in line:
      maxpri=int(dcomm[iline+1])
      break
  else:
    maxpri=20

  # write daltcomm file
  string="""**DALTONINPUT
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
**END OF INPUTS\n""" % maxpri
  writefile(daltcomm,string)

  # write daltaoin file
  dold=readfile(oldgeom)
  dnew=readfile(newgeom)
  string=''
  for i,line in enumerate(dold):
    if not i == 3:
      string+=line
    else:
      words = line.split()
      num = int(words[1]) * 2
      string+="s  %2i    0           0.10D-14\n"%num
  for i, line in enumerate(dnew):
    if i > 3:
      string+=line
  writefile(daltaoin,string)

# ======================================================================= #
def rundalton(QMin):
  workdir=QMin['scratchdir']+'/OVERLAP'
  string='%s/dalton.x -m %i > dalton.out 2> dalton.err' % (QMin['columbus'],QMin['colmem'])
  runerror=runProgram(string,workdir)
  if runerror!=0:
    print 'DALTON call not successful!'
    sys.exit(105)

# ======================================================================= #
def write_ciovin(QMin,imult):
  if 'ncore' in QMin and QMin['ncore']>=0:
    icore=QMin['ncore']
  else:
    icore=QMin['frozenmap'][imult]
  nintegrals={'seward':1,'dalton':2}[QMin['integrals']]
  string='''a_mo=mocoef1
b_mo=mocoef2
a_det=dets1
b_det=dets2
ao_read=%i
ncore=%i''' % (nintegrals, icore)
  writefile(os.path.join(QMin['scratchdir'],'OVERLAP','cioverlap.in'),string)

# ======================================================================= #
def runcioverlap(QMin,mult):
  workdir=QMin['scratchdir']+'/OVERLAP'
  os.environ['OMP_NUM_THREADS']=str(QMin['ncpu'])
  string='%s -f cioverlap.in -m %i &> cioverlap.out' % (QMin['wfoverlap'],QMin['colmem'])
  runerror=runProgram(string,workdir)
  if runerror!=0:
    print 'cioverlap call not successful!'
    sys.exit(106)

  # backup dyson.out
  if 'backup' in QMin:
    f=QMin['scratchdir']+'/OVERLAP/cioverlap.out'
    fdest=QMin['backup']+'/cioverlap.out.m%i' % (mult)
    shutil.copy(f,fdest)

# ======================================================================= #
def write_dysonin(QMin,imult,jmult):
  if 'ncore' in QMin and QMin['ncore']>=0:
    icore=QMin['ncore']
  else:
    icore=min(QMin['frozenmap'][imult],QMin['frozenmap'][jmult])
  nintegrals={'seward':1,'dalton':2}[QMin['integrals']]
  string='''a_mo=mocoef1
b_mo=mocoef2
a_det=dets1
b_det=dets2
same_aos
ao_read=%i
ncore=%i
''' % (nintegrals, icore)
  if 'ndocc' in QMin:
    string+='ndocc=%i\n' % QMin['ndocc']
  writefile(os.path.join(QMin['scratchdir'],'DYSON','dyson.in'),string)

# ======================================================================= #
def rundyson(QMin,imult,jmult):
  workdir=QMin['scratchdir']+'/DYSON'
  os.environ['OMP_NUM_THREADS']=str(QMin['ncpu'])
  string='%s -f dyson.in -m %i &> dyson.out' % (QMin['wfoverlap'],QMin['colmem'])
  runerror=runProgram(string,workdir)
  if runerror!=0:
    print 'cioverlap call not successful!'
    sys.exit(107)

  # backup dyson.out
  if 'backup' in QMin:
    f=QMin['scratchdir']+'/DYSON/dyson.out'
    fdest=QMin['backup']+'/dyson.out.m%i.m%i' % (imult,jmult)
    shutil.copy(f,fdest)

## ======================================================================= #
#def make_dets(job,QMin):
  ## find all multiplicities which are relevant for dyson calculations
  #mults=set([])
  #for i in QMin['ionlist']:
    #m1,s1,m2,s1=tuple(i)
    #mults.add(m1)
    #mults.add(m2)
  #for mult in mults:
    #job2,drt=tuple(QMin['multmap'][mult])
    #if job2==job:
      ## remove "dum" dummy files of civecconsolidate
      #ls=os.listdir(QMin['scratchdir']+'/JOB/WORK')
      #for ifile in ls:
        #if 'dum' in ifile:
          #os.remove(QMin['scratchdir']+'/JOB/WORK/'+ifile)
      ## run civecconsolidate in JOB/WORK for drt
      #workdir=QMin['scratchdir']+'/JOB/WORK'
      #eivecfile='eivectors.red.drt%i.sp' % (drt)
      #slaterfile='slaterfile.red.drt%i.sp' % (drt)
      #if 'civecconsolidate' in QMin:
        #string=QMin['civecconsolidate']
      #else:
        #string=QMin['columbus']+'/civecconsolidate'
      #string+=' -d %f %s %s dum1 dum2 dum3 < civecconsolidate.in >> civecconsolidate.out' % (QMin['civecthres'],eivecfile,slaterfile)
      #runerror=runProgram(string,workdir)
      #if runerror!=0:
        #print 'civecconsolidate call not successful!'
        #sys.exit(108)
      ## move determinant files to KEEP
      #for istate in range(1,QMin['states'][mult-1]+1):
        #f=QMin['scratchdir']+'/JOB/WORK/determinants.vec_%i' % (istate)
        #fdest=QMin['scratchdir']+'/KEEP/determinants.vec_s%i_m%i' % (istate,mult)
        #shutil.move(f,fdest)

# ======================================================================= #
def make_dets_new(job,QMin):
  # find all multiplicities and drts for current job
  mults=QMin['multmap'][job]
  # for all DRTs, call read_civfl
  os.chdir(os.path.join(QMin['scratchdir'],'JOB','WORK'))

  # delete cidrtfl, if it exists
  if not QMin['socimap'][job]==-1:
    filename='cidrtfl'
    if os.path.isfile(filename):
      os.remove(filename)

  for imult in mults:

    if QMin['socimap'][job]==1 or QMin['socimap'][job]==0:
      idrt=QMin['multmap'][imult][1]
      f='civfl.drt%i' % (idrt)
      fdest='civfl'
      link(f,fdest,force=True)
      f='cidrtfl.%i' % (idrt)
      fdest='cidrtfl'
      link(f,fdest,force=True)
      f='civout.drt%i' % (idrt)
      fdest='civout'
      link(f,fdest,force=True)
    elif QMin['socimap'][job]==-1:
      # in the case of single-drt calculations the files already have the correct names
      pass

    #print 'Running cipc.x for multiplicity: %i' % (imult)
    ca=civfl_ana(QMin['wfthres'],DEBUG,QMin['columbus'])
    for istate in range(1,1+QMin['states'][imult-1]):
      #ca.call_cipc(istate,ms=0.,mem=QMin['colmem'])
      ca.call_cipc(istate,ms=0.5*(imult-1),mem=QMin['colmem'])
    ca.write_det_file(QMin['states'][imult-1])
    f=os.path.join(QMin['scratchdir'],'JOB','WORK','dets')
    fdest=os.path.join(QMin['savedir'],'dets.%i' % (imult))
    shutil.move(f,fdest)
    if not 'frozenmap' in QMin:
      QMin['frozenmap']={}
    QMin['frozenmap'][imult]=ca.nfct
  print ''
  return QMin

# ======================================================================= #
def mo_convert(fromfile,tofile,directory,QMin):
  # write mo_convertin
  string='''&input
formatin=1
formatout=0
mo_in="%s"
mo_out="%s"
&end
''' % (fromfile,tofile)
  filename=os.path.join(directory,'mo_convertin')
  writefile(filename,string)

  # run mo_convert.x
  string='%s/mo_convert.x' % (QMin['columbus'])
  runerror=runProgram(string,directory)
  if runerror!=0:
    print 'mo_convert.x call not successful!'
    sys.exit(109)

  # remove in and ls files
  toremove=['mo_convertin','mo_convertls']
  for i in toremove:
    filename=os.path.join(directory,i)
    if os.path.isfile(filename):
      os.remove(filename)

# ======================================================================= #
def get_info_from_det(filename):
  lines=readfile(filename)
  line=lines[0].split()
  nstate=int(line[0])
  norb=int(line[1])
  ndet=int(line[2])

  nelec_map={'d':2, 'a':1, 'b':1, 'e':0}
  det=lines[1].split()[0]
  nelec=0
  for i in det:
    nelec+=nelec_map[i]
  return nelec,norb,ndet,nstate

# ======================================================================= #
#def do_dyson(QMin,QMout,pair):
  ## write input file
  #m1,s1,m2,s2=tuple(pair)
  #mocoef1=QMin['savedir']+'/mocoef.%s' % (QMin['multmap'][m1][0].replace('/','_'))
  #mocoef2=QMin['savedir']+'/mocoef.%s' % (QMin['multmap'][m2][0].replace('/','_'))
  #file1=QMin['savedir']+'/dets.%i' % (m1)
  #file2=QMin['savedir']+'/dets.%i' % (m2)
  #nelec1,norb1,ndet1,nstate1=get_info_from_det(file1)
  #nelec2,norb2,ndet2,nstate2=get_info_from_det(file2)

  #if s1>nstate1:
    #print 'Not enough states in %s to calculate Dyson pair %s!' % (file1,pair)
    #sys.exit(110)
  #if s2>nstate2:
    #print 'Not enough states in %s to calculate Dyson pair %s!' % (file2,pair)
    #sys.exit(111)

  #if nelec1>nelec2:
    ## m1,s1 is reference
    #infos=(mocoef1,mocoef2,file1,file2,QMin['dysonthres'],QMin['integrals'],s1,s2)
    #jobdir=QMin['multmap'][m1][0].replace('/','_')
  #else:
    ## m2,s2 is reference
    #infos=(mocoef2,mocoef1,file2,file1,QMin['dysonthres'],QMin['integrals'],s2,s1)
    #jobdir=QMin['multmap'][m2][0].replace('/','_')

  #string='''ref_mo=%s
#ion_mo=%s
#ref_det=%s
#ion_det=%s
#c2_threshold=%.12f
#integrals=%s
#refstate=%i
#ionstate=%i
#''' % infos
  #if 'ncore' in QMin:
    #string+='ncore=%i\n' % QMin['ncore']
  #elif 'frozenmap' in QMin:
    #string+='ncore=%i\n' % (min(QMin['frozenmap'][m1],QMin['frozenmap'][m2]))
  #if 'ndocc' in QMin:
    #string+='ndocc=%i\n' % QMin['ndocc']

  #writefile(QMin['scratchdir']+'/DYSON/dyson.input',string)



  #if QMin['integrals']=='seward':
    ## link RUNFILE, ONEINT, molcas.env
    #f=QMin['scratchdir']+'/KEEP/ONEINT.%s' % (jobdir)
    #fdest=QMin['scratchdir']+'/DYSON/ONEINT'
    #link(f,fdest)

    #f=QMin['scratchdir']+'/KEEP/RUNFILE.%s' % (jobdir)
    #fdest=QMin['scratchdir']+'/DYSON/RUNFILE'
    #link(f,fdest)

    #f=QMin['scratchdir']+'/KEEP/molcas.env.%s' % (jobdir)
    #fdest=QMin['scratchdir']+'/DYSON/molcas.env'
    #link(f,fdest)

  ## run dyson
  #os.environ['OMP_NUM_THREADS']=str(QMin['ncpu'])
  #workdir=QMin['scratchdir']+'/DYSON'
  #string=QMin['dyson']+' > dyson.output          # Pair: %s' % (pair)
  #runerror=runProgram(string,workdir)
  #if runerror!=0:
    #print 'dyson call not successful!'
    #sys.exit(112)

  ## extract result
  ##with open(QMin['scratchdir']+'/DYSON/dyson.output') as f:
    ##dysonout=f.readlines()
  #dysonout=readfile(QMin['scratchdir']+'/DYSON/dyson.output')
  #for line in dysonout:
    #if '<psid|psid>' in line:
      #dysonnorm=float(line.split()[2])
      #break

  ## backup dyson.out
  #if 'backup' in QMin:
    #f=QMin['scratchdir']+'/DYSON/dyson.output'
    #fdest=QMin['backup']+'/dyson.output.m%is%i.m%is%i' % tuple(pair)
    #shutil.move(f,fdest)

  ## update QMout
  #if not 'prop' in QMout:
    #QMout['prop']=makecmatrix(QMin['nmstates'],QMin['nmstates'])
  #for i in range(QMin['nmstates']):
    #for j in range(QMin['nmstates']):
      #m3,s3,ms3=tuple(QMin['statemap'][i+1])
      #m4,s4,ms4=tuple(QMin['statemap'][j+1])
      #if (m3,s3,m4,s4)==(m1,s1,m2,s2) or (m4,s4,m3,s3)==(m1,s1,m2,s2):
        #if abs(ms3-ms4)==0.5:
          #QMout['prop'][i][j]=complex(dysonnorm,0.)
          #QMout['prop'][j][i]=complex(dysonnorm,0.)

  #return QMout




# ======================================================================= #
def runeverything(tasks, QMin):

  if PRINT or DEBUG:
    print '=============> Entering RUN section <=============\n\n'

  QMout={}
  #states=QMin['states']
  #nstates=QMin['nstates']
  #nmstates=QMin['nmstates']
  for task in tasks:
    if DEBUG:
      print task
    if task[0]=='movetoold':
      movetoold(QMin)
    if task[0]=='mkdir':
      mkdir(task[1])
    if task[0]=='link':
      if len(task)==4:
        link(task[1],task[2],task[3])
      else:
        link(task[1],task[2])
    if task[0]=='getmo':
      getmo(task[1],QMin)
    if task[0]=='getinput':
      getinput(task[1],QMin)
    if task[0]=='getmcinput':
      getmcinput(task[1],QMin)
    if task[0]=='checkinput':
      checkinput(QMin)
    if task[0]=='writegeom':
      writegeom(QMin)
    if task[0]=='write':
      writefile(task[1],task[2])
    if task[0]=='runc':
      runCOLUMBUS(QMin)
    if task[0]=='save_data':
      save_data(task[1],QMin)
    if task[0]=='keep_data':
      keep_data(task[1],QMin)
    if task[0]=='backupdata':
      backupdata(task[1],QMin)
    if task[0]=='copymolden':
      copymolden(task[1],QMin)
    #if task[0]=='make_dets':
      #make_dets(task[1],QMin)
    if task[0]=='make_dets_new':
      QMin=make_dets_new(task[1],QMin)
    if task[0]=='cleanup':
      cleandir(task[1])
    if task[0]=='writemolcas':
      writemolcas(task[1],task[2],task[3])
    if task[0]=='runmolcas':
      runmolcas(QMin)
    if task[0]=='writedalton':
      writedalton(task[1],task[2],task[3],task[4],task[5])
    if task[0]=='rundalton':
      rundalton(QMin)
    if task[0]=='get_COLout':
      QMout=get_COLout(QMin,QMout,task[1])
    if task[0]=='mo_convert':
      mo_convert(task[1],task[2],task[3],QMin)

    if task[0]=='write_ciovin':
      write_ciovin(QMin,task[1])
    if task[0]=='runcioverlap':
      runcioverlap(QMin,task[1])
    if task[0]=='get_CIOoutput':
      QMout=get_CIOoutput(QMin,QMout,task[1],task[2])

    if task[0]=='write_dysonin':
      write_dysonin(QMin,task[1],task[2])
    if task[0]=='rundyson':
      rundyson(QMin,task[1],task[2])
    if task[0]=='get_DYSONoutput':
      QMout=get_DYSONoutput(QMin,QMout,task[1],task[2])

    #if task[0]=='dyson':
      #QMout=do_dyson(QMin,QMout,task[1])

  # if no dyson pairs were calculated because of selection rules, put an empty matrix
  if not 'prop' in QMout and 'ion' in QMin:
    QMout['prop']=makecmatrix(QMin['nmstates'],QMin['nmstates'])

  # Phases from overlaps
  if 'phases' in QMin:
    if not 'phases' in QMout:
      QMout['phases']=[ complex(1.,0.) for i in range(QMin['nmstates']) ]
    if 'overlap' in QMout:
      for i in range(QMin['nmstates']):
        if QMout['overlap'][i][i].real<0.:
          QMout['phases'][i]=complex(-1.,0.)

  return QMout





# ======================================================================= #
# ======================================================================= #
# ======================================================================= #
class civfl_ana:
    def __init__(self, maxsqnorm=1.0, debug=False, columbus=os.environ['COLUMBUS']):
        self.det_dict = {} # dictionary with determinant strings and cicoefficient information
        self.nmot = -1  # number of MOs
        self.niot = -1  # number of internal orbs
        self.nfct = -1  # number of frozen orbs
        self.nfvt = -1  # number of frozen virtuals
        self.ncsf = -1
        self.maxsqnorm = maxsqnorm
        self.sqcinorms = {} # CI-norms
        self.debug = debug
        self.columbus=columbus
# ================================================== #
    def read_cipcls(self, istate):
        """
        Read pre-generated cipcls files.
        """
        fname = 'cipcls.det%i'%istate
        #print "Reading %s ..."%fname
        cip = open(fname, 'r')
        fstring = cip.read()
        cip.close()
        self.read_cipinfo(istate, fstring)
# ================================================== #
    def call_cipc(self, istate, ms="0", csfbuf=50000, mem=1000):
        """
        Call cipc.x and analyse the information on the fly.
        The input is batched to limit the amount of memory used.
        The batch size is controlled by csfbuf.
        """
        command = ["%s/cipc.x"%self.columbus, "-m", "%i"%mem]
        istart = 1
        maxiter=20000
        for i in xrange(maxiter):
            iend = istart + csfbuf
            cipstr  = "2\n4\n1\n%s\n0\n"%ms # initialize determinant print out
            cipstr += "7\n5\n%i\n0\n"%istate # read the coefficients
            cipstr += "1\n%i %i\n0/\n"%(istart, iend) # first and last CSF to print
            cipstr += "0\n" # finish
            print "%s/cipc.x for state %i, CSFs %i to %i"%(self.columbus, istate, istart, iend)
            starttime=datetime.datetime.now()
            sys.stdout.write('\t%s' % (starttime))
            cipx = sp.Popen(command, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
            cipout, ciperr = cipx.communicate(cipstr)
            if self.debug:
                open('pycipcin.st%i.%i'%(istate,i), 'w').write(cipstr)
                open('pycipcls.st%i.%i'%(istate,i), 'w').write(cipout)
            if not 'end of cipc' in ciperr:
                print " ERROR in cipc.x during determinant generation!"
                #print "\n Standard output:\n", cipout
                print "\n Standard error:\n", ciperr, 'Exit code:',cipx.returncode
                print cipout
                sys.exit(113)
            if self.debug:
                print " cipc.x finished succesfully"
            self.read_cipinfo(istate, cipout)
            endtime=datetime.datetime.now()
            sys.stdout.write('\t%s\t\tRuntime: %s\n\n' % (endtime,endtime-starttime))
            if (iend > self.ncsf) or (self.sqcinorms[istate] > self.maxsqnorm):
                if self.debug:
                    print "Finished, iend = %i, sqcinorm = %.4f.\n"%(iend, self.sqcinorms[istate])
                break
            else:
                istart = iend + 1
        else:
            print '%i CSFs read, maxiter reached.' % (csfbuf*maxiter)
            sys.exit(114)
# ================================================== #
    def read_cipinfo(self, istate, fstring):
        if istate in self.sqcinorms:
            sqcinorm = self.sqcinorms[istate]
        else:
            sqcinorm = 0.
        flist = fstring.split("\n")
        flist.reverse()
        csfsec = False
        n_ext = 0; ext1 = 0; ext2 = 0
        while(True):
            line = flist.pop()
            if 'drt header information' in line:
                line = flist.pop()
                line = flist.pop()
                words = line.replace('=',' ').split()
                self.nmot = int(words[1])
                self.niot = int(words[3])
                self.nfct = int(words[5])
                self.nfvt = int(words[7])
                #print 'Orbital information parsed:'
                if self.debug:
                    print '  nmot = %i, niot = %i, nfct = %i, nfvt = %i'%(self.nmot, self.niot, self.nfct, self.nfvt)
            elif 'ncsft:' in line:
                self.ncsf = int(line.split()[-1])
                #print '  ncsf = %i'%self.ncsf
            if 'indcsf' in line:
                csfsec = True
                line = flist.pop()
                continue
            elif ('csfs were printed in this range' in line):
                if self.debug:
                    print "All CSFs of this batch read in.\n"
                break
            if not csfsec: continue
            if len(line) == 0:
                if self.debug:
                    print "All CSFs of this batch read in.\n"
                break
            words = line.split()
            if ('idet' in line):
                if sqcinorm > self.maxsqnorm:
                    if self.debug:
                        print "Stopping at sqcinorm = %.4f"%sqcinorm
                    break
                else:
                    continue
            elif (len(words) > 3):
                wtype = words[3]
                n_ext_el = {'z*':0, 'z':0, 'y':1, 'x':2, 'w':2}[wtype]
                n_exto = line.count(':') # number of external orbitals
                ext1 = ext2 = -1
                rwords = line.replace(':', ' ').split()
                if n_exto == 1:
                    ext1 = int(rwords[-2])
                elif n_exto == 2:
                    ext1 = int(rwords[-4])
                    ext2 = int(rwords[-2])
                CSFstr = words[-1]
                nel = 2*CSFstr.count('3') + CSFstr.count('1') + CSFstr.count('2')
                if self.debug:
                    print "-> %2s-CSF %s: nel = %i, n_ext_el = %i, n_exto = %i, ext1 = %2i, ext2 = %2i"%(wtype, CSFstr, nel, n_ext_el, n_exto, ext1, ext2)
            else:
                coeff = float(words[1])
                # For even electron systems there is phase change if one external
                #   electron is moved from the front (where it is in the DRT)
                #   to the back (where we expect it to be).
                # With two external electrons, this cancels out.
                # For odd electron systems, the phases should never change.
                if (nel%2==0) and (n_ext_el==1): coeff = -coeff
                #if n_ext_el==1: coeff = -coeff
                det = self.det_string(words[-1], n_exto, ext1, ext2)
                if self.debug:
                    print "%25s -> %s: % .10f"%(words[-1], det, coeff)
                if det in self.det_dict:
                    if istate in self.det_dict[det]:
                        self.det_dict[det][istate] += coeff
                    else:
                        self.det_dict[det][istate]  = coeff
                else:
                    self.det_dict[det] = {istate:coeff}
                sqcinorm += coeff **2
        self.sqcinorms[istate] = sqcinorm
# ================================================== #
    def det_string(self, cipstr, n_ext, ext1, ext2):
        retstr  = self.nfct * 'd'
        retstr += self.det_labels(cipstr[n_ext:])
        for iorb in xrange(self.nfct + self.niot + 1, self.nmot + 1):
            if   iorb == ext1:
                retstr += self.det_labels(cipstr[0])
            elif iorb == ext2:
                retstr += self.det_labels(cipstr[1])
            else:
                retstr += 'e'
        return retstr
# ================================================== #
    def det_labels(self, cipstr):
        return cipstr.replace('#','d').replace('+','a').replace('-','b').replace('.','e')
# ================================================== #
    def sort_key(self, key):
        """
        For specifying the sorting order of the determinants.
        """
        return key.replace('d', '0').replace('a', '1').replace('b', '1')
# ================================================== #
    def sort_key2(self, key):
        """
        For specifying the sorting order of the determinants.
        """
        return key.replace('d', '0').replace('a', '0').replace('b', '1').replace('e', '1')
# ================================================== #
    def write_det_file(self, nstate, wname='dets', wform=' % 14.10f'):
        wf = open(wname, 'w')
        wf.write("%i %i %i\n"%(nstate, self.nmot, len(self.det_dict)))
        for det in sorted(sorted(self.det_dict, key=self.sort_key2), key=self.sort_key):
            wf.write(det)
            for istate in xrange(1, nstate+1):
                try:
                    coeff = self.det_dict[det][istate]
                except KeyError:
                    coeff = 0.
                wf.write(wform%coeff)
            wf.write('\n')
        wf.close()
        if self.debug:
            print "File %s written."%wname






# ========================== Main Code =============================== #
def main():
  '''This script realises an interface between the semi-classical dynamics code SHARC and the quantum chemistry program MOLPRO 2012. It allows the automatised calculation of non-relativistic and spin-orbit Hamiltonians, Dipole moments, gradients and non-adiabatic couplings at the CASSCF level of theory for an arbitrary number of states of different multiplicities. It also includes a small number of MOLPRO error handling capabilities (restarting non-converged calculations etc.).

  Input is realised through two files and a number of environment variables.

  QM.in:
    This file contains all information which are known to SHARC and which are independent of the used quantum chemistry code. This includes the current geometry and velocity, the number of states/multiplicities, the time step and the kind of quantities to be calculated.

  MOLPRO.template:
    This file is a minimal MOLPRO input file containing all molecule-specific parameters, like memory requirement, basis set, Douglas-Kroll-Hess transformation, active space and state-averaging.

  Environment variables:
    Additional information, which are necessary to run MOLPRO, but which do not actually belong in a MOLPRO input file.
    The necessary variables are:
      * QMEXE: is the path to the MOLPRO executable
      * SCRATCHDIR: is the path to a scratch directory for fast I/O Operations.
    Some optional variables are concerned with MOLPRO error handling (defaults in parenthesis):
      * GRADACCUDEFAULT: default accuracy for MOLPRO CPMCSCF (1e-7)
      * GRADACCUMAX: loosest allowed accuracy for MOLPRO CPMCSCF (1e-2)
      * GRADACCUSTEP: factor for decreasing the accuracy for MOLPRO CPMCSCF (1e-1)
      * CHECKNACS: check whether non-adiabatic couplings are corrupted by intruder states (False)
      * CHECKNACS_MRCIO: threshold for intruder state detection, see mrcioverlapsok() (0.85)'''

  # Retrieve PRINT and DEBUG
  try:
    envPRINT=os.getenv('SH2COL_PRINT')
    if envPRINT and envPRINT.lower()=='false':
      global PRINT
      PRINT=False
    envDEBUG=os.getenv('SH2COL_DEBUG')
    if envDEBUG and envDEBUG.lower()=='true':
      global DEBUG
      DEBUG=True
  except ValueError:
    print 'PRINT or DEBUG environment variables do not evaluate to logical values!'
    sys.exit(115)

  # Process Command line arguments
  if len(sys.argv)!=2:
    print 'Usage:\n./SHARC_COLUMBUS.py <QMin>\n'
    print 'version:',version
    print 'date:',versiondate
    print 'changelog:\n',changelogstring
    sys.exit(116)
  QMinfilename=sys.argv[1]

  # Print header
  printheader()

  # Read QMinfile
  QMin=readQMin(QMinfilename)
  printQMin(QMin)

  # Process Tasks
  Tasks=gettasks(QMin)
  if DEBUG:
    printtasks(Tasks)

  # Do the COLUMBUS, cioverlaps and dyson runs, extract QMout
  QMout=runeverything(Tasks,QMin)

  printQMout(QMin,QMout)

  # Measure time
  runtime=measuretime()
  QMout['runtime']=runtime

  # Write QMout
  writeQMout(QMin,QMout,QMinfilename)

  if PRINT or DEBUG:
    print datetime.datetime.now()
    print '#================ END ================#'

if __name__ == '__main__':
    main()
