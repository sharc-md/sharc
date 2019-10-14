#!/usr/bin/env python2
"""

Standalone script for performing SHARC dynamics, reading a single QMout file at the beginning

@version: 0.1.0
@author: Maximilian F.S.J. Menger (slight modifications of SHARC_LVC)
@description: Python Module for SHARC Dynamics using a single QMout file

"""
from __future__ import print_function

import shutil
import sys
import os

import numpy
import numpy as np

import pprint

from sharc.pysharc.interface import SHARC_INTERFACE

# ******************************
#
# Helper functions
#
# ******************************

# QMout targets


def check_keys(lines, targets):
    keys = []
    for line in lines:
        for t in targets:
            if targets[t]['line'] in line: keys.append(t)
    return keys 
# ======================================================================= #
def readfile(filename):
  try:
    f=open(filename)
    out=f.readlines()
    f.close()
  except IOError:
    print('File %s does not exist!' % filename)
    sys.exit(12)
  return out
# ======================================================================= #
def read_QMout(path,nstates,natom,request=None):
  targets={'h':     {  'flag' : 1,
                         'type' : complex,
                         'dim'  :  (nstates,nstates),
                         'line' : '! 1 Hamiltonian Matrix',
                      },
           'dm':      {  'flag' : 2,
                         'type' : complex,
                         'dim'  :  (3,nstates,nstates),
                         'line' : '! 2 Dipole Moment Matrices',
                      },
           'grad':    {  'flag' : 3,
                         'type' : float,
                         'dim'  :  (nstates,natom,3),
                         'line' : '! 3 Gradient Vectors',
                      },
           'nacdr':   {  'flag': 5,
                         'type': float,
                         'dim' :  (nstates,nstates,natom,3),
                         'line': 'Not Implemented!',
                      },
           'overlap':   {  'flag': 6,
                         'type': complex,
                         'dim' :  (nstates,nstates),
                         'line': '! 6 Overlap matrix',
                          },

           }
  # read QM.out
  lines=readfile(path)
  if request is None:
    request = check_keys(lines, targets)

  # obtain all targets
  QMout={}
  for t in targets:
    if t in request:
      iline=-1
      while True:
        iline+=1
        if iline>=len(lines):
          print('Could not find target %s with flag %i in file %s!' % (t,targets[t]['flag'],path))
          sys.exit(11)
        line=lines[iline]
        if '! %i' % (targets[t]['flag']) in line:
          break
      values=[]
      # =========== single matrix
      if len(targets[t]['dim'])==2:
        iline+=1
        for irow in range(targets[t]['dim'][0]):
          iline+=1
          line=lines[iline].split()
          if targets[t]['type']==complex:
            row=[ complex(float(line[2*i]),float(line[2*i+1])) for i in range(targets[t]['dim'][1]) ]
          elif targets[t]['type']==float:
            row=[ float(line[i]) for i in range(targets[t]['dim'][1]) ]
          values.append(row)
      # =========== list of matrices
      elif len(targets[t]['dim'])==3:
        for iblocks in range(targets[t]['dim'][0]):
          iline+=1
          block=[]
          for irow in range(targets[t]['dim'][1]):
            iline+=1
            line=lines[iline].split()
            if targets[t]['type']==complex:
              row=[ complex(float(line[2*i]),float(line[2*i+1])) for i in range(targets[t]['dim'][2]) ]
            elif targets[t]['type']==float:
              row=[ float(line[i]) for i in range(targets[t]['dim'][2]) ]
            block.append(row)
          values.append(block)
      # =========== matrix of matrices
      elif len(targets[t]['dim'])==4:
        for iblocks in range(targets[t]['dim'][0]):
          sblock=[]
          for jblocks in range(targets[t]['dim'][1]):
            iline+=1
            block=[]
            for irow in range(targets[t]['dim'][2]):
              iline+=1
              line=lines[iline].split()
              if targets[t]['type']==complex:
                row=[ complex(float(line[2*i]),float(line[2*i+1])) for i in range(targets[t]['dim'][3]) ]
              elif targets[t]['type']==float:
                row=[ float(line[i]) for i in range(targets[t]['dim'][3]) ]
              block.append(row)
            sblock.append(block)
          values.append(sblock)
      QMout[t]=values

  #pprint.pprint(QMout)
  return QMout


class SHARC_QMOUT(SHARC_INTERFACE):
    """ Class for SHARC QMout """
    # Name of the interface
    interface = 'QMout'
    # store atom ids
    save_atids = False
    # store atom names
    save_atnames = False
    # accepted units:  0 : Bohr, 1 : Angstrom
    iunit = 0
    # not supported keys
    not_supported = ['nacdt', 'dmdr' ]

    def do_qm_job(self, tasks, Crd):
        """ Here you should perform all your qm calculations depending on the tasks, that were asked """
        return self.storage['QMout']


    def readParameter(self, fname, *args, **kwargs):
        """ read QMout file for the calculation here """
        fname='QM/QM.out'
        self.storage['QMout'] = read_QMout(fname,self.states['nmstates'],self.NAtoms,kwargs['QMout_keys'])


def getCommandoLine():
    """
        Get Commando line option with argpase
        
    """

    import argparse

    parser = argparse.ArgumentParser("Perform SHARC QMout calculations")
    parser.add_argument("input", metavar="FILE",type=str,
                        default="input", nargs='?',
                        help="input file")
    parser.add_argument("QMout", metavar="FILE",type=str,
                        default="QM/QM.out",
                        help="QMout file for the dynamics")

    args = parser.parse_args()

    return args.input, args.QMout


def main():
    """
        Main Function if program is called as standalone 

    """
    inp_file, qmout_file = getCommandoLine()

   # init SHARC_LVC class
    qmout = SHARC_QMOUT()
    # get 
#    # run sharc dynamics
    qmout.run_sharc(inp_file, qmout_file, QMout_keys = None)

if __name__ == "__main__":
    main()



