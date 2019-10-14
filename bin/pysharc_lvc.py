#!/usr/bin/env python2
"""

Standalone script for performing SHARC/LVC dynamics.

@version: 0.1.0
@author: Felix Plasser, Maximilian F.S.J. Menger (slight modifications of SHARC_LVC)
@description: Python Module for SHARC Dynamics using the LVC model for the molecular system.

"""
from __future__ import print_function

import shutil
import sys
import os

import numpy
import numpy as np

from sharc.pysharc.interface import SHARC_INTERFACE

# ******************************
#
# SHARC_LVC.py functions
#
# ******************************

def itnmstates(states):
    '''Takes an array of the number of states in each multiplicity and 
    generates an iterator over all states specified. 
    Iterates also over all MS values of all states.
    
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
    3 float: MS value'''
    
    for i in range(len(states)):
        if states[i]<1:
            continue
        for k in range(i+1):
            for j in range(states[i]):
                yield i+1,j+1,k-i/2.
    return

def checkscratch(SCRATCHDIR):
    '''Checks whether SCRATCHDIR is a file or directory. 
    If a file, it quits with exit code 1, if its a directory, it passes. 
    If SCRATCHDIR does not exist, tries to create it.

    Arguments:
    1 string: path to SCRATCHDIR
    '''

    exist=os.path.exists(SCRATCHDIR)
    if exist:
        isfile=os.path.isfile(SCRATCHDIR)
        if isfile:
            print('$SCRATCHDIR=%s exists and is a file!' % (SCRATCHDIR))
            sys.exit(16)
    else:
        try:
            os.makedirs(SCRATCHDIR)
        except OSError:
            print('Can not create SCRATCHDIR=%s\n' % (SCRATCHDIR))
            sys.exit(17)
# =========================================================
def find_lines(nlines, match, strings):
    smatch = match.lower().split()
    iline = -1
    while True:
        iline += 1
        if iline == len(strings):
            return []
        line = strings[iline].lower().split()
        if tuple(line) == tuple(smatch):
            return strings[iline+1:iline+nlines+1]

# =========================================================
def read_LVC_mat(nmstates, header, rfile):
    mat=[ [ complex(0.,0.) for i in range(nmstates) ] for j in range(nmstates) ]

    # real part
    tmp=find_lines(nmstates,header+' R',rfile)
    if not tmp==[]:
        for i, line in enumerate(tmp):
            for j, val in enumerate(line.split()):
                mat[i][j] += float(val)

    # imaginary part
    tmp=find_lines(nmstates,header+' I',rfile)
    if not tmp==[]:
        for i, line in enumerate(tmp):
            for j, val in enumerate(line.split()):
                mat[i][j] += float(val) * 1j

    return mat
# getQMout
# =========================================================
def diagonalize(A):
    """ diagonalize Hamiltonian """
    Hd, U = np.linalg.eigh(A)
    return Hd, U

# =========================================================
def transform(A, U):
    '''returns U^T.A.U'''
    return np.dot(np.array(U).T, np.dot(A, U))

# =========================================================
def getQMout(QMin,SH2LVC):
    '''Calculates the MCH Hamiltonian, SOC matrix ,overlap matrix, gradients, DM'''
  
    QMout={}
  
    nmult = len(QMin['states'])
    r3N = range(3*QMin['natom'])
  
    # Diagonalize Hamiltonian and expand to the full ms-basis
    U  = [[ 0. for i in range(QMin['nmstates']) ] for j in range(QMin['nmstates']) ]
    Hd = [[ 0. for i in range(QMin['nmstates']) ] for j in range(QMin['nmstates']) ]
    dHfull = [[[ 0. for i in range(QMin['nmstates']) ] for j in range(QMin['nmstates']) ] for iQ in r3N]
    offs = 0
    for imult in range(nmult):
        dim = QMin['states'][imult]
        if not dim == 0:
            Hdtmp,Utmp=diagonalize(SH2LVC['H'][imult])
            for ms in range(imult+1):
                for i in range(dim):
                    Hd[i+offs][i+offs] = Hdtmp[i]
                    U[i+offs][offs:offs+dim] = Utmp[i]
                    for iQ in r3N:
                        dHfull[iQ][i+offs][offs:offs+dim] = SH2LVC['dH'][iQ][imult][i]
                offs += dim
#  print "QMout1: CPU time: % .3f s, wall time: %.3f s"%(time.clock() - tc, time.time() - tt)
  
    # Transform the gradients to the MCH basis
    dE = [[[0. for iQ in range(3*QMin['natom'])] for istate in range(QMin['nmstates'])] for jstate in range(QMin['nmstates'])]
    for iQ in r3N:
        dEmat = transform(dHfull[iQ], U)
        for istate in range(QMin['nmstates']):
            for jstate in range(QMin['nmstates']):
                dE[istate][jstate][iQ] = dEmat[istate][jstate]
  #  print "QMout2: CPU time: % .3f s, wall time: %.3f s"%(time.clock() - tc, time.time() - tt)
  
    # Convert the gradient to Cartesian coordinates
    #   -> It would be more efficent to do this only for unique Ms values
    VOdE = [0. for i in r3N]
    grad = []
    for istate in range(QMin['nmstates']):
        OdE = [0. for iQ in r3N]
        for iQ in r3N:
            if abs(SH2LVC['Om'][iQ]) > 1.e-8:
                OdE[iQ] = dE[istate][istate][iQ] * SH2LVC['Om'][iQ]**0.5
        VOdE = numpy.dot(SH2LVC['V'], OdE)
  
        grad.append([])
        for iat in range(QMin['natom']):
            grad[-1].append([VOdE[3*iat] * SH2LVC['Ms'][3*iat], VOdE[3*iat+1] * SH2LVC['Ms'][3*iat+1], VOdE[3*iat+2] * SH2LVC['Ms'][3*iat+2]])
  #  print "QMout3: CPU time: % .3f s, wall time: %.3f s"%(time.clock() - tc, time.time() - tt)
  
    if 'nacdr' in QMin:
        nonac = [ [0., 0., 0.] for iat in range(QMin['natom'])]
        QMout['nacdr'] = [[ nonac for istate in range(QMin['nmstates'])] for jstate in range(QMin['nmstates'])]
        istate =- 1
        for imult,ist,ims in itnmstates(QMin['states']):
            istate += 1
  
            jstate =- 1
            for jmult,jst,jms in itnmstates(QMin['states']):
                jstate += 1
  
                if imult == jmult and ims == jms and istate < jstate:
                    OdE = [0. for iQ in r3N]
                    for iQ in r3N:
                        if abs(SH2LVC['Om'][iQ]) > 1.e-8:
                            OdE[iQ] = dE[istate][jstate][iQ] * SH2LVC['Om'][iQ]**0.5
                    VOdE = numpy.dot(SH2LVC['V'], OdE)
  
                    deriv = []
                    for iat in range(QMin['natom']):
                        deriv.append([VOdE[3*iat] * SH2LVC['Ms'][3*iat], VOdE[3*iat+1] * SH2LVC['Ms'][3*iat+1], VOdE[3*iat+2] * SH2LVC['Ms'][3*iat+2]])
  
                    Einv = (Hd[jstate][jstate] - Hd[istate][istate]) ** (-1.)
  
                    QMout['nacdr'][istate][jstate] = [ [ c * Einv for c in d] for d in deriv]
                    QMout['nacdr'][jstate][istate] = [ [-c * Einv for c in d] for d in deriv]
  
    # transform dipole matrices
    dipole=[]
    for idir in range(3):
        Dmatrix=transform(SH2LVC['dipole'][idir+1],U).tolist()
        dipole.append(Dmatrix)
  
    # get overlap matrix
    if 'overlap' in QMin:
        Uoldfile=os.path.join(QMin['savedir'],'Uold.out')
        if 'init' in QMin:
            overlap = [ [ float(i==j) for i in range(QMin['nmstates']) ] for j in range(QMin['nmstates']) ]
        else:
            Uold = [[float(v) for v in line.split()] for line in open(Uoldfile, 'r').readlines()]
            overlap = numpy.dot(numpy.array(Uold).T,U)
        QMout['overlap']=overlap
  
  
    Ufile=os.path.join(QMin['savedir'],'U.out')
    f = open(Ufile, 'w')
    for line in U:
        for c in line:
            f.write(str(c) + ' ')
        f.write('\n')
    f.close()
  
    # transform SOC matrix
    SO=transform(SH2LVC['soc'],U)
    for i in range(QMin['nmstates']):
        SO[i][i]=complex(0.,0.)
    Hfull=[ [ Hd[i][j]+SO[i][j] for i in range(QMin['nmstates']) ] for j in range(QMin['nmstates']) ]
  
    # assign QMout elements
    QMout['h']=Hfull
    QMout['dm']=dipole
    QMout['grad']=grad
    #QMout['dmdr']=dmdr
    QMout['runtime']=0.
  
    #pprint.pprint(QMout,width=192)
  
    return QMout


class SHARC_LVC(SHARC_INTERFACE):
    """
    Class for SHARC LVC
    """
    # Name of the interface
    interface = 'LVC'
    # store atom ids
    save_atids = True
    # store atom names
    save_atnames = False
    # accepted units:  0 : Bohr, 1 : Angstrom
    iunit = 0
    # not supported keys
    not_supported = ['nacdt', 'dmdr' ]

    def do_qm_job(self, tasks, Crd):
        """

        Here you should perform all your qm calculations

        depending on the tasks, that were asked

        """
        QMin = self.parseTasks(tasks)
        self.build_lvc_hamiltonian(Crd)
        QMout = getQMout(QMin, self.storage['SH2LVC'])
        return QMout


    def parseTasks(self, tasks):
        """
        these things should be interface dependent

        so write what you love, it covers basically everything
        after savedir information in QMin

        """

        # find init, samestep, restart
        QMin = dict( (key, value) for key, value in self.QMin.items() )
        QMin['natom'] = self.NAtoms

        key_tasks = tasks['tasks'].lower().split()


        if any( [self.not_supported in key_tasks ] ):
            print("not supported keys: ", self.not_supported)
            sys.exit(16)

        for key in key_tasks:
            QMin[key] = []

        for key in self.states:
            QMin[key] = self.states[key]


        if 'init' in QMin:
            checkscratch(QMin['savedir'])
        if not 'init' in QMin and not 'samestep' in QMin and not 'restart' in QMin:
            fromfile=os.path.join(QMin['savedir'],'U.out')
            if not os.path.isfile(fromfile):
                print('ERROR: savedir does not contain U.out! Maybe you need to add "init" to QM.in.')
                sys.exit(1)
            tofile=os.path.join(QMin['savedir'],'Uold.out')
            shutil.copy(fromfile,tofile)

        for key in ['grad', 'nacdr']:
            if tasks[key].strip() != "":
                QMin[key] = []

        QMin['pwd'] = os.getcwd()
        return QMin

    def compute_displacement(self, Crd):
        """
        Compute the difference between the Vectors
        """

        disp =  [ (Crd[i][j] - self.storage['SH2LVC']['CEq'][i*3+j])
                    for i in range(self.NAtoms) for j in range(3) ]
        return disp


    def readParameter(self, fname, *args, **kwargs):
        """
        read basic parameter files for the calculation

        here was before called read_SH2LVC

        add SH2LVC to storage
        """
        # reads LVC.template, deletes comments and blank lines

        #print fname

        try:
          fname = os.path.abspath(fname)
          f=open(fname.strip())
        except IOError:
          try:
            fname = os.path.abspath('SH2LVC.inp')
            f=open(fname)
          except IOError:
            print('Input file "LVC.template" not found.')
            sys.exit(1)
        sh2lvc=f.readlines()
        f.close()

        self.storage['V0'] = sh2lvc[0].strip()
        SH2LVC = self.read_V0(fname=os.path.join(os.path.dirname(fname),self.storage['V0']))
        # read EPSILON
        tmp = find_lines(1, 'epsilon',sh2lvc)
        eps = []
        if not tmp==[]:
            neps = int(tmp[0])
            tmp = find_lines(neps+1, 'epsilon', sh2lvc)
            for line in tmp[1:]:
                words = line.split()
                eps.append((int(words[0])-1, int(words[1])-1, float(words[-1])))
        SH2LVC['epsilon'] = eps
        # read KAPPA
        tmp = find_lines(1, 'kappa', sh2lvc)
        kappa = []
        if not tmp==[]:
            nkappa = int(tmp[0])
            tmp = find_lines(nkappa+1, 'kappa', sh2lvc)
            for line in tmp[1:]:
                words = line.split()
                kappa.append((int(words[0])-1, int(words[1])-1, int(words[2])-1, float(words[-1])))
        SH2LVC['kappa'] = kappa
        # read LAMBDA
        tmp = find_lines(1, 'lambda', sh2lvc)
        lam = []
        if not tmp==[]:
            nlam = int(tmp[0])
            tmp = find_lines(nlam+1, 'lambda', sh2lvc)
            for line in tmp[1:]:
                words = line.split()
                lam.append((int(words[0])-1, int(words[1])-1, int(words[2])-1, int(words[3])-1, float(words[-1])))
        SH2LVC['lambda'] = lam
        # read DIPOLE
        nmstates = self.states['nmstates']
        SH2LVC['dipole'] = {}
        SH2LVC['dipole'][1]= read_LVC_mat(nmstates, 'DMX', sh2lvc)
        SH2LVC['dipole'][2]= read_LVC_mat(nmstates, 'DMY', sh2lvc)
        SH2LVC['dipole'][3]= read_LVC_mat(nmstates, 'DMZ', sh2lvc)
        # obtain the SOC matrix
        SH2LVC['soc'] = read_LVC_mat(nmstates, 'SOC', sh2lvc)
        #  save SH2LVC
        self.storage['SH2LVC'] = SH2LVC
        return

    def build_lvc_hamiltonian(self, Crd):
        """
        does everything that was before done by **read_SH2LVC**
        except fileio  etc. everything there should be already
        done by readParameter?
        """
        states = self.states['states']
        nmult = len(states)
        r3N = range(3*self.NAtoms)
        # get access to SH2LVC
        SH2LVC = self.storage['SH2LVC']
        Om = SH2LVC['Om']
        # compute displacement compared to reference structure
        disp = self.compute_displacement(Crd)
        # Transform the coordinates to dimensionless mass-weighted normal modes
        MR = [SH2LVC['Ms'][i] * disp[i] for i in r3N]
        MRV = [ sum(MR[j] * SH2LVC['V'][j][i] for j in r3N) for i in r3N]
        Q =  [MRV[i] * Om[i]**0.5 for i in r3N]
        # Compute the ground state potential and gradient
        V0 = sum(0.5 * Om[i] * Q[i]*Q[i] for i in r3N)
        HMCH =  [[[
                    V0 if istate == jstate else 0.
                       for istate in range(states[imult])]
                       for jstate in range(states[imult])]
                       for imult in range(nmult)]
        #l
        dHMCH = [[[[
          Om[i] * Q[i] if istate == jstate else 0.
                       for istate in range(states[imult])]
                       for jstate in range(states[imult])]
                       for imult in range(nmult)]
                       for i in r3N]
        # Add the vertical energies (epsilon)
        # Enter in separate lines as:
        # <n_epsilon>
        # <mult> <state> <epsilon>
        # <mult> <state> <epsilon>
        for imult, istate, val in SH2LVC['epsilon']:
            HMCH[imult][istate][istate] += val
        # Add the intrastate LVC constants (kappa)
        # Enter in separate lines as:
        # <n_kappa>
        # <mult> <state> <mode> <kappa>
        # <mult> <state> <mode> <kappa>
        for imult, istate, i, val  in SH2LVC['kappa']:
            HMCH[imult][istate][istate]  += val * Q[i]
            dHMCH[i][imult][istate][istate] += val
        # Add the interstate LVC constants (lambda)
        # Enter in separate lines as:
        # <n_lambda>
        # <mult> <state1> <state2> <mode> <lambda>
        # <mult> <state1> <state2> <mode> <lambda>
        for imult, istate, jstate, i, val  in SH2LVC['lambda']:
            HMCH[imult][istate][jstate]  += val * Q[i]
            HMCH[imult][jstate][istate]  += val * Q[i]
            dHMCH[i][imult][istate][jstate] += val
            dHMCH[i][imult][jstate][istate] += val
        # save HMCH and dHMCH
        SH2LVC['H'] = HMCH
        SH2LVC['dH'] = dHMCH

    def read_V0(self, fname='V0.txt'):
        """"
        change original read_V0 routine

        to only parse the parameters and equilibrium geometry
        and save everything in memory.

        Reads information about the ground-state potential
        from V0.txt.
        Returns the displacement vector.
        """
        SH2LVC = {}
        try:
            f=open(fname)
        except IOError:
            print('Input file "%s" not found.'%fname)
            sys.exit(1)
        v0=f.readlines()
        f.close()

        U_TO_AMU = 1.0/self.constants['au2u']

        SH2LVC['CEq'] = [  float(ele)
                for line in find_lines(self.NAtoms, 'Geometry',v0)
                for ele in line.split()[2:5] ]

        SH2LVC['Ms'] =  [ (float(line.split()[5])*U_TO_AMU)**.5
                for line in find_lines(self.NAtoms, 'Geometry',v0)
                for i in range(3) ]

        # Frequencies (a.u.)
        tmp = find_lines(1, 'Frequencies',v0)
        if tmp==[]:
            print('No Frequencies defined in %s!'%fname)
            sys.exit(24)
        SH2LVC['Om'] = list(map(float,tmp[0].split()))
        # Normal modes in mass-weighted coordinates
        tmp = find_lines(len(SH2LVC['Om']), 'Mass-weighted normal modes', v0)
        if tmp==[]:
            print('No normal modes given in %s!'%fname)
            sys.exit(24)
        SH2LVC['V']  = [list(map(float,line.split())) for line in tmp] # transformation matrix
        return SH2LVC


def getCommandoLine():
    """
        Get Commando line option with argpase
        
    """

    import argparse

    parser = argparse.ArgumentParser("Perform SHARC LVC calculations")
    parser.add_argument("input", metavar="FILE",type=str,
                        default="input", nargs='?',
                        help="input file")
    parser.add_argument("param", metavar="FILE",type=str,
                        default="QM/LVC.template", nargs='?',
                        help="param file, LVC.template")

    args = parser.parse_args()

    return args.input, args.param


def main():
    """
        Main Function if program is called as standalone 

    """

    inp_file, param = getCommandoLine()
    # init SHARC_LVC class
    lvc = SHARC_LVC()
    # run sharc dynamics
    lvc.run_sharc(inp_file, param)

if __name__ == "__main__":
    main()



