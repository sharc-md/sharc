#!/usr/bin/env python2

#******************************************
#
#    SHARC Program Suite
#
#    Copyright (c) 2019 University of Vienna
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
import time
# relative packages
from .. import sharc
from .constants import IAn2AName
from .constants import IAn2AName
from .tools import writeQMout, lst2dct
from . import fileio 


class SHARC_INTERFACE(object):
    """
    Basic SHARC Interface class
    """
    # Name of the interface
    interface = 'None'
    # basic dict for permanently storing info
    storage = {}
    # qmout object, will be set in run_sharc
    QMout = None
    # store atom ids
    save_atids = False
    # store atom names
    save_atnames = True
    # accepted units:  0 : Bohr, 1 : Angstrom
    iunit = 0
    # not supported keys
    not_supported = []
    # if True the QMin object instead of the Crd is given to do_qm_job!
    use_qmin = False
    # if False It is do_qm_job has to return a dct, else it should return None
    set_qmout = False
    #
    iskip = 1   # has no effect anymore


    def __init__(self, *args, **kwargs):
        """
        Init your interface, best is you

        set parameter files etc. for read

        """
        pass

    def readParameter(self, *args, **kwargs):
        """
        read basic parameter files for the calculation

        everything that should be read only once during
        the calculation

        """
        pass

    def do_qm_job(self, tasks, Crd):
        """

        Here you should perform all your qm calculations

        depending on the tasks, that were asked

        needs to return a QMout object if set_qmout is False
        else you are responsible to set QMout in do_qm_job
        - only for advanced users!

        QMout is a dict with corresponding entries

        'h'     : Hamiltonian, lst[NStates][NStates] of float, complex
        'ovlap' : Overlap, lst[NStates][NStates] of float, complex
        'dm'    : Dipole Matrix, lst[3][NStates][NStates] of floats
        'grad'  : Dict, int IState : lst Grad_IState,
                        Grad_IStat, lst[NAtoms][3] of floats
        'nacdr' : Dict,   int IState : dct Dict_2,
                  Dict_2, int JState : lst NACDR_I_J,
                  NACDR_I_J, lst[NAtoms][3] of floats

        alternativly, the old sharc format can be used,

        which means that grad is given as a lst[NSTates][NAtoms][3]
        and nacdr is given als lst[NStates][NStates][NAtoms][3]

        """
        pass


    def do_single_job(self, QMin=None, tasks=None, Crd=None, basic_info=None, IAn=None, AtNames=None, *args, **kwargs):
        """ Run a single point charc Calculations starting from QMin """

        # setup_calculation
        self.initial_setup(**kwargs)
        if IAn is not None:
            self.IAn = IAn
        if AtNames is not None:
            self.AtNames = AtNames
        if QMin is not None:
            # not implimented yet!
            tasks, Crd, basic_info = self.sharc_readQMin(QMin)
        # set static objects
        self.readParameter(*args, **kwargs)

        QMout = self.do_qm_job(tasks, Crd)

        return QMout
    def final_print(self):
        """
        Called before sharc_finalize
        """
        pass

    def initial_setup(self, **kwargs):
        """
        Is called before run_sharc was setup
        """
        pass

    def crash_function(self):
        """
        This function is called in case do_qm_calculation crashes!
        """
        self.sharc_writeQMin()
        pass

    def sharc_initial_setup(self):
        """
        This function is used to handle all setup tasks 
        for the interface class before the dynamics starts!
        """
        if self.use_qmin is True:
            self._sharc_QMin = sharc.QMin(self.NAtoms)
    

    def sharc_readQMin(self, QMin):
        """
        Reads standard QMin file, returns tasks, Crd, and basic_info in PySHARC format

        sets also self.IAn and self.AtNames

        """
        pass


    def sharc_writeQMout(self, QMout, QMoutfile='QM.out'):
        """
        writes QMout file, based on the QMout dct!
        """
        QMin = {
            'natom' : self.NAtoms,
            'states' : self.states['states'],
            'nstates' : self.states['nstates'],
            'nmstates' : self.states['nmstates'],
            }

        writeQMout(QMin, QMout, QMoutfile)


    def sharc_writeQMin(self, QMout='QM.in'):
        """ Writes a QMin file, based on the current tasks, etc.  """
        tasks, Crd = self.sharc_get_sharc_tasks(1, True)
        basic_info = sharc.get_basic_info()

        if self.use_qmin is False:
            txt = "%d\nQM.in created with SHARC_INTERFACE\n" % self.NAtoms
            if self.AtNames is None:
                for i in range(self.NAtoms):
                    txt += "%s   %12.8f %12.8f %12.8f \n" % ( IAn2AName[self.IAn[i]], Crd[i][0], Crd[i][1], Crd[i][2] )
            else:
                for i in range(self.NAtoms):
                    txt += "%s   %12.8f %12.8f %12.8f \n" % ( self.AtNames[i], Crd[i][0], Crd[i][1], Crd[i][2] )
        else:
            txt=""
        
        txt += "# Basic Info \n"
        if self.iunit == 0:
            txt += 'Unit Bohr\n'
        else:
            txt += 'Unit Angstrom\n'


        txt += 'states %s\n' % basic_info['states']
        txt += 'savedir %s\n' % basic_info['savedir']

        for task in tasks['tasks'].split():
            txt += "%s\n" % task
        if tasks['grad'].strip() != "":
            txt += "grad %s\n" % tasks['grad']
        if tasks['nacdr'].strip() != "":
            txt += "%s\n" % tasks['nacdr']

        fileio.writeOutput(QMout, txt)


    def sharc_set_basic(self, basic_info):

        self.states = self.getStates(basic_info['states'])
        if self.save_atids is True:
            self.IAn = basic_info['IAn']    # basic
        else:
            self.IAn = None

        if self.save_atnames is True:
            self.AtNames = [ IAn2AName[ian] for ian in  basic_info['IAn'] ]
        else:
            self.AtNames = None

        self.savedir = basic_info['savedir']
        self.NAtoms = basic_info['NAtoms']
        self.nsteps = basic_info['NSteps']
        self.istep  = basic_info['istep']
        self.constants = sharc.get_constants()
        self.QMin = { 'savedir' : basic_info['savedir'] }

    def sharc_get_sharc_tasks(self, icall, getCrd):
        tasks = sharc.get_all_tasks(icall)
        if self.use_qmin is True:
            return tasks, self._sharc_QMin

        if getCrd is True:
            Crd = sharc.get_crd(self.iunit)
        else:
            Crd = None
        return tasks, Crd


    def sharc_set_QMout(self, QMout, icall):
        """ set QMout """

        if self.set_qmout is True:
            # assume that QMout is set already! only for advanced users!
            return

        # set hamiltonian, dm only in first call
        if icall == 1:
            if 'h' in QMout:
                self.QMout.set_hamiltonian(QMout['h'])
            if 'dm' in QMout:
                self.QMout.set_dipolemoment(QMout['dm'])

        if 'overlap' in QMout:
            if type(QMout['overlap']) != type([]):
                # assumes type is numpy array
                QMout['overlap'] = [ list(ele) for ele in QMout['overlap'] ]
            self.QMout.set_overlap(QMout['overlap'])

        if 'grad' in QMout:
            if type(QMout['grad']) == type([]):
                self.QMout.set_gradient(lst2dct(QMout['grad']), icall)
            else:
                if QMout['grad'] is None:
                    QMout['grad'] = {}
                self.QMout.set_gradient(QMout['grad'], icall)
        if 'nacdr' in QMout:
            if type(QMout['nacdr']) == type([]):
                nacdr = {}
                for i, ele in enumerate(QMout['nacdr']):
                    nacdr[i] = lst2dct(ele)
                self.QMout.set_nacdr(nacdr, icall)

            else:
                self.QMout.set_nacdr(QMout['nacdr'], icall)

        return

    def sharc_qm_failure_handle(self, tasks, Crd):
        """In case the qm run crashes for some reason """

        try:
            return self.do_qm_job(tasks, Crd)
        except:
            self.crash_function()
            sharc.finalize_sharc()
            print("Unexpected error:", sys.exc_info()[0])
            sys.exit(101)

    def sharc_do_qm_calculation(self):
        """Replacement of the do_qm_calculation subroutine in main.f90 of source sharc """

        icall = 1
        tasks, Crd = self.sharc_get_sharc_tasks(icall, True)
        # call do_qm_job
        self.sharc_set_QMout(self.sharc_qm_failure_handle(tasks, Crd), icall)
        isecond = sharc.set_qmout(self.QMout, icall)
        if isecond == 1:
            icall = 2
            tasks, _ = self.sharc_get_sharc_tasks(icall, False)
            self.sharc_set_QMout(self.sharc_qm_failure_handle(tasks, Crd), icall)
            sharc.set_qmout(self.QMout, icall)
        return Crd


    def sharc_redo_qm_gradients(self, Crd):
        """ Replacement of the redo_qm_gradients subroutine in main.f90 of source sharc """

        icall = 3
        tasks, _ = self.sharc_get_sharc_tasks(icall, False)
        self.sharc_set_QMout(self.do_qm_job(tasks, Crd), icall)
        sharc.set_qmout(self.QMout, icall)



    def run_sharc(self, inp_file, *args, **kwargs):
        """Main SHARC Routine

        call do_qm_job for what ever job there is
        """
        # setup_calculation
        self.initial_setup(**kwargs)                    # nothing (dummy)
        # setup sharc
        IRestart = sharc.setup_sharc(inp_file)          # in C code, is declared as function
        # set basic info
        self.sharc_set_basic(sharc.get_basic_info())    # natoms, atoms, states, ...
        # setup the calculation stuff
        self.sharc_initial_setup()                      # nothing if use_qmin=False
        # set static objects
        self.readParameter(*args, **kwargs)             # nothing (dummy)
        # set QMout as sharc.QMout class
        self.QMout = sharc.QMout(self.interface, self.NAtoms, self.states['nmstates'])
        # if not Restart, do first QM calculation!
        if IRestart == 0:
            # do initial QM job
            sharc.initial_qm_pre()
            self.sharc_do_qm_calculation()
            sharc.initial_qm_post()
            # finish setup if not restart
            sharc.initial_step(IRestart)
        #do main sharc loop
        for istep in range(self.istep+1, self.nsteps+1):
            sharc.verlet_xstep(istep)
            # call do_qm_job
            Crd = self.sharc_do_qm_calculation()
            #
            IRedo = sharc.verlet_vstep()
            # do missing gradient calculations
            if IRedo == 1:
                self.sharc_redo_qm_gradients(Crd)
            # Verlet last step
            iexit = sharc.verlet_finalize(self.iskip)
            if iexit == 1:
                break
        self.final_print()
        # finalize sharc
        sharc.finalize_sharc()

    @classmethod
    def change_iskip(cls, value):
        assert type(value) == int
        assert value > 0

        cls.iskip = value

    def getStates(self, txt):
        """ convert txt string in sharc, nstates data!  """

        states={}
        states['states'] = list(map(int,txt.split()))
        #states['states'] = [ int(i) for i in txt.split()]
        states['nstates'] = 0
        states['nmstates'] = 0
        for mult, i in enumerate(states['states']):
            states['nstates']  += i
            states['nmstates'] += (mult+1)*i
        states['nmult'] = 0
        statemap = {}
        i = 1
        for imult, nstates in enumerate(states['states']):
            if nstates == 0:
                continue
            states['nmult'] += 1
            for ims in range(imult+1):
                ms = ims-imult/2.0
                for istate in range(nstates):
                    statemap[i] = [imult+1, istate+1, ms]
                    i += 1
        states['statemap'] = statemap
        return states




