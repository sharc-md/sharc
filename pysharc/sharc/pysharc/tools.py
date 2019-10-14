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


import sys
from . import fileio

def check_version(major, minor):
    if major != sys.version_info[0]:
        return False
    if minor < sys.version_info[1]:
        return False

    return True

if check_version(3,0):
    def lst2dct(lst):
        return { i : value for i, value in enumerate(lst) }
else:
    def lst2dct(lst):
        return dict((i, value) for i, value in enumerate(lst) )

def writeQMout(QMin, QMout, QMoutfile='QM.out'):
    """

    """
    string=''
    if 'h' in QMout:
        string+=writeQMoutsoc(QMin,QMout)
    if 'dm' in QMout:
        string+=writeQMoutdm(QMin,QMout)
    if 'grad' in QMout:
        string+=writeQMoutgrad(QMin,QMout)
    if 'overlap' in QMout:
        string+=writeQMoutnacsmat(QMin,QMout)
    if 'socdr' in QMout:
        string+=writeQMoutsocdr(QMin,QMout)
    if 'dmdr' in QMout:
        string+=writeQMoutdmdr(QMin,QMout)
    if 'ion' in QMout:
        string+=writeQMoutprop(QMin,QMout)
    if 'phases' in QMout:
        string+=writeQmoutPhases(QMin,QMout)
    fileio.writeOutput(QMoutfile, string)

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
            string+='%s %s ' % (eformat(QMout['h'][i][j].real,9,3),eformat(QMout['h'][i][j].imag,9,3))
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
                string+='%s %s ' % (eformat(QMout['dm'][xyz][i][j].real,9,3),eformat(QMout['dm'][xyz][i][j].imag,9,3))
            string+='\n'
        #string+='\n'
    string+='\n'
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
        string+='%i %i ! %i %i %i\n' % (natom,3,imult,istate,ims)
        for atom in range(natom):
            for xyz in range(3):
                string+='%s ' % (eformat(QMout['grad'][i][atom][xyz],9,3))
            string+='\n'
        #string+='\n'
        i+=1
    string+='\n'
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
            string+='%s %s ' % (eformat(QMout['overlap'][j][i].real,9,3),eformat(QMout['overlap'][j][i].imag,9,3))
        string+='\n'
    string+='\n'
    return string

# ======================================================================= #
def writeQMoutdmdr(QMin,QMout):

  states=QMin['states']
  nmstates=QMin['nmstates']
  natom=QMin['natom']
  string=''
  string+='! %i Dipole moment derivatives (%ix%ix3x%ix3, real)\n' % (12,nmstates,nmstates,natom)
  i=0
  for imult,istate,ims in itnmstates(states):
    j=0
    for jmult,jstate,jms in itnmstates(states):
      for ipol in range(3):
        string+='%i %i ! m1 %i s1 %i ms1 %i   m2 %i s2 %i ms2 %i   pol %i\n' % (natom,3,imult,istate,ims,jmult,jstate,jms,ipol)
        for atom in range(natom):
          for xyz in range(3):
            string+='%s ' % (eformat(QMout['dmdr'][ipol][i][j][atom][xyz],12,3))
          string+='\n'
        string+=''
      j+=1
    i+=1
  string+='\n'
  return string

# ======================================================================= #
def writeQMoutsocdr(QMin,QMout):

  states=QMin['states']
  nmstates=QMin['nmstates']
  natom=QMin['natom']
  string=''
  string+='! %i Spin-Orbit coupling derivatives (%ix%ix3x%ix3, complex)\n' % (13,nmstates,nmstates,natom)
  i=0
  for imult,istate,ims in itnmstates(states):
    j=0
    for jmult,jstate,jms in itnmstates(states):
        string+='%i %i ! m1 %i s1 %i ms1 %i   m2 %i s2 %i ms2 %i\n' % (natom,3,imult,istate,ims,jmult,jstate,jms)
        for atom in range(natom):
            for xyz in range(3):
                string+='%s %s ' % (eformat(QMout['socdr'][i][j][atom][xyz].real,12,3),eformat(QMout['socdr'][i][j][atom][xyz].imag,12,3))
        string+='\n'
        string+=''
        j+=1
    i+=1
  string+='\n'
  return string

# ======================================================================= #
def writeQMoutprop(QMin,QMout):

    nmstates=QMin['nmstates']

    # print property matrix (flag 11) for backwards compatibility
    string=''
    string+='! %i Property Matrix (%ix%i, complex)\n' % (11,nmstates,nmstates)
    string+='%i %i\n' % (nmstates,nmstates)
    for i in range(nmstates):
        for j in range(nmstates):
            string+='%s %s ' % (eformat(QMout['prop'][i][j].real,12,3),eformat(QMout['prop'][i][j].imag,12,3))
        string+='\n'
    string+='\n'

    # print property matrices (flag 20) in new format
    string+='! %i Property Matrices\n' % (20)
    string+='%i    ! number of property matrices\n' % (1)

    string+='! Property Matrix Labels (%i strings)\n' % (1)
    string+='Dyson norms\n'

    string+='! Property Matrices (%ix%ix%i, complex)\n' % (1,nmstates,nmstates)
    string+='%i %i   ! Dyson norms\n' % (nmstates,nmstates)
    for i in range(nmstates):
        for j in range(nmstates):
            string+='%s %s ' % (eformat(QMout['prop'][i][j].real,12,3),eformat(QMout['prop'][i][j].imag,12,3))
        string+='\n'
    string+='\n'

    return string

# ======================================================================= #
def writeQMoutTHEODORE(QMin,QMout):

    nmstates=QMin['nmstates']
    nprop=QMin['template']['theodore_n']
    if QMin['template']['qmmm']:
        nprop+=7
    if nprop<=0:
        return '\n'

    string=''

    string+='! %i Property Vectors\n' % (21)
    string+='%i    ! number of property vectors\n' % (nprop)

    string+='! Property Vector Labels (%i strings)\n' % (nprop)
    descriptors=[]
    if 'theodore' in QMin:
        for i in QMin['template']['theodore_prop']:
            descriptors.append('%s' % i)
            string+=descriptors[-1]+'\n'
        for i in range(len(QMin['template']['theodore_fragment'])):
            for j in range(len(QMin['template']['theodore_fragment'])):
                descriptors.append('Om_{%i,%i}' % (i+1,j+1))
                string+=descriptors[-1]+'\n'
    if QMin['template']['qmmm']:
        for label in QMout['qmmm_energies']:
            descriptors.append(label)
            string+=label+'\n'

    string+='! Property Vectors (%ix%i, real)\n' % (nprop,nmstates)
    if 'theodore' in QMin:
        for i in range(QMin['template']['theodore_n']):
            string+='! TheoDORE descriptor %i (%s)\n' % (i+1,descriptors[i])
            for j in range(nmstates):
                string+='%s\n' % (eformat(QMout['theodore'][j][i].real,12,3))
    if QMin['template']['qmmm']:
        for label in QMout['qmmm_energies']:
            string+='! QM/MM energy contribution (%s)\n' % (label)
            for j in range(nmstates):
                string+='%s\n' % (eformat(QMout['qmmm_energies'][label],12,3))
    string+='\n'

    return string

# ======================================================================= #
def writeQmoutPhases(QMin,QMout):

    string='! 7 Phases\n%i ! for all nmstates\n' % (QMin['nmstates'])
    for i in range(QMin['nmstates']):
        string+='%s %s\n' % (eformat(QMout['phases'][i].real,9,3),eformat(QMout['phases'][i].imag,9,3))
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


# =============================================================================================== #

