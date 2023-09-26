"""
version 1.1
author: Felix Plasser, Hans Georg Gallmetzer
usage: Package for reading a molden input vibration file.
    This input is used for normal mode analysis by nma.py.
    Also a multiple xyz vibration file that can be read by jmol can be printed out.
changes: xrange() changed to range() for compatability with python 3.11.4
"""

import os
import numpy
import file_handler

def convert_molden2xyz(path='.', name='molden.input'):
    """
    Sub to convert a molden file to an xyz file.
    """
    vm2x = vib_molden()
    vm2x.read_molden_file(os.path.join(path, name))
    vm2x.print_xyz_file(os.path.join(path, name.rpartition('.')[0] + '.xyz'))

class vib_molden:
    """
    Main class that reads in a molden input file.
    """
    def read_molden_file(self, file_path):
        mfile = open(file_path, 'r')
        Atoms = False
        FREQ = False
        FRNORMCOORD = False

        self.atoms = []
        self.freqs = []
        self.vibs = []
        actvib = -1
        for line in mfile:
            # what section are we in
            if '[' in line or '--' in  line:
                Atoms = False
                FREQ = False
                FRNORMCOORD = False

            if '[Atoms]' in line:
                Atoms = True
            elif '[FREQ]' in line:
                FREQ = True
            elif '[FR-NORM-COORD]' in line:
                FRNORMCOORD = True
            # extract the information in that section
            elif Atoms:
                words = line.split()
                self.atoms += [atom(words[0],words[3:6])]
            elif FREQ:
                self.freqs += [eval(line)]
            elif FRNORMCOORD:
                if 'vibration' in line or 'Vibration' in line:
                    actvib += 1
                    self.vibs += [vibration(self.freqs[actvib])]
                else:
                    self.vibs[actvib].add_vector(line.split())

    def ret_vib_matrix(self):
        """
        Return a normalised matrix that contains the coordinates for all the vibrations.
        The modes are in the lines of this matrix!
        """
        return numpy.array([vib.ret_joined_vector() for vib in self.vibs])

    def ret_freqs(self):
        """
        Return the frequencies.
        """
        return self.freqs

    def ret_eff_masses(self, mol_calc, mass_wt_pw=1):
        """
        Returns a list with the effective masses for all the modes.
        <mol_calc> is a struc_linalg.mol_calc instance with the mass information.
        """
        return [vib.ret_eff_mass(mol_calc, mass_wt_pw=mass_wt_pw) for vib in self.vibs]

    def ret_nma_header(self):
        """
        Returns a list that can be used for a header for a normal mode analysis table.
        """
        nr_list = ['1']
        freq_list = ['Time'] # frequency
        T_list = [''] # period
        for nr, freq in enumerate(self.freqs):
            try:
                T = 1/(freq * 2.9979E-5)
            except ZeroDivisionError:
                T = 0

            nr_list += [str(nr+2)]
            freq_list += [str(freq)[:6]]
            T_list += [str(T)[:6]]

        return [nr_list, freq_list, T_list]

    def print_xyz_file(self, file_path):
        """
        Creates an xyz file with the vibrational information.
        """
        out_str = ''
        for vib in self.vibs:
            out_str += str(len(self.atoms)) + '\n'
            out_str += 'Frequency ' + str(vib.frequency) + '\n'

            tmaker = file_handler.table_maker([5] + 6 * [20])
            for at_ind, atom in enumerate(self.atoms):
                tmaker.write_line([' ' + atom.name] + atom.pos + vib.vector_list[at_ind])
            out_str += tmaker.return_table()

        xyzfile = open(file_path, 'w')
        xyzfile.write(out_str)
        xyzfile.close()

class atom:
    def __init__(self, name, pos, length_factor=.529177): # units are changed from Bohr into Angstrom
        " Position in A "
        self.name = name.capitalize()
        self.pos = [eval(coor) *  length_factor for coor in pos]

class vibration:
    def __init__(self, frequency):
        self.vector_list = []
        self.frequency = frequency

    def add_vector(self, vector, length_factor=1.): # in this case the units are not changed to keep orthonormality
        self.vector_list += [[eval(coor) * length_factor for coor in vector]]

    def ret_joined_vector(self):
        """
        Returns a joined 3N-dimensional numpy vector for the vibration.
        """
        t_list = []
        for vector in self.vector_list:
            t_list += vector

        return numpy.array(t_list)

    def ret_eff_mass(self, mol_calc, mass_wt_pw=1):
        """
        Return the effective mass of the vibration according to structure <struc>.
        """
        M = mol_calc.ret_mass_matrix(power=mass_wt_pw)

        vec = self.ret_joined_vector()
        return numpy.dot(numpy.dot(vec, M), vec) / numpy.dot(vec,vec)

def make_molden_file(struc, freqs, vibs, out_file, title='Essential dynamics', num_at=None):
    """
    Subroutine for making a molden file.
    """
    if num_at == None:
        num_at = struc.ret_num_at()

    out_str = '[Molden Format]\n[Title]\n'+title+'\n'

    if not struc == None:
       out_str += '[Atoms] AU\n'
       out_str += ret_Atoms_table(struc)

    # wavenumbers
    out_str += '[FREQ]\n'
    for freq in freqs:
        out_str += ' ' + str(freq) + '\n'

    if not struc == None:
       out_str += '[FR-COORD]\n'
       out_str += ret_FRCOORD_table(struc)

    # normal modes
    out_str += '[FR-NORM-COORD]\n'
    for ind in range(len(vibs)):
        out_str += ' vibration' + str(ind + 1).rjust(5) + '\n'
#        tblm = file_handler.table_maker([1] + 3*[21])
        for j in range(num_at):
            out_str+= '% 14.8f % 14.8f % 14.8f\n'%(vibs[ind][3*j],vibs[ind][3*j+1],vibs[ind][3*j+2])
#            tblm.write_line([' '] + [coor for coor in vibs[ind][3*j:(3*j+3)]])
#        out_str += tblm.return_table()
        #print out_str

    w_file = open(out_file, 'w')
    w_file.write(out_str)
    w_file.close()

b_in_a = 0.529178 # conversion factor between Bohr and Angstrom

def ret_Atoms_table(struc):
    """
    Create the atoms part in the molden file.
    """
    tblmaker = file_handler.table_maker([6, 4, 3, 21, 21, 21])
    for i in range(struc.ret_num_at()):
        atom = struc.mol.GetAtom(i+1)
        tblmaker.write_line([struc.ret_symbol(i+1),i+1,atom.GetAtomicNum()]+[atom.x()/b_in_a]+[atom.y()/b_in_a]+[atom.z()/b_in_a])

    return tblmaker.return_table()

def ret_FRCOORD_table(struc):
    """
    Create the FRCOORD part in the molden file.
    """
    tblmaker = file_handler.table_maker([4, 21, 21, 21])
    for i in range(struc.ret_num_at()):
        atom = struc.mol.GetAtom(i+1)
        vec = atom.GetVector()
        tblmaker.write_line([struc.ret_symbol(i+1)]+[atom.x()/b_in_a]+[atom.y()/b_in_a]+[atom.z()/b_in_a])

    return tblmaker.return_table()


if __name__ == '__main__':
    import struc_linalg

    mol_calc = struc_linalg.mol_calc(def_file_path='coord', file_type='tmol')
    vmol = vib_molden()
    vmol.read_molden_file('molden.input')
    print(vmol.ret_eff_masses(mol_calc))
