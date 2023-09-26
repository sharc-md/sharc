"""
version 1.1
author: Felix Plasser, Hans Georg Gallmetzer
usage: package to perform linear algebra operations on structures.
changes: xrange() changed to range() for compatability with python 3.11.4
"""

# improvment: use return to pass structures by value (?).

import os, shutil, locale
import numpy
import openbabel
import superposition
import file_handler

class structure:
    """
    Class to manipulate a structure.
    """
    def __init__(self, name=''):
        self.name = name
    
    def read_file(self, file_path, file_type='tmol'):
        """
        Read in the structure from a file.
        """
        self.file_path = file_path
        self.file_type = file_type
        
        obconversion = openbabel.OBConversion()
        obconversion.SetInFormat(file_type)
        self.mol = openbabel.OBMol()
        obconversion.ReadFile(self.mol, file_path)

    def get_mol(self, mol, file_path, file_type='xyz'):
        """
        Read in an openbabel mol that is passed from a different routine.
        Can be used for accessing multiple structure xyz files.
        """
        self.file_path = file_path
        self.file_type = file_type
        self.mol = mol

    def read_file_vector(self, def_file_path, file_type, vector):
        """
        Initialise the structure with a default .mol from a file and a vector.
        """
        self.read_file(def_file_path, file_type) # like this because copying objects doesn't work

        for i in range(self.mol.NumAtoms()):
            atom = self.mol.GetAtom(i+1)
            atom.SetVector(vector[3*i], vector[3*i+1], vector[3*i+2])

    def read_file_3xN_matrix(self, def_file_path, file_type, coor_mat):
        """
        Initialise the structure with a default .mol and a 3xN matrix.
        """
        self.read_file(def_file_path, file_type)

        #print self.mol.NumAtoms(), file_type

        #for i, atom in enumerate(openbabel.OBMolAtomIter(self.mol)): ## this only works with the new version

        self.read_3xN_matrix(coor_mat)

    def read_3xN_matrix(self, coor_mat):
        """
        Read in a 3xN matrix.
        """
        for i in range(self.mol.NumAtoms()):
            atom = self.mol.GetAtom(i+1)
            atom.SetVector(coor_mat[i][0], coor_mat[i][1], coor_mat[i][2])

    def ret_vector(self):
        " All the coordinates in one vector "
        vec_list = []
        for i in range(self.mol.NumAtoms()):
            atom = self.mol.GetAtom(i+1)
            vec_list += [atom.x(), atom.y(), atom.z()]
           
        return numpy.array(vec_list)

    def ret_3xN_matrix(self):
        " Coordinates in a 3 x N matrix "
        mat_list = []

        #for i, atom in enumerate(openbabel.OBMolAtomIter(self.mol)):

        for i in range(self.mol.NumAtoms()):
            atom = self.mol.GetAtom(i+1)
            mat_list += [[atom.x(), atom.y(), atom.z()]]
            
        return numpy.array(mat_list)
        
    def ret_rotated_structure(self, theta, vec, name=''):
        """
        Return the structure rotated by <theta> around <vec>.
        """
        if name == '': name = self.name
        
        coor_mat = self.ret_3xN_matrix()
        
        sup = superposition.superposition()
        sup.set_rotation(theta, vec)
        
        #print coor_mat
        #print sup.ret_rotation_matrix()
        coor_mat = numpy.dot(coor_mat, sup.ret_rotation_matrix().transpose())
        
        ret_struc = structure(name=name)
        ret_struc.read_file_3xN_matrix(self.file_path, self.file_type, coor_mat)
        return ret_struc
        
    def ret_moved_structure(self, add_vec, name=''):
        """
        Move the structure by <add_vec>.
        """
        if name == '': name = self.name
        
        coor_mat = self.ret_3xN_matrix()
        coor_mat += add_vec
        
        ret_struc = structure(name=name)
        ret_struc.read_file_3xN_matrix(self.file_path, self.file_type, coor_mat)
        return ret_struc

    def ret_superimposed_structure(self, struc, mass_wt_pw=1, name='', manual_wts=[]):
        """
        Return this structure superimposed on another structure.
        <self> is fitted onto <struc>.
        <mass_wt_pw>=1 means regular mass weighting.
        <manual_wts> is a nested list [[ind1,wt1],...] with weights that can be put in by hand. Indeces start with 1.
        """
        if name == '': name = self.name
        
        coor_mat = self.ret_3xN_matrix()
        
        mass_vect = self.ret_mass_vector(power=mass_wt_pw)
        for ind, weight in manual_wts:
            mass_vect[ind-1] = weight
        
        sup = superposition.superposition()
        sup.superimpose(ref_points=struc.ret_3xN_matrix(), mv_points=coor_mat, weights=mass_vect)
       

        coor_mat = coor_mat - sup.ret_mv_av()
        coor_mat = numpy.dot(coor_mat, sup.ret_rotation_matrix().transpose())
        coor_mat = coor_mat + sup.ret_ref_av()

        #print coor_mat

        ret_struc = structure(name=name)
        ret_struc.read_file_3xN_matrix(self.file_path, self.file_type, coor_mat)
        
#         print 'rmsd', sup.ret_rmsd()
#         print 'angle', sup.ret_rotation_angle()
#         print 'axis', sup.ret_rotation_axis()
#         print 'centers', sup.ret_ref_av(), sup.ret_mv_av()
#         print 'matrix'
#         print sup.ret_rotation_matrix()
#         print 'orthogonal matrix?'
#         print numpy.dot(sup.ret_rotation_matrix(), sup.ret_rotation_matrix().transpose())

        return ret_struc
        
    def ret_renumbered_structure_file(self, ren_file, name=''):
        """
        Renumber the structure with input from <ren_file>.
        """
        ren_lines = open(ren_file, 'r').readlines()
        ren_list = [[eval(word) for word in file_handler.line_to_words(line)] for line in ren_lines]
        return self.ret_renumbered_structure(ren_list, name)

    def ret_renumbered_structure(self, perm_list, name=''):
        """
        Do permutations of atoms according to cyclic permutations in <perm_list>.
        This is important for molecules with symmetry or sigma bond rotations. A superposition can be done afterwards.
        example <perm_list>: [[1,2],[3],[6,7,8]...]
            switches 1 and 2, leaves 3 unchanged, and cyclically exchanges 6,7,8
        All atoms must be included exactly once.
        """
        
        if name == '': name = self.name
        
        start_mat = self.ret_3xN_matrix()

        num_at = self.mol.NumAtoms()
        
        trans_mat = numpy.zeros((num_at, num_at), float) # matrix multiplication with this matrix leads to the structure with permutated atom numbering
        for cycle in perm_list:
            for i in range(len(cycle)-1):
                trans_mat[cycle[i]-1, cycle[i+1]-1] = 1
            trans_mat[cycle[len(cycle)-1]-1, cycle[0]-1] = 1

        new_vec = numpy.dot(trans_mat, start_mat)

        ret_struc = structure(name=name)
        ret_struc.read_file_3xN_matrix(self.file_path, self.file_type, new_vec)

        return ret_struc

    def ret_bond_length(self, i, j):
        """
        Return the distance between atoms indexed i and j.
        """
        
        OBAtom_i = self.mol.GetAtom(i)
        OBAtom_j = self.mol.GetAtom(j)
        
        pos_i = numpy.array([OBAtom_i.x(), OBAtom_i.y(), OBAtom_i.z()])
        pos_j = numpy.array([OBAtom_j.x(), OBAtom_j.y(), OBAtom_j.z()])

        return numpy.dot(pos_i - pos_j, pos_i - pos_j)**.5
    
    def ret_bend(self, i, j, k):
        """
        Return the bending angle between atoms indexed i and j.
        """
        
        OBAtom_i = self.mol.GetAtom(i)
        OBAtom_j = self.mol.GetAtom(j)
        OBAtom_k = self.mol.GetAtom(k)
        
        pos_i = numpy.array([OBAtom_i.x(), OBAtom_i.y(), OBAtom_i.z()])
        pos_j = numpy.array([OBAtom_j.x(), OBAtom_j.y(), OBAtom_j.z()])
        pos_k = numpy.array([OBAtom_k.x(), OBAtom_k.y(), OBAtom_k.z()])

        vec1 = pos_i - pos_j
        vec2 = pos_k - pos_j
        
        len_1 = numpy.sqrt(numpy.dot(vec1,vec1))
        len_2 = numpy.sqrt(numpy.dot(vec2,vec2))

        return numpy.arccos(numpy.dot(vec1, vec2) / (len_1*len_2)) / numpy.pi * 180
    
    def ret_symbol(self, i):
        """
        Returns the symbol of atom i.
        """
        Z_symbol_dict = {1: 'H', 2: 'He',
3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar',
19: 'K', 20: 'Ca',
21:'Sc', 22:'Ti', 23:'V', 24:'Cr', 25:'Mn', 26:'Fe', 27:'Co', 28:'Ni', 29:'Cd', 30:'Zn',
31:'Ge', 32:'Ga', 33:'As', 34:'Se', 35:'Br', 36:'Kr',
37:'Rb', 38:'Sr',
39:'Y',  40:'Zr', 41:'Nb', 42:'Mo', 43:'Tc', 44:'Ru', 45:'Rh', 46:'Pd', 47:'Ag', 48:'Cd',
49:'In', 50:'Sn', 51:'Sb', 52:'Te',  53:'I', 54:'Xe',
55:'Cs', 56:'Ba',
57:'La', 72:'Hf', 73:'Ta',  74:'W', 75:'Re', 76:'Os', 77:'Ir', 78:'Pt', 79:'Au', 80:'Hg',
81:'Tl', 82:'Pb', 83:'Bi', 84:'Po', 85:'At', 86:'Rn'
}
        
        return Z_symbol_dict[self.mol.GetAtom(i).GetAtomicNum()]

    def ret_num_at(self):
        """
        Returns the number of atoms in the file.
        """
        return self.mol.NumAtoms()
    
    def ret_mass_vector(self, power):
        """
        Returns a vector with the masses of the atoms (each 1 time) taken to the <power> power.
        """
        mass_list = []
        for i in range(self.mol.NumAtoms()):
            atom = self.mol.GetAtom(i+1)
            mass_list += [atom.GetAtomicMass()**power]

        return numpy.array(mass_list, float)
        
    def make_coord_file(self, file_path, file_type='tmol'):
        """
        Write the structure file.
        """
        obconversion = openbabel.OBConversion()
        obconversion.SetOutFormat(file_type)
        obconversion.WriteFile(self.mol, file_path)

class mol_calc:
    """
    Class for doing calculations.
    A default file path 
    """
    def __init__(self, def_file_path, file_type='tmol'):
        # def_file_path is separately read in because with passing arguments the original objects would also be changed
            # maybe this would not occur if "return" would be used
        self.def_file_path = def_file_path
        self.file_type = file_type

    # return a set of structures
    def gram_schmidt(self, structures):
        """
        Returns an orthogonal system that spans the same space as *structures*.
        """
        b_list = [] # list of the orthogonal vectors

        for struc in structures:
            cn1 = struc.ret_vector()

            bn1 = cn1
            for bk in b_list:
                bn1 -= numpy.dot(cn1, bk) / numpy.dot(bk,bk) * bk

            bn1 = bn1 / numpy.dot(bn1,bn1)**.5
            b_list += [bn1]

        #print numpy.array([[numpy.dot(bi, bj) for bi in b_list] for bj in b_list])

        return [self.make_structure(bi) for bi in b_list]

        
    # return a structure
    def make_structure(self, vector, name=''):
        """
        Return the structure that corresponds to a vector.
        """
        ret_struc = structure(name)        
        ret_struc.read_file_vector(self.def_file_path, self.file_type, vector)
        return ret_struc
    
    def add(self, struc1, struc2):        
        return self.make_structure(struc1.ret_vector() + struc2.ret_vector())
        
    def subtract(self, struc1, struc2):
        #print struc1.ret_vector()
        #print struc2.ret_vector()
        return self.make_structure(struc1.ret_vector() - struc2.ret_vector())

    def scalar_mult(self, scalar, struc):    
        return self.make_structure(scalar * struc.ret_vector())

    def mean_structure(self, struc_list):
        """
        Return the mean structure of <struc_list>.
        """
        factor = 1. / len(struc_list)
        sum_struc = reduce(lambda struc1, struc2: self.add(struc1, struc2), struc_list)

        return self.scalar_mult(factor, sum_struc)

    def projected_structure(self, struc, ref_struc, diff_strucs, mass_wt=False):
        " Project the structure <struc> into the (hyper)plane defined by structures <ref_struc> and <diff_strucs> "
        shifted_struc = self.subtract(struc, ref_struc)

        o_diff_strucs = self.gram_schmidt(diff_strucs)

        ret_struc = ref_struc
        for o_struc in o_diff_strucs:
            factor = self.inner_product(o_struc, shifted_struc)
            plus_struc = self.scalar_mult(factor, o_struc)
            ret_struc = self.add(ret_struc, plus_struc)

        return ret_struc

    # return a scalar        
    def inner_product(self, struc1, struc2, mass_wt_pw=1):
        if mass_wt_pw == 0:
            return numpy.dot(struc1.ret_vector(), struc2.ret_vector())
        else:
            ## divide by the trace ???
            vec1 = numpy.dot(struc1.ret_vector(), self.ret_mass_matrix(power=mass_wt_pw))
            return numpy.dot(vec1, struc2.ret_vector())

    def norm(self, struc, mass_wt_pw=1):
        return numpy.sqrt(self.inner_product(struc, struc, mass_wt_pw))
        
    def distance(self, struc1, struc2, mass_wt_pw=1):
        """
        Return the distance between to structures in amu**(-1/2)*A. The distance as defined here is the RMSD times the squareroot of the molecular mass.
        """
        return self.norm(self.subtract(struc1, struc2), mass_wt_pw)

    def angle(self, struc1, struc2, mass_wt_pw=False):
        """
        Return the angle between two structure vectors or normal modes in degrees.
        """
        return numpy.arccos(self.inner_product(struc1, struc2, mass_wt_pw=mass_wt_pw) / self.norm(struc1, mass_wt_pw=mass_wt_pw) / self.norm(struc2, mass_wt_pw=mass_wt_pw)) / numpy.pi * 180

    # return different output
    def ret_mass_vector(self, power=1):
        """
        Returns a vector with the masses of the atoms (each 1 time) taken to the <power> power.
        """
        def_struc = structure()        
        def_struc.read_file(self.def_file_path, self.file_type) # mass weighing in this case not needed
        return def_struc.ret_mass_vector(power=power)
        
    def ret_mass_matrix(self, power=1):
        """
        Returns a diagonal matrix that contains the masses of the atoms (each 3 times).
        <power=.5> gives the squareroots which is typical mass weighting.
        """
        def_struc = structure()        
        def_struc.read_file(self.def_file_path, self.file_type) # mass weighing in this case not needed
        mass_list = []
        #for atom in openbabel.OBMolAtomIter(def_struc.mol):
        for i in range(def_struc.mol.NumAtoms()):
            atom = def_struc.mol.GetAtom(i+1)
            mass_list += 3*[atom.GetAtomicMass()**power]

        ret_mat = numpy.zeros((len(mass_list), len(mass_list)), dtype=float)
        for i in range(len(mass_list)):
            ret_mat[i,i]=mass_list[i]

        return ret_mat
    
    def distance_table(self, struc_list, mass_wt_pw, digits=4):
        """
        Return a table with the RMSDs between the structures in struc_list.
        <digits> specifies how many digits after the decimal point are printed out.
        """
        tm = file_handler.table_maker([5] + (len(struc_list)-1) * [digits+4])
        tm.write_line([''] + [struc.name for struc in struc_list[1:]])
        for i,struc1 in enumerate(struc_list[:-1]):
            tm.write_line([struc1.name] + [''] * i + [locale.format("%.*f", (digits, self.distance(struc1, struc2, mass_wt_pw=mass_wt_pw))) for struc2 in struc_list[i+1:]])
        return tm.return_table()
                 
if __name__ == '__main__':
    test_struc = structure()
    test_struc.read_file(file_path = '/home2/plasserf/BIP/opt/gs/b3lyp/SVP/coord', file_type='tmol')
    print(test_struc.ret_3xN_matrix())
