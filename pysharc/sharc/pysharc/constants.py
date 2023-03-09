# ******************************************
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
# ******************************************


# Number of frozen core orbitals

FROZENS = {'H': 0, 'He': 0,
           'Li': 1, 'Be': 1, 'B': 1, 'C': 1, 'N': 1, 'O': 1, 'F': 1, 'Ne': 1,
           'Na': 1, 'Mg': 1, 'Al': 5, 'Si': 5, 'P': 5, 'S': 5, 'Cl': 5, 'Ar': 5,
           'K': 5, 'Ca': 5,
           'Sc': 5, 'Ti': 5, 'V': 5, 'Cr': 5, 'Mn': 5, 'Fe': 5, 'Co': 5, 'Ni': 5, 'Cu': 5, 'Zn': 5,
           'Ga': 9, 'Ge': 9, 'As': 9, 'Se': 9, 'Br': 9, 'Kr': 9,
           'Rb': 9, 'Sr': 9,
           'Y': 14, 'Zr': 14, 'Nb': 14, 'Mo': 14, 'Tc': 14, 'Ru': 14, 'Rh': 14, 'Pd': 14, 'Ag': 14, 'Cd': 14,
           'In': 18, 'Sn': 18, 'Sb': 18, 'Te': 18, 'I': 18, 'Xe': 18,
           'Cs': 18, 'Ba': 18,
           'La': 23,
           'Ce': 23, 'Pr': 23, 'Nd': 23, 'Pm': 23, 'Sm': 23, 'Eu': 23, 'Gd': 23, 'Tb': 23, 'Dy': 23, 'Ho': 23, 'Er': 23, 'Tm': 23, 'Yb': 23, 'Lu': 23,
           'Hf': 23, 'Ta': 23, 'W': 23, 'Re': 23, 'Os': 23, 'Ir': 23, 'Pt': 23, 'Au': 23, 'Hg': 23,
           'Tl': 23, 'Pb': 23, 'Bi': 23, 'Po': 23, 'At': 23, 'Rn': 23,
           'Fr': 30, 'Ra': 30,
           'Ac': 30,
           'Th': 30, 'Pa': 30, 'U': 30, 'Np': 30, 'Pu': 30, 'Am': 30, 'Cm': 30, 'Bk': 30, 'Cf': 30, 'Es': 30, 'Fm': 30, 'Md': 30, 'No': 30, 'Lr': 30,
           'Rf': 30, 'Db': 30, 'Sg': 30, 'Bh': 30, 'Hs': 30, 'Mt': 30, 'Ds': 30, 'Rg': 30, 'Cn': 30,
           'Nh': 39, 'Fl': 39, 'Mc': 39, 'Lv': 39, 'Ts': 39, 'Og': 39
           }

FROZEN_ID = {1: 0, 2: 0, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1, 8: 1, 9: 1, 10: 1, 11: 1, 12: 1, 13: 5, 14: 5, 15: 5, 16: 5, 17: 5, 18: 5, 19: 5, 20: 5, 21: 5, 22: 5, 23: 5, 24: 5, 25: 5, 26: 5, 27: 5, 28: 5, 29: 5, 30: 5, 31: 9, 32: 9, 33: 9, 34: 9, 35: 9, 36: 9, 37: 9, 38: 9, 39: 14, 40: 14, 41: 14, 42: 14, 43: 14, 44: 14, 45: 14, 46: 14, 47: 14, 48: 14, 49: 18, 50: 18, 51: 18, 52: 18, 53: 18, 54: 18, 55: 18, 56: 18, 57: 23, 58: 23, 59: 23, 60: 23, 61: 23, 62: 23, 63: 23, 64: 23, 65: 23, 66: 23, 67: 23, 68: 23, 69: 23, 70: 23, 71: 23, 72: 23, 73: 23, 74: 23, 75: 23, 76: 23, 77: 23, 78: 23, 79: 23, 80: 23, 81: 23, 82: 23, 83: 23, 84: 23, 85: 23, 86: 23, 87: 30, 88: 30, 89: 30, 90: 30, 91: 30, 92: 30, 93: 30, 94: 30, 95: 30, 96: 30, 97: 30, 98: 30, 99: 30, 100: 30, 101: 30, 102: 30, 103: 30, 104: 30, 105: 30, 106: 30, 107: 30, 108: 30, 109: 30, 110: 30, 111: 30, 112: 30, 113: 39, 114: 39, 115: 39, 116: 39, 117: 39, 118: 39}
#
ATOMCHARGE = {'H': 1, 'He': 2,
              'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
              'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18,
              'K': 19, 'Ca': 20,
              'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30,
              'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36,
              'Rb': 37, 'Sr': 38,
              'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48,
              'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54,
              'Cs': 55, 'Ba': 56,
              'La': 57,
              'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71,
              'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80,
              'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86,
              'Fr': 87, 'Ra': 88,
              'Ac': 89,
              'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103,
              'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112,
              'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118
              }

IAn2AName = {
    1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F",
    10: "Ne", 11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P", 16: "S", 17: "Cl", 18: "Ar",
    19: "K", 20: "Ca", 21: "Sc", 22: "Ti", 23: "V", 24: "Cr", 25: "Mn", 26: "Fe", 27: "Co",
    28: "Ni", 29: "Cu", 30: "Zn", 31: "Ga", 32: "Ge", 33: "As", 34: "Se", 35: "Br", 36: "Kr",
    37: "Rb", 38: "Sr", 39: "Y", 40: "Zr", 41: "Nb", 42: "Mo", 43: "Tc", 44: "Ru", 45: "Rh", 46: "Pd",
    47: "Ag", 48: "Cd", 49: "In", 50: "Sn", 51: "Sb", 52: "Te", 53: "I", 54: "Xe", 55: "Cs", 56: "Ba",
    57: "La", 58: "Ce", 59: "Pr", 60: "Nd", 61: "Pm", 62: "Sm", 63: "Eu", 64: "Gd", 65: "Tb", 66: "Dy",
    67: "Ho", 68: "Er", 69: "Tm", 70: "Yb", 71: "Lu", 72: "Hf", 73: "Ta", 74: "W", 75: "Re", 76: "Os", 77: "Ir",
    78: "Pt", 79: "Au", 80: "Hg", 81: "Tl", 82: "Pb", 83: "Bi", 84: "Po", 85: "At", 86: "Rn", 87: "Fr",
    88: "Ra", 89: "Ac", 90: "Th", 91: "Pa", 92: "U", 93: "Np", 94: "Pu", 95: "Am", 96: "Cm", 97: "Bk",
    98: "Cf", 99: "Es", 100: "Fm", 101: "Md", 102: "No", 103: "Lr", 104: "Rf", 105: "Db", 106: "Sg",
    107: "Bh", 108: "Hs", 109: "Mt", 110: "Ds", 111: "Rg", 112: "Cn", 113: "Nh", 114: "Fl",
    115: "Mc", 116: "Lv", 117: "Ts", 118: "Og"
}

# hash table for conversion of multiplicity to the keywords used in MOLCAS
IToMult = {
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

# hash table for conversion of polarisations to the keywords used in MOLCAS
IToPol = {
    0: 'X',
    1: 'Y',
    2: 'Z',
    'X': 0,
    'Y': 1,
    'Z': 2
}
