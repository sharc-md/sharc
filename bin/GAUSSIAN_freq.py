#!/usr/bin/env python3

# Crude script for the conversion of
# Gaussian output files to frequency files
# for molden. Only works with freq=hpmodes

from sys import argv
#from constants import NUMBERS, au2a

au2a = 0.529177211
NUMBERS = {
    'H': 1,
    'He': 2,
    'Li': 3,
    'Be': 4,
    'B': 5,
    'C': 6,
    'N': 7,
    'O': 8,
    'F': 9,
    'Ne': 10,
    'Na': 11,
    'Mg': 12,
    'Al': 13,
    'Si': 14,
    'P': 15,
    'S': 16,
    'Cl': 17,
    'Ar': 18,
    'K': 19,
    'Ca': 20,
    'Sc': 21,
    'Ti': 22,
    'V': 23,
    'Cr': 24,
    'Mn': 25,
    'Fe': 26,
    'Co': 27,
    'Ni': 28,
    'Cu': 29,
    'Zn': 30,
    'Ga': 31,
    'Ge': 32,
    'As': 33,
    'Se': 34,
    'Br': 35,
    'Kr': 36,
    'Rb': 37,
    'Sr': 38,
    'Y': 39,
    'Zr': 40,
    'Nb': 41,
    'Mo': 42,
    'Tc': 43,
    'Ru': 44,
    'Rh': 45,
    'Pd': 46,
    'Ag': 47,
    'Cd': 48,
    'In': 49,
    'Sn': 50,
    'Sb': 51,
    'Te': 52,
    'I': 53,
    'Xe': 54,
    'Cs': 55,
    'Ba': 56,
    'La': 57,
    'Ce': 58,
    'Pr': 59,
    'Nd': 60,
    'Pm': 61,
    'Sm': 62,
    'Eu': 63,
    'Gd': 64,
    'Tb': 65,
    'Dy': 66,
    'Ho': 67,
    'Er': 68,
    'Tm': 69,
    'Yb': 70,
    'Lu': 71,
    'Hf': 72,
    'Ta': 73,
    'W': 74,
    'Re': 75,
    'Os': 76,
    'Ir': 77,
    'Pt': 78,
    'Au': 79,
    'Hg': 80,
    'Tl': 81,
    'Pb': 82,
    'Bi': 83,
    'Po': 84,
    'At': 85,
    'Rn': 86,
    'Fr': 87,
    'Ra': 88,
    'Ac': 89,
    'Th': 90,
    'Pa': 91,
    'U': 92,
    'Np': 93,
    'Pu': 94,
    'Am': 95,
    'Cm': 96,
    'Bk': 97,
    'Cf': 98,
    'Es': 99,
    'Fm': 100,
    'Md': 101,
    'No': 102,
    'Lr': 103,
    'Rf': 104,
    'Db': 105,
    'Sg': 106,
    'Bh': 107,
    'Hs': 108,
    'Mt': 109,
    'Ds': 110,
    'Rg': 111,
    'Cn': 112,
    'Nh': 113,
    'Fl': 114,
    'Mc': 115,
    'Lv': 116,
    'Ts': 117,
    'Og': 118
}

INV_NUMBERS = {v: k for k, v in NUMBERS.items()}


# Check arguments
if len(argv) != 2:
    print("Usage: GAUSSIAN_freq.py <gaussian.log>\n:")
    print("Convert Gaussian output file to molden file for")
    print("normal mode visualisation.")
    exit()

name, gaussian_file = argv

try:
    lines = open(gaussian_file, 'r').readlines()
except IOError:
    print("Could not open %s." % gaussian_file)
    exit()

# check if file is sucessfully completed orca file:
is_gaussian = False
finished = False
if "Entering Gaussian System" in lines[0]:
    is_gaussian = True
for line in lines:
    if 'hpmodes' in line.lower():
        finished = True

if not is_gaussian:
    print("File %s is not in gaussian output format (probably)!" % gaussian_file)
    exit()
elif is_gaussian and not finished:
    print("Run the job with freq=hpmodes!")
    exit()
elif is_gaussian and finished:
    print("Reading data from file %s..." % gaussian_file)

# Standard orientation: (from bottom)
# get coordinates
for iline, line in enumerate(lines[::-1]):
    original_index = len(lines) - 1 - iline
    if "Standard orientation:" in line:
        break
iline = original_index + 4
coords = []
while True:
    iline += 1
    line = lines[iline]
    if '---' in line:
        break
    s = line.split()
    atom = [ INV_NUMBERS[int(s[1])], float(s[3]), float(s[4]), float(s[5]) ]
    coords.append(atom)
natom = len(coords)
nfreq = 3*natom-6


# freq, int, modes
for iline, line in enumerate(lines):
    if "Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering" in line:
        break
iline += 4
freqs = []
modes = []
ints = []
while True:
    if "Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering" in lines[iline]:
        break
    s = lines[iline].split()
    nhere = len(s)
    for imode in range(nhere):
        freq =   float(lines[iline+2  ].split()[2+imode])
        ir   =   float(lines[iline+5  ].split()[3+imode])
        mode = [ float(lines[iline+7+i].split()[3+imode]) for i in range(3*natom) ]
        freqs.append(freq)
        ints.append(ir)
        modes.append(mode)
    iline += 7+3*natom

# generate molden file
out_file = gaussian_file + '.molden'
out = open(out_file, 'w')
out.write("[MOLDEN FORMAT]\n")
# write frequencies
out.write("[FREQ]\n")
for freq in freqs:
    out.write(str(freq) + '\n')
# write coordinates block (A.U.)
out.write("[FR-COORD]\n")
for coord in coords:
    out.write(coord[0] + ' ' + ' '.join([str(i/au2a) for i in coord[1:4]]) + '\n')
# write normal modes:
out.write("[FR-NORM-COORD]\n")
for i in range(nfreq):
    out.write("vibration %d\n" % (i + 1))
    for j in range(len(modes[i])):
        out.write(str(modes[i][j]) + ' ')
        if (j + 1) % 3 == 0:
            out.write('\n')
out.write('[INT]\n')
for i in range(nfreq):
    out.write('%16.9f\n' % ints[i])
out.close()
print("Molden output written to %s" % out_file)




