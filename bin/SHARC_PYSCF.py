#!/usr/bin/env python3

# ******************************************
#
#    SHARC Program Suite
#
#    Copyright (c) 2023 University of Vienna
#    Copyright (c) 2023 University of Minnesota
#    Copyright (c) 2024 University of Chicago
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

import os
import sys
import re
import datetime

from socket import gethostname

_version_ = "3.0"
_versiondate_ = datetime.date(2024, 2, 2)
_change_log_ = """"""

START_TIME = datetime.datetime.now()

# Global Variables for printing.
DEBUG = False  # Raw output
PRINT = True  # Formatted output

# conversion factors
AU_TO_ANG = 0.529177211
rcm_to_Eh = 4.556335e-6

def print_header():
    """Prints the formatted header of the log file. Prints version number and version date"""
    print(START_TIME, gethostname(), os.getcwd())
    if not PRINT:
        return
    string = "\n"
    string += "  " + "=" * 80 + "\n"
    string += "||" + " " * 80 + "||\n"
    string += "||" + " " * 27 + "SHARC - PySCF - Interface" + " " * 28 + "||\n"
    string += "||" + " " * 80 + "||\n"
    string += "||" + " " * 25 + "Authors: Matthew R. Hennefarth" + " " * 25 + "||\n"
    string += "||" + " " * 80 + "||\n"
    string += (
        "||"
        + " " * (36 - (len(_version_) + 1) // 2)
        + f"Version: {_version_}"
        + " " * (35 - (len(_version_)) // 2)
        + "||\n"
    )
    lens = len(_versiondate_.strftime("%d.%m.%y"))
    string += (
        "||"
        + " " * (37 - lens // 2)
        + f"Date: {_versiondate_.strftime('%d.%m.%y')}"
        + " " * (37 - (lens + 1) // 2)
        + "||\n"
    )
    string += "||" + " " * 80 + "||\n"
    string += "  " + "=" * 80 + "\n\n"
    print(string)
    if DEBUG:
        print(_change_log_)

def itnmstates(states):

    for i in range(len(states)):
        if states[i] < 1:
            continue
        for k in range(i + 1):
            for j in range(states[i]):
                yield i + 1, j + 1, k - i / 2.
    return

def get_pairs(lines, i):
    nacpairs = []
    while True:
        i += 1
        try:
            line = lines[i].lower()
        except IndexError:
            print('"keyword select" has to be completed with an "end" on another line!')
            sys.exit(1)
        if "end" in line:
            break

        fields = line.split()
        try:
            nacpairs.append([int(fields[0]), int(fields[1])])

        except ValueError:
            print(
                '"nacdr select" is followed by pairs of state indices, each pair on a new line!'
            )
            sys.exit(1)

    return nacpairs, i


def readQMin(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    qmin = {}
    try:
        natom = int(lines[0])
    
    except ValueError as e:
        print("First line must contain the number of atoms!")
        raise e

    qmin["natom"] = natom
    if len(lines) < natom + 4:
        print(
            """Input file must contain at least:
natom
comment
geometry
keyword 'states'
at least one task"""
        )
        sys.exit(1)

    qmin["comment"] = lines[1]

    # Get geometry and possibly velocity (for backup-analytical nonadiabatic couplings)
    qmin["geo"] = []
    qmin["veloc"] = []
    has_veloc = True
    for line in lines[2 : natom + 2]:
        line = line.split()
        line[1:] = [float(x) for x in line[1:]]
        qmin["geo"].append(line[:4])
        if len(line) >= 7:
            qmin["veloc"].append(line[4:7])

        else:
            has_veloc = False

    if not has_veloc:
        del qmin["veloc"]

    # TODO: I hate this and this should be fixed up....ugh!
    i = natom + 1
    while i + 1 < len(lines):
        i += 1
        line = lines[i]
        line = re.sub("#.*$", "", line)
        if len(line.split()) == 0:
            continue

        key = line.lower().split()[0]
        if "savedir" in key:
            args = line.split()[1:]

        else:
            args = line.lower().split()[1:]

        if key in qmin:
            print(f"Repeated keyword {key} in line {i+1} in input file!")
            continue

        if len(args) >= 1 and "select" in args[0]:
            pairs, i = get_pairs(lines, i)
            qmin[key] = pairs
        else:
            qmin[key] = args

    if "unit" in qmin:
        if qmin['unit'][0] == 'angstrom':
            factor = 1.0/AU_TO_ANG
        
        elif qmin['unit'][0] == 'bohr':
            factor = 1.0

        else:
            print(f"Don't know input unit {qmin['unit'][0]}")
            sys.exit(1)

    else:
        factor = 1.0/AU_TO_ANG
   
    for atom in qmin['geo']:
        atom[1:] = [xyz*factor for xyz in atom[1:]]

    if "states" not in qmin:
        print("Keyword 'states' not given!")
        sys.exit(1)

    qmin["states"] = [int(state) for state in qmin["states"]]
    
    reduc = 0
    for i in reversed(qmin['states']):
        if i == 0:
            reduc += 1

        else: 
            break

    if reduc > 0:
        qmin["states"] = qmin["states"][:-reduc]

    nstates = 0
    nmstates = 0
    for index, state in enumerate(qmin["states"]):
        nstates += state
        nmstates += state * (index + 1)
    
    qmin['nstates'] = nstates
    qmin['nmstates'] = nmstates

    possible_tasks = ["h", "soc", "dm", "grad", "overlap", "dmdr", "socdr", "ion", "phases"]
    if not any([i in qmin for i in possible_tasks]):
        print(f"No tasks found! Tasks are {possible_tasks}")
        sys.exit(1)

    if 'samestep' in qmin and 'init' in qmin:
        print("'init' and 'samestep' cannot both be present in inputfile")
        sys.exit(1)

    if 'phases'in qmin:
        qmin["overlap"] = []

    if "overlap" in qmin and "init" in qmin:
        print("'overlap' and 'phases' cannot both be calculated in the first timestep")
        sys.exit(1)

    if 'init' not in qmin and 'samestep' not in qmin:
        qmin['newstep'] = []

    if not any([i in qmin for i in ["h", "soc", "dm", "grad"]]) and "overlap" in qmin:
        qmin["h"] = []

    if len(qmin["states"]) > 8:
        print("Higher multiplicities than octets are not supported!")
        sys.exit(1)

    not_implemented_tasks = ["soc", "overlap", "nacdt", "dmdr", "ion", "theodore"]
    for task in not_implemented_tasks:
        if task in qmin:
            print(f"Within the SHARC-PySCF interface, '{task}' is not supported")
            sys.exit(1)


    if "h" in qmin and "soc" in qmin:
        del qmin["h"]

    if "molden" in qmin and "samestep" in qmin:
        print("HINT: not producing Molden files in 'samestep' mode!")
        del qmin['molden']

    if 'grad' in qmin:
        if len(qmin['grad']) == 0 or qmin['grad'][0] == "all":
            qmin['grad'] = [i+1 for i in range(nmstates)]
        
        else:
            for i in range(len(qmin['grad'])):
                try:
                    qmin['grad'][i] = int(qmin['grad'][i])
                
                except ValueError:
                    print("Arguments to keyword 'grad' must be 'all' or a list of integers")
                    sys.exit(1)

                if qmin['grad'][i] > nmstates:
                    print("State for requested gradient does not correspond to any state in QM input file state list!")
                    sys.exit(1)

    if 'overlap' in qmin:
        if len(qmin['overlap']) >= 1:
            overlap_pairs = qmin['overlap']
            for pair in overlap_pairs:
                if pair[0] > nmstates or pair[1] > nmstates:
                    print("State for requested overlap does not correspond to any state in QM input file state list!")
                    sys.exit(1)

        else:
            qmin['overlap'] = [[j+i, i+1] for i in range(nmstates) for j in range(i+1)]

    if 'nacdr' in qmin:
        if len(qmin['nacdr']) >= 1:
            nac_pairs = qmin['nacdr']
            for pair in nac_pairs:
                if pair[0] > nmstates or pair[1] > nmstates:
                    print("State for requested nacdr does not correspond to any state in QM input file state list!")

        else:
            qmin['nacdr'] = [[j+1, i+1] for i in range(nmstates) for j in range(i)]

    # obtain the statemap
    statemap = {}
    i = 1
    for imult, istate, ims in itnmstates(qmin['states']):
        statemap[i] = [imult, istate, ims]
        i += 1
    
    qmin['statemap'] = statemap

    gradmap = set()
    if 'grad' in qmin:
        for state in qmin['grad']:
            gradmap.add(tuple(statemap[state][0:2]))

    gradmap = sorted(gradmap)
    qmin["gradmap"] = gradmap

    nacmap = set()
    if 'nacdr' in qmin:
        for pair in qmin['nacdr']:
            s1 = statemap[pair[0]][:-1]
            s2 = statemap[pair[1]][:-1]
            if s1[0] != s2[0] or s1 == s2:
                continue
            nacmap.add(tuple(s1 + s2))

    nacmap = list(nacmap)
    nacmap.sort()
    qmin['nacmap'] = nacmap

    # TODO from 2222 of SHARC_MOLCAS, the MOLCAS.resources file loading stuff...

    return qmin


def main():
    try:
        env_print = os.getenv("SH2CAS_PRINT")
        if env_print and env_print.lower() == "false":
            global PRINT
            PRINT = False

        env_debug = os.getenv("SH2CAS_DEBUG")
        if env_debug and env_debug.lower() == "true":
            global DEBUG
            DEBUG = True

    except ValueError:
        print(
            "SH2CAS_PRINT or SH2CAS_DEBUG environment variables do no evaluate to booleans!"
        )

    if len(sys.argv) != 2:
        print(
            f"""Usage:
./SHARC_PYSCF.py <QMin>
version: {_version_}
date: {_versiondate_}
changelog: {_change_log_}"""
        )
        sys.exit(1)

    qmin_filename = sys.argv[1]

    print_header()
    qmin = readQMin(qmin_filename)
    import json

    print(json.dumps(qmin, indent=2, sort_keys=False))


if __name__ == "__main__":
    main()
