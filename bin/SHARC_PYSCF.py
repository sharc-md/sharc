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

from multiprocessing import Pool
import os
import shutil
import sys
import re
import datetime
import pprint
import numpy as np
from copy import deepcopy

from socket import gethostname
import time

from pyscf import lib, gto, mcscf

_version_ = "3.0"
_versiondate_ = datetime.date(2024, 2, 2)
_change_log_ = """"""

START_TIME = datetime.datetime.now()

# Global Variables for printing.
DEBUG = False  # Raw output
PRINT = True  # Formatted output

IToMult = {
    1: "Singlet",
    2: "Doublet",
    3: "Triplet",
    4: "Quartet",
    5: "Quintet",
    6: "Sextet",
    7: "Septet",
    8: "Octet",
    "Singlet": 1,
    "Doublet": 2,
    "Triplet": 3,
    "Quartet": 4,
    "Quintet": 5,
    "Sextet": 6,
    "Septet": 7,
    "Octet": 8,
}

# conversion factors
AU_TO_ANG = 0.529177211
rcm_to_Eh = 4.556335e-6


def eformat(f, prec, exp_digits):
    """Formats a float f into scientific notation with prec number of decimals and exp_digits number of exponent digits.

    String looks like:
    [ -][0-9]\\.[0-9]*E[+-][0-9]*

    Arguments:
    1 float: Number to format
    2 integer: Number of decimals
    3 integer: Number of exponent digits

    Returns:
    1 string: formatted number"""

    s = "% .*e" % (prec, f)
    mantissa, exp = s.split("e")
    return "%sE%+0*d" % (mantissa, exp_digits + 1, int(exp))


def measure_time():
    """Calculates the time difference between global variable starttime and the time of the call of measuretime.

    Prints the Runtime, if PRINT or DEBUG are enabled.

    Arguments:
    none

    Returns:
    1 float: runtime in seconds"""

    endtime = datetime.datetime.now()
    runtime = endtime - START_TIME
    if PRINT or DEBUG:
        hours = runtime.seconds // 3600
        minutes = runtime.seconds // 60 - hours * 60
        seconds = runtime.seconds % 60
        print(
            "==> Runtime:\n%i Days\t%i Hours\t%i Minutes\t%i Seconds\n\n"
            % (runtime.days, hours, minutes, seconds)
        )
    total_seconds = (
        runtime.days * 24 * 3600 + runtime.seconds + runtime.microseconds // 1.0e6
    )
    return total_seconds


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


def print_qmin(qmin):
    if DEBUG:
        pprint.pprint(qmin)

    if not PRINT:
        return

    print(f"==> QMin Job description for:\n{qmin['comment']}")
    TASK_TO_STRING = {
        "h": "H",
        "soc": "SOC",
        "dm": "DM",
        "grad": "Grad",
        "nacdr": "Nac(ddr)",
        "nacdt": "Nac(ddt)",
        "overlap": "Overlaps",
        "angular": "Angular",
        "ion": "Dyson norms",
        "dmdr": "DM-Grad",
        "socdr": "SOC-Grad",
        "phases": "Phases",
    }
    output = "Tasks: "
    for key, string in TASK_TO_STRING.items():
        if key in qmin:
            output += f"\t{string}"

    print(output)

    output = "States: "
    for i in itmult(qmin["states"]):
        output += f"\t{qmin['states'][i-1]} {IToMult[i]}"

    print(output)

    output = "Method: \t"
    tmp = "|".join([str(r) for r in qmin["template"]["roots"]])
    if qmin["template"]["method"].upper() == "L-PDFT":
        output += f"L({tmp})-PDFT"

    elif qmin["template"]["method"].upper() == "MC-PDFT":
        output += f"SA({tmp})-PDFT"

    else:
        output += f"SA({tmp})-{qmin['template']['method'].upper()}"

    output += f"({qmin['template']['nelecas']}, {qmin['template']['ncas']})/{qmin['template']['basis']}"
    print(output)

    oddmults = False
    for i in qmin["statemap"].values():
        if (qmin["template"]["nelecas"] + i[0]) % 2 == 0:
            oddmults = True

    if oddmults:
        output = "\t\t" + ["Even ", "Odd "][qmin["template"]["nelecas"]] % 2 == 0
        output += f"numbers of electrons are treated wiht CAS({qmin['template']['nelecas']}, {qmin['template']['ncas']})"
        print(output)

    output = "Found Geo"
    if "veloc" in qmin:
        output += " and Veloc! "

    else:
        output += "! "

    output += f"NAtom is {qmin['natom']}.\n"
    print(output)

    output = "\nGeometry in Bohrs:\n"
    for atom in qmin["geo"]:
        element = atom[0]
        coords = atom[1:]
        output += f"{element} "
        for x in coords:
            output += f"{x:7.4f} "

        output += "\n"

    print(output)

    if "veloc" in qmin:
        output = ""
        for index, veloc in enumerate(qmin["veloc"]):
            element = qmin["geo"][index][0]
            output += f"{element} "
            for v in veloc:
                output += f"{v:7.4f} "

            output += "\n"

        print(output)

    if "grad" in qmin:
        output = "Gradients:   "
        for i in range(1, qmin["nmstates"] + 1):
            if i in qmin["grad"]:
                output += "X"

            else:
                output += "."

        output += "\n"
        print(output)

    if "nacdr" in qmin:
        output = "Nonadiabatic couplings:\n"
        for i in range(1, qmin["nmstates"] + 1):
            for j in range(1, qmin["nmstates"] + 1):
                if [i, j] in qmin["nacdr"] or [j, i] in qmin["nacdr"]:
                    output += "X"

                else:
                    output += "."
            output += "\n"
        print(output)

    if "overlap" in qmin:
        output = "Overlaps:\n"
        for i in range(1, qmin["nmstates"] + 1):
            for j in range(1, qmin["nmstates"] + 1):
                if (i, j) in qmin["overlap"] or (j, i) in qmin["overlap"]:
                    output += "X"

                else:
                    output += "."
            output += "\n"
        print(output)

    print("\n")
    sys.stdout.flush()


def itmult(states):
    for i in range(len(states)):
        if states[i] < 1:
            continue
        yield i + 1
    return


def itnmstates(states):
    for i in range(len(states)):
        if states[i] < 1:
            continue
        for k in range(i + 1):
            for j in range(states[i]):
                yield i + 1, j + 1, k - i / 2.0
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


def removequotes(string):
    if string.startswith("'") and string.endswith("'"):
        return string[1:-1]
    elif string.startswith('"') and string.endswith('"'):
        return string[1:-1]
    else:
        return string


def getsh2caskey(sh2cas, key):
    for line in sh2cas:
        line = re.sub("#.*$", "", line)
        line = line.split(None, 1)
        if line == []:
            continue

        if key.lower() in line[0].lower():
            return line

    return ["", ""]


def get_sh2cas_environ(sh2cas, key, environ=True, crucial=True):
    line = getsh2caskey(sh2cas, key)
    if line[0]:
        line = line[1]
        line = removequotes(line).strip()
    else:
        if environ:
            line = os.getenv(key.upper())
            if not line:
                if crucial:
                    print(
                        f"Either set ${key.upper()} or give path to {key.upper()} in PYSCF.resources"
                    )
                    sys.exit(1)

                else:
                    return ""

        else:
            if crucial:
                print(f"Give path to {key.upper()} in PYSCF.resources")
                sys.exit(1)

            else:
                return ""

    line = os.path.expandvars(line)
    line = os.path.expanduser(line)
    if ";" in line:
        print(
            f"${key.upper()} contains a semicolon. Do you probably want to execute another command after {key.upper()}? I can't do that for you..."
        )
        sys.exit(1)
    return line


def check_directory(dir):
    """Checks where dir is a file or directory. If a file, quits with exit code 1. If a directory, it passes. If does not exist, then we try and create the directory"""

    if os.path.exists(dir):
        if not os.path.isdir(dir):
            print(f"{dir} exists but is not a directory! Quiting...")
            sys.exit(1)

    else:
        os.makedirs(dir)

    return True


def get_version():
    from pyscf import __version__ as pyscf_version

    min_version = (2, 6, 0)
    line = pyscf_version.split(".")
    line = [int(i) for i in line]
    for min, actual in zip(min_version, line):
        if actual < min:
            print(f"PySCF version {pyscf_version} not supported!")
            sys.exit(1)

        if actual > min:
            break

    if DEBUG:
        print(f"PySCF version {pyscf_version}")

    return pyscf_version


def readqmin(filename):
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
        element = line[0]
        coords = [float(x) for x in line[1:]]
        qmin["geo"].append([element] + coords[:3])
        if len(coords) >= 6:
            qmin["veloc"].append(coords[3:6])

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
        if qmin["unit"][0] == "angstrom":
            factor = 1.0 / AU_TO_ANG

        elif qmin["unit"][0] == "bohr":
            factor = 1.0

        else:
            print(f"Don't know input unit {qmin['unit'][0]}")
            sys.exit(1)

    else:
        factor = 1.0 / AU_TO_ANG

    for atom in qmin["geo"]:
        atom[1:] = [xyz * factor for xyz in atom[1:]]

    if "states" not in qmin:
        print("Keyword 'states' not given!")
        sys.exit(1)

    qmin["states"] = [int(state) for state in qmin["states"]]

    reduc = 0
    for i in reversed(qmin["states"]):
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

    qmin["nstates"] = nstates
    qmin["nmstates"] = nmstates

    possible_tasks = [
        "h",
        "soc",
        "dm",
        "grad",
        "overlap",
        "dmdr",
        "socdr",
        "ion",
        "phases",
    ]
    if not any([i in qmin for i in possible_tasks]):
        print(f"No tasks found! Tasks are {possible_tasks}")
        sys.exit(1)

    if "samestep" in qmin and "init" in qmin:
        print("'init' and 'samestep' cannot both be present in inputfile")
        sys.exit(1)

    if "phases" in qmin:
        qmin["overlap"] = []

    if "overlap" in qmin and "init" in qmin:
        print("'overlap' and 'phases' cannot both be calculated in the first timestep")
        sys.exit(1)

    if "init" not in qmin and "samestep" not in qmin:
        qmin["newstep"] = []

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
        del qmin["molden"]

    if "grad" in qmin:
        if len(qmin["grad"]) == 0 or qmin["grad"][0] == "all":
            qmin["grad"] = [i + 1 for i in range(nmstates)]

        else:
            for i in range(len(qmin["grad"])):
                try:
                    qmin["grad"][i] = int(qmin["grad"][i])

                except ValueError:
                    print(
                        "Arguments to keyword 'grad' must be 'all' or a list of integers"
                    )
                    sys.exit(1)

                if qmin["grad"][i] > nmstates:
                    print(
                        "State for requested gradient does not correspond to any state in QM input file state list!"
                    )
                    sys.exit(1)

    if "overlap" in qmin:
        if len(qmin["overlap"]) >= 1:
            overlap_pairs = qmin["overlap"]
            for pair in overlap_pairs:
                if pair[0] > nmstates or pair[1] > nmstates:
                    print(
                        "State for requested overlap does not correspond to any state in QM input file state list!"
                    )
                    sys.exit(1)

        else:
            qmin["overlap"] = [
                [j + i, i + 1] for i in range(nmstates) for j in range(i + 1)
            ]

    if "nacdr" in qmin:
        if len(qmin["nacdr"]) >= 1:
            nac_pairs = qmin["nacdr"]
            for pair in nac_pairs:
                if pair[0] > nmstates or pair[1] > nmstates:
                    print(
                        "State for requested nacdr does not correspond to any state in QM input file state list!"
                    )

        else:
            qmin["nacdr"] = [[j + 1, i + 1] for i in range(nmstates) for j in range(i)]

    # obtain the statemap
    statemap = {}
    i = 1
    for imult, istate, ims in itnmstates(qmin["states"]):
        statemap[i] = [imult, istate, ims]
        i += 1

    qmin["statemap"] = statemap

    gradmap = set()
    if "grad" in qmin:
        for state in qmin["grad"]:
            gradmap.add(tuple(statemap[state][0:2]))

    gradmap = sorted(gradmap)
    qmin["gradmap"] = gradmap

    nacmap = set()
    if "nacdr" in qmin:
        for pair in qmin["nacdr"]:
            s1 = statemap[pair[0]][:-1]
            s2 = statemap[pair[1]][:-1]
            if s1[0] != s2[0] or s1 == s2:
                continue
            nacmap.add(tuple(s1 + s2))

    nacmap = list(nacmap)
    nacmap.sort()
    qmin["nacmap"] = nacmap

    # TODO from 2222 of SHARC_MOLCAS, the MOLCAS.resources file loading stuff...

    pyscf_resource_filename = "PYSCF.resources"
    with open(pyscf_resource_filename, "r") as f:
        sh2cas = f.readlines()

    qmin["pwd"] = os.getcwd()

    line = get_sh2cas_environ(sh2cas, "scratchdir", environ=False, crucial=False)
    if line is None:
        line = os.path.join(qmin["pwd"], "SCRATCHDIR")

    line = os.path.expandvars(line)
    line = os.path.expanduser(line)
    line = os.path.abspath(line)
    qmin["scratchdir"] = line

    # Set up savedir
    if "savedir" in qmin:
        line = qmin["savedir"][0]

    else:
        line = get_sh2cas_environ(sh2cas, "savedir", environ=False, crucial=False)
        if line is None or line == "":
            line = os.path.join(qmin["pwd"], "SAVEDIR")

    line = os.path.expandvars(line)
    line = os.path.expanduser(line)
    line = os.path.abspath(line)

    if "init" in qmin:
        check_directory(line)

    qmin["savedir"] = line

    line = getsh2caskey(sh2cas, "debug")
    if line[0]:
        if len(line) <= 1 or "true" in line[1].lower():
            global DEBUG
            DEBUG = True

    line = getsh2caskey(sh2cas, "no_print")
    if line[0]:
        if len(line) <= 1 or "true" in line[1].lower():
            global PRINT
            PRINT = False

    qmin["memory"] = 4000
    line = getsh2caskey(sh2cas, "memory")
    if line[0]:
        try:
            qmin["memory"] = int(line[1])

        except ValueError:
            print("PYSCF memory does not evaluate to integer value!")
            sys.exit(1)

    else:
        print(
            "WARNING: Please set memory for PySCF in PYSCF.resources (in MB)! Using 4000 MB default value!"
        )

    os.environ["PYSCF_MAX_MEMORY"] = str(qmin["memory"])

    qmin["ncpu"] = 1
    line = getsh2caskey(sh2cas, "ncpu")
    if line[0]:
        try:
            qmin["ncpu"] = int(line[1])

        except ValueError:
            print("Number of CPUs does not evaluate to integer value!")
            sys.exit(1)

    qmin["delay"] = 0.0
    line = getsh2caskey(sh2cas, "delay")
    if line[0]:
        try:
            qmin["delay"] = float(line[1])

        except ValueError:
            print("Submit delay does not evaluate to numerical value!")
            sys.exit(1)

    line = getsh2caskey(sh2cas, "always_orb_init")
    if line[0]:
        qmin["always_orb_init"] = []

    line = getsh2caskey(sh2cas, "always_guess")
    if line[0]:
        qmin["always_guess"] = []

    if "always_orb_init" in qmin and "always_guess" in qmin:
        print("Keywords 'always_orb_init' and 'always_guess' cannot be used together!")
        sys.exit(1)

    # open template
    with open("PYSCF.template", "r") as f:
        template = f.readlines()

    template_dict = {}
    INTEGERS_KEYS = ["ncas", "nelecas", "roots", "grids-level", "verbose", "max-cycle-macro", "max-cycle-micro", "ah-max-cycle", "ah-start-cycle", "grad-max-cycle", "charge"]
    STRING_KEYS = ["basis", "method", "pdft-functional"]
    FLOAT_KEYS = ["conv-tol", "conv-tol-grad", "max-stepsize", "ah-start-tol", "ah-level-shift", "ah-conv-tol", "ah-lindep", "fix-spin-shift"]
    BOOL_KEYS = []

    template_dict["roots"] = [0 for _ in range(8)]
    template_dict["charge"] = 0

    template_dict["method"] = "casscf"
    template_dict["verbose"] = 3
  

    template_dict["fix-spin-shift"] = 0.2

    template_dict["pdft-functional"] = "tpbe"
    template_dict["grids-level"] = 4
    
    # CASSCF solver defaults
    template_dict["conv-tol"] = 1e-7
    template_dict["conv-tol-grad"] = 1e-4

    template_dict["max-stepsize"] = 0.02
    template_dict["max-cycle-macro"] = 50
    template_dict["max-cycle-micro"] = 4

    template_dict["ah-level-shift"] = 1e-8
    template_dict["ah-conv-tol"] = 1e-12
    template_dict["ah-max-cycle"] = 30
    template_dict["ah-lindep"] = 1e-14
    template_dict["ah-start-tol"] = 2.5
    template_dict["ah-start-cycle"] = 3
    

    # gradient solver defaults
    template_dict["grad-max-cycle"] = 50

    for line in template:
        orig = re.sub("#.*$", "", line).split(None, 1)
        line = re.sub("#.*$", "", line).lower().split()

        if len(line) == 0:
            continue

        key = line[0]
        line = line[1:]

        if key.startswith("spin"):
            template_dict["roots"][int(line[0]) - 1] = int(line[2])

        elif "roots" in key:
            for i, n in enumerate(line):
                template_dict["roots"][i] = int(n)

        elif key in INTEGERS_KEYS:
            template_dict[key] = int(line[0])

        elif key in STRING_KEYS:
            template_dict[key] = line[0]

        elif key in FLOAT_KEYS:
            template_dict[key] = float(line[0])

        elif key in BOOL_KEYS:
            template_dict[key] = True

    # Roots must be larger or equal to states
    for i, n in enumerate(template_dict["roots"]):
        if i == len(qmin["states"]):
            break

        if not n >= qmin["states"][i]:
            print(
                f"Too few states in state-averaging in multiplicity {i+1}! {qmin['states'][i]} requested, but only {n} given."
            )
            sys.exit(1)

    # condense roots list
    for i in range(len(template_dict["roots"]) - 1, 0, -1):
        if template_dict["roots"][i] == 0:
            template_dict["roots"].pop(i)

        else:
            break

    NECESSARY_KEYS = ["basis", "nelecas", "ncas"]
    for key in NECESSARY_KEYS:
        if key not in template_dict:
            print(f"Key {key} missing in template file!")
            sys.exit(1)

    ALLOWED_METHODS = ["casscf", "l-pdft", "mc-pdft", "cms-pdft"]
    for index, method in enumerate(ALLOWED_METHODS):
        if template_dict["method"] == method:
            qmin["method"] = index
            break

    else:
        print(f"Unknown method {template_dict['method']}")
        sys.exit(1)

    # find functional if pdft
    if qmin["method"] == 2 or qmin["method"] == 3 or qmin["method"] == 4:
        ALLOWED_FUNCTIONALS = ["tpbe", "ftpbe"]
        for index, func in enumerate(ALLOWED_FUNCTIONALS):
            if template_dict["pdft-functional"] == func:
                qmin["pdft-functional"] = index
                break

        else:
            print(
                f"Warning! No analytical gradients for L-PDFT with {template_dict['pdft-functional']} given!"
            )
            print(f"Allowed functionals are: {', '.join(ALLOWED_FUNCTIONALS)}")
            sys.exit(1)

        if "nacdr" in qmin:
            print("NACdr not allowed with L-PDFT!")
            sys.exit(1)

    qmin["template"] = template_dict

    # decide which type of gradients to do..
    # 0 = analytical CASSCF gradients in 1 thread/pyscf object (serially)
    # 1 = analytical CASSCF gradients in separate threads/pyscf objects. Possibly distributed over several CPUs (parallel)
    if "grad" in qmin or "nacdr" in qmin:
        if qmin["ncpu"] > 1:
            qmin["gradmode"] = 1

        else:
            qmin["gradmode"] = 0

    else:
        qmin["gradmode"] = 0

    qmin["ncpu"] = max(1, qmin["ncpu"])

    # check the save directory
    if "samestep" in qmin:
        if not os.path.isfile(os.path.join(qmin["savedir"], "pyscf.chk")):
            print("File 'pyscf.chk' missing in SAVEDIR!")
            sys.exit(1)

        if "overlap" in qmin:
            if not os.path.isfile(os.path.join(qmin["savedir"], "pyscf.old.chk")):
                print("File 'pyscf.old.chk' missing in SAVEDIR!")
                sys.exit(1)

    elif "overlap" in qmin:
        if not os.path.isfile(os.path.join(qmin["savedir"], "pyscf.chk")):
            print("File 'pyscf.chk' missing in SAVEDIR")
            sys.exit(1)

    qmin["version"] = get_version()

    if PRINT:
        print_qmin(qmin)

    return qmin


def generate_joblist(qmin):
    """Split the full job into subtasks, each with a qmin dict, a WORKDIR
    structure: joblist = [{WORKDIR: QMin, ..}, {..}, ..]
    each element of the joblist is a est of jobs, and all jobs from the first
    set need to be completed before the second set can be processed."""
    joblist = []
    if qmin["gradmode"] == 0:
        # Serial case on 1 cpu
        qmin_1 = deepcopy(qmin)
        qmin_1["master"] = []
        qmin_1["ncpu"] = 1
        qmin["nslots_pool"] = [1]
        joblist.append({"master": qmin_1})

    elif qmin["gradmode"] == 1:
        # Analytical gradients for several states on several cpus
        # do wave function and dm, soc, overlap always first
        # afterwards do gradients and nacdr asynchronously
        qmin_1 = deepcopy(qmin)
        qmin_1["master"] = []
        qmin_1["gradmap"] = []
        qmin_1["nacmap"] = []
        qmin["nslots_pool"] = [1]
        joblist.append({"master": qmin_1})

        qmin_2 = deepcopy(qmin)
        remove = [
            "h",
            "soc",
            "dm",
            "always_guess",
            "always_orb_init",
            "comment",
            "ncpu",
            "init",
            "veloc",
            "overlap",
            "ion",
            "molden",
        ]
        for r in remove:
            if r in qmin_2:
                del qmin_2[r]

        qmin_2["gradmode"] = 0
        qmin_2["pargrad"] = []
        qmin_2["samestep"] = []
        ntasks = len(qmin["gradmap"]) + len(qmin["nacmap"])

        # Determine number of slots (processes) and number of cpus for each slot
        # right now, we just do this...
        nslots = qmin["ncpu"]
        cpu_per_run = [1] * ntasks

        joblist.append({})
        icount = 0
        for grad in qmin["gradmap"]:
            qmin_3 = deepcopy(qmin_2)
            qmin_3["gradmap"] = [grad]
            qmin_3["nacmap"] = []
            qmin_3["ncpu"] = cpu_per_run[icount]
            icount += 1
            joblist[-1]["grad_%i_%i" % grad] = qmin_3

        for nac in qmin["nacmap"]:
            qmin_3 = deepcopy(qmin_2)
            qmin_3["nacmap"] = [nac]
            qmin_3["gradmap"] = []
            qmin_3["overlap"] = [
                [j + 1, i + 1] for i in range(qmin["nmstates"]) for j in range(i + 1)
            ]
            qmin_3["overlap_nacs"] = []
            qmin_3["ncpu"] = cpu_per_run[icount]
            icount += 1
            joblist[-1]["nacdr_%i_%i_%i_%i" % nac] = qmin_3

        qmin["nslots_pool"].append(nslots)

    if DEBUG:
        pprint.pprint(joblist, depth=3)

    return qmin, joblist


def move_chk_file(qmin):
    """Moves all relevant chk files in the savedir to old-chk files"""
    source = os.path.join(qmin["savedir"], "pyscf.chk")
    target = os.path.join(qmin["savedir"], "pyscf.old.chk")
    if not os.path.isfile(source):
        print(f"File {source} not found, cannot move to old!")
        sys.exit(1)

    if DEBUG:
        print(f"Copy:\t{source}\t==>\t{target}")

    shutil.copy(source, target)


def setup_workdir(qmin):
    """Make the scratch directory, or clean it if it exists. Copy any necessary files."""
    work_dir = qmin["scratchdir"]
    save_dir = qmin["savedir"]

    if os.path.exists(work_dir):
        if not os.path.isdir(work_dir):
            print(f"{work_dir} exists and is not a directory!")
            sys.exit(1)

        else:
            if DEBUG:
                print(f"Remake\t{work_dir}")
            shutil.rmtree(work_dir)
            os.makedirs(work_dir)

    else:
        if DEBUG:
            print(f"Making\t{work_dir}")
        os.makedirs(work_dir)

    source_chk = None
    if "always_guess" not in qmin:
        if "init" in qmin or "always_orb_init" in qmin:
            source_chk = os.path.join(qmin["pwd"], "pyscf.init.chk")

        elif "samestep" in qmin:
            source_chk = os.path.join(save_dir, "pyscf.chk")

        else:
            source_chk = os.path.join(save_dir, "pyscf.old.chk")

    if source_chk is not None and os.path.isfile(source_chk):
        target_chk = os.path.join(work_dir, "pyscf.old.chk")
        if DEBUG:
            print(f"Copying\t{source_chk}\t==>\t{target_chk}", flush=True)
        shutil.copy(source_chk, target_chk)


def save_chk_file(qmin):
    work_dir = qmin["scratchdir"]
    save_dir = qmin["savedir"]

    source_chk = os.path.join(work_dir, "pyscf.chk.master")
    if os.path.isfile(source_chk):
        target_chk = os.path.join(save_dir, "pyscf.chk")
        if DEBUG:
            print(f"Copying\t{source_chk}\t==>\t{target_chk}", flush=True)
        shutil.copy(source_chk, target_chk)


def build_mol(qmin):
    log_file = f"PySCF_{os.path.basename(qmin['scratchdir'])}.log"
    previous_chk = os.path.join(qmin["scratchdir"], "pyscf.old.chk")
    verbose = qmin["template"]["verbose"]

    if os.path.isfile(previous_chk) and "samestep" in qmin:
        if DEBUG:
            print(f"Loading mol from chkfile {previous_chk}", flush=True)
        mol = lib.chkfile.load_mol(previous_chk)
        mol.output = log_file
        mol.verbose = verbose 
        mol.build()

    else:
        mol = gto.Mole(
            atom=qmin["geo"],
            unit="Bohr",
            basis=qmin["template"]["basis"],
            output=log_file,
            verbose=verbose,
            symmetry=False,
            charge=qmin["template"]["charge"]
        )
        mol.build()

    return mol


def gen_solver(mol, qmin):
    mf = mol.RHF()
    mf.max_cycle = 0
    mf.run()

    ncas = qmin["template"]["ncas"]
    nelecas = qmin["template"]["nelecas"]

    if len(qmin["states"]) != 1:
        raise NotImplementedError("Not implemented for states other than singlets!")

    else:
        nroots = qmin["template"]["roots"][0]
        weights = [1.0/nroots,]*nroots

        if qmin["method"] == 0:
            solver = mcscf.CASSCF(mf, ncas, nelecas)

        else:
            functional = qmin["template"]["pdft-functional"]
            grids_level = qmin["template"]["grids-level"]
            try:
                from pyscf import mcpdft
                solver = mcpdft.CASSCF(mf, functional, ncas, nelecas, grids_level=grids_level)

            except ImportError as e:
                print("MC-PDFT requested but pyscf-forge not installed")
                raise e

        try:
            from mrh.my_pyscf.fci import csf_solver
            solver.fcisolver = csf_solver(mol, smult=1)
        
        except ImportError:
            solver.fix_spin_(ss=0, shift=qmin["template"]["fix-spin-shift"])
    
        if qmin["method"] == 1:
            solver = solver.multi_state(weights, method="lin")

        elif qmin["method"] == 3:
            solver = solver.multi_state(weights, method="cms")

        elif qmin["method"] == 2 or qmin["method"] == 0:
            solver = solver.state_average(weights)

        solver.conv_tol = qmin["template"]["conv-tol"]
        solver.conv_tol_grad = qmin["template"]["conv-tol-grad"]
        
        solver.max_stepsize = qmin["template"]["max-stepsize"]
        solver.max_cycle_macro = qmin["template"]["max-cycle-macro"]
        solver.max_cycle_micro = qmin["template"]["max-cycle-micro"]

        solver.ah_level_shift = qmin["template"]["ah-level-shift"]
        solver.ah_conv_tol = qmin["template"]["ah-conv-tol"]
        solver.ah_max_cycle = qmin["template"]["ah-max-cycle"]
        solver.ah_lindep = qmin["template"]["ah-lindep"]
        solver.ah_start_tol = qmin["template"]["ah-start-tol"]
        solver.ah_start_cycle = qmin["template"]["ah-start-cycle"]

    if "master" in qmin:
        solver.chkfile = os.path.join(qmin["scratchdir"], "pyscf.chk.master")
        solver.chk_ci = True

    old_chk = os.path.join(qmin["scratchdir"], "pyscf.old.chk")
    if os.path.isfile(old_chk):
        print(f"Updating solver from chk: {old_chk}", flush=True)
        solver.update(old_chk)
        # This doesn't work with fix_spin so I remove the old CI vector...
        solver.ci = None
        solver.mo_coeff = mcscf.project_init_guess(solver, solver.mo_coeff)

    solver.kernel(solver.mo_coeff)#, solver.ci)
    return solver


def get_dipole_elements(solver):
    mol = solver.mol
    mo_core = solver.mo_coeff[:, : solver.ncore]
    mo_cas = solver.mo_coeff[:, solver.ncore : solver.ncore + solver.ncas]

    nroots = solver.fcisolver.nroots

    dip_matrix = np.ones(shape=(3, nroots, nroots))

    # TODO: decide on gauge???? Use same as OpenMolcas
    # charge_center = (
        # np.einsum("z,zx->x", mol.atom_charges(), mol.atom_coords())
        # / mol.atom_charges().sum()
    # )

    # OpenMolcas I think uses this gauge center...
    gauge_center = (0,0,0)

    dm_core = 2 * mo_core @ mo_core.conj().T

    charges = mol.atom_charges()
    coords = mol.atom_coords()
    coords -= gauge_center
    nucl_term = charges.dot(coords)

    with mol.with_common_origin(gauge_center):
        dipole_ints = mol.intor("int1e_r")

    for state in range(nroots):
        casdm1 = solver.fcisolver._base_class.make_rdm1(
            solver.fcisolver, solver.ci[state], solver.ncas, solver.nelecas
        )
        dm1 = dm_core + mo_cas @ casdm1 @ mo_cas.conj().T
        dip_matrix[:, state, state] = nucl_term -np.einsum(
            "xij,ji->x", dipole_ints, dm1
        )

    for bra in range(nroots):
        for ket in range(bra + 1, nroots):
            t_dm = solver.fcisolver.trans_rdm1(
                solver.ci[bra], solver.ci[ket], solver.ncas, solver.nelecas
            )
            t_dm = mo_cas @ t_dm @ mo_cas.conj().T
            t_dip = -np.einsum("xij, ji->x", dipole_ints, t_dm)
            dip_matrix[:, bra, ket] = t_dip
            dip_matrix[:, ket, bra] = t_dip

    return dip_matrix


def get_grad(solver, qmin):
    grad = []
    err = 0
    zerograd = np.zeros(shape=(qmin["natom"], 3))

    solver_grad = solver.nuc_grad_method()

    solver_grad.max_cycle = qmin["template"]["grad-max-cycle"]

    for i in sorted(qmin["statemap"]):
        mult, state, _ = tuple(qmin["statemap"][i])
        if (mult, state) in qmin["gradmap"]:
            state = state - 1
            de = solver_grad.kernel(state=state)
            if not solver_grad.converged:
                print(f"Gradient failed to converge: {qmin['statemap'][i]}", flush=True)
                err = 1

            grad.append(de)

        else:
            grad.append(zerograd)

    return grad, err


def get_nac(solver, qmin):
    nmstates = qmin["nmstates"]
    nac = np.zeros(shape=(nmstates, nmstates, qmin["natom"], 3))
    err = 0

    solver_nac = solver.nac_method()

    for i in sorted(qmin["statemap"]):
        for j in sorted(qmin["statemap"]):
            m1, s1, ms1 = tuple(qmin["statemap"][i])
            m2, s2, ms2 = tuple(qmin["statemap"][j])
            if m1 != m2:
                continue

            if ms1 != ms2:
                continue

            if s1 == s2:
                continue

            if (m1, s1, m2, s2) in qmin["nacmap"]:
                bra = i - 1
                ket = j - 1
                nacdr = solver_nac.kernel(state=(bra, ket))
                nac[bra][ket] = nacdr
                nac[ket][bra] = -nacdr
                if not solver_nac.converged:
                    print(
                        f"NACDR failed to converge: {qmin['statemap'][i]}, {qmin['statemap'][j]}"
                    )
                    err = 1
    return nac, err


def run_calc(qmin):
    err = 0
    result = {}

    setup_workdir(qmin)
    mol = build_mol(qmin)

    solver = gen_solver(mol, qmin)

    result = {}
    if "h" in qmin:
        if not solver.converged:
            print("Calculator failed to converge!")
            err += 1

        result["energies"] = solver.e_states

    if "dm" in qmin:
        result["dipole"] = get_dipole_elements(solver)

    if "molden" in qmin:
        from pyscf.tools import molden

        molden.from_mcscf(solver, os.path.join(qmin["savedir"], "pyscf.molden"))

    if qmin["gradmap"]:
        result["grad"], e = get_grad(solver, qmin)
        err += e

    if qmin["nacmap"]:
        result["nacdr"], e = get_nac(solver, qmin)
        err += e

    save_chk_file(qmin)

    return err, result


def run_jobs(joblist, qmin):
    """Runs all of the jobs specified in joblist"""
    if "newstep" in qmin:
        move_chk_file(qmin)

    lib.param.TMPDIR = qmin["scratchdir"]
    lib.param.MAX_MEMORY = qmin["memory"]

    print(">>>>>>>>>>>>> Starting the job execution")

    error_codes = {}
    result = {}
    outputs = {}
    for idx, jobset in enumerate(joblist):
        if not jobset:
            continue

        pool = Pool(processes=qmin["nslots_pool"][idx])
        for job in jobset:
            qmin_1 = jobset[job]

            qmin_1["scratchdir"] = os.path.join(qmin_1["scratchdir"], job)
            outputs[job] = pool.apply_async(run_calc, [qmin_1])

            time.sleep(qmin["delay"])

        pool.close()
        pool.join()

        print("")

    for i in outputs:
        error_codes[i], result[i] = outputs[i].get()

    if PRINT:
        string = "  " + "=" * 40 + "\n"
        string += "||" + " " * 40 + "||\n"
        string += "||" + " " * 10 + "All Tasks completed!" + " " * 10 + "||\n"
        string += "||" + " " * 40 + "||\n"
        string += "  " + "=" * 40 + "\n"
        print(string)
        j = 0
        string = "Error Codes:\n\n"
        for i in error_codes:
            string += "\t%s\t%i" % (i + " " * (10 - len(i)), error_codes[i])
            j += 1
            if j == 4:
                j = 0
                string += "\n"
        print(string)

    if any((i != 0 for i in error_codes.values())):
        print("Some subprocesses did not finish successfully!")
        sys.exit(1)

    return result


def combine_result(qmin, result):
    output = {}
    nmstates = qmin["nmstates"]
    natom = qmin["natom"]

    if "h" in qmin:
        output["energies"] = result["master"]["energies"]
    if "dm" in qmin:
        output["dipole"] = result["master"]["dipole"]

    if "grad" in qmin:
        output["grad"] = np.zeros(shape=(nmstates, natom, 3))
        for job in result:
            if "grad" in result[job]:
                output["grad"] += result[job]["grad"]

    if "nacdr" in qmin:
        output["nacdr"] = np.zeros(shape=(nmstates, nmstates, natom, 3))
        for job in result:
            if "nacdr" in result[job]:
                output["nacdr"] += result[job]["nacdr"]

    return output


def write_ham(qmin, result):
    nmstates = qmin["nmstates"]

    string = f"! 1 Hamiltonian Matrix ({nmstates}x{nmstates}, complex)\n"
    string += f"{nmstates} {nmstates}\n"
    for i in range(nmstates):
        for j in range(nmstates):
            if i != j:
                string += f"{eformat(0.0, 9, 3)} {eformat(0.0, 9, 3)} "
            else:
                string += f"{eformat(result['energies'][i].real, 9, 3)} {eformat(result['energies'][i].imag, 9, 3)} "
        string += "\n"
    string += "\n"

    return string


def write_dm(qmin, result):
    nmstates = qmin["nmstates"]
    string = f"! 2 Dipole Moment Matrices (3x{nmstates}x{nmstates}, complex)\n"
    for dipole_xyz in result["dipole"]:
        string += f"{nmstates} {nmstates}\n"
        for bra in dipole_xyz:
            for element in bra:
                string += (
                    f"{eformat(element.real, 9, 3)} {eformat(element.imag, 9, 3)} "
                )

            string += "\n"

    return string


def write_grad(qmin, result):
    states = qmin["states"]
    nmstates = qmin["nmstates"]
    natom = qmin["natom"]
    string = f"! 3 Gradient Vectors ({nmstates}x{natom}x3, real)\n"
    for idx, (imult, istate, ims) in enumerate(itnmstates(states)):
        string += f"{natom} 3 ! {imult} {istate} {ims}\n"
        for atom in result["grad"][idx]:
            for coord in atom:
                string += f"{eformat(coord, 9, 3)} "
            string += "\n"

    return string


def write_nac(qmin, result):
    states = qmin["states"]
    nmstates = qmin["nmstates"]
    natom = qmin["natom"]
    string = f"! 5 Nonadiabatic couplings (ddr) ({nmstates}x{nmstates}x{natom}x3)\n"
    i = 0
    for imult, istate, ims in itnmstates(states):
        j = 0
        for jmult, jstate, jms in itnmstates(states):
            string += f"{natom} {3} ! {imult} {istate} {ims} {jmult} {jstate} {jms}\n"
            for atom in result["nacdr"][i][j]:
                for coord in atom:
                    string += f"{eformat(coord, 12, 3)} "
                string += "\n"

            j += 1
        i += 1

    return string


def write_qmout_time(runtime):
    return f"! 8 Runtime\n{eformat(runtime, 9, 3)}\n"


def write_qmout(qmin, result, qmin_filename):
    """Writes the requested quantities to the file which SHARC reads in. The filename is qmin_filename with everything after the first dot replaced by 'out'."""
    if "." in qmin_filename:
        idx = qmin_filename.find(".")
        outfilename = qmin_filename[:idx] + ".out"

    else:
        outfilename = qmin_filename + ".out"

    if PRINT:
        print(f"===> Writing output to file {outfilename} in SHARC Format\n")

    string = ""
    if "h" in qmin:
        string += write_ham(qmin, result)

    if "dm" in qmin:
        string += write_dm(qmin, result)

    if "grad" in qmin:
        string += write_grad(qmin, result)

    if "nacdr" in qmin:
        string += write_nac(qmin, result)

    string += write_qmout_time(result["runtime"])
    with open(os.path.join(qmin["pwd"], outfilename), "w") as f:
        f.write(string)

    return


def cleanup_directory(dir):
    if PRINT:
        print(f"===> Removing directory {dir}\n")

    try:
        shutil.rmtree(dir)

    except OSError:
        print("fCould not remove directory {dir}")


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
./SHARC_PYSCF.py <qmin>
version: {_version_}
date: {_versiondate_}
changelog: {_change_log_}"""
        )
        sys.exit(1)

    qmin_filename = sys.argv[1]

    print_header()
    qmin = readqmin(qmin_filename)

    qmin, joblist = generate_joblist(qmin)

    result = run_jobs(joblist, qmin)
    result = combine_result(qmin, result)

    runtime = measure_time()
    result["runtime"] = runtime

    write_qmout(qmin, result, qmin_filename)

    if not DEBUG:
        cleanup_directory(qmin["scratchdir"])
        if "cleanup" in qmin:
            cleanup_directory(qmin["savedir"])

    if PRINT or DEBUG:
        print("#================ END ================#")


if __name__ == "__main__":
    main()
