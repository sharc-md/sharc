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
import pprint

from socket import gethostname


_version_ = "3.0"
_versiondate_ = datetime.date(2024, 2, 2)
_change_log_ = """"""

START_TIME = datetime.datetime.now()

# Global Variables for printing.
DEBUG = True  # Raw output
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
                output += "X "

            else:
                output += ". "

        output += "\n"
        print(output)

    if "nacdr" in qmin:
        output = "Nonadiabatic couplings:\n"
        for i in range(1, qmin["nmstates"] + 1):
            for j in range(1, qmin["nmstates"] + 1):
                if (i, j) in qmin["nacdr"] or (j, i) in qmin["nacrd"]:
                    output += "X "

                else:
                    output += ". "
            output += "\n"
        print(output)

    if "overlap" in qmin:
        output = "Overlaps:\n"
        for i in range(1, qmin["nmstates"] + 1):
            for j in range(1, qmin["nmstates"] + 1):
                if (i, j) in qmin["overlap"] or (j, i) in qmin["overlap"]:
                    output += "X "

                else:
                    output += ". "
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

    min_version = (2, 4, 0)
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

        except:
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
    INTEGERS_KEYS = ["ncas", "nelecas", "roots"]
    STRING_KEYS = ["basis", "method", "pdft-functional"]
    FLOAT_KEYS = []
    BOOL_KEYS = []

    template_dict["roots"] = [0 for i in range(8)]

    template_dict["method"] = "casscf"
    template_dict["pdft-functional"] = "tpbe"

    for line in template:
        orig = re.sub("#.*$", "", line).split(None, 1)
        line = re.sub("#.*$", "", line).lower().split()

        if len(line) == 0:
            continue

        key = line[0]
        line = line[1:]

        if "spin" in key:
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

    ALLOWED_METHODS = ["casscf", "l-pdft"]
    for index, method in enumerate(ALLOWED_METHODS):
        if template_dict["method"] == method:
            qmin["method"] = index
            break

    else:
        print(f"Unknown method {template_dict['method']}")
        sys.exit(1)

    # find functional if pdft
    if qmin["method"] == 2:
        ALLOWED_FUNCTIONALS = ["tpbe"]
        for index, func in enumerate(ALLOWED_FUNCTIONALS):
            if template_dict["pdft-functional"] == func:
                qmin["pdft-functional"] == index
                break

        else:
            print(
                f"Warning! No analytical gradients for L-PDFT with {template_dict['pdft-functional']} given!"
            )
            print(f"Allowed functionals are: {', '.join(ALLOWED_FUNCTIONALS)}")
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


    if PRINT or DEBUG:
        print("#================ END ================#")


if __name__ == "__main__":
    main()
