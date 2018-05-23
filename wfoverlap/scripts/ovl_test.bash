#!/bin/bash

#******************************************
#
#    SHARC Program Suite
#
#    Copyright (c) 2018 University of Vienna
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

#!/bin/bash

# Test of cioverlap programs.
# The program checks whether the COLUMBUS and MOLCAS variables are set
#     and chooses the tests accordingly.

export OVDIR=${1}

if [ -z "$OVDIR" ]
then
    echo "Please give the main directory as argument:"
    echo " $0 <OVDIR>"
    exit 1
fi

TDIR="$OVDIR/test_jobs"

echo "Starting ovl_test.bash ..."
echo "TDIR=$TDIR"
echo

tests="water CH2_triplet CH2_dyson Adenine_fc CH2S_ricc2"
colt="CH2_dalton"
colmct="CH2_doublet CH2_CISD" # water_molcas - requires Molcas dev. version

if [ ! -z "$COLUMBUS" ]; then
    echo "found COLUMBUS=$COLUMBUS"
    echo "  including Tests: $colt"
    tests="$tests $colt"
    if [ ! -z "$MOLCAS" ]; then
        echo "found MOLCAS=$MOLCAS"
        echo "  including Tests: $colmct"
        tests="$tests $colmct"
    fi
fi

rm -r OVL_TEST
mkdir OVL_TEST
cd    OVL_TEST
PDIR=`pwd`

for dir in $tests
do
    echo
    echo "================================================"
    echo
    echo "Starting test $dir ..."
    sdir="$TDIR/$dir"

    rdir="$PDIR/$dir"
    if [ -d $rdir ]
    then
        echo " ERROR:"
        echo "  $rdir already exists!"
        echo "  Please delete it or run in a different directory."
        exit 5
    fi

    cp -r "$sdir/IN_FILES" $rdir
    cd $rdir

    bash ./run_test.bash 2> run_test.err
    chk=$?
    if [ "$chk" -ne 0 ]; then
        echo "  ... failed!"
    else
        echo
        echo "Checking output files:"
        for rfile in `ls "$sdir/REF_FILES"`
        do
            echo "  -> $rfile"
            diff -b -I walltime -I thread -I "memory limit" -I ^MO -I "P_ovl matrix" "$sdir/REF_FILES/$rfile" $rfile
            chk=$((chk+$?))
        done
    fi

    echo
    echo " *** Test $dir finished (error code: $chk)."
    tchk=$((tchk+chk))
done

echo
echo "================================================"
echo
echo " *** All tests finished (number of errors: $tchk)"
echo
echo "================================================"

exit $tchk
