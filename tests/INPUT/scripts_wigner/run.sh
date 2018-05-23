#!/bin/bash

$SHARC/wigner.py -n 1 -o initconds.ADF ADF.molden
$SHARC/wigner.py -n 1 -o initconds.COLUMBUS COLUMBUS.molden
$SHARC/wigner.py -n 1 -o initconds.GAUSSIAN GAUSSIAN.molden
$SHARC/wigner.py -n 1 -o initconds.MOLCAS MOLCAS.molden
$SHARC/wigner.py -n 1 -o initconds.MOLPRO MOLPRO.molden
$SHARC/wigner.py -n 1 -o initconds.ORCA ORCA.molden

$SHARC/wigner.py -l MOLPRO.molden