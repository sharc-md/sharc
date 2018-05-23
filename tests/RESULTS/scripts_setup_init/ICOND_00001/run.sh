#!/bin/bash

#$-N init_00001



PRIMARY_DIR=/user/mai/Documents/NewSHARC/SHARC_2.0/tests_suite/RUNNING_TESTS/scripts_setup_init/ICOND_00001//

cd $PRIMARY_DIR


$SHARC/SHARC_Analytical.py QM.in >> QM.log 2>> QM.err
