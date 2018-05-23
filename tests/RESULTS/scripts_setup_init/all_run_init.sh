#/bin/bash

CWD=/user/mai/Documents/NewSHARC/SHARC_2.0/tests_suite/RUNNING_TESTS/scripts_setup_init

cd $CWD/ICOND_00000//
bash run.sh
cd $CWD
echo ICOND_00000/ >> DONE
cd $CWD/ICOND_00001//
bash run.sh
cd $CWD
echo ICOND_00001/ >> DONE
