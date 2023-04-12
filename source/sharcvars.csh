setenv SHARC=/global/mai/SHARC_3.0_tests/install/sharc/bin
setenv PYSHARC=/global/mai/SHARC_3.0_tests/install/sharc/pysharc
setenv SHARCLIB=/global/mai/SHARC_3.0_tests/install/sharc/lib
setenv ANACONDA=/user/mai/anaconda/miniconda3/envs/pysharc_3.0_test
setenv PYTHONPATH=$SHARCLIB:$PYSHARC:$PYTHONPATH
setenv LD_LIBRARY_PATH=$SHARCLIB:$ANACONDA/lib:$LD_LIBRARY_PATH
