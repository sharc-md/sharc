cd QM
~/Developer/sharc/bin/SHARC_PYSCF.py QM.in >> QM.log 2>> QM.err
err=$?

exit $err
