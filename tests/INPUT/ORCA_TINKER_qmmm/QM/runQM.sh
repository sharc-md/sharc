cd QM
$SHARC/SHARC_ORCA.py QM.in >> QM.log 2>> QM.err
err=$?
exit $err