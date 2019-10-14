cd QM
$SHARC/SHARC_LVC.py QM.in >> QM.log 2>> QM.err
err=$?

exit $err