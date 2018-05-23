cd QM
$SHARC/SHARC_GAUSSIAN.py QM.in >> QM.log 2>> QM.err
err=$?

exit $err