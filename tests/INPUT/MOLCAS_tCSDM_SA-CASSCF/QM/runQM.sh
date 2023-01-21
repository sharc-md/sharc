cd QM
$SHARC/SHARC_MOLCAS.py QM.in >> QM.log 2>> QM.err
err=$?

exit $err