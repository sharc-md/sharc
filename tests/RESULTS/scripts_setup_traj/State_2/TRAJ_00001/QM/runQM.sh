cd QM
$SHARC/SHARC_Analytical.py QM.in >> QM.log 2>> QM.err
err=$?

exit $err