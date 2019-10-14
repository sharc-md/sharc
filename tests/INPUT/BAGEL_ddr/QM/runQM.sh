cd QM
$SHARC/SHARC_BAGEL.py QM.in >> QM.log 2>> QM.err
err=$?

rm *.xml
exit $err