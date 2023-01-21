cd QM


echo -n "hello 1  " >> QM.err
date >> QM.err
$SHARC/SHARC_MOLCAS.py QM.in >> QM.log 2>> QM.err
echo -n  "hello 2" >> QM.err
date >>QM.err
err=$?

exit $err
