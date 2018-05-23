PDIR=`pwd`

cd $PDIR/geo1/WORK
${OVDIR}/scripts/read_civfl.py 2 0.998001 >  prep.out
cd $PDIR/geo2/WORK
${OVDIR}/scripts/read_civfl.py 2 0.998001 >> prep.out
cd $PDIR/geo1.6-31g/WORK
${OVDIR}/scripts/read_civfl.py 4 0.99 >> prep.out

cd $PDIR
rm -f ONEINT RUNFILE
${OVDIR}/scripts/seward_double-mol.py geo1 geo2 run >> prep.out
${OVDIR}/../bin/wfoverlap.x > ciov.out

rm -f ONEINT RUNFILE
${OVDIR}/scripts/seward_double-mol.py geo1 geo1.6-31g run >> prep.out
${OVDIR}/../bin/wfoverlap.x -f ciov_proj.in > ciov_proj.out
