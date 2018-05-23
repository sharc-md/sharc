PDIR=`pwd`

cd $PDIR/GEO1/WORK
${OVDIR}/scripts/read_civfl.py 2  >  prep.out
cd $PDIR/GEO2/WORK
${OVDIR}/scripts/read_civfl.py 2  >> prep.out

cd $PDIR
rm -f ONEINT RUNFILE
${OVDIR}/scripts/seward_double-mol.py GEO1 GEO2 run >> prep.out
${OVDIR}/../bin/wfoverlap.x -f ciov_bin.in > ciov_bin.out
