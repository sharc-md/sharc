$OVDIR/../bin/wfoverlap.x -f ciov_all.in > ciov_all.out
$OVDIR/../bin/wfoverlap_ascii.x -f ciov_fc.in  > ciov_fc.out
$OVDIR/../bin/wfoverlap.x -f dyson.in > dyson.out

export OMP_NUM_THREADS=4 # If there are more threads, then more memory would be needed!
$OVDIR/../bin/wfoverlap.x -f ciov_fc.in -m 9  > ciov_low-mem.out
$OVDIR/../bin/wfoverlap.x -f ciov_fc_direct.in -m 9  > ciov_low-mem_direct.out
$OVDIR/../bin/wfoverlap.x -f ciov_fc.in -m 8  > ciov_very-low-mem.out

export OMP_NUM_THREADS=2 # If there are more threads, then more memory would be needed!
$OVDIR/../bin/wfoverlap.x -f dyson.in -m 4 > dyson.low-mem.out
$OVDIR/../bin/wfoverlap.x -f ciov_fc.in -m 9  > ciov_very-low-mem2.out
