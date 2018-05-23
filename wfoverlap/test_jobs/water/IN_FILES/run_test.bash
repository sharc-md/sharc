export OMP_NUM_THREADS=1
$OVDIR/../bin/wfoverlap.x -f ciovl.in > ciovl.out.1thr || exit $?
$OVDIR/../bin/wfoverlap.x -f ciovl.in.lumorb > ciovl.out.lumorb || exit $?
export OMP_NUM_THREADS=2
$OVDIR/../bin/wfoverlap.x -f ciovl.in > ciovl.out.2thr || exit $?
$OVDIR/../bin/wfoverlap.x -f ciovl.in.direct > ciovl.out.direct || exit $?
