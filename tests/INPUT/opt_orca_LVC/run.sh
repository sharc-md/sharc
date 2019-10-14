#/bin/bash

# $PRIMARYDIR and SCRADIR are set
COPY_DIR=$SCRADIR/TRAJ
PRIMARYDIR=$(pwd)
export PATH=$SHARC:$PATH

mkdir -p $COPY_DIR
cp -r $PRIMARYDIR/* $COPY_DIR
cd $COPY_DIR

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$orca
$orca/orca orca.inp > orca.log
err=$?

cp $COPY_DIR/* $PRIMARYDIR
rm -r $COPY_DIR
exit $err








