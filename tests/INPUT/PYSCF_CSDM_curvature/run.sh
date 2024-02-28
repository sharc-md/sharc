#/bin/bash

# $PRIMARYDIR and SCRADIR are set
COPY_DIR=$SCRADIR/TRAJ
PRIMARYDIR=$(pwd)

mkdir -p $COPY_DIR
cp -r $PRIMARYDIR/* $COPY_DIR
cd $COPY_DIR

$SHARC/sharc.x input
err=$?

cp $COPY_DIR/output.* $PRIMARYDIR
if [ ! $err == 0 ];
then
  cp $COPY_DIR/QM/* $PRIMARYDIR/QM/
fi
rm -r $COPY_DIR
exit $err
