#!/bin/bash

if [ -d State_2 ];
then
  rm -r State_2
fi

$SHARC/setup_traj.py < KEYSTROKES.setup_traj