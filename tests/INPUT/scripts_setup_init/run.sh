#!/bin/bash

if [ -d ICOND_00000 ];
then
  rm -r ICOND_00000
fi
if [ -d ICOND_00001 ];
then
  rm -r ICOND_00001
fi
if [ -f all_run_init.sh ];
then
  rm -r all_run_init.sh
fi

$SHARC/setup_init.py < KEYSTROKES.setup_init