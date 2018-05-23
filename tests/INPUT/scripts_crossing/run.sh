#!/bin/bash

if [ -f crossing.xyz ]
then
  rm crossing.xyz
fi

$SHARC/crossing.py < KEYSTROKES.crossing