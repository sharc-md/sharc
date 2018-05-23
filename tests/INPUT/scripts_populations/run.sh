#!/bin/bash

if [ -f pop.out ]
then
  rm pop.out
fi

$SHARC/populations.py < KEYSTROKES.populations