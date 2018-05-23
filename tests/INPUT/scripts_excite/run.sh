#!/bin/bash

if [ -f initconds.excited ]
then
  rm initconds.excited
fi
$SHARC/excite.py < KEYSTROKES.excite