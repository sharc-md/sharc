#!/bin/bash

#******************************************
#
#    SHARC Program Suite
#
#    Copyright (c) 2018 University of Vienna
#
#    This file is part of SHARC.
#
#    SHARC is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    SHARC is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    inside the SHARC manual.  If not, see <http://www.gnu.org/licenses/>.
#
#******************************************



if [ ! -f host_info ];
then
  echo 'No host_info file...'
  exit 1
fi

host=$(awk 'NR==1{print $NF}' host_info)
cwd=$(awk 'NR==2{print $NF}' host_info)

echo $host
echo $cwd

if [ "$host" == "" ] || [ "$cwd" == "" ];
then
  echo 'Not sufficient host or path info...'
  exit 1
fi

echo "Copying..."
if [ "$1" == "-lis" ];
then
  scp $USER@$host:$cwd/output.lis .
  exit 0
fi

scp  $USER@$host:$cwd/output.* .

if [ "$1" == "-res" ];
then
  scp  $USER@$host:$cwd/restart.* .
  scp  $USER@$host:$cwd/restart/* ./restart/
fi