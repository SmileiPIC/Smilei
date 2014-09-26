#!/bin/bash

MPIEXEC=mpiexec

oldDir=`pwd`

isnumber() { test "$1" && printf '%d' "$1" >/dev/null 2>&1; }


SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIRSMILEI="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
   # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
  [[ $SOURCE != /* ]] && SOURCE="$DIRSMILEI/$SOURCE"
done
DIRSMILEI="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

smilei=$DIRSMILEI/../../src/smilei

if [ "$#" -lt 1 ]; then
    echo "Illegal number of parameters"
    echo "$0 <namelist>"
    exit 1
fi

base=`basename $1 .in`
rm -rf $base
mkdir -p $base
cp $1 $base
cd $base
$smilei $1

python $DIRSMILEI/smileiQt.py

cd $oldDir
