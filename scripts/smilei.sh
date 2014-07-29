#!/bin/bash

MPIEXEC=mpiexec

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIRSMILEI="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
   # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
  [[ $SOURCE != /* ]] && SOURCE="$DIRSMILEI/$SOURCE"
done
DIRSMILEI="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

smilei=$DIRSMILEI/../src/smilei

if [ "$#" -le 1 ]; then
    echo "Illegal number of parameters"
    echo "$0 <namelist> <proc number>"
    exit 1
fi


base=`basename $1 .in`


for i in ${@:2}
do 
    dir="$base-$i"
    rm -rf $dir
    mkdir -p $dir
    cp $1 $dir
    cd $dir
    $MPIEXEC -np $i $smilei $1
    cd ..
done
