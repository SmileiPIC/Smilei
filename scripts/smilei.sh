#!/bin/bash

MPIEXEC=mpiexec


SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIRSMILEISH="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
   # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
  [[ $SOURCE != /* ]] && SOURCE="$DIRSMILEISH/$SOURCE"
done
DIRSMILEISH="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

smilei=$DIRSMILEISH/../smilei

if [ "$#" -lt 1 ]; then
    echo "usage: $0 [proc numbers] namelist [namelist] ..."
    exit 1
fi

if [[ $1 = [[:digit:]] ]]; then
    proc=$1
    suffix="_$proc"
    shift
else
    suffix=""
    proc=1
fi

outDirs=""

for nml in $@
do
    base="`basename $nml .py`"    
    dir="${base}${suffix}" 
    outDirs="${outDirs} ${dir}"
    
    rm -rf $dir
    mkdir -p $dir
    cp $nml $dir
    cd $dir
    $MPIEXEC -np $proc $smilei `basename $nml`
    cd ..
done
# echo ${outDirs}
# $DIRSMILEISH/../scripts/TPUPMC/smileiQt.py ${outDirs}
