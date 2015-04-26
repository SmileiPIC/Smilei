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
    base="`basename $nml .in`"    
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
# $DIRSMILEI/../scripts/TPUPMC/smileiQt.py ${outDirs}
