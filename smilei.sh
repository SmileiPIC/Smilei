#!/bin/bash

MPIEXEC=mpiexec

H=$PWD
smilei=$PWD/src/smilei

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

#outDirs=""

for nml in $@
do
    base="`basename $nml .py`"
    dir=`dirname $nml`/"${base}${suffix}"
    #outDirs="${outDirs} ${dir}"
    
    if [ -d $dir ]; then
    	echo "*   "
    	echo "*   ERROR: Directory $dir already exists. Please remove it first."
    	echo "*   "
    	exit 1
    fi
    #rm -rf $dir
    
    mkdir -p $dir
    cp $nml $dir
    cd $dir
    $MPIEXEC -np $proc $smilei `basename $nml`
    cd $H
done
# echo ${outDirs}
# $DIRSMILEI/../scripts/TPUPMC/smileiQt.py ${outDirs}
