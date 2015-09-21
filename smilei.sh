#!/bin/bash

MPIEXEC=mpirun

H=$PWD # current dir
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ) # dir of this script
smilei=$DIR/smilei # path to the smilei executable
script=$0 # name of this script

# Function to which check that there are arguments remaining
check () {
  if [ "$#" -lt 1 ]; then
    echo "usage: $script [proc numbers] namelist [outdir]"
    exit 1
  fi
}

# Check that at least one argument provided
check $@

# If first argument is a number, this is the number of procs
if [[ $1 =~ ^[0-9]+$ ]]; then
    proc=$1
    shift
else
    proc=1
fi

# Check there are arguments remaining
check $@

# Next argument is the path to the namelist
nml=$1
base="`basename $nml .py`" # strip directory and extension
shift

# Next argument, if exists, is the outdir 
outdir=`dirname $nml`/$base # by default, the outdir is same as namelist
if [ "$#" -gt 0 ]; then
outdir=`dirname $nml`/$base # by default, the outdir is same as namelist
  outdir=$1/$base # otherwise, provided as argument
  shift
fi

# If outdir already exist, error
if [ -d $outdir ]; then
  echo "*   "
  echo "*   ERROR: Directory $outdir already exists. Please remove it first."
  echo "*   "
  exit 1
fi

# Make new directory, go there, and run
mkdir -p $outdir
cp $nml $outdir
cd $outdir
$MPIEXEC -mca btl tcp,sm,self -np $proc $smilei `basename $nml`
cd $H

