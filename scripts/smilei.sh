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

if [ "$#" -ne 3 ]; then
    echo "Illegal number of parameters"
    echo "$0 <namelist> <stub> <proc number>"
    exit 1
fi

dir="$2-$3"
rm -rf $dir
mkdir $dir
cp $1 $dir
cd $dir
nml=`basename $1`

$MPIEXEC -np $3 $smilei $nml

#  a=1
#  bind "[" 'a=a-1; set title system(sprintf("sed \"%dq;d\" pippo-1/scalars.txt",a)); repl'
#  bind "]" 'a=a+1; set title system(sprintf("sed \"%dq;d\" pippo-1/scalars.txt",a)); repl'