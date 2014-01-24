#!/bin/bash
procnum="1"
mpiexec="openmpirun"
outdir="."


args=`getopt n:m:d:w:h $*`
for i
do
	case "$i"
	in
		-h)
			echo "help"
			exit 2;
			shift;; 
		-n)
			procnum=$2; shift;
			shift;;
		-m)
			mpiexec=$2; shift;
			shift;;
		-d)
			outdir=$2; shift;
			shift;;
		--)
			shift; break;;
	esac
done

#mkdir -p "$outdir"

basename $@

echo $@