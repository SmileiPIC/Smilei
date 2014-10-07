#!/bin/sh 
# this is a helper file for makefile to create dependency files .d

DIR="$1" 
shift 1 
cmd="s@\(^.*\)\.o:@$DIR/\1.d $DIR/\1.o:@"
#  >&2  echo "-------------- $@ | sed -e \"$cmd\""
$@ | sed -e "$cmd"
