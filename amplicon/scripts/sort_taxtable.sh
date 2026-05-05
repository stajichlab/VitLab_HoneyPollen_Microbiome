#!/usr/bin/bash -l
INFILE=$1

if [ -z $INFILE ]; then
	echo "need an input file as an argument"
fi
(head -n 2 $INFILE && tail -n +3 $INFILE | sort -V) > $INFILE.sorted
