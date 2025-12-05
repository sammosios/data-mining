#!/bin/bash

#check for installation of parallel-rsync

if [ "$#" -ne 1 ] ; then
	echo "Data file not supplied."
	echo "Usage ./plot {data-file.txt}"
	exit
fi

gnuplot -e "filename='$1'" graph.gnuplot

if command -v xdg-open >/dev/null 2>&1; then
    xdg-open graph.png || open graph.png
else
    open graph.png
fi