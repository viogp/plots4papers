#!/bin/tcsh

# ifort reformat125.f90 -L/gpfs/data/violeta/f90lib  -lvio
set exec = /cosma6/data/dp004/dc-gonz3/CUTE/CUTE_box/CUTE_box

limit datasize unlimited
limit stacksize unlimited
limit coredumpsize 0k
limit vmemoryuse 5000m
###########################################################################
# Read in INPUT ARGUMENTS   ~/lines/desi_hod_o2/o2/cuts/xi
###########################################################################

echo "Submitting job, param_file =" $parfil

$exec $parfil


