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
set nargs = $#argv
#echo "number of script arguments = $nargs"
if( $nargs == 1 )then
    set parfil = $1
    echo 'Parameter file =' $parfil
else
    echo 'ERROR - USAGE: CUTE paramfile' 
    exit 1
endif

$exec $parfil


