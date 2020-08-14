#!/bin/tcsh -ef

#set exec = hod.py
set exec = medians.py
#set exec = sfr_m.py

python3 $exec

echo 'The end'

