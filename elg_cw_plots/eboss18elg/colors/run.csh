#!/bin/tcsh -f

set exec = decam.py
#set exec = decam_eboss.py
#set exec = decam_desi.py

python $exec

echo 'The end'

