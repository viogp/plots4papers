#!/bin/tcsh -ef

#set exec = o2_fsat.py
#set exec = mass_fsat.py
#set exec = sfr_fsat.py
set exec = ssfr_fsat.py
#set exec = g_fixedlim.py
#set exec = g_lowlim.py
#set exec = g_highlim.py

python $exec

echo 'The end'

