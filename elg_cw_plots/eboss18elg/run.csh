#!/bin/tcsh -f

#set exec = lf_K.py
#set exec = lf_bJ.py
#set exec = ssfr_m.py
set exec = lf2_cuts.py
#set exec = lf2_2eboss_cuts.py
#set exec = lf2_calz.py

python3 $exec

echo 'The end'

