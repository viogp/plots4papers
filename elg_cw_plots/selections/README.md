**mass_cum.py, sfr_cum.py, lo2_cum.py** Calculate the cumulative abundance of all galaxies and ELGs vs stellar mass, SFR, L[OII] at different redshifts and store them into files.

**mass_selections.py, sfr_selections.py, lo2_selections.py** Create the environment input for a fixed set of number densities (the same for all the surveys): ascii tables with x,y,z and some other galactic properties. To run these, first run the corresponding ones in the 'find_cuts' folders.

**elg_selections.py** Create the environment input for ELGs: ascii tables with x,y,z and some other galactic properties.

**shuffled_elgs.py** Create the shuffled environment input for ELGs, modifying the elg_selections.py to consider the files output from ../environ/shuffle/shuffleMSI.py

**gsmf_files.py** Generate files with the GSMF (unrelated to the environment study).

Directories:
 * **find_cuts**: with code to find the stellar mass and SFR cuts giving the same number densities.

