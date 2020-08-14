1. Generate input files for CUTE in *generate_input*, using the adequate model:
  - **cute_input_r.py** 
  - **cute_input_gflux_z.py**
  - **cute_input_Lcuts_z.py**
  - **selections_cute_input.py** For the mass/sfr number density selections.

2. Copy the input file to the *run_cute* directory:
   '''> cp [file] ../run_cute/in.txt'''

3. Run CUTE from the *run_cute* directory
   '''> ./qsub_cute.csh  '''

   with the adequate model and z/r-space option.

   This calls run_cute.csh and generates the adequate inpute parameter file within the *params* directory.

4. Code for plots in **r_space_plots** and **z_space_plots**.