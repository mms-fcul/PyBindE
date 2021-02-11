import pandas as pd
import pyjoaov as pj
import freesasa as fs
import os
import time
initial_time = time.perf_counter()

#paths
working_path   = "/home/joaov/github/mmpbsa/test_system"+"/"
databases_path = "/home/joaov/github/mmpbsa/test_system/saved_dataframes/"

txt_df_res      = databases_path + "residues_rtp.txt"
txt_df_nb       = databases_path + "nonbonded_itp.txt"
txt_df_charges  = databases_path + "charges.txt"

gro_file        = working_path + "frame00100.gro" #"short-gro.gro"

saving_path     = working_path + "saved_files/" #save location for new .pdb from .gro

file_basename   = os.path.basename(gro_file)
file_name, file_extension = os.path.splitext(file_basename)
new_filepath    = saving_path + file_name + ".pdb"

# number of the atom where the monomers start and end
mon1_start = 1
mon1_end   = 1149
mon2_start = 1150
mon2_end   = 2298

# Cutoff settings for distances
cutoff   = False
cutoff_n = 4.5

# Settings for sasa calculations
use_gromacs  = False
pdb_terminus = ['NT3','CT4'] 

# Verbose mode
verbose=False

##############################################################################
databases=[txt_df_res,txt_df_nb,txt_df_charges]
files=[gro_file,databases,saving_path,new_filepath]
monomers=[mon1_start,mon1_end,mon2_start,mon2_end]

cutoff_settings = [cutoff,cutoff_n]
settings        = [cutoff_settings,use_gromacs,pdb_terminus,verbose]
##############################################################################
pj.mmpbsa_protocol(files,monomers,settings)
#################################################################################
## --This protocol reads the gro_file and all databases, matches the required  ##
## info in databases to the atoms and residues on the gro_file, calculates     ##
## distances between the atoms of the 2 monomers, calculates VdW and Coulomb   ## 
## energies (EMM).------------------------------------------------------------ ##
## --After that, calculates the sasa for the dimer, mon1 and  mon2 that will   ##
## be used to calculate Gnonpolar. For this step, gromacs OR freesasa (with    ##
## a convertion from gro to pdb) may be used.--------------------------------- ##
## --Next, it sets up a few variables needed for the Gpolar calculation,       ##
## using Delphi4py, and generates a summary file for all the energy terms      ##
## that contribute to Gbinding.----------------------------------------------- ##
## --The Gpolar energy for the test gro_file was calculated using Delphi4py    ##
## and manually added to the energy summary file upon generation               ##                            
#################################################################################

final_time = time.perf_counter()
elapsed_time=final_time-initial_time
print("Elapsed time: ", elapsed_time,"seconds.")
print("Elapsed time: ", elapsed_time/60,"minutes.")

name='01d_gro_en_calc.py, scipy cdist, gsas'
pj.log_timers(elapsed_time, name)
