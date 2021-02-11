import pandas as pd
import pyjoaov as pj
import freesasa as fs
import os
import time
initial_time = time.perf_counter()

#paths
working_path   = "/home/joaov/github/mmpbsa/testing/"
databases_path = "/home/joaov/github/mmpbsa/testing/saved_dataframes/"

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
use_gromacs  = False # if False, freesasa module will be used
pdb_terminus = ['NT3','CT4'] 

# Verbose mode
verbose=True

######################################################################
databases       = [txt_df_res,txt_df_nb,txt_df_charges]
files           = [gro_file,databases,saving_path,new_filepath]
monomers        = [mon1_start,mon1_end,mon2_start,mon2_end]

cutoff_settings = [cutoff,cutoff_n]
settings        = [cutoff_settings,use_gromacs,pdb_terminus,verbose]
######################################################################
pj.mmpbsa_protocol(files,monomers,settings)
######################################################################

final_time = time.perf_counter()
elapsed_time=final_time-initial_time
print("Elapsed time: ", elapsed_time,"seconds.")
#print("Elapsed time: ", elapsed_time/60,"minutes.")

name='01d_gro_en_calc.py, scipy cdist, fsasa'
pj.log_timers(elapsed_time, name)

