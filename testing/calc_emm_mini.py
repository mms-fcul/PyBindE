import pandas as pd
import pyjoaov as pj
import time
initial_time = time.perf_counter()

#paths
working_path   = "/home/joaov/github/mmpbsa/testing/"
databases_path = "/home/joaov/github/mmpbsa/testing/saved_dataframes/"

txt_df_res      = databases_path + "residues_rtp.txt"
txt_df_nb       = databases_path + "nonbonded_itp.txt"
txt_df_charges  = databases_path + "charges.txt"

gro_file        = working_path + "frame00100.gro" #"short-gro.gro"

# number of the atom where the monomers start and end
mon1_start = 1
mon1_end   = 1149
mon2_start = 1150
mon2_end   = 2298

######################################################################
databases=[txt_df_res,txt_df_nb,txt_df_charges]
files=[gro_file,databases]
monomers=[mon1_start,mon1_end,mon2_start,mon2_end]
######################################################################
VdW_en_total, Coul_en_total = pj.calc_emm(files,monomers)
######################################################################
print("En_VdW = ", VdW_en_total)
print("En_Coul = ", Coul_en_total)
######################################################################
final_time = time.perf_counter()
elapsed_time=final_time-initial_time
print("Elapsed time: ", elapsed_time,"seconds.")
name='calc_emm_mini.py, scipy cdist'
pj.log_timers(elapsed_time, name)
