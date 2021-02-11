#!/usr/bin/python3
import pandas as pd
import numpy as np
import time
import sys
sys.path.insert(0, '/home/joaov/github/mmpbsa/testing')
import pyjoaov as pj

initial_time = time.perf_counter()

#paths
working_path    = "/home/joaov/github/mmpbsa/testing/"
databases_path = "/home/joaov/github/mmpbsa/testing/saved_dataframes/"
saving_path     = working_path + "saved_files/"      # save location for new .pdb from .gro
gro_file        = working_path + "frame00100.gro"    # INPUT FILE
pdb_file        = saving_path  + 'corrected_gro.pdb'     # OUTPUT FILE
txt_df_res      = databases_path + "residues_rtp.txt"
# number of the atom where the monomers start and end
mon1_start = 1
mon1_end   = 1149
mon2_start = 1150
mon2_end   = 2298

gro_df = pj.read_gro_minus_index(gro_file) #create df with gro info, creates monomer index

monA_range = [mon1_start, mon1_end]
monB_range = [mon2_start, mon2_end]

gro_df = pj.make_index_and_subset(gro_df, monA_range, monB_range)

gro_df['atom_name'] = gro_df['atom_name'].replace('O2', 'OT2')
gro_df['atom_name'] = gro_df['atom_name'].replace('O1', 'OT1')

#Import databases
df_res      = pd.read_table(txt_df_res,     delim_whitespace=True,header=None,names=['res_name','atom_name','atom_type']) # df w/ atom names/types from rtp

ntr_list=[]
for i in df_res.loc[(df_res["res_name"] == 'NTR'), "atom_name"]:
    ntr_list.append(i)

ctr_list=[]
for i in df_res.loc[(df_res["res_name"] == 'CTR'), "atom_name"]:
    ctr_list.append(i)
    
#terminus = ['NTR','CTR']
pj.fix_terminus_improved(gro_df,ntr_list,ctr_list)

pdb_terminus=['NT3','CT4']
pj.gro2pdb_simple(gro_df, pdb_file, pdb_terminus)

final_time = time.perf_counter()
elapsed_time=final_time-initial_time
print("Written new pdb:",pdb_file)
print("Elapsed time: ", round(elapsed_time,3),"seconds.")
