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
saving_path     = working_path + "saved_files/"      # save location for new .pdb from .gro
gro_file        = working_path + "frame00100.gro"    # INPUT FILE
pdb_file        = saving_path  + 'teste_gro.pdb'     # OUTPUT FILE

# number of the atom where the monomers start and end
mon1_start = 1
mon1_end   = 1149
mon2_start = 1150
mon2_end   = 2298

gro_df = pj.read_gro_minus_index(gro_file) #create df with gro info, creates monomer index

monA_range = [mon1_start, mon1_end]
monB_range = [mon2_start, mon2_end]

gro_df = pj.make_index_and_subset(gro_df, monA_range, monB_range)

gro_df['atom'] = 'ATOM  '
gro_df['occ'] = 1.0
gro_df['T_factor'] = 0.0
gro_df['x_coord'] = (gro_df['x_coord'].astype(float)*10)
gro_df['y_coord'] = (gro_df['y_coord'].astype(float)*10)
gro_df['z_coord'] = (gro_df['z_coord'].astype(float)*10)
gro_df['empty'] = ''

gro_df = pd.DataFrame(gro_df,columns=['atom','atom_num','empty','atom_name','empty','res_name','chain_id', 'res_num','empty','empty','x_coord','y_coord','z_coord', 'occ','T_factor'])
np.savetxt(pdb_file,gro_df,delimiter='',fmt=('%6.6s', '%5.5s','%1s', '%-4.4s','%1s', '%-4.4s', '%1.1s','%4.4s','%1s','%3s', '%8.3f', '%8.3f', '%8.3f', '%6.2f', '%6.2f'))
#                                              atom    a_num  empty   a_name alt_i   r_name   chain_id r_num  code_i empty  x_coord  y_coord  z_coord   occ    T_factor

fin  = open(pdb_file, "rt")
data = fin.read() #read file contents to string

to_replace = ['HE11','HE12','HE21','HE22','HH11','HH12','HH21','HH22','HD11','HD12','HD21','HD22']
replace_to = ['1HE1','2HE1','1HE2','2HE2','1HH1','2HH1','1HH2','2HH2','1HD1','2HD1','1HD2','2HD2']

for i, pos in enumerate(to_replace): #replace all occurrences of the required string
  data = data.replace(pos,replace_to[i])
  #print("replaced naming convention:", pos,"-->",replace_to[i])

fin.close() #close the input file

fin = open(pdb_file, "wt") #open the input file in write mode
fin.write(data) #overrite the input file with the resulting data
fin.close() #close the file

final_time = time.perf_counter()
elapsed_time=final_time-initial_time
print("Written new pdb:",pdb_file)
print("Elapsed time: ", round(elapsed_time,3),"seconds.")
