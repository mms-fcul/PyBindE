import pandas as pd
import pyjoaov as pj
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

#distances cutoff
cutoff = False
cutoff_n = 4.5

#monomers have N and C terminus
terminus = True

######################################################################
gro_df = pj.read_gro_minus_index(gro_file) #create df with gro info, creates monomer index

monA_range = [mon1_start, mon1_end]
monB_range = [mon2_start, mon2_end]

gro_df = pj.make_index_and_subset(gro_df, monA_range, monB_range)

gro_df['atom_name'] = gro_df['atom_name'].replace('O2', 'OT2')
gro_df['atom_name'] = gro_df['atom_name'].replace('O1', 'OT1')

#Import databases
df_res      = pd.read_table(txt_df_res,     delim_whitespace=True,header=None,names=['res_name','atom_name','atom_type']) # df w/ atom names/types from rtp
df_charges  = pd.read_table(txt_df_charges, delim_whitespace=True,header=None,names=['atom_name','res_name','charge'])
df_nb       = pd.read_table(txt_df_nb,      delim_whitespace=True,header=None,names=['i','j','c6','c12']) # df from nb.itp to grab c6 and c12

#print(df_nb[df_nb.duplicated()])
df_nb = df_nb.rename(columns = {"i": "type_Ai", "j": "type_Bj"})
#df_nb = df_nb.drop_duplicates()

#Change C and N terminus residue to CTR and NTR
if terminus:
  
  ntr_list=[]
  for i in df_res.loc[(df_res["res_name"] == 'NTR'), "atom_name"]:
      ntr_list.append(i)

  ctr_list=[]
  for i in df_res.loc[(df_res["res_name"] == 'CTR'), "atom_name"]:
      ctr_list.append(i)

  pj.fix_terminus_improved(gro_df,ntr_list,ctr_list)

#add information about atom_type
gro_res_df = pd.merge(gro_df, df_res,  how = 'left', on = ['res_name', 'atom_name'])
gro_res_df['atom_type'] = gro_res_df['atom_type'].str.replace('CH2R', 'CH2r')

# add charges do gro df
gro_res_crg_df = pd.merge(gro_res_df, df_charges,  how = 'left', on = ['res_name', 'atom_name'])

# see if there's empty rows in atom_type
nans = gro_res_crg_df.query("atom_type != atom_type")
if not nans.empty:
  print(nans)
  pj.prYellow("There are atoms with no attribution in atom_type!")

monomers_df = gro_res_crg_df

monA = monomers_df.query('chain_id == "A"')
trimmed_monA = pd.DataFrame(monA,columns=['atom_num','atom_type','charge','x_coord','y_coord','z_coord'])

monB = monomers_df.query('chain_id == "B"')
trimmed_monB = pd.DataFrame(monB,columns=['atom_num','atom_type','charge','x_coord','y_coord','z_coord'])

# calculate distances between 2 mons
#print("Calculating distances...")
dist_df = pj.calc_cdist(trimmed_monA, trimmed_monB)

# fix some names between gro and databases
complete_df = pd.merge(dist_df, df_nb,  how = 'left', on = ['type_Ai', 'type_Bj'])
complete_df['type_Ai'] = complete_df['type_Ai'].str.replace('CH2R', 'CH2r')
complete_df['type_Bj'] = complete_df['type_Bj'].str.replace('CH2R', 'CH2r')


#Vdw and coulomb calculations
complete_df['En_VdW'] = (complete_df['c12']/complete_df['distance']**12)-(complete_df['c6']/complete_df['distance']**6)

f = 138.935458
Er = 1
complete_df['En_Coul'] = f*(complete_df['charge_Ai']*complete_df['charge_Bj'])/(Er*complete_df['distance'])


VdW_en_total = complete_df['En_VdW'].sum()
print("En_VdW = ", VdW_en_total)
Coul_en_total = complete_df['En_Coul'].sum()
print("En_Coul = ", Coul_en_total)


final_time = time.perf_counter()
elapsed_time=final_time-initial_time
print("Elapsed time: ", elapsed_time,"seconds.")
name='calc_emm.py, scipy cdist'
pj.log_timers(elapsed_time, name)
