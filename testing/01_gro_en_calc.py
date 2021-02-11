import pandas as pd
import pyjoaov as pj
import freesasa as fs
import MDAnalysis as mda
import os
import time
initial_time = time.perf_counter()
#paths
working_path   = "/home/joaov/github/mmpbsa/testing/"
databases_path = "/home/joaov/github/mmpbsa/testing/saved_dataframes/"

saving_path     = working_path + "saved_files/" #save location for new .pdb from .gro

gro_file        = working_path + "frame00100.gro" #"short-gro.gro"
#pdb_file        = saving_path + "frame00100.pdb"

txt_df_res      = databases_path + "residues_rtp.txt"
txt_df_nb       = databases_path + "nonbonded_itp.txt"
txt_df_charges  = databases_path + "charges.txt"

# number of the atom where the monomers start and end
mon1_start = 1
mon1_end   = 1149
mon2_start = 1150
mon2_end   = 2298

#distances cutoff
cutoff = False
cutoff_n = 4.5

######################################################################
gro_df = pj.read_gro(gro_file, mon1_start, mon1_end, mon2_start, mon2_end) #create df with gro info, creates monomer index

monA = gro_df.query('chain_id == "A"') #creates a dataframe with only monA
monB = gro_df.query('chain_id == "B"') #creates a dataframe with only monB

df_res = pd.read_table(txt_df_res,delim_whitespace=True,header=None,names=['res_name','atom_name','atom_type']) # df w/ atom names/types from rtp

mega_monA = pj.terminus_fix(monA, df_res) # grabs and fixes terminus, CTR and NTR added
mega_monB = pj.terminus_fix(monB, df_res) # grabs and fixes terminus, CTR and NTR added

merged_df = pd.merge(gro_df, df_res,  how = 'left', on = ['res_name', 'atom_name'])

merged_df.update(mega_monA) # use previously fixed terminus to update merged_df
merged_df.update(mega_monB) # use previously fixed terminus to update merged_df

charges = pd.read_table(txt_df_charges,delim_whitespace=True,header=None,names=['atom_name','res_name','charge'])

# add charges do gro df
complete_df = pd.merge(merged_df, charges,  how = 'left', on = ['res_name', 'atom_name'])

#nan_complete_df = complete_df.query("charge != charge")

nb_df = pd.read_table(txt_df_nb,delim_whitespace=True,header=None,names=['i','j','c6','c12']) # df from nb.itp to grab c6 and c12
complete_nb_df = nb_df

complete_nb_df = complete_nb_df.rename(columns = {"i": "type_Ai", "j": "type_Bj"})
complete_nb_df = complete_nb_df.drop_duplicates()

complete_monA = complete_df.query('chain_id == "A"') 
complete_monB = complete_df.query('chain_id == "B"') 
#print(complete_df.head(20))
#print(complete_df.tail(20))
#pj.df_snapshot(complete_df,"complete_df")
# calculate distances between 2 mons
dist_df = pd.DataFrame(pj.calc_dists(complete_monA, complete_monB, cutoff, cutoff_n))

# fix some names between gro and databases
complete_nb_df['type_Ai'] = complete_nb_df['type_Ai'].str.replace('CH2R', 'CH2r')
complete_nb_df['type_Bj'] = complete_nb_df['type_Bj'].str.replace('CH2R', 'CH2r')
dist_df['type_Ai'] = dist_df['type_Ai'].str.replace('CH2R', 'CH2r')
dist_df['type_Bj'] = dist_df['type_Bj'].str.replace('CH2R', 'CH2r')

copy_nb = complete_nb_df
copy_dist = dist_df

merged_copies = pd.merge(copy_dist, copy_nb,  how = 'left', on = ['type_Ai', 'type_Bj'])

nan_merged_copies = merged_copies.query("charge_Bj != charge_Bj")

merged_copies['En_VdW'] = (merged_copies['c12']/merged_copies['distance']**12)-(merged_copies['c6']/merged_copies['distance']**6)

f = 138.935458
Er = 1
merged_copies['En_Coul'] = f*(merged_copies['charge_Ai']*merged_copies['charge_Bj'])/(Er*merged_copies['distance'])

merged_copies = merged_copies.drop_duplicates(subset = ['atom_Ai', 'atom_Bj', 'distance'])

VdW_en_total = merged_copies['En_VdW'].sum()
print("En_VdW = ", VdW_en_total)
Coul_en_total = merged_copies['En_Coul'].sum()
print("En_Coul = ", Coul_en_total)
#pj.df_snapshot(merged_copies,"complete_df_01")
#convert gro to pdb

file_basename = os.path.basename(gro_file)
file_name, file_extension = os.path.splitext(file_basename)

new_filepath = saving_path + file_name + ".pdb"

universe = mda.Universe(gro_file)

PDB = new_filepath
with mda.Writer(PDB) as pdb:
  pdb.write(universe)
  
pj.prLightGray("New PDB has been written")

fin = open(PDB, "rt") #read input file
data = fin.read() #read file contents to string

#replace all occurrences of the required string
to_replace = ['HE11','HE12','HE21','HE22','HH11','HH12','HH21','HH22','HD11','HD12','HD21','HD22']
replace_to = ['1HE1','2HE1','1HE2','2HE2','1HH1','2HH1','1HH2','2HH2','1HD1','2HD1','1HD2','2HD2']

for i, pos in enumerate(to_replace):
  #print(i,pos)
  data = data.replace(pos,replace_to[i])
  #print("replaced naming convention:", pos,"-->",replace_to[i])

fin.close() #close the input file

fin = open(PDB, "wt") #open the input file in write mode
fin.write(data) #overrite the input file with the resulting data
fin.close()#close the file

pj.prYellow("New file is: "+new_filepath)



#calculate SASA estimate

use_gromacs = True

if use_gromacs:
  
  gmx20_path  = "/usr/bin/gmx"
  file_path   = PDB
  tpr_path    = "/home/joaov/github/mmpbsa/gromacs-dependent/files/sas.tpr"
  ndx_path    = "/home/joaov/github/mmpbsa/gromacs-dependent/files/pb.ndx"
  
  sasa    = pj.run_g_sasa(gmx20_path,file_path,tpr_path,ndx_path)

  print(sasa)
else:
  fs.setVerbosity(1)
  structure = fs.Structure(PDB)
  result =fs.calc(structure,fs.Parameters({'algorithm' : fs.ShrakeRupley,
                                           'n-points' : 10000}))

  #print("Total : %.2f A2" % result.totalArea())
  selections = fs.selectArea(('dimer, resi 1-198','monA, resi 1-99', 'monB, resi 100-198'), structure, result)
  list_keys=[]
  for key in selections:
      print (key, ": %.2f A2" % selections[key])
      list_keys.append(selections[key])
  sasa_P  = list_keys[0]
  sasa_MA = list_keys[1]
  sasa_MB = list_keys[2]
  sasa      = [sasa_P,sasa_MA,sasa_MB]
  
  
pdb_coord = pj.read_pdb(PDB)
pdb_df = pd.DataFrame(pdb_coord)

gcenter=pj.geom_center(pdb_df)
print("Geometric center of dimer in pdb is:",gcenter)

box_size=pj.appropriate_box_size(pdb_df)
print("Appropriate box side size of pdb is:",box_size)
  


#generate energy summary file energies.txt
g = 0.00542 # gamma -> kcal/(mol‚A2)
b = 0.92    # beta  -> kcal/mol

nonbonded = [VdW_en_total,Coul_en_total]
polar     = [0,0,0]
#polar     = [solv_dimer, solv_mon1, solv_mon2]
UserNote= "MISSING POLAR ENERGY VALUES"
pj.gen_en_summary(gro_file,g,b,nonbonded,sasa,polar,UserNote)
pj.prYellow("New energy summary has been created!")

final_time = time.perf_counter()
#final_time = (final_time.strftime('%d-%m-%Y %H:%M:%S'))

elapsed_time=final_time-initial_time
print("Elapsed time: ", elapsed_time,"seconds.")
print("Elapsed time: ", elapsed_time/60,"minutes.")

name='01_gro_en_calc.py'
pj.log_timers(elapsed_time, name)