#!/usr/bin/env python3
import MDAnalysis as mda
import os 

working_path = "/home/joaov/github/mmpbsa/testing/"
origin_file = working_path + "frame00100.gro"
saving_path = working_path + "saved_files/"

file_basename = os.path.basename(origin_file)
file_name, file_extension = os.path.splitext(file_basename)

new_filepath = saving_path + file_name + ".pdb"

universe = mda.Universe(origin_file)

def prYellow(skk):    print("\033[93m {}\033[00m" .format(skk)) 
def prLightGray(skk): print("\033[97m {}\033[00m" .format(skk))

PDB = new_filepath
with mda.Writer(PDB) as pdb:
  pdb.write(universe)
  
prLightGray("New PDB has been written")
  #read input file
fin = open(PDB, "rt")
#read file contents to string
data = fin.read()
#replace all occurrences of the required string
to_replace = ['HE11','HE12','HE21','HE22','HH11','HH12','HH21','HH22','HD11','HD12','HD21','HD22']
replace_to = ['1HE1','2HE1','1HE2','2HE2','1HH1','2HH1','1HH2','2HH2','1HD1','2HD1','1HD2','2HD2']

for i, pos in enumerate(to_replace):
  #print(i,pos)
  data = data.replace(pos,replace_to[i])
  #print("replaced naming convention:", pos,"-->",replace_to[i])
#close the input file
fin.close()
#open the input file in write mode
fin = open(PDB, "wt")
#overrite the input file with the resulting data
fin.write(data)
#close the file
fin.close()
prYellow("New file is: "+new_filepath)


def geom_center(df):
    x = (float(max(df['x_coord']))+float(min(df['x_coord'])))/2
    y = (float(max(df['y_coord']))+float(min(df['y_coord'])))/2
    z = (float(max(df['z_coord']))+float(min(df['z_coord'])))/2
    
    with open('geom_center_xyz.pkl', 'wb') as f:
        pickle.dump([x,y,z], f)
    return x,y,z

center=geom_center(gro_df)


def appropriate_box_size(df):
    x = max(df['x_coord'])-min(df['x_coord'])
    y = max(df['y_coord'])-min(df['y_coord'])
    z = max(df['z_coord'])-min(df['z_coord'])
    max_coor = max(x,y,z)
    return max_coor/0.8
appropriate_box_size(gro_df)