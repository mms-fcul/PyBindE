#!/usr/bin/python3

from os import read
import re
import pandas as pd
from biopandas.pdb import PandasPdb
import numpy as np
import subprocess
import matplotlib.pyplot as plt
pdb_file = "/home/joaov/python-mmpbsa/mmpbsa/testing/dimer.pdb"
# %%
res_conv = "/home/joaov/python-mmpbsa/mmpbsa/testing/res_atoms_db.rtp"

def read_pdb(file):
    pdb_lines = []
    pdb_dict = []
    with open(file) as pdbfile:
        for line in pdbfile:
            if line[:4] == 'ATOM' or line[:6] == "HETATM":
                splitted_line = [line[:6], line[6:11], line[12:16], line[17:20], 
                                 line[21], line[22:26], line[30:38], line[38:46], line[46:54]]
                if splitted_line[3] != 'SOL':
                    pdb_lines.append(splitted_line)
                else:
                    continue
                #print (splitted_line)
        for line in pdb_lines:
            line_dict = {
                'atom_num':    line[1].strip(), 
                'atom_name':   line[2].strip(),
                'res_name':    line[3].strip(), 
                'chain_id':    line[4].strip(), 
                'x_coord':    float(line[6].strip()), 
                'y_coord':    float(line[7].strip()), 
                'z_coord':    float(line[8].strip())
            }
            pdb_dict.append(line_dict)
    return pdb_dict


pdb_coord = read_pdb(pdb_file)
# print(pdb_coord)
pdb_df = pd.DataFrame(pdb_coord)
# print(df)
#treat df
to_replace = ['1HE1','2HE1','1HE2','2HE2','1HH1','2HH1','1HH2','2HH2','1HD1','2HD1','1HD2','2HD2']
replace_to = ['HE11','HE12','HE21','HE22','HH11','HH12','HH21','HH22','HD11','HD12','HD21','HD22']
df = pdb_df.replace(to_replace, replace_to)

def read_resdb(file):
    residues = []
    with open(file) as res_db:
        for line in res_db:
            if line.startswith('['):
                res = str(line).rstrip('\n')
                res1 = re.sub(r'\[','',res).strip()
                res2 = re.sub(r'\]', '', res1).strip()

            elif not line.startswith(';'):
                cleanline = re.sub('\s+', ' ', line).strip()
                fields = cleanline.split(' ')
                res_dict = {
                    'res':              res2,
                    'atom_name':        fields[0].rstrip(),
                    'atom_type':        fields[1].rstrip(),
                    'charge':     float(fields[2])
                    }
                residues.append(res_dict)
    return residues

df_res=pd.DataFrame(read_resdb(res_conv))

merged_df = pd.merge(df, df_res,  how='left', left_on=['res_name','atom_name'], right_on = ['res_name','atom_name'])

df=merged_df
subset_df= df[['atom_num','atom_name','atom_type','res_name','chain_id']]
# %%
# def make_index(PDB_dictionary):
monomers = []
index = {}
for atom in pdb_coord:
    monomer = atom['chain_id']
    if monomer not in monomers:
        monomers.append(monomer)
for n in monomers:
    index[n] = []
    for atom in pdb_coord:
        atom_num = atom['atom_num']
        monomer = atom['chain_id']
        if monomer == n:
            index[n].append(int(atom_num))
        else:
            continue
  #return index

# monomer_index= make_index(pdb_coord)
# print(monomer_index)
  
# ppdb=PandasPdb().read_pdb(pdb_file)  

#df = pd.DataFrame(pdb_coord)
# df=df.set_index('atom_num')
  
monA = df.query('chain_id == "A"')
monB = df.query('chain_id == "B"')
# print(monA)
# print(monB)

log_dist = "/home/joaov/python-mmpbsa/mmpbsa/testing/pdb_dist.log"
c = subprocess.Popen(["/bin/rm", '-f', log_dist])
c.communicate(b'\n')

cuttoff = 10
#with open(log_dist, "a") as log:
    #log.write(('{0:5} {1:5} {2:5}'.format("num_a", "num_b", "distance"+"\n")))
def calc_dists():
    distances_db = []
    for a in monA.index:
        xa = monA['x_coord'][a]
        ya = monA['y_coord'][a]
        za = monA['z_coord'][a]
        num_a = monA['atom_num'][a]
        a_coord = np.array((xa, ya, za), dtype = float)
        for b in monB.index:
            xb = monB['x_coord'][b]
            yb = monB['y_coord'][b]
            zb = monB['z_coord'][b]
            num_b = monB['atom_num'][b]
            b_coord = np.array((xb, yb, zb), dtype = float)
            distance = round(np.linalg.norm(b_coord - a_coord), 3)
            sdistance = str(distance)
            if distance <= cuttoff:
                dist_dict = {
                    'atom_Ai':      num_a, 
                    'atom_Bj':      num_b,
                    'distance':     distance
                }
                distances_db.append(dist_dict)
            else:
                continue
    return distances_db
dist_df = pd.DataFrame(calc_dists())
np.savetxt(r'/home/joaov/python-mmpbsa/mmpbsa/testing/test_dists.log', dist_df.values, fmt='%s', delimiter='\t')
        #print(distance)

#with open(log_dist, "a") as log:
#    log.write(('{0:5} {1:5} {2:5}'.format("num_a", "num_b", "distance"+"\n")))
#    log.write('{0:5} {1:5} {2:5}'.format(num_a, num_b, sdistance+"\n"))
# %%
def dist(atom1, atom2):
    a = atom1
    b = atom2
    num_a = df.index[df['atom_num'] == str(a)].tolist()
    num_b = df.index[df['atom_num'] == str(b)].tolist()
    xa = df['x_coord'][num_a]
    ya = df['y_coord'][num_a]
    za = df['z_coord'][num_a]
    xb = df['x_coord'][num_b]
    yb = df['y_coord'][num_b]
    zb = df['z_coord'][num_b]

    a_coord = np.array((xa, ya, za), dtype = float)
    b_coord = np.array((xb, yb, zb), dtype = float)
    distance = round(np.linalg.norm(b_coord - a_coord), 3)
    
    return distance
# dist(1,1150)

