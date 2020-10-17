#!/usr/bin/python3

from os import read
import re
import pandas as pd
from biopandas.pdb import PandasPdb
import numpy as np
import subprocess

pdb_file = "/home/joaov/python-mmpbsa/mmpbsa/testing/dimer.pdb"
# %%
import matplotlib.pyplot as plt
#ppdb=PandasPdb().read_pdb('/home/joaov/python-mmpbsa/mmpbsa/testing/dimer.pdb')
pdb_file = "/home/joaov/python-mmpbsa/mmpbsa/testing/dimer.pdb"
#atom_coord=ppdb.df["atom_number","x_coord","y_coord","z_coord"]

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
                'chain_id':    line[4].strip(),
                'x_coord':    float(line[6].strip()),
                'y_coord':    float(line[7].strip()),
                'z_coord':    float(line[8].strip())
            }
            pdb_dict.append(line_dict)
    return pdb_dict


pdb_coord = read_pdb(pdb_file)
# print(pdb_coord)
df = pd.DataFrame(pdb_coord)
# print(df)

# %%
#def make_index(PDB_dictionary):
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

#monomer_index= make_index(pdb_coord)
#print(monomer_index)
  
#ppdb=PandasPdb().read_pdb(pdb_file)  

df=pd.DataFrame(pdb_coord)
#df=df.set_index('atom_num')
  
monA=df.query('chain_id =="A"')
monB=df.query('chain_id =="B"')
#print(monA)
#print(monB)

log_dist="/home/joaov/python-mmpbsa/mmpbsa/testing/pdb_dist.log"
c=subprocess.Popen(["/bin/rm", '-f', log_dist])
c.communicate(b'\n')

with open(log_dist,"a") as log:
    log.write(('{0:5} {1:5} {2:5}'.format("num_a", "num_b", "distance"+"\n")))
    for a in monA.index:
        xa=monA['x_coord'][a]
        ya=monA['y_coord'][a]
        za=monA['z_coord'][a]
        num_a=monA['atom_num'][a]
        a_coord=np.array((xa, ya, za), dtype=float)
        for b in monB.index:
            xb=monB['x_coord'][b]
            yb=monB['y_coord'][b]
            zb=monB['z_coord'][b]
            num_b=monB['atom_num'][b]
            b_coord=np.array((xb, yb, zb), dtype=float)
            distance = round(np.linalg.norm(b_coord - a_coord),3)
            sdistance=str(distance)
            log.write('{0:5} {1:5} {2:5}'.format(num_a, num_b, sdistance+"\n"))
            #print(distance)
# %%
