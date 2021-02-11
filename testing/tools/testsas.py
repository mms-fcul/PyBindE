
import subprocess
import sys
from pathlib import Path
import sys
sys.path.insert(0, '/home/joaov/github/mmpbsa/testing')
import pyjoaov as pj
import pickle

file_path ="/home/joaov/BACKUPS/mmpbsa/testing/saved_files/frame00100.pdb"
g = 0.00542  # gamma -> kcal/(mol‚A2)
b = 0.92    # beta  -> kcal/mol

#gmx20 = Path("/programs/gromacs-OpenMPI2/gromacs-AMD/gromacs-2020.2/bin/gmx")
gmx20_path = Path("/usr/bin/gmx")
g4_sas = Path("/programs/gromacs-OpenMPI2/gromacs-AMD/gromacs-4.0.7_pH_I/bin/g_sas")
g4_editconf = Path("/programs/gromacs-OpenMPI2/gromacs-AMD/gromacs-4.0.7_pH_I/bin/editconf")
# gro="frame00100.gro" ----> gro right now is given by sys.argv[1]

filepath = "/home/joaov/github/mmpbsa/gromacs-dependent/files/"
tpr_path = filepath + "sas.tpr"
ndx_path = filepath + "pb.ndx"

results=pj.run_g_sasa(gmx20_path,file_path,tpr_path,ndx_path)

P   =results[0]
MA  =results[1]
MB  =results[2]


print("mon1",float(MA))
print("mon2",float(MB))
print("complex",float(P))
