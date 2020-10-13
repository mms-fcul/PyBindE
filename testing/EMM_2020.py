#!/usr/bin/python3

import argparse
import os
import subprocess
from subprocess import PIPE
import sys
from pathlib import Path

# -----------------------variables ------------------------
syst="b2M"
curr="TMP"
cutoff="3.0"
#conv=4.184 # 1 kcal = 4.184 kJ
conv=1 # 1 kcal = 4.184 kJ

#gmx20 = Path("/programs/gromacs-OpenMPI2/gromacs-AMD/gromacs-2020.2/bin/gmx")
gmx20= Path("/usr/bin/gmx")
'''
# ------------------------gromacs---------------------------
log_editconf = open('log_editconf.out','a')
a = subprocess.Popen([gmx20, "editconf", '-f', sys.argv[1],'-o', "conf.gro", '-n', syst+".ndx"], stdin=subprocess.PIPE, stdout=log_editconf, stderr=log_editconf)
a.communicate(b'Dimer\n')
a.wait()

text_mdp = ["integrator          =  md", "\n",
            "nsteps              =  0", "\n",
            "coulombtype         = Reaction-Field", "\n",
            ";coulombtype         = PME", "\n",
            "epsilon_rf          = 54.0", "\n",
            "rlist               = "+cutoff, "\n",
            "rcoulomb            = "+cutoff, "\n",
            "rvdw                = "+cutoff, "\n",
            "energygrps          = MonomerA MonomerB"
            ]

mdp_name=curr+".mdp"
with open(mdp_name, "w") as out_file:
    out_file.writelines(text_mdp)

# for gromacs-4.0.7
#grommp
log_grompp = open('log_grompp.out','a')
subprocess.run([gmx20, "grompp", '-f', curr+".mdp",
                      "-po",curr+"_out.mdp",
                      "-pp",curr+"_processed.top",
                      '-o', syst+".tpr",
                      "-c","conf.gro",
                      '-n', syst+".ndx",
                      "-p",syst+".top",
                      "-maxwarn","1000"], stdin=subprocess.PIPE, stdout=log_grompp, stderr=log_grompp)
#mdrun

log_mdrun = open('log_mdrun.out','a')
subprocess.run([gmx20, "mdrun", '-s', syst+".tpr",
                      "-x",curr+".xtc",
                      "-c",curr+".gro",
                      '-e',curr+".edr",
                      "-g",curr+".log"], stdin=subprocess.PIPE, stdout=log_mdrun, stderr=log_mdrun)
'''
# ------------------------------------------------------------------------
# commands to run next
'''
subprocess.Popen("echo {} {} {} | awk -v g={} -v b={} '{{print (($3*g+b) - ( ($1*g+b) + ($2*g+b))) * 4.184,\"kJ/mol\"}}'".format(MA, MB, P, g, b), shell=True)
'''
#############################
#falta fazer esta parte:
log_g_energy= open('log_g_energy.out','a')
choices=["Coul-SR:MonomerA-MonomerB", "LJ-SR:MonomerA-MonomerB"]

for i in choices:
 # p = subprocess.run([gmx20, "energy",'-f', curr+'.edr', '-o', curr+"_"+i+".xvg"], input=i.encode("utf8"), stdout=log_g_energy, stderr=log_g_energy)

  file=curr+"_"+i+".xvg"
  
  with open(file, 'r') as f:
    lines = f.read().splitlines()
    energy = (lines[-1]).split()[1]
  
  original_stdout = sys.stdout
  with open('energy_vals','a') as en:
    sys.stdout = en
    print(energy)
    sys.stdout = original_stdout

with open('energy_vals', 'r') as f:
    lines = f.read().splitlines()
    energy_Coul = lines[0]
    energy_VdW = lines[1]
print("E^MM(Coul)= {} kJ/mol".format(energy_Coul))
print("E^MM(VdW)= {} kJ/mol".format(energy_VdW))

#if conv != 1:
print("Convertion: E^MM(Coul)= {} kcal/mol".format(float(energy_Coul)/4.184))
print("Convertion: E^MM(VdW)= {} kcal/mol".format(float(energy_VdW)/4.184))
#done | awk -v c=${conv} '{printf "E^MM(Coul)= %6.2f kJ/mol\nE^MM(VdW) = %6.2f kJ/mol\n", ($1)/c,($2)/c}'
#############################
#c=subprocess.Popen(["/bin/rm", '-f', "conf.gro","traj.trr","TMP* *~ *# .*~ .*#"])
#c.communicate(b'\n')
#############################