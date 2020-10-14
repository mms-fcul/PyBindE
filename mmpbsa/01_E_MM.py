#!/usr/bin/python3

import subprocess
import sys
from pathlib import Path

# -----------------------variables ------------------------
syst="b2M"
curr="TMP"
cutoff="3.0"
#conv=4.184 # 1 kcal = 4.184 kJ
conv=1 # 1 kcal = 4.184 kJ

#gmx20 = Path("/programs/gromacs-OpenMPI2/gromacs-AMD/gromacs-2020.2/bin/gmx")
g4 = Path("/programs/gromacs-OpenMPI2/gromacs-AMD/gromacs-4.0.7_pH_I/bin/")
g4_editconf = Path("/programs/gromacs-OpenMPI2/gromacs-AMD/gromacs-4.0.7_pH_I/bin/editconf")
g4_grompp = Path("/programs/gromacs-OpenMPI2/gromacs-AMD/gromacs-4.0.7_pH_I/bin/grompp")
g4_mdrun = Path("/programs/gromacs-OpenMPI2/gromacs-AMD/gromacs-4.0.7_pH_I/bin/mdrun")
g4_g_energy = Path("/programs/gromacs-OpenMPI2/gromacs-AMD/gromacs-4.0.7_pH_I/bin/g_energy")
# ------------------------gromacs---------------------------
# for gromacs-4.0.7:
# #editconf
log_editconf = open('log_editconf.out','a')
a = subprocess.Popen([g4_editconf, '-f', sys.argv[1],'-o', "conf.gro", '-n', syst+".ndx"], stdin=subprocess.PIPE, stdout=log_editconf, stderr=log_editconf)
a.communicate(b'Dimer\n')
a.wait()
#make mdp
text_mdp = """
integrator      =  md
nsteps          =  0
coulombtype     = Reaction-Field
;coulombtype    = PME
epsilon_rf      = 54.0
rlist           = {0}
rcoulomb        = {0}
rvdw            = {0}
energygrps      = MonomerA MonomerB""".format(cutoff)

mdp_name=curr+".mdp"
with open(mdp_name, "w") as out_file:
    out_file.write(text_mdp)

#grompp
log_grompp = open('log_grompp.out','a')
subprocess.run([g4_grompp, '-f', curr+".mdp",
                      "-po",curr+"_out.mdp",
                      "-pp",curr+"_processed.top",
                      '-o', syst+".tpr",
                      "-c","conf.gro",
                      '-n', syst+".ndx",
                      "-p",syst+".top",
                      "-maxwarn","1000"], stdin=subprocess.PIPE, stdout=log_grompp, stderr=log_grompp)
#mdrun
log_mdrun = open('log_mdrun.out','a')
subprocess.run([g4_mdrun, '-s', syst+".tpr",
                      "-x",curr+".xtc",
                      "-c",curr+".gro",
                      '-e',curr+".edr",
                      "-g",curr+".log",
                      '-nice', "19"], stdin=subprocess.PIPE, stdout=log_mdrun, stderr=log_mdrun)
# ------------------------------------------------------------------------
# getting energies
#############################
log_g_energy= open('log_g_energy.out','a')
choices=["Coul-SR:MonomerA-MonomerB", "LJ-SR:MonomerA-MonomerB"]

for i in choices:
  p = subprocess.run([g4_g_energy,'-f', curr+'.edr', '-o', curr+"_"+i+".xvg"], input=i.encode("utf8"), stdout=log_g_energy, stderr=log_g_energy)

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
#############################
c=subprocess.Popen(["/bin/rm", '-f', "conf.gro","traj.trr","energy_vals","TMP* *~ *# .*~ .*#"])
c.communicate(b'\n')
#############################