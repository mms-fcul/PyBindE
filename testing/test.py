#!/usr/bin/python3
# -*- coding: utf-8 -*-
import argparse
import os
import subprocess
from subprocess import run, PIPE
import sys
from pathlib import Path

'''

with open("out_filename", 'w') as out_file:
     parsed_line = out_file.readlines()
     out_file.write(parsed_line)
     
curr="bbbbbb"

text_mdp = ["integrator          =  md", "\n",
            "nsteps              =  0", "\n",
            "coulombtype         = Reaction-Field", "\n",
            ";coulombtype         = PME", "\n",
            "epsilon_rf          = 54.0", "\n",
            "rlist               = $cutoff", "\n",
            "rcoulomb            = $cutoff", "\n",
            "rvdw                = $cutoff", "\n",
            "energygrps          = MonomerA MonomerB"
            ]

mdp_name=curr+".mdp"
with open(mdp_name, "w") as out_file:
  out_file.writelines(text_mdp)
'''
'''
for i in Coul-SR:MonomerA-MonomerB LJ-SR:MonomerA-MonomerB
do
    ${grom}/g_energy -f ${curr}.edr -o ${curr}_${i}.xvg >> ${curr}.out 2>&1 <<EOF
${i}
EOF
'''
#g4 = Path("/programs/gromacs-OpenMPI2/gromacs-AMD/gromacs-4.0.7_pH_I/bin/")
'''
choices = ["Coul-SR:MonomerA-MonomerB", "LJ-SR:MonomerA-MonomerB"]
log = open('log.out','a')
for i in choices:
  print(i)
 # a = subprocess.Popen([g4+"g_energy", '-f', curr+".edr",'-o', curr+"_"+i+".xvg"],  stdin=subprocess.PIPE, stdout=log, stderr=log)
  a = subprocess.Popen(['read',"var"],stdin=subprocess.PIPE,shell=True)
  a.communicate(b'banana\n')
  #a.wait()
  '''
  
tpr = "sas.tpr"
ndx = "pb.ndx"

gmx20 = Path("/programs/gromacs-OpenMPI2/gromacs-AMD/gromacs-2020.2/bin/gmx")
choices = ["Dimer\n", "MonomerA\n","MonomerB\n"]
for i in choices:
  p = subprocess.run(["gmx", 'sasa', '-f', "frame00000.gro", '-s',tpr, '-n', ndx, '-o', i], input=i.encode("utf8"))
  
  #p.communicate(i)
  #p.stdin.write(str(i)+"\n")
  #p.stdin.flush()
  #p.wait()
  '''
subprocess.Popen("echo {} {} {} | awk -v g={} -v b={} '{{print (($3*g+b) - ( ($1*g+b) + ($2*g+b))) * 4.184,\"kJ/mol\"}}'".format(MA, MB, P, g, b), shell=True)
'''

for i in choices:
 # p = subprocess.run([gmx20, "energy",'-f', curr+'.edr', '-o', curr+"_"+i+".xvg"], input=i.encode("utf8"), stdout=log_g_energy, stderr=log_g_energy)
  #tail=["tail -n 1", curr+"_"+i+".xvg"," | awk","'{{printf \"%f \", $2}}'"]
  file=curr+"_"+i+".xvg"
  #tail=["tail", "-n", "1",file]
  #print(file)
  #p=subprocess.Popen(tail, stdout=tail_test, shell=True)
  #line = subprocess.check_output(['tail', '-1', "{}".format(file)])
  #out = check_output(tail)
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
#done | awk -v c=${conv} '{printf "E^MM(Coul)= %6.2f kJ/mol\nE^MM(VdW) = %6.2f kJ/mol\n", ($1)/c,($2)/c}'
#############################
#c=subprocess.Popen(["/bin/rm", '-f', "conf.gro","traj.trr","TMP* *~ *# .*~ .*#"])
#c.communicate(b'\n')
#############################