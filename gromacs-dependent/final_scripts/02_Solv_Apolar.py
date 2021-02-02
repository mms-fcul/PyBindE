#!/usr/bin/python3
# -*- coding: utf-8 -*-

import subprocess
import sys
from pathlib import Path

# --------------------------------------------------------
# bash command
subprocess.Popen(["/bin/rm", '-f', "aux*.xvg", "conf.gro"])

# -----------------------variables ------------------------
g = 0.00542  # gamma -> kcal/(mol‚A2)
b = 0.92    # beta  -> kcal/mol

#gmx20 = Path("/programs/gromacs-OpenMPI2/gromacs-AMD/gromacs-2020.2/bin/gmx")
g4_sas = Path("/programs/gromacs-OpenMPI2/gromacs-AMD/gromacs-4.0.7_pH_I/bin/g_sas")
g4_editconf = Path("/programs/gromacs-OpenMPI2/gromacs-AMD/gromacs-4.0.7_pH_I/bin/editconf")
# gro="frame00100.gro" ----> gro right now is given by sys.argv[1]
tpr = "sas.tpr"
ndx = "pb.ndx"
# ------------------------gromacs---------------------------
# for gromacs 2020.x, this works:
'''
def run_sasa(output, option):
    a = subprocess.Popen([gmx20, 'editconf', '-f', sys.argv[1],'-o', "conf.gro", '-n', ndx], stdin=subprocess.PIPE)
    a.communicate(b'Dimer\n')
    a.wait()
    p = subprocess.Popen([gmx20, 'sasa', '-f', "conf.gro", '-s',tpr, '-n', ndx, '-o', output], stdin=subprocess.PIPE)
    if option == "Dimer":
        p.communicate(b'Dimer\n')
        p.wait()
    elif option == "MonA":
        p.communicate(b'MonomerA\n')
        p.wait()
    elif option == "MonB":
        p.communicate(b'MonomerB\n')
        p.wait()

run_sasa("aux_Prot.xvg", "Dimer")
run_sasa("aux_MonA.xvg", "MonA")
run_sasa("aux_MonB.xvg", "MonB")
'''
# for gromacs-4.0.7, this works:
def run_sasa(output, option):
    a = subprocess.Popen([g4_editconf, '-f', sys.argv[1],'-o', "conf.gro", '-n', ndx], stdin=subprocess.PIPE)
    a.communicate(b'Dimer\n')
    a.wait()
    p = subprocess.Popen([g4_sas, '-f', "conf.gro", '-s', tpr, '-n',ndx, '-ndots', '50', '-o', output], stdin=subprocess.PIPE)
    if option == "Dimer":
        # how to change these in a function? they must be a string..use if?
        p.communicate(b'Dimer\nDimer\n')
        p.wait()
    elif option == "MonA":
        p.communicate(b'MonomerA\nMonomerA\n')
        p.wait()
    elif option == "MonB":
        p.communicate(b'MonomerB\nMonomerB\n')
        p.wait()

run_sasa("aux_Prot.xvg", "Dimer")
run_sasa("aux_MonA.xvg", "MonA")
run_sasa("aux_MonB.xvg", "MonB")

# ------------------------------------------------------------------------
# commands to run next

MA = "`tail -n 1 aux_MonA.xvg | awk '{print $4*100}'`"
MB = "`tail -n 1 aux_MonB.xvg | awk '{print $4*100}'`"
P = "` tail -n 1 aux_Prot.xvg | awk '{print $4*100}'`"
# convert kcal in kJ (factor) = 4.184
subprocess.Popen("echo {} {} {} | awk -v g={} -v b={} '{{print (($3*g+b) - ( ($1*g+b) + ($2*g+b))) * 4.184,\"kJ/mol\"}}'".format(MA, MB, P, g, b), shell=True)

c=subprocess.Popen(["/bin/rm", '-f', "aux*.xvg", "conf.gro"])
c.communicate(b'\n')
