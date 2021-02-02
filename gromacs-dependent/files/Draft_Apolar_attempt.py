#!/usr/bin/python3
# -*- coding: utf-8 -*-

import argparse
import os
import subprocess
import sys
from pathlib import Path

#--------------------------------------------------------
#bash commands ran in OG version:
subprocess.Popen(["/bin/rm", '-f',"aux*.xvg","conf.gro"])

#-----------------------variables ------------------------
g=0.00542 # gamma -> kcal/(mol‚A2)
b=0.92    # beta  -> kcal/mol

gmx20 = Path("/home/joaov/downloads/gromacs-2020.3/build/bin/")#tested in home pc
g4_sas = Path("/home/joaov/downloads/gromacs-2020.3/build/bin/g_sas")#to do: edit path and test on lab pc
g4_editconf = Path("/home/joaov/downloads/gromacs-2020.3/build/bin/editconf")
#gro="frame00100.gro" ----> gro right now is given by sys.argv[1]
tpr="sas.tpr"
ndx="pb.ndx"

#------------------------gromacs---------------------------
#for gromacs 2020.3, this works:

def run_sasa(output,option):
    a = subprocess.Popen(['gmx', 'editconf','-f',sys.argv[1],'-o',"conf.gro",'-n',ndx],stdin=subprocess.PIPE)
    a.communicate(b'Dimer\n')
    a.wait()
    p = subprocess.Popen(['gmx', 'sasa','-f',"conf.gro",'-s',tpr,'-n',ndx,'-o',output],stdin=subprocess.PIPE)
    if option == "Dimer":
        p.communicate(b'Dimer\n')
        p.wait()
    elif option == "MonA":
        p.communicate(b'MonomerA\n')
        p.wait()
    elif option == "MonB":
        p.communicate(b'MonomerB\n')
        p.wait()
        
run_sasa("aux_Prot.xvg","Dimer")
run_sasa("aux_MonA.xvg","MonA")
run_sasa("aux_MonB.xvg","MonB")

#for gromacs4, it should be something like this:
'''
def run_sasa(output,option):
    a = subprocess.Popen([g4_editconf,'-f',sys.argv[1],'-o',"conf.gro",'-n',ndx],stdin=subprocess.PIPE)
    a.communicate(b'Dimer\n')
    a.wait()
    p = subprocess.Popen([g4_sas,'-f',gro,'-s',tpr,'-n',ndx,'-ndots',50,'-o',output],stdin=subprocess.PIPE)
    if option == "Dimer":
        p.communicate(b'Dimer\nDimer\n') #how to change these in a function? they must be a string..use if?
        p.wait()
    elif option == "MonA":
        p.communicate(b'MonomerA\nMonomerA\n')
        p.wait()
    elif option == "MonB":
        p.communicate(b'MonomerB\nMonomerB\n')
        p.wait()
        
run_sasa(aux_Prot.xvg,Dimer)
run_sasa(aux_MonA.xvg,MonA)
run_sasa(aux_MonB.xvg,MonB)
'''
#------------------------------------------------------------------------
#commands to run next, in final version 

MA = "`tail -n 1 aux_MonA.xvg | awk '{print $4*100}'`"
MB = "`tail -n 1 aux_MonB.xvg | awk '{print $4*100}'`"
P = "` tail -n 1 aux_Prot.xvg | awk '{print $4*100}'`"
#convert kcal in kJ (factor) = 4.184
subprocess.Popen("echo {} {} {} | awk -v g={} -v b={} '{{print (($3*g+b) - ( ($1*g+b) + ($2*g+b))) * 4.184,\"kJ/mol\"}}'".format(MA,MB,P,g,b),shell=True)

c=subprocess.Popen(["/bin/rm", '-f',"aux*.xvg","conf.gro"])
c.communicate(b'\n')