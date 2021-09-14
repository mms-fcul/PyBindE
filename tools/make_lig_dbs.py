#!/usr/bin/python3

import pandas as pd
import sys
import os
import numpy as np
import copy

lib_dir=os.path.dirname(os.path.realpath(__file__))
databases_path = "/home/joaov/Projetos/MMPBSA/validation/PyEBind_validation/databases/"

txt_df_res      = databases_path + "residues_rtp.txt"
txt_df_nb       = "/home/joaov/Projetos/MMPBSA/validation/FF/gromos54a7_atb.ff/nonbonded_itp.txt"
txt_df_charges  = databases_path + "charges.txt"


def get_itp(system_name):
    #system_name=gro_file.split("/")[-3]
    itp_file="/home/joaov/Projetos/MMPBSA/validation/ligands/02_ATB_files/final_itps/"+system_name+"_final.itp"
    
    lines = [] 
    found_section=False
    itp_dict=[]
    with open(itp_file) as itp:
        for line in itp:
            if line.startswith('[ atoms ]'):
                found_section = True
                continue
            if found_section:
                if line.startswith('[ bonds ]'):
                    found_section = False
                elif '[ atoms ]' in line:
                    continue
                elif line.startswith(';'):
                    continue
                else:
                    #s_line = str(line).rstrip('\n')
                    splitted_line = [line[0:5], line[5:11], line[11:16], line[16:23], line[23:31], line[31:36], line[36:45], line[45:54]]
                    lines.append(splitted_line)
            
        for line in lines:
            line_dict = {
                'atom_num':     line[0].strip(), 
                'atom_type':    line[1].strip(), 
                'res_num':      line[2].strip(), 
                'res_name':     line[3].strip(), 
                'atom_name':    line[4].strip(), 
                'cgnr':         line[5].strip(), 
                'charge':       line[6].strip(),
                'mass':         line[7].strip()
            }
            itp_dict.append(line_dict)


    itp_df = pd.DataFrame(itp_dict)
    return itp_df


system_name = sys.argv[1] 
itp_df=get_itp(system_name)
lig_charges = itp_df[['res_name','atom_name','charge']]
lig_res = itp_df[['res_name','atom_name','atom_type']]

#print("THIS IS LIG_RES",lig_res)

LIG_RES = "/home/joaov/Projetos/MMPBSA/validation/systems/"+system_name+"/mmpbsa/lig_res_types.txt"
np.savetxt(LIG_RES,lig_res, fmt=('%5.5s', '%5.5s', '%5.5s'))


f_crg = lib_dir + "/delphi4github/delphi4py/example/DataBaseT.crg"
f_siz = lib_dir + "/delphi4github/delphi4py/example/DataBaseT.siz"

df_res      = pd.read_table(txt_df_res,     delim_whitespace=True,header=None,names=['res_name','atom_name','atom_type']) 
df_charges  = pd.read_table(txt_df_charges, delim_whitespace=True,header=None,names=['atom_name','res_name','charge'])
df_nb       = pd.read_table(txt_df_nb,      delim_whitespace=True,header=None,names=['i','j','c6','c12']) 

#print(df_res)
#print(df_charges)
#print(df_nb)

# =======================MAKE RADII==============================
lig_res=lig_res.rename(columns = {"atom_type": "i"})
lig_res['j']='OW'
#lig_res['i'] = lig_res['i'].str.replace('SDmso', 'SDms')
lig_res['i'] = lig_res['i'].str.replace('CH2R', 'CH2r')
lig_res_nb=pd.merge(lig_res, df_nb,  how = 'left', on = ['i', 'j'])
#print(lig_res_nb)


nans = lig_res_nb.query("c6 != c6")
if not nans.empty: print(nans); print("There are atoms with no attribution in atom_type!")


rt=2

kT = 2.49432

WATER_RADIUS= 1.4241233005024403

def radius_from_lennard(c6, c12, rt=2):
    rmin = (2 * c12 / c6) ** (1.0 / 6)
    
    kT = 2.49432
    radius = (rmin ** (-6) + (rt * kT / c12) ** 0.5) ** (-1 / 6)

    return 10 * radius - WATER_RADIUS

#radius = radius_from_lennard(c6, c12, rt=2)
#import numpy as np
#formula at /programs/meadTools2.0.1/makepqr

lig_res_nb['radius'] = radius_from_lennard(lig_res_nb['c6'], lig_res_nb['c12'], rt=2)
# 10 * ( (((2 * lig_res_nb['c12'] / lig_res_nb['c6']) ** (1/6)) ** (-6)+ (rt * kT / lig_res_nb['c12']) ** 0.5) ** (-1/6) ) - WATER_RADIUS

#print(lig_res_nb)
nans = lig_res_nb.query("radius != radius")
#if not nans.empty: print(nans); print("There are atoms with no attribution in atom_type!")

print(lig_res_nb)
lig_res_nb.fillna(0.001,inplace=True)
nans = lig_res_nb.query("radius != radius")
print(lig_res_nb)
lig_res_nb['radius'] = lig_res_nb['radius'].astype(float)
classifier_df = lig_res_nb[['res_name','atom_name','i']]
lig_res_nb = lig_res_nb[['atom_name','res_name','radius']]
lig_res_nb = lig_res_nb.drop_duplicates()
LIG_RADII = "/home/joaov/Projetos/MMPBSA/validation/systems/"+system_name+"/mmpbsa/lig_radii.txt"
np.savetxt(LIG_RADII,lig_res_nb, fmt=('%-4.4s', '%4.4s', '%7.3f'))
print("LIG_RADII saved")
CLASSIFIER_PART = "/home/joaov/Projetos/MMPBSA/validation/systems/"+system_name+"/mmpbsa/add_to_classifier.txt"
np.savetxt(CLASSIFIER_PART,classifier_df, delimiter=' ',fmt='%s')
print("LIG_classifier saved")
# ========================MAKE CHARGES=============================
lig_charges = itp_df[['atom_name','res_name','charge']]
lig_charges['charge'] = lig_charges['charge'].astype(float)
LIG_CHARGES = "/home/joaov/Projetos/MMPBSA/validation/systems/"+system_name+"/mmpbsa/lig_charges.txt"
np.savetxt(LIG_CHARGES,lig_charges, fmt=('%-4.4s', '%4.4s', '%11.3f'))
print("LIG_CHARGES saved")












# ================================================================
#WATER_RADIUS= 1.4241233005024403
#def radius_from_lennard(c6, c12, rt=2):
#    rmin = (2 * c12 / c6) ** (1.0 / 6)
#    kT = 2.49432
#    radius = (rmin ** (-6) + (rt * kT / c12) ** 0.5) ** (-1 / 6)
#    return 10 * radius - WATER_RADIUS
#radius = radius_from_lennard(c6, c12, rt=2)
# ================================================================


#inverted_df_nb = copy.deepcopy(df_nb)
#inverted_df_nb = inverted_df_nb.rename(columns = {"i": "j", "j": "i"})
#total_nb = pd.concat([df_nb, inverted_df_nb], ignore_index=True)
#total_nb_drop=total_nb.drop_duplicates()
#print(total_nb_drop)

#total_res = pd.concat([df_res, lig_res], ignore_index=True)
#df_res_save = "/home/joaov/Projetos/MMPBSA/validation/systems/"+system_name+"/mmpbsa/residues_rtp.#txt"
#np.savetxt(df_res_save,total_res, fmt=('%5.4s', '%5.4s', '%5.4s'))
#print(total_res)

#total_charges = pd.concat([df_charges, lig_charges], ignore_index=True)
#df_charges_save = "/home/joaov/Projetos/MMPBSA/validation/systems/"+system_name+"/mmpbsa/charges.#txt"
#np.savetxt(df_charges_save,total_charges, fmt=('%5.4s', '%5.4s', '%5.4s'))
#print(total_charges)
#
#


#total_res=total_res.rename(columns = {"atom_type": "i"})
#total_res['j']='OW'
#print(total_res)
#total_res['i'] = total_res['i'].str.replace('SDmso', 'SDms')
#total_res['i'] = total_res['i'].str.replace('CH2R', 'CH2r')
#total_nb = pd.merge(total_res, df_nb,  how = 'left', on = ['i', 'j'])
#print(total_nb)
