#!/usr/bin/python3

import re
import pandas as pd
import numpy as np

#extracting atom info for each residue inside an itp
nonbonded_itp = "/home/joaov/Projetos/MMPBSA/validation/FF/gromos54a7_atb.ff/ffnonbonded.itp"

#pickled_df = "/home/joaov/python-mmpbsa/mmpbsa/testing/saved_dataframes/nonbonded_itp.pkl" #file used to save df
nb_df_txt = "/home/joaov/Projetos/MMPBSA/validation/FF/gromos54a7_atb.ff/nonbonded_itp.txt" #file used to save df

#equal_atom_lines = []
#diff_atom_lines = []
equal_atom_dicts = []
diff_atom_dicts = []
exclusion_list = ('[ pairtypes ]','[ atomtypes ]','[ nonbond_params ]')
with open(nonbonded_itp) as itp:
    for line in itp:
        line = line.split(";")[0] #remove comments
        line = re.sub('\s+', ' ', line) #simplify column separator
        line = str(line).strip('\n').strip() #clean line string
        if line.startswith('[ atomtypes ]'):
            curr_section='atomtypes'
            #equal_atom_lines.append(line)
            continue
        elif line.startswith('[ nonbond_params ]'):
            curr_section='nonbond'
            #diff_atom_lines.append(line)
            continue
        if line != '' and not line.startswith('[ ') and curr_section == 'atomtypes':
            #equal_atom_lines.append(line)
            fields = line.split(' ')
            fields = list(filter(None, fields)) 
            equal_dict = {
                'i':          fields[0].rstrip(), 
                'j':          fields[0].rstrip(), 
                'c6':   float(fields[5]), 
                'c12':  float(fields[6])
            }
            equal_atom_dicts.append(equal_dict)
        elif line != '' and not line.startswith('[ ') and curr_section == 'nonbond':
            #diff_atom_lines.append(line)
            fields = line.split(' ')
            fields = list(filter(None, fields))
            diff_dict = {
                'i':          fields[0].rstrip(), 
                'j':          fields[1].rstrip(), 
                'c6':   float(fields[3]), 
                'c12':  float(fields[4])
            }
            diff_atom_dicts.append(diff_dict)
        elif line.startswith('[ pairtypes ]'):
            break
        
equal_df = pd.DataFrame(equal_atom_dicts)
diff_df = pd.DataFrame(diff_atom_dicts)
alt_diff = diff_df.rename(columns = {"i": "j", "j": "i"}) # create switched refs
complete_diff_df = diff_df.append(alt_diff,ignore_index=True) # now we have the reverse combinations too
complete_df= equal_df.append(complete_diff_df)

#complete_df.to_pickle(pickled_df)

np.savetxt(nb_df_txt,complete_df, fmt=('%5.4s', '%5.4s', '%14.7e', '%14.7e'))