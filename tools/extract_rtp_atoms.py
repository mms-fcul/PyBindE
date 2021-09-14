#!/usr/bin/python3

import re
import pandas as pd
import numpy as np

#extracting atom info for each residue inside an rtp
rtp_file = "/home/machuque/gromacs/Constant_pH_MD/FFs_5.0/G54a7pH.ff/aminoacids.rtp" #"/home/joaov/python-mmpbsa/mmpbsa/testing/forcefield_files/ffG54a7pHt.rtp"
lines = [] 
found_section=False
exclusion_list = ('[ bonds ]','[ angles ]','[ impropers ]','[ dihedrals ]','[ exclusions ]',
                ' [ bonds ]',' [ angles ]',' [ impropers ]',' [ dihedrals ]',' [ exclusions ]',
                '[ bondedtypes ]')
aa_list = []
with open(rtp_file) as rtp:
    #next(rtp)
    #next(rtp)
    #next(rtp)
    for line in rtp:
        if line.startswith('[ ') and not line.startswith(exclusion_list):
            aa_list.append(line)
            found_section = True
            cline = str(line).rstrip('\n')
            lines.append(cline)
            continue
        if found_section:
            if line.startswith(exclusion_list):
                found_section = False
            elif '[ atoms ]' in line:
                continue
            else:
                s_line = str(line).rstrip('\n')
                lines.append(s_line)

log_rtp = "/home/joaov/python-mmpbsa/mmpbsa/testing/log_rtp.rtp"
with open(log_rtp, "w") as log:
    log.writelines("%s\n" % place for place in lines)
    
#adding CTR and NTR terminus
terminus_entries ="""[ CTR ]
   CT     C       0.270
  OT1    OM      -0.635
  OT2    OM      -0.635
 HO11     H       0.000
 HO21     H       0.000
 HO12     H       0.000
 HO22     H       0.000
[ NTR ]
    N    NT       0.129
   H1     H       0.248
   H2     H       0.248
   H3     H       0.248
   CA   CH1       0.127"""

with open(log_rtp,"a") as logrtp:
    logrtp.write(terminus_entries)
 
#make dataframe   
def read_resdb(file):
    residues = []
    with open(file) as res_db:
        for line in res_db:
            if line.startswith('['):
                res = str(line).rstrip('\n')
                res1 = re.sub(r'\[', '', res).strip()
                res2 = re.sub(r'\]', '', res1).strip()

            elif not line.startswith(';'):
                cleanline = re.sub('\s+', ' ', line).strip()
                fields = cleanline.split(' ')
                res_dict = {
                    'res_name':              res2,
                    'atom_name':        fields[0].rstrip(),
                    'atom_type':        fields[1].rstrip()
                    #'charge':     float(fields[2])
                }
                residues.append(res_dict)
    return residues

df_res = pd.DataFrame(read_resdb(log_rtp)) # df w/ atom names/types from rtp

#save to pickle (can be called and used with read_pickle)
df_res_txt = "/home/joaov/python-mmpbsa/mmpbsa/testing/saved_dataframes/residues_rtp.txt"
#df_res.to_pickle(pickle_df_res)
np.savetxt(df_res_txt,df_res, fmt=('%5.4s', '%5.4s', '%5.4s'))