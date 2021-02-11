#idea: generate new databases, modified to have complete NTR and CTR entries
#IDEA: find correct protonation

import pandas as pd
import numpy as np

working_path ="/home/joaov/python-mmpbsa/mmpbsa/testing/"
#pickle_charges = "/home/joaov/python-mmpbsa/mmpbsa/testing/saved_dataframes/charges_crg.pkl"
charges_db      = working_path + "DataBaseT.crg"
charges_txt     = working_path + "saved_dataframes/charges_crg.txt" # saving file

charges=pd.read_table(charges_db,delim_whitespace=True,skip_blank_lines=True,header=1,names=['atom_name','res_name','charge'],skiprows=1)

#mask=charges.query('res_name starts with CT')
def slices(df,string):
  dfa=df['res_name'].str.startswith(string)
  DF=df.loc[dfa,'res_name']
  return DF

ct_terminus=slices(charges,'CT')
nt_terminus=slices(charges,'NT')
#print(charges.loc[charges['res_name'].str.contains('CT')],charges.loc[charges['res_name'].str.contains('NT')])
#print(ct_terminus.unique())
#print(nt_terminus.unique())

#new_df = new_df.append(new_row, ignore_index=True)
ctr_list= \
{'atom_name':'CT',   'res_name':'CTR', 'charge':0.270  }, \
{'atom_name':'OT1',  'res_name':'CTR', 'charge':-0.635 }, \
{'atom_name':'OT2',  'res_name':'CTR', 'charge':-0.635 }, \
{'atom_name':'HO11', 'res_name':'CTR', 'charge':0.000  }, \
{'atom_name':'HO21', 'res_name':'CTR', 'charge':0.000  }, \
{'atom_name':'HO12', 'res_name':'CTR', 'charge':0.000  }, \
{'atom_name':'HO22', 'res_name':'CTR', 'charge':0.000  }

ntr_list= \
{'atom_name':'N',    'res_name':'NTR', 'charge':0.1290 }, \
{'atom_name':'H1',   'res_name':'NTR', 'charge':0.2480 }, \
{'atom_name':'H2',   'res_name':'NTR', 'charge':0.2480 }, \
{'atom_name':'H3',   'res_name':'NTR', 'charge':0.2480 }, \
{'atom_name':'CA',   'res_name':'NTR', 'charge':0.1270 }, \
{'atom_name':'C',    'res_name':'NTR', 'charge':0.4500 }, \
{'atom_name':'O',    'res_name':'NTR', 'charge':-0.4500} \

ctr_df = pd.DataFrame(ctr_list)
ntr_df = pd.DataFrame(ntr_list)

if not 'CTR' in ct_terminus.unique():
  charges=charges.append(ctr_df, ignore_index=True)

if not 'NTR' in nt_terminus.unique():
  charges=charges.append(ntr_df, ignore_index=True)

print(charges)

#print(charges[charges['res_name'].str.startswith('CT')])
#print(charges[charges['res_name'].str.startswith('NT')])
#charges.loc[mask,'res_name']
#print(charges.head(20))
#charges.to_pickle(pickle_charges)

#np.savetxt(charges_txt,charges.values, fmt=('%5.4s', '%5.4s', '%7.4f'))

