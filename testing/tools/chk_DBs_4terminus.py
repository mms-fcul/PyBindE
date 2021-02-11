#idea: generate new databases, modified to have complete NTR and CTR entries
#IDEA: find correct protonation

import pandas as pd
import numpy as np
import sys

arg1, arg2 = str(sys.argv[1]), str(sys.argv[2])
print("Arguments:",arg1,arg2)
working_path ="/home/joaov/python-mmpbsa/mmpbsa/testing/"
databases_path = "/home/joaov/github/mmpbsa/testing/saved_dataframes/"
#pickle_charges = "/home/joaov/python-mmpbsa/mmpbsa/testing/saved_dataframes/charges_crg.pkl"
charges_db      = working_path + "DataBaseT.crg"
charges_txt     = working_path + "saved_dataframes/charges_crg.txt" # saving file

txt_df_res      = databases_path + "residues_rtp.txt"
txt_df_nb       = databases_path + "nonbonded_itp.txt"
txt_df_charges  = databases_path + "charges.txt"

#Import databases
df_res      = pd.read_table(txt_df_res,     delim_whitespace=True,header=None,names=['res_name','atom_name','atom_type']) # df w/ atom names/types from rtp
df_charges  = pd.read_table(txt_df_charges, delim_whitespace=True,header=None,names=['atom_name','res_name','charge'])
#df_nb       = pd.read_table(txt_df_nb,      delim_whitespace=True,header=None,names=['i','j','c6','c12']) # df from nb.itp to grab c6 and c12

def slices(df,string):
  dfa=df['res_name'].str.startswith(string)
  DF=df.loc[dfa,'res_name']
  return DF

nt_terminus=slices(df_charges,arg1)
ct_terminus=slices(df_charges,arg2)

#print(df_charges.loc[df_charges['res_name'].str.contains(arg1)],df_charges.loc[df_charges['res_name'].str.contains(arg2)])
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

ctr_crg_df = pd.DataFrame(ctr_list)
ntr_crg_df = pd.DataFrame(ntr_list)

if not arg2 in ct_terminus.unique():
  print("There's no {} in charges database (DatabaseT.crg)".format(arg2))
  #df_charges=df_charges.append(ctr_crg_df, ignore_index=True)
else:
  print("{} is present in charges database (DatabaseT.crg)".format(arg2))

if not arg1 in nt_terminus.unique():
  print("There's no {} in charges database (DatabaseT.crg)".format(arg1))
  #df_charges=df_charges.append(ntr_crg_df, ignore_index=True)
else:
  print("{} is present in charges database (DatabaseT.crg)".format(arg1))

#print(df_charges)

nt_terminus=slices(df_res,arg1)
ct_terminus=slices(df_res,arg2)

if not arg2 in ct_terminus.unique():
  print("There's no {} in residues' database (.rtp)".format(arg2))
  #df_charges=df_charges.append(ctr_crg_df, ignore_index=True)
else:
  print("{} is present in residues' database (.rtp)".format(arg2))

if not arg1 in nt_terminus.unique():
  print("There's no {} in residues' database (.rtp)".format(arg1))
  #df_charges=df_charges.append(ntr_crg_df, ignore_index=True)
else:
  print("{} is present in residues' database (.rtp)".format(arg1))


#print(charges[charges['res_name'].str.startswith('CT')])
#print(charges[charges['res_name'].str.startswith('NT')])
#charges.loc[mask,'res_name']
#print(charges.head(20))
#charges.to_pickle(pickle_charges)

#np.savetxt(charges_txt,charges.values, fmt=('%5.4s', '%5.4s', '%7.4f'))

