#!/usr/bin/python3

import pandas as pd
import numpy as np

working_path ="/home/joaov/python-mmpbsa/mmpbsa/testing/"
#pickle_charges = "/home/joaov/python-mmpbsa/mmpbsa/testing/saved_dataframes/charges_crg.pkl"
charges_db      = working_path + "DataBaseT.crg"
charges_txt     = working_path + "saved_dataframes/charges_crg.txt"

charges=pd.read_table(charges_db,delim_whitespace=True,skip_blank_lines=True,header=1,names=['atom_name','res_name','charge'],skiprows=1) 

#charges.to_pickle(pickle_charges)
np.savetxt(charges_txt,charges.values, fmt=('%5.4s', '%5.4s', '%7.4f'))