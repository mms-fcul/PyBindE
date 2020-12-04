#!/usr/bin/python3

from os import read
import re
import pandas as pd
from biopandas.pdb import PandasPdb

nonbonded_itp = "/home/joaov/python-mmpbsa/mmpbsa/testing/nb.itp"
chargesDB = "/home/joaov/python-mmpbsa/mmpbsa/testing/DataBaseT.crg"
pdb_file = "/home/joaov/python-mmpbsa/mmpbsa/testing/dimer.pdb"
# %%

def read_nb(file):
    NB = []
    lines_no_comments = []
    clean_lines = []
    with open(file) as nb:
        next(nb)
        next(nb)
        lines = nb.read().splitlines()

        for line in lines:
            sep_comment = line.split(";")
            # print(fields)
            lines_no_comments.append(sep_comment[0])
        for line in lines_no_comments:
            cleanline = re.sub('\s+', ' ', line)
            clean_lines.append(cleanline.strip())
        # print(clean_lines)
        clean_lines = list(filter(None, clean_lines))
        for line in clean_lines:
            fields = line.split(' ')
            # print(fields)
            nb_dict = {
                'i':          fields[0], 
                'j':          fields[1], 
                'c6':   float(fields[3]), 
                'c12':  float(fields[4])
            }
            NB.append(nb_dict)
    return NB


NB = read_nb(nonbonded_itp)
# %%
def pair_search(criteria, string):
    search_list = []
    for pair in NB:
        if pair[criteria] == string:
            search_list.append(pair)
    return search_list
# print(pair_search('i','CH3'))
# %%
def read_charges(file):
    clean_lines = []
    charges_dataset = []
    with open(chargesDB) as db:
        next(db)
        next(db)
        next(db)
        lines = db.read().splitlines()
        for line in lines:
            cleanline = re.sub('\s+', ' ', line)
            clean_lines.append(cleanline)

        for line in clean_lines:
            fields = line.split(' ')
            # print(fields)
            charges_dict = {
                'atom':         fields[0], 
                'resnumbc':     fields[1], 
                'charge': float(fields[2])
            }
            # print(charges_dict)
            charges_dataset.append(charges_dict)
    return charges_dataset


charges = read_charges(chargesDB)

# %%


