from os import read
import re
import pandas as pd
from biopandas.pdb import PandasPdb
import numpy as np
import subprocess
# files
pdb_file = "/home/joaov/python-mmpbsa/mmpbsa/testing/dimer.pdb"
gro_file = "/home/joaov/python-mmpbsa/mmpbsa/testing/frame00100.gro"
res_conv = "/home/joaov/python-mmpbsa/mmpbsa/testing/res_atoms_db.rtp"
nonbonded_itp = "/home/joaov/python-mmpbsa/mmpbsa/testing/nb.itp" #use new_ instead of nb.itp (incomplete)
atomtypes_itp = "/home/joaov/python-mmpbsa/mmpbsa/testing/atomtypes_nb.itp"
charges_db = "/home/joaov/python-mmpbsa/mmpbsa/testing/DataBaseT.crg"
# number of the atom where the monomer 2 starts
mon2_start = 1150
# functions


def read_gro(file):
    gro_lines = []
    gro_dict = []
    with open(file) as gro:
        next(gro)
        next(gro)
        for line in gro:
            splitted_line = [line[0:5], line[5:9], line[10:15], line[16:20],
                             line[21:28], line[29:36], line[37:45]]
            gro_lines.append(splitted_line)
        for line in gro_lines:
            line_dict = {
                'res_num':    line[0].strip(),
                'res_name':   line[1].strip(),
                'atom_type':    line[2].strip(),
                'atom_num':    line[3].strip(),
                'x_coord':    (line[4].strip()),
                'y_coord':    (line[5].strip()),
                'z_coord':    (line[6].strip())
            }
            gro_dict.append(line_dict)
    gro_df = pd.DataFrame(gro_dict)
    gro_df.drop(gro_df.tail(1).index, inplace=True)
    monX = []
    for i in range(1, int(mon2_start)):
        monX.append("A")
    for i in range(1, int(mon2_start)):
        monX.append("B")
    gro_df['chain_id'] = monX
    return gro_df


def terminus_fix(gro_dataframe):
    minval = str(gro_dataframe['res_num'].min())  # 1
    subset_mingro = gro_dataframe.query(
        "res_num == @minval")  # only entries for res 1
    first_res = subset_mingro['res_name'].iloc[0]  # ILE
    subset_mingro['res_name'] = subset_mingro['res_name'].replace(
        first_res, 'NTR')  # replace all ILE with NTR
    # Grab only the NTR entries from df_res
    nt_res = df_res.query("res_name == 'NTR'")
    subset_mingro['copy_index'] = subset_mingro.index
    mingro_merge = pd.merge(subset_mingro, nt_res,  how='left', left_on=['res_name', 'atom_type'], right_on=[
                            'res_name', 'atom_type']).set_index('copy_index')  # merge dataframes and output a new column with atom names
    min_nans = mingro_merge.query("atom_name != atom_name")
    min_nans.loc[min_nans.atom_name != 'atom_name', 'res_name'] = first_res
    mingro_merge.update(min_nans)

    maxval = str(gro_dataframe['res_num'].max())  # 99
    subset_maxgro = gro_dataframe.query("res_num == @maxval")
    last_res = subset_maxgro['res_name'].iloc[-1]  # MET
    subset_maxgro['res_name'] = subset_maxgro['res_name'].replace(
        last_res, 'CTR')
    subset_maxgro['atom_type'] = subset_maxgro['atom_type'].replace(
        'O2', 'OT2')
    subset_maxgro['atom_type'] = subset_maxgro['atom_type'].replace(
        'O1', 'OT1')
    subset_maxgro['atom_type'] = subset_maxgro['atom_type'].replace('C', 'CT')
    ct_res = df_res.query("res_name == 'CTR'")
    # saves a copy from og index
    subset_maxgro['copy_index'] = subset_maxgro.index
    maxgro_merge = pd.merge(subset_maxgro, ct_res,  how='left', left_on=[
                            'res_name', 'atom_type'], right_on=['res_name', 'atom_type']).set_index('copy_index')
    max_nans = maxgro_merge.query("atom_name != atom_name")
    max_nans.loc[max_nans.atom_name != 'atom_name', 'res_name'] = last_res
    maxgro_merge.update(max_nans)
    return mingro_merge.append(maxgro_merge)


# with open(log_dist, "a") as log:
    #log.write(('{0:5} {1:5} {2:5}'.format("num_a", "num_b", "distance"+"\n")))
def calc_dists(mon1,mon2):
    distances_db = []
    for a in mon1.index:
        xa = mon1['x_coord'][a]
        ya = mon1['y_coord'][a]
        za = mon1['z_coord'][a]
        num_a = mon1['atom_num'][a]
        name_a = mon1['atom_name'][a]
        a_coord = np.array((xa, ya, za), dtype=float)
        for b in mon2.index:
            xb = mon2['x_coord'][b]
            yb = mon2['y_coord'][b]
            zb = mon2['z_coord'][b]
            num_b = mon2['atom_num'][b]
            name_b = mon2['atom_name'][b]
            b_coord = np.array((xb, yb, zb), dtype=float)
            distance = round(np.linalg.norm(b_coord - a_coord), 3)
            if cutoff:
                if distance <= cutoff_n:
                    dist_dict = {
                        'atom_Ai':      num_a,
                        'atom_Bj':      num_b,
                        'name_Ai':      name_a,
                        'name_Bj':      name_b,
                        'distance':     distance
                    }
                    distances_db.append(dist_dict)
                else:
                    continue
            else:
                dist_dict = {
                    'atom_Ai':      num_a,
                    'atom_Bj':      num_b,
                    'name_Ai':      name_a,
                    'name_Bj':      name_b,
                    'distance':     distance
                }
                distances_db.append(dist_dict)
    return distances_db

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
                    'atom_type':        fields[0].rstrip(),
                    'atom_name':        fields[1].rstrip()
                    #'charge':     float(fields[2])
                }
                residues.append(res_dict)
    return residues

def read_nb(filenb,atomtypes):
    NB = []
    lines_no_comments = []
    clean_lines = []
    with open(atomtypes) as atypes:
        next(atypes)
        next(atypes)
        lines = atypes.read().splitlines()
        for line in lines:
            sep_comment = line.split(";")
            lines_no_comments.append(sep_comment[0])
        for line in lines_no_comments:
            cleanline = re.sub('\s+', ' ', line)
            clean_lines.append(cleanline.strip())
        clean_lines = list(filter(None, clean_lines))
        for line in clean_lines:
            fields = line.split(' ')
            nb_dict = {
                'i':          fields[0].rstrip(), 
                'j':          fields[0].rstrip(), 
                'c6':   float(fields[5]), 
                'c12':  float(fields[6])
            }
            NB.append(nb_dict)
    lines_no_comments = []
    clean_lines = []
    with open(filenb) as nb:
        next(nb)
        next(nb)
        lines = nb.read().splitlines()
        for line in lines:
            sep_comment = line.split(";")
            lines_no_comments.append(sep_comment[0])
        for line in lines_no_comments:
            cleanline = re.sub('\s+', ' ', line)
            clean_lines.append(cleanline.strip())
        clean_lines = list(filter(None, clean_lines))
        for line in clean_lines:
            fields = line.split(' ')
            nb_dict = {
                'i':          fields[0].rstrip(), 
                'j':          fields[1].rstrip(), 
                'c6':   float(fields[3]), 
                'c12':  float(fields[4])
            }
            NB.append(nb_dict)
    return NB

def read_charges(file):
    clean_lines = []
    charges_dataset = []
    with open(file) as db:
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

# start of pipeline
gro_df = read_gro(gro_file) #create df with gro info, creates monomer index

monA = gro_df.query('chain_id == "A"') #creates a dataframe with only monA
monB = gro_df.query('chain_id == "B"') #creates a dataframe with only monB

df_res = pd.DataFrame(read_resdb(res_conv)) # df w/ atom names/types from rtp

mega_monA = terminus_fix(monA) # grabs and fixes terminus, CTR and NTR added
mega_monB = terminus_fix(monB) # grabs and fixes terminus, CTR and NTR added

# merge to get a column with atom names
merged_df = pd.merge(gro_df, df_res,  how='left', left_on=[
                     'res_name', 'atom_type'], right_on=['res_name', 'atom_type'])

merged_df.update(mega_monA) # use previously fixed terminus to update merged_df
merged_df.update(mega_monB) # use previously fixed terminus to update merged_df

charges = pd.DataFrame(read_charges(charges_db)) #create charges df
charges = charges.rename(columns={"atom": "atom_type", "resnumbc": "res_name"})
#add chages do gro df
complete_df = pd.merge(merged_df, charges,  how='left', left_on=[
                     'res_name', 'atom_type'], right_on=['res_name', 'atom_type'])

nb_df=pd.DataFrame(read_nb(nonbonded_itp,atomtypes_itp)) # df from nb.itp to grab c6 and c12
alt_nb_df = nb_df.rename(columns={"i": "j", "j": "i"}) # create switched refs
complete_nb_df=nb_df.append(alt_nb_df) # now we have the reverse combinations too

complete_nb_df=complete_nb_df.rename(columns={"i": "name_Ai", "j": "name_Bj"})
complete_nb_df=complete_nb_df.drop_duplicates() # [3364 rows x 4 columns]

complete_monA = complete_df.query('chain_id == "A"') 
complete_monB = complete_df.query('chain_id == "B"') 

cutoff = True
cutoff_n = 5
dist_df = pd.DataFrame(calc_dists(complete_monA,complete_monB)) # calculate distances between 2 mons

'''log_dist = "/home/joaov/python-mmpbsa/mmpbsa/testing/pdb_dist.log"
c = subprocess.Popen(["/bin/rm", '-f', log_dist])
c.communicate(b'\n')
# np.savetxt(r'/home/joaov/python-mmpbsa/mmpbsa/testing/test_dists.log', dist_df.values, fmt='%s', delimiter='\t')'''

#fix some names between gro and databases
complete_nb_df['name_Ai'] = complete_nb_df['name_Ai'].str.replace('CH2R','CH2r')
complete_nb_df['name_Bj'] = complete_nb_df['name_Bj'].str.replace('CH2R','CH2r')
dist_df['name_Ai'] = dist_df['name_Ai'].str.replace('CH2R','CH2r')
dist_df['name_Bj'] = dist_df['name_Bj'].str.replace('CH2R','CH2r')

copy_nb = complete_nb_df
copy_dist = dist_df

merged_copies = pd.merge(copy_dist, copy_nb,  how='left', left_on=[
                     'name_Ai', 'name_Bj'], right_on=['name_Ai', 'name_Bj'])

nan_merged_copies = merged_copies.query("c6 != c6") #none

merged_copies['En_VdW'] = (merged_copies['c12']/merged_copies['distance']**12)-(merged_copies['c6']/merged_copies['distance']**6)

#Vcoulomb= f* (q1*q2)/(Er*rij)
#the charges come from databaseT.crg?
#where do I get Er?78 for water?
#rij "separation" is distance?
f=138.935458
Er=78
merged_copies['En_Coul'] = f*(merged_copies['charge_Ai']*merged_copies['charge_Bj'])/(Er*merged_copies['distance'])

#print(merged_copies)
VdW_en_total = merged_copies['En_VdW'].sum()
Coul_en_total = merged_copies['En_Coul'].sum()
print(VdW_en_total)
print(Coul_en_total)
#print(VdW_en_total)