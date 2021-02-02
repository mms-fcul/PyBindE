import pandas as pd
import numpy as np
#from Bio.PDB import PDBParser, ShrakeRupley
import freesasa as fs
import pickle

working_path = "/home/joaov/github/mmpbsa/testing/"
databases_path = "/home/joaov/github/mmpbsa/testing/saved_dataframes/"

pdb_file        = "/home/joaov/github/mmpbsa/testing/saved_files/frame00100.pdb"#working_path + "dimer.pdb"
gro_file        = working_path + "frame00100.gro" #"short-gro.gro"

txt_df_res      = databases_path + "residues_rtp.txt"
txt_df_nb       = databases_path + "nonbonded_itp.txt"
txt_df_charges  = databases_path + "charges.txt"


# number of the atom where the monomers start and end
mon1_start = 1
mon1_end   = 1149
mon2_start = 1150
mon2_end   = 2298

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
            #print (splitted_line)
        for line in gro_lines:
            line_dict = {
                'res_num':     line[0].strip(), 
                'res_name':    line[1].strip(), 
                'atom_name':   line[2].strip(), 
                'atom_num':    line[3].strip(), 
                'x_coord':    (line[4].strip()), 
                'y_coord':    (line[5].strip()), 
                'z_coord':    (line[6].strip())
            }
            gro_dict.append(line_dict)
    gro_df = pd.DataFrame(gro_dict)
    gro_df.drop(gro_df.tail(1).index, inplace = True)
    monX = []
    for i in range (int(mon1_start), int(mon1_end)+1):
        monX.append("A")
    for i in range (int(mon2_start), int(mon2_end)+1):
        monX.append("B")
    gro_df['chain_id'] = monX
    return gro_df

# start of pipeline
gro_df = read_gro(gro_file) #create df with gro info, creates monomer index

monA = gro_df.query('chain_id == "A"') #creates a dataframe with only monA
monB = gro_df.query('chain_id == "B"') #creates a dataframe with only monB
#print(gro_df)

#parser = PDBParser(QUIET=1)
#structure = parser.get_structure("dimer", pdb_file)
##print(structure)
#sr = ShrakeRupley()
#sr.compute(structure, level="S") 
#print(round(structure.sasa, 2))


structure = fs.Structure(pdb_file)
result =fs.calc(structure)
#area_classes = fs.classifyResults(result, structure)
print("Total : %.2f A2" % result.totalArea())
#for key in area_classes:
#    print(key, ": %.2f A2" % area_classes[key])
#    
selections = fs.selectArea(('all, resi 1-198','chainA, resi 1-99', 'chainB, resi 100-198'), structure, result)
list_keys=[]
for key in selections:
    print (key, ": %.2f A2" % selections[key])
    list_keys.append(selections[key])
    # Saving the objects:

with open('sasa_vals.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
    pickle.dump(list_keys, f)

