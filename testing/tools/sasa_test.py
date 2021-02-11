import freesasa as fs
import pyjoaov as pj
import os
import time
initial_time = time.perf_counter()
#initial_time = (initial_time.strftime('%d-%m-%Y %H:%M:%S'))

#paths
gro_file="/home/joaov/github/mmpbsa/testing/frame00100.gro"
pdb_file="/home/joaov/github/mmpbsa/testing/saved_files/frame00100.pdb" #frame00100.pdb"

working_path   = "/home/joaov/github/mmpbsa/testing/"
gro_file        = working_path + "frame00100.gro" #"short-gro.gro"
saving_path     = working_path + "saved_files/" #save location for new .pdb from .gro

file_basename   = os.path.basename(gro_file)
file_name, file_extension = os.path.splitext(file_basename)
new_filepath    = saving_path + file_name + ".pdb"

# number of the atom where the monomers start and end
mon1_start = 1
mon1_end   = 1149
mon2_start = 1150
mon2_end   = 2298

gro_df = pj.read_gro_minus_index(gro_file) #create df with gro info, creates monomer index

monA_range = [mon1_start, mon1_end]
monB_range = [mon2_start, mon2_end]

gro_df = pj.make_index_and_subset(gro_df, monA_range, monB_range)

monA_df = gro_df.query('chain_id == "A"') #creates a dataframe with only monA
monB_df = gro_df.query('chain_id == "B"') #creates a dataframe with only monB

prot = saving_path+"Prot.pdb"
monA = saving_path+"monA.pdb"
monB = saving_path+"monB.pdb"

pj.gro2pdb(gro_df,prot)
pj.gro2pdb(monA_df,monA)
pj.gro2pdb(monB_df,monB)

fs.setVerbosity(1)
#print(fs.getVerbosity())
#print(fs.probeRadius())

def run_freesasa(file):
    structure = fs.Structure(file)
    result =fs.calc(structure,fs.Parameters({'algorithm' : fs.ShrakeRupley,
                                             'probe-radius' : 1.4,
                                              'n-points' : 1000}))
    #result =fs.calc(structure,fs.Parameters({'algorithm' : fs.LeeRichards,
    #                                               'n-slices' : 100}))
    area =result.totalArea()
    #area_classes = fs.classifyResults(result, structure)
    #print("Total : %.2f A2" % result.totalArea())
    #for key in area_classes:
    #    print(key, ": %.2f A2" % area_classes[key])
    #    
    g = 0.00542 # gamma -> kcal/(mol‚A2)
    b = 0.92    # beta  -> kcal/mol
    energy=result.totalArea()*g+b
    return area, energy

area_prot, energy_prot = run_freesasa(prot)
area_monA, energy_monA = run_freesasa(monA)
area_monB, energy_Prot = run_freesasa(monB)

print("Area prot:",area_prot,"En:", energy_prot*4.184)
print("Area monA:",area_monA,"En:", energy_monA*4.184)
print("Area monB:",area_monB,"En:", energy_Prot*4.184)

print("Energy =",(energy_prot-(energy_monA+energy_Prot))*4.184)
final_time = time.perf_counter()
#final_time = (final_time.strftime('%d-%m-%Y %H:%M:%S'))

elapsed_time=final_time-initial_time
print("Elapsed time: ", elapsed_time,"seconds.")

