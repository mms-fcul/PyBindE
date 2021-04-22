#!/usr/bin/python3



#################################################################################
## --This protocol reads the gro_file and all databases, matches the required  ##
## info in databases to the atoms and residues on the gro_file, calculates     ##
## distances between the atoms of the 2 monomers, calculates VdW and Coulomb   ## 
## energies (EMM).------------------------------------------------------------ ##
## --After that, calculates the SASA for the dimer, mon1 and  mon2 that will   ##
## be used to calculate Gnonpolar. For this step, freesasa (with a convertion  ##
## from gro to pdb) will be used.--------------------------------------------- ##
## --Next, it sets up a few variables needed for the Gpolar calculation,       ##
## using DelPhi4Py, and generates a summary file for the 4 energy terms        ##
## that contribute to Gbinding.----------------------------------------------- ##                         
#################################################################################



import pandas as pd
import sys
sys.path.insert(0, './')
import mmpbsa_module as pj
import time
initial_time = time.perf_counter()

#paths
#working_path   = "/test/"
databases_path = "./databases/"

txt_df_res      = databases_path + "residues_rtp.txt"
txt_df_nb       = databases_path + "nonbonded_itp.txt"
txt_df_charges  = databases_path + "charges.txt"

# optional paths
delphi4py_path = "./delphi4github/delphi4py/" 
sasa_classif_path = "classifier.config"

gro_file        =  (sys.argv[1])

saving_path     = './saved_files/'  #save location for new .pdb from .gro
energy_summary_file_path='./saved_files/energy_summary.txt'
# number of the atom where the monomers start and end
mon1_start = 1
mon1_end   = 1149
mon2_start = 1150
mon2_end   = 2298

# termini to use
termini = ['NT3','CT4'] 

# distance cutoff
cutoff = False
cutoff_n = 4.5

# Verbose mode
verbose = True

# DelpPhi4Py Settings
f_crg = "./delphi4github/delphi4py/example/DataBaseT.crg"
f_siz = "./delphi4github/delphi4py/example/DataBaseT.siz"
fpdb = pj.make_pdb_path(gro_file,saving_path)
epsin=4

scale = 2
gsize = 181         # it's going to be calculated and overwritten
acent=[]            # it's going to be calculated and overwritten
nlit = 500
nonit = 50
convergence = 0.01

# From here onward, PLEASE don't modify unless you are sure of the implications in doing so.
######################################################################
databases       = [txt_df_res,txt_df_nb,txt_df_charges] # residues DB, Nonbonded DB, Charges DB
files           = [gro_file, databases, saving_path, delphi4py_path, sasa_classif_path] #gro_file, list of databases, pdb file to create, delphy4py path, configs file for sasa
monomers        = [mon1_start,mon1_end,mon2_start,mon2_end]

cutoff_settings = [cutoff,cutoff_n]
delphi_settings = [f_crg, f_siz, fpdb, scale, gsize, nlit, nonit, convergence, acent, epsin] # 10 arguments
settings        = [cutoff_settings,delphi_settings,termini,verbose]
######################################################################

def mmpbsa_protocol(files,monomers,settings):
    
    gro_file, databases, saving_path, delphi4py_path, sasa_classif_path                 = files
    txt_df_res, txt_df_nb, txt_df_charges                                               = databases
    mon1_start, mon1_end, mon2_start, mon2_end                                          = monomers

    cutoff_settings, delphi_settings, termini, verbose                                  = settings
    cutoff, cutoff_n                                                                    = cutoff_settings
    f_crg, f_siz, fpdb, scale, gsize, nlit, nonit, convergence, acent, epsin            = delphi_settings #9

    gro_df = pj.read_gro_df(gro_file) #create df with gro info, creates monomer index

    monA_range = [mon1_start, mon1_end]
    monB_range = [mon2_start, mon2_end]

    gro_df = pj.find_termini(gro_df, termini)
    gro_df = pj.make_index_and_subset(gro_df,monA_range,monB_range)

    gro_df['atom_name'] = gro_df['atom_name'].replace('O2', 'OT2')
    gro_df['atom_name'] = gro_df['atom_name'].replace('O1', 'OT1')

    #Import databases
    df_res, df_charges, df_nb = pj.import_databases(txt_df_res,txt_df_charges,txt_df_nb)

    #add information about atom_type
    gro_res_df = pd.merge(gro_df, df_res,  how = 'left', on = ['res_name', 'atom_name'])
    gro_res_df['atom_type'] = gro_res_df['atom_type'].str.replace('CH2R', 'CH2r')
    gro_res_df = pj.change_res(gro_res_df, 'NTR', termini[0])
    gro_res_df = pj.change_res(gro_res_df, 'CTR', termini[1])

    # add charges do gro df
    gro_res_crg_df = pd.merge(gro_res_df, df_charges,  how = 'left', on = ['res_name', 'atom_name'])

    nans = gro_res_crg_df.query("atom_type != atom_type")
    if not nans.empty: print(nans.head(20)); print(nans.tail(20)); pj.prYellow("There are atoms with no attribution in atom_type!")

    monomers_df = gro_res_crg_df

    monA_df = monomers_df.query('chain_id == "A"')
    trimmed_monA = pd.DataFrame(monA_df,columns=['atom_num','atom_type','charge','x_coord','y_coord','z_coord'])

    monB_df = monomers_df.query('chain_id == "B"')
    trimmed_monB = pd.DataFrame(monB_df,columns=['atom_num','atom_type','charge','x_coord','y_coord','z_coord'])

    # calculate distances between 2 mons
    if verbose: print("Calculating distances...")

    dist_df = pj.calc_cdist(trimmed_monA, trimmed_monB)
    if cutoff: dist_df = dist_df.query("distance <= {}".format(cutoff_n))

    # fix some names between gro and databases
    complete_df = pd.merge(dist_df, df_nb,  how = 'left', on = ['type_Ai', 'type_Bj'])

    #Vdw and coulomb calculations
    complete_df['En_VdW'] = (complete_df['c12']/complete_df['distance']**12)-(complete_df['c6']/complete_df['distance']**6)

    f = 138.935458
    complete_df['En_Coul'] = f*(complete_df['charge_Ai']*complete_df['charge_Bj'])/(epsin*complete_df['distance'])

    VdW_en_total = complete_df['En_VdW'].sum()
    Coul_en_total = complete_df['En_Coul'].sum()
    if verbose: print(complete_df); print("En_VdW = ", VdW_en_total); print("En_Coul = ", Coul_en_total)

    nans = complete_df.query("charge_Bj != charge_Bj")
    if not nans.empty:
        print("empty entries in complete_df:",nans)
    
    gro_df['atom_name'] = gro_df['atom_name'].replace('OT2', 'O2')
    gro_df['atom_name'] = gro_df['atom_name'].replace('OT1', 'O1')
    pj.gro2pdb_simple(gro_df,fpdb,termini)
    if verbose: pj.prYellow("New file is: "+fpdb)

    import freesasa as fs
    fs.setVerbosity(1)
    
    area_prot, area_monA, area_monB=pj.run_freesasa_custom(fpdb, verbose=True) # (fpdb,classifier_path)
    
    g = 0.00542 # gamma -> kcal/(mol‚A2)
    b = 0.92    # beta  -> kcal/mol 
    
    energy_prot = (area_prot*g+b)
    energy_monA = (area_monA*g+b)
    energy_monB = (area_monB*g+b)
    
    if verbose: 
        print("Area prot:",area_prot,"En:", energy_prot*4.184)
        print("Area monA:",area_monA,"En:", energy_monA*4.184)
        print("Area monB:",area_monB,"En:", energy_monB*4.184)

        print("Energy =",(energy_prot-(energy_monA+energy_monB))*4.184)
    
    pdb_coord = pj.read_pdb(fpdb)
    pdb_df = pd.DataFrame(pdb_coord)

    x,y,z=pj.geom_center(pdb_df)
    acent=[x,y,z]
    delphi_settings[8]=acent

    gsize=pj.max_gsize(pdb_df)
    delphi_settings[4]=gsize
    if verbose: print("Geometric center of dimer in pdb is:",acent)
    print(delphi_settings)
    box_size=pj.appropriate_box_size(pdb_df)
    if verbose: print("Appropriate box side size of pdb is:",box_size)
    print("gsize = ", gsize)
  
    # Run DelPhi4Py
    monomers = mon1_start, mon1_end, mon2_start, mon2_end

    #if __name__ == "__main__":
    solvation_dimer, solvation_mon1, solvation_mon2 = pj.d4p_run_all(monomers,delphi_settings,delphi4py_path)

    if verbose: print(solvation_dimer); print(solvation_mon1); print(solvation_mon2)
    
    nonbonded   = [VdW_en_total,Coul_en_total]
    sasa        = [area_prot,area_monA,area_monB]
    polar       = [solvation_dimer, solvation_mon1, solvation_mon2]
    
    return nonbonded, sasa, polar


def run_protocol(files,monomers,settings):
    nonbonded, sasa, polar = mmpbsa_protocol(files,monomers,settings)

    g = 0.00542 # gamma -> kcal/(mol‚A2)
    b = 0.92    # beta  -> kcal/mol
    #nonbonded = [VdW_en_total,Coul_en_total]
    #sasa      = [sasa_P,sasa_MA,sasa_MB]
    #polar     = [solv_dimer, solv_mon1, solv_mon2]

    #EMM --------------------------------------
    Ebonded = 0
    VdW_en_total  = nonbonded[0]
    Coul_en_total = nonbonded[1]
    EMM = Ebonded + VdW_en_total + Coul_en_total
    #Gnonpolar --------------------------------
    A_complex = float(sasa[0])
    A_mon1    = float(sasa[1])
    A_mon2    = float(sasa[2])

    Gnonpolar_complex = g*A_complex + b
    Gnonpolar_mon1    = g*A_mon1    + b
    Gnonpolar_mon2    = g*A_mon2    + b

    Gnonpolar = Gnonpolar_complex - (Gnonpolar_mon1 + Gnonpolar_mon2)
    Gnonpolar = Gnonpolar*4.184 # kJ/mol

    #Gpolar -----------------------------------
    t=310 #K
    kb=0.008314463
    Gpolar_complex = polar[0]*kb*t
    Gpolar_mon1    = polar[1]*kb*t
    Gpolar_mon2    = polar[2]*kb*t
    Gpolar = Gpolar_complex - (Gpolar_mon1 + Gpolar_mon2)
    if verbose: print("gpolar",Gpolar) # kJ/(mol.K)

    Gsolv = Gpolar + Gnonpolar
    if verbose: print("solv",Gsolv)
    G_binding = EMM + Gsolv

    return VdW_en_total, Coul_en_total, Gnonpolar, Gpolar

VdW_en_total, Coul_en_total, Gnonpolar, Gpolar = run_protocol(files,monomers,settings)
binding_energy= VdW_en_total + Coul_en_total + Gnonpolar + Gpolar
with open(energy_summary_file_path,'a+') as outfile:
    outfile.write('{:<20} {:<20} {:<20} {:<20} {:<20} {:<20}\n'.format(pj.get_file_name(sys.argv[1]), VdW_en_total, Coul_en_total, Gnonpolar, Gpolar, binding_energy)) 

if verbose:
    final_time = time.perf_counter()
    elapsed_time=final_time-initial_time
    print("Elapsed time: ", elapsed_time,"seconds.")


