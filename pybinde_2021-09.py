#!/usr/bin/python3
# version 2021.09.13
import pandas as pd
import sys
import os
import module_mmpbsa as bind
import time
import argparse
import freesasa as fs

initial_time = time.perf_counter()

#paths
#working_path   = "/test/"
lib_dir=os.path.dirname(os.path.realpath(__file__))
databases_path = lib_dir + "/databases/"

txt_df_res      = databases_path + "residues_rtp.txt"
txt_df_nb       = databases_path + "nonbonded_itp.txt"

# optional paths
delphi4py_path = os.path.dirname(os.path.realpath(__file__)) + "/delphi4py/"
#print("DELPHI PATH:",delphi4py_path)
sasa_classif_path = "classifier.config"


# distance cutoff
cutoff = False
cutoff_n = 4.5

# Verbose mode
verbose = False

# DelpPhi4Py Settings
# f_crg = lib_dir + "/delphi4py/example/DataBaseT.crg"
# f_siz = lib_dir + "/delphi4py/example/DataBaseT.siz"

scale = 2.5
nlit = 1000
nonit = 0
convergence = 0.001

# From here onward, PLEASE don't modify unless you are sure of the implications in doing so.
######################################################################

parser = argparse.ArgumentParser(description='====< Calculates MM-PBSA for a given .gro file >====')

parser.add_argument('gro',type=str,metavar='gro', help='input (.gro) file path')
parser.add_argument('atom_ranges',type=int, metavar='atom_ranges', help='Range of atom nrs that defines start and end of obj1 and obj2. Example: 1 1901 1902 1973', nargs=4)
parser.add_argument('-dbs','--databases',type=str, metavar='',required=False, help='(opt.) charges & radii databases path')
#parser.add_argument('-dbs','--databases',type=str, metavar='',required=False, default='./delphi4py/example/', help='(opt.) charges & radii databases path')
#parser.add_argument('-crg','--charges'  ,type=str, metavar='',required=False, default='./delphi4py/example/DataBaseT.crg', help='(opt.) charges database file path')
#parser.add_argument('-siz','--radii'    ,type=str, metavar='',required=False, default='./delphi4py/example/DataBaseT.siz', help='(opt.) radii database file path')
parser.add_argument('-sav','--saving_path',type=str, metavar='',required=False, default='./', help='(opt.) saving path for generated files')
parser.add_argument('-ensav','--energies_saving_path',type=str, metavar='',required=False, default='energies.txt', help='(opt.) saving path for generated energy files')
parser.add_argument('-ep','--epsin'    ,type=int, metavar='',required=False, default=4, help='(opt.) dielectric constant to be used')
parser.add_argument('-ter','--termini' ,type=str, metavar='',required=False, default=['NT3','CT4'], nargs=2, help='(opt.) termini to be used (Example: NT3 CT4)')

#group = parser.add_mutually_exclusive_group()
#group.add_argument('-q','--quiet', action='store_true',help='print quiet')
#group.add_argument('-v','--verbose', action='store_true',help='print verbose')
args = parser.parse_args()
#print(args)

if __name__ == '__main__':
    gro_file = args.gro
    epsin=args.epsin
    saving_path=args.saving_path
    if args.energies_saving_path:
      energy_summary_file_path=args.energies_saving_path
    else:
        energy_summary_file_path=args.energies_saving_path+'/energies_ep'+str(epsin)+'.txt'

    termini = args.termini

    monA_range, monB_range = [args.atom_ranges[0],args.atom_ranges[1]],[args.atom_ranges[2],args.atom_ranges[3]]

    if args.databases is None:
      databases='/delphi4py/example/'
      f_crg = lib_dir + databases + "DataBaseT.crg"
      f_siz = lib_dir + databases + "DataBaseT.siz"
      txt_df_res     =  "../databases/residues_rtp.txt"
      txt_df_nb      =  "../databases/nonbonded_itp.txt"
      classifier_path= "../databases/classifier.config"
    else:
      f_crg = args.databases + "/DataBaseT.crg"
      f_siz = args.databases + "/DataBaseT.siz"
      txt_df_res     = args.databases + "/residues_rtp.txt"
      txt_df_nb      = args.databases + "/nonbonded_itp.txt"
      classifier_path= args.databases + "/classifier.config"

    fpdb = bind.make_pdb_path(gro_file,saving_path)
    if verbose: print("fpdb:",fpdb)
    if verbose: print(os.getcwd())

    if verbose: print(os.path.dirname(os.path.realpath(__file__)))

    gro_df = bind.read_gro_df(gro_file,monA_range,monB_range) #create df with gro info, creates monomer index
    gro_df = bind.find_termini(gro_df)


    if verbose: print(gro_df.query("res_name == @termini[0]"))



    correct_protonation=True
    #bind.correct_protonation_state('LYSH','LY2','LY3','HZ3','HZ3',gro_df)
    #bind.correct_protonation_state('CYS2','CYS','CY0','HG','HG1',gro_df)
    
    #if termini are given (cpH) dont do
    if correct_protonation:
      # CORRECT_protonation_state(res_to_find,unprotonated_state, protonated_state,proton_name,new_proton_name,df)
      bind.correct_protonation_state('ASP','AS4','AS1','HD2','HD21',gro_df)
      bind.correct_protonation_state('ASPH','AS4','AS1','HD2','HD21',gro_df)
      bind.correct_protonation_state('CYS','CYS','CY0','HG','HG1',gro_df)
      bind.correct_protonation_state('CYSH','CYS','CY0','HG','HG1',gro_df)
      bind.correct_protonation_state('CYSC','CY3','CY3','HG','HG1',gro_df)
      bind.correct_protonation_state('CYS2','CYS','CYS','xx','xx',gro_df)

      bind.correct_protonation_state('GLU','GL4','GL1','HE2','HE21',gro_df)
      bind.correct_protonation_state('LYS','LY2','LY3','HZ3','HZ3',gro_df)
      bind.correct_protonation_state('LYSH','LY2','LY3','HZ3','HZ3',gro_df)
      bind.correct_protonation_state('TYR','TY2','TY0','HH','HH1',gro_df)

      bind.correct_protonation_state('CTR','CT4','CT1','HO','HO21',gro_df)
      bind.correct_protonation_state('NTR','NT0','NT3','H3','H3',gro_df)
      bind.rename_atom(gro_df,'NT5','H','H1')
      bind.correct_protonation_state('NTP','NT5','NT7','H2','H2',gro_df)

      #bind.correct_protonation_state('SER','SE3','SE0','HG','HG1',gro_df)
      #bind.correct_protonation_state('THR','TH3','TH0','HG1','HG1',gro_df)
      #bind.correct_protonation_state('HIS','HIS','AS2','HD2','HD21',gro_df)
      bind.correct_protonation_state_HIS(gro_df)
    #corr_atom_names=['HD11','HD21','HE12','HE22']
    #for x in ['HD1','HD2','HE1','HE2']: bind.rename_atom(gro_df,'HIS',x,'HD21')

    #Import databases
    df_res, df_charges, df_nb = bind.import_databases(txt_df_res,f_crg,txt_df_nb)

    #add information about atom_type
    gro_res_df = pd.merge(gro_df, df_res,  how = 'left', on = ['res_name', 'atom_name'])
    gro_res_df['atom_type'] = gro_res_df['atom_type'].str.replace('CH2R', 'CH2r')

    # add charges do gro df
    gro_res_crg_df = pd.merge(gro_res_df, df_charges,  how = 'left', on = ['res_name', 'atom_name'])
    if verbose: print(gro_res_crg_df)
    nans = gro_res_crg_df.query("atom_type != atom_type")
    nans_nr=gro_res_crg_df['atom_type'].isna().sum()
    if not nans.empty: print(nans.head(20)); print(nans.tail(20)); bind.prYellow("There are {} atoms with no attribution in atom_type!".format(nans_nr))


    #def guess_atom_type(res_num,og_df):
    #  mask1 = df["res_num"].startswith(og_df["res_num"])
    #  mask2 = df["atom_name"] == 'C'
    #  df.loc[mask1 & mask2,'atom_name'] = 'CT'



    monA_df = gro_res_crg_df.query('chain_id == "A"')
    trimmed_monA = pd.DataFrame(monA_df,columns=['atom_num','atom_type','charge','x_coord','y_coord','z_coord'])

    monB_df = gro_res_crg_df.query('chain_id == "B"')
    trimmed_monB = pd.DataFrame(monB_df,columns=['atom_num','atom_type','charge','x_coord','y_coord','z_coord'])

    # calculate distances between 2 mons
    if verbose: print("Calculating distances...")

    dist_df = bind.calc_cdist(trimmed_monA, trimmed_monB)
    if cutoff: dist_df = dist_df.query("distance <= {}".format(cutoff_n))

     # fix some names between gro and databases
    complete_df = pd.merge(dist_df, df_nb,  how = 'left', on = ['type_Ai', 'type_Bj'])
    if verbose: print(complete_df)
    nans = complete_df.query("charge_Ai != charge_Ai")
    nans_nr=complete_df['charge_Ai'].isna().sum()
    if not nans.empty: print(nans.head(20)); print(nans.tail(20)); print(pd.unique(nans['atom_Ai'])); bind.prYellow("There are {} atoms with no charge!".format(nans_nr))
    # Vdw and coulomb calculations
    complete_df['En_VdW'] = (complete_df['c12']/complete_df['distance']**12)-(complete_df['c6']/complete_df['distance']**6)

    f = 138.935458
    complete_df['En_Coul'] = f*(complete_df['charge_Ai'].astype(float)*complete_df['charge_Bj'].astype(float))/(epsin*complete_df['distance'].astype(float))

    VdW_en_total  = complete_df['En_VdW'].sum()
    Coul_en_total = complete_df['En_Coul'].sum()
    if verbose: print(complete_df); print("En_VdW = ", VdW_en_total); print("En_Coul = ", Coul_en_total)


    nans = complete_df.query("charge_Bj != charge_Bj")
    if not nans.empty:
        print("empty entries in complete_df:",nans)

    bind.gro2pdb_simple(gro_df,fpdb,termini)
    if verbose: bind.prYellow("New file is: "+fpdb)


    fs.setVerbosity(fs.nowarnings)

    area_prot, area_monA, area_monB=bind.run_freesasa_custom(fpdb,100,classifier_path,verbose=True) # (fpdb,classifier_path) "/home/joaov/Projetos/MMPBSA/validation/PyEBind_validation/classifier.config"

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

    pdb_coord = bind.read_pdb(fpdb)
    pdb_df = pd.DataFrame(pdb_coord)


    x,y,z=bind.geom_center(pdb_df)
    acent=[x,y,z]

    gsize=bind.max_gsize(pdb_df, scale)

    if verbose: print("Geometric center of dimer in pdb is:",acent)


    box_size=bind.appropriate_box_size(pdb_df)
    if verbose: print("Appropriate box side size of pdb is:",box_size)
    if verbose: print("gsize = ", gsize)

    # Run DelPhi4Py
    monomers = monA_range[0], monA_range[1], monB_range[0], monB_range[1]

    #if __name__ == "__main__":
    delphi_settings = [f_crg, f_siz, fpdb, scale, gsize, nlit, nonit, convergence, acent, epsin]

    solvation_dimer, solvation_mon1, solvation_mon2 = bind.d4p_run_all(monomers,delphi_settings,delphi4py_path)

    if verbose: print(solvation_dimer); print(solvation_mon1); print(solvation_mon2)

    nonbonded   = [VdW_en_total,Coul_en_total]
    sasa        = [area_prot,area_monA,area_monB]
    polar       = [solvation_dimer, solvation_mon1, solvation_mon2]

    if verbose: print(nonbonded, sasa, polar)

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
    
    final_time = time.perf_counter()
    elapsed_time=final_time-initial_time
    if verbose:
        print("Elapsed time: ", elapsed_time,"seconds.")

    #VdW_en_total, Coul_en_total, Gnonpolar, Gpolar = run_protocol(files,monomers,settings)
    binding_energy= VdW_en_total + Coul_en_total + Gnonpolar + Gpolar
    with open(energy_summary_file_path,'a+') as outfile:
        outfile.write('{:<20} {:<20} {:<20} {:<20} {:<20} {:<20} {:<20}\n'.format(bind.get_file_name(sys.argv[1]), VdW_en_total, Coul_en_total, Gnonpolar, Gpolar, binding_energy, round(elapsed_time,2)))
        
    print("File:                              ",bind.get_file_name(sys.argv[1]))
    print("VdW Energy (kJ/mol):               ", VdW_en_total)
    print("Coul Energy (kJ/mol):              ", Coul_en_total)
    print("Nonpolar Solvation Energy (kJ/mol):", Gnonpolar)
    print("Polar Solvation Energy (kJ/mol):   ", Gpolar)
    print("Binding Energy (kJ/mol):           ", binding_energy)
    print("Running Time (s):                  ", round(elapsed_time,2))
    print()
