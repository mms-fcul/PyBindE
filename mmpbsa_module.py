
import pandas as pd
import numpy as np
import freesasa as fs
import os
from copy import copy
import datetime
from scipy.spatial.distance import cdist

def files_in_path(path):
  files=[]
  #path="/home/joaov/Projetos/MMPBSA/dim005/tmp"
  for filename in os.listdir(path):
      if filename.endswith(".gro"):
          files.append(os.path.join(path, filename))
  return files

def read_gro_df(file):
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
  return gro_df

def find_termini(df, termini=['NT3','CT4']):
  basic_aa_list=[
  "ALA",
  "CYS",'CY0', 'CY1', 'CY2', 'CY3',
  "ASP",'AS0', 'AS1', 'AS2', 'AS3', 'AS4',
  "GLU",'GL0', 'GL1', 'GL2', 'GL3', 'GL4',
  "PHE",
  "GLY",
  "HIS",'HI0', 'HI1', 'HI2',
  "ILE",
  "LYS",'LY0', 'LY1', 'LY2', 'LY3',
  "LEU",
  "MET",
  "ASN",
  "PRO",
  "GLN",
  "ARG",
  "SER",'SE0', 'SE1', 'SE2', 'SE3'
  "THR",'TH0', 'TH1', 'TH2', 'TH3',
  "VAL",
  "TRP",
  "TYR",'TY0', 'TY1', 'TY2',
  "NTR",'NT0', 'NT1', 'NT2', 'NT3',
  "CTR",'CT0', 'CT1', 'CT2', 'CT3', 'CT4']


  df['res_num']=df['res_num'].astype(int)

  ntr_tags=['H1','H2','H3']
  
  mask1 = df["res_name"].isin(basic_aa_list)
  mask2 = df["atom_name"].isin(ntr_tags)
  ntrs=pd.unique(df.loc[mask1 & mask2,'res_num'])

  ctr_tags=['OXT','O1','O2','OT1','OT2']

  mask2 = df["atom_name"].isin(ctr_tags)
  ctrs=pd.unique(df.loc[mask1 & mask2,'res_num'])


  ntr_list=['N','H1','H2','H3','CA','C','O']
  mask1 = df["res_num"].isin(ntrs)
  mask2 = df["atom_name"].isin(ntr_list)
  df.loc[mask1 & mask2,'res_name'] = 'NTR'

  mask1 = df["res_num"].isin(ctrs)
  mask2 = df["atom_name"] == 'C'
  df.loc[mask1 & mask2,'atom_name'] = 'CT'

  ctr_list=['CT','OT1','OT2','O1','O2','HO11','HO21','HO12','HO22']
  mask2 = df["atom_name"].isin(ctr_list)
  df.loc[mask1 & mask2,'res_name'] = 'CTR'

  print(ntrs,ctrs)
  alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
  for num,terminus in enumerate(zip(ntrs,ctrs)):
    mask = df['res_num'].between(terminus[0],terminus[1])
    df.loc[mask, 'real_chain_id'] = alphabet[num]
    print(num,alphabet[num])
  return df

def make_index_and_subset(dataframe,monA_range,monB_range):
  df_copy = dataframe.copy()
  df_copy["atom_num"] = pd.to_numeric(df_copy["atom_num"])
  monA_start = monA_range[0]
  monA_end   = monA_range[1]
  mask = df_copy['atom_num'].between(monA_start,monA_end)
  df_copy.loc[mask, 'chain_id'] = 'A'

  monB_start = monB_range[0]
  monB_end   = monB_range[1]
  mask = df_copy['atom_num'].between(monB_start,monB_end)
  df_copy.loc[mask, 'chain_id'] = 'B'
  df_copy_subset = df_copy.query("chain_id == chain_id")
  return df_copy_subset

def import_databases(txt_df_res,txt_df_charges,txt_df_nb):
  df_res      = pd.read_table(txt_df_res,     delim_whitespace=True,header=None,names=['res_name','atom_name','atom_type']) # df w/ atom names/types from rtp
  df_charges  = pd.read_table(txt_df_charges, delim_whitespace=True,header=None,names=['atom_name','res_name','charge'])
  df_nb       = pd.read_table(txt_df_nb,      delim_whitespace=True,header=None,names=['i','j','c6','c12']) # df from nb.itp to grab c6 and c12

  #print(df_nb[df_nb.duplicated()])
  df_nb = df_nb.rename(columns = {"i": "type_Ai", "j": "type_Bj"})
  return df_res, df_charges, df_nb

def calc_cdist(set1, set2):
  a = set1[["x_coord", "y_coord", "z_coord"]]
  b = set2[["x_coord", "y_coord", "z_coord"]]

  df = set1[["atom_num", "atom_type","charge"]].merge(set2[["atom_num", "atom_type","charge"]], how="cross", suffixes=("_Ai", "_Bj"))
  df["distance"] = cdist(a, b).ravel()
  df.rename(columns = {'atom_num_Ai':'atom_Ai',
                       'atom_num_Bj':'atom_Bj',
                       'atom_type_Ai':'type_Ai',
                       'atom_type_Bj':'type_Bj'}, inplace = True)
  return df

def gro2pdb_simple(gro_df,name,pdb_terminus=['NT3','CT4']):
    
  gro = gro_df.copy()
  gro.loc[:,'atom'] = 'ATOM  '
  gro.loc[:,'occ'] = 1.0
  gro.loc[:,'T_factor'] = 0.0
  gro.loc[:,'x_coord'] = (gro['x_coord'].astype(float)*10)
  gro.loc[:,'y_coord'] = (gro['y_coord'].astype(float)*10)
  gro.loc[:,'z_coord'] = (gro['z_coord'].astype(float)*10)
  gro.loc[:,'empty'] = ''
  #print(gro)
  #gro = gro.replace('OT1','O1')
  #gro = gro.replace('OT2','O2')
  gro = gro.replace('NTR',pdb_terminus[0])
  gro = gro.replace('CTR',pdb_terminus[1])
  gro = pd.DataFrame(gro,columns=['atom','atom_num','empty','atom_name','empty','res_name','chain_id', 'res_num','empty','empty','x_coord','y_coord','z_coord', 'occ','T_factor'])
  #                                         atom     a_num  empty  a_name   alt_i  r_name   chain_id  r_num  code_i empty  x_coord  y_coord  z_coord   occ     T_factor
  np.savetxt(name,gro,delimiter='',fmt=('%6.6s', '%5.5s','%1s', '%-4.4s','%1s', '%-4.4s', '%1.1s','%4.4s','%1s','%3s', '%8.3f', '%8.3f', '%8.3f', '%6.2f', '%6.2f'))

def make_pdb_path(gro_file,saving_path):
  file_basename   = os.path.basename(gro_file)
  file_name, file_extension = os.path.splitext(file_basename)
  new_filepath    = saving_path + file_name + ".pdb"
  return new_filepath

def run_freesasa_custom(file,classifier_path="/home/joaov/Projetos/MMPBSA/test_system/classifier.config",verbose=False):
  c =fs.Classifier(classifier_path)

  structure = fs.Structure(file,c,({
  'hetatm' : False,          # False: skip HETATM
                              # True: include HETATM
  'hydrogen' : True,        # False: ignore hydrogens
                              # True: include hydrogens
  'join-models' : False,     # False: Only use the first MODEL
                              # True: Include all MODELs
  'skip-unknown' : False,    # False: Guess radius for unknown atoms
                              #     based on element
                              # True: Skip unknown atoms
  'halt-at-unknown' : False  # False: set radius for unknown atoms,
                              #    that can not be guessed to 0.
                              # True: Throw exception on unknown atoms.
  }))

  result =fs.calc(structure,fs.Parameters({'algorithm' : fs.ShrakeRupley,
                                              'probe-radius' : 1.4,
                                              'n-points' : 1000}))
  #result =fs.calc(structure,fs.Parameters({'algorithm' : fs.LeeRichards,
  #   
  # 'n-slices' : 100}))

  area_prot =result.totalArea()
  #energy_prot=result.totalArea()*g+b

  structureArray = fs.structureArray(file, {'separate-chains': True,'hydrogen' : True, 'separate-models': False},c)
  if verbose: print(structureArray)
  #en_list=[]
  area_list=[]
  for model in structureArray:
      #print(dir(model))
      result = fs.calc(model,fs.Parameters({'algorithm' : fs.ShrakeRupley,
                                              'probe-radius' : 1.4,
                                              'n-points' : 1000}))
      #energy=result.totalArea()*g+b
      area= result.totalArea()
      #print(model.chainLabel(1) ,area,'En:',energy)
      area_list.append(area)

  #area_monA, area_monB = area_list
  areas=[area_prot,area_list[0],area_list[1]]
  
  return areas

def read_pdb(file):
    pdb_lines = []
    pdb_dict = []
    with open(file) as pdbfile:
        for line in pdbfile:
            if line[:4] == 'ATOM' or line[:6] == "HETATM":
                splitted_line = [line[:6], line[6:11], line[12:16], line[17:20],
                                 line[21], line[22:26], line[30:38], line[38:46], line[46:54]]
                if splitted_line[3] != 'SOL':
                    pdb_lines.append(splitted_line)
                else:
                    continue
                #print (splitted_line)
        for line in pdb_lines:
            line_dict = {
                'atom_num':    line[1].strip(),
                'atom_name':   line[2].strip(),
                'chain_id':    line[4].strip(),
                'x_coord':    float(line[6].strip()),
                'y_coord':    float(line[7].strip()),
                'z_coord':    float(line[8].strip())
            }
            pdb_dict.append(line_dict)
    return pdb_dict

#functions to determine data needed for delphi
def geom_center(df):
    x = (float(max(df['x_coord']))+float(min(df['x_coord'])))/2
    y = (float(max(df['y_coord']))+float(min(df['y_coord'])))/2
    z = (float(max(df['z_coord']))+float(min(df['z_coord'])))/2
    
    #with open('geom_center_xyz.pkl', 'wb') as f:
        #pickle.dump([x,y,z], f)
    return x,y,z

def appropriate_box_size(df):
    x = max(df['x_coord'])-min(df['x_coord'])
    y = max(df['y_coord'])-min(df['y_coord'])
    z = max(df['z_coord'])-min(df['z_coord'])
    max_coor = max(x,y,z)
    return max_coor/0.8

def max_gsize(df):
    x = max(df['x_coord'])
    y = max(df['y_coord'])
    z = max(df['z_coord'])
    max_coor = max(x,y,z)
    if max_coor*2 % 2 == 0:
        max_coor = max_coor*2 +1
    else:
        max_coor = max_coor*2
    return max_coor

def df_snapshot(df, name):
    base_filename = name +'.txt'
    with open(base_filename,'w') as outfile:
        df.to_string(outfile)
    #Neatly allocate all columns and rows to a .txt file


def log_timers(elapsed_time, name):
    base_filename = 'speed_benchmarks.txt'
    line_to_write = '{0:<20} {1:<8}s'.format(name, round(elapsed_time,2))
    with open(base_filename,'a') as outfile:
        outfile.write(line_to_write)
        outfile.write("\n")

def d4p_run_all(monomers,delphi_settings,delphi4py_path="./delphi4github/delphi4py/"):
  
    import sys
    sys.path.insert(0, delphi4py_path)
    from delphi4py import DelPhi4py

    f_crg, f_siz, fpdb, scale, gsize, nlit, nonit, convergence, acent, epsin = delphi_settings
    
    mon1_start, mon1_end, mon2_start, mon2_end = monomers
    #indices
    i_mon1_start = mon1_start - mon1_start #0
    i_mon1_end   = mon1_end   - mon1_start #1148
    i_mon2_start = mon2_start - mon1_start #1149
    i_mon2_end   = mon2_end   - mon1_start #2297

    def d4p_read_dimer(fpdb):
        delphimol = DelPhi4py(
            f_crg,
            f_siz,
            fpdb,
            gsize,  # grid size #nodes in each axis has to be an odd number
            scale,  # scale = 1 / (distance between nodes)
            "single",  # precision
            epsin=epsin,
            conc=0.1,  # ionic strength
            ibctyp=4,  # boundary type: Coulombic
            res2=convergence,  # convergence criterion
            nlit=nlit,  # number of linear iterations
            outputfile="LOG_readFiles",
        )
        return delphimol

    def d4p_run_dimer(delphimol):
        natoms = delphimol.natoms
        p_atpos = delphimol.p_atpos

        delphimol.runDelPhi(
            nonit=nonit,
            nlit=nlit,
            acent=acent,
            relpar=0.2,
            relfac=0.2,
            outputfile="LOG_runDelPhi_dimer",
        )

        return delphimol.getSolvation()


    def d4p_run_monomer(delphimol,i_mon_start,i_mon_end,LOG):
        p_atpos = delphimol.p_atpos     # coordinates
        p_rad3 = delphimol.p_rad3       # raios
        p_chrgv4 = delphimol.p_chrgv4   # cargas
        atinf = delphimol.atinf         # info sobre o atomo p/ warnings

        #print(delphimol)
        for site_atom_position, atom_position in enumerate(range(i_mon_start, i_mon_end+1)):
            #print(i,pos)
            #p_atpos[i][0] = p_atpos[pos][0] # X coord of atom i
            #p_atpos[i][1] = p_atpos[pos][1] # Y coord of atom i
            #p_atpos[i][2] = p_atpos[pos][2] # Z coord of atom i
            
            p_atpos[site_atom_position] = delphimol.p_atpos[atom_position]
            p_rad3[site_atom_position] = delphimol.p_rad3[atom_position]
            p_chrgv4[site_atom_position] = delphimol.p_chrgv4[atom_position]
            atinf[site_atom_position].value = delphimol.atinf[atom_position].value


        natoms = i_mon_end - i_mon_start + 1
        delphimol.changeStructureSize(p_atpos, p_rad3, p_chrgv4, atinf, delphimol.p_iatmed, natoms=natoms)

        delphimol.runDelPhi(
            nonit=nonit,
            nlit=nlit,
            acent=acent,
            relpar=0.2,
            relfac=0.2,
            outputfile=LOG,
        )
        #print("successful exit")
        return delphimol.getSolvation()
#if __name__ == "__main__":


    delphimol = d4p_read_dimer(fpdb)

    original_delphimol = copy(delphimol)

    solvation_dimer = d4p_run_dimer(delphimol)

    solvation_mon1 = d4p_run_monomer(delphimol,i_mon1_start,i_mon1_end,"LOG_runDelPhi_monomer1")

    delphimol = original_delphimol
    solvation_mon2 = d4p_run_monomer(delphimol,i_mon2_start,i_mon2_end,"LOG_runDelPhi_monomer2")

    return solvation_dimer, solvation_mon1, solvation_mon2

def mmpbsa_fullprotocol(gro_file,other_files,monomers,settings):
    
    databases, saving_path, delphi4py_path, sasa_classif_path                           = other_files
    txt_df_res, txt_df_nb, txt_df_charges                                               = databases
    mon1_start, mon1_end, mon2_start, mon2_end                                          = monomers

    cutoff_settings, delphi_settings, termini, verbose                                  = settings
    cutoff, cutoff_n                                                                    = cutoff_settings
    f_crg, f_siz, fpdb, scale, gsize, nlit, nonit, convergence, acent, epsin            = delphi_settings #9

    fpdb=saving_path+(gro_file.split('/')[-1]).split('.')[0]+'.pdb'
    delphi_settings[2]=fpdb

    gro_df = read_gro_df(gro_file) #create df with gro info, creates monomer index

    monA_range = [mon1_start, mon1_end]
    monB_range = [mon2_start, mon2_end]

    gro_df = find_termini(gro_df, termini)
    gro_df = make_index_and_subset(gro_df,monA_range,monB_range)

    gro_df['atom_name'] = gro_df['atom_name'].replace('O2', 'OT2')
    gro_df['atom_name'] = gro_df['atom_name'].replace('O1', 'OT1')

    #Import databases
    df_res, df_charges, df_nb = import_databases(txt_df_res,txt_df_charges,txt_df_nb)

    #add information about atom_type
    gro_res_df = pd.merge(gro_df, df_res,  how = 'left', on = ['res_name', 'atom_name'])
    gro_res_df['atom_type'] = gro_res_df['atom_type'].str.replace('CH2R', 'CH2r')
    gro_res_df = change_res(gro_res_df, 'NTR', termini[0])
    gro_res_df = change_res(gro_res_df, 'CTR', termini[1])

    # add charges do gro df
    gro_res_crg_df = pd.merge(gro_res_df, df_charges,  how = 'left', on = ['res_name', 'atom_name'])

    nans = gro_res_crg_df.query("atom_type != atom_type")
    if not nans.empty: print(nans.head(20)); print(nans.tail(20)); prYellow("There are atoms with no attribution in atom_type!")

    monomers_df = gro_res_crg_df

    monA_df = monomers_df.query('chain_id == "A"')
    trimmed_monA = pd.DataFrame(monA_df,columns=['atom_num','atom_type','charge','x_coord','y_coord','z_coord'])

    monB_df = monomers_df.query('chain_id == "B"')
    trimmed_monB = pd.DataFrame(monB_df,columns=['atom_num','atom_type','charge','x_coord','y_coord','z_coord'])

    # calculate distances between 2 mons
    if verbose: print("Calculating distances...")

    dist_df = calc_cdist(trimmed_monA, trimmed_monB)
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
    gro2pdb_simple(gro_df,fpdb,termini)
    if verbose: prYellow("New file is: "+fpdb)

    import freesasa as fs
    fs.setVerbosity(1)
    
    area_prot, area_monA, area_monB=run_freesasa_custom(fpdb, verbose=True) # (fpdb,classifier_path)
    
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
    
    pdb_coord = read_pdb(fpdb)
    pdb_df = pd.DataFrame(pdb_coord)

    x,y,z=geom_center(pdb_df)
    acent=[x,y,z]
    delphi_settings[8]=acent

    gsize=max_gsize(pdb_df)
    delphi_settings[4]=gsize
    if verbose: print("Geometric center of dimer in pdb is:",acent)

    box_size=appropriate_box_size(pdb_df)
    if verbose: print("Appropriate box side size of pdb is:",box_size)
  
    # Run DelPhi4Py
    monomers = mon1_start, mon1_end, mon2_start, mon2_end

    #if __name__ == "__main__":
    solvation_dimer, solvation_mon1, solvation_mon2 = d4p_run_all(monomers,delphi_settings,delphi4py_path)

    if verbose: print(solvation_dimer); print(solvation_mon1); print(solvation_mon2)
    
    nonbonded   = [VdW_en_total,Coul_en_total]
    sasa        = [area_prot,area_monA,area_monB]
    polar       = [solvation_dimer, solvation_mon1, solvation_mon2]
    
    #return nonbonded, sasa, polar
    
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

    #return VdW_en_total, Coul_en_total, Gnonpolar, Gpolar
    binding_energy= VdW_en_total + Coul_en_total + Gnonpolar + Gpolar
    with open('/home/joaov/Projetos/MMPBSA/test_system/saved_files/01g_energy_summary.txt','a+') as outfile:
      outfile.write('{:<20} {:<20} {:<20} {:<20} {:<20} {:<20}\n'.format(get_file_name(gro_file), VdW_en_total, Coul_en_total, Gnonpolar, Gpolar, binding_energy)) 





# other functions
def change_res(df, to_replace, replace_to):
    df = df.replace(to_replace,replace_to)
    return df

def get_file_name(file):
  file_basename   = os.path.basename(file)
  return file_basename

#color prints
def prYellow(skk):    print("\033[93m {}\033[00m" .format(skk)) 
def prLightGray(skk): print("\033[97m {}\033[00m" .format(skk))