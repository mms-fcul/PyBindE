
import pandas as pd
import numpy as np
import pickle
import MDAnalysis as mda
import subprocess
import freesasa as fs
from copy import copy
import datetime
from scipy.spatial.distance import cdist

# This file's purpose is to serve as a custom python module to
# save useful functions related to the python implementation 
# of MMPBSA written by joaov.

#read convetional gro files, skipping the first 2 lines
def read_gro(file,mon1_start,mon1_end,mon2_start,mon2_end):
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

def read_gro_minus_index(file):
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

def fix_terminus_improved(df,ntr_list,ctr_list):
    terminus = ['NTR','CTR']
    monA = df.query('chain_id == "A"')
    monB = df.query('chain_id == "B"')
    
    #mask = gro_res_df['atom_type'].isnull()
    #list=gro_res_df.loc[ mask,'res_num'].unique()

    min_res_numA = monA['res_num'].min()  # 1
    max_res_numA = monA['res_num'].max()  # 99
    min_res_numB = monB['res_num'].min()  # 100
    max_res_numB = monB['res_num'].max()  # 198

    ntr_res=[min_res_numA,min_res_numB]
    ctr_res=[max_res_numA,max_res_numB]
    
    mask1 = df["res_num"].isin(ctr_res)
    mask2 = df["atom_name"] == 'C'
    df.loc[mask1 & mask2,'atom_name'] = 'CT'

    mask1 = df["res_num"].isin(ntr_res)
    mask2 = df["atom_name"].isin(ntr_list)
    df.loc[mask1 & mask2,'res_name'] = terminus[0]

    mask1 = df["res_num"].isin(ctr_res)
    mask2 = df["atom_name"].isin(ctr_list)
    df.loc[mask1 & mask2,'res_name'] = terminus[1]


#fix dataframe terminus for monomers' C and N terminus

def terminus_fix(gro_dataframe,df_res):
  minval = str(gro_dataframe['res_num'].min())  # 1
  subset_mingro = gro_dataframe.query(
      "res_num == @minval")  # only entries for res 1
  first_res = subset_mingro['res_name'].iloc[0]  # ILE
  subset_mingro['res_name'] = subset_mingro['res_name'].replace(first_res, 'NTR')  # replace all ILE with NTR
  # Grab only the NTR entries from df_res
  nt_res = df_res.query("res_name == 'NTR'")
  subset_mingro['copy_index'] = subset_mingro.index
  mingro_merge = pd.merge(subset_mingro, nt_res,  how = 'left', on = ['res_name', 'atom_name']).set_index('copy_index')  # merge dataframes and output a new column with atom names
  min_nans = mingro_merge.query("atom_type != atom_type")
  min_nans.loc[min_nans.atom_name != 'atom_type', 'res_name'] = first_res
  mingro_merge.update(min_nans)

  maxval = str(gro_dataframe['res_num'].max())  # 99
  subset_maxgro = gro_dataframe.query("res_num == @maxval")
  last_res = subset_maxgro['res_name'].iloc[-1]  # MET
  subset_maxgro['res_name'] = subset_maxgro['res_name'].replace(last_res, 'CTR')
  subset_maxgro['atom_name'] = subset_maxgro['atom_name'].replace('O2', 'OT2')
  subset_maxgro['atom_name'] = subset_maxgro['atom_name'].replace('O1', 'OT1')
  subset_maxgro['atom_name'] = subset_maxgro['atom_name'].replace('C', 'CT')
  ct_res = df_res.query("res_name == 'CTR'")
  # saves a copy from og index
  subset_maxgro['copy_index'] = subset_maxgro.index
  maxgro_merge = pd.merge(subset_maxgro, ct_res,  how = 'left', on = ['res_name', 'atom_name']).set_index('copy_index')
  max_nans = maxgro_merge.query("atom_type != atom_type")
  max_nans.loc[max_nans.atom_name != 'atom_type', 'res_name'] = last_res
  maxgro_merge.update(max_nans)
  return mingro_merge.append(maxgro_merge)

def make_index(dataframe,monA_range,monB_range):
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
    return df_copy

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



#calculate distances between all atoms of 2 monomers
def calc_dists(mon1, mon2, cutoff, cutoff_n):
  distances_db = []
  for a in mon1.index:
      xa = mon1['x_coord'][a]
      ya = mon1['y_coord'][a]
      za = mon1['z_coord'][a]
      num_a = mon1['atom_num'][a]
      type_a = mon1['atom_type'][a]
      charge_a = mon1['charge'][a]
      a_coord = np.array((xa, ya, za), dtype = float)
      for b in mon2.index:
          xb = mon2['x_coord'][b]
          yb = mon2['y_coord'][b]
          zb = mon2['z_coord'][b]
          num_b = mon2['atom_num'][b]
          type_b = mon2['atom_type'][b]
          charge_b = mon2['charge'][b]
          b_coord = np.array((xb, yb, zb), dtype = float)
          distance = np.linalg.norm(b_coord - a_coord)
          if cutoff:
              if distance <= cutoff_n:
                  dist_dict = {
                      'atom_Ai':      num_a, 
                      'atom_Bj':      num_b, 
                      'type_Ai':      type_a, 
                      'type_Bj':      type_b, 
                      'distance':     distance, 
                      'charge_Ai':    charge_a, 
                      'charge_Bj':    charge_b
                  }
                  distances_db.append(dist_dict)
              else:
                  continue
          else:
              dist_dict = {
                  'atom_Ai':      num_a, 
                  'atom_Bj':      num_b, 
                  'type_Ai':      type_a, 
                  'type_Bj':      type_b, 
                  'distance':     distance, 
                  'charge_Ai':    charge_a, 
                  'charge_Bj':    charge_b
              }
              distances_db.append(dist_dict)
  return distances_db

#color prints
def prYellow(skk):    print("\033[93m {}\033[00m" .format(skk)) 
def prLightGray(skk): print("\033[97m {}\033[00m" .format(skk))


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


def gen_en_summary(filename,g,b,nonbonded,sasa,polar,note,saving_path):

  #g = 0.00542 # gamma -> kcal/(mol‚A2)
  #b = 0.92    # beta  -> kcal/mol
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

  Gpolar = Gpolar_complex - (Gpolar_mon1 + Gpolar_mon2) # kJ/(mol.K)

  Gsolv = Gpolar + Gnonpolar

  G_binding = EMM + Gsolv
  
  now = datetime.datetime.now()
  now = (now.strftime('%d-%m-%Y %H:%M:%S'))
  #G_binding = G_complex - G_monomer1 - G_monomer2

  #G_binding = G_complex - G_monomer1 - G_monomer2
  #            G_x  = EMM_x - T*S + Gsolv_x
  #                                 Gsolv_x = Gpolar_x + Gnonpolar_x
  #                   EMM_x = Ebonded + E_VdW + E_Coul
  text_energies = """Generated on {}
filename            =  {}

E_VdW               =  {}
E_Coul              =  {}

Gnonpolar_complex   =  {}
Gnonpolar_mon1      =  {}
Gnonpolar_mon2      =  {}

Gnonpolar_total     =  {}

Gpolar_complex      =  {}
Gpolar_mon1         =  {}
Gpolar_mon2         =  {}

Gpolar_total        =  {}

Gsolv               =  {}

EMM                 =  {}

G_binding           =  {}

User Note : {}  """.format(now,
             filename,
      round(VdW_en_total,2),
      round(Coul_en_total,2),
      round(Gnonpolar_complex,2),
      round(Gnonpolar_mon1,2),   
      round(Gnonpolar_mon2,2),
      round(Gnonpolar,2),
      round(Gpolar_complex,2) ,  
      round(Gpolar_mon1,2)  ,    
      round(Gpolar_mon2,2)  ,
      round(Gpolar,2),
      round(Gsolv,2),
      round(EMM,2),
      round(G_binding,2),
      note
      )

  energies = saving_path
  with open(energies, "w") as out_file:
      out_file.write(text_energies)
      


def run_g_sasa(gmx20_path,file_path,tpr_path,ndx_path):
  
    subprocess.Popen(["/bin/rm", '-f', "aux*.xvg", "conf.gro"])
    subprocess.Popen(["mkdir", 'temp_local'])
    
    # for gromacs 2020.x, this works:
    #a = subprocess.Popen([gmx20_path, 'editconf', '-f', file_path,'-o', "conf.gro", '-n', ndx_path], stdin=subprocess.PIPE)
    #a.communicate(b'Dimer\n')
    #a.wait()
    p = subprocess.Popen([gmx20_path, 'sasa', '-f', file_path, '-s',tpr_path, '-n', ndx_path, '-ndots','1000', '-o', "temp_local/TMP_aux_Prot.xvg", '-surface','"Dimer"','-output','"Dimer"'], stdin=subprocess.PIPE)
    p.wait()
    p = subprocess.Popen([gmx20_path, 'sasa', '-f', file_path, '-s',tpr_path, '-n', ndx_path, '-ndots','1000', '-o', "temp_local/TMP_aux_MonA.xvg", '-surface','"MonomerA"','-output','"MonomerA"'], stdin=subprocess.PIPE)
    p.wait()
    p = subprocess.Popen([gmx20_path, 'sasa', '-f', file_path, '-s',tpr_path, '-n', ndx_path, '-ndots','1000', '-o', "temp_local/TMP_aux_MonB.xvg", '-surface','"MonomerB"','-output','"MonomerB"'], stdin=subprocess.PIPE)
    p.wait()

    P   = subprocess.check_output(["tail -n 1 temp_local/TMP_aux_Prot.xvg | awk '{print $2*100}'"], shell=True).decode("utf-8")
    MA  = subprocess.check_output(["tail -n 1 temp_local/TMP_aux_MonA.xvg | awk '{print $2*100}'"], shell=True).decode("utf-8")
    MB  = subprocess.check_output(["tail -n 1 temp_local/TMP_aux_MonB.xvg | awk '{print $2*100}'"], shell=True).decode("utf-8")

    c=subprocess.run(["/bin/rm -f aux*.xvg conf.gro #conf*#"], shell=True)
    #energy= (P*g+b)
    components=[P,MA,MB]
    return components

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

def gro2pdb_4pymol(gro_df,name):
    
    gro_df['atom'] = 'ATOM  '
    gro_df['occ'] = 1.0
    gro_df['T_factor'] = 0.0
    gro_df['x_coord'] = (gro_df['x_coord'].astype(float)*10)
    gro_df['y_coord'] = (gro_df['y_coord'].astype(float)*10)
    gro_df['z_coord'] = (gro_df['z_coord'].astype(float)*10)
    gro_df['empty'] = ''
    #print(gro_df)

    gro_df = pd.DataFrame(gro_df,columns=['atom','atom_num','empty','atom_name','empty','res_name','chain_id', 'res_num','empty','empty','x_coord','y_coord','z_coord', 'occ','T_factor'])
    #                                         atom     a_num  empty  a_name   alt_i  r_name   chain_id  r_num  code_i empty  x_coord  y_coord  z_coord   occ     T_factor
    np.savetxt(name,gro_df,delimiter='',fmt=('%6.6s', '%5.5s','%1s', '%-4.4s','%1s', '%-4.4s', '%1.1s','%4.4s','%1s','%3s', '%8.3f', '%8.3f', '%8.3f', '%6.2f', '%6.2f'))
    PDB = name
    fin = open(PDB, "rt")
    
    data = fin.read()
    #replace all occurrences of the required string
    to_replace = ['HE11','HE12','HE21','HE22','HH11','HH12','HH21','HH22','HD11','HD12','HD21','HD22']
    replace_to = ['1HE1','2HE1','1HE2','2HE2','1HH1','2HH1','1HH2','2HH2','1HD1','2HD1','1HD2','2HD2']

    for i, pos in enumerate(to_replace):
        data = data.replace(pos,replace_to[i])
    #close the input file
    file_head= ['TITLE     PDB Created by gro2pdb','MODEL        1']
    file_tail= ['ENDMDL','END']

    fin.close() #close the input file

    fin = open(PDB, "wt") #open the input file in write mode
    
    for line in file_head:
        fin.write(line+ '\n')
        
    fin.write(data) #overrite the input file with the resulting data
    
    for line in file_tail:
        fin.write(line+ '\n')
        
    fin.close() #close the file
    
def gro2pdb_simple(gro_df,name,pdb_terminus):
    
    gro = gro_df.copy()
    gro.loc[:,'atom'] = 'ATOM  '
    gro.loc[:,'occ'] = 1.0
    gro.loc[:,'T_factor'] = 0.0
    gro.loc[:,'x_coord'] = (gro['x_coord'].astype(float)*10)
    gro.loc[:,'y_coord'] = (gro['y_coord'].astype(float)*10)
    gro.loc[:,'z_coord'] = (gro['z_coord'].astype(float)*10)
    gro.loc[:,'empty'] = ''
    #print(gro)

    gro = gro.replace('NTR',pdb_terminus[0])
    gro = gro.replace('CTR',pdb_terminus[1])
    gro = pd.DataFrame(gro,columns=['atom','atom_num','empty','atom_name','empty','res_name','chain_id', 'res_num','empty','empty','x_coord','y_coord','z_coord', 'occ','T_factor'])
    #                                         atom     a_num  empty  a_name   alt_i  r_name   chain_id  r_num  code_i empty  x_coord  y_coord  z_coord   occ     T_factor
    np.savetxt(name,gro,delimiter='',fmt=('%6.6s', '%5.5s','%1s', '%-4.4s','%1s', '%-4.4s', '%1.1s','%4.4s','%1s','%3s', '%8.3f', '%8.3f', '%8.3f', '%6.2f', '%6.2f'))


#calculate distances between all atoms of 2 monomers
def calc_dists2(mon1, mon2, cutoff, cutoff_n):
  distances_db = []
  for a in mon1.index:
      xa = float(mon1['x_coord'][a])
      ya = float(mon1['y_coord'][a])
      za = float(mon1['z_coord'][a])
      num_a = mon1['atom_num'][a]
      type_a = mon1['atom_type'][a]
      charge_a = mon1['charge'][a]
      #a_coord = np.array((xa, ya, za), dtype = float)
      for b in mon2.index:
          xb = float(mon2['x_coord'][b])
          yb = float(mon2['y_coord'][b])
          zb = float(mon2['z_coord'][b])
          num_b = mon2['atom_num'][b]
          type_b = mon2['atom_type'][b]
          charge_b = mon2['charge'][b]
          #b_coord = np.array((xb, yb, zb), dtype = float)
          distance2 = (xb-xa)**2+(yb-ya)**2+(zb-za)**2
          #distance=np.linalg.norm(b_coord - a_coord)
          if cutoff:
              if distance <= cutoff_n:
                  dist_dict = {
                      'atom_Ai':      num_a, 
                      'atom_Bj':      num_b, 
                      'type_Ai':      type_a, 
                      'type_Bj':      type_b, 
                      'distance':     distance2, 
                      'charge_Ai':    charge_a, 
                      'charge_Bj':    charge_b
                  }
                  distances_db.append(dist_dict)
              else:
                  continue
          else:
              dist_dict = {
                  'atom_Ai':      num_a, 
                  'atom_Bj':      num_b, 
                  'type_Ai':      type_a, 
                  'type_Bj':      type_b, 
                  'distance':     distance2, 
                  'charge_Ai':    charge_a, 
                  'charge_Bj':    charge_b
              }
              distances_db.append(dist_dict)
  df=pd.DataFrame(distances_db)
  df['distance']=np.sqrt((df['distance']))
  return df

#calculate distances between all atoms of 2 monomers
def calc_dists_pandas(mon1, mon2):
  
  mon1.rename(columns = {'atom_num':'atom_Ai',
                         'atom_type':'type_Ai',
                         'charge':'charge_Ai',
                         'x_coord': 'xa',
                         'y_coord': 'ya',
                         'z_coord': 'za',}, inplace = True)
  
  mon2.rename(columns = {'atom_num':'atom_Bj',
                         'atom_type':'type_Bj',
                         'charge':'charge_Bj',
                         'x_coord': 'xb',
                         'y_coord': 'yb',
                         'z_coord': 'zb',}, inplace = True)
  mon2=mon2.reset_index(drop=True)

  df = pd.concat([mon1, mon2], axis=1)

  df['xa'] = df['xa'].astype(float)
  df['ya'] = df['ya'].astype(float)
  df['za'] = df['za'].astype(float)
  df['xb'] = df['xb'].astype(float)
  df['yb'] = df['yb'].astype(float)
  df['zb'] = df['zb'].astype(float)
    
  df['distance']=((df['xb'])-(df['xa']))**2+((df['yb'])-(df['ya']))**2+((df['zb'])-(df['za']))**2
  df['distance']=np.sqrt((df['distance']))
  
  dist_df = pd.DataFrame(df,columns=['atom_Ai','atom_Bj','type_Ai','type_Bj','charge_Ai','charge_Bj','distance'])
  return dist_df


def calc_cdist_alt(set1, set2):
  #setup
  a = set1[["xa", "ya", "za"]]
  b = set2[["xb", "yb", "zb"]]
  #o suminho
  distances = cdist(a, b)
  #assembly
  index = pd.MultiIndex.from_product([set1["atom_Ai"], set2["atom_Bj"]], names=["atom_Ai", "atom_Bj"])
  df = pd.DataFrame(distances.ravel(), index=index, columns=["distance"]).reset_index()
  #treatment
  df=pd.merge(df,set1,how = 'left', on = ['atom_Ai'])
  df=pd.merge(df,set2,how = 'left', on = ['atom_Bj'])
  df=df.drop(["xa", "ya", "za","xb", "yb", "zb"],axis=1)
  
  return df

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

def mmpbsa_protocol(files,monomers,settings):
    
    gro_file, databases, saving_path, new_filepath     = files
    txt_df_res,txt_df_nb,txt_df_charges                = databases
    mon1_start,mon1_end,mon2_start,mon2_end            = monomers
    
    cutoff_settings,use_gromacs,pdb_terminus,verbose   = settings
    cutoff,cutoff_n                                    = cutoff_settings
    
    
    gro_df = read_gro_minus_index(gro_file) #create df with gro info, creates monomer index

    monA_range = [mon1_start, mon1_end]
    monB_range = [mon2_start, mon2_end]

    gro_df = make_index_and_subset(gro_df, monA_range, monB_range)

    gro_df['atom_name'] = gro_df['atom_name'].replace('O2', 'OT2')
    gro_df['atom_name'] = gro_df['atom_name'].replace('O1', 'OT1')

    #Import databases
    df_res      = pd.read_table(txt_df_res,     delim_whitespace=True,header=None,names=['res_name','atom_name','atom_type']) # df w/ atom names/types from rtp
    df_charges  = pd.read_table(txt_df_charges, delim_whitespace=True,header=None,names=['atom_name','res_name','charge'])
    df_nb       = pd.read_table(txt_df_nb,      delim_whitespace=True,header=None,names=['i','j','c6','c12']) # df from nb.itp to grab c6 and c12

    #print(df_nb[df_nb.duplicated()])
    df_nb = df_nb.rename(columns = {"i": "type_Ai", "j": "type_Bj"})
    #df_nb = df_nb.drop_duplicates()

    #Change C and N terminus residue to CTR and NTR
    has_terminus=True
    if has_terminus:
    
        ntr_list=[]
        for i in df_res.loc[(df_res["res_name"] == 'NTR'), "atom_name"]:
            ntr_list.append(i)

        ctr_list=[]
        for i in df_res.loc[(df_res["res_name"] == 'CTR'), "atom_name"]:
            ctr_list.append(i)

    #terminus = ['NTR','CTR']
    fix_terminus_improved(gro_df,ntr_list,ctr_list)

    #add information about atom_type
    gro_res_df = pd.merge(gro_df, df_res,  how = 'left', on = ['res_name', 'atom_name'])
    gro_res_df['atom_type'] = gro_res_df['atom_type'].str.replace('CH2R', 'CH2r')

    # add charges do gro df
    gro_res_crg_df = pd.merge(gro_res_df, df_charges,  how = 'left', on = ['res_name', 'atom_name'])

    nans = gro_res_crg_df.query("atom_type != atom_type")
    if not nans.empty:
        print(nans)
        prYellow("There are atoms with no attribution in atom_type!")

    monomers_df = gro_res_crg_df

    monA_df = monomers_df.query('chain_id == "A"')
    trimmed_monA = pd.DataFrame(monA_df,columns=['atom_num','atom_type','charge','x_coord','y_coord','z_coord'])

    monB_df = monomers_df.query('chain_id == "B"')
    trimmed_monB = pd.DataFrame(monB_df,columns=['atom_num','atom_type','charge','x_coord','y_coord','z_coord'])

    #print(gro_res_crg_df.head(20))
    #print(gro_res_crg_df.tail(20))
    #df_snapshot(gro_res_crg_df,"gro_res_crg_df")
    # calculate distances between 2 mons
    if verbose: print("Calculating distances...")

    dist_df = calc_cdist(trimmed_monA, trimmed_monB)
    if cutoff: dist_df = dist_df.query("distance <= {}".format(cutoff_n))
    #if verbose: print(dist_df)
    # fix some names between gro and databases
    complete_df = pd.merge(dist_df, df_nb,  how = 'left', on = ['type_Ai', 'type_Bj'])
    complete_df['type_Ai'] = complete_df['type_Ai'].str.replace('CH2R', 'CH2r')
    complete_df['type_Bj'] = complete_df['type_Bj'].str.replace('CH2R', 'CH2r')


    #Vdw and coulomb calculations
    complete_df['En_VdW'] = (complete_df['c12']/complete_df['distance']**12)-(complete_df['c6']/complete_df['distance']**6)

    f = 138.935458
    Er = 1
    complete_df['En_Coul'] = f*(complete_df['charge_Ai']*complete_df['charge_Bj'])/(Er*complete_df['distance'])

    #print(complete_df[complete_df.duplicated()])
    #complete_df = complete_df.drop_duplicates(subset = ['atom_Ai', 'atom_Bj', 'distance'])
    VdW_en_total = complete_df['En_VdW'].sum()
    Coul_en_total = complete_df['En_Coul'].sum()
    if verbose: 
        print(complete_df)
        print("En_VdW = ", VdW_en_total)
        print("En_Coul = ", Coul_en_total)
    #df_snapshot(complete_df,"complete_df_01c")

    nans = complete_df.query("charge_Bj != charge_Bj")
    if not nans.empty:
        print("empty entries in complete_df:",nans)
    #convert gro to pdb
    #pdb_terminus=['NT3','CT4']
    gro2pdb_simple(gro_df,new_filepath,pdb_terminus)
    if verbose: prYellow("New file is: "+new_filepath)
    PDB = new_filepath
    #calculate SASA estimate

    #use_gromacs = False

    if use_gromacs:
    
        gmx20_path  = "/usr/bin/gmx"
        file_path   = gro_file
        tpr_path    = "/home/joaov/github/mmpbsa/gromacs-dependent/files/sas.tpr"
        ndx_path    = "/home/joaov/github/mmpbsa/gromacs-dependent/files/pb.ndx"
        
        sasa_components = run_g_sasa(gmx20_path,file_path,tpr_path,ndx_path)
        
        sasa_P  = sasa_components[0]
        sasa_MA = sasa_components[1]
        sasa_MB = sasa_components[2]
        print(sasa_components)
    
    else:
        
        prot = saving_path+"Prot.pdb"
        monA = saving_path+"monA.pdb"
        monB = saving_path+"monB.pdb"
        
        pdb_terminus=['NT3','CT4']
        
        gro2pdb_simple(gro_df,prot,pdb_terminus)
        gro2pdb_simple(monA_df,monA,pdb_terminus)
        gro2pdb_simple(monB_df,monB,pdb_terminus)

        fs.setVerbosity(1)

        area_prot, energy_prot = run_freesasa(prot)
        area_monA, energy_monA = run_freesasa(monA)
        area_monB, energy_Prot = run_freesasa(monB)

        if verbose: 
            print("Area prot:",area_prot,"En:", energy_prot*4.184)
            print("Area monA:",area_monA,"En:", energy_monA*4.184)
            print("Area monB:",area_monB,"En:", energy_Prot*4.184)

            print("Energy =",(energy_prot-(energy_monA+energy_Prot))*4.184)
        
        sasa_P  = area_prot
        sasa_MA = area_monA
        sasa_MB = area_monB

    pdb_coord = read_pdb(PDB)
    pdb_df = pd.DataFrame(pdb_coord)

    gcenter=geom_center(pdb_df)
    if verbose: print("Geometric center of dimer in pdb is:",gcenter)

    box_size=appropriate_box_size(pdb_df)
    if verbose: print("Appropriate box side size of pdb is:",box_size)
    sasa      = [sasa_P,sasa_MA,sasa_MB]


    #generate energy summary file energies.txt
    g = 0.00542 # gamma -> kcal/(mol‚A2)
    b = 0.92    # beta  -> kcal/mol

    nonbonded = [VdW_en_total,Coul_en_total]
    #polar     = [-8971.42,-3854.31,-3713.99]
    polar =pd.read_pickle("/home/joaov/github/mmpbsa/testing/saved_files/solvation_vals.pkl")
    #polar     = [solv_dimer, solv_mon1, solv_mon2]
    UserNote= ""
    gen_en_summary(gro_file,g,b,nonbonded,sasa,polar,UserNote,saving_path+"energies_01d.txt")
    if verbose: prYellow("New energy summary has been created!")

def calc_emm(files,monomers):

    gro_file, databases = files
    txt_df_res,txt_df_nb,txt_df_charges = databases
    mon1_start,mon1_end,mon2_start,mon2_end = monomers
    
    gro_df = read_gro_minus_index(gro_file) #create df with gro info, creates monomer index

    monA_range = [mon1_start, mon1_end]
    monB_range = [mon2_start, mon2_end]

    gro_df = make_index_and_subset(gro_df, monA_range, monB_range)

    gro_df['atom_name'] = gro_df['atom_name'].replace('O2', 'OT2')
    gro_df['atom_name'] = gro_df['atom_name'].replace('O1', 'OT1')

    #Import databases
    df_res      = pd.read_table(txt_df_res,     delim_whitespace=True,header=None,names=['res_name','atom_name','atom_type']) # df w/ atom names/types from rtp
    df_charges  = pd.read_table(txt_df_charges, delim_whitespace=True,header=None,names=['atom_name','res_name','charge'])
    df_nb       = pd.read_table(txt_df_nb,      delim_whitespace=True,header=None,names=['i','j','c6','c12']) # df from nb.itp to grab c6 and c12

    #print(df_nb[df_nb.duplicated()])
    df_nb = df_nb.rename(columns = {"i": "type_Ai", "j": "type_Bj"})
    #df_nb = df_nb.drop_duplicates()

    #Change C and N terminus residue to CTR and NTR
    has_terminus=True
    if has_terminus:
    
        ntr_list=[]
        for i in df_res.loc[(df_res["res_name"] == 'NTR'), "atom_name"]:
            ntr_list.append(i)

        ctr_list=[]
        for i in df_res.loc[(df_res["res_name"] == 'CTR'), "atom_name"]:
            ctr_list.append(i)

    #terminus = ['NTR','CTR']
    fix_terminus_improved(gro_df,ntr_list,ctr_list)

    #add information about atom_type
    gro_res_df = pd.merge(gro_df, df_res,  how = 'left', on = ['res_name', 'atom_name'])
    gro_res_df['atom_type'] = gro_res_df['atom_type'].str.replace('CH2R', 'CH2r')

    # add charges do gro df
    gro_res_crg_df = pd.merge(gro_res_df, df_charges,  how = 'left', on = ['res_name', 'atom_name'])

    nans = gro_res_crg_df.query("atom_type != atom_type")
    if not nans.empty:
        print(nans)
        prYellow("There are atoms with no attribution in atom_type!")

    monomers_df = gro_res_crg_df

    monA_df = monomers_df.query('chain_id == "A"')
    trimmed_monA = pd.DataFrame(monA_df,columns=['atom_num','atom_type','charge','x_coord','y_coord','z_coord'])

    monB_df = monomers_df.query('chain_id == "B"')
    trimmed_monB = pd.DataFrame(monB_df,columns=['atom_num','atom_type','charge','x_coord','y_coord','z_coord'])


    dist_df = calc_cdist(trimmed_monA, trimmed_monB)

    # fix some names between gro and databases
    complete_df = pd.merge(dist_df, df_nb,  how = 'left', on = ['type_Ai', 'type_Bj'])
    complete_df['type_Ai'] = complete_df['type_Ai'].str.replace('CH2R', 'CH2r')
    complete_df['type_Bj'] = complete_df['type_Bj'].str.replace('CH2R', 'CH2r')


    #Vdw and coulomb calculations
    complete_df['En_VdW'] = (complete_df['c12']/complete_df['distance']**12)-(complete_df['c6']/complete_df['distance']**6)

    f = 138.935458
    Er = 1
    complete_df['En_Coul'] = f*(complete_df['charge_Ai']*complete_df['charge_Bj'])/(Er*complete_df['distance'])

    #print(complete_df[complete_df.duplicated()])
    #complete_df = complete_df.drop_duplicates(subset = ['atom_Ai', 'atom_Bj', 'distance'])

    VdW_en_total = complete_df['En_VdW'].sum()
    Coul_en_total = complete_df['En_Coul'].sum()
    return VdW_en_total, Coul_en_total
