
import pandas as pd
import numpy as np
import pickle
import MDAnalysis as mda
import subprocess
from copy import copy
import datetime

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
    
    monA = df.query('chain_id == "A"')
    monB = df.query('chain_id == "B"')

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
    df.loc[mask1 & mask2,'res_name'] = 'NTR'

    mask1 = df["res_num"].isin(ctr_res)
    mask2 = df["atom_name"].isin(ctr_list)
    df.loc[mask1 & mask2,'res_name'] = 'CTR'


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
  Gnonpolar = Gnonpolar*4.184

  #Gpolar -----------------------------------
  Gpolar_complex = polar[0]
  Gpolar_mon1    = polar[1]
  Gpolar_mon2    = polar[2]

  Gpolar = Gpolar_complex - (Gpolar_mon1 + Gpolar_mon2)

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
  # --------------------------------------------------------
  # bash command
  subprocess.Popen(["/bin/rm", '-f', "aux*.xvg", "conf.gro"])
  # -----------------------variables ------------------------
  #gmx20 = Path("/programs/gromacs-OpenMPI2/gromacs-AMD/gromacs-2020.2/bin/gmx")
  #gmx20_path = Path("/usr/bin/gmx")
  #g4_sas = Path("/programs/gromacs-OpenMPI2/gromacs-AMD/gromacs-4.0.7_pH_I/bin/g_sas")
  #g4_editconf = Path("/programs/gromacs-OpenMPI2/gromacs-AMD/gromacs-4.0.7_pH_I/bin/editconf")

  #filepath = "/home/joaov/github/mmpbsa/gromacs-dependent/files/"
  tpr = tpr_path #filepath + "sas.tpr"
  ndx = ndx_path #filepath + "pb.ndx"
  # ------------------------gromacs---------------------------
  # for gromacs 2020.x, this works:

  def run_sasa(output, option):
      a = subprocess.Popen([gmx20_path, 'editconf', '-f', file_path,'-o', "conf.gro", '-n', ndx], stdin=subprocess.PIPE)
      a.communicate(b'Dimer\n')
      a.wait()
      p = subprocess.Popen([gmx20_path, 'sasa', '-f', "conf.gro", '-s',tpr, '-n', ndx, '-ndots','50', '-o', output], stdin=subprocess.PIPE)
      if option == "Dimer":
          p.communicate(b'Dimer\n')
          p.wait()
      elif option == "MonA":
          p.communicate(b'MonomerA\n')
          p.wait()
      elif option == "MonB":
          p.communicate(b'MonomerB\n')
          p.wait()

  subprocess.Popen(["mkdir", 'temp_local'])
  temps ="/home/joaov/github/mmpbsa/gromacs-dependent/final_scripts/TMP/"

  run_sasa("temp_local/TMP_aux_Prot.xvg", "Dimer")
  run_sasa("temp_local/TMP_aux_MonA.xvg", "MonA")
  run_sasa("temp_local/TMP_aux_MonB.xvg", "MonB")
  
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

def gro2pdb(gro_df,name):
    
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
    fin.close()
    #open the input file in write mode
    fin = open(PDB, "wt")
    fin.write(data)
    fin.close()