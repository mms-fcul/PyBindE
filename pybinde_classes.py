#!/usr/bin/python3
# version 2021.09.25
import pandas as pd
import os
import numpy as np
import freesasa as fs
from copy import copy
from scipy.spatial.distance import cdist

lib_dir=os.path.dirname(os.path.realpath(__file__))
databases_path = lib_dir + "/databases/"

class pybinde(object):
    """"""
    
    
    def __init__(
        self,
        gro_file,
        atom_ranges,
        epsin = 4,
        database_charges = os.path.dirname(os.path.realpath(__file__)) + "/databases/DataBaseT.crg",
        database_radii = os.path.dirname(os.path.realpath(__file__)) + "/databases/DataBaseT.siz",
        database_restypes = os.path.dirname(os.path.realpath(__file__)) + "/databases/residues_rtp.txt",
        database_nonbonded = os.path.dirname(os.path.realpath(__file__)) + "/databases/nonbonded_itp.txt",
        database_classifier = os.path.dirname(os.path.realpath(__file__)) + "/databases/classifier.config",
        charges_df = None,
        nonbonded_df = None,
        restypes_df = None,
        obj1 = None,
        obj2 = None,
        dist_df = None,
        vdw_energy=0,
        coul_energy=0,
        pdb_path='',
        SA_whole=0,
        SA_obj1=0,
        SA_obj2=0,
        nonpolar_energy=0,
        delphi4py_path = os.path.dirname(os.path.realpath(__file__)) + "/delphi4py/",
        scale = 2.5,
        nlit = 1000,
        nonit = 0,
        convergence = 0.001,
        polar_energy = 0,
        binding_energy= 0
        ):
        self.gro_file = gro_file
        self.atom_ranges = atom_ranges
        self.epsin = int(epsin)
        self.database_charges = database_charges
        self.database_radii = database_radii
        self.database_restypes = database_restypes
        self.database_nonbonded = database_nonbonded
        self.database_classifier = database_classifier
        self.charges_df = charges_df
        self.nonbonded_df = nonbonded_df
        self.restypes_df = restypes_df
        self.obj1 = obj1
        self.obj2 = obj2
        self.dist_df = dist_df
        self.vdw_energy = vdw_energy
        self.coul_energy = coul_energy
        self.pdb_path = str(pdb_path)
        self.SA_whole = SA_whole
        self.SA_obj1 = SA_obj1
        self.SA_obj2 = SA_obj2
        self.nonpolar_energy = nonpolar_energy
        self.delphi4py_path = delphi4py_path
        self.scale = float(scale)
        self.nlit = int(nlit)
        self.nonit = int(nonit)
        self.convergence = float(convergence)
        self.polar_energy = polar_energy
        self.binding_energy = binding_energy


    
    
    def delphi4py_params(
        scale = 2.5,
        nlit = 1000,
        nonit = 0,
        convergence = 0.001
        
    ):
        scale = float(scale)
        nlit = int(nlit)
        nonit = int(nonit)
        res2 = float(res2)

    def read_gro_df(self):
        
        obj1_range,obj2_range = [self.atom_ranges[0],self.atom_ranges[1]],[self.atom_ranges[2],self.atom_ranges[3]]
        gro_lines = []
        gro_dict = []
        with open(self.gro_file) as gro:
            next(gro)
            next(gro)
            for line in gro:
                splitted_line = [line[0:5], line[5:9], line[10:15], line[15:20],
                                line[21:28], line[29:36], line[37:45]]
                gro_lines.append(splitted_line)
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
        # subset and atributte chain_id
        gro_df["atom_num"] = pd.to_numeric(gro_df["atom_num"])
        monA_start = obj1_range[0]
        monA_end   = obj1_range[1]
        mask = gro_df['atom_num'].between(monA_start,monA_end)
        gro_df.loc[mask, 'chain_id'] = 'A'

        monB_start = obj2_range[0]
        monB_end   = obj2_range[1]
        mask = gro_df['atom_num'].between(monB_start,monB_end)
        gro_df.loc[mask, 'chain_id'] = 'B'
        gro_df_subset = gro_df.query("chain_id == chain_id")
        return gro_df_subset
    
    def import_databases(self):
        df_res      = pd.read_table(self.database_restypes,     delim_whitespace=True,header=None,  names=['res_name','atom_name','atom_type']) # df w/ atom names/types from rtp
        df_charges  = pd.read_table(self.database_charges, delim_whitespace=True,header=None,   names=  ['atom_name','res_name','charge'])
        df_nb       = pd.read_table(self.database_nonbonded,      delim_whitespace=True,header=None, names=['i','j', 'c6','c12']) # df from nb.itp to grab c6 and c12

        #print(df_nb[df_nb.duplicated()])
        df_nb = df_nb.rename(columns = {"i": "type_Ai", "j": "type_Bj"})
        return df_res, df_charges, df_nb
    
    def prYellow(self,skk):    print("\033[93m {}\033[00m" .format(skk))
    
    def find_termini(self,df, termini=['NT3','CT4']):
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
        # proline termini
        pro_ntr_tags= ['H1','H2']
        pro_ntr_list= ['N','H1','H2','CA','CD','C','O' ]
        mask1 = df["res_name"] == 'PRO'
        mask2 = df["atom_name"].isin(pro_ntr_tags)
        pro_ntrs=pd.unique(df.loc[mask1 & mask2,'res_num'])

        mask1 = df["res_num"].isin(pro_ntrs)
        mask2 = df["atom_name"].isin(pro_ntr_list)
        df.loc[mask1 & mask2,'res_name'] = 'NTP' # or NT6, NT7
        #print(df.loc[mask1 & mask2,'res_name'])
        # general
        mask1 = df["res_name"].isin(basic_aa_list)
        mask2 = df["atom_name"].isin(ntr_tags)
        ntrs=pd.unique(df.loc[mask1 & mask2,'res_num'])

        ctr_tags=['OXT','O1','O2','OT1','OT2']

        mask2 = df["atom_name"].isin(ctr_tags)
        ctrs=pd.unique(df.loc[mask1 & mask2,'res_num'])


        ntr_list=['N','H1','H2','H3','CA','C','O']
        mask1 = df["res_num"].isin(ntrs)
        mask2 = df["atom_name"].isin(ntr_list)
        df.loc[mask1 & mask2,'res_name'] = termini[0]

        mask1 = df["res_num"].isin(ctrs)
        mask2 = df["atom_name"] == 'C'
        df.loc[mask1 & mask2,'atom_name'] = 'CT'

        ctr_list=['CT','OT1','OT2','O1','O2','HO11','HO21','HO12','HO22']
        mask2 = df["atom_name"].isin(ctr_list)
        df.loc[mask1 & mask2,'res_name'] = termini[1]

        #print(ntrs,ctrs)
        alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        for num,terminus in enumerate(zip(ntrs,ctrs)):
          mask = df['res_num'].between(terminus[0],terminus[1])
          df.loc[mask, 'real_chain_id'] = alphabet[num]
          #print(num,alphabet[num])
        return df
    
    
    def correct_protonation_state(self,res_to_find,base_state, protonated_state,proton_name,new_proton_name,df):
  
        mask1 = df["res_name"] == res_to_find; RES_atom_nrs=list(df.loc[mask1,'atom_num'])
        mask2 = df["atom_num"].isin(RES_atom_nrs)
        df.loc[mask1 & mask2,'res_name'] = base_state
        #print(RES_atom_nrs)
        #print(RES_atom_nrs)

        mask1 = df["res_name"] == base_state
        mask2 = df["atom_name"] == proton_name
        RES_prot_nums=list(df.loc[mask1 & mask2,'res_num'])
        mask2 = df["res_num"].isin(RES_prot_nums)
        RES_H_atom_nrs=list(df.loc[mask1 & mask2,'atom_num'])
        #print(RES_prot_nums)
        #print(RES_H_atom_nrs)

        mask3 = df["atom_num"].isin(RES_H_atom_nrs)
        df.loc[mask1 & mask2 & mask3,'res_name'] = protonated_state

        mask1 = df["res_name"] == protonated_state; RESH_atom_nrs=list(df.loc[mask1,'atom_num'])
        mask2 = df["atom_name"] == proton_name
        mask3 = df["atom_num"].isin(RESH_atom_nrs)
        df.loc[mask1 & mask2 & mask3,'atom_name'] = new_proton_name


    def correct_protonation_state_HIS(self,df):
        mask1 = df["res_name"] == 'HIS'
        mask2 = df["atom_name"].isin(['HD1','HD2','HE1','HE2'])

        HIS_res_nrs = (df.loc[df['res_name'] == 'HIS','res_num'])

        HISH_identifiers=set(['HD1','HD2','HE1','HE2'])
        HISA_identifiers=set(['HD1','HD2','HE1'])
        HISB_identifiers=set(['HD2','HE1','HE2'])

        for i in pd.unique(HIS_res_nrs):
            each_HIS_df=pd.DataFrame({'res_num':i,'atoms':[list(df.loc[(df['res_num'] == i) &   mask2,'atom_name'])]})

            is_HISH=(each_HIS_df.atoms.map(HISH_identifiers.issubset)[0])
            is_HISA=(each_HIS_df.atoms.map(HISA_identifiers.issubset)[0])
            is_HISB=(each_HIS_df.atoms.map(HISB_identifiers.issubset)[0])

            if is_HISH:
                mask3 = df["res_num"] == i
                df.loc[mask3,'res_name'] = 'HI2' #'HISH' #hi2
            elif is_HISA:
                mask3 = df["res_num"] == i
                df.loc[mask3,'res_name'] = 'HI0' #'HISA' #hi0
            elif is_HISB:
                mask3 = df["res_num"] == i
                df.loc[mask3,'res_name'] = 'HI1' #'HISB' #hi1
            else:
                print("Error in HIS protonation identification!")

    def rename_atom(self,df,residue,old_atom_name,new_atom_name):
        mask1 = df["res_name"] == residue
        mask2 = df["atom_name"] == old_atom_name
        df.loc[mask1 & mask2,'atom_name'] = new_atom_name
    
    
    
    
    
    def construct_objects(self):
        print("⌬  Welcome to PyBindE! ⌬")
        print("Calculation started...")
        print()
        print("File:                              ",(self.gro_file))
        gro_df = self.read_gro_df()
        gro_df = self.find_termini(gro_df)
        
        self.correct_protonation_state('ASP','AS4','AS1','HD2','HD21',gro_df)
        self.correct_protonation_state('ASPH','AS4','AS1','HD2','HD21',gro_df)
        self.correct_protonation_state('CYS','CYS','CY0','HG','HG1',gro_df)
        self.correct_protonation_state('CYSH','CYS','CY0','HG','HG1',gro_df)
        self.correct_protonation_state('CYSC','CY3','CY3','HG','HG1',gro_df)
        self.correct_protonation_state('CYS2','CYS','CYS','xx','xx',gro_df)
        self.correct_protonation_state('GLU','GL4','GL1','HE2','HE21',gro_df)
        self.correct_protonation_state('LYS','LY2','LY3','HZ3','HZ3',gro_df)
        self.correct_protonation_state('LYSH','LY2','LY3','HZ3','HZ3',gro_df)
        self.correct_protonation_state('TYR','TY2','TY0','HH','HH1',gro_df)
        self.correct_protonation_state('CTR','CT4','CT1','HO','HO21',gro_df)
        self.correct_protonation_state('NTR','NT0','NT3','H3','H3',gro_df)
        self.rename_atom(gro_df,'NT5','H','H1')
        self.correct_protonation_state('NTP','NT5','NT7','H2','H2',gro_df)
        #self.correct_protonation_state('SER','SE3','SE0','HG','HG1',gro_df)
        #self.correct_protonation_state('THR','TH3','TH0','HG1','HG1',gro_df)
        #self.correct_protonation_state('HIS','HIS','AS2','HD2','HD21',gro_df)
        self.correct_protonation_state_HIS(gro_df)
        
        self.gro_df = gro_df
        
        df_res, df_charges, df_nb = self.import_databases()
        self.restypes_df, self.charges_df, self.nonbonded_df  = df_res, df_charges, df_nb
        #add information about atom_type
        gro_res_df = pd.merge(gro_df, df_res,  how = 'left', on = ['res_name', 'atom_name'])
        gro_res_df['atom_type'] = gro_res_df['atom_type'].str.replace('CH2R', 'CH2r')
        
        # add charges do gro df
        gro_res_crg_df = pd.merge(gro_res_df, df_charges,  how = 'left', on = ['res_name', 'atom_name'])
        #print(gro_res_crg_df)
        nans = gro_res_crg_df.query("atom_type != atom_type")
        nans_nr=gro_res_crg_df['atom_type'].isna().sum()
        if not nans.empty: print(nans.head(50)); self.prYellow("There are {} atoms with no attribution in atom_type!".format(nans_nr))
        
        nans = gro_res_crg_df.query("charge != charge")
        nans_nr=gro_res_crg_df['charge'].isna().sum()
        if not nans.empty: print(nans.head(50)); self.prYellow("There are {} atoms with no attribution in charge!".format(nans_nr))
        
       
        
        monA_df = gro_res_crg_df.query('chain_id == "A"')
        trimmed_monA = pd.DataFrame(monA_df,columns=['atom_num','atom_type','charge','x_coord','y_coord','z_coord'])

        monB_df = gro_res_crg_df.query('chain_id == "B"')
        trimmed_monB = pd.DataFrame(monB_df,columns=['atom_num','atom_type','charge','x_coord','y_coord','z_coord'])
        self.obj1, self.obj2 = trimmed_monA, trimmed_monB
        return trimmed_monA, trimmed_monB
    
    def calculate_pair_distances(self, cutoff_n=None):
        #print("calculating distances...")
        #obj1, obj2 = self.trimmed_monA, self.trimmed_monB
        a = self.obj1[["x_coord", "y_coord", "z_coord"]]
        b = self.obj2[["x_coord", "y_coord", "z_coord"]]

        dist_df = self.obj1[["atom_num", "atom_type","charge"]].merge(self.obj2[["atom_num", "atom_type","charge"]],       how="cross", suffixes=("_Ai", "_Bj"))
        dist_df["distance"] = cdist(a.astype(float), b.astype(float)).ravel()
        dist_df.rename(columns = {'atom_num_Ai':'atom_Ai',
                             'atom_num_Bj':'atom_Bj',
                             'atom_type_Ai':'type_Ai',
                             'atom_type_Bj':'type_Bj'}, inplace = True)
        
        if cutoff_n: dist_df = dist_df.query("distance <= {}".format(cutoff_n))
        self.dist_df = dist_df
        return dist_df
    
    def calculate_MM_energy(self):
        
        
        complete_df = pd.merge(self.dist_df, self.nonbonded_df,  how = 'left', on = ['type_Ai', 'type_Bj'])
        #print(complete_df)
        nans = complete_df.query("charge_Ai != charge_Ai")
        nans_nr=complete_df['charge_Ai'].isna().sum()
        if not nans.empty: print(nans.head(20)); print(nans.tail(20)); print(pd.unique(nans ['atom_Ai'])); self.prYellow("There are {} atom pairs with missing charge_Ai!".format(nans_nr))
        
        # Vdw and coulomb calculations
        complete_df['En_VdW'] = (complete_df['c12']/complete_df['distance']**12)-(complete_df['c6'] /complete_df['distance']**6)

        f = 138.935458
        complete_df['En_Coul'] = f*(complete_df['charge_Ai'].astype(float)*complete_df['charge_Bj'] .astype(float))/(self.epsin*complete_df['distance'].astype(float))

        self.vdw_energy  = complete_df['En_VdW'].sum()
        self.coul_energy = complete_df['En_Coul'].sum()
        #print(complete_df); print("En_VdW = ", VdW_en_total); print("En_Coul = ",   Coul_en_total)

        print("VdW Energy (kJ/mol):               ", self.vdw_energy )
        print("Coul Energy (kJ/mol):              ", self.coul_energy)
    
    
    
    
    def create_pdb(self,pdb_saving_dir):
        if not pdb_saving_dir.endswith("/"):
            pdb_saving_dir = pdb_saving_dir + "/"
        self.pdb_path = pdb_saving_dir + os.path.splitext(self.gro_file)[0] + ".pdb"
        #print(self.pdb_path)
        gro = self.gro_df.copy()
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
        #gro = gro.replace('NTR',pdb_terminus[0])
        #gro = gro.replace('CTR',pdb_terminus[1])
        gro = pd.DataFrame(gro,columns=['atom','atom_num','empty','atom_name','empty','res_name',   'chain_id', 'res_num','empty','empty','x_coord','y_coord','z_coord', 'occ','T_factor'])
        #                                         atom     a_num  empty  a_name   alt_i  r_name         chain_id  r_num  code_i empty  x_coord  y_coord  z_coord   occ     T_factor
        #print("THIS IS THE PDB:")
        #print(gro)
        np.savetxt(self.pdb_path,gro,delimiter='',fmt=('%6.6s', '%5.5s','%1s', '%-4.4s','%1s', '%-4.4s', '%1.1s','%4.4s','%1s','%3s', '%8.3f', '%8.3f', '%8.3f', '%6.2f', '%6.2f'))
        print("PDB file saved to:                 ", self.pdb_path)
    
    def run_freesasa_custom(self,npoints,verbose=False):
        c =fs.Classifier(self.database_classifier)
        #print(classifier_path)
        structure = fs.Structure(self.pdb_path,c,({
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

        #result =fs.calc(structure,fs.Parameters({'algorithm' : fs.ShrakeRupley,
        #                                            'probe-radius' : 1.4,
        #                                            'n-points' : 1000}))
        result =fs.calc(structure,fs.Parameters({'algorithm' : fs.LeeRichards,
                                                 'probe-radius' : 1.4,
                                                 'n-slices' : npoints}))

        area_prot =result.totalArea()
        #energy_prot=result.totalArea()*g+b

        structureArray = fs.structureArray(self.pdb_path, {'separate-chains': True,'hydrogen' : True,    'separate-models': False},c)
        #if verbose: print(structureArray)
        #en_list=[]
        area_list=[]
        for model in structureArray:
            #print(dir(model))
            #result = fs.calc(model,fs.Parameters({'algorithm' : fs.ShrakeRupley,
            #                                        'probe-radius' : 1.4,
            #                                        'n-points' : 1000}))
            result =fs.calc(model,fs.Parameters({'algorithm' : fs.LeeRichards,
                                                 'probe-radius' : 1.4,
                                                 'n-slices' : npoints}))
            #energy=result.totalArea()*g+b
            area= result.totalArea()
            #print(model.chainLabel(1) ,area,'En:',energy)
            area_list.append(area)

        #area_monA, area_monB = area_list
        areas=[area_prot,area_list[0],area_list[1]]

        return areas
    
    
    
    def calculate_SA_energy(self):
        fs.setVerbosity(fs.nowarnings)
        self.SA_whole, self.SA_obj1, self.SA_obj2 = self.run_freesasa_custom(100)
        
        g = 0.00542 # gamma -> kcal/(mol‚A2)
        b = 0.92    # beta  -> kcal/mol

        energy_whole = (self.SA_whole*g+b)
        energy_obj1 = (self.SA_obj1*g+b)
        energy_obj2 = (self.SA_obj2*g+b)
        
        self.nonpolar_energy = (energy_whole -(energy_obj1 + energy_obj2))*4.184 # kJ/mol
        print("Nonpolar Solvation Energy (kJ/mol):", self.nonpolar_energy)
        
        
        
    def read_pdb(self,file):
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
        
        pdb_df = pd.DataFrame(pdb_dict)
        return pdb_df
    
    #functions to determine data needed for delphi
    def geom_center(self,df):
        x = (float(max(df['x_coord']))+float(min(df['x_coord'])))/2
        y = (float(max(df['y_coord']))+float(min(df['y_coord'])))/2
        z = (float(max(df['z_coord']))+float(min(df['z_coord'])))/2

        #with open('geom_center_xyz.pkl', 'wb') as f:
            #pickle.dump([x,y,z], f)
        return x,y,z

    def appropriate_box_size(self,df):
        x = max(df['x_coord'])-min(df['x_coord'])
        y = max(df['y_coord'])-min(df['y_coord'])
        z = max(df['z_coord'])-min(df['z_coord'])
        max_coor = max(x,y,z)
        return max_coor/0.8

    def max_gsize(self, df, scale, perfil=0.8):
        x = max(df['x_coord'])- min(df['x_coord'])
        y = max(df['y_coord'])- min(df['y_coord'])
        z = max(df['z_coord'])- min(df['z_coord'])
        max_coor = max(x,y,z)
        gsize = int(max_coor * scale / perfil)
        if gsize % 2 == 0:
            gsize += 1

        return gsize

    def df_snapshot(self,df, name):
        base_filename = name +'.txt'
        with open(base_filename,'w') as outfile:
            df.to_string(outfile)
        #Neatly allocate all columns and rows to a .txt file
    
    def d4p_run_all(self,scale,nlit,nonit,conv):

        import sys
        sys.path.insert(0, self.delphi4py_path)
        from delphi4py import DelPhi4py


        mon1_start, mon1_end, mon2_start, mon2_end = self.atom_ranges
        #indices
        i_mon1_start = mon1_start - mon1_start #0
        i_mon1_end   = mon1_end   - mon1_start #1148
        i_mon2_start = mon2_start - mon1_start #1149
        i_mon2_end   = mon2_end   - mon1_start #2297

        def d4p_read_dimer(fpdb):
            delphimol = DelPhi4py(
                self.database_charges,
                self.database_radii,
                fpdb,
                self.gsize,  # grid size #nodes in each axis has to be an odd number
                scale,  # scale = 1 / (distance between nodes)
                "single",  # precision
                epsin=self.epsin,
                conc=0.1,  # ionic strength
                ibctyp=4,  # boundary type: Coulombic
                res2=conv,  # convergence criterion
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
                acent=self.acent,
                relpar=0.9,
                relfac=0.9,
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
            delphimol.changeStructureSize(p_atpos, p_rad3, p_chrgv4, atinf, delphimol.p_iatmed,     natoms=natoms)

            delphimol.runDelPhi(
                nonit=nonit,
                nlit=nlit,
                acent=self.acent,
                relpar=0.9,
                relfac=0.9,
                outputfile=LOG,
            )
            #print("successful exit")
            return delphimol.getSolvation()
#if __name__ == "__main__":


        delphimol = d4p_read_dimer(self.pdb_path)

        original_delphimol = copy(delphimol)

        solvation_dimer = d4p_run_dimer(delphimol)

        solvation_mon1 = d4p_run_monomer(delphimol,i_mon1_start,i_mon1_end, "LOG_runDelPhi_monomer1")

        delphimol = original_delphimol
        solvation_mon2 = d4p_run_monomer(delphimol,i_mon2_start,i_mon2_end, "LOG_runDelPhi_monomer2")

        return solvation_dimer, solvation_mon1, solvation_mon2
    
    def calculate_PB_energy(self,scale,nlit,nonit,conv):
        pdb_df = self.read_pdb(self.pdb_path)
        


        x,y,z=self.geom_center(pdb_df)
        self.acent=[x,y,z]

        self.gsize=self.max_gsize(pdb_df, self.scale)

        #print("Geometric center of dimer in pdb is:",self.acent)


        box_size=self.appropriate_box_size(pdb_df)
        #print("Appropriate box side size of pdb is:",box_size)
        #print("gsize = ", self.gsize)

        # Run DelPhi4Py

        solvation_whole, solvation_obj1, solvation_obj2 = self.d4p_run_all(scale,nlit,nonit,conv)

        t=310 #K
        kb=0.008314463
        Gpolar_whole = solvation_whole*kb*t
        Gpolar_obj1  = solvation_obj1*kb*t
        Gpolar_obj2  = solvation_obj2*kb*t
        self.polar_energy = Gpolar_whole - (Gpolar_obj1 + Gpolar_obj2)
        print("Polar Solvation Energy (kJ/mol):   ", self.polar_energy)
    
    def calculate_binding_energy(self):
        self.binding_energy = self.vdw_energy + self.coul_energy + self.nonpolar_energy + self.polar_energy
        print("Binding Energy (kJ/mol):           ", self.binding_energy)
        
    def save_results(self,energy_summary_file_path,elapsed_time):
        file_basename   = os.path.basename(self.gro_file)
        with open(energy_summary_file_path,'a+') as outfile:
            outfile.write('{:<20} {:<20} {:<20} {:<20} {:<20} {:<20} {:<20}\n'.format(file_basename, self.vdw_energy , self.coul_energy , self.nonpolar_energy , self.polar_energy, self.binding_energy, round(elapsed_time,2)))
        print()
        print("Results saved/appended to:         ",energy_summary_file_path)