import pickle

filename="frame00100.gro"

with open('/home/joaov/github/mmpbsa/testing/saved_pickles/energy_vars.pkl', 'rb') as f:
    [VdW_en_total, Coul_en_total] = pickle.load(f)

VdW_en_total  = float(VdW_en_total)
Coul_en_total = float(Coul_en_total)

T = 310
S = 0

#EMM --------------------------------------
Ebonded = 0

EMM = Ebonded + VdW_en_total + Coul_en_total


with open('/home/joaov/github/mmpbsa/testing/saved_pickles/solvation_vals.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
    [solvation_dimer, solvation_mon1, solvation_mon2] = pickle.load(f)

Gpolar_complex = solvation_dimer
Gpolar_mon1    = solvation_mon1
Gpolar_mon2    = solvation_mon2


g = 0.00542 # gamma -> kcal/(mol‚A2)
b = 0.92    # beta  -> kcal/mol

# gromacs 2020 SASA vals
with open('/home/joaov/github/mmpbsa/testing/saved_pickles/g_sas_vals.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
    [MA, MB, P] = pickle.load(f)
# freesasa SASA vals
#with open('/home/joaov/github/mmpbsa/testing/saved_pickles/sasa_vals.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
#    [P,MA,MB] = pickle.load(f)

A_complex = float(P)
A_mon1    = float(MA)
A_mon2    = float(MB)

Gnonpolar_complex = g*A_complex + b
Gnonpolar_mon1    = g*A_mon1    + b
Gnonpolar_mon2    = g*A_mon2    + b

Gnonpolar = Gnonpolar_complex - (Gnonpolar_mon1 + Gnonpolar_mon2)
Gnonpolar = Gnonpolar*4.184
#Gsolv_complex = Gpolar_complex + Gnonpolar_complex
#Gsolv_mon1    = Gpolar_mon1    + Gnonpolar_mon1
#Gsolv_mon2    = Gpolar_mon2    + Gnonpolar_mon2

Gpolar = Gpolar_complex - (Gpolar_mon1 + Gpolar_mon2)

Gsolv = Gpolar + Gnonpolar

#G_complex  = EMM_complex - T*S + Gsolv_complex
#G_monomer1 = EMM_mon1    - T*S + Gsolv_mon1
#G_monomer2 = EMM_mon2    - T*S + Gsolv_mon2


G_binding = EMM + Gsolv
#G_binding = G_complex - G_monomer1 - G_monomer2

#G_binding = G_complex - G_monomer1 - G_monomer2
#            G_x  = EMM_x - T*S + Gsolv_x
#                                 Gsolv_x = Gpolar_x + Gnonpolar_x
#                   EMM_x = Ebonded + E_VdW + E_Coul
text_energies = """
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
""".format(filename,
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
    round(G_binding,2)
    )

energies = "energies.txt"
with open(energies, "w") as out_file:
    out_file.write(text_energies)