###### Hypothesis generation of alternative carbon sources for LGG ######
#   This code does the following:
        # creates a simDM model of LGG 
        #Defines the lower bounds of all non essential amino acids as a carbon equvilent to glucose (6 carbons = -10)
        #Simulates growth rate via FBA the simDM with glucose as a control, then iteratively replaces glucose with each amino acid and calculates and prints the fold change vs glucose 

import cobra 
from set_dms import set_dm
from run_fba import fba

model = cobra.io.read_sbml_model('iGR632_v37.xml')
set_model = set_dm(model, 'model_min')
model2 = set_model[0]

amino_acids = {
    'EX_gly_e': -30,
    'EX_ala_L_e': -20,
    'EX_ser_L_e': -20,
    'EX_val_L_e': -12,
    'EX_thr_L_e': -15,
    'EX_leu_L_e': -10,
    'EX_pro_L_e': -12,
    'EX_met_L_e': -12,
    'EX_phe_L_e': -6.66666666666667,
    'EX_trp_L_e': -5.45454545454545,
    'EX_tyr_L_e': -6.66666666666667,
    'EX_asp_L_e': -15, 
    'EX_gln_L_e': -12,
    'EX_lys_L_e': -10,
    'EX_arg_L_e': -10,
    'EX_his_L_e': -10
}

# === Full glucose condition ===
glc_fba = fba(model2)
print("-10 Glc: " + str(glc_fba.objective_value))
model2.reactions.EX_glc_D_e.lower_bound = -0.05
changed_fba = fba(model2)
print("-1 Glc: " + str(changed_fba.objective_value))

# === Amino acid  ===
aa_results = {}
for amino_acid, concentration in amino_acids.items():
    flag=0
    og_concentration = model2.reactions.get_by_id(amino_acid).lower_bound
    model2.reactions.get_by_id(amino_acid).lower_bound = concentration
    aa_fba = fba(model2)
    if aa_fba.objective_value > changed_fba.objective_value:
        flag=1
        print(amino_acid + ": "+ str(aa_fba.objective_value/changed_fba.objective_value)+ " ; improved")
    model2.reactions.get_by_id(amino_acid).lower_bound = og_concentration
