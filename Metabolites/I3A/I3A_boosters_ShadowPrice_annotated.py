
import cobra
import os 
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis
os.chdir('/Users/grichmond/Desktop/Code Catch')
from run_fba import fba
os.chdir('/Users/grichmond/Desktop/reduced_shadow/Just 0/ATP')
ATPmodel = cobra.io.read_sbml_model('iGR_ATPmain.xml')
os.chdir('/Users/grichmond/Desktop')
from set_dms import set_dm
I3Amodel = cobra.io.read_sbml_model('iGR632_v37.xml')


########### Identification of potential boosters based on FVA ###########
#   This section of the code uses the function phase_minMax_sim  to constrain the biomass values within a range of preset biomass. Here, this is based on regions where the slope of the production envelope is changing
#   phase_minMax_sim works by running FVA on a range-locked biomass objective and generating a reduced model where only reactions carrying flux and their metabolites are retained. 
        # Then, if lactate dehydrogenase reaction is still retained in the model, this reaction is set to 0 to prevent any lactate related fermentation 
        # Shadow prices and reduced costs for mets and rxns respectively are calculated from FBA solution vector of the reduced model 
#   Then, the code saves the results of the reduced costs and shadow prices as an excel for continued manual curation to pick interesting boosters to explore 
def phase_minMax_sim(model, min_BM_val, max_BM_val, biomass_rxn, objective_rxn, lac_check=0):
    results={}
    model.reactions.get_by_id(biomass_rxn).bounds= min_BM_val, max_BM_val
    x=flux_variability_analysis(model)
    reduced_model=model.copy()
    to_remove=[]
    for rxn, results in x.iterrows():
        if results[0]<=1e-6 and results[1]<=1e-6:  
            to_remove.append(rxn)
    if lac_check == 1:
        if "LDH_L" not in to_remove: 
            print("lactate had non zero FVA flux")
            model.reactions.LDH_L.bounds=0,0
        else:
            print("lactate removed in FVA step") 
    to_remove = [reduced_model.reactions.get_by_id(rxn_id) for rxn_id in to_remove]
    reduced_model.remove_reactions(to_remove, remove_orphans=True) 
    y=fba(reduced_model, objective=objective_rxn)
    return y.shadow_prices, y.reduced_costs

dm25_i3a_bmVals=[0.15]
dm25_i3a_bm_Range={0.00:0.24}
dm57_i3a_bmVals=[0.143, 0.381]
dm57_i3a_bm_Range={0.00:0.237, 0.32:0.42}

I3A_tmodelc25, media25 = set_dm(I3Amodel, '25')
I3A_tmodelc57, media57 = set_dm(I3Amodel, '57')

phase_25_i3a_results={}
for bm_min, bm_max in dm25_i3a_bm_Range.items():
    phase_I3A_25_SP, phase_I3A_25_RC=phase_minMax_sim(I3A_tmodelc25, bm_min, bm_max, 'curated_biomass', 'IAA_I3A', lac_check=0)
    phase_25_i3a_results[str(str(bm_min) + "-" + str(bm_max))]= phase_I3A_25_SP, phase_I3A_25_RC

phase_57_i3a_results={}
for bm_min, bm_max in dm57_i3a_bm_Range.items():
    phase_I3A_57_SP, phase_I3A_57_RC=phase_minMax_sim(I3A_tmodelc57, bm_min, bm_max, 'curated_biomass', 'IAA_I3A', lac_check=0)
    phase_57_i3a_results[str(str(bm_min) + "-" + str(bm_max))]= phase_I3A_57_SP, phase_I3A_57_RC

def save_to_excel(filename, data_dict, kind="pair"):
    """
    Save shadow prices (SP) and reduced costs (RC) from result dicts to one Excel workbook.
    
    Parameters
    ----------
    filename : str
        Output Excel file name (no extension needed).
    data_dict : dict
        Dictionary where key -> BM_val or range, and value -> tuple of (SP_df, RC_df)
    kind : str
        'open', 'pair', or 'phase'
        'open' means flat outputs with 25_SP etc.
        'pair' means each BM_val or range gets two sheets (SP and RC)
    """
    with pd.ExcelWriter(f"{filename}.xlsx", engine="openpyxl") as writer:
        if kind == "open":
            # expects keys like "25_SP", "25_RC"
            for sheet_name, df in data_dict.items():
                df.to_excel(writer, sheet_name=str(sheet_name)[:31])
        else:
            # expects dict: {BM_value: (SP_df, RC_df)}
            for key, (SP, RC) in data_dict.items():
                sheet_base = str(key).replace(":", "-").replace(".", "_")
                SP.to_excel(writer, sheet_name=f"{sheet_base}_SP"[:31])
                RC.to_excel(writer, sheet_name=f"{sheet_base}_RC"[:31])
    print(f"✅ Saved: {filename}.xlsx")

save_to_excel("I3A_phase_25_results", phase_25_i3a_results, kind="pair")
save_to_excel("I3A_phase_57_results", phase_57_i3a_results, kind="pair")


############# Explore and Filter list of potential boosters #############
#   This code uses the standard approach to generating a production envelope using the PE utils. 
#   However, different boosters are included via the addition of an intracellular exchange constrained to literature-supported supplementation concentrations relative to glucose 
#   Then, a production envelope is generated
import cobra
from cobra.flux_analysis import flux_variability_analysis
import os 
import numpy as np

os.chdir('/Users/grichmond/Desktop/Code Catch')
from run_fba import fba
os.chdir('/Users/grichmond/Desktop')
from set_dms import set_dm
tmodelc = cobra.io.read_sbml_model('iGR632_v37.xml')

def pe_Data(model, objective='curated_biomass', target='EX_i3a_e', num_points=10):
    base_fba = fba(model, objective=objective)
    max_bm = float(base_fba.objective_value)
    bm_values = np.linspace(0, max_bm, num_points)
    results={}
    for b in bm_values:
        mcopy = model.copy()
        mcopy.reactions.get_by_id(objective).bounds = b,b
        sol= flux_variability_analysis(mcopy, processes=1)
        target_min=sol.loc[target, 'minimum']
        target_max=sol.loc[target, 'maximum']
        results[b]= target_min, target_max 
    return results 

testing_media = ["57", "25"]

boosters = {
    "rib_D_c": {"met": tmodelc.metabolites.rib_D_c, "lb": -2.5},   #25% glucose  
    "quln_c": {"met": tmodelc.metabolites.quln_c, "lb": -0.025},    #100uM
    "indole_c": {"met": tmodelc.metabolites.indole_c, "lb": -0.125}, #0.5mM 
}

overall_Results = {}

for media in testing_media:
    media_results = {}
    media_model, media_condition = set_dm(tmodelc, media)

    # --- 1. Unsupplemented condition ---
    unsupp_Data = pe_Data(media_model)
    media_results[f"{media}_unsupplemented"] = unsupp_Data

    # --- 2. Booster-supplemented conditions ---
    for name, info in boosters.items():
        booster_obj = info["met"]
        lb = info["lb"]

        boosted_model = media_model.copy()
        boosted_model.add_boundary(booster_obj, type="sink")
        sink_id = f"SK_{booster_obj.id}"

        boosted_model.reactions.get_by_id(sink_id).lower_bound = lb
        booster_data = pe_Data(boosted_model)
        media_results[name] = booster_data

    # --- 3. Ribose-only condition (ribose = -10, glucose = 0) ---
    ribose_only_model = media_model.copy()

    for rxn in ribose_only_model.exchanges:
        if "glc__D" in rxn.id.lower():
            rxn.lower_bound = 0
    ribose_met = ribose_only_model.metabolites.rib_D_c
    ribose_only_model.add_boundary(ribose_met, type="sink")
    ribose_sink = ribose_only_model.reactions.get_by_id(f"SK_{ribose_met.id}")
    ribose_sink.lower_bound = -10

    ribose_only_data = pe_Data(ribose_only_model)
    media_results["ribose_only"] = ribose_only_data

    overall_Results[media] = media_results

import matplotlib.pyplot as plt
import numpy as np
import os

def plot_production_envelopes(overall_Results, save_dir="/Users/grichmond/Desktop/"):
    """
    Plot all production envelopes on one combined figure.
    - Color encodes condition (unsupplemented, rib_D_c, ribose_only, quln_c, indole_c)
    - Line style encodes media (DM57 = solid, DM25 = dotted)
    - Adds slight jitter to biomass points to improve visibility when overlapping.
    """
    os.makedirs(save_dir, exist_ok=True)

    # --- Color + linestyle maps ---
    color_map = {
        "unsupplemented": "black",
        "rib_D_c": "red",
        "ribose_only": "orange",
        "quln_c": "blue",
        "indole_c": "purple"
    }

    linestyle_map = {
        "57": "solid",
        "25": "dotted"
    }

    plt.figure(figsize=(9, 7))
    plt.rcParams.update({
        "font.family": "Helvetica",
        "axes.linewidth": 1.2,
        "axes.edgecolor": "gray",
        "grid.color": "lightgray",
        "grid.linestyle": "--",
        "grid.linewidth": 0.8,
        "legend.frameon": False,
    })

    for media, media_data in overall_Results.items():
        for condition, pe_results in media_data.items():
            bm_vals = np.array(sorted(pe_results.keys()))
            maxs = np.array([pe_results[b][1] for b in bm_vals])

            bm_vals = np.append(bm_vals, bm_vals[-1])
            maxs = np.append(maxs, 0)

            jitter = np.random.uniform(-0.02, 0.02, size=bm_vals.shape) * bm_vals.max()
            bm_jittered = bm_vals + jitter

            color = "gray"
            for key in color_map:
                if key in condition.lower():
                    color = color_map[key]
                    break

            style = linestyle_map.get(media, "solid")

            label = f"DM{media} - {condition}"
            plt.plot(bm_jittered, maxs, 'o-', color=color, linestyle=style, label=label, alpha=0.85)

    plt.title("I3A Production Envelopes (DM57 vs DM25)", fontsize=15, fontweight="bold")
    plt.xlabel("Biomass Flux (curated_biomass)", fontsize=13)
    plt.ylabel("Target Flux (EX_i3a_e)", fontsize=13)
    plt.grid(alpha=0.3)
    plt.legend(title="Condition / Media", fontsize=9, loc="best")
    plt.tight_layout()

    svg_path = os.path.join(save_dir, "Combined_Production_Envelope.svg")
    plt.savefig(svg_path, format="svg", bbox_inches="tight")
    print(f"✅ Saved combined plot: {svg_path}")
plot_production_envelopes(overall_Results)  
