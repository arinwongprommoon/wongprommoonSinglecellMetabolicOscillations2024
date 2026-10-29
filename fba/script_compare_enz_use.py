#!/usr/bin/env python3
import re

import pandas as pd

from yeast8model import Yeast8Model

# If flux is below this number, treat it as zero
tol = 1.1e-10

# Define properties of sequential-faster and parallel-faster models to
# investigate
models = {
    # Sequential faster = high glucose
    "sequential_faster": {
        "glc_exch_rate": 16.89,
        "amm_exch_rate": None,
        # Create keys to store enzyme fluxes later
        "carbohydrate_fluxes": None,
        "protein_fluxes": None,
    },
    # Parallel faster = custom glucose & ammonium exchange
    "parallel_faster": {
        "glc_exch_rate": 0.194 * 8.6869,
        "amm_exch_rate": 0.71 * 1.4848,
        "carbohydrate_fluxes": None,
        "protein_fluxes": None,
    },
}

# Investigate shared enzyme usage between carbohydrate production and protein
# production
for model, model_props in models.items():
    wt = Yeast8Model("ecYeastGEM_batch_8-6-0.xml")
    # Set glucose exchange rate
    wt.model.reactions.get_by_id("r_1714").bounds = (
        -model_props["glc_exch_rate"],
        0,
    )
    wt.model.reactions.get_by_id("r_1714_REV").bounds = (
        0,
        model_props["glc_exch_rate"],
    )
    # Set ammonium exchange rate
    if model_props["amm_exch_rate"] is None:
        pass
    else:
        wt.model.reactions.get_by_id("r_1654").bounds = (
            -model_props["amm_exch_rate"],
            0,
        )
        wt.model.reactions.get_by_id("r_1654_REV").bounds = (
            0,
            model_props["amm_exch_rate"],
        )
    # Optimise and ablate, store fluxes in each round
    wt.solution = wt.optimize()
    wt.ablation_result = wt.ablate()
    ablation_enzyme_fluxes = wt.ablation_enzyme_fluxes

    # Store enzyme usage reactions that are non-zero during pure carbohydrate
    # synthesis
    nonzero_in_carbohydrate = ablation_enzyme_fluxes["carbohydrate"] > tol
    # Store enzyme usage reactions that are non-zero during pure protein synthesis
    nonzero_in_protein = ablation_enzyme_fluxes["protein"] > tol
    # Find enzyme usage reactions that are non-zero in both pure carbohydrate
    # synthesis and pure protein synthesis
    nonzero_in_both = nonzero_in_carbohydrate & nonzero_in_protein
    # Convert enzyme usage reaction IDs to metabolite IDs of corresponding enzymes
    nonzero_in_both_list = nonzero_in_both[nonzero_in_both].index.to_list()
    enzyme_list = [
        re.sub("draw_prot_(\w+)", lambda x: f"prot_{x.group(1)}[c]", rxn)
        for rxn in nonzero_in_both_list
    ]
    # Get flux of common reactions during carbohydrate synthesis
    carbohydrate_usgs = ablation_enzyme_fluxes["carbohydrate"]
    carbohydrate_fluxes = carbohydrate_usgs[
        carbohydrate_usgs.index.isin(nonzero_in_both_list)
    ]
    model_props["carbohydrate_fluxes"] = carbohydrate_fluxes
    # Get flux of common reactions during protein synthesis
    protein_usgs = ablation_enzyme_fluxes["protein"]
    protein_fluxes = protein_usgs[protein_usgs.index.isin(nonzero_in_both_list)]
    model_props["protein_fluxes"] = protein_fluxes

# Assemble fluxes into a table
fluxes_compare = pd.concat(
    [
        models["sequential_faster"]["carbohydrate_fluxes"],
        models["sequential_faster"]["protein_fluxes"],
        models["parallel_faster"]["carbohydrate_fluxes"],
        models["parallel_faster"]["protein_fluxes"],
    ],
    axis=1,
)
fluxes_compare.columns = [
    "carbohydrate_fluxes_highglucose",
    "protein_fluxes_highglucose",
    "carbohydrate_fluxes_parallel",
    "protein_fluxes_parallel",
]
fluxes_compare = fluxes_compare.reset_index()
fluxes_compare["index"] = fluxes_compare["index"].str.replace("draw_prot_", "")
fluxes_compare = fluxes_compare.rename(columns={"index": "UNIPROT ID"})

# Merge with Uniprot lookup table
uniprot = pd.read_csv("idmapping_2024_07_27.tsv", sep="\t")
uniprot_subset = uniprot[["Entry", "Entry Name", "Protein names", "Gene Names"]]
uniprot_subset = uniprot_subset.rename(columns={"Entry": "UNIPROT ID"})

merged = pd.merge(uniprot_subset, fluxes_compare)
merged = merged.sort_values(by="Entry Name")
# Export
merged.to_csv("fluxes.tsv", sep="\t")
