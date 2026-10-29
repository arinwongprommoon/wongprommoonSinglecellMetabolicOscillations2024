#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from yeast8model import Yeast8Model

# This model: high glucose uptake; pyruvate and ammonium unrestricted
glc_exch_rate = 16.89


def get_topn_list(series, ntop):
    """Get top N flux-carrying reactions from a Series."""
    return series.sort_values(ascending=False)[:ntop].index.to_list()


def rxns_to_hues(rxn_list, hue_lookup):
    """Convert reactions to hues"""
    return [hue_lookup.get(rxn_id, np.nan) for rxn_id in rxn_list]


# Construct and solve model
wt = Yeast8Model("ecYeastGEM_batch_8-6-0.xml")
wt.model.reactions.get_by_id("r_1714").bounds = (-glc_exch_rate, 0)
wt.model.reactions.get_by_id("r_1714_REV").bounds = (0, glc_exch_rate)
wt.solution = wt.optimize()
wt.ablation_result = wt.ablate()
ablation_enzyme_fluxes = wt.ablation_enzyme_fluxes

list_components = [
    component.replace("original", "parallel")
    for component in ablation_enzyme_fluxes.keys()
]

# Take all reactions with non-zero flux in the original (un-ablated) model
ntop = np.sum(ablation_enzyme_fluxes["original"] != 0)

# Top N reactions, original (un-ablated); rank in this list defines the hue
original_topn_list = get_topn_list(ablation_enzyme_fluxes["original"], ntop)
hue_lookup = dict(zip(original_topn_list, range(ntop)))

# Find hues for all components
hues_array = []
for series in ablation_enzyme_fluxes.values():
    topn_list = get_topn_list(series, ntop)
    hues = rxns_to_hues(topn_list, hue_lookup)
    hues_array.append(hues)

hues_array = np.array(hues_array).T


# Visualise
fig, ax = plt.subplots(figsize=(5, 8))
sns.heatmap(
    hues_array,
    xticklabels=list_components,
    cmap="plasma_r",
    cbar=False,
    ax=ax,
)
ax.set_xlabel("Biomass component")
ax.set_ylabel("Rank")

fig.savefig("CompareEnzUse_glc16p89_pyrUnres_ammUnres.pdf")
plt.close(fig)
