#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from matplotlib.backends.backend_pdf import PdfPages
from yeast8model import Yeast8Model

# "highglc", "glc_highratio", "highpyr", "pyrhighratio" : use that preset model
# "all" : generate PDFs for all preset models
# "custom" : use custom model, defined by model_options variable
preset_model_options = "all"

model_options = {
    "glc_exch_rate": 16.89,
    "pyr_exch_rate": None,
    "amm_exch_rate": None,
}

compute_options = {
    # Top N reactions (of original enzyme usage fluxes) for the topflux plot.
    # If None, it takes all the reactions with non-zero flux.
    # If 0, it takes all the reactions.
    "topflux/ntop": None,
}

preset_models = {
    "highglc": {
        "glc_exch_rate": 16.89,
        "pyr_exch_rate": None,
        "amm_exch_rate": None,
    },
    "glchighratio": {
        "glc_exch_rate": 0.194 * 8.6869,
        "pyr_exch_rate": None,
        "amm_exch_rate": 0.71 * 1.4848,
    },
    "highpyr": {
        "glc_exch_rate": 0,
        "pyr_exch_rate": 8.8888,
        "amm_exch_rate": None,
    },
    "pyrhighratio": {
        "glc_exch_rate": 0,
        "pyr_exch_rate": 0.839 * 4.4444,
        "amm_exch_rate": 0.903 * 1.0,
    },
}


def prettyfloat(x):
    if x is None:
        # Unrestricted
        out_str = "Unres"
    else:
        out_str = f"{x:05.2f}".replace(".", "p")
    return out_str


def get_topn_list(series, ntop):
    """Get top N flux-carrying reactions from a Series."""
    return series.sort_values(ascending=False)[:ntop].index.to_list()


def rxns_to_hues(rxn_list, hue_lookup):
    """Convert reactions to hues"""
    hues = []
    for rxn_id in rxn_list:
        try:
            hue = hue_lookup[rxn_id]
            hues.append(hue)
        except KeyError:
            hues.append(np.nan)
    return hues


def drawplots(model_options):
    if model_options["glc_exch_rate"] is None:
        wt.model.reactions.get_by_id("r_1714").bounds = (-16.89, 0)
        wt.model.reactions.get_by_id("r_1714_REV").bounds = (0, 16.89)
    else:
        wt.model.reactions.get_by_id("r_1714").bounds = (
            -model_options["glc_exch_rate"],
            0,
        )
        wt.model.reactions.get_by_id("r_1714_REV").bounds = (
            0,
            model_options["glc_exch_rate"],
        )

    if model_options["pyr_exch_rate"] is None:
        pass
    else:
        # wt.model.reactions.get_by_id("r_1714").bounds = (0, 0)
        wt.model.reactions.get_by_id("r_2033").bounds = (
            -model_options["pyr_exch_rate"],
            0,
        )
        wt.model.reactions.get_by_id("r_2033_REV").bounds = (
            0,
            model_options["pyr_exch_rate"],
        )

    if model_options["amm_exch_rate"] is None:
        pass
    else:
        wt.model.reactions.get_by_id("r_1654").bounds = (
            -model_options["amm_exch_rate"],
            0,
        )
        wt.model.reactions.get_by_id("r_1654_REV").bounds = (
            0,
            model_options["amm_exch_rate"],
        )

    wt.solution = wt.optimize()
    wt.ablation_result = wt.ablate()
    ablation_enzyme_fluxes = wt.ablation_enzyme_fluxes

    # Convert dictionary of pandas dataframes to numpy array for various inputs
    enz_use_array = np.stack([df.to_numpy() for df in ablation_enzyme_fluxes.values()])
    # Remove enzymes that have all-zeros across components
    # because (a) they're not informative,
    # (b) they cause problems in downstream functions
    enz_use_array = enz_use_array[:, np.any(enz_use_array, axis=0)]

    list_components = list(ablation_enzyme_fluxes.keys())
    list_components = [
        component.replace("original", "parallel") for component in list_components
    ]

    # Draw plot
    # take all reactions with non-zero flux
    if compute_options["topflux/ntop"] is None:
        ntop = np.sum(ablation_enzyme_fluxes["original"] != 0)
        print(f"topflux: number of reactions = {ntop}")
    # take all reactions
    elif compute_options["topflux/ntop"] == 0:
        ntop = len(ablation_enzyme_fluxes["original"])
        print(f"topflux: number of reactions = {ntop}")
    # take top N reactions
    else:
        ntop = np.sum(ablation_enzyme_fluxes["original"] != 0)
        print(f"topflux: number of reactions = {ntop}")

    # List of top N reactions, original (un-ablated)
    original_topn_list = get_topn_list(ablation_enzyme_fluxes["original"], ntop)

    # Assign 'hues' and create lookup table
    hue_lookup = dict((zip(original_topn_list, range(ntop))))

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
    )
    ax.set_xlabel("Biomass component")
    ax.set_ylabel("Rank")

    filename = (
        "CompareEnzUse"
        + "_glc"
        + prettyfloat(model_options["glc_exch_rate"])
        + "_pyr"
        + prettyfloat(model_options["pyr_exch_rate"])
        + "_amm"
        + prettyfloat(model_options["amm_exch_rate"])
    )

    # Save all open figures to PDF
    pdf_filename = filename + ".pdf"
    with PdfPages(pdf_filename) as pdf:
        for fig in range(1, plt.gcf().number + 1):
            pdf.savefig(fig)
    # Close all figures
    plt.close("all")


if preset_model_options == "custom":
    wt = Yeast8Model("ecYeastGEM_batch_8-6-0.xml")
    drawplots(model_options)
elif preset_model_options == "all":
    for _, val in preset_models.items():
        wt = Yeast8Model("ecYeastGEM_batch_8-6-0.xml")
        drawplots(val)
else:
    wt = Yeast8Model("ecYeastGEM_batch_8-6-0.xml")
    drawplots(preset_models[model_options])
