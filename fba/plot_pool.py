#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np

from ablation import get_ablation_ratio
from matplotlib.backends.backend_pdf import PdfPages
from yeast8model import Yeast8Model

# "highglc", "glc_highratio", "highpyr", "pyrhighratio" : use that preset model
# "all" : generate PDFs for all preset models
# "custom" : use custom model, defined by model_options variable
preset_model_options = "glchighratio"

model_options = {
    "glc_exch_rate": 16.89,
    "pyr_exch_rate": None,
    "amm_exch_rate": None,
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


def pool_vs_ablation(ymodel, fractions):
    # get original pool upper bound
    orig_epool = ymodel.model.reactions.get_by_id("prot_pool_exchange").bounds[1]

    abl_res_list = []
    glc_exchange_fluxes = []
    glc_exchange_rev_fluxes = []
    for fraction in fractions:
        # impose pool constraint
        ub = fraction * orig_epool
        ymodel.model.reactions.get_by_id("prot_pool_exchange").bounds = (0, ub)
        sol = ymodel.optimize()
        # ablate & ratio
        abl_res = ymodel.ablate()
        # append lists
        abl_res_list.append(abl_res)
        glc_exchange_fluxes.append(sol.fluxes["r_1714"])
        glc_exchange_rev_fluxes.append(sol.fluxes["r_1714_REV"])

    return abl_res_list, glc_exchange_fluxes, glc_exchange_rev_fluxes


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

    # Define steps
    fractions = np.linspace(2, 0, num=51)

    # Generate results -- normally should take a couple of minutes
    abl_list, glc_exchange_fluxes, glc_exchange_rev_fluxes = pool_vs_ablation(
        wt, fractions
    )

    # Define lists for plotting
    gr_list = [abl_res.ablated_flux.iloc[0] for abl_res in abl_list]
    ratio_list = [get_ablation_ratio(abl_res) for abl_res in abl_list]
    abl_flux_list = [abl_res.ablated_flux.iloc[1:].to_list() for abl_res in abl_list]
    abl_flux_array = np.array(abl_flux_list)
    biomass_components_list = abl_list[0].priority_component.iloc[1:].to_list()

    # Draw plot
    fig, ax_ratio = plt.subplots()
    color_ratio = "k"
    color_gr = "tab:blue"

    ax_ratio.set_xlabel("Relative size of enzyme-available proteome pool")

    ax_ratio.set_ylabel(r"$\tau_{\mathrm{seq}/\mathrm{par}}$", color=color_ratio)
    ax_ratio.plot(fractions[:-1], ratio_list[:-1], "o-", color=color_ratio)
    ax_ratio.tick_params(axis="y", color=color_ratio, labelcolor=color_ratio)
    ax_ratio.set_xlim((0, 2))
    # ax_ratio.set_ylim((0.65, 1))

    ax_gr = ax_ratio.twinx()

    ax_gr.set_ylabel("Growth rate", color=color_gr)
    ax_gr.plot(fractions, gr_list, "o-", color=color_gr)
    ax_gr.tick_params(axis="y", color=color_gr, labelcolor=color_gr)
    ax_gr.set_xlim((0, 2))
    ax_gr.set_ylim((0, 3.5))

    filename = (
        "epool"
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
    drawplots(preset_models[preset_model_options])
