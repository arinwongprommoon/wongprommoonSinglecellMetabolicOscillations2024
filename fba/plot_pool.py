#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np

from ablation import get_ablation_ratio
from yeast8model import Yeast8Model


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


# Construct models of cells of interest
glc_exch_rate = 16.89
wt_ec = Yeast8Model("ecYeastGEM_batch_8-6-0.xml")
wt_ec.model.reactions.get_by_id("r_1714").bounds = (-glc_exch_rate, 0)
wt_ec.model.reactions.get_by_id("r_1714_REV").bounds = (0, glc_exch_rate)

# Define steps
fractions = np.linspace(20, 0, num=51)

# Generate results -- normally should take a couple of minutes
abl_list, glc_exchange_fluxes, glc_exchange_rev_fluxes = pool_vs_ablation(
    wt_ec, fractions
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
ax_ratio.set_ylim((0.65, 1))

ax_gr = ax_ratio.twinx()

ax_gr.set_ylabel("Growth rate", color=color_gr)
ax_gr.plot(fractions, gr_list, "o-", color=color_gr)
ax_gr.tick_params(axis="y", color=color_gr, labelcolor=color_gr)
ax_gr.set_xlim((0, 20))
ax_gr.set_ylim((0, 3.5))

plt.savefig("epool.pdf")
