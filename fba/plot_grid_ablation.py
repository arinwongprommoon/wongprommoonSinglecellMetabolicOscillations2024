#!/usr/bin/env python3
import pickle

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from ablation import vget_ablation_ratio
from matrix import ArrayCollection
from grid import heatmap_ablation_grid

model_options = {
    # "glc" or "pyr"
    "carbon_source": "glc",
    # "carbon_source": "pyr",
}


@np.vectorize
def vget_gr(x):
    return x.ablated_flux[0]


saturation_glc = 8.6869
saturation_pyr = 4.4444

exch_rate_dict = {}

if model_options["carbon_source"] == "glc":
    # build exch_rate_dict
    saturation_amm = 1.4848
    exch_rate_dict["r_1714"] = np.linspace(0, 2 * saturation_glc, 32)
    # plot options
    x_axis = exch_rate_dict["r_1714"]
    saturation_carb = saturation_glc
    grid_xlabel_leader = "Glucose exchange"
elif model_options["carbon_source"] == "pyr":
    # build exch_rate_dict
    saturation_amm = 1.0
    exch_rate_dict["r_2033"] = np.linspace(0, 2 * saturation_pyr, 32)
    # plot options
    x_axis = exch_rate_dict["r_2033"]
    saturation_carb = saturation_pyr
    grid_xlabel_leader = "Pyruvate exchange"

exch_rate_dict["r_1654"] = np.linspace(0, 2 * saturation_amm, 32)
y_axis = exch_rate_dict["r_1654"]

# Set up axes parameters
grid_xlabel = f"{grid_xlabel_leader}\n (% growth rate saturation)"
grid_ylabel = f"Ammonium exchange\n (% growth rate saturation)"

# Load saved data
grid_filename = "ec_grid_" + model_options["carbon_source"] + "_amm"
grid_filepath = grid_filename + ".pkl"
with open(grid_filepath, "rb") as handle:
    ablation_result_array = pickle.load(handle)

# Compute data
ratio_array = vget_ablation_ratio(ablation_result_array)
log2ratio_array = np.log2(ratio_array)
gr = ArrayCollection(vget_gr(ablation_result_array), x_axis, y_axis)

# Masks
# ratio > 1 is equivalent to log2(ratio) > 0, so this contour is valid for the
# log2ratio heatmap.
ratio_mask = ratio_array > 1
# For growth rate gradients, non-zero as an 'epsilon' value to get sensible-
# looking contours.
gr_x_mask = gr.gradient.x > 0.01
gr_y_mask = gr.gradient.y > 0.01


# Draw the log2(ratio) heatmap with ratio and growth rate contours
fig, ax = plt.subplots()
heatmap_ablation_grid(
    ax,
    exch_rate_dict,
    log2ratio_array,
    percent_saturation=True,
    saturation_point=(saturation_carb, saturation_amm),
    saturation_grid=True,
    vmin=-0.52,
    vmax=+0.27,
    showticklabels=True,
    center=0,
    cmap="RdBu_r",
    cbar_label=r"$\log_{2}(Ratio)$",
)
# Ratio contour (ratio > 1)
ax.contour(np.rot90(ratio_mask), origin="lower", colors="k")
# Growth rate contours: carbon- and nitrogen-limiting regions
ax.contour(np.rot90(gr_x_mask), origin="lower", colors="C1", linestyles="dashed")
ax.contour(np.rot90(gr_y_mask), origin="lower", colors="C2", linestyles="dashed")
ax.set_xlabel(grid_xlabel)
ax.set_ylabel(grid_ylabel)
ax.set_title(r"$\log_{2}(Ratio)$")


with PdfPages(grid_filename + ".pdf") as pdf:
    for fig_num in range(1, plt.gcf().number + 1):
        pdf.savefig(fig_num)
