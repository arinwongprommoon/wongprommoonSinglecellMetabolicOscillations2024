#!/usr/bin/env python3
import operator
import pickle

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from ablation import vget_ablation_ratio
from matrix import ArrayCollection
from grid import heatmap_ablation_grid

model_options = {
    # "glc" or "pyr"
    "carbon_source": "pyr",
}


@np.vectorize
def vget_gr(x):
    return x.ablated_flux[0]


@np.vectorize
def vget_carb(x):
    carb = x.ablated_est_time[3]
    return carb


@np.vectorize
def vget_prot(x):
    prot = x.ablated_est_time[2]
    return prot


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
grid_ylabel_leader = "Ammonium exchange"

# Set up axes parameters
grid_xlabel = f"{grid_xlabel_leader}\n (% growth rate saturation)"
grid_ylabel = f"{grid_ylabel_leader}\n (% growth rate saturation)"
# For quiver
X, Y = np.meshgrid(np.linspace(0, 31, 32), np.linspace(0, 31, 32))

# Load saved data
grid_filename = "ec_grid_" + model_options["carbon_source"] + "_amm"
grid_filepath = grid_filename + ".pkl"
with open(grid_filepath, "rb") as handle:
    ablation_result_array = pickle.load(handle)

# Compute data
ratio_array = vget_ablation_ratio(ablation_result_array)
ratio = ArrayCollection(ratio_array, x_axis, y_axis)

log2ratio = ArrayCollection(np.log2(ratio_array), x_axis, y_axis)

gr = ArrayCollection(vget_gr(ablation_result_array), x_axis, y_axis)

# Masks
ratio_array_mask = ratio.array > 1
# For growth rate gradients, non-zero as an 'epsilon' value to get sensible-
# looking contours.
gr_x_mask = gr.gradient.x > 0.01
gr_y_mask = gr.gradient.y > 0.01


def riced_heatmap(
    ax,
    acoll,
    attribute="array",
    cbar_label=" ",
    title=" ",
    vmin=None,
    vmax=None,
    center=None,
    cmap="RdBu_r",
    ratio_contour=True,
    isratio=False,
    gr_contour=False,
    quiver=False,
):
    """Convenience function to draw heatmaps with quivers

    Parameters
    ----------
    ax : matplotlib.pyplot Axes
        axes to draw on
    acoll : ArrayCollection
        array collection object
    attribute : string
        attribute of acoll to access to draw on heatmap. default "array"
    cbar_label : string
        colour bar label
    title : string
        title of plot
    vmin : float
        min value to show on heatmap
    vmax : float
        max value to show on heatmap
    center : float
        centre value for heatmap
    cmap : string
        matplotlib colour palette to use for colours
    ratio_contour : bool
       if true, draw ratio contour.  further options in 'isratio' parameter.
    isratio : bool
       if true, treats the input array as a ratio, and define contour based on
       where values are less than or greater than 1.  if false, draws contour
       based on the regular definition of ratio.
    gr_contour : bool
       if true, draw contours based on carbon- and nitrogen-limiting regions
       with respect to growth rate.
    quiver : bool
        if true, draw quiver based on susceptibility

    """
    data = operator.attrgetter(attribute)(acoll)
    heatmap_ablation_grid(
        ax,
        exch_rate_dict,
        data,
        percent_saturation=True,
        saturation_point=(saturation_carb, saturation_amm),
        saturation_grid=True,
        vmin=vmin,
        vmax=vmax,
        showticklabels=True,
        center=center,
        cmap=cmap,
        cbar_label=cbar_label,
    )
    if ratio_contour:
        if isratio:
            mask = data > 1
            ax.contour(np.rot90(mask), origin="lower", colors="k", linestyles="dotted")
        else:
            ax.contour(
                np.rot90(ratio_array_mask),
                origin="lower",
                colors="k",
            )
    if gr_contour:
        ax.contour(
            np.rot90(gr_x_mask), origin="lower", colors="C1", linestyles="dashed"
        )
        ax.contour(
            np.rot90(gr_y_mask), origin="lower", colors="C2", linestyles="dashed"
        )
    if quiver:
        ax.quiver(
            X,
            Y,
            acoll.sus_sp.y,
            -acoll.sus_sp.x,
            acoll.sus_sp.magnitudes,
            cmap="autumn",
        )
    ax.set_xlabel(grid_xlabel)
    ax.set_ylabel(grid_ylabel)
    ax.set_title(title)


# Plot!
fig_heatmap_log2ratio, ax_heatmap_log2ratio = plt.subplots()
# This can safely re-use the contour computed on ratio because ratio > 1
# is equivalent to log2ratio > 0.
riced_heatmap(
    ax_heatmap_log2ratio,
    acoll=log2ratio,
    cbar_label=r"$\log_{2}(Ratio)$",
    title=r"$\log_{2}(Ratio)$",
    vmin=-0.52,
    vmax=+0.27,
    center=0,
    gr_contour=True,
)

pdf_filename = grid_filename + ".pdf"
with PdfPages(pdf_filename) as pdf:
    for fig in range(1, plt.gcf().number + 1):
        pdf.savefig(fig)
