#!/usr/bin/env python3
#!/usr/bin/env python3
import pickle

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from ablation import vget_ablation_ratio
from types import SimpleNamespace

model_options = {
    # "glc" or "pyr"
    # "carbon_source": "glc",
    "carbon_source": "pyr",
}


@np.vectorize
def vget_gr(x):
    return x.ablated_flux[0]


class ArrayCollection:
    """2D array of values, with the exchange-rate-0 edge masked, and its gradient.

    Attributes
    ----------
    array : 2d numpy.ndarray
        Copy of the input, with the first row and first column (corresponding
        to an exchange rate of 0) replaced with NaNs, because they are prone to
        weird values.
    gradient : SimpleNamespace
        Gradients of `array` along axis 0 (`gradient.x`) and axis 1
        (`gradient.y`).
    """

    def __init__(self, raw_array):
        self.array = np.array(raw_array, dtype=float)
        self.array[0, :] = np.nan
        self.array[:, 0] = np.nan

        grad_x, grad_y = np.gradient(self.array)
        self.gradient = SimpleNamespace(x=grad_x, y=grad_y)


saturation_glc = 8.6869
saturation_pyr = 4.4444

if model_options["carbon_source"] == "glc":
    saturation_amm = 1.4848
    # Exchange rates of glucose (r_1714)
    x_axis = np.linspace(0, 2 * saturation_glc, 32)
    saturation_carb = saturation_glc
    grid_xlabel_leader = "Glucose exchange"
elif model_options["carbon_source"] == "pyr":
    saturation_amm = 1.0
    # Exchange rates of pyruvate (r_2033)
    x_axis = np.linspace(0, 2 * saturation_pyr, 32)
    saturation_carb = saturation_pyr
    grid_xlabel_leader = "Pyruvate exchange"

# Exchange rates of ammonium (r_1654)
y_axis = np.linspace(0, 2 * saturation_amm, 32)

# Set up axes parameters
grid_xlabel = f"{grid_xlabel_leader}\n (% growth rate saturation)"
grid_ylabel = "Ammonium exchange\n (% growth rate saturation)"

# Load saved data
grid_filename = "ec_grid_" + model_options["carbon_source"] + "_amm"
grid_filepath = grid_filename + ".pkl"
with open(grid_filepath, "rb") as handle:
    ablation_result_array = pickle.load(handle)

# Compute data
ratio_array = vget_ablation_ratio(ablation_result_array)
ratio = ArrayCollection(ratio_array)
log2ratio = ArrayCollection(np.log2(ratio_array))
gr = ArrayCollection(vget_gr(ablation_result_array))

# Masks
ratio_array_mask = ratio.array > 1
# For growth rate gradients, non-zero as an 'epsilon' value to get sensible-
# looking contours.
gr_x_mask = gr.gradient.x > 0.01
gr_y_mask = gr.gradient.y > 0.01

# Axis tick labels as percent of saturation. The y axis is drawn 'in reverse'.
xticklabels = x_axis / saturation_carb * 100
yticklabels = y_axis[::-1] / saturation_amm * 100

# Draw the log2(ratio) heatmap with ratio and growth rate contours
fig, ax = plt.subplots()

# Grid lines to show where saturation (100%) is. The y position is found on the
# un-reversed axis to keep the line position consistent with the x axis.
ax.axvline(np.searchsorted(xticklabels, 100, side="left"), color="k", linewidth=1)
ax.axhline(
    np.searchsorted(yticklabels[::-1], 100, side="right"), color="k", linewidth=1
)

# Rounding directly on the x/yticklabels variables because of known
# matplotlib-seaborn bug:
# - https://github.com/mwaskom/seaborn/issues/1005
sns.heatmap(
    data=np.rot90(log2ratio.array),
    xticklabels=np.around(xticklabels, decimals=1),
    yticklabels=np.around(yticklabels, decimals=1),
    vmin=-0.52,
    vmax=+0.27,
    center=0,
    cmap="RdBu_r",
    cbar_kws={"label": r"$\log_{2}(Ratio)$"},
    ax=ax,
)
# Ratio contour (ratio > 1)
ax.contour(np.rot90(ratio_array_mask), origin="lower", colors="k")
# Growth rate contours: carbon- and nitrogen-limiting regions
ax.contour(np.rot90(gr_x_mask), origin="lower", colors="C1", linestyles="dashed")
ax.contour(np.rot90(gr_y_mask), origin="lower", colors="C2", linestyles="dashed")
ax.set_xlabel(grid_xlabel)
ax.set_ylabel(grid_ylabel)
# Hide tick labels
ax.set(xticklabels=[])
ax.set(yticklabels=[])
ax.set_title(r"$\log_{2}(Ratio)$")


with PdfPages(grid_filename + ".pdf") as pdf:
    for fig_num in range(1, plt.gcf().number + 1):
        pdf.savefig(fig_num)
