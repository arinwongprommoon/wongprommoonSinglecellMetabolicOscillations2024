#!/usr/bin/env python3
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from yeast8model import Yeast8Model

# Construct models of strains
glc_exch_rate = 16.89

wt = Yeast8Model("ecYeastGEM_batch_8-6-0.xml")
wt.model.reactions.get_by_id("r_1714").bounds = (-glc_exch_rate, 0)
wt.model.reactions.get_by_id("r_1714_REV").bounds = (0, glc_exch_rate)
wt.solution = wt.optimize()

# Ablate and save
ablation_results = wt.ablate()
component_list = ablation_results.priority_component.to_list()[1:]
ymodel_keys = ablation_results.keys()

# Construct data structures for bar plots
# Hard-coding positions because honestly, it is easier.
est_times = ablation_results.iloc[1:, 2].to_list()
Tpar = ablation_results.iloc[2, 3]

# Draw
# Set plot parameters
plt.rcParams.update(
    {
        "font.size": 12,
        "axes.prop_cycle": plt.cycler(color=plt.cm.Pastel2.colors),
    }
)
# Set colour palette
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

n = len(component_list)

fig, ax = plt.subplots(figsize=(3.0, 4.8))

# Stacked bar: estimated time per biomass component
bottom = 0
for time, color in zip(est_times, colors[:n]):
    ax.bar(-0.2, time, bottom=bottom, width=0.3, color=color)
    bottom += time

# Single bar: parallel time
ax.bar(0.2, Tpar, width=0.3, color=colors[n])

ax.set_xticks([0])
ax.set_xticklabels(["ecYeast8"])
ax.set_xlabel("Model")
ax.set_ylabel("Estimated time (h)")

# Legend: colour = biomass component (+ parallel time)
# assumes that colors (usually 10 elements) is longer than
# component_list (usually 7 elements)
# Using colour patches ref:
# https://matplotlib.org/stable/tutorials/intermediate/legend_guide.html#creating-artists-specifically-for-adding-to-the-legend-aka-proxy-artists
handles = [
    mpatches.Patch(color=c, label=label)
    # Does not render the subscript via LaTeX, but easily fixable in a vector
    # graphics editor
    for c, label in zip(colors[: n + 1], component_list + [r"T_{par}"])
]
# Move legend out of the way
ax.legend(
    handles=handles,
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),  # just outside the right edge of the axes
    frameon=False,
)

# bbox_inches="tight" expands the canvas so the outside legend is not clipped
fig.savefig("combined_barplot.pdf", bbox_inches="tight")
plt.close(fig)
