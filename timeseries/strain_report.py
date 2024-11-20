#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.backends.backend_pdf import PdfPages
from postprocessor.core.processes.butter import butter, butterParameters
from postprocessor.core.processes.standardscaler import (
    standardscaler,
    standardscalerParameters,
)
from postprocessor.routines.heatmap import heatmap
from postprocessor.routines.single_birth_plot import single_birth_plot
from wela.dataloader import dataloader

from process_pipeline import apply_postprocesses
from signalcollection import SignalCollection

# To draw the BY4741 time series in Fig 1A, define:
# data_options["dataset_name"] = "491_2022_11_21_flavin_by4741_by4742_tsa1tsa2morgan_02_00"
# plot_choices["ts/cell_index"] = [
#    "by4741_019;3;1",
#    "by4741_019;4;2",
#    "by4741_019;11;2",
# ]

# To draw the pHluorin time series in Fig 1A, define:
# data_options["dataset_name"] = "18617_2020_02_21_protAgg_downUpShift_2_0_2_pHluorin_Ura7HA_Ura7HR_00"
# data_options["interval_end_hour"] = 5
# plot_choices["ts/cell_index"] = [
#    "pHluorin_001;57;2",
#    "pHluorin_001;65;1",
#    "pHluorin_001;82;1",
# ]
# "sampling_freq" = 1/3 (line 81)

# For Fig S2C, use the above, but then define:
# data_options["interval_end_hour"] = 15
# plot_choices["ts/cell_index"] = ["pHluorin_001;42;12"]

# To draw the heatmap in Fig 1E, define:
# data_options["dataset_name"] = "19972_2021_05_28_flavin_htb2_glucose_limitation_hard_Delft_01_01_htb2mCherry"
# And define "interval_start_hour" and "interval_end_hour" in such ways to
# capture 0-7 h, 7-15 h, 15-22 h.


data_options = {
    # Name of dataset
    "dataset_name": "26643_2022_05_23_flavin_htb2_glucose_20gpL_01_00_htb2mCherry",
    # Start time, in hours
    "interval_start_hour": 0,
    # End time, in hours
    "interval_end_hour": 24,
}

process_options = {
    # Critical frequency for Butterworth filter
    "critical_freqs": (1 / 350),
    # Signal-to-noise ratio cutoff
    "snr_cutoff": 0.01766784,
}

plot_options = {
    # Font size.  I believe matplotlib default is 10.  Good one for presentations is 16.
    "fontsize": 10,
}

plot_choices = {
    # Plot a few time series
    "ts": True,
    "ts/cell_index": [
        "htb2mCherry_001;1;10",
        "htb2mCherry_001;1;13",
        "htb2mCherry_001;1;19",
    ],
    "ts/legend": True,
    # Plot heatmap
    "heatmap": True,
}

# For core post-processes
# Structure: appending string, process object, parameters object, operates on, gives out
process_dict = {
    "/butter": {
        "runner": butter(
            butterParameters.from_dict(
                {
                    "order": 2,
                    "critical_freqs": process_options["critical_freqs"],
                    "filter_type": "highpass",
                    "sampling_freq": 1 / 5,
                }
            )
        ),
        "signame_endswith": "",
        "input_sigtype": "continuous",
        "output_sigtype": "continuous",
    },
    "_scale": {
        "runner": standardscaler(standardscalerParameters.default()),
        "signame_endswith": "/butter",
        "input_sigtype": "continuous",
        "output_sigtype": "continuous",
    },
}

### IMPORT AND SELECT DATA
# assumes that the experiment name begins with the experiment ID and the ID is
# separated from the rest by an underscore
experimentID = data_options["dataset_name"].split("_")[0]
# Load
dl = dataloader(wdir="tsv_data/individual", ls=False)
dl.load(data_options["dataset_name"] + ".tsv")
# Grab flavin and buddings
flavin_df = dl.wide_df(signal="flavin")
buddings_df = dl.wide_df(signal="buddings")
# Slice by interval
interval_start = data_options["interval_start_hour"] * 12
interval_end = data_options["interval_end_hour"] * 12
flavin_sliced = flavin_df.iloc[:, interval_start:interval_end]
buddings_sliced = buddings_df.iloc[:, interval_start:interval_end]
# Package into SignalCollection object
strain = SignalCollection()
strain.signals["flavin"] = flavin_sliced
strain.signals["buddings"] = buddings_sliced
strain.sigtypes["flavin"] = "continuous"
strain.sigtypes["buddings"] = "binary"

# Apply all processes to all signal collections
apply_postprocesses(strain, process_dict)

### PLOTTING

# Set font size
fontsize = plot_options["fontsize"]
plt.rcParams.update({"font.size": fontsize})

# TIME SERIES -- INDIVIDUAL

# Plot a few time series
if plot_choices["ts"]:
    cell_index_list = plot_choices["ts/cell_index"]
    num_cells = len(cell_index_list)
    tick_spacing = 2

    fig_ts, ax_ts = plt.subplots(
        ncols=1, nrows=num_cells, figsize=(8, 2 * num_cells), sharex=True
    )
    timepoints = strain.signals["buddings"].columns
    for ax_index, cell_index in enumerate(cell_index_list):
        flavin_ts = strain.signals["flavin/butter_scale"].loc[cell_index, :]
        birth_mask = strain.signals["buddings"].loc[cell_index, :]
        plot_title = ""
        ylabel = "Normalised fluorescence (AU)"
        single_birth_plot(
            trace_timepoints=timepoints,
            trace_values=flavin_ts,
            trace_name="flavin",
            birth_mask=birth_mask,
            trace_color="#103080",
            plot_title=plot_title,
            ylabel=ylabel,
            ax=ax_ts[ax_index],
        )
        # Set tick spacing
        ax_ts[ax_index].xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
        # Override axis labelling behaviour of single_birth_plot(),
        # in favour of common label specified by plt
        ax_ts[ax_index].set_xlabel("")
        ax_ts[ax_index].set_ylabel("")
        # Remove legend
        if not plot_choices["ts/legend"]:
            ax_ts[ax_index].get_legend().remove()
    fig_ts.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor="none", top=False, bottom=False, left=False, right=False)
    plt.grid(False)
    plt.xlabel("Time (hr)")
    plt.ylabel("Normalised fluorescence (AU)")

# TIME SERIES -- POPULATION

# Plot heatmap
if plot_choices["heatmap"]:
    fig_heatmap, ax_heatmap = plt.subplots()
    heatmap(
        strain.signals["flavin/butter"].dropna(how="all"),
        trace_name="flavin",
        buddings_df=strain.signals["buddings"].loc[
            strain.signals["flavin/butter"]
            .dropna()
            .index.intersection(strain.signals["buddings"].index)
        ],
        plot_title="Flavin fluorescence",
        unit_scaling=1 / 12,
        xtick_step=2,
        xlabel="Time (hours)",
        ax=ax_heatmap,
    )

# Save figures
pdf_filename = f"{experimentID}_plots.pdf"
with PdfPages(pdf_filename) as pdf:
    for fig in range(1, plt.gcf().number + 1):
        pdf.savefig(fig)
# Close all figures
plt.close("all")
