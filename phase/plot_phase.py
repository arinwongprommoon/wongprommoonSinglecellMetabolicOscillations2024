"""Create phase dataframes and plot results."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from wela.dataloader import dataloader

from all_experiments_plots import (
    augment_data,
    catplot_from_period_dict,
    plot_periods_ratios,
    plot_single_cell_data,
)
from data_dict import data_dict, get_oscillators
from find_phases import get_signals
from single_experiment_plots import (
    compare_periods_for_single_cells,
    find_2dmax,
    plot_flavin_phase,
    plot_periods,
    plot_time_series,
    plot_torus,
    raster_plot,
)

use_oscillators = False
no_plots_for_peaks = 3

rng = np.random.default_rng()


def add_experiment(summary_dict, data_dict, data_id, group):
    """Add details of experiment to dict."""
    summary_dict["omero_no"].append(data_dict[data_id][0].split("_")[0])
    summary_dict["condition"].append(data_dict[data_id][2])
    summary_dict["group"].append(group)


def add_to_summary_dict(summary_dict, data_dict, data_id, group, phases_df):
    """Add results for all experiments to a dict to generate a dataframe ."""
    for i, (xv, yv) in enumerate(
        zip(
            ["flavin_period", "htb2_period", "htb2_period"],
            ["buddings_period", "buddings_period", "flavin_period"],
        )
    ):
        # periods
        twodg = find_2dmax(phases_df, xv, yv)
        if i < 2:
            # flavin and htb2
            for dpt in augment_data(
                twodg.get("x_mean", np.nan), twodg.get("x_std", np.nan)
            ):
                summary_dict["period_variable"].append(xv.split("_")[0])
                summary_dict["period"].append(dpt)
                add_experiment(summary_dict, data_dict, data_id, group)
        if i == 0:
            # buddings
            for dpt in augment_data(
                twodg.get("y_mean", np.nan), twodg.get("y_std", np.nan)
            ):
                summary_dict["period_variable"].append(yv.split("_")[0])
                summary_dict["period"].append(dpt)
        _, ratio = compare_periods_for_single_cells(phases_df, xv, yv)
        # period-to-period
        for dpt in augment_data(np.mean(ratio), np.std(ratio)):
            summary_dict["ratio_type"].append(f"{xv}_to_{yv}")
            summary_dict["ratio"].append(dpt)
            if i == 0:
                add_experiment(summary_dict, data_dict, data_id, group)


def add_to_period_dict(period_dict, period_res):
    """Add results of plot_period to dict."""
    for field in ["flavin", "buddings", "htb2"]:
        for key in [f"{field}_period", f"{field}_period_std"]:
            if key in period_res:
                period_dict[key].append(period_res[key])
            else:
                period_dict[key].append(np.nan)


def make_single_experiment_plots(
    condition, phases_df, flavin_phase_df, period_dict, title
):
    """Plot results for a single experiment."""
    plot_torus(phases_df, title)
    period_res = plot_periods(phases_df, title)
    add_to_period_dict(period_dict, period_res)
    plot_flavin_phase(flavin_phase_df, title)
    if "pyr" in condition:
        tstar = 12
        xlim = [-6.1, 0.1]
    elif "glu K" in condition:
        tstar = 14
        xlim = [-4.1, 0.1]
    else:
        tstar = 8
        xlim = [-4.1, 0.1]
    if "buddings_phase" in phases_df.columns and "htb2_phase" in phases_df.columns:
        rsignals = ["buddings", "htb2"]
    else:
        rsignals = ["buddings"]
    lin_order = raster_plot(
        phases_df,
        signals=rsignals,
        tstar=tstar,
        title=f"{title}",
        xlim=xlim,
    )
    if False and "1649" in data_dict[data_id][0]:
        for cell_id in lin_order[50:56]:
            plot_time_series(
                cell_id,
                phases_df,
                title,
                tlim=[tstar + xlim[0], tstar + 3],
            )


#########################
# analyse all experiments
#########################
summary_dict = {
    "omero_no": [],
    "condition": [],
    "group": [],
    "period": [],
    "period_variable": [],
    "ratio": [],
    "ratio_type": [],
}
period_dict = {"omero_no": [], "condition": [], "group": []}
for field in ["flavin", "buddings", "htb2"]:
    period_dict[f"{field}_period"] = []
    period_dict[f"{field}_period_std"] = []
all_phases_dfs = []

# run through each experiment
for data_id in data_dict:
    if "19972" in data_dict[data_id][0] or "20212" in data_dict[data_id][0]:
        continue
    elif "613" in data_dict[data_id][0]:
        data_dict[data_id][2] = "2% glu K"
    dl = dataloader(wdir="tsv_data/individual", ls=False)
    dl.load(data_dict[data_id][0], use_tsv=True)
    if use_oscillators:
        osc = get_oscillators(data_dict[data_id][1])
        groups = list(osc.keys())
    else:
        groups = list(dl.df.group.unique())
    condition = data_dict[data_id][2]
    omero_no = data_dict[data_id][0].split("_")[0]
    for group in groups:
        if use_oscillators:
            odf = dl.df[dl.df.id.isin(osc[group][0])]
        else:
            odf = dl.df
        odf = odf.rename(columns={"mean_mCherry": "mCherry"})
        if not odf.empty:
            signal_names = ["flavin", "buddings"]
            if "mCherry" in odf.columns:
                signal_names.append("mCherry")
            title = f"{omero_no} - {condition}\n{group}"
            # randomly choose cells to plot
            i_to_plot = rng.integers(
                low=0,
                high=odf.id.unique().shape[0],
                size=no_plots_for_peaks,
            )
            # get results
            if "pyr" in condition:
                res = get_signals(
                    signal_names,
                    dl,
                    group,
                    odf,
                    title,
                    i_to_plot,
                    flavin_prominence=2.5,
                    flavin_lowpass_frac=0.25,
                    htb2_filter=False,
                )
            elif "0.001% glu" in condition:
                # no detectable oscillations
                res = get_signals(
                    signal_names,
                    dl,
                    group,
                    odf,
                    title,
                    i_to_plot,
                    flavin_prominence=0.5,
                    flavin_lowpass_frac=0.25,
                    htb2_filter=True,
                )
            else:
                res = get_signals(
                    signal_names,
                    dl,
                    group,
                    odf,
                    title,
                    i_to_plot,
                    flavin_prominence=0.5,
                    flavin_lowpass_frac=0.25,
                    htb2_filter=False,
                )
            # convert res to two dataframes
            # ignore flavin_phase_at_ ...
            phases_df = pd.DataFrame.from_dict(
                {key: res[key] for key in res if "_at_" not in key}
            )
            phases_df.columns = phases_df.columns.str.replace("mCherry", "htb2")
            for key in ["omero_no", "condition", "group"]:
                phases_df[key] = locals()[key]
                period_dict[key].append(locals()[key])
            all_phases_dfs.append(phases_df)
            # dataframe for flavin_phase_at_ ... fields
            endpt = np.min([res[key].size for key in res if "_at_" in key])
            flavin_phase_df = pd.DataFrame.from_dict(
                {key: res[key][:endpt] for key in res if "_at_" in key}
            )
            # for a dataframe of all experiments
            add_to_summary_dict(summary_dict, data_dict, data_id, group, phases_df)
            make_single_experiment_plots(
                condition, phases_df, flavin_phase_df, period_dict, title
            )
            plot_periods(phases_df, title, xv="flavin_period", yv="buddings_period")


# combine summary statistics from all experiments
summary_df = pd.DataFrame.from_dict(summary_dict)
summary_df["id"] = summary_df.omero_no + " " + summary_df.group
plot_periods_ratios(summary_df)
plot_single_cell_data(all_phases_dfs)
catplot_from_period_dict(period_dict)

# Print drawings
pdf_filename = f"all_plots.pdf"
with PdfPages(pdf_filename) as pdf:
    for fig in range(1, plt.gcf().number + 1):
        pdf.savefig(fig)
# Close all figures
plt.close("all")
