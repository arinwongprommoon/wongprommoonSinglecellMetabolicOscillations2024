import matplotlib.pylab as plt
import numpy as np
import seaborn as sns
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar
from scipy.stats import pearsonr


def get_subdf(df, xv, yv):
    """Specialise to one entry per cell in the dataframe."""
    data = (
        df[
            [
                xv,
                yv,
                xv.split("_")[0] + "_oscillation_id",
                yv.split("_")[0] + "_oscillation_id",
                "cell_id",
            ]
        ]
        .drop_duplicates()
        .dropna()
    )
    return data


def find_1dmax(x, prob):
    """Find maximum of a 1D probability."""
    spl = CubicSpline(x, prob)
    res = minimize_scalar(lambda x: -spl(x), bounds=(x[0], x[-1]), method="bounded")
    if res.success:
        return res.x
    else:
        return -1


def find_2dmax(phases_df, xv="flavin_period", yv="htb2_period", nobins=20):
    """Estimate maximum of a period-period plot."""
    if xv in phases_df.columns and yv in phases_df.columns:
        data = (
            phases_df[
                [
                    xv,
                    yv,
                    xv.split("_")[0] + "_oscillation_id",
                    yv.split("_")[0] + "_oscillation_id",
                ]
            ]
            .drop_duplicates()
            .dropna()
        )
        h, xedges, yedges = np.histogram2d(
            data[xv].values,
            data[yv].values,
            range=[[0, 8], [0, 8]],
            bins=nobins,
        )
        probx = np.sum(h, axis=1) / np.sum(h)
        proby = np.sum(h, axis=0) / np.sum(h)
        xmidpoint = (xedges[:-1] + xedges[1:]) / 2
        ymidpoint = (yedges[:-1] + yedges[1:]) / 2
        x0 = find_1dmax(xmidpoint, probx)
        y0 = find_1dmax(ymidpoint, proby)
        varx = np.trapz((xmidpoint - x0) ** 2 * probx, xmidpoint)
        vary = np.trapz((ymidpoint - y0) ** 2 * proby, ymidpoint)
        res_dict = {
            "x_mean": x0,
            "y_mean": y0,
            "x_std": np.sqrt(varx),
            "y_std": np.sqrt(vary),
        }
    else:
        res_dict = {}
    return res_dict


def compare_periods_for_single_cells(phases_df, xv="flavin_period", yv="htb2_period"):
    """Find statistics for single-cell time series."""
    if xv in phases_df.columns and yv in phases_df.columns:
        data = get_subdf(phases_df, xv, yv)
        corr = []
        ratio = []
        for i in data.cell_id.unique():
            x_values = data[data.cell_id == i][xv].values
            y_values = data[data.cell_id == i][yv].values
            if np.unique(x_values).size > 1 and np.unique(y_values).size > 1:
                corr.append(pearsonr(x_values, y_values).statistic)
            ratio.append(np.log2(x_values / y_values))
        corr = np.array(corr)
        ratio = np.concatenate(ratio)
        return corr, ratio
    else:
        return np.nan, np.nan


def plot_torus(phases_df, title, xv="flavin_phase", yv="htb2_phase"):
    """Plot one phase versus another."""
    if "htb2_phase" in phases_df.columns:
        plt.figure()
        ax = sns.kdeplot(
            data=phases_df,
            x=xv,
            y=yv,
            fill=True,
            thresh=0,
            levels=100,
        )
        ax.set_xlim((0, 1))
        ax.set_ylim((0, 1))
        plt.title(title)
        plt.show(block=False)


def plot_periods(phases_df, title, xv=None, yv=None):
    """Plot one period versus another."""
    # plot periods
    if xv is None and yv is None:
        if "htb2_period" in phases_df.columns:
            xvyv = zip(
                ["flavin_period", "flavin_period", "htb2_period"],
                ["htb2_period", "buddings_period", "buddings_period"],
            )
        else:
            xvyv = zip(["flavin_period"], ["buddings_period"])
    else:
        xvyv = zip([xv], [yv])
    res = {}
    for xv, yv in xvyv:
        #
        # hist plot
        #
        data = get_subdf(phases_df, xv, yv)
        g = sns.jointplot(data, x=xv, y=yv, kind="hist")
        # set limits for both axes to be the same
        max_limit = max(g.ax_joint.get_xlim()[1], g.ax_joint.get_ylim()[1])
        g.ax_joint.set_xlim(0, max_limit)
        g.ax_joint.set_ylim(0, max_limit)
        # find mode
        twodg = find_2dmax(phases_df, xv, yv)
        res[f"{xv}"] = twodg["x_mean"]
        res[f"{xv}_std"] = twodg["x_std"]
        res[f"{yv}"] = twodg["y_mean"]
        res[f"{yv}_std"] = twodg["y_std"]
        # add mode to plot
        plt.plot(twodg["x_mean"], twodg["y_mean"], "ro")
        plt.gca().set_aspect("equal", adjustable="box")
        plt.suptitle(title)
        plt.show(block=False)
        #
        # kde plot
        #
        plt.figure()
        g = sns.kdeplot(data, x=xv, y=yv, fill=True)
        _, x_max = g.get_xlim()
        diag = np.linspace(0, x_max, 100)
        plt.plot(diag, diag, "k--", alpha=0.4)
        plt.suptitle(title)
        plt.tight_layout()
        plt.show(block=False)
    return res


def plot_flavin_phase(flavin_phase_df, title):
    """Plot flavin phase at peaks."""
    cols = sns.color_palette("tab10")
    plt.figure()
    ax = sns.histplot(data=flavin_phase_df, element="step", palette=cols[1:])
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    plt.title(title)
    plt.tight_layout()
    plt.show(block=False)


def raster_plot(
    df,
    signals,
    tstar,
    anchor_signal="flavin",
    signal_colour="r",
    title=None,
    xlim=None,
    cutoff=None,
):
    """
    Plot signal peaks relative to two consecutive peaks in the anchor_signal.

    Plot all cells following Bieler ... Naef, Mol Syst Biol, 2014, Fig 1B.

    Parameters
    ----------
    df: phases_df
    signal: List[str]
        Names of signal to plot.
    tstar: float
        Time at which to find flavin oscillations, which must be large
        enough for one full oscillation to have occurred.
    anchor_signal: str
        Name of anchor_signal used to find peaks with which to compare the
        peaks of signal.
    signal_colour : str
        Colour in which to plot signal peaks.
    title: str
        Title for plot.
    xlim: (float, float)
        Limits for the x-axis.
    cutoff: int, optional
        Display all cells from index cutoff onwards.
    """
    for signal in signals:
        if f"{signal}_period" in df.columns:
            df = df.dropna(subset=[f"{anchor_signal}_period", f"{signal}_period"])
        else:
            print(f"{signal} is unavailable.")
            return
    # specialise to particular anchor_signal oscillations
    peak_dict = {}
    # run through each cell time-series to check anchor_signal
    anchor_signal_period = np.empty(0)
    for lin in df.cell_id.unique():
        ldf = df[df.cell_id == lin]
        # find times of peaks in anchor_signal
        peak_times = np.sort(ldf[ldf[f"{anchor_signal}_phase"] == 0].time.values)
        if np.any(peak_times):
            tstar_i = np.argmin((peak_times - tstar) ** 2)
            if tstar_i > 0:
                peak_dict[lin] = {}
                peaks_tstar = peak_times[tstar_i - 1 : tstar_i + 1]
                anchor_signal_period = np.append(
                    anchor_signal_period, peaks_tstar[1] - peaks_tstar[0]
                )
                peak_dict[lin][f"{anchor_signal}_peaks"] = peaks_tstar
    longest_period = np.max(anchor_signal_period)
    # run through each cell time-series for the signals
    find_peaks_between_peaks(df, peak_dict, signals, anchor_signal, longest_period)
    # order by longest period
    lin_order = np.array(list(peak_dict.keys()))[np.argsort(anchor_signal_period)[::-1]]
    # lin_order = filter_cells_for_raster(
    #     lin_order, df, peak_dict, signals, anchor_signal
    # )
    if cutoff:
        lin_order = lin_order[cutoff:]
    # raster plot
    print(f" Plotting {lin_order.size} cells.")
    cols = sns.color_palette("tab10")
    plt.figure()
    y = 0
    for lin in lin_order:
        plt.plot(
            peak_dict[lin][f"{anchor_signal}_peaks"],
            y * np.ones(peak_dict[lin][f"{anchor_signal}_peaks"].size),
            ".",
            color=cols[0],
        )
        for c_i, signal in enumerate(signals):
            if np.any(peak_dict[lin][f"{signal}_peaks"]):
                plt.plot(
                    peak_dict[lin][f"{signal}_peaks"],
                    y * np.ones(peak_dict[lin][f"{signal}_peaks"].size),
                    ".",
                    color=cols[c_i + 1],
                    alpha=0.6,
                )
        y += 1
    plt.ylabel("cell index")
    plt.xlabel("time")
    if xlim is not None:
        plt.xlim(xlim)
    plt.ylim(bottom=-1)
    plt.grid()
    if title is None:
        plt.title(signal)
    else:
        plt.title(title)
    plt.show(block=False)
    # histograms of time of closest budding and htb2 max
    for i, key in enumerate(["buddings_peaks", "htb2_peaks"]):
        dist = [
            -peak_dict[lin][key][-1]
            for lin in peak_dict
            if key in peak_dict[lin] and np.any(peak_dict[lin][key])
        ]
        if np.any(dist):
            plt.figure()
            plt.hist(dist, color=cols[i + 1], alpha=0.4)
            plt.xlabel(f"time of {key[:-1]}")
            if title is None:
                plt.title(key)
            else:
                plt.title(title)
            plt.tight_layout()
            plt.show(block=False)
    return lin_order


def filter_cells_for_raster(
    lin_order, df, peak_dict, signals, anchor_signal, snr_quantile=None
):
    """Select cells based on qualities of cell-cycle and flavin signals."""
    # select cells using cell-cycle characteristics
    if "buddings" in signals and "htb2" in signals:
        selected_cells = select_cells_by_cell_cycle(peak_dict, anchor_signal)
        original_size = lin_order.size
        lin_order = np.array([cid for cid in lin_order if cid in selected_cells])
        print(
            f" Discarding {original_size - lin_order.size} cells for poor "
            "cell cycle."
        )
    # select cells using signal-to-noise of flavin signal
    if snr_quantile:
        snr = [df[df.cell_id == cid].flavin_snr.mean() for cid in df.cell_id.unique()]
        snr_min = np.quantile(snr, snr_quantile)
        selected_cells = [
            cid
            for cid in df.cell_id.unique()
            if df[df.cell_id == cid].flavin_snr.mean() > snr_min
        ]
        original_size = lin_order.size
        lin_order = np.array([cid for cid in lin_order if cid in selected_cells])
        print(
            f" Discarding {original_size - lin_order.size} cells for poor "
            "flavin signal."
        )
    return lin_order


def find_peaks_between_peaks(df, peak_dict, signals, anchor_signal, longest_period):
    """
    Find times of signal peaks given anchor_signal.

    Required for raster_plot.
    """
    for lin in peak_dict:
        ldf = df[df.cell_id == lin]
        l_ldf = ldf[
            (ldf.time > peak_dict[lin][f"{anchor_signal}_peaks"][1] - longest_period)
            & (ldf.time < peak_dict[lin][f"{anchor_signal}_peaks"][1])
        ]
        for signal in signals:
            # find times of peaks in signal
            peaks_tstar = l_ldf[l_ldf[f"{signal}_phase"] == 0].time.values
            # zero time relative to last peak in anchor_signal
            peak_dict[lin][f"{signal}_peaks"] = (
                peaks_tstar - peak_dict[lin][f"{anchor_signal}_peaks"][1]
            )
        # zero time for anchor_signal relative to last peak in anchor_signal
        peak_dict[lin][f"{anchor_signal}_peaks"] -= peak_dict[lin][
            f"{anchor_signal}_peaks"
        ][1]


def select_cells_by_cell_cycle(peak_dict, anchor_signal):
    """
    Select cells with one htb2 peak between buddings.

    Buddings are those closest to the anchor_signal.
    These cells likely have correctly identified peaks.
    """
    selected = []
    for lin in peak_dict:
        anchor_pks = peak_dict[lin][f"{anchor_signal}_peaks"]
        bud_pks = peak_dict[lin]["buddings_peaks"]
        htb2_pks = peak_dict[lin]["htb2_peaks"]
        if np.any(bud_pks) and np.any(htb2_pks):
            bud_r = np.argmin((bud_pks - anchor_pks.max()) ** 2)
            bud_l = np.argmin((bud_pks - anchor_pks.min()) ** 2)
            no_htb2_pks = np.array(
                [
                    (True if pk > bud_pks[bud_l] and pk < bud_pks[bud_r] else False)
                    for pk in htb2_pks
                ]
            )
            no_htb2_pks = no_htb2_pks[no_htb2_pks].size
            # select cells if there is one htb2 peak between buddings
            if no_htb2_pks == 1:
                selected.append(lin)
    return selected


def plot_time_series(cell_id, phases_df, title, tlim=None):
    """Plot peaks on signals."""
    df = phases_df[phases_df.cell_id == cell_id]
    signal_names = [field.split("_")[0] for field in df.columns if "_phase" in field]
    plt.figure(figsize=(10, 4))
    for i, signal_name in enumerate(signal_names):
        if signal_name in df:
            plt.subplot(len(signal_names), 1, i + 1)
            plt.plot(df.time.values, df[signal_name].values, ".-")
            pk_pts = np.where(df[f"{signal_name}_peaks"].values)
            plt.plot(
                df.time.values[pk_pts],
                df[signal_name].values[pk_pts],
                "ro",
            )
            plt.ylabel(signal_name)
            if tlim:
                plt.xlim(tlim)
    plt.suptitle(
        f"{title}: {int(cell_id):3d}" f" (snr = {np.mean(df.flavin_snr.values):.2f})"
    )
    plt.tight_layout()
    plt.show(block=False)
