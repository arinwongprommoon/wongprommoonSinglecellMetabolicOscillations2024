"""Use peaks in signals to define phases."""

from dataclasses import dataclass

import matplotlib.pylab as plt
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from wela.fft_filtering import flavin_filter

from peaks import bud_to_bud_distance, get_flavin_peaks, get_htb2_peaks

# samples taken every 5 minutes
fs = 1 / 5


def get_signals(
    signal_names,
    dl,
    group,
    odf,
    title,
    i_to_plot,
    flavin_prominence,
    flavin_lowpass_frac,
    htb2_filter=False,
):
    """Define results dict and call get_all_phases."""
    res = {"cell_id": [], "time": [], "flavin_snr": []}
    signals = {}
    for signal_name in signal_names:
        t, signals[signal_name] = dl.get_time_series(signal_name, group=group, df=odf)
        # convert to hours
        if np.max(t) > 48:
            t = t / 60
        # filter signal
        if signal_name != "buddings":
            signals[signal_name] = flavin_filter(signals[signal_name])
        # interpolate NaN between data points - necessary for find_peaks
        signals[signal_name] = (
            pd.DataFrame(signals[signal_name]).interpolate(axis=1).to_numpy()
        )
        res[f"{signal_name}"] = []
        res[f"{signal_name}_peaks"] = []
        res[f"{signal_name}_phase"] = []
        res[f"{signal_name}_period"] = []
        res[f"{signal_name}_oscillation_id"] = []
    res["flavin_phase_at_budding"] = []
    if "mCherry" in signal_names:
        res["flavin_phase_at_htb2"] = []
    distance = bud_to_bud_distance(t, dl, group, odf)
    print(f" bud-to-bud distance = {distance} units")
    # fill res dict
    get_all_phases(
        res,
        t,
        signals,
        distance,
        i_to_plot,
        title,
        flavin_prominence,
        flavin_lowpass_frac,
        htb2_filter,
    )
    return res


def get_all_phases(
    res,
    t,
    signals,
    distance,
    i_to_plot,
    title,
    flavin_prominence,
    flavin_lowpass_frac,
    htb2_filter,
):
    """Find peaks and phase and period for all signals for all cells."""
    nocells = signals["flavin"].shape[0]
    for i in range(nocells):
        lres = {}
        for signal_name in signals:
            # single cell data
            data = signals[signal_name][i, :]
            # remove any remaining NaN at the beginning of the time series
            cdata = data[~np.isnan(data)]
            # find peaks
            if signal_name == "flavin":
                peaks, snr = get_flavin_peaks(
                    cdata,
                    distance,
                    fs,
                    flavin_prominence,
                    flavin_lowpass_frac,
                    # snr_cutoff_freq=_freq=1 / (2 * distance * 1 / fs),
                )
            elif signal_name == "buddings":
                peaks, _ = find_peaks(cdata, distance=distance / 3)
            elif signal_name == "mCherry":
                peaks = get_htb2_peaks(cdata, distance, fs=fs, htb2_filter=htb2_filter)
            if np.any(peaks):
                if np.any(np.isnan(data)):
                    # increment indices because of initial NaNs
                    peaks += 1 + np.where(np.diff(np.isnan(data)))[0][0]
                # get phase and period for a lineage
                phase_res = get_phase(t, peaks, cell_index=i)
                lres[f"{signal_name}_phase"] = phase_res.phase
                lres[f"{signal_name}_period"] = phase_res.period
                lres[f"{signal_name}_oscillation_id"] = phase_res.oscillation_id
                lres[f"{signal_name}_peaks"] = peaks
                lres[signal_name] = data
                lres["cell_id"] = phase_res.cell_id
                lres["flavin_snr"] = snr * np.ones(data.shape)
        # check peaks in all signals
        if np.all(
            [np.any(lres.get(f"{signal_name}_peaks", False)) for signal_name in signals]
        ):
            # add lres to res
            for signal_name in signals:
                res[f"{signal_name}"].append(lres[f"{signal_name}"])
                # make peaks binary
                s_peaks = np.zeros(len(lres[f"{signal_name}"]))
                s_peaks[lres[f"{signal_name}_peaks"]] = 1
                res[f"{signal_name}_peaks"].append(s_peaks)
                res[f"{signal_name}_phase"].append(lres[f"{signal_name}_phase"])
                res[f"{signal_name}_period"].append(lres[f"{signal_name}_period"])
                res[f"{signal_name}_oscillation_id"].append(
                    lres[f"{signal_name}_oscillation_id"]
                )
            res["cell_id"].append(lres["cell_id"])
            res["time"].append(t)
            res["flavin_snr"].append(lres["flavin_snr"])
            # find flavin phases at budding
            if "flavin_phase" in lres and "buddings_peaks" in lres:
                res["flavin_phase_at_budding"].append(
                    lres["flavin_phase"][lres["buddings_peaks"]]
                )
            # find flavin phase at Htb2 maxima
            if "flavin_phase" in lres and "mCherry_peaks" in lres:
                res["flavin_phase_at_htb2"].append(
                    lres["flavin_phase"][lres["mCherry_peaks"]]
                )
            # plot peaks
            if i in i_to_plot:
                plot_time_series_from_lres(t, lres, title, i)
    for key in res:
        res[key] = np.concatenate(res[key])


def get_phase(t, peaks, cell_index):
    """
    Find phase and period from peaks in a single time series.

    Both are 1D arrays of the same shape as input t to allow
    comparison at particular time points

    The period is therefore repeated multiple times.
    """

    @dataclass
    class phase_result:
        phase: np.ndarray
        period: np.ndarray
        oscillation_id: np.ndarray
        cell_id: np.ndarray

    dt = np.median(np.diff(np.round(t, 2)))
    phase = np.empty(0)
    period = np.empty(0)
    oscillation_id = np.empty(0)
    for ind, i in enumerate(range(len(peaks) - 1)):
        local_t = t[peaks[i] : peaks[i + 1]]
        # period
        local_period = len(local_t) * np.ones(len(local_t))
        period = np.append(period, local_period)
        # phase in units of pi
        local_phase = local_t - local_t[0]
        local_phase = local_phase / (np.max(local_phase) + dt)
        phase = np.append(phase, local_phase)
        # define a unique ID for each oscillation
        # prevents dropduplicates dropping multiple oscillations
        local_id = ind * np.ones(len(local_t))
        oscillation_id = np.append(oscillation_id, local_id)
    # extend to length of t
    full_phase = np.nan * np.ones(t.size)
    full_phase[peaks[0] : peaks[-1]] = phase
    full_period = np.nan * np.ones(t.size)
    full_period[peaks[0] : peaks[-1]] = period * dt
    full_oscillation_id = np.nan * np.ones(t.size)
    full_oscillation_id[peaks[0] : peaks[-1]] = oscillation_id
    cell_id = np.ones(t.size) * cell_index
    return phase_result(
        phase=full_phase,
        period=full_period,
        oscillation_id=full_oscillation_id,
        cell_id=cell_id,
    )


def plot_time_series_from_lres(t, lres, title, i_timeseries):
    """Plot peaks on signals."""
    if "mCherry" in lres:
        nosubplots = 3
        signal_names = ["flavin", "buddings", "mCherry"]
    else:
        nosubplots = 2
        signal_names = ["flavin", "buddings"]
    plt.figure(figsize=(10, 4))
    for i, signal_name in enumerate(signal_names):
        if signal_name in lres:
            plt.subplot(nosubplots, 1, i + 1)
            plt.plot(t, lres[signal_name])
            plt.plot(
                t[lres[f"{signal_name}_peaks"]],
                lres[signal_name][lres[f"{signal_name}_peaks"]],
                "ro",
            )
            if signal_name == "mCherry":
                plt.ylabel("htb2")
            else:
                plt.ylabel(signal_name)
    plt.suptitle(f"{title}: {i_timeseries}")
    plt.tight_layout()
    plt.show(block=False)
