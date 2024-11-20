import matplotlib.pylab as plt
import numpy as np
from scipy.signal import find_peaks
from wela.fft_filtering import classic_periodogram, lowpass_filter, snr


def bud_to_bud_distance(t, dl, group, odf):
    """Estimate cell-cycle duration using buddings."""
    _, buddings = dl.get_time_series("buddings", group=group, df=odf)
    dt = np.median(np.round(np.diff(t), 2))
    tm = np.repeat(
        t.reshape(
            1,
            t.size,
        ),
        buddings.shape[0],
        axis=0,
    )
    b2b_dist = np.diff(tm[buddings.astype(bool)]) / dt
    b2b_dist = b2b_dist.astype(int)
    distance = np.median(b2b_dist)
    return distance


def get_flavin_peaks(
    cdata, distance, fs, prominence=1, lowpass_frac=0.3, snr_cutoff_freq=None
):
    """Find peaks in a single time series of flavin fluorescence."""
    nyquist = fs / 2
    filt_data = lowpass_filter(
        cdata - np.mean(cdata), lowpass_frac * nyquist, fs=fs, N=4
    )
    distance = np.max([4, distance / 12])
    peaks, _ = find_peaks(filt_data, distance=distance, prominence=prominence)
    if False:
        # for debugging
        plt.figure(figsize=(10, 4))
        plt.plot(cdata, ".-", alpha=0.5)
        plt.plot(filt_data)
        plt.plot(peaks, cdata[peaks], "ro")
        plt.show(block=False)
        breakpoint()
    # noise in number of predicted peaks - works poorly
    pks = [
        find_peaks(filt_data, distance=distance, prominence=prom)[0].size
        for prom in np.linspace(0.8, 1.2, 10) * prominence
    ]
    if np.any(pks):
        eta = np.std(pks) / np.mean(pks)
    else:
        eta = np.nan
    # use signal-to-noise - works poorly too
    if snr_cutoff_freq:
        freqs, power = classic_periodogram(cdata, 1 / fs)
        eta = snr(freqs, power, snr_cutoff_freq)
    else:
        eta = np.nan
    return peaks, eta


def get_htb2_peaks(cdata, distance, htb2_filter=False, fs=None):
    """Find peaks in a single time series of Htb2 fluorescence."""
    if htb2_filter and fs is not None:
        cdata = lowpass_filter(cdata - np.mean(cdata), 0.05 * fs, fs=fs, N=4)
    nsig = cdata - np.min(cdata)
    nsig /= np.max(nsig)
    # find all peaks
    # no distance because we want the right of two plateau peaks
    peaks, props = find_peaks(nsig, prominence=0.01)
    # find peaks with a prominent drop on the right
    right_prominence = nsig[peaks] - nsig[props["right_bases"]]
    # find drop from peak to minimum before next peak in time
    drop_to_next_min = np.array(
        [
            nsig[peak1] - np.min(nsig[peak1:peak2])
            for peak1, peak2 in zip(peaks[:-1], peaks[1:])
        ]
    )
    drop_to_next_min = np.append(drop_to_next_min, 1)
    # peaks have sharp right drop of sufficient height
    peaks = peaks[(right_prominence > 0.3) & (drop_to_next_min > 0.3)]
    # sort by distance
    new_sig = np.zeros(nsig.size)
    new_sig[peaks] = 1
    if distance is not None:
        new_peaks, _ = find_peaks(new_sig, distance=distance / 4)
    else:
        new_peaks, _ = find_peaks(new_sig)
    return new_peaks
