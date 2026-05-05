"""
ae_picker.pickers
=================
First-arrival picking algorithms for acoustic emission (AE) waveforms.

All functions operate on raw NumPy arrays and are completely independent of
the file format — any reader that yields ``(time, amp)`` arrays is compatible.

Algorithms
----------
aic_picker          Maeda (1985) AIC criterion  — recommended primary picker
energy_onset_picker Energy-ratio onset detector  — robust fallback for short records
stalta_picker       Classical STA/LTA ratio      — configurable trigger picker

Picking strategy (Maria dataset)
---------------------------------
* **Unfiltered signals** (record starts at onset, P arrives in first ~20 µs):
  Use ``aic_picker`` with ``search_end ≈ 10 samples`` (50 µs at 200 kHz).
  Expanding the window causes the global AIC minimum to drift to later
  high-energy coda features.  Use ``energy_onset_picker`` as a cross-check.

* **Filtered signals** (record centred at t = 0, ±250 µs window):
  Use ``aic_picker`` over the central 40 % of the record.
"""

import numpy as np


def aic_picker(amp, search_start=0, search_end=None):
    """
    Maeda (1985) AIC first-arrival picker.

    The AIC of a split point *k* measures how well the signal before and
    after *k* are each described by a stationary Gaussian model::

        AIC(k) = k · log(var(x[0 : k+1])) + (N−k−1) · log(var(x[k+1 : N]))

    The sample index that minimises AIC within the search window is the
    estimated onset of the first arrival.

    Parameters
    ----------
    amp : array-like
        Single-channel amplitude time series (any units).
    search_start : int
        First sample index to consider as a candidate onset.
        Skip at least 1 sample so that the pre-onset variance estimate
        is based on at least one data point.
    search_end : int or None
        Last sample index to consider (inclusive).  If ``None``, the
        search extends to ``N − 2``.

        **Tuning guidance** — restrict this window carefully:

        * For records pre-windowed to start at the P onset
          (e.g. unfiltered format), use ``search_end ≈ 10`` samples
          (50 µs at 200 kHz).  Larger windows cause the minimum to
          migrate to later coda peaks.
        * For centred records (filtered format), use
          ``search_end = int(0.70 * N)`` with
          ``search_start = int(0.30 * N)``.

    Returns
    -------
    pick_idx : int
        Sample index of the estimated first arrival.
    aic : np.ndarray, shape (N,)
        Full AIC function; ``np.nan`` outside ``[search_start, search_end]``.

    References
    ----------
    Maeda, N. (1985). A method for reading and checking phase times in
    auto-processing system of seismic wave data.
    *Zisin (Journal of the Seismological Society of Japan)*, 38, 365–379.
    """
    x = np.asarray(amp, dtype=float)
    N = len(x)

    if search_end is None:
        search_end = N - 2
    search_end = min(search_end, N - 2)

    aic = np.full(N, np.nan)
    for k in range(search_start, search_end + 1):
        var_pre  = np.var(x[: k + 1])
        var_post = np.var(x[k + 1 :])
        if var_pre > 0 and var_post > 0:
            aic[k] = (k + 1) * np.log(var_pre) + (N - k - 2) * np.log(var_post)

    seg = aic[search_start : search_end + 1]
    if not np.any(~np.isnan(seg)):
        return search_start, aic

    pick_idx = search_start + int(np.nanargmin(seg))
    return pick_idx, aic


def energy_onset_picker(amp, time, noise_window_s=25e-6, threshold=3.0,
                         search_start=1):
    """
    Energy-ratio first-arrival picker.

    Returns the first sample index at which the instantaneous energy
    ``amp[i]²`` exceeds ``threshold`` times the background energy level
    estimated from the pre-onset noise window.

    This picker is robust when the record begins immediately at the onset
    (leaving too few samples for reliable AIC variance estimation) or when
    the signal-to-noise ratio is low.

    Parameters
    ----------
    amp : array-like
        Single-channel amplitude time series.
    time : array-like
        Corresponding time axis (seconds).  Used only to convert
        ``noise_window_s`` to a sample count.
    noise_window_s : float
        Duration of the pre-onset noise window used to estimate background
        energy (seconds).  Default: 25 µs.
    threshold : float
        Energy trigger ratio: pick fires when ``amp[i]² > threshold × bg``.
        Default: 3.0 (energy must be 3× the background level).
    search_start : int
        First sample index to evaluate.  Default: 1 (skip index 0 to
        avoid triggering on the very first sample).

    Returns
    -------
    pick_idx : int or None
        Sample index of the first detected onset, or ``None`` if the
        threshold is never exceeded within the record.
    """
    x  = np.asarray(amp, dtype=float)
    t  = np.asarray(time, dtype=float)
    dt = float(np.median(np.diff(t)))

    n_noise = max(2, int(round(noise_window_s / dt)))
    bg = np.mean(x[:n_noise] ** 2)
    if bg == 0:
        return None

    for i in range(search_start, len(x)):
        if x[i] ** 2 > threshold * bg:
            return i
    return None


def stalta_picker(amp, time, sta_s, lta_s, threshold=3.0, search_start=0):
    """
    Classical STA/LTA (Short-Term Average / Long-Term Average) first-arrival picker.

    Computes the ratio of short-term to long-term absolute amplitude averages
    at each sample and returns the first index where the ratio exceeds the
    trigger threshold.

    Parameters
    ----------
    amp : array-like
        Single-channel amplitude time series.
    time : array-like
        Corresponding time axis (seconds).  Used to convert window durations
        to sample counts.
    sta_s : float
        Short-term average window length (seconds).
        Typical values: 5–20 µs for AE data at 200 kHz.
    lta_s : float
        Long-term average window length (seconds).
        Typical values: 50–200 µs for AE data at 200 kHz.
    threshold : float
        STA/LTA ratio that triggers a pick.  Default: 3.0.
    search_start : int
        First sample index to search for a trigger.

    Returns
    -------
    pick_idx : int or None
        First sample where ``ratio ≥ threshold``, or ``None`` if the
        threshold is never reached.
    ratio : np.ndarray, shape (N,)
        Full STA/LTA ratio array (zero before the first valid LTA window).
    """
    x   = np.abs(np.asarray(amp, dtype=float))
    t   = np.asarray(time, dtype=float)
    dt  = float(np.median(np.diff(t)))

    sta_n = max(1, int(round(sta_s / dt)))
    lta_n = max(sta_n + 1, int(round(lta_s / dt)))
    N     = len(x)

    ratio = np.zeros(N)
    for i in range(lta_n, N):
        sta = np.mean(x[i - sta_n : i])
        lta = np.mean(x[i - lta_n : i - sta_n])
        ratio[i] = sta / lta if lta > 0 else 0.0

    pick_idx = None
    for i in range(max(search_start, lta_n), N):
        if ratio[i] >= threshold:
            pick_idx = i
            break

    return pick_idx, ratio
