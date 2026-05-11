"""
ae_picker.pickers
=================
First-arrival picking algorithms for acoustic emission (AE) waveforms.

All functions operate on raw NumPy arrays and are completely independent of
the file format — any reader that yields ``(time, amp)`` arrays is compatible.

Algorithms
----------
aic_picker                Maeda (1985) AIC criterion  — recommended primary picker
prepend_noise_aic_picker  AIC with synthetic pre-noise window prepended
envelope_onset_picker     Hilbert-envelope onset with MAD noise stats + persistence
envelope_offset_picker    Hilbert-envelope offset (end-of-event) detector
stalta_picker             Classical STA/LTA ratio picker
refined_stalta_picker     Recursive STA/LTA trigger → envelope-refined onset

Picking strategy (Maria dataset)
---------------------------------
* **Unfiltered signals** (record starts at onset, P arrives in first ~20 µs):
  Use ``aic_picker`` with ``search_end ≈ 10 samples`` (50 µs at 200 kHz).

* **Filtered signals** (record centred at t = 0, ±250 µs window):
  Use ``aic_picker`` over the central 40 % of the record.
"""

import numpy as np
from scipy.signal import hilbert


# ── private helpers ───────────────────────────────────────────────────────────

def _moving_average(x, n):
    """Centered moving average with edge-padding; returns same length as x."""
    n = int(max(1, n))
    if n == 1:
        return x.copy()
    pad  = n // 2
    xpad = np.pad(x, (pad, pad), mode="edge")
    w    = np.ones(n, dtype=float) / n
    y    = np.convolve(xpad, w, mode="valid")
    return y[:len(x)]


# ── public pickers ────────────────────────────────────────────────────────────

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


def prepend_noise_aic_picker(amp, time, n_prepend=None, prepend_duration_s=50e-6,
                              search_start=None, search_end=None, rng_seed=0):
    """
    AIC picker with a synthetic pre-onset noise window prepended.

    For records that start at the P onset (e.g. unfiltered format), the AIC
    picker has almost no pre-onset samples to anchor the noise variance
    estimate, causing the minimum to drift to the search boundary.  This
    function prepends ``n_prepend`` samples of Gaussian noise whose amplitude
    is matched to the RMS of the first few real samples, giving AIC a stable
    noise baseline before the onset.

    The returned ``pick_idx`` is corrected back to the original (non-prepended)
    sample index, so it can be used directly with the original ``amp`` and
    ``time`` arrays.

    Parameters
    ----------
    amp : array-like
        Single-channel amplitude time series.
    time : array-like
        Corresponding time axis (seconds).
    n_prepend : int or None
        Number of synthetic noise samples to prepend.  If ``None``,
        ``prepend_duration_s / dt`` is used.
    prepend_duration_s : float
        Duration of the prepended noise window (seconds).  Used only when
        ``n_prepend`` is ``None``.  Default: 50 µs.
    search_start : int or None
        First candidate sample in the *padded* array.  Defaults to
        ``n_prepend`` so the search starts where real signal begins.
    search_end : int or None
        Last candidate sample in the *padded* array (inclusive).  Defaults
        to ``n_prepend + int(0.005 * N_orig)`` — same 0.5 % fraction used
        by the batch runner for unfiltered data.
    rng_seed : int
        Seed for the random number generator (ensures reproducibility).
        Default: 0.

    Returns
    -------
    pick_idx : int
        Estimated onset sample index in the **original** (unpadded) array.
        Returns 0 if the AIC minimum falls inside the prepended region.
    aic_padded : np.ndarray
        Full AIC function on the padded array (length ``n_prepend + N_orig``).
    n_prepend : int
        Number of samples that were prepended (useful for interpreting the
        AIC curve).
    """
    x  = np.asarray(amp,  dtype=float)
    t  = np.asarray(time, dtype=float)
    dt = float(np.median(np.diff(t)))
    N  = len(x)

    if n_prepend is None:
        n_prepend = max(10, int(round(prepend_duration_s / dt)))

    n_ref     = max(2, min(int(0.05 * N), 20))
    noise_rms = float(np.std(x[:n_ref]))
    if noise_rms == 0:
        noise_rms = float(np.std(x)) * 0.01

    rng   = np.random.default_rng(rng_seed)
    noise = rng.normal(0.0, noise_rms, size=n_prepend)
    x_pad = np.concatenate([noise, x])

    if search_start is None:
        search_start = n_prepend
    if search_end is None:
        search_end = n_prepend + max(2, int(0.005 * N))

    _, aic_padded = aic_picker(x_pad, search_start=search_start,
                               search_end=search_end)

    seg = aic_padded[search_start : search_end + 1]
    if not np.any(~np.isnan(seg)):
        return 0, aic_padded, n_prepend

    padded_idx = search_start + int(np.nanargmin(seg))
    pick_idx   = max(0, padded_idx - n_prepend)
    return pick_idx, aic_padded, n_prepend



def envelope_onset_picker(amp, time, search_start_s=None, search_end_s=None,
                           k_high=3.0, k_low=3.0,
                           min_persist_ms=0.3, noise_lookback_s=40e-6):
    """
    Hilbert-envelope onset picker with MAD noise statistics and persistence gating.

    Detects the first arrival by:

    1. Computing the Hilbert envelope of ``amp``.
    2. Estimating noise level from a pre-search window using the median
       absolute deviation (MAD) — robust to outliers.
    3. Confirming an event exists by finding the first sample where the
       envelope exceeds ``med + k_high × σ``.
    4. Walking backward from that sample to find the first point that is
       continuously above ``med + k_low × σ`` for ``min_persist_ms``.

    Parameters
    ----------
    amp : array-like
        Single-channel amplitude time series.
    time : array-like
        Corresponding time axis (seconds).
    search_start_s : float or None
        Search window start time (seconds).  Defaults to ``time[0]``.
    search_end_s : float or None
        Search window end time (seconds).  Defaults to ``time[-1]``.
    k_high : float
        High-threshold multiplier (noise σ above median).  Used to confirm
        a detectable event exists before walking backward to the onset.
        Default: 3.0.
    k_low : float
        Low-threshold multiplier.  The onset is the earliest sample at which
        the envelope is continuously above ``med + k_low × σ`` for at least
        ``min_persist_ms``.  Default: 3.0.
    min_persist_ms : float
        Minimum duration (ms) the envelope must stay above ``thr_low`` to
        confirm the onset.  Prevents triggering on transient noise spikes.
        Default: 0.3 ms.
    noise_lookback_s : float
        Duration of the noise estimation window immediately preceding
        ``search_start_s``.  Default: 40 µs.

    Returns
    -------
    pick_idx : int or None
        Sample index of the estimated onset, or ``None`` if no event is
        detected within the search window.
    """
    x  = np.asarray(amp, dtype=float)
    t  = np.asarray(time, dtype=float)
    dt = float(np.median(np.diff(t)))
    fs = 1.0 / dt
    N  = len(x)

    i0 = int(round((search_start_s - t[0]) / dt)) if search_start_s is not None else 0
    i1 = int(round((search_end_s   - t[0]) / dt)) if search_end_s   is not None else N
    i0 = max(0, min(i0, N - 1))
    i1 = max(i0 + 2, min(i1, N))

    env = np.abs(hilbert(x))

    pre0  = max(0, i0 - int(noise_lookback_s * fs))
    noise = env[pre0:i0] if i0 > pre0 else env[:max(1, i0)]
    med   = np.median(noise)
    mad   = np.median(np.abs(noise - med)) + 1e-12
    sigma = 1.4826 * mad

    thr_high = med + k_high * sigma
    thr_low  = med + k_low  * sigma

    seg        = env[i0:i1]
    above_high = seg > thr_high
    if not above_high.any():
        return None
    i_hi = i0 + int(np.argmax(above_high))

    persist = max(1, int((min_persist_ms / 1000.0) * fs))
    j = i_hi
    while j > i0 + persist:
        if np.all(env[j - persist : j] > thr_low):
            j -= 1
        else:
            break

    return j


def envelope_offset_picker(amp, time, onset_idx, search_end_s=None,
                            k_low=3.0, end_persist_ms=1.0,
                            noise_lookback_s=30e-6):
    """
    Hilbert-envelope offset (end-of-event) picker.

    Starting from ``onset_idx``, moves forward until the envelope stays
    continuously below ``med + k_low × σ`` for ``end_persist_ms``.

    Parameters
    ----------
    amp : array-like
        Single-channel amplitude time series.
    time : array-like
        Corresponding time axis (seconds).
    onset_idx : int
        Sample index of the event onset (from e.g. ``envelope_onset_picker``).
    search_end_s : float or None
        Latest time (seconds) to search for the offset.  Defaults to the
        end of the record.
    k_low : float
        Low-threshold multiplier for offset detection.  Default: 3.0.
    end_persist_ms : float
        Duration (ms) the envelope must remain below ``thr_low`` to confirm
        the offset.  Default: 1.0 ms.
    noise_lookback_s : float
        Duration of the noise estimation window before the onset.  Default: 30 µs.

    Returns
    -------
    offset_idx : int
        Sample index of the estimated event offset.  Returns the end of the
        search window if the envelope never drops below threshold.
    """
    x  = np.asarray(amp,  dtype=float)
    t  = np.asarray(time, dtype=float)
    dt = float(np.median(np.diff(t)))
    fs = 1.0 / dt
    N  = len(x)

    i0 = int(onset_idx)
    i1 = int(round((search_end_s - t[0]) / dt)) if search_end_s is not None else N
    i1 = max(i0 + 2, min(i1, N))

    env  = np.abs(hilbert(x))
    pre0 = max(0, i0 - int(noise_lookback_s * fs))
    noise = env[pre0:i0] if i0 > pre0 else env[:max(1, i0)]

    med   = np.median(noise)
    mad   = np.median(np.abs(noise - med)) + 1e-12
    sigma = 1.4826 * mad
    thr_low = med + k_low * sigma

    persist = max(1, int((end_persist_ms / 1000.0) * fs))
    for j in range(i0, i1 - persist):
        if np.all(env[j : j + persist] < thr_low):
            return j

    return i1


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


def refined_stalta_picker(amp, time,
                           sta_s=5e-6, lta_s=50e-6,
                           sta_thresh=2.5, lta_thresh=1.5,
                           k_high=6.0, k_low=2.5,
                           min_persist_ms=0.3,
                           pre_pick_s=50e-6, post_pick_s=200e-6,
                           post_offset_s=500e-6, end_persist_ms=0.05,
                           search_start=0):
    """
    Recursive STA/LTA trigger followed by envelope-refined onset picking.

    This two-stage approach combines the sensitivity of the recursive
    (causal, exponentially weighted) STA/LTA characteristic function with
    the precision of envelope-based hysteresis/persistence onset detection:

    1. **STA/LTA stage** — computes the recursive characteristic function
       (using ObsPy's ``recursive_sta_lta`` when available; falls back to a
       simple causal moving-average ratio otherwise) and extracts trigger
       on/off indices via threshold crossing.
    2. **Refinement stage** — for each trigger, calls
       ``envelope_onset_picker`` over a tight window around the STA/LTA
       onset to walk back to the true first motion.

    Parameters
    ----------
    amp : array-like
        Single-channel amplitude time series.
    time : array-like
        Corresponding time axis (seconds).
    sta_s : float
        Short-term average window length (seconds).  Default: 5 µs.
    lta_s : float
        Long-term average window length (seconds).  Default: 50 µs.
    sta_thresh : float
        STA/LTA ratio to trigger ON.  Default: 2.5.
    lta_thresh : float
        STA/LTA ratio to trigger OFF.  Default: 1.5.
    k_high : float
        High-threshold multiplier for ``envelope_onset_picker``.  Default: 6.0.
    k_low : float
        Low-threshold multiplier for ``envelope_onset_picker``.  Default: 2.5.
    min_persist_ms : float
        Minimum persistence duration (ms) for envelope onset.  Default: 0.3 ms.
    pre_pick_s : float
        Time (seconds) to look before the STA/LTA trigger when setting the
        envelope search window.  Default: 50 µs.
    post_pick_s : float
        Time (seconds) to look after the STA/LTA trigger.  Default: 200 µs.
    post_offset_s : float
        Maximum duration (seconds) to search for the event offset after the
        refined onset.  Default: 500 µs.
    end_persist_ms : float
        Persistence duration (ms) for offset detection.  Default: 0.05 ms.
    search_start : int
        First sample index at which a trigger is allowed to fire.  Triggers
        whose onset falls before this index are discarded.  Useful for
        filtered records centred at t=0 where arrivals cannot precede a known
        index (e.g. pass the index of t=0 to prevent pre-onset false triggers).
        Default: 0 (no restriction).

    Returns
    -------
    pick_idx : int or None
        Sample index of the refined first arrival, or ``None`` if no trigger
        was found.
    cft : np.ndarray, shape (N,)
        STA/LTA characteristic function used for triggering.
    triggers : list of [int, int]
        Raw ``[onset, offset]`` trigger pairs from the STA/LTA stage (sample
        indices).  Empty list if no triggers were found.
    """
    x  = np.asarray(amp,  dtype=float)
    t  = np.asarray(time, dtype=float)
    dt = float(np.median(np.diff(t)))
    fs = 1.0 / dt
    N  = len(x)

    nsta = max(1, int(round(sta_s / dt)))
    nlta = max(nsta + 1, int(round(lta_s / dt)))

    # ── STA/LTA characteristic function ──────────────────────────────────
    try:
        from obspy.signal.trigger import recursive_sta_lta, trigger_onset
        cft          = recursive_sta_lta(x, nsta, nlta)
        raw_triggers = [list(pair) for pair in trigger_onset(cft, sta_thresh, lta_thresh)]
    except ImportError:
        # Fallback: causal moving-average ratio (same as stalta_picker)
        cft = np.zeros(N)
        for i in range(nlta, N):
            sta_val = np.mean(np.abs(x[i - nsta : i]))
            lta_val = np.mean(np.abs(x[i - nlta : i - nsta]))
            cft[i]  = sta_val / lta_val if lta_val > 0 else 0.0

        in_trig, raw_triggers, onset = False, [], None
        for i in range(N):
            if not in_trig and cft[i] >= sta_thresh:
                in_trig, onset = True, i
            elif in_trig and cft[i] <= lta_thresh:
                raw_triggers.append([onset, i])
                in_trig = False
        if in_trig and onset is not None:
            raw_triggers.append([onset, N - 1])

    # Discard triggers that fired before the allowed search window
    raw_triggers = [tr for tr in raw_triggers if tr[0] >= search_start]

    if not raw_triggers:
        return None, cft, raw_triggers

    # ── Envelope refinement of the first trigger ──────────────────────────
    onset_raw = raw_triggers[0][0]
    i_start   = max(0, onset_raw - int(pre_pick_s  * fs))
    i_end     = min(N, onset_raw + int(post_pick_s * fs))

    refined = envelope_onset_picker(
        amp, time,
        search_start_s = t[i_start],
        search_end_s   = t[i_end - 1],
        k_high=k_high, k_low=k_low,
        min_persist_ms=min_persist_ms,
    )
    pick_idx = refined if refined is not None else onset_raw

    return pick_idx, cft, raw_triggers
