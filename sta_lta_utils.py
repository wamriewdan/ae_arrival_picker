import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from obspy.signal.trigger import recursive_sta_lta, trigger_onset

def pick_onset_from_envelope(
    x, fs,
    search_start, search_end,
    k_high=3.0, k_low=3.0,
    min_persist_ms=0.3,
    noise_lookback_s=0.04,
):
    """
    Envelope-based onset picker using hysteresis thresholds derived from robust noise stats (MAD)
    + persistence.

    k_high : float
        High-threshold multiplier (noise sigma above noise median) used to confirm an event
        exists in the search window before walking backward to the onset.
        Lower values (e.g. 3.0) allow the first confirmed crossing to occur earlier on the
        rising edge, giving more accurate onset times. This is safe because the STA/LTA stage
        already gates the search window to a known-event region — k_high only needs to beat
        within-window noise fluctuations, not reject unrelated noise bursts.

    Returns onset sample index or None.
    """
    x = np.asarray(x, dtype=float)
    env = np.abs(hilbert(x))

    i0 = max(0, int(search_start))
    i1 = min(len(x), int(search_end))
    if i1 <= i0 + 2:
        return None

    # Noise stats from a pre-event window before i0
    pre0 = max(0, i0 - int(noise_lookback_s * fs))
    pre1 = i0
    noise = env[pre0:pre1] if pre1 > pre0 else env[:max(1, i0)]

    med = np.median(noise)
    mad = np.median(np.abs(noise - med)) + 1e-12
    sigma = 1.4826 * mad  # robust std estimate

    seg = env[i0:i1]

    if k_high is None:
        k_high = 3.0
    thr_high = med + k_high * sigma
    thr_low  = med + k_low  * sigma

    # 1) confirm event: first exceed high threshold
    above_high = seg > thr_high
    if not above_high.any():
        return None
    idx_hi = np.argmax(above_high)
    i_hi = i0 + idx_hi

    # 2) walk backward to find first sustained low-threshold crossing
    persist = max(1, int((min_persist_ms / 1000) * fs))
    j = i_hi
    while j > i0 + persist:
        if np.all(env[j - persist:j] > thr_low):
            j -= 1
        else:
            break

    return j

def pick_offset_from_envelope(
    x, fs,
    start_from,
    search_end,
    k_low=3.0,
    end_persist_ms=1.0,
    noise_lookback_s=0.03
):
    """
    Simple offset picker: from start_from, move forward until envelope stays below thr_low
    for end_persist_ms. Returns offset sample index.
    """
    x = np.asarray(x, dtype=float)
    env = np.abs(hilbert(x))

    i0 = max(0, int(start_from))
    i1 = min(len(x), int(search_end))
    if i1 <= i0 + 2:
        return i1

    # Estimate noise from a window just before i0
    pre0 = max(0, i0 - int(noise_lookback_s * fs))
    pre1 = i0
    noise = env[pre0:pre1] if pre1 > pre0 else env[:max(1, i0)]

    med = np.median(noise)
    mad = np.median(np.abs(noise - med)) + 1e-12
    sigma = 1.4826 * mad
    thr_low = med + k_low * sigma

    persist = max(1, int((end_persist_ms / 1000) * fs))

    # Find first index where env remains below thr_low for 'persist' samples
    for j in range(i0, i1 - persist):
        if np.all(env[j:j + persist] < thr_low):
            return j

    return i1

def refined_sta_lta_detection(
    tr, fs,
    sta_thresh=2.5, lta_thresh=1.5,
    sta_window=0.0005, lta_window=0.005,
    # envelope picker params
    k_high=6.0, k_low=2.5,
    min_persist_ms=0.3,
    # search window around trigger onset
    pre_pick_s=0.002, post_pick_s=0.008,
    # offset control
    post_offset_s=0.10, end_persist_ms=1.0
):
    """
    STA/LTA trigger -> refined onset using envelope hysteresis/persistence.
    Returns list of [onset_idx, offset_idx].
    """
    nsta = int(sta_window * fs)
    nlta = int(lta_window * fs)
    if nsta < 1 or nlta < 1:
        raise ValueError("sta_window/lta_window too small for the given fs.")

    cft = recursive_sta_lta(tr.data, nsta, nlta)

    # Initial triggers from STA/LTA ratio
    raw_triggers = trigger_onset(cft, sta_thresh, lta_thresh)

    refined_triggers = []
    n = len(tr.data)

    for onset, offset in raw_triggers:
        # Focus picking window near trigger onset (prevents late bias)
        i_start = max(0, onset - int(pre_pick_s * fs))
        i_end   = min(n, onset + int(post_pick_s * fs))

        i_ref_start = pick_onset_from_envelope(
            tr.data, fs,
            i_start, i_end,
            k_high=k_high, k_low=k_low,
            min_persist_ms=min_persist_ms
        )

        # Fallback: if envelope picker fails, use original onset
        if i_ref_start is None:
            i_ref_start = onset

        # Offset: search forward from refined onset (or original offset) up to a cap
        off_search_end = min(n, i_ref_start + int(post_offset_s * fs))
        i_ref_end = pick_offset_from_envelope(
            tr.data, fs,
            start_from=i_ref_start,
            search_end=off_search_end,
            k_low=k_low,
            end_persist_ms=end_persist_ms
        )

        # Ensure sane ordering
        if i_ref_end <= i_ref_start:
            i_ref_end = min(n, i_ref_start + int(0.005 * fs))  # at least 5 ms

        refined_triggers.append([int(i_ref_start), int(i_ref_end)])

    return refined_triggers

def _moving_average(x, n):
    """Centered moving average with edge-padding, returns same length as x."""
    n = int(max(1, n))
    if n == 1:
        return x.copy()
    pad = n // 2
    xpad = np.pad(x, (pad, pad), mode="edge")
    w = np.ones(n, dtype=float) / n
    y = np.convolve(xpad, w, mode="valid")
    # If n is even, valid length is len(x)+1; trim to match.
    return y[:len(x)]

def plot_sta_lta_overlay(
    tr,
    fs,
    triggers=None,
    sta_window=0.001,
    lta_window=0.05,
    sta_thresh=2.5,   # on-threshold for ratio (trigger_on)
    lta_thresh=1.5,   # off-threshold for ratio (trigger_off)
    t_shift=0.0,      # start time of the trace in seconds (e.g. tr.stats.starttime offset)
    manual_picks=None,   # list of pick times (in seconds, global time like t_shift)
    manual_pick_color='red',
    manual_pick_style='-',
    manual_pick_label='Manual Pick',
    manual_pick_alpha=0.8,
    manual_picks_samples=None,  # list of sample indices
    fig_size=(12, 8),
    trace_color='black',
    sta_color='#0072B2',
    lta_color='#E69F00',
    ratio_color='#555555',
    ton_color="#009E73",
    toff_color='#CC79A7',
    onset_vline_color='#D55E00',
    offset_vline_color='#56B4E9',
    legend_loc="upper right",
    legend_fontsize=11,
    dpi=300,
    tlim=None,        # (t0, t1) in global seconds (same coordinate as x-axis); None = full record
    use_abs=True,
    save=False,
    save_dir='./',
    save_name='sta_lta_overlay.jpg'
):
    """
    3-panel plot:
      (1) waveform
      (2) STA and LTA overlays (computed as moving averages of |x|)
      (3) STA/LTA ratio with thresholds and trigger markers

    triggers : list of [i_start, i_end] sample indices relative to the start of tr.
    t_shift  : absolute start time of the trace (seconds). The x-axis runs from
               t_shift to t_shift + len(tr)/fs. tlim and trigger lines are all
               expressed in the same global time coordinate.
    """
    x = np.asarray(tr.data, dtype=float)
    n = len(x)
    # Global time axis: starts at t_shift
    t = t_shift + np.arange(n) / fs

    # Use absolute amplitude for STA/LTA (typical in AE / microseismic detection)
    a = np.abs(x) if use_abs else x

    nsta = int(sta_window * fs)
    nlta = int(lta_window * fs)
    if nsta < 1 or nlta < 1:
        raise ValueError("sta_window and lta_window must produce at least 1 sample each.")

    sta = _moving_average(a, nsta)
    lta = _moving_average(a, nlta)

    # Avoid division by zero
    eps = np.finfo(float).eps
    ratio = sta / np.maximum(lta, eps)

    # Time windowing: tlim is in global seconds, convert to local sample indices
    if tlim is not None:
        t0, t1 = tlim
        i0 = max(0, int((t0 - t_shift) * fs))
        i1 = min(n, int((t1 - t_shift) * fs))
    else:
        i0, i1 = 0, n

    # Plot
    fig, axes = plt.subplots(3, 1, figsize=fig_size, sharex=True, dpi=dpi)

    # (1) Waveform
    axes[0].plot(t[i0:i1], x[i0:i1], linewidth=1, color=trace_color)
    axes[0].set_ylabel("Amplitude")
    axes[0].set_title("Waveform")

    # (2) STA & LTA overlays
    axes[1].plot(t[i0:i1], sta[i0:i1], linewidth=1, color=sta_color, label=f"STA ({sta_window*1e3:.1f} ms)")
    axes[1].plot(t[i0:i1], lta[i0:i1], linewidth=1, color=lta_color, label=f"LTA ({lta_window*1e3:.1f} ms)")
    axes[1].set_ylabel("Avg |x|")
    axes[1].set_title("STA & LTA (STA rises sharply at event onset)")
    # axes[1].legend(loc="upper right")

    # (3) Recursive STA/LTA characteristic function — same algorithm used by refined_sta_lta_detection,
    #     so the threshold lines correspond exactly to what triggered (or didn't trigger) the picks.
    cft = recursive_sta_lta(x, nsta, nlta)
    axes[2].plot(t[i0:i1], cft[i0:i1], linewidth=1, color=ratio_color, label="STA/LTA Ratio")
    axes[2].axhline(sta_thresh, linestyle="--", linewidth=1, color=ton_color, label=f"Trigger on = {sta_thresh}")
    axes[2].axhline(lta_thresh, linestyle=":", linewidth=1, color=toff_color, label=f"Trigger off = {lta_thresh}")
    axes[2].set_ylabel("STA/LTA")
    axes[2].set_xlabel("Time (s)")
    axes[2].set_title("STA/LTA Ratio + Thresholds (recursive, causal)")
    # axes[2].legend(loc="upper right")
    
    for ax in axes:
        ax.get_xaxis().get_major_formatter().set_useOffset(False)


    # Manual picks (in seconds, global time)
    if manual_picks is not None:
        for i, tpick in enumerate(manual_picks):
            for ax in axes:
                ax.axvline(
                    tpick,
                    linestyle=manual_pick_style,
                    color=manual_pick_color,
                    linewidth=1.0,
                    alpha=manual_pick_alpha,
                    label=manual_pick_label if i == 0 else None  # avoid duplicate legend entries
                )

    # Manual picks from sample indices
    if manual_picks_samples is not None:
        for i, ipick in enumerate(manual_picks_samples):
            tpick = t_shift + ipick / fs
            for ax in axes:
                ax.axvline(
                    tpick,
                    linestyle=manual_pick_style,
                    color=manual_pick_color,
                    linewidth=1.0,
                    label=manual_pick_label if i == 0 else None
                )
                
    # Trigger markers: sample indices are local to tr, convert to global time via t_shift
    if triggers is not None:
        for j, (onset, offset) in enumerate(triggers):
            ton  = t_shift + onset  / fs
            toff = t_shift + offset / fs
            for ax in axes:
                ax.axvline(ton,  linestyle="-", color=onset_vline_color, linewidth=1.0, 
                           label="STA/LTA Pick" if j==0 else None)
                # ax.axvline(toff, linestyle=":",  color=offset_vline_color, linewidth=1, label="STA/LTA Pick")

    if triggers:
        print(f"Onset times (s):  {[t_shift + onset/fs  for onset, _      in triggers]}")
        print(f"Offset times (s): {[t_shift + offset/fs for _,     offset in triggers]}")

    from collections import OrderedDict

    for ax in axes:
        handles, labels = ax.get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))  # removes duplicates
        if by_label:
            ax.legend(by_label.values(), by_label.keys(), loc=legend_loc, fontsize=legend_fontsize)
    
    plt.tight_layout()
    if save:
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        plt.savefig(os.path.join(save_dir, save_name), dpi=dpi, bbox_inches='tight')
        print(f"Saved figure to {os.path.join(save_dir, save_name)}")
    plt.show()
    return fig, axes


import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.ticker import ScalarFormatter


def plot_selected_ae_events_from_stream(
    stream,
    event_ids,
    event_times,
    pad=0.001,
    receiver_order="descending",
    normalize_per_trace=True,
    scale=0.75,
    line_color="black",
    line_width=0.7,
    fig_size=(10, 6),
    dpi=300,
    ncols=2,
    title=None,
    save=False,
    save_path="selected_ae_events.png",
):
    """
    Plot selected AE events from a continuous ObsPy Stream in a multi-panel layout.

    Parameters
    ----------
    stream : obspy.Stream
        Continuous stream with one Trace per receiver.
    event_ids : list
        Event IDs to plot, e.g. [1, 4, 75, 76].
    event_times : dict
        Dictionary mapping event_id -> (event_start, event_end), in seconds
        relative to the trace start.
    pad : float
        Padding before and after each event window, in seconds.
    receiver_order : {"descending", "ascending"}
        "descending" puts receiver 8 at top and 1 at bottom.
    normalize_per_trace : bool
        If True, normalize each receiver trace independently.
    scale : float
        Vertical scaling of traces around each receiver center.
    """

    if len(stream) == 0:
        raise ValueError("Input stream is empty.")
    if not isinstance(event_times, dict):
        raise TypeError("event_times must be a dict: {event_id: (start, end)}")
    if len(event_ids) == 0:
        raise ValueError("event_ids is empty.")

    fs = float(stream[0].stats.sampling_rate)
    npts = int(stream[0].stats.npts)

    for tr in stream:
        if float(tr.stats.sampling_rate) != fs:
            raise ValueError("All traces must have the same sampling rate.")
        if int(tr.stats.npts) != npts:
            raise ValueError("All traces must have the same number of samples.")

    data = np.vstack([np.asarray(tr.data, dtype=float) for tr in stream])
    n_receivers = data.shape[0]

    if receiver_order.lower() == "descending":
        order = np.arange(n_receivers - 1, -1, -1)
        ytick_labels = [str(i) for i in range(n_receivers, 0, -1)]
    elif receiver_order.lower() == "ascending":
        order = np.arange(n_receivers)
        ytick_labels = [str(i) for i in range(1, n_receivers + 1)]
    else:
        raise ValueError("receiver_order must be 'descending' or 'ascending'.")

    data = data[order, :]

    n_events = len(event_ids)
    ncols = min(max(1, ncols), n_events)
    nrows = int(np.ceil(n_events / ncols))

    fig, axes = plt.subplots(
        nrows, ncols, figsize=fig_size, dpi=dpi, squeeze=False
    )
    axes = axes.ravel()

    # Trace centers exactly at tick positions
    y_positions = np.arange(n_receivers, 0, -1)

    for iax, event_id in enumerate(event_ids):
        ax = axes[iax]

        if event_id not in event_times:
            raise KeyError(f"Event ID {event_id} not found in event_times.")

        t_start, t_end = event_times[event_id]
        if t_end <= t_start:
            raise ValueError(f"Event {event_id}: end time must be greater than start time.")

        t0 = max(0.0, t_start - pad)
        t1 = min(npts / fs, t_end + pad)

        i0 = max(0, int(np.floor(t0 * fs)))
        i1 = min(npts, int(np.ceil(t1 * fs)))

        if i1 <= i0:
            raise ValueError(f"Event {event_id}: invalid extraction window.")

        win = data[:, i0:i1]
        t = np.arange(i0, i1) / fs  # absolute time

        if normalize_per_trace:
            denom = np.max(np.abs(win), axis=1, keepdims=True)
            denom[denom == 0] = 1.0
            win_plot = win / denom
        else:
            denom = np.max(np.abs(win))
            denom = 1.0 if denom == 0 else denom
            win_plot = win / denom

        # Plot each trace centered exactly on its y tick
        for ir in range(n_receivers):
            # y = y_positions[ir] + scale * win_plot[ir]
            # ax.plot(t, y, color=line_color, lw=line_width)
            trace = scale * win_plot[ir]
            y = y_positions[ir] + (trace - np.mean(trace))
            ax.plot(t, y, color=line_color, lw=line_width)

        # Shade event duration
        ax.axvspan(t_start, t_end, color="red", alpha=0.05, zorder=0)

        # Y-axis formatting
        ax.set_ylim(0.5, n_receivers + 0.5)
        ax.set_yticks(y_positions)
        ax.set_yticklabels(ytick_labels, fontsize=9)
        ax.set_ylabel("Receiver", fontsize=10)

        # X-axis formatting: exact time, no +offset
        formatter = ScalarFormatter(useOffset=False)
        formatter.set_scientific(False)
        ax.xaxis.set_major_formatter(formatter)
        ax.ticklabel_format(axis="x", style="plain", useOffset=False)

        ax.set_xlabel("Time (s)", fontsize=10)
        ax.tick_params(axis="x", labelsize=9)
        ax.tick_params(axis="y", labelsize=9)
        ax.tick_params(direction="out", length=4, width=0.8)

        for spine in ax.spines.values():
            spine.set_linewidth(0.8)
            spine.set_color("0.5")

        # Event ID in red oval
        ellipse = Ellipse(
            xy=(0.90, 0.88),
            width=0.12,
            height=0.12,
            transform=ax.transAxes,
            facecolor="none",
            edgecolor="red",
            linewidth=1.5,
            clip_on=False,
        )
        ax.add_patch(ellipse)
        ax.text(
            0.90, 0.88, f"{event_id}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=10,
            color="black",
        )

    for j in range(n_events, len(axes)):
        axes[j].axis("off")

    if title:
        fig.suptitle(title, fontsize=12, y=0.98)
        fig.tight_layout(rect=(0, 0, 1, 0.965))
    else:
        fig.tight_layout()

    if save:
        save_dir = os.path.dirname(save_path)
        if save_dir and not os.path.exists(save_dir):
            os.makedirs(save_dir)
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
        print(f"Saved figure to {save_path}")

    plt.show()
    # return fig, axes


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec

from scipy import signal
from scipy.signal import spectrogram, get_window
from obspy import Stream, Trace, UTCDateTime


# ─────────────────────────────────────────────────────────────────────────────
# Harmonic notch filters
# ─────────────────────────────────────────────────────────────────────────────

def design_harmonic_notches(fs, f0=8000.0, pad_hz=20.0, max_k=None):
    """
    Build a cascaded SOS notch filter that removes f0 and its harmonics,
    each with a +/- pad_hz half-width (i.e., total bandwidth = 2*pad_hz).
    """
    nyq = 0.5 * fs
    if max_k is None:
        max_k = int(np.floor(nyq // f0))
    harmonics = [k * f0 for k in range(1, max_k + 1) if (k * f0 + pad_hz) < nyq]
    if not harmonics:
        raise ValueError("No valid harmonics below Nyquist with the given f0/pad_hz.")

    sos_all = []
    for f in harmonics:
        bw = 2.0 * pad_hz
        Q  = max(1.0, f / bw)
        w0 = f / nyq
        b, a = signal.iirnotch(w0=w0, Q=Q)
        sos = signal.tf2sos(b, a)
        sos_all.append(sos)

    sos_all = np.vstack(sos_all)
    return sos_all, harmonics


def apply_notches(x, fs, f0=8000.0, pad_hz=20.0):
    """Apply harmonic notch filters to a 1-D NumPy array (zero-phase)."""
    sos, harm = design_harmonic_notches(fs, f0, pad_hz)
    y = signal.sosfiltfilt(sos, x)
    return y, harm


def notch_harmonics_stream(st, f0=8000.0, pad_hz=20.0):
    """
    Apply harmonic notch filters to every trace in an ObsPy Stream (in-place
    on a copy).  A 1 % cosine taper is applied first to suppress filtfilt
    edge transients.
    """
    st = st.copy()
    for tr in st:
        fs = float(tr.stats.sampling_rate)
        sos, harm = design_harmonic_notches(fs, f0, pad_hz)
        tr.taper(max_percentage=0.01, type="cosine")
        tr.data = signal.sosfiltfilt(sos, tr.data).astype(tr.data.dtype, copy=False)
        tr.stats.processing = tr.stats.get("processing", []) + [
            f"harmonic_notch: f0={f0} Hz, pad={pad_hz} Hz, harmonics={harm}"
        ]
    return st


# ─────────────────────────────────────────────────────────────────────────────
# Preprocessing
# ─────────────────────────────────────────────────────────────────────────────

def trim_traces(st, tcut_start=40.0, tcut_end=None):
    """
    Remove the first `tcut_start` seconds (and optionally beyond `tcut_end`
    seconds) from each Trace. Returns a new Stream.
    """
    st_trimmed = st.copy()
    for tr in st_trimmed:
        fs = tr.stats.sampling_rate
        ncut_start = int(tcut_start * fs)
        ncut_end   = int(tcut_end  * fs) if tcut_end is not None else None
        tr.data = tr.data[ncut_start:ncut_end]
        tr.stats.starttime += tcut_start
        tr.stats.npts = tr.data.size
    return st_trimmed


def preprocess_stream(
    st,
    freqmin,
    freqmax=None,
    *,
    filter_type="bandpass",
    taper_max_percentage=0.02,
    corners=4,
    zerophase=True,
):
    """
    Return a *copy* of `st` after:
      1. Demean + linear detrend
      2. Cosine taper
      3. Bandpass or highpass filter

    Normalization is intentionally omitted here so that subsequent steps
    (e.g. notch_harmonics_stream) operate on physical-unit amplitudes.
    Normalization for display is handled inside plot_wave_and_spectrogram.

    Always pass the raw (unfiltered) stream to avoid double-processing.
    """
    if filter_type not in {"bandpass", "highpass"}:
        raise ValueError("filter_type must be 'bandpass' or 'highpass'.")

    st = st.copy()
    eps = 1e-6

    for tr in st:
        fs  = float(tr.stats.sampling_rate)
        nyq = 0.5 * fs

        tr.detrend("demean")
        tr.detrend("linear")
        tr.taper(max_percentage=taper_max_percentage, type="cosine")

        if filter_type == "bandpass":
            if freqmax is None:
                raise ValueError("freqmax is required for bandpass.")
            f1 = max(eps, min(float(freqmin), nyq - eps))
            f2 = max(f1 + eps, min(float(freqmax), nyq - eps))
            tr.filter("bandpass", freqmin=f1, freqmax=f2,
                      corners=corners, zerophase=zerophase)
        else:
            f1 = max(eps, min(float(freqmin), nyq - eps))
            tr.filter("highpass", freq=f1, corners=corners, zerophase=zerophase)

    return st


# ─────────────────────────────────────────────────────────────────────────────
# Visualisation
# ─────────────────────────────────────────────────────────────────────────────

# Journal-grade rcParams applied once at import time; override freely.
plt.rcParams.update({
    "font.family":       "serif",
    "font.size":         10,
    "axes.labelsize":    11,
    "axes.titlesize":    11,
    "xtick.labelsize":    9,
    "ytick.labelsize":    9,
    "xtick.direction":  "in",
    "ytick.direction":  "in",
    "xtick.top":         False,
    "ytick.right":       False,
    "axes.linewidth":    0.8,
    "lines.linewidth":   0.6,
})


def plot_wave_and_spectrogram(tr, fmax=None, nperseg=None, overlap=0.90,
                              dyn_range_db=80, norm_percentile=99.9,
                              dpi=300, title=None, tmin=None, tmax=None,
                              ypad=0.1, fig_size=(7.0, 4.5), save_path=None):
    """
    Journal-ready waveform + PSD spectrogram with exactly the same time axis.

    Spectrogram quantity
    --------------------
    Power Spectral Density (PSD) computed via Welch's method (Hann window,
    `scaling='density'`) and expressed in dB:

        PSD_dB = 10 · log₁₀( Sxx / (1 [unit²/Hz]) )

    Values are negative because Sxx [unit²/Hz] << 1 for any physically
    normalised signal — the reference is 1 [unit²/Hz], analogous to dB SPL
    being referenced to 20 μPa.  The dB *differences* (dynamic range) are
    what matter for visualisation.

    Parameters
    ----------
    tr              : obspy.Trace  — preprocessed and notch-filtered
    fmax            : float | None — upper frequency bound [Hz]
    nperseg         : int   | None — STFT window length [samples]; auto → ~2 ms
    overlap         : float        — fractional window overlap in [0, 1)
    dyn_range_db    : float        — colour dynamic range in dB
    norm_percentile : float        — percentile used for waveform normalisation
    dpi             : int          — 300 for print, 150 for screen
    title           : str   | None — panel title (omit for journal figures)
    tmin            : float | None — display window start [s]
    tmax            : float | None — display window end   [s]
    ypad            : float        — fractional y-axis padding (default 0.1)
    save_path       : str   | None — if given, save figure to this path

    Returns
    -------
    fig, (ax_wave, ax_spec)
    """
    x  = np.asarray(tr.data, dtype=float)
    fs = float(tr.stats.sampling_rate)
    N  = x.size

    # ── Time axis & view window ───────────────────────────────────────────
    t_wave = np.arange(N) / fs
    t0 = float(tmin) if tmin is not None else t_wave[0]
    t1 = float(tmax) if tmax is not None else t_wave[-1]
    if t0 >= t1:
        raise ValueError(f"tmin ({t0}) must be less than tmax ({t1}).")

    # ── Normalise for display ─────────────────────────────────────────────
    scale = np.nanpercentile(np.abs(x), norm_percentile)
    if not np.isfinite(scale) or scale == 0:
        scale = np.max(np.abs(x)) or 1.0
    x_norm = x / scale

    # ── Auto y-limits from the visible window ─────────────────────────────
    view_mask = (t_wave >= t0) & (t_wave <= t1)
    x_view    = x_norm[view_mask]
    y_abs     = np.max(np.abs(x_view)) if x_view.size else 1.0
    y_lim     = y_abs * (1.0 + ypad)

    # ── STFT / PSD ────────────────────────────────────────────────────────
    # mode='psd'  → Sxx in [unit²/Hz] → 10·log₁₀ (power dB)
    if nperseg is None:
        nperseg = int(0.002 * fs)
    nperseg  = int(2 ** np.ceil(np.log2(max(256, nperseg))))
    noverlap = int(overlap * nperseg)

    f_ax, t_spec, Sxx = spectrogram(
        x, fs=fs,
        window=get_window("hann", nperseg),
        nperseg=nperseg, noverlap=noverlap, nfft=nperseg,
        detrend=False, scaling="density", mode="psd",
    )
    Sxx_db = 10.0 * np.log10(np.maximum(Sxx, np.finfo(float).eps))

    if fmax is not None:
        fmask  = f_ax <= fmax
        f_ax   = f_ax[fmask]
        Sxx_db = Sxx_db[fmask, :]
    f_khz = f_ax * 1e-3

    tmask    = (t_spec >= t0) & (t_spec <= t1)
    Sxx_view = Sxx_db[:, tmask] if tmask.any() else Sxx_db
    vmax     = np.percentile(Sxx_view, 99.5)
    vmin     = vmax - dyn_range_db

    # ── Figure layout ─────────────────────────────────────────────────────
    # layout="constrained" is designed for GridSpec + colorbars and handles
    # spacing automatically — no tight_layout() call needed or wanted.
    # tight_layout() raises a UserWarning for asymmetric grids (gs[0,1]
    # empty, gs[1,1] occupied by the colorbar); constrained_layout does not.
    fig = plt.figure(figsize=fig_size, dpi=dpi)
    gs  = GridSpec(2, 2, figure=fig,
                   height_ratios=[1, 2.0],
                   width_ratios=[1, 0.03],
                   hspace=0.00, wspace=0.05)

    ax_spec = fig.add_subplot(gs[1, 0])
    ax_wave = fig.add_subplot(gs[0, 0], sharex=ax_spec)
    ax_cbar = fig.add_subplot(gs[1, 1])

    # — Waveform ──────────────────────────────────────────────────────────
    ax_wave.plot(t_wave, x_norm, color="k", lw=0.5, rasterized=True)
    ax_wave.set_ylabel("Amplitude (norm.)")
    ax_wave.set_xlim(t0, t1)
    ax_wave.set_ylim(-y_lim, y_lim)
    # Hide x-axis on waveform panel
    ax_wave.tick_params(axis='x', which='both', bottom=False, labelbottom=False)

    if title:
        ax_wave.set_title(title, fontsize=11)

    # — Spectrogram ───────────────────────────────────────────────────────
    im = ax_spec.pcolormesh(
        t_spec, f_khz, Sxx_db,
        shading="auto", cmap="magma", vmin=vmin, vmax=vmax,
        rasterized=True,
    )
    ax_spec.set_ylim(0, f_khz[-1])
    ax_spec.set_xlabel("Time (s)")
    ax_spec.set_ylabel("Frequency (kHz)")
    ax_spec.tick_params(top=False)
    ax_spec.spines[["top", "right"]].set_visible(True)    # Ensure spectrogram keeps labels
    ax_spec.tick_params(axis='x', which='both', bottom=True, labelbottom=True)

    # — Colorbar ──────────────────────────────────────────────────────────
    cbar = fig.colorbar(im, cax=ax_cbar)
    cbar.set_label("PSD (dB re. 1 unit²/Hz)", fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_ticks(np.linspace(vmin, vmax, 5))
    cbar.set_ticklabels([f"{v:.0f}" for v in np.linspace(vmin, vmax, 5)])
    
    # # Add annotations and highlights
    # event_times = [58, 98, 137]  # example

    # for t_event in event_times:
    #     ax_wave.axvline(t_event, color='cyan', linestyle='--', lw=0.5, alpha=0.7)
    #     ax_spec.axvline(t_event, color='cyan', linestyle='--', lw=0.5, alpha=0.7)
        
    # ax_wave.annotate(
    #     "Main AE event",
    #     xy=(98, 0.6*y_lim),
    #     xytext=(108, 1.0*y_lim),
    #     arrowprops=dict(arrowstyle="->", lw=1),
    #     fontsize=9
    # )
    
    # ax_spec.axhspan(2, 48, color='white', alpha=0.1)
    
    # ax_spec.text(
    #     100, 1.0,
    #     "Low-frequency noise",
    #     fontsize=8,
    #     color="white"
    # )

    # ax_wave.text(0.01, 0.9, "(a)", transform=ax_wave.transAxes, fontweight="bold")
    # ax_spec.text(0.01, 0.9, "(b)", transform=ax_spec.transAxes, fontweight="bold", color="white")
    
    
    fig.tight_layout()
    if save_path is not None:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
        print(f"Saved → {save_path}")

    plt.show()
    # return fig, (ax_wave, ax_spec)

