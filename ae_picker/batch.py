"""
ae_picker.batch
===============
Batch processing helpers for acoustic emission (AE) first-arrival picking.

Public API
----------
pick_file(filepath, signal_type, ...)
    Pick all channels in a single file using any combination of pickers.
run(data_dir, ...)
    Scan a directory tree and batch-process all waveform files found.
"""

import numpy as np
import pandas as pd
from pathlib import Path

from .io      import read_unfiltered, read_filtered
from .pickers import (aic_picker, prepend_noise_aic_picker,
                      envelope_onset_picker,
                      stalta_picker, refined_stalta_picker)
from .plot    import plot_channels


# ── internal helpers ──────────────────────────────────────────────────────────

def _nan_time(channels, idx):
    """Return time value (s) for index *idx*, or NaN if idx is None."""
    if idx is None:
        return np.nan
    return float(channels[0]["time"][0])   # placeholder — overridden below


def _t(ch, idx):
    """Time value at sample *idx* in channel dict *ch*, or NaN."""
    if idx is None:
        return np.nan
    return float(ch["time"][idx])


# ── public API ────────────────────────────────────────────────────────────────

def pick_file(filepath, signal_type,
              pickers=("aic", "refined_stalta"),
              aic_search_frac_unfiltered=(0.0005, 0.005),
              aic_search_frac_filtered=(0.30, 0.70),
              sta_s=5e-6, lta_s=50e-6, stalta_threshold=3.0,
              refined_sta_thresh=2.5, refined_lta_thresh=1.5,
              refined_search_start_frac=0.0,
              envelope_k_high=6.0, envelope_k_low=2.5,
              prepend_duration_s=50e-6,
              plot=True, save_plots=True, show_plots=True, output_dir=None):
    """
    Pick all channels in a single AE waveform file.

    Parameters
    ----------
    filepath : str or Path
    signal_type : {'unfiltered', 'filtered'}
        Determines which reader and default search windows are used.
    pickers : tuple of str
        Any subset of
        ``('aic', 'aic_prepend', 'envelope', 'stalta', 'refined_stalta')``.
        ``'aic'`` is the original baseline AIC picker.
        ``'aic_prepend'`` prepends synthetic Gaussian noise before AIC.
        ``'envelope'`` uses Hilbert-envelope hysteresis with MAD noise stats.
        ``'refined_stalta'`` runs recursive STA/LTA then refines the onset
        with the envelope picker.
        All listed pickers are run; unlisted pickers produce NaN columns.
    aic_search_frac_unfiltered : tuple (start_frac, end_frac)
        AIC search window expressed as fractions of the record length, used
        for **unfiltered** files.  Default ``(0.0005, 0.005)`` ≈ 1–50 µs
        at 200 kHz — keeps the search in the onset region only.
    aic_search_frac_filtered : tuple (start_frac, end_frac)
        AIC search window for **filtered** files.  Default ``(0.30, 0.70)``
        searches the central 40 % of the ±250 µs record.
    sta_s, lta_s : float
        STA and LTA window lengths for ``stalta_picker`` (seconds).
    stalta_threshold : float
        STA/LTA trigger ratio for the classical ``stalta_picker``.
    refined_sta_thresh, refined_lta_thresh : float
        Trigger-on / trigger-off thresholds for ``refined_stalta_picker``.
    envelope_k_high, envelope_k_low : float
        MAD multipliers for ``envelope_onset_picker`` and
        ``refined_stalta_picker``'s refinement stage.
    prepend_duration_s : float
        Duration of the synthetic noise window prepended before running
        ``prepend_noise_aic_picker`` (seconds).  Default: 50 µs.  Only used
        when ``'aic_prepend'`` is in *pickers*.
    plot : bool
        If ``True``, a channel plot is produced for each file.
    save_plots : bool
        If ``True`` and *plot* is ``True``, the figure is saved to
        *output_dir* as ``<stem>_picks.png``.
    show_plots : bool
        If ``True``, display the generated figure interactively.  Set this
        to ``False`` for command-line or batch workflows that only need
        files written to disk.
    output_dir : Path or None
        Directory for saved figures.  Required when *save_plots* is ``True``.

    Returns
    -------
    results : list of dict
        One dict per channel with keys::

            file, signal_type, channel,
            aic_pick_s,            aic_pick_us,
            aic_prepend_pick_s,    aic_prepend_pick_us,
            stalta_pick_s,         stalta_pick_us,
            refined_stalta_pick_s, refined_stalta_pick_us,
            header_P_s,            header_S_s
    """
    filepath = Path(filepath)

    if signal_type == "unfiltered":
        meta, channels = read_unfiltered(filepath)
        aic_frac = aic_search_frac_unfiltered
    else:
        meta, channels = read_filtered(filepath)
        aic_frac = aic_search_frac_filtered

    N = len(channels[0]["time"])
    ss = max(1, int(aic_frac[0] * N))
    se = max(ss + 2, int(aic_frac[1] * N))

    aic_idxs              = []
    aic_prepend_idxs      = []
    envelope_idxs         = []
    stalta_idxs           = []
    refined_stalta_idxs   = []
    results               = []

    for ch in channels:
        t, amp = ch["time"], ch["amp"]

        # ── AIC (baseline) ────────────────────────────────────────────────
        if "aic" in pickers:
            p_aic, _ = aic_picker(amp, search_start=ss, search_end=se)
        else:
            p_aic = None
        aic_idxs.append(p_aic)

        # ── AIC with prepended noise ──────────────────────────────────────
        if "aic_prepend" in pickers:
            p_pre, _, _ = prepend_noise_aic_picker(
                amp, t, prepend_duration_s=prepend_duration_s)
        else:
            p_pre = None
        aic_prepend_idxs.append(p_pre)

        # ── Envelope onset (Hilbert + MAD) ────────────────────────────────
        if "envelope" in pickers:
            p_env = envelope_onset_picker(
                amp, t,
                search_start_s=t[ss],
                search_end_s=t[se],
                k_high=envelope_k_high,
                k_low=envelope_k_low)
        else:
            p_env = None
        envelope_idxs.append(p_env)

        # ── Classical STA/LTA ─────────────────────────────────────────────
        if "stalta" in pickers:
            p_sl, _ = stalta_picker(amp, t, sta_s=sta_s, lta_s=lta_s,
                                    threshold=stalta_threshold,
                                    search_start=ss)
        else:
            p_sl = None
        stalta_idxs.append(p_sl)

        # ── Refined STA/LTA ───────────────────────────────────────────────
        if "refined_stalta" in pickers:
            p_rs, _, _ = refined_stalta_picker(
                amp, t, sta_s=sta_s, lta_s=lta_s,
                sta_thresh=refined_sta_thresh,
                lta_thresh=refined_lta_thresh,
                k_high=envelope_k_high,
                k_low=envelope_k_low,
                search_start=int(refined_search_start_frac * N))
        else:
            p_rs = None
        refined_stalta_idxs.append(p_rs)

        t_aic = _t(ch, p_aic)
        t_pre = _t(ch, p_pre)
        t_env = _t(ch, p_env)
        t_sl  = _t(ch, p_sl)
        t_rs  = _t(ch, p_rs)

        def _us(val):
            return val * 1e6 if not np.isnan(val) else np.nan

        results.append({
            "file":                      filepath.name,
            "signal_type":               signal_type,
            "channel":                   len(results) + 1,
            "aic_pick_s":                t_aic,
            "aic_pick_us":               _us(t_aic),
            "aic_prepend_pick_s":        t_pre,
            "aic_prepend_pick_us":       _us(t_pre),
            "envelope_pick_s":           t_env,
            "envelope_pick_us":          _us(t_env),
            "stalta_pick_s":             t_sl,
            "stalta_pick_us":            _us(t_sl),
            "refined_stalta_pick_s":     t_rs,
            "refined_stalta_pick_us":    _us(t_rs),
            "header_P_s":                meta.get("pick_P_s", np.nan),
            "header_S_s":                meta.get("pick_S_s", np.nan),
        })

    if plot:
        picks_for_plot = {}
        if "aic"             in pickers:
            picks_for_plot["AIC"]             = aic_idxs
        if "aic_prepend"     in pickers:
            picks_for_plot["AIC+prepend"]     = aic_prepend_idxs
        if "envelope"        in pickers:
            picks_for_plot["Envelope"]        = envelope_idxs
        if "stalta"          in pickers:
            picks_for_plot["STA/LTA"]         = stalta_idxs
        if "refined_stalta"  in pickers:
            picks_for_plot["Refined STA/LTA"] = refined_stalta_idxs

        savepath = None
        if save_plots and output_dir is not None:
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            savepath = Path(output_dir) / (filepath.stem + "_picks.png")

        fig = plot_channels(meta, channels, picks_for_plot,
                            signal_type=signal_type, savepath=savepath)
        import matplotlib.pyplot as plt
        if show_plots:
            plt.show()
        plt.close(fig)

    return results


def run(data_dir, signal_dirs=("unfiltered_signals", "filtered_signals"),
        output_dir=None, pickers=("aic", "refined_stalta"),
        prepend_duration_s=50e-6,
        plot=True, save_plots=True, show_plots=True):
    """
    Batch-process all ``.txt`` waveform files found under *data_dir*.

    The function scans each directory listed in *signal_dirs* (relative to
    *data_dir*), infers the signal type from the directory name, and calls
    ``pick_file`` on every file found.

    Parameters
    ----------
    data_dir : str or Path
        Root directory containing the signal sub-directories.
    signal_dirs : tuple of str
        Sub-directory names to process.  Each name must contain either
        ``'unfiltered'`` or ``'filtered'`` to determine the reader and
        default AIC window.
    output_dir : str, Path, or None
        Where to write ``picks_summary.csv`` and PNG figures.
        Defaults to ``data_dir / 'picks_output'``.
    pickers : tuple of str
        Any subset of
        ``('aic', 'aic_prepend', 'envelope', 'stalta', 'refined_stalta')``.
    prepend_duration_s : float
        Noise prepend duration passed to ``prepend_noise_aic_picker``.
    plot : bool
        Produce a per-file waveform plot.
    save_plots : bool
        Save each figure to *output_dir*.
    show_plots : bool
        Display plots interactively while processing.

    Returns
    -------
    df : pandas.DataFrame
        Combined results for all files and channels.  Also saved as
        ``picks_summary.csv`` inside *output_dir*.
    """
    data_dir   = Path(data_dir)
    output_dir = Path(output_dir) if output_dir else data_dir / "picks_output"
    output_dir.mkdir(parents=True, exist_ok=True)

    _type_map = {"unfiltered": "unfiltered", "filtered": "filtered"}

    all_results = []

    for sub in signal_dirs:
        folder = data_dir / sub
        if not folder.exists():
            print(f"[ae_picker.batch] Skipping '{folder}' — directory not found.")
            continue

        signal_type = None
        for key, val in _type_map.items():
            if key in sub.lower():
                signal_type = val
                break
        if signal_type is None:
            print(f"[ae_picker.batch] Cannot infer signal type from '{sub}', skipping.")
            continue

        files = sorted(folder.glob("*.txt"))
        print(f"\n=== {sub} [{signal_type}] ({len(files)} files) ===")

        for fp in files:
            print(f"  {fp.name}")
            rows = pick_file(fp, signal_type,
                             pickers=pickers,
                             prepend_duration_s=prepend_duration_s,
                             plot=plot, save_plots=save_plots,
                             show_plots=show_plots,
                             output_dir=output_dir)
            all_results.extend(rows)

    df = pd.DataFrame(all_results)
    csv_path = output_dir / "picks_summary.csv"
    df.to_csv(csv_path, index=False)
    print(f"\nSummary saved: {csv_path}")
    return df
