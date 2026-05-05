"""
ae_picker.plot
==============
Plotting helpers for acoustic emission (AE) waveform visualisation.

All functions return a ``matplotlib.figure.Figure`` object so callers can
further customise the plot or embed it in a GUI / notebook.

Public API
----------
plot_channels(meta, channels, picks, ...)
    Multi-channel waveform plot with named pick lines.
plot_comparison(df, ref_col, pick_cols, ...)
    Scatter + residual chart comparing computed picks to a reference.
plot_onset_zoom(files, reader_fn, ...)
    Side-by-side channel zoom browser for rapid visual QC.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path


# ── colour palette for picks (cycles if more than 6 pickers) ─────────────────
_PICK_COLOURS = ["crimson", "darkorange", "royalblue", "seagreen", "purple", "saddlebrown"]
_PICK_STYLES  = ["-", "--", "-.", ":", "-", "--"]


def plot_channels(meta, channels, picks, signal_type="unfiltered",
                  zoom_us=None, savepath=None):
    """
    Plot all channels of one AE record with named pick overlays.

    Parameters
    ----------
    meta : dict
        Metadata dict returned by a reader (``read_unfiltered`` /
        ``read_filtered``).  Used for the figure title and for drawing
        header P/S picks on unfiltered signals.
    channels : list of dict
        Channel list returned by a reader.  Each dict has ``'time'`` and
        ``'amp'`` as numpy arrays.
    picks : dict[str, list[int | None]]
        Named picks to overlay.  Keys are labels (e.g. ``'AIC'``,
        ``'Energy'``); values are lists of sample indices, one per channel.
        ``None`` entries are silently skipped.

        Example::

            picks = {
                'AIC':    [12, 14, 10, 15],   # one index per channel
                'Energy': [15, None, 18, 20],
            }

    signal_type : {'unfiltered', 'filtered'}
        Controls whether header P/S reference picks from ``meta`` are drawn.
    zoom_us : float or None
        If set, the x-axis of every subplot is restricted to
        ``[t[0], t[0] + zoom_us]`` microseconds.  Useful for inspecting
        the onset region.  Default: ``None`` (full record shown).
    savepath : str, Path, or None
        If provided, the figure is saved to this path (PNG, PDF, etc.) at
        150 dpi before being returned.

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    n_ch = len(channels)
    fig  = plt.figure(figsize=(13, 2.8 * n_ch))
    gs   = gridspec.GridSpec(n_ch, 1, hspace=0.55)

    title = meta.get("id", meta.get("filename", "Signal"))
    fig.suptitle(f"{title}  [{signal_type}]", fontsize=11, fontweight="bold")

    pick_names   = list(picks.keys())
    pick_colours = [_PICK_COLOURS[i % len(_PICK_COLOURS)] for i in range(len(pick_names))]
    pick_styles  = [_PICK_STYLES[i % len(_PICK_STYLES)]   for i in range(len(pick_names))]

    for ch_idx, ch in enumerate(channels):
        ax = fig.add_subplot(gs[ch_idx])
        t  = ch["time"]
        a  = ch["amp"]

        x_us = t * 1e6
        ax.plot(x_us, a, color="steelblue", lw=0.7, label="Waveform")

        # ── computed picks ────────────────────────────────────────────────
        for name, colour, style in zip(pick_names, pick_colours, pick_styles):
            idx_list = picks[name]
            if ch_idx >= len(idx_list):
                continue
            idx = idx_list[ch_idx]
            if idx is None:
                continue
            t_pick = t[idx] * 1e6
            ax.axvline(t_pick, color=colour, lw=1.6, ls=style,
                       label=f"{name}  {t_pick:.1f} µs")

        # ── header picks (unfiltered only) ────────────────────────────────
        if signal_type == "unfiltered":
            p_P = meta.get("pick_P_s", np.nan)
            p_S = meta.get("pick_S_s", np.nan)
            if not np.isnan(p_P):
                ax.axvline(p_P * 1e6, color="limegreen", lw=1.2, ls=":",
                           label=f"Hdr P  {p_P*1e6:.1f} µs")
            if not np.isnan(p_S):
                ax.axvline(p_S * 1e6, color="mediumpurple", lw=1.2, ls=":",
                           label=f"Hdr S  {p_S*1e6:.1f} µs")

        # ── x-axis zoom ───────────────────────────────────────────────────
        if zoom_us is not None:
            ax.set_xlim(x_us[0], x_us[0] + zoom_us)

        ax.set_title(f"Channel {ch_idx + 1}", fontsize=9)
        ax.set_xlabel("Time [µs]", fontsize=8)
        ax.set_ylabel("Amplitude [V]", fontsize=8)
        ax.tick_params(labelsize=8)
        ax.grid(True, ls=":", alpha=0.5)
        ax.legend(fontsize=7, loc="upper right", ncol=2)

    if savepath is not None:
        fig.savefig(savepath, dpi=150, bbox_inches="tight")

    return fig


def plot_comparison(df, ref_col="header_P_us", pick_cols=("aic_pick_us", "energy_pick_us"),
                    savepath=None):
    """
    Compare computed picks against a reference (e.g. header picks).

    Produces two side-by-side panels:

    * **Left** — scatter plot: computed pick vs reference pick, with a 1:1
      diagonal.  Ideal picks fall on the diagonal.
    * **Right** — residual bar chart: (computed pick − reference) per
      channel entry, one bar group per picker.

    Parameters
    ----------
    df : pandas.DataFrame
        Output of ``batch.run()`` or ``batch.pick_file()``.
    ref_col : str
        Column in *df* to use as the reference pick (µs).
        Default: ``'header_P_us'``.
    pick_cols : tuple of str
        Columns in *df* to compare against the reference.
        Default: ``('aic_pick_us', 'energy_pick_us')``.
    savepath : str, Path, or None
        If provided, the figure is saved to this path.

    Returns
    -------
    fig : matplotlib.figure.Figure
        ``None`` if *df* has no rows with a valid reference pick.
    """
    sub = df.dropna(subset=[ref_col]).copy()
    if sub.empty:
        print(f"[plot_comparison] No rows with a valid '{ref_col}' — nothing to plot.")
        return None

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # ── left: scatter ─────────────────────────────────────────────────────
    ax = axes[0]
    colours = _PICK_COLOURS[: len(pick_cols)]
    markers = ["o", "^", "s", "D"]
    all_vals = [sub[ref_col].values]

    for col, colour, marker in zip(pick_cols, colours, markers):
        if col not in sub.columns:
            continue
        vals = sub[col].values
        all_vals.append(vals[~np.isnan(vals)])
        label = col.replace("_pick_us", "").replace("_", " ").title()
        ax.scatter(sub[ref_col], vals, c=colour, marker=marker,
                   s=55, zorder=3, label=label)

    if all_vals:
        flat = np.concatenate([v for v in all_vals if len(v)])
        lim  = [np.nanmin(flat) * 0.9, np.nanmax(flat) * 1.1]
        ax.plot(lim, lim, "k--", lw=0.8, label="1 : 1")
        ax.set_xlim(lim)
        ax.set_ylim(lim)

    ref_label = ref_col.replace("_us", "").replace("_", " ").title()
    ax.set_xlabel(f"{ref_label} [µs]")
    ax.set_ylabel("Computed pick [µs]")
    ax.set_title("Picks vs reference")
    ax.legend(fontsize=8)
    ax.grid(True, ls=":", alpha=0.5)

    # ── right: residuals ──────────────────────────────────────────────────
    ax  = axes[1]
    x   = np.arange(len(sub))
    w   = 0.8 / max(len(pick_cols), 1)
    offsets = np.linspace(-(len(pick_cols) - 1) / 2, (len(pick_cols) - 1) / 2,
                          len(pick_cols)) * w

    for col, colour, offset in zip(pick_cols, colours, offsets):
        if col not in sub.columns:
            continue
        residuals = (sub[col] - sub[ref_col]).values
        label = col.replace("_pick_us", "").replace("_", " ").title()
        ax.bar(x + offset, residuals, w * 0.9, color=colour,
               alpha=0.8, label=label)

    ax.axhline(0, color="k", lw=0.8)
    tick_labels = [
        f"{r['file'][:10]}\nCh{r['channel']}"
        for _, r in sub.iterrows()
    ]
    ax.set_xticks(x)
    ax.set_xticklabels(tick_labels, fontsize=6, rotation=45, ha="right")
    ax.set_ylabel(f"Pick − {ref_label} [µs]")
    ax.set_title("Residuals vs reference")
    ax.legend(fontsize=8)
    ax.grid(True, ls=":", alpha=0.5)

    plt.tight_layout()

    if savepath is not None:
        fig.savefig(savepath, dpi=150, bbox_inches="tight")

    return fig


def plot_onset_zoom(files, reader_fn, zoom_us=200,
                    noise_window_s=25e-6, energy_threshold=3.0,
                    aic_search_end=10, savepath_dir=None):
    """
    Side-by-side channel zoom browser — one figure per file.

    For each file, plots the first ``zoom_us`` microseconds of every channel
    with AIC and energy-onset picks overlaid.  Designed for rapid visual QC
    of a batch of unfiltered records.

    Parameters
    ----------
    files : list of Path
        Files to browse.
    reader_fn : callable
        Reader function (e.g. ``read_unfiltered`` or ``read_filtered``) that
        accepts a file path and returns ``(meta, channels)``.
    zoom_us : float
        X-axis width in microseconds.  Default: 200 µs.
    noise_window_s : float
        Noise window for ``energy_onset_picker``.  Default: 25 µs.
    energy_threshold : float
        Energy trigger ratio.  Default: 3.0.
    aic_search_end : int
        ``search_end`` passed to ``aic_picker``.  Default: 10 samples.
    savepath_dir : Path or None
        If provided, each figure is saved as
        ``<savepath_dir>/<stem>_zoom.png``.

    Returns
    -------
    figures : list of matplotlib.figure.Figure
    """
    from .pickers import aic_picker, energy_onset_picker

    figures = []
    for fp in files:
        meta, channels = reader_fn(fp)
        n_ch = len(channels)

        fig, axes = plt.subplots(1, n_ch, figsize=(4 * n_ch, 3), sharey=False)
        if n_ch == 1:
            axes = [axes]

        stem  = meta.get("id", meta.get("filename", Path(fp).stem))
        fig.suptitle(f"{stem}  (first {zoom_us} µs)", fontsize=9, fontweight="bold")

        for i, ax in enumerate(axes):
            t, amp = channels[i]["time"], channels[i]["amp"]
            mask   = (t * 1e6) <= (t[0] * 1e6 + zoom_us)

            ax.plot(t[mask] * 1e6, amp[mask], color="steelblue", lw=1.0)

            p_aic, _ = aic_picker(amp, search_start=1,
                                  search_end=aic_search_end)
            p_en     = energy_onset_picker(amp, t,
                                           noise_window_s=noise_window_s,
                                           threshold=energy_threshold)

            ax.axvline(t[p_aic] * 1e6, color="crimson", lw=1.5, ls="-",
                       label=f"AIC {t[p_aic]*1e6:.1f} µs")
            if p_en is not None:
                ax.axvline(t[p_en] * 1e6, color="darkorange", lw=1.5, ls="--",
                           label=f"Energy {t[p_en]*1e6:.1f} µs")

            pP = meta.get("pick_P_s", np.nan)
            if not np.isnan(pP):
                ax.axvline(pP * 1e6, color="limegreen", lw=1.2, ls=":",
                           label=f"Hdr P {pP*1e6:.1f}")

            ax.set_title(f"Ch {i + 1}", fontsize=8)
            ax.set_xlabel("Time [µs]", fontsize=7)
            if i == 0:
                ax.set_ylabel("Amplitude [V]", fontsize=7)
            ax.tick_params(labelsize=7)
            ax.grid(True, ls=":", alpha=0.5)
            ax.legend(fontsize=6, loc="upper right")

        plt.tight_layout()

        if savepath_dir is not None:
            out = Path(savepath_dir) / (Path(fp).stem + "_zoom.png")
            fig.savefig(out, dpi=150, bbox_inches="tight")

        figures.append(fig)

    return figures
