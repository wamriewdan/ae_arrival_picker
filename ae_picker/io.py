"""
ae_picker.io
============
File readers for acoustic emission (AE) waveform data.

Both readers return the same ``(meta, channels)`` structure, making all
downstream pickers and plotters format-agnostic.  To support a different
dataset, write a new reader that returns this same structure::

    meta     : dict          – metadata (keys vary by format)
    channels : list of dict  – each dict has
                  'time' : np.ndarray  [seconds]
                  'amp'  : np.ndarray  [volts]

Public API
----------
read_unfiltered(filepath)
read_filtered(filepath)
"""

import re
import numpy as np
import pandas as pd
from pathlib import Path


def read_unfiltered(filepath):
    """
    Read a 4-channel unfiltered AE waveform file.

    Expected file format
    --------------------
    Line 1 : project / channel ID  (semicolon-separated)
    Line 2 : equipment string
    Line 3 : ``Flight dist [m]: <val>  Offset [s]: <val>``
    Line 4 : ``Picking time P [s]: <val>  Picking time S [s]: <val>``
    Line 5 : ``Measure depth [m] <val>``
    Line 6 : column headers — ``Time [s]  Amplitude [V]`` repeated × n_channels
    Lines 7+: tab-separated numeric data (one row per time sample)

    The file may have trailing tab characters on each row; these are handled
    transparently via ``index_col=False`` and ``dropna``.

    Parameters
    ----------
    filepath : str or Path

    Returns
    -------
    meta : dict
        Keys: ``id``, ``equip``, ``flight_dist_m``, ``offset_s``,
        ``pick_P_s``, ``pick_S_s``, ``depth_m``.
        ``pick_P_s`` / ``pick_S_s`` are ``np.nan`` when not present in the
        header (some files store only a P pick or neither pick).
    channels : list of dict
        One entry per channel pair in the file, each with keys
        ``'time'`` (np.ndarray, seconds) and ``'amp'`` (np.ndarray, volts).
    """
    filepath = Path(filepath)
    meta = {}

    with open(filepath, encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()

    meta["id"]    = lines[0].strip()
    meta["equip"] = lines[1].strip()

    m = re.search(
        r"Flight dist \[m\]:\s*([\d.eE+\-]+)\s+Offset \[s\]:\s*([\d.eE+\-]+)",
        lines[2],
    )
    if m:
        meta["flight_dist_m"] = float(m.group(1))
        meta["offset_s"]      = float(m.group(2))

    m = re.search(
        r"Picking time P \[s\]:\s*([\d.eE+\-]+)\s+Picking time S \[s\]:\s*([\d.eE+\-]+)",
        lines[3],
    )
    if m:
        meta["pick_P_s"] = float(m.group(1))
        meta["pick_S_s"] = float(m.group(2))
    else:
        meta["pick_P_s"] = np.nan
        meta["pick_S_s"] = np.nan

    m = re.search(r"Measure depth \[m\]\s+([\d.eE+\-]+)", lines[4])
    if m:
        meta["depth_m"] = float(m.group(1))

    # index_col=False prevents trailing-tab phantom columns from being
    # absorbed as the row index, which would shift every column by one.
    df = pd.read_csv(filepath, sep="\t", skiprows=5, header=0, index_col=False)
    df = df.dropna(axis=1, how="all")

    channels = []
    for i in range(0, df.shape[1] - 1, 2):
        channels.append({
            "time": pd.to_numeric(df.iloc[:, i],     errors="coerce").values,
            "amp":  pd.to_numeric(df.iloc[:, i + 1], errors="coerce").values,
        })

    return meta, channels


def read_filtered(filepath):
    """
    Read a 2-channel filtered AE waveform file centred at t = 0.

    Expected file format
    --------------------
    Line 1 : column headers —
              ``Time [s]  Amplitude [V]  Time [s].2  Amplitude [V].2``
    Lines 2+: tab-separated numeric data; time axis spans approximately
              ±250 µs around the nominal arrival time.

    Column pairing is done by name prefix (``Time`` → time axis,
    ``Amplitude`` → amplitude axis), so additional channel pairs are
    picked up automatically.

    Parameters
    ----------
    filepath : str or Path

    Returns
    -------
    meta : dict
        Keys: ``filename`` (Path stem of the file).
    channels : list of dict
        One entry per channel pair, each with keys
        ``'time'`` (np.ndarray, seconds) and ``'amp'`` (np.ndarray, volts).
    """
    filepath = Path(filepath)
    meta     = {"filename": filepath.stem}

    df   = pd.read_csv(filepath, sep="\t", header=0)
    cols = df.columns.tolist()

    time_cols = [c for c in cols if c.lower().startswith("time")]
    amp_cols  = [c for c in cols if c.lower().startswith("amplitude")]

    channels = []
    for t_col, a_col in zip(time_cols, amp_cols):
        channels.append({
            "time": df[t_col].values.astype(float),
            "amp":  df[a_col].values.astype(float),
        })

    return meta, channels


def read_segy(filepath, component=None):
    """
    Read a SEG-Y file and return the standard ``(meta, channels)`` structure.

    Each trace in the file becomes one entry in *channels*.  The time axis
    is reconstructed from the trace header sample interval.

    Parameters
    ----------
    filepath : str or Path
    component : str or None
        Optional label stored in ``meta['component']`` (e.g. ``'Z'``).
        Has no effect on the data returned.

    Returns
    -------
    meta : dict
        Keys: ``filename``, ``format`` (``'segy'``), ``n_traces``,
        ``dt_s`` (sample interval in seconds), ``n_samples``,
        and ``component`` (when provided).
    channels : list of dict
        One entry per trace, each with ``'time'`` (np.ndarray, seconds)
        and ``'amp'`` (np.ndarray, counts or field units).

    Raises
    ------
    ImportError
        When ObsPy is not installed.
    """
    try:
        from obspy import read as _obspy_read
    except ImportError:
        raise ImportError(
            "ObsPy is required to read SEG-Y files. "
            "Install it with: pip install ae-picker[obspy]"
        )

    filepath = Path(filepath)
    st = _obspy_read(str(filepath), format="SEGY")

    dt_s      = st[0].stats.delta
    n_samples = st[0].stats.npts
    t_axis    = np.arange(n_samples, dtype=float) * dt_s

    channels = []
    for tr in st:
        channels.append({
            "time": t_axis.copy(),
            "amp":  tr.data.astype(float),
        })

    meta = {
        "filename": filepath.stem,
        "format":   "segy",
        "n_traces": len(st),
        "dt_s":     dt_s,
        "n_samples": n_samples,
    }
    if component is not None:
        meta["component"] = component

    return meta, channels


def read_mseed(filepath):
    """
    Read a miniSEED file and return the standard ``(meta, channels)`` structure.

    Each trace in the stream becomes one entry in *channels*.

    Parameters
    ----------
    filepath : str or Path

    Returns
    -------
    meta : dict
        Keys: ``filename``, ``format`` (``'mseed'``), ``n_traces``,
        and per-trace ``network``, ``station``, ``channel`` from the
        first trace's stats.
    channels : list of dict
        One entry per trace, each with ``'time'`` (np.ndarray, seconds
        from trace start) and ``'amp'`` (np.ndarray).

    Raises
    ------
    ImportError
        When ObsPy is not installed.
    """
    try:
        from obspy import read as _obspy_read
    except ImportError:
        raise ImportError(
            "ObsPy is required to read miniSEED files. "
            "Install it with: pip install ae-picker[obspy]"
        )

    filepath = Path(filepath)
    st = _obspy_read(str(filepath), format="MSEED")

    channels = []
    for tr in st:
        dt_s = tr.stats.delta
        n    = tr.stats.npts
        channels.append({
            "time": np.arange(n, dtype=float) * dt_s,
            "amp":  tr.data.astype(float),
        })

    first = st[0].stats
    meta = {
        "filename":  filepath.stem,
        "format":    "mseed",
        "n_traces":  len(st),
        "network":   first.network,
        "station":   first.station,
        "channel":   first.channel,
    }

    return meta, channels
