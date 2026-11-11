"""
ae_picker.presets
=================
Frequency-adaptive parameter presets for first-arrival picking.

Removes the manual scaling burden when moving from AE data (~8 MHz) to
seismic data (1–20 kHz).  One call to ``suggest_parameters`` returns a
complete set of STA/LTA and onset-picker parameters appropriate for the
recording's sampling rate.

Public API
----------
suggest_parameters(fs, f_dominant=None, profile=None) -> dict
"""

__all__ = ["suggest_parameters"]

# Per-profile defaults.  All time values in seconds, persist in ms.
_PROFILE_DEFAULTS = {
    "ae": {
        "sta_s":                5e-6,
        "lta_s":                100e-6,
        "min_persist_ms":       0.1,
        "noise_lookback_s":     40e-6,
        "pre_pick_s":           50e-6,
        "post_pick_s":          200e-6,
        "search_window_hint_s": 50e-6,
    },
    "ps_log": {
        "sta_s":                1e-3,
        "lta_s":                10e-3,
        "min_persist_ms":       1.0,
        "noise_lookback_s":     2e-3,
        "pre_pick_s":           2e-3,
        "post_pick_s":          10e-3,
        "search_window_hint_s": 5e-3,
    },
    "refraction": {
        "sta_s":                10e-3,
        "lta_s":                100e-3,
        "min_persist_ms":       5.0,
        "noise_lookback_s":     20e-3,
        "pre_pick_s":           20e-3,
        "post_pick_s":          100e-3,
        "search_window_hint_s": 50e-3,
    },
}


def suggest_parameters(fs, f_dominant=None, profile=None):
    """
    Return recommended picker parameters scaled to the recording's sampling rate.

    Parameters
    ----------
    fs : float
        Sampling rate in Hz.  Examples: ``8e6`` (AE), ``4000`` (PS log),
        ``500`` (refraction survey).
    f_dominant : float or None
        Expected dominant signal frequency in Hz.  When provided,
        ``sta_s = 2 / f_dominant`` and ``lta_s = 20 / f_dominant``
        override the profile defaults (all other parameters remain at
        their profile values).
    profile : {'ae', 'ps_log', 'refraction', None}
        Explicit profile selection.  When ``None`` the profile is inferred
        automatically from ``fs``:

        * ``fs > 100_000`` Hz → ``'ae'``
        * ``fs > 1_000``   Hz → ``'ps_log'``
        * otherwise        → ``'refraction'``

    Returns
    -------
    params : dict
        Keys:

        ``sta_s``
            STA window length (seconds).
        ``lta_s``
            LTA window length (seconds).
        ``min_persist_ms``
            Minimum persistence for onset confirmation (milliseconds).
        ``noise_lookback_s``
            Pre-search noise estimation window (seconds).
        ``pre_pick_s``
            Look-before-trigger for envelope refinement (seconds).
        ``post_pick_s``
            Look-after-trigger for envelope refinement (seconds).
        ``search_window_hint_s``
            Recommended ± window around an expected t0 (seconds).
        ``profile``
            Resolved profile name (``'ae'``, ``'ps_log'``, or
            ``'refraction'``).

    Raises
    ------
    ValueError
        If *profile* is not one of the recognised strings.
    ValueError
        If *f_dominant* is non-positive.

    Examples
    --------
    >>> p = suggest_parameters(fs=8e6)
    >>> p['profile']
    'ae'
    >>> p['sta_s']
    5e-06
    >>> p = suggest_parameters(fs=1000, profile='refraction')
    >>> p['sta_s']
    0.01
    >>> p = suggest_parameters(fs=4000, f_dominant=500)
    >>> abs(p['sta_s'] - 2 / 500) < 1e-15
    True
    >>> abs(p['lta_s'] - 20 / 500) < 1e-15
    True
    """
    # ── auto-detect profile ───────────────────────────────────────────────────
    if profile is None:
        if fs > 100_000:
            profile = "ae"
        elif fs > 1_000:
            profile = "ps_log"
        else:
            profile = "refraction"

    if profile not in _PROFILE_DEFAULTS:
        raise ValueError(
            f"Unknown profile {profile!r}. "
            f"Valid choices: {sorted(_PROFILE_DEFAULTS)}."
        )

    params = dict(_PROFILE_DEFAULTS[profile])
    params["profile"] = profile

    # ── f_dominant override ───────────────────────────────────────────────────
    if f_dominant is not None:
        f_dominant = float(f_dominant)
        if f_dominant <= 0:
            raise ValueError(
                f"f_dominant must be a positive frequency in Hz, got {f_dominant!r}."
            )
        params["sta_s"] = 2.0 / f_dominant
        params["lta_s"] = 20.0 / f_dominant

    return params
