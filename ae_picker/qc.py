"""
ae_picker.qc
============
Inter-trace quality-control utilities for first-arrival picks.

Public API
----------
flag_outlier_picks(df, offsets, pick_col, ...) -> pd.Series
    Flag picks that deviate significantly from a linear (or robust) moveout
    model fitted through all traces.
"""

import numpy as np
import pandas as pd

__all__ = ["flag_outlier_picks"]


def flag_outlier_picks(df, offsets, pick_col,
                       method="mad", n_mad=3.0,
                       model="linear",
                       return_model=False):
    """
    Flag first-arrival picks that are inconsistent with the inter-trace trend.

    A simple 1-D travel-time model (linear or robust) is fitted through
    ``(offset, pick_time)`` pairs.  Picks whose residual from the model
    exceeds ``n_mad`` times the Median Absolute Deviation (MAD) of all
    residuals are flagged as outliers.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame produced by ``batch.run`` or ``batch.pick_file``.
    offsets : array-like
        Source-receiver offsets or depths (metres), one per row of *df*.
    pick_col : str
        Name of the column in *df* to evaluate
        (e.g. ``'aic_pick_s'``, ``'refined_stalta_pick_s'``).
    method : {'mad'}
        Outlier-detection method.  Currently only ``'mad'`` is supported
        (Median Absolute Deviation).
    n_mad : float
        Number of MAD units beyond which a pick is flagged.  Default: 3.0.
    model : {'linear', 'robust'}
        Travel-time model used to fit the picks:

        ``'linear'``
            Ordinary least-squares line via ``numpy.polyfit`` (degree 1).
        ``'robust'``
            Theil–Sen estimator via ``scipy.stats.theilslopes``.
            Requires SciPy; raises ``ImportError`` with an install hint
            if not available.
    return_model : bool
        When ``True``, also return a dict describing the fitted model.
        Default: ``False``.

    Returns
    -------
    outlier_mask : pandas.Series (bool)
        Index-aligned with *df*; ``True`` where the pick is flagged as an
        outlier.  Rows with ``NaN`` picks are always flagged ``True``.
    model_info : dict
        Only returned when *return_model* is ``True``::

            {
                'intercept':   float,   # intercept in seconds
                'slope':       float,   # slope in s/m
                'velocity_ms': float,   # apparent velocity = 1/slope (m/s),
                                        # np.inf when slope == 0
            }

    Raises
    ------
    ValueError
        If *pick_col* is not in *df*.
    ValueError
        If *offsets* length does not match the number of rows in *df*.
    ImportError
        If ``model='robust'`` and SciPy is not installed.

    Notes
    -----
    ``model='linear'`` uses ordinary least squares, which can be influenced by
    extreme outliers when the dataset is small or when more than ~25 % of
    picks are bad.  In those cases, prefer ``model='robust'`` (Theil–Sen),
    which is breakdown-resistant up to 29 % outliers.

    Examples
    --------
    >>> import pandas as pd
    >>> # 10 traces on a 1 ms / 10 m moveout; last trace is a cycle-skip
    >>> df = pd.DataFrame({'aic_pick_s': [0.010, 0.020, 0.030, 0.040, 0.050,
    ...                                    0.060, 0.070, 0.080, 0.090, 0.200]})
    >>> offsets = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    >>> mask = flag_outlier_picks(df, offsets, 'aic_pick_s', n_mad=3.0)
    >>> mask.iloc[-1]   # cycle-skip flagged
    True
    >>> mask.iloc[0]    # clean trace not flagged
    False
    """
    if pick_col not in df.columns:
        raise ValueError(f"Column {pick_col!r} not found in df.")

    offsets = np.asarray(offsets, dtype=float)
    if len(offsets) != len(df):
        raise ValueError(
            f"offsets has {len(offsets)} elements but df has {len(df)} rows."
        )

    if method != "mad":
        raise ValueError(f"method must be 'mad', got {method!r}.")

    picks = df[pick_col].values.astype(float)

    # Rows with NaN picks are immediately outliers; exclude from fitting
    valid = ~np.isnan(picks)

    if valid.sum() < 2:
        # Not enough data to fit a model — flag everything NaN as outlier,
        # mark valids as non-outliers (insufficient context)
        outlier = ~valid
        if return_model:
            return (pd.Series(outlier, index=df.index),
                    {"intercept": np.nan, "slope": np.nan, "velocity_ms": np.nan})
        return pd.Series(outlier, index=df.index)

    x_fit = offsets[valid]
    y_fit = picks[valid]

    # ── fit travel-time model ──────────────────────────────────────────────────
    if model == "linear":
        slope, intercept = np.polyfit(x_fit, y_fit, 1)
    elif model == "robust":
        try:
            from scipy.stats import theilslopes
        except ImportError:
            raise ImportError(
                "SciPy is required for model='robust'. "
                "Install it with: pip install scipy"
            )
        result = theilslopes(y_fit, x_fit)
        slope     = float(result.slope)
        intercept = float(result.intercept)
    else:
        raise ValueError(f"model must be 'linear' or 'robust', got {model!r}.")

    # ── residuals and MAD flagging ────────────────────────────────────────────
    predicted = slope * offsets + intercept
    residuals = picks - predicted   # NaN where picks are NaN

    res_valid = residuals[valid]
    mad = float(np.median(np.abs(res_valid - np.median(res_valid))))

    if mad == 0:
        # All valid picks lie exactly on the model line; flag nothing
        outlier = ~valid
    else:
        outlier = np.abs(residuals) > n_mad * mad
        outlier[~valid] = True   # NaN picks always flagged

    velocity_ms = (1.0 / slope) if slope != 0 else np.inf

    if return_model:
        model_info = {
            "intercept":   float(intercept),
            "slope":       float(slope),
            "velocity_ms": float(velocity_ms),
        }
        return pd.Series(outlier, index=df.index), model_info

    return pd.Series(outlier, index=df.index)
