# ae-picker v2 — Changes

This document describes every change introduced on the `ae-picker-v2` branch relative to the last `main` release.  All changes are backward-compatible: no existing call sites break.

---

## New modules

### `ae_picker/presets.py` — frequency-adaptive parameter presets

**`suggest_parameters(fs, f_dominant=None, profile=None) -> dict`**

Returns a complete set of picker parameters scaled to the recording's sampling rate, removing the need to manually convert window lengths between AE and seismic scales.

| Parameter | Role |
|---|---|
| `fs` | Sampling rate in Hz |
| `f_dominant` | Optional dominant frequency; overrides profile defaults with `sta_s = 2/f_dom`, `lta_s = 20/f_dom` |
| `profile` | `'ae'`, `'ps_log'`, `'refraction'`, or `None` (auto-detected from `fs`) |

Auto-detection rules: `fs > 100 000` Hz → `'ae'`; `fs > 1 000` Hz → `'ps_log'`; otherwise → `'refraction'`.

Returned dict keys: `sta_s`, `lta_s`, `min_persist_ms`, `noise_lookback_s`, `pre_pick_s`, `post_pick_s`, `search_window_hint_s`, `profile`.

Profile defaults:

| Profile | `sta_s` | `lta_s` | `min_persist_ms` |
|---|---|---|---|
| `ae` | 5 µs | 100 µs | 0.1 |
| `ps_log` | 1 ms | 10 ms | 1.0 |
| `refraction` | 10 ms | 100 ms | 5.0 |

---

### `ae_picker/qc.py` — inter-trace outlier flagging

**`flag_outlier_picks(df, offsets, pick_col, method='mad', n_mad=3.0, model='linear', return_model=False) -> pd.Series`**

Fits a 1-D travel-time model across traces and flags picks whose residual exceeds `n_mad` times the Median Absolute Deviation (MAD).  Rows with `NaN` picks are always flagged.

| Parameter | Description |
|---|---|
| `df` | DataFrame from `batch.run` / `batch.pick_file` |
| `offsets` | Source-receiver offsets or depths (m), one per row |
| `pick_col` | Column to evaluate, e.g. `'aic_pick_s'` |
| `model` | `'linear'` (numpy OLS) or `'robust'` (scipy Theil–Sen; requires SciPy) |
| `return_model` | When `True`, also returns `{'intercept', 'slope', 'velocity_ms'}` |

Returns an index-aligned boolean `pd.Series`; `True` = outlier.

---

## Modified modules

### `ae_picker/io.py` — two new seismic readers

Both require ObsPy (`pip install "ae-picker[obspy]"`) and raise `ImportError` with an install hint if it is absent.  Both return the same `(meta, channels)` contract as the existing readers.

**`read_segy(filepath, component=None) -> (meta, channels)`**

Reads a SEG-Y file via ObsPy.  `component` optionally selects a single trace by component code.

**`read_mseed(filepath) -> (meta, channels)`**

Reads a miniSEED file via ObsPy.

Both are exported from `ae_picker/__init__.py` inside a `try/except ImportError` guard (import-time failure is silenced; the functions raise at call time).

---

### `ae_picker/pickers.py` — moveout-guided search windows

**`aic_picker`** gains three new optional keyword arguments:

```python
aic_picker(amp, search_start=0, search_end=None,
           expected_t0_s=None, search_window_s=None, time=None)
```

When all three new kwargs are provided, `search_start` and `search_end` are derived from the time array as `[expected_t0_s − search_window_s, expected_t0_s + search_window_s]`.  When any are absent the existing `search_start` / `search_end` index behaviour is unchanged.

**`refined_stalta_picker`** gains two new optional keyword arguments:

```python
refined_stalta_picker(..., expected_t0_s=None, search_window_s=None)
```

When both are provided, triggers whose onset falls outside `[expected_t0_s − search_window_s, expected_t0_s + search_window_s]` are discarded before the envelope-refinement stage.

All new kwargs default to `None`; no existing calls break.

---

### `ae_picker/plot.py` — `plot_section`

**`plot_section(channels, picks=None, offsets_m=None, time_axis='x', wiggle_scale=None, normalise='trace', clip=0.95, savepath=None) -> fig`**

Produces a standard wiggle-trace section.  Positive lobes are filled (SEG convention).  Pick times are overlaid as scatter markers.

| Parameter | Description | Default |
|---|---|---|
| `offsets_m` | Source-receiver offsets or depths (m) | `None` → sequential index |
| `time_axis` | `'x'` → time horizontal (refraction); `'y'` → time vertical (PS log / VSP) | `'x'` |
| `wiggle_scale` | Max excursion as fraction of trace spacing; `None` → auto (80% of spacing) | `None` |
| `normalise` | `'trace'`, `'global'`, or `'none'` | `'trace'` |
| `clip` | Amplitude clip percentile | `0.95` |
| `savepath` | Output file path | `None` |

---

### `ae_picker/batch.py` — generic reader and profile support

**`pick_file`** gains one new keyword argument:

```python
pick_file(filepath, signal_type=None, reader_fn=None, ...)
```

- `reader_fn`: any callable with the signature `meta, channels = reader_fn(filepath)`.  When provided, the built-in readers are bypassed; `signal_type` is used only as a display label (defaults to `'custom'`).

**`run`** gains three new keyword arguments:

```python
run(data_dir, ..., reader_fn=None, file_glob='*.txt', profile=None, sta_s=None, lta_s=None)
```

- `reader_fn`: forwarded to `pick_file`; bypasses built-in readers.
- `file_glob`: glob pattern for file selection inside each sub-directory.  Default `'*.txt'` preserves existing behaviour; pass `'*.segy'` for refraction surveys.
- `profile`: when set, `suggest_parameters` is called on the first file's measured sampling rate and its `sta_s` / `lta_s` values become the defaults.  Explicitly supplied `sta_s` / `lta_s` always take precedence.

The output DataFrame gains an `envelope_pick_s` / `envelope_pick_us` column pair (the envelope onset picker was already run but not recorded in v1).

---

### `ae_picker/__init__.py` — updated exports

All new public names added to `__all__`:

| Name | Source |
|---|---|
| `suggest_parameters` | `presets.py` |
| `flag_outlier_picks` | `qc.py` |
| `read_segy` | `io.py` (ObsPy-guarded) |
| `read_mseed` | `io.py` (ObsPy-guarded) |
| `plot_section` | `plot.py` |

---

## Compatibility

- All existing calls to `run()`, `pick_file()`, `aic_picker()`, and `refined_stalta_picker()` with no new kwargs produce identical output.
- The `envelope_pick_s` / `envelope_pick_us` columns are new in the CSV; downstream code that selects specific columns by name is unaffected.
- `read_segy` and `read_mseed` are listed in `__all__` but only imported at module load if ObsPy is present.  Code that checks `hasattr(ae_picker, 'read_segy')` before calling will work correctly regardless of whether ObsPy is installed.
