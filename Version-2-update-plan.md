# ae_picker v2 — Seismic Use Case Extension

## Context

A potential user asked about applying ae_picker to P-wave refraction surveys and PS/sonic logs. The conversation exposed two structural gaps:

1. **Parameter scaling**: All STA/LTA defaults target AE data at ~8 MHz. Seismic data runs at 1–20 kHz with dominant frequencies of 50 Hz–5 kHz — requiring windows 100–1000× larger. Users must manually compute every parameter.
2. **No trace-to-trace consistency**: Each channel is picked independently with no awareness of moveout. Refraction and PS log QC depends on comparing picks across traces.

v2 adds a preset/scaling layer, a generic reader interface, standard seismic file readers, a section plot, moveout-guided search windows, and an outlier-flagging utility — without breaking any existing calls.

---

## Files to Modify / Create

| File | Action |
|---|---|
| `ae_picker/presets.py` | **New** — `suggest_parameters()` |
| `ae_picker/qc.py` | **New** — `flag_outlier_picks()` |
| `ae_picker/pickers.py` | **Modify** — add `expected_t0_s`, `search_window_s`, `time` to `aic_picker` and `refined_stalta_picker` |
| `ae_picker/io.py` | **Modify** — add `read_segy()` and `read_mseed()` |
| `ae_picker/batch.py` | **Modify** — add `reader_fn`, `file_glob`, `profile` kwargs to `pick_file` and `run` |
| `ae_picker/plot.py` | **Modify** — add `plot_section()` |
| `ae_picker/__init__.py` | **Modify** — export all new public names |

---

## 1. `presets.py` — Frequency-Adaptive Presets (Must Have)

**Why**: removes the manual scaling burden; one parameter unlocks the whole picker.

```python
suggest_parameters(fs, f_dominant=None, profile=None) -> dict
```

- `profile`: `"ae"` | `"refraction"` | `"ps_log"` | `None` (auto from `fs`)
- Auto-detect rule: `fs > 100_000` → `"ae"`, `fs > 1_000` → `"ps_log"`, else → `"refraction"`
- When `f_dominant` is given: `sta_s = 2 / f_dominant`, `lta_s = 20 / f_dominant` (overrides profile defaults)
- Returns dict with: `sta_s`, `lta_s`, `min_persist_ms`, `noise_lookback_s`, `pre_pick_s`, `post_pick_s`, `search_window_hint_s`

Target values per profile:

| Profile | fs | f_dom | sta_s | lta_s | min_persist_ms |
|---|---|---|---|---|---|
| `ae` | 200 kHz–8 MHz | 100 kHz–1 MHz | 2–10 µs | 50–200 µs | 0.05–0.3 |
| `ps_log` | 4–20 kHz | 500 Hz–5 kHz | 0.2–2 ms | 2–20 ms | 0.5–2 |
| `refraction` | 1–10 kHz | 50–500 Hz | 2–20 ms | 20–200 ms | 2–10 |

Export: add `suggest_parameters` to `__init__.py` imports and `__all__`.

---

## 2. Generic Batch Interface (Must Have)

**Why**: `pick_file` and `run` are locked to `.txt` Maria format. Any other data source requires a full reimplementation.

Pattern already exists in `plot_onset_zoom` (accepts `reader_fn` callable) — apply the same to `batch.py`.

**`pick_file` change:**
```python
pick_file(filepath, signal_type=None, reader_fn=None, ...)
```
- If `reader_fn` is provided: `meta, channels = reader_fn(filepath)`; `signal_type` used only as a label (defaults to `"custom"`).
- If `reader_fn=None`: existing `signal_type in {"unfiltered","filtered"}` branch unchanged.

**`run` change:**
```python
run(data_dir, signal_dirs=..., reader_fn=None, file_glob="*.txt", ...)
```
- `file_glob`: default `"*.txt"` preserves existing behaviour; pass `"*.segy"` for refraction.
- `profile` kwarg: calls `suggest_parameters` and uses result as defaults for `sta_s`, `lta_s`, etc. (caller-supplied values still override).

No existing call sites break.

---

## 3. ObsPy-Backed Readers (Should Have)

**Why**: SEG-Y is the dominant refraction format; miniSEED is common for PS logs. ObsPy is already an optional dep.

Add to bottom of `ae_picker/io.py`:

```python
read_segy(filepath, component=None) -> (meta, channels)
read_mseed(filepath) -> (meta, channels)
```

Both return the same `(meta, channels)` contract as existing readers:
- `channels`: list of `{"time": ndarray[s], "amp": ndarray[float]}`
- `meta`: dict with `filename`, `format`, and format-specific fields

Both raise `ImportError` with install hint (`pip install ae-picker[obspy]`) if ObsPy absent — same pattern as `refined_stalta_picker` in `pickers.py`.

Export: add both to `__init__.py` inside a `try/except ImportError` guard on the import; add to `__all__` unconditionally (they raise at call time, not import time).

---

## 4. `plot_section` — Multi-Trace Wiggle Section (Should Have)

**Why**: standard QC display for refraction and PS logs. `plot_channels` stacks traces by channel but doesn't show moveout.

Add to bottom of `ae_picker/plot.py`:

```python
plot_section(channels, picks=None,
             offsets_m=None,
             time_axis="x",
             wiggle_scale=None,
             normalise="trace",
             clip=0.95,
             savepath=None) -> fig
```

- `offsets_m`: list of floats (source-receiver offsets or depths). `None` → sequential index.
- `time_axis="x"`: time horizontal (refraction convention); `"y"`: time vertical (PS log / VSP convention).
- `wiggle_scale`: `None` → auto (max excursion = 80% of trace spacing).
- `normalise`: `"trace"` | `"global"` | `"none"`.
- `clip`: amplitude clip percentile (default 0.95) to suppress spikes.
- Fills positive wiggle lobes (standard SEG convention); marks picks as scatter markers at pick time.

Export: add `plot_section` to `__init__.py`.

---

## 5. Moveout-Guided Search Window (Nice to Have)

**Why**: constrains picks to physically expected arrival times; prevents false triggers in pre-arrival noise.

**`aic_picker` change** (add optional time-domain window):
```python
aic_picker(amp, search_start=0, search_end=None,
           expected_t0_s=None, search_window_s=None, time=None)
```
- When `expected_t0_s` + `search_window_s` + `time` are all provided, derive `search_start`/`search_end` from time array. Otherwise: unchanged behaviour.

**`refined_stalta_picker` change:**
```python
refined_stalta_picker(..., expected_t0_s=None, search_window_s=None, search_start=0)
```
- When provided, discards triggers outside `[expected_t0_s - search_window_s, expected_t0_s + search_window_s]` — consistent with existing `search_start` filter at `pickers.py:510`.

All new kwargs default to `None`; no existing calls break.

---

## 6. `qc.py` — Inter-Trace Consistency Check (Nice to Have)

**Why**: even with windowed pickers, single-trace outliers (noise bursts, cycle-skipping) need flagging before writing CSV.

```python
flag_outlier_picks(df, offsets, pick_col,
                   method="mad", n_mad=3.0,
                   model="linear",
                   return_model=False) -> pd.Series (bool, True = outlier)
```

- `df`: DataFrame from `batch.run` / `batch.pick_file`.
- `offsets`: array-like of offsets (m) or depths (m), same length as `df`.
- `pick_col`: column name to evaluate (e.g. `"aic_pick_s"`).
- `model`: `"linear"` (numpy polyfit) | `"robust"` (scipy Theil–Sen, import-guarded).
- Returns index-aligned boolean Series; does not mutate `df`.
- `return_model=True` → `(mask, {"intercept", "slope", "velocity_ms"})`.

Export: add `flag_outlier_picks` to `__init__.py`.

---

## Implementation Order

1. `presets.py` (standalone, no internal deps)
2. `pickers.py` moveout kwargs (no new imports)
3. `io.py` new readers (ObsPy guard)
4. `batch.py` `reader_fn` + `profile` (needs `presets.py`)
5. `plot.py` `plot_section` (matplotlib/numpy only)
6. `qc.py` (numpy; optional scipy guard)
7. `__init__.py` (last — ties all together)

---

## Verification

- `suggest_parameters(fs=1000, profile="refraction")` returns `sta_s` in ms range, `lta_s` ~10× larger.
- `suggest_parameters(fs=8e6)` returns values matching current AE defaults.
- `pick_file(fp, reader_fn=my_reader)` processes a non-.txt file without error.
- `run(dir, file_glob="*.segy", reader_fn=read_segy)` batches a SEG-Y dataset.
- `plot_section(channels, picks, offsets_m=[...], time_axis="x")` renders a wiggle section with pick markers.
- `flag_outlier_picks(df, offsets, "aic_pick_s")` returns a boolean Series of correct length with outliers flagged.
- All existing notebook cells and `batch.run(...)` calls with no new kwargs produce identical output.
