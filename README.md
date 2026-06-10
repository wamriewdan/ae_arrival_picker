# ae-picker

[DOI: 10.5281/zenodo.20126965](https://doi.org/10.5281/zenodo.20126965)


`ae-picker` is a Python package for automatic first-arrival time picking in
acoustic emission (AE), microseismic, and seismic waveform data.

It includes waveform readers for AE and standard seismic formats, multiple
picking algorithms, frequency-adaptive presets, batch processing, QC utilities,
and plotting helpers for review and comparison.

## Features

- AIC, envelope-based, and STA/LTA picking workflows
- **Frequency-adaptive presets** for AE, PS log, and refraction surveys
- Readers for AE text formats, SEG-Y, and miniSEED (ObsPy-backed)
- Batch processing with custom reader and glob support; writes `picks_summary.csv`
- **Inter-trace QC** with outlier flagging via linear or robust moveout models
- Plot helpers: waveform review, picker comparison, wiggle section, STA/LTA diagnostic, spectrogram
- A command-line interface exposed as `ae-picker`

## Algorithm notes

### AIC picker

The AIC picker treats the waveform as two segments split at sample `k` and
selects the split that minimizes the Akaike Information Criterion:

```text
AIC(k) = k * log(var(x[0:k])) + (N - k - 1) * log(var(x[k+1:N]))
```

Short explanation:

- Before the arrival, the signal is assumed to behave like background noise.
- After the arrival, the variance changes because the waveform contains the
  first motion and subsequent energy.
- The best onset is the sample where this two-segment model fits best.

### Refined STA/LTA picker

The refined STA/LTA picker in this package is a two-stage method:

1. A recursive STA/LTA trigger finds a coarse onset window.
2. A Hilbert-envelope refinement step walks back to the earliest persistent
   onset inside that window.

The STA/LTA trigger is based on the ratio

```text
CFT(i) = STA(i) / LTA(i)
```

where the fallback implementation in this package uses

```text
STA(i) = mean(|x[i-nsta:i]|)
LTA(i) = mean(|x[i-nlta:i]|)
```

and triggers when `CFT(i)` rises above the on-threshold and ends when it falls
below the off-threshold.

The envelope-refinement stage uses the analytic-signal envelope

```text
e(i) = |H{x}(i)|
```

and robust noise statistics

```text
sigma ~= 1.4826 * MAD(e_noise)
thr_high = median(e_noise) + k_high * sigma
thr_low  = median(e_noise) + k_low  * sigma
```

Short explanation:

- STA/LTA provides a stable first trigger for emergent arrivals.
- The envelope stage then refines that trigger by requiring a persistent rise
  above a lower threshold, which improves onset timing.
- In this repository, the "refined STA/LTA" picker is this package-specific
  combination of recursive STA/LTA and envelope-based onset refinement.

Both `aic_picker` and `refined_stalta_picker` accept optional moveout-guided
search window parameters (`expected_t0_s`, `search_window_s`, `time`) to
constrain picking to physically expected arrival times.


## Frequency-Adaptive Presets

`suggest_parameters` returns sensible STA/LTA window lengths and thresholds
calibrated to the data's sampling rate and dominant frequency.  Passing the
result as kwargs to the pickers removes all manual unit-conversion.

```python
from ae_picker import suggest_parameters, refined_stalta_picker

params = suggest_parameters(fs=1000, profile="refraction")
# {'sta_s': 0.01, 'lta_s': 0.1, 'min_persist_ms': 5.0, ...}

pick_idx, cft, triggers = refined_stalta_picker(
    amp, time,
    sta_s=params["sta_s"],
    lta_s=params["lta_s"],
)
```

Supported profiles:

| Profile | Typical `fs` | `f_dominant` | `sta_s` | `lta_s` |
|---|---|---|---|---|
| `"ae"` | 200 kHz – 8 MHz | 100 kHz – 1 MHz | 2 – 10 µs | 50 – 200 µs |
| `"ps_log"` | 4 – 20 kHz | 500 Hz – 5 kHz | 0.2 – 2 ms | 2 – 20 ms |
| `"refraction"` | 1 – 10 kHz | 50 – 500 Hz | 2 – 20 ms | 20 – 200 ms |

When `profile=None`, the profile is auto-detected from `fs`.  Supplying
`f_dominant` overrides the profile defaults with `sta_s = 2 / f_dominant`
and `lta_s = 20 / f_dominant`.


## Plotting

### `plot_section` — multi-trace wiggle section

Produces a standard wiggle-trace section for refraction surveys, PS logs, or
any multi-channel dataset where moveout is visible.  Positive lobes are filled
(SEG convention); pick times are overlaid as scatter markers.

```python
from ae_picker import plot_section

fig = plot_section(
    channels,
    picks={"AIC": aic_picks},
    offsets_m=[0, 5, 10, 15, 20],   # source-receiver offsets (m); None → index
    time_axis="x",                   # "x" (refraction) or "y" (PS log / VSP)
    normalise="trace",               # "trace" | "global" | "none"
    clip=0.95,
    savepath="section.png",
)
```

Key parameters:

| Parameter | Description | Default |
|---|---|---|
| `offsets_m` | Source-receiver offsets or depths (m) | `None` (sequential index) |
| `time_axis` | `"x"` → time horizontal; `"y"` → time vertical | `"x"` |
| `wiggle_scale` | Max excursion as fraction of trace spacing; `None` → auto | `None` |
| `normalise` | Amplitude normalisation scope | `"trace"` |
| `clip` | Amplitude clip percentile to suppress spikes | `0.95` |
| `savepath` | File path to save the figure | `None` |

---

### `plot_sta_lta_overlay` — STA/LTA diagnostic plot

Produces a three-panel figure for diagnosing the STA/LTA picker on a single
channel:

1. **Waveform** — raw amplitude time series.
2. **STA & LTA envelopes** — centered moving averages of `|amp|` showing how
   the short- and long-term averages evolve around the arrival.
3. **STA/LTA characteristic function** — the recursive ratio with trigger-on
   and trigger-off threshold lines, and optional onset marker lines.

```python
from ae_picker import refined_stalta_picker, plot_sta_lta_overlay

amp, time = channels[0]["amp"], channels[0]["time"]
N = len(time)

pick_idx, cft, triggers = refined_stalta_picker(
    amp, time,
    sta_s=2e-6, lta_s=100e-6,
    sta_thresh=2.5, lta_thresh=1.5,
    search_start=int(0.5 * N),
)

fig, axes = plot_sta_lta_overlay(
    amp, time,
    sta_s=2e-6, lta_s=100e-6,
    sta_thresh=2.5, lta_thresh=1.5,
    triggers=triggers,
    tlim_us=(-100, 200),
)
```

Key parameters:

| Parameter | Description | Default |
|---|---|---|
| `sta_s`, `lta_s` | STA and LTA window lengths (seconds) | `5e-6`, `50e-6` |
| `sta_thresh`, `lta_thresh` | Trigger on/off threshold lines | `2.5`, `1.5` |
| `triggers` | `[[onset, offset], ...]` sample pairs from `refined_stalta_picker` | `None` |
| `manual_picks_s` | Additional pick times (seconds) to overlay | `None` |
| `tlim_us` | X-axis limits in microseconds relative to `time[0]` | `None` (full record) |
| `savepath` | File path to save the figure | `None` |

---

### `plot_wave_and_spectrogram` — waveform + PSD spectrogram

Produces a two-panel figure with the normalised waveform above and the
Power Spectral Density spectrogram (Hann window, Welch method) below, sharing
the same time axis.  Suitable for journal figures.

```python
from ae_picker import plot_wave_and_spectrogram

fig, (ax_wave, ax_spec) = plot_wave_and_spectrogram(
    amp, time,
    fmax=500e3,
    tmin=time[0],
    tmax=time[-1],
    dyn_range_db=80,
    title="Channel 1",
    savepath="ch1_spec.png",
)
```

Key parameters:

| Parameter | Description | Default |
|---|---|---|
| `fmax` | Upper frequency limit for the spectrogram (Hz) | `None` (Nyquist) |
| `nperseg` | STFT window length (samples); `None` → auto (~2 ms) | `None` |
| `overlap` | Fractional window overlap `[0, 1)` | `0.90` |
| `dyn_range_db` | Colour dynamic range in dB | `80` |
| `norm_percentile` | Percentile used to normalise the waveform | `99.9` |
| `tmin`, `tmax` | Display window start/end (seconds) | `None` (full record) |
| `savepath` | File path to save the figure | `None` |


## QC — inter-trace outlier flagging

`flag_outlier_picks` fits a linear (or robust Theil–Sen) moveout model to a
column of pick times and flags traces whose residuals exceed a MAD-based
threshold.

```python
from ae_picker import flag_outlier_picks

mask = flag_outlier_picks(
    df,
    offsets=[0, 5, 10, 15, 20],
    pick_col="aic_pick_s",
    method="mad",
    n_mad=3.0,
    model="linear",        # "linear" | "robust" (scipy Theil–Sen)
)
# mask is a boolean Series; True = outlier
clean_df = df[~mask]
```

Pass `return_model=True` to also retrieve the fitted velocity, intercept, and slope.


## Installation

Install from PyPI after the package is published:

```bash
pip install ae-picker
```

Install from a local checkout:

```bash
git clone https://github.com/wamriewdan/ae_arrival_picker
cd ae_arrival_picker
pip install .
```

Optional extras:

```bash
pip install -e ".[dev]"
pip install ".[obspy]"   # enables read_segy() and read_mseed()
```

The package is installed as `ae-picker` but imported as `ae_picker`.

## Quick start

Batch process a dataset with the Python API:

```python
from ae_picker import run

df = run("path/to/your/data_dir", plot=False, save_plots=False)
print(df.head())
```

Run the CLI:

```bash
ae-picker path/to/your/data_dir --save-plots
```

Pick a single AE file:

```python
from ae_picker import aic_picker, plot_channels, read_unfiltered

meta, channels = read_unfiltered("my_ae_file.txt")

aic_picks = []
for ch in channels:
    idx, _ = aic_picker(ch["amp"], search_start=1, search_end=10)
    aic_picks.append(idx)
    print(f"P arrival at {ch['time'][idx] * 1e6:.1f} us")

fig = plot_channels(meta, channels, picks={"AIC": aic_picks})
fig.savefig("picks.png", dpi=150)
```

Process a SEG-Y refraction dataset:

```python
from ae_picker import read_segy, suggest_parameters, refined_stalta_picker, run

params = suggest_parameters(fs=1000, profile="refraction")

df = run(
    "refraction_data/",
    reader_fn=read_segy,
    file_glob="*.segy",
    sta_s=params["sta_s"],
    lta_s=params["lta_s"],
)
```


## Batch input layout

The default batch runner expects a root directory containing one or both of these
subdirectories:

```text
your_data_dir/
|-- unfiltered_signals/
`-- filtered_signals/
```

For custom formats, pass `reader_fn` and `file_glob` to `run()`:

```python
run("my_dir/", reader_fn=my_reader, file_glob="*.dat")
```

`run()` and `ae-picker` scan the target folder, infer the signal type from the
folder name or the `reader_fn` label, and write results to
`picks_output/picks_summary.csv` by default.

## Seismic readers (ObsPy)

`read_segy` and `read_mseed` require the optional ObsPy dependency
(`pip install "ae-picker[obspy]"`).  Both return the same
`(meta, channels)` contract as the built-in readers.

```python
from ae_picker import read_segy, read_mseed

meta, channels = read_segy("line01.segy")
meta, channels = read_mseed("station.mseed")
```

Both raise `ImportError` with an install hint if ObsPy is not available.

## References

1. Maeda, N. (1985). *A Method for Reading and Checking Phase Time in
   Auto-Processing System of Seismic Wave Data*. Zisin, 38(3), 365-379.
   https://doi.org/10.4294/zisin1948.38.3_365
2. Allen, R. V. (1978). *Automatic Earthquake Recognition and Timing from
   Single Traces*. Bulletin of the Seismological Society of America, 68(5),
   1521-1532. https://doi.org/10.1785/BSSA0680051521
3. Beyreuther, M., Barsch, R., Krischer, L., Megies, T., Behr, Y., and
   Wassermann, J. (2010). *ObsPy: A Python Toolbox for Seismology*.
   Seismological Research Letters, 81(3), 530-533.
   https://doi.org/10.1785/gssrl.81.3.530

## Citation

If you use `ae-picker`, please cite:

1. Wamriew, D. (2026). ae-picker. Zenodo. https://doi.org/10.5281/zenodo.20126965

## License

MIT
