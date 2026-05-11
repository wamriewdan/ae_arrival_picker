# ae-picker

`ae-picker` is a Python package for automatic first-arrival time picking in
acoustic emission (AE) and microseismic waveform data.

It includes waveform readers, multiple picking algorithms, batch processing,
and plotting helpers for review and comparison.

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
LTA(i) = mean(|x[i-nlta:i-nsta]|)
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


## Features

- AIC, envelope-based, and STA/LTA picking workflows
- Readers for the repository's filtered and unfiltered waveform formats
- Batch processing that writes a single `picks_summary.csv`
- Plot helpers for waveform review and picker comparison
- A command-line interface exposed as `ae-picker`

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
pip install ".[obspy]"
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

Pick a single file:

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


## Batch input layout

The batch runner expects a root directory containing one or both of these
subdirectories:

```text
your_data_dir/
|-- unfiltered_signals/
`-- filtered_signals/
```

`run()` and `ae-picker` scan those folders, infer the signal type from the
folder name, and write results to `picks_output/picks_summary.csv` by default.

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

## License

MIT
