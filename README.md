# ae-picker

`ae-picker` is a Python package for automatic first-arrival time picking in
acoustic emission (AE) and microseismic waveform data.

It includes waveform readers, multiple picking algorithms, batch processing,
and plotting helpers for review and comparison.

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

## License

MIT
