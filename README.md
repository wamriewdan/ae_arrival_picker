# ae-picker

`ae-picker` is a Python package for automatic first-arrival time picking in
acoustic emission (AE) and microseismic waveform data.

It provides multiple picking algorithms, format-agnostic readers, batch
processing helpers, and plotting utilities in a package that can be installed
from PyPI or directly from source.

## Features

- AIC, envelope-based, and STA/LTA picking workflows
- Readers for the repository's filtered and unfiltered waveform formats
- Batch processing that writes a single `picks_summary.csv`
- Plot helpers for waveform review and picker comparison
- A command-line interface exposed as `ae-picker`

## Installation

Install from PyPI after the package has been published:

```bash
pip install ae-picker
```

Install from a local checkout:

```bash
git clone https://github.com/wamriewdan/ae_arrival_picker
cd ae_arrival_picker
pip install .
```

For development:

```bash
pip install -e ".[dev]"
```

If you want the optional ObsPy-backed recursive STA/LTA trigger:

```bash
pip install "ae-picker[obspy]"
```

## Quick start

Use the Python API:

```python
from ae_picker import run

df = run("path/to/your/data_dir", plot=False, save_plots=False)
print(df.head())
```

Use the CLI:

```bash
ae-picker path/to/your/data_dir --save-plots
```

Single-file usage:

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

## Project layout

```text
ae_arrival_picker/
|-- ae_picker/
|   |-- __init__.py
|   |-- __main__.py
|   |-- batch.py
|   |-- cli.py
|   |-- io.py
|   |-- pickers.py
|   `-- plot.py
|-- notebooks/
|-- pyproject.toml
|-- README.md
`-- LICENSE
```

## Publishing to PyPI

`pip install ae-picker` only works for end users after the project is
published under that distribution name on PyPI.

Build the distribution files:

```bash
python -m build
```

Upload them:

```bash
python -m twine upload dist/*
```

If `ae-picker` is already taken on PyPI, you will need to choose a different
distribution name. The Python import can still remain `ae_picker`.

For the repository's exact GitHub Actions Trusted Publishing setup, see
`RELEASING.md`.

For a safe pre-release check, the repo also includes a separate manual
TestPyPI workflow.

Routine validation is handled by a separate CI workflow that runs on every
push and pull request.

## License

MIT
