# ae_picker

**ae_picker** is an open-source Python package for automatic first-arrival time picking in acoustic emission (AE) and microseismic waveform data.

It provides three complementary picking algorithms, format-agnostic readers, and publication-quality plots — all in a clean, pip-installable package designed to be adapted to any AE dataset.

---

## Features

- **Three picking algorithms** — AIC (Maeda 1985), energy-onset ratio, and STA/LTA
- **Format-agnostic design** — pickers operate on NumPy arrays; plug in your own reader for any file format
- **Batch processing** — scan a directory tree and export all picks to a single CSV
- **Publication-quality plots** — per-file channel plots, scatter comparison, and onset zoom browser
- **Demo notebook** — end-to-end walkthrough using real AE data from borehole monitoring experiments

---

## Installation

```bash
git clone https://github.com/wamriewdan/ae_arrival_picker
cd ae_arrival_picker
pip install -e .
```

For running the demo notebook, also install the optional dependencies:

```bash
pip install -e ".[dev]"
```

---

## Quick start

```python
from ae_picker import run

# Batch-pick all .txt files in a data directory
df = run("path/to/your/data_dir")
# → saves picks_output/picks_summary.csv and one PNG per file
```

Single-file usage:

```python
from ae_picker import read_unfiltered, aic_picker, plot_channels

meta, channels = read_unfiltered("my_ae_file.txt")

aic_picks = []
for ch in channels:
    idx, _ = aic_picker(ch["amp"], search_start=1, search_end=10)
    aic_picks.append(idx)
    print(f"P arrival at {ch['time'][idx]*1e6:.1f} µs")

fig = plot_channels(meta, channels, picks={"AIC": aic_picks})
fig.savefig("picks.png", dpi=150)
```

---

## Supported file formats

| Format | Channels | Description |
|---|---|---|
| `unfiltered` | 4 | 5-line metadata header; tab-separated time/amplitude pairs; record starts at P onset |
| `filtered` | 2 | 1-line column header; centred at t = 0; ±250 µs window |

Both formats return the same `(meta, channels)` structure — see `ae_picker/io.py` for the full specification.

**Extending to your own dataset:** Write a reader function that returns `(meta_dict, channels_list)` where each channel is `{'time': ndarray, 'amp': ndarray}`.  All pickers and plotting functions will work immediately with your data.

---

## Algorithms

### AIC (Akaike Information Criterion)
*Maeda (1985)* — treats the waveform as two concatenated stationary Gaussian processes and finds the split point that minimises:

```
AIC(k) = k · log(var(x[0:k])) + (N−k−1) · log(var(x[k+1:N]))
```

**Best for:** Records where the P onset is the dominant variance change.
**Key parameter:** `search_end` — restrict to the onset region to prevent drift to later coda features.

### Energy-onset ratio
Finds the first sample where instantaneous energy exceeds a multiple of the pre-onset background energy:

```
amp[i]² > threshold × mean(amp[:n_noise]²)
```

**Best for:** Short records that start immediately at the onset (leaving too few samples for AIC).
**Key parameters:** `noise_window_s`, `threshold`.

### STA/LTA
Classical Short-Term Average / Long-Term Average ratio trigger:

```
ratio[i] = mean(|amp[i-sta_n : i]|) / mean(|amp[i-lta_n : i-sta_n]|)
```

Pick fires at the first sample where `ratio ≥ threshold`.

**Best for:** Longer continuous records; widely used in seismology.
**Key parameters:** `sta_s`, `lta_s`, `threshold`.

---

## Demo notebook

`notebooks/demo_maria_project.ipynb` walks through the full workflow using AE data from borehole monitoring experiments (Maria dataset, KAUST RockGem group):

1. Load and inspect both file formats
2. Pick a single file with each algorithm
3. Batch-process all files
4. Export results to CSV
5. Compare computed picks against reference header picks

Set `DATA_DIR` in the second cell to the folder containing your `unfiltered_signals/` and `filtered_signals/` subdirectories before running.

---

## Project structure

```
ae_picker/
├── ae_picker/
│   ├── __init__.py    # public API exports
│   ├── io.py          # file readers
│   ├── pickers.py     # AIC, energy-onset, STA/LTA
│   ├── plot.py        # visualisation helpers
│   └── batch.py       # pick_file() and run()
├── notebooks/
│   └── demo_maria_project.ipynb
├── pyproject.toml
├── README.md
└── .gitignore
```

---

## Citation

If you use **ae_picker** in your research, please cite:

> Wamriew, D. (2026). *ae_picker: first-arrival time picking for acoustic emission waveforms*. GitHub. https://github.com/wamriewdan/ae_arrival_picker

the underlying AIC algorithm:

> Maeda, N. (1985). A method for reading and checking phase times in auto-processing system of seismic wave data. *Zisin (Journal of the Seismological Society of Japan)*, 38, 365–379.

and the STA/LTA method:

> Allen, R. (1978). Automatic earthquake recognition and timing from single traces. *Bulletin of the Seismological Society of America*, 68(5), 1521–1532.

and its modern implementation context:

> Beyreuther, M., Barsch, R., Krischer, L., Megies, T., Behr, Y., & Wassermann, J. (2010). ObsPy: A Python toolbox for seismology. *Seismological Research Letters*, 81(3), 530–533.

---

## License

MIT — see `LICENSE` for details.
