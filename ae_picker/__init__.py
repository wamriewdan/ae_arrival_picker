"""
ae_picker
=========
A Python package for first-arrival time picking in acoustic emission (AE)
and microseismic waveform data.

Modules
-------
io       : file readers  (``read_unfiltered``, ``read_filtered``)
pickers  : picking algorithms  (AIC, energy-onset, STA/LTA)
plot     : visualisation helpers
batch    : single-file and batch processing

Quick start
-----------
>>> from ae_picker import run
>>> df = run("path/to/data_dir")          # picks all .txt files, saves CSV + PNGs

>>> from ae_picker import read_unfiltered, aic_picker
>>> meta, channels = read_unfiltered("my_file.txt")
>>> pick_idx, aic = aic_picker(channels[0]["amp"], search_start=1, search_end=10)
>>> print(f"P arrival at {channels[0]['time'][pick_idx]*1e6:.1f} µs")
"""

from .io      import read_unfiltered, read_filtered
from .pickers import (aic_picker, prepend_noise_aic_picker,
                      energy_onset_picker,
                      envelope_onset_picker, envelope_offset_picker,
                      stalta_picker, refined_stalta_picker)
from .plot    import plot_channels, plot_comparison, plot_onset_zoom, plot_sta_lta_overlay
from .batch   import pick_file, run

__all__ = [
    # io
    "read_unfiltered",
    "read_filtered",
    # pickers
    "aic_picker",
    "prepend_noise_aic_picker",
    "energy_onset_picker",
    "stalta_picker",
    # plot
    "plot_channels",
    "plot_comparison",
    "plot_onset_zoom",
    # batch
    "pick_file",
    "run",
]

__version__ = "0.1.0"
__author__  = "RockGem Group, KAUST"
