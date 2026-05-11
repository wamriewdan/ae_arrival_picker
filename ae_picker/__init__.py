"""Public package exports for ae_picker."""

from .batch import pick_file, run
from .io import read_filtered, read_unfiltered
from .pickers import (
    aic_picker,
    envelope_offset_picker,
    envelope_onset_picker,
    prepend_noise_aic_picker,
    refined_stalta_picker,
    stalta_picker,
)
from .plot import (
    plot_channels,
    plot_comparison,
    plot_onset_zoom,
    plot_sta_lta_overlay,
)

__all__ = [
    "read_unfiltered",
    "read_filtered",
    "aic_picker",
    "prepend_noise_aic_picker",
    "envelope_onset_picker",
    "envelope_offset_picker",
    "stalta_picker",
    "refined_stalta_picker",
    "plot_channels",
    "plot_comparison",
    "plot_onset_zoom",
    "plot_sta_lta_overlay",
    "pick_file",
    "run",
]

__version__ = "0.2.0"
__author__ = "Daniel Wamriew, wamriewdan@gmail.com"
