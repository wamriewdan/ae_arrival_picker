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
    plot_section,
    plot_sta_lta_overlay,
    plot_wave_and_spectrogram,
)
from .presets import suggest_parameters
from .qc import flag_outlier_picks

# ObsPy-backed readers are optional; they raise ImportError at call time
# if ObsPy is not installed.
try:
    from .io import read_mseed, read_segy
except ImportError:
    pass

__all__ = [
    # io — built-in readers
    "read_unfiltered",
    "read_filtered",
    # io — optional seismic readers (raise ImportError at call time if ObsPy absent)
    "read_segy",
    "read_mseed",
    # pickers
    "aic_picker",
    "prepend_noise_aic_picker",
    "envelope_onset_picker",
    "envelope_offset_picker",
    "stalta_picker",
    "refined_stalta_picker",
    # plot
    "plot_channels",
    "plot_comparison",
    "plot_onset_zoom",
    "plot_section",
    "plot_sta_lta_overlay",
    "plot_wave_and_spectrogram",
    # batch
    "pick_file",
    "run",
    # presets
    "suggest_parameters",
    # qc
    "flag_outlier_picks",
]

__version__ = "0.2.0"
__author__ = "Daniel Wamriew, wamriewdan@gmail.com"
