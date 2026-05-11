"""Command-line interface for ae_picker."""

from __future__ import annotations

import argparse
from pathlib import Path

from .batch import run


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ae-picker",
        description=(
            "Batch-pick first arrivals from acoustic emission waveform folders "
            "and write picks_summary.csv."
        ),
    )
    parser.add_argument(
        "data_dir",
        help="Root directory containing waveform folders such as unfiltered_signals/.",
    )
    parser.add_argument(
        "--signal-dir",
        dest="signal_dirs",
        action="append",
        default=None,
        help=(
            "Signal subdirectory to process. Repeat the flag to pass multiple "
            "directories. Defaults to unfiltered_signals and filtered_signals."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory for picks_summary.csv and any saved figures.",
    )
    parser.add_argument(
        "--picker",
        dest="pickers",
        action="append",
        default=None,
        choices=("aic", "aic_prepend", "envelope", "stalta", "refined_stalta"),
        help="Picker to enable. Repeat the flag to run multiple pickers.",
    )
    parser.add_argument(
        "--prepend-duration-us",
        type=float,
        default=50.0,
        help="Synthetic noise prepend duration in microseconds for aic_prepend.",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Display plots interactively while processing.",
    )
    parser.add_argument(
        "--save-plots",
        action="store_true",
        help="Save one pick plot per waveform into the output directory.",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    signal_dirs = tuple(args.signal_dirs) if args.signal_dirs else (
        "unfiltered_signals",
        "filtered_signals",
    )
    pickers = tuple(args.pickers) if args.pickers else ("aic", "refined_stalta")
    output_dir = Path(args.output_dir) if args.output_dir else None
    should_plot = args.plot or args.save_plots

    df = run(
        data_dir=Path(args.data_dir),
        signal_dirs=signal_dirs,
        output_dir=output_dir,
        pickers=pickers,
        prepend_duration_s=args.prepend_duration_us * 1e-6,
        plot=should_plot,
        save_plots=args.save_plots,
        show_plots=args.plot,
    )

    print(f"Processed {len(df)} channel rows.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
