#!/usr/bin/env python3
"""
qPCR Analysis Pipeline - Main Entry Point

This script processes qPCR data from multiple instrument platforms,
merges with plate maps, computes statistics, performs quality control,
and generates comprehensive summary reports.

Usage:
    python qpcr_pipeline_main.py [workdir] [--internal-control TARGET_NAME]

    workdir: Optional path to working directory (defaults to current directory)
    --internal-control: Optional target name (e.g. RNaseP) to correct SQ Mean
        against. Must match a Target value produced by target_mapping.csv.
        Correction is disabled unless this is given.
"""
import argparse
from pathlib import Path

from src.config import PipelineConfig
from src.pipeline import QPCRPipeline


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Process qPCR data: merge with plate maps, QC, and summarize."
    )
    parser.add_argument(
        "workdir",
        nargs="?",
        default=None,
        help="Working directory containing plate_maps.csv and raw/ (defaults to the current directory).",
    )
    parser.add_argument(
        "--internal-control",
        metavar="TARGET_NAME",
        default=None,
        help=(
            "Enable internal control correction, normalizing SQ Mean against this target "
            "(e.g. RNaseP). Must match a Target value assigned via target_mapping.csv "
            "(case-insensitive). Disabled unless given."
        ),
    )
    return parser.parse_args()


def main():
    """Main entry point for the pipeline."""
    args = parse_args()
    workdir = Path(args.workdir).resolve() if args.workdir else None

    config = None
    if args.internal_control:
        # workdir here is just a placeholder to satisfy PipelineConfig's
        # required field — QPCRPipeline.__init__ re-targets the config (workdir
        # and every path derived from it) at the actually-resolved working
        # directory right after this.
        config = PipelineConfig(workdir=Path.cwd(), internal_control_target=args.internal_control)

    pipeline = QPCRPipeline(workdir=workdir, config=config)
    pipeline.run()


if __name__ == "__main__":
    main()
