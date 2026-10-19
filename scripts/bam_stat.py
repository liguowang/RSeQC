#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Summarize mapping statistics for a BAM or SAM alignment file."""

from __future__ import annotations
import argparse
import sys
from pathlib import Path
from typing import Sequence
from qcmodule import SAM
from importlib.metadata import version

__author__ = "Liguo Wang"
__version__ = version("RSeQC")

def build_parser() -> argparse.ArgumentParser:
    """Create and configure the command-line argument parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-i",
        "--input-file",
        required=True,
        type=Path,
        help="Alignment file in BAM or SAM format.",
    )
    parser.add_argument(
        "-q",
        "--mapq",
        dest="map_qual",
        type=int,
        default=30,
        metavar="INT",
        help=(
            "Minimum mapping quality (Phred-scaled) used to identify "
            "uniquely mapped reads. Default: %(default)s"
        ),
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    return parser


def validate_args(
    parser: argparse.ArgumentParser,
    input_file: Path,
    map_qual: int,
) -> None:
    """Validate command-line arguments."""
    if not input_file.is_file():
        parser.error(f"input file does not exist or is not a file: {input_file}")

    if map_qual < 0:
        parser.error("--mapq must be zero or greater")


def main(argv: Sequence[str] | None = None) -> int:
    """Run the BAM/SAM statistics command."""
    parser = build_parser()
    args = parser.parse_args(argv)

    validate_args(parser, args.input_file, args.map_qual)

    try:
        alignment = SAM.ParseBAM(str(args.input_file))
        alignment.stat(q_cut=args.map_qual)
    except (OSError, ValueError) as exc:
        parser.exit(1, f"{parser.prog}: error: {exc}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
