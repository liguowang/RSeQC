#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Calculate per-cycle Phred quality-score distributions for aligned reads."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Sequence
from qcmodule import SAM
from importlib.metadata import version

__version__ = version("RSeQC")
__author__ = "Liguo Wang"


DESCRIPTION = (
    "Calculate the Phred quality-score distribution at each read position. "
    "Input reads should have a fixed, consistent length."
)

EPILOG = """
Example
-------
read_quality.py -i sample.bam -o sample -r 1 -q 30

Notes
-----
* --reduce affects only the generated boxplot data.
* A value of 1 retains maximum precision.
* Larger values reduce the size of vectors passed to R, using less memory
  at the cost of reduced plotting precision.
"""


def build_parser() -> argparse.ArgumentParser:
    """Create and configure the command-line parser."""
    parser = argparse.ArgumentParser(
        description=DESCRIPTION,
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--input-file",
        required=True,
        type=Path,
        help="Input alignment file in BAM or SAM format.",
    )
    parser.add_argument(
        "-o",
        "--out-prefix",
        dest="output_prefix",
        required=True,
        type=Path,
        help="Prefix for generated output files.",
    )
    parser.add_argument(
        "-r",
        "--reduce",
        dest="reduce_fold",
        type=int,
        default=1,
        metavar="INT",
        help=(
            "Ignore quality-score observations occurring fewer than this "
            "number of times when preparing boxplot data. "
            "Default: %(default)s"
        ),
    )
    parser.add_argument(
        "-q",
        "--mapq",
        dest="map_qual",
        type=int,
        default=30,
        metavar="INT",
        help=(
            "Minimum mapping quality for an alignment to be treated as "
            "uniquely mapped. Default: %(default)s"
        ),
    )
    parser.add_argument(
        "--skip-plot",
        action="store_true",
        help="Generate quality-profile data but do not execute the R script.",
    )
    parser.add_argument(
        "--rscript",
        default="Rscript",
        metavar="PATH",
        help="Rscript executable to use. Default: %(default)s",
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
    output_prefix: Path,
    reduce_fold: int,
    map_qual: int,
) -> None:
    """Validate command-line arguments."""
    if not input_file.is_file():
        parser.error(
            f"input alignment file does not exist or is not a file: "
            f"{input_file}"
        )

    if input_file.suffix.lower() not in {".bam", ".sam", ".cram"}:
        parser.error(
            "input alignment file must have a .bam, .sam, or .cram "
            f"extension: {input_file}"
        )

    if reduce_fold <= 0:
        parser.error("--reduce must be greater than zero")

    if map_qual < 0:
        parser.error("--mapq must be zero or greater")

    output_parent = output_prefix.parent
    if not output_parent.exists():
        parser.error(f"output directory does not exist: {output_parent}")

    if not output_parent.is_dir():
        parser.error(f"output parent is not a directory: {output_parent}")


def run_plot_script(
    parser: argparse.ArgumentParser,
    output_prefix: Path,
    rscript_executable: str,
) -> None:
    """Execute the generated R plotting script safely."""
    script_path = Path(f"{output_prefix}.qual.r")

    if not script_path.is_file():
        parser.exit(
            1,
            f"{parser.prog}: error: expected R script was not created: "
            f"{script_path}\n",
        )

    resolved_rscript = shutil.which(rscript_executable)
    if resolved_rscript is None:
        parser.exit(
            1,
            f"{parser.prog}: error: Rscript executable not found: "
            f"{rscript_executable}\n",
        )

    try:
        subprocess.run(
            [resolved_rscript, str(script_path)],
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        parser.exit(
            exc.returncode or 1,
            f"{parser.prog}: error: R plotting failed for {script_path}\n",
        )


def main(argv: Sequence[str] | None = None) -> int:
    """Run per-cycle read-quality analysis."""
    parser = build_parser()
    args = parser.parse_args(argv)

    validate_args(
        parser=parser,
        input_file=args.input_file,
        output_prefix=args.output_prefix,
        reduce_fold=args.reduce_fold,
        map_qual=args.map_qual,
    )

    try:
        alignment = SAM.ParseBAM(str(args.input_file))
        alignment.readsQual_boxplot(
            outfile=str(args.output_prefix),
            q_cut=args.map_qual,
            shrink=args.reduce_fold,
        )
    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(1, f"{parser.prog}: error: {exc}\n")

    if not args.skip_plot:
        run_plot_script(
            parser=parser,
            output_prefix=args.output_prefix,
            rscript_executable=args.rscript,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
