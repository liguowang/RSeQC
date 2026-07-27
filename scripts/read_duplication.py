#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Calculate sequence-based and mapping-based read duplication rates."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Sequence

from qcmodule import SAM


__author__ = "Liguo Wang"
__version__ = "5.05"


DESCRIPTION = (
    "Estimate read duplication using both sequence identity and alignment "
    "position."
)

EPILOG = """
Definitions
-----------
Sequence-based duplication
    Reads with identical sequences are considered duplicates.

Mapping-based duplication
    Reads mapped to the same genomic location are considered duplicates.

Example
-------
read_duplication.py \
    -i sample.bam \
    -o sample \
    -u 500 \
    -q 30
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
        "-u",
        "--up-limit",
        dest="upper_limit",
        type=int,
        default=500,
        metavar="INT",
        help=(
            "Upper occurrence limit used only for plotting. "
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
        help="Generate duplication data but do not execute the R plot script.",
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
    upper_limit: int,
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

    if upper_limit <= 0:
        parser.error("--up-limit must be greater than zero")

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
    script_path = Path(f"{output_prefix}.DupRate_plot.r")

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
    """Run read-duplication analysis."""
    parser = build_parser()
    args = parser.parse_args(argv)

    validate_args(
        parser=parser,
        input_file=args.input_file,
        output_prefix=args.output_prefix,
        upper_limit=args.upper_limit,
        map_qual=args.map_qual,
    )

    try:
        alignment = SAM.ParseBAM(str(args.input_file))
        alignment.readDupRate(
            outfile=str(args.output_prefix),
            up_bound=args.upper_limit,
            q_cut=args.map_qual,
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
