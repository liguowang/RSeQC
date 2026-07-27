#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Calculate the distribution of mismatches across aligned reads."""

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
    "Calculate the positional distribution of mismatches across aligned "
    "reads. The BAM file must contain MD tags."
)

EPILOG = """
Example
-------
mismatch_profile.py -i sample.bam -l 101 -o sample -n 1000000 -q 30

Notes
-----
* --read-align-length is usually the original read length.
* CIGAR strings such as 101M, 68M140N33M, and 53M1D48M all correspond
  to a read alignment length of 101.
* The BAM file must contain MD tags.
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
        "--input",
        dest="input_bam",
        required=True,
        type=Path,
        help="Input BAM file.",
    )
    parser.add_argument(
        "-l",
        "--read-align-length",
        dest="read_alignment_length",
        required=True,
        type=int,
        metavar="INT",
        help="Aligned read length, usually the original read length.",
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
        "-n",
        "--read-num",
        dest="read_number",
        type=int,
        default=1_000_000,
        metavar="INT",
        help=(
            "Maximum number of aligned reads with mismatches to sample. "
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
        help="Minimum mapping quality. Default: %(default)s",
    )
    parser.add_argument(
        "--skip-plot",
        action="store_true",
        help="Generate profile files but do not execute the R plotting script.",
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
    input_bam: Path,
    output_prefix: Path,
    read_alignment_length: int,
    read_number: int,
    map_qual: int,
) -> None:
    """Validate command-line arguments."""
    if not input_bam.is_file():
        parser.error(
            f"input BAM file does not exist or is not a file: {input_bam}"
        )

    if input_bam.suffix.lower() != ".bam":
        parser.error(f"input file must have a .bam extension: {input_bam}")

    if read_alignment_length <= 0:
        parser.error("--read-align-length must be greater than zero")

    if read_number <= 0:
        parser.error("--read-num must be greater than zero")

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
    script_path = Path(f"{output_prefix}.mismatch_profile.r")

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
    """Generate the mismatch profile and optional plot."""
    parser = build_parser()
    args = parser.parse_args(argv)

    validate_args(
        parser=parser,
        input_bam=args.input_bam,
        output_prefix=args.output_prefix,
        read_alignment_length=args.read_alignment_length,
        read_number=args.read_number,
        map_qual=args.map_qual,
    )

    try:
        alignment = SAM.ParseBAM(str(args.input_bam))
        alignment.mismatchProfile(
            read_length=args.read_alignment_length,
            read_num=args.read_number,
            q_cut=args.map_qual,
            outfile=str(args.output_prefix),
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
