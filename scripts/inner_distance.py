#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Estimate the inner distance between paired-end RNA-seq reads."""

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
    "Estimate the inner distance (insert size) between paired-end RNA-seq "
    "reads using a BED gene model."
)

EPILOG = """
Inner-distance definition
-------------------------
fragment size = read-1 length + inner distance + read-2 length

The inner distance may be negative when paired reads overlap.

Example
-------
inner_distance.py -i sample.bam -r genes.bed12 -o sample -k 1000000 -l -250 \\
    -u 250 -s 5
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
        required=True,
        type=Path,
        help="Prefix for generated output files.",
    )
    parser.add_argument(
        "-r",
        "--refgene",
        dest="ref_gene",
        required=True,
        type=Path,
        help="Reference gene model in BED12 format.",
    )
    parser.add_argument(
        "-k",
        "--sample-size",
        dest="sample_size",
        type=int,
        default=1_000_000,
        metavar="INT",
        help=(
            "Number of read pairs sampled to estimate inner distance. "
            "Default: %(default)s"
        ),
    )
    parser.add_argument(
        "-l",
        "--lower-bound",
        dest="lower_bound_size",
        type=int,
        default=-250,
        metavar="INT",
        help=(
            "Lower histogram bound in bp. "
            "Default: %(default)s"
        ),
    )
    parser.add_argument(
        "-u",
        "--upper-bound",
        dest="upper_bound_size",
        type=int,
        default=250,
        metavar="INT",
        help=(
            "Upper histogram bound in bp. "
            "Default: %(default)s"
        ),
    )
    parser.add_argument(
        "-s",
        "--step",
        dest="step_size",
        type=int,
        default=5,
        metavar="INT",
        help=(
            "Histogram bin width in bp. "
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
            "Minimum mapping quality for an alignment to be considered "
            "uniquely mapped. Default: %(default)s"
        ),
    )
    parser.add_argument(
        "--skip-plot",
        action="store_true",
        help="Generate result files but do not execute the R plotting script.",
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
    ref_gene: Path,
    output_prefix: Path,
    sample_size: int,
    lower_bound: int,
    upper_bound: int,
    step_size: int,
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

    if not ref_gene.is_file():
        parser.error(
            f"reference BED12 file does not exist or is not a file: "
            f"{ref_gene}"
        )

    if sample_size <= 0:
        parser.error("--sample-size must be greater than zero")

    if step_size <= 0:
        parser.error("--step must be greater than zero")

    if lower_bound >= upper_bound:
        parser.error("--lower-bound must be smaller than --upper-bound")

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
    script_path = Path(f"{output_prefix}.inner_distance_plot.r")

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
    """Run inner-distance estimation."""
    parser = build_parser()
    args = parser.parse_args(argv)

    validate_args(
        parser=parser,
        input_file=args.input_file,
        ref_gene=args.ref_gene,
        output_prefix=args.out_prefix,
        sample_size=args.sample_size,
        lower_bound=args.lower_bound_size,
        upper_bound=args.upper_bound_size,
        step_size=args.step_size,
        map_qual=args.map_qual,
    )

    try:
        alignment = SAM.ParseBAM(str(args.input_file))
        alignment.mRNA_inner_distance(
            outfile=str(args.out_prefix),
            low_bound=args.lower_bound_size,
            up_bound=args.upper_bound_size,
            step=args.step_size,
            refbed=str(args.ref_gene),
            sample_size=args.sample_size,
            q_cut=args.map_qual,
        )
    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(1, f"{parser.prog}: error: {exc}\n")

    if not args.skip_plot:
        run_plot_script(
            parser=parser,
            output_prefix=args.out_prefix,
            rscript_executable=args.rscript,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
