#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Assess whether splice-junction discovery has reached sequencing saturation."""

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
    "Estimate splice-junction saturation by repeatedly subsampling alignments "
    "and measuring the number of discovered junctions."
)

EPILOG = """
Interpretation
--------------
As sequencing depth approaches saturation, progressively fewer new splice
junctions should be discovered.

Example
-------
junction_saturation.py -i sample.bam -r genes.bed12 -o sample -l 5 -u 100 \\
    -s 5 -m 50 -v 1 -q 30
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
        dest="refgene_bed",
        required=True,
        type=Path,
        help=(
            "Reference gene model in BED format, used to classify known "
            "splice junctions."
        ),
    )
    parser.add_argument(
        "-l",
        "--percentile-floor",
        dest="percentile_low_bound",
        type=int,
        default=5,
        metavar="INT",
        help="Starting sampling percentile. Default: %(default)s",
    )
    parser.add_argument(
        "-u",
        "--percentile-ceiling",
        dest="percentile_up_bound",
        type=int,
        default=100,
        metavar="INT",
        help="Ending sampling percentile. Default: %(default)s",
    )
    parser.add_argument(
        "-s",
        "--percentile-step",
        dest="percentile_step",
        type=int,
        default=5,
        metavar="INT",
        help=(
            "Sampling percentile step. Smaller values perform more "
            "subsampling iterations. Default: %(default)s"
        ),
    )
    parser.add_argument(
        "-m",
        "--min-intron",
        dest="minimum_intron_size",
        type=int,
        default=50,
        metavar="INT",
        help="Minimum intron length in bp. Default: %(default)s",
    )
    parser.add_argument(
        "-v",
        "--min-coverage",
        dest="minimum_splice_read",
        type=int,
        default=1,
        metavar="INT",
        help=(
            "Minimum number of supporting reads required to call a junction. "
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
        help="Generate saturation data but do not execute the R plotting script.",
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
    refgene_bed: Path,
    output_prefix: Path,
    percentile_low_bound: int,
    percentile_up_bound: int,
    percentile_step: int,
    minimum_intron_size: int,
    minimum_splice_read: int,
    map_qual: int,
) -> None:
    """Validate command-line arguments."""
    if not input_file.is_file():
        parser.error(
            f"input alignment file does not exist or is not a file: {input_file}"
        )

    if input_file.suffix.lower() not in {".bam", ".sam", ".cram"}:
        parser.error(
            "input alignment file must have a .bam, .sam, or .cram "
            f"extension: {input_file}"
        )

    if not refgene_bed.is_file():
        parser.error(
            f"reference BED file does not exist or is not a file: "
            f"{refgene_bed}"
        )

    if not 0 <= percentile_low_bound <= 100:
        parser.error("--percentile-floor must be between 0 and 100")

    if not 0 <= percentile_up_bound <= 100:
        parser.error("--percentile-ceiling must be between 0 and 100")

    if percentile_up_bound < percentile_low_bound:
        parser.error(
            "--percentile-ceiling must be greater than or equal to "
            "--percentile-floor"
        )

    if percentile_step <= 0:
        parser.error("--percentile-step must be greater than zero")

    if percentile_step > percentile_up_bound:
        parser.error(
            "--percentile-step cannot be greater than "
            "--percentile-ceiling"
        )

    if minimum_intron_size <= 0:
        parser.error("--min-intron must be greater than zero")

    if minimum_splice_read <= 0:
        parser.error("--min-coverage must be greater than zero")

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
    script_path = Path(f"{output_prefix}.junctionSaturation_plot.r")

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
    """Run splice-junction saturation analysis."""
    parser = build_parser()
    args = parser.parse_args(argv)

    validate_args(
        parser=parser,
        input_file=args.input_file,
        refgene_bed=args.refgene_bed,
        output_prefix=args.out_prefix,
        percentile_low_bound=args.percentile_low_bound,
        percentile_up_bound=args.percentile_up_bound,
        percentile_step=args.percentile_step,
        minimum_intron_size=args.minimum_intron_size,
        minimum_splice_read=args.minimum_splice_read,
        map_qual=args.map_qual,
    )

    try:
        alignment = SAM.ParseBAM(str(args.input_file))
        alignment.saturation_junction(
            outfile=str(args.out_prefix),
            refgene=str(args.refgene_bed),
            sample_start=args.percentile_low_bound,
            sample_end=args.percentile_up_bound,
            sample_step=args.percentile_step,
            min_intron=args.minimum_intron_size,
            recur=args.minimum_splice_read,
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
