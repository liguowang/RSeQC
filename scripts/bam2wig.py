#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Convert a sorted, indexed BAM file into WIG coverage files.

SAM input is not supported.
"""

from __future__ import annotations
import argparse
import sys
from pathlib import Path
from typing import Sequence
from qcmodule import SAM
from importlib.metadata import version

__version__ = version("RSeQC")
__author__ = "Liguo Wang"

def build_parser() -> argparse.ArgumentParser:
    """Create and configure the command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "-i",
        "--input-file",
        required=True,
        type=Path,
        help=(
            "Coordinate-sorted BAM file. The BAM index must be available as "
            "<file>.bai or alongside the BAM as <stem>.bai."
        ),
    )
    parser.add_argument(
        "-s",
        "--chrom-size",
        "--chromSize",
        dest="chrom_size",
        required=True,
        type=Path,
        help=(
            "Two-column chromosome-size file containing chromosome name and "
            "length. Names must match those in the BAM file."
        ),
    )
    parser.add_argument(
        "-o",
        "--out-prefix",
        required=True,
        type=Path,
        help=(
            "Output prefix. Unstranded data produce one WIG file; stranded "
            "data produce forward and reverse WIG files."
        ),
    )
    parser.add_argument(
        "-t",
        "--wigsum",
        dest="total_wigsum",
        type=float,
        default=None,
        metavar="FLOAT",
        help=(
            "Target total WIG sum used for normalization. For example, "
            "1,000,000,000 corresponds to the coverage from ten million "
            "100-nt reads. Omit to disable normalization."
        ),
    )
    parser.add_argument(
        "-u",
        "--skip-multi-hits",
        dest="skip_multi",
        action="store_true",
        help="Exclude non-unique alignments.",
    )
    parser.add_argument(
        "-d",
        "--strand",
        dest="strand_rule",
        default=None,
        help=(
            "Strand rule, for example '1++,1--,2+-,2-+'. Use "
            "infer_experiment.py if the library protocol is unknown."
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
            "Minimum mapping quality used to identify uniquely mapped reads. "
            "Default: %(default)s"
        ),
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def load_chrom_sizes(path: Path) -> dict[str, int]:
    """Read a two-column chromosome-size file."""
    chrom_sizes: dict[str, int] = {}

    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()

            if not stripped or stripped.startswith("#"):
                continue

            fields = stripped.split()
            if len(fields) < 2:
                raise ValueError(
                    f"{path}:{line_number}: expected at least two columns"
                )

            chromosome = fields[0]

            try:
                size = int(fields[1])
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{line_number}: chromosome size must be an integer"
                ) from exc

            if size <= 0:
                raise ValueError(
                    f"{path}:{line_number}: chromosome size must be positive"
                )

            if chromosome in chrom_sizes:
                raise ValueError(
                    f"{path}:{line_number}: duplicate chromosome {chromosome!r}"
                )

            chrom_sizes[chromosome] = size

    if not chrom_sizes:
        raise ValueError(f"no chromosome sizes were found in {path}")

    return chrom_sizes


def find_bam_index(bam_path: Path) -> Path | None:
    """Return an existing BAM index path, if present."""
    candidates = (
        Path(f"{bam_path}.bai"),
        bam_path.with_suffix(".bai"),
    )

    for candidate in candidates:
        if candidate.is_file():
            return candidate

    return None


def validate_args(
    parser: argparse.ArgumentParser,
    bam_path: Path,
    chrom_size_path: Path,
    map_qual: int,
    total_wigsum: float | None,
) -> None:
    """Validate command-line arguments."""
    if not bam_path.is_file():
        parser.error(f"input BAM file does not exist or is not a file: {bam_path}")

    if bam_path.suffix.lower() != ".bam":
        parser.error(f"input must be a BAM file: {bam_path}")

    if not chrom_size_path.is_file():
        parser.error(
            "chromosome-size file does not exist or is not a file: "
            f"{chrom_size_path}"
        )

    if find_bam_index(bam_path) is None:
        parser.error(
            "BAM index not found; expected either "
            f"{bam_path}.bai or {bam_path.with_suffix('.bai')}"
        )

    if map_qual < 0:
        parser.error("--mapq must be zero or greater")

    if total_wigsum is not None and total_wigsum <= 0:
        parser.error("--wigsum must be greater than zero")


def calculate_normalization_factor(
    bam_path: Path,
    chrom_sizes: dict[str, int],
    target_wigsum: float | None,
    skip_multi: bool,
) -> float | None:
    """Calculate the WIG normalization factor, when requested."""
    if target_wigsum is None:
        return None

    alignment = SAM.ParseBAM(str(bam_path))
    wig_sum = alignment.calWigSum(
        chrom_sizes=chrom_sizes,
        skip_multi=skip_multi,
    )

    print(f"Total WIG sum: {wig_sum}", file=sys.stderr)

    if wig_sum is None or wig_sum <= 0:
        raise ValueError(
            "normalization cannot be calculated because the observed "
            f"WIG sum is {wig_sum!r}"
        )

    return target_wigsum / wig_sum


def main(argv: Sequence[str] | None = None) -> int:
    """Run BAM-to-WIG conversion."""
    parser = build_parser()
    args = parser.parse_args(argv)

    validate_args(
        parser,
        args.input_file,
        args.chrom_size,
        args.map_qual,
        args.total_wigsum,
    )

    print(f"Skip multi-hits: {args.skip_multi}", file=sys.stderr)

    try:
        chrom_sizes = load_chrom_sizes(args.chrom_size)

        normalization_factor = calculate_normalization_factor(
            bam_path=args.input_file,
            chrom_sizes=chrom_sizes,
            target_wigsum=args.total_wigsum,
            skip_multi=args.skip_multi,
        )

        alignment = SAM.ParseBAM(str(args.input_file))
        alignment.bamTowig(
            outfile=str(args.out_prefix),
            chrom_sizes=chrom_sizes,
            chrom_file=str(args.chrom_size),
            q_cut=args.map_qual,
            skip_multi=args.skip_multi,
            strand_rule=args.strand_rule,
            WigSumFactor=normalization_factor,
        )

    except (OSError, ValueError, RuntimeError, ZeroDivisionError) as exc:
        parser.exit(1, f"{parser.prog}: error: {exc}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
