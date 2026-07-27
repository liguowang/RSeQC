#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

from __future__ import annotations
import argparse
import sys
from pathlib import Path
from typing import Sequence
from qcmodule import SAM

DESCRIPTION = (
    "Infer RNA-seq library layout and strandedness from a SAM/BAM file."
)

EPILOG = """
Examples
--------
Pair-end strand-specific
  1++,1--,2+-,2-+ (dUTP))
     * read1 mapped to '+' strand indicates parental gene on '+' strand
     * read1 mapped to '-' strand indicates parental gene on '-' strand
     * read2 mapped to '+' strand indicates parental gene on '-' strand
     * read2 mapped to '-' strand indicates parental gene on '+' strand
     
  1+-,1-+,2++,2-- (Illumina TruSeq):
     * read1 mapped to '+' strand indicates parental gene on '-' strand
     * read1 mapped to '-' strand indicates parental gene on '+' strand
     * read2 mapped to '+' strand indicates parental gene on '+' strand
     * read2 mapped to '-' strand indicates parental gene on '-' strand
     
Single-end:
  ++,--     (same strand)
     * read mapped to '+' strand indicates parental gene on '+' strand
     * read mapped to '-' strand indicates parental gene on '-' strand
  +-,-+     (opposite strand)
     * read mapped to '+' strand indicates parental gene on '-' strand
     * read mapped to '-' strand indicates parental gene on '+' strand
Notes
-----
* A sample size of at least 200,000 alignments is recommended.
* Reads should be aligned as non-strand-specific before running this tool.
"""

__author__ = "Liguo Wang"
__version__ = "5.05"


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
        help="Input alignment file in SAM or BAM format.",
    )
    parser.add_argument(
        "-r",
        "--refgene",
        dest="refgene_bed",
        required=True,
        type=Path,
        help="Reference gene model in BED format.",
    )
    parser.add_argument(
        "-s",
        "--sample-size",
        dest="sample_size",
        type=int,
        default=200_000,
        metavar="INT",
        help=(
            "Number of alignments sampled from the SAM/BAM file. "
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
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def validate_args(
    parser: argparse.ArgumentParser,
    input_file: Path,
    refgene_bed: Path,
    sample_size: int,
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

    if not refgene_bed.is_file():
        parser.error(
            f"reference BED file does not exist or is not a file: "
            f"{refgene_bed}"
        )

    if sample_size <= 0:
        parser.error("--sample-size must be greater than zero")

    if map_qual < 0:
        parser.error("--mapq must be zero or greater")

    if sample_size < 1_000:
        print(
            "Warning: sample size is below 1,000; the inferred protocol "
            "may be unreliable.",
            file=sys.stderr,
        )


def print_results(
    protocol: str,
    strand_pattern_1: float,
    strand_pattern_2: float,
    undetermined: float,
) -> None:
    """Print inferred experiment design in the original output format."""
    undetermined = max(float(undetermined), 0.0)

    if protocol == "PairEnd":
        print("\nThis is PairEnd Data")
        print(f'Fraction of reads failed to determine: {undetermined:.4f}')
        print(
            'Fraction of reads explained by '
            '"1++,1--,2+-,2-+": '
            f"{strand_pattern_1:.4f}"
        )
        print(
            'Fraction of reads explained by '
            '"1+-,1-+,2++,2--": '
            f"{strand_pattern_2:.4f}"
        )
        return

    if protocol == "SingleEnd":
        print("\nThis is SingleEnd Data")
        print(f'Fraction of reads failed to determine: {undetermined:.4f}')
        print(
            'Fraction of reads explained by "++,--": '
            f"{strand_pattern_1:.4f}"
        )
        print(
            'Fraction of reads explained by "+-,-+": '
            f"{strand_pattern_2:.4f}"
        )
        return

    print(f"Unknown data type: {protocol}")


def main(argv: Sequence[str] | None = None) -> int:
    """Infer and report RNA-seq experiment design."""
    parser = build_parser()
    args = parser.parse_args(argv)

    validate_args(
        parser=parser,
        input_file=args.input_file,
        refgene_bed=args.refgene_bed,
        sample_size=args.sample_size,
        map_qual=args.map_qual,
    )

    try:
        alignment = SAM.ParseBAM(str(args.input_file))
        protocol, strand_pattern_1, strand_pattern_2, undetermined = (
            alignment.configure_experiment(
                refbed=str(args.refgene_bed),
                sample_size=args.sample_size,
                q_cut=args.map_qual,
            )
        )
    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(1, f"{parser.prog}: error: {exc}\n")

    print_results(
        protocol=protocol,
        strand_pattern_1=strand_pattern_1,
        strand_pattern_2=strand_pattern_2,
        undetermined=undetermined,
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
