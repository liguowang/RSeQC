#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Calculate fragment-size statistics for each transcript or gene.

For each BED12 record, report:
1. Number of fragments used
2. Mean fragment size
3. Median fragment size
4. Standard deviation of fragment size

The original fragment-length calculation and read-filtering logic are
preserved.
"""

from __future__ import annotations
import argparse
import sys
from pathlib import Path
from typing import Iterator, Sequence
import numpy as np
import pysam
from importlib.metadata import version

__version__ = version("RSeQC")
__author__ = "Liguo Wang"

DESCRIPTION = (
    "Calculate per-transcript RNA fragment-size statistics from paired-end "
    "alignments and a BED12 gene model."
)

EPILOG = """
Example
-------
RNA_fragment_size.py -i sample.bam -r genes.bed12 -q 30 -n 3

Notes
-----
* The BAM file must be coordinate-sorted and indexed.
* Only read 1 from paired-end alignments is used.
* QC-failed, duplicate, secondary, mate-unmapped, and low-MAPQ alignments
  are excluded.
* Fragment sizes are calculated in transcript space using the original
  exon-overlap method.
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
        dest="input_file",
        required=True,
        type=Path,
        help="Input coordinate-sorted BAM file.",
    )
    parser.add_argument(
        "-r",
        "--refgene",
        dest="refgene_bed",
        required=True,
        type=Path,
        help="Reference gene model in standard 12-column BED format.",
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
        "-n",
        "--frag-num",
        dest="fragment_num",
        type=int,
        default=3,
        metavar="INT",
        help=(
            "Minimum number of fragments required to report statistics. "
            "Default: %(default)s"
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Write the table to this file instead of standard output.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def find_bam_index(bam_path: Path) -> Path | None:
    """Return an existing BAM index path, supporting both common names."""
    candidates = (
        Path(f"{bam_path}.bai"),
        bam_path.with_suffix(".bai"),
    )

    for candidate in candidates:
        if candidate.is_file():
            return candidate

    return None


def validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Validate command-line arguments."""
    if not args.input_file.is_file():
        parser.error(f"input BAM file does not exist: {args.input_file}")

    if args.input_file.suffix.lower() != ".bam":
        parser.error(f"input alignment must be a BAM file: {args.input_file}")

    if find_bam_index(args.input_file) is None:
        parser.error(
            "cannot find BAM index; expected either "
            f"{args.input_file}.bai or {args.input_file.with_suffix('.bai')}"
        )

    if not args.refgene_bed.is_file():
        parser.error(f"reference BED file does not exist: {args.refgene_bed}")

    if args.map_qual < 0:
        parser.error("--mapq must be zero or greater")

    if args.fragment_num <= 0:
        parser.error("--frag-num must be greater than zero")

    if args.output is not None:
        output_parent = args.output.parent

        if not output_parent.exists():
            parser.error(f"output directory does not exist: {output_parent}")

        if not output_parent.is_dir():
            parser.error(f"output parent is not a directory: {output_parent}")


def overlap_length2(
    intervals1: Sequence[Sequence[int]],
    intervals2: Sequence[Sequence[int]],
) -> int:
    """Calculate pairwise inclusive overlap length.

    This preserves the original inclusive-coordinate implementation.
    """
    overlap = 0

    for first in intervals1:
        for second in intervals2:
            overlap += len(
                range(
                    max(first[0], second[0]),
                    min(first[-1], second[-1]) + 1,
                )
            )

    return overlap


def parse_bed12_record(
    line: str,
    line_number: int,
) -> tuple[
    str,
    int,
    int,
    str,
    list[list[int]],
]:
    """Parse one BED12 record and build inclusive exon ranges."""
    fields = line.split()

    if len(fields) < 12:
        raise ValueError(
            f"BED line {line_number} has {len(fields)} columns; expected 12"
        )

    chromosome = fields[0]
    transcript_start = int(fields[1])
    transcript_end = int(fields[2])
    gene_name = fields[3]

    block_count = int(fields[9])
    block_sizes = [
        int(value)
        for value in fields[10].rstrip(",").split(",")
        if value
    ]
    relative_starts = [
        int(value)
        for value in fields[11].rstrip(",").split(",")
        if value
    ]

    if len(block_sizes) != block_count:
        raise ValueError(
            f"BED line {line_number}: blockCount={block_count}, "
            f"but found {len(block_sizes)} block sizes"
        )

    if len(relative_starts) != block_count:
        raise ValueError(
            f"BED line {line_number}: blockCount={block_count}, "
            f"but found {len(relative_starts)} block starts"
        )

    exon_starts = [
        transcript_start + relative_start
        for relative_start in relative_starts
    ]
    exon_ends = [
        exon_start + block_size
        for exon_start, block_size in zip(exon_starts, block_sizes)
    ]

    exon_ranges = [
        [start + 1, end + 1]
        for start, end in zip(exon_starts, exon_ends)
    ]

    return (
        chromosome,
        transcript_start,
        transcript_end,
        gene_name,
        exon_ranges,
    )


def format_result(
    chromosome: str,
    transcript_start: int,
    transcript_end: int,
    gene_name: str,
    fragment_sizes: Sequence[int],
    minimum_fragments: int,
) -> str:
    """Format one output row."""
    count = len(fragment_sizes)

    if count < minimum_fragments:
        mean_value = 0
        median_value = 0
        std_value = 0
    else:
        mean_value = float(np.mean(fragment_sizes))
        median_value = float(np.median(fragment_sizes))
        std_value = float(np.std(fragment_sizes))

    return "\t".join(
        str(value)
        for value in (
            chromosome,
            transcript_start,
            transcript_end,
            gene_name,
            count,
            mean_value,
            median_value,
            std_value,
        )
    )


def fragment_size(
    bedfile: Path,
    samfile: pysam.AlignmentFile,
    qcut: int = 30,
    ncut: int = 5,
) -> Iterator[str]:
    """Calculate fragment-size statistics for each BED12 record."""
    with bedfile.open("r", encoding="utf-8") as bed_handle:
        for line_number, line in enumerate(bed_handle, start=1):
            if not line.strip():
                continue

            if line.startswith(("#", "track", "browser")):
                continue

            (
                chromosome,
                transcript_start,
                transcript_end,
                gene_name,
                exon_ranges,
            ) = parse_bed12_record(line, line_number)

            try:
                aligned_reads = samfile.fetch(
                    chromosome,
                    transcript_start,
                    transcript_end,
                )
            except (ValueError, KeyError):
                yield format_result(
                    chromosome,
                    transcript_start,
                    transcript_end,
                    gene_name,
                    [],
                    ncut,
                )
                continue

            fragment_sizes: list[int] = []

            for aligned_read in aligned_reads:
                if not aligned_read.is_paired:
                    continue

                if aligned_read.is_read2:
                    continue

                if aligned_read.mate_is_unmapped:
                    continue

                if aligned_read.is_qcfail:
                    continue

                if aligned_read.is_duplicate:
                    continue

                if aligned_read.is_secondary:
                    continue

                if aligned_read.mapping_quality < qcut:
                    continue

                read_start = aligned_read.reference_start
                mate_start = aligned_read.next_reference_start

                if read_start > mate_start:
                    read_start, mate_start = mate_start, read_start

                if read_start < transcript_start or mate_start > transcript_end:
                    continue

                read_length = aligned_read.query_length
                if read_length is None:
                    continue

                mapped_range = [[read_start + 1, mate_start]]
                fragment_length = (
                    overlap_length2(exon_ranges, mapped_range) + read_length
                )
                fragment_sizes.append(fragment_length)

            yield format_result(
                chromosome,
                transcript_start,
                transcript_end,
                gene_name,
                fragment_sizes,
                ncut,
            )


def write_report(
    output,
    rows: Iterator[str],
) -> None:
    """Write the fragment-size report."""
    print(
        "\t".join(
            (
                "chrom",
                "tx_start",
                "tx_end",
                "symbol",
                "frag_count",
                "frag_mean",
                "frag_median",
                "frag_std",
            )
        ),
        file=output,
    )

    for row in rows:
        print(row, file=output)


def main(argv: Sequence[str] | None = None) -> int:
    """Run RNA fragment-size analysis."""
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(parser, args)

    try:
        with pysam.AlignmentFile(str(args.input_file), "rb") as samfile:
            rows = fragment_size(
                bedfile=args.refgene_bed,
                samfile=samfile,
                qcut=args.map_qual,
                ncut=args.fragment_num,
            )

            if args.output is None:
                write_report(sys.stdout, rows)
            else:
                with args.output.open("w", encoding="utf-8") as output:
                    write_report(output, rows)

                print(f"Created: {args.output}", file=sys.stderr)

    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(1, f"{parser.prog}: error: {exc}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
