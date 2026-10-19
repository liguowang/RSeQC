#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Split a paired-end BAM into read-1, read-2, and unmapped BAM files.

Mapped records are rewritten as independent single-end alignments. Unmapped
records are written unchanged to ``PREFIX.unmap.bam``.

The historical flag-conversion behavior is preserved for mapped reads:
paired-end and mate-related flags are removed, while reverse-strand,
secondary, QC-failure, and duplicate flags are retained.
"""

from __future__ import annotations

import argparse
import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence
import pysam
from importlib.metadata import version

__version__ = version("RSeQC")
__author__ = "Liguo Wang"


DESCRIPTION = (
    "Split a paired-end BAM into separate single-end BAM files for read 1, "
    "read 2, and unmapped records."
)

EPILOG = """
Outputs
-------
PREFIX.R1.bam
    Mapped records marked as read 1.

PREFIX.R2.bam
    Other mapped records, historically including read-2 records.

PREFIX.unmap.bam
    Unmapped records written without modification.

Example
-------
split_paired_bam.py \
    -i sample.bam \
    -o sample_split \
    --index-output

Notes
-----
* Mapped records are rewritten as unpaired single-end alignments.
* Reference position, CIGAR, MAPQ, sequence, qualities, and auxiliary tags
  are retained.
* Mate-reference, mate-position, template-length, and paired-end flags are
  intentionally removed from mapped output records.
"""


@dataclass
class SplitCounts:
    """Record counts for the three output BAM files."""

    total: int = 0
    read1: int = 0
    read2: int = 0
    unmapped: int = 0


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
        help="Input paired-end BAM file.",
    )
    parser.add_argument(
        "-o",
        "--out-prefix",
        dest="output_prefix",
        required=True,
        type=Path,
        help="Prefix for the R1, R2, and unmapped BAM files.",
    )
    parser.add_argument(
        "--index-output",
        action="store_true",
        help="Create BAI indexes for the output BAM files.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow existing output BAM and index files to be replaced.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable detailed progress logging.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def configure_logging(verbose: bool) -> None:
    """Configure command-line logging."""
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.DEBUG if verbose else logging.INFO,
    )


def output_paths(prefix: Path) -> tuple[Path, Path, Path]:
    """Return the R1, R2, and unmapped output BAM paths."""
    return (
        Path(f"{prefix}.R1.bam"),
        Path(f"{prefix}.R2.bam"),
        Path(f"{prefix}.unmap.bam"),
    )


def possible_index_paths(bam_path: Path) -> tuple[Path, Path]:
    """Return common BAI paths for a BAM file."""
    return (
        Path(f"{bam_path}.bai"),
        bam_path.with_suffix(".bai"),
    )


def validate_args(
    parser: argparse.ArgumentParser,
    args: argparse.Namespace,
) -> tuple[Path, Path, Path]:
    """Validate command-line arguments."""
    if not args.input_file.is_file():
        parser.error(f"input BAM file does not exist: {args.input_file}")

    if args.input_file.suffix.lower() != ".bam":
        parser.error(f"input alignment must be a BAM file: {args.input_file}")

    output_parent = args.output_prefix.parent

    if not output_parent.exists():
        parser.error(f"output directory does not exist: {output_parent}")

    if not output_parent.is_dir():
        parser.error(f"output parent is not a directory: {output_parent}")

    paths = output_paths(args.output_prefix)

    if not args.overwrite:
        candidates = list(paths)

        for path in paths:
            candidates.extend(possible_index_paths(path))

        existing = [path for path in candidates if path.exists()]

        if existing:
            parser.error(
                "output file already exists; use --overwrite to replace: "
                + ", ".join(str(path) for path in existing)
            )

    return paths


def single_end_flag(alignment: pysam.AlignedSegment) -> int:
    """Build the historical single-end flag for a mapped record."""
    flag = 0

    if alignment.is_reverse:
        flag |= 0x0010

    if alignment.is_secondary:
        flag |= 0x0100

    if alignment.is_qcfail:
        flag |= 0x0200

    if alignment.is_duplicate:
        flag |= 0x0400

    return flag


def to_single_end_alignment(
    old_alignment: pysam.AlignedSegment,
    header: pysam.AlignmentHeader,
) -> pysam.AlignedSegment:
    """Copy a mapped alignment into a new unpaired record."""
    new_alignment = pysam.AlignedSegment(header)

    new_alignment.query_name = old_alignment.query_name
    new_alignment.flag = single_end_flag(old_alignment)
    new_alignment.reference_id = old_alignment.reference_id
    new_alignment.reference_start = old_alignment.reference_start
    new_alignment.mapping_quality = old_alignment.mapping_quality
    new_alignment.cigartuples = old_alignment.cigartuples
    new_alignment.query_sequence = old_alignment.query_sequence

    if old_alignment.query_qualities is not None:
        new_alignment.query_qualities = old_alignment.query_qualities

    new_alignment.set_tags(old_alignment.get_tags(with_value_type=True))

    return new_alignment


def split_paired_bam(
    input_file: Path,
    read1_bam: Path,
    read2_bam: Path,
    unmapped_bam: Path,
) -> SplitCounts:
    """Split paired-end records into three BAM files."""
    counts = SplitCounts()

    logging.info("Splitting %s", input_file)

    with pysam.AlignmentFile(str(input_file), "rb") as input_bam:
        with (
            pysam.AlignmentFile(
                str(read1_bam),
                "wb",
                template=input_bam,
            ) as read1_output,
            pysam.AlignmentFile(
                str(read2_bam),
                "wb",
                template=input_bam,
            ) as read2_output,
            pysam.AlignmentFile(
                str(unmapped_bam),
                "wb",
                template=input_bam,
            ) as unmapped_output,
        ):
            for old_alignment in input_bam:
                counts.total += 1

                if old_alignment.is_unmapped:
                    unmapped_output.write(old_alignment)
                    counts.unmapped += 1
                    continue

                new_alignment = to_single_end_alignment(
                    old_alignment,
                    input_bam.header,
                )

                if old_alignment.is_read1:
                    read1_output.write(new_alignment)
                    counts.read1 += 1
                else:
                    # Preserve the original behavior: any mapped record not
                    # marked read 1 is written to the R2 output.
                    read2_output.write(new_alignment)
                    counts.read2 += 1

    return counts


def index_bam(path: Path) -> None:
    """Create a BAI index for one BAM file."""
    logging.info("Indexing %s", path)

    try:
        pysam.index("-f", str(path))
    except pysam.utils.SamtoolsError as exc:
        raise RuntimeError(
            f"could not index {path}; the output BAM may not be "
            f"coordinate-sorted: {exc}"
        ) from exc


def print_report(
    read1_bam: Path,
    read2_bam: Path,
    unmapped_bam: Path,
    counts: SplitCounts,
) -> None:
    """Print a concise split summary."""
    print(f"{'Total records:':<55}{counts.total}")
    print(f"{str(read1_bam) + ' (Read 1):':<55}{counts.read1}")
    print(f"{str(read2_bam) + ' (Read 2):':<55}{counts.read2}")
    print(f"{str(unmapped_bam) + ' (Unmapped):':<55}{counts.unmapped}")


def main(argv: Sequence[str] | None = None) -> int:
    """Run paired-BAM splitting."""
    parser = build_parser()
    args = parser.parse_args(argv)

    configure_logging(args.verbose)
    read1_bam, read2_bam, unmapped_bam = validate_args(parser, args)

    try:
        counts = split_paired_bam(
            input_file=args.input_file,
            read1_bam=read1_bam,
            read2_bam=read2_bam,
            unmapped_bam=unmapped_bam,
        )

        if counts.total != counts.read1 + counts.read2 + counts.unmapped:
            raise RuntimeError("internal count mismatch after BAM splitting")

        for path in (read1_bam, read2_bam, unmapped_bam):
            if not path.is_file():
                raise OSError(f"expected output BAM was not created: {path}")

        if args.index_output:
            for path in (read1_bam, read2_bam, unmapped_bam):
                index_bam(path)

    except (
        OSError,
        ValueError,
        RuntimeError,
        KeyError,
        IndexError,
        pysam.utils.SamtoolsError,
    ) as exc:
        logging.error("%s", exc)
        return 1

    print_report(
        read1_bam=read1_bam,
        read2_bam=read2_bam,
        unmapped_bam=unmapped_bam,
        counts=counts,
    )
    logging.info("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

