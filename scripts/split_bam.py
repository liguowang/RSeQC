#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Split a BAM file according to exon regions from a BED gene list.

Three BAM files are produced:

* ``PREFIX.in.bam`` contains alignments consumed by the input gene list.
* ``PREFIX.ex.bam`` contains alignments not consumed by the gene list.
* ``PREFIX.junk.bam`` contains QC-failed or unmapped alignments.

The original read-start and mate-start classification algorithm is preserved.
"""

from __future__ import annotations

import argparse
import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence
import pysam
from bx.intervals.intersection import Intersecter, Interval
from qcmodule import BED
from importlib.metadata import version

__version__ = version("RSeQC")
__author__ = "Liguo Wang"


DESCRIPTION = (
    "Split alignments into exon-overlapping, non-overlapping, and junk BAM "
    "files using exon regions from a BED gene model."
)

EPILOG = """
Classification
--------------
in.bam
    A mapped alignment whose read start, or mapped mate start, overlaps an
    exon from the supplied BED gene model.

ex.bam
    A mapped alignment whose tested start positions do not overlap an exon.

junk.bam
    An unmapped or QC-failed alignment.

Example
-------
split_bam.py \
    -i sample.bam \
    -r genes.bed12 \
    -o sample_split \
    --index-output

Notes
-----
* The classification intentionally uses read-start and mate-start positions,
  matching the historical implementation.
* Secondary and duplicate alignments are not filtered.
* Output BAMs preserve the input BAM header.
"""


@dataclass
class SplitCounts:
    """Alignment counts for each output category."""

    total: int = 0
    consumed: int = 0
    excluded: int = 0
    junk: int = 0


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
        help="Input BAM file.",
    )
    parser.add_argument(
        "-r",
        "--genelist",
        dest="gene_list",
        required=True,
        type=Path,
        help="Reference gene list in BED format.",
    )
    parser.add_argument(
        "-o",
        "--out-prefix",
        dest="output_prefix",
        required=True,
        type=Path,
        help="Prefix for the three output BAM files.",
    )
    parser.add_argument(
        "--index-output",
        action="store_true",
        help="Create BAI indexes for the output BAM files after splitting.",
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
    """Return the three output BAM paths."""
    return (
        Path(f"{prefix}.in.bam"),
        Path(f"{prefix}.ex.bam"),
        Path(f"{prefix}.junk.bam"),
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

    if not args.gene_list.is_file():
        parser.error(f"gene-list BED file does not exist: {args.gene_list}")

    output_parent = args.output_prefix.parent
    if not output_parent.exists():
        parser.error(f"output directory does not exist: {output_parent}")

    if not output_parent.is_dir():
        parser.error(f"output parent is not a directory: {output_parent}")

    paths = output_paths(args.output_prefix)

    if not args.overwrite:
        existing = [
            path
            for path in (
                *paths,
                *(Path(f"{path}.bai") for path in paths),
                *(path.with_suffix(".bai") for path in paths),
            )
            if path.exists()
        ]

        if existing:
            parser.error(
                "output file already exists; use --overwrite to replace: "
                + ", ".join(str(path) for path in existing)
            )

    return paths


def build_interval_trees(exons) -> dict[str, Intersecter]:
    """Build chromosome-specific exon interval trees."""
    ranges: dict[str, Intersecter] = {}

    for exon in exons:
        chromosome = str(exon[0]).upper()
        start = int(exon[1])
        end = int(exon[2])

        tree = ranges.setdefault(chromosome, Intersecter())
        tree.add_interval(Interval(start, end))

    return ranges


def load_exon_ranges(gene_list: Path) -> dict[str, Intersecter]:
    """Read BED exons and build interval trees."""
    logging.info("Reading gene model %s", gene_list)

    bed = BED.ParseBED(str(gene_list))
    exons = bed.getExon()

    if not exons:
        raise ValueError(f"no exon intervals were found in {gene_list}")

    exon_ranges = build_interval_trees(exons)
    logging.info(
        "Loaded exon intervals for %d chromosome(s)",
        len(exon_ranges),
    )
    return exon_ranges


def overlaps_exon_start(
    exon_ranges: dict[str, Intersecter],
    chromosome: str,
    position: int,
) -> bool:
    """Return whether a one-base start interval overlaps an exon."""
    tree = exon_ranges.get(chromosome)
    if tree is None:
        return False

    return bool(tree.find(position, position + 1))


def classify_alignment(
    aligned_read: pysam.AlignedSegment,
    samfile: pysam.AlignmentFile,
    exon_ranges: dict[str, Intersecter],
) -> str:
    """Classify one alignment using the historical algorithm."""
    if aligned_read.is_qcfail or aligned_read.is_unmapped:
        return "junk"

    chromosome = samfile.get_reference_name(
        aligned_read.reference_id
    ).upper()
    read_start = aligned_read.reference_start

    if aligned_read.mate_is_unmapped:
        if overlaps_exon_start(exon_ranges, chromosome, read_start):
            return "in"
        return "ex"

    mate_start = aligned_read.next_reference_start

    if (
        overlaps_exon_start(exon_ranges, chromosome, read_start)
        or overlaps_exon_start(exon_ranges, chromosome, mate_start)
    ):
        return "in"

    return "ex"


def split_bam(
    input_file: Path,
    gene_list: Path,
    in_bam: Path,
    ex_bam: Path,
    junk_bam: Path,
) -> SplitCounts:
    """Split the BAM file and return output counts."""
    exon_ranges = load_exon_ranges(gene_list)
    counts = SplitCounts()

    logging.info("Splitting %s", input_file)

    with pysam.AlignmentFile(str(input_file), "rb") as samfile:
        with (
            pysam.AlignmentFile(
                str(in_bam),
                "wb",
                template=samfile,
            ) as consumed_output,
            pysam.AlignmentFile(
                str(ex_bam),
                "wb",
                template=samfile,
            ) as excluded_output,
            pysam.AlignmentFile(
                str(junk_bam),
                "wb",
                template=samfile,
            ) as junk_output,
        ):
            outputs = {
                "in": consumed_output,
                "ex": excluded_output,
                "junk": junk_output,
            }

            for aligned_read in samfile:
                counts.total += 1

                category = classify_alignment(
                    aligned_read,
                    samfile,
                    exon_ranges,
                )
                outputs[category].write(aligned_read)

                if category == "in":
                    counts.consumed += 1
                elif category == "ex":
                    counts.excluded += 1
                else:
                    counts.junk += 1

    return counts


def index_bam(path: Path) -> None:
    """Create a BAI index for one BAM file."""
    logging.info("Indexing %s", path)

    try:
        pysam.index("-f", str(path))
    except pysam.utils.SamtoolsError as exc:
        raise RuntimeError(
            f"could not index {path}; output BAM may not be "
            "coordinate-sorted: {exc}"
        ) from exc


def print_report(
    in_bam: Path,
    ex_bam: Path,
    junk_bam: Path,
    counts: SplitCounts,
) -> None:
    """Print the historical-style summary report."""
    print(f"{'Total records:':<55}{counts.total}")
    print(
        f"{str(in_bam) + ' (Alignments consumed by input gene list):':<55}"
        f"{counts.consumed}"
    )
    print(
        f"{str(ex_bam) + ' (Alignments not consumed by input gene list):':<55}"
        f"{counts.excluded}"
    )
    print(
        f"{str(junk_bam) + ' (QC-failed, unmapped reads):':<55}"
        f"{counts.junk}"
    )


def main(argv: Sequence[str] | None = None) -> int:
    """Run BAM splitting."""
    parser = build_parser()
    args = parser.parse_args(argv)

    configure_logging(args.verbose)
    in_bam, ex_bam, junk_bam = validate_args(parser, args)

    try:
        counts = split_bam(
            input_file=args.input_file,
            gene_list=args.gene_list,
            in_bam=in_bam,
            ex_bam=ex_bam,
            junk_bam=junk_bam,
        )

        if counts.total != (
            counts.consumed + counts.excluded + counts.junk
        ):
            raise RuntimeError(
                "internal count mismatch after BAM splitting"
            )

        for path in (in_bam, ex_bam, junk_bam):
            if not path.is_file():
                raise OSError(f"expected output BAM was not created: {path}")

        if args.index_output:
            for path in (in_bam, ex_bam, junk_bam):
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
        in_bam=in_bam,
        ex_bam=ex_bam,
        junk_bam=junk_bam,
        counts=counts,
    )
    logging.info("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
