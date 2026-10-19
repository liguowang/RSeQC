#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Randomly divide a BAM file into approximately equal subsets.

All records sharing the same query name are assigned to the same output file,
so paired mates remain together. Output files are named
``<prefix>_<index>.bam``.
"""

from __future__ import annotations
import argparse
import random
import sys
from contextlib import ExitStack
from pathlib import Path
from typing import Sequence
import pysam
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
        help="Input BAM file.",
    )
    parser.add_argument(
        "-n",
        "--subset-num",
        required=True,
        type=int,
        metavar="INT",
        help="Number of output BAM subsets.",
    )
    parser.add_argument(
        "-o",
        "--out-prefix",
        required=True,
        type=Path,
        help="Output prefix. Files are named <prefix>_<index>.bam.",
    )
    parser.add_argument(
        "-s",
        "--skip-unmap",
        action="store_true",
        help="Exclude unmapped alignments.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        metavar="INT",
        help="Random seed for reproducible assignment.",
    )
    parser.add_argument(
        "--index",
        action="store_true",
        help="Create a BAM index for each output file after writing.",
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
    subset_num: int,
    output_prefix: Path,
) -> None:
    """Validate command-line arguments."""
    if not input_file.is_file():
        parser.error(f"input BAM file does not exist or is not a file: {input_file}")

    if input_file.suffix.lower() != ".bam":
        parser.error(f"input file must have a .bam extension: {input_file}")

    if subset_num <= 0:
        parser.error("--subset-num must be greater than zero")

    output_parent = output_prefix.parent
    if not output_parent.exists():
        parser.error(f"output directory does not exist: {output_parent}")

    if not output_parent.is_dir():
        parser.error(f"output parent is not a directory: {output_parent}")


def output_paths(prefix: Path, subset_num: int) -> list[Path]:
    """Return the output BAM paths."""
    return [
        Path(f"{prefix}_{index}.bam")
        for index in range(subset_num)
    ]


def divide_bam(
    input_file: Path,
    paths: Sequence[Path],
    skip_unmapped: bool,
    seed: int | None,
) -> tuple[int, int, list[int]]:
    """Assign query-name groups to output BAM files.

    Every alignment sharing the same query name is written to the same output
    file. This keeps paired mates, secondary alignments, and supplementary
    alignments together.

    Returns
    -------
    retained_count
        Number of alignments written.
    skipped_count
        Number of unmapped alignments skipped.
    per_file_counts
        Number of alignments written to each output file.
    """
    rng = random.Random(seed)
    per_file_counts = [0] * len(paths)
    retained_count = 0
    skipped_count = 0
    query_assignments: dict[str, int] = {}

    with pysam.AlignmentFile(str(input_file), "rb") as source:
        with ExitStack() as stack:
            outputs = [
                stack.enter_context(
                    pysam.AlignmentFile(
                        str(path),
                        "wb",
                        template=source,
                    )
                )
                for path in paths
            ]

            print(f"Dividing {input_file} ...", file=sys.stderr)

            for alignment in source.fetch(until_eof=True):
                if skip_unmapped and alignment.is_unmapped:
                    skipped_count += 1
                    continue

                query_name = alignment.query_name
                if query_name is None:
                    raise ValueError(
                        "encountered an alignment without a query name"
                    )

                output_index = query_assignments.get(query_name)
                if output_index is None:
                    output_index = rng.randrange(len(outputs))
                    query_assignments[query_name] = output_index

                outputs[output_index].write(alignment)

                per_file_counts[output_index] += 1
                retained_count += 1

    print("Done", file=sys.stderr)
    return retained_count, skipped_count, per_file_counts


def create_indexes(paths: Sequence[Path]) -> None:
    """Create BAM indexes for output files."""
    for path in paths:
        try:
            pysam.index(str(path))
        except pysam.SamtoolsError as exc:
            raise RuntimeError(
                f"could not index {path}; output BAM may not be coordinate-sorted"
            ) from exc


def main(argv: Sequence[str] | None = None) -> int:
    """Run BAM subdivision."""
    parser = build_parser()
    args = parser.parse_args(argv)

    validate_args(
        parser=parser,
        input_file=args.input_file,
        subset_num=args.subset_num,
        output_prefix=args.out_prefix,
    )

    paths = output_paths(args.out_prefix, args.subset_num)

    existing = [path for path in paths if path.exists()]
    if existing:
        parser.error(
            "refusing to overwrite existing output file(s): "
            + ", ".join(str(path) for path in existing)
        )

    try:
        retained_count, skipped_count, counts = divide_bam(
            input_file=args.input_file,
            paths=paths,
            skip_unmapped=args.skip_unmap,
            seed=args.seed,
        )

        if args.index:
            create_indexes(paths)

    except (OSError, ValueError, RuntimeError, pysam.SamtoolsError) as exc:
        parser.exit(1, f"{parser.prog}: error: {exc}\n")

    for path, count in zip(paths, counts):
        print(f"{path}\t{count}")

    print(f"Total alignments written: {retained_count}", file=sys.stderr)

    if args.skip_unmap:
        print(f"Unmapped alignments skipped: {skipped_count}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())

