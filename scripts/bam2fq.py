#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Convert alignments in BAM or SAM format to FASTQ."""

from __future__ import annotations
import argparse
import gzip
import shutil
import sys
from pathlib import Path
from typing import Iterable, Sequence
from qcmodule import SAM

__author__ = "Liguo Wang"
__version__ = "5.0.5"

def build_parser() -> argparse.ArgumentParser:
    """Create and configure the command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-i",
        "--input-file",
        required=True,
        type=Path,
        help="Alignment file in BAM or SAM format.",
    )
    parser.add_argument(
        "-o",
        "--out-prefix",
        required=True,
        type=Path,
        help="Prefix for the output FASTQ file or files. It does not \
            synchronize mates or guarantee that the R1 and R2 FASTQ files \
            contain matching records in the same order.",
    )
    parser.add_argument(
        "-s",
        "--single-end",
        action="store_true",
        help="Treat the input as single-end sequencing data.",
    )
    parser.add_argument(
        "-c",
        "--compress",
        action="store_true",
        help="Compress output FASTQ files with gzip.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    return parser


def validate_input(parser: argparse.ArgumentParser, input_file: Path) -> None:
    """Validate the input alignment file."""
    if not input_file.is_file():
        parser.error(f"input file does not exist or is not a file: {input_file}")


def output_fastq_paths(prefix: Path, single_end: bool) -> list[Path]:
    """Return the FASTQ paths produced by SAM.ParseBAM.bam2fq."""
    if single_end:
        return [Path(f"{prefix}.fastq")]

    return [
        Path(f"{prefix}.R1.fastq"),
        Path(f"{prefix}.R2.fastq"),
    ]


def gzip_file(path: Path, *, remove_original: bool = True) -> Path:
    """Compress one file using Python's gzip implementation."""
    if not path.is_file():
        raise FileNotFoundError(f"expected output file was not created: {path}")

    compressed_path = Path(f"{path}.gz")

    with path.open("rb") as source, gzip.open(compressed_path, "wb") as destination:
        shutil.copyfileobj(source, destination)

    if remove_original:
        path.unlink()

    return compressed_path


def gzip_outputs(paths: Iterable[Path]) -> list[Path]:
    """Compress all generated FASTQ files."""
    compressed_paths = []

    for path in paths:
        print(f"Compressing {path} ...", file=sys.stderr)
        compressed_paths.append(gzip_file(path))

    return compressed_paths


def main(argv: Sequence[str] | None = None) -> int:
    """Run the BAM/SAM-to-FASTQ conversion."""
    parser = build_parser()
    args = parser.parse_args(argv)

    validate_input(parser, args.input_file)

    try:
        alignment = SAM.ParseBAM(str(args.input_file))
        alignment.bam2fq(
            prefix=str(args.out_prefix),
            paired=not args.single_end,
        )

        outputs = output_fastq_paths(args.out_prefix, args.single_end)

        if args.compress:
            gzip_outputs(outputs)
            print("Compression complete.", file=sys.stderr)

    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(1, f"{parser.prog}: error: {exc}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
