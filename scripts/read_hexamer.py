#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Calculate normalized hexamer frequencies from FASTA or FASTQ files."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence
from qcmodule import FrameKmer
from importlib.metadata import version

__version__ = version("RSeQC")
__author__ = "Liguo Wang"


DESCRIPTION = (
    "Calculate normalized hexamer frequencies for one or more read files, "
    "with optional reference-genome and reference-transcriptome comparisons."
)

EPILOG = """
Example
-------
read_hexamer.py \
    -i reads_1.fastq,reads_2.fastq \
    -r genome.fa \
    -g transcripts.fa

Notes
-----
* Input sequence files may be in FASTA or FASTQ format.
* Multiple read files are supplied as a comma-separated list.
* Hexamers containing N are omitted from the output.
* Output is written to standard output.
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
        dest="input_reads",
        required=True,
        help=(
            "Comma-separated FASTA/FASTQ read files, for example "
            "'reads_1.fq,reads_2.fa'."
        ),
    )
    parser.add_argument(
        "-r",
        "--refgenome",
        dest="ref_genome",
        type=Path,
        help="Optional reference-genome FASTA file.",
    )
    parser.add_argument(
        "-g",
        "--refgene",
        dest="ref_gene",
        type=Path,
        help="Optional reference mRNA/transcript FASTA file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Write the frequency table to this file instead of standard output.",
    )
    parser.add_argument(
        "--skip-missing",
        action="store_true",
        help=(
            "Skip missing read files with a warning. By default, any missing "
            "input file is treated as an error."
        ),
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def parse_read_files(
    parser: argparse.ArgumentParser,
    value: str,
    *,
    skip_missing: bool,
) -> list[Path]:
    """Parse and validate the comma-separated read-file list."""
    files = [Path(item.strip()) for item in value.split(",") if item.strip()]

    if not files:
        parser.error("--input must contain at least one file")

    valid_files: list[Path] = []

    for path in files:
        if path.is_file():
            valid_files.append(path)
            continue

        if skip_missing:
            print(
                f"Warning: input file does not exist and will be skipped: {path}",
                file=sys.stderr,
            )
        else:
            parser.error(f"input file does not exist: {path}")

    if not valid_files:
        parser.error("none of the requested input read files exists")

    return valid_files


def validate_optional_file(
    parser: argparse.ArgumentParser,
    path: Path | None,
    label: str,
) -> None:
    """Validate an optional reference file."""
    if path is not None and not path.is_file():
        parser.error(f"{label} file does not exist: {path}")


def unique_display_name(path: Path, existing: set[str]) -> str:
    """Create a stable, unique column name for an input file."""
    candidate = path.name

    if candidate not in existing:
        return candidate

    candidate = str(path)

    if candidate not in existing:
        return candidate

    suffix = 2
    while f"{candidate}#{suffix}" in existing:
        suffix += 1

    return f"{candidate}#{suffix}"


def calculate_hexamer_frequencies(
    path: Path,
) -> tuple[dict[str, float], float]:
    """Calculate raw hexamer counts and their total."""
    print(
        f"Calculate hexamer frequencies for {path} ...",
        end=" ",
        file=sys.stderr,
    )

    counts = FrameKmer.kmer_freq_file(
        fastafile=str(path),
        word_size=6,
        step_size=1,
        frame=0,
    )
    total = float(sum(counts.values()))

    print("Done", file=sys.stderr)
    return counts, total


def collect_inputs(
    read_files: Sequence[Path],
    ref_genome: Path | None,
    ref_gene: Path | None,
) -> tuple[list[str], dict[str, dict[str, float]], dict[str, float]]:
    """Calculate hexamer tables for all requested inputs."""
    names: list[str] = []
    tables: dict[str, dict[str, float]] = {}
    totals: dict[str, float] = {}
    existing_names: set[str] = set()

    ordered_paths = list(read_files)

    if ref_genome is not None:
        ordered_paths.append(ref_genome)

    if ref_gene is not None:
        ordered_paths.append(ref_gene)

    for path in ordered_paths:
        name = unique_display_name(path, existing_names)
        existing_names.add(name)

        counts, total = calculate_hexamer_frequencies(path)

        names.append(name)
        tables[name] = counts
        totals[name] = total

    return names, tables, totals


def normalized_frequency(
    counts: dict[str, float],
    total: float,
    kmer: str,
) -> float:
    """Return a normalized k-mer frequency without masking zero totals."""
    if total <= 0:
        return 0.0

    return float(counts.get(kmer, 0.0)) / total


def write_report(
    output,
    names: Sequence[str],
    tables: dict[str, dict[str, float]],
    totals: dict[str, float],
) -> None:
    """Write the normalized hexamer-frequency table."""
    print("Hexamer\t" + "\t".join(names), file=output)

    for kmer in FrameKmer.all_possible_kmer(6):
        if "N" in kmer:
            continue

        values = [
            normalized_frequency(tables[name], totals[name], kmer)
            for name in names
        ]
        formatted_values = "\t".join(f"{value:.12g}" for value in values)
        print(f"{kmer}\t{formatted_values}", file=output)


def main(argv: Sequence[str] | None = None) -> int:
    """Run hexamer-frequency analysis."""
    parser = build_parser()
    args = parser.parse_args(argv)

    read_files = parse_read_files(
        parser,
        args.input_reads,
        skip_missing=args.skip_missing,
    )
    validate_optional_file(parser, args.ref_genome, "reference genome")
    validate_optional_file(parser, args.ref_gene, "reference transcript")

    try:
        names, tables, totals = collect_inputs(
            read_files=read_files,
            ref_genome=args.ref_genome,
            ref_gene=args.ref_gene,
        )

        if args.output is None:
            write_report(sys.stdout, names, tables, totals)
        else:
            output_parent = args.output.parent

            if not output_parent.exists():
                parser.error(
                    f"output directory does not exist: {output_parent}"
                )

            if not output_parent.is_dir():
                parser.error(
                    f"output parent is not a directory: {output_parent}"
                )

            with args.output.open("w", encoding="utf-8") as output:
                write_report(output, names, tables, totals)

            print(f"Created: {args.output}", file=sys.stderr)

    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(1, f"{parser.prog}: error: {exc}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())

