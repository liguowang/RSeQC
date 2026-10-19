#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Generate a DNA sequence logo from FASTA, FASTQ, or sequence-only input.

The command is useful for visualizing nucleotide composition across sample
barcodes, cellular barcodes, molecular barcodes, and other fixed-length DNA
sequences.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Sequence
from qcmodule import fastq
from importlib.metadata import version

__version__ = version("RSeQC")
__author__ = "Liguo Wang"


DESCRIPTION = (
    "Generate a nucleotide count matrix and DNA sequence logo from "
    "fixed-length FASTA or FASTQ sequences."
)

EPILOG = """
Examples
--------
Generate a PDF logo from FASTQ:
    sc_seqLogo.py \
        -i barcodes.fastq.gz \
        -o barcodes \
        --iformat fq \
        --oformat pdf

Generate an SVG logo from FASTA and remove sequences containing N:
    sc_seqLogo.py \
        -i barcodes.fa \
        -o barcodes \
        --iformat fa \
        --oformat svg \
        --exclude-N

Highlight positions 4 through 8:
    sc_seqLogo.py \
        -i barcodes.fastq \
        -o barcodes \
        --highlight-start 4 \
        --highlight-end 8

Notes
-----
* All input sequences must have the same length.
* Input may be plain text or compressed using gzip, compress, or bzip2.
* Highlight positions are zero-based.
"""


INPUT_FORMATS = ("fq", "fa")
OUTPUT_FORMATS = ("pdf", "png", "svg")
STACK_ORDERS = ("big_on_top", "small_on_top", "fixed")


def build_parser() -> argparse.ArgumentParser:
    """Create and configure the command-line parser."""
    parser = argparse.ArgumentParser(
        description=DESCRIPTION,
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--infile",
        dest="in_file",
        required=True,
        type=Path,
        help=(
            "Input FASTA or FASTQ file. Compressed input is supported by "
            "qcmodule.fastq."
        ),
    )
    parser.add_argument(
        "-o",
        "--outfile",
        dest="out_file",
        required=True,
        type=Path,
        help="Prefix for the count matrix and sequence-logo files.",
    )
    parser.add_argument(
        "--iformat",
        dest="in_format",
        type=str.lower,
        choices=INPUT_FORMATS,
        default="fq",
        help="Input format. Default: %(default)s",
    )
    parser.add_argument(
        "--oformat",
        dest="out_format",
        type=str.lower,
        choices=OUTPUT_FORMATS,
        default="pdf",
        help="Sequence-logo output format. Default: %(default)s",
    )
    parser.add_argument(
        "-n",
        "--nseq-limit",
        dest="max_seq",
        type=int,
        default=None,
        metavar="INT",
        help="Maximum number of sequences to process. Default: all",
    )
    parser.add_argument(
        "--font-name",
        default="sans",
        metavar="NAME",
        help=(
            "Font used for logo characters. Valid names may be listed with "
            "logomaker.list_font_names(). Default: %(default)s"
        ),
    )
    parser.add_argument(
        "--stack-order",
        choices=STACK_ORDERS,
        default="big_on_top",
        help="Vertical ordering of nucleotides within each stack. Default: %(default)s",
    )
    parser.add_argument(
        "--flip-below",
        action="store_true",
        help="Flip characters drawn below the x-axis upside down.",
    )
    parser.add_argument(
        "--shade-below",
        type=float,
        default=0.0,
        metavar="FLOAT",
        help=(
            "Shading applied to characters below the x-axis, from 0 to 1. "
            "Default: %(default)s"
        ),
    )
    parser.add_argument(
        "--fade-below",
        type=float,
        default=0.0,
        metavar="FLOAT",
        help=(
            "Fading applied to characters below the x-axis, from 0 to 1. "
            "Default: %(default)s"
        ),
    )
    parser.add_argument(
        "--exclude-N",
        "--excludeN",
        dest="exclude_n",
        action="store_true",
        help="Exclude sequences containing the ambiguous nucleotide N.",
    )
    parser.add_argument(
        "--highlight-start",
        dest="highlight_start",
        type=int,
        default=None,
        metavar="INT",
        help="Zero-based first highlighted logo position.",
    )
    parser.add_argument(
        "--highlight-end",
        dest="highlight_end",
        type=int,
        default=None,
        metavar="INT",
        help="Zero-based last highlighted logo position.",
    )
    parser.add_argument(
        "--step-size",
        type=int,
        default=10_000,
        metavar="INT",
        help=(
            "Number of sequences processed per matrix-update batch. "
            "Default: %(default)s"
        ),
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


def validate_args(
    parser: argparse.ArgumentParser,
    args: argparse.Namespace,
) -> None:
    """Validate command-line arguments."""
    if not args.in_file.is_file():
        parser.error(f"input sequence file does not exist: {args.in_file}")

    if args.max_seq is not None and args.max_seq <= 0:
        parser.error("--nseq-limit must be greater than zero")

    if args.step_size <= 0:
        parser.error("--step-size must be greater than zero")

    if not 0.0 <= args.shade_below <= 1.0:
        parser.error("--shade-below must be between 0 and 1")

    if not 0.0 <= args.fade_below <= 1.0:
        parser.error("--fade-below must be between 0 and 1")

    if not args.font_name.strip():
        parser.error("--font-name cannot be empty")

    if args.highlight_start is None and args.highlight_end is not None:
        parser.error(
            "--highlight-start is required when --highlight-end is supplied"
        )

    if args.highlight_start is not None and args.highlight_end is None:
        parser.error(
            "--highlight-end is required when --highlight-start is supplied"
        )

    if args.highlight_start is not None:
        if args.highlight_start < 0:
            parser.error("--highlight-start must be zero or greater")

        if args.highlight_end < 0:
            parser.error("--highlight-end must be zero or greater")

        if args.highlight_end < args.highlight_start:
            parser.error(
                "--highlight-end must be greater than or equal to "
                "--highlight-start"
            )

    output_parent = args.out_file.parent
    if not output_parent.exists():
        parser.error(f"output directory does not exist: {output_parent}")

    if not output_parent.is_dir():
        parser.error(f"output parent is not a directory: {output_parent}")


def sequence_iterator(input_file: Path, input_format: str):
    """Return the requested FASTA or FASTQ sequence iterator."""
    if input_format == "fq":
        return fastq.fastq_iter(str(input_file), mode="seq")

    return fastq.fasta_iter(str(input_file))


def expected_output_paths(
    output_prefix: Path,
    output_format: str,
) -> tuple[Path, Path]:
    """Return the expected count-matrix and logo paths."""
    count_matrix = Path(f"{output_prefix}.count_matrix.csv")
    logo_file = Path(f"{output_prefix}.logo.{output_format}")
    return count_matrix, logo_file


def main(argv: Sequence[str] | None = None) -> int:
    """Generate a sequence count matrix and DNA logo."""
    parser = build_parser()
    args = parser.parse_args(argv)

    configure_logging(args.verbose)
    validate_args(parser, args)

    count_matrix_path, logo_path = expected_output_paths(
        args.out_file,
        args.out_format,
    )

    logging.debug("Input file: %s", args.in_file)
    logging.debug("Input format: %s", args.in_format)
    logging.debug("Output prefix: %s", args.out_file)
    logging.debug("Output format: %s", args.out_format)

    try:
        file_iter = sequence_iterator(args.in_file, args.in_format)

        matrix = fastq.seq2countMat(
            file_iter,
            step_size=args.step_size,
            exclude_N=args.exclude_n,
            limit=args.max_seq,
        )

        if matrix is None or matrix.empty:
            raise ValueError(
                "no usable sequences were found; check the input format, "
                "sequence content, and --exclude-N setting"
            )

        sequence_length = len(matrix.index)

        if args.highlight_start is not None:
            if args.highlight_start >= sequence_length:
                raise ValueError(
                    f"--highlight-start {args.highlight_start} is outside "
                    f"the valid range 0..{sequence_length - 1}"
                )

            if args.highlight_end >= sequence_length:
                raise ValueError(
                    f"--highlight-end {args.highlight_end} is outside "
                    f"the valid range 0..{sequence_length - 1}"
                )

        matrix.to_csv(
            count_matrix_path,
            index=True,
            index_label="Index",
        )

        if not count_matrix_path.is_file():
            raise OSError(
                f"count matrix was not created: {count_matrix_path}"
            )

        fastq.make_logo(
            matrix,
            outfile=str(args.out_file),
            exclude_N=args.exclude_n,
            font_name=args.font_name,
            stack_order=args.stack_order,
            flip_below=args.flip_below,
            shade_below=args.shade_below,
            fade_below=args.fade_below,
            highlight_start=args.highlight_start,
            highlight_end=args.highlight_end,
            oformat=args.out_format,
        )

        if not logo_path.is_file():
            raise RuntimeError(
                f"sequence logo was not created: {logo_path}"
            )

    except (
        OSError,
        ValueError,
        RuntimeError,
        KeyError,
        IndexError,
    ) as exc:
        logging.error("%s", exc)
        return 1

    logging.info("Created %s", count_matrix_path)
    logging.info("Created %s", logo_path)
    logging.info("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
