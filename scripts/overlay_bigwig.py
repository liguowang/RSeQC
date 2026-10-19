#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Apply an arithmetic operation to two BigWig signal tracks."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Callable, Sequence
import numpy as np
import pyBigWig
from qcmodule import BED
from qcmodule import twoList
from importlib.metadata import version

__version__ = version("RSeQC")
__author__ = "Liguo Wang"



DESCRIPTION = (
    "Combine two BigWig signal tracks position by position and write the "
    "result as a variableStep WIG file."
)

EPILOG = """
Supported actions
-----------------
Add
    Add corresponding values.

Average
    Calculate the arithmetic mean.

Division
    Divide BigWig 1 by BigWig 2 using the implementation in qcmodule.twoList.

Max
    Select the larger value.

Min
    Select the smaller value.

Product
    Multiply corresponding values.

Subtract
    Subtract BigWig 2 from BigWig 1.

geometricMean
    Calculate the geometric mean.

Example
-------
overlay_bigwig.py -i sample1.bw -j sample2.bw -a Subtract -o difference.wig -c 100000
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
        "--bwfile1",
        dest="bigwig_file1",
        required=True,
        type=Path,
        help="First BigWig file.",
    )
    parser.add_argument(
        "-j",
        "--bwfile2",
        dest="bigwig_file2",
        required=True,
        type=Path,
        help=(
            "Second BigWig file. Both files should use the same reference "
            "genome."
        ),
    )
    parser.add_argument(
        "-a",
        "--action",
        required=True,
        choices=(
            "Add",
            "Average",
            "Division",
            "Max",
            "Min",
            "Product",
            "Subtract",
            "geometricMean",
        ),
        help="Arithmetic operation applied to corresponding signal values.",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output_wig",
        required=True,
        type=Path,
        help="Output variableStep WIG file.",
    )
    parser.add_argument(
        "-c",
        "--chunk",
        dest="chunk_size",
        type=int,
        default=100_000,
        metavar="INT",
        help=(
            "Chromosome chunk size in bp. Smaller chunks use less memory. "
            "Default: %(default)s"
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow an existing output file to be replaced.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Validate command-line arguments."""
    for path, label in (
        (args.bigwig_file1, "first BigWig"),
        (args.bigwig_file2, "second BigWig"),
    ):
        if not path.is_file():
            parser.error(f"{label} file does not exist: {path}")

    if args.chunk_size <= 0:
        parser.error("--chunk must be greater than zero")

    output_parent = args.output_wig.parent
    if not output_parent.exists():
        parser.error(f"output directory does not exist: {output_parent}")

    if not output_parent.is_dir():
        parser.error(f"output parent is not a directory: {output_parent}")

    if args.output_wig.exists() and not args.overwrite:
        parser.error(
            f"output file already exists: {args.output_wig}; "
            "use --overwrite to replace it"
        )


def open_bigwig(path: Path):
    """Open and validate a BigWig file."""
    bigwig = pyBigWig.open(str(path))

    if bigwig is None:
        raise OSError(f"could not open BigWig file: {path}")

    if not bigwig.isBigWig():
        bigwig.close()
        raise ValueError(f"file is not a BigWig: {path}")

    return bigwig


def resolve_action(name: str) -> Callable:
    """Resolve a supported qcmodule.twoList operation."""
    callback = getattr(twoList, name, None)

    if callback is None or not callable(callback):
        raise ValueError(f"unsupported overlay action: {name}")

    return callback


def combined_chromosome_sizes(
    first: dict[str, int],
    second: dict[str, int],
) -> dict[str, int]:
    """Combine chromosome headers while checking shared lengths."""
    combined: dict[str, int] = {}

    for chromosome, size in first.items():
        combined[chromosome] = int(size)

    for chromosome, size in second.items():
        size = int(size)

        if chromosome in combined and combined[chromosome] != size:
            raise ValueError(
                f"chromosome size mismatch for {chromosome!r}: "
                f"{combined[chromosome]} versus {size}"
            )

        combined[chromosome] = size

    return combined


def interval_has_signal(
    bigwig,
    chromosome: str,
    start: int,
    end: int,
) -> bool:
    """Return whether a BigWig contains signal in an interval."""
    try:
        statistics = bigwig.stats(chromosome, start, end)
    except RuntimeError:
        return False

    return bool(statistics) and statistics[0] is not None


def interval_values(
    bigwig,
    chromosome: str,
    start: int,
    end: int,
) -> np.ndarray:
    """Return interval values, using zeros when the chromosome is absent."""
    length = end - start

    try:
        values = bigwig.values(chromosome, start, end)
    except RuntimeError:
        return np.zeros(length, dtype=float)

    if values is None:
        return np.zeros(length, dtype=float)

    array = np.asarray(values, dtype=float)

    if array.size != length:
        padded = np.zeros(length, dtype=float)
        padded[: min(array.size, length)] = array[:length]
        array = padded

    return np.nan_to_num(array)


def overlay_bigwigs(
    bigwig_file1: Path,
    bigwig_file2: Path,
    output_wig: Path,
    *,
    action: str,
    chunk_size: int,
) -> None:
    """Overlay two BigWigs using the original chunked WIG algorithm."""
    callback = resolve_action(action)
    bigwig1 = open_bigwig(bigwig_file1)

    try:
        bigwig2 = open_bigwig(bigwig_file2)

        try:
            print(
                "Get chromosome sizes from BigWig headers ...",
                file=sys.stderr,
            )

            chromosome_sizes = combined_chromosome_sizes(
                bigwig1.chroms(),
                bigwig2.chroms(),
            )

            with output_wig.open("w", encoding="utf-8") as output:
                for chromosome, chromosome_size in chromosome_sizes.items():
                    print(
                        f"Processing {chromosome} ...",
                        file=sys.stderr,
                    )
                    output.write(f"variableStep chrom={chromosome}\n")

                    for _, start, end in BED.tillingBed(
                        chrName=chromosome,
                        chrSize=chromosome_size,
                        stepSize=chunk_size,
                    ):
                        first_has_signal = interval_has_signal(
                            bigwig1,
                            chromosome,
                            start,
                            end,
                        )
                        second_has_signal = interval_has_signal(
                            bigwig2,
                            chromosome,
                            start,
                            end,
                        )

                        if not first_has_signal and not second_has_signal:
                            continue

                        signal1 = interval_values(
                            bigwig1,
                            chromosome,
                            start,
                            end,
                        )
                        signal2 = interval_values(
                            bigwig2,
                            chromosome,
                            start,
                            end,
                        )

                        result = callback(signal1, signal2)
                        coordinate = start

                        for value in result:
                            coordinate += 1
                            if value != 0:
                                print(
                                    f"{coordinate}\t{value:.2f}",
                                    file=output,
                                )

        finally:
            bigwig2.close()

    finally:
        bigwig1.close()


def main(argv: Sequence[str] | None = None) -> int:
    """Run BigWig overlay."""
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(parser, args)

    try:
        overlay_bigwigs(
            bigwig_file1=args.bigwig_file1,
            bigwig_file2=args.bigwig_file2,
            output_wig=args.output_wig,
            action=args.action,
            chunk_size=args.chunk_size,
        )
    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(1, f"{parser.prog}: error: {exc}\n")

    print(f"Created: {args.output_wig}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

