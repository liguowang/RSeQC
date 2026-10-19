#!/usr/bin/env python3
"""Collect valid BAM files from several forms of input.

Input may be:
1. A directory containing one or more BAM files.
2. A plain-text file containing BAM file paths, one per line.
3. A single BAM file.
4. A comma-separated list of BAM files.

A BAM file is considered valid here when it exists, is non-empty, and has
an accompanying BAI index. Both common BAI naming conventions are accepted:
``sample.bam.bai`` and ``sample.bai``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

__author__ = "Liguo Wang"

def _find_bai(bam_path):
    """Return the BAI path for *bam_path*, or ``None`` if no index exists."""
    bam = Path(bam_path)

    # samtools commonly uses sample.bam.bai, while sample.bai is also seen.
    candidates = (
        Path(str(bam) + ".bai"),
        bam.with_suffix(".bai"),
    )
    for bai in candidates:
        if bai.is_file():
            return bai
    return None


def isbamfile(infile):
    """Return True if *infile* is a non-empty BAM with an existing BAI index."""
    bam = Path(infile).expanduser()

    if not bam.is_file() or bam.suffix.lower() != ".bam":
        return False

    try:
        if bam.stat().st_size == 0:
            print(f"The size of {bam} is 0! Skip it.", file=sys.stderr)
            return False
    except OSError as exc:
        print(f"Cannot access {bam}: {exc}. Skip it.", file=sys.stderr)
        return False

    if _find_bai(bam) is None:
        print(
            f"Warning: BAM index does not exist for {bam}! Skip it.",
            file=sys.stderr,
        )
        return False

    return True


def _iter_bams_from_list_file(list_file):
    """Yield BAM paths listed in a plain-text file.

    Relative paths are first interpreted relative to the current working
    directory, preserving historical behavior. If that path is not a valid
    BAM, a path relative to the list file's directory is also tried.
    """
    list_path = Path(list_file).expanduser()

    with list_path.open("r") as handle:
        for line_number, line in enumerate(handle, start=1):
            entry = line.strip()
            if not entry or entry.startswith("#"):
                continue

            candidate = Path(entry).expanduser()
            if isbamfile(candidate):
                yield candidate
                continue

            if not candidate.is_absolute():
                relative_candidate = list_path.parent / candidate
                if relative_candidate != candidate and isbamfile(relative_candidate):
                    yield relative_candidate


def get_bam_files(input, printit=False):
    """Return valid BAM files specified by *input*.

    Parameters
    ----------
    input : str or os.PathLike
        Directory, BAM file, text list file, or comma-separated BAM paths.
    printit : bool, default=False
        If True, print each returned BAM path.

    Returns
    -------
    list of str
        Valid BAM file paths. Duplicate paths are removed while preserving
        discovery order.
    """
    input_text = str(input).strip()
    input_path = Path(input_text).expanduser()
    bam_files = []

    if input_path.is_dir():
        # Preserve the original recursive directory search and symlink
        # following behavior, but make output deterministic.
        for path in sorted(input_path.rglob("*")):
            if path.is_file() and isbamfile(path):
                bam_files.append(str(path.resolve()))

    elif isbamfile(input_path):
        bam_files.append(str(input_path))

    elif input_path.is_file():
        try:
            bam_files.extend(
                str(path) for path in _iter_bams_from_list_file(input_path)
            )
        except (OSError, UnicodeError) as exc:
            print(
                f"Cannot read BAM list file {input_path}: {exc}",
                file=sys.stderr,
            )

    else:
        # Treat a non-file input as a comma-separated list. Whitespace around
        # entries is ignored.
        for item in input_text.split(","):
            item = item.strip()
            if item and isbamfile(item):
                bam_files.append(item)

    # Remove duplicates without changing the first-seen order.
    bam_files = list(dict.fromkeys(bam_files))

    if printit:
        for bam_file in bam_files:
            print(bam_file)

    return bam_files


def main(argv=None):
    """Command-line entry point."""
    parser = argparse.ArgumentParser(
        description="Collect valid BAM files and verify their BAI indexes."
    )
    parser.add_argument(
        "input",
        help=(
            "BAM file, directory, text file containing BAM paths, "
            "or comma-separated BAM files"
        ),
    )
    args = parser.parse_args(argv)

    get_bam_files(args.input, printit=True)


if __name__ == "__main__":
    main()
