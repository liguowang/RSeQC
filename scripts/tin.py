#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Calculate transcript integrity number (TIN) for each transcript or gene.

TIN is conceptually similar to RIN but is calculated per transcript. Scores
range from 0 to 100, where higher values indicate more uniform transcript
coverage.

The historical TIN algorithm, read filters, sampling logic, and optional
intronic-background subtraction are preserved.
"""

from __future__ import annotations

import argparse
import logging
import math
import sys
from pathlib import Path
from typing import Iterator, Sequence
import numpy as np
import pysam
from bx.intervals.intersection import Intersecter, Interval
from qcmodule import BED
from qcmodule import getBamFiles
from importlib.metadata import version

__version__ = version("RSeQC")
__author__ = "Liguo Wang"


DESCRIPTION = (
    "Calculate transcript-level TIN scores and sample-level mean, median, "
    "and standard deviation from one or more sorted, indexed BAM files."
)

EPILOG = """
Input forms accepted by --input
-------------------------------
* One BAM file
* A comma-separated BAM list
* A directory containing BAM files
* A text file listing one BAM path per line

Example
-------
tin.py \
    -i sample.bam \
    -r genes.bed12 \
    -c 10 \
    -n 100

With intronic-background subtraction:
    tin.py \
        -i sample.bam \
        -r genes.bed12 \
        -s

Notes
-----
* BAM files must be coordinate-sorted and indexed.
* The BED file must contain 12 columns.
* TIN is reported as 0 when coverage is insufficient.
"""


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=DESCRIPTION,
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input_files",
        required=True,
        help="BAM input specification; see the accepted forms below.",
    )
    parser.add_argument(
        "-r",
        "--refgene",
        dest="ref_gene_model",
        required=True,
        type=Path,
        help="Reference gene model in standard BED12 format.",
    )
    parser.add_argument(
        "-c",
        "--minCov",
        dest="minimum_coverage",
        type=int,
        default=10,
        metavar="INT",
        help="Minimum number of distinct read starts required. Default: %(default)s",
    )
    parser.add_argument(
        "-n",
        "--sample-size",
        dest="sample_size",
        type=int,
        default=100,
        metavar="INT",
        help="Number of approximately equally spaced transcript positions. Default: %(default)s",
    )
    parser.add_argument(
        "-s",
        "--subtract-background",
        dest="subtract_bg",
        action="store_true",
        help="Subtract background estimated from intronic signal.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Directory for TIN output files. Default: current directory",
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
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.DEBUG if verbose else logging.INFO,
    )


def find_bam_index(bam_path: Path) -> Path | None:
    candidates = (
        Path(f"{bam_path}.bai"),
        bam_path.with_suffix(".bai"),
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


def validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    if not args.ref_gene_model.is_file():
        parser.error(f"reference BED file does not exist: {args.ref_gene_model}")

    if args.minimum_coverage < 0:
        parser.error("--minCov must be zero or greater")

    if args.sample_size <= 0:
        parser.error("--sample-size must be greater than zero")

    if args.sample_size > 1000:
        logging.warning(
            "--sample-size is greater than 1000; reduce it if performance is poor"
        )

    if not args.output_dir.exists():
        parser.error(f"output directory does not exist: {args.output_dir}")

    if not args.output_dir.is_dir():
        parser.error(f"output path is not a directory: {args.output_dir}")


def uniqify(values: Sequence[int]) -> list[int]:
    """Remove duplicates while preserving order."""
    return list(dict.fromkeys(values))


def shannon_entropy(values: Sequence[float]) -> float:
    """Calculate Shannon entropy using natural logarithms."""
    total = sum(values)
    if total <= 0:
        return 0.0

    entropy = 0.0
    for value in values:
        probability = value / total
        entropy += probability * math.log(probability)

    return 0.0 if entropy == 0 else -entropy


def build_interval_trees(intervals) -> dict[str, Intersecter]:
    ranges: dict[str, Intersecter] = {}

    for chromosome, start, end in intervals:
        tree = ranges.setdefault(str(chromosome), Intersecter())
        tree.add_interval(Interval(int(start), int(end)))

    return ranges


def union_exons(refbed: Path) -> dict[str, Intersecter]:
    bed = BED.ParseBED(str(refbed))
    all_exons = bed.getExon()
    unioned_exons = BED.unionBed3(all_exons)
    return build_interval_trees(unioned_exons)


def estimate_bg_noise(
    chromosome: str,
    tx_start: int,
    tx_end: int,
    samfile: pysam.AlignmentFile,
    exon_ranges: dict[str, Intersecter],
) -> float:
    """Estimate transcript background from non-exonic aligned read bases."""
    intron_signal = 0.0
    tree = exon_ranges.get(chromosome)

    for aligned_read in samfile.fetch(chromosome, tx_start, tx_end):
        if aligned_read.is_qcfail:
            continue
        if aligned_read.is_unmapped:
            continue
        if aligned_read.is_secondary:
            continue

        read_start = aligned_read.reference_start
        if read_start < tx_start or read_start >= tx_end:
            continue

        read_length = aligned_read.query_length
        if read_length is None:
            continue

        if tree is not None and tree.find(read_start, read_start + read_length):
            continue

        intron_signal += read_length

    return intron_signal


def genomic_positions(
    refbed: Path,
    sample_size: int,
) -> Iterator[tuple[str, str, int, int, int, list[int]]]:
    """Yield sampled genomic transcript positions using the original logic."""
    with refbed.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            if line.startswith(("#", "track", "browser")):
                continue

            fields = line.split()

            try:
                if len(fields) < 12:
                    raise ValueError("fewer than 12 columns")

                chromosome = fields[0]
                tx_start = int(fields[1])
                tx_end = int(fields[2])
                gene_name = fields[3]

                block_sizes = [
                    int(value)
                    for value in fields[10].rstrip(",").split(",")
                    if value
                ]
                block_starts = [
                    int(value)
                    for value in fields[11].rstrip(",").split(",")
                    if value
                ]

                exon_starts = [
                    tx_start + start
                    for start in block_starts
                ]
                exon_ends = [
                    start + size
                    for start, size in zip(exon_starts, block_sizes)
                ]

                mrna_size = sum(block_sizes)
                intron_size = max(0, tx_end - tx_start - mrna_size)

            except (ValueError, IndexError) as exc:
                logging.warning(
                    "Skipping malformed BED12 line %d: %s",
                    line_number,
                    exc,
                )
                continue

            chosen_bases = [tx_start + 1, tx_end]

            if mrna_size <= sample_size:
                for start, end in zip(exon_starts, exon_ends):
                    chosen_bases.extend(range(start + 1, end + 1))
            else:
                step_size = max(1, int(mrna_size / sample_size))
                all_bases: list[int] = []
                exon_bounds: list[int] = []

                for start, end in zip(exon_starts, exon_ends):
                    all_bases.extend(range(start + 1, end + 1))
                    exon_bounds.extend((start + 1, end))

                sampled_indexes = range(0, len(all_bases), step_size)
                sampled_bases = [all_bases[index] for index in sampled_indexes]
                chosen_bases = uniqify(exon_bounds + sampled_bases)

            yield (
                gene_name,
                chromosome,
                tx_start,
                tx_end,
                intron_size,
                chosen_bases,
            )


def check_min_reads(
    samfile: pysam.AlignmentFile,
    chromosome: str,
    tx_start: int,
    tx_end: int,
    cutoff: int,
) -> bool:
    """Check the original distinct-read-start coverage criterion."""
    read_starts: set[int] = set()

    try:
        for aligned_read in samfile.fetch(chromosome, tx_start, tx_end):
            if aligned_read.is_qcfail:
                continue
            if aligned_read.is_unmapped:
                continue
            if aligned_read.is_secondary:
                continue

            read_start = aligned_read.reference_start
            if read_start < tx_start or read_start >= tx_end:
                continue

            read_starts.add(read_start)

            # Preserve historical strict-greater-than behavior.
            if len(read_starts) > cutoff:
                return True

    except (ValueError, KeyError):
        return False

    return False


def genebody_coverage(
    samfile: pysam.AlignmentFile,
    chromosome: str,
    positions: Sequence[int],
    bg_level: float = 0,
) -> list[float]:
    """Calculate coverage at sampled positions using the original filters."""
    coverage: list[float] = []
    position_set = set(positions)
    start = positions[0] - 1
    end = positions[-1]

    try:
        for pileup_column in samfile.pileup(
            chromosome,
            start,
            end,
            truncate=True,
        ):
            reference_position = pileup_column.reference_pos + 1

            if reference_position not in position_set:
                continue

            if pileup_column.n == 0:
                coverage.append(0.0)
                continue

            covered_reads = 0.0

            for pileup_read in pileup_column.pileups:
                if pileup_read.is_del:
                    continue
                if pileup_read.alignment.is_qcfail:
                    continue
                if pileup_read.alignment.is_secondary:
                    continue
                if pileup_read.alignment.is_unmapped:
                    continue

                covered_reads += 1.0

            coverage.append(covered_reads)

    except (ValueError, KeyError):
        return []

    if bg_level <= 0:
        return coverage

    corrected: list[float] = []
    for value in coverage:
        subtracted = int(value - bg_level)
        corrected.append(subtracted if subtracted > 0 else 0)

    return corrected


def tin_score(coverage: Sequence[float], sampled_length: int) -> float:
    """Calculate TIN using the historical entropy formula."""
    if not coverage or sampled_length <= 0:
        return 0.0

    effective = [float(value) for value in coverage if float(value) > 0]
    if not effective:
        return 0.0

    entropy = shannon_entropy(effective)
    return 100.0 * math.exp(entropy) / sampled_length


def output_paths(bam_file: Path, output_dir: Path) -> tuple[Path, Path]:
    """Return transcript-level and summary output paths."""
    stem = bam_file.name
    if stem.lower().endswith(".bam"):
        stem = stem[:-4]

    return (
        output_dir / f"{stem}.tin.xls",
        output_dir / f"{stem}.summary.txt",
    )


def process_bam(
    bam_file: Path,
    ref_gene_model: Path,
    minimum_coverage: int,
    sample_size: int,
    subtract_background: bool,
    exon_ranges: dict[str, Intersecter] | None,
    output_dir: Path,
) -> None:
    """Calculate TIN for one BAM file."""
    tin_path, summary_path = output_paths(bam_file, output_dir)
    sample_tins: list[float] = []

    with pysam.AlignmentFile(str(bam_file), "rb") as samfile:
        with (
            tin_path.open("w", encoding="utf-8") as output,
            summary_path.open("w", encoding="utf-8") as summary,
        ):
            print(
                "\t".join(("geneID", "chrom", "tx_start", "tx_end", "TIN")),
                file=output,
            )
            print(
                "\t".join(("Bam_file", "TIN(mean)", "TIN(median)", "TIN(stdev)")),
                file=summary,
            )

            finished = 0

            for (
                gene_name,
                chromosome,
                tx_start,
                tx_end,
                intron_size,
                picked_positions,
            ) in genomic_positions(ref_gene_model, sample_size):
                finished += 1

                if not check_min_reads(
                    samfile,
                    chromosome,
                    tx_start,
                    tx_end,
                    minimum_coverage,
                ):
                    score = 0.0
                else:
                    noise_level = 0.0

                    if subtract_background and exon_ranges is not None:
                        intron_signal = estimate_bg_noise(
                            chromosome,
                            tx_start,
                            tx_end,
                            samfile,
                            exon_ranges,
                        )
                        if intron_size > 0:
                            noise_level = intron_signal / intron_size

                    positions = sorted(picked_positions)
                    coverage = genebody_coverage(
                        samfile,
                        chromosome,
                        positions,
                        noise_level,
                    )
                    score = tin_score(coverage, len(picked_positions))
                    sample_tins.append(score)

                print(
                    "\t".join(
                        map(
                            str,
                            (
                                gene_name,
                                chromosome,
                                tx_start,
                                tx_end,
                                score,
                            ),
                        )
                    ),
                    file=output,
                )

                if finished % 100 == 0:
                    logging.info("%d transcripts finished", finished)

            if sample_tins:
                mean_tin = float(np.mean(sample_tins))
                median_tin = float(np.median(sample_tins))
                std_tin = float(np.std(sample_tins))
            else:
                mean_tin = median_tin = std_tin = 0.0

            print(
                "\t".join(
                    map(
                        str,
                        (
                            bam_file.name,
                            mean_tin,
                            median_tin,
                            std_tin,
                        ),
                    )
                ),
                file=summary,
            )

    logging.info("Created %s", tin_path)
    logging.info("Created %s", summary_path)


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    configure_logging(args.verbose)
    validate_args(parser, args)

    logging.info("Get BAM file(s) ...")
    bam_files = sorted(
        Path(path)
        for path in getBamFiles.get_bam_files(args.input_files)
    )

    if not bam_files:
        logging.error("No BAM files found")
        return 1

    for bam_file in bam_files:
        if not bam_file.is_file():
            logging.error("BAM file does not exist: %s", bam_file)
            return 1

        if find_bam_index(bam_file) is None:
            logging.error(
                "BAM index not found for %s; expected %s or %s",
                bam_file,
                f"{bam_file}.bai",
                bam_file.with_suffix(".bai"),
            )
            return 1

    logging.info("Total %d BAM file(s)", len(bam_files))
    for bam_file in bam_files:
        logging.info("  %s", bam_file)

    try:
        exon_ranges = (
            union_exons(args.ref_gene_model)
            if args.subtract_bg
            else None
        )

        for bam_file in bam_files:
            logging.info("Processing %s", bam_file)
            process_bam(
                bam_file=bam_file,
                ref_gene_model=args.ref_gene_model,
                minimum_coverage=args.minimum_coverage,
                sample_size=args.sample_size,
                subtract_background=args.subtract_bg,
                exon_ranges=exon_ranges,
                output_dir=args.output_dir,
            )

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

    logging.info("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
