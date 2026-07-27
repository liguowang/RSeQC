#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Summarize read distribution across genomic annotation categories.

The following alignments are skipped:
- QC-failed
- PCR/optical duplicate
- Unmapped
- Secondary/non-primary

The original category precedence and midpoint-based assignment algorithm are
preserved.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from bx.intervals.intersection import Intersecter, Interval

from qcmodule import BED
from qcmodule import SAM
from qcmodule import bam_cigar


__author__ = "Liguo Wang"
__version__ = "5.05"


DESCRIPTION = (
    "Report how aligned read blocks are distributed across CDS exons, UTRs, "
    "introns, and flanking intergenic regions."
)

EPILOG = """
Assignment priority
-------------------
Each aligned block is assigned by its midpoint using this order:

1. CDS exon
2. 5' UTR
3. 3' UTR
4. Intron
5. Upstream/downstream intergenic bins
6. Unassigned

Example
-------
read_distribution.py -i sample.bam -r genes.bed12
"""


@dataclass(frozen=True)
class RegionModel:
    """Processed genomic regions and their total sizes."""

    cds_exon: dict[str, Intersecter]
    intron: dict[str, Intersecter]
    utr_5: dict[str, Intersecter]
    utr_3: dict[str, Intersecter]
    upstream_1kb: dict[str, Intersecter]
    upstream_5kb: dict[str, Intersecter]
    upstream_10kb: dict[str, Intersecter]
    downstream_1kb: dict[str, Intersecter]
    downstream_5kb: dict[str, Intersecter]
    downstream_10kb: dict[str, Intersecter]

    cds_exon_bases: int
    intron_bases: int
    utr_5_bases: int
    utr_3_bases: int
    upstream_1kb_bases: int
    upstream_5kb_bases: int
    upstream_10kb_bases: int
    downstream_1kb_bases: int
    downstream_5kb_bases: int
    downstream_10kb_bases: int


@dataclass
class DistributionCounts:
    """Read-block counts for each annotation category."""

    total_reads: int = 0
    total_tags: int = 0
    unassigned_tags: int = 0

    cds_exon: int = 0
    intron: int = 0
    utr_5: int = 0
    utr_3: int = 0

    upstream_1kb: int = 0
    upstream_5kb: int = 0
    upstream_10kb: int = 0
    downstream_1kb: int = 0
    downstream_5kb: int = 0
    downstream_10kb: int = 0

    qc_failed_reads: int = 0
    duplicate_reads: int = 0
    secondary_reads: int = 0
    unmapped_reads: int = 0


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
        help="Input alignment file in BAM or SAM format.",
    )
    parser.add_argument(
        "-r",
        "--refgene",
        dest="ref_gene_model",
        required=True,
        type=Path,
        help="Reference gene model in BED format.",
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
    ref_gene_model: Path,
) -> None:
    """Validate command-line arguments."""
    if not input_file.is_file():
        parser.error(
            f"input alignment file does not exist or is not a file: {input_file}"
        )

    if input_file.suffix.lower() not in {".bam", ".sam", ".cram"}:
        parser.error(
            "input alignment file must have a .bam, .sam, or .cram "
            f"extension: {input_file}"
        )

    if not ref_gene_model.is_file():
        parser.error(
            f"reference gene model does not exist or is not a file: "
            f"{ref_gene_model}"
        )


def calculate_size(intervals: Sequence[Sequence[object]]) -> int:
    """Calculate total length of BED-like intervals."""
    return sum(int(interval[2]) - int(interval[1]) for interval in intervals)


def build_interval_trees(
    intervals: Sequence[Sequence[object]],
) -> dict[str, Intersecter]:
    """Build chromosome-specific interval trees."""
    ranges: dict[str, Intersecter] = {}

    for interval in intervals:
        chromosome = str(interval[0]).upper()
        start = int(interval[1])
        end = int(interval[2])

        tree = ranges.setdefault(chromosome, Intersecter())
        tree.add_interval(Interval(start, end))

    return ranges


def overlaps_point(
    chromosome: str,
    ranges: dict[str, Intersecter],
    position: int,
) -> bool:
    """Return whether a 1-bp interval at position overlaps a region."""
    tree = ranges.get(chromosome)
    if tree is None:
        return False

    return bool(tree.find(position, position))


def purify_regions(
    regions,
    cds_exon,
    utr_5,
    utr_3,
    intron,
):
    """Remove annotated genic regions from a flanking-region set."""
    result = BED.subtractBed3(regions, cds_exon)
    result = BED.subtractBed3(result, utr_5)
    result = BED.subtractBed3(result, utr_3)
    result = BED.subtractBed3(result, intron)
    return result


def process_gene_model(gene_model: Path) -> RegionModel:
    """Build processed annotation regions using the original algorithm."""
    print(f"Processing {gene_model} ...", end=" ", file=sys.stderr)

    bed = BED.ParseBED(str(gene_model))

    utr_3 = BED.unionBed3(bed.getUTR(utr=3))
    utr_5 = BED.unionBed3(bed.getUTR(utr=5))
    cds_exon = BED.unionBed3(bed.getCDSExon())
    intron = BED.unionBed3(bed.getIntron())

    utr_5 = BED.subtractBed3(utr_5, cds_exon)
    utr_3 = BED.subtractBed3(utr_3, cds_exon)

    intron = BED.subtractBed3(intron, cds_exon)
    intron = BED.subtractBed3(intron, utr_5)
    intron = BED.subtractBed3(intron, utr_3)

    upstream_1kb = BED.unionBed3(
        bed.getIntergenic(direction="up", size=1_000)
    )
    upstream_5kb = BED.unionBed3(
        bed.getIntergenic(direction="up", size=5_000)
    )
    upstream_10kb = BED.unionBed3(
        bed.getIntergenic(direction="up", size=10_000)
    )
    downstream_1kb = BED.unionBed3(
        bed.getIntergenic(direction="down", size=1_000)
    )
    downstream_5kb = BED.unionBed3(
        bed.getIntergenic(direction="down", size=5_000)
    )
    downstream_10kb = BED.unionBed3(
        bed.getIntergenic(direction="down", size=10_000)
    )

    upstream_1kb = purify_regions(
        upstream_1kb, cds_exon, utr_5, utr_3, intron
    )
    upstream_5kb = purify_regions(
        upstream_5kb, cds_exon, utr_5, utr_3, intron
    )
    upstream_10kb = purify_regions(
        upstream_10kb, cds_exon, utr_5, utr_3, intron
    )
    downstream_1kb = purify_regions(
        downstream_1kb, cds_exon, utr_5, utr_3, intron
    )
    downstream_5kb = purify_regions(
        downstream_5kb, cds_exon, utr_5, utr_3, intron
    )
    downstream_10kb = purify_regions(
        downstream_10kb, cds_exon, utr_5, utr_3, intron
    )

    model = RegionModel(
        cds_exon=build_interval_trees(cds_exon),
        intron=build_interval_trees(intron),
        utr_5=build_interval_trees(utr_5),
        utr_3=build_interval_trees(utr_3),
        upstream_1kb=build_interval_trees(upstream_1kb),
        upstream_5kb=build_interval_trees(upstream_5kb),
        upstream_10kb=build_interval_trees(upstream_10kb),
        downstream_1kb=build_interval_trees(downstream_1kb),
        downstream_5kb=build_interval_trees(downstream_5kb),
        downstream_10kb=build_interval_trees(downstream_10kb),
        cds_exon_bases=calculate_size(cds_exon),
        intron_bases=calculate_size(intron),
        utr_5_bases=calculate_size(utr_5),
        utr_3_bases=calculate_size(utr_3),
        upstream_1kb_bases=calculate_size(upstream_1kb),
        upstream_5kb_bases=calculate_size(upstream_5kb),
        upstream_10kb_bases=calculate_size(upstream_10kb),
        downstream_1kb_bases=calculate_size(downstream_1kb),
        downstream_5kb_bases=calculate_size(downstream_5kb),
        downstream_10kb_bases=calculate_size(downstream_10kb),
    )

    print("Done", file=sys.stderr)
    return model


def assign_tag(
    chromosome: str,
    midpoint: int,
    model: RegionModel,
    counts: DistributionCounts,
) -> None:
    """Assign one aligned block using the original category precedence."""
    if overlaps_point(chromosome, model.cds_exon, midpoint):
        counts.cds_exon += 1
        return

    overlaps_5utr = overlaps_point(chromosome, model.utr_5, midpoint)
    overlaps_3utr = overlaps_point(chromosome, model.utr_3, midpoint)

    if overlaps_5utr and not overlaps_3utr:
        counts.utr_5 += 1
        return

    if overlaps_3utr and not overlaps_5utr:
        counts.utr_3 += 1
        return

    if overlaps_5utr and overlaps_3utr:
        counts.unassigned_tags += 1
        return

    if overlaps_point(chromosome, model.intron, midpoint):
        counts.intron += 1
        return

    overlaps_up_10kb = overlaps_point(
        chromosome, model.upstream_10kb, midpoint
    )
    overlaps_down_10kb = overlaps_point(
        chromosome, model.downstream_10kb, midpoint
    )

    if overlaps_up_10kb and overlaps_down_10kb:
        counts.unassigned_tags += 1
        return

    if overlaps_point(chromosome, model.upstream_1kb, midpoint):
        counts.upstream_1kb += 1
        counts.upstream_5kb += 1
        counts.upstream_10kb += 1
        return

    if overlaps_point(chromosome, model.upstream_5kb, midpoint):
        counts.upstream_5kb += 1
        counts.upstream_10kb += 1
        return

    if overlaps_up_10kb:
        counts.upstream_10kb += 1
        return

    if overlaps_point(chromosome, model.downstream_1kb, midpoint):
        counts.downstream_1kb += 1
        counts.downstream_5kb += 1
        counts.downstream_10kb += 1
        return

    if overlaps_point(chromosome, model.downstream_5kb, midpoint):
        counts.downstream_5kb += 1
        counts.downstream_10kb += 1
        return

    if overlaps_down_10kb:
        counts.downstream_10kb += 1
        return

    counts.unassigned_tags += 1


def count_read_distribution(
    input_file: Path,
    model: RegionModel,
) -> DistributionCounts:
    """Count read-block distribution using the original algorithm."""
    counts = DistributionCounts()
    alignment_file = SAM.ParseBAM(str(input_file)).samfile

    print(f"Processing {input_file} ...", end=" ", file=sys.stderr)

    for aligned_read in alignment_file:
        if aligned_read.is_qcfail:
            counts.qc_failed_reads += 1
            continue

        if aligned_read.is_duplicate:
            counts.duplicate_reads += 1
            continue

        if aligned_read.is_secondary:
            counts.secondary_reads += 1
            continue

        if aligned_read.is_unmapped:
            counts.unmapped_reads += 1
            continue

        counts.total_reads += 1

        chromosome = alignment_file.getrname(aligned_read.tid).upper()
        exons = bam_cigar.fetch_exon(
            chromosome,
            aligned_read.pos,
            aligned_read.cigar,
        )
        counts.total_tags += len(exons)

        for exon in exons:
            start = int(exon[1])
            end = int(exon[2])
            midpoint = start + int((end - start) / 2)

            assign_tag(
                chromosome=chromosome,
                midpoint=midpoint,
                model=model,
                counts=counts,
            )

    print("Finished\n", file=sys.stderr)
    return counts


def tags_per_kb(tag_count: int, base_count: int) -> float:
    """Calculate tags per kb using the original +1 denominator."""
    return tag_count * 1_000.0 / (base_count + 1)


def print_row(
    label: str,
    bases: int,
    tags: int,
) -> None:
    """Print one formatted distribution row."""
    print(
        f"{label:<20}{bases:<20d}{tags:<20d}"
        f"{tags_per_kb(tags, bases):<18.2f}"
    )


def print_report(
    model: RegionModel,
    counts: DistributionCounts,
) -> None:
    """Print the original read-distribution report."""
    print(f"{'Total Reads':<30}{counts.total_reads}")
    print(f"{'Total Tags':<30}{counts.total_tags}")
    print(
        f"{'Total Assigned Tags':<30}"
        f"{counts.total_tags - counts.unassigned_tags}"
    )

    print("=" * 69)
    print(
        f"{'Group':<20}{'Total_bases':<20}"
        f"{'Tag_count':<20}{'Tags/Kb':<20}"
    )

    print_row("CDS_Exons", model.cds_exon_bases, counts.cds_exon)
    print_row("5'UTR_Exons", model.utr_5_bases, counts.utr_5)
    print_row("3'UTR_Exons", model.utr_3_bases, counts.utr_3)
    print_row("Introns", model.intron_bases, counts.intron)

    print_row(
        "TSS_up_1kb",
        model.upstream_1kb_bases,
        counts.upstream_1kb,
    )
    print_row(
        "TSS_up_5kb",
        model.upstream_5kb_bases,
        counts.upstream_5kb,
    )
    print_row(
        "TSS_up_10kb",
        model.upstream_10kb_bases,
        counts.upstream_10kb,
    )
    print_row(
        "TES_down_1kb",
        model.downstream_1kb_bases,
        counts.downstream_1kb,
    )
    print_row(
        "TES_down_5kb",
        model.downstream_5kb_bases,
        counts.downstream_5kb,
    )
    print_row(
        "TES_down_10kb",
        model.downstream_10kb_bases,
        counts.downstream_10kb,
    )

    print("=" * 69)


def main(argv: Sequence[str] | None = None) -> int:
    """Run read-distribution analysis."""
    parser = build_parser()
    args = parser.parse_args(argv)

    validate_args(
        parser=parser,
        input_file=args.input_file,
        ref_gene_model=args.ref_gene_model,
    )

    try:
        model = process_gene_model(args.ref_gene_model)
        counts = count_read_distribution(args.input_file, model)
        print_report(model, counts)

    except (OSError, ValueError, RuntimeError, IndexError) as exc:
        parser.exit(1, f"{parser.prog}: error: {exc}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
