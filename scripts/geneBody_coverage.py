#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Calculate RNA-seq read coverage across the gene body.

Notes
-----
1. Input BAM files must be coordinate-sorted and indexed.
2. SAM input is not supported.
3. Transcripts shorter than ``--minimum-length`` are skipped.
4. The minimum allowed transcript length is 100 bp.

The coverage, percentile, normalization, and skewness algorithms are preserved
from the original implementation.
"""

from __future__ import annotations
import argparse
import collections
import operator
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Iterable, Sequence
import numpy as np
import pysam
from qcmodule import getBamFiles
from qcmodule import mystat

__author__ = "Liguo Wang"
__version__ = "5.05"

def valid_name(value: str) -> str:
    """Convert a string into a valid R variable name."""
    symbols = set(
        "abcdefghijklmnopqrstuvwxyz"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "0123456789_."
    )
    digits = set("0123456789")

    result = "_".join(value.split())

    if not result:
        return "sample"

    if result[0] in digits:
        result = "V" + result

    return "".join(character if character in symbols else "_" for character in result)


def printlog(message: str, log_file: Path = Path("log.txt")) -> None:
    """Write a timestamped message to stderr and the log file."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    formatted = f"@ {timestamp}: {message}"

    print(formatted, file=sys.stderr)

    with log_file.open("a", encoding="utf-8") as handle:
        print(formatted, file=handle)


def pearson_moment_coefficient(values: Sequence[float]) -> float:
    """Calculate skewness using the original Pearson-moment algorithm."""
    middle_value = values[int(len(values) / 2)]
    sigma = np.std(values, ddof=1)

    standardized_cubes = [
        ((value - middle_value) / sigma) ** 3
        for value in values
    ]
    return float(np.mean(standardized_cubes))


def genebody_percentile(
    refbed: Path,
    mrna_length_cutoff: int = 100,
) -> dict[str, tuple[str, str, list[int]]]:
    """Return 100 percentile positions for each eligible BED12 transcript."""
    gene_percentiles: dict[str, tuple[str, str, list[int]]] = {}
    transcript_count = 0

    with refbed.open("r", encoding="utf-8") as handle:
        for line in handle:
            try:
                if line.startswith(("#", "track", "browser")):
                    continue

                fields = line.split()
                chrom = fields[0]
                tx_start = int(fields[1])
                tx_end = int(fields[2])
                gene_name = fields[3]
                strand = fields[5]

                gene_id = "_".join(
                    str(value)
                    for value in (
                        chrom,
                        tx_start,
                        tx_end,
                        gene_name,
                        strand,
                    )
                )

                exon_starts = [
                    int(value)
                    for value in fields[11].rstrip(",\n").split(",")
                ]
                exon_starts = [
                    start + tx_start
                    for start in exon_starts
                ]

                exon_sizes = [
                    int(value)
                    for value in fields[10].rstrip(",\n").split(",")
                ]
                exon_ends = [
                    start + size
                    for start, size in zip(exon_starts, exon_sizes)
                ]

                transcript_count += 1

            except (IndexError, ValueError):
                print(
                    "[NOTE: input BED must be 12-column] skipped this line: "
                    + line,
                    end=" ",
                    file=sys.stderr,
                )
                continue

            gene_all_bases: list[int] = []

            for start, end in zip(exon_starts, exon_ends):
                # Preserve the original 1-based coordinate construction.
                gene_all_bases.extend(range(start + 1, end + 1))

            if len(gene_all_bases) < mrna_length_cutoff:
                continue

            gene_percentiles[gene_id] = (
                chrom,
                strand,
                mystat.percentile_list(gene_all_bases),
            )

    printlog(f"Total {transcript_count} transcripts loaded")
    return gene_percentiles


def genebody_coverage(
    bam_path: Path,
    position_list: dict[str, tuple[str, str, list[int]]],
) -> collections.defaultdict[int, int]:
    """Calculate aggregated coverage using the original pileup algorithm."""
    aggregated_coverage: collections.defaultdict[int, int] = (
        collections.defaultdict(int)
    )

    gene_finished = 0

    with pysam.AlignmentFile(str(bam_path), "rb") as samfile:
        for chrom, strand, positions in position_list.values():
            coverage = {position: 0.0 for position in positions}

            chrom_start = positions[0] - 1
            if chrom_start < 0:
                chrom_start = 0

            chrom_end = positions[-1]

            try:
                # Preserve the original reference-existence check.
                next(samfile.pileup(chrom, 1, 2), None)
            except (ValueError, OSError):
                continue

            for pileup_column in samfile.pileup(
                chrom,
                chrom_start,
                chrom_end,
                truncate=True,
            ):
                reference_position = pileup_column.pos + 1

                if reference_position not in positions:
                    continue

                if pileup_column.n == 0:
                    coverage[reference_position] = 0
                    continue

                covered_reads = 0

                for pileup_read in pileup_column.pileups:
                    if pileup_read.is_del:
                        continue

                    alignment = pileup_read.alignment

                    if alignment.is_qcfail:
                        continue
                    if alignment.is_secondary:
                        continue
                    if alignment.is_unmapped:
                        continue
                    if alignment.is_duplicate:
                        continue

                    covered_reads += 1

                coverage[reference_position] = covered_reads

            transcript_coverage = [
                coverage[position]
                for position in sorted(coverage)
            ]

            if strand == "-":
                transcript_coverage.reverse()

            for index, value in enumerate(transcript_coverage):
                aggregated_coverage[index] += value

            gene_finished += 1

            if gene_finished % 100 == 0:
                print(
                    f"\t{gene_finished} transcripts finished\r",
                    end=" ",
                    file=sys.stderr,
                )

    return aggregated_coverage


def write_r_code(
    dataset: Sequence[tuple[str, Sequence[float], float]],
    file_prefix: Path,
    *,
    output_format: str = "pdf",
    column_count: int = 100,
) -> Path:
    """Generate the original R plotting script."""
    r_script = Path(f"{file_prefix}.r")
    names: list[str] = []

    with r_script.open("w", encoding="utf-8") as output:
        for name, data, _skewness in dataset:
            names.append(name)
            print(
                f"{name} <- c({','.join(str(value) for value in data)})",
                file=output,
            )

        tick_positions = [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        tick_labels = [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

        if len(names) >= 3:
            print(
                "data_matrix <- matrix(c("
                + ",".join(names)
                + f"), byrow=T, ncol={column_count})",
                file=output,
            )
            print(
                "rowLabel <- c("
                + ",".join(f'"{name}"' for name in names)
                + ")",
                file=output,
            )
            print(file=output)
            print(
                f'{output_format.lower()}("{file_prefix}.heatMap.'
                f'{output_format.lower()}")',
                file=output,
            )
            print("rc <- cm.colors(ncol(data_matrix))", file=output)
            print(
                "heatmap(data_matrix, scale=c(\"none\"),keep.dendro=F, "
                "labRow=rowLabel,Colv=NA,Rowv=NA,labCol=NA,"
                "col=cm.colors(256),margins=c(6,8),ColSideColors=rc,"
                "cexRow=1,cexCol=1,"
                "xlab=\"Gene body percentile (5'->3')\","
                "add.expr=x_axis_expr <- axis("
                f"side=1,at=c({','.join(str(i) for i in tick_positions)}),"
                "labels=c("
                + ",".join(f'"{i}"' for i in tick_labels)
                + ")))",
                file=output,
            )
            print("dev.off()", file=output)

        print(file=output)
        print(
            f'{output_format.lower()}("{file_prefix}.curves.'
            f'{output_format.lower()}")',
            file=output,
        )
        print(f"x=1:{column_count}", file=output)
        print(
            'icolor = colorRampPalette(c('
            '"#7fc97f","#beaed4","#fdc086","#ffff99",'
            '"#386cb0","#f0027f"))'
            f"({len(names)})",
            file=output,
        )

        if len(names) == 1:
            print(
                f"plot(x,{names[0]},type='l',"
                "xlab=\"Gene body percentile (5'->3')\","
                "ylab=\"Coverage\",lwd=0.8,col=icolor[1])",
                file=output,
            )

        elif 2 <= len(names) <= 6:
            print(
                f"plot(x,{names[0]},type='l',"
                "xlab=\"Gene body percentile (5'->3')\","
                "ylab=\"Coverage\",lwd=0.8,col=icolor[1])",
                file=output,
            )

            for index, name in enumerate(names[1:], start=2):
                print(
                    f"lines(x,{name},type='l',col=icolor[{index}])",
                    file=output,
                )

            print(
                f"legend(0,1,fill=icolor[1:{len(names)}],legend=c("
                + ",".join(f"'{name}'" for name in names)
                + "))",
                file=output,
            )

        elif len(names) > 6:
            print(
                "layout(matrix(c(1,1,1,2,1,1,1,2,1,1,1,2),"
                "4,4,byrow=TRUE))",
                file=output,
            )
            print(
                f"plot(x,{names[0]},type='l',"
                "xlab=\"Gene body percentile (5'->3')\","
                "ylab=\"Coverage\",lwd=0.8,col=icolor[1])",
                file=output,
            )

            for index, name in enumerate(names[1:], start=2):
                print(
                    f"lines(x,{name},type='l',col=icolor[{index}])",
                    file=output,
                )

            print("par(mar=c(1,0,2,1))", file=output)
            print("plot.new()", file=output)
            print(
                f"legend(0,1,fill=icolor[1:{len(names)}],legend=c("
                + ",".join(f"'{name}'" for name in names)
                + "))",
                file=output,
            )

        print("dev.off()", file=output)

    return r_script


def build_parser() -> argparse.ArgumentParser:
    """Create and configure the command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "-i",
        "--input",
        dest="input_files",
        required=True,
        help=(
            "A BAM file, comma-separated BAM files, a directory containing "
            "BAM files, or a text file listing BAM paths."
        ),
    )
    parser.add_argument(
        "-r",
        "--refgene",
        dest="ref_gene_model",
        required=True,
        type=Path,
        help="Reference gene model in BED12 format.",
    )
    parser.add_argument(
        "-l",
        "--minimum-length",
        "--minimum_length",
        dest="min_mrna_length",
        type=int,
        default=100,
        metavar="INT",
        help=(
            "Minimum transcript length in bp. Values below 100 are not "
            "allowed. Default: %(default)s"
        ),
    )
    parser.add_argument(
        "-f",
        "--format",
        dest="output_format",
        choices=("pdf", "png", "jpeg"),
        default="pdf",
        help="Plot output format. Default: %(default)s",
    )
    parser.add_argument(
        "-o",
        "--out-prefix",
        dest="output_prefix",
        required=True,
        type=Path,
        help="Prefix for output files.",
    )
    parser.add_argument(
        "--skip-plot",
        action="store_true",
        help="Generate data and R code but do not execute the R script.",
    )
    parser.add_argument(
        "--rscript",
        default="Rscript",
        metavar="PATH",
        help="Rscript executable to use. Default: %(default)s",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Validate command-line arguments."""
    if not args.ref_gene_model.is_file():
        parser.error(
            f"reference BED12 file does not exist: {args.ref_gene_model}"
        )

    if args.min_mrna_length < 100:
        parser.error("--minimum-length cannot be smaller than 100")

    output_parent = args.output_prefix.parent
    if not output_parent.exists():
        parser.error(f"output directory does not exist: {output_parent}")

    if not output_parent.is_dir():
        parser.error(f"output parent is not a directory: {output_parent}")


def make_unique_sample_name(name: str, seen: list[str]) -> str:
    """Return a unique R-safe sample name using the original naming scheme."""
    count = seen.count(name)
    return name if count == 0 else f"{name}.{count}"


def load_dataset(coverage_file: Path) -> list[tuple[str, list[float], float]]:
    """Load, normalize, and sort coverage using the original algorithm."""
    dataset: list[tuple[str, list[float], float]] = []

    with coverage_file.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()

            if not stripped or stripped.startswith("Percentile"):
                continue

            fields = stripped.split()
            name = fields[0]
            values = [float(value) for value in fields[1:]]

            skewness = pearson_moment_coefficient(values)
            minimum = min(values)
            maximum = max(values)

            # Preserve the original min-max normalization formula.
            normalized = [
                (value - minimum) / (maximum - minimum)
                for value in values
            ]

            dataset.append((name, normalized, skewness))

    dataset.sort(key=operator.itemgetter(2), reverse=True)
    return dataset


def run_r_script(
    parser: argparse.ArgumentParser,
    r_script: Path,
    executable: str,
) -> None:
    """Execute the generated R script safely."""
    resolved = shutil.which(executable)

    if resolved is None:
        parser.exit(
            1,
            f'{parser.prog}: error: Rscript executable not found: '
            f'"{executable}"\n',
        )

    try:
        subprocess.run(
            [resolved, str(r_script)],
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        parser.exit(
            exc.returncode or 1,
            f"{parser.prog}: error: R plotting failed for {r_script}\n",
        )


def main(argv: Sequence[str] | None = None) -> int:
    """Run gene-body coverage analysis."""
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(parser, args)

    coverage_file = Path(f"{args.output_prefix}.geneBodyCoverage.txt")
    plot_prefix = Path(f"{args.output_prefix}.geneBodyCoverage")

    try:
        printlog("Read BED file (reference gene model) ...")
        gene_percentiles = genebody_percentile(
            refbed=args.ref_gene_model,
            mrna_length_cutoff=args.min_mrna_length,
        )

        printlog("Get BAM file(s) ...")
        bam_files = [
            Path(path)
            for path in getBamFiles.get_bam_files(args.input_files)
        ]

        if not bam_files:
            raise ValueError("no BAM files were found")

        for bam_file in bam_files:
            print(f"\t{bam_file}", file=sys.stderr)

        sample_names: list[str] = []

        with coverage_file.open("w", encoding="utf-8") as output:
            print(
                "Percentile\t"
                + "\t".join(str(value) for value in range(1, 101)),
                file=output,
            )

            for bam_file in bam_files:
                if not bam_file.is_file():
                    print(
                        f"Warning: BAM file does not exist; skipped: {bam_file}",
                        file=sys.stderr,
                    )
                    continue

                printlog(f"Processing {bam_file.name} ...")
                coverage = genebody_coverage(bam_file, gene_percentiles)

                if not coverage:
                    print(
                        f"\nCannot get coverage signal from "
                        f"{bam_file.name}! Skip",
                        file=sys.stderr,
                    )
                    continue

                base_name = valid_name(bam_file.stem)
                sample_name = make_unique_sample_name(
                    base_name,
                    sample_names,
                )
                sample_names.append(base_name)

                print(
                    sample_name
                    + "\t"
                    + "\t".join(
                        str(coverage[index])
                        for index in sorted(coverage)
                    ),
                    file=output,
                )

        dataset = load_dataset(coverage_file)

        if not dataset:
            raise ValueError("no valid coverage profiles were generated")

        print("\n\n", file=sys.stderr)
        print("\tSample\tSkewness", file=sys.stderr)

        for name, _data, skewness in dataset:
            print(f"\t{name}\t{skewness}", file=sys.stderr)

        r_script = write_r_code(
            dataset,
            plot_prefix,
            output_format=args.output_format,
        )

        if not args.skip_plot:
            printlog("Running R script ...")
            run_r_script(parser, r_script, args.rscript)

        printlog(f"Created: {coverage_file}")
        printlog(f"Created: {r_script}")

    except (OSError, ValueError, ZeroDivisionError) as exc:
        parser.exit(1, f"{parser.prog}: error: {exc}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
