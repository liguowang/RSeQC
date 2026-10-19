#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Calculate RNA-seq coverage across the gene body from a BigWig file.

The transcript-coordinate construction, 100-point percentile sampling, and
coverage aggregation algorithm are preserved from the original implementation.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Sequence
import numpy as np
import pyBigWig
from qcmodule import mystat
from importlib.metadata import version

__version__ = version("RSeQC")
__author__ = "Liguo Wang, Santiago Revale"


def build_parser() -> argparse.ArgumentParser:
    """Create and configure the command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "-i",
        "--input-file",
        required=True,
        type=Path,
        help="Coverage signal file in BigWig format.",
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
        "-o",
        "--out-prefix",
        dest="output_prefix",
        required=True,
        type=Path,
        help="Prefix for generated output files.",
    )
    parser.add_argument(
        "-t",
        "--graph-type",
        dest="graph_type",
        choices=("pdf", "png", "bmp", "jpeg", "tiff"),
        default="pdf",
        help="Plot file type. Default: %(default)s",
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
    """Validate input files and output location."""
    if not args.input_file.is_file():
        parser.error(f"BigWig file does not exist: {args.input_file}")

    if not args.ref_gene_model.is_file():
        parser.error(
            f"reference BED12 file does not exist: {args.ref_gene_model}"
        )

    output_parent = args.output_prefix.parent
    if not output_parent.exists():
        parser.error(f"output directory does not exist: {output_parent}")

    if not output_parent.is_dir():
        parser.error(f"output parent is not a directory: {output_parent}")


def parse_bed12_line(
    line: str,
    *,
    path: Path,
    line_number: int,
) -> tuple[str, int, int, str, str, list[int], list[int]]:
    """Parse one BED12 record using the original coordinate calculations."""
    fields = line.split()

    if len(fields) < 12:
        raise ValueError(
            f"{path}:{line_number}: expected 12 BED columns, "
            f"found {len(fields)}"
        )

    chrom = fields[0]
    tx_start = int(fields[1])
    tx_end = int(fields[2])
    gene_name = fields[3]
    strand = fields[5]

    exon_starts = [
        int(value)
        for value in fields[11].rstrip(",\n").split(",")
        if value
    ]
    exon_starts = [
        exon_start + tx_start
        for exon_start in exon_starts
    ]

    exon_sizes = [
        int(value)
        for value in fields[10].rstrip(",\n").split(",")
        if value
    ]
    exon_ends = [
        exon_start + exon_size
        for exon_start, exon_size in zip(exon_starts, exon_sizes)
    ]

    return (
        chrom,
        tx_start,
        tx_end,
        gene_name,
        strand,
        exon_starts,
        exon_ends,
    )


def coverage_gene_body_bigwig(
    bigwig_file: Path,
    refbed: Path,
    output_prefix: Path,
    *,
    graph_type: str = "png",
) -> tuple[Path, Path, int]:
    """Calculate gene-body coverage with the original sampling algorithm."""
    coverage: defaultdict[int, float] = defaultdict(float)
    gene_count = 0

    bigwig = pyBigWig.open(str(bigwig_file), "r")
    if bigwig is None:
        raise OSError(f"could not open BigWig file: {bigwig_file}")

    try:
        chromosomes = set(bigwig.chroms())

        print("Calculating coverage over gene body ...", file=sys.stderr)

        with refbed.open("r", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                if line.startswith(("#", "track", "browser")):
                    continue

                try:
                    (
                        chrom,
                        _tx_start,
                        _tx_end,
                        _gene_name,
                        strand,
                        exon_starts,
                        exon_ends,
                    ) = parse_bed12_line(
                        line,
                        path=refbed,
                        line_number=line_number,
                    )
                except (ValueError, IndexError):
                    print(
                        "[NOTE: input BED must be 12-column] "
                        f"skipped line {line_number}: {line.rstrip()}",
                        file=sys.stderr,
                    )
                    continue

                if chrom not in chromosomes:
                    continue

                gene_count += 1
                gene_all_bases: list[int] = []
                skip_gene = False

                for start, end in zip(exon_starts, exon_ends):
                    # Preserve the original 1-based genomic coordinates.
                    gene_all_bases.extend(range(start + 1, end + 1))

                    # Preserve the original early length check.
                    if len(gene_all_bases) < 100:
                        skip_gene = True
                        break

                if skip_gene:
                    continue

                gene_all_bases.sort(reverse=(strand == "-"))
                percentile_bases = mystat.percentile_list(gene_all_bases)

                for index, genomic_position in enumerate(percentile_bases):
                    signal = bigwig.values(
                        chrom,
                        genomic_position - 1,
                        genomic_position,
                    )
                    coverage[index] += float(np.nan_to_num(signal[0]))

                print(
                    f"\t{gene_count} genes finished\r",
                    end=" ",
                    file=sys.stderr,
                )

    finally:
        bigwig.close()

    print(file=sys.stderr)

    data_file = Path(f"{output_prefix}.geneBodyCoverage.txt")
    r_script = Path(f"{output_prefix}.geneBodyCoverage_plot.r")

    y_coordinates: list[str] = []

    with data_file.open("w", encoding="utf-8") as output:
        output.write("percentile\tcount\n")

        for index in coverage:
            value = coverage[index]
            y_coordinates.append(str(value))
            output.write(f"{index}\t{value}\n")

    with r_script.open("w", encoding="utf-8") as output:
        plot_file = f"{output_prefix}.geneBodyCoverage.{graph_type}"
        output.write(f"{graph_type}('{plot_file}')\n")
        output.write("x=1:100\n")
        output.write(f"y=c({','.join(y_coordinates)})\n")
        output.write(
            "plot("
            f"x, y/{gene_count}, "
            "xlab=\"percentile of gene body (5'->3')\", "
            "ylab='average wigsum', "
            "type='s'"
            ")\n"
        )
        output.write("dev.off()\n")

    return data_file, r_script, gene_count


def run_r_script(
    parser: argparse.ArgumentParser,
    script_path: Path,
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
            [resolved, str(script_path)],
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        parser.exit(
            exc.returncode or 1,
            f"{parser.prog}: error: R plotting failed for {script_path}\n",
        )


def main(argv: Sequence[str] | None = None) -> int:
    """Run BigWig-based gene-body coverage analysis."""
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(parser, args)

    try:
        data_file, r_script, gene_count = coverage_gene_body_bigwig(
            bigwig_file=args.input_file,
            refbed=args.ref_gene_model,
            output_prefix=args.output_prefix,
            graph_type=args.graph_type,
        )

        if gene_count == 0:
            raise ValueError(
                "no valid BED12 records matched chromosomes in the BigWig"
            )

        if not args.skip_plot:
            run_r_script(
                parser=parser,
                script_path=r_script,
                executable=args.rscript,
            )

        print(f"Created: {data_file}", file=sys.stderr)
        print(f"Created: {r_script}", file=sys.stderr)

    except (OSError, ValueError, RuntimeError, ZeroDivisionError) as exc:
        parser.exit(1, f"{parser.prog}: error: {exc}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
