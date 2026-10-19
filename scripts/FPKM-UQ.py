#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Calculate raw counts, FPKM, and upper-quartile normalized FPKM.

The FPKM and FPKM-UQ calculations intentionally preserve the original
RSeQC/TCGA-compatible algorithm.
"""

from __future__ import annotations
import argparse
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Sequence
import numpy as np
from importlib.metadata import version

__version__ = version("RSeQC")
__author__ = "Liguo Wang"


def printlog(message: str) -> None:
    """Write a timestamped progress message to standard error."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"@ {timestamp}: {message}", file=sys.stderr)


def build_parser() -> argparse.ArgumentParser:
    """Create and configure the command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--bam",
        dest="bam_file",
        required=True,
        type=Path,
        help=(
            "Coordinate-sorted BAM file. For TCGA-compatible results, use "
            "alignments produced by the corresponding TCGA RNA-seq workflow."
        ),
    )
    parser.add_argument(
        "--gtf",
        dest="gtf_file",
        required=True,
        type=Path,
        help="Gene model in GTF format.",
    )
    parser.add_argument(
        "--info",
        dest="info_file",
        required=True,
        type=Path,
        help="Gene information file containing exon lengths and gene types.",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="out_prefix",
        required=True,
        type=Path,
        help="Prefix for output files.",
    )
    parser.add_argument(
        "--log2",
        dest="log_scale",
        action="store_true",
        help="Report log2(FPKM + 1) and log2(FPKM-UQ + 1).",
    )
    parser.add_argument(
        "--htseq-count",
        dest="htseq_count",
        default="htseq-count",
        metavar="PATH",
        help="htseq-count executable to use. Default: %(default)s",
    )
    parser.add_argument(
        "--print-htseq-command",
        action="store_true",
        help="Print the htseq-count command and exit without running it.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Validate input files and output location."""
    if not args.bam_file.is_file():
        parser.error(f"BAM file does not exist: {args.bam_file}")

    if args.bam_file.suffix.lower() != ".bam":
        parser.error(f"input alignment file must be BAM: {args.bam_file}")

    if not args.gtf_file.is_file():
        parser.error(f"GTF file does not exist: {args.gtf_file}")

    if not args.info_file.is_file():
        parser.error(f"gene information file does not exist: {args.info_file}")

    output_parent = args.out_prefix.parent
    if not output_parent.exists():
        parser.error(f"output directory does not exist: {output_parent}")

    if not output_parent.is_dir():
        parser.error(f"output parent is not a directory: {output_parent}")


def resolve_executable(executable: str) -> str:
    """Resolve an executable name or path."""
    resolved = shutil.which(executable)
    if resolved is None:
        raise FileNotFoundError(
            f'cannot find htseq-count executable "{executable}"'
        )
    return resolved


def htseq_command(
    executable: str,
    bam_file: Path,
    gtf_file: Path,
) -> list[str]:
    """Build the original TCGA-compatible htseq-count command."""
    return [
        executable,
        "-f",
        "bam",
        "-r",
        "pos",
        "-s",
        "no",
        "-a",
        "10",
        "-t",
        "exon",
        "-i",
        "gene_id",
        "-m",
        "intersection-nonempty",
        str(bam_file.resolve()),
        str(gtf_file.resolve()),
    ]


def format_command(command: Sequence[str]) -> str:
    """Return a shell-readable command for display only."""
    import shlex

    return " ".join(shlex.quote(part) for part in command)


def run_htseq(
    bam_file: Path,
    gtf_file: Path,
    out_file: Path,
    *,
    executable: str = "htseq-count",
    print_command_only: bool = False,
) -> None:
    """Run htseq-count using the original parameter set."""
    resolved = resolve_executable(executable)
    command = htseq_command(resolved, bam_file, gtf_file)

    if print_command_only:
        print(format_command(command))
        return

    printlog(f"Running: {format_command(command)}")

    with out_file.open("w", encoding="utf-8") as output:
        try:
            subprocess.run(
                command,
                stdout=output,
                check=True,
            )
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                f"htseq-count failed with exit status {exc.returncode}"
            ) from exc


def read_gene_information(
    info_file: Path,
) -> tuple[dict[str, int], dict[str, str], set[str]]:
    """Read gene sizes, display information, and protein-coding gene IDs."""
    printlog(f"Read gene information file: {info_file}")

    gene_sizes: dict[str, int] = {}
    gene_info: dict[str, str] = {}
    protein_coding: set[str] = set()

    with info_file.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()

            if not stripped or stripped.startswith("gene_id"):
                continue

            fields = stripped.split()
            if len(fields) < 11:
                raise ValueError(
                    f"{info_file}:{line_number}: expected at least 11 columns"
                )

            gene_id = fields[0]

            try:
                exon_length = int(fields[10])
            except ValueError as exc:
                raise ValueError(
                    f"{info_file}:{line_number}: exon_length must be an integer"
                ) from exc

            gene_sizes[gene_id] = exon_length
            gene_info[gene_id] = "\t".join(fields[1:6])

            if fields[6] == "protein_coding":
                protein_coding.add(gene_id)

    print(f"\tTotal genes: {len(gene_sizes)}", file=sys.stderr)
    print(
        f"\tTotal protein-coding genes: {len(protein_coding)}",
        file=sys.stderr,
    )

    return gene_sizes, gene_info, protein_coding


def read_htseq_counts(count_file: Path) -> list[tuple[str, int]]:
    """Read non-summary records from an htseq-count output file."""
    counts: list[tuple[str, int]] = []

    with count_file.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()

            if not stripped or stripped.startswith("__"):
                continue

            fields = stripped.split()
            if len(fields) < 2:
                raise ValueError(
                    f"{count_file}:{line_number}: expected two columns"
                )

            try:
                count = int(fields[1])
            except ValueError as exc:
                raise ValueError(
                    f"{count_file}:{line_number}: count must be an integer"
                ) from exc

            counts.append((fields[0], count))

    return counts


def calculate_fpkm(
    count_file: Path,
    info_file: Path,
    out_file: Path,
    *,
    log2_flag: bool = False,
) -> None:
    """Calculate FPKM and FPKM-UQ with the original algorithm."""
    gene_sizes, gene_info, protein_coding = read_gene_information(info_file)

    printlog(
        "Read gene count file to calculate 75 percentile count and total "
        f"count: {count_file}"
    )
    all_counts = read_htseq_counts(count_file)

    # Preserve the original algorithm: calculate both normalization values
    # from counts belonging only to protein-coding genes.
    protein_coding_counts = [
        count
        for gene_id, count in all_counts
        if gene_id in protein_coding
    ]

    if not protein_coding_counts:
        raise ValueError("no protein-coding gene counts were found")

    uq_count = np.percentile(sorted(protein_coding_counts), 75)
    total_count = sum(protein_coding_counts)

    print(
        f"\tTotal protein-coding genes: {len(protein_coding_counts)}",
        file=sys.stderr,
    )
    print(
        "\tThe 75 percentile count of protein-coding genes: "
        f"{uq_count:f}",
        file=sys.stderr,
    )
    print(
        f"\tThe total count of protein-coding genes: {total_count:f}",
        file=sys.stderr,
    )

    print(
        f"Read gene count file to calculate FPKM and FPKM-UQ: {count_file}",
        file=sys.stderr,
    )

    with out_file.open("w", encoding="utf-8") as output:
        if log2_flag:
            header = (
                "gene_ID",
                "symbol",
                "chrom",
                "start",
                "end",
                "strand",
                "raw_count",
                "FPKM(log2(x+1))",
                "FPKM-UQ(log2(x+1))",
            )
        else:
            header = (
                "gene_ID",
                "symbol",
                "chrom",
                "start",
                "end",
                "strand",
                "raw_count",
                "FPKM",
                "FPKM-UQ",
            )

        print("\t".join(header), file=output)

        for gene_id, count in all_counts:
            if gene_id not in gene_sizes:
                # The original formula cannot be evaluated without exon length
                # and display metadata. Skip safely rather than raising a
                # secondary KeyError.
                print(
                    f"Warning: {gene_id} is absent from {info_file}; skipped",
                    file=sys.stderr,
                )
                continue

            gene_size = gene_sizes[gene_id]

            try:
                if log2_flag:
                    fpkm_uq = np.log2(
                        (count * 1_000_000_000)
                        / (gene_size * uq_count)
                        + 1
                    )
                    fpkm = np.log2(
                        (count * 1_000_000_000)
                        / (gene_size * total_count)
                        + 1
                    )
                else:
                    fpkm_uq = (
                        count * 1_000_000_000
                        / (gene_size * uq_count)
                    )
                    fpkm = (
                        count * 1_000_000_000
                        / (gene_size * total_count)
                    )
            except (ZeroDivisionError, FloatingPointError):
                fpkm_uq = "NA"
                fpkm = "NA"

            print(
                "\t".join(
                    (
                        gene_id,
                        gene_info[gene_id],
                        str(count),
                        str(fpkm),
                        str(fpkm_uq),
                    )
                ),
                file=output,
            )


def main(argv: Sequence[str] | None = None) -> int:
    """Run htseq-count and calculate FPKM/FPKM-UQ."""
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(parser, args)

    count_file = Path(f"{args.out_prefix}.htseq.counts.txt")
    fpkm_file = Path(f"{args.out_prefix}.FPKM-UQ.txt")

    try:
        if args.print_htseq_command:
            run_htseq(
                bam_file=args.bam_file,
                gtf_file=args.gtf_file,
                out_file=count_file,
                executable=args.htseq_count,
                print_command_only=True,
            )
            return 0

        printlog("Running htseq-count ...")
        run_htseq(
            bam_file=args.bam_file,
            gtf_file=args.gtf_file,
            out_file=count_file,
            executable=args.htseq_count,
        )

        if args.log_scale:
            printlog("Calculate log2(FPKM + 1) and log2(FPKM-UQ + 1) ...")
        else:
            printlog("Calculate FPKM and FPKM-UQ ...")

        calculate_fpkm(
            count_file=count_file,
            info_file=args.info_file,
            out_file=fpkm_file,
            log2_flag=args.log_scale,
        )

        printlog(f"Created: {count_file}")
        printlog(f"Created: {fpkm_file}")

    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(1, f"{parser.prog}: error: {exc}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
