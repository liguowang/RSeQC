#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Annotate splice junctions against a reference gene model."""

from __future__ import annotations
import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Sequence
from qcmodule import SAM


__author__ = "Liguo Wang"
__version__ = "5.05"


DESCRIPTION = (
    "Annotate RNA splicings at both read and consolidated junction levels "
    "using a gene model in BED format."
)

EPILOG = """
Notes
-----
* A long read may contain multiple splice junctions.
* Multiple reads spanning the same intron are consolidated into one junction.
* The BED and Interact files are generated from <prefix>.junction.xls.

Example
-------
junction_annotation.py -i sample.bam -r genes.bed12 -o sample -m 50 -q 30
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
        help=(
            "Reference gene model in BED format. A pooled annotation is "
            "recommended for junction classification."
        ),
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
        "-m",
        "--min-intron",
        dest="min_intron",
        type=int,
        default=50,
        metavar="INT",
        help="Minimum intron length in bp. Default: %(default)s",
    )
    parser.add_argument(
        "-q",
        "--mapq",
        dest="map_qual",
        type=int,
        default=30,
        metavar="INT",
        help=(
            "Minimum mapping quality for an alignment to be considered "
            "uniquely mapped. Default: %(default)s"
        ),
    )
    parser.add_argument(
        "--skip-plot",
        action="store_true",
        help="Generate annotation files but do not execute the R plotting script.",
    )
    parser.add_argument(
        "--rscript",
        default="Rscript",
        metavar="PATH",
        help="Rscript executable to use. Default: %(default)s",
    )
    parser.add_argument(
        "--skip-bed",
        action="store_true",
        help="Do not generate the BED12 junction file.",
    )
    parser.add_argument(
        "--skip-interact",
        action="store_true",
        help="Do not generate the UCSC Interact junction file.",
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
    output_prefix: Path,
    min_intron: int,
    map_qual: int,
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

    if min_intron <= 0:
        parser.error("--min-intron must be greater than zero")

    if map_qual < 0:
        parser.error("--mapq must be zero or greater")

    output_parent = output_prefix.parent
    if not output_parent.exists():
        parser.error(f"output directory does not exist: {output_parent}")

    if not output_parent.is_dir():
        parser.error(f"output parent is not a directory: {output_parent}")


def junction_color(annotation: str) -> str:
    """Return the legacy RGB color for a junction class."""
    return {
        "annotated": "205,0,0",
        "partial_novel": "0,205,0",
        "complete_novel": "0,0,205",
    }.get(annotation, "0,0,0")


def iter_junction_records(infile: Path):
    """Yield parsed records from the junction summary file."""
    with infile.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()

            if not stripped or stripped.startswith("chrom"):
                continue

            fields = stripped.split()
            if len(fields) != 5:
                print(
                    f"Warning: skipped malformed junction record at "
                    f"{infile}:{line_number}",
                    file=sys.stderr,
                )
                continue

            chrom, start_text, end_text, score_text, annotation = fields

            try:
                start = int(start_text)
                end = int(end_text)
                score = int(score_text)
            except ValueError:
                print(
                    f"Warning: skipped non-numeric junction record at "
                    f"{infile}:{line_number}",
                    file=sys.stderr,
                )
                continue

            yield chrom, start, end, score, annotation


def generate_bed12(infile: Path, size: int = 1) -> Path:
    """Generate a BED12 file from the junction summary."""
    if size <= 0:
        raise ValueError("BED block size must be greater than zero")

    outfile = infile.with_suffix(".bed")

    with outfile.open("w", encoding="utf-8") as output:
        for chrom, start_value, end_value, score, annotation in iter_junction_records(infile):
            start = start_value - size
            end = end_value + size
            thick_start = start
            thick_end = end
            block_count = 2
            block_sizes = f"{size},{size}"
            block_starts = f"0,{end - size - start}"

            values = (
                chrom,
                start,
                end,
                annotation,
                score,
                ".",
                thick_start,
                thick_end,
                junction_color(annotation),
                block_count,
                block_sizes,
                block_starts,
            )
            print("\t".join(str(value) for value in values), file=output)

    return outfile


def generate_interact(
    infile: Path,
    bam_file: Path,
    size: int = 1,
) -> Path:
    """Generate a UCSC Interact BED file from the junction summary."""
    if size <= 0:
        raise ValueError("Interact anchor size must be greater than zero")

    outfile = infile.with_name(infile.stem + ".Interact.bed")

    with outfile.open("w", encoding="utf-8") as output:
        print(
            'track type=interact name="Splice junctions" '
            f'description="Splice junctions detected from {bam_file}" '
            "maxHeightPixels=200:200:50 visibility=full",
            file=output,
        )

        for chrom, start_value, end_value, score, annotation in iter_junction_records(infile):
            chrom_start = start_value - size
            chrom_end = end_value + size

            name = (
                f"{chrom}:{chrom_start}-{chrom_end}_{annotation}"
            )
            value = float(score)
            experiment = "RNAseq_junction"

            source_start = chrom_start
            source_end = chrom_start + size
            source_name = f"{chrom}:{source_start}-{source_end}"

            target_start = chrom_end - size
            target_end = chrom_end
            target_name = f"{chrom}:{target_start}-{target_end}"

            values = (
                chrom,
                chrom_start,
                chrom_end,
                name,
                score,
                value,
                experiment,
                junction_color(annotation),
                chrom,
                source_start,
                source_end,
                source_name,
                ".",
                chrom,
                target_start,
                target_end,
                target_name,
                ".",
            )
            print("\t".join(str(value) for value in values), file=output)

    return outfile


def run_plot_script(
    parser: argparse.ArgumentParser,
    output_prefix: Path,
    rscript_executable: str,
) -> None:
    """Execute the generated R plotting script safely."""
    script_path = Path(f"{output_prefix}.junction_plot.r")

    if not script_path.is_file():
        parser.exit(
            1,
            f"{parser.prog}: error: expected R script was not created: "
            f"{script_path}\n",
        )

    resolved_rscript = shutil.which(rscript_executable)
    if resolved_rscript is None:
        parser.exit(
            1,
            f"{parser.prog}: error: Rscript executable not found: "
            f"{rscript_executable}\n",
        )

    try:
        subprocess.run(
            [resolved_rscript, str(script_path)],
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        parser.exit(
            exc.returncode or 1,
            f"{parser.prog}: error: R plotting failed for {script_path}\n",
        )


def main(argv: Sequence[str] | None = None) -> int:
    """Run splice-junction annotation."""
    parser = build_parser()
    args = parser.parse_args(argv)

    validate_args(
        parser=parser,
        input_file=args.input_file,
        ref_gene_model=args.ref_gene_model,
        output_prefix=args.output_prefix,
        min_intron=args.min_intron,
        map_qual=args.map_qual,
    )

    try:
        alignment = SAM.ParseBAM(str(args.input_file))
        alignment.annotate_junction(
            outfile=str(args.output_prefix),
            refgene=str(args.ref_gene_model),
            min_intron=args.min_intron,
            q_cut=args.map_qual,
        )

        junction_file = Path(f"{args.output_prefix}.junction.xls")
        if not junction_file.is_file():
            raise FileNotFoundError(
                f"expected junction summary was not created: {junction_file}"
            )

        if not args.skip_plot:
            run_plot_script(
                parser=parser,
                output_prefix=args.output_prefix,
                rscript_executable=args.rscript,
            )

        if not args.skip_bed:
            print("Create BED file ...", file=sys.stderr)
            bed_file = generate_bed12(junction_file)
            print(f"Created: {bed_file}", file=sys.stderr)

        if not args.skip_interact:
            print("Create Interact file ...", file=sys.stderr)
            interact_file = generate_interact(
                junction_file,
                args.input_file,
            )
            print(f"Created: {interact_file}", file=sys.stderr)

    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(1, f"{parser.prog}: error: {exc}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
