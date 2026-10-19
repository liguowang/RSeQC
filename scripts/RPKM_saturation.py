#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Assess whether transcript RPKM estimates have reached sequencing saturation.

The original resampling analysis, percent-relative-error calculation, quartile
grouping, and R plotting logic are preserved.
"""

from __future__ import annotations

import argparse
import collections
import operator
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Sequence
import numpy as np
from qcmodule import SAM
from importlib.metadata import version

__version__ = version("RSeQC")
__author__ = "Liguo Wang"


DESCRIPTION = (
    "Estimate transcript-level RPKM saturation by repeatedly subsampling "
    "mapped reads and comparing each estimate with the full-depth value."
)

EPILOG = """
Example
-------
RPKM_saturation.py \
    -i sample.bam \
    -r genes.bed12 \
    -o sample \
    -l 5 \
    -u 100 \
    -s 5 \
    -c 0.01 \
    -q 30

Strand-specific example
-----------------------
RPKM_saturation.py \
    -i sample.bam \
    -r genes.bed12 \
    -o sample \
    -d '1++,1--,2+-,2-+'

Notes
-----
* Use infer_experiment.py when the strand rule is unknown.
* Transcripts with mean RPKM below --rpkm-cutoff are omitted from the plot.
* The generated boxplots show percent relative error by expression quartile.
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
        "-o",
        "--out-prefix",
        dest="output_prefix",
        required=True,
        type=Path,
        help="Prefix for generated output files.",
    )
    parser.add_argument(
        "-r",
        "--refgene",
        dest="refgene_bed",
        required=True,
        type=Path,
        help="Reference gene model in BED format.",
    )
    parser.add_argument(
        "-d",
        "--strand",
        dest="strand_rule",
        default=None,
        help=(
            "Sequencing strand rule, for example "
            "'1++,1--,2+-,2-+'. Omit for unstranded data."
        ),
    )
    parser.add_argument(
        "-l",
        "--percentile-floor",
        dest="percentile_low_bound",
        type=int,
        default=5,
        metavar="INT",
        help="Starting sampling percentile. Default: %(default)s",
    )
    parser.add_argument(
        "-u",
        "--percentile-ceiling",
        dest="percentile_up_bound",
        type=int,
        default=100,
        metavar="INT",
        help="Ending sampling percentile. Default: %(default)s",
    )
    parser.add_argument(
        "-s",
        "--percentile-step",
        dest="percentile_step",
        type=int,
        default=5,
        metavar="INT",
        help=(
            "Sampling percentile step. Smaller values perform more "
            "resampling iterations. Default: %(default)s"
        ),
    )
    parser.add_argument(
        "-c",
        "--rpkm-cutoff",
        dest="rpkm_cutoff",
        type=float,
        default=0.01,
        metavar="FLOAT",
        help=(
            "Omit transcripts with mean RPKM below this value from the plot. "
            "Default: %(default)s"
        ),
    )
    parser.add_argument(
        "-q",
        "--mapq",
        dest="map_qual",
        type=int,
        default=30,
        metavar="INT",
        help=(
            "Minimum mapping quality for an alignment to be treated as "
            "uniquely mapped. Default: %(default)s"
        ),
    )
    parser.add_argument(
        "--skip-plot",
        action="store_true",
        help="Generate saturation data and R script without executing R.",
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
    if not args.input_file.is_file():
        parser.error(
            f"input alignment file does not exist or is not a file: "
            f"{args.input_file}"
        )

    if args.input_file.suffix.lower() not in {".bam", ".sam", ".cram"}:
        parser.error(
            "input alignment file must have a .bam, .sam, or .cram "
            f"extension: {args.input_file}"
        )

    if not args.refgene_bed.is_file():
        parser.error(
            f"reference BED file does not exist or is not a file: "
            f"{args.refgene_bed}"
        )

    if not 0 <= args.percentile_low_bound <= 100:
        parser.error("--percentile-floor must be between 0 and 100")

    if not 0 <= args.percentile_up_bound <= 100:
        parser.error("--percentile-ceiling must be between 0 and 100")

    if args.percentile_up_bound < args.percentile_low_bound:
        parser.error(
            "--percentile-ceiling must be greater than or equal to "
            "--percentile-floor"
        )

    if args.percentile_step <= 0:
        parser.error("--percentile-step must be greater than zero")

    if args.percentile_step > args.percentile_up_bound:
        parser.error(
            "--percentile-step cannot be greater than "
            "--percentile-ceiling"
        )

    if args.rpkm_cutoff < 0:
        parser.error("--rpkm-cutoff must be zero or greater")

    if args.map_qual < 0:
        parser.error("--mapq must be zero or greater")

    output_parent = args.output_prefix.parent
    if not output_parent.exists():
        parser.error(f"output directory does not exist: {output_parent}")

    if not output_parent.is_dir():
        parser.error(f"output parent is not a directory: {output_parent}")


def square_error(values: Sequence[float]) -> list[float] | None:
    """Calculate relative error against the final, full-depth RPKM value."""
    true_rpkm = values[-1]
    value_range = max(values) - min(values)

    if true_rpkm == 0 or value_range == 0:
        return None

    return [
        abs(value - true_rpkm) / true_rpkm
        for value in values
    ]


def parse_rpkm_table(
    infile: Path,
    rpkm_cutoff: float,
) -> tuple[
    list[str],
    dict[str, list[float]],
    dict[str, float],
]:
    """Parse the saturation table using the original filtering rules."""
    rpkm_errors: dict[str, list[float]] = {}
    rpkm_means: dict[str, float] = {}
    header: list[str] | None = None

    with infile.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()

            if not stripped:
                continue

            fields = stripped.split()

            if fields[0].startswith("#"):
                if len(fields) < 7:
                    raise ValueError(
                        f"header on line {line_number} has fewer than 7 fields"
                    )
                header = [value.replace("%", "") for value in fields[6:]]
                continue

            if header is None:
                raise ValueError("no header line was found in the RPKM table")

            if len(fields) < 7:
                raise ValueError(
                    f"data line {line_number} has fewer than 7 fields"
                )

            key = "\t".join(fields[:6])
            values = [float(value) for value in fields[6:]]

            if len(values) != len(header):
                raise ValueError(
                    f"line {line_number} contains {len(values)} RPKM values, "
                    f"but the header contains {len(header)} sampling points"
                )

            if max(values) == 0:
                continue

            if max(values) - min(values) == 0:
                continue

            if float(np.mean(values)) < rpkm_cutoff:
                continue

            errors = square_error(values)
            if errors is None:
                continue

            rpkm_errors[key] = errors
            rpkm_means[key] = float(np.mean(values))

    if header is None:
        raise ValueError("no header line was found in the RPKM table")

    return header, rpkm_errors, rpkm_means


def r_vector(values: Sequence[str]) -> str:
    """Format values as an R vector."""
    return ",".join(values)


def show_saturation(
    infile: Path,
    outfile: Path,
    rpkm_cutoff: float = 0.01,
) -> None:
    """Generate the original quartile-based R saturation plot script."""
    header, rpkm_errors, rpkm_means = parse_rpkm_table(
        infile,
        rpkm_cutoff,
    )
    gene_count = len(rpkm_means)

    if gene_count == 0:
        raise ValueError(
            "no transcripts passed the saturation-plot filtering criteria"
        )

    quantiles = {
        "Q1": (0.00, 0.25),
        "Q2": (0.25, 0.50),
        "Q3": (0.50, 0.75),
        "Q4": (0.75, 1.00),
    }

    pdf_path = outfile.with_suffix(".pdf")

    with outfile.open("w", encoding="utf-8") as rout:
        print(f"pdf({str(pdf_path)!r})", file=rout)
        print("par(mfrow=c(2,2))", file=rout)

        sorted_genes = sorted(
            rpkm_means.items(),
            key=operator.itemgetter(1),
        )

        for quantile in sorted(quantiles):
            lower_bound, upper_bound = quantiles[quantile]
            normalized_rpkm: collections.defaultdict[
                str, list[str]
            ] = collections.defaultdict(list)

            for line_count, (gene_key, _) in enumerate(
                sorted_genes,
                start=1,
            ):
                if (
                    line_count > gene_count * lower_bound
                    and line_count <= gene_count * upper_bound
                ):
                    for index, error in enumerate(rpkm_errors[gene_key]):
                        normalized_rpkm[header[index]].append(str(error))

            plot_headers = header[:-1]

            if not plot_headers:
                raise ValueError(
                    "at least two sampling percentages are required"
                )

            print(
                "name=c(%s)"
                % r_vector([repr(value) for value in plot_headers]),
                file=rout,
            )

            for sample in plot_headers:
                print(
                    f"S{sample}=c({r_vector(normalized_rpkm[sample])})",
                    file=rout,
                )

            series = ",".join(f"100*S{sample}" for sample in plot_headers)
            print(
                "boxplot("
                f"{series},"
                "names=name,"
                "outline=FALSE,"
                "ylab='Percent Relative Error',"
                f"main={quantile!r},"
                "xlab='Resampling percentage'"
                ")",
                file=rout,
            )

        print("dev.off()", file=rout)


def run_plot_script(
    parser: argparse.ArgumentParser,
    script_path: Path,
    rscript_executable: str,
) -> None:
    """Execute the generated R plotting script safely."""
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
    """Run RPKM-saturation analysis."""
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(parser, args)

    rpkm_table = Path(f"{args.output_prefix}.eRPKM.xls")
    plot_script = Path(f"{args.output_prefix}.saturation.r")

    try:
        alignment = SAM.ParseBAM(str(args.input_file))
        alignment.saturation_RPKM(
            outfile=str(args.output_prefix),
            refbed=str(args.refgene_bed),
            sample_start=args.percentile_low_bound,
            sample_end=args.percentile_up_bound,
            sample_step=args.percentile_step,
            strand_rule=args.strand_rule,
            q_cut=args.map_qual,
        )

        if not rpkm_table.is_file():
            raise OSError(
                f"expected saturation table was not created: {rpkm_table}"
            )

        show_saturation(
            infile=rpkm_table,
            outfile=plot_script,
            rpkm_cutoff=args.rpkm_cutoff,
        )

    except (OSError, ValueError, RuntimeError, IndexError) as exc:
        parser.exit(1, f"{parser.prog}: error: {exc}\n")

    if not args.skip_plot:
        run_plot_script(
            parser=parser,
            script_path=plot_script,
            rscript_executable=args.rscript,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
