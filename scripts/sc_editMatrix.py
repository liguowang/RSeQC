#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Visualize error-correction edits in cellular barcodes and UMIs.

The command extracts raw-versus-corrected barcode and UMI edits from a tagged
single-cell RNA-seq BAM file, writes edit-count matrices, and optionally
generates heatmaps showing edit position, substitution type, and frequency.
"""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Sequence
from qcmodule import heatmap
from qcmodule import scbam
from importlib.metadata import version

__version__ = version("RSeQC")
__author__ = "Liguo Wang"


DESCRIPTION = (
    "Generate edit-count matrices and heatmaps for error-corrected cellular "
    "barcodes and UMIs in a tagged single-cell RNA-seq BAM file."
)

EPILOG = """
Default BAM tags
----------------
CR
    Cellular barcode reported by the sequencer.

CB
    Error-corrected cellular barcode.

UR
    UMI reported by the sequencer.

UB
    Error-corrected UMI.

Examples
--------
Generate matrices and heatmaps:
    sc_editMatrix.py -i possorted_genome_bam.bam -o sample

Generate matrices only:
    sc_editMatrix.py -i possorted_genome_bam.bam -o sample --skip-heatmap

Install the required R package automatically:
    sc_editMatrix.py \
        -i possorted_genome_bam.bam \
        -o sample \
        --install-r-deps
"""


SUPPORTED_FILE_TYPES = ("pdf", "png", "tiff", "bmp", "jpeg")
SUPPORTED_ANGLES = (0, 45, 90, 270, 315)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=DESCRIPTION,
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-i", "--infile", dest="in_file", required=True, type=Path,
        help="Input BAM file.",
    )
    parser.add_argument(
        "-o", "--outfile", dest="out_file", required=True, type=Path,
        help="Prefix for generated matrix and heatmap files.",
    )
    parser.add_argument(
        "--limit", dest="reads_num", type=int, default=None, metavar="INT",
        help="Maximum number of alignments to process. Default: all",
    )
    parser.add_argument(
        "--cr-tag", dest="cr_tag", default="CR", metavar="TAG",
        help="Sequencer-reported cellular-barcode tag. Default: %(default)s",
    )
    parser.add_argument(
        "--cb-tag", dest="cb_tag", default="CB", metavar="TAG",
        help="Error-corrected cellular-barcode tag. Default: %(default)s",
    )
    parser.add_argument(
        "--ur-tag", dest="ur_tag", default="UR", metavar="TAG",
        help="Sequencer-reported UMI tag. Default: %(default)s",
    )
    parser.add_argument(
        "--ub-tag", dest="ub_tag", default="UB", metavar="TAG",
        help="Error-corrected UMI tag. Default: %(default)s",
    )
    parser.add_argument(
        "--cell-width", type=int, default=15, metavar="INT",
        help="Heatmap cell width in points. Default: %(default)s",
    )
    parser.add_argument(
        "--cell-height", type=int, default=10, metavar="INT",
        help="Heatmap cell height in points. Default: %(default)s",
    )
    parser.add_argument(
        "--font-size", type=int, default=8, metavar="INT",
        help="Heatmap font size. Default: %(default)s",
    )
    parser.add_argument(
        "--angle", dest="col_angle", type=int, choices=SUPPORTED_ANGLES,
        default=45, help="Column-label angle. Default: %(default)s",
    )
    parser.add_argument(
        "--text-color", default="black", metavar="COLOR",
        help="Color of numeric labels in heatmap cells. Default: %(default)s",
    )
    parser.add_argument(
        "--file-type", type=str.lower, choices=SUPPORTED_FILE_TYPES,
        default="pdf", help="Heatmap output format. Default: %(default)s",
    )
    parser.add_argument(
        "--no-num", action="store_true",
        help="Do not print numeric values inside heatmap cells.",
    )
    parser.add_argument(
        "--skip-heatmap", action="store_true",
        help="Generate edit-count matrices without creating heatmaps.",
    )
    parser.add_argument(
        "--rscript", default="Rscript", metavar="PATH",
        help="Rscript executable used for dependency checks. Default: %(default)s",
    )
    parser.add_argument(
        "--install-r-deps", action="store_true",
        help="Install the R package pheatmap when it is missing.",
    )
    parser.add_argument(
        "--cran-mirror",
        default="https://cloud.r-project.org",
        metavar="URL",
        help="CRAN mirror used with --install-r-deps. Default: %(default)s",
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Enable detailed progress logging.",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}",
    )
    return parser


def configure_logging(verbose: bool) -> None:
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.DEBUG if verbose else logging.INFO,
    )


def validate_tag(parser: argparse.ArgumentParser, option_name: str, tag: str) -> None:
    if len(tag) != 2:
        parser.error(f"{option_name} must be exactly two characters")
    if any(character.isspace() for character in tag):
        parser.error(f"{option_name} cannot contain whitespace")


def find_bam_index(bam_file: Path) -> Path | None:
    candidates = (Path(f"{bam_file}.bai"), bam_file.with_suffix(".bai"))
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


def validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    if not args.in_file.is_file():
        parser.error(f"input BAM file does not exist: {args.in_file}")
    if args.in_file.suffix.lower() != ".bam":
        parser.error(f"input file must have a .bam extension: {args.in_file}")
    if args.reads_num is not None and args.reads_num <= 0:
        parser.error("--limit must be greater than zero")

    for option_name, tag in (
        ("--cr-tag", args.cr_tag),
        ("--cb-tag", args.cb_tag),
        ("--ur-tag", args.ur_tag),
        ("--ub-tag", args.ub_tag),
    ):
        validate_tag(parser, option_name, tag)

    if args.cell_width <= 0:
        parser.error("--cell-width must be greater than zero")
    if args.cell_height <= 0:
        parser.error("--cell-height must be greater than zero")
    if args.font_size <= 0:
        parser.error("--font-size must be greater than zero")
    if not args.text_color.strip():
        parser.error("--text-color cannot be empty")
    if not args.cran_mirror.strip():
        parser.error("--cran-mirror cannot be empty")

    output_parent = args.out_file.parent
    if not output_parent.exists():
        parser.error(f"output directory does not exist: {output_parent}")
    if not output_parent.is_dir():
        parser.error(f"output parent is not a directory: {output_parent}")

    bam_index = find_bam_index(args.in_file)
    if bam_index is not None:
        try:
            if bam_index.stat().st_mtime < args.in_file.stat().st_mtime:
                logging.warning(
                    "BAM index is older than the BAM file: %s. "
                    "Rebuild it with 'samtools index -f %s'.",
                    bam_index,
                    args.in_file,
                )
        except OSError:
            pass


def resolve_rscript(executable: str) -> str:
    resolved = shutil.which(executable)
    if resolved is None:
        raise OSError(f"Rscript executable not found: {executable}")
    return resolved


def r_package_available(rscript: str, package: str) -> bool:
    result = subprocess.run(
        [
            rscript,
            "-e",
            (
                f"quit(status=ifelse(requireNamespace('{package}', "
                "quietly=TRUE), 0, 1))"
            ),
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return result.returncode == 0


def install_r_package(
    rscript: str,
    package: str,
    cran_mirror: str,
) -> None:
    logging.info("Installing R package %s from %s", package, cran_mirror)
    subprocess.run(
        [
            rscript,
            "-e",
            (
                f"install.packages('{package}', "
                f"repos={cran_mirror!r})"
            ),
        ],
        check=True,
    )


def ensure_r_dependencies(
    rscript_executable: str,
    install_missing: bool,
    cran_mirror: str,
) -> None:
    rscript = resolve_rscript(rscript_executable)

    if r_package_available(rscript, "pheatmap"):
        return

    if not install_missing:
        raise RuntimeError(
            "R package 'pheatmap' is not installed. Install it with "
            "\"Rscript -e \\\"install.packages('pheatmap', "
            "repos='https://cloud.r-project.org')\\\"\", rerun with "
            "--install-r-deps, or use --skip-heatmap."
        )

    install_r_package(rscript, "pheatmap", cran_mirror)

    if not r_package_available(rscript, "pheatmap"):
        raise RuntimeError("R package 'pheatmap' is still unavailable after installation")


def expected_matrix_paths(output_prefix: Path) -> tuple[Path, Path]:
    return (
        Path(f"{output_prefix}.CB_edits_count.csv"),
        Path(f"{output_prefix}.UMI_edits_count.csv"),
    )


def expected_heatmap_path(output_prefix: Path, file_type: str) -> Path:
    return Path(f"{output_prefix}.{file_type}")


def generate_heatmap(
    matrix_file: Path,
    output_prefix: Path,
    *,
    file_type: str,
    cell_width: int,
    cell_height: int,
    col_angle: int,
    font_size: int,
    text_color: str,
    no_numbers: bool,
) -> Path:
    if not matrix_file.is_file():
        raise OSError(f"expected edit-count matrix was not created: {matrix_file}")

    output_path = expected_heatmap_path(output_prefix, file_type)
    if output_path.exists():
        output_path.unlink()

    heatmap.make_heatmap(
        infile=str(matrix_file),
        outfile=str(output_prefix),
        filetype=file_type,
        cell_width=cell_width,
        cell_height=cell_height,
        col_angle=col_angle,
        font_size=font_size,
        text_color=text_color,
        no_numbers=no_numbers,
        log2_scale=True,
    )

    if not output_path.is_file():
        raise RuntimeError(
            f"heatmap generation failed; expected output was not created: {output_path}"
        )

    return output_path


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    configure_logging(args.verbose)
    validate_args(parser, args)

    try:
        scbam.barcode_edits(
            infile=str(args.in_file),
            outfile=str(args.out_file),
            limit=args.reads_num,
            CR_tag=args.cr_tag,
            CB_tag=args.cb_tag,
            UR_tag=args.ur_tag,
            UB_tag=args.ub_tag,
        )

        cb_matrix, umi_matrix = expected_matrix_paths(args.out_file)

        if not cb_matrix.is_file():
            raise OSError(
                f"expected cellular-barcode edit matrix was not created: {cb_matrix}"
            )
        if not umi_matrix.is_file():
            raise OSError(f"expected UMI edit matrix was not created: {umi_matrix}")

        logging.info("Created %s", cb_matrix)
        logging.info("Created %s", umi_matrix)

        if not args.skip_heatmap:
            ensure_r_dependencies(
                rscript_executable=args.rscript,
                install_missing=args.install_r_deps,
                cran_mirror=args.cran_mirror,
            )

            cb_heatmap = generate_heatmap(
                cb_matrix,
                Path(f"{args.out_file}.CB_edits_heatmap"),
                file_type=args.file_type,
                cell_width=args.cell_width,
                cell_height=args.cell_height,
                col_angle=args.col_angle,
                font_size=args.font_size,
                text_color=args.text_color,
                no_numbers=args.no_num,
            )
            umi_heatmap = generate_heatmap(
                umi_matrix,
                Path(f"{args.out_file}.UMI_edits_heatmap"),
                file_type=args.file_type,
                cell_width=args.cell_width,
                cell_height=args.cell_height,
                col_angle=args.col_angle,
                font_size=args.font_size,
                text_color=args.text_color,
                no_numbers=args.no_num,
            )

            logging.info("Created %s", cb_heatmap)
            logging.info("Created %s", umi_heatmap)

    except subprocess.CalledProcessError as exc:
        logging.error("R dependency installation failed with exit code %s", exc.returncode)
        return exc.returncode or 1
    except (OSError, ValueError, RuntimeError, KeyError) as exc:
        logging.error("%s", exc)
        return 1

    logging.info("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
