#!/usr/bin/env python3
"""Generate heatmaps with the R package pheatmap."""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path

SUPPORTED_FILETYPES = {"png", "pdf", "tiff", "bmp", "jpeg"}
VALID_COLUMN_ANGLES = {0, 45, 90, 270, 315}


def _r_quote(value):
    """Quote a value for use as an R single-quoted string."""
    value = str(value).replace("\\", "\\\\").replace("'", "\\'")
    return "'" + value + "'"


def make_heatmap(
    infile,
    outfile,
    filetype,
    cell_width,
    cell_height,
    col_angle,
    font_size,
    text_color,
    no_numbers=False,
    log2_scale=False,
):
    """Generate a heatmap using the R package pheatmap.

    The input is a comma-separated table with a header; its first column is
    used as row names. Clustering remains disabled, matching the historical
    RSeQC behavior.
    """
    infile = Path(infile)
    outfile = Path(outfile)

    if not infile.is_file():
        raise FileNotFoundError(f"Input file does not exist: {infile}")

    filetype = str(filetype).lower().lstrip(".")
    if filetype == "jpg":
        filetype = "jpeg"
    if filetype not in SUPPORTED_FILETYPES:
        supported = ", ".join(sorted(SUPPORTED_FILETYPES))
        raise ValueError(
            f"Unsupported filetype {filetype!r}. Supported formats: {supported}"
        )

    if col_angle not in VALID_COLUMN_ANGLES:
        valid = ", ".join(str(x) for x in sorted(VALID_COLUMN_ANGLES))
        raise ValueError(f"col_angle must be one of: {valid}")
    if cell_width <= 0 or cell_height <= 0:
        raise ValueError("cell_width and cell_height must be greater than zero")
    if font_size <= 0:
        raise ValueError("font_size must be greater than zero")

    r_file = Path(str(outfile) + ".r")
    heatmap_file = Path(f"{outfile}.{filetype}")
    logging.info('Writing R code to "%s"', r_file)

    matrix_expr = "log2(as.matrix(dat) + 1)" if log2_scale else "as.matrix(dat)"
    display_numbers = "FALSE" if no_numbers else "TRUE"

    if no_numbers:
        logging.info("Numerical values will not be displayed on the heatmap")
    else:
        logging.info("Numerical values will be displayed on the heatmap")
        if log2_scale:
            logging.info("Numbers will be displayed on the log2 scale")

    args = [
        matrix_expr,
        f"filename={_r_quote(heatmap_file)}",
        f"cellwidth={cell_width}",
        f"cellheight={cell_height}",
        f"display_numbers={display_numbers}",
        f"angle_col={col_angle}",
        f"fontsize={font_size}",
    ]
    if not no_numbers:
        args.append(f"number_color={_r_quote(text_color)}")

    args.extend([
        "cluster_rows=FALSE",
        "cluster_cols=FALSE",
        "scale='none'",
        "color=colorRampPalette(c('#f1eef6','#d7b5d8','#df65b0','#ce1256'))(50)",
    ])

    r_code = [
        'if (!requireNamespace("pheatmap", quietly=TRUE)) {',
        '    stop("The R package pheatmap is required but is not installed.")',
        '}',
        'library(pheatmap)',
        (
            "dat <- read.table(file=" + _r_quote(infile)
            + ", sep=',', header=TRUE, row.names=1, check.names=FALSE)"
        ),
        "pheatmap(",
        "    " + ",\\n    ".join(args),
        ")",
        "",
    ]

    r_file.write_text("\\n".join(r_code), encoding="utf-8")

    logging.info('Running R script file "%s"', r_file)
    try:
        subprocess.run(["Rscript", str(r_file)], check=True)
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            "Rscript was not found. Install R or ensure Rscript is on PATH."
        ) from exc
    except subprocess.CalledProcessError:
        logging.exception('R script failed: "%s"', r_file)
        raise

    logging.info('Heatmap saved to "%s"', heatmap_file)
    return heatmap_file
