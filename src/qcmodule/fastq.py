#!/usr/bin/env python3
"""Utilities for reading FASTA/FASTQ files and generating sequence logos."""

from __future__ import annotations

import collections
import logging
import sys

import logomaker
import matplotlib.pyplot as plt
import pandas as pd

from qcmodule import ireader


def fasta_iter(infile):
    """Yield sequence lines from a FASTA file.

    Parameters
    ----------
    infile : str
        Input FASTA file.

    Yields
    ------
    str
        Sequence lines.
    """
    logging.info('Reading FASTA file "%s" ...', infile)
    for line in ireader.reader(infile):
        line = line.strip()
        if not line or line.startswith(">"):
            continue
        yield line


def fastq_iter(infile, mode="seq"):
    """Yield sequence or quality strings from a four-line FASTQ file.

    Parameters
    ----------
    infile : str
        Input FASTQ file.
    mode : {"seq", "qual"}
        Select whether sequence or quality strings are yielded.

    Yields
    ------
    str
        Sequence or quality string.

    Raises
    ------
    ValueError
        If ``mode`` is invalid or the FASTQ file contains an incomplete
        four-line record.
    """
    if mode not in {"seq", "qual"}:
        raise ValueError("mode must be 'seq' or 'qual'")

    logging.info('Reading FASTQ file "%s" ...', infile)

    record = []
    for line in ireader.reader(infile):
        # Do not discard blank lines here: FASTQ is record-oriented and
        # silently skipping a line would shift all subsequent records.
        record.append(line.rstrip("\r\n"))

        if len(record) == 4:
            header, seq, plus, qual = record

            if not header.startswith("@"):
                raise ValueError(
                    f"Invalid FASTQ record: header does not start with '@': {header!r}"
                )
            if not plus.startswith("+"):
                raise ValueError(
                    f"Invalid FASTQ record for {header!r}: third line does not start with '+'"
                )
            if len(seq) != len(qual):
                raise ValueError(
                    f"Invalid FASTQ record for {header!r}: sequence and quality "
                    f"lengths differ ({len(seq)} != {len(qual)})"
                )

            yield seq if mode == "seq" else qual
            record = []

    if record:
        raise ValueError(
            f"Incomplete FASTQ record at end of {infile!r}: "
            f"found {len(record)} of 4 lines"
        )


def qual2countMat(q_obj, limit, step_size=100000):
    """Generate a positional Phred-quality count data frame.

    Parameters
    ----------
    q_obj : iterable
        Generator returned by ``fastq_iter``.
    limit : int or None
        Read at most this many quality strings. ``None`` means no limit.
    step_size : int
        Report progress every ``step_size`` records.

    Returns
    -------
    pandas.DataFrame
        Rows are read positions and columns are Phred scores.
    """
    if limit is not None and limit < 0:
        raise ValueError("limit must be non-negative or None")
    if step_size <= 0:
        raise ValueError("step_size must be greater than zero")

    dat = collections.defaultdict(dict)
    count = 0

    if limit == 0:
        return pd.DataFrame()

    for qstr in q_obj:
        count += 1
        for pos, char in enumerate(qstr):
            q_score = ord(char) - 33
            dat[pos][q_score] = dat[pos].get(q_score, 0) + 1

        if count % step_size == 0:
            print(
                f"{count} quality sequences finished\r",
                end=" ",
                file=sys.stderr,
            )
        if limit is not None and count >= limit:
            break

    logging.info("%d quality sequences finished", count)
    logging.info("Make data frame from dict of dict ...")
    qual_mat = pd.DataFrame.from_dict(dat)
    logging.info("Filling NA as zero ...")
    qual_mat = qual_mat.fillna(0)
    qual_mat.index.name = "pos"
    return qual_mat.T


def seq2countMat(s_obj, limit, step_size=100000, exclude_N=False):
    """Generate a positional nucleotide-count data frame.

    Parameters
    ----------
    s_obj : iterable
        Generator returned by ``fastq_iter`` or ``fasta_iter``.
    limit : int or None
        Examine at most this many input sequences. ``None`` means no limit.
    step_size : int
        Report progress every ``step_size`` input sequences.
    exclude_N : bool
        If True, sequences containing uppercase or lowercase ``N`` are skipped
        from the count matrix.

    Returns
    -------
    pandas.DataFrame
        Rows are read positions and columns are observed bases.
    """
    if limit is not None and limit < 0:
        raise ValueError("limit must be non-negative or None")
    if step_size <= 0:
        raise ValueError("step_size must be greater than zero")

    mat = collections.defaultdict(dict)
    count = 0

    if limit == 0:
        return pd.DataFrame()

    for seq in s_obj:
        count += 1

        if not (exclude_N and "N" in seq.upper()):
            for pos, base in enumerate(seq):
                mat[pos][base] = mat[pos].get(base, 0) + 1

        if count % step_size == 0:
            print(f"{count} sequences finished\r", end=" ", file=sys.stderr)

        # Keep the original interpretation of limit: number of input
        # sequences examined, including sequences skipped by exclude_N.
        if limit is not None and count >= limit:
            break

    logging.info("%d sequences finished", count)
    logging.info("Make data frame from dict of dict ...")
    count_mat = pd.DataFrame.from_dict(mat)
    logging.info("Filling NA as zero ...")
    count_mat = count_mat.fillna(0)
    count_mat.index.name = "pos"
    return count_mat.T


def make_logo(
    mat,
    outfile,
    exclude_N=False,
    font_name="sans",
    stack_order="big_on_top",
    flip_below=True,
    shade_below=0.0,
    fade_below=0.0,
    highlight_start=0,
    highlight_end=15,
    oformat="pdf",
):
    """Create mean-centered and standard nucleotide sequence logos.

    Parameters
    ----------
    mat : pandas.DataFrame
        Nucleotide frequencies table.
    outfile : str
        Output filename prefix.
    exclude_N : bool
        If True, use colors for A/C/G/T only.
    font_name : str
        Character font used by logomaker.
    stack_order : {"big_on_top", "small_on_top", "fixed"}
        Character stacking order.
    flip_below : bool
        Flip characters below the x-axis.
    shade_below : float
        Shading amount for characters below the x-axis, in [0, 1].
    fade_below : float
        Fading amount for characters below the x-axis, in [0, 1].
    highlight_start, highlight_end : int or None
        Inclusive positions to highlight. Set either to None to disable.
    oformat : str
        Matplotlib output format, for example ``pdf`` or ``png``.
    """
    if not 0.0 <= shade_below <= 1.0:
        raise ValueError("shade_below must be between 0.0 and 1.0")
    if not 0.0 <= fade_below <= 1.0:
        raise ValueError("fade_below must be between 0.0 and 1.0")
    if stack_order not in {"big_on_top", "small_on_top", "fixed"}:
        raise ValueError(
            "stack_order must be 'big_on_top', 'small_on_top', or 'fixed'"
        )

    highlight = highlight_start is not None and highlight_end is not None
    if highlight:
        if not isinstance(highlight_start, int) or not isinstance(highlight_end, int):
            raise TypeError("highlight_start and highlight_end must be integers or None")
        if highlight_start < 0 or highlight_end < 0 or highlight_start > highlight_end:
            raise ValueError("Incorrect highlight positions")
        if len(mat) and highlight_end >= len(mat):
            raise ValueError(
                f"highlight_end ({highlight_end}) is outside the logo "
                f"range 0..{len(mat) - 1}"
            )

    logging.info("Making logo ...")
    logging.debug("Font name is: %s", font_name)
    logging.debug("Stack order is: %s", stack_order)
    logging.debug("Flip below flag: %s", flip_below)
    logging.debug("Shade below score: %f", shade_below)
    logging.debug("Fade below score: %f", fade_below)

    if exclude_N:
        logging.info("'N' will be excluded.")
        color = {"A": "green", "C": "blue", "G": "orange", "T": "red"}
    else:
        logging.info("'N' will be kept.")
        color = {
            "A": "green",
            "C": "blue",
            "G": "orange",
            "T": "red",
            "N": "grey",
        }

    fmt = oformat.lower()

    mean_centered_file = f"{outfile}.logo.mean_centered.{fmt}"
    logging.info('Mean-centered logo saved to "%s".', mean_centered_file)
    fig, ax = plt.subplots()
    logo = logomaker.Logo(
        mat,
        ax=ax,
        center_values=True,
        color_scheme=color,
        font_name=font_name,
        stack_order=stack_order,
        flip_below=flip_below,
        shade_below=shade_below,
        fade_below=fade_below,
    )
    if highlight:
        logging.info(
            "Highlight logo from %d to %d", highlight_start, highlight_end
        )
        logo.highlight_position_range(
            pmin=highlight_start, pmax=highlight_end
        )
    fig.savefig(mean_centered_file)
    plt.close(fig)

    logo_file = f"{outfile}.logo.{fmt}"
    logging.info('Logo saved to "%s".', logo_file)
    fig, ax = plt.subplots()
    logo = logomaker.Logo(
        mat,
        ax=ax,
        center_values=False,
        color_scheme=color,
        font_name=font_name,
        stack_order=stack_order,
        flip_below=flip_below,
        shade_below=shade_below,
        fade_below=fade_below,
    )
    if highlight:
        logging.info(
            "Highlight logo from %d to %d", highlight_start, highlight_end
        )
        logo.highlight_position_range(
            pmin=highlight_start, pmax=highlight_end
        )
    fig.savefig(logo_file)
    plt.close(fig)
