#!/usr/bin/env python3
"""Utilities for generating and analyzing DNA k-mers.

DNA sequences are expected to contain A, C, G, T, and optionally N.
"""

from __future__ import annotations

import itertools
import math
import re
from collections import Counter


_DNA_PATTERN = re.compile(r"^[ACGTN]+$")


def word_generator(seq, word_size, step_size, frame=0):
    """Yield fixed-length words from a sequence.

    Parameters
    ----------
    seq : str
        Input sequence.
    word_size : int
        Length of each word.
    step_size : int
        Distance between successive word starts.
    frame : int, default=0
        Starting offset in the sequence.
    """
    if word_size <= 0:
        raise ValueError("word_size must be greater than zero")
    if step_size <= 0:
        raise ValueError("step_size must be greater than zero")
    if frame < 0:
        raise ValueError("frame must be non-negative")

    for i in range(frame, len(seq), step_size):
        word = seq[i:i + word_size]
        if len(word) == word_size:
            yield word


def seq_generator(fastafile):
    """Yield ``[name, sequence]`` records from a FASTA-like file.

    Sequence lines containing characters other than A, C, G, T, or N are
    ignored, preserving the behavior of the original implementation.
    Header lines may begin with ``>`` or ``@``.
    """
    name = None
    seq_parts = []

    with open(fastafile, "r") as handle:
        for line in handle:
            line = line.strip().upper()

            if not line or line.startswith("#"):
                continue

            if line.startswith((">", "@")):
                if name is not None:
                    yield [name, "".join(seq_parts)]
                name = line.split()[0][1:]
                seq_parts = []
            elif _DNA_PATTERN.fullmatch(line):
                # Preserve the old behavior: valid sequence before a header
                # is accumulated and ultimately associated with an empty name.
                seq_parts.append(line)

    if name is not None:
        yield [name, "".join(seq_parts)]
    elif seq_parts:
        yield ["", "".join(seq_parts)]


def all_possible_kmer(l):
    """Yield all length-``l`` combinations of A, C, G, T, and N."""
    if l < 0:
        raise ValueError("k-mer length must be non-negative")

    for bases in itertools.product("ACGTN", repeat=l):
        yield "".join(bases)


def kmer_freq_file(fastafile, word_size, step_size=1, frame=0, min_count=0):
    """Return k-mer counts from sequences in a FASTA-like file.

    K-mers containing ``N`` are omitted from the returned dictionary.
    """
    count_table = Counter()

    for _name, seq in seq_generator(fastafile):
        count_table.update(
            word_generator(
                seq,
                word_size=word_size,
                step_size=step_size,
                frame=frame,
            )
        )

    ret_dict = {}
    for kmer in all_possible_kmer(word_size):
        count = count_table.get(kmer, 0)
        if count >= min_count and "N" not in kmer:
            ret_dict[kmer] = count

    return ret_dict


def kmer_freq_seq(seq, word_size, step_size=1, frame=0, min_count=0):
    """Print k-mer counts for a DNA sequence."""
    count_table = Counter(
        word_generator(
            seq,
            word_size=word_size,
            step_size=step_size,
            frame=frame,
        )
    )

    for kmer in all_possible_kmer(word_size):
        count = count_table.get(kmer, 0)
        if count >= min_count:
            print(f"{kmer}\t{count}")


def kmer_ratio(seq, word_size, step_size, coding, noncoding):
    """Calculate the frame-0 coding/noncoding k-mer log-ratio score.

    This preserves the active algorithm in the original implementation:
    only frame 0 is evaluated. K-mers absent from either reference mapping
    are skipped.

    Returns
    -------
    float
        Mean frame-0 score. Returns 0.0 if the sequence is shorter than
        ``word_size`` or if no usable k-mers are available.
    """
    if len(seq) < word_size:
        return 0.0

    sum_of_log_ratio_0 = 0.0
    frame0_count = 0

    for kmer in word_generator(
        seq=seq,
        word_size=word_size,
        step_size=step_size,
        frame=0,
    ):
        if kmer not in coding or kmer not in noncoding:
            continue

        coding_count = coding[kmer]
        noncoding_count = noncoding[kmer]

        if coding_count > 0 and noncoding_count > 0:
            sum_of_log_ratio_0 += math.log(coding_count / noncoding_count)
        elif coding_count > 0 and noncoding_count == 0:
            sum_of_log_ratio_0 += 1
        elif coding_count == 0 and noncoding_count == 0:
            continue
        elif coding_count == 0 and noncoding_count > 0:
            sum_of_log_ratio_0 -= 1
        else:
            continue

        frame0_count += 1

    if frame0_count == 0:
        return 0.0

    return sum_of_log_ratio_0 / frame0_count
