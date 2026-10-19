#!/usr/bin/env python3
"""Utilities for finding open reading frames (ORFs) in transcript sequences."""

from __future__ import annotations

__author__ = "Liguo Wang"

DEFAULT_START_CODONS = ("ATG",)
DEFAULT_STOP_CODONS = ("TAG", "TAA", "TGA")

start_codon = list(DEFAULT_START_CODONS)
stop_codon = list(DEFAULT_STOP_CODONS)

_COMPLEMENT = str.maketrans(
    {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "N": "N",
        "X": "X",
    }
)


def _reverse_comp(seq):
    """Return the reverse complement of a DNA sequence."""
    try:
        return seq.upper().translate(_COMPLEMENT)[::-1]
    except KeyError as exc:
        raise ValueError(f"Unsupported nucleotide in sequence: {exc.args[0]!r}") from exc


def _parse_codons(value, defaults, label):
    """Parse and validate a comma-separated codon specification."""
    if value is None:
        return tuple(defaults)

    codons = tuple(
        codon.strip().upper()
        for codon in value.split(",")
        if codon.strip()
    )
    if not codons:
        raise ValueError(f"No {label} codons were specified")

    for codon in codons:
        if len(codon) != 3:
            raise ValueError(f"Invalid {label} codon {codon!r}: codons must have length 3")
        if any(base not in "ACGTNX" for base in codon):
            raise ValueError(
                f"Invalid {label} codon {codon!r}: unsupported nucleotide"
            )
    return codons


def _find_longest_orf_range(seq, start_codons, stop_codons):
    """Return ``(start, stop_start)`` for the longest complete ORF.

    The stop coordinate is the position at which the stop codon begins.
    This preserves the coordinate convention of the historical implementation.
    If multiple ORFs have the same maximum length, the last one in genomic
    sequence order is returned, matching the old dictionary iteration behavior.
    """
    start_positions = []
    stop_positions = []

    for codon in start_codons:
        start = seq.find(codon)
        while start != -1:
            start_positions.append(start)
            start = seq.find(codon, start + 1)

    for codon in stop_codons:
        stop = seq.find(codon)
        while stop != -1:
            stop_positions.append(stop)
            stop = seq.find(codon, stop + 1)

    # For each start codon, the historical implementation chose the first
    # in-frame downstream stop codon.
    candidates = []
    for start in sorted(set(start_positions)):
        in_frame_stops = [
            stop
            for stop in stop_positions
            if stop > start and (stop - start) % 3 == 0
        ]
        if in_frame_stops:
            stop = min(in_frame_stops)
            candidates.append((start, stop))

    if not candidates:
        return None

    max_length = max(stop - start for start, stop in candidates)
    longest = [
        (start, stop)
        for start, stop in candidates
        if stop - start == max_length
    ]
    return longest[-1]


def longest_orf(seq, strandness, sc=None, tc=None):
    """Return the longest complete ORF on the requested strand.

    Parameters
    ----------
    seq : str
        Transcript sequence.
    strandness : {"+", "-"}
        Strand to search.
    sc : str or None
        Optional comma-separated start codons.
    tc : str or None
        Optional comma-separated stop codons.

    Returns
    -------
    str
        ORF sequence excluding the terminal stop codon, preserving historical
        RSeQC behavior. Returns an empty string if no complete ORF is found.
    """
    if strandness not in {"+", "-"}:
        raise ValueError("strandness must be '+' or '-'")

    start_codons = _parse_codons(sc, DEFAULT_START_CODONS, "start")
    stop_codons = _parse_codons(tc, DEFAULT_STOP_CODONS, "stop")

    dna_seq = seq.upper()
    if strandness == "-":
        dna_seq = _reverse_comp(dna_seq)

    orf_range = _find_longest_orf_range(dna_seq, start_codons, stop_codons)
    if orf_range is None:
        return ""

    start, stop = orf_range
    return dna_seq[start:stop]


def _transcript_offset_to_genome_plus(offset, exon_sizes, exon_starts):
    """Map a 0-based transcript offset to a genomic coordinate on '+'."""
    remaining = offset
    for exon_size, exon_start in zip(exon_sizes, exon_starts):
        if remaining > exon_size:
            remaining -= exon_size
            continue
        return exon_start + remaining
    raise ValueError("ORF coordinate exceeds transcript exon length")


def _transcript_offset_to_genome_minus(offset, exon_sizes, exon_ends):
    """Map a 0-based transcript offset to a genomic coordinate on '-'."""
    remaining = offset
    for exon_size, exon_end in zip(reversed(exon_sizes), reversed(exon_ends)):
        if remaining > exon_size:
            remaining -= exon_size
            continue
        return exon_end - remaining
    raise ValueError("ORF coordinate exceeds transcript exon length")


def longest_orf_bed(seq, bedline, sc=None, tc=None):
    """Set BED12 thickStart/thickEnd to the longest complete ORF.

    The ORF search and coordinate conventions intentionally match
    ``longest_orf`` and the historical RSeQC implementation.

    Returns
    -------
    str or None
        Updated BED12 line, or ``None`` when no complete ORF is found.
    """
    start_codons = _parse_codons(sc, DEFAULT_START_CODONS, "start")
    stop_codons = _parse_codons(tc, DEFAULT_STOP_CODONS, "stop")

    fields = bedline.split()
    if len(fields) < 12:
        raise ValueError("bedline must contain at least 12 BED fields")

    tx_start = int(fields[1])
    strand = fields[5]
    if strand not in {"+", "-"}:
        raise ValueError("BED strand must be '+' or '-'")

    exon_sizes = [int(value) for value in fields[10].rstrip(",").split(",") if value]
    exon_offsets = [int(value) for value in fields[11].rstrip(",").split(",") if value]

    if len(exon_sizes) != len(exon_offsets):
        raise ValueError("BED exon size and exon start counts do not match")
    if not exon_sizes:
        raise ValueError("BED record contains no exons")

    exon_starts = [tx_start + offset for offset in exon_offsets]
    exon_ends = [
        exon_start + exon_size
        for exon_start, exon_size in zip(exon_starts, exon_sizes)
    ]

    dna_seq = seq.upper()
    if strand == "-":
        dna_seq = _reverse_comp(dna_seq)

    orf_range = _find_longest_orf_range(dna_seq, start_codons, stop_codons)
    if orf_range is None:
        return None

    cds_start, cds_stop = orf_range

    if strand == "+":
        cds_pos1 = _transcript_offset_to_genome_plus(
            cds_start, exon_sizes, exon_starts
        )
        cds_pos2 = _transcript_offset_to_genome_plus(
            cds_stop, exon_sizes, exon_starts
        ) + 3
        fields[6] = str(cds_pos1)
        fields[7] = str(cds_pos2)
    else:
        cds_pos1 = _transcript_offset_to_genome_minus(
            cds_start, exon_sizes, exon_ends
        )
        cds_pos2 = _transcript_offset_to_genome_minus(
            cds_stop, exon_sizes, exon_ends
        ) - 3
        fields[7] = str(cds_pos1)
        fields[6] = str(cds_pos2)

    return "\t".join(fields)
