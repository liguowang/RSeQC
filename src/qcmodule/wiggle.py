#!/usr/bin/env python3
"""Utilities for reading and querying UCSC WIG (wiggle) files."""

from __future__ import annotations

import math
import sys

import bx.wiggle
from bx.binned_array import BinnedArray


__author__ = "Liguo Wang"


def _validate_range(start, end):
    """Validate a 0-based, half-open interval."""
    if start < 0:
        raise ValueError("start must be non-negative")
    if end < start:
        raise ValueError("end must be greater than or equal to start")


def _valid_scores(scores):
    """Return finite/non-NaN numeric scores from an iterable."""
    return [float(value) for value in scores if not math.isnan(float(value))]


class ParseWig:
    """Read and query a WIG file.

    If a coordinate occurs more than once, the last value encountered is
    retained, matching the historical behavior of this class.
    """

    def __init__(self, wigFile):
        """Read a WIG file into chromosome-specific ``BinnedArray`` objects."""
        self.scores = {}
        point_count = 0

        with open(wigFile, "r") as handle:
            for chrom, pos, value in bx.wiggle.Reader(handle):
                chrom = chrom.upper()
                if chrom not in self.scores:
                    self.scores[chrom] = BinnedArray()
                self.scores[chrom][pos] = value
                point_count += 1

                if point_count % 100000 == 0:
                    print(
                        f"{point_count} datapoints loaded \r",
                        end="",
                        file=sys.stderr,
                    )

        print(f"total {point_count} points loaded")

    def _scores_by_range(self, chrom, start, end):
        """Return raw scores for a genomic interval."""
        _validate_range(start, end)
        chrom = chrom.upper()
        if chrom not in self.scores:
            raise KeyError(f"Chromosome {chrom!r} is not present in WIG file")
        return [self.scores[chrom][pos] for pos in range(start, end)]

    def fetch_all_scores(self, chr, st, end):
        """Return all scores in the 0-based, half-open interval ``[st, end)``.

        Positions without a WIG score are represented by NaN.
        """
        return self._scores_by_range(chr, st, end)

    def fetch_max_scores(self, chr, st, end):
        """Return the maximum raw score in ``[st, end)``.

        This preserves the historical behavior, including NaN values supplied
        by ``BinnedArray``.
        """
        return max(self._scores_by_range(chr, st, end))

    def fetch_min_scores(self, chr, st, end):
        """Return the minimum raw score in ``[st, end)``."""
        return min(self._scores_by_range(chr, st, end))

    def fetch_avg_scores(self, chr, st, end):
        """Return the historical interval average for ``[st, end)``.

        Missing positions contribute zero to the numerator because NaN scores
        are omitted, but the denominator remains the full interval length.
        This intentionally preserves the original RSeQC behavior.
        """
        scores = self._scores_by_range(chr, st, end)
        if not scores:
            raise ValueError("cannot calculate an average for an empty interval")
        valid = _valid_scores(scores)
        return sum(valid) / len(scores)

    def fetch_sum_scores(self, chr, st, end):
        """Return the sum of defined scores in ``[st, end)``."""
        return sum(_valid_scores(self._scores_by_range(chr, st, end)))


class ParseWig2:
    """Read and query a WIG file while summing duplicate coordinates.

    This class is slower than :class:`ParseWig`, but if the same genomic
    coordinate occurs multiple times, its values are added together.
    """

    def __init__(self, wigFile):
        """Read a WIG file and sum values at duplicate coordinates."""
        self.scores = {}
        point_count = 0

        with open(wigFile, "r") as handle:
            for chrom, pos, value in bx.wiggle.Reader(handle):
                chrom = chrom.upper()
                if chrom not in self.scores:
                    self.scores[chrom] = BinnedArray()

                current = self.scores[chrom][pos]
                if math.isnan(float(current)):
                    self.scores[chrom][pos] = value
                else:
                    self.scores[chrom][pos] += value

                point_count += 1
                if point_count % 100000 == 0:
                    print(
                        f"{point_count} datapoints loaded \r",
                        end="",
                        file=sys.stderr,
                    )

        print(f"total {point_count} points loaded")

    def _scores_by_range(self, chrom, start, end):
        _validate_range(start, end)
        chrom = chrom.upper()
        if chrom not in self.scores:
            raise KeyError(f"Chromosome {chrom!r} is not present in WIG file")
        return [self.scores[chrom][pos] for pos in range(start, end)]

    def _scores_by_positions(self, chrom, positions):
        chrom = chrom.upper()
        if chrom not in self.scores:
            raise KeyError(f"Chromosome {chrom!r} is not present in WIG file")
        return [self.scores[chrom][pos] for pos in positions]

    def fetch_all_scores_by_range(self, chr, st, end):
        """Return all scores in the 0-based, half-open interval ``[st, end)``."""
        return self._scores_by_range(chr, st, end)

    def fetch_all_scores_by_positions(self, chr, lst):
        """Return scores at the specified positions."""
        return self._scores_by_positions(chr, lst)

    def fetch_max_scores_by_range(self, chr, st, end):
        """Return the maximum raw score in ``[st, end)``."""
        return max(self._scores_by_range(chr, st, end))

    def fetch_max_scores_by_positions(self, chr, lst):
        """Return the maximum raw score at the specified positions."""
        return max(self._scores_by_positions(chr, lst))

    def fetch_min_scores_by_range(self, chr, st, end):
        """Return the minimum raw score in ``[st, end)``."""
        return min(self._scores_by_range(chr, st, end))

    def fetch_min_scores_by_positions(self, chr, lst):
        """Return the minimum raw score at the specified positions."""
        return min(self._scores_by_positions(chr, lst))

    def fetch_avg_scores_by_range(self, chr, st, end):
        """Return the historical interval average for ``[st, end)``.

        As in the original implementation, missing positions are omitted from
        the numerator but retained in the full interval-length denominator.
        """
        scores = self._scores_by_range(chr, st, end)
        if not scores:
            raise ValueError("cannot calculate an average for an empty interval")
        return sum(_valid_scores(scores)) / len(scores)

    def fetch_avg_scores_by_positions(self, chr, lst):
        """Return the average of defined scores at specified positions."""
        scores = _valid_scores(self._scores_by_positions(chr, lst))
        if not scores:
            raise ValueError("no defined scores are available at these positions")
        return sum(scores) / len(scores)

    def fetch_sum_scores_by_range(self, chr, st, end):
        """Return the sum of defined scores in ``[st, end)``."""
        return sum(_valid_scores(self._scores_by_range(chr, st, end)))

    def fetch_sum_scores_by_positions(self, chr, lst):
        """Return the sum of defined scores at specified positions."""
        return sum(_valid_scores(self._scores_by_positions(chr, lst)))

    def distriub_wig(self, bed, till_count=100):
        """Parse BED12 records for the historical coverage-distribution method.

        The original implementation stops after parsing exon coordinates and
        never performs the documented distribution calculation. That behavior
        is intentionally not invented here; this method currently validates
        and parses the BED12 input, as the original code did.
        """
        if till_count <= 0:
            raise ValueError("till_count must be greater than zero")

        print(f"Reading {bed} ...", file=sys.stderr)
        with open(bed, "r") as handle:
            for line in handle:
                if line.startswith(("#", "track", "browser")):
                    continue

                try:
                    fields = line.rstrip("\r\n").split()
                    if len(fields) < 12:
                        raise ValueError("BED record has fewer than 12 fields")

                    tx_start = int(fields[1])
                    exon_starts = [
                        int(value)
                        for value in fields[11].rstrip(",").split(",")
                    ]
                    exon_starts = [
                        start + tx_start for start in exon_starts
                    ]
                    exon_sizes = [
                        int(value)
                        for value in fields[10].rstrip(",").split(",")
                    ]
                    if len(exon_starts) != len(exon_sizes):
                        raise ValueError(
                            "BED exon start and exon size counts differ"
                        )

                    # Preserve parsing performed by the historical method.
                    _exon_ends = [
                        start + size
                        for start, size in zip(exon_starts, exon_sizes)
                    ]
                except (ValueError, IndexError) as exc:
                    print(
                        "[NOTE: input BED must be valid BED12] "
                        f"skipped this line: {line.rstrip()} ({exc})",
                        file=sys.stderr,
                    )
                    continue
