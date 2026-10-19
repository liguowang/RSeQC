#!/usr/bin/env python3
"""Utilities for manipulating FASTA-format sequence files.

This module was originally written for Python 2.7 and later converted to
Python 3.
"""

from __future__ import annotations

import collections
import sys
from importlib.metadata import PackageNotFoundError, version
from optparse import OptionParser

import numpy

from qcmodule import FrameKmer

try:
    __version__ = version("RSeQC")
except PackageNotFoundError:
    __version__ = "unknown"

__author__ = "Liguo Wang"


class Fasta:
    """Manipulate FASTA-format sequence files."""

    def __init__(self, fastafile=None):
        """Initialize object; lowercase sequence is converted to uppercase."""
        self.seqs = {}
        self.IDs = []
        self.transtab = str.maketrans("ACGTNX", "TGCANX")
        self.filename = fastafile

        if fastafile is None:
            return

        name = None
        seq_parts = []

        with open(fastafile, "r") as fasta_handle:
            for line in fasta_handle:
                line = line.strip()
                if not line:
                    continue

                if line.startswith(">"):
                    if name is not None:
                        self.seqs[name] = "".join(seq_parts)
                    name = line[1:]
                    seq_parts = []
                    self.IDs.append(name)
                    print("\tloading " + name + " ...", file=sys.stderr)
                else:
                    if name is None:
                        raise ValueError(
                            "Invalid FASTA file: sequence data found before first header"
                        )
                    seq_parts.append(line.upper())

        if name is not None:
            self.seqs[name] = "".join(seq_parts)

    def addSeq(self, id, seq):
        """Add a sequence to the current data."""
        if id in self.seqs:
            print(id + " already exists!", file=sys.stderr)
            return
        self.seqs[id] = seq.upper()
        self.IDs.append(id)

    def getNames(self, file=None):
        """Return all sequence IDs."""
        return self.IDs

    def getSeq(self, seqID=None):
        """Return the specified sequence, or all sequences if seqID is None."""
        if seqID is None:
            return list(self.seqs.values())
        return self.seqs[seqID]

    def printSeqs(self, n=50):
        """Print all sequences in FASTA format."""
        for name, seq in self.seqs.items():
            print(">" + name)
            for i in range(0, len(seq), n):
                print(seq[i:i + n])

    def getSeqLen(self, seqID=None):
        """Return sequence lengths."""
        seqlen = collections.defaultdict(dict)
        if seqID is None:
            for name, seq in self.seqs.items():
                seqlen[name] = len(seq)
        else:
            try:
                seqlen[seqID] = len(self.seqs[seqID])
            except KeyError:
                print("Not found", file=sys.stderr)
        return seqlen

    def countBase(self, pattern=None):
        """Count bases, or non-overlapping occurrences of a supplied pattern."""
        if pattern is None:
            print("ID\tTotal\tA\tC\tG\tT\tN\tX")
            for name, seq in self.seqs.items():
                print(
                    f"{name}\t{len(seq)}\t{seq.count('A')}\t{seq.count('C')}\t"
                    f"{seq.count('G')}\t{seq.count('T')}\t{seq.count('N')}\t"
                    f"{seq.count('X')}"
                )
        else:
            for name, seq in self.seqs.items():
                print(f"{name}\t{len(seq)}\t{seq.count(pattern)}")

    def cal_entropy(self, l=3):
        """Calculate entropy for each sequence."""
        for seq_id, seq in self.seqs.items():
            entropy = 0
            dna_chars_uniq = FrameKmer.all_possible_kmer(l)
            dna_len = len(seq)

            if dna_len == 0:
                yield seq_id, entropy
                continue

            for kmer in dna_chars_uniq:
                if "N" in kmer:
                    continue
                prop = seq.count(kmer) / float(dna_len)
                if prop == 0:
                    continue
                information = numpy.log2(1.0 / prop)
                entropy += prop * information
            yield seq_id, entropy

    def revComp(self, seqID=None):
        """Return reverse-complemented sequence(s).

        If seqID is specified, return that sequence's reverse complement.
        Otherwise print all reverse-complemented sequences in FASTA format.
        """
        if seqID is None:
            for name, seq in self.seqs.items():
                print(">" + name + "_rev")
                print(seq.translate(self.transtab)[::-1])
            return None
        return self.seqs[seqID].translate(self.transtab)[::-1]

    def getUniqSeqs(self):
        """Print non-redundant sequences, reporting each sequence once."""
        seq2Name = {}
        seq2Count = {}
        for key, value in self.seqs.items():
            seq2Name[value] = key
            if value in seq2Count:
                seq2Count[value] += 1
            else:
                seq2Count[value] = 1

        for value in seq2Name:
            print(">" + str(seq2Name[value]) + "_" + str(seq2Count[value]))
            print(value)

    def findPattern(self, pat, outfile, seqID=None, rev=True):
        """Find pattern occurrences and write coordinates in BED format."""
        length = len(pat)
        pat = pat.upper()

        if seqID is None:
            sequences = self.seqs.items()
        else:
            sequences = [(seqID, self.seqs[seqID])]

        with open(outfile, "w") as fout:
            for name, seq in sequences:
                start = 0
                while True:
                    hit = seq.find(pat, start)
                    if hit == -1:
                        break
                    print(
                        f"{name}\t{hit}\t{hit + length}\t{pat}\t0\t+",
                        file=fout,
                    )
                    start = hit + 1

            if rev:
                pat_rev = pat.translate(self.transtab)[::-1]
                if seqID is None:
                    sequences = self.seqs.items()
                else:
                    sequences = [(seqID, self.seqs[seqID])]

                for name, seq in sequences:
                    start = 0
                    while True:
                        hit = seq.find(pat_rev, start)
                        if hit == -1:
                            break
                        print(
                            f"{name}\t{hit}\t{hit + length}\t{pat}\t0\t-",
                            file=fout,
                        )
                        start = hit + 1

    def fetchSeq(self, chr=None, st=None, end=None, infile=None, outfile=None):
        """Fetch sequence by coordinates.

        Coordinates are 0-based, half-open. ``infile`` may contain BED3,
        BED6, or BED12 records.
        """
        if infile is not None and outfile is not None:
            with open(outfile, "w") as fout, open(infile, "r") as fin:
                for line in fin:
                    fields = line.strip().split()
                    if not fields:
                        continue

                    if len(fields) == 3:
                        print(
                            fields[0] + ":" + fields[1] + "-" + fields[2]
                            + "\tstrand=+",
                            file=fout,
                        )
                        print(
                            self.seqs[fields[0]][
                                int(fields[1]):int(fields[2])
                            ].upper(),
                            file=fout,
                        )
                    elif len(fields) >= 6:
                        if fields[5] == "-":
                            print(
                                fields[0] + ":" + fields[1] + "-" + fields[2]
                                + "\tstrand=-",
                                file=fout,
                            )
                            print(
                                self.seqs[fields[0]][
                                    int(fields[1]):int(fields[2])
                                ].translate(self.transtab)[::-1].upper(),
                                file=fout,
                            )
                        else:
                            print(
                                fields[0] + ":" + fields[1] + "-" + fields[2]
                                + "\tstrand=+",
                                file=fout,
                            )
                            print(
                                self.seqs[fields[0]][
                                    int(fields[1]):int(fields[2])
                                ].upper(),
                                file=fout,
                            )
                    else:
                        raise ValueError(
                            "Input coordinate file must contain BED3, BED6, "
                            "or BED12 records"
                        )
            return None

        try:
            return self.seqs[chr][st:end].upper()
        except (KeyError, TypeError):
            print(
                "cannot fetch sequence from "
                + str(self.filename)
                + " for "
                + str(chr)
                + ":"
                + str(st)
                + "-"
                + str(end),
                file=sys.stderr,
            )
            return ""


def main():
    parser = OptionParser()
    parser.add_option("-i", "--input_file", dest="in_file", help="input file name")
    options, _args = parser.parse_args()
    obj = Fasta(options.in_file)
    obj.printSeqs(n=80)


if __name__ == "__main__":
    main()
