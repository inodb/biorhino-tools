#!/usr/bin/env python
"""Cut up fasta file in non-overlapping or overlapping parts of equal length.
"""
import argparse
from Bio import SeqIO


def filter_fasta(fastfiles, min_length):
    for ff in fastfiles:
        for record in SeqIO.parse(ff, "fasta"):
            if len(record.seq) >= min_length:
                print ">%s\n%s" % (record.id, record.seq)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "contigs", nargs="+", help="Fasta files with contigs\n")
    parser.add_argument("-m", "--min_length", default=0, type=int, help="Minimum length\n")
    args = parser.parse_args()
    filter_fasta(args.contigs, args.min_length)
