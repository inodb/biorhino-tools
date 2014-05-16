#!/usr/bin/env python
"""Remove fasta ids from a fasta file
"""
import argparse
from Bio import SeqIO


def remove_fasta_ids(fastfiles, ids_file):
    rem_ids = dict([(s.strip(), 0) for s in open(ids_file).readlines()])
    count_removed = 0

    for ff in fastfiles:
        for record in SeqIO.parse(ff, "fasta"):
            try:
                rem_ids[record.id] += 1
                count_removed += 1
            except KeyError:
                print ">%s\n%s" % (record.id, record.seq)

    assert(count_removed == len(rem_ids))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "contigs", nargs="+", help="Fasta files with contigs\n")
    parser.add_argument("-r", "--remove_fasta_ids", default=0, help="File with one fasta id per line. Remove those ids.\n")
    args = parser.parse_args()
    remove_fasta_ids(args.contigs, args.remove_fasta_ids)
