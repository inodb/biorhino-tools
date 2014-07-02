#!/usr/bin/env python
"""Remove fasta sequences from a fasta file by specifying a list of fasta ids to remove
"""
import argparse
from Bio import SeqIO


def get_fasta_ids(fastfiles, ids_file):
    get_ids = dict([(s.strip(), 0) for s in open(ids_file).readlines()])
    count_retrieved = 0
    count_removed = 0

    for ff in fastfiles:
        for record in SeqIO.parse(ff, "fasta"):
            try:
                get_ids[record.id] += 1
                print ">%s\n%s" % (record.id, record.seq)
                count_retrieved += 1
            except KeyError:
                count_removed += 1

    assert(count_retrieved == len(get_ids))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "contigs", nargs="+", help="Fasta files with contigs\n")
    parser.add_argument("-g", "--get_fasta_ids", default=0, help="File with one fasta id per line. Get those ids.\n")
    args = parser.parse_args()
    get_fasta_ids(args.contigs, args.get_fasta_ids)
