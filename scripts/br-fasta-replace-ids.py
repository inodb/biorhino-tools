#!/usr/bin/env python
"""Replace ids in fasta file with given file with id-newid pairs.
"""
import argparse
from Bio import SeqIO


def replace_ids(fastfiles, replace_ids_file, sep):
    replace_dict = {}
    with open(replace_ids_file) as rpfh:
        for line in rpfh:
            splits = line.split(sep)
            assert(len(splits) == 2)
            replace_dict[splits[0]] = splits[1].rstrip('\n')

    for ff in fastfiles:
        for record in SeqIO.parse(ff, "fasta"):
            print ">%s\n%s" % (replace_dict[record.id], record.seq)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "replace_ids_file", help="Tab-separated file with key-value pairs that is to be replaced\n")
    parser.add_argument(
        "contigs", nargs="+", help="Fasta files with contigs\n")
    parser.add_argument("-d", "--delimiter", default="\t", help="Delimiter in replace_ids_file\n")
    args = parser.parse_args()
    replace_ids(args.contigs, args.replace_ids_file, args.delimiter)
