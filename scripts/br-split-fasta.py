#!/usr/bin/env python
"""
Splits each record of a fasta file into a separate fasta file where the
filename is the name of the accession id.
"""
import os
import errno
import argparse
from Bio import SeqIO


def make_dir(directory):
    """Make directory unless existing. Ignore error in the latter case."""
    try:
        os.makedirs(directory)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def split_fasta(fastfiles, output_folder):
    make_dir(output_folder)

    for ff in fastfiles:
        for record in SeqIO.parse(ff, "fasta"):
            with open(output_folder + "/" + record.id.replace(" ", "_").replace("/", "_") + ".fa", "w") as fh:
                fh.write(">%s\n%s" % (record.id, record.seq))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "contigs", nargs="+", help="Fasta files with contigs\n")
    parser.add_argument(
        "output_folder", help="Output folder with splitted records\n")
    args = parser.parse_args()
    split_fasta(args.contigs, args.output_folder.rstrip('/'))
