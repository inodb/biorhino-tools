#!/usr/bin/env python
"""Get pairs that properly align to given references. Bam file should be sorted by read name.
"""
from __future__ import print_function

import argparse
import sys
from signal import signal, SIGPIPE, SIG_DFL

import pysam


def rc(s):
    """Reverse complement DNA string s"""
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([basecomplement[c] for c in s][::-1])


def get_pairs_aligning_to_refs(bamfiles, ref_names_file, pair1_file, pair2_file, read_format):
    ref_names = dict([(s.strip(), 0) for s in open(ref_names_file).readlines()])

    for i, bamfile in enumerate(bamfiles):
        with pysam.Samfile(bamfile, 'rb') as bamh, open(pair1_file, 'w') as ph1, open(pair2_file, 'w') as ph2:
            ref_ids = [bamh.gettid(r) for r in ref_names]

            prev_read1 = ""
            prev_read2 = ""
            for read in bamh:
                # check if read is aligned
                try:
                    ref = bamh.getrname(read.tid)
                except ValueError:
                    continue

                # check if the read is aligned to given ref and if both reads
                # of the pair properly align to the same reference
                if read.tid in ref_ids and read.is_proper_pair and read.tid == read.mrnm:
                    if read.is_read1:
                        # make sure read is only printed once in case of
                        # multiple alignment
                        if read.qname != prev_read1:
                            print(read_format.format(qname=read.qname,
                                pairid=1,
                                seq=read.seq if not read.is_reverse else rc(read.seq),
                                qual=read.qual if not read.is_reverse else read.qual[::-1]), file=ph1)
                            prev_read1 = read.qname
                    else:
                        # make sure read is only printed once in case of
                        # multiple alignment
                        if read.qname != prev_read2:
                            print(read_format.format(qname=read.qname,
                                pairid=2,
                                seq=read.seq if not read.is_reverse else rc(read.seq),
                                qual=read.qual if not read.is_reverse else read.qual[::-1]), file=ph2)
                            prev_read2 = read.qname


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "bamfiles", nargs="+", help="Bam files\n")
    parser.add_argument("-r", "--ref_names", help="File with one ref name per line. Get those ids.\n")
    parser.add_argument("-1", "--pair1", help="File with pair1 reads\n")
    parser.add_argument("-2", "--pair2", help="File with pair2 reads\n")
    parser.add_argument("-f", "--read_format", default="@{qname}/{pairid}\n{seq}\n+\n{qual}",
    help="Python format string for read, you can use qname, pairid, seq and qual (default '@{qname}/{pairid}\\n{seq}\\n+\\n{qual}'")
    args = parser.parse_args()

    # ignore broken pipe error when piping output
    # http://newbebweb.blogspot.pt/2012/02/python-head-ioerror-errno-32-broken.html
    signal(SIGPIPE, SIG_DFL)

    get_pairs_aligning_to_refs(args.bamfiles, args.ref_names, args.pair1, args.pair2, args.read_format)
