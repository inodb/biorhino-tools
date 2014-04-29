#!/usr/bin/env python
"""Cut up fasta file in non-overlapping or overlapping parts.

Usage:
cut-up-fasta.py <contigs.fasta> [contigs.fasta ...]
Options:
    -l INT   Cut contigs up in chunks of given parameter [1999]
    -o INT   Overlap between successive chunks. The last o bases of a chunk
             overlap with the first o bases of the next chunk. Has to be
             smaller than l. [1900]
"""
import sys
import getopt
from Bio import SeqIO


def cut_up_fasta(fastfiles, chunk_size, overlap):
    for ff in fastfiles:
        for record in SeqIO.parse(ff, "fasta"):
            if (len(record.seq) > chunk_size):
                i = 0
                for split_seq in chunks(record.seq, chunk_size, overlap):
                    print ">%s.%i\n%s" % (record.id, i, split_seq)
                    i = i + 1
            else:
                print ">%s\n%s" % (record.id, record.seq)

    return 0


def chunks(l, n=1999, o=1900, merge_last=True):
    """ Yield successive n-sized chunks from l with given overlap o between the
    chunks.
    """
    assert n > o

    if not merge_last:
        for i in xrange(0, len(l), n - o):
            yield l[i:i + n]
    else:
        for i in xrange(0, len(l) - n + 1, n - o):
            yield l[i:i + n] if i + n + n - o <= len(l) else l[i:]


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hl:o:", ["help"])
        except getopt.error, msg:
            raise Usage(msg)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "for help use --help"
        return 2
    # process options
    chunk_size = 1000
    overlap = 900
    for o, a in opts:
        if o in ("-h", "--help"):
            print __doc__
            return 0
        if o in ("-l"):
            chunk_size = int(a)
        if o in ("-o"):
            overlap = int(a)
    if overlap >= chunk_size:
        print >>sys.stderr, "Overlap not smaller than chunk size"
        return 2
    # process arguments
    if (len(args) >= 1):
        return cut_up_fasta(args, chunk_size, overlap)
    else:
        print >>sys.stderr, "At least one argument required"
        print __doc__
        return 2

if __name__ == "__main__":
    sys.exit(main())
