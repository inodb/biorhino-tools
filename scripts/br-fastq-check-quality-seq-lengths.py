#!/usr/bin/env python
"""Checks if fastq quality lengths are equal to the sequence lengths.

Usage:
fastq-check-quality-seq-lengths <fastqfile.fastq>
"""
import sys
import getopt
import screed

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def check_quality_seq_lengths(fastqfile):
    for record in screed.open(fastqfile):
        if len(record.sequence) != len(record.accuracy):
            print ">%s\n%s\n+\n%s" % (record.name, record.sequence,
                    record.accuracy) 

    return 0

def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help"])
        except getopt.error, msg:
             raise Usage(msg)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "for help use --help"
        return 2
    # process options
    for o, a in opts:
        if o in ("-h", "--help"):
            print __doc__
            return 0
    # process arguments
    if (len(args) == 1):
        return check_quality_seq_lengths(fastqfile=args[0])
    else:
        print >>sys.stderr, "One argument required"
        print __doc__
        return 2

if __name__ == "__main__":
    sys.exit(main())
