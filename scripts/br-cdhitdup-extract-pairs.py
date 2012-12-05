#!/usr/bin/env python
"""The sequence duplicate removal tool cd-hit-dup concatenates pairs and uses
the concatenated pairs to cluster. It outputs the concatenated pairs instead of
the actual pairs. To extract the actualy pairs, this script can be used. Supply
the original pairs and the fasta file that cd-hit-dup outputs.

Usage:
cdhitdup-extract-pairs <pairstoextract.fastq> <idsfromcdhit.fastq>
"""
import sys
import getopt
import screed

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def get_ids(fastfile, remove_pair_info=True):
    ids = {}

    for record in screed.open(fastfile):
        if remove_pair_info:
            ids[record.name[:-2]] = 1
        else:
            ids[record.name] = 1
    
    return ids

def print_ids(fastqfile, ids, remove_pair_info=True):
    for record in screed.open(fastqfile):
        if remove_pair_info:
            try:
                ids[record.name[:-2]]
            except KeyError:
                continue
            print "@%s\n%s\n+\n%s" % (record.name, record.sequence,
                    record.accuracy)
        if not remove_pair_info:
            try:
                ids[record.name]
            except KeyError:
                continue
            print "@%s\n%s\n+\n%s" % (record.name, record.sequence,
                    record.accuracy)
            
    return 0

def extract_fast_ids(fastqfile, idsfastqfile):
    ids = get_ids(idsfastqfile)
    if len(ids) == 0:
        print >>sys.stderr, "No ids found"
        return 2
    print_ids(fastqfile, ids)
    
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
    if (len(args) == 2):
        return extract_fast_ids(fastqfile=args[0], idsfastqfile=args[1])
    else:
        print >>sys.stderr, "Two arguments required"
        print __doc__
        return 2

if __name__ == "__main__":
    sys.exit(main())
