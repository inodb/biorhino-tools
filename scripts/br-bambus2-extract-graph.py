#!/usr/bin/env python
"""Get the bambus2 dot graph representation from the contigs in given fasta or
fastq file(s). The entire dot graph output is usually to large to view with a
dot viewer which makes subsetting of the file a necessity. Obviously the
provided contigs should be a subset of the contigs used to construct the
bambus2 scaffold graph. The subsetted dotfile is printed to stdout.

Usage:
bambus2-extract-graph <bambus2output.agp> <bambus2output.dot> <contigs.fast(q|a)> [contigs.fast(q|a) ...]
"""
import sys
import getopt
import screed

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def get_contig_ids(fastqafile):
    contig_ids = []

    for record in screed.open(fastqafile):
        contig_ids.append(record.name)
    
    return contig_ids

def subset_dot_file(ids, agpfile, bambusdotfile, inputtype="contigfile"):
    if inputtype == "contigfile":
        contig_ids = []
        for c in ids:
            contig_ids.extend(get_contig_ids(c))
        if len(contig_ids) == 0:
            print >>sys.stderr, "No contigs found"
            return 2

        # Get corresponding scaffolds
        scaffold_ids = []
        for line in open(agpfile):   
            fields = line.split()                                        
            if (fields[5] in contig_ids):
                scaffold_ids.append(fields[0])                          
    # Print first 9 lines, contains heading of dotfile
    bdfd = open(bambusdotfile)
    for n in xrange(9):
        try:
            line = bdfd.next()
        except StopIteration:
            print >>sys.stderr, "Dot file seems to be empty"
            return 2
        sys.stdout.write(line)

    
    # Extract subgraphs that contain requested contigs
    print_lines = False
    while True:
        try:
            line = bdfd.next()
        except StopIteration:
            break
        if line.startswith("subgraph"):
            if line[line.find("_") + 1:line.find("{")] in scaffold_ids:
                print_lines = True
        if line.startswith("}") and print_lines:
            sys.stdout.write(line)
            print_lines = False
        if print_lines:
            sys.stdout.write(line)

    # Close opening graph bracket
    sys.stdout.write("}\n")

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
    if (len(args) >= 3):
        return subset_dot_file(agpfile=args[0], bambusdotfile=args[1], ids=args[2:])
    else:
        print >>sys.stderr, "Three arguments required"
        print __doc__
        return 2

if __name__ == "__main__":
    sys.exit(main())
