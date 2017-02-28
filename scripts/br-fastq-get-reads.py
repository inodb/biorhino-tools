#!/usr/bin/env python
"""
Use paste to make single line

E.g. like:

time for f in ext_files/CleanData/J538_MCF7/J538_MCF7_USD16070002_HT3FKCCXX_L2_1.clean.fq.gz; do 
    paste <(zcat $f | paste - - - -) <(zcat ${f/_1.clean/_2.clean} | paste - - - - ) | \
        python ~/git/biorhino-tools/scripts/br-fastq-get-reads.py <(cat asm/unaligned_reads.txt asm/read_names_pwb15.txt | sort -u) - | \
    tee >(cut -f 1-4 | tr "\t" "\n" > asm/$(basename $f .fq.gz).fastq) | cut -f 5-8 | tr "\t" "\n" > asm/$(basename ${f/_1.clean/_2.clean} .fq.gz).fastq;
done
"""
import argparse
import os
import sys

def get_reads(key_file, pair_record_file):
    keys = {}
    with open(key_file) as key_fileh:
        for line in key_fileh:
            key = line.rstrip('\n')
            keys[key] = 1
    keys = set(keys.keys())

    for line in sys.stdin:
        if line.split()[0][1:] in keys:
            print line.rsplit('\n')[0]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("keys", help="File with one read name per line. Get those ids.\n")
    parser.add_argument("pair_record", help="Single line with all reads info using paste")
    args = parser.parse_args()

    get_reads(args.keys, args.pair_record)
