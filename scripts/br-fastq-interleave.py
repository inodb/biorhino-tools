#!/usr/bin/env python
import argparse
import sys


def main(pair1, pair2):
    fp1 = open(pair1)
    fp2 = open(pair2)

    for line1 in fp1:
        sys.stdout.write(line1)
        for i in range(3):
            sys.stdout.write(next(fp1))
        for i in range(4):
            sys.stdout.write(next(fp2))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("pair1", help="First of pair of fastq reads\n")
    parser.add_argument("pair2", help="Second of pair of fastq reads\n")
    args = parser.parse_args()
    main(args.pair1, args.pair2)
