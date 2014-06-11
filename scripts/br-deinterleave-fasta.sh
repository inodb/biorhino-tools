#!/bin/bash
# Author: nathanhaigh 
# Edit to fasta: inodb
#
# Usage: deinterleave_fasta.sh < interleaved.fasta f.fastq r.fastq
# 
# Deinterleaves a FASTA file of paired reads into two FASTA
# files specified on the command line.
# 
# Original code for fastq: https://gist.github.com/3521724

paste - - - -  | tee >(cut -f 1-2 | tr "\t" "\n" > $1) | cut -f 3-4 | tr "\t" "\n" > $2
