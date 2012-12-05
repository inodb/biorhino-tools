#!/bin/bash
HELPDOC=$( cat <<EOF
Converts sam file to sorted bam file for given sam file and reference fasta.

Usage:
    bash [options] indob-sam2sortedbam.sh <reference fasta> <sam file>
Options:
    -h      This help documentation.
EOF
) 
 
# Halt on error
set -e

# Parse options
while getopts ":h" opt; do
    case $opt in
        h)
            echo "$HELPDOC"
            exit 0
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            echo "$HELPDOC"
            exit 1
            ;;
    esac
done

# Parse arguments
if [ "$#" -ne "2" ]
then
    echo "Invalid number of arguments: 2 needed but $# supplied" >&2
    echo "$HELPDOC"
    exit 1
fi
REF=$1
SAM=$2
echo "Reference file: $REF"
echo "Sam file: $SAM"

samtools faidx $REF
BASE=${SAM: 0: ${#SAM} - 4}
samtools view -bt $REF.fai $SAM > ${BASE}.bam
samtools sort ${BASE}.bam ${BASE}-sorted
samtools index ${BASE}-sorted.bam
rm ${BASE}.bam
