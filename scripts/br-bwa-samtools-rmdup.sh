#!/bin/bash
set -e
HELPDOC=$( cat <<EOF
Maps given paired library to given reference with bwa and uses samtools
to removes duplicates.

Usage:
    bash [options] bwa-rmdup.sh <reads1> <reads2> <qname> <ref> <rname>
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
if [ "$#" -ne "5" ]
then
    echo "Invalid number of arguments: 5 needed but $# supplied" >&2
    echo "$HELPDOC"
    exit 1
fi
Q1=$1
Q2=$2
QNAME=$3
REF=$4
RNAME=$5
echo "Reads 1: $Q1"
echo "Reads 2: $Q2"
echo "Query name: $QNAME"
echo "Reference file: $REF"
echo "Reference name: $RNAME"

if [ ! -e ${REF}.bwt ]
then
    bwa index $REF
fi
bwa aln $REF -1 $Q1 > ${QNAME}1.sai
bwa aln $REF -2 $Q2 > ${QNAME}2.sai
bwa sampe -a 3000 $REF ${QNAME}1.sai ${QNAME}2.sai $Q1 $Q2 > ${RNAME}_${QNAME}.sam
samtools faidx $REF
samtools view -bt $REF.fai ${RNAME}_${QNAME}.sam > ${RNAME}_${QNAME}.bam
samtools sort ${RNAME}_${QNAME}.bam ${RNAME}_${QNAME}-sorted
samtools index ${RNAME}_${QNAME}-sorted.bam
samtools rmdup -S ${RNAME}_${QNAME}-sorted.bam ${RNAME}_${QNAME}-sorted-rmdup.bam
samtools index ${RNAME}_${QNAME}-sorted-rmdup.bam
