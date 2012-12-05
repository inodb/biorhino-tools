#!/bin/bash
HELPDOC=$( cat <<EOF
Maps given single-end library to given reference with bwa and uses picard to
removes duplicates.

Usage:
    bash inodb-bwa-se-picard-mrkdup.sh [options] <single reads> <qname> <ref> <rname>
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
if [ "$#" -ne "4" ]
then
    echo "Invalid number of arguments: 4 needed but $# supplied" >&2
    echo "$HELPDOC"
    exit 1
fi
Q1=$1
QNAME=$2
REF=$3
RNAME=$4
echo "SE Reads: $Q1"
echo "Query name: $QNAME"
echo "Reference file: $REF"
echo "Reference name: $RNAME"

if [ ! -e ${REF}.bwt ]
then
    bwa index $REF
fi
bwa aln $REF $Q1 > ${QNAME}.sai
bwa samse -n 10 $REF ${QNAME}.sai $Q1 > ${RNAME}_${QNAME}.sam
samtools faidx $REF
samtools view -bt $REF.fai ${RNAME}_${QNAME}.sam > ${RNAME}_${QNAME}.bam
samtools sort ${RNAME}_${QNAME}.bam ${RNAME}_${QNAME}-s
samtools index ${RNAME}_${QNAME}-s.bam
java -jar /bubo/home/h16/inod/glob/downloaded_software/picard-tools-1.66/MarkDuplicates.jar \
    INPUT=${RNAME}_${QNAME}-s.bam \
    OUTPUT=${RNAME}_${QNAME}-smd.bam \
    METRICS_FILE=${RNAME}_${QNAME}-smd.metrics \
    AS=TRUE \
    VALIDATION_STRINGENCY=LENIENT \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
    REMOVE_DUPLICATES=TRUE
samtools sort ${RNAME}_${QNAME}-smd.bam ${RNAME}_${QNAME}-smds
samtools index ${RNAME}_${QNAME}-smds.bam
rm ${RNAME}_${QNAME}.sam ${RNAME}_${QNAME}.bam ${RNAME}_${QNAME}-smd.bam ${RNAME}_${QNAME}-s.bam ${RNAME}_${QNAME}-s.bam.bai ${QNAME}.sai
genomeCoverageBed -ibam ${RNAME}_${QNAME}-smds.bam > ${RNAME}_${QNAME}-smds.coverage
inodb-plot-genomeCoverageBed.R ${RNAME}_${QNAME}-smds.coverage
