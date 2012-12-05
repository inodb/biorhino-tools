#!/bin/bash
HELPDOC=$( cat <<EOF
Maps given paired library to given reference with bwa and uses samtools
to removes duplicates.

Usage:
    bash `basename $0` [options] <reads1> <reads2> <qname> <ref> <rname>
Options:
    -k      Keep all output from intermediate steps.
    -h      This help documentation.
EOF
) 
 
# Halt on error
set -e

# Default parameters
RMTMPFILES=true

# Parse options
while getopts "kh" opt; do
    case $opt in
        k)
            RMTMPFILES=false
            ;;
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

# Index reference, Burrows-Wheeler Transform
if [ ! -e ${REF}.bwt ]
then
    bwa index $REF
fi

# Find SA coordinates in the reads (whatever that may be)
bwa aln $REF -1 $Q1 > ${QNAME}1.sai
bwa aln $REF -2 $Q2 > ${QNAME}2.sai

# Align Paired end and bam it
bwa sampe -a 3000 $REF ${QNAME}1.sai ${QNAME}2.sai $Q1 $Q2 > ${RNAME}_${QNAME}.sam
samtools faidx $REF
samtools view -bt $REF.fai ${RNAME}_${QNAME}.sam > ${RNAME}_${QNAME}.bam
samtools sort ${RNAME}_${QNAME}.bam ${RNAME}_${QNAME}-s
samtools index ${RNAME}_${QNAME}-s.bam

# Mark duplicates and sort
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

# Determine Genome Coverage
genomeCoverageBed -ibam ${RNAME}_${QNAME}-smds.bam > ${RNAME}_${QNAME}-smds.coverage

# Determine insert sizes
java -jar /bubo/home/h16/inod/glob/downloaded_software/picard-tools-1.66/CollectInsertSizeMetrics.jar \
    HISTOGRAM_FILE=${RNAME}_${QNAME}-smds.inserthist \
    INPUT=${RNAME}_${QNAME}-smds.bam \
    OUTPUT=${RNAME}_${QNAME}-smds.insertmetrics \
    AS=TRUE \
    VALIDATION_STRINGENCY=LENIENT

# Plot insert sizes and coverage
inodb-plot-CollectInsertSizeMetrics.R ${RNAME}_${QNAME}-smds.insertmetrics
inodb-plot-genomeCoverageBed.R ${RNAME}_${QNAME}-smds.coverage

# Remove temp files
if $RMTMPFILES
then
    rm ${RNAME}_${QNAME}.sam \
       ${RNAME}_${QNAME}.bam \
       ${RNAME}_${QNAME}-smd.bam \
       ${RNAME}_${QNAME}-s.bam \
       ${RNAME}_${QNAME}-s.bam.bai \
       ${QNAME}1.sai \
       ${QNAME}2.sai
fi
