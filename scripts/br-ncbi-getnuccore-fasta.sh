#!/bin/bash
set -e
HELPDOC=$( cat <<EOF
Downloads fasta files of given NCBI IDs from the nuccore db and pipes them to
stdout.

Usage:
    bash [options] inodb-ncbi-getnuccore-fasta.sh [nuccore idi ...]
Options:
    -h      This help documentation.
EOF
) 

urls=''
for id in ${@:1}
do
    urls="${urls}http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&id=${id} "
done

curl $urls
