#!/usr/bin/env bash

set -e -u

if [ "$#" -ne 6 ]
then
    echo "Usage: splicing-analysis <REF_DIR> <BAM_FILE> <JUNCTIONS_BED> <FPKM_FILE> <GENOME_VERSION> <OUT_SAMPLE_NAME>

        REF_DIR             Directory containing reference files
        BAM_FILE            BAM file containing reads
        JUNCTIONS_BED       BED file containing junctions
        FPKM_FILE           genes.fpkm file obtained from cufflinks
        GENOME_VERSION      Human Genome version prefix used in REF_DIR files
        OUT_SAMPLE_NAME     Name of output directory

eg: splicing-analysis /home/data/hg38_ref/ ./accepted_hits.bam ./junctions.bed ./genes.fpkm hg38 MCF7
    "
    exit
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

echo "======> Computing coverage..."
$DIR/src/compute_coverage.sh $1 $2 $3 $5 $6
echo "======> Analyzing splice sites..."
python $DIR/src/analyze.py $1 $4 $5 $6
