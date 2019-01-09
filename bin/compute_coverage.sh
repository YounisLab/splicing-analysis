#!/usr/bin/env bash

# The following files are used from ref_dir:
# - only_exons_refGene_longest.bed
# - only_introns_refGene_longest.bed
# - exon_exon_junctions.bed
# - exon_intron_junctions.bed
# - refGene_full

# The original `making_all_reference_files.sh`
# saves files in the GFF format.
# If the corresponding .bed files do not exist,
# use the `gff2bed` utility to convert.

set -e -u

if [ "$#" -ne "5" ]
then
    echo "Usage: compute_coverage.sh <ref_dir> <aligned BAM file> <junctions BED file> <genome_version> <sample_name>"
    echo "eg: compute_coverage.sh /home/data/hg38_ref/ ./accepted_hits.bam ./junctions.bed hg38 MCF7"
    exit 1
fi

sample=$5
wd=${5}_cvg
mkdir -p $wd

sort_using_bam () {
    # Args: BED file to sort

    bname=`basename $1 .bed`
    echo ">> Sorting BED file $1..."
    bedtools sort -faidx $wd/${sample}_chr_order.txt -i $1 > $wd/${bname}_sorted.bed
}

echo ">> Operating on sample $5"
echo ">> Generating BAM index..."
samtools index $2 ${2}.bai
echo ">> Generating chromosomal order from BAM file..."
chr=$wd/${5}_chr_order.txt
samtools idxstats $2 | cut -f 1-2 > $chr
sort_using_bam $1/${4}_only_exons_refGene_longest.bed
sort_using_bam $1/${4}_only_introns_refGene_longest.bed
sort_using_bam $1/${4}_exon_exon_junctions.bed
sort_using_bam $1/${4}.refGene_full.bed

echo ">> Running intersectBed on exons..."
intersectBed -a $wd/${4}_only_exons_refGene_longest_sorted.bed -b $2 -sorted -g $chr -F 1 -c > $wd/${5}_exon_cvg.bed
echo ">> Running intersectBed on introns..."
intersectBed -a $wd/${4}_only_introns_refGene_longest_sorted.bed -b $2 -sorted -g $chr -F 1 -c > $wd/${5}_intron_cvg.bed
echo ">> Running intersectBed on exon-exon junctions..."
intersectBed -a $wd/${4}_exon_exon_junctions_sorted.bed -b $3 -F 0.8 -loj> $wd/${5}_exon_exon_cvg.bed
echo ">> Running intersectBed on exons with -wao..."
intersectBed -a $wd/${4}_only_exons_refGene_longest_sorted.bed -b $2 -sorted -g $chr -wao> $wd/${5}_exon_wao_cvg.bed
echo ">> Running intersectBed on introns with -wao..."
intersectBed -a $wd/${4}_only_introns_refGene_longest_sorted.bed -b $2 -sorted -g $chr -wao> $wd/${5}_intron_wao_cvg.bed
echo ">> Running intersectBed on refGene_full"
intersectBed -a $wd/${4}.refGene_full_sorted.bed -b $2 -sorted -g $chr -c > $wd/${5}_full_cvg.bed

echo ">> Coverage on $5 completed."
