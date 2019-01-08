#!/usr/bin/env python

"""
Performs intron_analysis of RNA-seq data given
appropriate coverage files.

NOTE:
 To be used exclusively with reference files created by
 the `making_all_reference_files.sh` and coverage files
 created by `compute_coverage.sh`.

NOTE:
 Use build_*_junctions_ref() in utils.py to generate
 exon_exon_junctions or exon_intron_junctions.

NOTE:
 .gtf & .gff are 1-based while
 .bed & .bam are 0-based. This matters when
 using bedtools with different formats.
"""

import subprocess as sp

import argparse
import os
import os.path
import util

# Define script arguments

parser = argparse.ArgumentParser(
    description="Performs intron_analysis given appropriate coverage files.",
    epilog="Usage eg: python intron_analysis.py <ref_dir> <genes.fpkm> <genome_version> <sample_name>")

parser.add_argument("ref_dir", help="Directory containing reference files.")
parser.add_argument("fpkm_arg", help="FPKM file.")
parser.add_argument("hg_arg", help="Genome version.")
parser.add_argument("sample_name", help="Sample name used in compute_coverage.sh. \
                     Directory of coverage files must be in the same dir as this  \
                     script and must be named ${sample_name}_cvg")

args = parser.parse_args()

# Open coverage files
cvg_base = args.sample_name + "_cvg" + "/" + args.sample_name
exons_cvg = open(cvg_base + "_exon_cvg.bed")
introns_cvg = open(cvg_base + "_intron_cvg.bed")
exon_exon_cvg = open(cvg_base + "_exon_exon_cvg.bed")
exon_wao_cvg = open(cvg_base + "_exon_wao_cvg.bed")
intron_wao_cvg = open(cvg_base + "_intron_wao_cvg.bed")
total_cvg = open(cvg_base + "_full_cvg.bed")
chr_file = open(args.sample_name + "_cvg" + "/" + args.hg_arg + "_only_exons_refGene_longest_sorted.bed")

print ">> Processing sample " + args.sample_name
print ">> Running gene-id mapper..."

geneid_map = open(args.ref_dir + "/" + args.hg_arg + ".gene_name_refseqid")

# Uncomment below for ID'ing exons and exon-exon junctions.
# util.append_id(geneid_map, exons_f, 9)
# util.append_id(geneid_map, exon_exon_2junc_f, 21)
introns_id_f = util.append_id(geneid_map, introns_cvg, 9)

print ">> Building exon tables..."
exon_t = util.build_cvg_table(exons_cvg)
print ">> Building exon-exon junction tables..."
exon_exon_t = util.build_junc_table(exon_exon_cvg)
print ">> Building exon-intron junction tables..."
exon_intron_t = util.build_exon_intron_table(exon_wao_cvg, intron_wao_cvg, chr_file)
print ">> Building gene_fpkm tables..."
fpkm_t = util.build_fpkm_table_from_stringtie(open(args.fpkm_arg))
print ">> Building gene length tables..."
gene_length_t = util.build_gene_length_table(open("%s/%s_gene_length.txt" % (args.ref_dir, args.hg_arg)))
print ">> Building mrna tables..."
mrna_t = util.build_mrna_table(open("%s/%s_select_feature_lengths" % (args.ref_dir, args.hg_arg)))
print ">> Merging tables into output file..."
out_file = open(args.sample_name + "_intron_analysis.txt", "w+")
util.do_analysis(introns_id_f, exon_t, exon_exon_t, exon_intron_t, fpkm_t, gene_length_t, mrna_t, out_file, args.sample_name)
print ">> Computing total coverage..."
total_output = open(args.sample_name + "_total_cvg.txt", "w+")
util.compute_total_cvg(total_cvg, total_output)

print ">> %s intron analysis finished. Come again! (*^_^)/" % args.sample_name
