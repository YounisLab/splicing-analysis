#!/usr/bin/env nextflow

params.help = null
params.ref_dir = null
params.bam_file = null
params.junc_bed = null
params.fpkm = null
params.genome = null
params.sample_name = null

def helpMessage() {
    log.info """
    Usage: nextflow run splicing-analysis.nf --ref_dir <REF_DIR> --bam_file <BAM_FILE> --junc_bed <JUNCTIONS_BED> ---fpkm <FPKM_FILE> --genome <GENOME_VERSION> --sample_name <OUT_SAMPLE_NAME>

        REF_DIR             Directory containing reference files
        BAM_FILE            BAM file containing reads
        JUNCTIONS_BED       BED file containing junctions
        FPKM_FILE           gene_abund.tab file obtained from stringtie
        GENOME_VERSION      Human Genome version prefix used in REF_DIR files
        OUT_SAMPLE_NAME     Name of output directory

eg: nextflow run splicing-analysis.nf --ref_dir /home/data/hg38_ref/ --bam_file ./accepted_hits.bam --junc_bed ./junctions.bed --fpkm ./gene_abund.tab --genome hg38 --sample_name MCF7
    """
}

if (params.help) {
    helpMessage()
    exit 0
}

if (!params.ref_dir) {
    helpMessage()
    exit 1, "REF_DIR not specified."
}

if (!params.bam_file || !params.junc_bed || !params.fpkm) {
    helpMessage()
    exit 1, "Missing one or more input files."
}

if (!params.genome) {
    helpMessage()
    exit 1, "Please specify GENOME_VERSION."
}

if (!params.sample_name) {
    helpMessage()
    exit 1, "Please specify OUT_SAMPLE_NAME."
}

process analysis {

    input:
     file ref_dir from Channel.fromPath(params.ref_dir).collect()
     file bam_file from Channel.fromPath(params.bam_file)
     file junc_bed from Channel.fromPath(params.junc_bed)
     file fpkm from Channel.fromPath(params.fpkm)

    output:
    file("${params.sample_name}_intron_analysis.txt") into ANALYSIS_DIR_1
    file("${params.sample_name}_total_cvg.txt") into ANALYSIS_DIR_2

    publishDir "$baseDir", mode: 'copy'

    script:
    """
    echo "===> Computing coverage..."
    compute_coverage.sh $ref_dir $bam_file $junc_bed $params.genome $params.sample_name
    echo "===> Performing analysis..."
    analyze.py $ref_dir $fpkm $params.genome $params.sample_name
    """
}

workflow.onComplete {
    log.info "Analysis completed."
}
