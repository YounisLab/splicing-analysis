# splicing-analysis

## Prerequisites

1. Install [docker CE](https://docs.docker.com/install/linux/docker-ce/ubuntu/).

## Building

```
docker image build -t splicing-analysis .
```

## Running

```
Usage: nextflow run splicing-analysis.nf --ref_dir <REF_DIR> --bam_file <BAM_FILE> --junc_bed <JUNCTIONS_BED> ---fpkm <FPKM_FILE> --genome <GENOME_VERSION> --sample_name <OUT_SAMPLE_NAME>

        REF_DIR             Directory containing reference files
        BAM_FILE            BAM file containing reads
        JUNCTIONS_BED       BED file containing junctions
        FPKM_FILE           gene_abund.tab file obtained from stringtie
        GENOME_VERSION      Human Genome version prefix used in REF_DIR files
        OUT_SAMPLE_NAME     Name of output directory

eg: nextflow run splicing-analysis.nf --ref_dir /home/data/hg38_ref/ --bam_file ./accepted_hits.bam --junc_bed ./junctions.bed --fpkm ./gene_abund.tab --genome hg38 --sample_name MCF7

```
