# splicing-analysis

## Prerequisites

1. Install [docker CE](https://docs.docker.com/install/linux/docker-ce/ubuntu/).

## Building

```
docker image build -t splicing-analysis .
```

## Running

### Without Docker:

```
Usage: splicing-analysis <REF_DIR> <BAM_FILE> <JUNCTIONS_BED> <FPKM_FILE> <GENOME_VERSION> <OUT_SAMPLE_NAME>

        REF_DIR             Directory containing reference files
        BAM_FILE            BAM file containing reads
        JUNCTIONS_BED       BED file containing junctions
        FPKM_FILE           genes.fpkm file obtained from cufflinks
        GENOME_VERSION      Human Genome version prefix used in REF_DIR files
        OUT_SAMPLE_NAME     Name of output directory

```

## Development

```
docker run -it --rm -v $PWD:/home/ splicing-analysis
```

This maps the current directory to a folder inside the container. This makes it
so that there is no need to re-build the image upon making source code changes.
