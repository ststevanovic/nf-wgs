# FASTQC

**FASTQC** aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines, [see documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

## Requirements

- <a href="https://www.nextflow.io/">NextFlow</a>
- <a href="https://www.docker.com/">Docker</a>

## Running unnamed workflow (double-end test)

```
 nextflow -log logs/fastqc.log run modules/fastqc/fastqc.nf -c modules/fastqc/fastqc.config
```

## Running named workflow (single-end test)

```
 nextflow -log logs/fastqc.log run modules/fastqc/fastqc.nf -entry test_fastqc_single_end -c modules/fastqc/fastqc.config
```

## Pipeline parameters

```
nextflow run modules/fastqc/fastqc.nf --help
```

```
USAGE:

nextflow \
-log logs/fastqc.log \
run modules/fastqc/fastqc.nf -c modules/fastqc/fastqc.config \
--test_data read_r1.fq.gz

Arguments:
    --test_data     FILE(s)            Fastq(.gz) file(s) for read1 [, read2] (default: data/reads/*_dna_ercc_{1,2}.fq.gz)
    --outdir        DIR                Output directory(default: results/TEST/FastQC/)


====================================================================
```
