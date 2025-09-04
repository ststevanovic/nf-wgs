# Samtools
Samtools is a suite of programs for interacting with high-throughput sequencing data. This module
uses Samtools for indexing reference genome file  and sorting BAM file. It produces comprehensive statistics from alignment file. Samtools use HTSlib internally, but these source packages contain their own copies of htslib so they can be built independently.
## Requirements

- <a href="https://www.nextflow.io/">NextFlow</a>
- <a href="https://www.docker.com/">Docker</a>

## Running module's named workflow

```
nextflow \
        -log logs/samtools.log \                                  
        run modules/Samtools/samtools.nf \                                 
        -entry samtools_fa_index_report \
        -c modules/Samtools/samtools.config \
        [--resume]
```

## Pipeline parameters

```
nextflow run modules/samtools/samtools.nf --help
```

```
Arguments:
    --fasta         FILE               Fasta reference genome ( default: data/reference_genome/test.fa )
    --outdir        DIR                Output directory ( default: results/TEST/Samtools/ )
    --bam           FILE               Bam input for reporting stats (default: results/TEST/Samtools/bam_sored.aligned.bam)

====================================================================
```
