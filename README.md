# WGS Nextflow Pipeline

A modular and scalable whole genome sequencing (WGS) pipeline implemented in Nextflow DSL2.  
This pipeline integrates standard tools for alignment, variant calling, structural variation detection, quality control, and annotation.

**Status**: Prototype  
**Nextflow version**: DSL2  
**Pipeline version**: 0.1.0

## Table of Contents

- [Introduction](#introduction)
- [Input](#input)
- [Output](#output)
- [Usage](#usage)
- [Parameters](#parameters)
- [Pipeline Steps](#pipeline-steps)
- [Software Tools](#software-tools)
- [Acknowledgements](#acknowledgements)
- [TODO](#todo)

## Introduction

This pipeline processes raw FASTQ reads and performs:

- Alignment with BWA-MEM  
- BAM sorting, duplicate marking, and BQSR with GATK  
- Germline variant calling with HaplotypeCaller  
- Structural variant calling with Manta  
- Quality metrics via FastQC and Qualimap  
- Annotation using SnpEff and VEP  
- Summary reporting with MultiQC

The pipeline is written using Nextflow DSL2 and follows a modular design for reusability and scalability.

## Input

| Input Type       | Description                                  | Example Path                                             |
|------------------|----------------------------------------------|----------------------------------------------------------|
| Paired-end FASTQ | Input reads directory                        | data/reads/*_{1,2}.fastq.gz                              |
| Reference Genome | FASTA + .fai + .dict                         | data/reference_genome/test.fa                            |
| Known Indels     | VCF + index (gold standard)                  | Mills_and_1000G_gold_standard.indels.hg19.vcf.gz         |
| dbSNP            | dbSNP VCF + index                            | data/dbsnp/test.dbsnp.vcf.gz                             |
| BAM/BAI          | BAM and index files for Manta and downstream tools | data/novak-bam-bai/*{.bam,.bam.bai}               |

## Output

| Output Directory         | Contents                                                 |
|--------------------------|----------------------------------------------------------|
| results/FASTQC/          | Raw read QC reports                                      |
| results/BWAmem/          | Aligned BAM files, sorted and deduplicated              |
| results/GATK/            | Recalibrated BAMs, BQSR reports, variant VCFs           |
| results/Manta/           | Structural variants (SVs) in VCF and BED                |
| results/Haplotyper/      | Germline SNPs/INDELs (GATK HaplotypeCaller)             |
| results/Qualimap/        | Coverage and alignment QC metrics                        |
| results/SnpEff/          | Annotated VCF using SnpEff                               |
| results/VEP/             | Annotated VCF using VEP                                  |
| results/MultiQC/         | Summary QC report (multiqc_report.html)                 |

## Usage

1. Clone the pipeline

    git clone https://github.com/your-org/wgs-nextflow-pipeline.git
    cd wgs-nextflow-pipeline

2. Prepare your inputs

    - Place FASTQ files in `data/reads/`
    - Reference genome in `data/reference_genome/`
    - dbSNP and known indels in appropriate subdirectories

3. Run the pipeline

    nextflow run main.nf -profile standard

## Parameters

These parameters can be set in `nextflow.config` or passed at runtime:

| Parameter             | Description                                        | Default                                      |
|-----------------------|----------------------------------------------------|----------------------------------------------|
| --reads_dir           | Path to paired-end reads                           | data/reads                                   |
| --fasta               | Reference genome FASTA                             | data/reference_genome/test.fa                |
| --fai                 | FASTA index file                                   | data/reference_genome/test.fa.fai            |
| --dict                | Reference genome dictionary file                   | results/TEST/GATK/test.dict                  |
| --dbsnp_gz            | dbSNP VCF file                                     | data/dbsnp/test.dbsnp.vcf.gz                 |
| --dbsnp_idx_gz        | Index for dbSNP VCF                                | data/dbsnp/test.dbsnp.vcf.gz.tbi             |
| --bambai              | BAM + BAI files for Manta input                    | data/novak-bam-bai/*{.bam,.bam.bai}          |
| --inpdir              | Directory containing BAMs for Qualimap             | results/TEST/Haplotyper                      |
| --vcf                 | Input VCF files for annotation                     | results/TEST/*/*.vcf                         |
| --snpeff_db           | SnpEff DB version                                  | GRCh38.86                                    |
| --snpeff_cache        | Path to SnpEff cache                               | null                                         |
| --genome              | Genome build name (used by VEP, SnpEff)            | GRCh38                                       |
| --with_tumor          | Enable tumor-normal Manta mode                     | false                                        |
| --publish_dir_mode    | Output publishing mode                             | copy                                         |

## Pipeline Steps

- Quality Control
    - FastQC
- Alignment
    - BWA-MEM
    - Samtools Sort
- Post-alignment Processing
    - GATK MarkDuplicates
    - GATK BaseRecalibrator
    - GATK ApplyBQSR
- Variant Calling
    - GATK HaplotypeCaller
    - Manta (for structural variation)
- Quality Assessment
    - Qualimap
- Annotation
    - SnpEff
    - VEP
- Reporting
    - MultiQC

## Software Tools

| Tool         | Function                                 |
|--------------|------------------------------------------|
| FastQC       | Raw read QC                              |
| BWA-MEM      | Read alignment                           |
| Samtools     | BAM sorting and indexing                 |
| GATK         | Duplicate marking, BQSR, variant calling |
| Manta        | Structural variant detection             |
| Qualimap     | QC metrics on BAM files                  |
| SnpEff       | Variant effect annotation                |
| VEP          | Variant annotation (Ensembl)             |
| MultiQC      | Aggregated summary report                |

## Acknowledgements

- Inspired by nf-core community best practices  
- Uses open-source bioinformatics tools  
- Contributions welcome

## TODO

- Emit `.dict` and `.bai` automatically if missing  
- Add tumor-normal pairing logic for Manta  
- Auto-generate Tabix indexes if absent  
- Add optional SnpEff cache support  
