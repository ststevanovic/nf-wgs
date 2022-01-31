# Introduction

Bwamem is a pipeline which combines BWA, Samtools and Genome Analysis Toolkit (GATK).
BWA is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome. It consists of three algorithms: BWA-backtrack, BWA-SW and BWA-MEM. In this pipeline, first step is BWA index takes genome reference as an input and make index files. The output is sorted using Samtools.
It then takes these files for

It deduplicates...

## Additional description:

    GATK verison: v1.0
    BWA version:
    Software requirements: { checkout docker file } // read from yml ?

## Inputs

Requried reads_ch, fasta_ch, BWAINDEX.out.bwa_index
Optional: dbsnp_gz dbsnp_idx_gz golden_indel_gz golden_indel_idx_gz

## Outputs

aligned_sorted_MarkDup_bqsr.bam, ${name}\_recal_data.table \*data.table

## Usage

nextflow -log logs/bwa.log run modules/bwamem/bwa.nf \
-entry gatk*workflow \
-c modules/bwamem/bwa.config \
-genome path/to/genome \
-reads path/to/reads \
[--dbsnp_gz /path/dbsnp_gz ] \
[--dbsnp_idx_gz = /path/to/dbsnp_idx_gz ] \
[ --golden_indel* = /path/to/golden_indel ] \
[ --gz golden_indel_idx_gz = /path/to/golden_indel_idx_gz ] \
[--resume ]

# Build and test

Git clone the wgs repository here (development):
cd wgs/
Run the usage command to test the bwamem workflow.
