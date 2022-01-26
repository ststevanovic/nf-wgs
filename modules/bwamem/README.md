# Introduction

Bwamem is a pipeline which combines BWA and GATK tools.
It takes genome reference sequence and make index.
Steps...
It deduplicates...

## Additional description:

    GATK verison: v1.0
    BWA version:
    Software requirements: { checkout docker file } // read from yarn ?

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
