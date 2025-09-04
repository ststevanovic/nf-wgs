# Introduction

BWAmem is a workflow with BWA, Samtools and Genome Analysis Toolkit (GATK) combined.
Burrows-Wheeler Aligner (BWA) is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome. It consists of three algorithms: BWA-backtrack, BWA-SW and BWA-MEM. Samtools is a set of utilities that manipulate alignments in the SAM (Sequence Alignment/Map), BAM, and CRAM formats. It converts between the formats, does sorting, merging and indexing, and can retrieve reads in any regions swiftly. The Genome Analysis Toolkit (GATK) is a set of bioinformatic tools for analyzing high-throughput sequencing (HTS) and variant call format (VCF) data. The toolkit is well established for germline short variant discovery from whole genome and exome sequencing data.

*Steps in this pipeline:*

1. Create index from input reference genome ( Samtools fai )
2. Map reads to reference genome ( BWA mem align ) and sort aligned reads ( Samtools sort )
3. Markduplicates ( GATK MarkDuplicates )
4. Recalibrate ( GATK BaseRecalibrator )
5. Create BQSR ( GATK BQSR )

### Additional description:

GATK verison: v4.2.4.0
[GATK MarkDuplicates docs](https://gatk.broadinstitute.org/hc/en-us/community/posts/360076320032-MarkDuplicates)  
[GATK BaseRecalibrator docs](https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator)

BWA version:  0.7.17  
[BWA docs](http://bio-bwa.sourceforge.net/bwa.shtml)

Samtools version: 1.14  
[Samtools docs](http://www.htslib.org/doc/samtools.html)

### Inputs:

- Required: reads, fasta, fasta_index
- Optional: dbsnp_gz, dbsnp_idx_gz, golden_indel_gz, golden_indel_idx_gz,

### Outputs:

- aligned.sorted.MarkedDupicates.BQSR.bam
- ${name}_reacal_data.table 
- reports

## Usage

The **bwamem** Nextflow script is located in /modules/bwa directory and config file is bwamem.config directory. 

If you don't specify --output_dir, Nextflow will create results dir in main directory under results/TEST and place your results there. If you want to place your results somewhere else, provide additional path by using --output_dir parameter in your Nextflow command line.

In case when the run is unsuccessful, you can fix your issues and then run it from the failure point, by using -resume parameter in your Nextflow command line.
# Build and test

1. Git clone the wgs repository [here](https://github.com/strahinja08/nf-wgs.git)
2. cd wgs/
3. Run the bellow command to test the BWAmem workflow.

nextflow \   
*-log* logs/bwa.log \  
*run* modules/bwamem/bwa.nf \  
*-entry* gatk_workflow \  
*-c* modules/bwamem/bwa.config \  
*-genome* path/to/genome \  
*-reads* path/to/reads \  
[ *-\-dbsnp_gz* /path/dbsnp_gz ] \  
[ *-\-dbsnp_idx_gz* /path/to/dbsnp_idx_gz ] \  
[ *-\-golden_indel* /path/to/golden_indel ] \  
[ *-\-gz golden_indel_idx_gz* /path/to/golden_indel_idx_gz ] \  
[ *-\-resume* ]
