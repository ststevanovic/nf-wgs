# Introduction

Haplotyper module contains Genome Analysis Toolkit (GATK) HaplotypeCaller. 
It is used to call germline SNPs and indels via local re-assembly of haplotypes.

*Steps in this pipeline:*

1. Create index from input BAM file ( Samtools fai )
2. Create VCF and VCF-index file ( GATK HaplotypeCaller )

### Additional description:

GATK verison: v4.2.4.0 
[GATK HaplotypCaller docs](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller)  

Samtools version: 1.14  
[Samtools docs](http://www.htslib.org/doc/samtools.html)

### Inputs

- Required: 
* BAM - Input bam file(s) from which to make variant calls
* fasta - Reference genome 
* fai - indexed reference genome 
* dict - dictionary of the reference genome

- Optional:  
* outdir - Output directory
* input_intervals - intervals file

### Outputs

-   VCF file
-   VCF index file
-   reports

## Usage

The **haplotyper** nextflow script is located in /modules/haplotyper directory and config file is haplotyper.config in the same directory. 

If you don't specify --outdir, Nextflow will create results dir in main directory under results/TEST and place your results there. If you want to place your results somewhere else, provide additional path by using --outdir parameter in your Nextflow command line.

In case when the run is unsuccessful, you can fix your issues and then run it from the failure point, by using --resume parameter in your Nextflow command line.

# Build and test

1. Git clone the wgs repository [here](https://github.com/strahinja08/nf-wgs.git)
2. cd wgs/
3. Run the bellow command to test the Haplotyper (haplotyper_module workflow).

nextflow \   
*-log* logs/hap.log \  
*run* modules/haplotyper/haplotyper.nf \  
*-entry* haplotyper_module \  
*-c* modules/haplotyper/haplotyper.config \  
[*-fasta* path/to/fasta \  ]
[*-fai* path/to/reads \  ]
[*-dict* path/to/dict \ ]
[ *-\-resume* ]