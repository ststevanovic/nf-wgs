# Manta

## Introduction
Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads. It is optimized for analysis of germline variation in small sets of individuals and somatic variation in tumor/normal sample pairs. Manta discovers, assembles and scores large-scale SVs, medium-sized indels and large insertions within a single efficient workflow.

## Input params
    
    --fasta             = Reference genome (default: "$rootDir/data/reference_genome/test.fa")
    --fai               = Reference fasta index (default: "$rootDir/data/reference_genome/test.fa.fai")
    
    --normal_bam          = Normal BAM file for sample (default: "$rootDir/results/TEST/Haplotyper/*somatic*_bqsr.bam")  
    --tumor_bam           = Tumor BAM file (default: "$rootDir/results/TEST/Haplotyper/  *germline_bqsr.bam")  
  
    --normal_bai        = Bam index for normal sample (default: "$rootDir/results/TEST/Haplotyper/*somatic*_bqsr.bam.bai")  
    --tumor_bai         = Bam index for tumor sample (default: "$rootDir/results/TEST/Haplotyper/*germline_bqsr.bam.bai")  
    --outdir            = Output dir (default: "$rootDir/results/TEST/VariantCalling/")  

## Outputs
    Manta_${tumor_sample_id}_vs_${normal_sample_id}.candidateSmallIndels.vcf.gz
    Manta_${tumor_sample_id}_vs_${normal_sample_id}.candidateSV.vcf.gz
    Manta_${tumor_sample_id}_vs_${normal_sample_id}.diploidSV.vcf.gz
    Manta_${tumor_sample_id}_vs_${normal_sample_id}.somaticSV.vcf.gz


## Errs

[2022-03-15T10:35:51.903040Z] [ea3e2b3cd77a] [15_1] [TaskManager] [ERROR] Failed to complete command task: 'getAlignmentStats_generateStats_001' launched from master workflow, error code: 1, command: '/manta/libexec/GetAlignmentStats --ref test.fa --output-file Manta/workspace/alignmentStats.xml.tmpdir/alignmentStats.xml.001.xml --align-file patient_3_germline_bqsr.bam'

2022-03-15T10:36:07.666831Z] [ea3e2b3cd77a] [15_1] [WorkflowRunner] [ERROR] [2022-03-15T10:35:51.812844Z] [ea3e2b3cd77a] [15_1] [getAlignmentStats_generateStats_001]  At least 100 high-confidence read pairs are required to determine pair orientation.

Sol1 
I have found the reason why I failed running some samples, and shared to anyone who might meet the same error issue. Please check the fastq quatlity format, if you use phred64 fastq as input and mapped to the reference, it might be failed, transferred to phred33 by seqtk will solve the problem（https://github.com/lh3/seqtk）

Sol2
