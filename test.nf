#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
                        Testing, other notes
========================================================================================
 #### Homepage / Documentation
 https://github.com/lifebit-ai/fred-hutch-gatk
----------------------------------------------------------------------------------------
 cmd : nextflow run test.nf -entry desired-workflow1 

 cmdASK: nextflow run path/to/module.nf

*/

include { BWA; BWAMEM; BWAMEM_SAMTOOLS_SORT } from './modules/bwamem/bwa.nf'
include { MARKDUPLICATES; BASERECALIBRATOR; APPLYBQSR } from './modules/bwamem/bwa.nf'


// TODO/ASK: if path/to/module.nf nextflow config inside a module ...
workflow bwa_only {
    Channel
        .fromPath(params.reference_genome)
        .set { fasta_ch }

    Channel
        .fromFilePairs(params.reads)
        .set { reads_ch }
    
    BWA ( fasta_ch )    
    BWAMEM ( reads_ch, fasta_ch, BWA.out.bwa_index )
}

workflow bwa_sam_combined_test {
    Channel
        .fromPath(params.reference_genome)
        .set { fasta_ch }

    fasta_ch.view()

    Channel
        .fromFilePairs(params.reads)
        .set { reads_ch }
    
    BWA ( fasta_ch )    
    BWAMEM_SAMTOOLS_SORT ( reads_ch, fasta_ch, BWA.out.bwa_index )

}

workflow bwa_gatk {
    // # 1
    MARKDUPLICATES(name) // from path (bam_sort) file
    // # 2
    Channel.fromPath(params.fasta)
            .ifEmpty { exit 1, "fasta annotation file not found: ${params.fasta}" }
            .set { fasta_baserecalibrator_ch}
    // # 3
    fasta_baserecalibrator_ch
        .mix( fai_baserecalibrator, dict_baserecalibrator, dbsnp, dbsnp_idx, golden_indel, golden_indel_idx )
        .set{ baserecalibrator_index_ch }
    // # 4
    bam_markdup_baserecalibrator_ch
        .combine( baserecalibrator_index_ch )
        .set{ baserecalibrator_ch }
    // # 5
    BASERECALIBRATOR( baserecalibrator_ch )
    // # 6
    baserecalibrator_table_ch
        .join( bam_markdup_applybqsr_ch )
        .set{ applybqsr_ch }
    // # 6
    APPLYBQSR( applybqsr_ch )
}
