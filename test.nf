#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
                        Testing, other notes
========================================================================================
 #### Homepage / Documentation
 https://github.com/lifebit-ai/fred-hutch-gatk
----------------------------------------------------------------------------------------
 cmd : nextflow run test.nf -entry desired-workflow1 -C module.config

 cmdASK: nextflow run path/to/module.nf 

*/

// 
// ASK: Validate inputs, this might be handy, meaning
//
// params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
// if (params.fasta) {
//   Channel.fromPath(params.fasta)
//         .ifEmpty { exit 1, "fasta annotation file not found: ${params.fasta}" }
//         .into { fasta_scatter_intervals; fasta_bwa; fasta_baserecalibrator; fasta_haplotypecaller; fasta_genotypeVcfs; fasta_variantrecalibrator_snps; fasta_variantrecalibrator_tranches; fasta_variant_eval; fasta_structural_variantcaller }
//   }
// }
// ************

include { BWAINDEX; BWAMEM; BWAMEM_SAMTOOLS_SORT } from './modules/bwamem/bwa.nf'
include { MARKDUPLICATES; BASERECALIBRATOR; APPLYBQSR } from './modules/bwamem/bwa.nf'


// TODO/ASK: if path/to/module.nf nextflow config inside a module ...
// workflow bwa_only {
//     Channel
//         .fromPath(params.reference_genome)
//         .set { fasta_ch }

//     Channel
//         .fromFilePairs(params.reads)
//         .set { reads_ch }
    
//     BWAINDEX ( fasta_ch )    
//     BWAMEM ( reads_ch, fasta_ch, BWA.out.bwa_index )
// }

// workflow bwa_sam_combined_test {
//     Channel
//         .fromPath(params.reference_genome)
//         .set { fasta_ch }

//     fasta_ch.view()

//     Channel
//         .fromFilePairs(params.reads)
//         .set { reads_ch }
    
//     BWAINDEX ( fasta_ch )    
//     BWAMEM_SAMTOOLS_SORT ( reads_ch, fasta_ch, BWA.out.bwa_index )

// }

// workflow bwa_gatk {
//     // # 1
//     MARKDUPLICATES(name) // from path (bam_sort) file
//     // # 2
//     Channel.fromPath(params.fasta)
//             .ifEmpty { exit 1, "fasta annotation file not found: ${params.fasta}" }
//             .set { fasta_baserecalibrator_ch}
//     // # 3
//     fasta_baserecalibrator_ch
//         .mix( fai_baserecalibrator, dict_baserecalibrator, dbsnp, dbsnp_idx, golden_indel, golden_indel_idx )
//         .set{ baserecalibrator_index_ch }
//     // # 4
//     bam_markdup_baserecalibrator_ch
//         .combine( baserecalibrator_index_ch )
//         .set{ baserecalibrator_ch }
//     // # 5
//     BASERECALIBRATOR( baserecalibrator_ch )
//     // # 6
//     baserecalibrator_table_ch
//         .join( bam_markdup_applybqsr_ch )
//         .set{ applybqsr_ch }
//     // # 6
//     APPLYBQSR( applybqsr_ch )
// }

workflow test_this {    
    
    
    Channel
        .fromPath(params.reference_genome)
        .set { fasta_ch }

    BWAINDEX ( fasta_ch )

    Channel
        .of ( BWAINDEX.out.bwa_index )
        .branch{
            branch_1: it
            branch_2: it
        }.set{ results_ch }

    results_ch.branch_1.view()

    results_ch.branch_2.view()

}