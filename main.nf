nextflow.enable.dsl=2

log.info """\
         W G S - N F   P I P E L I N E    
         ===================================
         reads              : ${params.reads}
         multiqc            : ${params.outmultiqc}
         outdir             : ${params.outdir}
         reference_genome   : ${params.reference_genome}
         """
         .stripIndent()

include { FASTQC } from './modules/fastqc/fastqc.nf'
include { BWAINDEX; BWAMEM_SAMTOOLS_SORT } from './modules/bwamem/bwa.nf'
// include { MARKDUPLICATES } from './modules/bwamem/bwa.nf'
// include { BASERECALIBRATOR } from './modules/bwamem/bwa.nf'
// include { SAMTOOLS_REPORT } from 'modules/samtools/samtools.nf'
// include { QUALIMAP } from 'modules/qualimap/qualimap.nf'
// include { HAPOLOTYPER } from 'modules/haplotyper/haplotyper.nf'
include { MULTIQC } from './modules/multiqc/multiqc.nf'

// unnamed 
workflow { 
    // channels
    Channel
        .fromFilePairs( params.reads, checkIfExists: true )
        .set { reads_ch }

    Channel
        .fromPath( params.reference_genome )
        .set{ fasta_ch }

    // reads
    FASTQC ( reads_ch )

    // TODO: FASTQC.out.fastqc where to use this ?

    // bwa index (preprocess)
    BWAINDEX ( fasta_ch )
    // bwa mem & samtools for sorting
    bwamem_sam_ch = BWAMEM_SAMTOOLS_SORT ( reads_ch, fasta_ch, BWAINDEX.out.bwa_index )
    // 

    // TODO: qualimap

    // TODO: haplotyper

    // TODO: ...

    MULTIQC( FASTQC.out.fastqc.collect().ifEmpty([]) )
}

// TODO: named 
// workflow EXAMPLE {
//   take:
//     transcriptome
//     read_pairs_ch

//   main:
//     INDEX( transcriptome )
//     QUANT( INDEX.out, read_pairs_ch )
//     FASTQC( read_pairs_ch )
//     MULTIQC(
//             QUANT.out.mix(FASTQC.out).collect(),
//             file(params.multiqc) )
// }

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : null)
}

workflow.onError {
    log.info "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}