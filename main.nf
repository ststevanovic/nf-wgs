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

include { FASTQC } from './modules/FastQC/fastqc.nf'
include { BWAINDEX; BWAMEM_SAMTOOLS_SORT } from './modules/BWAmem/bwa.nf'
// include { MARKDUPLICATES } from './modules/BWAmem/bwa.nf'
// include { BASERECALIBRATOR } from './modules/BWAmem/bwa.nf'
// include { SAMTOOLS_REPORT } from 'modules/Samtools/samtools.nf'
// include { QUALIMAP } from 'modules/Qualimap/qualimap.nf'
// include { HAPOLOTYPER } from 'modules/Haplotyper/haplotyper.nf'
include { MULTIQC } from './modules/MultiQC/MultiQC.nf'

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
    //     MULTIQC(
    //     QUANT.out.mix(FASTQC.out).collect(),
    //     file(params.multiqc) )

}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : null)
}

workflow.onError {
    log.info "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}