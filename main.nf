nextflow.enable.dsl=2

params.reads = "$baseDir/data/reads/*_ercc_{1,2}.fq.gz"
params.multiqc = "$baseDir/results/multiqc"
params.outdir = "$baseDir/results"

log.info """\
         W G S - N F   P I P E L I N E    
         ===================================
         reads        : ${params.reads}
         multiqc      : ${params.multiqc}
         outdir       : ${params.outdir}
         """
         .stripIndent()

include { FASTQC } from 'modules/fastqc/fastqc.nf'
// include { BWA } from 'modules/bwa/bwa.nf'
// include { SAMTOOLS } from 'modules/samtools/samtools.nf'
// include { QUALIMAP } from 'modules/qualimap/qualimap.nf'
// include { HAPOLOTYPER } from 'modules/haplotyper/haplotyper.nf'
// include { MULTIQC } from 'modules/multiqc/multiqc.nf'

// unnamed 
workflow { 
    channel
    .fromFilePairs( params.reads, checkIfExists: true )
    .set { read_pairs_ch }

    fastqc_ch = FASTQC( read_pairs_ch )

    // bwa x3 processes

    // samtools

    // qualimap

    // haplotyper

    // ...

    // multiqc
}

// // named 
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
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}