#!/usr/bin/env nextflow
nextflow.enable.dsl=2

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "                     USAGE                              "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "--------------------  QC  ------------------------------"
    log.info ""
    log.info "nextflow -log logs/qualimap.log run modules/Qualimap/qualimap.nf -c modules/Qualimap/qualimap.config"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--inpdir               FOLDER               Folder containing bam files"
    log.info ""
    log.info "Optional arguments:"
    log.info "--outdir               PATH                 Output directory for html and zip files (default='results/TEST/Qualimap')"
    log.info "--output_format        STRING               Output format for individual qualimap files, html or pdf (default=html)"
    log.info ""
    log.info "Flags:"
    log.info "--help                                      Display this message"
    log.info ""
    exit 0
} else {
  log.info "Input           = ${params.inpdir}"
  log.info "Cpus            = ${params.cpu}"
  log.info "Memory          = ${params.mem}"
  log.info "Output          = ${params.outdir}"
  log.info "Format          = ${params.output_format}"
  log.info "help            = ${params.help}"
}

assert (params.inpdir != null) : "please provide the --inpdir option" 

process QUALIMAP {
    label "qualimap_run_settings"
    tag "$sample_id"

    input:
    tuple val(sample_id), file(bam)

    output:
    path ("${sample_id}"), emit: qualimap_results

    shell:
    mem_qm  = params.mem -2 //params.mem.intdiv(2)
    '''

    qualimap bamqc -nt !{params.cpu} --skip-duplicated -bam !{bam} --java-mem-size=!{mem_qm}G -outdir !{bam_tag} -outformat !{params.output_format}
    '''
}

workflow {

  Channel.fromPath( params.inpdir+'/*.bam' )
       .ifEmpty { error "Cannot find any bam file in: ${params.inpdir}" }
       .map{ file -> tuple(file.baseName, file) } 
       .groupTuple(by: 0)
       .set{ bam_init }

  QUALIMAP( bam_init )
}