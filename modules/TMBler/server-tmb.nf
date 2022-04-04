nextflow.enable.dsl=2
log.info """
nextflow -log logs/tmbler.log run TMBler/server-tmb.nf -c TMBler/server-tmb.config
"""

process TMBler {
    label "tmbler_module"
    publishDir "${outdir}", mode: params.publish_dir_mode

    input:
    tuple val(port), val(vcfs)
    val genome
    val outdir
    
    output:
    stdout
    path "*.tsv", emit: tmbout

    script:
    """
    Rscript /scripts/tmbler.R --args ${vcfs} ${genome}
    """
}

workflow {
    Channel
        .fromPath( params.vcfs )
        .map( file -> "list('${file.baseName}', '${file}')" )
        .collect()
        .map( file -> "list(${file.join(', ')})")
        .set{ ch_vcfs }
    
    Channel
        .of( params.genome )
        .set{ ch_genome }

    Channel
        .value( params.port )
        .set{ ch_port }

    ch_port.combine(ch_vcfs).set{ ch_group_port }

    TMBler( ch_group_port, ch_genome, params.outdir )
}