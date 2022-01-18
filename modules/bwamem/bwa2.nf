#!/usr/bin/env nextflow

// *** NOT USED IN THIS PIPELINE... ***

// use dsl2

params.reads = "$baseDir/data/reads/*_ercc_{1,2}.fq.gz"
channel.fromPath(params.reads)
    .println { file -> "$file.name \tdirectory: $file.parent"}
params.outdir = "$baseDir/tests/bwa/expected"

log.info """\
         BWA MEM 2 - N F   P I P E L I N E
         ===================================
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()


process index {
    
    input:
    path reads from params.transcriptome_file
     
    output:
    path 'index' into index_ch

    script:       
    """
    bwa-mem2 index --threads $task.cpus -t $transcriptome -i index
    """
}