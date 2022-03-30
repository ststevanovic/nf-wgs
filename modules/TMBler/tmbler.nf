nextflow.enable.dsl=2
log.info """
nextflow -log logs/tmbler.log run modules/TMBler/tmbler.nf -c modules/TMBler/tmbler.config
"""

process TMB {
    label "tmbler_module"
    containerOptions "--name tmbler \
                  --rm \
                  -e USER=\$(whoami) \
                  -e PASSWORD=helloworld \
                  -e USERID=\$UID \
                  -p 8787:${port}"
    
    input:
    tuple val(port), file(vcf)
    
    output:
    stdout

    script:
    
    genome = params.genome
    """
    #!/usr/bin/env Rscript

    genome_assembly <- "${genome}"
    vcf_file <- "${vcf}"

    print(vcf_file)    
    print("Done.")
    """
}
    // vcfs <- readVcfFiles(vcfFiles = vcf_file, assembly = genome_assembly)
//     TMB_res <- applyTMB(inputForTMB = vcfs_all, assembly = genome_assembly)
//     DT::datatable(TMB_res)

process test_inpt {
    label "tmbler_module"

    input:
    tuple val(port), file(vcf)

    output:
    stdout

    script:
    """
    #!/usr/bin/env Rscript

    print("${port}")
    print("${vcf}")
    """

}

workflow {
    Channel.fromPath(params.vcfs).set{ch_vcfs}

        ch_vcfs
        .collect()
        .withIndex()
        .flatMap()
        .map {
            it -> [it[1]+8787, it[0] ] 
        }
        .set{ ch_vcfs_ports }
    
    // ch_vcfs_ports.view()
    TMB( ch_vcfs_ports )

    // test_inpt(ch_vcfs_ports)
    // test_inpt.out.view()

    // TMB.out.view()
}

