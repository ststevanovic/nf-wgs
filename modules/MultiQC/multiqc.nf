/*
================================================================================
                                     MultiQC
================================================================================
*/

process MULTIQC {
    // publishDir "${params.outdir}/Reports/MultiQC", mode: 'copy' should be defined inside config
    // output dir from environemnt.config used...
    container 'multiqc:v1.11'

    input:
    path files_path
    // file ("FastQC/*") 
    // file ('GATK/*') 
    // file ('Haplotyper/*') 
    // file ('Manta/*') 
    // file ('Qualimap/*') 
    // file ('GATK/*.recal_data.table') 
    // file ('SNPeff/*') 
    // file ('VEP/*') 
    
    output:
    path "*multiqc_report.html"                ,emit: report
    path "*_data"                              ,emit: data
    path "*_plots"              ,optional:true ,emit: plots
    path "versions.yml"                        ,emit: versions

    shell:
    '''
    multiqc -v -f .
    '''
}

process MULTIQC_VERSION {
    shell:
    '''
    multiqc --version
    '''
}