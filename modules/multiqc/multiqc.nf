process MULTIQC {
    container 'multiqc:latest' // ???
    publishDir params.outdir, mode:'copy'
    errorStrategy params.errStrategy // 'ignore'
    cpus params.cpu
    memory params.mem+'G'
   
    input:
    // apply mix
    path '*' from fastqc_ch.collect()
    file qualimap_results from qualimap_results.collect()
    
    // path or file ?
    output:
    // results/multiqc/multiqc_report.html
    // path 'multiqc_report.html'
    file("multiqc_report.html") into params.outdir
    file("multiqc_data/") into final_output_data

    shell:
    config = config_file.name != 'NO_FILE' ? "--config $config_file" : ''
    '''
    multiqc !{config} -v .
    '''
}

/*
================================================================================
                                     MultiQC
================================================================================
*/

// STEP MULTIQC

// process MultiQC {
//     publishDir "${params.outdir}/Reports/MultiQC", mode: params.publish_dir_mode

//     input:
//         file (multiqcConfig) from ch_multiqc_config
//         file (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
//         file (versions) from ch_software_versions_yaml.collect()
//         file workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")
//         file ('bamQC/*') from bamQCReport.collect().ifEmpty([])
//         file ('BCFToolsStats/*') from bcftoolsReport.collect().ifEmpty([])
//         file ('FastQC/*') from fastQCReport.collect().ifEmpty([])
//         file ('TrimmedFastQC/*') from trimGaloreReport.collect().ifEmpty([])
//         file ('MarkDuplicates/*') from duplicates_marked_report.collect().ifEmpty([])
//         file ('DuplicatesMarked/*.recal.table') from baseRecalibratorReport.collect().ifEmpty([])
//         file ('SamToolsStats/*') from samtoolsStatsReport.collect().ifEmpty([])
//         file ('snpEff/*') from snpeffReport.collect().ifEmpty([])
//         file ('VCFTools/*') from vcftoolsReport.collect().ifEmpty([])

//     output:
//         file "*multiqc_report.html" into ch_multiqc_report
//         file "*_data"
//         file "multiqc_plots"

//     when: !('multiqc' in skipQC)

//     script:
//     rtitle = ''
//     rfilename = ''
//     if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
//         rtitle = "--title \"${workflow.runName}\""
//         rfilename = "--filename " + workflow.runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report"
//     }
//     custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
//     """
//     multiqc -f ${rtitle} ${rfilename} ${custom_config_file} .
//     """
// }

// ch_multiqc_report.dump(tag:'MultiQC')

// MULTIQC W CONFIG

// params.multiqc_config = "$baseDir/assets/multiqc_config.yaml"
// Channel.fromPath(params.multiqc_config, checkIfExists: true).set { ch_config_for_multiqc }

// process multiqc {
//     input:
//     file multiqc_config from ch_config_for_multiqc
//     file ('fastqc/*') from fastqc_results.collect().ifEmpty([])

//     output:
//     file "multiqc_report.html" into multiqc_report
//     file "multiqc_data"

//     script:
//     """
//     multiqc --config $multiqc_config .
//     """
// }