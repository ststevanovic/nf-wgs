/*
================================================================================
                                     MultiQC
================================================================================
*/

process MULTIQC {
    input:
    // path '*' 
    // path '*' from fastqc_ch.collect()
    file ('FastQC/*') // from fastqc_ch.collect().ifEmpty([])

    // file ('QualimapResults/*') from qualimapResults.collect()
        // file ('bamQC/*') from bamQCReport.collect().ifEmpty([])
        // file ('TrimmedFastQC/*') from trimGaloreReport.collect().ifEmpty([])
        // file ('MarkDuplicates/*') from duplicates_marked_report.collect().ifEmpty([])
//         file ('DuplicatesMarked/*.recal.table') from baseRecalibratorReport.collect().ifEmpty([])
        // file ('SamToolsStats/*') from samtoolsStatsReport.collect().ifEmpty([])
//         file ('snpEff/*') from snpeffReport.collect().ifEmpty([])
//         file ('VCFTools/*') from vcftoolsReport.collect().ifEmpty([])
    
    output:
    path "*multiqc_report.html"                ,emit: report
    path "*_data"                              ,emit: data
    path "*_plots"              ,optional:true ,emit: plots
    // path "versions.yml"                        ,emit: versions

    shell:
    '''
    multiqc -v .
    '''
}

// process MULTIQC_VERSION {
//     // if params version specified, run this

//     shell:
//     '''
//     multiqc --version
//     '''
// }