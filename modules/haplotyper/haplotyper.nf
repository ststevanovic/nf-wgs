applybqsr = baserecalibrator_table.join(bam_markdup_applybqsr)

process ApplyBQSR {
    tag "$baserecalibrator_table"
    container 'gatk:latest'

    input:
    set val(name), file(baserecalibrator_table), file(bam_markdup) from applybqsr

    output:
    set val(name), file("${name}_bqsr.bam") into bam_bqsr

    script:
    """
    gatk ApplyBQSR -I $bam_markdup -bqsr $baserecalibrator_table -O ${name}_bqsr.bam
    """
}
