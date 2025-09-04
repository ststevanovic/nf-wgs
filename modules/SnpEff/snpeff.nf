nextflow.enable.dsl=2

log.info """
    nextflow -log logs/snpeff.log run modules/SnpEff/snpeff.nf -c modules/SnpEff/snpeff.config
"""



process SNPEFF {
    tag "${sample_id}"

    input:
    tuple val(sample_id), file(vcf) 
    file(data_dir) 
    val snpeff_db

    output:
    tuple file("${sample_id_reduced}.snpEff.genes.txt"), file("${sample_id_reduced}.snpEff.html"), file("${sample_id_reduced}.snpEff.csv"), emit: snpeff_report
    tuple val(sample_id_reduced), file("${sample_id_reduced}.snpEff.ann.vcf"), emit: snpeff_vcf
    
    // Note: anything from Haplotyper ?
    script:
    sample_id_reduced = sample_id.minus("_bqsr")
    cache = (params.snpeff_cache) ? "-dataDir \${PWD}/${data_dir}" : ""
    """
    snpEff -Xmx8g \
        ${snpeff_db} \
        -csvStats ${sample_id_reduced}.snpEff.csv \
        -nodownload \
        ${cache} \
        -canon \
        -v \
        ${vcf} \
        > ${sample_id_reduced}.snpEff.ann.vcf
    mv snpEff_summary.html ${sample_id_reduced}.snpEff.html
    """
}

process SNPEFF_VERSION {
    tag "SnpEff version"
    publishDir $params.outdir 

    output:
    path "snpeff_version.yml", emit: snpeff_version

    script:
    """
    snpEff -version &> snpeff_version.yml 2>&1 || true
    """
}

workflow {
    ch_snpeff_cache = params.snpeff_cache ? Channel.value(file(params.snpeff_cache)) : "null"
    ch_snpeff_db = params.snpeff_db ? Channel.value(params.snpeff_db) : "null"

    Channel 
        .fromPath(params.vcf)
        .map{file -> tuple(file.baseName, file)}
        .set{ ch_vcf }

 
    
    SNPEFF(
        ch_vcf,
        ch_snpeff_cache,
        ch_snpeff_db
    )
}

// as now have the list of VCFs to annotate, the first step is to annotate with allele frequencies, if there are any
// (vcfSnpeff, vcfVep) = vcfAnnotation.into(2)

