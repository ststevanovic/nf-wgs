nextflow.enable.dsl=2
log.info"""
nextflow -log logs/vep.log run modules/VEP/vep.nf -c modules/VEP/vep.config
"""
process VEP_VERSION {
    publishDir params.outdir, mode: params.publish_dir_mode, saveAs: {
        "vep_version/${it}"
    }

    output:
    file "vep_version.yaml"

    script:
    """
    vep --help &> v_vep.txt 2>&1 || true
    cat v_vep.txt >> vep_version.yaml
    """
}

process VEP {
    tag "${sample_id}"

    publishDir params.outdir, mode: params.publish_dir_mode, saveAs: {
        if (it == "${sample_id}.VEP.summary.html") "Reports/${it}"
        else it
    }

    input:
        tuple val(sample_id), file(vcf)

    output:
        tuple val(sample_id), file("${sample_id}.VEP.ann.vcf"), emit: vepVCF
        path "${sample_id}.VEP.summary.html", emit: vepReport

    script:
    // TODO: pass in channel..
    genome = params.genome
    """
    mkdir ${sample_id}
    vep \
        -i ${vcf} \
        -o ${sample_id}.VEP.ann.vcf \
        --assembly ${genome} \
        --everything \
        --filter_common \
        --fork ${task.cpus} \
        --format vcf \
        --per_gene \
        --stats_file ${sample_id}.VEP.summary.html \
        --total_length \
        --offline   \
        --dir_cache "/.vep/" \
        --vcf
    rm -rf ${sample_id}
    """
}

workflow {
    Channel 
        .fromPath(params.vcf)
        .map{file -> tuple(file.baseName, file)}
        .set{ ch_vcf }

    VEP_VERSION()
    
    VEP( ch_vcf )
}   

