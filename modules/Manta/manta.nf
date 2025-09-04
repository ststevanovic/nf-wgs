nextflow.enable.dsl=2

log.info"""
nextflow -log logs/manta.log run modules/Manta/manta.nf -c modules/Manta/manta.config
"""

process MANTA {
    label "powerup"

    tag "${sample_id}"


    input:
        tuple val(sample_id), file(bam), file(bai)
        file(fasta) 
        file(fasta_fai) 

    output:
        path "*"


    script:
    prefix="${sample_id}"
    """
    configManta.py \
        --bam ${bam} \
        --reference ${fasta} \
        --runDir ${prefix}
    python ${prefix}/runWorkflow.py -m local 
    """
}

process MANTAGERMLINE {
    label "powerup"

    tag "${sample_id}_vs_${tumor_sample_id}"

    input:
        tuple val(sample_id), file(bam), file(bai)
        tuple val(tumor_sample_id), file(tumor_bam), file(tumor_bai) 
        file(fasta) 
        file(fasta_fai) 

    output:
        path "*", emit: manta_vcf

    script:

    prefix="${sample_id}_vs_${tumor_sample_id}"
    """
    configManta.py \
        --bam ${bam} \
        --tumorBam ${tumor_bam} \
        --reference ${fasta} \
        --runDir ${prefix}
    python ${prefix}/runWorkflow.py -m local 
    """
}

workflow {
    Channel.fromPath( params.fasta ).set{ ch_fasta }
    Channel.fromPath( params.fai ).set{ ch_fai }
    
    Channel
        .fromFilePairs(params.bambai)
        .map{ it -> [it[0], [it[1]]].flatten() }
        .set{ch_sample}
        
    MANTA( ch_sample, ch_fasta.collect(), ch_fai.collect() )

    if (params.with_tumor) {
        Channel.fromFilePairs(params.bambai_tumor).set{ch_tumor_sample}
        // TODO: found no workflow with this name
        manta_germline_workflow( ch_sample, ch_tumor_sample, ch_fasta, ch_fai )
    }

}