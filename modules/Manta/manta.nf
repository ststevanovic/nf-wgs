nextflow.enable.dsl=2

process MantaSingle {
    // label 'cpus_max'
    // label 'memory_max'
    tag "${sample_id}"

    input:
        set patient_id, sample_id, file(bam), file(bai)  
        file(fasta) 
        file(fasta_fai) 
        file(target_bed) 

    output:
        set val("Manta"), patient_id, sample_id, file("*.vcf.gz"), file("*.vcf.gz.tbi"), emit: manta_single_vcf
    
    script:
    beforeScript = params.target_bed ? "bgzip --threads ${task.cpus} -c ${target_bed} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
    options = params.target_bed ? "--exome --callRegions call_targets.bed.gz" : ""
    status = statusMap[idPatient, sample_id]
    inputbam = status == 0 ? "--bam" : "--tumorBam"
    vcftype = status == 0 ? "diploid" : "tumor"
    """
    ${beforeScript}
    configManta.py \
        ${inputbam} ${bam} \
        --reference ${fasta} \
        ${options} \
        --runDir Manta
    python Manta/runWorkflow.py -m local -j ${task.cpus}
    mv Manta/results/variants/candidateSmallIndels.vcf.gz \
        Manta_${sample_id}.candidateSmallIndels.vcf.gz
    mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
        Manta_${sample_id}.candidateSmallIndels.vcf.gz.tbi
    mv Manta/results/variants/candidateSV.vcf.gz \
        Manta_${sample_id}.candidateSV.vcf.gz
    mv Manta/results/variants/candidateSV.vcf.gz.tbi \
        Manta_${sample_id}.candidateSV.vcf.gz.tbi
    mv Manta/results/variants/${vcftype}SV.vcf.gz \
        Manta_${sample_id}.${vcftype}SV.vcf.gz
    mv Manta/results/variants/${vcftype}SV.vcf.gz.tbi \
        Manta_${sample_id}.${vcftype}SV.vcf.gz.tbi
    """
}

manta_single_vcf = manta_single_vcf.dump(tag:'Single Manta')

// Kako se ovo interpretira ? 
// pairBam (pairBamManta, pairBamStrelka, pairBamStrelkaBP, pairBamMsisensor, pairBamCNVkit, pairBam) = pairBam.into(6)

process Manta {
    // label 'cpus_max'
    // label 'memory_max'

    tag "${tumor_sample_id}_vs_${normal_sample_id}"


    input:
        set patient_id, normal_sample_id, file(normal_bam), file(normal_bai), tumor_sample_id, file(tumor_bam), file(tumor_bai) // from ch_pairBamManta
        file(fasta) 
        file(fasta_fai) 
        file(target_bed) 

    output:
        set val("Manta"), patient_id, val("${tumor_sample_id}_vs_${normal_sample_id}"), file("*.vcf.gz"), file("*.vcf.gz.tbi"), emit: manta_vcf
        set patient_id, normal_sample_id, tumor_sample_id, file("*.candidateSmallIndels.vcf.gz"), file("*.candidateSmallIndels.vcf.gz.tbi"), emit: manta_to_strelka


    script:
    beforeScript = params.target_bed ? "bgzip --threads ${task.cpus} -c ${target_bed} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
    options = params.target_bed ? "--exome --callRegions call_targets.bed.gz" : ""
    """
    ${beforeScript}
    configManta.py \
        --normalBam ${normal_bam} \
        --tumorBam ${tumor_bam} \
        --reference ${fasta} \
        ${options} \
        --runDir Manta
    python Manta/runWorkflow.py -m local -j ${task.cpus}
    mv Manta/results/variants/candidateSmallIndels.vcf.gz \
        Manta_${tumor_sample_id}_vs_${normal_sample_id}.candidateSmallIndels.vcf.gz
    mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
        Manta_${tumor_sample_id}_vs_${normal_sample_id}.candidateSmallIndels.vcf.gz.tbi
    mv Manta/results/variants/candidateSV.vcf.gz \
        Manta_${tumor_sample_id}_vs_${normal_sample_id}.candidateSV.vcf.gz
    mv Manta/results/variants/candidateSV.vcf.gz.tbi \
        Manta_${tumor_sample_id}_vs_${normal_sample_id}.candidateSV.vcf.gz.tbi
    mv Manta/results/variants/diploidSV.vcf.gz \
        Manta_${tumor_sample_id}_vs_${normal_sample_id}.diploidSV.vcf.gz
    mv Manta/results/variants/diploidSV.vcf.gz.tbi \
        Manta_${tumor_sample_id}_vs_${normal_sample_id}.diploidSV.vcf.gz.tbi
    mv Manta/results/variants/somaticSV.vcf.gz \
        Manta_${tumor_sample_id}_vs_${normal_sample_id}.somaticSV.vcf.gz
    mv Manta/results/variants/somaticSV.vcf.gz.tbi \
        Manta_${tumor_sample_id}_vs_${normal_sample_id}.somaticSV.vcf.gz.tbi
    """
}

manta_vcf = manta_vcf.dump(tag:'Manta')

workflow single_workflow {
    // call ch here
    take:
        ch_singleBam
        ch_fasta
        ch_fai
        ch_target_bed
}

workflow paired_workflow {
    // call ch here
    take:
        ch_fasta
        ch_fai
        ch_target_bed
    main:
        Manta()
    emit:
    // where my_data goes ? 
        my_data = process.out
}
// 
workflow {
    // define ch here
    Channel.fromPath( params.fasta ).set{ ch_fasta }
    Channel.fromPath( params.fai ).set{ ch_fai }
    Channel.fromPath( params.vcf ).set{ ch_target_bed }
    // ch_singleBam
    // ch_fasta
    // ch_fai
    // ch_tar get_bed
   
}