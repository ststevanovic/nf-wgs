nextflow.enable.dsl=2

process BWA {
  // ASK: output dir modular selection (env.config)
  tag { genome.baseName }

  input:
  file genome

  output:
  path "${genome.baseName}.{amb,ann,bwt,pac,sa}", emit: bwa_index

//   [ ! -d "bwa_index" ] && mkdir bwa_index ]
  """
  bwa \\
      index \\
      -a bwtsw \\
      -p ${genome.baseName} \\
      ${genome} 
  """
}

process BWAMEM {
    tag "${name}"
    
    
    input:
    tuple val(name), path(reads)
    path fasta
    path index

    
    output:
    tuple val(name), file("${name}.aligned.bam")

    script:

    """
    bwa mem -p -v 3 -t 16 -M ${fasta.baseName} ${reads} >> ${name}.aligned.bam
    """
}

process BWAMEM_SAMTOOLS_SORT {
    tag "${name}"

    input:
    tuple val(name), path(reads)
    path fasta
    path index

    output:
    tuple val(name), file("${name}_sorted.aligned.bam")

    script:

    """
    bwa mem -p -v 3 -t 16 -M ${fasta.baseName} ${reads} | samtools sort -o ${name}_sorted.aligned.bam
    """
}


// gatk tools
process MARKDUPLICATES {
  tag "$name"
  tag "$bam_sort" // ??

  input:
  set val(name), file(bam_sort) from bam_sort

  output:
  set val(name), file("${name}_MarkDup.bam"), emit: bam_markdup_baserecalibrator_ch, bam_markdup_applybqsr_ch
  file "metrics.txt", emit: markdup_multiqc_ch

  """
  gatk MarkDuplicates -I $bam_sort -M metrics.txt -O ${name}_MarkDup.bam
  """

}

process BASERECALIBRATOR {
  tag "$bam_markdup"

  input:
  set val(name), file(bam_markdup), file(fasta), file(fai), file(dict), file(dbsnp), file(dbsnp_idx), file(golden_indel), file(golden_indel_idx) 

  output:
  set val(name), file("${name}_recal_data.table"), emit: baserecalibrator_table_ch
  file ("*data.table"), emit: baseRecalibratorReport_ch

  """
  gatk BaseRecalibrator \
  -I $bam_markdup \
  --known-sites $dbsnp \
  --known-sites $golden_indel \
  -O ${name}_recal_data.table \
  -R $fasta
  """
}

process APPLYBQSR {
    tag "$baserecalibrator_table"

    input:
    set val(name), file(baserecalibrator_table), file(bam_markdup) 

    output:
    set val(name), file("${name}_bqsr.bam"), emit: bam_bqsr_ch // will be needed for Haplotyper, (IndexBam?)

    script:
    """
    gatk ApplyBQSR -I $bam_markdup -bqsr $baserecalibrator_table -O ${name}_bqsr.bam
    """
}

workflow {

  // 
  // ASK: Validate inputs, this might be handy, meaning
  //
  // params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
  // if (params.fasta) {
  //   Channel.fromPath(params.fasta)
  //         .ifEmpty { exit 1, "fasta annotation file not found: ${params.fasta}" }
  //         .into { fasta_scatter_intervals; fasta_bwa; fasta_baserecalibrator; fasta_haplotypecaller; fasta_genotypeVcfs; fasta_variantrecalibrator_snps; fasta_variantrecalibrator_tranches; fasta_variant_eval; fasta_structural_variantcaller }
  //   }
  // }
  // ************
  
  
  // # 1
  MARKDUPLICATES(name) // from path (bam_sort) file
  // # 2
  Channel.fromPath(params.fasta)
        .ifEmpty { exit 1, "fasta annotation file not found: ${params.fasta}" }
        .set { fasta_baserecalibrator_ch}
  // # 3
  fasta_baserecalibrator_ch
    .mix( fai_baserecalibrator, dict_baserecalibrator, dbsnp, dbsnp_idx, golden_indel, golden_indel_idx )
    .set{ baserecalibrator_index_ch }
  // # 4
  bam_markdup_baserecalibrator_ch
    .combine( baserecalibrator_index_ch )
    .set{ baserecalibrator_ch }
  // # 5
  BASERECALIBRATOR( baserecalibrator_ch )
  // # 6
  baserecalibrator_table_ch
    .join( bam_markdup_applybqsr_ch )
    .set{ applybqsr_ch }
  // # 6
  APPLYBQSR( applybqsr_ch )
}
