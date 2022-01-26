nextflow.enable.dsl=2

log.info """
  nextflow -log logs/bwa.log run modules/bwamem/bwa.nf -entry gatk_workflow -c modules/bwamem/bwa.config
"""
  // TODO: check dump, use emit

process BWAINDEX {
  // TODO: ASK: output dir modular selection (env.config)
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
    tuple val(name), file("${name}_sorted.aligned.bam"), emit: bam

    script:

    """
    bwa mem -p -v 3 -t 16 -M ${fasta.baseName} ${reads} | samtools sort -o ${name}_sorted.aligned.bam
    """
}


// gatk tools
process MARKDUPLICATES {
  tag "$bam_sort" 

  input:
  tuple val(name), file(bam_sort) 

  output:
  tuple val(name), file("${name}_MarkDup.bam"), emit: bam_markdup
  // TODO: can't emit two var -- is this done by splitChannels ?
  // bam_markdup_baserecalibrator, bam_markdup_applybqsr
  // TODO: fix for multiqc needed
  // tuple val("marked_dup_metrics") file("marked_dup_metrics.txt"), emit: markdup_multiqc

  """
  gatk MarkDuplicates -M marked_dup_metrics.txt -I $bam_sort -O ${name}_MarkDup.bam
  """

}

process BASERECALIBRATOR {
  tag "$bam_markdup"

  input:
  tuple val(name), file(bam_markdup), file(fasta), file(fai), file(dict), file(dbsnp), file(dbsnp_idx), file(golden_indel), file(golden_indel_idx) 

  output:
  tuple val(name), file("${name}_recal_data.table"), emit: baserecalibrator_table
  path "*data.table", emit: baseRecalibratorReport // -> multiqc

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
    tuple val(name), file(baserecalibrator_table), file(bam_markdup) 

    output:
    tuple val(name), file("${name}_bqsr.bam"), emit: bam_bqsr_ch // will be needed for Haplotyper, (IndexBam?)

    script:
    """
    gatk ApplyBQSR -I $bam_markdup -bqsr $baserecalibrator_table -O ${name}_bqsr.bam
    """
}

// util
 process gunzip_dbsnp {
    tag "$dbsnp_gz"

    input:
    file dbsnp_gz 
    file dbsnp_idx_gz 

    output:
    // into dbsnp, dbsnp_variantrecalibrator_snps, dbsnp_variantrecalibrator_indels
    path "*.vcf", emit: dbsnp
    // into dbsnp_idx, dbsnp_idx_variantrecalibrator_snps, dbsnp_idx_variantrecalibrator_indels
    path "*.vcf.idx", emit: dbsnp_idx

    script:
    if ( "${dbsnp_gz}".endsWith(".gz") ) {
     """
     gunzip -d --force $dbsnp_gz
     gunzip -d --force $dbsnp_idx_gz
     """
   } else {
     """
     cp $dbsnp_gz dbsnp.vcf
     cp $dbsnp_idx_gz dbsnp.vcf.idx
     """
   }
  }


  process gunzip_golden_indel {
    tag "$golden_indel_gz"

    input:
    file golden_indel_gz
    file golden_indel_idx_gz

    output:
    // golden_indel, golden_indel_variantrecalibrator_indels
    path "*.vcf", emit: golden_indel
    // golden_indel_idx, golden_indel_idx_variantrecalibrator_indels
    path "*.vcf.idx", emit: golden_indel_idx

    script:
    if ( "${golden_indel_gz}".endsWith(".gz") ) {
     """
     gunzip -d --force $golden_indel_gz
     gunzip -d --force $golden_indel_idx_gz
     """
    } else {
     """
     cp $golden_indel_gz golden_indel.vcf
     cp $golden_indel_idx_gz golden_indel.vcf.idx
     """
   }
  }

  process DICTIONARY {
    input:
    file genome  // .fasta

    output:
    path "${genome.baseName}.dict", emit: dict

     """
    gatk CreateSequenceDictionary \
    R=${genome} \
    O=${genome.baseName}.dict
    """
  }
 
workflow gatk_workflow {
  // if (params.bwasam) {
  // Channel
  //     .fromPath(params.bwasam)
  //     .map { file -> tuple(file.baseName, file) }
  //     .ifEmpty { exit 1, "Cannot find any reads matching: ${params.bwasam}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
  //     // .set { reads_fastqc; reads_files }
  //     .set { reads_fastqc }
      
  //     // .dump(tag:'input')
  // }  
  // reads_fastqc.view()

  
  // TODO: TEMP:, from emit
  
  if (params.bwasam) {
  Channel
    .fromPath( params.bwasam )
    .map { file -> tuple(file.baseName, file) }
    .ifEmpty { exit 1, "Cannot find any bams matching: ${params.bwasam}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .view{ "BAMfile: $it" }
    .set{ bam_sort_ch }
  }  

  // // //
  
  // # 2 - channels
  if (params.genome) {
  Channel.fromPath( params.genomes_base )
        .ifEmpty { exit 1, "fasta annotation file not found: ${params.genome}" }
        .set { fasta_baserecalibrator_ch }
  }

  // # 2 - utils
  if (params.dbsnp_gz) {
      Channel.fromPath( params.dbsnp_gz )
            .ifEmpty { exit 1, "dbsnp annotation file not found: ${params.dbsnp_gz}" }
            .set { dbsnp_gz_ch }
  }

  if (params.dbsnp_idx_gz) {
      Channel.fromPath( params.dbsnp_idx_gz )
            .ifEmpty { exit 1, "dbsnp_idx_gz annotation file not found: ${params.dbsnp_idx_gz}" }
            .set { dbsnp_idx_gz_ch }
  }

  if (params.golden_indel_gz) {
      Channel.fromPath( params.golden_indel_gz )
            .ifEmpty { exit 1, "golden_indel_gz annotation file not found: ${params.golden_indel_gz}" }
            .set { golden_indel_gz_ch }
  }
  
  if (params.golden_indel_idx_gz) {
      Channel.fromPath(params.golden_indel_idx_gz)
            .ifEmpty { exit 1, "golden_indel_idx_gz annotation file not found: ${params.golden_indel_idx_gz}" }
            .set { golden_indel_idx_gz_ch }
  }

  Channel
    .fromPath( params.fastaFai )
    .view { "Fasta FAI : $it" }
    .set{ fai_baserecalibrator_ch }


  MARKDUPLICATES( bam_sort_ch ) // emits bam_markdup
  
  // ------ utils -----
  gunzip_dbsnp( dbsnp_gz_ch, dbsnp_idx_gz_ch )                      // snips
  gunzip_golden_indel( golden_indel_gz_ch, golden_indel_idx_gz_ch)  // indels
  dict_baserecalibrator = DICTIONARY( fasta_baserecalibrator_ch )   // reference
  // ------------------
  
  // TODO: no space left on device...
  gunzip_dbsnp.out
  
  // fasta_baserecalibrator_ch
  //   .mix( 
  //     fasta_baserecalibrator_ch,
  //     fai_baserecalibrator_ch,
  //     dict_baserecalibrator, 
  //     util_dbsnp.out.dbsnp, 
  //     util_dbsnp.out.dbsnp_idx, 
  //     util_dbsnp_gi.out.golden_indel, 
  //     util_dbsnp_gi.out.golden_indel_idx )
  //   .set{ baserecalibrator_index_ch }
  
  // MARKDUPLICATES.out.bam_markdup // here !!! format: MARKDUPLICATES .. markduplicates // MARKDUPLICATES.out.out ? naming convention and docs ?
  //   .combine( baserecalibrator_index_ch )
  //   .set{ baserecalibrator_ch }

  // // // # 2-processes
  // BASERECALIBRATOR( baserecalibrator_ch )
  
  // // // # 3
  // BASERECALIBRATOR.out.baserecalibrator_table
  //   .join( MARKDUPLICATES.out.bam_markdup )
  //   .set{ applybqsr_ch }
  
  // APPLYBQSR( applybqsr_ch )
}
