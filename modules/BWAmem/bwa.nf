nextflow.enable.dsl=2

log.info """
  nextflow -log logs/bwa.log run modules/BWAmem/bwa.nf -c modules/BWAmem/bwa.config
"""

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
    tuple val(name), file("${name}.sorted.aligned.bam"), emit: bam
    
    script:
    cn = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ""
    readGroup = "@RG\\tID:${name}\\t${cn}PU:${name}\\tSM:${name}\\tLB:${name}\\tPL:illumina"
    """
    bwa mem -p -v 3 -t 16 -M ${fasta.baseName} ${reads} -R \"${readGroup}\" | samtools sort -o ${name}.sorted.aligned.bam
    """
}


// gatk tools
process MARKDUPLICATES {
  tag "$bam_sort" 

  input:
  tuple val(name), file(bam_sort) 

  output:
  tuple val(name), file("${name}_MarkDup.bam"), emit: bam_markdup

  """
  gatk MarkDuplicates -M marked_dup_metrics.txt -I $bam_sort -O ${name}.markdup.bam
  """

}

process BASERECALIBRATOR {
  tag "$bam_markdup"

  input:
  tuple val(name), file(bam_markdup)
  file(fasta) 
  file(index) 
  file(dbsnp)
  file(dbsnp_idx)

  // file(golden_indel)

  output:
  tuple val(name), file("${name}_recal_data.table"), emit: baserecalibrator_table
  path "*data.table", emit: baseRecalibratorReport // -> multiqc

  """
  gatk BaseRecalibrator \
  -I $bam_markdup \
  --known-sites $dbsnp \
  -O ${name}_recal_data.table \
  -R $fasta
  """
}

process APPLYBQSR {
    tag "$baserecalibrator_table"

    input:
    tuple val(name), file(baserecalibrator_table), file(bam_markdup)
    file(fasta)
    file(index)

    output:
    tuple val(name), file("${name}.bqsr.bam"), emit: bam_bqsr_ch // will be needed for Haplotyper, Manta
    tuple val(name), file("${name}.bqsr.bam.bai"), emit: indexed_bam_bqsr

    script:
    """
    gatk ApplyBQSR -I $bam_markdup -R $fasta -bqsr $baserecalibrator_table -O ${name}.bqsr.bam
    samtools index ${name}_bqsr.bam
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

  process dictionary {
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
 
workflow {
  
  if (params.dbsnp_gz) {
      Channel.fromPath( params.dbsnp_gz )
            .ifEmpty { exit 1, "dbsnp annotation file not found: ${params.dbsnp_gz}" }
            .set { ch_dbsnp }
  }

  if (params.dbsnp_idx_gz) {
      Channel.fromPath( params.dbsnp_idx_gz )
            .ifEmpty { exit 1, "dbsnp_idx_gz annotation file not found: ${params.dbsnp_idx_gz}" }
            .set { ch_dbsnp_idx }
  }

  Channel
    .fromFilePairs( params.reads, checkIfExists: true )
    .set{ ch_reads }

  Channel
    .fromPath( params.genomes_base )
    .ifEmpty { exit 1, "fasta file not found: ${params.genome}" }
    .set{ ch_fasta }

  dict_baserecalibrator = dictionary( ch_fasta )   

  Channel 
          .fromPath( params.fastaFai )
          .mix( dict_baserecalibrator )
          .set { ch_index }
    
 
  BWAINDEX ( ch_fasta )
  BWAMEM_SAMTOOLS_SORT ( ch_reads, ch_fasta.collect(), BWAINDEX.out.bwa_index.collect() )
  
  MARKDUPLICATES( BWAMEM_SAMTOOLS_SORT.out.bam ) 
  ch_bam_dedup = MARKDUPLICATES.out.bam_markdup

  BASERECALIBRATOR( 
    ch_bam_dedup,
    ch_fasta.collect(),
    ch_index.collect(),
    ch_dbsnp.collect(),
    ch_dbsnp_idx.collect(),
  )
  
  BASERECALIBRATOR.out.baserecalibrator_table
    .join( ch_bam_dedup )
    .set{ ch_table_bam }

  Channel
    .of( BASERECALIBRATOR.out.baserecalibrator_table )

  APPLYBQSR( 
    ch_table_bam,
    ch_fasta.collect(), 
    ch_index.collect() 
    )
}
