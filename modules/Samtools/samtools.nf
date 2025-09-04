nextflow.enable.dsl=2

log.info "===================================================================="
log.info " Indexing Fasta fasta file, sorting BAM file                        "
log.info " nextflow -log logs/samtools.log  \\                                "
log.info " run modules/Samtools/samtools.nf \\                                "
log.info " -entry samtools_fai -c modules/Samtools/samtools.config            "
log.info "===================================================================="


process SAMTOOLS_FAI { 
  tag { "samtools FAI on ${fasta.baseName}" }

  input:
  file fasta

  output:
  path "${fasta.baseName}.fai", emit: fasta_index_fai

  script:
  """
  samtools faidx ${fasta} -o ${fasta.baseName}.fai
  """
}

process SAMTOOLS_SORT_BAM {
  input:
  file sam
  
  output:
  path "${sam}_sorted.bam", emit: sorted_bam

  script:
  """
  samtools sort ${sam}.sam
  """
}

process SAMTOOLS_REPORT { 
  input:
  path bam

  output:
  path "${bam.baseName}_reports.txt" 

  script:
  """
  samtools stats ${bam} > "${bam.baseName}_reports.txt"
  """
}

workflow samtools_fa_index_report {

  if (params.fasta)  {
    Channel
          .fromPath( params.fasta )
          .view { "Fasta : $it" }
          .set{ ch_fasta }
  }

  Channel
    .fromPath( params.bam )
    .view{ "Bam file : $it" }
    .set{ ch_bam }

  SAMTOOLS_FAI( ch_fasta )

  SAMTOOLS_REPORT( ch_bam )

}