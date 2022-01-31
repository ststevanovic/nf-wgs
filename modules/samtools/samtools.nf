/*
 * Process : Create a FASTA genome index with samtools
 * nextflow -log logs/samtools.log run modules/samtools/samtools.nf -entry samtools_fai -c modules/samtools/samtools.config
 */
nextflow.enable.dsl=2

process SAMTOOLS_FAI { 
  tag { "samtools FAI on ${genome.baseName}" }

  input:
  file genome

  output:
  //  genome_index_ch 
  path "${genome.baseName}.fai", emit: genome_index_fai

  script:
  """
  samtools faidx ${genome}
  """
}

process SAMTOOLS_SORT_BAM {
  input:
  
  output:

  script:
  """
  samtools sort ${input}.sam -o ${input}_sorted.bam
  """
}

process SAMTOOLS_REPORT { 
  input:
  path genome

  output:
  // into genome_index_ch 
  path "${genome}.fai" 

  script:
  """
  samtools faidx ${genome} 
  """
}

workflow samtools_fa_index {
  
  if (params.genome)  {
    Channel
      .fromPath( params.genome )
      .view( "Fasta genome : -> ${it}" )
      .set{ fasta_ch }
  }

  SAMTOOLS_FAI( fasta_ch )

}