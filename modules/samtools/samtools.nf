/*
 * Process : Create a FASTA genome index with samtools
 */

process SAMTOOLS { 
  container 'samtools:latest'

  input:
    path genome from BLANK 

  output:
    path "${genome}.fai" into genome_index_ch 

  script:
  """
  samtools faidx ${genome} 
  """
}