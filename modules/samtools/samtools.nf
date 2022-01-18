/*
 * Process : Create a FASTA genome index with samtools
 */

process SAMTOOLS_FAIDX { 
  
  input:
  
  file genome 

  output:
  //  genome_index_ch 
  path "${genome}.fai", emit genome_index_fai

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
  path genome from BLANK 

  output:
  path "${genome}.fai" into genome_index_ch 

  script:
  """
  samtools faidx ${genome} 
  """
}

// TODO: named 
// workflow preprocess {
//   take:
//     params.genome

//   main:
//     SAMTOOLS_PREPARE_GENOME()
  // output:
    // emit: genome_index_ch

// }