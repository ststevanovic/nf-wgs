/*
 * Process : Create a FASTA genome index with samtools
 */

process SAMTOOLS_PREPARE_GENOME { 
  input:
  // ???
  file genome from params.genome

  output:
  path "${genome}.fai" into genome_index_ch 
  //??? bwa_index_amb, bwa_index_ann, bwa_index_bwt, bwa_index_pac, bwa_index_sa

  script:
  """
  samtools faidx ${genome}
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