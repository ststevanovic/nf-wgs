process BWA {
  tag { genome.baseName }

  input:
  file genome

  output:
  path "bwa_index", emit: bwa_index

  """
  [ ! -d "bwa_index" ] && mkdir bwa_index ]
  bwa \\
      index \\
      -a bwtsw \\
      -p bwa_index/${genome.baseName} \\
      ${genome} 
  """
}

// index files called from bwa_index directory
// bwa_index = bwa_index_amb.merge(bwa_index_ann, bwa_index_bwt, bwa_index_pac, bwa_index_sa)
// bwa = reads_bwa.combine(bwa_index)

process BWA_MEM {
    tag "$reads"

    // bwa_index

    input:
    set val(name), file(params.reads), file(fasta), file(amb), file(ann), file(bwt), file(pac), file(sa)

    output:
    set val(name), file("${name}.sam")
    
    emit: sam

    """
    bwa mem -M -R '@RG\\tID:${name}\\tSM:${name}\\tPL:Illumina' $fasta $reads > ${name}.sam
    """
}