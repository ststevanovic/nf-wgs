// TODO: define
bwa_index = bwa_index_amb.merge(bwa_index_ann, bwa_index_bwt, bwa_index_pac, bwa_index_sa)
bwa = reads_bwa.combine(bwa_index)

process BWA {
    tag "$reads"

    input:
    set val(name), file(reads), file(fasta), file(amb), file(ann), file(bwt), file(pac), file(sa) from bwa

    output:
    set val(name), file("${name}.sam") into sam

    """
    bwa mem -M -R '@RG\\tID:${name}\\tSM:${name}\\tPL:Illumina' $fasta $reads > ${name}.sam
    """
}

// TODO: into params

process MARKDUPLICATES {
  tag "$bam_sort"

  input:
  set val(name), file(bam_sort) from bam_sort

  output:
  set val(name), file("${name}_MarkDup.bam") into bam_markdup_baserecalibrator, bam_markdup_applybqsr
  file "metrics.txt" into markdup_multiqc

  """
  gatk MarkDuplicates -I $bam_sort -M metrics.txt -O ${name}_MarkDup.bam
  """

}

baserecalibrator_index = fasta_baserecalibrator.merge(fai_baserecalibrator, dict_baserecalibrator, dbsnp, dbsnp_idx, golden_indel, golden_indel_idx)
baserecalibrator = bam_markdup_baserecalibrator.combine(baserecalibrator_index)

process BASERECALIBRATOR {
  tag "$bam_markdup"

  input:
  set val(name), file(bam_markdup), file(fasta), file(fai), file(dict), file(dbsnp), file(dbsnp_idx), file(golden_indel), file(golden_indel_idx) from baserecalibrator

  output:
  set val(name), file("${name}_recal_data.table") into baserecalibrator_table
  file ("*data.table") into baseRecalibratorReport

  """
  gatk BaseRecalibrator \
  -I $bam_markdup \
  --known-sites $dbsnp \
  --known-sites $golden_indel \
  -O ${name}_recal_data.table \
  -R $fasta
  """
}