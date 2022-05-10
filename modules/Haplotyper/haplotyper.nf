nextflow.enable.dsl=2

log.info """
  nextflow -log logs/hap.log run modules/Haplotyper/haplotyper.nf -entry haplotyper_module -c modules/Haplotyper/haplotyper.config
"""

process HAPLOTYPCALLER {
  tag "$sample"

  input:
  tuple val(sample), file(bam_bqsr), file(bai)
  file(fasta)
  file(fai)
  file(dict)

  output:
  path("${sample}.vcf"), emit: vcf
  path("${sample}.vcf.idx"), emit: vcfindex
  val(sample), emit: sample

  script:
  """
    gatk HaplotypeCaller \
    --java-options "-Xmx4g" \
    -R $fasta \
    -O ${sample}.vcf \
    -I $bam_bqsr \
    --QUIET

  """
}


workflow {
        
        Channel
                .fromPath(params.fasta)
                .ifEmpty { exit 1, "fasta annotation file not found: ${params.fasta}" }
                .view()
                .set { ch_fasta }

        Channel
                .fromPath(params.fai)
                .ifEmpty { exit 1, "fai annotation file not found: ${params.fai}" }
                .set { ch_fai }
        
        Channel
                .fromPath(params.dict)
                .ifEmpty { exit 1, "dict annotation file not found: ${params.dict}" }
                .set { ch_dict }
        
        Channel
                .fromPath(params.bai)
                .ifEmpty { exit 1, "bai annotation file not found: ${params.bai}" }
                .set { ch_bai }

        
        HAPLOTYPCALLER( ch_bai, ch_fasta.collect(), ch_fai.collect(), ch_dict.collect() )
        
}
    