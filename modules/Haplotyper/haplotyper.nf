nextflow.enable.dsl=2

log.info """
  nextflow -log logs/hap.log run modules/Haplotyper/haplotyper.nf -entry haplotyper_module -c modules/Haplotyper/haplotyper.config
"""

process ScatterIntervals {
    tag "$fasta"

    input:
    file fasta
    file index
    file dict


    output:
    path("skip_Ns.interval_list"), emit: scattered_intervals

    script:
    """
    gatk ScatterIntervalsByNs \
    -O skip_Ns.interval_list \
    -R ${fasta}
    """
}

process SkipIntervals {
    tag "$intervals"

    input:
    file(intervals)

    output:
    path("skip_Ns.bed"), emit: skipped_intervals

    script:
    """
    grep -v @ ${intervals}
    #awk '{print \$1":"\$2"-"\$3}'
    awk 'BEGIN { OFS = "" }{ print \$1,":",\$2,"-",\$3 }' ${intervals} > intervals.interval_list
    gatk IntervalListToBed \
    -I ${intervals} \
    -O skip_Ns.bed
    """
}

process MakeIntervals {
    tag "$intervals"

    input:
    file(bed) 

    output:
    path("ACGTmers_interval_size.interval_list"), emit: expanded_intervals

    // TODO: awk coursera
    script:
    """
    cut -f 1-3 $bed > ACGTmers.per-interval.bed
    awk '{print \$1, \$2, \$3, \$3-\$2}' $bed | sort -k4nr | cut -f 1-3 > ACGTmers_interval_size.bed
    awk '{print \$1":"\$2+1"-"\$3}' ACGTmers.per-interval.bed > ACGTmers.per-interval.txt
    bedtools intersect -a ACGTmers.per-interval.bed -wa | awk '{print \$1":"\$2+1"-"\$3}' | sort -u > ACGTmers_interval_target_intersect.txt
    cp ACGTmers_interval_size.bed ACGTmers_interval_size_temp.txt
    awk 'FNR==NR{a[\$1] = \$1;next}{if (\$1 in a) print \$1}' ACGTmers_interval_size_temp.txt ACGTmers_interval_target_intersect.txt > ACGTmers.per-interval_temp.txt
    awk 'FNR==NR{a[\$1] = \$1;next}{if (\$1 in a) print \$1}' ACGTmers_interval_target_intersect.txt ACGTmers.per-interval.txt > ACGTmers_interval_size.interval_list
    mv ACGTmers.per-interval_temp.txt ACGTmers.per-interval.txt
    """
}

 process IndexBam {
  tag "$bam"

  input:
  tuple val(name), file(bam)

  output:
  tuple val(name), file(bam), file("${bam}.bai"), emit: indexed_bam_bqsr

  script:
  """
  samtools index ${bam};
  """
  }

process HaplotypeCaller {
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
  int mem = (Runtime.getRuntime().totalMemory()) >> 30
  """
    gatk HaplotypeCaller \
    --java-options "-Xmx4g" \
    -R $fasta \
    -O ${sample}.vcf \
    -I $bam_bqsr \
    --QUIET

  """
}

workflow haplotyper_module { 

    if (params.fasta) {
        Channel
                .fromPath(params.fasta)
                .ifEmpty { exit 1, "fasta annotation file not found: ${params.fasta}" }
                .view()
                .set { ch_fasta }
    }

        Channel
                .fromPath(params.fai)
                .ifEmpty { exit 1, "fai annotation file not found: ${params.fai}" }
                .set { ch_fai }

        Channel
                .fromPath(params.dict)
                .ifEmpty { exit 1, "dict annotation file not found: ${params.dict}" }
                .set { ch_dict }

        Channel
                .fromPath( params.bambqsr )
                .map { file -> tuple(file.baseName, file) }
                .ifEmpty { exit 1, "bambqsr list file for IndexBam not found: ${params.bambqsr}" }
                .view{ "bam_bqsr: $it" } 
                .set { ch_bam_bqsr }

    IndexBam( ch_bam_bqsr )
    indexed_bam_bqsr = IndexBam.out.indexed_bam_bqsr

    HaplotypeCaller( indexed_bam_bqsr, ch_fasta, ch_fai, ch_dict )
}

// example, not to be used
workflow haplotyper_with_intervals {

     if (params.fasta) {
        Channel
                .fromPath(params.fasta)
                .ifEmpty { exit 1, "fasta annotation file not found: ${params.fasta}" }
                .view()
                .set { ch_fasta }
    }

        Channel
                .fromPath(params.fai)
                .ifEmpty { exit 1, "fai annotation file not found: ${params.fai}" }
                .set { ch_fai }

        Channel
                .fromPath(params.dict)
                .ifEmpty { exit 1, "dict annotation file not found: ${params.dict}" }
                .set { ch_dict }

    if (params.intervals) {
        Channel
                .fromPath(params.intervals)
                .ifEmpty { exit 1, "Interval list file for HaplotypeCaller not found: ${params.intervals}" }
                .set { ch_intervals }
    }
        Channel
                .fromPath( params.bambqsr )
                .map { file -> tuple(file.baseName, file) }
                .ifEmpty { exit 1, "bambqsr list file for IndexBam not found: ${params.bambqsr}" }
                .view{ "bam_bqsr: $it" } 
                .set { ch_bam_bqsr }

    IndexBam( ch_bam_bqsr )
    indexed_bam_bqsr = IndexBam.out.indexed_bam_bqsr

    ScatterIntervals( ch_fasta, ch_fai, ch_dict )
    scattered_intervals = ScatterIntervals.out.scattered_intervals

    SkipIntervals( scattered_intervals )
    skipped_intervals = SkipIntervals.out.skipped_intervals

    MakeIntervals(SkipIntervals.out.skipped_intervals)           
    intervals_file = MakeIntervals.out.expanded_intervals 

    //  params intervals to be True
    intervals_file
            .splitText()
            .map { it -> it.trim() }
            .set { ch_interval }

    ch_fasta
            .mix(ch_fai, ch_dict, indexed_bam_bqsr, intervals_file)
            .set { ch_haplotypcaller_index }
    
    ch_interval
            .combine( ch_haplotypcaller_index )
            .set{ ch_haplotypcaller }

    haplotypecaller( ch_haplotypcaller )
}
    