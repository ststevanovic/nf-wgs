nextflow.enable.dsl=2

// input params for nextflow.config
// bwa mem
params.publish_dir_mode       = "copy"
params.reads_dir              = "$rootDir/data/reads"
params.outdir                 = "$rootDir/results/STANDARD"
params.sequencing_center      = null
params.fasta                  = "$rootDir/data/reference_genome/test.fa"
params.fai                    = "$rootDir/data/reference_genome/test.fa.fai"
params.golden_indel_gz        = 's3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz'
params.golden_indel_idx_gz    = 's3://ngi-igenomes/igenomes/Homo_sapiens/GATK/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx.gz'
    // TODO: --> http://www.htslib.org/doc/tabix.html
    // TODO: snpdb by date and version
params.dbsnp_gz               = "$rootDir/data/dbsnp/test.dbsnp.vcf.gz"
params.dbsnp_idx_gz           = "$rootDir/data/dbsnp/test.dbsnp.vcf.gz.tbi" 
// haplotyper 
params.fai                    = "$rootDir/data/reference_genome/test.fa.fai"
// TODO: should emit
params.dict                   = "$rootDir/results/TEST/GATK/test.dict"
// TODO: should emit
params.bai                    = "$rootDir/results/TEST/GATK/*.bqsr.bam.bai"
// manta
// TODO: should emit
    // bambai testing: "$rootDir/results/TEST/GATK/*germline*_bqsr{.bam,.bam.bai}"
params.bambai                 = "$rootDir/data/novak-bam-bai/*{.bam,.bam.bai}"
    // Germline/tumor setup 
params.with_tumor             = false 
// TODO: should channel from emited with germline_option
params.bambai_tumor           = "$rootDir/results/TEST/GATK/*germline_bqsr{.bam,.bam.bai}"
// qualimap
// TODO: should emit
params.inpdir                 = "$baseTestDir/results/TEST/Haplotyper" // or Manta
params.output_format          = 'html'
params.help                   = false

// Annotations| snpeff, vep
params.vcf                    = "$rootDir/results/TEST/*/*.vcf"
params.snpeff_db              = "GRCh38.86"
params.snpeff_cache           = null
params.genome                 = "GRCh38" 



log.info """\
         W G S - N F   P I P E L I N E    
         ===================================
         reads              : ${params.reads}
         """
         .stripIndent()

include { FASTQC } from './modules/FastQC/fastqc.nf'
include { dictionary, BWAINDEX, BWAMEM_SAMTOOLS_SORT } from './modules/BWAmem/bwa.nf'
include { MARKDUPLICATES, BASERECALIBRATOR, MARKDUPLICATES } from './modules/BWAmem/bwa.nf'
include { HAPLOTYPECALLER } from './modules/Haplotyper/haplotyper.nf'
include { MANTA } from './modules/Manta/manta.nf'
include { QUALIMAP } from './modules/Qualimap/qualimap.nf'
include { SNPEFF } from './modules/SNPeff/snpeff.nf'
include { VEP } from './modules/VEP/vep.nf'
include { MULTIQC } from './modules/MultiQC/MultiQC.nf'


// Check this out!
// import ParamsChecker
// include {printHeader; helpMessage} from './help' params(params)
// include { BuildBWAFastaIndexes } from './bwa/buildIndex' params(params)


workflow { 
    // Channels
    reads_dir = params.reads_dir
    Channel
        .fromFilePairs(             reads_dir + '/*_{1,2}.fq',       size: -1  )
        .mix( Channel.fromFilePairs(reads_dir + '/*_{1,2}.fq.gz',    size: -1) )
        .mix( Channel.fromFilePairs(reads_dir + '/*_{1,2}.fastq',    size: -1) )
        .mix( Channel.fromFilePairs(reads_dir + '/*_{1,2}.fastq.gz', size: -1) )
        .ifEmpty { exit 1, "Reads not found " }
        .set{ ch_reads }

    Channel
        .fromPath( params.fasta )
        .ifEmpty { exit 1, "fasta file not found: ${params.fasta}" }
        .set{ ch_fasta }
    
    Channel
        .fromPath(params.fai)
        .ifEmpty { exit 1, "fai annotation file not found: ${params.fai}" }
        .set { ch_fai }

    dict_baserecalibrator = dictionary( ch_fasta )   

    ch_fai
        .mix( dict_baserecalibrator )
        .set { ch_index }

    Channel
        .fromPath( params.dbsnp_gz )
        .ifEmpty { exit 1, "dbsnp annotation file not found: ${params.dbsnp_gz}" }
        .set { ch_dbsnp }

    Channel
        .fromPath( params.dbsnp_idx_gz )
        .ifEmpty { exit 1, "dbsnp_idx_gz annotation file not found: ${params.dbsnp_idx_gz}" }
        .set { ch_dbsnp_idx }

    Channel
        .fromPath(params.bai)
        .ifEmpty { exit 1, "bai annotation file not found: ${params.bai}" }
        .set { ch_bai }
    
    // manta config
    // bambai              = "$rootDir/data/novak-bam-bai/*{.bam,.bam.bai}"
    Channel
        .fromFilePairs(params.bambai)
        .ifEmpty { exit 1, "bam and bam index files not found: ${params.bambai}" }
        .map{ it -> [it[0], [it[1]]].flatten() }
        // .set{ch_sample} update
        .set{ ch_bam_bai }

    // qualimap inpdir          = "$baseTestDir/results/TEST/Haplotyper"
    Channel
        .fromPath( params.inpdir + '/*.bam' )
        .ifEmpty { error "Cannot find any bam file: ${params.inpdir}" }
        .map{ file -> tuple(file.baseName, file) } 
        .groupTuple(by: 0)
        .collect() // if getting all bams from haplotyper and manta CHECK!!!
        .set{ ch_bams }

    // snpeff 
    // snpeff_db       = "GRCh38.86"
    // snpeff_cache    = null
    ch_snpeff_cache = params.snpeff_cache ? Channel.value(file(params.snpeff_cache)) : "null"
    ch_snpeff_db = params.snpeff_db ? Channel.value(params.snpeff_db) : "null"

    //  snpeff vcf                 = "$rootDir/results/TEST/Haplotyper/*.vcf"
    // vep vcf                 = "$rootDir/results/TEST/Haplotyper/*.vcf"
    Channel 
        .fromPath(params.vcf)
        .ifEmpty { error "Cannot find any vcf file: ${params.vcf}" }
        .map{file -> tuple(file.baseName, file)}
        .set{ ch_vcf }

    // Preprocess
    FASTQC ( ch_reads )
    
    // Alignment
    BWAINDEX ( ch_fasta )
    BWAMEM_SAMTOOLS_SORT ( ch_reads, ch_fasta.collect(), BWAINDEX.out.bwa_index.collect() )
    
    // Deduplication and Recalibration
    MARKDUPLICATES( BWAMEM_SAMTOOLS_SORT.out.bam ) 
    BASERECALIBRATOR( 
        MARKDUPLICATES.out.bam_markdup,
        ch_fasta.collect(),
        ch_index.collect(),
        ch_dbsnp.collect(),
        ch_dbsnp_idx.collect(),
    )
    BASERECALIBRATOR.out.baserecalibrator_table
        .join( ch_bam_dedup )
        .set{ ch_table_bam }
    APPLYBQSR( 
        ch_table_bam,
        ch_fasta.collect(), 
        ch_index.collect() 
        )

    // Variant Calling
    HAPLOTYPECALLER( ch_bai, ch_fasta.collect(), ch_fai.collect(), ch_dict.collect() )
    MANTA( ch_bam_bai, ch_fasta.collect(), ch_fai.collect() )
    
    // variant calls quality ???
    QUALIMAP( ch_bams )

    // Annotation
    SNPEFF(
        ch_vcf,
        ch_snpeff_cache,
        ch_snpeff_db
    )
    VEP( ch_vcf )

    // Reports
    Channel 
        .of  ( FASTQC.out )
        .mix ( BWAINDEX.out )
        .mix ( BWAMEM_SAMTOOLS_SORT.out )
        .mix ( MARKDUPLICATES.out )
        .mix ( BASERECALIBRATOR.out )
        .mix ( APPLYBQSR.out )
        .mix ( HAPLOTYPECALLER.out )
        .mix ( MANTA.out )
        .mix ( QUALIMAP.out )
        .mix ( SNPEFF.out )
        .mix ( VEP.out ) 
        .collect()
        .set{ ch_reports }
    
    MULTIQC( ch_reports )

}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : null)
}

workflow.onError {
    log.info "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}


