#!/usr/bin/env nextflow

/*
================================================================================
                                  nf-core/sarek
================================================================================
Started March 2016.
Ported to nf-core May 2019.
--------------------------------------------------------------------------------
nf-core/sarek:
  An open-source analysis pipeline to detect germline or somatic variants
  from whole genome or targeted sequencing
--------------------------------------------------------------------------------
 @Homepage
 https://nf-co.re/sarek
--------------------------------------------------------------------------------
 @Documentation
 https://nf-co.re/sarek/docs
--------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run nf-core/sarek/download_cache.nf -profile docker --genome <genome> --help
                                      [--snpeff_cache <pathToSNPEFFcache> --snpeff_db_version <snpEff DB version>]
                                      [--vep_cache <pathToVEPcache> --vep_cache_version <VEP cache version> --species <species>]
                                      [--cadd_cache <pathToCADDcache> --cadd_version <CADD Version>]
    Options:
      --help                   [bool] You're reading it
      --snpeff_cache           [file] Path to snpEff cache
      --snpeff_db_version       [str] snpEff DB version
                                      Default: ${params.genomes[params.genome].snpeff_db}
      --vep_cache              [file] Path to VEP cache
      --vep_cache_version       [int] VEP cache version
                                      Default: ${params.genomes[params.genome].vep_cache_version}
      --species                 [str] Species
                                      Default: ${params.genomes[params.genome].species}
      --cadd_cache             [file] Path to CADD cache
      --cadd_version            [str] CADD version to download
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help) exit 0, helpMessage()