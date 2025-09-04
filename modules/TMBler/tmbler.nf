nextflow.enable.dsl=2
log.info """
nextflow -log logs/tmbler.log run modules/TMBler/tmbler.nf -c modules/TMBler/tmbler.config
"""

// Note 2 combine collect <-
// Note 3 fastq propagira <-
process TMBler {
    label "tmbler_module"
    //   Note 1: port host needs to mirror the guest 
    containerOptions "--name tmbler \
                  --rm \
                  -e USER=\$(whoami) \
                  -e PASSWORD=helloworld \
                  -e USERID=\$UID \
                  -p ${port}:${port} \
                  -v ${PWD}:${PWD}"

    publishDir "${outdir}", mode: params.publish_dir_mode

    input:
    tuple val(port), val(vcfs)
    val genome
    val outdir
    
    output:
    stdout
    path "*.tsv", emit: tmbout

    script:
    genome = params.genome
    """
    #!/usr/bin/env Rscript
    
    library("devtools")
    library("TMBleR")

    genome_assembly <- "${genome}"
    
    v_keys = sapply(${vcfs}, "[[", 1)
    v_values = sapply(${vcfs}, "[[", 2)
    vcf_files <- setNames(as.list(v_values), v_keys)
    
    vcfs <- readVcfFiles(vcfFiles = vcf_files, assembly = genome_assembly)

    design_file_path = system.file( "extdata"
                            , "ExamplePanel_GeneIDs.txt"
                            , package = "TMBleR"
                            , mustWork = TRUE)
    
    design <- readDesign(
                            filename = design_file_path,
                            assembly = genome_assembly,
                            ids = "entrezgene_id"
                            )
        

    vcfs_NoCancer_NoSynonymous <- applyFilters(  
                            vcfs = vcfs,
                            assembly = genome_assembly,
                            design = design,
                            remove.cancer = T,
                            variantType = c("synonymous")
                            )
    
    vcfs_nonfiltered <- applyFilters(
                        vcfs = vcfs,
                        assembly = genome_assembly, 
                        design = design
                        ) 


    vcfs_all <- c(
                vcfs_NoCancer_NoSynonymous,
                vcfs_nonfiltered
                )

    TMB_res <- applyTMB(
                inputForTMB = vcfs_all, 
                assembly = genome_assembly
                )
           
    DT::datatable(TMB_res)
    write.table(TMB_res, file='TMB_res.tsv', quote=FALSE, sep='\t')  
    
    """
}



process test_inputs {
  label "tmbler_module"
  
  input:
    tuple val(port), val(vcfs)
  
  output:
    stdout
  
  script:
  """
  #!/usr/bin/env Rscript
  library("devtools")
  library("TMBleR")

  v_keys = sapply(${vcfs}, "[[", 1)
  v_values = sapply(${vcfs}, "[[", 2)
  hh <- setNames(as.list(v_values), v_keys)
  vcf_files <- list(hh)
  
  print(vcf_files)
  print(${port})
  """
}

workflow {
    // Note 2: The input VCFs should belong to patient group 
    // i.e. (single) Somatic mutational profile
    Channel
        .fromPath( params.vcfs )
        .map( file -> "list('${file.baseName}', '${file}')" )
        .collect()
        .map( file -> "list(${file.join(', ')})")
        .set{ ch_vcfs }
    
    Channel
        .of( params.genome )
        .set{ ch_genome }

    Channel
        .value( params.port )
        .set{ ch_port }

    // Note 3: In case of parallel `TMBling`,
    // dynamically allocated ports per patient group
    // it[0] list of tuples, i.g. [ [patient1, patient1.vcf], [patient2, patient2.vcf] ]
    if ( params.vcfs_collections ) {
         Channel
            .fromPath( params.vcfs_collections )
            .set{ ch_vcfs_collections }
       
        ch_vcf_collections
            .collect()
            .withIndex()
            .flatMap()
            .map {
                it -> [it[1]+8787, it[0] ] 
            }
        .set{ ch_groups_ports }
        
        TMBler( ch_groups_ports, ch_genome )
    }

    // ch_vcfs.view()
    ch_port.combine(ch_vcfs).set{ ch_group_port }

    // Regular TMBling
    TMBler( ch_group_port, ch_genome, params.outdir )
    // TMBler.out.view()
}