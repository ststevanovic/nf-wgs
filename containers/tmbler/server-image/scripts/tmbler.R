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
    

vcfs_NoCancer <- applyFilters(  
                        vcfs = vcfs,
                        assembly = genome_assembly,
                        design = design,
                        remove.cancer = T
                        )

vcfs_NoSynonymous <- applyFilters(
                        vcfs = vcfs, 
                        assembly = genome_assembly, 
                        design = design, 
                        variantType = c("synonymous")
                        )

vcfs_NoCancer_NoSynonymous <- applyFilters(
                        vcfs = vcfs, 
                        assembly = genome_assembly,
                        design = design, 
                        remove.cancer = T, 
                        variantType = c("synonymous")
                        )

vcfs_NoSynonymous_VAFFilter <- applyFilters(
                        vcfs = vcfs, 
                        assembly = genome_assembly, 
                        design = design, 
                        vaf.cutoff = 0.05, 
                        variantType = c("synonymous")
                        )

vcfs_VAFFilter <- applyFilters(
                        vcfs = vcfs, 
                        assembly = genome_assembly, 
                        design = design, 
                        vaf.cutoff = 0.05,
                        tsList = NULL, 
                        remove.cancer = T
                        )

vcfs_NoCancer_VAFFilter <- applyFilters(
                        vcfs = vcfs, 
                        assembly = genome_assembly, 
                        design = design, 
                        vaf.cutoff = 0.05, 
                        remove.cancer = T
                        )

vcfs_NoCancer_VAFFilter_NoSynonymous <- applyFilters(
                        vcfs = vcfs, 
                        assembly = genome_assembly, 
                        design = design, 
                        vaf.cutoff = 0.05, 
                        remove.cancer = T, 
                        variantType = c("synonymous")
                        )


vcfs_nonfiltered <- applyFilters(
                        vcfs = vcfs,
                        assembly = genome_assembly, 
                        design = design
                        ) 


vcfs_all <- c(
            vcfs_NoCancer, 
            vcfs_NoCancer_NoSynonymous, 
            vcfs_NoSynonymous, 
            vcfs_NoSynonymous_VAFFilter, 
            vcfs_VAFFilter, 
            vcfs_NoCancer_VAFFilter, 
            vcfs_NoCancer_VAFFilter_NoSynonymous, 
            vcfs_nonfiltered 
            )

TMB_res <- applyTMB(
            inputForTMB = vcfs_all, 
            assembly = genome_assembly
            )
        
DT::datatable(TMB_res)
write.table(TMB_res, file='TMB_res.tsv', quote=FALSE, sep='\t')  