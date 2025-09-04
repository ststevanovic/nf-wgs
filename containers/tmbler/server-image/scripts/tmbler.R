 #!/usr/bin/env Rscript

library("devtools")
library("TMBleR")

args <- commandArgs(trailingOnly=TRUE)

vcfs <- eval(parse(text=(args[1])))
genome_assembly <- as.character(args[2])
design_file <- as.character(args[3])


v_keys = sapply(vcfs, "[[", 1)
v_values = sapply(vcfs, "[[", 2)
vcf_files <- setNames(as.list(v_values), v_keys)

vcfs <- readVcfFiles(vcfFiles = vcf_files, assembly = genome_assembly)

design_file_path = system.file( "extdata"
                        , design_file
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
        
write.table(TMB_res, file='TMBquant.tsv', quote=FALSE, sep='\t')  