# TMBler

TMB has been proposed as a predictive biomarker for immunotherapy response in cancer patients, as it is thought to enrich for tumors with high neoantigen load. TMB assessed by Whole Exome Sequencing (WES) is considered the gold standard but remains confined to research settings. Targeted enrichment panels of various genomic sizes are emerging as a more sustainable methodology for assessing TMB in the clinical setting.

# Usage

```
    nextflow -log logs/tmbler.log run modules/TMBler/tmbler.nf -c modules/TMBler/tmbler.config
```

Note: VAFFilter is not implemented due to vcf-content dependencies.  
This filter might break the pipeline.  
Thus, parameters which can be included in applyFilters but are not required:  
```
vaf.cutoff = 0.05 
tsList = NULL

```

## Input params
1. Required  
- genome : "hg19" or "hg38" (default: hg19)) 
- vcfs : .vcf Path to vcfs obtained in WGS/WES   
- outdir : Results output directory (default "results/TBMler/")
2. Optional:  
- gene_panels: [.bam, .txt ], default  "ExamplePanel_GeneIDs.txt"  

## Outputs

1. TMB_out.tsv : Tumor Mutationl Budren (TMB) quantification file