# SnpEff module

SnpEff is a variant annotation and effect prediction tool. It annotates and predicts the effects of genetic variants (such as amino acid changes).

A typical SnpEff use case would be:

1. Input: The inputs are predicted variants (SNPs, insertions, deletions and MNPs). The input file is usually obtained as a result of a sequencing experiment, and it is usually in variant call format (VCF).  
2. Output: SnpEff analyzes the input variants. It annotates the variants and calculates the effects they produce on known genes (e.g. amino acid changes). A list of effects and annotations that SnpEff can calculate can be found [here](http://pcingola.github.io/SnpEff/se_inputoutput/#effect-prediction-details).  

The following flags are applied for SnpEff in this pipeline:

    -dataDir <path>    = Override data_dir parameter from config file.
    -nodownload        = Do not download a SnpEff database, if not available locally.
    -csvStats          = Create CSV summary file instead of HTML
    -canon             = Only use canonical transcripts.
    -v                 = verbose mode
# Usage
To test the module cd to project / and run:
```
nextflow -log logs/snpeff.log run modules/SnpEff/snpeff.nf -c modules/SnpEff/snpeff.config 
```
See the arguments section to apply changes.

# Arguments
### Required: 
    vcf          = Variant call file (default: "$rootDir/results/TEST/Haplotyper/*.vcf")
    snpeff_cache = Path to SnpEff data directory (default: null)
    snpeff_db    = SnpEff database, currently set to dock.image (default: "GRCh38.86")

### Optional: 
    publish_dir_mode  = "copy"
    outdir            = Output directory (default: "$rootDir/results/TEST/SNPeff/")

# Outputs

    $sample_id.snpEff.genes.txt     = Reports in txt format  
    $sample_id.snpEff.html          = Reports in html format  
    $sample_id.snpEff.csv           = Reports in csv format  
    $sample_id.snpEff.ann.vcf       = Annotated VCF file  
