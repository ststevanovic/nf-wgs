# Qualimap bamqc module

Qualimap is a platform-independent application written in Java and R that provides both a Graphical User Inteface (GUI) and a command-line interface to facilitate the quality control of alignment sequencing data and its derivatives like feature counts. Bam QC reports information for the evaluation of the quality of the provided alignment data (a BAM file). In short, the basic statistics of the alignment (number of reads, coverage, GC-content, etc.) are summarized and a number of useful graphs are produced. This analysis can be performed with any kind of sequencing data, e.g. whole-genome sequencing, exome sequencing, RNA-seq, ChIP-seq, etc.

This module workflow collects BAM file from results directory and do quality mapping upon these files.
# Usage
To test the module use:
```
nextflow -log logs/qualimap.log run modules/Qualimap/qualimap.nf -c modules/Qualimap/qualimap.config -entry qualimap_workflow 
```
See the arguments section to apply changes.
# Arguments:
### Required: 
    --inpcrd            "/path/to/bam" (default: results/TEST/*/*.bam)
### Optional: 
    --outdir            "/path/to/individual_results" (default: results/TEST/Qualimap/individual_results)
    --output_format     "PDF / html " (default: html)
    --help              Print out help (default: false)

### Outputs:
    report files in html/PDF format
