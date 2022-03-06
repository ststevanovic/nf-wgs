# Manta module

This module uses Manta tool. Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads. It is optimized for analysis of germline variation in small sets of individuals and somatic variation in tumor/normal sample pairs. Manta discovers, assembles and scores large-scale SVs, medium-sized indels and large insertions within a single efficient workflow. Manta divides the SV and indel discovery process into two primary steps: (1) scanning the genome to find SV associated regions and (2) analysis, scoring and output of SVs found in such regions.

## Inputs

--paireads        pair-end sequencing assay (or reads) in BAM/CRAM format
--reads        non-pair reads in BAM/CRAM format
--bai        index file matching the input BAM/CRAM file
--refbam     reference BAM/CRAM file

## Outputs
### primary
diploidSV.vcf.gz                 VCF 1.4 files, genotyped under a diploid model
somaticSV.vcf.gz                 VCF 1.4 files, genotyped under a somatic variant model
candidateSV.vcf.gz               Unscored SV and indel candidates
candidateSmallIndels.vcf.gz      simple insertion and deletion variants
#### tumor only analysis
tumorSV.vcf.gz  subset of candidateSV with score values above minimum scored variant size (default: 50)

### secondary
alignmentStatsSummary.txt       fragment length quantiles for each input alignment file
svLocusGraphStats.tsv           statistics and runtime information pertaining to the SV locus graph
svCandidateGenerationStats.tsv  statistics and runtime information pertaining to the SV candidate generation
svCandidateGenerationStats.xml  xml data backing the svCandidateGenerationStats.tsv report


## Usage

To use Manta module, Single Manta process run the following code from the root directory
```
nextflow -log logs/mantaSingle.log run modules/Manta/manta.nf -entry single_workflow -c modules/Manta/manta.config
```
To use Manta module, Somatic Pair Manta process run the following code from the root directory
```
nextflow -log logs/mantaPaired.log run modules/Manta/manta.nf -entry paired_workflow -c modules/Manta/manta.config
```

The module will produce {sample_name} directory containing primary and secondary manda outputs
The system first configures with
```
configManta.py \
        ${inputbam} ${bam} \
        --reference ${fasta} \
        ${options} \
        --runDir Manta
```
then it produce outputs
```
python Manta/runWorkflow.py -m local -j ${task.cpus}
```

