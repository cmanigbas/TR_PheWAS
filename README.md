# PheWAS in Tandem Repeats (TRs)
This repository contains scripts for performing PheWAS using TR average copy number with various phenotypes in the UK Biobank cohort, including data generation, data normalization, QC steps, PheWAS, and fine-mapping. All coordinates used in this pipeline were in UCSC Human Genome hg38.

# Catalog of TRs 
We derived a catalog of  TRs with motif sizes ranging from 2-20bp that we hypothesized would be enriched for
functional effects and which are either highly polymorphic or were observed to undergo rare expansion in the
human population
  * [ExpansionHunter](https://github.com/Illumina/ExpansionHunter) and [hipSTR](https://github.com/HipSTR-Tool/HipSTR) tools to genotype TRs that overlap regulatory regions
  * [TRFinder annotations from UCSC Genome browser](https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=simpleRepeat&hgta_table=simpleRepeat&hgta_doSchema=describe+table+schema) for TRs composed of purely G and C bases
  * [ExpansionHunter Denovo]() and [STRetch]() tools to identify TRs showing long outlier allele sizes

# TR QC Filtering 
Using All of Us short read GS and long read GS, we performed genotyping of the set of TRs in both the Illumina data using ExpansionHunter and in the PacBio data using [TRGT](https://github.com/PacificBiosciences/trgt). For quality filtering, we removed:
  * TRs with no variation in TRGT genotypes
  * TRs with no variation in ExpansionHunter genotypes
  * TRs with Spearman correlation between TRGT and ExpansionHunter genotype < 0.5
  * TRs with median purity scores from TRGT < 0.75

# PheWAS
For each sequencing center, SC and deCODE, we ran association testing using [REGENIE](https://rgcgithub.github.io/regenie/), incorporating covariates of sex, age, age squared, GS insert size, and top five PCs. For binary traits, we utilized the Saddle Point Approximation function to reduce the type I error rate, and for quantitative traits, we applied a rank based inverse normal transformation. Results for each sub-cohort for combined using [METAL](https://github.com/statgen/METAL). 

# Conditional Analysis
We performed conditional analysis in each TR:trait pair that was significant in PheWAS. We performed association analysis between that trait and all SNVs located within ±500kb using REGENIE. This was performed separately in each of the two sequencing sub-cohorts, before combining p-values together using METAL. We removed SNVs that had either minor allele count<100, missingness rate >5%, HWE<10-300, quality score (QUAL) <30, mapping quality (MQ) <40, read depth (DP)
<10, genotype quality (GQ) <20, strand bias (SB) <0.25 or >0.75 and allelic balance for heterozygous calls (ABhet) <0.2. To condition the TR, we separated individuals by genotype state based on the lead associated SNV with the trait. 

# Fine-mapping
We performed statistical fine mapping using [CAVIAR](http://genetics.cs.ucla.edu/caviar/) on the Sanger Center sub-cohort using the
top 100 most significantly associated SNVs per locus, as defined above. We calculated LD between SNVs and
TRs by pairwise correlation. P-values from REGENIE were converted into z-scores and CAVIAR run using parameters ρ=0.95, γ=0.01 and the maximum number of causal variants was set to 2


# Identification of TR QTLs
