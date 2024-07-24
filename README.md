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


# Conditional Analysis

# Fine-mapping

# Identification of TR QTLs
