# RUNX1 ChIP-seq Analysis BF528
## Overview
RUNX1 is a transcription factor with essential roles in hematopoiesis and leukemia. This project aims to reproduce and analyze RUNX1 binding profiles using publicly available ChIP-seq datasets, integrating the results with RNA-seq data. Analyses replicate Figure 2 and relevant supplementary materials from the original publication:[Gao et al., PMC5071180](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5071180/).

* Perform quality control and preprocessing of ChIP-seq reads.
* Align reads to the human genome (GRCh38) and call peaks.
* Identify reproducible peaks and annotate them to genomic features.
* Conduct motif enrichment analysis.
* Visualize coverage profiles across genes of interest.
* Compare ChIP-seq results with RNA-seq to explore regulatory effects.

## Dependencies
* Package management: Miniconda, Bioconda, Conda-Forge
* Workflow & Programming Dependencies:
  * snakemake v8.5.2
  * snakemake-executor-plugin-cluster-generic v1.0.7
  * pandas v2.2.1
  * jupyterlab v4.1.2
* High-performance computing: Boston University Shared Computing Cluster (SCC)
* Tools:
  * BEDTools v2.31.1
  * Bowtie2 v2.5.3
  * Deeptools v3.5.4
  * FastQC v0.12.1-0
  * HOMER v4.11
  * MultiQC v1.20
  * SAMtools v1.19.2
  * Trimmomatic v0.39
 
## Workflow Summary
### 1. Data Acquisition & QC
* Download ChIP-seq reads from EMBL-ENA using wget.
* Perform initial quality control with FastQC.
* Trimm adapters and low-quality bases using Trimmomatic.
* Generate Bowtie2 genome index for GRCh38.
* QC Insights
  * IP samples showed higher duplication rates due to enrichment of RUNX1-bound fragments.
  * GC content differed between IP and Input, reflecting sequence bias from transcription factor binding.

### 2. Alignment & Processing
* Align reads with Bowtie2.
* Convert and sort alignments using SAMtools.
* Aggregate QC reports with MultiQC.
* Generate bigWig files for coverage visualization (deeptools bamCoverage).
* Assess sample correlations (multiBigwigSummary and plotCorrelation).

### 3. Peak Calling & Annotation
* Create tag directories and call peaks using HOMER.
* Merge replicates and filtered reproducible peaks using BEDTools.
* Annotate peaks to genomic features and remove blacklisted regions.
* Perform motif enrichment analysis using HOMER.

### 4. Data Visualization & Integration
* Generate gene body signal profiles using deeptools plotProfile.
* Visualize peaks in key loci (MALAT1, NEAT1) using Integrative Genomics Viewer (IGV).
* Integrate ChIP-seq and RNA-seq data to reproduce figures from the publication.
  
## Results
* Replicate peaks: Replicate 1 = 106,072, Replicate 2 = 50,952
* Reproducible peaks (intersect): 10,924
* Filtered reproducible peaks: 10,785

## Aknowledgements
* Data from
  
  Gao, Y. et al., PMC5071180. "RUNX1 ChIP-seq study in hematopoietic cells." Nucleic Acids Research, 2016.

## Author
Mabel Chang

Developed as the final project for BF528: Applications in Translational Bioinformatics at Boston University (S24)
