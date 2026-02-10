Population Genetic Analysis of Plasmodium falciparum Isolates Across Transmission Gradients in Cameroon

Project Overview:

This repository contains a reproducible Nextflow DSL2 pipeline for automated processing of Plasmodium falciparum Illumina paired-end sequencing data, starting from raw FASTQ files and producing high-quality, annotated genomic VCF (gVCF) files of SNPs and INDELs for downstream population genetic analyses.
The workflow is designed for scalability, reproducibility, and compatibility with large public malaria genomics datasets (e.g., Pf3k, Pf7, ENA).

## Pipeline Workflow:

The pipeline integrates established bioinformatics tools to ensure robust variant discovery, including tools like:
- FASTP for adapter trimming and removal of low quality bases (Quality control & trimmingof raw reads)
- KRAKEN2 for deconvolution and removal of human host reads from parasite reads
- BWA-MEM for mapping reads to the P. falciparum reference genome (Read Alignment)
- PICARD and VCFtools for marking and removing duplicates, BAM sorting and indexing (Post-alignment processing of BAM files)
- GATK (pamgen best practice) for SNP and INDEL variants discovery & quality filtering (Variant Calling)
- BCFtools + custom BED files for masking hypervariable and low-confidence genomic regions Masking problematic regions
- SNPEFF with custom database for functional annotation of variants
- Generation of chromosome-specific VCF files (14 nuclear + apicoplast + mitochondrial) with BCFtools
- It also provides quality control steps after key processes of the pipeline with multiQC reporting

FASTQ
  │
  ├── FastQC → MultiQC
  │
  ├── fastp (trim adapters)
  │
  ├── FastQC (post-trim)
  │
  ├── Kraken2 (host removal)
  │
  ├── BWA-MEM (alignment)
  │
  ├── MarkDuplicates (GATK)
  │
  ├── HaplotypeCaller (gVCF)
  │
  ├── GenotypeGVCFs (cohort VCF)
  │
  ├── SelectVariants (SNPs)
  │
  └── SnpEff (annotation)

## Pipeline Outputs:

The pipeline produces:
* Parasite-classified reads
* BAM files with quality assessment of the alignment output
* High-quality gVCF files
* Masked and filtered snp-based and indel-based VCF files
* Annotated chromosome-specific snp-based VCF datasets

These outputs are optimized for downstream population genetic analyses such as:
Population structure inference
Transmission dynamics studies
Evolutionary history reconstruction
Drug resistance surveillance
Sample clustering and relatedness analyses

## Requirements

- Nextflow ≥ 23
- Conda / Mamba
- 32 GB RAM recommended

## HOW TO RUN THE PIPELINE
Use nextflow run main.nf -profile conda 
(N.B. You do not need to manually create environments as Nextflow does it automatically. The first run takes time as the envs build but after that, it's cached).

Relevance:

This workflow is particularly suited for studies investigating:
Malaria parasite evolution
Genomic variation across transmission intensity gradients
Regional and global parasite population structure
It enables standardized and reproducible analysis across both local datasets and large public resources.
