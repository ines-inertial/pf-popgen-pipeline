## All results or outputs are written to the directory specified by '--outdir'

results/
│
├── qc/
│ ├── fastqc_raw/
│ ├── fastqc_trim/
│ └── multiqc_raw/
│
├── trimmed/
│ └── *.trimmed.fq.gz
│
├── host_filter/
│ └── clean.fq.gz
│
├── alignment/
│ └── *.sorted.bam
│
├── dedup/
│ └── *.dedup.bam
│
├── gvcf/
│ └── *.g.vcf.gz
│
├── vcf/
│ └── cohort.vcf.gz
│
├── vcf_filtered/
│ └── cohort.snps.vcf.gz
│
└── annotation/
  └── annotated.vcf 
