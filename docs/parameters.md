# Pipeline Parameters

This document describes all user-configurable parameters in the *Pf Population Genomics Pipeline*.

---

## Required Parameters

| Parameter | Description | Example |
|----------|-------------|---------|
| `--reads` | Input paired-end FASTQ files (must follow pattern *_1.fastq.gz and *_2.fastq.gz) | `"data/*_{1,2}.fastq.gz"` |
| `--reference_fasta` | Reference genome FASTA file | `Pf3D7.fasta` |
| `--kraken_db` | Kraken2 database path | `/path/to/kraken_db` |
| `--snpeff_data_dir` | Directory containing SnpEff databases | `/path/to/snpeff/data` |

---

## Optional Parameters

| Parameter | Default | Description |
|----------|--------|-------------|
| `--outdir` | `results/` | Output directory |
| `--snpeff_db` | `Plasmodium_falciparum` | SnpEff database name |

---

## Example Run

```bash
nextflow run main.nf \
  --reads "data/*_{1,2}.fastq.gz" \
  --reference_fasta Pf3D7.fasta \
  --kraken_db /db/kraken \
  --snpeff_data_dir /db/snpeff \
  --outdir results

