# Zi-Mendelian (Nextflow)

Zi-Mendelian is a Nextflow DSL2 pipeline for Mendelian variant
prioritisation in rare disease genomics, integrating ANNOVAR
annotation and a Zimbabwe/African population frequency layer.

## Features
- Nextflow DSL2
- Optional VEP annotation
- ANNOVAR-based multi-annotation
- Population-aware prioritisation using Zimbabwe cohort data
- Designed for low-resource compute environments

## Requirements
- Nextflow (>= 22.10)
- Python 3
- ANNOVAR
- bgzip, tabix

## Usage

```bash
nextflow run main.nf \
  --vcf_dir /path/to/vcfs \
  --outdir results
