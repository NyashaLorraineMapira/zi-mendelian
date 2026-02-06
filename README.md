# Zi-Mendelian

Zi-Mendelian is a Nextflow DSL2 pipeline for Mendelian variant
prioritisation in rare disease genomics. The pipeline integrates
ANNOVAR-based variant annotation with a population-aware filtering
layer informed by a Zimbabwean cohort.

## Pipeline overview
Zi-Mendelian performs the following steps:
1. Optional VEP annotation of input VCF files
2. Conversion of VCF files to ANNOVAR input format
3. Multi-database ANNOVAR annotation
4. Variant prioritisation using Zimbabwe/African population frequencies

## Requirements
- Nextflow (>= 22.10, DSL2 enabled)
- Python 3
- ANNOVAR
- bgzip and tabix

## Input
- Hard-filtered VCF files (e.g. `S*.hard-filtered.vcf`)

## Usage

```bash
nextflow run main.nf \
  --vcf_dir /path/to/vcfs \
  --outdir results
Enable VEP annotation:

nextflow run main.nf --run_vep true
Output
ANNOVAR multianno annotation files

Zimbabwe-prioritised variant tables (*.zim_prioritised.tsv)

Reproducibility
Zi-Mendelian is implemented using Nextflow DSL2 to ensure reproducible
and portable execution across computational environments.

License
MIT
