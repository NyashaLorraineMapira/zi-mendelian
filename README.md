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

## Parameters

### Input / Output
| Parameter | Description | Default |
|---------|-------------|---------|
| `--vcf_dir` | Directory containing input VCF files | Project directory |
| `--vcf_pattern` | Pattern used to match VCF files | `S*.hard-filtered.vcf` |
| `--outdir` | Output directory | `results/` |

### Annotation
| Parameter | Description | Default |
|---------|-------------|---------|
| `--run_vep` | Enable VEP annotation step | `false` |
| `--vep_exe` | Path to VEP executable | `vep` |
| `--vep_cache` | Path to VEP cache directory | Not set |

### ANNOVAR
| Parameter | Description | Default |
|---------|-------------|---------|
| `--annovar_dir` | Path to ANNOVAR installation | `~/annovar` |
| `--humandb` | ANNOVAR human database directory | `annovar/humandb` |
| `--build` | Genome build | `hg38` |
| `--annovar_protocol` | ANNOVAR databases used | `refGene,ensGene,gnomad,clinvar` |
| `--annovar_operation` | Operations for ANNOVAR protocols | `g,g,f,f,f` |

### Population Prioritisation
| Parameter | Description | Default |
|---------|-------------|---------|
| `--run_zim_layer` | Enable Zimbabwe/African prioritisation layer | `true` |
| `--zim_db` | Zimbabwe cohort frequency file | `refs/zimbabwe_cohort_freq.tsv` |

### Resources
| Parameter | Description | Default |
|---------|-------------|---------|
| `--max_cpus` | Maximum CPUs per process | `4` |
| `--max_mem` | Maximum memory per process | `7 GB` |

### Low-resource execution
For machines with limited memory and CPUs:

```bash
nextflow run main.nf -profile low_resource
