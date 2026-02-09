#!/usr/bin/env bash
set -euo pipefail

# Run from project root: bash tests/run_pipeline_test.sh

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

VCF_DIR="${ROOT_DIR}/tests/test_data"
OUTDIR="${ROOT_DIR}/tests/results/pipeline_test"
VCF_PATTERN="S_test.hard-filtered.vcf"
ZIM_DB="${ROOT_DIR}/tests/test_data/test_zim_db.tsv"

mkdir -p "${OUTDIR}"

echo "[INFO] Running Zi-Mendelian annotation pipeline on synthetic VCF..."
nextflow run "${ROOT_DIR}/main.nf" \
  --vcf_dir "${VCF_DIR}" \
  --vcf_pattern "${VCF_PATTERN}" \
  --outdir "${OUTDIR}" \
  --run_zim_layer true \
  --run_vep false \
  --zim_db "${ZIM_DB}"

echo "[INFO] Comparing pipeline outputs to expected files..."

diff -u \
  "${ROOT_DIR}/tests/expected/S_test.hard-filtered.hg38_multianno.txt" \
  "${OUTDIR}/annovar/annotated/S_test.hard-filtered.hg38_multianno.txt"

diff -u \
  "${ROOT_DIR}/tests/expected/S_test.hard-filtered.zim_prioritised.tsv" \
  "${OUTDIR}/prioritisation/S_test.hard-filtered.zim_prioritised.tsv"

echo "[OK] full Zi-Mendelian annotation curated test passed"
