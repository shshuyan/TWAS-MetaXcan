#!/bin/bash
# Step 5: summary imputation for one chromosome ($1) and sub-batch ($2)
set -euo pipefail; source "$(dirname "$0")/config.sh"
CHR=$1; SB=$2
OUT=$OUTPUT/summary_imputation/${TRAIT}_chr${CHR}_sb${SB}_reg0.1_ff0.01_by_region.txt.gz
[ -s "$OUT" ] && { echo "skip chr$CHR sb$SB"; exit 0; }
mkdir -p $OUTPUT/summary_imputation $OUTPUT/logs
$PY_IMP $GWAS_TOOLS/gwas_summary_imputation.py \
 -by_region_file $DATA/eur_ld.bed.gz \
 -gwas_file $OUTPUT/harmonized_gwas/$TRAIT.txt.gz \
 -parquet_genotype $DATA/reference_panel_1000G/chr${CHR}.variants.parquet \
 -parquet_genotype_metadata $DATA/reference_panel_1000G/variant_metadata.parquet \
 -window 100000 -parsimony 7 -chromosome $CHR \
 -regularization 0.1 -frequency_filter 0.01 \
 -sub_batches 10 -sub_batch $SB --standardise_dosages \
 -output $OUT > $OUTPUT/logs/impute_chr${CHR}_sb${SB}.log 2>&1
echo "done chr$CHR sb$SB"
