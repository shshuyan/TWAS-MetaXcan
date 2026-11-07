#!/bin/bash
# Step 7: S-PrediXcan in every GTEx v8 tissue with MASHR models
set -euo pipefail; source "$(dirname "$0")/config.sh"
mkdir -p $OUTPUT/spredixcan/eqtl $OUTPUT/logs
run() {
  T=$1; OUT=$OUTPUT/spredixcan/eqtl/${TRAIT}__PM__${T}.csv
  [ -s "$OUT" ] && return 0
  $PY_MX $METAXCAN/SPrediXcan.py \
   --gwas_file $OUTPUT/processed_summary_imputation/imputed_$TRAIT.txt.gz \
   --snp_column panel_variant_id --effect_allele_column effect_allele \
   --non_effect_allele_column non_effect_allele --zscore_column zscore \
   --model_db_path $DATA/models/eqtl/mashr/mashr_${T}.db \
   --covariance $DATA/models/eqtl/mashr/mashr_${T}.txt.gz \
   --keep_non_rsid --additional_output --model_db_snp_key varID --throw \
   --output_file $OUT > $OUTPUT/logs/spredixcan_${T}.log 2>&1
  echo "done $T"
}
export -f run; export OUTPUT TRAIT PY_MX METAXCAN DATA
ls $DATA/models/eqtl/mashr/mashr_*.db | sed 's#.*/mashr_##; s#\.db$##' | xargs -P $NJOBS -I{} bash -c 'run {}'
