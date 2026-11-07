#!/bin/bash
# Step 8: S-MultiXcan across all tissues
set -euo pipefail; source "$(dirname "$0")/config.sh"
mkdir -p $OUTPUT/smultixcan/eqtl
$PY_MX $METAXCAN/SMulTiXcan.py \
 --models_folder $DATA/models/eqtl/mashr \
 --models_name_pattern "mashr_(.*).db" \
 --snp_covariance $DATA/models/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz \
 --metaxcan_folder $OUTPUT/spredixcan/eqtl/ \
 --metaxcan_filter "${TRAIT}__PM__(.*).csv" \
 --metaxcan_file_name_parse_pattern "(.*)__PM__(.*).csv" \
 --gwas_file $OUTPUT/processed_summary_imputation/imputed_$TRAIT.txt.gz \
 --snp_column panel_variant_id --effect_allele_column effect_allele \
 --non_effect_allele_column non_effect_allele --zscore_column zscore --keep_non_rsid \
 --model_db_snp_key varID --cutoff_condition_number 30 --verbosity 7 --throw \
 --output $OUTPUT/smultixcan/eqtl/${TRAIT}_smultixcan.txt
