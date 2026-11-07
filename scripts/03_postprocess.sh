#!/bin/bash
# Step 6: merge the 220 imputation batches
set -euo pipefail; source "$(dirname "$0")/config.sh"
mkdir -p $OUTPUT/processed_summary_imputation
$PY_IMP $GWAS_TOOLS/gwas_summary_imputation_postprocess.py \
 -gwas_file $OUTPUT/harmonized_gwas/$TRAIT.txt.gz \
 -folder $OUTPUT/summary_imputation \
 -pattern "$TRAIT.*" -parsimony 7 \
 -output $OUTPUT/processed_summary_imputation/imputed_$TRAIT.txt.gz
