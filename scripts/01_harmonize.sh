#!/bin/bash
# Step 4: harmonize GWAS, liftover hg19 -> hg38
set -euo pipefail; source "$(dirname "$0")/config.sh"
mkdir -p $OUTPUT/harmonized_gwas
$PY_IMP $GWAS_TOOLS/gwas_parsing.py \
 -gwas_file $DATA/gwas/cad.add.160614.website.txt.gz \
 -liftover $DATA/liftover/hg19ToHg38.over.chain.gz \
 -snp_reference_metadata $DATA/reference_panel_1000G/variant_metadata.txt.gz METADATA \
 -output_column_map markername variant_id \
 -output_column_map noneffect_allele non_effect_allele \
 -output_column_map effect_allele effect_allele \
 -output_column_map beta effect_size \
 -output_column_map p_dgc pvalue \
 -output_column_map chr chromosome \
 --chromosome_format \
 -output_column_map bp_hg19 position \
 -output_column_map effect_allele_freq frequency \
 --insert_value sample_size 184305 --insert_value n_cases 60801 \
 -output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
 -output $OUTPUT/harmonized_gwas/$TRAIT.txt.gz
