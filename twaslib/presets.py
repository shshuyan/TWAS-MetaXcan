"""Column presets for common GWAS summary-statistics layouts.

Each preset maps our canonical field names to the column names in the raw file.
Canonical fields:
  variant_id  rsID or any variant label (optional)
  chrom, pos  chromosome and base-pair position (required)
  ea, nea     effect allele and other allele (required)
  eaf         effect-allele frequency
  beta        effect size (log-odds for case-control)
  odds_ratio  used to derive beta when beta is absent
  se          standard error
  pval        p-value
  mlogp       -log10 p-value (used instead of pval when present)
  zscore      z-score
  n           per-variant sample size
  n_cases     per-variant number of cases

A trait config may give `input.columns` to override or extend any preset.
"""

PRESETS = {
    # CARDIoGRAMplusC4D 1000G additive meta-analysis (the MetaXcan tutorial file)
    "cardiogram": {
        "variant_id": "markername", "chrom": "chr", "pos": "bp_hg19",
        "ea": "effect_allele", "nea": "noneffect_allele", "eaf": "effect_allele_freq",
        "beta": "beta", "se": "se_dgc", "pval": "p_dgc",
    },
    # FinnGen public release files (R5+), build hg38, effect allele = alt
    "finngen": {
        "variant_id": "rsids", "chrom": "#chrom", "pos": "pos",
        "ea": "alt", "nea": "ref", "eaf": "af_alt",
        "beta": "beta", "se": "sebeta", "pval": "pval", "mlogp": "mlogp",
    },
    # GWAS Catalog harmonised summary statistics (GWAS-SSF, *.h.tsv.gz), build hg38
    "gwas_catalog": {
        "variant_id": "rsid", "chrom": "chromosome", "pos": "base_pair_location",
        "ea": "effect_allele", "nea": "other_allele", "eaf": "effect_allele_frequency",
        "beta": "beta", "odds_ratio": "odds_ratio", "se": "standard_error",
        "pval": "p_value", "n": "n",
    },
    # Pan-UK Biobank per-ancestry files (e.g. *_EAS columns are renamed via input.columns)
    "pan_ukbb": {
        "chrom": "chr", "pos": "pos", "ea": "alt", "nea": "ref",
        "eaf": "af_meta", "beta": "beta_meta", "se": "se_meta", "mlogp": "neglog10_pval_meta",
    },
    # BioBank Japan / pheweb.jp downloads
    "bbj": {
        "variant_id": "SNPID", "chrom": "CHR", "pos": "POS",
        "ea": "Allele2", "nea": "Allele1", "eaf": "AF_Allele2",
        "beta": "BETA", "se": "SE", "pval": "p.value", "n": "N",
    },
    # Plain/generic layout; anything else goes in input.columns
    "generic": {
        "variant_id": "SNP", "chrom": "CHR", "pos": "BP",
        "ea": "A1", "nea": "A2", "eaf": "FRQ",
        "beta": "BETA", "se": "SE", "pval": "P", "n": "N",
    },
}


def resolve_columns(fmt, overrides=None):
    if fmt not in PRESETS:
        raise SystemExit(f"Unknown input.format '{fmt}'. Known: {sorted(PRESETS)} or 'gwasvcf'.")
    cols = dict(PRESETS[fmt])
    for k, v in (overrides or {}).items():
        if v is None:
            cols.pop(k, None)
        else:
            cols[k] = v
    return cols
