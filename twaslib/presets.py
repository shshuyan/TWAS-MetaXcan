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
    if fmt == "custom":
        return {k: v for k, v in (overrides or {}).items() if v}
    if fmt not in PRESETS:
        raise SystemExit(f"Unknown input.format '{fmt}'. Known: {sorted(PRESETS)} or 'gwasvcf'.")
    cols = dict(PRESETS[fmt])
    for k, v in (overrides or {}).items():
        if v is None:
            cols.pop(k, None)
        else:
            cols[k] = v
    return cols


# Synonyms for auto-detecting columns (case-insensitive). Earlier entries win, so
# GWAS Catalog harmonised `hm_*` columns are preferred over the submitted ones.
SYNONYMS = {
    "chrom": ["hm_chrom", "chromosome", "#chrom", "chrom", "chr", "#chr", "chr_name", "#chromosome"],
    "pos": ["hm_pos", "base_pair_location", "pos", "position", "bp", "bp_hg38", "bp_hg19", "pos_b38", "pos_b37", "genpos"],
    "ea": ["hm_effect_allele", "effect_allele", "alt", "ea", "tested_allele", "allele_b", "allele2", "a1"],
    "nea": ["hm_other_allele", "other_allele", "non_effect_allele", "noneffect_allele", "ref", "nea", "allele_a", "allele1", "a2"],
    "beta": ["hm_beta", "beta", "b", "effect", "effect_size", "log_odds", "beta_meta"],
    "odds_ratio": ["hm_odds_ratio", "odds_ratio", "or"],
    "se": ["standard_error", "se", "sebeta", "stderr", "se_beta", "se_meta"],
    "pval": ["p_value", "pval", "p", "pvalue", "p.value", "p-value", "p_bolt_lmm", "p_dgc"],
    "mlogp": ["mlogp", "neg_log_10_p_value", "neglog10_pval", "log10p", "neglog10_pval_meta"],
    "eaf": ["hm_effect_allele_frequency", "effect_allele_frequency", "af_alt", "eaf", "a1freq", "frq", "freq",
            "af", "maf_effect", "af_controls", "a1_freq"],
    "n": ["n", "n_total", "sample_size", "neff", "obs_ct", "n_obs"],
    "n_cases": ["n_cases", "num_cases", "ncase", "n_case"],
    "variant_id": ["rsid", "rsids", "snp", "markername", "snpid", "id", "rs_id", "variant_id", "variant"],
}
# Allele1/Allele2 mean different things in different tools; flag them for the user to check.
AMBIGUOUS = {"allele1", "allele2", "a1", "a2"}


def autodetect(header):
    """Guess a column mapping from a list of header names. Returns (mapping, warnings)."""
    low = {h.lower(): h for h in header}
    cols, warn = {}, []
    for field, names in SYNONYMS.items():
        for n in names:
            if n in low and low[n] not in cols.values():
                cols[field] = low[n]
                if n in AMBIGUOUS:
                    warn.append(f"'{low[n]}' was used as {field}; check which allele the effect refers to")
                break
    if "beta" in cols and "odds_ratio" in cols:
        cols.pop("odds_ratio")
    for req in ("chrom", "pos", "ea", "nea"):
        if req not in cols:
            warn.append(f"no column found for required field '{req}'")
    if not any(k in cols for k in ("pval", "mlogp", "zscore")):
        warn.append("no p-value column found")
    if not any(k in cols for k in ("beta", "odds_ratio", "zscore")):
        warn.append("no effect-size column found (beta / odds ratio / z)")
    return cols, warn
