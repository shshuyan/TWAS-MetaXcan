"""Convert GWAS summary statistics to and from GWAS-VCF.

GWAS-VCF (Lyon et al. 2021, Genome Biology; MRC-IEU OpenGWAS format) stores one study per
sample column. FORMAT fields used here:
  ES effect size of ALT, SE standard error, LP -log10 p, AF ALT frequency,
  SS sample size, NC number of cases, EZ z-score.
ALT is always the effect allele. Before harmonization REF is just the non-effect allele
reported by the study; the harmonize step aligns alleles to the reference panel.

All sorting, indexing, filtering, extraction and merging is done with bcftools / tabix.
"""
import gzip
import io
import logging
import math
import os
import subprocess

import numpy as np
import pandas as pd
from scipy import stats

from .presets import resolve_columns

FORMAT_FIELDS = ["ES", "SE", "LP", "AF", "SS", "NC", "EZ"]
HEADER_FORMAT = [
    '##FORMAT=<ID=ES,Number=A,Type=Float,Description="Effect size estimate relative to the alternative allele">',
    '##FORMAT=<ID=SE,Number=A,Type=Float,Description="Standard error of effect size estimate">',
    '##FORMAT=<ID=LP,Number=A,Type=Float,Description="-log10 p-value for effect estimate">',
    '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency in the association study">',
    '##FORMAT=<ID=SS,Number=A,Type=Float,Description="Sample size used to estimate genetic effect">',
    '##FORMAT=<ID=NC,Number=A,Type=Float,Description="Number of cases used to estimate genetic effect">',
    '##FORMAT=<ID=EZ,Number=A,Type=Float,Description="Z-score provided if it was used to derive the EFFECT and SE fields">',
]
CONTIGS = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
CHROM_FIX = {"23": "X", "24": "Y", "25": "X", "26": "MT", "M": "MT"}
LP_CAP = 320.0


def bcftools(cfg):
    return cfg.pipe["tools"]["bcftools"]


def run(cmd, **kw):
    logging.info("$ %s", cmd if isinstance(cmd, str) else " ".join(cmd))
    subprocess.run(cmd, check=True, shell=isinstance(cmd, str), executable="/bin/bash" if isinstance(cmd, str) else None, **kw)


def _norm_chrom(s):
    s = s.astype(str).str.replace(r"^chr", "", regex=True, case=False).str.upper()
    return s.replace(CHROM_FIX)


def _fmt(x, spec):
    return "." if x is None or (isinstance(x, float) and not math.isfinite(x)) else format(x, spec)


def standardize_chunk(df, cols, inp):
    """Map one raw chunk to canonical columns and derive missing statistics."""
    out = pd.DataFrame(index=df.index)
    get = lambda k: df[cols[k]] if k in cols and cols[k] in df.columns else None
    for k in ("chrom", "pos", "ea", "nea"):
        if get(k) is None:
            raise SystemExit(f"Required column for '{k}' not found (looked for '{cols.get(k)}'). "
                             f"Columns in file: {list(df.columns)}")
    out["chrom"] = _norm_chrom(get("chrom"))
    out["pos"] = pd.to_numeric(get("pos"), errors="coerce")
    out["ea"] = get("ea").astype(str).str.upper()
    out["nea"] = get("nea").astype(str).str.upper()
    vid = get("variant_id")
    # VCF IDs: missing -> ".", several rsIDs (FinnGen "rs1,rs2") -> first one
    out["id"] = (vid.fillna(".").astype(str).str.split(r"[,;]").str[0].replace({"": ".", "nan": "."})
                 if vid is not None else ".")
    num = lambda k: pd.to_numeric(get(k), errors="coerce") if get(k) is not None else pd.Series(np.nan, index=df.index)

    beta, se, eaf = num("beta"), num("se"), num("eaf")
    if get("beta") is None and get("odds_ratio") is not None:
        beta = np.log(num("odds_ratio"))
    z = num("zscore")
    if get("mlogp") is not None:
        lp = num("mlogp")
    else:
        p = num("pval")
        lp = -np.log10(p.clip(lower=1e-320))
    # Fill gaps: z from beta/se, else from p and the sign of beta
    z = z.where(z.notna(), beta / se)
    z_from_p = stats.norm.isf(np.power(10.0, -lp.clip(upper=300)) / 2) * np.sign(beta)
    z = z.where(z.notna(), z_from_p)
    lp = lp.where(lp.notna(), -(stats.norm.logsf(np.abs(z)) + np.log(2)) / np.log(10))
    se = se.where(se.notna(), np.abs(beta / z))
    if inp.get("eaf_is_nea_freq"):
        eaf = 1 - eaf

    out["ES"], out["SE"], out["LP"], out["AF"], out["EZ"] = beta, se, lp.clip(upper=LP_CAP), eaf, z
    n = num("n")
    if get("n") is None and get("n_cases") is not None and get("n_controls") is not None:
        n = num("n_cases") + num("n_controls")  # e.g. EADB / GWAS Catalog n_cas + n_con
    out["SS"] = n.where(n.notna(), inp.get("sample_size") or np.nan)
    nc = num("n_cases")
    out["NC"] = nc.where(nc.notna(), inp.get("n_cases") or np.nan)

    flag = inp.get("exclude_flag")  # e.g. Pan-UKBB low_confidence_EUR == true
    if flag and flag in df.columns:
        out = out[~df[flag].astype(str).str.lower().isin(["true", "1", "yes"])]
    ok = (out.pos.notna() & out.chrom.isin(CONTIGS) & out.ea.str.fullmatch(r"[ACGT]+")
          & out.nea.str.fullmatch(r"[ACGT]+") & (out.ea != out.nea) & out.LP.notna())
    return out[ok]


def vcf_header(trait_id, build, meta):
    h = ["##fileformat=VCFv4.2", f"##genomeBuild={build}", "##source=twaslib.gwasvcf"]
    h += [f"##contig=<ID={c}>" for c in CONTIGS]
    h += HEADER_FORMAT
    kv = ",".join(f'{k}="{v}"' for k, v in meta.items() if v not in (None, ""))
    h.append(f"##SAMPLE=<ID={trait_id},{kv}>" if kv else f"##SAMPLE=<ID={trait_id}>")
    h.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", trait_id]))
    return "\n".join(h) + "\n"


def _vcf_lines(d):
    fmt = ":".join(FORMAT_FIELDS)
    for r in d.itertuples(index=False):
        sample = ":".join([_fmt(r.ES, ".6g"), _fmt(r.SE, ".6g"), _fmt(r.LP, ".5g"), _fmt(r.AF, ".5g"),
                           _fmt(r.SS, ".0f"), _fmt(r.NC, ".0f"), _fmt(r.EZ, ".5g")])
        yield f"{r.chrom}\t{int(r.pos)}\t{r.id}\t{r.nea}\t{r.ea}\t.\tPASS\t.\t{fmt}\t{sample}\n"


def write_vcf(cfg, chunks, out_path, build, meta):
    """Stream canonical chunks into a sorted, bgzipped, indexed GWAS-VCF."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    tmp = out_path + ".tmp.vcf.gz"
    tmpdir = os.path.dirname(out_path) + "/sort_tmp"
    sort = subprocess.Popen(f"{bcftools(cfg)} sort -m 1G -T {tmpdir} -Oz -o {tmp} -W=tbi",
                            shell=True, stdin=subprocess.PIPE, text=True, executable="/bin/bash")
    sort.stdin.write(vcf_header(cfg.trait_id, build, meta))
    n = 0
    for d in chunks:
        sort.stdin.writelines(_vcf_lines(d))
        n += len(d)
    sort.stdin.close()
    if sort.wait() != 0:
        raise SystemExit("bcftools sort failed")
    os.replace(tmp, out_path)
    os.replace(tmp + ".tbi", out_path + ".tbi")
    logging.info("Wrote %d variants to %s", n, out_path)
    return n


def raw_to_vcf(cfg):
    """Step: raw summary statistics (any preset) -> GWAS-VCF in the input build."""
    inp = cfg.trait["input"]
    path = inp["path"]
    if inp.get("format") == "gwasvcf":
        # Already GWAS-VCF (e.g. from OpenGWAS): rename the sample to the trait id, sort and index
        os.makedirs(os.path.dirname(cfg.raw_vcf), exist_ok=True)
        run(f"{bcftools(cfg)} view {path} | sed '/^#CHROM/s/\\t[^\\t]*$/\\t{cfg.trait_id}/' "
            f"| sed 's/^chr//' | {bcftools(cfg)} sort -Oz -o {cfg.raw_vcf} -W=tbi")
        return
    cols = resolve_columns(inp.get("format", "generic"), inp.get("columns"))
    sep = inp.get("separator", "\t")
    usecols = list(dict.fromkeys(list(cols.values()) + ([inp["exclude_flag"]] if inp.get("exclude_flag") else [])))
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        header = fh.readline().rstrip("\n")
    present = pd.read_csv(io.StringIO(header + "\n"), sep=sep, engine="python").columns
    usecols = [c for c in usecols if c in present]
    reader = pd.read_csv(path, sep=sep, usecols=usecols, chunksize=1_000_000, dtype=str,
                         na_values=["NA", "nan", ".", ""])
    meta = {"StudyType": "CaseControl" if inp.get("n_cases") else "Continuous",
            "TotalCases": inp.get("n_cases"), "TotalControls": inp.get("n_controls"),
            "SampleSize": inp.get("sample_size"), "Ancestry": cfg.trait.get("ancestry"),
            "Source": inp.get("source"), "Description": cfg.trait.get("description", "").replace('"', "'")}
    write_vcf(cfg, (standardize_chunk(c, cols, inp) for c in reader), cfg.raw_vcf, cfg.build, meta)


EXPORT_COLUMNS = ["variant_id", "chromosome", "position", "non_effect_allele", "effect_allele",
                  "frequency", "pvalue", "zscore", "effect_size", "standard_error", "sample_size", "n_cases"]


def vcf_to_tsv(cfg, vcf, out_tsv, include=None, regions=None):
    """Step: GWAS-VCF -> the tab-separated layout expected by gwas_parsing.py.

    `include` is a bcftools -i expression (e.g. 'FORMAT/AF>=0.01 && FORMAT/AF<=0.99'),
    `regions` a bcftools -r string, both applied with bcftools before export.
    """
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    q = f"{bcftools(cfg)} query"
    if include:
        q += f" -i '{include}'"
    if regions:
        q += f" -r {regions}"
    q += r" -f '%ID\t%CHROM\t%POS\t%REF\t%ALT[\t%AF\t%LP\t%EZ\t%ES\t%SE\t%SS\t%NC]\n' " + vcf
    inp = cfg.trait.get("input", {})
    proc = subprocess.Popen(q, shell=True, stdout=subprocess.PIPE, text=True, executable="/bin/bash")
    names = ["variant_id", "chromosome", "position", "non_effect_allele", "effect_allele",
             "frequency", "LP", "zscore", "effect_size", "standard_error", "sample_size", "n_cases"]
    tmp = out_tsv + ".tmp"
    n = 0
    with gzip.open(tmp, "wt") as out:
        out.write("\t".join(EXPORT_COLUMNS) + "\n")
        for d in pd.read_csv(proc.stdout, sep="\t", names=names, na_values=["."], chunksize=1_000_000,
                             dtype={"chromosome": str, "variant_id": str}):
            d["pvalue"] = np.power(10.0, -d.LP)
            # Very small p underflows; use the z-score for those
            tiny = d.pvalue <= 1e-300
            d.loc[tiny, "pvalue"] = 1e-300
            if inp.get("sample_size"):
                d["sample_size"] = d.sample_size.fillna(inp["sample_size"])
            if inp.get("n_cases"):
                d["n_cases"] = d.n_cases.fillna(inp["n_cases"])
            d["variant_id"] = d.variant_id.where(d.variant_id.notna(),
                                                 "chr" + d.chromosome + "_" + d.position.astype(str) + "_"
                                                 + d.non_effect_allele + "_" + d.effect_allele)
            d[EXPORT_COLUMNS].to_csv(out, sep="\t", index=False, header=False, na_rep="NA")
            n += len(d)
    if proc.wait() != 0:
        raise SystemExit("bcftools query failed")
    os.replace(tmp, out_tsv)
    logging.info("Exported %d variants to %s", n, out_tsv)


def harmonized_to_vcf(cfg):
    """Write the harmonized (hg38, reference-aligned) GWAS as GWAS-VCF for cross-trait merging."""
    def chunks():
        for d in pd.read_csv(cfg.harmonized, sep="\t", chunksize=1_000_000, na_values=["NA"],
                             dtype={"chromosome": str}):
            d = d[d.panel_variant_id.notna()]
            # panel_variant_id = chr{c}_{pos}_{ref}_{alt}_b38; effect allele is the panel alt
            parts = d.panel_variant_id.str.split("_", expand=True)
            out = pd.DataFrame({"chrom": parts[0].str.replace("chr", "", regex=False), "pos": parts[1].astype(int),
                                "nea": parts[2], "ea": parts[3], "id": d.variant_id.astype(str),
                                "ES": d.effect_size, "SE": d.standard_error,
                                "LP": -np.log10(d.pvalue.clip(lower=1e-320)), "AF": d.frequency,
                                "SS": d.sample_size, "NC": d.n_cases, "EZ": d.zscore})
            yield out
    meta = {"Ancestry": cfg.trait.get("ancestry"), "Harmonized": "GTEx v8 / 1000G panel, hg38"}
    write_vcf(cfg, chunks(), cfg.harmonized_vcf, "hg38", meta)
