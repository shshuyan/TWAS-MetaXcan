"""Pipeline steps. Each step is skipped when its output already exists (use --force to redo)."""
import glob
import logging
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

from . import gwasvcf
from .config import LIFTOVER_CHAINS

STEPS = ["download", "standardize", "export", "harmonize", "impute", "postprocess",
         "spredixcan", "smultixcan", "summarize"]


def _done(path, force):
    if os.path.exists(path) and os.path.getsize(path) > 0 and not force:
        logging.info("  exists, skipping: %s", path)
        return True
    return False


def _sh(cmd, log=None):
    logging.info("$ %s", " ".join(cmd) if isinstance(cmd, list) else cmd)
    if log:
        os.makedirs(os.path.dirname(log), exist_ok=True)
        with open(log, "w") as fh:
            r = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
        if r.returncode != 0:
            tail = open(log).read()[-3000:]
            raise SystemExit(f"Command failed (log: {log}):\n{tail}")
    else:
        subprocess.run(cmd, check=True)


def step_download(cfg, force=False):
    inp = cfg.trait["input"]
    if os.path.exists(inp["path"]):
        logging.info("  input present: %s", inp["path"])
        return
    url = inp.get("url")
    if not url:
        raise SystemExit(f"Input file {inp['path']} is missing and input.url is not set. "
                         f"Download it and put it at that path (see the comments in the trait config).")
    os.makedirs(os.path.dirname(inp["path"]), exist_ok=True)
    tmp = inp["path"] + ".part"
    _sh(["curl", "-fSL", "--retry", "5", "-o", tmp, url])
    os.replace(tmp, inp["path"])


def step_standardize(cfg, force=False):
    if not _done(cfg.raw_vcf, force):
        gwasvcf.raw_to_vcf(cfg)


def step_export(cfg, force=False):
    if not _done(cfg.export_tsv, force):
        gwasvcf.vcf_to_tsv(cfg, cfg.raw_vcf, cfg.export_tsv,
                           include=cfg.trait.get("vcf_filter"),
                           # reference panel and models are autosomal: keep only the trait's chromosomes
                           regions=cfg.trait.get("vcf_regions")
                           or ",".join(str(c) for c in cfg.trait.get("chromosomes", range(1, 23))))


def step_harmonize(cfg, force=False):
    out = cfg.harmonized
    if not _done(out, force):
        if force and os.path.exists(out):
            os.remove(out)  # gwas_parsing.py refuses to overwrite
        os.makedirs(os.path.dirname(out), exist_ok=True)
        ref, py = cfg.ref, cfg.pipe["python"]["gwasimputation"]
        cmd = [py, os.path.join(cfg.pipe["paths"]["gwas_tools"], "gwas_parsing.py"),
               "-gwas_file", cfg.export_tsv,
               "-snp_reference_metadata", ref["variant_metadata_txt"], "METADATA",
               "--chromosome_format"]
        for c in gwasvcf.EXPORT_COLUMNS:
            cmd += ["-output_column_map", c, c]
        if cfg.build != "hg38":
            cmd += ["-liftover", os.path.join(ref["liftover_dir"], LIFTOVER_CHAINS[cfg.build])]
        cmd += ["-output_order", "variant_id", "panel_variant_id", "chromosome", "position",
                "effect_allele", "non_effect_allele", "frequency", "pvalue", "zscore",
                "effect_size", "standard_error", "sample_size", "n_cases", "-output", out]
        _sh(cmd, log=cfg.p("logs", "harmonize.log"))
    if cfg.pipe.get("write_harmonized_vcf", True) and not _done(cfg.harmonized_vcf, force):
        gwasvcf.harmonized_to_vcf(cfg)


def _impute_job(cfg, chrom, sb, force):
    imp = cfg.pipe["imputation"]
    out = cfg.p("03_summary_imputation", f"{cfg.trait_id}_chr{chrom}_sb{sb}_by_region.txt.gz")
    if _done(out, force):
        return out
    ref = cfg.ref
    # write into a subfolder so a half-written file never matches the postprocess pattern
    tmp = cfg.p("03_summary_imputation", "partial", os.path.basename(out))
    os.makedirs(os.path.dirname(tmp), exist_ok=True)
    if os.path.exists(tmp):
        os.remove(tmp)
    cmd = [cfg.pipe["python"]["gwasimputation"],
           os.path.join(cfg.pipe["paths"]["gwas_tools"], "gwas_summary_imputation.py"),
           "-by_region_file", ref["ld_blocks"],
           "-gwas_file", cfg.harmonized,
           "-parquet_genotype", ref["genotype_parquet"].format(chrom=chrom),
           "-parquet_genotype_metadata", ref["variant_metadata_parquet"],
           "-window", str(imp["window"]), "-parsimony", "7", "-chromosome", str(chrom),
           "-regularization", str(imp["regularization"]), "-frequency_filter", str(imp["frequency_filter"]),
           "-sub_batches", str(imp["sub_batches"]), "-sub_batch", str(sb),
           "--standardise_dosages", "-output", tmp]
    _sh(cmd, log=cfg.p("logs", "impute", f"chr{chrom}_sb{sb}.log"))
    os.replace(tmp, out)
    return out


def step_impute(cfg, force=False):
    if cfg.trait.get("skip_imputation"):
        logging.info("  skip_imputation: true, S-PrediXcan will use the harmonized GWAS")
        return
    if _done(cfg.imputed, force):
        return
    os.makedirs(cfg.p("03_summary_imputation"), exist_ok=True)
    nb = cfg.pipe["imputation"]["sub_batches"]
    chroms = cfg.trait.get("chromosomes", list(range(1, 23)))
    jobs = [(c, b) for c in chroms for b in range(nb)]
    logging.info("  %d imputation jobs, %d in parallel", len(jobs), cfg.pipe["njobs"])
    with ThreadPoolExecutor(cfg.pipe["njobs"]) as ex:
        futs = {ex.submit(_impute_job, cfg, c, b, force): (c, b) for c, b in jobs}
        for i, f in enumerate(as_completed(futs), 1):
            f.result()
            logging.info("  imputation %d/%d done (chr%s sb%s)", i, len(jobs), *futs[f])


def step_postprocess(cfg, force=False):
    if cfg.trait.get("skip_imputation") or _done(cfg.imputed, force):
        return
    os.makedirs(os.path.dirname(cfg.imputed), exist_ok=True)
    _sh([cfg.pipe["python"]["gwasimputation"],
         os.path.join(cfg.pipe["paths"]["gwas_tools"], "gwas_summary_imputation_postprocess.py"),
         "-gwas_file", cfg.harmonized, "-folder", cfg.p("03_summary_imputation"),
         "-pattern", f"{re.escape(cfg.trait_id)}_chr.*", "-parsimony", "7",
         "-output", cfg.imputed], log=cfg.p("logs", "postprocess.log"))


def _tissues(cfg, mtype):
    mdir = cfg.ref["models"][mtype]["dir"]
    pat = cfg.ref["models"][mtype].get("prefix", "mashr_")
    have = sorted(re.sub(rf"^{pat}|\.db$", "", os.path.basename(f)) for f in glob.glob(f"{mdir}/{pat}*.db"))
    want = cfg.trait.get("tissues", "all")
    if want in (None, "all"):
        return have
    missing = set(want) - set(have)
    if missing:
        raise SystemExit(f"Tissues without {mtype} models: {sorted(missing)}")
    return list(want)


def _spredixcan_job(cfg, mtype, tissue, force):
    m = cfg.ref["models"][mtype]
    pre = m.get("prefix", "mashr_")
    out = os.path.join(cfg.spredixcan_dir(mtype), f"{cfg.trait_id}__PM__{tissue}.csv")
    if _done(out, force):
        return
    _sh([cfg.pipe["python"]["metaxcan"], os.path.join(cfg.pipe["paths"]["metaxcan"], "SPrediXcan.py"),
         "--gwas_file", cfg.twas_input,
         "--snp_column", "panel_variant_id", "--effect_allele_column", "effect_allele",
         "--non_effect_allele_column", "non_effect_allele", "--zscore_column", "zscore",
         "--model_db_path", os.path.join(m["dir"], f"{pre}{tissue}.db"),
         "--covariance", os.path.join(m["dir"], f"{pre}{tissue}.txt.gz"),
         "--keep_non_rsid", "--additional_output", "--model_db_snp_key", "varID", "--throw",
         "--output_file", out], log=cfg.p("logs", "spredixcan", mtype, f"{tissue}.log"))


def step_spredixcan(cfg, force=False):
    for mtype in cfg.trait.get("model_types", ["eqtl"]):
        os.makedirs(cfg.spredixcan_dir(mtype), exist_ok=True)
        tissues = _tissues(cfg, mtype)
        logging.info("  S-PrediXcan %s: %d tissues", mtype, len(tissues))
        with ThreadPoolExecutor(cfg.pipe["njobs"]) as ex:
            for f in as_completed([ex.submit(_spredixcan_job, cfg, mtype, t, force) for t in tissues]):
                f.result()


def step_smultixcan(cfg, force=False):
    for mtype in cfg.trait.get("model_types", ["eqtl"]):
        m = cfg.ref["models"][mtype]
        out = cfg.smultixcan_out(mtype)
        if _done(out, force):
            continue
        os.makedirs(os.path.dirname(out), exist_ok=True)
        pre = m.get("prefix", "mashr_")
        _sh([cfg.pipe["python"]["metaxcan"], os.path.join(cfg.pipe["paths"]["metaxcan"], "SMulTiXcan.py"),
             "--models_folder", m["dir"], "--models_name_pattern", f"{pre}(.*).db",
             "--snp_covariance", m["smultixcan_covariance"],
             "--metaxcan_folder", cfg.spredixcan_dir(mtype) + "/",
             "--metaxcan_filter", f"{re.escape(cfg.trait_id)}__PM__(.*).csv",
             "--metaxcan_file_name_parse_pattern", "(.*)__PM__(.*).csv",
             "--gwas_file", cfg.twas_input,
             "--snp_column", "panel_variant_id", "--effect_allele_column", "effect_allele",
             "--non_effect_allele_column", "non_effect_allele", "--zscore_column", "zscore",
             "--keep_non_rsid", "--model_db_snp_key", "varID",
             "--cutoff_condition_number", "30", "--verbosity", "7", "--throw",
             "--output", out], log=cfg.p("logs", f"smultixcan_{mtype}.log"))


def step_summarize(cfg, force=False):
    from .summarize import summarize
    summarize(cfg)


def run_steps(cfg, steps, force=False):
    for s in steps:
        logging.info("[%s] step: %s", cfg.trait_id, s)
        globals()[f"step_{s}"](cfg, force)
