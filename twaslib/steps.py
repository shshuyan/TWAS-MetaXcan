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
         "spredixcan", "smultixcan", "summarize", "candidate_recheck"]


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
    _sh(["curl", "-fsSL", "--show-error", "--retry", "5", "-o", tmp, url])
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
           "-gwas_file", _chrom_gwas(cfg, chrom),
           "-parquet_genotype", ref["genotype_parquet"].format(chrom=chrom),
           "-parquet_genotype_metadata", ref["variant_metadata_parquet"],
           "-window", str(imp["window"]), "-parsimony", "7", "-chromosome", str(chrom),
           "-regularization", str(imp["regularization"]), "-frequency_filter", str(imp["frequency_filter"]),
           "-sub_batches", str(imp["sub_batches"]), "-sub_batch", str(sb),
           "--standardise_dosages", "-output", tmp]
    _sh(cmd, log=cfg.p("logs", "impute", f"chr{chrom}_sb{sb}.log"))
    os.replace(tmp, out)
    return out


def _chrom_gwas(cfg, chrom):
    return cfg.p("03_summary_imputation", "by_chrom", f"chr{chrom}.txt.gz")


def _split_by_chrom(cfg, chroms):
    """Write one harmonized file per chromosome so each imputation job loads only its chromosome
    (less to parse per job; peak memory is dominated by the reference genotypes)."""
    import gzip
    import pandas as pd
    if all(os.path.exists(_chrom_gwas(cfg, c)) for c in chroms):
        return
    os.makedirs(os.path.dirname(_chrom_gwas(cfg, chroms[0])), exist_ok=True)
    logging.info("  splitting harmonized GWAS by chromosome")
    handles = {}
    try:
        for d in pd.read_csv(cfg.harmonized, sep="\t", chunksize=1_000_000, dtype=str):
            for c, g in d.groupby("chromosome"):
                key = c.replace("chr", "")
                if key not in {str(x) for x in chroms}:
                    continue
                if key not in handles:
                    handles[key] = gzip.open(_chrom_gwas(cfg, key) + ".tmp", "wt", compresslevel=3)
                    g.iloc[:0].to_csv(handles[key], sep="\t", index=False)
                g.to_csv(handles[key], sep="\t", index=False, header=False, na_rep="NA")
    finally:
        for h in handles.values():
            h.close()
    for key in handles:
        os.replace(_chrom_gwas(cfg, key) + ".tmp", _chrom_gwas(cfg, key))
    for c in chroms:  # chromosomes with no variants still need a (header-only) file
        if not os.path.exists(_chrom_gwas(cfg, c)):
            with gzip.open(_chrom_gwas(cfg, c), "wt") as fh:
                fh.write("\t".join(pd.read_csv(cfg.harmonized, sep="\t", nrows=0).columns) + "\n")


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
    _split_by_chrom(cfg, list(chroms))
    logging.info("  %d imputation jobs, %d in parallel", len(jobs), cfg.pipe["njobs"])
    failed = _run_pool(cfg, jobs, lambda c, b: _impute_job(cfg, c, b, force), "imputation")
    if failed:
        # Most failures are out-of-memory kills from running in parallel: retry one at a time
        logging.info("  retrying %d failed imputation jobs one at a time", len(failed))
        failed = _run_pool(cfg, failed, lambda c, b: _impute_job(cfg, c, b, force), "imputation retry", njobs=1)
    if failed:
        raise SystemExit(f"{len(failed)} imputation jobs failed: {failed[:5]} (see {cfg.p('logs', 'impute')})")


def _run_pool(cfg, jobs, fn, label, njobs=None):
    """Run fn(*job) for each job in parallel; return the jobs that failed."""
    njobs = njobs or cfg.pipe["njobs"]
    failed = []
    with ThreadPoolExecutor(njobs) as ex:
        futs = {ex.submit(fn, *j): j for j in jobs}
        for i, f in enumerate(as_completed(futs), 1):
            try:
                f.result()
                logging.info("  %s %d/%d done %s", label, i, len(jobs), futs[f])
            except (SystemExit, Exception) as e:  # noqa: BLE001
                logging.warning("  %s %s failed: %s", label, futs[f], str(e).splitlines()[0][:200])
                failed.append(futs[f])
    return failed


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
        jobs = [(mtype, t) for t in tissues]
        failed = _run_pool(cfg, jobs, lambda m, t: _spredixcan_job(cfg, m, t, force), f"S-PrediXcan {mtype}")
        if failed:
            failed = _run_pool(cfg, failed, lambda m, t: _spredixcan_job(cfg, m, t, force),
                               f"S-PrediXcan {mtype} retry", njobs=1)
        if failed:
            raise SystemExit(f"S-PrediXcan failed for {failed} (see {cfg.p('logs', 'spredixcan')})")


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


def step_candidate_recheck(cfg, force=False):
    """Opt-in post-screen evidence checks. Never silently replace a genetic source."""
    setting = cfg.trait.get("candidate_recheck")
    if not setting or not setting.get("enabled", False):
        logging.info("  candidate_recheck not enabled")
        return
    from .candidate_recheck import run_recheck
    from .recheck_inputs import load_settings, sha256
    evidence_config = setting["config"]
    variables = {"trait_id": cfg.trait_id, "results": cfg.results, "work": cfg.work}
    evidence, resolved_path = load_settings(evidence_config, variables)
    out = evidence["output"]
    if not os.path.isabs(out):
        out = os.path.join(str(resolved_path.parent), out)
    marker = os.path.join(out, "summary.json")
    if os.path.exists(marker):
        # Reusing stale evidence after input/config changes would invalidate review.
        import json
        manifest_path = os.path.join(out, "input_manifest.json")
        if not os.path.exists(manifest_path):
            raise SystemExit("Recheck completion marker lacks input manifest: " + out)
        with open(manifest_path) as fh:
            manifest = json.load(fh)
        changed = [r["path"] for r in manifest if not os.path.isfile(r["path"]) or sha256(r["path"]) != r["sha256"]]
        with open(os.path.join(out, "resolved_config.json")) as fh:
            old_config = json.load(fh)
        if changed or old_config != evidence or force:
            raise SystemExit("Recheck inputs/config/code changed or --force requested. Set a NEW evidence output directory; old evidence is preserved.")
        logging.info("  unchanged completed candidate recheck: %s", out)
        return
    run_recheck(evidence_config, variables=variables)


def run_steps(cfg, steps, force=False):
    for s in steps:
        logging.info("[%s] step: %s", cfg.trait_id, s)
        globals()[f"step_{s}"](cfg, force)
