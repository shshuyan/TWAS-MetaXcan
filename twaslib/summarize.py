"""Per-trait result tables, and cross-trait comparison."""
import glob
import logging
import os
import re

import pandas as pd


def _load_spredixcan(cfg, mtype):
    frames = []
    for f in sorted(glob.glob(os.path.join(cfg.spredixcan_dir(mtype), f"{cfg.trait_id}__PM__*.csv"))):
        d = pd.read_csv(f)
        d.insert(0, "tissue", re.search(r"__PM__(.*)\.csv$", f).group(1))
        frames.append(d)
    if not frames:
        raise SystemExit(f"No S-PrediXcan output in {cfg.spredixcan_dir(mtype)}")
    return pd.concat(frames, ignore_index=True).dropna(subset=["pvalue"]).sort_values("pvalue")


def summarize(cfg):
    os.makedirs(cfg.results, exist_ok=True)
    alpha = cfg.pipe.get("alpha", 0.05)
    lines = [f"# TWAS results: {cfg.trait_id}", "", cfg.trait.get("description", ""), "",
             f"- Input: `{os.path.basename(cfg.trait['input']['path'])}` "
             f"(format `{cfg.trait['input'].get('format')}`, build {cfg.build}, ancestry {cfg.trait.get('ancestry', 'EUR')})",
             f"- Imputation: {'skipped' if cfg.trait.get('skip_imputation') else 'summary-gwas-imputation, 1000G EUR'}",
             ""]
    for mtype in cfg.trait.get("model_types", ["eqtl"]):
        unit = "gene" if mtype == "eqtl" else "intron"
        sp = _load_spredixcan(cfg, mtype)
        # Same layout as the MetaXcan tutorial's "Gene level results" (one tissue, all S-PrediXcan columns)
        for t in cfg.pipe.get("example_tissues", ["Whole_Blood"]):
            src = os.path.join(cfg.spredixcan_dir(mtype), f"{cfg.trait_id}__PM__{t}.csv")
            if os.path.exists(src):
                pd.read_csv(src).sort_values("pvalue").to_csv(
                    os.path.join(cfg.results, f"spredixcan_{mtype}_{t}.csv"), index=False)
        thr = alpha / len(sp)
        sig = sp[sp.pvalue < thr]
        sp.to_csv(os.path.join(cfg.results, f"spredixcan_{mtype}_all.tsv.gz"), sep="\t", index=False)
        sig.to_csv(os.path.join(cfg.results, f"spredixcan_{mtype}_significant.tsv"), sep="\t", index=False)
        name = "gene_name" if "gene_name" in sp else "gene"
        top = (sig.groupby(["gene", name] if name != "gene" else ["gene"])
                  .agg(n_tissues=("tissue", "nunique"), best_p=("pvalue", "min"))
                  .reset_index().sort_values("best_p"))
        best_t = sig.drop_duplicates("gene").set_index("gene").tissue
        top["best_tissue"] = top.gene.map(best_t)

        sm_path = cfg.smultixcan_out(mtype)
        sm = pd.read_csv(sm_path, sep="\t").dropna(subset=["pvalue"]).sort_values("pvalue")
        sm_thr = alpha / len(sm)
        sm_sig = sm[sm.pvalue < sm_thr]
        sm.to_csv(os.path.join(cfg.results, f"smultixcan_{mtype}_all.tsv.gz"), sep="\t", index=False)
        sm_sig.to_csv(os.path.join(cfg.results, f"smultixcan_{mtype}_significant.tsv"), sep="\t", index=False)

        lines += [f"## {mtype}: S-PrediXcan (per tissue)", "",
                  f"- {sp.tissue.nunique()} tissues, {len(sp):,} {unit}-tissue tests, Bonferroni p < {thr:.2e}",
                  f"- Significant pairs: {len(sig):,} ({len(top)} unique {unit}s)", "",
                  f"| {unit} | tissues significant | best tissue | best p |", "|---|---|---|---|"]
        for r in top.head(30).itertuples():
            lines.append(f"| {getattr(r, name)} | {r.n_tissues} | {r.best_tissue} | {r.best_p:.2e} |")
        sm_name = "gene_name" if "gene_name" in sm else "gene"
        lines += ["", f"## {mtype}: S-MultiXcan (joint across tissues)", "",
                  f"- {len(sm):,} {unit}s tested, Bonferroni p < {sm_thr:.2e}; significant: {len(sm_sig)}", "",
                  f"| {unit} | p | tissues used | best single-tissue p |", "|---|---|---|---|"]
        for r in sm.head(30).itertuples():
            lines.append(f"| {getattr(r, sm_name)} | {r.pvalue:.2e} | {int(r.n)} | {r.p_i_best:.2e} |")
        lines.append("")
    with open(os.path.join(cfg.results, "SUMMARY.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    logging.info("Wrote %s", os.path.join(cfg.results, "SUMMARY.md"))


def compare(cfgs, out_dir, mtype="eqtl"):
    """Cross-trait gene table from S-MultiXcan, plus a merged multi-trait GWAS-VCF."""
    from . import gwasvcf
    os.makedirs(out_dir, exist_ok=True)
    alpha = cfgs[0].pipe.get("alpha", 0.05)
    wide = None
    for c in cfgs:
        sm = pd.read_csv(c.smultixcan_out(mtype), sep="\t").dropna(subset=["pvalue"])
        sm = sm[["gene", "gene_name", "pvalue"]].rename(columns={"pvalue": c.trait_id})
        wide = sm if wide is None else wide.merge(sm, on=["gene", "gene_name"], how="outer")
    traits = [c.trait_id for c in cfgs]
    sig = pd.DataFrame({t: wide[t] < alpha / wide[t].notna().sum() for t in traits})
    wide["n_traits_significant"] = sig.sum(axis=1)
    wide = wide.sort_values(["n_traits_significant"] + traits, ascending=[False] + [True] * len(traits))
    wide.to_csv(os.path.join(out_dir, f"smultixcan_{mtype}_cross_trait.tsv"), sep="\t", index=False)
    shared = wide[wide.n_traits_significant >= 2]
    shared.to_csv(os.path.join(out_dir, f"smultixcan_{mtype}_shared_genes.tsv"), sep="\t", index=False)

    vcfs = [c.harmonized_vcf for c in cfgs if os.path.exists(c.harmonized_vcf)]
    if len(vcfs) >= 2:
        merged = os.path.join(out_dir, "harmonized_multi_trait.hg38.vcf.gz")
        gwasvcf.run(f"{gwasvcf.bcftools(cfgs[0])} merge -m none -Oz -o {merged} -W=tbi " + " ".join(vcfs))
    logging.info("Cross-trait: %d genes significant in >= 2 traits", len(shared))
    return shared
