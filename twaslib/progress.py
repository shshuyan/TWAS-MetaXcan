"""Work out a trait's progress from the files on disk (works whoever started the run)."""
import glob
import os


def _exists(p):
    return os.path.exists(p) and os.path.getsize(p) > 0


def progress(cfg):
    steps = []
    inp = cfg.trait["input"]
    part = inp["path"] + ".part"
    dl = {"name": "download", "done": _exists(inp["path"])}
    if not dl["done"] and os.path.exists(part):
        dl["detail"] = f"{os.path.getsize(part) / 1e6:.0f} MB"
    steps.append(dl)
    steps.append({"name": "standardize", "done": _exists(cfg.raw_vcf)})
    steps.append({"name": "export", "done": _exists(cfg.export_tsv)})
    steps.append({"name": "harmonize", "done": _exists(cfg.harmonized)})
    skip = cfg.trait.get("skip_imputation")
    nchr = len(cfg.trait.get("chromosomes", range(1, 23)))
    total = nchr * cfg.pipe["imputation"]["sub_batches"]
    n_imp = len(glob.glob(cfg.p("03_summary_imputation", "*_by_region.txt.gz")))
    steps.append({"name": "impute", "done": skip or n_imp >= total or _exists(cfg.imputed),
                  "detail": "skipped" if skip else f"{n_imp}/{total}"})
    steps.append({"name": "postprocess", "done": skip or _exists(cfg.imputed), "detail": "skipped" if skip else None})
    for m in cfg.trait.get("model_types", ["eqtl"]):
        mdir = cfg.ref["models"][m]["dir"]
        n_t = len(glob.glob(os.path.join(mdir, "*.db")))
        want = cfg.trait.get("tissues", "all")
        n_t = n_t if want in (None, "all") else len(want)
        n_done = len(glob.glob(os.path.join(cfg.spredixcan_dir(m), "*.csv")))
        steps.append({"name": f"spredixcan ({m})", "done": n_done >= n_t and n_t > 0, "detail": f"{n_done}/{n_t}"})
    for m in cfg.trait.get("model_types", ["eqtl"]):
        steps.append({"name": f"smultixcan ({m})", "done": _exists(cfg.smultixcan_out(m))})
    steps.append({"name": "summarize", "done": _exists(os.path.join(cfg.results, "SUMMARY.md"))})
    screen = cfg.trait.get("drug_screen", {})
    if screen.get("enabled"):
        from .drug_screen import reuse_screen
        try:
            done = reuse_screen(screen["config"], {"trait_id": cfg.trait_id, "results": cfg.results, "work": cfg.work})
            steps.append({"name": "drug_screen", "done": done})
        except (SystemExit, ValueError, OSError, KeyError) as exc:
            steps.append({"name": "drug_screen", "done": False, "detail": str(exc)})
    n = sum(bool(s["done"]) for s in steps)
    current = next((s["name"] for s in steps if not s["done"]), None)
    return {"steps": steps, "done": n, "total": len(steps), "current": current, "finished": current is None}
