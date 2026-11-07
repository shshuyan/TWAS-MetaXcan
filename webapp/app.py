"""TWAS web app: search public GWAS, download, run the pipeline, browse gene-level results.

  /opt/conda/envs/MetaXcan/bin/python -m uvicorn webapp.app:app --host 0.0.0.0 --port 8000
  (run from the repository root, or use `bash webapp/start.sh`)
"""
import glob
import logging
import os
import queue
import re
import signal
import subprocess
import sys
import threading
import time

import pandas as pd
import yaml
from fastapi import Body, FastAPI, HTTPException
from fastapi.responses import FileResponse, HTMLResponse, PlainTextResponse

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from twaslib import config as tconfig  # noqa: E402
from twaslib import sources  # noqa: E402
from twaslib.progress import progress  # noqa: E402

TRAITS_DIR = os.path.join(ROOT, "config", "traits")
LOG_DIR = os.path.join(ROOT, "logs")
PY = sys.executable
_pipe = yaml.safe_load(open(os.path.join(ROOT, "config", "pipeline.yaml")))
_paths = tconfig._subst(_pipe["paths"], {"root": ROOT})
CATALOG_DIR = os.path.join(_paths["raw_dir"], "_catalogs")
RAW_DIR = _paths["raw_dir"]

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
app = FastAPI(title="TWAS pipeline")

# ----------------------------------------------------------------------------- job runner
_jobs = queue.Queue()
_state = {}          # trait_id -> {"status": queued|running|done|failed|cancelled, ...}
_procs = {}


def _worker():
    while True:
        tid, cfg_path, args = _jobs.get()
        if _state.get(tid, {}).get("status") == "cancelled":
            continue
        # One pipeline at a time: wait for runs started outside the app (memory-bound)
        while _external_running(""):
            _state[tid]["waiting"] = "waiting for another pipeline run to finish"
            time.sleep(30)
            if _state.get(tid, {}).get("status") == "cancelled":
                break
        if _state.get(tid, {}).get("status") == "cancelled":
            continue
        os.makedirs(LOG_DIR, exist_ok=True)
        log = os.path.join(LOG_DIR, f"{tid}.log")
        _state[tid] = {"status": "running", "started": time.time(), "log": log}
        with open(log, "a") as fh:
            fh.write(f"\n==== {time.ctime()} twas.py run {' '.join(args)}\n")
            fh.flush()
            p = subprocess.Popen([PY, os.path.join(ROOT, "twas.py"), "run", cfg_path] + args, cwd=ROOT,
                                 stdout=fh, stderr=subprocess.STDOUT, start_new_session=True)
            _procs[tid] = p
            rc = p.wait()
        _procs.pop(tid, None)
        if _state[tid]["status"] != "cancelled":
            _state[tid].update(status="done" if rc == 0 else "failed", rc=rc, ended=time.time())


threading.Thread(target=_worker, daemon=True).start()


def _external_running(cfg_path):
    """A run started outside the app (e.g. run_detached.sh) for this config."""
    try:
        out = subprocess.run(["pgrep", "-af", "twas.py run"], capture_output=True, text=True).stdout
    except FileNotFoundError:
        return False
    name = os.path.basename(cfg_path)
    mine = {str(p.pid) for p in _procs.values()}
    return any(name in line and line.split()[0] not in mine for line in out.splitlines()
               if "pgrep" not in line)


def _trait_files():
    return sorted(f for f in glob.glob(os.path.join(TRAITS_DIR, "*.yaml")) if not f.endswith("TEMPLATE.yaml"))


def _load(tid):
    for f in _trait_files():
        with open(f) as fh:
            if (yaml.safe_load(fh) or {}).get("id") == tid:
                return f, tconfig.load(f)
    raise HTTPException(404, f"trait {tid} not found")


# ----------------------------------------------------------------------------- API
@app.get("/", response_class=HTMLResponse)
def index():
    return open(os.path.join(ROOT, "webapp", "static", "index.html"), encoding="utf-8").read()


@app.get("/api/search")
def api_search(q: str, ancestry: str = "", src: str = "finngen,pan_ukbb,gwas_catalog"):
    if len(q.strip()) < 2:
        raise HTTPException(400, "query too short")
    return sources.search(CATALOG_DIR, q, ancestry or None, src.split(","))


@app.post("/api/catalogs/refresh")
def api_refresh():
    sources.load_catalogs(CATALOG_DIR, refresh=True)
    return {"ok": True}


@app.get("/api/resolve")
def api_resolve(source: str, id: str, pop: str = ""):
    try:
        return sources.resolve(CATALOG_DIR, source, id, pop or None)
    except KeyError:
        raise HTTPException(404, f"{source} {id} not found")
    except Exception as e:  # noqa: BLE001
        raise HTTPException(502, f"could not resolve {source} {id}: {e}")


@app.get("/api/tissues")
def api_tissues():
    cfg_any = _pipe["reference"]["models"]["eqtl"]["dir"]
    d = tconfig._subst(cfg_any, {**_paths, "root": ROOT, "data_dir": _paths["data_dir"]})
    return sorted(re.sub(r"^mashr_|\.db$", "", os.path.basename(f)) for f in glob.glob(os.path.join(d, "*.db")))


@app.post("/api/jobs")
def api_create(spec: dict = Body(...)):
    tid = spec.get("trait_id", "")
    if not re.fullmatch(r"[A-Za-z0-9_.-]{1,80}", tid):
        raise HTTPException(400, "trait id may only contain letters, digits, _ . -")
    if spec.get("build") not in tconfig.BUILD_ALIASES:
        raise HTTPException(400, "set the genome build (hg19 or hg38)")
    path = os.path.join(TRAITS_DIR, f"{tid}.yaml")
    if os.path.exists(path) and not spec.get("overwrite"):
        raise HTTPException(409, f"config/traits/{tid}.yaml already exists")
    sources.write_trait_config(path, spec, RAW_DIR)
    try:
        tconfig.load(path)
    except SystemExit as e:
        os.remove(path)
        raise HTTPException(400, str(e))
    return _enqueue(tid, path, [])


def _enqueue(tid, path, args):
    if _state.get(tid, {}).get("status") in ("queued", "running") or _external_running(path):
        raise HTTPException(409, f"{tid} is already queued or running")
    _state[tid] = {"status": "queued", "queued": time.time()}
    _jobs.put((tid, path, args))
    return {"trait_id": tid, "status": "queued", "config": os.path.relpath(path, ROOT)}


@app.get("/api/traits")
def api_traits():
    out = []
    for f in _trait_files():
        try:
            cfg = tconfig.load(f)
        except SystemExit as e:
            out.append({"trait_id": os.path.basename(f)[:-5], "error": str(e), "config": os.path.relpath(f, ROOT)})
            continue
        st = dict(_state.get(cfg.trait_id, {}))
        if not st.get("status") and _external_running(f):
            st["status"], st["external"] = "running", True
        pr = progress(cfg)
        st.setdefault("status", "done" if pr["finished"] else "idle")
        out.append({"trait_id": cfg.trait_id, "description": cfg.trait.get("description", ""),
                    "ancestry": cfg.trait.get("ancestry"), "source": cfg.trait["input"].get("source"),
                    "skip_imputation": bool(cfg.trait.get("skip_imputation")),
                    "model_types": cfg.trait.get("model_types", ["eqtl"]),
                    "config": os.path.relpath(f, ROOT), "progress": pr, **st})
    return out


@app.post("/api/traits/{tid}/run")
def api_run(tid: str, body: dict = Body(default={})):
    path, _ = _load(tid)
    args = []
    if body.get("from_step"):
        args += ["--from", body["from_step"]]
    return _enqueue(tid, path, args)


@app.post("/api/traits/{tid}/cancel")
def api_cancel(tid: str):
    st = _state.get(tid)
    if not st or st.get("status") not in ("queued", "running"):
        raise HTTPException(409, "not queued or running in this app")
    st["status"] = "cancelled"
    p = _procs.get(tid)
    if p:
        os.killpg(p.pid, signal.SIGTERM)
    return {"ok": True}


@app.get("/api/traits/{tid}/log", response_class=PlainTextResponse)
def api_log(tid: str, lines: int = 80):
    f = os.path.join(LOG_DIR, f"{tid}.log")
    if not os.path.exists(f):
        cand = sorted(glob.glob(os.path.join(LOG_DIR, "*.log")), key=os.path.getmtime)
        cand = [c for c in cand if tid in open(c, errors="replace").read(20000)]
        if not cand:
            return "(no log yet)"
        f = cand[-1]
    with open(f, errors="replace") as fh:
        return "".join(fh.readlines()[-lines:])


def _result_frame(cfg, kind, mtype):
    if kind == "smultixcan":
        p = cfg.smultixcan_out(mtype)
        if not os.path.exists(p):
            raise HTTPException(404, "S-MultiXcan results not available yet")
        d = pd.read_csv(p, sep="\t")
        keep = ["gene", "gene_name", "pvalue", "n", "n_indep", "p_i_best", "t_i_best", "z_mean", "z_sd"]
    else:
        files = glob.glob(os.path.join(cfg.spredixcan_dir(mtype), "*__PM__*.csv"))
        if not files:
            raise HTTPException(404, "S-PrediXcan results not available yet")
        fr = []
        for f in files:
            x = pd.read_csv(f, usecols=lambda c: c in {"gene", "gene_name", "zscore", "effect_size", "pvalue",
                                                       "pred_perf_r2", "n_snps_used", "n_snps_in_model"})
            x.insert(0, "tissue", re.search(r"__PM__(.*)\.csv$", f).group(1))
            fr.append(x)
        d = pd.concat(fr, ignore_index=True)
        keep = ["gene", "gene_name", "tissue", "zscore", "effect_size", "pvalue", "pred_perf_r2",
                "n_snps_used", "n_snps_in_model"]
    d = d[[c for c in keep if c in d.columns]].dropna(subset=["pvalue"]).sort_values("pvalue")
    alpha = cfg.pipe.get("alpha", 0.05)
    d["bonferroni_threshold"] = alpha / len(d)
    d["significant"] = d.pvalue < d.bonferroni_threshold
    return d


@app.get("/api/traits/{tid}/results")
def api_results(tid: str, kind: str = "smultixcan", model_type: str = "eqtl", q: str = "", limit: int = 200,
                significant_only: bool = False):
    _, cfg = _load(tid)
    d = _result_frame(cfg, kind, model_type)
    n_sig, n_tests = int(d.significant.sum()), len(d)
    thr = float(d.bonferroni_threshold.iloc[0]) if n_tests else None
    if q:
        d = d[d.gene_name.fillna("").str.contains(q, case=False) | d.gene.str.contains(q, case=False)]
    if significant_only:
        d = d[d.significant]
    top = d.head(limit)
    rows = top.astype(object).where(top.notna(), None).to_dict("records")
    return {"n_tests": n_tests, "n_significant": n_sig, "threshold": thr, "rows": rows}


@app.get("/api/traits/{tid}/download")
def api_download(tid: str, kind: str = "smultixcan", model_type: str = "eqtl"):
    _, cfg = _load(tid)
    d = _result_frame(cfg, kind, model_type)
    out = os.path.join(cfg.work, f"{tid}_{kind}_{model_type}_gene_level.tsv")
    d.to_csv(out, sep="\t", index=False)
    return FileResponse(out, filename=os.path.basename(out), media_type="text/tab-separated-values")


@app.get("/api/traits/{tid}/summary", response_class=PlainTextResponse)
def api_summary(tid: str):
    _, cfg = _load(tid)
    p = os.path.join(cfg.results, "SUMMARY.md")
    return open(p).read() if os.path.exists(p) else "(summary not available yet)"


@app.get("/api/traits/{tid}/config", response_class=PlainTextResponse)
def api_config(tid: str):
    path, _ = _load(tid)
    return open(path).read()
