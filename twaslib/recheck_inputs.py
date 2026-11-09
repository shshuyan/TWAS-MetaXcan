"""Validated local evidence inputs. No downloads, LLM, or implicit fallback."""
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd
import yaml


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def read_table(path):
    return pd.read_csv(path, sep="," if str(path).endswith((".csv", ".csv.gz")) else "\t")


def symbol(value):
    if pd.isna(value):
        return ""
    return str(value).strip().upper()


def load_settings(path, variables=None):
    path = Path(path).resolve()
    settings = yaml.safe_load(path.read_text())
    if not isinstance(settings, dict):
        raise ValueError("Recheck config must be a mapping")
    env = {"root": str(Path(__file__).resolve().parents[1]), "config_dir": str(path.parent), **(variables or {})}
    def expand(value):
        if isinstance(value, str):
            for k, v in env.items():
                value = value.replace("${" + k + "}", str(v))
            if "${" in value:
                raise ValueError("Unresolved config variable: " + value)
        elif isinstance(value, dict):
            return {k: expand(v) for k, v in value.items()}
        elif isinstance(value, list):
            return [expand(v) for v in value]
        return value
    return expand(settings), path


def input_path(value, base):
    p = Path(value).expanduser()
    p = p if p.is_absolute() else base / p
    p = p.resolve()
    if not p.is_file():
        raise ValueError(f"Missing input file: {p}")
    return p


def load_candidates(path, aliases=None, target_path=None):
    if str(path).endswith(".json"):
        records = json.loads(Path(path).read_text())
        if isinstance(records, dict):
            records = records.get("candidates")
    else:
        records = read_table(path).replace({np.nan: None}).to_dict("records")
    if not isinstance(records, list) or not records:
        raise ValueError("Candidates must be a nonempty list/table")
    targets = {}
    if target_path:
        data = json.loads(Path(target_path).read_text())
        for r in data.get("candidate_targets", []):
            targets.setdefault(r["pert_id"], set()).update(map(symbol, r.get("gene_symbols", [])))
    result = []
    for r in records:
        cid = r.get("candidate_id") or r.get("pert_id")
        name = r.get("drug") or r.get("pert_name")
        if not cid or not name:
            raise ValueError("Each candidate needs candidate_id/pert_id and drug/pert_name")
        raw_targets = r.get("targets") or []
        if isinstance(raw_targets, str):
            raw_targets = raw_targets.split(";")
        t = set(map(symbol, raw_targets)) | targets.get(cid, set())
        result.append({"candidate_id": str(cid), "drug": str(name), "aliases": (aliases or {}).get(cid, []),
                       "targets": sorted(t - {""}), "upstream_record": r})
    if len({r["candidate_id"] for r in result}) != len(result):
        raise ValueError("Duplicate candidate IDs")
    identities = {}
    for r in result:
        for name in [r["drug"], *r["aliases"]]:
            key = name.strip().casefold()
            if key in identities and identities[key] != r["candidate_id"]:
                raise ValueError(f"Ambiguous candidate name/alias: {name}")
            identities[key] = r["candidate_id"]
    return result


def load_genetics(spec, path):
    d = read_table(path)
    mode = spec.get("mode")
    column = spec.get("gene_column", "gene_name" if mode == "twas" else "mapped_genes")
    if column not in d:
        raise ValueError(f"Genetics input lacks {column}")
    if mode == "twas":
        if spec.get("complete_tested_results") is not True:
            raise ValueError("TWAS requires complete_tested_results: true; significant-only tables cannot define background")
        pc = spec.get("pvalue_column", "pvalue")
        if pc not in d:
            raise ValueError(f"TWAS input lacks {pc}")
        d = d.assign(_gene=d[column].map(symbol), _p=pd.to_numeric(d[pc], errors="coerce"))
        d = d[(d._gene != "") & d._p.notna()]
        if d._gene.duplicated().any():
            raise ValueError("TWAS input must have one test per gene: use S-MultiXcan or one preselected tissue")
        if not len(d) or (~d._p.between(0, 1)).any():
            raise ValueError("Invalid/empty TWAS p-values")
        threshold = float(spec.get("alpha", .05)) / len(d)
        hits = set(d.loc[d._p < threshold, "_gene"])
        return hits, set(d._gene), {"mode": mode, "threshold": threshold, "tests_before_assay_intersection": len(d),
                                   "mapping_strength": "predicted-expression association; not colocalization/causality"}
    if mode != "gwas_mapped":
        raise ValueError("genetics.mode must be twas or gwas_mapped; no automatic fallback")
    if not spec.get("study"):
        raise ValueError("GWAS mapped enrichment requires an explicit study identifier")
    pc = spec.get("pvalue_column", "pvalue")
    p = pd.to_numeric(d[pc], errors="coerce")
    if p.isna().any() or (~p.between(0, 1)).any():
        raise ValueError("Invalid GWAS p-values")
    hits = {symbol(g) for value in d.loc[p <= float(spec.get("threshold", 5e-8)), column].dropna()
            for g in re.split(r";", str(value))} - {""}
    return hits, None, {"mode": mode, "study": spec["study"], "threshold": spec.get("threshold", 5e-8),
                        "associations": int((p <= spec.get("threshold", 5e-8)).sum()),
                        "mapping_strength": "locus-mapped genes only; not causal/eQTL direction or LD-aware enrichment"}


def load_de(path):
    d = read_table(path).rename(columns={"logFC": "log_fc", "P.Value": "p_value", "adj.P.Val": "fdr_bh"})
    for c in ["gene", "log_fc", "p_value", "fdr_bh"]:
        if c not in d:
            raise ValueError(f"DE input lacks {c}")
    d["gene"] = d.gene.map(symbol)
    for c in ["log_fc", "p_value", "fdr_bh"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.dropna(subset=["log_fc", "p_value", "fdr_bh"])
    d = d[(d.gene != "") & np.isfinite(d.log_fc)]
    if not len(d) or (~d.p_value.between(0, 1)).any() or (~d.fdr_bh.between(0, 1)).any():
        raise ValueError("Invalid DE input")
    return d.sort_values(["p_value", "gene"]).drop_duplicates("gene").set_index("gene")


def disease_query(d, background, n=100, fdr=.05, minimum=10):
    eligible = d.loc[d.index.intersection(sorted(background))]
    sig = eligible[eligible.fdr_bh < fdr]
    up = sig[sig.log_fc > 0].sort_values(["p_value"]).head(n).index.tolist()
    down = sig[sig.log_fc < 0].sort_values(["p_value"]).head(n).index.tolist()
    if min(len(up), len(down)) < minimum:
        return up, down, "insufficient_query"
    return up, down, "ready"


def load_lincs(path, candidates, cells):
    wanted = {c["candidate_id"]: c for c in candidates}
    signatures, seen = [], set()
    with Path(path).open() as f:
        for line in f:
            s = json.loads(line)
            if s.get("pert_id") not in wanted or s.get("cell_id") not in cells:
                continue
            sid = s["sig_id"]
            if sid in seen:
                raise ValueError("Duplicate LINCS signature ID: " + sid)
            seen.add(sid)
            u, d = set(map(symbol, s["up"])) - {""}, set(map(symbol, s["down"])) - {""}
            if u & d:
                raise ValueError("LINCS signature contains contradictory gene directions: " + sid)
            signatures.append({**s, "up": u, "down": d})
    return signatures
