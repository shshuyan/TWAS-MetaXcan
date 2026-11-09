"""Optional post-screen AD genetic-enrichment and perturbation-reversal checks.

All inputs are explicit local files. Frozen candidates are never replaced with new
hits. The two evidence axes are kept separate; their p-values are not pooled.
"""
from collections import Counter, defaultdict
from datetime import datetime, timezone
import json
import logging
from pathlib import Path
import platform
import re

import h5py
import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata

from .recheck_inputs import (disease_query, input_path, load_candidates, load_de,
                            load_genetics, load_lincs, load_settings, read_table, sha256, symbol)
from .recheck_stats import bh, cmap_calibration, cmap_connectivity, enrichment, lincs_calibration


def json_write(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def table_write(path, records):
    pd.DataFrame(records).to_csv(path, sep="\t", index=False)


def name_key(name):
    return str(name).strip().casefold()


class CMapReference:
    def __init__(self, matrix, mapping):
        with h5py.File(matrix) as f:
            raw = np.asarray(f["assay"], dtype=np.float32)
            ids = [x.decode() if isinstance(x, bytes) else str(x) for x in f["rownames"][:].ravel()]
            self.names = [x.decode() if isinstance(x, bytes) else str(x) for x in f["colnames"][:].ravel()]
        if raw.shape != (len(self.names), len(ids)):
            if raw.shape == (len(ids), len(self.names)):
                raw = raw.T
            else:
                raise ValueError("CMap dimensions do not match row/column names")
        if not np.isfinite(raw).all():
            raise ValueError("Nonfinite CMap ranks")
        # Input contract: legacy build02, lower ranks denote higher expression.
        self.ranks = rankdata(raw, axis=1, method="average").astype(np.float32)
        d = read_table(mapping).fillna("")
        alias = dict(zip(d.cmap_gene_id.astype(str), d.symbol.map(symbol)))
        genes = [alias.get(g, "") for g in ids]
        counts = Counter(genes)
        self.gi = {g: i for i, g in enumerate(genes) if g and counts[g] == 1}
        self.genes = genes
        group = defaultdict(list)
        for i, name in enumerate(self.names):
            if len(name.split("__")) < 2:
                raise ValueError("Malformed CMap signature ID: " + name)
            group[name_key(name.split("__")[0])].append(i)
        self.drugs = sorted(group)
        self.groups = [group[d] for d in self.drugs]
        self.by_drug = group

    def match(self, candidate):
        names = [candidate["drug"], *candidate["aliases"]]
        found = {name_key(x) for x in names if name_key(x) in self.by_drug}
        if len(found) > 1:
            raise ValueError("Multiple CMap identities for " + candidate["candidate_id"])
        return next(iter(found)) if found else None

    def footprint(self, compound, size):
        hit = np.any((self.ranks[self.by_drug[compound]] <= size) |
                     (self.ranks[self.by_drug[compound]] > self.ranks.shape[1] - size), axis=0)
        return {g for g, ix in self.gi.items() if hit[ix]}


def finalize_genetic_tests(records, alpha):
    tested = [r for r in records if r["status"] == "tested"]
    for r, q in zip(tested, bh([r["pvalue"] for r in tested])):
        r["q_bh"] = float(q)
        r["support"] = bool(q < alpha and r["fold_enrichment"] is not None and r["fold_enrichment"] > 1)
    for r in records:
        r.setdefault("q_bh", None); r.setdefault("support", False)
    return records


def evidence_decision(genetic_support, reversal_states, genetics_mode):
    """Missing/insufficient data cannot turn into a negative biological claim."""
    if not reversal_states or any(s in {"not_covered", "insufficient_query"} for s in reversal_states):
        return "incomplete_reversal_evidence"
    reversal = all(s == "supported" for s in reversal_states)
    if genetic_support and reversal:
        return "both_signals_exploratory" if genetics_mode == "gwas_mapped" else "both_association_signals"
    if reversal:
        return "reversal_only_genetics_not_supported"
    if genetic_support:
        return "genetics_only_reversal_not_reproduced"
    return "neither_check_supported"


def run_recheck(config_path, output=None, variables=None):
    cfg, config_path = load_settings(config_path, variables)
    base = config_path.parent
    if cfg.get("disease") != "AD":
        raise ValueError("This recheck is scoped to disease: AD")
    alpha = float(cfg.get("alpha", .05))
    b = int(cfg.get("permutations", 10000)); seed = int(cfg.get("seed", 20261108))
    if not 0 < alpha < 1 or b < 100:
        raise ValueError("Require 0<alpha<1 and permutations>=100 (10000 recommended)")
    minimum_cohorts = int(cfg.get("minimum_cohorts", 2))
    if minimum_cohorts < 2:
        raise ValueError("Candidate verification requires at least two disease cohorts")
    manifest = []
    def track(value, role):
        p = input_path(value, base)
        manifest.append({"role": role, "path": str(p), "bytes": p.stat().st_size, "sha256": sha256(p)})
        return p
    candidate_path = track(cfg["candidates"], "frozen_candidates")
    target_path = track(cfg["target_context"], "target_context") if cfg.get("target_context") else None
    candidates = load_candidates(candidate_path, cfg.get("aliases"), target_path)
    genetic_path = track(cfg["genetics"]["path"], "genetic_input")
    hits, tested_genes, genetic_meta = load_genetics(cfg["genetics"], genetic_path)
    cohorts = []
    for s in cfg["cohorts"]:
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", s.get("id", "")) or s.get("role") not in {"discovery", "validation"}:
            raise ValueError("Cohorts require id and role")
        p = track(s["differential_expression"], "DE:" + s["id"])
        cohorts.append((s, load_de(p)))
    if len(cohorts) < minimum_cohorts or len({s["id"] for s, d in cohorts}) != len(cohorts):
        raise ValueError("Insufficient or duplicate cohort IDs")
    if {s["role"] for s, d in cohorts} != {"discovery", "validation"}:
        raise ValueError("Both discovery and validation roles are required")
    de_hashes = [r["sha256"] for r in manifest if r["role"].startswith("DE:")]
    if len(set(de_hashes)) != len(de_hashes):
        raise ValueError("Same DE file content cannot be counted as two cohorts")
    ls = cfg["lincs"]
    if not isinstance(ls.get("reused_from_upstream_screen"), bool):
        raise ValueError("lincs.reused_from_upstream_screen must explicitly be true or false")
    cells = ls.get("cells", ["NEU", "NPC"]); primary = ls.get("primary_cell", "NEU")
    if primary not in cells:
        raise ValueError("primary_cell must be included in lincs.cells")
    lp = track(ls["signatures"], "lincs_signatures")
    ap = track(ls["gene_universe"], "lincs_assay_annotation")
    annotation = read_table(ap)
    assay = set(annotation[ls.get("gene_column", "gene")].map(symbol)) - {""}
    if len(assay) < 20:
        raise ValueError("LINCS assay annotation has insufficient genes")
    signatures = load_lincs(lp, candidates, cells)
    by_id_cell = defaultdict(list)
    for i, s in enumerate(signatures):
        by_id_cell[(s["pert_id"], s["cell_id"])].append(i)
    cs = cfg.get("cmap")
    cmap = None
    if cs:
        if cs.get("rank_direction") != "lower_is_up":
            raise ValueError("CMap rank_direction must explicitly be lower_is_up for the legacy matrix")
        cmap = CMapReference(track(cs["matrix"], "cmap_matrix"), track(cs["gene_mapping"], "cmap_mapping"))
    out = Path(output or cfg["output"])
    out = (out if out.is_absolute() else base / out).resolve()
    # A fresh directory prevents partial results from being mistaken for an older complete run.
    out.mkdir(parents=True, exist_ok=False)
    json_write(out / "resolved_config.json", cfg)
    json_write(out / "frozen_candidates.json", candidates)
    raw_genes = sorted(hits)
    table_write(out / "genetic_hit_genes.tsv", [{"gene": g} for g in raw_genes])
    genetic_records, target_records, signature_records, reversal_records, query_records = [], [], [], [], []
    assay_meta = []
    for source, universe in [("lincs:" + primary, assay), *([("cmap:legacy", set(cmap.gi))] if cmap else [])]:
        background = universe if tested_genes is None else universe & tested_genes
        assay_meta.append({"source": source, "background_genes": len(background), "genetic_hits_in_background": len(hits & background)})
        for c in candidates:
            if source.startswith("lincs:"):
                ix = by_id_cell.get((c["candidate_id"], primary), [])
                footprint = set().union(*(signatures[i]["up"] | signatures[i]["down"] for i in ix)) if ix else set()
                count, identity = len(ix), "exact_pert_id"
            else:
                match = cmap.match(c)
                footprint = cmap.footprint(match, int(cs.get("tail_size", 250))) if match else set()
                count, identity = (len(cmap.by_drug[match]) if match else 0), "name_or_explicit_alias_unverified_structure"
            r = enrichment(footprint, hits, background)
            r.update(candidate_id=c["candidate_id"], drug=c["drug"], source=source,
                     reference_signatures=count, identity_basis=identity,
                     footprint_genes_outside_background=len(footprint - background),
                     genetics_mode=genetic_meta["mode"])
            if not count:
                r["status"] = "not_covered"
            genetic_records.append(r)
    finalize_genetic_tests(genetic_records, alpha)
    for c in candidates:
        target_records.append({"candidate_id": c["candidate_id"], "drug": c["drug"], "known_targets": c["targets"],
                               "genetic_hit_targets": sorted(set(c["targets"]) & hits),
                               "status": "annotation_only_no_direction" if c["targets"] else "targets_not_provided"})
    query_n = int(cfg.get("query_size", 100)); fdr = float(cfg.get("query_fdr", .05))
    if query_n < 10 or not 0 < fdr < 1:
        raise ValueError("query_size must be >=10 and query_fdr in (0,1)")
    for cohort, d in cohorts:
        acc = cohort["id"]
        for source in ["lincs", *(["cmap"] if cmap else [])]:
            background = (assay if source == "lincs" else set(cmap.gi)) & set(d.index)
            up, down, status = disease_query(d, background, query_n, fdr)
            query_records.append({"cohort": acc, "source": source, "status": status, "up": len(up), "down": len(down),
                                  "background_genes": len(background), "all_de_genes": len(d), "gene_sets": {"up": up, "down": down}})
            if status != "ready":
                for c in candidates:
                    for cell in (cells if source == "lincs" else ["legacy"]):
                        reversal_records.append({"candidate_id": c["candidate_id"], "drug": c["drug"], "cohort": acc,
                                                 "source": source, "cell": cell, "status": "insufficient_query"})
                continue
            if source == "lincs":
                genes = sorted(background); gi = {g: i for i, g in enumerate(genes)}
                signs = np.zeros((len(signatures), len(genes)), dtype=np.int8)
                for i, s in enumerate(signatures):
                    for g in s["up"] & background: signs[i, gi[g]] = 1
                    for g in s["down"] & background: signs[i, gi[g]] = -1
                keys = sorted(by_id_cell)
                if keys:
                    indiv, observed, p, q = lincs_calibration(signs, [gi[g] for g in up], [gi[g] for g in down],
                                                            [by_id_cell[k] for k in keys], b, seed)
                    for i, s in enumerate(signatures):
                        opposite = (set(up) & s["down"]) | (set(down) & s["up"])
                        same = (set(up) & s["up"]) | (set(down) & s["down"])
                        signature_records.append({"cohort": acc, "candidate_id": s["pert_id"], "source": source,
                                                  "signature_id": s["sig_id"], "cell": s["cell_id"], "dose": s.get("dose"),
                                                  "time": s.get("time"), "score": float(indiv[i]),
                                                  "opposite_genes": sorted(opposite), "same_genes": sorted(same),
                                                  "opposite_genetic_genes": sorted(opposite & hits),
                                                  "url": s.get("signature_url")})
                else:
                    p = q = observed = []
                for c in candidates:
                    for cell in cells:
                        key = c["candidate_id"], cell
                        r = {"candidate_id": c["candidate_id"], "drug": c["drug"], "cohort": acc,
                             "source": source, "cell": cell, "identity_basis": "exact_pert_id", "status": "not_covered"}
                        if key in by_id_cell:
                            j = keys.index(key); ix = by_id_cell[key]
                            r.update(score=float(observed[j]), empirical_p=float(p[j]), q_bh=float(q[j]),
                                     signatures=len(ix), negative_signature_fraction=float((indiv[ix] < 0).mean()),
                                     multiplicity_tests=len(keys), status="supported" if observed[j] < 0 and q[j] < alpha else "not_supported")
                        reversal_records.append(r)
            else:
                genes = sorted(background); ix = np.array([cmap.gi[g] for g in genes])
                logging.info("%s: calibrating all %d CMap compounds, B=%d", acc, len(cmap.drugs), b)
                observed, p, q = cmap_calibration(cmap.ranks, [cmap.gi[g] for g in up], [cmap.gi[g] for g in down],
                                                 ix, cmap.groups, b, seed,
                                                 lambda n: logging.info("%s CMap %d/%d", acc, n, b))
                table_write(out / (acc + "_all_cmap_compounds.tsv"),
                            [{"compound": drug, "score": float(observed[i]), "empirical_p": float(p[i]), "q_bh": float(q[i])}
                             for i, drug in enumerate(cmap.drugs)])
                indiv = cmap_connectivity(cmap.ranks, [cmap.gi[g] for g in up], [cmap.gi[g] for g in down])
                for c in candidates:
                    match = cmap.match(c)
                    r = {"candidate_id": c["candidate_id"], "drug": c["drug"], "cohort": acc, "source": source,
                         "cell": "legacy", "status": "not_covered", "identity_basis": "name_or_explicit_alias_unverified_structure"}
                    if match:
                        j = cmap.drugs.index(match); ref_ix = cmap.by_drug[match]
                        r.update(matched_compound=match, score=float(observed[j]), empirical_p=float(p[j]), q_bh=float(q[j]),
                                 signatures=len(ref_ix), multiplicity_tests=len(cmap.drugs),
                                 status="supported" if observed[j] < 0 and q[j] < alpha else "not_supported")
                        for i in ref_ix:
                            signature_records.append({"cohort": acc, "candidate_id": c["candidate_id"], "source": source,
                                                      "signature_id": cmap.names[i], "cell": cmap.names[i].split("__")[1],
                                                      "score": float(indiv[i])})
                    reversal_records.append(r)
    candidate_rows = []
    for c in candidates:
        gr = [r for r in genetic_records if r["candidate_id"] == c["candidate_id"]]
        rr = [r for r in reversal_records if r["candidate_id"] == c["candidate_id"]]
        primary_rows = [r for r in rr if r["source"] == "lincs" and r["cell"] == primary]
        primary_genetics = next(r for r in gr if r["source"] == "lincs:" + primary)
        row = {"candidate_id": c["candidate_id"], "drug": c["drug"], "genetics_mode": genetic_meta["mode"],
               "genetic_status": ("supported" if primary_genetics["support"] else
                                  "not_supported" if primary_genetics["status"] == "tested" else primary_genetics["status"]),
               "genetic_q": primary_genetics["q_bh"], "genetic_overlap": primary_genetics["overlap_genes"],
               "decision": evidence_decision(primary_genetics["support"], [r["status"] for r in primary_rows], genetic_meta["mode"])
                           if primary_genetics["status"] == "tested" else "incomplete_genetic_evidence",
               "genetic_support_sources": [r["source"] for r in gr if r["support"]],
               "reference_reused_from_upstream": ls["reused_from_upstream_screen"],
               "independent_drug_experiment": False,
               "wet_lab_ready": False}
        covered_cmap = [r for r in rr if r["source"] == "cmap" and r["status"] not in {"not_covered", "insufficient_query"}]
        row["cmap_cross_cohort_status"] = ("not_covered" if not covered_cmap else
                                            "supported" if all(r["status"] == "supported" for r in covered_cmap) and len(covered_cmap) == len(cohorts)
                                            else "not_consistently_supported")
        for r in rr:
            prefix = r["source"] + "_" + r["cell"] + "_" + r["cohort"]
            row[prefix + "_status"] = r["status"]
            row[prefix + "_score"] = r.get("score")
            row[prefix + "_q"] = r.get("q_bh")
        candidate_rows.append(row)
    methods = {"genetics": genetic_meta, "genetic_backgrounds": assay_meta,
               "genetic_test": "one-sided Fisher overrepresentation; BH across all covered frozen candidate x primary-LINCS/legacy-CMap footprints",
               "genetic_limitations": "locus/gene dependence and perturbation ascertainment not corrected; exploratory, not a heritability enrichment or causal test",
               "cmap": "legacy lower-is-up unweighted KS; median across reference cells; empirical BH across all reference compounds; not CLUE tau",
               "lincs": "signed tail-overlap median per candidate/cell; empirical BH across all covered frozen candidate/cell groups per cohort",
               "calibration": "random gene queries preserving up/down sizes in assay-measured disease background; not sample permutations; no coexpression preservation",
               "permutations": b, "seed": seed, "alpha": alpha, "query_size": query_n, "query_fdr": fdr,
               "primary_cell": primary, "cohorts": [s for s, d in cohorts],
               "upstream_reference_reused": ls["reused_from_upstream_screen"],
               "independent_donors_verified": False, "pvalues_combined": False, "llm_calls": 0}
    metadata = {"status": "complete", "completed_utc": datetime.now(timezone.utc).isoformat(),
                "candidate_count": len(candidates), "candidates_sha256": sha256(candidate_path),
                "decisions": dict(Counter(r["decision"] for r in candidate_rows)),
                "wet_lab_ready_count": 0, "methods": methods,
                "software": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, "pandas": pd.__version__}}
    for filename, records in [("genetic_enrichment", genetic_records), ("target_genetic_context", target_records),
                              ("reversal_by_cohort", reversal_records), ("signature_evidence", signature_records),
                              ("candidate_recheck", candidate_rows), ("disease_queries", query_records)]:
        json_write(out / (filename + ".json"), records)
        table_write(out / (filename + ".tsv"), records)
    for p, role in [(config_path, "recheck_config"), (Path(__file__), "recheck_code"),
                    (Path(__file__).with_name("recheck_stats.py"), "statistics_code"),
                    (Path(__file__).with_name("recheck_inputs.py"), "input_code"),
                    (Path(__file__).with_name("recheck_report.py"), "report_code")]:
        manifest.append({"role": role, "path": str(p), "bytes": p.stat().st_size, "sha256": sha256(p)})
    json_write(out / "input_manifest.json", manifest)
    from .recheck_report import write_report
    write_report(out, candidate_rows, genetic_records, reversal_records, metadata)
    # Summary written last is the completion marker.
    json_write(out / "summary.json", metadata)
    return out
