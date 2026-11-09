"""Discover drugs from directional TWAS results across the supplied perturbation library."""
from collections import defaultdict
from datetime import datetime, timezone
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from .perturbation_reference import CMapReference
from .screen_io import json_write
from .drug_screen_inputs import load_screen_lincs, load_twas_queries, tissue_query
from .drug_screen_report import write_screen_report
from .screen_io import input_path, load_settings, read_table, sha256, symbol
from .screen_stats import bh, cmap_calibration, cmap_connectivity, lincs_calibration


def screen_output(cfg, config_path, output=None):
    path = Path(output or cfg["output"]).expanduser()
    return (path if path.is_absolute() else config_path.parent / path).resolve()


def reuse_screen(config_path, variables=None, force=False):
    """Reuse only complete outputs with the same config, sources, and implementation."""
    import json
    cfg, path = load_settings(config_path, variables)
    out = screen_output(cfg, path)
    if not out.exists():
        return False
    try:
        summary = json.loads((out / "summary.json").read_text())
        old = json.loads((out / "resolved_config.json").read_text())
        manifest = json.loads((out / "input_manifest.json").read_text())
        outputs = json.loads((out / "output_manifest.json").read_text())
        unchanged = summary["status"] == "complete" and old == cfg and bool(manifest) and bool(outputs)
        for record in manifest:
            p = Path(record["path"])
            unchanged &= p.is_file() and sha256(p) == record["sha256"]
        for record in outputs:
            p = out / record["file"]
            unchanged &= p.is_file() and sha256(p) == record["sha256"]
        if unchanged and not force:
            return True
    except (OSError, ValueError, KeyError, TypeError):
        pass
    raise SystemExit("Drug-screen output incomplete, changed, or --force requested. Use a NEW output directory: " + str(out))


def run_drug_screen(config_path, output=None, variables=None):
    cfg, config_path = load_settings(config_path, variables)
    if cfg.get("disease") != "AD":
        raise ValueError("This drug-screen implementation currently requires disease: AD")
    if "candidates" in cfg or "cohorts" in cfg:
        raise ValueError("drug-screen starts from TWAS: do not provide a candidate shortlist or GEO cohorts")
    alpha = float(cfg.get("alpha", .05))
    permutations, seed = cfg.get("permutations", 10000), cfg.get("seed", 20261109)
    if not 0 < alpha < 1 or not isinstance(permutations, int) or permutations < 100 or not isinstance(seed, int):
        raise ValueError("Require 0<alpha<1, integer permutations>=100 and an integer seed")
    primary = cfg.get("primary", {})
    if primary.get("source") not in {"lincs", "cmap"} or not primary.get("cell"):
        raise ValueError("Predeclare primary.source (lincs/cmap) and primary.cell")
    manifest = []
    def track(value, role):
        p = input_path(value, config_path.parent)
        manifest.append({"role": role, "path": str(p), "bytes": p.stat().st_size, "sha256": sha256(p)})
        return p
    spec = cfg["twas"]
    selected, full_sp, full_sm, twas_meta = load_twas_queries(
        spec, track(spec["spredixcan"], "full_spredixcan"), track(spec["smultixcan"], "full_smultixcan"))
    refs = {}
    if cfg.get("lincs"):
        ls = cfg["lincs"]
        if ls.get("compound_only") is not True:
            raise ValueError("LINCS input must be declared compound_only: true; curated drug signatures only")
        signatures = load_screen_lincs(track(ls["signatures"], "lincs_signatures"), ls["cells"])
        annotation = read_table(track(ls["gene_universe"], "lincs_assay_annotation"))
        universe = set(annotation[ls.get("gene_column", "gene")].map(symbol)) - {""}
        if len(universe) < 20:
            raise ValueError("LINCS assay universe must include at least 20 genes")
        groups = defaultdict(list)
        for i, s in enumerate(signatures):
            groups[(s["pert_id"], s["cell_id"])].append(i)
        keys = sorted(groups)
        refs["lincs"] = {"signatures": signatures, "universe": universe, "keys": keys,
                         "groups": [groups[k] for k in keys],
                         "names": {s["pert_id"]: s["pert_name"] for s in signatures}}
    if cfg.get("cmap"):
        cs = cfg["cmap"]
        if cs.get("rank_direction") != "lower_is_up":
            raise ValueError("Legacy CMap requires rank_direction: lower_is_up")
        cmap = CMapReference(track(cs["matrix"], "cmap_matrix"), track(cs["gene_mapping"], "cmap_mapping"))
        groups = defaultdict(list)
        for i, name in enumerate(cmap.names):
            compound, cell = name.split("__")[:2]
            groups[(compound.casefold().strip(), cell)].append(i)
        keys = sorted(groups)
        refs["cmap"] = {"reference": cmap, "universe": set(cmap.gi), "keys": keys,
                        "groups": [groups[k] for k in keys], "names": {drug: drug for drug, cell in keys}}
    if primary["source"] not in refs or not any(cell == primary["cell"] for _, cell in refs[primary["source"]]["keys"]):
        raise ValueError("Primary source/cell has no reference coverage; choose it before screening")
    out = screen_output(cfg, config_path, output)
    out.mkdir(parents=True, exist_ok=False)
    json_write(out / "resolved_config.json", cfg)
    full_sp.to_csv(out / "spredixcan_all_annotated.tsv.gz", sep="\t", index=False)
    full_sm.to_csv(out / "smultixcan_all_annotated.tsv.gz", sep="\t", index=False)
    selected.to_csv(out / "twas_signed_genes.tsv", sep="\t", index=False)
    results, conditions, queries = [], [], []
    for tissue in spec["tissues"]:
        frame = selected[selected.tissue == tissue]
        for source, ref in refs.items():
            query = tissue_query(frame, ref["universe"], spec)
            queries.append({"tissue": tissue, "source": source, **query})
            keys, groups = ref["keys"], ref["groups"]
            if query["status"] != "ready":
                results.extend({"source": source, "compound_id": cid, "drug": ref["names"][cid],
                                "cell": cell, "tissue": tissue, "signatures": len(groups[j]),
                                "status": "insufficient_query", "score": None, "empirical_p": None}
                               for j, (cid, cell) in enumerate(keys))
                continue
            genes = query["background"]
            up, down = query["risk_positive"], query["risk_negative"]
            logging.info("Screening %s/%s: %d drug-cell groups; %d positive, %d negative TWAS genes; B=%d",
                         tissue, source, len(keys), len(up), len(down), permutations)
            if source == "lincs":
                gi = {gene: i for i, gene in enumerate(genes)}
                signs = np.zeros((len(ref["signatures"]), len(genes)), dtype=np.int8)
                background = set(genes)
                for i, s in enumerate(ref["signatures"]):
                    for gene in s["up"] & background: signs[i, gi[gene]] = 1
                    for gene in s["down"] & background: signs[i, gi[gene]] = -1
                individual, observed, pvalues, _ = lincs_calibration(
                    signs, [gi[g] for g in up], [gi[g] for g in down], groups, permutations, seed)
                for i, s in enumerate(ref["signatures"]):
                    opposite = (set(up) & s["down"]) | (set(down) & s["up"])
                    same = (set(up) & s["up"]) | (set(down) & s["down"])
                    conditions.append({"source": source, "compound_id": s["pert_id"], "drug": s["pert_name"],
                                       "tissue": tissue, "cell": s["cell_id"], "signature_id": s["sig_id"],
                                       "dose": s.get("dose"), "time": s.get("time"), "score": float(individual[i]),
                                       "opposite_genes": sorted(opposite), "same_genes": sorted(same),
                                       "url": s.get("signature_url")})
            else:
                cmap = ref["reference"]
                ui, di = [cmap.gi[g] for g in up], [cmap.gi[g] for g in down]
                observed, pvalues, _ = cmap_calibration(cmap.ranks, ui, di, np.array([cmap.gi[g] for g in genes]),
                                                       groups, permutations, seed,
                                                       lambda n: logging.info("%s CMap %d/%d", tissue, n, permutations))
                individual = cmap_connectivity(cmap.ranks, ui, di)
                for i, sid in enumerate(cmap.names):
                    drug, cell, *_ = sid.split("__")
                    conditions.append({"source": source, "compound_id": drug.casefold().strip(), "drug": drug,
                                       "tissue": tissue, "cell": cell, "signature_id": sid,
                                       "score": float(individual[i]), "condition_metadata": "encoded_in_signature_id"})
            for j, (cid, cell) in enumerate(keys):
                results.append({"source": source, "compound_id": cid, "drug": ref["names"][cid], "cell": cell,
                                "tissue": tissue, "signatures": len(groups[j]), "score": float(observed[j]),
                                "empirical_p": float(pvalues[j]), "status": "tested",
                                "negative_signature_fraction": float((individual[groups[j]] < 0).mean())})
    # One family across every screened drug x cell x tissue x source; no shortlist filtering.
    tested = [r for r in results if r["status"] == "tested"]
    for record, q in zip(tested, bh([r["empirical_p"] for r in tested])):
        record["q_bh_global"] = float(q)
        record["status"] = "supported" if record["score"] < 0 and q < alpha else "not_supported"
    for record in results:
        record.setdefault("q_bh_global", None)
        record["primary_context"] = (record["source"] == primary["source"] and record["cell"] == primary["cell"]
                                     and record["tissue"] == spec["primary_tissue"])
        record["candidate"] = record["primary_context"] and record["status"] == "supported"
        record["wet_lab_ready"] = False
    results.sort(key=lambda r: (r["q_bh_global"] if r["q_bh_global"] is not None else 2,
                                r["score"] if r["score"] is not None else 2,
                                r["source"], r["compound_id"], r["cell"], r["tissue"]))
    candidates = [r for r in results if r["candidate"]]
    for rank, record in enumerate(candidates, 1):
        record["candidate_rank"] = rank
        record["other_contexts"] = [{k: r[k] for k in ["cell", "tissue", "status", "score", "q_bh_global"]}
                                    for r in results if r["source"] == record["source"]
                                    and r["compound_id"] == record["compound_id"] and not r["primary_context"]]
    def table(name, records, columns=None):
        json_write(out / (name + ".json"), records)
        pd.DataFrame(records, columns=columns).to_csv(out / (name + ".tsv"), sep="\t", index=False)
    table("drug_screen", results)
    table("candidate_drugs", candidates, list(results[0]) + ["candidate_rank", "other_contexts"] if not candidates else None)
    json_write(out / "twas_queries.json", queries)
    with (out / "signature_evidence.jsonl").open("w") as handle:
        import json
        for record in conditions:
            handle.write(json.dumps(record, ensure_ascii=False, allow_nan=False) + "\n")
    experiments = [{"question": "Does this drug reverse the tissue-specific AD-risk TWAS expression signature?",
                    "hypothesis": "Drug-induced expression opposes the prespecified signed genetic association query.",
                    "experiment": {"method": r["source"], "tissue": r["tissue"], "cell": r["cell"],
                                   "permutations": permutations, "seed": seed},
                    "result": r, "interpretation": r["status"],
                    "next_question": "Check colocalization, independent neural perturbations, target direction and exposure."}
                   for r in results]
    with (out / "virtual_experiments.jsonl").open("w") as handle:
        import json
        for record in experiments:
            handle.write(json.dumps(record, ensure_ascii=False, allow_nan=False) + "\n")
    summary = {"status": "complete", "completed_utc": datetime.now(timezone.utc).isoformat(),
               "input_kind": cfg.get("input_kind", "unspecified"),
               "mode": "twas_drug_discovery", "primary": {**primary, "tissue": spec["primary_tissue"]},
               "twas": twas_meta, "permutations": permutations, "seed": seed, "alpha": alpha,
               "correction": "BH across all tested source x tissue x drug x cell groups in this run",
               "tested_groups": len(tested), "total_groups": len(results), "candidate_count": len(candidates),
               "minimum_empirical_p": 1 / (permutations + 1),
               "isolated_hit_bh_resolution": min(1., len(tested) / (permutations + 1)),
               "isolated_hit_resolution_limited": len(tested) / (permutations + 1) >= alpha,
               "query_status": "ready" if any(q["status"] == "ready" for q in queries) else "insufficient_query",
               "references": {source: {"compounds": len({cid for cid, cell in ref["keys"]}),
                                       "drug_cell_groups": len(ref["keys"]),
                                       "signatures": len(ref["signatures"] if source == "lincs" else ref["reference"].names)}
                              for source, ref in refs.items()},
               "geo_required": False, "candidate_list_required": False, "llm_calls": 0,
               "independent_validation": False, "wet_lab_ready_count": 0,
               "limitations": ["TWAS association is not causality or measured patient differential expression",
                               "random-gene null does not preserve LD, gene correlation or model ascertainment",
                               "drug-cell condition medians are not independent biological replication",
                               "LINCS tail overlap and legacy CMap KS are not official CLUE tau",
                               "only the supplied reference library is screened; no claim of complete global LINCS coverage",
                               "CMap and LINCS compound identities are not merged without verified structure mapping"]}
    write_screen_report(out, summary, candidates, queries)
    code_root = Path(__file__).resolve().parent
    for name in ["drug_screen.py", "drug_screen_inputs.py", "drug_screen_report.py", "perturbation_reference.py",
                 "screen_io.py", "screen_stats.py"]:
        track(str(code_root / name), "implementation:" + name)
    json_write(out / "input_manifest.json", manifest)
    json_write(out / "output_manifest.json", [{"file": p.name, "sha256": sha256(p)} for p in sorted(out.iterdir()) if p.is_file()])
    json_write(out / "summary.json", summary)  # completion marker is always last
    return out
