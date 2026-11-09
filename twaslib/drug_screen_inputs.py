"""Build tissue-specific, signed TWAS queries without GEO or a drug shortlist."""
import json

import numpy as np
import pandas as pd

from .screen_io import read_table, symbol
from .screen_stats import bh


def load_twas_queries(spec, sp_path, sm_path):
    if spec.get("complete_tested_results") is not True:
        raise ValueError("Screening requires complete_tested_results: true, not significant-only tables")
    if spec.get("model_type") != "eqtl":
        raise ValueError("Directional expression screening currently requires model_type: eqtl")
    if spec.get("trait_direction") != "higher_is_AD_risk":
        raise ValueError("Confirm GWAS trait_direction: higher_is_AD_risk before interpreting Z signs")
    tissues = spec.get("tissues")
    primary = spec.get("primary_tissue")
    if not isinstance(tissues, list) or not tissues or len(set(tissues)) != len(tissues) or primary not in tissues:
        raise ValueError("Specify unique tissues and a primary_tissue within them before screening")
    alpha = float(spec.get("alpha", .05))
    if not 0 < alpha < 1:
        raise ValueError("TWAS alpha must lie in (0,1)")
    for key in ["require_joint_significance", "exclude_direction_conflicts"]:
        if not isinstance(spec.get(key, True), bool):
            raise ValueError(key + " must be boolean")
    sp, sm = read_table(sp_path), read_table(sm_path)
    for frame, columns in [(sp, ["gene_name", "tissue", "pvalue", "zscore"]), (sm, ["gene_name", "pvalue"])]:
        missing = set(columns) - set(frame)
        if missing:
            raise ValueError("TWAS input missing columns: " + ", ".join(sorted(missing)))
        frame["gene_name"] = frame.gene_name.map(symbol)
        frame["pvalue"] = pd.to_numeric(frame.pvalue, errors="coerce")
        frame["valid_test"] = (frame.gene_name != "") & frame.pvalue.between(0, 1)
        if "status" in frame:
            frame["valid_test"] &= frame.status.isna() | frame.status.astype(str).str.lower().isin(["", "ok", "success"])
    if sp.duplicated(["gene_name", "tissue"]).any() or sm.gene_name.duplicated().any():
        raise ValueError("Ambiguous duplicate genes: require one row per gene/tissue and one joint result per gene")
    sp["zscore"] = pd.to_numeric(sp.zscore, errors="coerce")
    sp["valid_test"] &= np.isfinite(sp.zscore)
    sm_valid = sm[sm.valid_test]
    chosen = sp[sp.tissue.isin(tissues) & sp.valid_test].copy()
    if not len(sm_valid) or set(tissues) - set(chosen.tissue):
        raise ValueError("Missing valid joint results or a requested S-PrediXcan tissue")
    # Correct before joint-hit selection or perturbation-background intersection.
    chosen["sp_q_bh"] = bh(chosen.pvalue)
    sm_threshold = alpha / len(sm_valid)
    sm["joint_significant"] = sm.valid_test & (sm.pvalue < sm_threshold)
    joint = sm_valid.set_index("gene_name").pvalue
    chosen["joint_pvalue"] = chosen.gene_name.map(joint)
    chosen["joint_significant"] = chosen.joint_pvalue < sm_threshold
    chosen["risk_direction"] = np.sign(chosen.zscore).astype(int)
    chosen["eligible_before_conflicts"] = (chosen.sp_q_bh < alpha) & (chosen.risk_direction != 0)
    if spec.get("require_joint_significance", True):
        chosen["eligible_before_conflicts"] &= chosen.joint_significant
    eligible = chosen[chosen.eligible_before_conflicts]
    conflicts = set(eligible.groupby("gene_name").risk_direction.nunique().loc[lambda s: s > 1].index)
    chosen["direction_conflict"] = chosen.gene_name.isin(conflicts)
    chosen["query_eligible"] = chosen.eligible_before_conflicts
    if spec.get("exclude_direction_conflicts", True):
        chosen["query_eligible"] &= ~chosen.direction_conflict
    chosen = chosen.sort_values(["tissue", "pvalue", "gene_name"])
    sp = sp.merge(chosen[["gene_name", "tissue", "sp_q_bh", "joint_pvalue", "joint_significant",
                          "risk_direction", "direction_conflict", "query_eligible"]],
                  on=["gene_name", "tissue"], how="left", validate="one_to_one")
    metadata = {"sp_tests_corrected": len(chosen), "joint_tests_corrected": len(sm_valid),
                "sp_invalid_rows": int((~sp.valid_test).sum()), "joint_invalid_rows": int((~sm.valid_test).sum()),
                "joint_bonferroni_threshold": sm_threshold, "sp_correction": "BH across all prespecified gene-tissue tests",
                "conflicting_genes": sorted(conflicts), "tissues": tissues, "primary_tissue": primary,
                "direction_source": "tissue-specific S-PrediXcan zscore; never S-MultiXcan z_mean",
                "signature_meaning": "genetically predicted expression association with AD risk, not measured disease DE"}
    return chosen, sp, sm, metadata


def tissue_query(frame, universe, spec):
    background = set(frame.gene_name)
    if spec.get("require_joint_significance", True):
        background &= set(frame.loc[frame.joint_pvalue.notna(), "gene_name"])
    background &= set(universe)
    eligible = frame[frame.query_eligible & frame.gene_name.isin(background)]
    limit = spec.get("max_genes_per_direction")
    minimum = spec.get("minimum_genes_per_direction", 10)
    if not isinstance(minimum, int) or minimum < 1 or (limit is not None and (not isinstance(limit, int) or limit < minimum)):
        raise ValueError("Require positive minimum_genes_per_direction and optional max >= minimum")
    up = eligible[eligible.risk_direction > 0].gene_name.tolist()
    down = eligible[eligible.risk_direction < 0].gene_name.tolist()
    return {"risk_positive": up[:limit], "risk_negative": down[:limit], "background": sorted(background),
            "eligible_positive_before_limit": len(up), "eligible_negative_before_limit": len(down),
            "status": "ready" if min(len(up[:limit]), len(down[:limit])) >= minimum else "insufficient_query"}


def load_screen_lincs(path, cells):
    if not isinstance(cells, list) or not cells or len(set(cells)) != len(cells):
        raise ValueError("LINCS cells must be an explicit, unique nonempty list")
    records, seen, names = [], set(), {}
    with path.open() as handle:
        for line in handle:
            record = json.loads(line)
            if record.get("cell_id") not in cells:
                continue
            if record.get("pert_type", "trt_cp") != "trt_cp":
                continue
            for field in ["sig_id", "pert_id", "pert_name"]:
                if not record.get(field):
                    raise ValueError("LINCS signature lacks " + field)
            sid, cid = record["sig_id"], record["pert_id"]
            if sid in seen:
                raise ValueError("Duplicate LINCS signature ID: " + sid)
            seen.add(sid)
            if cid in names and names[cid] != record["pert_name"]:
                raise ValueError("Conflicting compound names for LINCS pert_id " + cid)
            names[cid] = record["pert_name"]
            up, down = {symbol(g) for g in record["up"]} - {""}, {symbol(g) for g in record["down"]} - {""}
            if up & down:
                raise ValueError("Contradictory LINCS gene directions: " + sid)
            records.append({**record, "up": up, "down": down})
    if not records:
        raise ValueError("No compound signatures in the configured LINCS cells")
    return records
