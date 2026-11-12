"""Continuous TWAS–CMap reversal screen with primary-family error control."""
from collections import defaultdict
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .screen_io import symbol
from .screen_stats import bh


RESULT_COLUMNS = [
    "disease", "compound_id", "drug", "cell", "score", "empirical_p",
    "q_bh_primary_family", "signature_count", "negative_signature_fraction",
    "primary_context", "candidate",
]


def weighted_cosine_scores(matrix, disease_z, z_clip=5.0):
    """Score signature rows against a clipped, magnitude-weighted disease vector."""
    matrix = np.asarray(matrix, dtype=float)
    disease_z = np.asarray(disease_z, dtype=float)
    if matrix.ndim != 2 or disease_z.ndim != 1 or matrix.shape[1] != len(disease_z):
        raise ValueError("matrix must be signatures × genes aligned to disease_z")
    if not np.isfinite(matrix).all() or not np.isfinite(disease_z).all() or z_clip <= 0:
        raise ValueError("Continuous score inputs must be finite and z_clip positive")
    clipped = np.clip(disease_z, -z_clip, z_clip)
    weights = np.abs(clipped) / z_clip
    disease_norm = np.sqrt(np.sum(weights * clipped ** 2))
    signature_norm = np.sqrt(np.sum(weights[None, :] * matrix ** 2, axis=1))
    denom = disease_norm * signature_norm
    numerator = np.sum(weights[None, :] * matrix * clipped[None, :], axis=1)
    return np.divide(numerator, denom, out=np.zeros_like(numerator), where=denom > 0)


def _group_medians(values, groups):
    return np.array([np.median(values[group]) for group in groups], dtype=float)


def permutation_calibration(matrix, disease_z, groups, permutations, seed, z_clip=5.0, strata=None):
    """Calibrate compound/cell medians against deterministic gene-label permutations."""
    if not isinstance(permutations, int) or permutations < 1:
        raise ValueError("permutations must be a positive integer")
    observed_individual = weighted_cosine_scores(matrix, disease_z, z_clip)
    observed = _group_medians(observed_individual, groups)
    counts = np.zeros(len(groups), dtype=int)
    rng = np.random.default_rng(seed)
    strata = np.asarray(["all"] * len(disease_z) if strata is None else strata)
    if strata.ndim != 1 or len(strata) != len(disease_z):
        raise ValueError("strata must align one-to-one with disease_z")
    stratum_indices = [np.flatnonzero(strata == value) for value in np.unique(strata)]
    for _ in range(permutations):
        permuted = np.asarray(disease_z, dtype=float).copy()
        for indices in stratum_indices:
            permuted[indices] = disease_z[indices][rng.permutation(len(indices))]
        null = _group_medians(weighted_cosine_scores(matrix, permuted, z_clip), groups)
        counts += null <= observed + 1e-12
    return observed, (counts + 1) / (permutations + 1), observed_individual


def _write_table(path, frame):
    frame.to_csv(path, sep="\t", index=False)


def _empty_results():
    return pd.DataFrame(columns=RESULT_COLUMNS)


def screen_disease(config, reference, output):
    """Run one disease independently and write tidy primary/sensitivity outputs."""
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    required = {"disease", "twas", "tissue", "primary_cells", "secondary_cells", "z_clip",
                "minimum_genes", "permutations", "seed", "alpha"}
    if required - set(config):
        raise ValueError("Disease screen config missing: " + ", ".join(sorted(required - set(config))))
    if set(config["primary_cells"]) & set(config["secondary_cells"]):
        raise ValueError("Primary and secondary cells must be disjoint")
    if not 0 < float(config["alpha"]) < 1:
        raise ValueError("alpha must lie in (0,1)")
    twas = pd.read_csv(config["twas"], sep="\t")
    needed = {"gene_name", "tissue", "zscore"}
    if needed - set(twas):
        raise ValueError("TWAS input missing: " + ", ".join(sorted(needed - set(twas))))
    twas = twas[twas.tissue == config["tissue"]].copy()
    twas["gene_name"] = twas.gene_name.map(symbol)
    twas["zscore"] = pd.to_numeric(twas.zscore, errors="coerce")
    twas = twas[(twas.gene_name != "") & np.isfinite(twas.zscore)]
    if twas.gene_name.duplicated().any():
        raise ValueError("TWAS input has duplicate genes in the selected tissue")
    twas = twas[twas.gene_name.isin(reference.gene_index)].sort_values("gene_name").reset_index(drop=True)
    twas["zscore_clipped"] = twas.zscore.clip(-float(config["z_clip"]), float(config["z_clip"]))
    twas["weight"] = twas.zscore_clipped.abs() / float(config["z_clip"])
    _write_table(output / "aligned_twas_genes.tsv.gz", twas)

    summary = {
        "disease": config["disease"], "input_kind": config.get("input_kind", "unspecified"),
        "tissue": config["tissue"], "primary_cells": config["primary_cells"],
        "secondary_cells": config["secondary_cells"], "shared_genes": len(twas),
        "minimum_genes": int(config["minimum_genes"]), "z_clip": float(config["z_clip"]),
        "permutations": int(config["permutations"]), "seed": int(config["seed"]),
        "alpha": float(config["alpha"]), "score": "weighted_cosine; negative is reversal",
        "null": "gene-label permutation within CMap landmark/inferred strata, preserving TWAS Z and weight distribution",
    }
    if len(twas) < int(config["minimum_genes"]):
        empty = _empty_results()
        _write_table(output / "all_primary_tests.tsv", empty)
        _write_table(output / "secondary_contexts.tsv", empty)
        _write_table(output / "candidate_compounds.tsv", empty)
        _write_table(output / "signature_evidence.tsv.gz", pd.DataFrame(columns=["sig_id", "score"]))
        _write_table(output / "permutation_summary.tsv", pd.DataFrame(columns=["compound_id", "cell", "empirical_p"]))
        summary.update({"status": "insufficient_gene_overlap", "primary_hypotheses": 0, "candidate_count": 0})
        (output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
        return output

    metadata, matrices = [], []
    for block_meta, block in reference.iter_signature_blocks(twas.gene_name.tolist()):
        metadata.append(block_meta)
        matrices.append(block)
    sig = pd.concat(metadata, ignore_index=True)
    matrix = np.vstack(matrices)
    groups = defaultdict(list)
    for index, row in sig.iterrows():
        groups[(str(row.pert_id), str(row.cell_id))].append(index)
    keys = sorted(groups)
    group_indices = [groups[key] for key in keys]
    observed, pvalues, individual = permutation_calibration(
        matrix, twas.zscore.to_numpy(float), group_indices,
        int(config["permutations"]), int(config["seed"]), float(config["z_clip"]),
        strata=np.array([reference.gene_strata[gene] for gene in twas.gene_name]),
    )
    sig = sig.copy()
    sig["score"] = individual
    sig.insert(0, "disease", config["disease"])
    _write_table(output / "signature_evidence.tsv.gz", sig)

    rows = []
    names = dict(zip(sig.pert_id.astype(str), sig.pert_iname.astype(str)))
    primary_cells = set(config["primary_cells"])
    secondary_cells = set(config["secondary_cells"])
    for (compound, cell), group, score, pvalue in zip(keys, group_indices, observed, pvalues):
        if cell not in primary_cells | secondary_cells:
            continue
        rows.append({
            "disease": config["disease"], "compound_id": compound, "drug": names[compound], "cell": cell,
            "score": float(score), "empirical_p": float(pvalue), "q_bh_primary_family": np.nan,
            "signature_count": len(group), "negative_signature_fraction": float((individual[group] < 0).mean()),
            "primary_context": cell in primary_cells, "candidate": False,
        })
    results = pd.DataFrame(rows, columns=RESULT_COLUMNS)
    primary = results[results.primary_context].copy().reset_index(drop=True)
    if len(primary):
        primary["q_bh_primary_family"] = bh(primary.empirical_p.to_numpy(float))
        primary["candidate"] = (primary.score < 0) & (primary.q_bh_primary_family < float(config["alpha"]))
    secondary = results[~results.primary_context].copy().reset_index(drop=True)
    primary = primary.sort_values(["q_bh_primary_family", "score", "compound_id"]).reset_index(drop=True)
    secondary = secondary.sort_values(["empirical_p", "score", "compound_id"]).reset_index(drop=True)
    candidate_rows = []
    for compound, frame in primary.groupby("compound_id", sort=True):
        supported = frame.candidate
        # A compound appears once: at least one prespecified primary context must pass,
        # and no tested primary context may point in the disease-mimicking direction.
        if supported.any() and (frame.score < 0).all():
            candidate_rows.append({
                "compound_id": compound, "drug": frame.drug.iloc[0],
                "primary_cells_tested": len(frame), "primary_cells_supported": int(supported.sum()),
                "median_primary_score": float(frame.score.median()),
                "best_q_bh_primary_family": float(frame.q_bh_primary_family.min()),
                "worst_primary_score": float(frame.score.max()),
            })
    candidates = pd.DataFrame(candidate_rows, columns=[
        "compound_id", "drug", "primary_cells_tested", "primary_cells_supported",
        "median_primary_score", "best_q_bh_primary_family", "worst_primary_score",
    ]).sort_values(["best_q_bh_primary_family", "median_primary_score", "compound_id"]).reset_index(drop=True)
    candidates.insert(0, "candidate_rank", np.arange(1, len(candidates) + 1))
    _write_table(output / "all_primary_tests.tsv", primary)
    _write_table(output / "secondary_contexts.tsv", secondary)
    _write_table(output / "candidate_compounds.tsv", candidates)
    _write_table(output / "permutation_summary.tsv", primary[["compound_id", "cell", "score", "empirical_p", "q_bh_primary_family"]])
    summary.update({
        "status": "complete", "primary_hypotheses": len(primary), "secondary_hypotheses": len(secondary),
        "candidate_count": len(candidates), "minimum_empirical_p": 1 / (int(config["permutations"]) + 1),
        "isolated_hit_bh_resolution": min(1.0, len(primary) / (int(config["permutations"]) + 1)),
        "fdr_family": "all compound x primary-cell tests within this disease",
    })
    (output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return output
