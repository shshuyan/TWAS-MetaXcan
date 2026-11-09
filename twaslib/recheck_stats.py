"""Statistics for exploratory AD candidate checks; no causal/efficacy inference."""
import numpy as np
from scipy.stats import fisher_exact


def bh(pvalues):
    p = np.asarray(pvalues, dtype=float)
    if not np.isfinite(p).all() or ((p < 0) | (p > 1)).any():
        raise ValueError("BH requires finite p-values in [0,1]")
    if not len(p):
        return p.copy()
    order = np.argsort(p, kind="stable")
    adjusted = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty(len(p)); out[order] = np.minimum(adjusted, 1)
    return out


def enrichment(footprint, hits, universe):
    """One-sided overrepresentation, conditional on declared gene-set sizes."""
    universe = set(universe)
    f, h = set(footprint) & universe, set(hits) & universe
    overlap = f & h
    a, b, c = len(overlap), len(f - h), len(h - f)
    d = len(universe) - a - b - c
    if not universe or not f:
        return {"status": "not_testable", "pvalue": None, "odds_ratio": None,
                "fold_enrichment": None, "overlap_genes": [], "table": [[a, b], [c, d]]}
    if not h or h == universe or f == universe:
        odds, p = None, 1.0
    else:
        odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        odds = float(odds) if np.isfinite(odds) else None
    expected = len(f) * len(h) / len(universe)
    return {"status": "tested", "pvalue": float(p), "odds_ratio": odds,
            "fold_enrichment": len(overlap) / expected if expected else None,
            "overlap_genes": sorted(overlap), "table": [[a, b], [c, d]],
            "background_genes": len(universe), "genetic_hit_genes": len(h),
            "footprint_genes": len(f), "expected_overlap": expected}


def ks_enrichment(query_ranks, universe_size):
    """Signed unweighted KS: reference rank 1 denotes strongest upregulation."""
    x = np.sort(np.asarray(query_ranks), axis=1) / universe_size
    n = x.shape[1]
    if not n:
        raise ValueError("Empty KS query")
    above = (np.arange(1, n + 1) / n - x).max(axis=1)
    below = (x - np.arange(n) / n).max(axis=1)
    return np.where(above > below, above, -below)


def cmap_connectivity(ranks, up, down):
    eu = ks_enrichment(ranks[:, up], ranks.shape[1])
    ed = ks_enrichment(ranks[:, down], ranks.shape[1])
    return np.where(np.sign(eu) != np.sign(ed), (eu - ed) / 2, 0)


def group_medians(values, groups):
    if isinstance(groups, tuple):
        ix, valid = groups
        return np.nanmedian(np.where(valid, values[ix], np.nan), axis=1)
    return np.array([np.median(values[g]) for g in groups])


def median_index(groups):
    width = max(map(len, groups))
    ix = np.zeros((len(groups), width), dtype=int)
    valid = np.zeros_like(ix, dtype=bool)
    for i, group in enumerate(groups):
        ix[i, :len(group)] = group
        valid[i, :len(group)] = True
    return ix, valid


def cmap_calibration(ranks, up, down, background, groups, permutations, seed, progress=None):
    groups = median_index(groups)
    observed = group_medians(cmap_connectivity(ranks, up, down), groups)
    rng = np.random.default_rng(seed)
    count = np.zeros(len(observed), dtype=int)
    for i in range(permutations):
        ix = rng.choice(background, len(up) + len(down), replace=False)
        null = group_medians(cmap_connectivity(ranks, ix[:len(up)], ix[len(up):]), groups)
        count += null <= observed + 1e-12
        if progress and (i + 1) % 1000 == 0:
            progress(i + 1)
    p = (count + 1) / (permutations + 1)
    return observed, p, bh(p)


def lincs_calibration(signs, up, down, groups, permutations, seed):
    """Tail-overlap sign score, calibrated at the candidate/cell median level."""
    groups = median_index(groups)
    n = len(up) + len(down)
    def score(u, d):
        return (signs[:, u].sum(axis=1) - signs[:, d].sum(axis=1)) / n
    individual = score(up, down)
    observed = group_medians(individual, groups)
    rng = np.random.default_rng(seed)
    count = np.zeros(len(observed), dtype=int)
    for _ in range(permutations):
        ix = rng.choice(signs.shape[1], n, replace=False)
        null = group_medians(score(ix[:len(up)], ix[len(up):]), groups)
        count += null <= observed + 1e-12
    p = (count + 1) / (permutations + 1)
    return individual, observed, p, bh(p)
