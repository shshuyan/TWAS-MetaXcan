"""Verify a completed real/synthetic output without rerunning Monte Carlo.

python tests/verify_recheck_output.py /path/to/completed-output
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.stats import fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from twaslib.recheck_stats import bh


def verify(out):
    s = json.loads((out / "summary.json").read_text())
    assert s["status"] == "complete"
    manifest = json.loads((out / "input_manifest.json").read_text())
    for r in manifest:
        p = Path(r["path"]); h = hashlib.sha256()
        with p.open("rb") as f:
            for b in iter(lambda: f.read(1024 * 1024), b""):
                h.update(b)
        assert h.hexdigest() == r["sha256"], str(p)
    genes = json.loads((out / "genetic_enrichment.json").read_text())
    tested = [r for r in genes if r["status"] == "tested"]
    for r in tested:
        a = np.array(r["table"])
        assert (a >= 0).all() and a.sum() == r["background_genes"]
        assert a[0, 0] == len(r["overlap_genes"])
        assert np.isclose(fisher_exact(a, alternative="greater").pvalue, r["pvalue"])
    assert np.allclose(bh([r["pvalue"] for r in tested]), [r["q_bh"] for r in tested])
    reversal = json.loads((out / "reversal_by_cohort.json").read_text())
    signatures = json.loads((out / "signature_evidence.json").read_text())
    for r in reversal:
        if r["status"] not in {"supported", "not_supported"}:
            continue
        assert r["status"] == ("supported" if r["score"] < 0 and r["q_bh"] < s["methods"]["alpha"] else "not_supported")
        raw = [x["score"] for x in signatures if x["candidate_id"] == r["candidate_id"] and x["cohort"] == r["cohort"]
               and x["source"] == r["source"] and (r["source"] == "cmap" or x["cell"] == r["cell"])]
        assert len(raw) == r["signatures"]
        assert np.isclose(np.median(raw), r["score"])
        assert np.isclose(r["empirical_p"] * (s["methods"]["permutations"] + 1),
                          round(r["empirical_p"] * (s["methods"]["permutations"] + 1)))
    for cohort in [c["id"] for c in s["methods"]["cohorts"]]:
        lin = [r for r in reversal if r["cohort"] == cohort and r["source"] == "lincs" and "empirical_p" in r]
        assert np.allclose(bh([r["empirical_p"] for r in lin]), [r["q_bh"] for r in lin])
        p = out / (cohort + "_all_cmap_compounds.tsv")
        if p.exists():
            import pandas as pd
            d = pd.read_csv(p, sep="\t")
            assert np.allclose(bh(d.empirical_p), d.q_bh)
    candidates = json.loads((out / "candidate_recheck.json").read_text())
    assert len(candidates) == s["candidate_count"]
    assert all(not r["wet_lab_ready"] for r in candidates)
    experiment_count = sum(1 for _ in (out / "virtual_experiments.jsonl").open())
    assert experiment_count == len(genes) + len(reversal)
    result = {"status": "passed", "input_hashes_verified": len(manifest), "candidate_count": len(candidates),
              "genetic_tests": len(tested), "reversal_records": len(reversal),
              "signature_records": len(signatures), "virtual_experiment_records": experiment_count,
              "verification_scope": "hashes, Fisher, BH, score medians, empirical-p granularity, decisions and trace records; Monte Carlo not repeated"}
    (out / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("output", type=Path)
    print(json.dumps(verify(parser.parse_args().output), indent=2))
