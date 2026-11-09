"""Scientific edge cases + CLI/pipeline integration; synthetic inputs only here."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import h5py
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
import yaml

from twaslib.candidate_recheck import evidence_decision, run_recheck
from twaslib.recheck_inputs import disease_query, load_candidates, load_de, load_genetics
from twaslib.recheck_stats import bh, cmap_connectivity, enrichment, lincs_calibration, cmap_calibration
from twaslib import config, steps, progress


def fixture(root):
    genes = [f"G{i}" for i in range(100)]
    c = [{"pert_id": "A", "pert_name": "reverse"}, {"pert_id": "B", "pert_name": "mimic"},
         {"pert_id": "C", "pert_name": "uncovered"}]
    (root / "candidates.json").write_text(json.dumps(c))
    pd.DataFrame({"gene": genes}).to_csv(root / "assay.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_name": genes, "pvalue": [1e-8] * 10 + [.8] * 90}).to_csv(root / "twas.tsv", sep="\t", index=False)
    records = []
    for cid, u, d in [("A", genes[20:40], genes[:20]), ("B", genes[:20], genes[20:40])]:
        records.append({"pert_id": cid, "pert_name": cid, "sig_id": cid + "_NEU", "cell_id": "NEU", "up": u, "down": d})
    (root / "lincs.jsonl").write_text("\n".join(json.dumps(r) for r in records))
    for k in [1, 2]:
        pd.DataFrame({"gene": genes, "log_fc": [1] * 20 + [-1] * 20 + [0] * 60,
                      "p_value": [1e-5 * k] * 40 + [.7] * 60,
                      "fdr_bh": [.001] * 40 + [.8] * 60}).to_csv(root / f"de{k}.tsv", sep="\t", index=False)
    # Rank 1 means drug-up. A reverses disease; B mimics it.
    reverse = np.r_[np.arange(81, 101), np.arange(1, 21), np.arange(21, 81)]
    mimic = np.r_[np.arange(1, 21), np.arange(81, 101), np.arange(21, 81)]
    with h5py.File(root / "cmap.h5", "w") as f:
        f["assay"] = np.array([reverse, mimic], dtype=float)
        f["rownames"] = np.array([str(i).encode() for i in range(100)])[None, :]
        f["colnames"] = np.array([b"reverse__MCF7__trt_cp", b"mimic__PC3__trt_cp"])[None, :]
    pd.DataFrame({"cmap_gene_id": range(100), "symbol": genes}).to_csv(root / "mapping.tsv", sep="\t", index=False)
    cfg = {"disease": "AD", "candidates": "candidates.json", "output": "out", "permutations": 499,
           "query_size": 20, "genetics": {"mode": "twas", "path": "twas.tsv", "complete_tested_results": True},
           "cohorts": [{"id": "D1", "role": "discovery", "differential_expression": "de1.tsv"},
                       {"id": "D2", "role": "validation", "differential_expression": "de2.tsv"}],
           "lincs": {"signatures": "lincs.jsonl", "gene_universe": "assay.tsv", "cells": ["NEU"],
                     "reused_from_upstream_screen": True},
           "cmap": {"matrix": "cmap.h5", "gene_mapping": "mapping.tsv", "rank_direction": "lower_is_up", "tail_size": 20}}
    (root / "evidence.yaml").write_text(yaml.safe_dump(cfg))
    return cfg


class ScientificTests(unittest.TestCase):
    def test_fisher_uses_full_declared_background(self):
        r = enrichment(set(range(10)), set(range(5)), set(range(100)))
        self.assertEqual(r["table"], [[5, 5], [0, 90]])
        self.assertAlmostEqual(r["pvalue"], fisher_exact(r["table"], alternative="greater").pvalue)
        self.assertEqual(r["fold_enrichment"], 10)
        self.assertEqual(enrichment({1}, set(), set(range(100)))["pvalue"], 1)

    def test_ks_rank_orientation_and_zero_rule(self):
        ranks = np.arange(1, 101, dtype=float)[None, :]
        self.assertGreater(cmap_connectivity(ranks, list(range(10)), list(range(90, 100)))[0], .8)
        self.assertLess(cmap_connectivity(ranks, list(range(90, 100)), list(range(10)))[0], -.8)
        self.assertEqual(cmap_connectivity(ranks, list(range(10)), list(range(10, 20)))[0], 0)

    def test_bh_family_and_missing_evidence(self):
        self.assertTrue(np.allclose(bh([.01, .04, .03, 1]), [.04, .053333333, .053333333, 1]))
        self.assertEqual(evidence_decision(True, ["supported", "not_covered"], "twas"), "incomplete_reversal_evidence")
        self.assertEqual(evidence_decision(True, ["supported", "not_supported"], "twas"), "genetics_only_reversal_not_reproduced")

    def test_calibration_preserves_three_or_more_group_dimensions(self):
        signs = np.array([[1, -1, 0, 0], [-1, 1, 0, 0], [0, 0, 1, -1]])
        groups = [[0], [1], [2]]
        _, obs, p, q = lincs_calibration(signs, [0], [1], groups, 100, 1)
        self.assertEqual(len(obs), 3); self.assertEqual(len(p), 3); self.assertEqual(len(q), 3)
        ranks = np.array([[1, 4, 2, 3], [4, 1, 2, 3], [2, 3, 1, 4]])
        obs, p, q = cmap_calibration(ranks, [0], [1], np.arange(4), groups, 100, 1)
        self.assertEqual(len(obs), 3); self.assertEqual(len(p), 3); self.assertEqual(len(q), 3)

    def test_duplicate_twas_genes_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "bad.tsv"
            pd.DataFrame({"gene_name": ["A", "A"], "pvalue": [.01, .02]}).to_csv(p, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "one test per gene"):
                load_genetics({"mode": "twas", "complete_tested_results": True}, p)

    def test_gwas_mapping_is_explicit_and_thresholded(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "gwas.tsv"
            pd.DataFrame({"mapped_genes": ["A;B", "C"], "pvalue": [1e-9, .1]}).to_csv(p, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "study"):
                load_genetics({"mode": "gwas_mapped"}, p)
            hits, background, meta = load_genetics({"mode": "gwas_mapped", "study": "test"}, p)
            self.assertEqual(hits, {"A", "B"}); self.assertIsNone(background)
            self.assertIn("not causal", meta["mapping_strength"])

    def test_query_does_not_relax_fdr(self):
        d = pd.DataFrame({"log_fc": [1, -1], "fdr_bh": [.9, .9], "p_value": [.8, .8]}, index=["A", "B"])
        self.assertEqual(disease_query(d, {"A", "B"})[2], "insufficient_query")

    def test_end_to_end_sign_identity_and_traceability(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fixture(root)
            out = run_recheck(root / "evidence.yaml")
            result = {r["candidate_id"]: r for r in json.loads((out / "candidate_recheck.json").read_text())}
            self.assertEqual(result["A"]["lincs_NEU_D1_status"], "supported")
            self.assertEqual(result["B"]["lincs_NEU_D1_status"], "not_supported")
            self.assertEqual(result["C"]["lincs_NEU_D1_status"], "not_covered")
            self.assertEqual(result["C"]["cmap_legacy_D1_status"], "not_covered")
            self.assertLess(result["A"]["cmap_legacy_D1_score"], 0)
            self.assertGreater(result["B"]["cmap_legacy_D1_score"], 0)
            self.assertFalse(any(r["wet_lab_ready"] for r in result.values()))
            self.assertTrue((out / "input_manifest.json").is_file())
            self.assertEqual(json.loads((out / "summary.json").read_text())["status"], "complete")
            with self.assertRaises(FileExistsError): run_recheck(root / "evidence.yaml")

    def test_same_de_content_is_not_replication(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fixture(root)
            (root / "de2.tsv").write_bytes((root / "de1.tsv").read_bytes())
            with self.assertRaisesRegex(ValueError, "Same DE"):
                run_recheck(root / "evidence.yaml")

    def test_pipeline_hook_resume_and_stale_input_rejection(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fixture(root)
            cfg = config.Config("fixture", {"candidate_recheck": {"enabled": True, "config": str(root / "evidence.yaml")}},
                                {}, {}, str(root / "work"), str(root / "results"))
            steps.step_candidate_recheck(cfg)
            steps.step_candidate_recheck(cfg)  # identical inputs can resume
            (root / "de1.tsv").write_text((root / "de1.tsv").read_text() + "\n")
            with self.assertRaisesRegex(SystemExit, "changed"):
                steps.step_candidate_recheck(cfg)

    def test_cli_is_callable(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fixture(root)
            r = subprocess.run([sys.executable, str(Path(__file__).resolve().parents[1] / "twas.py"), "recheck", str(root / "evidence.yaml")],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("Completed AD candidate recheck", r.stdout)


if __name__ == "__main__":
    unittest.main()
