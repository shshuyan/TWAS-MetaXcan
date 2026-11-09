"""TWAS-led discovery: scientific direction, whole-library multiplicity, and CLI tests."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import h5py
import numpy as np
import pandas as pd
import yaml

from twaslib import config, steps
from twaslib.drug_screen import run_drug_screen
from twaslib.drug_screen_inputs import load_twas_queries, tissue_query
from twaslib.screen_stats import bh


def fixture(root):
    genes = [f"G{i:03}" for i in range(200)]
    sp = pd.DataFrame({"gene_name": genes, "tissue": "Brain_Cortex",
                       "pvalue": [1e-8] * 40 + [.8] * 160,
                       "zscore": [6.] * 20 + [-6.] * 20 + [.1] * 160})
    sp.to_csv(root / "sp.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_name": genes, "pvalue": [1e-9] * 40 + [.9] * 160,
                  "z_mean": [-99.] * 200}).to_csv(root / "sm.tsv", sep="\t", index=False)
    pd.DataFrame({"gene": genes}).to_csv(root / "universe.tsv", sep="\t", index=False)
    records = []
    for drug, cell, up, down in [("R", "NEU", genes[20:40], genes[:20]),
                                 ("M", "NEU", genes[:20], genes[20:40]),
                                 ("R", "NPC", genes[:20], genes[20:40]),
                                 ("N", "NPC", genes[20:40], genes[:20])]:
        records.append({"pert_id": drug, "pert_name": "drug_" + drug, "cell_id": cell,
                        "sig_id": drug + "_" + cell, "up": up, "down": down, "dose": "1 uM", "time": "24 h"})
    (root / "lincs.jsonl").write_text("".join(json.dumps(r) + "\n" for r in records))
    reverse = np.r_[np.arange(181, 201), np.arange(1, 21), np.arange(21, 181)]
    mimic = np.r_[np.arange(1, 21), np.arange(181, 201), np.arange(21, 181)]
    with h5py.File(root / "cmap.h5", "w") as f:
        f["assay"] = np.array([reverse, mimic, mimic], dtype=float)
        f["rownames"] = np.array([str(i).encode() for i in range(200)])[None, :]
        f["colnames"] = np.array([b"other_drug__MCF7__1uM__24h", b"mimic__MCF7__1uM__24h", b"other_drug__PC3__1uM__24h"])[None, :]
    pd.DataFrame({"cmap_gene_id": range(200), "symbol": genes}).to_csv(root / "mapping.tsv", sep="\t", index=False)
    cfg = {"disease": "AD", "input_kind": "synthetic_test", "output": "screen", "permutations": 499, "seed": 17,
           "primary": {"source": "lincs", "cell": "NEU"},
           "twas": {"model_type": "eqtl", "complete_tested_results": True, "trait_direction": "higher_is_AD_risk",
                    "spredixcan": "sp.tsv", "smultixcan": "sm.tsv", "tissues": ["Brain_Cortex"],
                    "primary_tissue": "Brain_Cortex", "minimum_genes_per_direction": 10},
           "lincs": {"signatures": "lincs.jsonl", "gene_universe": "universe.tsv", "cells": ["NEU", "NPC"], "compound_only": True},
           "cmap": {"matrix": "cmap.h5", "gene_mapping": "mapping.tsv", "rank_direction": "lower_is_up"}}
    write_cfg(root, cfg)
    return cfg


def write_cfg(root, cfg):
    (root / "screen.yaml").write_text(yaml.safe_dump(cfg))


class DrugScreenTests(unittest.TestCase):
    def test_whole_library_discovery_signs_cells_and_global_correction(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fixture(root)
            out = run_drug_screen(root / "screen.yaml")
            rows = json.loads((out / "drug_screen.json").read_text())
            self.assertEqual(len(rows), 7)  # all four LINCS groups + all three CMap groups
            by_key = {(r["source"], r["compound_id"], r["cell"]): r for r in rows}
            self.assertTrue(by_key["lincs", "R", "NEU"]["candidate"])
            self.assertFalse(by_key["lincs", "M", "NEU"]["candidate"])
            self.assertGreater(by_key["lincs", "R", "NPC"]["score"], 0)
            self.assertEqual(by_key["lincs", "N", "NPC"]["status"], "supported")
            self.assertFalse(by_key["lincs", "N", "NPC"]["candidate"])  # no promotion from a secondary cell
            self.assertLess(by_key["cmap", "other_drug", "MCF7"]["score"], 0)
            self.assertGreater(by_key["cmap", "other_drug", "PC3"]["score"], 0)
            self.assertTrue(np.allclose(bh([r["empirical_p"] for r in rows]), [r["q_bh_global"] for r in rows]))
            self.assertEqual([r["compound_id"] for r in json.loads((out / "candidate_drugs.json").read_text())], ["R"])
            self.assertEqual(len(pd.read_csv(out / "spredixcan_all_annotated.tsv.gz", sep="\t")), 200)
            self.assertTrue(all(not r["wet_lab_ready"] for r in rows))
            summary = json.loads((out / "summary.json").read_text())
            self.assertEqual(summary["minimum_empirical_p"], 1 / 500)
            self.assertEqual(summary["isolated_hit_bh_resolution"], 7 / 500)
            self.assertFalse(summary["isolated_hit_resolution_limited"])
            with self.assertRaises(FileExistsError): run_drug_screen(root / "screen.yaml")

    def test_direction_comes_from_single_tissue_not_joint_zmean(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); cfg = fixture(root)
            chosen, _, _, meta = load_twas_queries(cfg["twas"], root / "sp.tsv", root / "sm.tsv")
            q = tissue_query(chosen, set(chosen.gene_name), cfg["twas"])
            self.assertEqual(len(q["risk_positive"]), 20); self.assertEqual(len(q["risk_negative"]), 20)
            self.assertEqual(meta["joint_tests_corrected"], 200)
            sp = pd.read_csv(root / "sp.tsv", sep="\t"); sp.zscore *= -1
            sp.to_csv(root / "sp.tsv", sep="\t", index=False)
            out = run_drug_screen(root / "screen.yaml")
            self.assertEqual([r["compound_id"] for r in json.loads((out / "candidate_drugs.json").read_text())], ["M"])

    def test_missing_single_tissue_direction_never_falls_back_to_joint_mean(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); cfg = fixture(root)
            pd.read_csv(root / "sp.tsv", sep="\t").drop(columns="zscore").to_csv(root / "sp.tsv", sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "zscore"):
                run_drug_screen(root / "screen.yaml")

    def test_conflicts_are_retained_and_default_excluded(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); cfg = fixture(root)
            sp = pd.read_csv(root / "sp.tsv", sep="\t")
            second = sp.assign(tissue="Brain_Hippocampus").copy(); second.loc[0, "zscore"] = -6
            pd.concat([sp, second]).to_csv(root / "sp.tsv", sep="\t", index=False)
            cfg["twas"]["tissues"].append("Brain_Hippocampus")
            chosen, _, _, meta = load_twas_queries(cfg["twas"], root / "sp.tsv", root / "sm.tsv")
            self.assertEqual(meta["conflicting_genes"], ["G000"])
            self.assertFalse(chosen.loc[chosen.gene_name == "G000", "query_eligible"].any())
            self.assertEqual(meta["sp_tests_corrected"], 400)
            cfg["twas"]["exclude_direction_conflicts"] = False
            chosen, _, _, _ = load_twas_queries(cfg["twas"], root / "sp.tsv", root / "sm.tsv")
            self.assertTrue(chosen.loc[chosen.gene_name == "G000", "query_eligible"].all())

    def test_empty_signature_does_not_relax_thresholds(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fixture(root)
            sp = pd.read_csv(root / "sp.tsv", sep="\t"); sp.pvalue = .8
            sp.to_csv(root / "sp.tsv", sep="\t", index=False)
            out = run_drug_screen(root / "screen.yaml")
            summary = json.loads((out / "summary.json").read_text())
            self.assertEqual(summary["tested_groups"], 0)
            self.assertEqual(summary["query_status"], "insufficient_query")
            self.assertEqual(json.loads((out / "candidate_drugs.json").read_text()), [])
            self.assertEqual(len(pd.read_csv(out / "candidate_drugs.tsv", sep="\t")), 0)
            self.assertTrue(all(r["status"] == "insufficient_query" for r in json.loads((out / "drug_screen.json").read_text())))

    def test_background_and_significance_are_before_assay_intersection(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); cfg = fixture(root)
            cfg["twas"]["max_genes_per_direction"] = 10
            chosen, _, _, meta = load_twas_queries(cfg["twas"], root / "sp.tsv", root / "sm.tsv")
            q = tissue_query(chosen, set(chosen.gene_name), cfg["twas"])
            self.assertEqual(len(q["risk_positive"]), 10)
            self.assertEqual(q["eligible_positive_before_limit"], 20)
            self.assertEqual(len(q["background"]), 200)
            self.assertEqual(meta["joint_bonferroni_threshold"], .05 / 200)

    def test_joint_association_gate_is_explicit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); cfg = fixture(root)
            sm = pd.read_csv(root / "sm.tsv", sep="\t"); sm.pvalue = .9
            sm.to_csv(root / "sm.tsv", sep="\t", index=False)
            chosen, _, _, _ = load_twas_queries(cfg["twas"], root / "sp.tsv", root / "sm.tsv")
            self.assertFalse(chosen.query_eligible.any())
            cfg["twas"]["require_joint_significance"] = False
            chosen, _, _, _ = load_twas_queries(cfg["twas"], root / "sp.tsv", root / "sm.tsv")
            self.assertEqual(int(chosen.query_eligible.sum()), 40)

    def test_invalid_orientation_duplicate_tissue_and_splicing_rejected(self):
        for field, value in [("trait_direction", "unknown"), ("model_type", "sqtl"),
                             ("complete_tested_results", False), ("tissues", ["Brain_Cortex", "Brain_Cortex"])]:
            with self.subTest(field=field), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp); cfg = fixture(root); cfg["twas"][field] = value; write_cfg(root, cfg)
                with self.assertRaises(ValueError): run_drug_screen(root / "screen.yaml")

    def test_pipeline_resume_rejects_input_and_output_changes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fixture(root)
            cfg = config.Config("AD", {"drug_screen": {"enabled": True, "config": str(root / "screen.yaml")}},
                                {}, {}, str(root / "work"), str(root / "results"))
            steps.step_drug_screen(cfg); steps.step_drug_screen(cfg)
            (root / "screen" / "candidate_drugs.tsv").write_text("changed\n")
            with self.assertRaisesRegex(SystemExit, "NEW output"):
                steps.step_drug_screen(cfg)

    def test_conditions_are_aggregated_by_median_not_best_score(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); cfg = fixture(root)
            records = [json.loads(line) for line in (root / "lincs.jsonl").read_text().splitlines()]
            reverse = records[0]
            records.append({**reverse, "sig_id": "R_NEU_opposite", "up": reverse["down"], "down": reverse["up"]})
            (root / "lincs.jsonl").write_text("".join(json.dumps(r) + "\n" for r in records))
            out = run_drug_screen(root / "screen.yaml")
            rows = json.loads((out / "drug_screen.json").read_text())
            record = next(r for r in rows if r["source"] == "lincs" and r["compound_id"] == "R" and r["cell"] == "NEU")
            self.assertEqual(record["score"], 0)
            self.assertEqual(record["signatures"], 2)
            self.assertFalse(record["candidate"])

    def test_duplicate_input_and_stale_source_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); cfg = fixture(root)
            sp = pd.read_csv(root / "sp.tsv", sep="\t")
            pd.concat([sp, sp.iloc[:1]]).to_csv(root / "sp.tsv", sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "duplicate"):
                run_drug_screen(root / "screen.yaml")
            sp.to_csv(root / "sp.tsv", sep="\t", index=False)
            pipeline = config.Config("AD", {"drug_screen": {"enabled": True, "config": str(root / "screen.yaml")}},
                                     {}, {}, str(root / "work"), str(root / "results"))
            steps.step_drug_screen(pipeline)
            (root / "sm.tsv").write_text((root / "sm.tsv").read_text() + "\n")
            with self.assertRaisesRegex(SystemExit, "NEW output"):
                steps.step_drug_screen(pipeline)

    def test_cli_runs_without_geo_or_candidate_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); cfg = fixture(root)
            del cfg["lincs"]; cfg["primary"] = {"source": "cmap", "cell": "MCF7"}
            write_cfg(root, cfg)
            result = subprocess.run([sys.executable, str(Path(__file__).resolve().parents[1] / "twas.py"),
                                     "drug-screen", str(root / "screen.yaml")], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            candidates = json.loads((root / "screen" / "candidate_drugs.json").read_text())
            self.assertEqual([r["compound_id"] for r in candidates], ["other_drug"])


if __name__ == "__main__":
    unittest.main()
