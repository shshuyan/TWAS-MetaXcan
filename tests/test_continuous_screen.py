"""Hand-derived statistical and disease-screen contracts."""
import json
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np
import pandas as pd

from twaslib.cmap_level5 import Level5Reference
from twaslib.continuous_screen import (matched_signature_calibration, permutation_calibration,
                                       screen_disease, weighted_cosine_scores)


def screening_fixture(root):
    genes = [str(i) for i in range(1, 7)]
    symbols = [f"G{i}" for i in range(1, 7)]
    sigs = ["A1", "A2", "B1", "S1"]
    disease = np.array([5, 4, 3, -3, -4, -5], dtype=float)
    matrix = np.array([-disease, -disease, disease, -disease], dtype=np.float32).T
    with h5py.File(root / "level5.gctx", "w") as handle:
        handle["/0/DATA/0/matrix"] = matrix
        handle["/0/META/ROW/id"] = np.array(genes, dtype="S")
        handle["/0/META/COL/id"] = np.array(sigs, dtype="S")
    pd.DataFrame([
        {"sig_id": "A1", "pert_id": "A", "pert_iname": "reverse", "pert_type": "trt_cp", "cell_id": "NEU", "pert_dose": 1, "pert_time": 6},
        {"sig_id": "A2", "pert_id": "A", "pert_iname": "reverse", "pert_type": "trt_cp", "cell_id": "NEU", "pert_dose": 10, "pert_time": 24},
        {"sig_id": "B1", "pert_id": "B", "pert_iname": "mimic", "pert_type": "trt_cp", "cell_id": "NEU", "pert_dose": 1, "pert_time": 24},
        {"sig_id": "S1", "pert_id": "S", "pert_iname": "secondary-only", "pert_type": "trt_cp", "cell_id": "A549", "pert_dose": 1, "pert_time": 24},
    ]).to_csv(root / "sig.tsv", sep="\t", index=False)
    pd.DataFrame({"pr_gene_id": genes, "pr_gene_symbol": symbols, "pr_is_lm": 1}).to_csv(root / "gene.tsv", sep="\t", index=False)
    pd.DataFrame({"pert_id": ["A", "B", "S"], "pert_iname": ["reverse", "mimic", "secondary-only"]}).to_csv(root / "pert.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_name": symbols, "tissue": "Brain_Cortex", "zscore": disease, "pvalue": .01}).to_csv(root / "twas.tsv", sep="\t", index=False)
    return Level5Reference(root / "level5.gctx", root / "sig.tsv", root / "gene.tsv", root / "pert.tsv", cells=["NEU", "A549"])


class ContinuousStatisticsTests(unittest.TestCase):
    def test_weighted_cosine_has_hand_computed_direction_and_zero(self):
        matrix = np.array([[-1, 1], [1, -1], [1, 1]], dtype=float)
        scores = weighted_cosine_scores(matrix, np.array([5, -5], dtype=float), z_clip=5)
        self.assertTrue(np.allclose(scores, [-1, 1, 0]))

    def test_permutation_calibration_uses_group_median_and_is_seeded(self):
        disease = np.array([5, 4, -4, -5], dtype=float)
        matrix = np.array([-disease, disease, -disease], dtype=float)
        groups = [[0, 2], [1]]
        first = permutation_calibration(matrix, disease, groups, permutations=99, seed=13, z_clip=5)
        second = permutation_calibration(matrix, disease, groups, permutations=99, seed=13, z_clip=5)
        self.assertTrue(np.array_equal(first[0], [-1, 1]))
        self.assertTrue(np.array_equal(first[0], second[0]))
        self.assertTrue(np.array_equal(first[1], second[1]))
        self.assertTrue(np.all((first[1] >= .01) & (first[1] <= 1)))

    def test_permutations_shuffle_only_within_prespecified_gene_strata(self):
        disease = np.array([5, 3, -3, -5], dtype=float)
        matrix = np.array([-disease], dtype=float)
        observed, empirical_p, _ = permutation_calibration(
            matrix, disease, [[0]], permutations=19, seed=3, z_clip=5,
            strata=np.array(["s1", "s2", "s3", "s4"]),
        )
        self.assertTrue(np.allclose(observed, [-1]))
        self.assertEqual(empirical_p.tolist(), [1])

    def test_primary_calibration_matches_signature_context_and_group_size(self):
        scores = np.array([-1, -.9, 1, .9], dtype=float)
        groups = [[0, 1], [2, 3]]
        strata = np.array(["NEU|1|24"] * 4)
        observed, empirical_p = matched_signature_calibration(
            scores, groups, strata, permutations=99, seed=11,
        )
        self.assertTrue(np.allclose(observed, [-.95, .95]))
        self.assertEqual(empirical_p.tolist(), [.01, 1])

    def test_primary_family_fdr_and_secondary_nonpromotion(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reference = screening_fixture(root)
            out = screen_disease({
                "disease": "AD", "input_kind": "synthetic_test", "twas": str(root / "twas.tsv"),
                "tissue": "Brain_Cortex", "primary_cells": ["NEU"], "secondary_cells": ["A549"],
                "z_clip": 5, "minimum_genes": 6, "permutations": 199, "seed": 7, "alpha": .05,
            }, reference, root / "out")
            primary = pd.read_csv(out / "all_primary_tests.tsv", sep="\t")
            secondary = pd.read_csv(out / "secondary_contexts.tsv", sep="\t")
            candidates = pd.read_csv(out / "candidate_compounds.tsv", sep="\t")
            evidence = pd.read_csv(out / "signature_evidence.tsv.gz", sep="\t")
            self.assertEqual(primary.compound_id.tolist(), ["A", "B"])
            self.assertEqual(candidates.compound_id.tolist(), ["A"])
            self.assertEqual(primary.loc[primary.compound_id == "A", "signature_count"].item(), 2)
            self.assertAlmostEqual(primary.loc[primary.compound_id == "A", "score"].item(), -1)
            self.assertEqual(secondary.compound_id.tolist(), ["S"])
            self.assertTrue(secondary.q_bh_primary_family.isna().all())
            self.assertFalse(secondary.candidate.any())
            self.assertEqual(len(evidence), 4)
            summary = json.loads((out / "summary.json").read_text())
            self.assertEqual(summary["primary_hypotheses"], 2)
            self.assertIn("signature-label", summary["null"])

    def test_too_few_shared_genes_is_auditable_not_relaxed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reference = screening_fixture(root)
            twas = pd.read_csv(root / "twas.tsv", sep="\t").iloc[:2]
            twas.to_csv(root / "twas.tsv", sep="\t", index=False)
            out = screen_disease({
                "disease": "AD", "input_kind": "synthetic_test", "twas": str(root / "twas.tsv"),
                "tissue": "Brain_Cortex", "primary_cells": ["NEU"], "secondary_cells": ["A549"],
                "z_clip": 5, "minimum_genes": 6, "permutations": 199, "seed": 7, "alpha": .05,
            }, reference, root / "out")
            summary = json.loads((out / "summary.json").read_text())
            self.assertEqual(summary["status"], "insufficient_gene_overlap")
            self.assertEqual(summary["shared_genes"], 2)
            self.assertEqual(len(pd.read_csv(out / "candidate_compounds.tsv", sep="\t")), 0)


if __name__ == "__main__":
    unittest.main()
