"""Integration contract for four independent disease correction families."""
import json
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np
import pandas as pd
import yaml

from twaslib.four_disease_screen import run_four_disease_screen


DISEASES = ["AD", "LSC", "atopic_dermatitis", "IPF"]


def batch_fixture(root):
    genes = [str(i) for i in range(1, 7)]
    symbols = [f"G{i}" for i in range(1, 7)]
    disease_z = np.array([5, 4, 3, -3, -4, -5], dtype=float)
    cells = ["NEU", "NPC", "MCH58", "HCC515", "A549"]
    signatures, columns, vectors = [], [], []
    for cell in cells:
        for compound, vector in [("A", -disease_z), ("B", disease_z)]:
            sid = f"{compound}_{cell}"
            columns.append(sid)
            vectors.append(vector)
            signatures.append({"sig_id": sid, "pert_id": compound, "pert_iname": "reverse" if compound == "A" else "mimic",
                               "pert_type": "trt_cp", "cell_id": cell, "pert_dose": 1, "pert_time": 24})
    with h5py.File(root / "level5.gctx", "w") as handle:
        handle["/0/DATA/0/matrix"] = np.asarray(vectors, dtype=np.float32).T
        handle["/0/META/ROW/id"] = np.array(genes, dtype="S")
        handle["/0/META/COL/id"] = np.array(columns, dtype="S")
    pd.DataFrame(signatures).to_csv(root / "sig.tsv", sep="\t", index=False)
    pd.DataFrame({"pr_gene_id": genes, "pr_gene_symbol": symbols, "pr_is_lm": 1}).to_csv(root / "gene.tsv", sep="\t", index=False)
    pd.DataFrame({"pert_id": ["A", "B"], "pert_iname": ["reverse", "mimic"]}).to_csv(root / "pert.tsv", sep="\t", index=False)
    disease_specs = []
    contexts = {"AD": (["NEU", "NPC"], ["A549"]), "LSC": (["MCH58"], ["A549"]),
                "atopic_dermatitis": (["MCH58"], ["A549"]), "IPF": (["HCC515"], ["A549"])}
    for disease in DISEASES:
        twas = root / f"{disease}.tsv"
        pd.DataFrame({"gene_name": symbols, "tissue": "Tissue", "zscore": disease_z, "pvalue": .01}).to_csv(twas, sep="\t", index=False)
        primary, secondary = contexts[disease]
        disease_specs.append({"disease": disease, "twas": twas.name, "tissue": "Tissue",
                              "primary_cells": primary, "secondary_cells": secondary})
    config = {
        "schema_version": 1, "input_kind": "synthetic_test", "output": "batch",
        "analysis": {"z_clip": 5, "minimum_genes": 6, "permutations": 199, "seed": 19, "alpha": .05},
        "cmap": {"gctx": "level5.gctx", "sig_info": "sig.tsv", "gene_info": "gene.tsv", "pert_info": "pert.tsv",
                 "cells": cells},
        "diseases": disease_specs,
    }
    path = root / "screen.yaml"
    path.write_text(yaml.safe_dump(config, sort_keys=False))
    return path


class FourDiseaseScreenTests(unittest.TestCase):
    def test_four_diseases_have_isolated_fdr_and_auditable_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = batch_fixture(root)
            out = run_four_disease_screen(config)

            index = json.loads((out / "summary.json").read_text())
            self.assertEqual(index["diseases"], DISEASES)
            self.assertEqual(index["correction_families"], 4)
            expected_primary_counts = {"AD": 2, "LSC": 1, "atopic_dermatitis": 1, "IPF": 1}
            for disease in DISEASES:
                disease_out = out / disease
                primary = pd.read_csv(disease_out / "all_primary_tests.tsv", sep="\t")
                candidates = pd.read_csv(disease_out / "candidate_compounds.tsv", sep="\t")
                self.assertTrue((primary.disease == disease).all())
                self.assertEqual(candidates.compound_id.tolist(), ["A"])
                self.assertEqual(candidates.primary_cells_supported.iloc[0], expected_primary_counts[disease])
                self.assertNotIn("pooled_rank", candidates.columns)
                self.assertTrue((disease_out / "resolved_config.yaml").is_file())
                self.assertTrue((disease_out / "input_manifest.json").is_file())
                self.assertTrue((disease_out / "output_manifest.json").is_file())
                self.assertIn("synthetic", (disease_out / "report_zh.md").read_text())
            cross = pd.read_csv(out / "cross_disease_descriptive.tsv", sep="\t")
            self.assertEqual(set(cross.disease), set(DISEASES))
            self.assertNotIn("pooled_rank", cross.columns)
            input_manifest = json.loads((out / "input_manifest.json").read_text())
            self.assertTrue(all(len(record["sha256"]) == 64 for record in input_manifest))
            self.assertTrue(any(record["role"].startswith("implementation:") for record in input_manifest))
            runtime = json.loads((out / "runtime_versions.json").read_text())
            self.assertIn("python", runtime)
            self.assertIn("numpy", runtime)
            self.assertIn("h5py", runtime)
            output_manifest = json.loads((out / "output_manifest.json").read_text())
            self.assertTrue(any(record["file"] == "AD/candidate_compounds.tsv" for record in output_manifest))

    def test_existing_batch_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = batch_fixture(Path(tmp))
            run_four_disease_screen(config)
            with self.assertRaises(FileExistsError):
                run_four_disease_screen(config)


if __name__ == "__main__":
    unittest.main()
