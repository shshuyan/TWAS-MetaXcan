"""GCTX boundary tests using the official Level 5 group names."""
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np
import pandas as pd

from twaslib.cmap_level5 import Level5Reference


def write_reference(root, transpose=False):
    genes = ["101", "102", "103"]
    signatures = ["SIG_A", "SIG_B", "SIG_KD", "SIG_KD2"]
    values = np.array([[1, 10, 100, 1000], [2, 20, 200, 2000], [3, 30, 300, 3000]], dtype=np.float32)
    with h5py.File(root / "level5.gctx", "w") as handle:
        handle["/0/DATA/0/matrix"] = values.T if transpose else values
        handle["/0/META/ROW/id"] = np.array(genes, dtype="S")
        handle["/0/META/COL/id"] = np.array(signatures, dtype="S")
    pd.DataFrame([
        {"sig_id": "SIG_A", "pert_id": "P1", "pert_iname": "drug-a", "pert_type": "trt_cp", "cell_id": "NEU", "pert_dose": 1, "pert_time": 6},
        {"sig_id": "SIG_B", "pert_id": "P2", "pert_iname": "drug-b", "pert_type": "trt_cp", "cell_id": "A549", "pert_dose": 2, "pert_time": 24},
        {"sig_id": "SIG_KD", "pert_id": "KD", "pert_iname": "knockdown", "pert_type": "trt_sh", "cell_id": "NEU", "pert_dose": 0, "pert_time": 24},
        {"sig_id": "SIG_KD2", "pert_id": "KD2", "pert_iname": "knockdown-2", "pert_type": "trt_sh", "cell_id": "NEU", "pert_dose": 0, "pert_time": 24},
    ]).to_csv(root / "sig.tsv", sep="\t", index=False)
    pd.DataFrame({"pr_gene_id": genes, "pr_gene_symbol": ["G1", "G2", "G3"], "pr_is_lm": [1, 1, 0]}).to_csv(root / "gene.tsv", sep="\t", index=False)
    pd.DataFrame({"pert_id": ["P1", "P2"], "pert_iname": ["drug-a", "drug-b"]}).to_csv(root / "pert.tsv", sep="\t", index=False)
    return root / "level5.gctx", root / "sig.tsv", root / "gene.tsv", root / "pert.tsv"


class Level5ReferenceTests(unittest.TestCase):
    def test_detects_both_matrix_orientations_and_filters_compounds(self):
        for transpose in [False, True]:
            with self.subTest(transpose=transpose), tempfile.TemporaryDirectory() as tmp:
                paths = write_reference(Path(tmp), transpose)
                reference = Level5Reference(*paths, cells=["NEU", "A549"])
                blocks = list(reference.iter_signature_blocks(["G3", "G1"], block_size=1))
                self.assertEqual(reference.genes, ["G1", "G2", "G3"])
                self.assertEqual(reference.signatures.sig_id.tolist(), ["SIG_A", "SIG_B"])
                self.assertEqual([b[0].sig_id.tolist() for b in blocks], [["SIG_A"], ["SIG_B"]])
                self.assertTrue(np.array_equal(np.vstack([b[1] for b in blocks]), [[3, 1], [30, 10]]))

    def test_unknown_cell_and_duplicate_gene_symbol_fail_explicitly(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = write_reference(root)
            with self.assertRaisesRegex(ValueError, "no compound signatures"):
                Level5Reference(*paths, cells=["MISSING"])
            genes = pd.read_csv(root / "gene.tsv", sep="\t")
            genes.loc[1, "pr_gene_symbol"] = "G1"
            genes.to_csv(root / "gene.tsv", sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "duplicate gene symbol"):
                Level5Reference(*paths, cells=["NEU"])


if __name__ == "__main__":
    unittest.main()
