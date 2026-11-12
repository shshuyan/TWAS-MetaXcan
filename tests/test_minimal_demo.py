"""The published demo command proves the complete software loop."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pandas as pd


class MinimalDemoTests(unittest.TestCase):
    def test_demo_command_emits_deterministic_four_disease_results_and_figures(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            script = Path(__file__).resolve().parents[1] / "scripts" / "run_minimal_four_disease_demo.py"
            outputs = []
            for name in ["first", "second"]:
                target = root / name
                result = subprocess.run([sys.executable, str(script), "--output", str(target), "--formats", "png", "svg"],
                                        capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                batch = target / "results"
                outputs.append(batch)
                summary = json.loads((batch / "summary.json").read_text())
                self.assertEqual(summary["input_kind"], "synthetic_test")
                self.assertEqual(summary["diseases"], ["AD", "LSC", "atopic_dermatitis", "IPF"])
                for disease in summary["diseases"]:
                    candidates = pd.read_csv(batch / disease / "candidate_compounds.tsv", sep="\t")
                    self.assertGreaterEqual(len(candidates), 1)
                    self.assertIn("DEMO / synthetic", (batch / disease / "report_zh.md").read_text())
                figures = list((batch / "figures").glob("figure_*.png"))
                self.assertEqual(len(figures), 8)
                manifest = json.loads((batch / "output_manifest.json").read_text())
                self.assertTrue(any(item["file"] == "figures/figure_01_pipeline.png" for item in manifest))
            for disease in ["AD", "LSC", "atopic_dermatitis", "IPF"]:
                first = (outputs[0] / disease / "all_primary_tests.tsv").read_bytes()
                second = (outputs[1] / disease / "all_primary_tests.tsv").read_bytes()
                self.assertEqual(first, second)

    def test_demo_refuses_to_overwrite_existing_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target = root / "existing"
            target.mkdir()
            marker = target / "keep.txt"
            marker.write_text("user data")
            script = Path(__file__).resolve().parents[1] / "scripts" / "run_minimal_four_disease_demo.py"
            result = subprocess.run([sys.executable, str(script), "--output", str(target)], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(marker.read_text(), "user data")


if __name__ == "__main__":
    unittest.main()
